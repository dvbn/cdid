# SPDX-License-Identifier: GPL-3.0-only
#
# CP08 common-boundary performance: reusable public input preparation through
# the panel/frame/CP04-reader/occupancy boundary and exact native cp04_payload
# accounting. The request layer forms the canonical row permutation, sorted unit
# order and outcome counts once, validates direct plain numeric-column formulas
# on raw values, and (on the exact first eligibility) seals an invocation-local
# record bound to its own arrays; the pipeline binds it by object identity and
# the prepared constructors assemble the unchanged sealed representations.
#
# UNEXECUTED at authoring time (no R in the authoring job). Every expected value
# comes from independent references written here (the frozen recursive R payload
# walk, brute-force pair counts over raw rows, closed-form bounds of the
# unchanged preflights) or from the unchanged standalone/general route of the
# same installed build, forced in disclosed test-local namespace clones. The
# prepared route never supplies its own truth.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

# ---------------------------------------------------------------------------
# Independent references and fixtures
# ---------------------------------------------------------------------------

# Frozen old R payload walk (R/chained_cp04_corollary1.R before the native kernel).
pc_ref_payload <- function(x, scan_attrs = TRUE) {
  if (is.environment(x) || is.null(x)) return(c(entries = 0, char_bytes = 0))
  ans <- c(entries = length(x) + 8, char_bytes = 0)
  if (is.list(x)) {
    for (v in x) ans <- ans + pc_ref_payload(v)
  } else if (is.character(x)) {
    bytes <- nchar(x, type = "bytes", allowNA = TRUE)
    bytes[is.na(bytes)] <- 0
    ans <- ans + c(entries = sum(ceiling(bytes / 2)), char_bytes = sum(bytes))
  }
  if (scan_attrs && !is.null(attributes(x))) ans <- ans + pc_ref_payload(attributes(x), FALSE)
  ans
}

pc_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 2e7, work = 1e10, output = 1e7) api$ipw_owner_control(
  max_chunk_rows = chunk, max_retained_entries = retained, max_matrix_entries = matrix,
  max_work_ops = work, max_output_entries = output)

# Benchmark-shaped complete balanced panel in public coding (G = 0 never
# treated): one treated cohort 3, unit covariates x (integer-valued double) and
# z, optional non-contiguous ids, row shuffle, second treated cohort, a
# singleton unit (one observed outcome) or a unit weight column.
pc_panel <- function(N = 40L, T = 6L, seed = 17L, ids = seq_len(N), rows_perm = FALSE, second = FALSE,
                     singleton = FALSE, integer_ids = FALSE) {
  set.seed(seed)
  g <- rep(0, N); g[seq(1L, N, by = 2L)] <- 3
  if (second) g[seq(2L, N, by = 4L)] <- 4
  d <- data.frame(id = rep(ids, each = T), period = rep(seq_len(T), N))
  u <- rep(seq_len(N), each = T)
  d$G <- g[u]
  x <- stats::rbinom(N, 3, 0.5); z <- round(stats::rnorm(N), 3)
  d$x <- as.double(x[u]); d$z <- z[u]
  d$f <- factor(c("a", "b")[1 + (u %% 2L)])
  d$w <- 1 + (u %% 3L) / 4
  d$one <- 1
  d$y <- 0.3 * u + 0.1 * d$period + 0.2 * d$x * d$period + ifelse(d$G > 0 & d$period >= d$G, 1.5, 0) +
    0.35 * sin(0.7 * u * d$period) + stats::rnorm(nrow(d))
  if (singleton) d$y[d$id == ids[N] & d$period > 1L] <- NA_real_
  if (integer_ids) d$id <- as.integer(d$id)
  if (rows_perm) d <- d[sample.int(nrow(d)), ]
  rownames(d) <- NULL
  d
}

pc_resolved <- function(chunk = 7L) api$chained_public_resolve_control(if (is.null(chunk)) NULL else list(chunk_rows = chunk))

pc_prepare <- function(d, est = "ipw", xformla = NULL, frame = NULL, weightsname = NULL, clustervars = NULL,
                       resolved = pc_resolved(), ns = api) {
  ns$chained_request_prepare_inputs(data = d, frame = frame, yname = "y", tname = "period", idname = "id", gname = "G",
    xformla = xformla, weightsname = weightsname, clustervars = clustervars, anticipation = 0,
    outcome_regression = identical(est, "reg"), doubly_robust = identical(est, "dr"), or_resource = resolved$profile)
}

# Sealed owners from a prepared input: the standalone constructors and the
# prepared-intake constructors on the very same arrays.
pc_owners_original <- function(prep, ctrl, ns = api) {
  po <- ns$construct_ipw_panel_owner(prep$data_internal, "id", "period", "y", "G", ctrl)
  fo <- ns$construct_ipw_frame_owner(po, prep$frame_ids, prep$frame_Z, weights = prep$frame_w, control = ctrl,
    frame_g = prep$frame_g_arg, allow_frame_only = FALSE, design_on_A_only = !is.null(prep$cp04))
  rd <- if (is.null(prep$cp04)) NULL else ns$construct_cp04_frame_reader(fo, prep$frame_ids, prep$cp04$Z,
    prep$cp04$frame_g, prep$frame_w, ctrl, U = prep$cp04$U, rho = prep$cp04$rho, universe = prep$cp04$universe,
    justification = prep$cp04$justification, allow_unknown_off_A = TRUE)
  list(po = po, fo = fo, rd = rd)
}
pc_owners_prepared <- function(prep, ctrl, ns = api) {
  po <- ns$construct_ipw_panel_owner_prepared(prep$data_internal, "id", "period", "y", "G", ctrl, prepared = prep$prepared)
  fo <- ns$construct_ipw_frame_owner_prepared(po, prep$frame_ids, prep$frame_Z, prep$frame_w, ctrl, prep$frame_g_arg,
    prepared = prep$prepared, design_on_A_only = !is.null(prep$cp04))
  rd <- if (is.null(prep$cp04)) NULL else ns$construct_cp04_frame_reader_prepared(fo, prep$frame_ids, prep$cp04$Z,
    prep$cp04$frame_g, prep$frame_w, ctrl, prepared = prep$prepared, U = prep$cp04$U, rho = prep$cp04$rho,
    universe = prep$cp04$universe, justification = prep$cp04$justification, allow_unknown_off_A = TRUE)
  list(po = po, fo = fo, rd = rd)
}
# Field comparison of sealed objects: the seal environments differ by
# construction; every other retained field (arrays, counts, postings, keys,
# accounting, work diagnostics, control, preflight) must be identical.
pc_strip <- function(x) {
  if (is.null(x)) return(NULL)
  x <- unclass(x)
  x[setdiff(names(x), c("ref", "panel_ref", "frame_ref"))]
}
pc_expect_same_owners <- function(a, b) {
  expect_identical(pc_strip(a$po), pc_strip(b$po))
  expect_identical(pc_strip(a$fo), pc_strip(b$fo))
  expect_identical(pc_strip(a$rd), pc_strip(b$rd))
  expect_true(api$ipw_owner_seal_ok(b$po, "panel_owner"))
  expect_true(api$ipw_owner_seal_ok(b$fo, "frame_owner"))
  if (!is.null(b$rd)) expect_true(api$ipw_owner_seal_ok(b$rd, "cp04_reader"))
  expect_identical(b$fo$panel_ref, b$po$ref)
  if (!is.null(b$rd)) expect_identical(b$rd$frame_ref, b$fo$ref)
}

# Disclosed test-local clone of the namespace (the dc_counted pattern): every
# internal function re-bound into `e`, optional counters, optional overrides.
pc_clone <- function(wrap = character(), overrides = list()) {
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
  make_wrap <- function(nm, inner) {
    force(nm); force(inner)
    function(...) {
      assign(nm, get(nm, envir = counts, inherits = FALSE) + 1, envir = counts)
      inner(...)
    }
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  for (nm in names(overrides)) {
    f <- overrides[[nm]]
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  e$.counts <- counts
  e
}
pc_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0)

# The reader counts its own payload with a NULL `ref` slot before sealing.
pc_strip_ref_slot <- function(rd) { rd["ref"] <- list(NULL); rd }

pc_route_counters <- c("construct_ipw_panel_owner", "construct_ipw_panel_owner_prepared",
  "construct_ipw_frame_owner", "construct_ipw_frame_owner_prepared",
  "construct_cp04_frame_reader", "construct_cp04_frame_reader_prepared",
  "pipeline_occupancy_accumulate_complete", "chained_request_prepared_bind")

# Public-shaped pipeline run (the dc_run pattern) in namespace `ns`, with or
# without the prepared record of `prep`.
pc_run <- function(prep, resolved = pc_resolved(), ns = api, with_prepared = TRUE, cluster = NULL) {
  if (!is.null(prep$cp04)) prep$cp04$automatic_chunk_rows <- is.null(resolved$chunk_rows)
  dims <- api$chained_public_dimensions(prep, biters = 49L, clustered = !is.null(prep$cluster_full))
  built <- api$chained_public_pipeline_control(dims, resolved, gmm_weighting = "optimal", gmm_covariance = "tikhonov")
  ns$run_chained_ipw_pipeline(
    data = prep$data_internal, idname = "id", tname = "period", yname = "y", gname = "G",
    control = built$control, xformla = prep$xformla, weightsname = NULL, cluster_vec = prep$cluster_full,
    target_cohorts = prep$target_cohorts, anticipation = 0L, control_group = "nevertreated", base_period = "varying",
    weighting = "optimal", reduce_structural_redundancy = TRUE, run_numerical_transfer = FALSE, root_witness = NULL,
    aggregations = character(0L), cband = FALSE, alp = 0.05, biters = 49L, seed = NULL,
    frame_id = prep$frame_ids, frame_Z = prep$frame_Z, frame_weights = prep$frame_w, frame_g = prep$frame_g_arg,
    allow_frame_only = prep$n_frame_only > 0L, na.rm = FALSE, bstrap = FALSE, caller_rng = TRUE,
    gmm_covariance = "tikhonov", pt_assumption = "post", pt_start = NULL, cp04 = prep$cp04,
    prepared = if (with_prepared) prep$prepared else NULL)
}

# Brute-force (cohort, s, t, count) occupancy over raw panel rows: every
# within-unit pair of observed outcomes of a pair-eligible unit.
pc_brute_occupancy <- function(d) {
  d <- d[!is.na(d$y), c("id", "period", "G")]
  d$G[d$G == 0] <- Inf
  rows <- list()
  for (u in unique(d$id)) {
    tt <- sort(d$period[d$id == u])
    if (length(tt) < 2L) next
    g <- d$G[d$id == u][1L]
    for (i in seq_len(length(tt) - 1L)) for (j in seq.int(i + 1L, length(tt)))
      rows[[length(rows) + 1L]] <- c(g, tt[i], tt[j])
  }
  m <- do.call(rbind, rows)
  key <- paste(m[, 1], m[, 2], m[, 3], sep = "|")
  agg <- as.data.frame(table(key), stringsAsFactors = FALSE)
  parts <- do.call(rbind, strsplit(agg$key, "|", fixed = TRUE))
  out <- data.frame(cohort = as.double(parts[, 1]), s = as.double(parts[, 2]), t = as.double(parts[, 3]),
    count = as.double(agg$Freq), stringsAsFactors = FALSE)
  out <- out[order(out$cohort, out$s, out$t), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# 1. Native payload accounting equals the frozen recursive R walk
# ---------------------------------------------------------------------------

test_that("C_gmm_cp08_payload reproduces the frozen recursive cp04_payload on attributes, NA, UTF-8, aliases, functions, formulas and pairlists", {
  v <- c(1.5, 2, NA)
  s <- c("plain", NA_character_, "", "été", "qß", strrep("x", 1001))
  shared <- list(v, v, list(v, s), s)
  df <- data.frame(a = 1:4, b = c("u", "v", NA, "ü"), stringsAsFactors = FALSE)
  df_rn <- df; rownames(df_rn) <- c("r1", "r2", "r3", "r4")
  df_compact_explicit <- df; attr(df_compact_explicit, "row.names") <- c(NA_integer_, 4L)
  df_compact_automatic <- df; attr(df_compact_automatic, "row.names") <- c(NA_integer_, -4L)
  f_closure <- function(x, y = 2) { x + y }
  arr1 <- array(1:3, 3L, dimnames = list(c("a", "b", "c")))
  pl <- pairlist(a = 1, 2, b = "three")
  pl_untagged <- pairlist(1, "two", NULL)
  fx <- pc_panel(N = 12L, T = 4L)
  ctrl <- pc_ctrl()
  prep <- pc_prepare(fx, est = "dr", xformla = ~x)
  own <- pc_owners_original(prep, ctrl)
  lv <- api$construct_ipw_link_view(own$fo, 1, 2, g = 3, K = Inf, control = ctrl, link_id = "L1_2",
    functional_key = "cmpcat:pc:1-2")
  posix <- as.POSIXlt(c("2020-01-01", "2020-06-30"), tz = "UTC")
  cases <- list(
    NULL, new.env(), structure(new.env(), class = "classed_env"), 1:7, v, s, c(a = 1, b = 2), TRUE, NA,
    list(), list(NULL), shared, list(x = shared, y = list(NULL, NULL)), df, df_rn, df_compact_explicit, df_compact_automatic, data.frame(),
    matrix(1:6, 2, dimnames = list(c("r1", "r2"), c("c1", "c2", "c3"))), arr1,
    factor(c("lo", "hi", "lo")), expression(x + 1, y), quote(f(x, y = 2)), as.name("sym"),
    ~ x + z, y ~ x, f_closure, sum, `+`, pl, pl_untagged, structure(list(1, "a"), class = "x_pair", note = "attribute"),
    structure(1:3, class = "myint"), posix, list(posix, posix),
    own$po, own$fo, own$rd, lv, attributes(own$rd), own$rd$state, prep$cp04, prep$units,
    structure(list(a = 1), names = c("")), list(a = NULL, b = list(c = NULL)))
  for (k in seq_along(cases)) {
    x <- cases[[k]]
    expect_identical(api$cp04_payload(x), pc_ref_payload(x), info = paste("case", k))
    expect_identical(api$cp04_payload(x, scan_attrs = FALSE), pc_ref_payload(x, FALSE), info = paste("case (no attrs)", k))
  }
  # Return contract: named double vector; every alias occurrence is charged again.
  pay <- api$cp04_payload(shared)
  expect_identical(names(pay), c("entries", "char_bytes"))
  expect_type(pay, "double")
  Pv <- unname(pc_ref_payload(v)["entries"]); Ps <- unname(pc_ref_payload(s)["entries"])
  expect_identical(unname(pay["entries"]), (4 + 8) + Pv + Pv + ((2 + 8) + Pv + Ps) + Ps)
  expect_identical(Pv, 3 + 8)
  # The private replacement preserves the old R condition and length dispatch.
  # In particular numeric flags are not a new stricter helper contract.
  for (flag in list(TRUE, FALSE, 1, 0, NA, NULL, c(TRUE, FALSE))) {
    for (x in list(v, list(a = v))) {
      old <- tryCatch(pc_ref_payload(x, flag), error = identity)
      now <- tryCatch(api$cp04_payload(x, flag), error = identity)
      if (inherits(old, "error")) {
        expect_s3_class(now, "error")
        expect_identical(conditionMessage(now), conditionMessage(old))
      } else expect_identical(now, old)
    }
  }
  had_method <- exists("length.pc_payload", globalenv(), inherits = FALSE)
  old_method <- get0("length.pc_payload", globalenv(), inherits = FALSE)
  withr::defer(if (had_method) assign("length.pc_payload", old_method, globalenv()) else
    rm("length.pc_payload", envir = globalenv()))
  assign("length.pc_payload", function(x) 2L, globalenv())
  classed <- structure(list(1, 2, 3), class = "pc_payload")
  expect_identical(api$cp04_payload(classed), pc_ref_payload(classed))
  expect_identical(api$cp04_payload(classed, FALSE), pc_ref_payload(classed, FALSE))
  # The retained incremental catalog charge still equals the walk on every prefix.
  records <- list(list(pair = list(s = 1, t = 2), link_id = "Lé"), list(), list(w = list(NULL, "qß")))
  st <- list(n = 0L, records = c(entries = 0, char_bytes = 0)); recs <- list()
  for (k in seq_along(records)) {
    recs[[k]] <- records[[k]]
    st <- api$or_baseline_payload_append(st, recs, api$or_baseline_abort)
    expect_identical(st$total, pc_ref_payload(recs))
  }
})

# ---------------------------------------------------------------------------
# 2. Prepared owners equal the standalone owners field for field
# ---------------------------------------------------------------------------

test_that("prepared panel/frame/reader equal the standalone constructors on the eligible benchmark shapes and ID permutations", {
  set.seed(101L)
  shuffled_ids <- sample(1000 + 7 * seq_len(40L))
  shapes <- list(
    list(d = pc_panel(), est = "ipw", xf = NULL),
    list(d = pc_panel(), est = "ipw", xf = ~x),
    list(d = pc_panel(seed = 5L), est = "reg", xf = NULL),
    list(d = pc_panel(seed = 5L), est = "reg", xf = ~x + z),
    list(d = pc_panel(seed = 9L), est = "dr", xf = ~1),
    list(d = pc_panel(seed = 9L), est = "dr", xf = ~x),
    list(d = pc_panel(ids = shuffled_ids, rows_perm = TRUE, seed = 3L), est = "ipw", xf = ~x),
    list(d = pc_panel(ids = shuffled_ids, rows_perm = TRUE, seed = 3L), est = "dr", xf = ~x + z),
    list(d = pc_panel(integer_ids = TRUE, rows_perm = TRUE, seed = 21L), est = "reg", xf = ~x),
    list(d = pc_panel(N = 7L, T = 3L, seed = 2L), est = "ipw", xf = NULL))
  for (k in seq_along(shapes)) {
    sh <- shapes[[k]]
    for (chunk in c(7L, 256L)) {
      prep <- pc_prepare(sh$d, est = sh$est, xformla = sh$xf)
      expect_s3_class(prep$prepared, "cdid_request_prepared_intake")
      expect_true(api$ipw_owner_seal_ok(prep$prepared, "prepared_intake"))
      # The record references, never copies, the arrays it describes.
      expect_true(api$gmm_sexp_aliases(prep$prepared$source, prep$data_internal))
      expect_true(api$gmm_sexp_aliases(prep$prepared$frame_Z, prep$frame_Z))
      expect_true(api$gmm_sexp_aliases(prep$prepared$frame_ids, prep$frame_ids))
      if (is.null(sh$xf)) {
        pe <- environment(prep$xformla)
        # The package default has no retained preparation frame at all.
        expect_identical(pe, baseenv())
        expect_false(exists("prepared", pe, inherits = FALSE))
        expect_false(exists("row_order", pe, inherits = FALSE))
        expect_false(exists("data_internal", pe, inherits = FALSE))
      }
      ctrl <- pc_ctrl(chunk)
      a <- pc_owners_original(prep, ctrl)
      b <- pc_owners_prepared(prep, ctrl)
      pc_expect_same_owners(a, b)
      # The prepared owner's sorted keys ARE the record's vector; the sorted
      # order matches an independent sort of the unit identifiers.
      expect_true(api$gmm_sexp_aliases(b$po$unit_id_key, prep$prepared$unit_key))
      expect_identical(b$po$unit_id_key, sort(as.double(unique(sh$d$id))))
      expect_identical(b$po$units$n_obs, as.integer(table(factor(sh$d$id[!is.na(sh$d$y)], levels = sort(unique(sh$d$id))))))
    }
  }
  # A supplied frame in a permuted order with F = A is eligible and aligned by identity.
  d <- pc_panel(seed = 11L)
  fr <- data.frame(id = rev(unique(d$id)), G = rev(d$G[!duplicated(d$id)]), x = rev(d$x[!duplicated(d$id)]))
  for (est in c("ipw", "dr")) {
    prep <- pc_prepare(d, est = est, xformla = ~x, frame = fr)
    expect_s3_class(prep$prepared, "cdid_request_prepared_intake")
    expect_identical(prep$frame_ids, fr$id)
    ctrl <- pc_ctrl(16L)
    pc_expect_same_owners(pc_owners_original(prep, ctrl), pc_owners_prepared(prep, ctrl))
  }
})

# ---------------------------------------------------------------------------
# 3. Direct-formula validation equals the original NT-row route; frame design unchanged
# ---------------------------------------------------------------------------

test_that("direct plain numeric-column formulas keep the original frame design, errors and dated routing; context formulas keep the NT route", {
  d <- pc_panel(seed = 13L)
  first <- d[!duplicated(d$id), ]
  ordinary <- pc_clone(overrides = list(cp04_design_direct_columns = function(tf, vars, frame, data) FALSE))
  for (est in c("ipw", "reg", "dr")) for (xf in list(NULL, ~1, ~x, ~x + z)) {
    a <- pc_prepare(d, est = est, xformla = xf)
    b <- pc_prepare(d, est = est, xformla = xf, ns = ordinary)
    expect_identical(a$frame_Z, b$frame_Z)
    expect_identical(a$cp04, b$cp04)
    expect_identical(a$units, b$units)
    expect_identical(a$data_internal, b$data_internal)
    expect_identical(a$or_basis, b$or_basis)
    # A default ~1 is created in each distinct namespace-clone call frame.
    # Compare its exact formula AST and all other fields, not those environments.
    aa <- a[setdiff(names(a), "prepared")]; bb <- b[setdiff(names(b), "prepared")]
    if (is.null(xf)) { environment(aa$xformla) <- environment(bb$xformla) <- environment() }
    expect_identical(aa, bb)
    expect_null(b$prepared) # the original NT route is never the prepared route
  }
  # Independent reference of the ~x IPW frame design: the first-observed frame's model matrix.
  a <- pc_prepare(d, est = "ipw", xformla = ~x)
  fz <- stats::model.matrix(~x, data = first); rownames(fz) <- NULL
  attr(fz, "assign") <- NULL
  expect_identical(a$frame_Z, fz)
  # Factors, interactions and row-wise transforms decline the direct route (and the record) and run the NT route.
  for (xf in list(~f, ~x:z, ~I(x^2), ~log(1 + x))) {
    p <- pc_prepare(d, est = "ipw", xformla = xf)
    q <- pc_prepare(d, est = "ipw", xformla = xf, ns = ordinary)
    expect_null(p$prepared)
    expect_identical(p$frame_Z, q$frame_Z)
    expect_identical(p$frame_Z, { m <- stats::model.matrix(xf, data = first); rownames(m) <- NULL; attr(m, "assign") <- NULL; m })
  }
  # Context-dependent calls (scale/poly evaluate differently on N frame rows and
  # NT panel rows) keep the original typed refusal on every route.
  for (xf in list(~scale(z), ~poly(z, 2))) for (est in c("ipw", "reg", "dr")) {
    e1 <- tryCatch(pc_prepare(d, est = est, xformla = xf), error = identity)
    e2 <- tryCatch(pc_prepare(d, est = est, xformla = xf, ns = ordinary), error = identity)
    expect_s3_class(e1, "cdid_chained_request_covariate_error")
    expect_identical(conditionMessage(e1), conditionMessage(e2))
  }
  # Same typed refusals on both routes: NA in a panel covariate, frame/data disagreement.
  dna <- d; dna$x[3L] <- NA_real_
  for (est in c("ipw", "dr")) {
    e1 <- tryCatch(pc_prepare(dna, est = est, xformla = ~x), error = identity)
    e2 <- tryCatch(pc_prepare(dna, est = est, xformla = ~x, ns = ordinary), error = identity)
    expect_s3_class(e1, "cdid_chained_request_covariate_error")
    expect_identical(conditionMessage(e1), conditionMessage(e2))
    expect_identical(class(e1), class(e2))
  }
  fr_bad <- first[, c("id", "G", "x")]; fr_bad$x[2L] <- fr_bad$x[2L] + 1
  for (est in c("ipw", "dr")) {
    e1 <- tryCatch(pc_prepare(d, est = est, xformla = ~x, frame = fr_bad), error = identity)
    e2 <- tryCatch(pc_prepare(d, est = est, xformla = ~x, frame = fr_bad, ns = ordinary), error = identity)
    expect_s3_class(e1, "cdid_chained_request_covariate_error")
    expect_identical(conditionMessage(e1), conditionMessage(e2))
  }
  # A genuinely changing covariate still becomes the dated reference route on both, and is never eligible.
  dv <- d; dv$x[dv$id == d$id[1L] & dv$period == 2L] <- 9
  for (est in c("ipw", "reg")) {
    p <- pc_prepare(dv, est = est, xformla = ~x); q <- pc_prepare(dv, est = est, xformla = ~x, ns = ordinary)
    expect_identical(p$weight_reduction, "reference_recipes")
    expect_identical(p$reference_route, q$reference_route)
    expect_null(p$prepared)
  }
  # Duplicate (id, period) rows are refused on the values themselves; identifiers
  # that only a 15-digit string rendering would merge stay distinct units.
  expect_error(pc_prepare(rbind(d, d[1L, ]), est = "ipw"), class = "cdid_chained_request_input_error")
  dz <- d[d$id %in% unique(d$id)[1:2], ]; dz$period[1L] <- 0; dz$period[2L] <- -0
  expect_error(pc_prepare(dz, est = "ipw"), class = "cdid_chained_request_input_error")
  near <- d[d$id %in% unique(d$id)[1:2], ]
  near$id <- ifelse(near$id == unique(near$id)[1L], 1, 1 + 2^-50)
  pn <- pc_prepare(near, est = "ipw")
  expect_identical(nrow(pn$units), 2L)
  expect_identical(as.character(unique(near$id)), c("1", "1")) # the old formatted key could not tell them apart
})

# ---------------------------------------------------------------------------
# 4. Pipeline binding: the prepared route is taken and yields the original route's outputs exactly
# ---------------------------------------------------------------------------

test_that("the pipeline binds the record by identity, takes the prepared constructors and reproduces the original route's estimates, SEs, IFs and covariance", {
  for (case in list(list(est = "ipw", xf = NULL), list(est = "ipw", xf = ~x), list(est = "reg", xf = NULL),
                    list(est = "reg", xf = ~x), list(est = "dr", xf = ~1), list(est = "dr", xf = ~x + z))) {
    d <- pc_panel(N = 36L, T = 5L, seed = 29L, ids = sample(2000 + seq_len(36L)), rows_perm = TRUE)
    prep <- pc_prepare(d, est = case$est, xformla = case$xf)
    expect_s3_class(prep$prepared, "cdid_request_prepared_intake")
    counted <- pc_clone(wrap = pc_route_counters)
    res_new <- pc_run(prep, ns = counted, with_prepared = TRUE)
    expect_identical(pc_n(counted, "chained_request_prepared_bind"), 1)
    expect_identical(pc_n(counted, "construct_ipw_panel_owner_prepared"), 1)
    expect_identical(pc_n(counted, "construct_ipw_frame_owner_prepared"), 1)
    expect_identical(pc_n(counted, "construct_cp04_frame_reader_prepared"), if (identical(case$est, "ipw")) 0 else 1)
    expect_identical(pc_n(counted, "construct_ipw_panel_owner"), 0)
    expect_identical(pc_n(counted, "construct_ipw_frame_owner"), 0)
    expect_identical(pc_n(counted, "construct_cp04_frame_reader"), 0)
    expect_gte(pc_n(counted, "pipeline_occupancy_accumulate_complete"), 1)
    original <- pc_clone(wrap = pc_route_counters)
    res_old <- pc_run(prep, ns = original, with_prepared = FALSE)
    expect_identical(pc_n(original, "construct_ipw_panel_owner_prepared"), 0)
    expect_identical(pc_n(original, "construct_ipw_panel_owner"), 1)
    expect_identical(res_new$reporting$estimates, res_old$reporting$estimates)
    expect_identical(res_new$reporting$se, res_old$reporting$se)
    expect_identical(res_new$reporting$if_full, res_old$reporting$if_full)
    expect_identical(res_new$reporting$vcov, res_old$reporting$vcov)
    expect_identical(res_new$reporting$target_meta, res_old$reporting$target_meta)
    expect_identical(pc_strip(res_new$panel_owner), pc_strip(res_old$panel_owner))
    expect_identical(pc_strip(res_new$frame_owner), pc_strip(res_old$frame_owner))
    expect_true(all(is.finite(res_new$reporting$estimates[res_new$reporting$available])))
  }
  # Genuine public requests: the exported entry with the record bound versus the
  # same entry in a clone whose binding declines (original constructors).
  skip_if_not_installed("did")
  d <- pc_panel(N = 36L, T = 5L, seed = 31L)
  declined <- pc_clone(wrap = pc_route_counters,
    overrides = list(chained_request_prepared_bind = function(prepared, ...) NULL))
  taken <- pc_clone(wrap = pc_route_counters)
  for (case in list(list(est = "ipw", xf = ~x), list(est = "reg", xf = ~1), list(est = "dr", xf = ~x))) {
    args <- list(yname = "y", tname = "period", idname = "id", gname = "G", data = d, xformla = case$xf,
      est_method = case$est, control_group = "nevertreated", bstrap = FALSE, cband = FALSE, base_period = "varying",
      mode = "chained", pt_assumption = "post")
    mp_new <- do.call(taken$att_gt, args)
    mp_old <- do.call(declined$att_gt, args)
    expect_s3_class(mp_new, "cdid_MP")
    expect_identical(mp_new$att, mp_old$att)
    expect_identical(mp_new$se, mp_old$se)
    expect_identical(mp_new$inffunc, mp_old$inffunc)
    expect_identical(mp_new$V_analytical, mp_old$V_analytical)
    expect_true(all(is.finite(mp_new$att)))
  }
  expect_identical(pc_n(taken, "construct_ipw_panel_owner_prepared"), 3)
  expect_identical(pc_n(taken, "construct_ipw_panel_owner"), 0)
  expect_identical(pc_n(declined, "construct_ipw_panel_owner_prepared"), 0)
  expect_identical(pc_n(declined, "construct_ipw_panel_owner"), 3)
  # The public MP record is identical apart from the call object of each clone.
  expect_identical(mp_new$chained$recoding, mp_old$chained$recoding)
  expect_identical(mp_new$chained$frame, mp_old$chained$frame)
})

# ---------------------------------------------------------------------------
# 5. Negative eligibility keeps the original route (no record is built)
# ---------------------------------------------------------------------------

test_that("character or decorated IDs, context formulas, non-unit or dated weights, incomplete or extra rosters, several treated cohorts, nonfinite outcomes and reference routes build no record", {
  d <- pc_panel(seed = 41L)
  expect_s3_class(pc_prepare(d, est = "ipw", xformla = ~x)$prepared, "cdid_request_prepared_intake")
  expect_s3_class(pc_prepare(d, est = "reg", xformla = ~x, weightsname = "one")$prepared, "cdid_request_prepared_intake") # unit weights
  dc <- d; dc$id <- sprintf("unit-%04d", dc$id)
  expect_null(pc_prepare(dc, est = "ipw")$prepared)
  expect_null(pc_prepare(d, est = "ipw", xformla = ~I(x^2))$prepared)
  expect_null(pc_prepare(d, est = "dr", xformla = ~f)$prepared)
  expect_null(pc_prepare(d, est = "reg", xformla = ~x, weightsname = "w")$prepared)
  dw <- d; dw$w[dw$period == 3L] <- dw$w[dw$period == 3L] + 1
  expect_null(api$chained_request_prepare_inputs(data = dw, yname = "y", tname = "period", idname = "id",
    gname = "G", weightsname = "w", fix_weights = "first_period")$prepared)
  expect_null(pc_prepare(pc_panel(singleton = TRUE, seed = 41L), est = "ipw")$prepared)
  # Still pair-eligible but not calendar-complete: keep the original route.
  expect_null(pc_prepare(d[-3L, ], est = "ipw")$prepared)
  dn <- d; dn$y[3L] <- NA_real_
  expect_null(pc_prepare(dn, est = "dr", xformla = ~x)$prepared)
  fr_extra <- rbind(d[!duplicated(d$id), c("id", "G", "x")], data.frame(id = 99999L, G = 0, x = 1))
  expect_null(pc_prepare(d, est = "ipw", xformla = ~x, frame = fr_extra)$prepared)
  expect_null(pc_prepare(pc_panel(second = TRUE, seed = 41L), est = "dr", xformla = ~x)$prepared)
  dinf <- d; dinf$y[5L] <- Inf
  expect_null(pc_prepare(dinf, est = "ipw")$prepared)
  dv <- d; dv$x[dv$id == d$id[1L] & dv$period == 2L] <- 9
  expect_null(pc_prepare(dv, est = "reg", xformla = ~x)$prepared)
  # The ineligible inputs still run the original constructors through the pipeline unchanged.
  counted <- pc_clone(wrap = pc_route_counters)
  prep <- pc_prepare(pc_panel(second = TRUE, seed = 41L), est = "ipw", xformla = ~x)
  res <- pc_run(prep, ns = counted)
  expect_identical(pc_n(counted, "construct_ipw_panel_owner"), 1)
  expect_identical(pc_n(counted, "construct_ipw_panel_owner_prepared"), 0)
  expect_identical(pc_n(counted, "pipeline_occupancy_accumulate_complete"), 1) # the complete-roster closed form is owner-level
  expect_s3_class(res, "cdid_ipw_pipeline_result")
})

# ---------------------------------------------------------------------------
# 6. Binding and mutation: the record is accepted only for its own objects
# ---------------------------------------------------------------------------

test_that("the binding declines any other object, mutated data, tampered records or unlicensed arguments, and the prepared constructors refuse foreign inputs", {
  d <- pc_panel(seed = 47L)
  prep <- pc_prepare(d, est = "dr", xformla = ~x)
  bind <- function(prepared = prep$prepared, data = prep$data_internal, frame_id = prep$frame_ids, frame_Z = prep$frame_Z,
                   frame_weights = prep$frame_w, frame_g = prep$frame_g_arg, allow_frame_only = FALSE, weightsname = NULL,
                   clustername = NULL, cp04 = prep$cp04, custom = NULL)
    api$chained_request_prepared_bind(prepared, data, "id", "period", "y", "G", frame_id, frame_Z, frame_weights,
      frame_g, allow_frame_only, weightsname, clustername, cp04, custom)
  expect_true(api$gmm_sexp_aliases(bind(), prep$prepared))
  expect_null(bind(prepared = NULL))
  mutated <- prep$data_internal; mutated$y[1L] <- mutated$y[1L] + 1
  expect_null(bind(data = mutated))
  expect_null(bind(data = prep$data_internal[seq_len(nrow(prep$data_internal)), ]))
  expect_null(bind(frame_Z = prep$frame_Z + 0))
  expect_null(bind(frame_id = prep$frame_ids + 0))
  expect_null(bind(frame_weights = rep(1, length(prep$frame_ids))))
  expect_null(bind(frame_g = prep$frame_g_arg + 0))
  expect_null(bind(allow_frame_only = TRUE))
  expect_null(bind(weightsname = "w"))
  expect_null(bind(clustername = "cl"))
  expect_null(bind(custom = function(...) NULL))
  cp_other <- prep$cp04; cp_other$Z <- prep$cp04$Z + 0
  expect_null(bind(cp04 = cp_other))
  cp_ref <- prep$cp04; cp_ref$reference_source <- list(V = prep$frame_w)
  expect_null(bind(cp04 = cp_ref))
  expect_null(bind(cp04 = NULL)) # an IPW call cannot bind an OR/DR record
  tampered <- prep$prepared; tampered$row_order <- rev(tampered$row_order)
  expect_null(bind(prepared = tampered))
  relabeled <- prep$prepared; class(relabeled) <- "cdid_request_prepared_intake"
  expect_null(bind(prepared = relabeled)) # the copied list is no longer the sealed snapshot
  # The IPW record binds only an IPW call.
  prep_ipw <- pc_prepare(d, est = "ipw", xformla = ~x)
  expect_null(api$chained_request_prepared_bind(prep_ipw$prepared, prep_ipw$data_internal, "id", "period", "y", "G",
    prep_ipw$frame_ids, prep_ipw$frame_Z, prep_ipw$frame_w, prep_ipw$frame_g_arg, FALSE, NULL, NULL, prep$cp04, NULL))
  # Prepared constructors refuse inputs the record does not describe.
  ctrl <- pc_ctrl()
  expect_error(api$construct_ipw_panel_owner_prepared(mutated, "id", "period", "y", "G", ctrl, prepared = prep$prepared),
    class = "cdid_ipw_owner_input_error")
  po <- api$construct_ipw_panel_owner_prepared(prep$data_internal, "id", "period", "y", "G", ctrl, prepared = prep$prepared)
  po_alien <- api$construct_ipw_panel_owner(prep$data_internal, "id", "period", "y", "G", ctrl)
  expect_error(api$construct_ipw_frame_owner_prepared(po_alien, prep$frame_ids, prep$frame_Z, prep$frame_w, ctrl,
    prep$frame_g_arg, prepared = prep$prepared, design_on_A_only = TRUE), class = "cdid_ipw_owner_input_error")
  expect_error(api$construct_ipw_frame_owner_prepared(po, prep$frame_ids, prep$frame_Z + 0, prep$frame_w, ctrl,
    prep$frame_g_arg, prepared = prep$prepared, design_on_A_only = TRUE), class = "cdid_ipw_owner_input_error")
  fo <- api$construct_ipw_frame_owner_prepared(po, prep$frame_ids, prep$frame_Z, prep$frame_w, ctrl, prep$frame_g_arg,
    prepared = prep$prepared, design_on_A_only = TRUE)
  expect_error(api$construct_cp04_frame_reader_prepared(fo, prep$frame_ids, prep$cp04$Z, prep$cp04$frame_g,
    rep(1, length(prep$frame_ids)), ctrl, prepared = prep$prepared, U = prep$cp04$U, rho = prep$cp04$rho,
    allow_unknown_off_A = TRUE), class = "cdid_cp04_integrity_error")
  expect_error(api$construct_cp04_frame_reader_prepared(fo, prep$frame_ids, prep$cp04$Z, prep$cp04$frame_g,
    prep$frame_w, ctrl, prepared = tampered, U = prep$cp04$U, rho = prep$cp04$rho, allow_unknown_off_A = TRUE),
    class = "cdid_cp04_integrity_error")
  # Caller data are unchanged by preparation and construction.
  expect_identical(d, pc_panel(seed = 47L))
})

# ---------------------------------------------------------------------------
# 7. Exact adjacent-bound refusals are shared with the standalone constructors
# ---------------------------------------------------------------------------

test_that("prepared constructors admit and refuse at exactly the standalone preflight and budget bounds", {
  d <- pc_panel(N = 23L, T = 5L, seed = 53L)
  prep <- pc_prepare(d, est = "dr", xformla = ~x)
  n <- nrow(prep$data_internal); N <- length(prep$frame_ids); p <- ncol(prep$frame_Z)
  same_refusal <- function(f_old, f_new, cls = "cdid_ipw_owner_bound_error") {
    e1 <- tryCatch(f_old(), error = identity); e2 <- tryCatch(f_new(), error = identity)
    expect_s3_class(e1, cls); expect_s3_class(e2, cls)
    expect_identical(conditionMessage(e1), conditionMessage(e2))
    expect_identical(class(e1), class(e2))
  }
  # Panel owner: work bound 128 (1 + ceil(log2(n + 1))) n + 256 and retained bound 64 n + 256.
  W <- 128 * (1 + ceiling(log2(n + 1))) * n + 256
  R0 <- 64 * n + 256
  po_old <- function(ctrl) api$construct_ipw_panel_owner(prep$data_internal, "id", "period", "y", "G", ctrl)
  po_new <- function(ctrl) api$construct_ipw_panel_owner_prepared(prep$data_internal, "id", "period", "y", "G", ctrl,
    prepared = prep$prepared)
  same_refusal(function() po_old(pc_ctrl(work = W - 1)), function() po_new(pc_ctrl(work = W - 1)))
  expect_identical(pc_strip(po_old(pc_ctrl(work = W))), pc_strip(po_new(pc_ctrl(work = W))))
  same_refusal(function() po_old(pc_ctrl(retained = R0 - 1)), function() po_new(pc_ctrl(retained = R0 - 1)))
  expect_identical(pc_strip(po_old(pc_ctrl(retained = R0))), pc_strip(po_new(pc_ctrl(retained = R0))))
  # Frame owner: retained bound n_a p + n_a + n_u + 2 n_u + 256 with F = A.
  ctrl <- pc_ctrl()
  po <- po_new(ctrl); po_s <- po_old(ctrl)
  Rf <- N * p + N + N + 2 * N + 256
  fo_old <- function(c2) api$construct_ipw_frame_owner(po_s, prep$frame_ids, prep$frame_Z, weights = prep$frame_w,
    control = c2, frame_g = prep$frame_g_arg, design_on_A_only = TRUE)
  fo_new <- function(c2) api$construct_ipw_frame_owner_prepared(po, prep$frame_ids, prep$frame_Z, prep$frame_w, c2,
    prep$frame_g_arg, prepared = prep$prepared, design_on_A_only = TRUE)
  same_refusal(function() fo_old(pc_ctrl(retained = Rf - 1)), function() fo_new(pc_ctrl(retained = Rf - 1)))
  expect_identical(pc_strip(fo_old(pc_ctrl(retained = Rf))), pc_strip(fo_new(pc_ctrl(retained = Rf))))
  # Reader: retained budget borrowed + reader bound (both constructors charge it before any scan).
  fo <- fo_new(ctrl); fo_s <- fo_old(ctrl)
  borrowed <- po$accounting$actual_retained_entries + fo$accounting$retained_matrix_entries + fo$accounting$retained_vector_entries
  dsel <- ncol(prep$cp04$Z)
  meta <- pc_ref_payload(list(dimnames(prep$cp04$Z), attributes(prep$frame_ids)))
  Rr <- borrowed + N * (dsel + 14) + 8 * N + 4096 + unname(meta["entries"])
  rd_old <- function(c2) api$construct_cp04_frame_reader(fo_s, prep$frame_ids, prep$cp04$Z, prep$cp04$frame_g, prep$frame_w,
    c2, U = prep$cp04$U, rho = prep$cp04$rho, allow_unknown_off_A = TRUE)
  rd_new <- function(c2) api$construct_cp04_frame_reader_prepared(fo, prep$frame_ids, prep$cp04$Z, prep$cp04$frame_g,
    prep$frame_w, c2, prepared = prep$prepared, U = prep$cp04$U, rho = prep$cp04$rho, allow_unknown_off_A = TRUE)
  same_refusal(function() rd_old(pc_ctrl(retained = Rr - 1)), function() rd_new(pc_ctrl(retained = Rr - 1)))
  expect_identical(pc_strip(rd_old(pc_ctrl(retained = Rr))), pc_strip(rd_new(pc_ctrl(retained = Rr))))
  expect_identical(rd_new(ctrl)$accounting$owned_entries, unname(pc_ref_payload(pc_strip_ref_slot(rd_new(ctrl)))["entries"]))
})

# ---------------------------------------------------------------------------
# 8. Complete-roster occupancy: closed form equals the general algorithm and brute force
# ---------------------------------------------------------------------------

test_that("the complete-roster occupancy closed form equals the general accumulation (cells, types and diagnostics) and brute-force pair counts", {
  general <- pc_clone(overrides = list(pipeline_occupancy_complete_roster = function(p_own) FALSE))
  ctrl <- pc_ctrl(5L)
  panels <- list(
    complete = pc_panel(N = 23L, T = 5L, seed = 59L),
    complete_two_cohorts = pc_panel(N = 23L, T = 5L, seed = 59L, second = TRUE),
    complete_T2 = pc_panel(N = 9L, T = 2L, seed = 61L),
    one_unit = pc_panel(N = 1L, T = 4L, seed = 61L),
    own_cohorts = local({ d <- pc_panel(N = 11L, T = 4L, seed = 67L); d$G <- 10 + d$id; d }),
    incomplete = local({ d <- pc_panel(N = 23L, T = 5L, seed = 59L); d[-c(3L, 17L, 40L), ] }),
    missing_outcome = local({ d <- pc_panel(N = 23L, T = 5L, seed = 59L); d$y[c(2L, 50L)] <- NA_real_; d }))
  for (nm in names(panels)) {
    d <- panels[[nm]]
    po <- api$construct_ipw_panel_owner(transform(d, G = ifelse(G == 0, Inf, G)), "id", "period", "y", "G", ctrl)
    complete <- !nm %in% c("incomplete", "missing_outcome")
    expect_identical(api$pipeline_occupancy_complete_roster(po), complete, info = nm)
    occ <- api$derive_owner_sparse_occupancy(po, ctrl)
    occ_general <- general$derive_owner_sparse_occupancy(po, ctrl)
    expect_identical(pc_strip(occ), pc_strip(occ_general), info = nm)
    expect_identical(occ$occupancy, pc_brute_occupancy(d), info = nm)
    expect_equal(as.double(occ$diagnostics$total_occupancy_count), as.double(po$counts$n_pairs), info = nm)
    pre <- api$preflight_owner_sparse_occupancy(po, ctrl)
    acc <- api$pipeline_occupancy_accumulate(po, pre)
    acc_general <- general$pipeline_occupancy_accumulate(po, pre)
    expect_identical(acc, acc_general, info = nm)
    expect_type(acc$cohort_a_index, "integer"); expect_type(acc$s_index, "integer")
    expect_type(acc$t_index, "integer"); expect_type(acc$count, "double")
  }
  # The accumulation plan refusal is raised at the same reservation on both paths.
  d <- panels$complete
  po <- api$construct_ipw_panel_owner(transform(d, G = ifelse(G == 0, Inf, G)), "id", "period", "y", "G", ctrl)
  pre <- api$preflight_owner_sparse_occupancy(po, ctrl)
  tight <- pre; tight$accumulation_scratch_entries <- 0
  e1 <- tryCatch(api$pipeline_occupancy_accumulate(po, tight), error = identity)
  e2 <- tryCatch(general$pipeline_occupancy_accumulate(po, tight), error = identity)
  expect_s3_class(e1, "cdid_pipeline_metadata_internal_error")
  expect_identical(conditionMessage(e1), conditionMessage(e2))
})
