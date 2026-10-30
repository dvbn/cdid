# CP08 validated catalog-to-design handoff and stage-local dependency demand
# (R/chained_pipeline_metadata.R, R/chained_comparison_catalog.R,
# R/chained_design.R): one adapter preflight, one declared-qualification
# validation, one complete link validation and the unchanged graph core on the
# normal same-owner sealed-source route; helper groups demanded at the stage
# that calls them.
#
# UNEXECUTED at authoring time (no R in the authoring job). Every expected value
# comes from the independent packaged references (hash-bound bit-register
# comparator, Cartesian catalog enumeration, brute-force occupancy), from the
# unchanged standalone constructors of the same installed build, or from
# disclosed test-local namespace clones that force the generic bridge. Counter
# clones are supplementary: the genuine installed namespace runs every parity
# check and the public six-method family.
library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
if (!requireNamespace("digest", quietly = TRUE)) stop("Declared digest test dependency required.", call. = FALSE)
api <- asNamespace("cdid")
stopifnot(isNamespace(api), environmentIsLocked(api))

# ---------------------------------------------------------------------------
# Independent references (base R only, hash-bound) and fixtures
# ---------------------------------------------------------------------------

ch_load_reference <- function(dir, file, sha) {
  path <- testthat::test_path(dir, file)
  if (!file.exists(path)) stop("missing packaged reference: ", file, call. = FALSE)
  if (!identical(digest::digest(file = path, algo = "sha256"), sha)) stop("reference hash mismatch: ", file, call. = FALSE)
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  env
}
ref_cal <- ch_load_reference("reference", "cp02_calendar_design_reference.R",
  "ed66dddd0b84ebe4f0e4bb1ef64eb0866d92eab0f1be2b697762279dcdda3031")
ref_cat <- ch_load_reference("reference", "cp02_comparison_catalog_reference.R",
  "2307f19615bc3bc98d664a072a036c5a7c9828b5b93782a53c62cbcbd6b9ac58")
ref_fix <- ch_load_reference("fixtures", "pipeline_metadata_reference.R",
  "fcdb1169e06eb4c142d94278f69fd95fcd502ea5cc3d5036871fcade6f39a6ed")
ch_cmp <- ref_cal$ref_exact_cmp_sum

ch_ctrl <- function(chunk = 7L, retained = 1e7, matrix = 1e7, work = 1e10, output = 1e7) api$ipw_owner_control(
  max_chunk_rows = chunk, max_retained_entries = retained, max_matrix_entries = matrix,
  max_work_ops = work, max_output_entries = output)
ch_bounds <- list(max_candidates = 1e6, max_rows = 1e6, max_pool_entries = 1e7, max_scratch = 1e8, max_work = 1e9)

# Complete balanced long panel in owner coding (Inf never treated): cohort 3 for
# id %% 3 == 0, deterministic nonadditive outcome, one unit covariate; no RNG.
ch_panel <- function(N = 41L, T = 20L) {
  d <- expand.grid(time = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "time")]
  d$time <- as.double(d$time)
  d$g <- ifelse(d$id %% 3L == 0L, 3, Inf)
  d$y <- 0.3 * d$id + 0.2 * d$time + ifelse(d$time >= d$g, 1, 0) + ((7 * d$id + 3 * d$time) %% 5) / 10
  d$x <- (d$id %% 5) / 4
  rownames(d) <- NULL
  d
}
# The same shape in public coding (G = 0 never treated) for `cdid::att_gt`.
ch_public_panel <- function(N = 41L, T = 8L) {
  d <- expand.grid(period = seq_len(T), id = seq_len(N), KEEP.OUT.ATTRS = FALSE)[, c("id", "period")]
  d$G <- ifelse(d$id %% 3L == 0L, 3, 0)
  d$x <- (d$id %% 5) / 4
  d$y <- 0.3 * d$id + 0.2 * d$period + 0.15 * d$x * d$period + ifelse(d$G > 0 & d$period >= d$G, 1.5, 0) +
    ((7 * d$id + 3 * d$period) %% 5) / 10
  rownames(d) <- NULL
  d
}

ch_owners <- function(d, ctrl, with_x = FALSE) {
  po <- api$construct_ipw_panel_owner(d, "id", "time", "y", "g", ctrl)
  u <- d[!duplicated(d$id), ]
  u <- u[order(u$id), ]
  Z <- if (with_x) cbind(1, u$x) else matrix(1, nrow(u), 1L)
  fo <- api$construct_ipw_frame_owner(po, u$id, Z, weights = NULL, control = ctrl)
  list(po = po, fo = fo, units = u)
}

ch_build <- function(ns, own, targets, control_group, conventions, ctrl = ch_ctrl(), bounds = ch_bounds,
                     anticipation = 0, requested_post = NULL, ...) {
  ns$build_owner_comparison_catalog(panel_owner = own$po, targets = targets, anticipation = anticipation,
    control_group = control_group, conventions = conventions, requested_post = requested_post,
    control = ctrl, catalog_bounds = bounds, frame_owner = own$fo, ...)
}

# The design inputs exactly as the normal pipeline derives them from the first
# catalog: every conditional causal endpoint is declared; every calendar date
# at or after a row-bearing cohort is a requested post cell.
ch_design_inputs <- function(meta) {
  rows <- meta$catalog$rows
  s_c <- rows$s_role %in% c("anticipation_conditional", "post_conditional")
  t_c <- rows$t_role %in% c("anticipation_conditional", "post_conditional")
  dq <- unique(data.frame(g = as.double(c(rows$g[s_c], rows$g[t_c])), t = as.double(c(rows$s[s_c], rows$t[t_c])),
    stringsAsFactors = FALSE))
  dq <- dq[order(dq$g, dq$t), , drop = FALSE]
  rownames(dq) <- NULL
  if (nrow(dq) == 0L) dq <- NULL
  cal <- meta$catalog$spec$calendar
  row_g <- sort(unique(as.double(rows$g)))
  tpt <- do.call(rbind, lapply(row_g, function(g) data.frame(g = g, t = cal[cal >= g], stringsAsFactors = FALSE)))
  rownames(tpt) <- NULL
  list(dq = dq, tpt = tpt)
}

ch_strip <- function(x) {
  if (is.environment(x)) return(NULL)
  if (is.list(x)) {
    a <- attributes(x)
    y <- lapply(x, ch_strip)
    attributes(y) <- a
    return(y)
  }
  x
}
ch_err <- function(expr) tryCatch({ expr; NULL }, error = identity)
ch_same_refusal <- function(a, b, info) {
  expect_true(inherits(a, "error"), info = info)
  expect_true(inherits(b, "error"), info = info)
  expect_identical(class(a), class(b), info = info)
  expect_identical(conditionMessage(a), conditionMessage(b), info = info)
}
ch_key <- function(...) paste(vapply(list(...), function(v) paste(sprintf("%a", as.double(v)), collapse = ","), character(1L)), collapse = "|")

# Disclosed test-local clone of the namespace: every internal function re-bound
# into `e` (counters, overrides). With `parent = api` (the counting clone) the
# genuine namespace stays below the clone; with the imports environment as
# parent (the sandbox) the genuine namespace is NOT in the lookup chain, so a
# helper removed from the sandbox is genuinely absent for its dependency checks.
ch_clone <- function(wrap = character(), overrides = list(), parent = api) {
  e <- new.env(parent = parent)
  counts <- new.env(parent = emptyenv())
  genuine <- "gmm_require_native_dll"
  skip <- c(".__NAMESPACE__.", ".__S3MethodsTable__.")
  for (nm in setdiff(ls(api, all.names = TRUE), skip)) {
    f <- get(nm, envir = api, inherits = FALSE)
    if (is.function(f) && !is.primitive(f) && identical(environment(f), api) && !nm %in% genuine) {
      environment(f) <- e
      assign(nm, f, envir = e)
    } else if (!identical(parent, api)) {
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
  # Overrides first, so a counted name counts the override that replaces it.
  for (nm in names(overrides)) {
    f <- overrides[[nm]]
    environment(f) <- e
    assign(nm, f, envir = e)
  }
  for (nm in wrap) {
    if (!exists(nm, envir = e, inherits = FALSE)) stop("No internal function `", nm, "` to count.", call. = FALSE)
    assign(nm, 0, envir = counts)
    assign(nm, make_wrap(nm, get(nm, envir = e, inherits = FALSE)), envir = e)
  }
  e$.counts <- counts
  e
}
ch_n <- function(e, nm) get0(nm, envir = e$.counts, inherits = FALSE, ifnotfound = 0)
ch_snapshot <- function(e) mget(ls(e$.counts, all.names = TRUE), envir = e$.counts)
ch_delta <- function(before, after, nm) after[[nm]] - before[[nm]]
# Records every helper-group request of the clone's stage demand.
ch_stage_recorder <- function(e) {
  rec <- new.env(parent = emptyenv())
  rec$groups <- list()
  inner <- get("pipeline_metadata_require_stage", envir = e, inherits = FALSE)
  assign("pipeline_metadata_require_stage", function(groups) {
    rec$groups[[length(rec$groups) + 1L]] <- groups
    inner(groups)
  }, envir = e)
  rec
}
ch_generic_bridge <- function(handoff, target_post_times = NULL, pt_assumption = "post", pt_start = NULL) {
  build_calendar_chained_design(handoff$links, calendar = handoff$calendar, cohorts = handoff$cohorts,
    anticipation = handoff$anticipation, target_post_times = target_post_times,
    pt_assumption = pt_assumption, pt_start = pt_start)
}

ch_wrap_names <- c("catalog_adapter_preflight", "catalog_adapter_dimensions", "catalog_design_envelope",
  "comparison_catalog_design_links", "build_comparison_catalog_design", "catalog_design_handoff",
  "build_catalog_handoff_design", "adapt_calendar_design_links", "adapt_calendar_design_links_validation",
  "validate_design_links", "design_links_validation", "chained_design_core", "build_chained_design",
  "build_calendar_chained_design", "calendar_reference_date", "parse_control_pool_string",
  "normalize_target_post_times", "exact_cmp_sum", "build_comparison_catalog",
  "pipeline_metadata_resolve_catalog_mappings", "derive_owner_sparse_occupancy", "build_owner_comparison_catalog")

ch_fixtures <- function() {
  fx <- ref_fix$make_pipeline_metadata_fixtures()
  g3 <- 3.0
  g3_adj <- 3.0 + 2 * .Machine$double.eps
  list(
    complete_190 = list(d = ch_panel(41L, 20L), targets = 3, cg = "nevertreated", anticipation = 0),
    rotating_nyt = list(d = fx$rotating_gapped, targets = c(3, 4), cg = "notyettreated", anticipation = 0,
      requested_post = data.frame(g = c(3, 4), t = c(3, 4))),
    rotating_nyt_delta1 = list(d = fx$rotating_gapped, targets = c(3, 4), cg = "notyettreated", anticipation = 1),
    singletons_nyt = list(d = fx$with_singletons, targets = c(3, 4, 5), cg = "notyettreated", anticipation = 0),
    adjacent_nyt = list(d = fx$adjacent_binary64, targets = c(g3, g3_adj), cg = "notyettreated", anticipation = 0),
    balanced_nt = list(d = fx$balanced, targets = c(3, 4), cg = "nevertreated", anticipation = 0)
  )
}

# ---------------------------------------------------------------------------
# 1. Parity: the admitted same-source route, the normal fallback and the
#    standalone wrappers produce identical complete metadata/design objects;
#    independent references fix references, qualification, pools and roles
# ---------------------------------------------------------------------------

test_that("fused handoff equals the generic bridge and standalone constructors on every metadata/design field, with independent references", {
  ctrl <- ch_ctrl()
  for (nm in names(ch_fixtures())) {
    fx <- ch_fixtures()[[nm]]
    own <- ch_owners(fx$d, ctrl)
    conv <- paste0("conv:handoff:", nm)
    meta0 <- ch_build(api, own, fx$targets, fx$cg, conv, anticipation = fx$anticipation,
      requested_post = fx$requested_post, build_design = FALSE)
    expect_true(api$ipw_owner_seal_ok(meta0, "pipeline_metadata"), info = nm)
    expect_null(meta0$design, info = nm)
    expect_null(meta0$design_links, info = nm)
    inp <- ch_design_inputs(meta0)
    cat <- meta0$catalog
    cal <- cat$spec$calendar; coh <- cat$spec$cohorts; delta <- cat$spec$anticipation
    expect_identical(delta, as.double(fx$anticipation), info = nm)

    # Admitted same-source route and the normal fallback (no source).
    meta_reuse <- ch_build(api, own, fx$targets, fx$cg, conv, anticipation = fx$anticipation,
      requested_post = fx$requested_post, occupancy = meta0$occupancy_object, declared_qualified = inp$dq,
      target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0)
    meta_plain <- ch_build(api, own, fx$targets, fx$cg, conv, anticipation = fx$anticipation,
      requested_post = fx$requested_post, occupancy = meta0$occupancy_object, declared_qualified = inp$dq,
      target_post_times = inp$tpt, build_design = TRUE)
    expect_true(api$ipw_owner_seal_ok(meta_reuse, "pipeline_metadata"), info = nm)
    expect_true(isTRUE(meta_reuse$diagnostics$catalog_reuse$reused), info = nm)
    expect_false(isTRUE(meta_plain$diagnostics$catalog_reuse$reused), info = nm)
    expect_identical(meta_reuse$catalog_source_ref, meta0$ref, info = nm)
    expect_identical(meta_reuse$catalog, meta0$catalog, info = nm)
    expect_identical(meta_reuse$row_mapping, meta0$row_mapping, info = nm)
    expect_identical(meta_reuse$operator_mapping, meta0$operator_mapping, info = nm)
    expect_identical(meta_reuse$occupancy_ref, meta0$occupancy_ref, info = nm)
    for (f in c("catalog", "row_mapping", "operator_mapping", "design_links", "design", "occupancy", "cell_indices",
      "owner_summary", "control", "construction")) {
      expect_identical(ch_strip(meta_reuse[[f]]), ch_strip(meta_plain[[f]]), info = paste(nm, f))
    }
    expect_identical(meta_reuse$diagnostics$catalog_reuse$owned_design_entries,
      meta_plain$diagnostics$catalog_reuse$owned_design_entries, info = nm)
    expect_true(meta_reuse$diagnostics$all_catalog_rows_retained, info = nm)
    expect_true(meta_reuse$diagnostics$all_catalog_operators_retained, info = nm)
    expect_false(meta_reuse$diagnostics$formatted_float_keys_used, info = nm)
    expect_identical(meta_reuse$diagnostics$endpoints_read, 0L, info = nm)

    # Standalone wrappers (strict preflights) and the former two-pass
    # composition give the identical links table and design object.
    links_ref <- api$comparison_catalog_design_links(cat, inp$dq)
    design_ref <- api$build_comparison_catalog_design(cat, inp$dq, inp$tpt)
    design_cal <- api$build_calendar_chained_design(links_ref, cal, coh, delta, target_post_times = inp$tpt)
    design_two_pass <- api$build_chained_design(api$adapt_calendar_design_links(links_ref, cal, coh, delta), inp$tpt)
    expect_identical(meta_reuse$design_links, design_ref$links, info = nm)
    for (col in names(links_ref)) expect_identical(meta_reuse$design$links[[col]], links_ref[[col]], info = paste(nm, col))
    expect_identical(meta_reuse$design, design_ref, info = nm)
    expect_identical(design_ref, design_cal, info = nm)
    expect_identical(design_ref, design_two_pass, info = nm)
    expect_s3_class(meta_reuse$design, "cdid_chained_design")
    # The handoff packet itself: the bridge table, canonical calendar fields and typed pools.
    h <- api$catalog_design_handoff(cat, inp$dq)
    expect_identical(names(h), c("links", "calendar", "cohorts", "anticipation", "pools"), info = nm)
    expect_identical(h$links, links_ref, info = nm)
    expect_identical(h$calendar, cal, info = nm); expect_identical(h$cohorts, coh, info = nm)
    expect_identical(h$anticipation, delta, info = nm)
    expect_identical(h$pools$index, cat$rows$pool_id, info = nm)
    expect_type(h$pools$index, "integer")
    expect_identical(length(h$pools$values), cat$pools$n_pools, info = nm)
    expect_identical(length(h$pools$text), cat$pools$n_pools, info = nm)
    expect_false(anyDuplicated(h$pools$text) > 0L, info = nm)
    expect_identical(h$pools$text[h$pools$index], h$links$control_pool, info = nm)
    pl <- cat$pools
    for (p in seq_len(pl$n_pools)) {
      K <- coh[pl$entries[seq.int(pl$start[p], length.out = pl$size[p])]]
      expect_identical(h$pools$values[[p]], as.double(K), info = paste(nm, "pool", p))
      # The typed pool IS the exact parse of the text it was rendered into.
      expect_identical(api$parse_control_pool_string(h$pools$text[p]), h$pools$values[[p]], info = paste(nm, "pool", p))
    }
    expect_identical(api$build_catalog_handoff_design(h, inp$tpt), design_ref, info = nm)

    # Independent references: exact reference dates, declared qualification
    # from the independent role comparator, pool text from the typed set, and
    # the six endpoint roles of the design vertices.
    L <- meta_reuse$design$links
    expect_identical(L$ref, ref_cal$ref_calendar_reference(L$g, cal, delta), info = nm)
    dq_keys <- if (is.null(inp$dq)) character() else paste(sprintf("%a", inp$dq$g), sprintf("%a", inp$dq$t))
    for (i in seq_len(nrow(L))) {
      K <- h$pools$values[[h$pools$index[i]]]
      for (side in c("s", "t")) {
        u <- L[[side]][i]
        role <- ref_cat$ref_cat_role(u, L$g[i], L$ref[i], K, delta, ch_cmp)
        want_q <- role == "reference" || (role %in% c("anticipation_conditional", "post_conditional") &&
          paste(sprintf("%a", L$g[i]), sprintf("%a", u)) %in% dq_keys)
        expect_identical(L[[paste0(side, "_qualified")]][i], want_q, info = paste(nm, side, i))
        v <- L[[paste0(if (side == "s") "start" else "end", "_vertex_index")]][i]
        expect_identical(meta_reuse$design$vertices$role[v],
          ref_cal$ref_endpoint_role(u, L$g[i], L$ref[i], delta, L[[paste0(side, "_qualified")]][i]), info = paste(nm, side, i))
        expect_identical(meta_reuse$design$vertices$time[v], u, info = paste(nm, side, i))
        expect_identical(meta_reuse$design$vertices$g[v], L$g[i], info = paste(nm, side, i))
      }
      expect_identical(L$control_pool[i], paste0("{", paste(sprintf("%.17g", K), collapse = ","), "}"), info = paste(nm, i))
    }
    # Charges: the fused route charges exactly what the original two adapter
    # preflights charged (links envelope and complete design envelope).
    ad <- api$catalog_adapter_preflight(cat, inp$dq, inp$tpt, design = TRUE)
    ad_links <- api$catalog_adapter_preflight(cat, inp$dq)
    expect_identical(meta_reuse$diagnostics$catalog_reuse$owned_design_entries, ad$scratch + ad_links$scratch, info = nm)
    expect_identical(meta_reuse$diagnostics$catalog_reuse$owned_work_ops, ad$work + ad_links$work, info = nm)
    expect_gt(ad$scratch, 0)
  }
})

test_that("the fused design links are the independent Cartesian catalog enumeration (operators, pools, provenance)", {
  ctrl <- ch_ctrl()
  for (nm in c("rotating_nyt", "rotating_nyt_delta1", "balanced_nt")) {
    fx <- ch_fixtures()[[nm]]
    own <- ch_owners(fx$d, ctrl)
    meta0 <- ch_build(api, own, fx$targets, fx$cg, paste0("conv:enum:", nm), anticipation = fx$anticipation,
      requested_post = fx$requested_post, build_design = FALSE)
    inp <- ch_design_inputs(meta0)
    meta <- ch_build(api, own, fx$targets, fx$cg, paste0("conv:enum:", nm), anticipation = fx$anticipation,
      requested_post = fx$requested_post, occupancy = meta0$occupancy_object, declared_qualified = inp$dq,
      target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0)
    bf <- ref_fix$brute_force_panel_occupancy(fx$d, "id", "time", "y", "g")
    want <- ref_cat$ref_catalog_enumerate(bf$calendar, bf$cohorts, fx$targets, fx$anticipation, bf$sparse_occupancy, fx$cg, ch_cmp)
    L <- meta$design$links
    h <- api$catalog_design_handoff(meta$catalog, inp$dq)
    got_keys <- vapply(seq_len(nrow(L)), function(i) ref_cat$ref_cat_key(L$g[i], L$s[i], L$t[i], h$pools$values[[h$pools$index[i]]]), character(1L))
    want_keys <- vapply(want$rows, function(x) ref_cat$ref_cat_key(x$g, x$s, x$t, x$K), character(1L))
    expect_false(anyDuplicated(got_keys) > 0L, info = nm)
    expect_identical(sort(got_keys, method = "radix"), sort(want_keys, method = "radix"), info = nm)
    expect_identical(nrow(L), length(want$rows), info = nm)
    expect_identical(nrow(meta$catalog$rows), length(want$rows), info = nm)
    # One operator per (g, K); equal pools of one target share the operator and
    # therefore the functional key carried into the design.
    op_key <- vapply(seq_len(nrow(L)), function(i) ch_key(L$g[i], h$pools$values[[h$pools$index[i]]]), character(1L))
    expect_identical(match(op_key, unique(op_key)), match(L$functional_key, unique(L$functional_key)), info = nm)
    expect_true(all(L$catalog_row == meta$catalog$rows$row_id), info = nm)
  }
})

# ---------------------------------------------------------------------------
# 2. Exact structure of the 190-row complete-roster design: rows, ordering,
#    identifiers, incidence, types and attributes
# ---------------------------------------------------------------------------

test_that("the 190-row complete-roster design keeps every row, the original link/vertex ordering, exact identifiers, incidence and attributes", {
  ctrl <- ch_ctrl()
  own <- ch_owners(ch_panel(41L, 20L), ctrl)
  meta0 <- ch_build(api, own, 3, "nevertreated", "conv:190", build_design = FALSE)
  inp <- ch_design_inputs(meta0)
  expect_identical(inp$tpt, data.frame(g = 3, t = as.double(3:20), stringsAsFactors = FALSE))
  meta <- ch_build(api, own, 3, "nevertreated", "conv:190", occupancy = meta0$occupancy_object,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0)
  des <- meta$design
  L <- des$links
  V <- des$vertices
  expect_identical(nrow(meta$catalog$rows), 190L)
  expect_identical(nrow(meta$row_mapping), 190L)
  expect_identical(nrow(meta$design_links), 190L)
  expect_identical(nrow(L), 190L)
  expect_identical(L$catalog_row, 1:190)
  expect_identical(meta$catalog$pools$n_pools, 1L)
  expect_identical(unique(L$control_pool), "{Inf}")
  raw_names <- c("g", "s", "t", "ref", "functional_key", "control_pool", "s_qualified", "t_qualified", "catalog_row")
  expect_identical(meta$design_links, L)
  expect_identical(names(L), c(raw_names, "link_index", "link_id", "start_vertex_index", "end_vertex_index",
    "start_param_index", "end_param_index", "link_type", "component_id"))
  expect_identical(attr(L, "row.names"), 1:190)
  # Link order is the catalog row order; identifiers are the canonical ranks.
  canonical <- order(L$g, L$s, L$t, L$functional_key, method = "radix")
  expect_identical(L$link_index, 1:190)
  expect_identical(L$link_id, paste0("link_", order(canonical)))
  expect_false(anyDuplicated(L$link_id) > 0L)
  expect_identical(L$s, meta$catalog$rows$s); expect_identical(L$t, meta$catalog$rows$t)
  expect_identical(L$ref, rep(2, 190L))
  # Vertex order: (cohort, kind, time, private key); indices and labels unique.
  kind_map <- c(reference = 1L, pre = 2L, anticipation_nuisance = 3L, anticipation = 4L, post_nuisance = 5L, post = 6L)
  fn <- sort(unique(L$functional_key), method = "radix")
  fv <- ifelse(V$role %in% c("pre", "anticipation_nuisance", "post_nuisance"), match(V$functional_key, fn), 0L)
  expect_identical(order(match(V$g, sort(unique(L$g))), unname(kind_map[V$role]), V$time, fv, method = "radix"), seq_len(nrow(V)))
  expect_identical(V$vertex_index, seq_len(nrow(V)))
  expect_false(anyDuplicated(V$vertex_id) > 0L)
  expect_identical(sum(V$role == "reference"), 1L)
  expect_identical(V$time[V$role == "reference"], 2)
  expect_identical(sort(V$time[V$role == "pre"]), 1)
  expect_identical(sort(V$time[V$role == "post"]), as.double(3:20))
  expect_identical(V$vertex_id[V$role == "reference"], "ref(g=3,t=2)")
  expect_true(all(paste0("ATT(g=3,t=", 3:20, ")") %in% V$vertex_id))
  # PT-post: every non-reference vertex is a parameter; incidence is the exact map.
  expect_identical(des$parameters$vertex_index, which(V$role != "reference"))
  expect_identical(des$parameters$param_id, V$vertex_id[V$role != "reference"])
  expect_identical(des$incidence$start, ifelse(is.na(L$start_param_index), 0L, L$start_param_index))
  expect_identical(des$incidence$end, ifelse(is.na(L$end_param_index), 0L, L$end_param_index))
  expect_identical(des$incidence$n_links, 190L)
  expect_identical(des$incidence$link_ids, L$link_id)
  expect_identical(des$incidence$param_ids, des$parameters$param_id)
  expect_identical(L$start_vertex_index, match(paste(L$g, L$s, ifelse(L$s_qualified, "", L$functional_key)),
    paste(V$g, V$time, ifelse(V$role %in% c("pre", "post_nuisance", "anticipation_nuisance"), V$functional_key, ""))))
  expect_identical(L$end_vertex_index, match(paste(L$g, L$t, ifelse(L$t_qualified, "", L$functional_key)),
    paste(V$g, V$time, ifelse(V$role %in% c("pre", "post_nuisance", "anticipation_nuisance"), V$functional_key, ""))))
  expect_identical(L$component_id, des$forest$component_id[L$start_vertex_index])
  expect_identical(des$post_targets$status, rep("estimable", 18L))
  expect_identical(des$rank_diagnostics$n_links, 190L)
  expect_true(des$rank_diagnostics$has_full_column_rank)
  # Types and attributes.
  expect_type(L$g, "double"); expect_type(L$s, "double"); expect_type(L$t, "double"); expect_type(L$ref, "double")
  expect_type(L$functional_key, "character"); expect_type(L$control_pool, "character")
  expect_type(L$s_qualified, "logical"); expect_type(L$t_qualified, "logical")
  expect_type(L$catalog_row, "integer"); expect_type(L$link_index, "integer"); expect_type(L$link_id, "character")
  expect_type(L$start_vertex_index, "integer"); expect_type(L$end_vertex_index, "integer")
  expect_type(L$start_param_index, "integer"); expect_type(L$end_param_index, "integer")
  expect_identical(attr(L, "calendar_spec"), des$calendar_spec)
  expect_s3_class(des$calendar_spec, "cdid_calendar_spec")
  expect_null(des$calendar_spec$pt_assumption)
  expect_identical(des$calendar_spec, api$validate_calendar_spec(meta$catalog$spec$calendar, meta$catalog$spec$cohorts, 0))
  expect_setequal(names(attributes(meta$design_links)), c("names", "class", "row.names", "calendar_spec"))
  # Every original row keeps its own retained record and mapping (no fewer local comparisons).
  expect_identical(meta$row_mapping$row_id, 1:190)
  expect_identical(meta$row_mapping$functional_key, L$functional_key)
  expect_identical(length(unique(L$link_id)), 190L)
})

# ---------------------------------------------------------------------------
# 3. Declared qualifications, PT-post/PT-all, references, targets, collisions
# ---------------------------------------------------------------------------

test_that("PT-all handoff equals the standalone PT-all design; adjacent binary64 dates and cohorts never collide", {
  ctrl <- ch_ctrl()
  own <- ch_owners(ch_panel(41L, 12L), ctrl)
  for (pt in list(list(pt_assumption = "post", pt_start = NULL), list(pt_assumption = "all", pt_start = 1))) {
    conv <- paste0("conv:pt:", pt$pt_assumption)
    meta0 <- ch_build(api, own, 3, "nevertreated", conv, build_design = FALSE,
      pt_assumption = pt$pt_assumption, pt_start = pt$pt_start)
    inp <- ch_design_inputs(meta0)
    meta <- ch_build(api, own, 3, "nevertreated", conv, occupancy = meta0$occupancy_object,
      declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0,
      pt_assumption = pt$pt_assumption, pt_start = pt$pt_start)
    ref <- api$build_comparison_catalog_design(meta0$catalog, inp$dq, inp$tpt,
      pt_assumption = pt$pt_assumption, pt_start = pt$pt_start)
    expect_identical(meta$design, ref, info = pt$pt_assumption)
    expect_identical(meta$design$calendar_spec$pt_assumption, if (pt$pt_assumption == "all") "all" else NULL)
    if (pt$pt_assumption == "all") {
      expect_false(is.null(meta$design$pt_all_spec))
      expect_identical(meta$design$pt_all_spec$pt_start, 1)
      expect_identical(meta$design$calendar_spec$pt_start, 1)
      expect_true(all(meta$design$vertices$anchored_value[meta$design$vertices$role == "pre" & meta$design$vertices$time >= 1] == 0))
    } else {
      expect_null(meta$design$pt_all_spec)
    }
    expect_identical(meta$construction$pt_assumption, pt$pt_assumption)
    expect_identical(meta$construction$pt_start, pt$pt_start)
  }
  # A different PT request is a different construction: the source is not reused.
  meta_post <- ch_build(api, own, 3, "nevertreated", "conv:pt:x", build_design = FALSE)
  inp <- ch_design_inputs(meta_post)
  e <- ch_err(ch_build(api, own, 3, "nevertreated", "conv:pt:x", occupancy = meta_post$occupancy_object,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta_post,
    pt_assumption = "all", pt_start = 1))
  expect_true(inherits(e, "cdid_pipeline_metadata_reuse_error"))

  # Adjacent binary64 values: identical in default display, distinct everywhere here.
  fx <- ch_fixtures()$adjacent_nyt
  own_adj <- ch_owners(fx$d, ctrl)
  g3 <- 3.0; g3_adj <- 3.0 + 2 * .Machine$double.eps
  t1 <- 1.0; t1_adj <- 1.0 + .Machine$double.eps
  expect_identical(format(g3), format(g3_adj)); expect_identical(as.character(t1), as.character(t1_adj))
  meta0 <- ch_build(api, own_adj, fx$targets, fx$cg, "conv:adj", build_design = FALSE)
  inp <- ch_design_inputs(meta0)
  meta <- ch_build(api, own_adj, fx$targets, fx$cg, "conv:adj", occupancy = meta0$occupancy_object,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0)
  h <- api$catalog_design_handoff(meta$catalog, inp$dq)
  expect_true(all(c(g3, g3_adj) %in% meta$design$links$g))
  expect_true(all(c(t1, t1_adj) %in% meta$design$vertices$time))
  expect_true(g3_adj %in% unlist(h$pools$values))
  expect_true(any(grepl("3.0000000000000004", h$pools$text, fixed = TRUE)))
  for (p in seq_along(h$pools$text)) expect_identical(api$parse_control_pool_string(h$pools$text[p]), h$pools$values[[p]])
  expect_false(anyDuplicated(h$pools$text) > 0L)
  expect_false(anyDuplicated(meta$design$vertices$vertex_id) > 0L)
  expect_identical(meta$design, api$build_comparison_catalog_design(meta$catalog, inp$dq, inp$tpt))
  expect_identical(length(unique(meta$design$vertices$time[meta$design$vertices$role == "reference"])), 2L)
})

test_that("malformed qualifications, targets, references, duplicates and bounds refuse identically on the admitted route, the fallback and the standalone wrappers, before any design construction", {
  ctrl <- ch_ctrl()
  fx <- ch_fixtures()$rotating_nyt
  own <- ch_owners(fx$d, ctrl)
  meta0 <- ch_build(api, own, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post, build_design = FALSE)
  inp <- ch_design_inputs(meta0)
  cat <- meta0$catalog
  expect_identical(cat$spec$calendar, c(1, 2, 3, 4)); expect_identical(cat$spec$cohorts, c(3, 4, Inf))
  routes <- function(dq, tpt) list(
    reuse = ch_err(ch_build(api, own, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post,
      occupancy = meta0$occupancy_object, declared_qualified = dq, target_post_times = tpt, build_design = TRUE,
      metadata_source = meta0)),
    fallback = ch_err(ch_build(api, own, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post,
      occupancy = meta0$occupancy_object, declared_qualified = dq, target_post_times = tpt, build_design = TRUE)),
    standalone = ch_err(api$build_comparison_catalog_design(cat, dq, tpt)))
  bad_dq <- list(
    pre = list(dq = data.frame(g = 3, t = 1), cls = "cdid_catalog_qualification_error"),
    off_calendar = list(dq = data.frame(g = 3, t = 3.5), cls = "cdid_catalog_qualification_error"),
    not_a_target = list(dq = data.frame(g = 5, t = 4), cls = "cdid_catalog_qualification_error"),
    duplicate = list(dq = data.frame(g = c(3, 3), t = c(4, 4)), cls = "cdid_catalog_duplicate_error"),
    wrong_columns = list(dq = data.frame(g = 3, time = 4), cls = "cdid_catalog_qualification_error"),
    character = list(dq = data.frame(g = "3", t = "4", stringsAsFactors = FALSE), cls = "cdid_catalog_qualification_error"))
  for (nm in names(bad_dq)) {
    r <- routes(bad_dq[[nm]]$dq, inp$tpt)
    for (route in names(r)) expect_true(inherits(r[[route]], bad_dq[[nm]]$cls), info = paste(nm, route))
    ch_same_refusal(r$reuse, r$fallback, nm)
    ch_same_refusal(r$reuse, r$standalone, nm)
    expect_true(inherits(ch_err(api$comparison_catalog_design_links(cat, bad_dq[[nm]]$dq)), bad_dq[[nm]]$cls), info = nm)
    expect_true(inherits(ch_err(api$catalog_design_handoff(cat, bad_dq[[nm]]$dq)), bad_dq[[nm]]$cls), info = nm)
  }
  bad_tpt <- list(
    before_g = data.frame(g = 3, t = 2), off_calendar = data.frame(g = 3, t = 3.5), unknown_cohort = data.frame(g = 7, t = 4),
    at_reference = data.frame(g = 4, t = 3))
  for (nm in names(bad_tpt)) {
    r <- routes(inp$dq, bad_tpt[[nm]])
    for (route in names(r)) expect_true(inherits(r[[route]], "cdid_design_target_error"), info = paste(nm, route))
    ch_same_refusal(r$reuse, r$fallback, nm)
    ch_same_refusal(r$reuse, r$standalone, nm)
  }
  # Owner caps: the design charge is checked against the unchanged caps before
  # the handoff or the design core is entered (zero constructions).
  e <- ch_clone(c("catalog_design_handoff", "build_catalog_handoff_design", "chained_design_core"))
  meta_ok <- ch_build(api, own, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post,
    occupancy = meta0$occupancy_object, declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE)
  need <- meta_ok$diagnostics$mapping_accounting$retained_entries + meta_ok$diagnostics$catalog_reuse$owned_design_entries
  tight <- ch_ctrl(retained = need - 1, output = need - 1)
  r <- ch_err(ch_build(e, own, fx$targets, fx$cg, "conv:neg", ctrl = tight, requested_post = fx$requested_post,
    occupancy = meta0$occupancy_object, declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE))
  expect_true(inherits(r, "cdid_ipw_owner_bound_error"))
  for (nm in c("catalog_design_handoff", "build_catalog_handoff_design", "chained_design_core")) expect_identical(ch_n(e, nm), 0, info = nm)
  # Catalog scratch cap just below the complete design envelope: the catalog
  # bound refusal precedes any handoff construction.
  ad <- api$catalog_adapter_preflight(cat, inp$dq, inp$tpt, design = TRUE)
  low <- ch_bounds; low$max_scratch <- ad$scratch - 1
  r <- ch_err(ch_build(e, own, fx$targets, fx$cg, "conv:neg", bounds = low, requested_post = fx$requested_post,
    occupancy = meta0$occupancy_object, declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE))
  expect_true(inherits(r, "cdid_catalog_bound_error"))
  for (nm in c("catalog_design_handoff", "build_catalog_handoff_design", "chained_design_core")) expect_identical(ch_n(e, nm), 0, info = nm)
  small <- cat; small$resources$bounds$max_scratch <- 0
  expect_true(inherits(ch_err(api$build_comparison_catalog_design(small, inp$dq, inp$tpt)), "cdid_catalog_bound_error"))
  expect_true(inherits(ch_err(api$comparison_catalog_design_links(small, inp$dq)), "cdid_catalog_bound_error"))

  # Foreign or mutated owner/source: the admitted route keeps the original refusals.
  stale <- meta0
  stale$catalog$rows$t[1L] <- stale$catalog$rows$t[1L] + 0.5
  expect_true(inherits(ch_err(ch_build(api, own, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = stale)),
    "cdid_pipeline_metadata_integrity_error"))
  own_x <- ch_owners(fx$d, ctrl, with_x = TRUE)
  foreign <- ch_err(api$build_owner_comparison_catalog(panel_owner = own$po, targets = fx$targets, anticipation = 0,
    control_group = fx$cg, conventions = "conv:neg", requested_post = fx$requested_post, control = ctrl,
    catalog_bounds = ch_bounds, frame_owner = own_x$fo, declared_qualified = inp$dq, target_post_times = inp$tpt,
    build_design = TRUE, metadata_source = meta0))
  expect_true(inherits(foreign, "cdid_ipw_owner_mismatch_error"))
  expect_identical(conditionMessage(foreign), "Frame owner belongs to a different `cdid_ipw_panel_owner` snapshot.")
  other <- ch_owners(ch_fixtures()$balanced_nt$d, ctrl)
  expect_true(inherits(ch_err(ch_build(api, other, fx$targets, fx$cg, "conv:neg", requested_post = fx$requested_post,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE, metadata_source = meta0)),
    "cdid_pipeline_metadata_owner_mismatch_error"))
  # A mutated catalog handed to the standalone wrappers is refused by their preflight.
  mut <- cat; mut$rows$pool_id[1L] <- mut$pools$n_pools + 1L
  expect_true(inherits(ch_err(api$comparison_catalog_design_links(mut, inp$dq)), "cdid_catalog_input_error"))
  expect_true(inherits(ch_err(api$build_comparison_catalog_design(mut, inp$dq, inp$tpt)), "cdid_catalog_input_error"))
  expect_true(inherits(ch_err(api$comparison_catalog_design_links(list(rows = cat$rows), inp$dq)), "cdid_catalog_input_error"))
  expect_true(inherits(ch_err(api$catalog_design_handoff(list(rows = cat$rows), inp$dq)), "cdid_catalog_input_error"))
})

test_that("the private design core admits only a packet bound to its own links and still validates the typed pools semantically", {
  ctrl <- ch_ctrl()
  fx <- ch_fixtures()$rotating_nyt
  own <- ch_owners(fx$d, ctrl)
  meta0 <- ch_build(api, own, fx$targets, fx$cg, "conv:packet", build_design = FALSE)
  inp <- ch_design_inputs(meta0)
  cat <- meta0$catalog
  h <- api$catalog_design_handoff(cat, inp$dq)
  expect_gte(length(h$pools$values), 2L)
  ref <- api$build_comparison_catalog_design(cat, inp$dq, inp$tpt)
  expect_identical(api$build_catalog_handoff_design(h, inp$tpt), ref)
  internal <- "cdid_design_internal_error"
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(list(a = 1), inp$tpt)), internal))
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(h$links, inp$tpt)), internal))
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(h[c("links", "calendar")], inp$tpt)), internal))
  # Pools not bound to the actual column: permuted index, edited text, foreign pool.
  perm <- h; perm$pools$index <- rev(perm$pools$index)
  expect_false(identical(perm$pools$text[perm$pools$index], perm$links$control_pool))
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(perm, inp$tpt)), internal))
  edited <- h; edited$links$control_pool[1L] <- h$pools$text[h$pools$index[1L] %% length(h$pools$text) + 1L]
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(edited, inp$tpt)), internal))
  short <- h; short$pools$values <- short$pools$values[-1L]
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(short, inp$tpt)), internal))
  dbl <- h; dbl$pools$index <- as.double(dbl$pools$index)
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(dbl, inp$tpt)), internal))
  # Unchanged text must not authenticate different, otherwise legal values.
  # A never-treated singleton is clean at every date and belongs to this universe.
  k_diff <- which(vapply(h$pools$values, function(z) !identical(z, Inf), logical(1L)))[1L]
  expect_false(is.na(k_diff))
  wrong_values <- h; wrong_values$pools$values[[k_diff]] <- Inf
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(wrong_values, inp$tpt)), internal))
  # For semantic negatives, make the typed values and text agree first; these
  # are invalid cohort sets rather than malformed private binding packets.
  bind_pool <- function(packet, k) {
    packet$pools$text[k] <- paste0("{", paste(sprintf("%.17g", packet$pools$values[[k]]), collapse = ","), "}")
    packet$links$control_pool <- packet$pools$text[packet$pools$index]
    packet
  }
  k1 <- h$pools$index[1L]
  bad_g <- h; bad_g$pools$values[[k1]] <- sort(c(bad_g$pools$values[[k1]], h$links$g[1L]))
  bad_g <- bind_pool(bad_g, k1)
  e_g <- ch_err(api$build_catalog_handoff_design(bad_g, inp$tpt))
  expect_true(inherits(e_g, "cdid_design_qualification_error")); expect_false(inherits(e_g, internal))
  bad_u <- h; bad_u$pools$values[[k1]] <- sort(c(bad_u$pools$values[[k1]], 99))
  bad_u <- bind_pool(bad_u, k1)
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(bad_u, inp$tpt)), "cdid_design_cohort_error"))
  early <- h; early$pools$values[[k1]] <- sort(c(early$pools$values[[k1]], 3))
  early <- bind_pool(early, k1)
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(early, inp$tpt)), "cdid_design_qualification_error"))
  # Duplicated links and a false reference are caught by the one validation.
  dup <- h; dup$links <- rbind(h$links, h$links[1L, ]); dup$pools$index <- c(h$pools$index, h$pools$index[1L])
  rownames(dup$links) <- NULL
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(dup, inp$tpt)), "cdid_design_duplicate_link_error"))
  wrong_ref <- h; wrong_ref$links$ref[h$links$g == 3] <- 1
  expect_true(inherits(ch_err(api$build_catalog_handoff_design(wrong_ref, inp$tpt)), "cdid_design_reference_error"))
  # The same binding refusal through the validation body directly.
  cs <- api$validate_calendar_spec(h$calendar, h$cohorts, h$anticipation)
  l <- h$links; attr(l, "calendar_spec") <- cs
  expect_true(inherits(ch_err(api$design_links_validation(l, cs, pools = perm$pools)), internal))
  v <- api$design_links_validation(l, cs, pools = h$pools)
  expect_identical(v$links, api$validate_design_links(l, cs))
  expect_identical(v, api$design_links_validation(l, cs))
})

# ---------------------------------------------------------------------------
# 4. Standalone wrappers stay strict; one validation per constructor
# ---------------------------------------------------------------------------

test_that("standalone design wrappers keep their contracts: validation packet, adaptation, graph core and legacy no-calendar path", {
  ctrl <- ch_ctrl()
  fx <- ch_fixtures()$rotating_nyt_delta1
  own <- ch_owners(fx$d, ctrl)
  meta0 <- ch_build(api, own, fx$targets, fx$cg, "conv:wrap", anticipation = 1, build_design = FALSE)
  inp <- ch_design_inputs(meta0)
  cat <- meta0$catalog
  cal <- cat$spec$calendar; coh <- cat$spec$cohorts
  links <- api$comparison_catalog_design_links(cat, inp$dq)
  cs <- api$validate_calendar_spec(cal, coh, 1)
  l <- links; attr(l, "calendar_spec") <- cs
  v <- api$design_links_validation(l, cs)
  expect_identical(names(v), c("links", "calendar_spec", "s_cmp", "t_cmp"))
  expect_identical(v$calendar_spec, cs)
  expect_identical(api$validate_design_links(l, cs), v$links)
  expect_identical(api$validate_design_links(l), v$links)
  expect_type(v$s_cmp, "integer"); expect_length(v$s_cmp, nrow(links))
  expect_identical(v$s_cmp, vapply(seq_len(nrow(links)), function(i) ch_cmp(links$s[i], 1, links$g[i]), integer(1L)))
  expect_identical(v$t_cmp, vapply(seq_len(nrow(links)), function(i) ch_cmp(links$t[i], 1, links$g[i]), integer(1L)))
  a <- api$adapt_calendar_design_links(links, cal, coh, 1)
  av <- api$adapt_calendar_design_links_validation(links, cal, coh, 1)
  expect_identical(a, av$links)
  expect_identical(av, v)
  expect_identical(api$build_chained_design(a, inp$tpt), api$chained_design_core(av, inp$tpt))
  expect_identical(api$build_calendar_chained_design(links, cal, coh, 1, target_post_times = inp$tpt), api$build_chained_design(a, inp$tpt))
  expect_identical(formals(api$build_chained_design), formals(function(links, target_post_times = NULL) NULL))
  expect_identical(formals(api$validate_design_links), formals(function(links, calendar_spec = NULL) NULL))
  expect_identical(formals(api$build_calendar_chained_design), formals(function(links, calendar, cohorts, anticipation = 0,
    target_post_times = NULL, pt_assumption = "post", pt_start = NULL) NULL))
  # Anticipation roles are reached and agree with the independent role function.
  des <- api$build_calendar_chained_design(links, cal, coh, 1, target_post_times = inp$tpt)
  expect_true(any(des$vertices$role %in% c("anticipation", "anticipation_nuisance")))
  for (i in seq_len(nrow(links))) {
    expect_identical(des$vertices$role[des$links$start_vertex_index[i]],
      ref_cal$ref_endpoint_role(links$s[i], links$g[i], links$ref[i], 1, links$s_qualified[i]))
    expect_identical(des$vertices$role[des$links$end_vertex_index[i]],
      ref_cal$ref_endpoint_role(links$t[i], links$g[i], links$ref[i], 1, links$t_qualified[i]))
  }
  # Legacy no-calendar path: the same wrapper, no cutoff signs, unchanged output through the core.
  legacy <- data.frame(g = c(3, 3, 4, 4), s = c(1, 2, 1, 1), t = c(2, 3, 3, 4), ref = c(2, 2, 3, 3),
    functional_key = c("op:A", "op:A", "op:B1", "op:B2"), stringsAsFactors = FALSE)
  tgt <- data.frame(g = c(3, 3, 4), time = c(3, 4, 4))
  lv <- api$design_links_validation(legacy)
  expect_null(lv$s_cmp); expect_null(lv$t_cmp); expect_null(lv$calendar_spec)
  expect_identical(lv$links, api$validate_design_links(legacy))
  expect_identical(api$build_chained_design(legacy, tgt), api$chained_design_core(lv, tgt))
  expect_null(api$build_chained_design(legacy, tgt)$calendar_spec)
  expect_identical(api$build_chained_design(legacy, tgt)$parameters$role, c("pre", "post", "post", "pre", "pre", "post"))
  # Strict wrapper refusals are unchanged.
  undeclared <- data.frame(g = 6, s = 3, t = 7, functional_key = "opaque", stringsAsFactors = FALSE)
  expect_error(api$build_calendar_chained_design(undeclared, c(0, 3, 7, 9, 12, 14), c(6, 10, 14, Inf), 2), class = "cdid_design_qualification_error")
  bad_ref <- links; bad_ref$ref <- 4
  expect_error(api$build_calendar_chained_design(bad_ref, cal, coh, 1), class = "cdid_design_reference_error")
  expect_error(api$build_chained_design(rbind(legacy, legacy[1L, ]), tgt), class = "cdid_design_duplicate_link_error")
})

# ---------------------------------------------------------------------------
# 5. Disclosed constructor/preflight counts of the fused route versus the
#    generic bridge, and the helper groups each stage demands
# ---------------------------------------------------------------------------

test_that("the admitted route runs one adapter preflight, one qualification validation, one link validation and no pool parse; the generic bridge keeps its own counts", {
  ctrl <- ch_ctrl()
  fx <- ch_fixtures()$rotating_nyt
  own <- ch_owners(fx$d, ctrl)
  e <- ch_clone(ch_wrap_names)
  rec <- ch_stage_recorder(e)
  meta0 <- ch_build(e, own, fx$targets, fx$cg, "conv:count", requested_post = fx$requested_post, build_design = FALSE)
  s0 <- ch_snapshot(e)
  g0 <- unique(unlist(rec$groups))
  expect_true(length(rec$groups) > 0)
  expect_true(all(g0 %in% c("owners", "occupancy", "catalog")))
  expect_false(any(c("design", "views", "readers") %in% g0))
  expect_identical(s0$catalog_adapter_preflight, 0); expect_identical(s0$catalog_design_handoff, 0)
  expect_identical(s0$build_catalog_handoff_design, 0); expect_identical(s0$design_links_validation, 0)
  expect_identical(s0$build_comparison_catalog, 1); expect_identical(s0$derive_owner_sparse_occupancy, 1)
  expect_identical(s0$pipeline_metadata_resolve_catalog_mappings, 1)
  inp <- ch_design_inputs(meta0)
  rec$groups <- list()
  meta <- ch_build(e, own, fx$targets, fx$cg, "conv:count", requested_post = fx$requested_post,
    occupancy = meta0$occupancy_object, declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE,
    metadata_source = meta0)
  s1 <- ch_snapshot(e)
  g1 <- unique(unlist(rec$groups))
  expect_true("design" %in% g1)
  expect_false(any(c("views", "readers") %in% g1))
  d <- function(nm) ch_delta(s0, s1, nm)
  expect_identical(d("build_owner_comparison_catalog"), 1)
  expect_identical(d("catalog_adapter_preflight"), 1)
  expect_identical(d("catalog_adapter_dimensions"), 1)
  # One lower-bound mirror in the stage-1 preflight, one inside the single
  # adapter preflight, one links envelope on the same actual dimensions.
  expect_identical(d("catalog_design_envelope"), 3)
  expect_identical(d("catalog_design_handoff"), 1)
  expect_identical(d("build_catalog_handoff_design"), 1)
  expect_identical(d("design_links_validation"), 1)
  expect_identical(d("chained_design_core"), 1)
  expect_identical(d("normalize_target_post_times"), 1)
  expect_identical(d("calendar_reference_date"), 1)
  expect_identical(d("parse_control_pool_string"), 0)
  for (nm in c("comparison_catalog_design_links", "build_comparison_catalog_design", "adapt_calendar_design_links",
    "adapt_calendar_design_links_validation", "validate_design_links", "build_chained_design", "build_calendar_chained_design",
    "build_comparison_catalog", "pipeline_metadata_resolve_catalog_mappings", "derive_owner_sparse_occupancy")) {
    expect_identical(d(nm), 0, info = nm)
  }
  fused_cmp <- d("exact_cmp_sum")
  expect_identical(meta$design, api$build_comparison_catalog_design(meta0$catalog, inp$dq, inp$tpt))

  # The generic standalone bridge in the same clone: two preflights, the
  # adaptation plus its validation, one pool parse per interned pool.
  s2 <- ch_snapshot(e)
  gen <- e$build_comparison_catalog_design(meta0$catalog, inp$dq, inp$tpt)
  s3 <- ch_snapshot(e)
  d2 <- function(nm) ch_delta(s2, s3, nm)
  expect_identical(gen, meta$design)
  expect_identical(d2("catalog_adapter_preflight"), 2); expect_identical(d2("catalog_adapter_dimensions"), 2)
  expect_identical(d2("catalog_design_envelope"), 2)
  expect_identical(d2("comparison_catalog_design_links"), 1); expect_identical(d2("catalog_design_handoff"), 1)
  expect_identical(d2("build_catalog_handoff_design"), 0)
  expect_identical(d2("build_calendar_chained_design"), 1); expect_identical(d2("adapt_calendar_design_links_validation"), 1)
  expect_identical(d2("design_links_validation"), 1); expect_identical(d2("chained_design_core"), 1)
  expect_identical(d2("calendar_reference_date"), 2)
  expect_identical(d2("parse_control_pool_string"), as.double(meta0$catalog$pools$n_pools))
  expect_gt(d2("exact_cmp_sum"), fused_cmp)
  # The former two-pass composition validates twice.
  s4 <- ch_snapshot(e)
  two <- e$build_chained_design(e$adapt_calendar_design_links(meta$design_links, meta0$catalog$spec$calendar,
    meta0$catalog$spec$cohorts, 0), inp$tpt)
  s5 <- ch_snapshot(e)
  expect_identical(two, meta$design)
  expect_identical(ch_delta(s4, s5, "design_links_validation"), 2)
  expect_identical(ch_delta(s4, s5, "calendar_reference_date"), 3)
  expect_identical(ch_delta(s4, s5, "chained_design_core"), 1)

  # Links-only design request: one preflight (links envelope), the handoff, no design core.
  s6 <- ch_snapshot(e)
  only_links <- ch_build(e, own, fx$targets, fx$cg, "conv:count", requested_post = fx$requested_post,
    occupancy = meta0$occupancy_object, declared_qualified = inp$dq, build_design = TRUE, metadata_source = meta0)
  s7 <- ch_snapshot(e)
  expect_identical(ch_delta(s6, s7, "catalog_adapter_preflight"), 1)
  expect_identical(ch_delta(s6, s7, "catalog_design_envelope"), 2)
  expect_identical(ch_delta(s6, s7, "catalog_design_handoff"), 1)
  expect_identical(ch_delta(s6, s7, "build_catalog_handoff_design"), 0)
  expect_identical(ch_delta(s6, s7, "design_links_validation"), 0)
  expect_null(only_links$design)
  expect_identical(only_links$design_links, api$comparison_catalog_design_links(meta0$catalog, inp$dq))
  expect_identical(only_links$design_links, meta$design_links[names(only_links$design_links)] |> (function(x) {attr(x, "calendar_spec") <- NULL; x})())
  expect_identical(only_links$diagnostics$catalog_reuse$owned_design_entries, api$catalog_adapter_preflight(meta0$catalog, inp$dq)$scratch)
})

test_that("helper groups partition the exhaustive list; each entry demands only its stage; sandbox clones re-check live bindings and refuse missing or non-function helpers", {
  original_exhaustive <- c(
    "cdid_abort", "pair_count_safe", "pair_scalar_vector",
    "gmm_sexp_aliases", "gmm_scalar_numeric", "gmm_scalar_logical", "gmm_scalar_character",
    "logit_cert_count_output", "design_vector", "exact_cmp_sum", "validate_calendar_spec",
    "cdid_ipw_owner_abort", "ipw_owner_check_interrupt", "ipw_owner_seal", "ipw_owner_seal_ok",
    "ipw_owner_refuse_decorative_attrs", "ipw_owner_int_bound", "ipw_owner_safe_product",
    "ipw_owner_safe_sum", "ipw_owner_control", "validate_ipw_owner_control",
    "ipw_owner_match_sorted_double", "ipw_owner_upper_bound", "validate_ipw_panel_owner", "validate_ipw_frame_owner",
    "traverse_ipw_panel_pair_positions", "query_ipw_panel_pair_positions",
    "construct_ipw_nuisance_view", "construct_ipw_link_view", "ipw_owner_reset_trace",
    "catalog_abort", "validate_catalog_bounds", "validate_catalog_occupancy",
    "validate_catalog_requests", "catalog_size", "catalog_product",
    "build_comparison_catalog", "catalog_adapter_preflight", "catalog_adapter_dimensions",
    "catalog_design_envelope", "comparison_catalog_design_links", "build_comparison_catalog_design")
  groups <- api$pipeline_metadata_dependency_groups
  expect_type(groups, "list")
  expect_identical(names(groups), c("owners", "occupancy", "catalog", "design", "views", "readers"))
  all_fns <- unlist(groups, use.names = FALSE)
  expect_false(anyDuplicated(all_fns) > 0L)
  expect_setequal(all_fns, c(original_exhaustive, "catalog_design_handoff", "build_catalog_handoff_design"))
  for (fn in all_fns) expect_true(is.function(get0(fn, envir = api, mode = "function", inherits = FALSE)), info = fn)
  expect_identical(groups$views, c("construct_ipw_nuisance_view", "construct_ipw_link_view"))
  expect_identical(groups$readers, c("traverse_ipw_panel_pair_positions", "query_ipw_panel_pair_positions"))
  expect_true(all(c("catalog_adapter_preflight", "catalog_design_handoff", "build_catalog_handoff_design") %in% groups$design))
  expect_true(all(c("build_comparison_catalog", "exact_cmp_sum", "validate_calendar_spec") %in% groups$catalog))
  expect_true(all(c("validate_ipw_panel_owner", "validate_ipw_owner_control", "validate_catalog_bounds", "ipw_owner_seal") %in% groups$owners))
  expect_true(all(c("pair_count_safe", "ipw_owner_upper_bound", "validate_catalog_occupancy") %in% groups$occupancy))
  expect_identical(formals(api$pipeline_metadata_require_deps), NULL)
  expect_true(api$pipeline_metadata_require_deps())
  expect_true(api$pipeline_metadata_require_stage("owners"))
  expect_error(api$pipeline_metadata_require_stage("bogus"), class = "cdid_pipeline_metadata_dependency_error")
  expect_error(api$pipeline_metadata_require_stage(character()), class = "cdid_pipeline_metadata_dependency_error")
  expect_error(api$pipeline_metadata_require_stage(NULL), class = "cdid_pipeline_metadata_dependency_error")
  # The genuine locked namespace records each group once, under its own identity only.
  for (grp in names(groups)) expect_true(api$cdid_dependency_check_done(paste0("pipeline_metadata:", grp), api), info = grp)

  ctrl <- ch_ctrl()
  fx <- ch_fixtures()$rotating_nyt
  own <- ch_owners(fx$d, ctrl)
  meta_api <- ch_build(api, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE)
  inp <- ch_design_inputs(meta_api)
  dep <- "cdid_pipeline_metadata_dependency_error"

  # Sandbox: the genuine namespace is not in the lookup chain, so removals are real.
  s <- ch_clone(parent = parent.env(api))
  expect_false(isNamespace(s))
  expect_true(s$pipeline_metadata_require_deps())
  m <- ch_build(s, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE)
  expect_identical(m$catalog, meta_api$catalog)
  expect_identical(m$row_mapping, meta_api$row_mapping)
  expect_false(api$cdid_dependency_check_done("pipeline_metadata:owners", s))
  # Views removed: the exhaustive check and the view stage refuse; catalog and
  # design stages, which never call a view constructor, do not.
  rm("construct_ipw_link_view", envir = s)
  expect_error(s$pipeline_metadata_require_deps(), class = dep)
  expect_error(s$pipeline_metadata_require_stage("views"), class = dep)
  expect_true(s$pipeline_metadata_require_stage(c("owners", "occupancy", "catalog", "design")))
  m1 <- ch_build(s, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE)
  expect_identical(m1$catalog, meta_api$catalog)
  m2 <- ch_build(s, own, fx$targets, fx$cg, "conv:sandbox", occupancy = m1$occupancy_object, declared_qualified = inp$dq,
    target_post_times = inp$tpt, build_design = TRUE, metadata_source = m1)
  expect_identical(m2$design, api$build_comparison_catalog_design(meta_api$catalog, inp$dq, inp$tpt))
  expect_error(s$construct_pipeline_catalog_row_views(m1, 1L), class = dep)
  spec <- s$extract_pipeline_catalog_row_spec(m1, 1L)
  expect_identical(spec$row_id, 1L)
  expect_error(s$pipeline_metadata_row_views_from_spec(spec, own$fo, ctrl), class = dep)
  # The genuine namespace constructs the views for the same row.
  rv <- api$construct_pipeline_catalog_row_views(meta_api, 1L)
  expect_s3_class(rv$link_view, "cdid_ipw_link_view"); expect_identical(rv$endpoints_read, 0L)
  # Design helper removed: only a requested design refuses, at entry, before any owner work.
  s2 <- ch_clone(wrap = c("derive_owner_sparse_occupancy", "validate_ipw_panel_owner"), parent = parent.env(api))
  rm("catalog_design_handoff", envir = s2)
  expect_error(ch_build(s2, own, fx$targets, fx$cg, "conv:sandbox", occupancy = meta_api$occupancy_object,
    declared_qualified = inp$dq, target_post_times = inp$tpt, build_design = TRUE), class = dep)
  expect_identical(ch_n(s2, "validate_ipw_panel_owner"), 0)
  m3 <- ch_build(s2, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE)
  expect_identical(m3$catalog, meta_api$catalog)
  expect_identical(ch_n(s2, "derive_owner_sparse_occupancy"), 1)
  # Catalog helper removed: the catalog stage refuses at entry.
  s3 <- ch_clone(wrap = "derive_owner_sparse_occupancy", parent = parent.env(api))
  rm("build_comparison_catalog", envir = s3)
  expect_error(ch_build(s3, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE), class = dep)
  expect_identical(ch_n(s3, "derive_owner_sparse_occupancy"), 0)
  expect_error(s3$preflight_owner_comparison_catalog(own$po, fx$targets, 0, fx$cg, "conv:sandbox", control = ctrl,
    catalog_bounds = ch_bounds), class = dep)
  # Owner-only entries keep working without the catalog helper.
  expect_s3_class(s3$pipeline_metadata_control(ctrl, ch_bounds), "cdid_pipeline_metadata_control")
  expect_identical(s3$derive_owner_sparse_occupancy(own$po, ctrl)$occupancy, meta_api$occupancy)
  # Function-type check (never mode = "any"): a non-function binding is not a helper.
  s4 <- ch_clone(parent = parent.env(api))
  assign("exact_cmp_sum", "not a function", envir = s4)
  expect_error(s4$pipeline_metadata_require_stage("catalog"), class = dep)
  expect_error(ch_build(s4, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE), class = dep)
  expect_true(s4$pipeline_metadata_require_stage("owners"))
  # Live binding: a cloned environment is never memoized, so a later removal is seen.
  s5 <- ch_clone(parent = parent.env(api))
  m5 <- ch_build(s5, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE)
  expect_identical(m5$catalog, meta_api$catalog)
  rm("validate_catalog_requests", envir = s5)
  expect_error(ch_build(s5, own, fx$targets, fx$cg, "conv:sandbox", build_design = FALSE), class = dep)
  expect_false(api$cdid_dependency_check_done("pipeline_metadata:catalog", s5))
  expect_true(api$cdid_dependency_check_done("pipeline_metadata:catalog", api))
})

# ---------------------------------------------------------------------------
# 6. Public six-method family: genuine fused route versus a disclosed clone
#    forcing the generic bridge, at the existing tolerance
# ---------------------------------------------------------------------------

test_that("public IPW/OR/DR with and without covariates: the fused handoff and the generic bridge agree on IDs, points, IF, covariance and tests", {
  d <- ch_public_panel(41L, 8L)
  # The public entry must be a namespace closure for the clones to instrument it.
  expect_true(is.function(api$att_gt)); expect_identical(environment(api$att_gt), api)
  generic <- ch_clone(wrap = c("build_catalog_handoff_design", "catalog_design_handoff", "comparison_catalog_design_links",
    "build_comparison_catalog_design", "build_owner_comparison_catalog", "build_comparison_catalog"),
    overrides = list(build_catalog_handoff_design = ch_generic_bridge))
  counted <- ch_clone(wrap = c("build_catalog_handoff_design", "catalog_design_handoff", "comparison_catalog_design_links",
    "build_comparison_catalog_design", "build_owner_comparison_catalog", "build_comparison_catalog", "catalog_adapter_preflight"))
  run <- function(ns, est, form) ns$att_gt(yname = "y", tname = "period", idname = "id", gname = "G", data = d,
    est_method = est, xformla = if (form == "1") ~1 else ~x, control_group = "nevertreated", bstrap = FALSE, cband = FALSE,
    base_period = "varying", mode = "chained", gmm_weighting = "optimal", gmm_covariance = "tikhonov",
    pt_assumption = "post", control = list(chunk_rows = 7))
  for (form in c("1", "x")) for (est in c("ipw", "reg", "dr")) {
    info <- paste(est, form)
    before <- ch_snapshot(counted)
    mp <- run(api, est, form)
    mp_c <- run(counted, est, form)
    after <- ch_snapshot(counted)
    mp_g <- run(generic, est, form)
    expect_s3_class(mp, "cdid_MP")
    expect_true(is.matrix(mp$inffunc), info = info)
    expect_identical(nrow(mp$inffunc), 41L, info = info)
    expect_gt(sum(is.finite(mp$att) & is.finite(mp$se)), 0)
    for (f in c("group", "t")) expect_identical(mp[[f]], mp_g[[f]], info = paste(info, f))
    for (f in c("att", "se", "V_analytical", "inffunc", "W", "Wpval")) {
      expect_equal(mp[[f]], mp_g[[f]], tolerance = 1e-12, info = paste(info, f))
      expect_equal(mp[[f]], mp_c[[f]], tolerance = 1e-12, info = paste(info, f))
    }
    # Genuine route: two catalog calls, one catalog, one handoff, one fused design, no generic bridge.
    expect_identical(ch_delta(before, after, "build_owner_comparison_catalog"), 2, info = info)
    expect_identical(ch_delta(before, after, "build_comparison_catalog"), 1, info = info)
    expect_identical(ch_delta(before, after, "catalog_design_handoff"), 1, info = info)
    expect_identical(ch_delta(before, after, "build_catalog_handoff_design"), 1, info = info)
    expect_identical(ch_delta(before, after, "catalog_adapter_preflight"), 1, info = info)
    expect_identical(ch_delta(before, after, "comparison_catalog_design_links"), 0, info = info)
    expect_identical(ch_delta(before, after, "build_comparison_catalog_design"), 0, info = info)
    expect_identical(ch_n(generic, "build_catalog_handoff_design"), ch_n(generic, "catalog_design_handoff"), info = info)
  }
})
