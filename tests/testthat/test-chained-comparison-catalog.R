# Strict tests for the private CP02 approved comparison catalog (D10 retention,
# D13 fixed-horizon reuse); proofs in docs/theory/cp02_approved_comparison_catalog.tex.
# These tests use the loaded cdid namespace only. They never source package
# files or dyn.load code; a missing or foreign namespace is a hard error.

if (!("cdid" %in% loadedNamespaces())) {
  stop("Load the candidate cdid namespace before running this test file.", call. = FALSE)
}
local({
  ns <- asNamespace("cdid")
  desc <- file.path(getNamespaceInfo(ns, "path"), "DESCRIPTION")
  if (!identical(unname(getNamespaceName(ns)), "cdid") || !file.exists(desc) ||
      !identical(unname(read.dcf(desc, fields = "Package")[1L, 1L]), "cdid")) {
    stop("Foreign `cdid` namespace; refusing to test it.", call. = FALSE)
  }
  for (f in c("build_comparison_catalog", "comparison_catalog_design_links", "build_comparison_catalog_design",
              "validate_calendar_spec", "exact_cmp_sum", "calendar_reference_date", "eligible_control_cohorts",
              "build_calendar_chained_design", "materialize_chained_design", "check_level_estimable")) {
    fn <- get0(f, envir = ns, mode = "function", inherits = FALSE)
    if (is.null(fn) || !identical(environment(fn), ns)) {
      stop("`", f, "` is not bound in the loaded cdid namespace; no source/dyn.load fallback.", call. = FALSE)
    }
  }
})

catalog_reference_sha <- c(
  cp02_calendar_design_reference.R = "ed66dddd0b84ebe4f0e4bb1ef64eb0866d92eab0f1be2b697762279dcdda3031",
  cp02_comparison_catalog_reference.R = "2307f19615bc3bc98d664a072a036c5a7c9828b5b93782a53c62cbcbd6b9ac58"
)

load_catalog_reference <- function() {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to bind packaged references.", call. = FALSE)
  env <- new.env(parent = baseenv())
  for (file in names(catalog_reference_sha)) {
    path <- testthat::test_path("reference", file)
    if (!file.exists(path)) stop("missing packaged reference: ", file, call. = FALSE)
    expect_identical(digest::digest(file = path, algo = "sha256"), catalog_reference_sha[[file]])
    sys.source(path, envir = env)
  }
  env
}

# Test-local explicit caps (not package defaults).
cat_caps <- list(max_candidates = 1e6, max_rows = 1e6, max_pool_entries = 1e7, max_scratch = 1e8, max_work = 1e8)

cat_build <- function(fx, control_group = fx$control_group, requested_post = NULL, bounds = cat_caps,
                      conventions = "conv:test") {
  cdid:::build_comparison_catalog(fx$calendar, fx$cohorts, fx$targets, fx$anticipation, fx$occupancy,
    control_group, conventions, requested_post, bounds)
}

cat_grid_occupancy <- function(calendar, cohorts, rule) {
  out <- list()
  for (si in seq_along(calendar)) for (ti in seq_along(calendar)) if (si < ti) {
    for (ki in seq_along(cohorts)) {
      out[[length(out) + 1L]] <- data.frame(cohort = cohorts[ki], s = calendar[si], t = calendar[ti],
        count = rule(si, ti, ki))
    }
  }
  do.call(rbind, out)
}

cat_witness <- function() {
  # plan/evidence/cp02_nyt_horizon_choice (script sha256 bf80c18c...): six
  # equal-mass units observed exactly twice on pairs 12 and 14, cohorts 3/4/5.
  list(calendar = 1:5, cohorts = c(3, 4, 5), targets = 3, anticipation = 0, control_group = "notyettreated",
       occupancy = data.frame(cohort = rep(c(3, 4, 5), 2L), s = 1, t = rep(c(2, 4), each = 3L), count = 1))
}

cat_decode <- function(cat) {
  pl <- cat$pools; coh <- cat$spec$cohorts; cal <- cat$spec$calendar
  lapply(seq_len(nrow(cat$rows)), function(r) {
    x <- cat$rows[r, ]
    K <- coh[pl$entries[seq.int(pl$start[x$pool_id], length.out = pl$size[x$pool_id])]]
    hz <- if (x$n_horizons > 0L) cal[seq.int(x$horizon_first_index, x$horizon_last_index)] else numeric()
    list(g = x$g, s = x$s, t = x$t, ref = cal[x$ref_index], K = K, original = x$original,
         horizons = hz, n_horizons = x$n_horizons)
  })
}

# Full independent comparison with the Cartesian reference.
cat_check_reference <- function(cat, fx, ref, control_group = fx$control_group) {
  cmp <- ref$ref_exact_cmp_sum
  want <- ref$ref_catalog_enumerate(fx$calendar, fx$cohorts, fx$targets, fx$anticipation, fx$occupancy,
                                    control_group, cmp)
  got <- cat_decode(cat)
  kg <- vapply(got, function(x) ref$ref_cat_key(x$g, x$s, x$t, x$K), character(1L))
  kw <- vapply(want$rows, function(x) ref$ref_cat_key(x$g, x$s, x$t, x$K), character(1L))
  expect_false(anyDuplicated(kg) > 0L)
  expect_identical(sort(kg, method = "radix"), sort(kw, method = "radix"))
  for (i in seq_along(got)) {
    w <- want$rows[[match(kg[i], kw)]]
    expect_identical(got[[i]]$K, w$K)
    expect_identical(got[[i]]$ref, w$ref)
    expect_identical(got[[i]]$original, w$original)
    expect_identical(got[[i]]$horizons, w$horizons)
    expect_identical(got[[i]]$n_horizons, length(w$horizons))
    r <- cat$rows[i, ]
    expect_identical(r$s_role, ref$ref_cat_role(w$s, w$g, w$ref, w$K, fx$anticipation, cmp))
    expect_identical(r$t_role, ref$ref_cat_role(w$t, w$g, w$ref, w$K, fx$anticipation, cmp))
  }
  res <- cat$resources
  expect_identical(res$candidates, as.double(want$n_candidates))
  expect_identical(res$nonempty_candidates, as.integer(want$n_nonempty))
  expect_identical(res$rows, length(want$rows))
  expect_identical(res$pool_entries_preintern, as.double(sum(vapply(want$rows, function(x) length(x$K), numeric(1L)))))
  # Operator identity: same (g, K) under the one declared convention.
  okey <- vapply(got, function(x) paste(sprintf("%a", x$g), paste(sprintf("%a", x$K), collapse = ",")), character(1L))
  expect_identical(match(okey, unique(okey)), match(cat$rows$operator_id, unique(cat$rows$operator_id)))
  # Coordinate identity: shared per (g, u); private per (g, K, u).
  ckey <- function(u, role, x) {
    if (role == "reference") paste("ref", sprintf("%a", x$g))
    else if (grepl("_private$", role)) paste("priv", sprintf("%a", x$g), sprintf("%a", u), paste(sprintf("%a", x$K), collapse = ","))
    else paste("shared", sprintf("%a", x$g), sprintf("%a", u))
  }
  ek <- c(vapply(seq_along(got), function(i) ckey(got[[i]]$s, cat$rows$s_role[i], got[[i]]), character(1L)),
          vapply(seq_along(got), function(i) ckey(got[[i]]$t, cat$rows$t_role[i], got[[i]]), character(1L)))
  eid <- c(cat$rows$s_coordinate_id, cat$rows$t_coordinate_id)
  expect_identical(match(ek, unique(ek)), match(eid, unique(eid)))
  # References agree with the accepted helper and the independent comparator.
  for (j in seq_len(nrow(cat$targets))) {
    g <- cat$targets$g[j]
    rr <- ref$ref_cat_reference(g, sort(as.double(fx$calendar)), fx$anticipation, cmp)
    expect_identical(cat$targets$ref[j], rr)
    if (is.na(rr)) {
      expect_error(cdid:::calendar_reference_date(g, fx$calendar, fx$anticipation), class = "cdid_design_reference_error")
    } else {
      expect_identical(cdid:::calendar_reference_date(g, fx$calendar, fx$anticipation), rr)
    }
  }
  nT <- length(fx$calendar); nC <- length(fx$cohorts)
  # Independent conservative visit envelope, replacing the preserved failed
  # assumption that only four operations per candidate are ever performed.
  H <- nT + nC + nrow(fx$occupancy) + length(fx$targets) + res$requested_rows + 1
  D <- want$n_candidates + res$target_pairs + 1
  B <- as.double(nchar(cat$conventions, type = "bytes"))
  expect_identical(res$work_bound, 64 * H * (1 + ceiling(log2(H + 1))) + 16 * B +
    128 * D * (nC + 1 + ceiling(log2(D + 1))) + 16 * B * want$n_candidates)
  expect_lte(sum(res$work_performed), res$work_bound)
  invisible(want)
}

cat_fixtures <- function() {
  list(
    regular_nt_inf = list(calendar = 1:6, cohorts = c(2, 3, 4, 6, Inf), targets = c(2, 3, 4), anticipation = 0,
      control_group = "notyettreated",
      occupancy = cat_grid_occupancy(1:6, c(2, 3, 4, 6, Inf), function(si, ti, ki) (si * 7 + ti * 3 + ki) %% 4)),
    regular_delta1 = list(calendar = 1:6, cohorts = c(2, 3, 4, 6, Inf), targets = c(2, 3, 4), anticipation = 1,
      control_group = "notyettreated",
      occupancy = cat_grid_occupancy(1:6, c(2, 3, 4, 6, Inf), function(si, ti, ki) (si + ti * 5 + ki * 2) %% 3)),
    irregular_offgrid_d0 = list(calendar = c(2000, 2002, 2005, 2009, 2011), cohorts = c(2003, 2007, 2009, 2010, Inf),
      targets = c(2007, 2009, 2010), anticipation = 0, control_group = "notyettreated",
      occupancy = cat_grid_occupancy(c(2000, 2002, 2005, 2009, 2011), c(2003, 2007, 2009, 2010, Inf),
                                     function(si, ti, ki) (si + 2 * ti + 3 * ki) %% 3)),
    irregular_offgrid_d4 = list(calendar = c(2000, 2002, 2005, 2009, 2011), cohorts = c(2003, 2007, 2009, 2010, Inf),
      targets = c(2003, 2007, 2009, 2010), anticipation = 4, control_group = "notyettreated",
      occupancy = cat_grid_occupancy(c(2000, 2002, 2005, 2009, 2011), c(2003, 2007, 2009, 2010, Inf),
                                     function(si, ti, ki) (si * ki + ti) %% 2)),
    no_nt = list(calendar = c(-3, 0, 1.5, 4, 8), cohorts = c(0.5, 1.5, 4, 9), targets = c(0.5, 1.5, 4),
      anticipation = 2, control_group = "notyettreated",
      occupancy = cat_grid_occupancy(c(-3, 0, 1.5, 4, 8), c(0.5, 1.5, 4, 9), function(si, ti, ki) (ti - si + ki) %% 3))
  )
}

test_that("the approved five-period exact bridge is recovered only with operator-specific early coordinates", {
  ref <- load_catalog_reference()
  fx <- cat_witness()
  cat <- cat_build(fx, requested_post = data.frame(g = 3, t = c(3, 4, 5)))
  cat_check_reference(cat, fx, ref)
  dec <- cat_decode(cat)
  expect_identical(lapply(dec, `[[`, "K"), list(c(4, 5), 5, 5))
  expect_identical(cat$rows$s, c(1, 1, 1)); expect_identical(cat$rows$t, c(2, 2, 4))
  expect_identical(cat$rows$original, c(TRUE, FALSE, TRUE))
  expect_identical(lapply(dec, `[[`, "horizons"), list(3, 4, 4))
  expect_identical(cat$rows$s_role, rep("pre_private", 3L))
  expect_identical(cat$rows$t_role, c("reference", "reference", "post_conditional"))
  expect_identical(cat$rows$operator_id, c(1L, 2L, 2L))
  # Pair 12 under {4,5} and under {5} keep distinct private date-1 coordinates.
  expect_identical(length(unique(cat$rows$s_coordinate_id)), 2L)
  expect_identical(cat$rows$s_coordinate_id[2L], cat$rows$s_coordinate_id[3L])
  expect_identical(unlist(cat$resources[c("candidates", "pool_entries_preintern")]),
                   c(candidates = 7, pool_entries_preintern = 4))
  expect_identical(unlist(cat$resources[c("nonempty_candidates", "rows", "pool_entries_interned", "pools", "operators", "coordinates")]),
                   c(nonempty_candidates = 5L, rows = 3L, pool_entries_interned = 3L, pools = 2L, operators = 2L, coordinates = 4L))
  expect_identical(cat$requested_supported$t, 4)
  expect_identical(cat$requested_unsupported$t, c(3, 5))
  expect_identical(cat$requested_unsupported$reason, rep("no_target_pair_occupancy", 2L))

  # Witness exact means for rows (12,{4,5}), (12,{5}), (14,{5}).
  mu <- c(2, 3, 4)
  des <- cdid:::materialize_chained_design(
    cdid:::build_comparison_catalog_design(cat, data.frame(g = 3, t = 4), target_post_times = 4), 1000)
  expect_identical(des$parameters$role, c("pre", "pre", "post"))
  W <- unname(des$W)
  expect_identical(W, matrix(c(-1L, 0L, 0L, 0L, -1L, 0L, 0L, -1L, 1L), 3L, byrow = TRUE))
  expect_identical(ref$ref_rowspace_check(W, numeric(3L))$rank, 3L)
  expect_identical(drop(W %*% c(-2, -3, 1)), mu)
  cert <- cdid:::check_level_estimable(des, "ATT(g=3,t=4)")
  expect_true(cert$estimable)
  expect_identical(sum(cert$witness_links$sign * mu[cert$witness_links$link_index]), 1)

  # Original-only rows (largest per-pair pools): rank 2, ATT(3,4) not estimable.
  links <- cdid:::comparison_catalog_design_links(cat, data.frame(g = 3, t = 4))
  orig <- cdid:::build_calendar_chained_design(links[cat$rows$original, ], 1:5, c(3, 4, 5), 0, 4)
  expect_identical(orig$rank_diagnostics$rank, 2L)
  expect_false(orig$post_targets$level_estimable)

  # Wrong-operator discriminator. Declared pools make a merged early operator a
  # refused contradiction; hiding them, the 12 rows collide as duplicate links
  # with means 2 and 3 (an inconsistent false restriction), and the false merger
  # of the two original operators reports ATT(3,4) = 4 - 2 = 2 instead of 1.
  merged <- links; merged$functional_key <- "merged_operator"
  expect_error(cdid:::build_calendar_chained_design(merged, 1:5, c(3, 4, 5), 0, 4),
               class = "cdid_design_profile_merge_error")
  merged$control_pool <- NULL
  expect_error(cdid:::build_calendar_chained_design(merged, 1:5, c(3, 4, 5), 0, 4),
               class = "cdid_design_duplicate_link_error")
  false_W <- rbind(c(-1, 0), c(-1, 0), c(-1, 1))
  expect_gt(max(abs(false_W %*% qr.solve(false_W, mu) - mu)), 0.4)
  mdes <- cdid:::build_calendar_chained_design(merged[cat$rows$original, ], 1:5, c(3, 4, 5), 0, 4)
  expect_identical(mdes$rank_diagnostics$rank, 2L)
  mw <- cdid:::check_level_estimable(mdes, "ATT(g=3,t=4)")$witness_links
  expect_identical(sum(mw$sign * mu[cat$rows$original][mw$link_index]), 2)
})

test_that("independent Cartesian enumeration matches exactly on regular, sparse, irregular, anticipation and no-NT fixtures", {
  ref <- load_catalog_reference()
  for (nm in names(cat_fixtures())) {
    fx <- cat_fixtures()[[nm]]
    for (cg in c("notyettreated", "nevertreated")) {
      cat <- cat_build(fx, control_group = cg)
      want <- cat_check_reference(cat, fx, ref, cg)
      if (cg == "nevertreated") expect_true(all(cat$rows$original & cat$rows$n_horizons == 0L))
      if (!any(is.infinite(fx$cohorts)) && cg == "nevertreated") expect_identical(nrow(cat$rows), 0L)
    }
  }
  # The approved family produces only conditional causal endpoints and private pre coordinates.
  cat <- cat_build(cat_fixtures()$irregular_offgrid_d4)
  expect_false(any(grepl("^(anticipation|post)_private$", c(cat$rows$s_role, cat$rows$t_role))))
  expect_true(all(c("pre_private", "reference", "anticipation_conditional", "post_conditional") %in%
                  c(cat$rows$s_role, cat$rows$t_role)))
})

test_that("every nonempty original pool is retained; identical pools are deduplicated with full provenance", {
  ref <- load_catalog_reference()
  fx <- cat_fixtures()$regular_nt_inf
  cat <- cat_build(fx)
  want <- cat_check_reference(cat, fx, ref)
  tp <- cat$target_pairs
  # Original NYT pool for every referenced supported target pair.
  for (i in seq_len(nrow(tp))) {
    occ <- fx$occupancy
    Ce <- sort(occ$cohort[occ$s == tp$s[i] & occ$t == tp$t[i] & occ$count > 0])
    K <- Ce[is.infinite(Ce) | (Ce > tp$g[i] & Ce > tp$t[i])]
    hit <- which(cat$rows$g == tp$g[i] & cat$rows$s == tp$s[i] & cat$rows$t == tp$t[i] & cat$rows$original)
    if (length(K)) {
      expect_identical(tp$original_status[i], "retained")
      expect_length(hit, 1L)
      expect_identical(cat_decode(cat)[[hit]]$K, K)
    } else {
      expect_identical(tp$original_status[i], "empty_pool")
      expect_length(hit, 0L)
    }
  }
  # For t >= g the original pool equals the h = t reevaluation: one row, both provenances.
  post <- cat$rows$t >= cat$rows$g & cat$rows$original
  expect_true(any(post))
  expect_true(all(cat$rows$horizon_first_index[post] == cat$rows$t_index[post]))
  # Horizon provenance is a contiguous interval, the row key is unique and rows per
  # target pair never exceed the suffix bound |F| + 1{Inf}.
  expect_true(all(is.na(cat$rows$horizon_first_index) |
    cat$rows$horizon_last_index - cat$rows$horizon_first_index + 1L == cat$rows$n_horizons))
  key <- paste(cat$rows$target_index, cat$rows$s_index, cat$rows$t_index, cat$rows$pool_id)
  expect_false(anyDuplicated(key) > 0L)
  expect_identical(sum(tp$n_rows), nrow(cat$rows))
  expect_true(all(tp$n_rows <= tp$n_other_supplied_cohorts))
  # Equal pools across pairs share one operator; pools are interned once.
  expect_identical(nrow(cat$operators), length(unique(paste(cat$rows$target_index, cat$rows$pool_id))))
  expect_identical(cat$pools$n_pools, length(unique(cat$rows$pool_id)))
  expect_identical(cat$resources$pool_entries_interned, sum(cat$pools$size))
  expect_lt(cat$resources$pool_entries_interned, cat$resources$pool_entries_preintern)
  # Different declared conventions give different operator keys but identical rows.
  cat2 <- cat_build(fx, conventions = "conv:other")
  expect_identical(cat2$rows, cat$rows)
  expect_false(any(cat2$operators$functional_key %in% cat$operators$functional_key))
})

test_that("sparse, zero, missing-Inf and earlier-cohort controls are handled exactly and reported honestly", {
  ref <- load_catalog_reference()
  # Occupancy shrinkage: cohort 5 absent at pair 12 gives K={4} there, distinct from {5}.
  fx <- cat_witness()
  fx$occupancy$count[fx$occupancy$cohort == 5 & fx$occupancy$t == 2] <- 0
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref)
  expect_identical(lapply(cat_decode(cat), `[[`, "K"), list(4, 5))
  expect_false(cat$rows$s_coordinate_id[1L] == cat$rows$s_coordinate_id[2L])
  # Equal shrunken pools on different pairs share the operator and its private coordinate.
  fx$occupancy$count <- 1
  fx$occupancy$count[fx$occupancy$cohort == 4 & fx$occupancy$t == 2] <- 0
  cat <- cat_build(fx)
  expect_identical(lapply(cat_decode(cat), `[[`, "K"), list(5, 5))
  expect_identical(cat$rows$operator_id, c(1L, 1L))
  expect_identical(cat$rows$s_coordinate_id[1L], cat$rows$s_coordinate_id[2L])
  expect_identical(cat$rows$original, c(TRUE, TRUE))
  expect_identical(cat$rows$n_horizons, c(2L, 1L))

  # Earlier-treated clean cohorts are not added even when the general helper allows them.
  fx <- list(calendar = 1:6, cohorts = c(4, 5, 6), targets = 5, anticipation = 0, control_group = "notyettreated",
             occupancy = data.frame(cohort = c(4, 5, 6, 4, 5, 6), s = c(1, 1, 1, 2, 2, 2), t = c(2, 2, 2, 5, 5, 5), count = 3))
  expect_identical(cdid:::eligible_control_cohorts(5, 2, fx$cohorts, 0), c(4, 6))
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref)
  expect_identical(lapply(cat_decode(cat), `[[`, "K"), list(6, 6))
  expect_false(any(cat$pools$entries == 1L))

  # No never-treated cohort under NT: no rows, honest empty-pool and request diagnostics.
  fx <- cat_witness()
  fx$occupancy <- rbind(fx$occupancy, data.frame(cohort = c(3, 4, 5), s = 1, t = 5, count = c(2, 1, 1)))
  nt <- cat_build(fx, control_group = "nevertreated", requested_post = data.frame(g = 3, t = c(4, 5)))
  expect_identical(nrow(nt$rows), 0L)
  expect_identical(nt$target_pairs$original_status, rep("empty_pool", 3L))
  expect_identical(nt$requested_unsupported$reason, rep("no_nonempty_control_pool", 2L))
  expect_error(cdid:::comparison_catalog_design_links(nt, NULL), class = "cdid_catalog_input_error")
  # Under NYT the 15 pair has no clean later cohort and no NT: reported, not dropped.
  nyt <- cat_build(fx, requested_post = data.frame(g = 3, t = c(4, 5)))
  expect_identical(nyt$target_pairs$original_status[nyt$target_pairs$t == 5], "empty_pool")
  expect_identical(nyt$requested_unsupported$t, 5)
  expect_identical(nyt$requested_unsupported$reason, "no_nonempty_control_pool")

  # A never-treated cohort with zero counts is not support; with positive counts it is.
  fx <- list(calendar = 1:4, cohorts = c(2, Inf), targets = 2, anticipation = 0, control_group = "nevertreated",
             occupancy = data.frame(cohort = c(2, Inf, 2, Inf), s = c(1, 1, 2, 2), t = c(3, 3, 4, 4), count = c(5, 0, 1, 7)))
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref)
  expect_identical(cat$rows$s, 2); expect_identical(cat$pools$entries, 2L)
  expect_identical(cat$target_pairs$original_status, c("empty_pool", "retained"))
})

test_that("missing references and requested targets are diagnosed explicitly, never invented or dropped", {
  fx <- list(calendar = 1:5, cohorts = c(1, 3, 4, Inf), targets = c(1, 3), anticipation = 0, control_group = "notyettreated",
             occupancy = data.frame(cohort = c(1, 3, 4, Inf, 1, Inf), s = c(1, 1, 1, 1, 2, 2), t = c(2, 2, 2, 2, 3, 3), count = 1))
  cat <- cat_build(fx, requested_post = data.frame(g = c(1, 1, 3, 3), t = c(2, 3, 4, 5)))
  expect_identical(cat$targets$status, c("no_admissible_reference", "reference_available"))
  expect_identical(cat$targets$ref, c(NA_real_, 2))
  expect_identical(cat$targets$n_candidates, c(0, 4))
  expect_identical(cat$target_pairs$original_status[cat$target_pairs$g == 1], rep("no_admissible_reference", 2L))
  expect_true(all(cat$rows$g == 3))
  expect_identical(cat$requested_unsupported$reason,
                   c("no_admissible_reference", "no_admissible_reference", "no_target_pair_occupancy", "no_target_pair_occupancy"))
  expect_identical(nrow(cat$requested_supported), 0L)
  expect_error(cdid:::calendar_reference_date(1, fx$calendar, 0), class = "cdid_design_reference_error")
  # The reference is the actual full-calendar date, never an ordinal g - delta - 1.
  cal <- c(2000, 2002, 2005, 2009)
  fx2 <- list(calendar = cal, cohorts = c(2009, Inf), targets = 2009, anticipation = 1, control_group = "notyettreated",
              occupancy = data.frame(cohort = c(2009, Inf), s = 2002, t = 2009, count = 1))
  c1 <- cat_build(fx2)
  expect_identical(c1$targets$ref, 2005)
  expect_identical(c(c1$rows$s_role, c1$rows$t_role), c("pre_private", "post_conditional"))
  fx2$anticipation <- 4
  c4 <- cat_build(fx2)
  expect_identical(c4$targets$ref, 2002)
  expect_identical(c(c4$rows$s_role, c4$rows$t_role), c("reference", "post_conditional"))
  fx2$anticipation <- 9
  c9 <- cat_build(fx2, requested_post = data.frame(g = 2009, t = 2009))
  expect_identical(c9$targets$status, "no_admissible_reference")
  expect_identical(c9$requested_unsupported$reason, "no_admissible_reference")
  expect_identical(nrow(c9$rows), 0L)
})

test_that("extreme stored dates use exact cutoffs, distinct stored values stay distinct, and overflow refuses", {
  ref <- load_catalog_reference()
  big <- 2^53 + 2; big4 <- 2^53 + 4
  # Exactly big + 1 < big4 although the rounded sum equals big4.
  fx <- list(calendar = c(2^53 - 2, 2^53, big, big4), cohorts = c(big, big4, Inf), targets = big, anticipation = 1,
             control_group = "notyettreated",
             occupancy = data.frame(cohort = c(big, big4, Inf), s = 2^53 - 2, t = big, count = 1))
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref)
  expect_identical(cat_decode(cat)[[1L]]$K, c(big4, Inf))
  expect_identical(cat$targets$ref, 2^53)
  # Declared design links carry the exact stored pool (round-tripped text).
  des <- cdid:::build_comparison_catalog_design(cat, data.frame(g = big, t = big), target_post_times = big)
  expect_identical(des$links$control_pool[1L], "{9007199254740996,Inf}")
  # The only pair starts before the reference: a private pre coordinate, no level.
  expect_identical(des$post_targets$status, "disconnected_component")
  # Cancellation: g - delta = 0 exactly, reference -1.
  fx <- list(calendar = c(-2, -1, 0, big4), cohorts = c(big4, Inf), targets = big4, anticipation = big4,
             control_group = "nevertreated",
             occupancy = cat_grid_occupancy(c(-2, -1, 0, big4), c(big4, Inf), function(si, ti, ki) 1))
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref, "nevertreated")
  expect_identical(cat$targets$ref, -1)
  # Values equal in default display remain distinct dictionary entries.
  eps <- 1 + 2^-52
  fx <- list(calendar = c(1, eps, 2, 3), cohorts = c(2, 3, Inf), targets = 2, anticipation = 0,
             control_group = "notyettreated",
             occupancy = cat_grid_occupancy(c(1, eps, 2, 3), c(2, 3, Inf), function(si, ti, ki) 1))
  cat <- cat_build(fx)
  cat_check_reference(cat, fx, ref)
  expect_identical(format(1), format(eps))
  expect_identical(cat$targets$ref, eps)
  expect_identical(sort(unique(cat$rows$s[cat$rows$s < 2])), c(1, eps))
  links <- cdid:::comparison_catalog_design_links(cat, data.frame(g = 2, t = c(2, 3)))
  expect_identical(cdid:::build_calendar_chained_design(links, fx$calendar, fx$cohorts, 0)$links$ref, rep(eps, nrow(links)))
  # Signed zero is one stored date; duplicates are refused by the accepted helper.
  expect_error(cat_build(list(calendar = c(-0, 0, 1), cohorts = c(1, Inf), targets = 1, anticipation = 0,
    control_group = "nevertreated", occupancy = data.frame(cohort = 1, s = 0, t = 1, count = 1))),
    class = "cdid_design_calendar_error")
  # Overflowing cutoff arithmetic is refused, not certified, including a late date
  # beyond every cohort and a low cohort whose cutoff c - delta overflows.
  expect_error(cat_build(list(calendar = c(-1e308, 0, 1e308), cohorts = c(1.5e308, Inf), targets = 1.5e308,
    anticipation = 1e308, control_group = "notyettreated",
    occupancy = data.frame(cohort = 1.5e308, s = 0, t = 1e308, count = 1))), class = "cdid_design_arithmetic_error")
  no_occ <- data.frame(cohort = numeric(), s = numeric(), t = numeric(), count = numeric())
  expect_error(cdid:::calendar_reference_date(1, c(0, 1e308), 1e308), class = "cdid_design_arithmetic_error")
  expect_error(cat_build(list(calendar = c(0, 1e308), cohorts = 1, targets = 1, anticipation = 1e308,
    control_group = "notyettreated", occupancy = no_occ)), class = "cdid_design_arithmetic_error")
  expect_error(cat_build(list(calendar = c(0, 1), cohorts = c(-1e308, 5), targets = 5, anticipation = 1e308,
    control_group = "notyettreated", occupancy = no_occ)), class = "cdid_design_arithmetic_error")
  # Without overflow an empty occupancy is a valid, honest empty catalog.
  empty <- cat_build(list(calendar = c(0, 1), cohorts = c(1, Inf), targets = 1, anticipation = 0,
    control_group = "notyettreated", occupancy = no_occ), requested_post = data.frame(g = 1, t = 1))
  expect_identical(nrow(empty$rows), 0L)
  expect_identical(empty$targets$ref, 0)
  expect_identical(empty$requested_unsupported$reason, "no_target_pair_occupancy")
  # d=0: 0<1 (1 comparison); d=1: 1<1 fails, advance, 1<Inf (2 comparisons).
  expect_identical(empty$resources$work_performed[["cutoff_comparisons"]], 3)
})

test_that("invalid catalog inputs are refused deterministically", {
  fx <- cat_witness()
  bad <- function(..., class) {
    args <- list(calendar = fx$calendar, cohorts = fx$cohorts, targets = fx$targets, anticipation = 0,
      occupancy = fx$occupancy, control_group = "notyettreated", conventions = "conv", requested_post = NULL,
      bounds = cat_caps)
    over <- list(...)
    for (nm in names(over)) args[nm] <- list(over[[nm]])
    expect_error(do.call(cdid:::build_comparison_catalog, args), class = class)
  }
  occ <- fx$occupancy
  bad(occupancy = rbind(occ, occ[1L, ]), class = "cdid_catalog_duplicate_error")
  bad(occupancy = transform(occ, count = replace(count, 1L, -1)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, count = replace(count, 1L, 0.5)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, count = replace(count, 1L, NA)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, count = replace(count, 1L, 2^53)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, t = replace(t, 1L, 2 + 1e-12)), class = "cdid_catalog_calendar_error")
  bad(occupancy = transform(occ, s = replace(s, 1L, 2)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, cohort = replace(cohort, 1L, 7)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, cohort = replace(cohort, 1L, NaN)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = transform(occ, s = as.character(s)), class = "cdid_catalog_occupancy_error")
  bad(occupancy = cbind(occ, extra = 1), class = "cdid_catalog_occupancy_error")
  bad(occupancy = occ[, 1:3], class = "cdid_catalog_occupancy_error")
  bad(targets = c(3, 3), class = "cdid_catalog_target_error")
  bad(targets = 7, class = "cdid_catalog_target_error")
  bad(targets = Inf, class = "cdid_catalog_target_error")
  bad(targets = numeric(), class = "cdid_catalog_target_error")
  bad(control_group = "all", class = "cdid_catalog_input_error")
  bad(control_group = NA_character_, class = "cdid_catalog_input_error")
  bad(conventions = "", class = "cdid_catalog_input_error")
  bad(conventions = c("a", "b"), class = "cdid_catalog_input_error")
  bad(anticipation = 0.5, class = "cdid_design_anticipation_error")
  bad(cohorts = c(3, 4, 4), class = "cdid_design_cohort_error")
  bad(requested_post = data.frame(g = 3, t = 2), class = "cdid_catalog_request_error")
  bad(requested_post = data.frame(g = 3, t = 4.5), class = "cdid_catalog_request_error")
  bad(requested_post = data.frame(g = 4, t = 4), class = "cdid_catalog_request_error")
  bad(requested_post = data.frame(g = c(3, 3), t = c(4, 4)), class = "cdid_catalog_duplicate_error")
  bad(requested_post = data.frame(g = 3, time = 4), class = "cdid_catalog_request_error")
  bad(bounds = cat_caps[-1L], class = "cdid_catalog_bound_error")
  bad(bounds = c(cat_caps, list(max_other = 1)), class = "cdid_catalog_bound_error")
  # Invalid caps are refused by validation itself, not by a later cap breach.
  for (v in list(-1, 1e6 + 0.5, Inf, 2^53, "1", c(1, 2), NA_real_)) {
    expect_error(cat_build(fx, bounds = modifyList(cat_caps, list(max_rows = v))),
                 regexp = "max_rows` must be an explicit finite whole number", class = "cdid_catalog_bound_error")
  }
  expect_error(cdid:::build_comparison_catalog(fx$calendar, fx$cohorts, 3, 0, fx$occupancy, "notyettreated", "conv", NULL),
               "bounds")
  # Declared qualifications: only referenced targets' causal stored dates.
  cat <- cat_build(fx)
  expect_error(cdid:::comparison_catalog_design_links(cat, data.frame(g = 3, t = 1)), class = "cdid_catalog_qualification_error")
  expect_error(cdid:::comparison_catalog_design_links(cat, data.frame(g = 3, t = 3.5)), class = "cdid_catalog_qualification_error")
  expect_error(cdid:::comparison_catalog_design_links(cat, data.frame(g = 4, t = 4)), class = "cdid_catalog_qualification_error")
  expect_error(cdid:::comparison_catalog_design_links(cat, data.frame(g = c(3, 3), t = c(4, 4))), class = "cdid_catalog_duplicate_error")
  expect_error(cdid:::comparison_catalog_design_links(list(rows = cat$rows), NULL), class = "cdid_catalog_input_error")
  # Undeclared conditional endpoints stay private: no post level is estimable.
  des <- cdid:::build_comparison_catalog_design(cat, NULL, target_post_times = 4)
  expect_identical(des$post_targets$status, "isolated_unsupported")
  expect_identical(sum(des$parameters$role == "post_nuisance"), 1L)
})

test_that("resource caps are exact, preflighted and refuse before the corresponding allocation", {
  ref <- load_catalog_reference()
  fx <- cat_fixtures()$regular_nt_inf
  req <- data.frame(g = c(2, 3), t = c(5, 6))
  cat <- cat_build(fx, requested_post = req)
  want <- cat_check_reference(cat, fx, ref)
  res <- cat$resources
  # Independent scratch formula (entries charged before each transient allocation).
  occ <- fx$occupancy
  n_pos <- sum(occ$count > 0)
  n_tp <- sum(occ$count > 0 & occ$cohort %in% fx$targets)
  n_ref <- sum(!is.na(cat$targets$ref))
  R <- res$rows; P <- res$pool_entries_preintern; Cn <- res$candidates
  expect_identical(res$scratch_charged,
    64 * (length(fx$calendar) + length(fx$cohorts) + length(fx$targets) + nrow(occ) + nrow(req) + 1) +
      16 * nchar(cat$conventions, type = "bytes") +
      256 * (Cn + n_tp + 1) + 16 * nchar(cat$conventions, type = "bytes") * Cn +
      5 * nrow(occ) + length(fx$calendar) + length(fx$targets) + 4 * n_pos + 8 * n_tp + 5 * Cn + 64 * P + R + 4 * R +
      4 * (n_ref + 2 * R) + 4 * nrow(req))
  need <- c(max_candidates = Cn, max_rows = R, max_pool_entries = P, max_scratch = res$scratch_charged,
            max_work = res$work_bound)
  for (nm in names(need)) {
    ok <- cat_build(fx, requested_post = req, bounds = modifyList(cat_caps, as.list(need[nm])))
    expect_identical(ok$rows, cat$rows)
    lo <- modifyList(cat_caps, setNames(list(need[[nm]] - 1), nm))
    expect_error(cat_build(fx, requested_post = req, bounds = lo), class = "cdid_catalog_bound_error")
  }
  exact <- cat_build(fx, requested_post = req, bounds = as.list(need))
  expect_identical(exact$rows, cat$rows)
  expect_identical(exact$resources$bounds, as.list(need))

  # A huge O(J T^3) candidate family is refused from its exact preflight count.
  nT <- 3000
  big <- list(calendar = as.double(seq_len(nT)), cohorts = c(2, Inf), targets = 2, anticipation = 0,
              control_group = "notyettreated",
              occupancy = data.frame(cohort = rep(c(2, Inf), each = nT - 1), s = 1, t = rep(2:nT, 2L), count = 1))
  n_cand <- sum(1 + (nT - pmax(2, 2:nT) + 1))
  msg <- tryCatch(cat_build(big, bounds = modifyList(cat_caps, list(max_candidates = 1e5))),
                  cdid_catalog_bound_error = function(e) conditionMessage(e))
  expect_identical(msg, sprintf("Exact candidate count %.0f exceeds explicit `max_candidates`.", n_cand))
  # Each refusal happens at its own preflight, before the work/allocation it guards.
  expect_error(cat_build(big, bounds = modifyList(cat_caps, list(max_work = 3000))),
               regexp = "Header validation/sort work", class = "cdid_catalog_bound_error")
  expect_error(cat_build(big, bounds = modifyList(cat_caps, list(max_candidates = 1e7, max_work = 1e7))),
               regexp = "Work bound .* exceeds", class = "cdid_catalog_bound_error")
  expect_error(cat_build(big, bounds = modifyList(cat_caps, list(max_scratch = 1000))),
               regexp = "Scratch for header validation/sort buffers", class = "cdid_catalog_bound_error")
  expect_error(cat_build(big, bounds = modifyList(cat_caps, list(max_candidates = 1e7, max_scratch = 5e6))),
               regexp = "Work bound .* exceeds|Scratch for candidate/dedup/output buffers", class = "cdid_catalog_bound_error")
})


test_that("header and byte budgets precede validation and adapter allocations", {
  fx <- cat_witness()
  bad <- fx; bad$calendar <- c(NaN,2:5)
  expect_error(cat_build(bad,bounds=modifyList(cat_caps,list(max_work=0))),
               regexp="before scans",class="cdid_catalog_bound_error")
  expect_error(cat_build(bad,bounds=modifyList(cat_caps,list(max_scratch=0))),
               regexp="header validation/sort",class="cdid_catalog_bound_error")
  expect_error(cat_build(fx,conventions=strrep("x",10000),bounds=modifyList(cat_caps,list(max_scratch=1000))),
               class="cdid_catalog_bound_error")
  cat <- cat_build(fx)
  expect_gt(cat$resources$work_bound,sum(cat$resources$work_performed))
  expect_match(cat$resources$work_definition,"repeated pair-pool scans")
  small <- cat;small$resources$bounds$max_scratch<-0
  expect_error(cdid:::comparison_catalog_design_links(small,NULL),class="cdid_catalog_bound_error")
  expect_error(cdid:::build_comparison_catalog_design(small,NULL,4),class="cdid_catalog_bound_error")
})


test_that("pool offsets and composite identities do not truncate exact long counts", {
  expect_identical(cdid:::catalog_pool_offsets(c(2,3)),c(1L,3L))
  expect_identical(cdid:::catalog_pool_offsets(c(2^31-1,2)),c(1,2^31))
  expect_error(cdid:::catalog_pool_offsets(c(2^52,1)),class="cdid_catalog_bound_error")
  expect_error(cdid:::catalog_product(2^27,2^27,"operator key"),class="cdid_catalog_bound_error")
})
