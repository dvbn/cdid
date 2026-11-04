# Private counterfactual outcome fits over caller-owned unique comparison blocks.
# The source owns support summaries/scales, calendar L, fit identities and row
# replay. No catalog, upstream forcing, link IF or backend budget lives here.

# Metadata-only preflight, also available before an owner's budget reservation.
# Dimension bounds describe retained operators; they are NOT budget charges.
or_outcome_group_layout <- function(recipe, supports, blocks, links, mode,
                                    route = "pooled") {
  or_pooled_check_recipe(recipe, recipe$columns, route)
  if (missing(mode) || !is.character(mode) || is.object(mode) ||
      !is.null(dim(mode)) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("pooled", "separate"))
    or_pooled_abort("Supply explicit private outcome mode: pooled or separate.", "input", route)
  if (!is.list(supports) || !length(supports) || !is.list(blocks) ||
      !length(blocks) || !is.list(links))
    or_pooled_abort("Supply compact supports, nonempty unique blocks and consuming links.",
      "input", route)
  ids_ok <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n)
  for (b in blocks) {
    if (!is.list(b) || !is.character(b$name) || length(b$name) != 1L ||
        is.na(b$name) || !nzchar(b$name) ||
        !ids_ok(b$support_ids, length(supports)) || !length(b$support_ids) ||
        anyDuplicated(b$support_ids) || !is.numeric(b$L) || is.object(b$L) ||
        !is.null(dim(b$L)) || length(b$L) != 1L || !is.finite(b$L) ||
        b$L <= 0 || b$L != floor(b$L))
      or_pooled_abort("Each unique block needs a name, distinct support IDs and positive calendar L.",
        "input", route)
  }
  block_names <- vapply(blocks, `[[`, character(1), "name")
  if (anyDuplicated(block_names))
    or_pooled_abort("Comparison block names must be unique; reuse caller block IDs for identical fits.",
      "input", route)
  for (l in links) {
    if (!is.list(l) || !ids_ok(l$regression, length(blocks)) ||
        length(l$regression) != 1L || !ids_ok(l$support_id, length(supports)) ||
        length(l$support_id) != 1L)
      or_pooled_abort("Each consuming link needs an original block ID and support ID.", "input", route)
  }
  B <- length(blocks)
  group_id <- if (mode == "pooled") rep.int(1L, B) else seq_len(B)
  block_map <- data.frame(block_id = seq_len(B), name = block_names,
    group_id = group_id, local_block_id = if (mode == "pooled") seq_len(B) else rep.int(1L, B),
    stringsAsFactors = FALSE)
  group_map <- lapply(seq_len(max(group_id)), function(g) {
    bids <- which(group_id == g)
    lids <- which(vapply(links, function(l) l$regression %in% bids, logical(1)))
    # In separate mode include precisely this block's controls and ALL of its
    # consuming prediction supports. Unrelated treated domains must stay out.
    sids <- if (mode == "pooled") seq_along(supports) else sort(unique(c(
      unlist(lapply(blocks[bids], `[[`, "support_ids"), use.names = FALSE),
      vapply(links[lids], function(l) as.integer(l$support_id), integer(1)))))
    list(block_ids = bids, support_ids = as.integer(sids), link_ids = lids,
      slope_dimension_bound = recipe$p,
      coefficient_dimension_bound = recipe$p + recipe$intercept * length(bids),
      operator_dimension_bound = c(recipe$p, recipe$p),
      summary_dimension_bound = c(recipe$p + 1L, recipe$p + 1L))
  })
  list(block_map = block_map, group_map = group_map)
}

# Same compact contract as or_pooled_fit_blocks, with an explicit private mode.
# The pooled call and its result are unchanged. Separate calls fit one unique
# caller block, never one link; no numerical deduplication or cross-response QR
# cache is inferred here. Owners must reserve scans/work/retention before use.
or_outcome_fit_blocks <- function(recipe, supports, blocks, links, rank_tol, N,
                                  visit_rows, mode, route = "pooled") {
  layout <- or_outcome_group_layout(recipe, supports, blocks, links, mode, route)
  if (!is.numeric(N) || is.object(N) || !is.null(dim(N)) || length(N) != 1L ||
      !is.finite(N) || N <= 0 || N != floor(N) || !is.function(visit_rows))
    or_pooled_abort("Supply the original positive integer roster N and a visit_rows callback.",
      "input", route)
  if (mode == "pooled") {
    groups <- list(or_pooled_fit_blocks(recipe, supports, blocks, links,
      rank_tol, N, visit_rows, route))
  } else {
    groups <- lapply(layout$group_map, function(m) {
      local_blocks <- lapply(blocks[m$block_ids], function(b) {
        b$support_ids <- match(b$support_ids, m$support_ids); b
      })
      local_links <- lapply(links[m$link_ids], function(l) {
        l$regression <- match(l$regression, m$block_ids)
        l$support_id <- match(l$support_id, m$support_ids); l
      })
      visit_local <- function(sid, visit) visit_rows(m$support_ids[sid], visit)
      or_pooled_fit_blocks(recipe, supports[m$support_ids], local_blocks, local_links,
        rank_tol, N, visit_local, route)
    })
  }
  rss <- sum(vapply(groups, `[[`, numeric(1), "rss"))
  retained <- lapply(groups, function(g) list(slope_dimension = length(g$keep),
    operator_dimension = dim(g$fit$R),
    summary_dimension = dim(if (recipe$intercept) g$summary$residual$factor else g$summary$factor),
    alias_map_dimension = dim(g$basis$alias_map),
    block_constants_dimension = dim(g$basis$block_constants),
    verification_scan = g$verification_scan, alias_rows_checked = g$basis$rows_checked))
  structure(c(list(mode = mode, recipe = recipe, N = N, route = route), layout,
    list(groups = groups, rss = rss, objective = rss / N, retained = retained)),
    class = "cdid_or_outcome_groups")
}

or_outcome_block_index <- function(object, block) {
  if (!inherits(object, "cdid_or_outcome_groups"))
    or_pooled_abort("Expected private counterfactual outcome groups.")
  if (is.character(block) && !is.object(block) && is.null(dim(block)) &&
      length(block) == 1L && !is.na(block)) {
    block <- match(block, object$block_map$name)
  }
  if (!is.numeric(block) || is.object(block) || !is.null(dim(block)) ||
      length(block) != 1L || !is.finite(block) || block != floor(block) ||
      block < 1L || block > nrow(object$block_map))
    or_pooled_abort("Unknown original outcome block ID or name.", "input", object$route)
  as.integer(block)
}

# One block's intercept/means and its group's retained slope basis/operator.
# In pooled mode fit/basis describe the joint group, local_block_id selects the
# intercept record. There is no standalone block slope solve in that mode.
or_outcome_block_fit <- function(object, block) {
  b <- or_outcome_block_index(object, block); m <- object$block_map[b, ]
  g <- object$groups[[m$group_id]]
  c(list(block_id = b, group_id = m$group_id, local_block_id = m$local_block_id),
    g$regressions[[m$local_block_id]], list(p = length(g$keep),
      coefficient_dimension = length(g$keep) + object$recipe$intercept,
      keep = g$keep, beta = g$beta, basis = g$basis, fit = g$fit))
}

or_outcome_block_p <- function(object, block) {
  b <- or_outcome_block_index(object, block)
  length(object$groups[[object$block_map$group_id[b]]]$keep)
}

# Z has the source recipe's columns, in order, on the declared fit/prediction
# domain already checked by visit_rows. Retained d is always L * z, even when
# separate slopes could equivalently be expressed on unscaled z.
or_outcome_block_predict <- function(object, block, Z) {
  f <- or_outcome_block_fit(object, block); recipe <- object$recipe
  if (!is.matrix(Z) || !is.numeric(Z) || is.object(Z) ||
      ncol(Z) != length(recipe$columns) || any(!is.finite(Z)) ||
      (!is.null(colnames(Z)) && !identical(colnames(Z), recipe$columns)) ||
      (recipe$intercept && any(Z[, recipe$intercept_index] != 1)))
    or_pooled_abort("Prediction rows must have the finite source outcome recipe columns.",
      "input", object$route)
  d <- f$L * Z[, recipe$slope_index[f$keep], drop = FALSE]
  as.numeric(f$alpha + if (f$p) as.vector(d %*% f$beta) else numeric(nrow(Z)))
}

# Normalized retained slope Gram solve, partialled only with an intercept.
# This is an operator accessor, not a coefficient or downstream link IF.
or_outcome_block_gram_solve <- function(object, block, b) {
  f <- or_outcome_block_fit(object, block)
  if (object$recipe$intercept) or_pooled_wls_gram_solve(f$fit, b, normalizer = object$N) else
    or_wls_gram_solve(f$fit, b, normalizer = object$N)
}
