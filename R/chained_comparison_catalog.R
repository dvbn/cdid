# Private CP02 approved comparison catalog: D10 retention of every valid
# original comparison plus D13 fixed-horizon not-yet-treated reevaluations.
# Metadata only: no outcomes, fits, weights, estimator, IF or inference. A row
# declares one population comparison (target g, pair e, actual pool K) under a
# single caller-declared opaque convention. Equal K is the same population
# operator under those conventions, never a certificate that fitted nuisances
# or samples coincide. Endpoint roles are conditional eligibilities; counts and
# control cleanliness never verify PT/MTAR/overlap/model/integrability premises.
# Proofs: docs/theory/cp02_approved_comparison_catalog.tex.

catalog_abort <- function(message, class) {
  cdid_abort(message, c(class, "cdid_catalog_error"))
}

# The same condition as `catalog_abort()`, with typed diagnostics (when they
# could be built) for the restatement of resource refusals.
catalog_abort_recorded <- function(message, class, diagnostics = NULL) {
  if (!is.list(diagnostics)) catalog_abort(message, class)
  stop(structure(class = c(class, "cdid_catalog_error", "cdid_error", "error", "condition"),
    list(message = message, call = NULL, diagnostics = diagnostics)))
}

# Typed counts of the adapter envelope check: the work (`max_work`, from the
# work limit) when it fails, otherwise the scratch (`max_scratch`, from the
# workspace limit).
catalog_bound_diagnostics <- function(env, b) {
  if (env$work > b$max_work) list(label = "catalog adapter output, string and graph work", requested = env$work,
    bound = b$max_work, limit_kind = "work") else list(label = "catalog adapter output, string and graph scratch",
    requested = env$scratch, bound = b$max_scratch, limit_kind = "workspace")
}

catalog_bound_names <- c("max_candidates", "max_rows", "max_pool_entries", "max_scratch", "max_work")

# Caller caps are explicit private controls; there are no defaults.
validate_catalog_bounds <- function(bounds) {
  if (!is.list(bounds) || is.object(bounds) || is.null(names(bounds)) ||
      length(bounds) != length(catalog_bound_names) || anyDuplicated(names(bounds)) ||
      !setequal(names(bounds), catalog_bound_names)) {
    catalog_abort("Supply explicit `bounds` naming exactly max_candidates, max_rows, max_pool_entries, max_scratch and max_work.",
      "cdid_catalog_bound_error")
  }
  out <- lapply(catalog_bound_names, function(nm) {
    x <- bounds[[nm]]
    if (!design_vector(x, 1L) || x < 0 || x != floor(x) || x > 2^53 - 1) {
      catalog_abort(paste0("`bounds$", nm, "` must be an explicit finite whole number in [0, 2^53 - 1]."),
        "cdid_catalog_bound_error")
    }
    as.double(x)
  })
  names(out) <- catalog_bound_names
  out
}

# Sparse caller-supplied cohort/pair occupancy: exactly columns cohort, s, t,
# count. Dates/cohorts are matched exactly to dictionary indices later.
validate_catalog_occupancy <- function(occupancy) {
  cols <- c("cohort", "s", "t", "count")
  if (!is.data.frame(occupancy) || ncol(occupancy) != 4L || anyDuplicated(names(occupancy)) ||
      !setequal(names(occupancy), cols)) {
    catalog_abort("`occupancy` must be a data.frame with exactly the columns cohort, s, t and count.",
      "cdid_catalog_occupancy_error")
  }
  n <- nrow(occupancy)
  k <- occupancy$cohort
  if (!design_vector(k, n, finite = FALSE) || any(is.nan(k) | k == -Inf)) {
    catalog_abort("Occupancy cohorts must be numeric; only +Inf (never treated) is non-finite.",
      "cdid_catalog_occupancy_error")
  }
  if (!design_vector(occupancy$s, n) || !design_vector(occupancy$t, n)) {
    catalog_abort("Occupancy dates must be finite numeric stored dates.", "cdid_catalog_occupancy_error")
  }
  cnt <- occupancy$count
  if (!design_vector(cnt, n) || any(cnt < 0 | cnt != floor(cnt) | cnt > 2^53 - 1)) {
    catalog_abort("Occupancy counts must be finite whole numbers in [0, 2^53 - 1].", "cdid_catalog_occupancy_error")
  }
  list(cohort = as.double(k), s = as.double(occupancy$s), t = as.double(occupancy$t), count = as.double(cnt))
}

validate_catalog_requests <- function(requested_post, cal, coh, tg) {
  if (is.null(requested_post)) return(list(gi = integer(), ti = integer()))
  if (!is.data.frame(requested_post) || ncol(requested_post) != 2L ||
      anyDuplicated(names(requested_post)) || !setequal(names(requested_post), c("g", "t"))) {
    catalog_abort("`requested_post` must be NULL or a data.frame with exactly columns g and t.", "cdid_catalog_request_error")
  }
  n <- nrow(requested_post)
  g <- requested_post$g; t <- requested_post$t
  if (!design_vector(g, n) || !design_vector(t, n)) {
    catalog_abort("Requested post targets require finite numeric g and t.", "cdid_catalog_request_error")
  }
  gi <- match(as.double(g), coh); ti <- match(as.double(t), cal)
  if (anyNA(gi) || any(!gi %in% tg)) {
    catalog_abort("A requested post target cohort must be one of `targets`.", "cdid_catalog_request_error")
  }
  if (anyNA(ti)) {
    catalog_abort("A requested post date must be a stored date of the full calendar.", "cdid_catalog_request_error")
  }
  if (any(as.double(t) < as.double(g))) {
    catalog_abort("Requested post targets require t >= g.", "cdid_catalog_request_error")
  }
  if (anyDuplicated(data.frame(gi, ti))) {
    catalog_abort("Duplicate requested post target.", c("cdid_catalog_duplicate_error", "cdid_catalog_request_error"))
  }
  list(gi = gi, ti = ti)
}

# Conservative explicit payload/work envelopes; not CPU instruction or RSS caps.
# Count products/sums refuse before inexact aggregate arithmetic or allocation.
catalog_size <- function(x, label) {
  if (is.numeric(x) && length(x) > 0L && all(is.finite(x)) && all(x >= 0) && all(x == floor(x))) {
    s <- sum(as.double(x))
    if (is.finite(s) && s <= (2^53 - 1)) return(s)
  }
  total <- 0
  for (v in x) {
    if (!is.finite(v) || v < 0 || v != floor(v) || v > (2^53 - 1) - total)
      catalog_abort(paste0("Inexact or unsupported size for ", label, "."), "cdid_catalog_bound_error")
    total <- total + v
  }
  total
}
catalog_product <- function(a, b, label) {
  if (length(a) == 1L && length(b) == 1L && is.finite(a) && is.finite(b) &&
      a >= 0 && b >= 0 && a == floor(a) && b == floor(b)) {
    p <- as.double(a) * as.double(b)
    if (is.finite(p) && p <= (2^53 - 1)) return(p)
  }
  if (a < 0 || b < 0 || !is.finite(a) || !is.finite(b) ||
      (a != 0 && b > (2^53 - 1) / a))
    catalog_abort(paste0("Inexact or unsupported size for ", label, "."), "cdid_catalog_bound_error")
  catalog_size(a * b, label)
}

# Offsets can exceed INT_MAX even though each cohort dictionary index cannot.
# Keep the historical integer representation for small dictionaries; use exact
# double offsets for long entry vectors. Check the actual 64-bit R vector limit.
catalog_pool_offsets <- function(pool_size) {
  total <- catalog_size(pool_size, "interned pool entries")
  if (total > 2^52 || .Machine$sizeof.pointer != 8L)
    catalog_abort("Interned pool entry vector exceeds the supported 64-bit R vector length.", "cdid_catalog_bound_error")
  starts <- cumsum(c(1, as.double(pool_size)))[seq_along(pool_size)]
  if (total < .Machine$integer.max) as.integer(starts) else starts
}

build_comparison_catalog <- function(calendar, cohorts, targets, anticipation, occupancy,
                                     control_group, conventions, requested_post, bounds,
                                     eligibility_occupancy = NULL) {
  bnd <- validate_catalog_bounds(bounds)
  # Header-only dimensions precede validation/coercion/sort/match. The accepted
  # helper domain uses integer dictionary and data-frame indices. C+1 is a real
  # sentinel, so C must be strictly below INT_MAX; counts remain separate.
  nT0 <- length(calendar); nC0 <- length(cohorts); nJ0 <- length(targets)
  nO0 <- if (is.data.frame(occupancy)) nrow(occupancy) else 0
  nQ0 <- if (is.data.frame(requested_post)) nrow(requested_post) else 0
  if (nT0 >= .Machine$integer.max || nC0 >= .Machine$integer.max ||
      nJ0 > .Machine$integer.max || nO0 >= .Machine$integer.max || nQ0 > .Machine$integer.max)
    catalog_abort("Catalog dictionary/data-frame index header exceeds supported integer range.", "cdid_catalog_bound_error")
  H <- catalog_size(c(nT0, nC0, nJ0, nO0, nQ0, 1), "header dimensions")
  byte_conventions <- if (is.character(conventions) && !is.object(conventions) && length(conventions) == 1L && !anyNA(conventions))
    as.double(nchar(conventions, type = "bytes")) else 0
  scratch <- 0
  charge <- function(n, what) {
    n <- catalog_size(n, what)
    if (n > bnd$max_scratch - scratch)
      catalog_abort(paste0("Scratch for ", what, " exceeds explicit `max_scratch` before allocation."), "cdid_catalog_bound_error")
    scratch <<- scratch + n
  }
  header_scratch <- catalog_size(c(catalog_product(64, H, "header payload"),
    catalog_product(16, byte_conventions, "convention bytes")), "header scratch")
  header_work <- catalog_size(c(catalog_product(64 * (1 + ceiling(log2(H + 1))), H, "header validation/sort work"),
    catalog_product(16, byte_conventions, "convention byte work")), "header work")
  if (header_work > bnd$max_work)
    catalog_abort("Header validation/sort work exceeds explicit `max_work` before scans.", "cdid_catalog_bound_error")
  charge(header_scratch, "header validation/sort buffers")
  spec <- validate_calendar_spec(calendar, cohorts, anticipation)
  cal <- spec$calendar; coh <- spec$cohorts; delta <- spec$anticipation
  nT <- length(cal); nC <- length(coh)
  if (as.double(nT) * nC > 2^53 - 1) {
    catalog_abort("Calendar/cohort dictionary product exceeds exact index arithmetic.", "cdid_catalog_bound_error")
  }
  inf_idx <- if (spec$has_never_treated) nC else 0L # +Inf sorts last.
  if (!design_vector(control_group, 1L, "character") ||
      !control_group %in% c("nevertreated", "notyettreated")) {
    catalog_abort("`control_group` must be explicitly \"nevertreated\" or \"notyettreated\".", "cdid_catalog_input_error")
  }
  nyt <- identical(control_group, "notyettreated")
  if (!design_vector(conventions, 1L, "character") || !nzchar(conventions)) {
    catalog_abort("Declare the common population universe/mixture/covariate/weight conventions as one nonempty opaque string.",
      "cdid_catalog_input_error")
  }
  if (!design_vector(targets, length(targets)) || length(targets) < 1L || anyDuplicated(as.double(targets))) {
    catalog_abort("`targets` must be a nonempty vector of unique finite cohorts.", "cdid_catalog_target_error")
  }
  tg <- match(as.double(targets), coh)
  if (anyNA(tg)) catalog_abort("Every target must belong to the supplied cohort universe.", "cdid_catalog_target_error")
  tg <- sort(tg)
  nJ <- length(tg)
  occ <- validate_catalog_occupancy(occupancy)
  n_occ <- length(occ$count)
  req <- validate_catalog_requests(requested_post, cal, coh, tg)
  n_req <- length(req$gi)

  base_work <- header_work

  # Exact dictionary indices; refuse off-grid dates, unknown cohorts, duplicates.
  charge(5 * n_occ, "occupancy indices")
  ki <- match(occ$cohort, coh); si <- match(occ$s, cal); ti <- match(occ$t, cal)
  if (anyNA(ki)) catalog_abort("Occupancy cohort outside the supplied cohort universe.", "cdid_catalog_occupancy_error")
  if (anyNA(si) || anyNA(ti)) {
    catalog_abort("Occupancy dates must be exactly stored dates of the full calendar.",
      c("cdid_catalog_occupancy_error", "cdid_catalog_calendar_error"))
  }
  if (any(si >= ti)) catalog_abort("Occupancy pairs require s < t.", "cdid_catalog_occupancy_error")
  ord <- order(si, ti, ki, method = "radix")
  if (n_occ > 1L) {
    a <- ord[-n_occ]; b <- ord[-1L]
    if (any(si[a] == si[b] & ti[a] == ti[b] & ki[a] == ki[b])) {
      catalog_abort("Duplicate occupancy (cohort, s, t) index; supply one explicit count per cell.",
        c("cdid_catalog_duplicate_error", "cdid_catalog_occupancy_error"))
    }
  }

  # cstar[d] = first cohort index k with d + delta < k, exactly (nC + 1 if none).
  # Monotone in both d and k, so a two-pointer scan needs at most nT + nC
  # comparisons. Every stored date and the smallest cohort are compared at least
  # once, so overflow refuses exactly when exact_cmp_sum refuses on some pair.
  charge(nT + nJ, "cutoff table and references")
  cstar <- integer(nT); ci <- 1L; n_cmp <- 0
  for (d in seq_len(nT)) {
    repeat {
      n_cmp <- n_cmp + 1
      below <- exact_cmp_sum(cal[d], delta, coh[min(ci, nC)]) < 0L
      if (ci > nC || below) break
      ci <- ci + 1L
    }
    cstar[d] <- ci
  }
  # r_g = max{d : d + delta < g} and d + delta < g iff cstar[d] <= index(g).
  ref_i <- findInterval(tg, cstar)
  ref_i[ref_i == 0L] <- NA_integer_
  hfirst <- findInterval(coh[tg], cal, left.open = TRUE) + 1L # first stored date >= g

  # Issue #16: restricted eligibility contexts. Each restricted context c is a
  # copy of the owner occupancy minus the pair counts of its excluded units;
  # a "pair" is then (context, s, t), and a target reads only its own
  # context's cells (unrestricted targets keep context 0, the owner cells).
  cell_cnt <- occ$count; cell_ctx <- NULL; tg_ctx <- NULL
  if (!is.null(eligibility_occupancy)) {
    ex <- catalog_eligibility_cells(eligibility_occupancy, occ, ki, si, ti, coh, cal, targets, tg)
    charge(6 * length(ex$count), "eligibility occupancy cells")
    base_work <- catalog_size(c(base_work, catalog_product(64 * (1 + ceiling(log2(length(ex$count) + 2))),
      length(ex$count) + 1, "eligibility occupancy cell work")), "header plus eligibility occupancy work")
    ki <- c(ki, ex$ki); si <- c(si, ex$si); ti <- c(ti, ex$ti)
    cell_cnt <- c(cell_cnt, ex$count); cell_ctx <- c(integer(n_occ), ex$context); tg_ctx <- ex$target_context
    ord <- order(cell_ctx, si, ti, ki, method = "radix")
  }

  # Positive cells grouped by pair (CSR); cohorts ascending within each pair.
  pos <- ord[cell_cnt[ord] > 0]
  n_pos <- length(pos)
  charge(4 * n_pos, "positive support index")
  cell_k <- ki[pos]; cs <- si[pos]; ct <- ti[pos]
  new_pair <- if (!n_pos) logical() else if (is.null(cell_ctx)) c(TRUE, cs[-1L] != cs[-n_pos] | ct[-1L] != ct[-n_pos]) else {
    cx <- cell_ctx[pos]
    c(TRUE, cx[-1L] != cx[-n_pos] | cs[-1L] != cs[-n_pos] | ct[-1L] != ct[-n_pos])
  }
  pair_of <- cumsum(new_pair)
  pair_start <- which(new_pair)
  pair_len <- diff(c(pair_start, n_pos + 1L))
  pair_s <- cs[pair_start]; pair_t <- ct[pair_start]
  pair_inf <- logical(length(pair_start))
  if (inf_idx > 0L) pair_inf[pair_of[cell_k == inf_idx]] <- TRUE

  # Supported target pairs: positive target occupancy at that same pair.
  tj <- match(cell_k, tg)
  tc <- if (is.null(cell_ctx)) which(!is.na(tj)) else which(!is.na(tj) & cell_ctx[pos] == tg_ctx[pmax(1L, tj)])
  n_tp <- length(tc)
  charge(8 * n_tp, "target-pair index")
  tc <- tc[order(tj[tc], pair_of[tc], method = "radix")]
  tp_j <- tj[tc]; tp_pair <- pair_of[tc]; tp_count <- cell_cnt[pos[tc]]
  tp_ref <- ref_i[tp_j]
  tp_h0 <- pmax(hfirst[tp_j], pair_t[tp_pair])
  tp_nh <- if (nyt) pmax(0L, nT - tp_h0 + 1L) else integer(n_tp)
  tp_nc <- (1 + as.double(tp_nh)) * !is.na(tp_ref)
  tp_inf <- pair_inf[tp_pair]

  # Exact candidate preflight before any candidate allocation.
  n_cand <- sum(tp_nc)
  if (n_cand > bnd$max_candidates) {
    catalog_abort(sprintf("Exact candidate count %.0f exceeds explicit `max_candidates`.", n_cand), "cdid_catalog_bound_error")
  }
  if (n_cand > .Machine$integer.max) catalog_abort("Candidate count exceeds supported index range.", "cdid_catalog_bound_error")
  # Every supported target pair and retained row may traverse its entire
  # actual cohort pool. Include sort/match/format/split visits and newly created
  # convention/key bytes, not merely candidate-loop iterations.
  stage_size <- catalog_size(c(n_cand, n_tp, 1), "candidate stage")
  stage_work <- catalog_size(c(catalog_product(128 * (nC + 1 + ceiling(log2(stage_size + 1))),
    stage_size, "candidate/pool/sort work"), catalog_product(catalog_product(16, byte_conventions, "key bytes"), n_cand, "operator key work")), "candidate work")
  work_bound <- catalog_size(c(base_work, stage_work), "whole catalog work")
  if (work_bound > bnd$max_work) {
    catalog_abort(sprintf("Work bound %.0f exceeds explicit `max_work`.", work_bound), "cdid_catalog_bound_error")
  }
  charge(catalog_size(c(catalog_product(256, stage_size, "candidate/output buffers"),
    catalog_product(catalog_product(16, byte_conventions, "operator key bytes"), n_cand, "operator key payload")), "candidate scratch"), "candidate/dedup/output buffers")
  charge(3 * n_cand, "candidate arrays")
  cand_tp <- integer(n_cand); cand_h <- integer(n_cand); cand_m <- integer(n_cand)
  later_finite <- function(q) {
    ce <- cell_k[seq.int(pair_start[tp_pair[q]], length.out = pair_len[tp_pair[q]])]
    ce[ce > tg[tp_j[q]] & ce != inf_idx]
  }
  at <- 0L
  for (q in which(tp_nc > 0)) {
    hs <- c(0L, if (tp_nh[q] > 0L) seq.int(tp_h0[q], nT))
    m <- 0L
    if (nyt) {
      # K = (C_e and Inf) plus {k in C_e : k > g, cutoff + delta < k}: an upper
      # suffix of the sorted finite later cohorts, determined by its size m.
      fin <- later_finite(q)
      lo <- cstar[c(pair_t[tp_pair[q]], hs[-1L])]
      m <- length(fin) - findInterval(lo - 1L, fin)
    }
    k <- at + seq_along(hs)
    cand_tp[k] <- q; cand_h[k] <- hs; cand_m[k] <- m
    at <- at + length(hs)
  }

  # Deduplicate only identical (target, pair, actual K) under the one declared
  # convention; provenance is kept as the original flag and a horizon interval.
  charge(2 * n_cand, "candidate deduplication")
  nonempty <- cand_m > 0L | tp_inf[cand_tp]
  o <- order(cand_tp, -cand_m, cand_h, method = "radix")
  o <- o[nonempty[o]]
  n_ne <- length(o)
  new_row <- if (n_ne) c(TRUE, cand_tp[o[-1L]] != cand_tp[o[-n_ne]] | cand_m[o[-1L]] != cand_m[o[-n_ne]]) else logical()
  n_rows <- sum(new_row)
  if (n_rows > bnd$max_rows) {
    catalog_abort(sprintf("Exact retained row count %d exceeds explicit `max_rows`.", n_rows), "cdid_catalog_bound_error")
  }
  grp <- cumsum(new_row); hh <- cand_h[o]
  r_tp <- cand_tp[o[new_row]]; r_m <- cand_m[o[new_row]]; r_inf <- tp_inf[r_tp]
  r_orig <- logical(n_rows); r_orig[grp[hh == 0L]] <- TRUE
  gH <- grp[hh > 0L]; vH <- hh[hh > 0L]
  r_nh <- tabulate(gH, nbins = n_rows)
  r_h1 <- r_h2 <- rep(NA_integer_, n_rows)
  r_h1[gH[!duplicated(gH)]] <- vH[!duplicated(gH)]
  r_h2[gH[!duplicated(gH, fromLast = TRUE)]] <- vH[!duplicated(gH, fromLast = TRUE)]

  # Intern pools as integer cohort-index sets; rows store only pool ids.
  n_entries <- sum(as.double(r_m)) + sum(r_inf)
  if (n_entries > bnd$max_pool_entries) {
    catalog_abort(sprintf("Pool entries %.0f exceed explicit `max_pool_entries` before interning.", n_entries),
      "cdid_catalog_bound_error")
  }
  # Cohort dictionary indices need at most ten digits plus delimiters. Reserve
  # raw character bytes, split strings, integer entries and temporary vectors.
  charge(catalog_size(c(catalog_product(64, n_entries, "pool/key payload"), n_rows), "pool interning"), "pool interning and character buffers")
  if (!nyt && n_rows > 0L && all(r_m == 0L) && all(r_inf)) {
    r_kmin <- rep(NA_integer_, n_rows)
    r_pool <- rep.int(1L, n_rows)
    pools <- list(n_pools = 1L, size = 1L, start = 1L, entries = as.integer(inf_idx), never_treated = TRUE)
    n_ukeys <- 1L
  } else {
    width <- nchar(as.character(nC))
    coh_fmt <- formatC(seq_len(nC), width = width, flag = "0")
    keys <- character(n_rows); r_kmin <- rep(NA_integer_, n_rows)
    for (r in seq_len(n_rows)) {
      fin <- later_finite(r_tp[r])
      fin <- fin[seq.int(length(fin) - r_m[r] + 1L, length.out = r_m[r])]
      if (r_m[r] > 0L) r_kmin[r] <- fin[1L]
      keys[r] <- paste(coh_fmt[c(fin, if (r_inf[r]) inf_idx)], collapse = ",")
    }
    ukeys <- sort(unique(keys), method = "radix")
    r_pool <- match(keys, ukeys)
    first_row <- match(ukeys, keys)
    pool_size <- r_m[first_row] + r_inf[first_row]
    pools <- list(n_pools = length(ukeys), size = as.integer(pool_size),
      start = catalog_pool_offsets(pool_size),
      entries = as.integer(unlist(strsplit(ukeys, ",", fixed = TRUE))),
      never_treated = r_inf[first_row])
    n_ukeys <- length(ukeys)
  }

  # Operators: one per (target, pool) under the declared conventions.
  r_gi <- tg[tp_j[r_tp]]
  catalog_product(nC, n_ukeys, "exact operator dictionary product")
  okey <- (as.double(r_gi) - 1) * n_ukeys + r_pool
  uo <- sort(unique(okey)); r_op <- match(okey, uo)
  op_first <- match(uo, okey)
  operators <- cdid_fast_df(list(operator_id = seq_along(uo), target_index = r_gi[op_first], g = coh[r_gi[op_first]],
    pool_id = r_pool[op_first],
    functional_key = paste0("cmpcat:", conventions, ":g", r_gi[op_first], ":K", r_pool[op_first], recycle0 = TRUE)))

  # Conditional endpoint roles; reference and pre regions follow r_g exactly.
  charge(4 * n_rows, "endpoint roles")
  r_ref <- tp_ref[r_tp]; r_s <- pair_s[tp_pair[r_tp]]; r_t <- pair_t[tp_pair[r_tp]]
  r_coh <- coh[r_gi]
  r_kmin_na <- is.na(r_kmin)
  role_labels <- c("post_private", "post_conditional", "anticipation_private", "anticipation_conditional")
  role_of <- function(u) {
    elig <- r_kmin_na | (r_kmin >= cstar[u]) # K inside E_u(delta)
    role <- role_labels[1L + elig + 2L * (cal[u] < r_coh)]
    role[u < r_ref] <- "pre_private"
    role[u == r_ref] <- "reference"
    role
  }
  s_role <- role_of(r_s); t_role <- role_of(r_t)

  # Coordinate dictionary: shared reference/conditional vertices per (g, u);
  # private vertices per (g, operator, u). Every referenced target is anchored.
  has_ref <- which(!is.na(ref_i))
  n_ref <- length(has_ref)
  charge(4 * (n_ref + 2 * n_rows), "coordinate interning")
  kind_of <- function(role) 1L + (role != "reference") + (role %in% c("pre_private", "anticipation_private", "post_private"))
  cj <- c(tg[has_ref], r_gi, r_gi)
  cu <- c(ref_i[has_ref], r_s, r_t)
  ck <- c(rep(1L, n_ref), kind_of(s_role), kind_of(t_role))
  co <- c(integer(n_ref), r_op, r_op) * (ck == 3L)
  oo <- order(cj, cu, ck, co, method = "radix"); nn <- length(oo)
  chg <- if (nn) c(TRUE, cj[oo[-1L]] != cj[oo[-nn]] | cu[oo[-1L]] != cu[oo[-nn]] |
                         ck[oo[-1L]] != ck[oo[-nn]] | co[oo[-1L]] != co[oo[-nn]]) else logical()
  cid <- integer(nn); cid[oo] <- cumsum(chg)
  keep <- oo[chg]
  ck_keep <- ck[keep]
  op_id_keep <- co[keep]
  op_id_keep[ck_keep != 3L] <- NA_integer_
  coordinates <- cdid_fast_df(list(coordinate_id = seq_along(keep), target_index = cj[keep], g = coh[cj[keep]],
    date_index = cu[keep], time = cal[cu[keep]],
    kind = c("reference", "shared_conditional", "private")[ck_keep],
    operator_id = op_id_keep))

  rows <- cdid_fast_df(list(row_id = seq_len(n_rows), target_index = r_gi, g = r_coh,
    s_index = r_s, t_index = r_t, s = cal[r_s], t = cal[r_t], ref_index = r_ref,
    operator_id = r_op, pool_id = r_pool, pool_size = as.integer(r_m + r_inf),
    original = r_orig, n_horizons = r_nh, horizon_first_index = r_h1, horizon_last_index = r_h2,
    s_role = s_role, t_role = t_role,
    s_coordinate_id = cid[n_ref + seq_len(n_rows)], t_coordinate_id = cid[n_ref + n_rows + seq_len(n_rows)]))

  tp_ne <- tabulate(cand_tp[nonempty], nbins = n_tp)
  tp_rows <- tabulate(r_tp, nbins = n_tp)
  orig_ne <- logical(n_tp); orig_ne[cand_tp[cand_h == 0L]] <- nonempty[cand_h == 0L]
  tp_status <- rep("retained", n_tp)
  tp_status[!orig_ne] <- "empty_pool"
  tp_status[is.na(tp_ref)] <- "no_admissible_reference"
  target_pairs <- cdid_fast_df(list(target_index = tg[tp_j], g = coh[tg[tp_j]],
    s_index = pair_s[tp_pair], t_index = pair_t[tp_pair], s = cal[pair_s[tp_pair]], t = cal[pair_t[tp_pair]],
    target_count = tp_count, n_other_supplied_cohorts = pair_len[tp_pair] - 1L,
    n_candidates = tp_nc, n_nonempty_candidates = tp_ne, n_rows = tp_rows,
    original_status = tp_status))
  j_status <- rep("reference_available", nJ)
  # Issue #16: a target read in a restricted eligibility context says so.
  if (!is.null(tg_ctx)) j_status[tg_ctx > 0L] <- "reference_available_eligibility_restricted"
  j_status[is.na(ref_i)] <- "no_admissible_reference"
  j_ncand <- numeric(nJ)
  for (ji in seq_len(nJ)) j_ncand[ji] <- sum(tp_nc[tp_j == ji])
  targets_out <- cdid_fast_df(list(target_index = tg, g = coh[tg], ref_index = ref_i, ref = cal[ref_i],
    status = j_status, n_target_pairs = tabulate(tp_j, nbins = nJ),
    n_candidates = j_ncand,
    n_nonempty_candidates = tabulate(tp_j[cand_tp[nonempty]], nbins = nJ),
    n_rows = tabulate(tp_j[r_tp], nbins = nJ)))

  # Requested post targets: catalog support only, never population identification.
  charge(4 * n_req, "requested target classification")
  qkey <- (as.double(req$gi) - 1) * nT + req$ti
  hit <- match(c((as.double(r_gi) - 1) * nT + r_s, (as.double(r_gi) - 1) * nT + r_t), qkey)
  erole <- c(s_role, t_role)
  n_touch <- tabulate(hit, nbins = n_req)
  n_cond <- tabulate(hit[erole == "post_conditional"], nbins = n_req)
  if (any(n_touch != n_cond)) {
    catalog_abort("Internal role invariant failed: a retained row touches a post date privately.", "cdid_catalog_internal_error")
  }
  tkey <- c((as.double(tg[tp_j]) - 1) * nT + pair_s[tp_pair], (as.double(tg[tp_j]) - 1) * nT + pair_t[tp_pair])
  reason <- rep("no_nonempty_control_pool", n_req)
  reason[!qkey %in% tkey] <- "no_target_pair_occupancy"
  if (!is.null(tg_ctx)) reason[!qkey %in% tkey & tg_ctx[match(req$gi, tg)] > 0L] <-
    "no_target_pair_occupancy_in_eligibility_context"
  reason[is.na(ref_i[match(req$gi, tg)])] <- "no_admissible_reference"
  sup <- n_cond > 0
  requested_supported <- cdid_fast_df(list(target_index = req$gi[sup], g = coh[req$gi[sup]], date_index = req$ti[sup],
    t = cal[req$ti[sup]], n_rows_touching = n_cond[sup]))
  requested_unsupported <- cdid_fast_df(list(target_index = req$gi[!sup], g = coh[req$gi[!sup]], date_index = req$ti[!sup],
    t = cal[req$ti[!sup]], reason = reason[!sup]))

  work_performed <- c(cutoff_comparisons = n_cmp, occupancy_entries = n_occ, reference_searches = nJ,
    candidate_evaluations = n_cand, dedup_entries = n_cand, role_evaluations = 2 * n_rows)
  if (n_cmp > nT + nC || sum(work_performed) > work_bound) {
    catalog_abort("Internal work accounting invariant failed.", "cdid_catalog_internal_error")
  }
  structure(list(
    spec = spec, control_group = control_group, conventions = conventions,
    targets = targets_out, target_pairs = target_pairs, rows = rows, pools = pools,
    operators = operators, coordinates = coordinates,
    requested_supported = requested_supported, requested_unsupported = requested_unsupported,
    resources = list(bounds = bnd, occupancy_rows = n_occ, requested_rows = n_req, positive_cells = n_pos, target_pairs = n_tp,
      candidates = n_cand, nonempty_candidates = n_ne, rows = n_rows,
      pool_entries_preintern = n_entries, pool_entries_interned = length(pools$entries),
      pools = pools$n_pools, operators = nrow(operators), coordinates = nrow(coordinates),
      scratch_charged = scratch, work_bound = work_bound, work_performed = work_performed,
      header_scratch = header_scratch, header_work = header_work, candidate_work = stage_work,
      work_definition = "Conservative scalar/vector-slot and character-byte visit envelope including validation, sort/search, repeated pair-pool scans, interning, roles and outputs; work_performed is only the listed completed loop counters, not total CPU instructions.",
      scratch_definition = "Conservative cumulative reserved atomic slots plus newly created string bytes; excludes caller-owned inputs, SEXP/list headers, allocator overhead and process RSS."),
    scope = paste("Private metadata catalog only. Conditional roles are not verified statistical assumptions;",
      "equal pools denote one population operator under the declared conventions, not equal fits or samples;",
      "zero occupancy is not population nonidentification; no estimator, IF, covariance or inference.")),
    class = "cdid_comparison_catalog")
}

# Issue #16: the restricted-context cells of `eligibility_occupancy`
# (list(target_context = context per supplied target, 0 = unrestricted;
# subtract = data.frame(context, cohort, s, t, count) of excluded units' pairs)).
# Every restricted context gets every owner cell with its excluded counts
# removed; a subtraction outside the owner occupancy or below zero is refused.
catalog_eligibility_cells <- function(eo, occ, ki, si, ti, coh, cal, targets, tg) {
  bad <- function(msg) catalog_abort(msg, c("cdid_catalog_occupancy_error", "cdid_catalog_input_error"))
  if (!is.list(eo) || !is.numeric(eo$target_context) || length(eo$target_context) != length(targets) ||
      anyNA(eo$target_context) || any(eo$target_context < 0 | eo$target_context != floor(eo$target_context)) ||
      !is.data.frame(eo$subtract) || !all(c("context", "cohort", "s", "t", "count") %in% names(eo$subtract)))
    bad("Invalid eligibility occupancy (target contexts and excluded-unit pair counts).")
  sub <- eo$subtract
  n_occ <- length(ki)
  n_ctx <- max(c(0L, as.integer(eo$target_context), as.integer(sub$context)))
  ctx_sorted <- as.integer(eo$target_context)[match(tg, match(as.double(targets), coh))]
  sk <- match(as.double(sub$cohort), coh); ss <- match(as.double(sub$s), cal); st <- match(as.double(sub$t), cal)
  if (anyNA(sk) || anyNA(ss) || anyNA(st) || anyNA(sub$count) || any(sub$count < 0) ||
      any(sub$context < 1L | sub$context > n_ctx))
    bad("Eligibility occupancy subtraction outside the cohort/calendar dictionaries.")
  base_key <- (as.double(si) * (length(cal) + 1) + ti) * (length(coh) + 1) + ki
  sub_cell <- match((as.double(ss) * (length(cal) + 1) + st) * (length(coh) + 1) + sk, base_key)
  if (anyNA(sub_cell)) bad("Eligibility occupancy subtracts a pair cell absent from the owner occupancy.")
  out_cnt <- rep(occ$count, n_ctx)
  hit <- (as.integer(sub$context) - 1L) * n_occ + sub_cell
  dec <- numeric(n_ctx * n_occ)
  dec_sum <- rowsum(as.double(sub$count), hit)
  dec[as.integer(rownames(dec_sum))] <- dec_sum[, 1L]
  out_cnt <- out_cnt - dec
  if (any(out_cnt < 0)) bad("Eligibility occupancy subtracts more units than the owner occupancy holds.")
  list(ki = rep(ki, n_ctx), si = rep(si, n_ctx), ti = rep(ti, n_ctx), count = out_cnt,
    context = rep(seq_len(n_ctx), each = n_occ), target_context = ctx_sorted)
}

# Reuse the explicit catalog controls for newly created bridge/design buffers.
# A class label is not a sampling or immutable-identity certificate. The helper
# validates finite dictionary/row/pool headers and derives the adapter's own
# work from the catalog's actual typed dimensions (CP08 scaling action 1).
#
# Typed dimensions: R rows, C cohorts, T dates, P interned pool entries in NP
# pools, PR = sum of the actual pool sizes of the R rows (row/pool visits),
# J targets, D declared qualifications, Q requested design grid cells, O
# operators, B convention bytes and Bf the longest operator key in bytes.
catalog_adapter_dimensions <- function(catalog, declared_qualified, target_post_times) {
  if (!inherits(catalog, "cdid_comparison_catalog") || !is.list(catalog) ||
      !is.list(catalog$resources) || !is.list(catalog$spec) ||
      !is.data.frame(catalog$rows) || !is.list(catalog$pools) ||
      !is.data.frame(catalog$targets) || !is.data.frame(catalog$operators))
    catalog_abort("Expected a catalog with actual dictionary, row, pool and resource headers.", "cdid_catalog_input_error")
  rows <- catalog$rows; pl <- catalog$pools
  R <- nrow(rows); C <- length(catalog$spec$cohorts); T <- length(catalog$spec$calendar)
  NP <- pl$n_pools; P <- length(pl$entries)
  if (!design_vector(NP, 1L) || NP < 0 || NP != floor(NP) || NP > .Machine$integer.max ||
      !is.numeric(pl$size) || length(pl$size) != NP || anyNA(pl$size) ||
      any(pl$size < 0 | pl$size != floor(pl$size) | pl$size > C) || length(pl$start) != NP ||
      !is.numeric(rows$pool_id) || length(rows$pool_id) != R || anyNA(rows$pool_id) ||
      any(rows$pool_id < 1 | rows$pool_id > NP | rows$pool_id != floor(rows$pool_id)))
    catalog_abort("Catalog pool/row index headers are malformed.", "cdid_catalog_input_error")
  catalog_product(NP, C, "interned pool size bound")
  if (catalog_size(c(sum(as.double(pl$size)), 0), "interned pool entries") != P)
    catalog_abort("Catalog pool sizes disagree with the interned entry vector.", "cdid_catalog_input_error")
  # Each actual row visits its own interned pool; R * C bounds the exact sum.
  catalog_product(R, C, "row/pool visit bound")
  PR <- sum(as.double(pl$size[rows$pool_id]))
  J <- nrow(catalog$targets); O <- nrow(catalog$operators)
  D <- if (is.data.frame(declared_qualified)) nrow(declared_qualified) else 0
  Q <- if (is.data.frame(target_post_times)) nrow(target_post_times) else
    catalog_product(J, length(target_post_times), "requested design grid")
  B <- if (is.character(catalog$conventions) && length(catalog$conventions) == 1L && !anyNA(catalog$conventions))
    as.double(nchar(catalog$conventions, type = "bytes")) else 0
  fk <- catalog$operators$functional_key
  if (O > 0L && (!is.character(fk) || length(fk) != O || anyNA(fk)))
    catalog_abort("Catalog operator keys are malformed.", "cdid_catalog_input_error")
  Bf <- if (O > 0L) max(as.double(nchar(fk, type = "bytes"))) else 0
  list(R = as.double(R), C = as.double(C), T = as.double(T), P = as.double(P), NP = as.double(NP),
    PR = PR, J = as.double(J), D = as.double(D), Q = as.double(Q), O = as.double(O), B = B, Bf = Bf)
}

# Conservative slot/byte envelope of the links bridge and, when `design`, of
# the accepted calendar design built from it (adaptation, validation twice,
# vertices, sorted CSR adjacency, both forests, labels and outputs), in the
# catalog's own units: work counts slot/byte visits, scratch counts reserved
# vector slots plus newly created string bytes plus eight header slots per
# created string. V = J + 2R + Q bounds vertices, E = R edges and A = J + 2R
# anchors (PT-all anchors every represented pre vertex); incidence_forest is
# linear in V + E + A (sorted adjacency, anchors grouped by component once).
# Itemized per typed element, the links stage allocates fewer than 20 slots
# per row (qualification keys and the 9-column table) and the design fewer
# than 200 per row and 90 per vertex (two validation passes, keys, vertex and
# parameter tables, two forests, post targets); the 64 slots charged per link
# stage element and 128 per design element of S2 (about 6 R plus V-scale
# terms) cover both at least twice.
# Created strings: pool text (at most 25 bytes per formatted entry plus
# delimiters), declared/qualification/tuple/vertex/time keys (at most 40 bytes
# each), labels (operator key plus at most 112 bytes) and link identifiers.
# Only actual per-cohort products enter: reference searches (J T), PT-all
# per-cohort vertex scans (J V), row-pool cleanliness (PR) and pool universe
# checks (R C), each exact comparison charged 64 visits. Applies equally to
# lower-bound dimensions, as in the early owner preflight before rows exist.
catalog_design_envelope <- function(d, design = FALSE) {
  lg <- function(n) 1 + ceiling(log2(n + 1))
  P2 <- function(a, b, label) catalog_product(a, b, label)
  pool_bytes <- catalog_size(c(P2(25, d$P, "formatted pool numbers"), P2(2, d$NP, "pool delimiters")), "pool text bytes")
  link_bytes <- catalog_size(c(pool_bytes, P2(24, d$D, "declared keys"), P2(48, d$R, "qualification keys")),
    "link string bytes")
  link_strings <- catalog_size(c(d$NP, d$D, P2(2, d$R, "qualification key strings")), "link strings")
  S1 <- catalog_size(c(d$R, d$C, d$T, d$P, d$NP, d$J, d$D, 1), "link stage dimensions")
  work <- catalog_size(c(P2(512 * lg(S1), S1, "link validation/key/pool work"), P2(64, d$D, "declared cutoff checks"),
    P2(16, link_bytes, "link string work")), "adapter link work")
  scratch <- catalog_size(c(P2(64, S1, "link vectors and table"), link_bytes, P2(8, link_strings, "link string headers")),
    "adapter link scratch")
  if (!isTRUE(design)) return(list(scratch = scratch, work = work, string_bytes = link_bytes, linear = S1))
  V <- catalog_size(c(d$J, P2(2, d$R, "link endpoints"), d$Q), "design vertex bound")
  A <- catalog_size(c(d$J, P2(2, d$R, "pre anchors")), "design anchor bound")
  S2 <- catalog_size(c(S1, d$PR, d$Q, V, d$R, A), "design dimensions")
  design_bytes <- catalog_size(c(P2(2, pool_bytes, "two validation parses"),
    P2(40, catalog_size(c(P2(3, V, "vertex/time/label keys"), P2(5, d$R, "endpoint and tuple keys")), "key strings"), "key bytes"),
    P2(catalog_size(c(d$Bf, 112), "label width"), V, "vertex labels"), P2(16, d$R, "link identifiers")),
    "design string bytes")
  design_strings <- catalog_size(c(P2(4, V, "vertex strings"), P2(6, d$R, "link strings"), P2(2, d$NP, "pool parses")),
    "design strings")
  comparisons <- catalog_size(c(P2(3, P2(d$J, d$T, "reference searches"), "three reference passes"),
    P2(d$J, V, "PT-all per-cohort vertex scans"), P2(8, d$PR, "row-pool cleanliness checks"),
    P2(12, d$R, "endpoint cutoff checks"), P2(d$R, d$C, "pool universe checks")), "design exact comparisons")
  work <- catalog_size(c(work, P2(512 * lg(S2), S2, "design validation/sort/adjacency/forest/output work"),
    P2(64, comparisons, "design cutoff comparison work"), P2(16, design_bytes, "design string work")),
    "adapter design work")
  scratch <- catalog_size(c(scratch, P2(128, S2, "design vectors, CSR, forests and outputs"),
    design_bytes, P2(8, design_strings, "design string headers")), "adapter design scratch")
  list(scratch = scratch, work = work, string_bytes = catalog_size(c(link_bytes, design_bytes), "adapter string bytes"),
    linear = S2, vertices = V, anchors = A, comparisons = comparisons)
}

catalog_adapter_preflight <- function(catalog, declared_qualified, target_post_times = NULL, design = FALSE) {
  d <- catalog_adapter_dimensions(catalog, declared_qualified, target_post_times)
  b <- validate_catalog_bounds(catalog$resources$bounds)
  env <- catalog_design_envelope(d, design = design)
  if (env$scratch > b$max_scratch || env$work > b$max_work)
    catalog_abort_recorded("Catalog adapter output/string/graph work exceeds its explicit stored controls before allocation.", "cdid_catalog_bound_error",
      diagnostics = tryCatch(catalog_bound_diagnostics(env, b), error = function(e) NULL))
  invisible(list(scratch = env$scratch, work = env$work, dimensions = d, envelope = env))
}

# Declared-qualification bridge to the accepted calendar design. A declared
# (g, u) is the caller's statement that the population premises hold from the
# actual reference; only conditionally eligible endpoints can use it.
comparison_catalog_design_links <- function(catalog, declared_qualified) {
  catalog_adapter_preflight(catalog, declared_qualified)
  catalog_design_handoff(catalog, declared_qualified)$links
}

# CP08 validated catalog-to-design handoff (the body of the strict bridge
# above). One complete semantic validation of the declared qualification
# packet, then the links table exactly as the bridge returns it, together with
# the catalog's canonical calendar/cohorts/anticipation and the typed pools:
# every interned pool is expanded once into its cohort set, rendered once into
# the `control_pool` text, and both are handed to the private design core
# (`build_catalog_handoff_design`) so the design parses no pool text and
# derives no reference/qualification a second time. The packet is bound to
# this invocation's catalog object; it certifies nothing and is never stored.
# Callers run `catalog_adapter_preflight` on the same catalog first.
catalog_design_handoff <- function(catalog, declared_qualified) {
  if (!inherits(catalog, "cdid_comparison_catalog")) {
    catalog_abort("Expected a private comparison catalog.", "cdid_catalog_input_error")
  }
  spec <- catalog$spec; cal <- spec$calendar; coh <- spec$cohorts; delta <- spec$anticipation
  rows <- catalog$rows
  if (nrow(rows) == 0L) catalog_abort("The catalog retains no rows to declare as design links.", "cdid_catalog_input_error")
  decl <- numeric()
  n_cal <- length(cal)
  if (!is.null(declared_qualified)) {
    dq <- declared_qualified
    if (!is.data.frame(dq) || ncol(dq) != 2L || anyDuplicated(names(dq)) || !setequal(names(dq), c("g", "t")) ||
        !design_vector(dq$g, nrow(dq)) || !design_vector(dq$t, nrow(dq))) {
      catalog_abort("`declared_qualified` must be NULL or a data.frame with finite numeric columns g and t.",
        "cdid_catalog_qualification_error")
    }
    gi <- match(as.double(dq$g), coh); ui <- match(as.double(dq$t), cal)
    j <- match(gi, catalog$targets$target_index)
    if (anyNA(j) || anyNA(catalog$targets$ref_index[j]) || anyNA(ui)) {
      catalog_abort("Declared qualifications need a referenced catalog target and a stored calendar date.",
        "cdid_catalog_qualification_error")
    }
    can_fast_dq <- exists("cdid_are_default_fns", mode = "function", inherits = TRUE) &&
      cdid_are_default_fns("exact_cmp_sum") &&
      is.finite(delta) && abs(delta) <= 1e14 &&
      all(is.finite(cal[ui])) && all(abs(cal[ui]) <= 1e14) &&
      all(is.finite(coh[gi])) && all(abs(coh[gi]) <= 1e14)
    bad_causal <- if (can_fast_dq) {
      any(cal[ui] + delta < coh[gi])
    } else {
      any(vapply(seq_along(gi), function(i) exact_cmp_sum(cal[ui[i]], delta, coh[gi[i]]) < 0L, logical(1L)))
    }
    if (bad_causal) {
      catalog_abort("Only causal endpoints (u + anticipation >= g) can be declared qualified; earlier coordinates stay private.",
        "cdid_catalog_qualification_error")
    }
    decl <- (as.double(gi) - 1) * n_cal + ui
    if (anyDuplicated(decl)) {
      catalog_abort("Duplicate declared qualification.", c("cdid_catalog_duplicate_error", "cdid_catalog_qualification_error"))
    }
  }
  qual <- function(role, u) {
    if (!length(decl)) return(role == "reference")
    role == "reference" |
      (role %in% c("anticipation_conditional", "post_conditional") &
        ((as.double(rows$target_index) - 1) * n_cal + u) %in% decl)
  }
  pl <- catalog$pools
  # One expansion per interned pool: the typed cohort set and its rendering.
  pool_values <- lapply(seq_len(pl$n_pools), function(p) {
    as.double(coh[pl$entries[seq.int(pl$start[p], length.out = pl$size[p])]])
  })
  pool_text <- vapply(pool_values, function(vals) {
    # %.17g round-trips doubles exactly; the standalone design path parses it back.
    paste0("{", paste(sprintf("%.17g", vals), collapse = ","), "}")
  }, character(1L))
  pool_index <- as.integer(rows$pool_id)
  links <- cdid_fast_df(list(g = rows$g, s = rows$s, t = rows$t, ref = cal[rows$ref_index],
    functional_key = catalog$operators$functional_key[rows$operator_id],
    control_pool = pool_text[pool_index],
    s_qualified = qual(rows$s_role, rows$s_index), t_qualified = qual(rows$t_role, rows$t_index),
    catalog_row = rows$row_id))
  list(links = links, calendar = cal, cohorts = coh, anticipation = delta,
    pools = list(index = pool_index, values = pool_values, text = pool_text))
}

build_comparison_catalog_design <- function(catalog, declared_qualified, target_post_times,
                                             pt_assumption = "post", pt_start = NULL) {
  catalog_adapter_preflight(catalog, declared_qualified, target_post_times, design = TRUE)
  links <- comparison_catalog_design_links(catalog, declared_qualified)
  build_calendar_chained_design(links, calendar = catalog$spec$calendar, cohorts = catalog$spec$cohorts,
    anticipation = catalog$spec$anticipation, target_post_times = target_post_times,
    pt_assumption = pt_assumption, pt_start = pt_start)
}
