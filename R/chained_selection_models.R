# Private additive probability-model backend for CP04 / D48 Step 2 Section A.
#
# Supports:
#   * selection_model in c("pair_logit", "independent_period_logit")
#   * selection_pooling in c("separate", "within_cohort", "across_cohorts")
#     ("separate" is the default for BOTH models)
#   * Full-B marginal date indicator reading (cp04_date_rows) on the sealed
#     cdid_cp04_frame_reader owner (preserving singletons in B \ A with a_index == 0
#     and keeping linked zero-observation panel units and F \ B units at 0).
#   * Formula intercept metadata handling:
#       - standard intercept + slopes (~ 1 + z1 + ...): local cell intercepts + shared slopes
#       - intercept-only (~ 1): local cell intercepts only (d_s = 0)
#       - homogeneous no-intercept (~ 0 + z1 + ...): shared slopes only, uncentered, no cell intercepts
#       - zero-column (~ 0 or ~ -1): eta = 0, specified marginal probability p = 1/2
#   * Compact block IRLS and profiled Fisher information via shared-WLS QR/Cholesky
#     primitives (or_pooled_wls_* / or_wls_*), preserving cp04_selection_fit / cp04_logit
#     on the separate positive-column pair_logit path.
#   * Original-unit score summation before Fisher-information inversion, pair-probability
#     and log-weight moment contractions, and bounded chunked IF replay / covariance streaming.

selection_model_design <- function(formula, data) {
  if (!inherits(formula, "formula") || !is.data.frame(data) || nrow(data) < 1L)
    cp04_abort("Supply a valid selection formula and nonempty data.frame.", "input")
  tt <- tryCatch(stats::terms(formula, data = data), error = function(e) {
    cp04_abort(paste0("Invalid selection formula: ", conditionMessage(e)), "input")
  })
  has_intercept <- identical(as.integer(attr(tt, "intercept")), 1L)
  mm <- tryCatch(stats::model.matrix(tt, data = data), error = function(e) {
    cp04_abort(paste0("Failed to construct selection design matrix: ", conditionMessage(e)), "input")
  })
  if (!is.matrix(mm) || !is.numeric(mm) || nrow(mm) != nrow(data))
    cp04_abort("Selection model.matrix row count disagrees with supplied data.", "input")
  Z <- matrix(as.double(mm), nrow = nrow(mm), ncol = ncol(mm),
              dimnames = list(NULL, colnames(mm)))
  attr(Z, "has_intercept") <- has_intercept
  if (!is.null(attr(mm, "assign"))) {
    attr(Z, "assign") <- as.integer(attr(mm, "assign"))
  }
  Z
}

extract_selection_design_spec <- function(Z, U = NULL, has_intercept = NULL,
                                           chunk_rows = nrow(Z)) {
  if (!is.matrix(Z) || is.object(Z) || !is.numeric(Z) || nrow(Z) < 1L || ncol(Z) < 0L)
    cp04_abort("Supply a numeric selection design matrix Z with nrow >= 1 and ncol >= 0.", "input")
  N <- nrow(Z)
  d <- ncol(Z)
  if (is.null(U)) U <- rep.int(TRUE, N)
  if (!is.logical(U) || is.object(U) || !is.null(dim(U)) || length(U) != N || anyNA(U) || !any(U))
    cp04_abort("Invalid fitting roster U for selection design.", "input")
  chunk_rows <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  for (start in seq.int(1L, N, by = chunk_rows)) {
    rr <- seq.int(start, min(N, start + chunk_rows - 1L))
    rr <- rr[U[rr]]
    cp04_finite(Z[rr, , drop = FALSE], "selection design on U")
  }

  is_intercept_column <- function(j) {
    for (start in seq.int(1L, N, by = chunk_rows)) {
      rr <- seq.int(start, min(N, start + chunk_rows - 1L))
      rr <- rr[U[rr]]
      if (any(Z[rr, j] != 1)) return(FALSE)
    }
    TRUE
  }
  hi <- NULL
  if (!is.null(has_intercept)) {
    if (!is.logical(has_intercept) || is.object(has_intercept) || !is.null(dim(has_intercept)) ||
        length(has_intercept) != 1L || is.na(has_intercept))
      cp04_abort("`has_intercept` must be a non-NA logical scalar.", "input")
    hi <- has_intercept
  } else if (!is.null(attr(Z, "has_intercept"))) {
    attr_hi <- attr(Z, "has_intercept")
    if (!is.logical(attr_hi) || is.object(attr_hi) || !is.null(dim(attr_hi)) ||
        length(attr_hi) != 1L || is.na(attr_hi))
      cp04_abort("Attribute `has_intercept` on Z must be a non-NA logical scalar.", "input")
    hi <- attr_hi
  } else if (!is.null(attr(Z, "intercept"))) {
    attr_int <- attr(Z, "intercept")
    if (!(is.numeric(attr_int) || is.logical(attr_int)) || is.object(attr_int) ||
        !is.null(dim(attr_int)) || length(attr_int) != 1L || is.na(attr_int) || !attr_int %in% c(0, 1))
      cp04_abort("Attribute `intercept` on Z must be 0 or 1.", "input")
    hi <- (as.integer(attr_int) == 1L)
  } else if (!is.null(attr(Z, "assign"))) {
    asg <- attr(Z, "assign")
    if (!is.numeric(asg) || length(asg) != d || anyNA(asg))
      cp04_abort("Attribute `assign` on Z must match ncol(Z).", "input")
    hi <- (0L %in% as.integer(asg))
  } else {
    hi <- if (d == 0L) FALSE else is_intercept_column(1L)
  }

  if (isTRUE(hi)) {
    if (d < 1L)
      cp04_abort("Intercept selection design requires at least one column.", "input")
    asg <- attr(Z, "assign")
    int_col <- if (!is.null(asg) && any(as.integer(asg) == 0L)) {
      which(as.integer(asg) == 0L)[1L]
    } else {
      1L
    }
    if (!is_intercept_column(int_col))
      cp04_abort("Declared intercept column of selection design must equal 1 on all U rows.", "input")
    slope_cols <- setdiff(seq_len(d), int_col)
    d_s <- length(slope_cols)
    branch <- if (d_s > 0L) "intercept_slope" else "intercept_only"
  } else {
    int_col <- integer(0)
    slope_cols <- seq_len(d)
    d_s <- d
    branch <- if (d_s > 0L) "no_intercept" else "zero_col"
  }

  cn <- colnames(Z)
  slope_names <- if (d_s > 0L && !is.null(cn) && all(nzchar(cn[slope_cols]))) {
    cn[slope_cols]
  } else if (d_s > 0L) {
    paste0("z", seq_len(d_s))
  } else {
    character(0)
  }

  list(
    N = N,
    d = d,
    has_intercept = isTRUE(hi),
    int_col = as.integer(int_col),
    slope_cols = as.integer(slope_cols),
    d_s = as.integer(d_s),
    branch = branch,
    slope_names = slope_names
  )
}

# Full-B marginal date indicator reader on the sealed CP04 full-F reader.
# Unlike cp04_pair_rows (which filters a_index > 0), cp04_date_rows also queries
# singleton panel units in B \ A (where a_index == 0 and panel$units$n_obs == 1)
# via sorted unit-ID lookup, while linked zero-observer panel units (n_obs == 0)
# and frame-only units (F \ B) remain 0. Completely outcome-free (never reads y).
cp04_date_rows <- function(reader, rows, date_spec) {
  validate_cp04_reader(reader)
  st <- reader$state
  ctrl <- reader$control
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows != floor(rows)) || any(rows < 1L | rows > length(st$a_index)) ||
      length(rows) > ctrl$max_chunk_rows)
    cp04_abort("Invalid row chunk for cp04_date_rows.", "input")
  rows <- as.integer(rows)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)

  dt <- NULL
  t_val <- NULL
  if (is.list(date_spec)) {
    dt <- date_spec$calendar_index
    if (!is.null(date_spec$t)) t_val <- as.double(date_spec$t)
    if (is.null(dt) && !is.null(t_val) && length(t_val) == 1L && is.finite(t_val)) {
      dt <- match(t_val, panel$all_times)
    }
  } else if (is.numeric(date_spec) && !is.object(date_spec) && is.null(dim(date_spec)) &&
             length(date_spec) == 1L && is.finite(date_spec)) {
    t_val <- as.double(date_spec)
    dt <- match(t_val, panel$all_times)
  } else {
    cp04_abort("Invalid date specification for cp04_date_rows.", "input")
  }
  if (is.null(dt) || length(dt) != 1L || is.na(dt) || !is.finite(dt) ||
      dt != floor(dt) || dt < 1L || dt > length(panel$all_times))
    cp04_abort("Requested selection date is not on the panel calendar.", "support")
  dt <- as.integer(dt)
  if (is.null(t_val) || length(t_val) != 1L || !is.finite(t_val)) {
    t_val <- panel$all_times[dt]
  } else if (panel$all_times[dt] != t_val) {
    cp04_abort("Date calendar_index disagrees with date value t.", "identity")
  }

  n <- length(rows)
  R <- numeric(n)
  ai <- st$a_index[rows]

  # 1. Repeat-observation units in A (n_obs >= 2): binary search on date postings.
  keep_a <- which(ai > 0L)
  post <- panel$postings
  if (length(keep_a) > 0L && post$date_count[dt] > 0L) {
    mt <- ipw_owner_match_sorted_double(
      ai[keep_a], post$post_a_index, post$date_start[dt], post$date_end[dt]
    )
    hit_a <- keep_a[!is.na(mt)]
    if (length(hit_a) > 0L) R[hit_a] <- 1
  }

  # 2. Off-A units (ai == 0): singletons in B \ A have n_obs == 1 at obs_start[u].
  # Linked zero-observation panel units (n_obs == 0) and frame-only units (is.na(u)) stay 0.
  off_a <- which(ai == 0L)
  if (length(off_a) > 0L && isTRUE(panel$n_singleton_units > 0L)) {
    id_sub <- st$id[rows[off_a]]
    id_key <- if (is.numeric(id_sub)) {
      id_d <- as.double(id_sub)
      id_d[id_d == 0.0] <- 0.0
      id_d
    } else {
      as.character(id_sub)
    }
    u <- ipw_owner_match_sorted_id(id_key, panel$unit_id_key, panel$id_type)
    is_sing <- which(!is.na(u) & panel$units$n_obs[u] == 1L)
    if (length(is_sing) > 0L) {
      obs_pos <- panel$units$obs_start[u[is_sing]]
      cal_idx <- panel$observations$calendar_index[obs_pos]
      hit_s <- off_a[is_sing[cal_idx == dt]]
      if (length(hit_s) > 0L) R[hit_s] <- 1
    }
  }

  list(R = R, calendar_index = dt, t = t_val)
}

validate_selection_mode_args <- function(selection_model, selection_pooling, universe = "F", U = NULL, rho = NULL) {
  if (!is.character(selection_model) || is.object(selection_model) || !is.null(dim(selection_model)) ||
      length(selection_model) != 1L || is.na(selection_model) ||
      !selection_model %in% c("pair_logit", "independent_period_logit"))
    cp04_abort("`selection_model` must be \"pair_logit\" or \"independent_period_logit\".", "unsupported")
  if (!is.character(selection_pooling) || is.object(selection_pooling) || !is.null(dim(selection_pooling)) ||
      length(selection_pooling) != 1L || is.na(selection_pooling) ||
      !selection_pooling %in% c("separate", "within_cohort", "across_cohorts"))
    cp04_abort("`selection_pooling` must be \"separate\", \"within_cohort\", or \"across_cohorts\".", "unsupported")
  if (identical(selection_model, "independent_period_logit")) {
    if (!identical(universe, "F") || (!is.null(U) && !all(U)) || (!is.null(rho) && any(rho != 1))) {
      cp04_abort(
        "`independent_period_logit` is supported only on the full frame U = F with rho = 1; Conditioning on A or B destroys marginal date independence.",
        "unsupported"
      )
    }
  } else {
    if (!is.character(universe) || length(universe) != 1L || is.na(universe) || !universe %in% c("F", "B", "A"))
      cp04_abort("`pair_logit` requires fitting universe in c(\"F\", \"B\", \"A\").", "unsupported")
  }
  invisible(TRUE)
}

# Pure deterministic catalog builder: deduplicates and canonically orders all
# active pair-cohort targets (e, k) and active selection cells (date cells (k, t)
# or pair cells (e, k)), then partitions cells into blocks B_1, ..., B_M.
# `families` (D60 recipe-cohort layout): recipe of each layout cohort; with it,
# across_cohorts pooling forms one block per recipe.
partition_selection_cells <- function(pair_specs, cohorts, all_times = NULL,
                                      selection_model = "pair_logit",
                                      selection_pooling = "separate", families = NULL) {
  validate_selection_mode_args(selection_model, selection_pooling)
  if (!is.list(pair_specs) || !length(pair_specs))
    cp04_abort("Supply a nonempty list of pair-cohort specifications.", "input")

  n_req <- length(pair_specs)
  req_k <- integer(n_req)
  req_ds <- integer(n_req)
  req_dt <- integer(n_req)
  req_s <- numeric(n_req)
  req_t <- numeric(n_req)
  # D60: a pair spec may carry its recipe. Cells, pair-cohorts and pooling
  # blocks are then keyed by recipe; pooling partitions cells inside one recipe
  # (no cross-recipe slope). Without recipes every key and order is unchanged.
  has_ref <- any(vapply(pair_specs, function(sp) is.list(sp) && !is.null(sp$reference), logical(1)))
  req_r <- integer(n_req)
  rtag <- function(r) if (has_ref) paste0(":r", r) else ""

  for (i in seq_len(n_req)) {
    sp <- pair_specs[[i]]
    if (!is.list(sp) || is.null(sp$cohort) || !is.list(sp$pair))
      cp04_abort("Each pair_spec must provide `cohort` index and `pair` descriptor.", "input")
    k <- as.integer(sp$cohort)
    if (length(k) != 1L || is.na(k) || k < 1L || k > length(cohorts))
      cp04_abort("Invalid cohort index in pair_spec.", "support")
    pr <- sp$pair
    s_val <- as.double(pr$s)
    t_val <- as.double(pr$t)
    ds <- pr$s_calendar_index
    dt <- pr$t_calendar_index
    if (is.null(ds) && !is.null(all_times)) ds <- match(s_val, all_times)
    if (is.null(dt) && !is.null(all_times)) dt <- match(t_val, all_times)
    if (is.null(ds)) ds <- as.integer(s_val)
    if (is.null(dt)) dt <- as.integer(t_val)
    ds <- as.integer(ds)
    dt <- as.integer(dt)
    if (length(s_val) != 1L || length(t_val) != 1L || !is.finite(s_val) || !is.finite(t_val) ||
        !(s_val < t_val) || length(ds) != 1L || length(dt) != 1L || is.na(ds) || is.na(dt) ||
        ds < 1L || dt < 1L || ds >= dt)
      cp04_abort("Invalid pair endpoints (s, t) in selection catalog.", "input")
    req_k[i] <- k
    req_ds[i] <- ds
    req_dt[i] <- dt
    req_s[i] <- s_val
    req_t[i] <- t_val
    if (has_ref) {
      r <- as.integer(sp$reference)
      if (length(r) != 1L || is.na(r) || r < 1L) cp04_abort("Invalid recipe in pair_spec.", "input")
      req_r[i] <- r
    }
  }

  # Canonical deduplication of active pair-cohorts (k, ds, dt[, recipe])
  pc_keys <- paste0(req_k, ":", req_ds, ":", req_dt, vapply(req_r, rtag, character(1)))
  first_idx <- which(!duplicated(pc_keys))
  ord_pc <- first_idx[order(req_k[first_idx], req_ds[first_idx], req_dt[first_idx], req_r[first_idx], method = "radix")]
  n_pc <- length(ord_pc)
  pair_cohorts <- vector("list", n_pc)
  canon_pc_keys <- character(n_pc)
  for (q in seq_len(n_pc)) {
    ii <- ord_pc[q]
    k <- req_k[ii]
    ds <- req_ds[ii]
    dt <- req_dt[ii]
    s_val <- req_s[ii]
    t_val <- req_t[ii]
    canon_pc_keys[q] <- paste0(k, ":", ds, ":", dt, rtag(req_r[ii]))
    pair_cohorts[[q]] <- list(
      pair_cohort_id = q,
      cohort = k,
      cohort_value = cohorts[k],
      pair = list(s = s_val, t = t_val, s_calendar_index = ds, t_calendar_index = dt),
      label = paste0("k", k, "_s", ds, "_t", dt, if (has_ref) paste0("_r", req_r[ii]) else "")
    )
    if (has_ref) pair_cohorts[[q]]$reference <- req_r[ii]
  }
  request_to_pc <- match(pc_keys, canon_pc_keys)

  # Build deduplicated, canonically ordered active selection cells
  pc_ref <- if (has_ref) vapply(pair_cohorts, `[[`, integer(1), "reference") else integer(n_pc)
  if (identical(selection_model, "independent_period_logit")) {
    raw_ck <- c(vapply(pair_cohorts, `[[`, integer(1), "cohort"),
                vapply(pair_cohorts, `[[`, integer(1), "cohort"))
    raw_dt <- c(vapply(pair_cohorts, function(x) x$pair$s_calendar_index, integer(1)),
                vapply(pair_cohorts, function(x) x$pair$t_calendar_index, integer(1)))
    raw_t <- c(vapply(pair_cohorts, function(x) x$pair$s, numeric(1)),
               vapply(pair_cohorts, function(x) x$pair$t, numeric(1)))
    raw_r <- c(pc_ref, pc_ref)
    ck_keys <- paste0(raw_ck, ":", raw_dt, vapply(raw_r, rtag, character(1)))
    keep_c <- which(!duplicated(ck_keys))
    ord_c <- keep_c[order(raw_ck[keep_c], raw_dt[keep_c], raw_r[keep_c], method = "radix")]
    C <- length(ord_c)
    cells <- vector("list", C)
    canon_c_keys <- character(C)
    for (c_id in seq_len(C)) {
      jj <- ord_c[c_id]
      k <- raw_ck[jj]
      dt <- raw_dt[jj]
      t_val <- raw_t[jj]
      canon_c_keys[c_id] <- paste0(k, ":", dt, rtag(raw_r[jj]))
      cells[[c_id]] <- list(
        cell_id = c_id,
        type = "date",
        cohort = k,
        cohort_value = cohorts[k],
        t = t_val,
        calendar_index = dt,
        date_spec = list(t = t_val, calendar_index = dt),
        label = paste0("k", k, "_t", dt, if (has_ref) paste0("_r", raw_r[jj]) else "")
      )
      if (has_ref) cells[[c_id]]$reference <- raw_r[jj]
    }
  } else {
    C <- n_pc
    cells <- vector("list", C)
    canon_c_keys <- canon_pc_keys
    for (c_id in seq_len(C)) {
      pc <- pair_cohorts[[c_id]]
      cells[[c_id]] <- list(
        cell_id = c_id,
        type = "pair",
        cohort = pc$cohort,
        cohort_value = pc$cohort_value,
        pair = pc$pair,
        label = pc$label
      )
      if (has_ref) cells[[c_id]]$reference <- pc$reference
    }
  }

  # Partition cells into blocks B_1, ..., B_M (inside one recipe when recipes exist)
  cell_cohorts <- vapply(cells, `[[`, integer(1), "cohort")
  cell_ref <- if (has_ref) vapply(cells, `[[`, integer(1), "reference") else integer(C)
  if (identical(selection_pooling, "separate")) {
    M <- C
    block_specs <- lapply(seq_len(M), function(m) {
      list(
        block_id = m,
        cell_ids = m,
        cohorts = cell_cohorts[m],
        B_m = 1L,
        label = paste0("separate_", cells[[m]]$label)
      )
    })
  } else if (identical(selection_pooling, "within_cohort")) {
    kr <- unique(data.frame(k = cell_cohorts, r = cell_ref))
    kr <- kr[order(kr$k, kr$r, method = "radix"), , drop = FALSE]
    M <- nrow(kr)
    block_specs <- lapply(seq_len(M), function(m) {
      k <- kr$k[m]
      cids <- which(cell_cohorts == k & cell_ref == kr$r[m])
      list(
        block_id = m,
        cell_ids = as.integer(cids),
        cohorts = k,
        B_m = length(cids),
        label = paste0("within_cohort_k", k, if (has_ref) paste0("_r", kr$r[m]) else "")
      )
    })
  } else {
    fam <- if (is.null(families)) cell_ref else as.integer(families[cell_cohorts])
    rr <- sort(unique(fam))
    M <- length(rr)
    block_specs <- lapply(seq_len(M), function(m) {
      cids <- which(fam == rr[m])
      list(
        block_id = m,
        cell_ids = as.integer(cids),
        cohorts = sort(unique(cell_cohorts[cids]), method = "radix"),
        B_m = length(cids),
        label = if (has_ref || !is.null(families)) paste0("across_cohorts_r", rr[m]) else "across_cohorts_all"
      )
    })
  }
  if (has_ref) for (m in seq_len(M)) block_specs[[m]]$reference <- cell_ref[block_specs[[m]]$cell_ids[1L]]

  for (m in seq_len(M)) {
    cids <- block_specs[[m]]$cell_ids
    for (b in seq_along(cids)) {
      cid <- cids[b]
      cells[[cid]]$block_id <- m
      cells[[cid]]$local_b <- as.integer(b)
    }
  }

  for (q in seq_len(n_pc)) {
    pc <- pair_cohorts[[q]]
    if (identical(selection_model, "independent_period_logit")) {
      sc <- match(paste0(pc$cohort, ":", pc$pair$s_calendar_index, rtag(pc_ref[q])), canon_c_keys)
      tc <- match(paste0(pc$cohort, ":", pc$pair$t_calendar_index, rtag(pc_ref[q])), canon_c_keys)
      pc$s_cell_id <- as.integer(sc)
      pc$t_cell_id <- as.integer(tc)
      pc$s_block_id <- cells[[sc]]$block_id
      pc$s_local_b <- cells[[sc]]$local_b
      pc$t_block_id <- cells[[tc]]$block_id
      pc$t_local_b <- cells[[tc]]$local_b
      pc$cell_ids <- c(as.integer(sc), as.integer(tc))
      pc$block_ids <- unique(c(cells[[sc]]$block_id, cells[[tc]]$block_id))
    } else {
      pc$pair_cell_id <- as.integer(q)
      pc$block_id <- cells[[q]]$block_id
      pc$local_b <- cells[[q]]$local_b
      pc$cell_ids <- as.integer(q)
      pc$block_ids <- cells[[q]]$block_id
    }
    pair_cohorts[[q]] <- pc
  }

  list(
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    cells = cells,
    block_specs = block_specs,
    pair_cohorts = pair_cohorts,
    request_to_pc = as.integer(request_to_pc)
  )
}

build_selection_block_catalog <- function(reader, link_views,
                                          selection_model = "pair_logit",
                                          selection_pooling = "separate",
                                          control = reader$control, reference = NULL,
                                          view_cohorts = NULL, families = NULL) {
  validate_cp04_reader(reader)
  ctrl <- validate_ipw_owner_control(control)
  st <- reader$state
  validate_selection_mode_args(selection_model, selection_pooling, st$universe, st$U, st$rho)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  if (!is.list(link_views) || !length(link_views))
    cp04_abort("Supply actual private link views.", "input")

  L <- length(link_views)
  raw_specs <- list()
  link_cohorts <- vector("list", L)
  link_pair <- vector("list", L)
  link_ids <- character(L)
  source_keys <- character(L)
  n_obs_vec <- integer(L)

  for (ell in seq_len(L)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    sk <- v$functional_key
    if (!gmm_scalar_character(sk) || !startsWith(sk, "cmpcat:"))
      cp04_abort("Selection models require the supplied comparison-catalog functional key.", "unsupported")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      cp04_abort("Legacy nuisance candidates are unsupported for corrected selection.", "unsupported")
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
        cp04_abort("Cross-pair or unverified legacy nuisance view sharing is unsupported.", "unsupported")
    }
    ks <- sort(v$cohorts$K, method = "radix")
    # A recipe-cohort layout supplies each view's layout cohort indices.
    cohorts <- if (is.null(view_cohorts)) match(c(v$cohorts$g, ks), st$cohorts) else as.integer(view_cohorts[[ell]])
    if (anyNA(cohorts) || any(v$support$by_cohort_pair_count < 1L) || v$support$n_observed_risk_units < 1L)
      cp04_abort("Unsupported full-F cohort or pair support.", "support")
    r <- if (is.null(reference)) NULL else reference$view_recipe[ell]
    for (k in cohorts) {
      a_k <- if (is.null(reference)) st$a[k] else chained_reference_stat(reference, r, k, "a")
      if (!is.finite(a_k) || a_k <= 0)
        cp04_abort("A relevant cohort has no repeat-observation mass.", "support")
      raw_specs[[length(raw_specs) + 1L]] <- if (is.null(reference)) list(cohort = k, pair = v$pair) else
        list(cohort = k, pair = v$pair, reference = r)
    }
    link_cohorts[[ell]] <- as.integer(cohorts)
    link_pair[[ell]] <- v$pair
    link_ids[ell] <- if (is.null(v$link_id)) "" else as.character(v$link_id)
    source_keys[ell] <- sk
    n_obs_vec[ell] <- as.integer(v$support$n_observed_risk_units)
  }
  if (any(!nzchar(link_ids)) || anyDuplicated(link_ids))
    cp04_abort("Distinct nonempty link IDs required.", "input")

  cat_obj <- partition_selection_cells(
    pair_specs = raw_specs,
    cohorts = st$cohorts,
    all_times = panel$all_times,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    families = families
  )

  # Map each link view to its ordered pair-cohort IDs
  pos <- 0L
  links <- vector("list", L)
  law_tag <- paste0(
    "cmpcat:CP04-", st$universe, "-weighted-",
    selection_model, "-", selection_pooling, "-a/p-etaQG|"
  )
  for (ell in seq_len(L)) {
    m_c <- length(link_cohorts[[ell]])
    pc_ids <- cat_obj$request_to_pc[seq.int(pos + 1L, pos + m_c)]
    pos <- pos + m_c
    links[[ell]] <- list(
      link_id = link_ids[ell],
      cohorts = link_cohorts[[ell]],
      pair = link_pair[[ell]],
      pair_cohort_ids = as.integer(pc_ids),
      n_observers = n_obs_vec[ell],
      functional_key = paste0(law_tag, substring(source_keys[ell], 8L)),
      source_functional_key = source_keys[ell]
    )
    if (!is.null(reference)) links[[ell]]$reference <- reference$view_recipe[ell]
  }

  cat_obj$links <- links
  cat_obj
}

# Scan a block through a cell/chunk provider. Only compact cell statistics and
# coordinate support facts survive the scan; responses/designs are never cached.
selection_block_scan <- function(provider, B_m, design_spec, N, spec_id) {
  d_s <- design_spec$d_s
  has_int <- design_spec$has_intercept
  n_fit <- 0
  pi <- r <- mass <- numeric(B_m)
  pos <- neg <- rep.int(FALSE, d_s)
  sep_cell <- NULL
  for (b in seq_len(B_m)) {
    nb <- wy <- wb <- 0
    all_one <- all_zero <- TRUE
    a_k <- mass_U_k <- NA_real_
    provider(b, function(sl) {
      nr <- length(sl$y)
      if (!is.numeric(sl$y) || anyNA(sl$y) || any(!sl$y %in% c(0, 1)))
        cp04_abort("Only binary selection responses in {0, 1} are supported.", "unsupported",
                   list(model = spec_id, cell = b))
      if (!is.matrix(sl$D) || nrow(sl$D) != nr || ncol(sl$D) != d_s || length(sl$w) != nr)
        cp04_abort("Selection chunk dimensions disagree.", "input")
      cp04_finite(sl$D, "selection block design")
      cp04_finite(sl$w, "selection block fitting weights", TRUE)
      if (length(sl$a_k) != 1L || !is.finite(sl$a_k) || sl$a_k <= 0 ||
          length(sl$mass_U_k) != 1L || !is.finite(sl$mass_U_k) || sl$mass_U_k <= 0)
        cp04_abort("A relevant cohort has no repeat-observation mass.", "support",
                   list(model = spec_id, cell = b))
      a_k <<- sl$a_k
      mass_U_k <<- sl$mass_U_k
      nb <<- nb + nr
      wy <<- wy + sum(sl$w * sl$y)
      wb <<- wb + sum(sl$w)
      all_one <<- all_one && all(sl$y == 1)
      all_zero <<- all_zero && all(sl$y == 0)
      for (j in seq_len(d_s)) {
        signed <- ifelse(sl$y == 1, sl$D[, j], -sl$D[, j])
        pos[j] <<- pos[j] || any(signed > 0)
        neg[j] <<- neg[j] || any(signed < 0)
      }
    })
    if (nb < 1L) cp04_abort("Selection cell has no fitting rows.", "support", list(model = spec_id, cell = b))
    n_fit <- n_fit + nb
    pi[b] <- wy / (mass_U_k * N)
    r[b] <- pi[b] / a_k
    mass[b] <- wb / N
    if (has_int && is.null(sep_cell) && (all_one || all_zero))
      sep_cell <- list(coordinate = b, cell = b, sign = if (all_one) 1 else -1)
  }
  sep <- sep_cell
  if (is.null(sep)) {
    j <- which(xor(pos, neg))[1L]
    if (!is.na(j)) sep <- list(coordinate = (if (has_int) B_m else 0L) + j,
                              slope = j, sign = if (pos[j]) 1 else -1)
  }
  list(n_fit = n_fit, pi = pi, r = r, mass = mass, separation = sep)
}

# Map the unchanged legacy CP04 fit/factor in the supplied design coordinates
# to the backend's canonical (local intercept, slopes) coordinates.
selection_cp04_block <- function(fit, design_spec, m_b, z_bar) {
  has_int <- design_spec$has_intercept
  list(branch = "cp04_separate_pair", has_intercept = has_int,
       B_m = 1L, d_s = design_spec$d_s,
       alpha = if (has_int) fit$coefficients[design_spec$int_col] else numeric(0),
       beta = fit$coefficients[design_spec$slope_cols],
       coefficients = fit$coefficients, H_chol = fit$H_chol,
       coordinate_order = c(design_spec$int_col, design_spec$slope_cols),
       m_b = m_b, z_bar = z_bar, slope_fit = NULL,
       pi = fit$pi, r = fit$r, diagnostics = fit$diagnostics)
}

# The small explicit-table adapter and production owner provider share ONE IRLS
# implementation. A provider invokes consume(sl) in canonical cell/row order,
# with at most max_chunk_rows rows and no retained per-cell panels.
fit_selection_block_slices <- function(slices, design_spec, glm_control, N,
                                       max_chunk_rows = N,
                                       allow_cp04_separate = FALSE,
                                       spec_id = NULL) {
  B_m <- length(slices)
  if (isTRUE(allow_cp04_separate) && B_m == 1L && design_spec$d > 0L) {
    sl <- slices[[1L]]
    canonical <- if (design_spec$has_intercept) cbind(1, sl$D) else sl$D
    coord <- c(design_spec$int_col, design_spec$slope_cols)
    X <- canonical[, order(coord), drop = FALSE]
    colnames(X) <- NULL
    fit <- cp04_logit(X, sl$y, sl$w, glm_control, N, spec_id)
    fit$pi <- sum(sl$w * sl$y) / (sl$mass_U_k * N)
    fit$r <- fit$pi / sl$a_k
    pr <- ipw_logit_values(X, fit$coefficients)$prob
    q <- sl$w * pr * (1 - pr)
    zb <- if (design_spec$has_intercept && design_spec$d_s > 0L)
      matrix(colSums(q * sl$D) / sum(q), 1L) else matrix(numeric(0), 1L, 0L)
    return(selection_cp04_block(fit, design_spec, sum(q) / N, zb))
  }
  provider <- function(b, consume) {
    sl <- slices[[b]]
    for (start in seq.int(1L, length(sl$y), by = max_chunk_rows)) {
      rr <- seq.int(start, min(length(sl$y), start + max_chunk_rows - 1L))
      consume(list(D = sl$D[rr, , drop = FALSE], y = sl$y[rr], w = sl$w[rr],
                   a_k = sl$a_k, mass_U_k = sl$mass_U_k))
    }
  }
  fit_selection_block_provider(provider, B_m, design_spec, glm_control, N, spec_id)
}

fit_selection_block_provider <- function(provider, B_m, design_spec, glm_control, N,
                                         spec_id = NULL) {
  gc <- validate_ipw_glm_control(glm_control)
  if (B_m < 1L) cp04_abort("Selection block requires at least one active cell.", "input")
  has_int <- design_spec$has_intercept
  d_s <- design_spec$d_s
  p_m <- (if (has_int) B_m else 0L) + d_s
  scan <- selection_block_scan(provider, B_m, design_spec, N, spec_id)
  fail <- function(message, type, extra = list()) {
    cp04_abort(message, type, c(list(model = spec_id, n_fit = scan$n_fit, p = p_m), extra))
  }
  if (scan$n_fit < p_m) fail("Fitting rows smaller than logit dimension.", "rank")
  if (!has_int && d_s == 0L) {
    return(list(branch = "zero_col", has_intercept = FALSE, B_m = B_m, d_s = 0L,
                alpha = numeric(0), beta = numeric(0), coefficients = numeric(0),
                H_chol = matrix(numeric(0), 0L, 0L), m_b = scan$mass * 0.25,
                z_bar = matrix(numeric(0), B_m, 0L), slope_fit = NULL,
                pi = scan$pi, r = scan$r,
                diagnostics = list(converged = TRUE, boundary = FALSE, rank = 0L,
                                   iter = 0L, warnings = character(0))))
  }
  rank_tol <- max(1e-12, min(1e-7, gc$epsilon / 1000))
  eps_prob <- 10 * .Machine$double.eps
  sep <- scan$separation
  binom_dev <- stats::binomial(link = "logit")$dev.resids
  alpha <- if (has_int) numeric(B_m) else numeric(0)
  beta <- numeric(d_s)
  boundary <- FALSE
  warnings_vec <- character(0)

  values <- function(sl, b, final = FALSE) {
    eta <- (if (has_int) alpha[b] else 0) +
      (if (d_s > 0L) drop(sl$D %*% beta) else rep.int(0, length(sl$y)))
    if (final) {
      odds <- exp(eta)
      mu <- odds / (1 + odds)
      cp04_finite(eta, "logit predictors")
      cp04_finite(odds, "logit odds", TRUE)
      cp04_finite(mu, "logit probabilities")
      if (any(mu <= eps_prob | mu >= 1 - eps_prob | mu * (1 - mu) <= 0))
        fail("Logit probabilities reached numerical boundary.", "boundary", diag)
    } else {
      if (any(!is.finite(eta))) return(NULL)
      if (any(abs(eta) > 30)) {
        boundary <<- TRUE
        warnings_vec <<- unique(c(warnings_vec, "fitted probabilities numerically 0 or 1 occurred"))
      }
      odds <- exp(eta)
      mu <- ifelse(eta >= 0, 1 / (1 + exp(-eta)), odds / (1 + odds))
      if (any(mu <= eps_prob | mu >= 1 - eps_prob)) {
        boundary <<- TRUE
        warnings_vec <<- unique(c(warnings_vec, "fitted probabilities numerically 0 or 1 occurred"))
        mu <- pmin(pmax(mu, eps_prob), 1 - eps_prob)
      }
    }
    list(eta = eta, mu = mu)
  }

  # Left-fold the same QR summaries as before, releasing each row chunk at once.
  # Cell partials are small factors; no lists of row chunks/eta/mu are retained.
  build_block_wls <- function(final = FALSE) {
    partials <- if (has_int) vector("list", B_m) else NULL
    pooled <- NULL
    m_b <- numeric(B_m)
    ok <- TRUE
    for (b in seq_len(B_m)) {
      cell_summary <- NULL
      provider(b, function(sl) {
        if (!ok) return(invisible(NULL))
        lv <- values(sl, b, final)
        if (is.null(lv)) { ok <<- FALSE; return(invisible(NULL)) }
        mu_eta <- lv$mu * (1 - lv$mu)
        q <- sl$w * mu_eta
        z <- if (final) rep.int(0, length(sl$y)) else lv$eta + (sl$y - lv$mu) / mu_eta
        if (any(!is.finite(q) | q <= 0) || any(!is.finite(z))) {
          ok <<- FALSE; return(invisible(NULL))
        }
        m_b[b] <<- m_b[b] + sum(q) / N
        if (has_int) {
          one <- or_pooled_wls_block_summary(sl$D, z, q)
          cell_summary <<- if (is.null(cell_summary)) one else
            or_pooled_wls_merge_chunks(list(cell_summary, one))
        } else {
          one <- or_wls_summary(sl$D, z, q)
          pooled <<- if (is.null(pooled)) one else or_wls_merge(list(pooled, one))
        }
      })
      if (!ok) return(NULL)
      if (has_int) partials[[b]] <- or_pooled_wls_partial_block(cell_summary)
    }
    list(summary = if (has_int) or_pooled_wls_merge(partials) else pooled, m_b = m_b)
  }
  deviance <- function() {
    dev <- 0
    for (b in seq_len(B_m)) provider(b, function(sl) {
      lv <- values(sl, b)
      dev <<- if (is.null(lv)) NaN else dev + sum(binom_dev(sl$y, lv$mu, sl$w))
    })
    dev
  }
  dev_old <- deviance()
  converged <- FALSE
  iter_used <- 0L
  for (iter in seq_len(gc$maxit)) {
    iter_used <- iter
    summary_iter <- tryCatch(build_block_wls(), error = function(e) e)
    if (is.null(summary_iter) || inherits(summary_iter, "error")) {
      if (iter == 1L) fail("Logit has deficient model rank.", "rank", list(iter = iter))
      boundary <- TRUE
      break
    }
    wls_fit <- tryCatch(
      if (has_int) or_pooled_wls_fit(summary_iter$summary, rank_tol = rank_tol) else
        or_wls_fit(summary_iter$summary, rank_tol = rank_tol), error = function(e) e)
    if (inherits(wls_fit, "error")) {
      if (iter == 1L || is.null(sep))
        fail("Logit has deficient model rank.", "rank", list(iter = iter, solver_message = conditionMessage(wls_fit)))
      fail("Logit has no finite root: signed-coordinate separation witness.", "separation",
           list(iter = iter, separation_witness = sep))
    }
    alpha <- if (has_int) as.double(wls_fit$intercepts) else numeric(0)
    beta <- as.double(wls_fit$coefficients)
    if (any(!is.finite(c(alpha, beta)))) {
      if (!is.null(sep)) fail("Logit has no finite root: signed-coordinate separation witness.", "separation",
                             list(iter = iter, separation_witness = sep))
      fail("Nonfinite logit coefficients during IRLS.", "numeric", list(iter = iter))
    }
    dev <- deviance()
    if (!is.finite(dev)) { boundary <- TRUE; break }
    if (abs(dev - dev_old) / (0.1 + abs(dev)) < gc$epsilon) { converged <- TRUE; break }
    dev_old <- dev
  }
  if (!converged) warnings_vec <- unique(c(warnings_vec, "glm.fit: algorithm did not converge"))
  diag <- list(converged = converged, boundary = boundary, rank = p_m,
               iter = iter_used, warnings = warnings_vec)
  if (!is.null(sep)) fail("Logit has no finite root: signed-coordinate separation witness.", "separation",
                          c(diag, list(separation_witness = sep)))
  if (boundary) fail("Logit encountered boundary probabilities.", "boundary", diag)
  if (!converged) fail("Logit did not converge under supplied controls.", "convergence", diag)
  cp04_finite(c(alpha, beta), "logit coefficients")
  info <- tryCatch(build_block_wls(final = TRUE), error = function(e) {
    if (inherits(e, "cdid_cp04_boundary_error")) stop(e)
    NULL
  })
  if (is.null(info)) fail("Logit information is singular/nonfinite.", "rank", diag)
  info_fit <- tryCatch(if (has_int) or_pooled_wls_fit(info$summary, rank_tol = rank_tol) else
    or_wls_fit(info$summary, rank_tol = rank_tol), error = function(e) NULL)
  if (is.null(info_fit)) fail("Logit information is singular/nonfinite.", "rank", diag)
  sf <- if (has_int) info_fit$slope_fit else info_fit
  if (d_s > 0L && (any(!is.finite(sf$R)) || any(abs(diag(sf$R)) <= 0)))
    fail("Logit information is singular/nonfinite.", "rank", diag)
  m_b <- if (has_int) vapply(info_fit$blocks, function(bl) bl$weight_sum / N, numeric(1)) else info$m_b
  cp04_finite(m_b, "cell Fisher mass", TRUE)
  z_bar <- if (has_int && d_s > 0L) do.call(rbind, lapply(info_fit$blocks, function(bl) as.double(bl$mean_d))) else
    matrix(numeric(0), B_m, 0L)
  cp04_finite(z_bar, "cell Fisher-weighted covariate means")
  slope_fit <- if (d_s > 0L) structure(
    list(coefficients = numeric(d_s), rss = 0, R = sf$R, pivot = sf$pivot,
         n = sf$n, p = sf$p, rank = sf$rank, rank_tol = sf$rank_tol, cutoff = sf$cutoff),
    class = "cdid_or_wls_fit") else NULL
  list(branch = design_spec$branch, has_intercept = has_int, B_m = B_m, d_s = d_s,
       alpha = alpha, beta = beta, coefficients = c(alpha, beta), H_chol = NULL,
       m_b = m_b, z_bar = z_bar, slope_fit = slope_fit,
       pi = scan$pi, r = scan$r, diagnostics = diag)
}

# Production scans immutable cohort-order ranges and filters U inside each
# bounded chunk. Even row-index vectors are not replicated across cells.
# A D60 cell reads its recipe's selection design, weights W_g and repeat rate
# a_k^{g,U} on the declared roster U of its cohort.
selection_owner_block_provider <- function(reader, cids, cells, design_spec, ctrl, selection_model, reference = NULL) {
  st <- reader$state
  function(b, consume) {
    cell <- cells[[cids[b]]]
    k <- cell$cohort
    r <- cell$reference
    a_k <- if (is.null(reference)) st$a[k] else chained_reference_stat(reference, r, k, "a")
    mass_U_k <- if (is.null(reference)) st$mass_U[k] else chained_reference_stat(reference, r, k, "mass_U")
    for (start in seq.int(st$cohort_start[k], st$cohort_end[k], by = ctrl$max_chunk_rows)) {
      idx <- seq.int(start, min(st$cohort_end[k], start + ctrl$max_chunk_rows - 1L))
      rows <- st$cohort_order[idx]
      rows <- rows[st$U[rows]]
      if (!length(rows)) next
      D <- if (is.null(reference)) st$Z[rows, design_spec$slope_cols, drop = FALSE] else
        chained_reference_s(reference, r, rows)[, design_spec$slope_cols, drop = FALSE]
      colnames(D) <- NULL
      y <- if (identical(selection_model, "independent_period_logit"))
        cp04_date_rows(reader, rows, cell$date_spec)$R else cp04_pair_rows(reader, rows, cell$pair)$S
      w <- if (is.null(reference)) st$w[rows] else chained_reference_weight(reference, r, rows)
      consume(list(D = D, y = y, w = w, a_k = a_k, mass_U_k = mass_U_k))
    }
    invisible(NULL)
  }
}

fit_pooled_selection_block <- function(reader, block_spec, cells, design_spec,
                                       glm_control, ctrl,
                                       selection_model = "pair_logit",
                                       selection_pooling = "separate", reference = NULL) {
  st <- reader$state
  cids <- block_spec$cell_ids
  provider <- selection_owner_block_provider(reader, cids, cells, design_spec, ctrl, selection_model, reference)
  if (identical(selection_model, "pair_logit") && identical(selection_pooling, "separate") &&
      length(cids) == 1L && design_spec$d > 0L) {
    cell <- cells[[cids[1L]]]
    sp <- list(pair = cell$pair, cohort = cell$cohort)
    if (!is.null(cell$reference)) sp$reference <- cell$reference
    fit <- cp04_selection_fit(reader, sp, glm_control, ctrl, reference)
    m_b <- 0
    z_sum <- numeric(design_spec$d_s)
    provider(1L, function(sl) {
      eta <- (if (design_spec$has_intercept) fit$coefficients[design_spec$int_col] else 0) +
        (if (design_spec$d_s > 0L) drop(sl$D %*% fit$coefficients[design_spec$slope_cols]) else rep.int(0, length(sl$y)))
      odds <- exp(eta)
      pr <- odds / (1 + odds)
      q <- sl$w * pr * (1 - pr)
      m_b <<- m_b + sum(q) / st$N
      if (design_spec$has_intercept && design_spec$d_s > 0L) z_sum <<- z_sum + colSums(q * sl$D)
    })
    zb <- if (design_spec$has_intercept && design_spec$d_s > 0L)
      matrix(z_sum / (m_b * st$N), 1L) else matrix(numeric(0), 1L, 0L)
    res <- selection_cp04_block(fit, design_spec, m_b, zb)
  } else {
    res <- fit_selection_block_provider(provider, length(cids), design_spec, glm_control, st$N, block_spec)
  }
  res$block_id <- block_spec$block_id
  res$cell_ids <- cids
  res
}

validate_selection_models <- function(backend) {
  if (!inherits(backend, "cdid_selection_models") ||
      !ipw_owner_seal_ok(backend, "selection_models"))
    cp04_abort("Expected the immutable private selection-models backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) && exists("cp04_reader", backend$reader_ref, inherits = FALSE)) {
    get("cp04_reader", backend$reader_ref, inherits = FALSE)
  } else {
    NULL
  }
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

# Conservative bounds for the actual implementation: one row chunk, compact
# cell QR partials, and (only on the preserved CP04 route) one full cohort fit.
# There is no C-by-C fitting allocation.
selection_fit_resources <- function(st, block_specs, cells, design_spec, gc, ctrl,
                                    selection_model, selection_pooling) {
  d <- design_spec$d
  d_s <- design_spec$d_s
  B <- min(ctrl$max_chunk_rows, st$N)
  C_max <- max(vapply(block_specs, function(bs) bs$B_m, numeric(1)))
  n_scan <- sum(vapply(cells, function(cell) {
    k <- cell$cohort
    st$cohort_end[k] - st$cohort_start[k] + 1
  }, numeric(1)))
  legacy_rows <- if (identical(selection_model, "pair_logit") &&
                    identical(selection_pooling, "separate") && d > 0L)
    max(st$cohort_end - st$cohort_start + 1) else 0
  scratch <- 128 * (B * (d + 8) + (C_max + 1) * (d_s + 2)^2 + legacy_rows * (d + 8))
  log_b <- 1 + ceiling(log2(st$N + 1))
  # Initial support/deviance, two scans per IRLS iteration, final information;
  # every repeated response lookup and every chunk/cell QR merge is charged.
  work <- 64 * (2 * gc$maxit + 4) *
    (n_scan * ((d_s + 2)^2 + log_b) + length(cells) * (d_s + 2)^3)
  list(scratch = scratch, work = work, scan_rows = n_scan,
       max_chunk_rows = B, legacy_rows = legacy_rows)
}

# `reference` (D60): cells keyed (pair or date, cohort, recipe) under each
# recipe's law P^g and information C^r; pooling inside one recipe.
fit_selection_models <- function(reader, link_views,
                                 selection_model = "pair_logit",
                                 selection_pooling = "separate",
                                 control = reader$control,
                                 glm_control,
                                 has_intercept = NULL,
                                 retained_external_entries = 0,
                                 external_work_ops = 0,
                                 reference = NULL, view_cohorts = NULL, families = NULL) {
  validate_cp04_reader(reader)
  ctrl <- validate_ipw_owner_control(control)
  gc <- validate_ipw_glm_control(glm_control)
  st <- reader$state
  validate_selection_mode_args(selection_model, selection_pooling, st$universe, st$U, st$rho)

  if (!is.numeric(retained_external_entries) || length(retained_external_entries) != 1L ||
      !is.finite(retained_external_entries) || retained_external_entries < 0 ||
      !is.numeric(external_work_ops) || length(external_work_ops) != 1L ||
      !is.finite(external_work_ops) || external_work_ops < 0)
    cp04_abort("Invalid external pipeline reservation.", "input")

  L <- if (is.list(link_views)) length(link_views) else 0L
  if (L < 1L) cp04_abort("Supply actual private link views.", "input")

  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries + retained_external_entries
  if (!is.null(reference)) {
    validate_chained_reference_reader(reference)
    if (!is.list(link_views) || length(reference$view_recipe) != length(link_views))
      cp04_abort("The reference reader must index exactly the supplied link views.", "integrity")
    base <- base + reference$accounting$owned_entries
  }
  d_cols <- if (is.null(reference)) ncol(st$Z) else length(reference$designs[[reference$role_slot["s"]]]$columns)
  # Preflight design validation before its bounded scans/copies.
  design_B <- min(ctrl$max_chunk_rows, st$N)
  cp04_budget(ctrl, base, 32 * design_B * (d_cols + 4),
              c(32 * st$N * (d_cols + 4), external_work_ops), 0,
              "selection design scan")
  design_spec <- if (is.null(reference)) extract_selection_design_spec(st$Z, st$U, has_intercept = has_intercept,
                                                chunk_rows = ctrl$max_chunk_rows) else
    chained_reference_selection_spec(reference, st$U, ctrl$max_chunk_rows)
  header_bound <- 2048 + 1024 * L + sum(vapply(link_views, function(v) {
    if (is.list(v) && is.list(v$cohorts)) 128 * length(v$cohorts$cohorts) else 0
  }, numeric(1)))
  cp04_budget(ctrl, c(base, header_bound), header_bound, 64 * L * L,
              header_bound, "selection model catalog metadata")

  cat_obj <- build_selection_block_catalog(
    reader = reader,
    link_views = link_views,
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    control = ctrl,
    reference = reference,
    view_cohorts = view_cohorts,
    families = families
  )

  cells <- cat_obj$cells
  block_specs <- cat_obj$block_specs
  pair_cohorts <- cat_obj$pair_cohorts
  links <- cat_obj$links
  C <- length(cells)
  M <- length(block_specs)
  d_s <- design_spec$d_s
  d_full <- design_spec$d

  metadata_payload <- cp04_payload(list(cells, block_specs, pair_cohorts, links, design_spec))
  metadata_entries <- unname(metadata_payload["entries"])

  # Retained bound for M blocks: each block stores alpha (B_m), beta (d_s), m_b (B_m),
  # z_bar (B_m x d_s), slope_fit R (d_s x d_s) or H_chol ((d_s+1)^2), plus header slots.
  model_bound <- sum(vapply(block_specs, function(bs) {
    B_m <- bs$B_m
    B_m * (d_s + 8) + (d_s + 2)^2 + 4096
  }, numeric(1))) + C * 64 + 4096

  resources <- selection_fit_resources(st, block_specs, cells, design_spec, gc, ctrl,
                                       selection_model, selection_pooling)
  scratch <- resources$scratch + header_bound
  sel_work <- resources$work + 32 * st$N * (d_full + 4)
  work <- sel_work + 64 * (L * L + unname(metadata_payload["char_bytes"]))
  retained_bound <- base + metadata_entries + model_bound + 4096

  cp04_budget(ctrl, retained_bound, scratch, c(work, external_work_ops),
              metadata_entries + model_bound + 4096, "selection models fit")

  blocks <- lapply(block_specs, function(bs) {
    fit_pooled_selection_block(
      reader = reader,
      block_spec = bs,
      cells = cells,
      design_spec = design_spec,
      glm_control = gc,
      ctrl = ctrl,
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      reference = reference
    )
  })

  # Stacked parameter summary in canonical order: all local cell intercepts (if has_intercept)
  # followed by all block slope vectors in canonical block order.
  alpha_all <- if (design_spec$has_intercept) {
    unlist(lapply(blocks, `[[`, "alpha"), use.names = FALSE)
  } else {
    numeric(0)
  }
  beta_all <- unlist(lapply(blocks, `[[`, "beta"), use.names = FALSE)
  theta_all <- c(alpha_all, beta_all)

  law_str <- paste0(
    "CP04:", st$universe, ":weighted-", selection_model, ":", selection_pooling, ":a_U/(rho*p_U)"
  )
  obj <- structure(
    list(
      selection_model = selection_model,
      selection_pooling = selection_pooling,
      design_spec = design_spec,
      cells = cells,
      block_specs = block_specs,
      blocks = blocks,
      pair_cohorts = pair_cohorts,
      links = links,
      alpha = alpha_all,
      beta = beta_all,
      coefficients = theta_all,
      glm_control = gc,
      law = law_str,
      universe = st$universe,
      transport_justification = st$justification,
      reader_ref = reader$ref,
      control = ctrl,
      accounting = list(
        owned_entries = 0,
        character_bytes = 0,
        borrowed_entries = base,
        retained_bound = retained_bound,
        fit_scratch_bound = scratch,
        work_ops_bound = work,
        n_active_cells = C,
        n_selection_blocks = M,
        n_pair_cohorts = length(pair_cohorts),
        n_parameters = length(theta_all),
        full_F_masks_per_model = 0L,
        full_F_by_cell_IF_retained = FALSE,
        global_dummy_matrix_constructed = FALSE,
        inverse_matrices_retained = FALSE,
        peak_rss_measured = FALSE
      )
    ),
    class = "cdid_selection_models"
  )
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$law <- paste0(obj$law, ":recipe-indexed-W_g-C^r")
  }
  obj["ref"] <- list(NULL)
  payload <- cp04_payload(obj)
  actual <- unname(payload["entries"])
  if (actual > metadata_entries + model_bound + 4096)
    cp04_abort("Retained selection model bound underestimated.", "bound")
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  ipw_owner_seal(obj, "selection_models")
}

resolve_pair_cohort_id <- function(backend, pair_cohort_id) {
  if (is.numeric(pair_cohort_id) && !is.object(pair_cohort_id) && is.null(dim(pair_cohort_id)) &&
      length(pair_cohort_id) == 1L && !is.na(pair_cohort_id) && pair_cohort_id == floor(pair_cohort_id) &&
      pair_cohort_id >= 1L && pair_cohort_id <= length(backend$pair_cohorts)) {
    return(as.integer(pair_cohort_id))
  }
  if (is.list(pair_cohort_id) && !is.null(pair_cohort_id$cohort) && !is.null(pair_cohort_id$pair)) {
    k <- as.integer(pair_cohort_id$cohort)
    s_val <- as.double(pair_cohort_id$pair$s)
    t_val <- as.double(pair_cohort_id$pair$t)
    r <- pair_cohort_id$reference
    for (q in seq_along(backend$pair_cohorts)) {
      pc <- backend$pair_cohorts[[q]]
      if (pc$cohort == k && pc$pair$s == s_val && pc$pair$t == t_val &&
          identical(as.integer(pc$reference), as.integer(r))) return(q)
    }
  }
  cp04_abort("Unknown or invalid `pair_cohort_id`.", "input")
}

predict_selection_cell <- function(backend, rows, cell_id) {
  validate_selection_models(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows != floor(rows)) || any(rows < 1L | rows > length(st$a_index)))
    cp04_abort("Invalid row indices for predict_selection_cell.", "input")
  if (!is.numeric(cell_id) || length(cell_id) != 1L || is.na(cell_id) ||
      cell_id < 1L || cell_id > length(backend$cells))
    cp04_abort("Invalid `cell_id`.", "input")

  cell <- backend$cells[[as.integer(cell_id)]]
  block <- backend$blocks[[cell$block_id]]
  b <- cell$local_b
  ref <- chained_backend_reference(backend)

  if (identical(block$branch, "zero_col")) {
    return(rep.int(0.5, length(rows)))
  }
  if (identical(block$branch, "cp04_separate_pair")) {
    sp <- list(pair = cell$pair, cohort = cell$cohort)
    if (!is.null(cell$reference)) sp$reference <- cell$reference
    return(cp04_predict_selection(reader, rows, sp, block, ref))
  }
  d_s <- block$d_s
  Z_slope <- if (d_s > 0L) {
    if (is.null(ref)) st$Z[rows, backend$design_spec$slope_cols, drop = FALSE] else
      chained_reference_s(ref, cell$reference, rows)[, backend$design_spec$slope_cols, drop = FALSE]
  } else matrix(0, length(rows), 0L)
  eta <- (if (block$has_intercept) block$alpha[b] else 0) +
    (if (d_s > 0L) drop(Z_slope %*% block$beta) else rep.int(0, length(rows)))
  p <- ipw_clip_logit(eta)$prob
  cp04_finite(p, "replayed selection probabilities", TRUE)
  if (any(p >= 1)) cp04_abort("Replayed selection probability at boundary.", "boundary", cell)
  p
}

predict_period_pair_selection <- function(backend, rows, pair_cohort_id) {
  validate_selection_models(backend)
  if (!identical(backend$selection_model, "independent_period_logit"))
    cp04_abort("`predict_period_pair_selection` requires selection_model == \"independent_period_logit\".", "unsupported")
  q <- resolve_pair_cohort_id(backend, pair_cohort_id)
  pc <- backend$pair_cohorts[[q]]
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state

  ps <- predict_selection_cell(backend, rows, pc$s_cell_id)
  pt <- predict_selection_cell(backend, rows, pc$t_cell_id)
  p_pair <- ps * pt
  cp04_finite(p_pair, "replayed product pair selection probabilities", TRUE)
  if (any(p_pair >= 1)) cp04_abort("Replayed product pair selection probability at boundary.", "boundary", pc)
  ref <- chained_backend_reference(backend)
  a_pc <- if (is.null(ref)) st$a[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
  w_rows <- if (is.null(ref)) st$w[rows] else chained_reference_weight(ref, pc$reference, rows)
  # The final raw product is floored once; raw endpoints stay for scores.
  p_weight <- ipw_clip_floor(p_pair)
  w_sel <- w_rows * a_pc / (st$rho[rows] * p_weight)
  cp04_finite(w_sel, "replayed selection weights", TRUE)
  list(prob = p_pair, prob_s = ps, prob_t = pt, prob_weight = p_weight,
       active = ipw_clip_floor_active(p_pair), weight = w_sel, a = a_pc, pair_cohort_id = q)
}

predict_pair_selection <- function(backend, rows, pair_cohort_id) {
  validate_selection_models(backend)
  if (identical(backend$selection_model, "independent_period_logit")) {
    return(predict_period_pair_selection(backend, rows, pair_cohort_id))
  }
  q <- resolve_pair_cohort_id(backend, pair_cohort_id)
  pc <- backend$pair_cohorts[[q]]
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state

  p_pair <- predict_selection_cell(backend, rows, pc$pair_cell_id)
  ref <- chained_backend_reference(backend)
  a_pc <- if (is.null(ref)) st$a[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
  w_rows <- if (is.null(ref)) st$w[rows] else chained_reference_weight(ref, pc$reference, rows)
  p_weight <- ipw_clip_floor(p_pair)
  w_sel <- w_rows * a_pc / (st$rho[rows] * p_weight)
  cp04_finite(w_sel, "replayed selection weights", TRUE)
  list(prob = p_pair, prob_weight = p_weight, active = ipw_clip_floor_active(p_pair),
       weight = w_sel, a = a_pc, pair_cohort_id = q)
}

# Diagnostic/test entry points below read the static roster design and rates;
# a recipe-indexed (D60) backend is consumed only through the recipe-keyed
# predictions, log-weight derivatives and upstream rows of the consumers.
selection_static_only <- function(backend, what) {
  if (!is.null(backend$reference_ref))
    cp04_abort(paste0("`", what, "` reads the static selection design; recipe-indexed (D60) selection models are ",
      "replayed only through their consumers' recipe-keyed channels."), "unsupported")
  invisible(TRUE)
}

# Accumulates selection forcing moments (t_alpha, t_beta) and repeat-rate forcing
# moment t_a for one or more consumer pair-cohorts (e, k).
# Default sign = -1 corresponds to the inverse-probability derivative
# d/dtheta [1 / p_{e,k}^{sel}(X_i)] = -(1 / p_{e,k}^{sel}(X_i)) * d/dtheta log p_{e,k}^{sel}(X_i)
# when observer_weights already include 1 / p_{e,k}^{sel}(X_i).
add_selection_moments <- function(backend, pair_cohort_id, rows, observer_weights,
                                  acc = NULL, sign = -1) {
  selection_static_only(backend, "add_selection_moments")
  validate_selection_models(backend)
  q <- resolve_pair_cohort_id(backend, pair_cohort_id)
  pc <- backend$pair_cohorts[[q]]
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state

  if (!is.numeric(observer_weights) || is.object(observer_weights) || !is.null(dim(observer_weights)) ||
      length(observer_weights) != length(rows) || any(!is.finite(observer_weights)) ||
      !is.numeric(sign) || length(sign) != 1L || !is.finite(sign))
    cp04_abort("Invalid `observer_weights` or `sign` in add_selection_moments.", "input")

  selection_if_budget(backend, max(1L, length(rows)), q, details = FALSE)
  if (is.null(acc)) {
    acc <- list(
      t_a = numeric(length(st$cohorts)),
      blocks = lapply(backend$blocks, function(bl) {
        list(
          t_alpha = if (bl$has_intercept) numeric(bl$B_m) else numeric(0),
          t_beta = numeric(bl$d_s)
        )
      })
    )
  }

  pred <- predict_pair_selection(backend, rows, q)
  acc$t_a[pc$cohort] <- acc$t_a[pc$cohort] + sum(observer_weights / st$a[pc$cohort]) / st$N
  d_s <- backend$design_spec$d_s
  Z_slope <- if (d_s > 0L) st$Z[rows, backend$design_spec$slope_cols, drop = FALSE] else matrix(0, length(rows), 0L)

  if (identical(backend$selection_model, "independent_period_logit")) {
    # Endpoint s: cell c_s = (k, s)
    ms <- pc$s_block_id
    bs <- pc$s_local_b
    bl_s <- backend$blocks[[ms]]
    ws <- sign * observer_weights * ((1 - pred$prob_s) * pred$active)
    if (bl_s$has_intercept) {
      acc$blocks[[ms]]$t_alpha[bs] <- acc$blocks[[ms]]$t_alpha[bs] + sum(ws) / st$N
    }
    if (bl_s$d_s > 0L) {
      acc$blocks[[ms]]$t_beta <- acc$blocks[[ms]]$t_beta + colSums(ws * Z_slope) / st$N
    }
    # Endpoint t: cell c_t = (k, t)
    mt <- pc$t_block_id
    bt <- pc$t_local_b
    bl_t <- backend$blocks[[mt]]
    wt <- sign * observer_weights * ((1 - pred$prob_t) * pred$active)
    if (bl_t$has_intercept) {
      acc$blocks[[mt]]$t_alpha[bt] <- acc$blocks[[mt]]$t_alpha[bt] + sum(wt) / st$N
    }
    if (bl_t$d_s > 0L) {
      acc$blocks[[mt]]$t_beta <- acc$blocks[[mt]]$t_beta + colSums(wt * Z_slope) / st$N
    }
  } else {
    m <- pc$block_id
    b <- pc$local_b
    bl <- backend$blocks[[m]]
    wp <- sign * observer_weights * ((1 - pred$prob) * pred$active)
    if (bl$has_intercept) {
      acc$blocks[[m]]$t_alpha[b] <- acc$blocks[[m]]$t_alpha[b] + sum(wp) / st$N
    }
    if (bl$d_s > 0L) {
      acc$blocks[[m]]$t_beta <- acc$blocks[[m]]$t_beta + colSums(wp * Z_slope) / st$N
    }
  }

  acc
}

# Pure block IF evaluator: computes xi_alpha (B x B_m) and xi_beta (B x d_s)
# for a chunk of original units from their B x B_m raw cell scores s_{ib}^{sel}.
# All repeated cell scores for original unit i within block m are summed BEFORE
# Fisher-information inversion.
evaluate_selection_block_if <- function(block_fit, Z_slope_rows, cell_scores, N) {
  B <- nrow(cell_scores)
  B_m <- block_fit$B_m
  d_s <- block_fit$d_s

  if (identical(block_fit$branch, "zero_col")) {
    return(list(
      xi_alpha = matrix(0, B, 0L),
      xi_beta = matrix(0, B, 0L),
      xi_theta = matrix(0, B, 0L)
    ))
  }
  if (identical(block_fit$branch, "cp04_separate_pair")) {
    s1 <- cell_scores[, 1L]
    canonical <- if (block_fit$has_intercept) cbind(1, Z_slope_rows) else Z_slope_rows
    coord <- block_fit$coordinate_order
    score_mat <- s1 * canonical[, order(coord), drop = FALSE]
    xi_raw <- t(cp04_solve(block_fit$H_chol, t(score_mat)))
    xi_theta <- xi_raw[, coord, drop = FALSE]
    return(list(
      xi_alpha = if (block_fit$has_intercept) xi_theta[, 1L, drop = FALSE] else matrix(0, B, 0L),
      xi_beta = if (block_fit$has_intercept) xi_theta[, -1L, drop = FALSE] else xi_theta,
      xi_theta = xi_theta
    ))
  }

  if (block_fit$has_intercept) {
    if (d_s > 0L) {
      sum_s <- rowSums(cell_scores)
      s_tilde_beta <- (sum_s * Z_slope_rows) - (cell_scores %*% block_fit$z_bar)
      xi_beta <- t(or_wls_gram_solve(block_fit$slope_fit, t(s_tilde_beta), normalizer = N))
      xi_alpha <- sweep(cell_scores, 2L, block_fit$m_b, "/") - (xi_beta %*% t(block_fit$z_bar))
    } else {
      xi_beta <- matrix(0, B, 0L)
      xi_alpha <- sweep(cell_scores, 2L, block_fit$m_b, "/")
    }
  } else {
    sum_s <- rowSums(cell_scores)
    s_beta0 <- sum_s * Z_slope_rows
    xi_beta <- t(or_wls_gram_solve(block_fit$slope_fit, t(s_beta0), normalizer = N))
    xi_alpha <- matrix(0, B, 0L)
  }

  list(
    xi_alpha = xi_alpha,
    xi_beta = xi_beta,
    xi_theta = cbind(xi_alpha, xi_beta)
  )
}

# Evaluates the original-unit influence contribution on `rows` for an accumulated
# selection moment object `acc` (from add_selection_moments) using the profiled
# contraction identity (DATE_SELECTION_SPEC.md Section 4.1 Eq. 7):
#   sum_b (t_b^alpha / m_b) s_{ib} + (t_m^beta - sum_b t_b^alpha z_bar_b)' xi_{i,beta_m}.
evaluate_pair_selection_loading <- function(backend, rows, acc, include_repeat_rate = TRUE) {
  selection_static_only(backend, "evaluate_pair_selection_loading")
  validate_selection_models(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows != floor(rows)) || any(rows < 1L | rows > st$N) ||
      length(rows) > backend$control$max_chunk_rows)
    cp04_abort("Invalid row chunk for evaluate_pair_selection_loading.", "input")
  rows <- as.integer(rows)
  B <- length(rows)
  selection_if_budget(backend, B, seq_along(backend$pair_cohorts), details = FALSE)
  ci <- st$cohort_index[rows]
  w <- st$w[rows]
  A <- as.double(st$a_index[rows] > 0L)
  U_rows <- st$U[rows]

  loading <- numeric(B)
  if (isTRUE(include_repeat_rate) && !is.null(acc$t_a)) {
    for (k in which(acc$t_a != 0)) {
      rr <- which(ci == k & U_rows)
      if (length(rr)) {
        a_if_k <- w[rr] * (A[rr] - st$a[k]) / st$mass_U[k]
        loading[rr] <- loading[rr] + acc$t_a[k] * a_if_k
      }
    }
  }

  d_s <- backend$design_spec$d_s
  Z_slope <- if (d_s > 0L) st$Z[rows, backend$design_spec$slope_cols, drop = FALSE] else matrix(0, B, 0L)

  for (m in seq_along(backend$blocks)) {
    bm <- acc$blocks[[m]]
    if (is.null(bm)) next
    if (all(bm$t_alpha == 0) && all(bm$t_beta == 0)) next
    bl <- backend$blocks[[m]]
    if (identical(bl$branch, "zero_col")) next

    cids <- bl$cell_ids
    B_m <- bl$B_m
    cell_scores <- matrix(0, B, B_m)
    for (b in seq_len(B_m)) {
      cell <- backend$cells[[cids[b]]]
      rr <- which(ci == cell$cohort & U_rows)
      if (!length(rr)) next
      pr <- predict_selection_cell(backend, rows[rr], cids[b])
      Y <- if (identical(backend$selection_model, "independent_period_logit")) {
        cp04_date_rows(reader, rows[rr], cell$date_spec)$R
      } else {
        cp04_pair_rows(reader, rows[rr], cell$pair)$S
      }
      cell_scores[rr, b] <- w[rr] * (Y - pr)
    }

    if (identical(bl$branch, "cp04_separate_pair")) {
      bif <- evaluate_selection_block_if(bl, Z_slope, cell_scores, st$N)
      t_vec <- c(bm$t_alpha, bm$t_beta)
      loading <- loading + drop(bif$xi_theta %*% t_vec)
    } else if (bl$has_intercept) {
      # Profiled contraction (Eq. 7): intercept term + centered slope term
      loading <- loading + drop(cell_scores %*% (bm$t_alpha / bl$m_b))
      if (bl$d_s > 0L) {
        t_beta_eff <- bm$t_beta - drop(crossprod(bl$z_bar, bm$t_alpha))
        if (any(t_beta_eff != 0)) {
          bif <- evaluate_selection_block_if(bl, Z_slope, cell_scores, st$N)
          loading <- loading + drop(bif$xi_beta %*% t_beta_eff)
        }
      }
    } else if (bl$d_s > 0L) {
      bif <- evaluate_selection_block_if(bl, Z_slope, cell_scores, st$N)
      loading <- loading + drop(bif$xi_beta %*% bm$t_beta)
    }
  }

  cp04_finite(loading, "selection loading IF")
  loading
}

# Build probe evaluation targets: either from caller-supplied `probes`
# (list(cohort = ..., Z = ..., rho = ...)) or default cohort-mean probes for requested pair_cohorts.
build_selection_probe_targets <- function(backend, pair_cohort_ids, probes = NULL) {
  selection_static_only(backend, "build_selection_probe_targets")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  d_s <- backend$design_spec$d_s
  slope_cols <- backend$design_spec$slope_cols
  selection_if_budget(backend, 1L, pair_cohort_ids, probes = probes, details = FALSE)

  if (is.null(probes)) {
    # One evaluation target per requested pair_cohort_id at its cohort U-weighted mean covariate
    n_targets <- length(pair_cohort_ids)
    targets <- vector("list", n_targets)
    means <- vector("list", length(st$cohorts))
    for (j in seq_len(n_targets)) {
      q <- pair_cohort_ids[j]
      pc <- backend$pair_cohorts[[q]]
      k <- pc$cohort
      if (is.null(means[[k]])) {
        sw <- sr <- 0
        sz <- numeric(d_s)
        for (start in seq.int(st$cohort_start[k], st$cohort_end[k], by = backend$control$max_chunk_rows)) {
          idx <- seq.int(start, min(st$cohort_end[k], start + backend$control$max_chunk_rows - 1L))
          rr <- st$cohort_order[idx]
          rr <- rr[st$U[rr]]
          wk <- st$w[rr]
          sw <- sw + sum(wk)
          sr <- sr + sum(wk * st$rho[rr])
          if (d_s > 0L) sz <- sz + colSums(wk * st$Z[rr, slope_cols, drop = FALSE])
        }
        means[[k]] <- list(z = sz / sw, rho = sr / sw)
      }
      z_eval <- means[[k]]$z
      rho_eval <- means[[k]]$rho
      targets[[j]] <- list(
        pair_cohort_id = q,
        cohort = k,
        z_slope = as.double(z_eval),
        rho = as.double(rho_eval),
        label = pc$label
      )
    }
    return(targets)
  }

  # Explicit probe list: list(cohort = ..., Z = ..., rho = ..., pair_cohort_ids = ...)
  if (!is.list(probes) || is.null(probes$cohort) || is.null(probes$Z))
    cp04_abort("`probes` must be a list with `cohort` and `Z`.", "input")
  Z_p <- probes$Z
  if (!is.matrix(Z_p) || !is.numeric(Z_p))
    cp04_abort("`probes$Z` must be a numeric matrix.", "input")
  n_p <- nrow(Z_p)
  if (length(probes$cohort) != n_p)
    cp04_abort("`probes$cohort` length must match nrow(probes$Z).", "input")
  coh_idx <- if (all(probes$cohort %in% st$cohorts)) {
    match(as.double(probes$cohort), st$cohorts)
  } else {
    as.integer(probes$cohort)
  }
  if (anyNA(coh_idx) || any(coh_idx < 1L | coh_idx > length(st$cohorts)))
    cp04_abort("Invalid cohort in `probes$cohort`.", "input")
  rho_p <- if (is.null(probes$rho)) rep.int(1, n_p) else as.double(probes$rho)
  if (length(rho_p) != n_p || any(!is.finite(rho_p) | rho_p <= 0))
    cp04_abort("Invalid `probes$rho`.", "input")

  Z_p_slope <- if (ncol(Z_p) == d_s) {
    Z_p
  } else if (ncol(Z_p) == backend$design_spec$d) {
    Z_p[, slope_cols, drop = FALSE]
  } else {
    cp04_abort("`probes$Z` column count must match selection design.", "input")
  }

  targets <- list()
  for (p_idx in seq_len(n_p)) {
    k <- coh_idx[p_idx]
    z_row <- as.double(Z_p_slope[p_idx, ])
    # All requested pair_cohort_ids belonging to cohort k
    matching_q <- pair_cohort_ids[vapply(backend$pair_cohorts[pair_cohort_ids], function(pc) pc$cohort == k, logical(1))]
    for (q in matching_q) {
      pc <- backend$pair_cohorts[[q]]
      targets[[length(targets) + 1L]] <- list(
        pair_cohort_id = q,
        cohort = k,
        z_slope = z_row,
        rho = rho_p[p_idx],
        label = paste0("p", p_idx, "_", pc$label)
      )
    }
  }
  if (!length(targets))
    cp04_abort("Supplied `probes` did not match any requested pair-cohorts.", "input")
  targets
}

predict_selection_probes <- function(backend, probes = NULL, pair_cohort_ids = NULL) {
  validate_selection_models(backend)
  if (is.null(pair_cohort_ids)) pair_cohort_ids <- seq_along(backend$pair_cohorts)
  targets <- build_selection_probe_targets(backend, pair_cohort_ids, probes)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  n_t <- length(targets)
  pair_prob <- numeric(n_t)
  multiplier <- numeric(n_t)

  eval_cell <- function(cell_id, z_slope) {
    cell <- backend$cells[[cell_id]]
    bl <- backend$blocks[[cell$block_id]]
    b <- cell$local_b
    if (identical(bl$branch, "zero_col")) return(0.5)
    eta <- (if (bl$has_intercept) bl$alpha[b] else 0) +
      (if (bl$d_s > 0L) sum(z_slope * bl$beta) else 0)
    p <- ipw_clip_logit(eta)$prob
    cp04_finite(p, "probe selection probability", TRUE)
    if (p >= 1) cp04_abort("Probe selection probability at boundary.", "boundary", cell)
    p
  }

  for (j in seq_len(n_t)) {
    tg <- targets[[j]]
    pc <- backend$pair_cohorts[[tg$pair_cohort_id]]
    if (identical(backend$selection_model, "independent_period_logit")) {
      ps <- eval_cell(pc$s_cell_id, tg$z_slope)
      pt <- eval_cell(pc$t_cell_id, tg$z_slope)
      pp <- ps * pt
      tg$prob_s <- ps
      tg$prob_t <- pt
    } else {
      pp <- eval_cell(pc$pair_cell_id, tg$z_slope)
    }
    cp04_finite(pp, "probe pair probability", TRUE)
    if (pp >= 1) cp04_abort("Probe pair probability at boundary.", "boundary", pc)
    pair_prob[j] <- pp
    multiplier[j] <- st$a[pc$cohort] / (tg$rho * ipw_clip_floor(pp))
    tg$pair_prob <- pp
    tg$multiplier <- multiplier[j]
    targets[[j]] <- tg
  }
  names(pair_prob) <- names(multiplier) <- vapply(targets, `[[`, character(1), "label")
  list(pair_probabilities = pair_prob, multipliers = multiplier, targets = targets)
}

selection_if_budget <- function(backend, B, pair_cohort_ids, probes = NULL, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  ctrl <- backend$control
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)

  active_blocks <- if (details) seq_along(backend$blocks) else
    unique(unlist(lapply(backend$pair_cohorts[pair_cohort_ids], `[[`, "block_ids")))
  active_cells <- unique(unlist(lapply(backend$block_specs[active_blocks], `[[`, "cell_ids")))
  n_theta <- length(backend$coefficients)
  n_a <- length(st$cohorts)
  n_probes <- if (is.null(probes)) length(pair_cohort_ids) else nrow(probes$Z) * length(pair_cohort_ids)
  widths <- (if (details) 3 else 1) * n_theta + (if (details) 2 else 1) * n_a +
    (if (details) 5 else 3) * n_probes + length(active_cells) * (if (details) 3 else 1) + 32

  max_bytes <- frame$id_byte_metadata$max
  metadata <- 8192 + 1024 * (length(backend$blocks) + length(backend$cells) + n_a) +
    n_probes * (backend$design_spec$d_s + 1024) + B * max_bytes
  output <- B * (widths + 16) + metadata
  scratch <- 32 * B * (widths + backend$design_spec$d + 8) +
    64 * (backend$design_spec$d_s + 2)^2 + metadata
  retained <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  log_b <- 1 + ceiling(log2(st$N + 1))
  per_row_work <- 32 * length(active_cells) * log_b +
    32 * sum(vapply(backend$blocks[active_blocks], function(bl) bl$B_m * (bl$d_s + 2)^2, numeric(1))) +
    32 * (n_theta + n_a + 2 * n_probes + 8)
  target_scan <- if (is.null(probes)) 32 * st$N * (backend$design_spec$d + 4) else 0
  work <- B * per_row_work + target_scan + 64 * metadata
  cp04_budget(ctrl, retained, scratch, work, output, "selection IF chunk")
  list(output = output, scratch = scratch, work = work,
       active_blocks = active_blocks, active_cells = active_cells)
}

evaluate_selection_if_chunk <- function(backend, rows,
                                        pair_cohort_ids = NULL,
                                        probes = NULL,
                                        scale = "frame",
                                        details = FALSE) {
  selection_static_only(backend, "evaluate_selection_if_chunk")
  validate_selection_models(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  ctrl <- backend$control

  if (is.null(pair_cohort_ids)) pair_cohort_ids <- seq_along(backend$pair_cohorts)
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
      anyNA(rows) || any(rows != floor(rows)) || any(rows < 1L | rows > st$N) ||
      anyDuplicated(rows) || length(rows) > ctrl$max_chunk_rows ||
      !is.numeric(pair_cohort_ids) || is.object(pair_cohort_ids) || !is.null(dim(pair_cohort_ids)) ||
      !length(pair_cohort_ids) || anyNA(pair_cohort_ids) || any(pair_cohort_ids != floor(pair_cohort_ids)) ||
      any(pair_cohort_ids < 1L | pair_cohort_ids > length(backend$pair_cohorts)) ||
      anyDuplicated(pair_cohort_ids) ||
      !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    cp04_abort("Invalid explicit IF row/pair_cohort batch, scale or details flag.", "input")

  rows <- as.integer(rows)
  pair_cohort_ids <- as.integer(pair_cohort_ids)
  B <- length(rows)
  pre <- selection_if_budget(backend, B, pair_cohort_ids, probes = probes, details = details)

  ci <- st$cohort_index[rows]
  w <- st$w[rows]
  A <- as.double(st$a_index[rows] > 0L)
  U_rows <- st$U[rows]

  # 1. Empirical repeat-rate IFs xi_{i, a_k^U}
  K <- length(st$cohorts)
  a_if <- matrix(0, B, K)
  for (k in unique(ci)) {
    rr <- which(ci == k & U_rows)
    if (length(rr)) {
      a_if[rr, k] <- w[rr] * (A[rr] - st$a[k]) / st$mass_U[k]
    }
  }

  # 2. Block-by-block cell score aggregation across original units and Fisher solve
  d_s <- backend$design_spec$d_s
  Z_slope <- if (d_s > 0L) st$Z[rows, backend$design_spec$slope_cols, drop = FALSE] else matrix(0, B, 0L)
  C <- length(backend$cells)
  M <- length(backend$blocks)

  xi_alpha_by_cell <- if (backend$design_spec$has_intercept) matrix(0, B, C) else matrix(0, B, 0L)
  xi_beta_by_block <- vector("list", M)
  cell_scores_all <- if (details) matrix(0, B, C) else NULL
  cell_probs_all <- if (details) matrix(0, B, C) else NULL
  cell_resp_all <- if (details) matrix(0, B, C) else NULL

  for (m in seq_len(M)) {
    bl <- backend$blocks[[m]]
    cids <- bl$cell_ids
    B_m <- bl$B_m
    if (m %in% pre$active_blocks || details) {
      cell_scores <- matrix(0, B, B_m)
      for (b in seq_len(B_m)) {
        cid <- cids[b]
        cell <- backend$cells[[cid]]
        rr <- which(ci == cell$cohort & U_rows)
        if (!length(rr)) next
        pr <- predict_selection_cell(backend, rows[rr], cid)
        Y <- if (identical(backend$selection_model, "independent_period_logit")) {
          cp04_date_rows(reader, rows[rr], cell$date_spec)$R
        } else {
          cp04_pair_rows(reader, rows[rr], cell$pair)$S
        }
        sc <- w[rr] * (Y - pr)
        cell_scores[rr, b] <- sc
        if (details) {
          cell_scores_all[rr, cid] <- sc
          cell_probs_all[rr, cid] <- pr
          cell_resp_all[rr, cid] <- Y
        }
      }
      bif <- evaluate_selection_block_if(bl, Z_slope, cell_scores, st$N)
      if (bl$has_intercept) {
        xi_alpha_by_cell[, cids] <- bif$xi_alpha
      }
      xi_beta_by_block[[m]] <- bif$xi_beta
    } else {
      xi_beta_by_block[[m]] <- matrix(0, B, bl$d_s)
    }
  }

  xi_beta_all <- if (d_s > 0L) do.call(cbind, xi_beta_by_block) else matrix(0, B, 0L)
  theta_if <- cbind(xi_alpha_by_cell, xi_beta_all)

  # 3. Evaluate pair-probability IFs and multiplier IFs on probe targets
  probe_pred <- predict_selection_probes(backend, probes = probes, pair_cohort_ids = pair_cohort_ids)
  targets <- probe_pred$targets
  n_t <- length(targets)
  pair_prob_if <- matrix(0, B, n_t)
  multiplier_if <- matrix(0, B, n_t)
  colnames(pair_prob_if) <- colnames(multiplier_if) <- names(probe_pred$pair_probabilities)

  cell_eta_if <- function(cid, z_slope) {
    cell <- backend$cells[[cid]]
    bl <- backend$blocks[[cell$block_id]]
    if (identical(bl$branch, "zero_col")) return(numeric(B))
    ans <- if (bl$has_intercept) xi_alpha_by_cell[, cid] else numeric(B)
    if (bl$d_s > 0L) {
      ans <- ans + drop(xi_beta_by_block[[cell$block_id]] %*% z_slope)
    }
    ans
  }

  for (j in seq_len(n_t)) {
    tg <- targets[[j]]
    pc <- backend$pair_cohorts[[tg$pair_cohort_id]]
    pp <- tg$pair_prob
    if (identical(backend$selection_model, "independent_period_logit")) {
      eta_s_if <- cell_eta_if(pc$s_cell_id, tg$z_slope)
      eta_t_if <- cell_eta_if(pc$t_cell_id, tg$z_slope)
      log_p_if <- (1 - tg$prob_s) * eta_s_if + (1 - tg$prob_t) * eta_t_if
    } else {
      eta_if <- cell_eta_if(pc$pair_cell_id, tg$z_slope)
      log_p_if <- (1 - pp) * eta_if
    }
    # Raw pair-probability IF; the multiplier a / (rho max(p, eps)) loses its
    # probability path when the floor binds and keeps its a_k path.
    xp <- pp * log_p_if
    xa <- a_if[, pc$cohort]
    qc <- ipw_clip_floor(pp)
    xm <- xa / (tg$rho * qc) - (st$a[pc$cohort] / (tg$rho * qc^2)) * (ipw_clip_floor_active(pp) * xp)
    pair_prob_if[, j] <- xp
    multiplier_if[, j] <- xm
  }

  cp04_finite(c(theta_if, a_if, pair_prob_if, multiplier_if), "selection chunk IF")

  scale_fac <- if (identical(scale, "core")) st$lambda else 1
  out <- list(
    rows = rows,
    id = st$id[rows],
    pair_cohort_ids = pair_cohort_ids,
    theta_if = scale_fac * theta_if,
    a_if = scale_fac * a_if,
    pair_prob_if = scale_fac * pair_prob_if,
    multiplier_if = scale_fac * multiplier_if,
    phi = scale_fac * multiplier_if,
    scale = scale,
    normalization_count = if (identical(scale, "core")) st$n_A else st$N,
    diagnostics = list(
      n_active_cells_evaluated = length(pre$active_cells),
      n_active_blocks_evaluated = length(pre$active_blocks),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  if (details) {
    out$details <- list(
      theta_if_frame = theta_if,
      a_if_frame = a_if,
      pair_prob_if_frame = pair_prob_if,
      multiplier_if_frame = multiplier_if,
      xi_alpha = xi_alpha_by_cell,
      xi_beta_by_block = xi_beta_by_block,
      cell_scores = cell_scores_all,
      cell_probabilities = cell_probs_all,
      cell_responses = cell_resp_all
    )
  }
  actual_output <- unname(cp04_payload(out)["entries"])
  if (actual_output > pre$output) cp04_abort("Selection IF output bound underestimated.", "bound")
  out
}

replay_selection_if <- function(backend, chunk_rows, pair_batch_size = NULL,
                                callback, probes = NULL,
                                scale = "frame", details = FALSE) {
  selection_static_only(backend, "replay_selection_if")
  validate_selection_models(backend)
  ctrl <- backend$control
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  n_pc <- length(backend$pair_cohorts)
  if (is.null(pair_batch_size)) pair_batch_size <- n_pc
  PB <- ipw_owner_int_bound(pair_batch_size, "pair_batch_size")
  if (B > ctrl$max_chunk_rows || !is.function(callback) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core") ||
      !identical(details, isTRUE(details)))
    cp04_abort("Invalid selection replay controls/callback.", "input")

  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N

  # Preflight cumulative work BEFORE any callback runs
  total_work <- 0
  for (ps in seq.int(1L, n_pc, by = PB)) {
    pp_ids <- seq.int(ps, min(n_pc, ps + PB - 1L))
    pre <- selection_if_budget(backend, min(B, N), pp_ids, probes = probes, details = details)
    total_work <- ipw_owner_safe_sum(
      c(total_work, ceiling(N / B) * pre$work),
      ctrl$max_work_ops,
      "total selection replay work"
    )
  }
  ipw_owner_safe_sum(total_work, ctrl$max_work_ops, "total selection replay work")

  batches <- 0L
  for (ps in seq.int(1L, n_pc, by = PB)) {
    pp_ids <- seq.int(ps, min(n_pc, ps + PB - 1L))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_selection_if_chunk(
        backend, rr, pair_cohort_ids = pp_ids,
        probes = probes, scale = scale, details = details
      )
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(
    n_frame_rows = N,
    normalization_count = if (identical(scale, "core")) reader$state$n_A else N,
    n_pair_cohorts = n_pc,
    batches = batches,
    scale = scale,
    work_ops_bound = total_work,
    all_F_identities = TRUE
  ))
}

stream_selection_covariance <- function(backend, chunk_rows,
                                        pair_cohort_ids = NULL,
                                        probes = NULL,
                                        scale = "frame",
                                        include_theta = TRUE) {
  selection_static_only(backend, "stream_selection_covariance")
  validate_selection_models(backend)
  ctrl <- backend$control
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  if (B > ctrl$max_chunk_rows || !identical(include_theta, isTRUE(include_theta)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    cp04_abort("Invalid covariance streaming controls.", "input")

  if (is.null(pair_cohort_ids)) pair_cohort_ids <- seq_along(backend$pair_cohorts)
  pair_cohort_ids <- as.integer(pair_cohort_ids)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  N <- st$N

  pre <- selection_if_budget(backend, min(B, N), pair_cohort_ids, probes = probes, details = FALSE)
  # Explicit probes can match at most every requested pair-cohort. Reserve this
  # upper bound before building targets or scanning default cohort means.
  target_bound <- if (is.null(probes)) length(pair_cohort_ids) else nrow(probes$Z) * length(pair_cohort_ids)
  n_theta <- if (include_theta) length(backend$coefficients) + length(st$cohorts) else 0L
  W_bound <- n_theta + 2L * target_bound
  cov_output_bound <- W_bound^2 + 2 * target_bound^2 + n_theta^2 + 4096
  cov_work_bound <- pre$work + ceiling(N / B) * (pre$work + 32 * min(B, N) * W_bound^2)
  retained <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  cp04_budget(ctrl, retained, pre$scratch + 8 * W_bound^2 + 8192,
              cov_work_bound, cov_output_bound, "selection covariance stream")
  probe_pred <- predict_selection_probes(backend, probes = probes, pair_cohort_ids = pair_cohort_ids)
  n_t <- length(probe_pred$targets)
  W <- n_theta + 2L * n_t

  cov_full <- matrix(0, W, W)
  cov_pair_prob <- matrix(0, n_t, n_t)
  cov_multiplier <- matrix(0, n_t, n_t)
  cov_theta <- if (include_theta && length(backend$coefficients) > 0L) {
    matrix(0, length(backend$coefficients), length(backend$coefficients))
  } else {
    matrix(0, 0L, 0L)
  }

  for (start in seq.int(1L, N, by = B)) {
    rr <- seq.int(start, min(N, start + B - 1L))
    ch <- evaluate_selection_if_chunk(
      backend, rr, pair_cohort_ids = pair_cohort_ids,
      probes = probes, scale = scale, details = FALSE
    )
    stack_ch <- if (include_theta) {
      cbind(ch$theta_if, ch$a_if, ch$pair_prob_if, ch$multiplier_if)
    } else {
      cbind(ch$pair_prob_if, ch$multiplier_if)
    }
    cov_full <- cov_full + crossprod(stack_ch)
    cov_pair_prob <- cov_pair_prob + crossprod(ch$pair_prob_if)
    cov_multiplier <- cov_multiplier + crossprod(ch$multiplier_if)
    if (include_theta && ncol(ch$theta_if) > 0L) {
      cov_theta <- cov_theta + crossprod(ch$theta_if)
    }
  }

  norm_n <- if (identical(scale, "core")) st$n_A else N
  denom <- as.double(norm_n)^2
  out <- list(
    covariance = cov_full / denom,
    pair_prob_covariance = cov_pair_prob / denom,
    multiplier_covariance = cov_multiplier / denom,
    theta_covariance = cov_theta / denom,
    scale = scale,
    normalization_count = norm_n
  )
  actual_out <- unname(cp04_payload(out)["entries"])
  if (actual_out > cov_output_bound) cp04_abort("Selection covariance output bound underestimated.", "bound")
  out
}

# Pure-helper table driver for standalone scientific verification without fake owners.
# Executes the exact same partition_selection_cells, fit_selection_block_slices,
# and evaluate_selection_block_if pipeline on explicit numeric arrays in chunks.
evaluate_selection_table_if <- function(Z, cohort_index, w, A,
                                        U = rep.int(TRUE, nrow(Z)),
                                        rho = rep.int(1, nrow(Z)),
                                        universe = "F",
                                        cohorts = sort(unique(cohort_index)),
                                        all_times = NULL,
                                        date_responses = NULL,
                                        pair_responses = NULL,
                                        pair_specs,
                                        selection_model = "pair_logit",
                                        selection_pooling = "separate",
                                        glm_control,
                                        has_intercept = NULL,
                                        probes = NULL,
                                        chunk_rows = nrow(Z),
                                        scale = "frame",
                                        details = FALSE) {
  validate_selection_mode_args(selection_model, selection_pooling, universe, U, rho)
  design_spec <- extract_selection_design_spec(Z, U, has_intercept = has_intercept)
  gc <- validate_ipw_glm_control(glm_control)
  N <- nrow(Z)
  if (length(cohort_index) != N || length(w) != N || length(A) != N || length(U) != N || length(rho) != N)
    cp04_abort("Table input lengths must match nrow(Z).", "input")
  cp04_finite(w, "table weights", TRUE)
  cp04_finite(rho[U], "table rho on U", TRUE)

  K <- length(cohorts)
  ci <- match(cohort_index, cohorts)
  if (anyNA(ci)) cp04_abort("Unknown cohort in cohort_index.", "input")
  if (any(A > 0 & !U)) cp04_abort("U must contain all A structurally.", "support")

  mass_F <- mass_U <- repeat_mass <- numeric(K)
  for (k in seq_len(K)) {
    ik <- which(ci == k)
    ik_u <- ik[U[ik]]
    ik_a <- ik[A[ik] > 0]
    mass_F[k] <- sum(w[ik]) / N
    mass_U[k] <- sum(w[ik_u]) / N
    repeat_mass[k] <- sum(w[ik_a]) / N
  }
  a_vec <- ifelse(mass_U > 0, repeat_mass / mass_U, 0)

  cat_obj <- partition_selection_cells(
    pair_specs = pair_specs,
    cohorts = cohorts,
    all_times = all_times,
    selection_model = selection_model,
    selection_pooling = selection_pooling
  )
  cells <- cat_obj$cells
  block_specs <- cat_obj$block_specs
  pair_cohorts <- cat_obj$pair_cohorts
  C <- length(cells)
  M <- length(block_specs)
  d_s <- design_spec$d_s
  slope_cols <- design_spec$slope_cols

  get_cell_response <- function(cid, rr) {
    cell <- cells[[cid]]
    if (identical(selection_model, "independent_period_logit")) {
      if (is.matrix(date_responses)) {
        return(as.double(date_responses[rr, cell$calendar_index]))
      } else if (is.list(date_responses)) {
        return(as.double(date_responses[[cell$calendar_index]][rr]))
      }
      cp04_abort("Supply `date_responses` for independent_period_logit.", "input")
    } else {
      key <- paste0(cell$pair$s_calendar_index, ":", cell$pair$t_calendar_index)
      if (is.list(pair_responses) && !is.null(pair_responses[[key]])) {
        return(as.double(pair_responses[[key]][rr]))
      } else if (is.matrix(pair_responses)) {
        col_nm <- paste0("S", cell$pair$s_calendar_index, cell$pair$t_calendar_index)
        if (!is.null(colnames(pair_responses)) && col_nm %in% colnames(pair_responses)) {
          return(as.double(pair_responses[rr, col_nm]))
        }
      }
      cp04_abort("Supply matching `pair_responses` for pair_logit.", "input")
    }
  }

  blocks <- vector("list", M)
  for (m in seq_len(M)) {
    bs <- block_specs[[m]]
    cids <- bs$cell_ids
    B_m <- bs$B_m
    slices <- vector("list", B_m)
    for (b in seq_len(B_m)) {
      cid <- cids[b]
      cell <- cells[[cid]]
      k <- cell$cohort
      rows_k <- which(ci == k & U)
      D_b <- if (d_s > 0L) Z[rows_k, slope_cols, drop = FALSE] else matrix(numeric(0), length(rows_k), 0L)
      colnames(D_b) <- NULL
      slices[[b]] <- list(
        D = D_b,
        y = get_cell_response(cid, rows_k),
        w = w[rows_k],
        a_k = a_vec[k],
        mass_U_k = mass_U[k]
      )
    }
    allow_sep <- identical(selection_model, "pair_logit") &&
      identical(selection_pooling, "separate") &&
      B_m == 1L && design_spec$d >= 1L
    blocks[[m]] <- fit_selection_block_slices(
      slices = slices,
      design_spec = design_spec,
      glm_control = gc,
      N = N,
      max_chunk_rows = chunk_rows,
      allow_cp04_separate = allow_sep,
      spec_id = bs
    )
    blocks[[m]]$block_id <- m
    blocks[[m]]$cell_ids <- cids
  }

  alpha_all <- if (design_spec$has_intercept) unlist(lapply(blocks, `[[`, "alpha"), use.names = FALSE) else numeric(0)
  beta_all <- unlist(lapply(blocks, `[[`, "beta"), use.names = FALSE)
  theta_all <- c(alpha_all, beta_all)

  # Evaluate chunked IF across 1..N
  a_if <- matrix(0, N, K)
  xi_alpha_by_cell <- if (design_spec$has_intercept) matrix(0, N, C) else matrix(0, N, 0L)
  xi_beta_by_block <- lapply(blocks, function(bl) matrix(0, N, bl$d_s))
  cell_scores_all <- matrix(0, N, C)

  eval_cell_rows <- function(cid, rr) {
    cell <- cells[[cid]]
    bl <- blocks[[cell$block_id]]
    b <- cell$local_b
    if (identical(bl$branch, "zero_col")) return(rep.int(0.5, length(rr)))
    D_r <- if (bl$d_s > 0L) Z[rr, slope_cols, drop = FALSE] else matrix(0, length(rr), 0L)
    eta <- (if (bl$has_intercept) bl$alpha[b] else 0) +
      (if (bl$d_s > 0L) drop(D_r %*% bl$beta) else rep.int(0, length(rr)))
    ipw_clip_logit(eta)$prob
  }

  for (start in seq.int(1L, N, by = chunk_rows)) {
    rr <- seq.int(start, min(N, start + chunk_rows - 1L))
    B_ch <- length(rr)
    ci_ch <- ci[rr]
    w_ch <- w[rr]
    A_ch <- as.double(A[rr] > 0)
    U_ch <- U[rr]
    Z_ch_slope <- if (d_s > 0L) Z[rr, slope_cols, drop = FALSE] else matrix(0, B_ch, 0L)

    for (k in unique(ci_ch)) {
      idx_k <- which(ci_ch == k & U_ch)
      if (length(idx_k)) {
        a_if[rr[idx_k], k] <- w_ch[idx_k] * (A_ch[idx_k] - a_vec[k]) / mass_U[k]
      }
    }

    for (m in seq_len(M)) {
      bl <- blocks[[m]]
      cids <- bl$cell_ids
      B_m <- bl$B_m
      cs_ch <- matrix(0, B_ch, B_m)
      for (b in seq_len(B_m)) {
        cid <- cids[b]
        cell <- cells[[cid]]
        idx_k <- which(ci_ch == cell$cohort & U_ch)
        if (!length(idx_k)) next
        pr <- eval_cell_rows(cid, rr[idx_k])
        Y <- get_cell_response(cid, rr[idx_k])
        sc <- w_ch[idx_k] * (Y - pr)
        cs_ch[idx_k, b] <- sc
        cell_scores_all[rr[idx_k], cid] <- sc
      }
      bif <- evaluate_selection_block_if(bl, Z_ch_slope, cs_ch, N)
      if (bl$has_intercept) xi_alpha_by_cell[rr, cids] <- bif$xi_alpha
      xi_beta_by_block[[m]][rr, ] <- bif$xi_beta
    }
  }

  xi_beta_all <- if (d_s > 0L) do.call(cbind, xi_beta_by_block) else matrix(0, N, 0L)
  theta_if <- cbind(xi_alpha_by_cell, xi_beta_all)

  # Probe predictions and IFs
  pair_cohort_ids <- seq_along(pair_cohorts)
  targets <- list()
  if (is.null(probes)) {
    for (q in pair_cohort_ids) {
      pc <- pair_cohorts[[q]]
      k <- pc$cohort
      rk <- which(ci == k & U)
      z_eval <- if (d_s > 0L) colSums(w[rk] * Z[rk, slope_cols, drop = FALSE]) / sum(w[rk]) else numeric(0)
      rho_eval <- sum(w[rk] * rho[rk]) / sum(w[rk])
      targets[[length(targets) + 1L]] <- list(
        pair_cohort_id = q, cohort = k, z_slope = as.double(z_eval), rho = rho_eval, label = pc$label
      )
    }
  } else {
    Z_p <- probes$Z
    coh_p <- if (all(probes$cohort %in% cohorts)) match(probes$cohort, cohorts) else as.integer(probes$cohort)
    rho_p <- if (is.null(probes$rho)) rep.int(1, nrow(Z_p)) else as.double(probes$rho)
    Z_p_slope <- if (ncol(Z_p) == d_s) Z_p else Z_p[, slope_cols, drop = FALSE]
    for (p_idx in seq_len(nrow(Z_p))) {
      k <- coh_p[p_idx]
      z_row <- as.double(Z_p_slope[p_idx, ])
      matching_q <- pair_cohort_ids[vapply(pair_cohorts, function(pc) pc$cohort == k, logical(1))]
      for (q in matching_q) {
        pc <- pair_cohorts[[q]]
        targets[[length(targets) + 1L]] <- list(
          pair_cohort_id = q, cohort = k, z_slope = z_row, rho = rho_p[p_idx],
          label = paste0("p", p_idx, "_", pc$label)
        )
      }
    }
  }

  eval_cell_probe <- function(cid, z_slope) {
    cell <- cells[[cid]]
    bl <- blocks[[cell$block_id]]
    b <- cell$local_b
    if (identical(bl$branch, "zero_col")) return(0.5)
    eta <- (if (bl$has_intercept) bl$alpha[b] else 0) + (if (bl$d_s > 0L) sum(z_slope * bl$beta) else 0)
    ipw_clip_logit(eta)$prob
  }
  cell_eta_if <- function(cid, z_slope) {
    cell <- cells[[cid]]
    bl <- blocks[[cell$block_id]]
    if (identical(bl$branch, "zero_col")) return(numeric(N))
    ans <- if (bl$has_intercept) xi_alpha_by_cell[, cid] else numeric(N)
    if (bl$d_s > 0L) ans <- ans + drop(xi_beta_by_block[[cell$block_id]] %*% z_slope)
    ans
  }

  n_t <- length(targets)
  pair_prob <- numeric(n_t)
  multiplier <- numeric(n_t)
  pair_prob_if <- matrix(0, N, n_t)
  multiplier_if <- matrix(0, N, n_t)
  for (j in seq_len(n_t)) {
    tg <- targets[[j]]
    pc <- pair_cohorts[[tg$pair_cohort_id]]
    if (identical(selection_model, "independent_period_logit")) {
      ps <- eval_cell_probe(pc$s_cell_id, tg$z_slope)
      pt <- eval_cell_probe(pc$t_cell_id, tg$z_slope)
      pp <- ps * pt
      log_p_if <- (1 - ps) * cell_eta_if(pc$s_cell_id, tg$z_slope) +
        (1 - pt) * cell_eta_if(pc$t_cell_id, tg$z_slope)
    } else {
      pp <- eval_cell_probe(pc$pair_cell_id, tg$z_slope)
      log_p_if <- (1 - pp) * cell_eta_if(pc$pair_cell_id, tg$z_slope)
    }
    xp <- pp * log_p_if
    xa <- a_if[, pc$cohort]
    qc <- ipw_clip_floor(pp)
    xm <- xa / (tg$rho * qc) - (a_vec[pc$cohort] / (tg$rho * qc^2)) * (ipw_clip_floor_active(pp) * xp)
    pair_prob[j] <- pp
    multiplier[j] <- a_vec[pc$cohort] / (tg$rho * qc)
    pair_prob_if[, j] <- xp
    multiplier_if[, j] <- xm
  }

  n_A <- sum(A > 0)
  lambda <- n_A / N
  scale_fac <- if (identical(scale, "core")) lambda else 1

  list(
    selection_model = selection_model,
    selection_pooling = selection_pooling,
    design_spec = design_spec,
    cells = cells,
    block_specs = block_specs,
    blocks = blocks,
    pair_cohorts = pair_cohorts,
    alpha = alpha_all,
    beta = beta_all,
    coefficients = theta_all,
    a = a_vec,
    pair_probabilities = pair_prob,
    multipliers = multiplier,
    theta_if = scale_fac * theta_if,
    a_if = scale_fac * a_if,
    pair_prob_if = scale_fac * pair_prob_if,
    multiplier_if = scale_fac * multiplier_if,
    analytic_stack = scale_fac * cbind(theta_if, a_if, pair_prob_if, multiplier_if),
    cell_scores = if (details) cell_scores_all else NULL,
    xi_alpha = if (details) xi_alpha_by_cell else NULL,
    xi_beta_by_block = if (details) xi_beta_by_block else NULL,
    scale = scale,
    normalization_count = if (identical(scale, "core")) n_A else N
  )
}

accumulate_selection_table_moments <- function(table_fit, Z, pair_cohort_id, rows,
                                               observer_weights, acc = NULL, sign = -1) {
  N <- nrow(Z)
  pc <- table_fit$pair_cohorts[[pair_cohort_id]]
  d_s <- table_fit$design_spec$d_s
  slope_cols <- table_fit$design_spec$slope_cols
  Z_slope <- if (d_s > 0L) Z[rows, slope_cols, drop = FALSE] else matrix(0, length(rows), 0L)

  if (is.null(acc)) {
    acc <- list(
      t_a = numeric(length(table_fit$a)),
      blocks = lapply(table_fit$blocks, function(bl) {
        list(
          t_alpha = if (bl$has_intercept) numeric(bl$B_m) else numeric(0),
          t_beta = numeric(bl$d_s)
        )
      })
    )
  }

  eval_c <- function(cid) {
    cell <- table_fit$cells[[cid]]
    bl <- table_fit$blocks[[cell$block_id]]
    b <- cell$local_b
    if (identical(bl$branch, "zero_col")) return(rep.int(0.5, length(rows)))
    eta <- (if (bl$has_intercept) bl$alpha[b] else 0) +
      (if (bl$d_s > 0L) drop(Z_slope %*% bl$beta) else rep.int(0, length(rows)))
    ipw_clip_logit(eta)$prob
  }

  acc$t_a[pc$cohort] <- acc$t_a[pc$cohort] + sum(observer_weights / table_fit$a[pc$cohort]) / N
  if (identical(table_fit$selection_model, "independent_period_logit")) {
    ps <- eval_c(pc$s_cell_id)
    pt <- eval_c(pc$t_cell_id)
    act <- ipw_clip_floor_active(ps * pt)
    ms <- pc$s_block_id; bs <- pc$s_local_b; bl_s <- table_fit$blocks[[ms]]
    ws <- sign * observer_weights * ((1 - ps) * act)
    if (bl_s$has_intercept) acc$blocks[[ms]]$t_alpha[bs] <- acc$blocks[[ms]]$t_alpha[bs] + sum(ws) / N
    if (bl_s$d_s > 0L) acc$blocks[[ms]]$t_beta <- acc$blocks[[ms]]$t_beta + colSums(ws * Z_slope) / N

    mt <- pc$t_block_id; bt <- pc$t_local_b; bl_t <- table_fit$blocks[[mt]]
    wt <- sign * observer_weights * ((1 - pt) * act)
    if (bl_t$has_intercept) acc$blocks[[mt]]$t_alpha[bt] <- acc$blocks[[mt]]$t_alpha[bt] + sum(wt) / N
    if (bl_t$d_s > 0L) acc$blocks[[mt]]$t_beta <- acc$blocks[[mt]]$t_beta + colSums(wt * Z_slope) / N
  } else {
    pp <- eval_c(pc$pair_cell_id)
    m <- pc$block_id; b <- pc$local_b; bl <- table_fit$blocks[[m]]
    wp <- sign * observer_weights * ((1 - pp) * ipw_clip_floor_active(pp))
    if (bl$has_intercept) acc$blocks[[m]]$t_alpha[b] <- acc$blocks[[m]]$t_alpha[b] + sum(wp) / N
    if (bl$d_s > 0L) acc$blocks[[m]]$t_beta <- acc$blocks[[m]]$t_beta + colSums(wp * Z_slope) / N
  }
  acc
}

evaluate_selection_table_loading <- function(table_fit, Z, acc, include_repeat_rate = TRUE) {
  N <- nrow(Z)
  loading <- numeric(N)
  if (isTRUE(include_repeat_rate) && !is.null(acc$t_a)) {
    loading <- loading + drop(table_fit$a_if %*% acc$t_a)
  }
  d_s <- table_fit$design_spec$d_s
  slope_cols <- table_fit$design_spec$slope_cols
  Z_slope <- if (d_s > 0L) Z[, slope_cols, drop = FALSE] else matrix(0, N, 0L)

  for (m in seq_along(table_fit$blocks)) {
    bm <- acc$blocks[[m]]
    bl <- table_fit$blocks[[m]]
    if (identical(bl$branch, "zero_col")) next
    cids <- bl$cell_ids
    cs <- table_fit$cell_scores[, cids, drop = FALSE]
    if (identical(bl$branch, "cp04_separate_pair")) {
      bif <- evaluate_selection_block_if(bl, Z_slope, cs, N)
      loading <- loading + drop(bif$xi_theta %*% c(bm$t_alpha, bm$t_beta))
    } else if (bl$has_intercept) {
      loading <- loading + drop(cs %*% (bm$t_alpha / bl$m_b))
      if (bl$d_s > 0L) {
        t_beta_eff <- bm$t_beta - drop(crossprod(bl$z_bar, bm$t_alpha))
        if (any(t_beta_eff != 0)) {
          bif <- evaluate_selection_block_if(bl, Z_slope, cs, N)
          loading <- loading + drop(bif$xi_beta %*% t_beta_eff)
        }
      }
    } else if (bl$d_s > 0L) {
      bif <- evaluate_selection_block_if(bl, Z_slope, cs, N)
      loading <- loading + drop(bif$xi_beta %*% bm$t_beta)
    }
  }
  loading
}

