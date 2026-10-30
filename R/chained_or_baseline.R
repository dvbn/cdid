# Private selection=none control-only shared-slope OR (D47). One borrowed common
# frame, explicit outcome recipe, actual pairs, fixed weights w/r and Z_a column
# roles. Each (pair, cohort) support is scanned once; ONE joint blockwise fit has
# a shared one-period slope and one intercept per unique (pair, pool) block.
# No reader selection predictors, fitting roster, transport or logits are used.

or_baseline_abort <- function(message, class = "input") {
  stop(structure(list(message = message, call = NULL),
    class = c(paste0("cdid_or_baseline_", class, "_error"),
      "cdid_ipw_fit_error", "cdid_error", "error", "condition")))
}

or_baseline_rank_tol <- function(rank_tol) {
  if (missing(rank_tol) || !is.numeric(rank_tol) || is.object(rank_tol) ||
      !is.null(dim(rank_tol)) || length(rank_tol) != 1L ||
      !is.finite(rank_tol) || rank_tol <= 0 || rank_tol >= 1)
    or_baseline_abort("Supply explicit finite private rank_tol in (0, 1).")
  rank_tol
}

# Bounded cohort chunks of ONE actual (pair, cohort) support. Z_a rows are read
# only on its observed pair rows (S == 1 implies A); nothing is read off A.
# A D60 support reads its recipe's frozen reference design instead of Z_a.
# `complete` (CP08 stage 3) reads the same pair rows by complete-date position.
or_baseline_visit <- function(reader, spec, ctrl, visit, reference = NULL, complete = NULL) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- if (is.null(complete)) NULL else or_baseline_complete_panel(reader, complete)
  first <- st$cohort_start[spec$cohort]; last <- st$cohort_end[spec$cohort]
  for (start in seq.int(first, last, by = ctrl$max_chunk_rows)) {
    rows <- st$cohort_order[seq.int(start, min(last, start + ctrl$max_chunk_rows - 1))]
    hit <- if (is.null(panel)) cp04_pair_rows(reader, rows, spec$pair, TRUE) else
      or_baseline_complete_pair_rows(reader, panel, rows, spec$pair, TRUE)
    obs <- which(hit$S == 1)
    visit(rows, hit$S, obs, hit$dy[obs], if (is.null(reference)) frame$Z_a[st$a_index[rows[obs]], , drop = FALSE] else
      chained_reference_design(reference, spec$reference, rows[obs]))
  }
  invisible(NULL)
}

# One scan: the rate among A and the raw-w design summary. Raw-w factors are
# rescaled by 1/sqrt(r) after the rate is known, without a second scan.
# A D60 support reads its context weight W_g on A rows only.
or_baseline_support <- function(reader, spec, recipe, ctrl, reference = NULL, complete = NULL) {
  st <- reader$state
  acc <- or_pooled_support_start(recipe, spec$transitions)
  repeat_mass <- 0; n_repeat <- 0
  or_baseline_visit(reader, spec, ctrl, function(rows, S, obs, dy, Z) {
    A <- st$a_index[rows] > 0L
    w <- if (is.null(reference)) st$w[rows] else {
      x <- rep(NA_real_, length(rows))
      x[A] <- chained_reference_weight(reference, spec$reference, rows[A])
      x
    }
    repeat_mass <<- repeat_mass + sum(w[A]); n_repeat <<- n_repeat + sum(A)
    acc <<- or_pooled_support_add(acc, recipe, Z, dy, w[obs], "baseline")
  }, reference, complete)
  or_baseline_support_rate(acc, repeat_mass, n_repeat, st$N)
}

# CP08 stage 3: complete-support execution of the eligible baseline OR/DR
# route, selected before any support scan from actual facts of the sealed
# owners: no dated reference recipe; F = A (every frame row a pair-eligible
# A-unit); constant unit weights; the intercept-only design (one Z_a column,
# every entry exactly 1) read by an intercept-only outcome recipe and, for DR,
# as the whole treatment design (finite fixed-X columns are also supported);
# finite target cohorts, with never-treated controls for OR or approved risk
# sets containing the never-treated pool for DR; and every required date complete (its
# sealed posting slice is exactly the A-units 1..n_A). NULL keeps the general
# route. No fit or outcome read occurs here; cheapest facts first. The new
# fact scans reserve their work and bounded scratch before reading arrays.
or_baseline_complete_facts <- function(reader, recipe, records, specs, ctrl, reference = NULL,
                                       treatment = NULL, retained_external_entries = 0,
                                       external_work_ops = 0, metadata_entries = 0) {
  st <- reader$state
  if (!is.null(reference) || isTRUE(st$reference_carrier) || !length(records) || !length(specs)) return(NULL)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  po <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  Z <- frame$Z_a; z_names <- colnames(Z)
  if (!inherits(po, "cdid_ipw_panel_owner") || !identical(po$n_a_units, frame$n_a_units) ||
      st$n_A != po$n_a_units || st$N != st$n_A || !is.matrix(Z) || !is.double(Z) || ncol(Z) < 1L ||
      is.null(z_names) || anyNA(z_names) || !all(nzchar(z_names)) || anyDuplicated(z_names) ||
      !inherits(recipe, "cdid_or_pooled_recipe") || !isTRUE(recipe$intercept) ||
      !is.numeric(recipe$p) || length(recipe$p) != 1L || !is.finite(recipe$p) ||
      recipe$p < 0L || recipe$p != floor(recipe$p) || !all(recipe$columns %in% z_names))
    return(NULL)
  or_cols <- match(recipe$columns, z_names)
  icpt_col <- or_cols[recipe$intercept_index]
  trt_cols <- integer()
  if (!is.null(treatment)) {
    if (!inherits(treatment, "cdid_dr_treatment_recipe") || !is.numeric(treatment$p) ||
        length(treatment$p) != 1L || !is.finite(treatment$p) || treatment$p < 1L ||
        treatment$p != floor(treatment$p) || !is.integer(treatment$index) ||
        length(treatment$index) != treatment$p || anyNA(treatment$index) ||
        any(treatment$index < 1L | treatment$index > ncol(Z)) || anyDuplicated(treatment$index) ||
        !identical(treatment$columns, z_names[treatment$index]))
      return(NULL)
    trt_cols <- treatment$index
  }
  used_cols <- sort(unique(as.integer(c(or_cols, trt_cols))))
  n_z <- length(used_cols)
  never <- which(st$cohorts == Inf)
  n_coh <- length(st$cohorts)
  if (length(never) != 1L || !all(vapply(records, function(x) {
    is.integer(x$target) && length(x$target) == 1L && !is.na(x$target) &&
      x$target >= 1L && x$target <= n_coh && is.finite(st$cohorts[x$target]) &&
      is.integer(x$controls) && length(x$controls) >= 1L && !anyNA(x$controls) &&
      all(x$controls >= 1L & x$controls <= n_coh) && !anyDuplicated(x$controls) &&
      !is.unsorted(x$controls) && !(x$target %in% x$controls) &&
      (if (is.null(treatment)) identical(x$controls, never) else never %in% x$controls)
  }, logical(1))))
    return(NULL)
  target <- unique(vapply(records, function(x) as.integer(x$target), integer(1)))
  control <- sort(unique(as.integer(unlist(lapply(records, `[[`, "controls")))))
  w <- st$w
  if (!is.double(w)) return(NULL)
  step <- min(as.double(po$n_a_units), as.double(ctrl$max_chunk_rows))
  dates <- sort(unique(as.integer(unlist(lapply(specs, function(s)
    c(s$pair$s_calendar_index, s$pair$t_calendar_index))))))
  facts <- length(dates) * po$n_a_units + (1 + n_z) * st$N + 8 * length(records)
  scratch <- 8 * (1 + n_z) * step + 8 * length(dates) + 1024
  cp04_budget(ctrl, c(reader$accounting$owned_entries, reader$accounting$borrowed_owner_entries,
    retained_external_entries, metadata_entries), scratch, c(external_work_ops, 64 * facts),
    8 * length(dates) + 1024, "complete-support facts")
  other_cols <- setdiff(used_cols, icpt_col)
  start <- 1
  while (start <= st$N) {
    end <- min(st$N, start + step - 1)
    rows <- seq.int(start, end)
    wi <- w[rows]; ai <- st$a_index[rows]
    if (anyNA(wi) || !all(wi == w[1L]) || anyNA(ai) || !all(ai >= 1L & ai <= po$n_a_units)) return(NULL)
    zi <- Z[ai, icpt_col]
    if (anyNA(zi) || !all(zi == 1)) return(NULL)
    if (length(other_cols)) {
      zx <- Z[ai, other_cols, drop = FALSE]
      if (anyNA(zx) || !all(is.finite(zx))) return(NULL)
    }
    start <- end + 1
  }
  for (d in dates) if (!isTRUE(ipw_owner_posting_date_complete(po$postings, po$n_a_units, d, step))) return(NULL)
  post <- po$postings
  out <- list(route = "complete_support", p = as.integer(recipe$p), n_z = n_z, panel_ref = reader$panel_ref, n_a = po$n_a_units, step = step,
    target = target, control = control, dates = dates, date_start = post$date_start[dates],
    date_end = post$date_end[dates], date_count = post$date_count[dates])
  if (n_z > 1L) out$fact_reads <- facts
  out
}

# The facts bind the sealed postings the reader reaches now; any other binding
# refuses at this invocation/replay boundary and is never re-planned elsewhere.
or_baseline_complete_panel <- function(reader, complete, abort = or_baseline_abort) {
  po <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  post <- if (is.list(po)) po$postings else NULL
  if (!is.list(complete) || !identical(complete$route, "complete_support") ||
      !identical(complete$panel_ref, reader$panel_ref) || !inherits(po, "cdid_ipw_panel_owner") ||
      !identical(po$n_a_units, complete$n_a) || !is.integer(post$post_a_index) ||
      !identical(post$date_start[complete$dates], complete$date_start) ||
      !identical(post$date_end[complete$dates], complete$date_end) ||
      !identical(post$date_count[complete$dates], complete$date_count))
    abort("Complete-support facts do not describe the reader's sealed postings.", "integrity")
  po
}

# `cp04_pair_rows` on complete dates of an F = A reader, by position: A-unit a
# occupies slot date_start + a - 1 of a complete date. Both slots are checked to
# hold a before any endpoint is read; the raw right-minus-left difference, row
# order and returned fields are those of `cp04_pair_rows`.
or_baseline_complete_pair_rows <- function(reader, panel, rows, pair, read_outcomes = FALSE,
                                           abort = or_baseline_abort) {
  post <- panel$postings
  a <- reader$state$a_index[rows]
  ps <- post$date_start[pair$s_calendar_index] - 1L + a
  pt <- post$date_start[pair$t_calendar_index] - 1L + a
  if (any(a < 1L) || !isTRUE(all(post$post_a_index[ps] == a)) || !isTRUE(all(post$post_a_index[pt] == a)))
    abort("A complete-date posting slot does not hold its A-unit.", "integrity")
  out <- list(S = rep(1, length(rows)))
  if (read_outcomes) {
    yy <- ipw_owner_read_y_positions(panel, post$post_obs_index[ps], post$post_obs_index[pt])
    out$dy <- yy$y_right_obs - yy$y_left_obs
    cp04_finite(out$dy, "observed pair differences")
  }
  out
}

# Eligible-route fit charges: Pass 1 visits each support's cohort rows once in
# B-row chunks; when p > 0, ONE common QR on the never-treated cohort (plus one
# conditional alias-reconstruction scan) replaces the S + G block QRs, while
# p = 0 uses no QR; Pass 2 visits the G control blocks on the never-treated
# cohort in B-row chunks to accumulate forcing moments and K as squares.
or_baseline_complete_fit_charges <- function(st, specs, S, G, L, p, B, metadata, header, name_bytes,
                                             complete) {
  P2 <- p + 2
  generic <- or_baseline_fit_charges(S, G, L, p, st$N, st$n_A, B, metadata, header, name_bytes)
  n_k <- vapply(specs, function(s) st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1, numeric(1))
  rows <- sum(n_k); chunks <- sum(ceiling(n_k / B))
  n_c_vec <- st$cohort_end[complete$control] - st$cohort_start[complete$control] + 1
  n_c <- sum(n_c_vec)
  c_chunks <- sum(ceiling(n_c_vec / B))
  pass2_rows <- G * n_c; pass2_chunks <- G * c_chunks
  facts <- if (!is.null(complete$fact_reads)) complete$fact_reads else
    length(complete$dates) * complete$n_a + 2 * st$N + 8 * L
  multi_ctrl <- length(complete$control) > 1L
  common_qr <- (if (p > 0L) 256 * n_c * P2^2 + 256 * c_chunks * P2^3 + 512 * P2^3 + 256 * rows * p^2 else 0) +
    if (multi_ctrl) 256 * pass2_rows * P2^2 + 256 * (pass2_chunks + G) * P2^3 else 0
  work <- 256 * (rows + pass2_rows) * P2 + common_qr +
    256 * (G + S) * P2 + 128 * (metadata + header + L * P2) + 64 * facts
  scratch <- header + metadata + 128 * B * (P2 + 8) +
    (if (p > 0L) 64 * n_c + 64 * 12 * P2^2 else 0) + 64 * (S + G) * P2 + 4096 +
    if (multi_ctrl) 64 * (2 * S + 2 * G + 8) * P2^2 + 64 * G * p * (p + 1) else 0
  list(models = generic$models + unname(cp04_payload(complete)["entries"]) + 1024,
    scratch = scratch, work = work, support_row_visits = rows, support_chunks = chunks,
    pass2_row_visits = pass2_rows, pass2_chunks = pass2_chunks,
    common_qr_chunks = if (p > 0L) c_chunks else 0,
    fact_reads = facts, generic_work = generic$work, generic_scratch = generic$scratch)
}

or_baseline_complete_record <- function(charges, complete, verification_scan = FALSE,
                                        alias_rows = 0, common_qr_compressions = 0L) {
  list(route = "complete_support", n_support_row_visits = charges$support_row_visits,
    n_support_chunks = charges$support_chunks, n_pass2_row_visits = charges$pass2_row_visits,
    n_pass2_chunks = charges$pass2_chunks, n_common_qr_chunks = charges$common_qr_chunks,
    n_common_qr_compressions = as.integer(common_qr_compressions),
    n_endpoint_searches = 0, n_alias_rescans = if (isTRUE(verification_scan)) 1 else 0,
    n_alias_rows_checked = as.double(alias_rows),
    n_complete_dates = length(complete$dates), fact_reads = charges$fact_reads,
    replaced_generic_fit_work_ops = charges$generic_work,
    replaced_generic_fit_scratch_entries = charges$generic_scratch)
}

# CP08 stage 4: native complete-support batches. Kernels are resolved at each
# call boundary as for the IPW batch plan: the existing resolver verifies the
# current `cdid` DLL (identity, its nine routines and arities, live arithmetic)
# and the kernel's registered symbol and arity are taken from that same DLL; a
# cached address is reused only while the resolver returns the very same
# symbol table. No R fallback exists.
or_baseline_native_cache <- new.env(parent = emptyenv())

or_baseline_native_kernel <- function(name) {
  arity <- c(C_gmm_cp08_support_chunk = 12L, C_gmm_cp08_pair_rows = 9L, C_gmm_cp08_fit_pass1 = 15L,
    C_gmm_cp08_fit_pass2 = 15L, C_gmm_cp08_or_if_chunk = 17L, C_gmm_cp08_dr_if_prepare = 14L,
    C_gmm_cp08_dr_if_assemble = 17L, C_gmm_cp08_dr_treatment_batch = 19L)
  syms <- gmm_require_native_dll()
  owner <- environment(sys.function())
  cache <- get0("or_baseline_native_cache", envir = owner, mode = "environment", inherits = FALSE)
  hit <- if (is.environment(cache)) get0(name, envir = cache, inherits = FALSE) else NULL
  if (is.list(hit) && identical(hit$owner, owner) && gmm_sexp_aliases(hit$syms, syms)) return(hit$address)
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  si <- if (name %in% names(arity)) tryCatch(base::getNativeSymbolInfo(name, PACKAGE = dlls[["cdid"]]),
    error = function(e) NULL)
  if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
      !identical(as.integer(si$numParameters), unname(arity[name])))
    cdid_gmm_abort(paste0("Native symbol registration/forcing verification failed for `", name, "`."),
      c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
  if (is.environment(cache)) assign(name, list(owner = owner, syms = syms, address = si$address), envir = cache)
  si$address
}

# A kernel status (code, stage, index, row) raises the original route's own
# condition for that first failure; anything else refuses as integrity.
or_baseline_native_refuse <- function(status, abort) {
  code <- status[1L]
  if (identical(code, 1L)) abort("A complete-date posting slot does not hold its A-unit.", "integrity")
  if (identical(code, 2L)) cdid_ipw_owner_abort("Observed risk-set outcome endpoints contain NA, NaN, or Inf.",
    c("cdid_ipw_owner_endpoint_error", "cdid_ipw_enclosure_endpoint_error"))
  if (identical(code, 3L)) cp04_abort("Nonfinite observed pair differences.", "numeric")
  if (identical(code, 4L)) or_pooled_abort("The declared intercept column must equal one on every observed pair row.",
    "input", "baseline")
  if (identical(code, 5L)) cp04_abort("Nonfinite complete baseline link IF.", "numeric")
  if (identical(code, 7L)) cp04_abort("Nonfinite complete baseline DR frame IF.", "numeric")
  if (identical(code, 8L)) cp04_abort("Nonfinite complete baseline DR link IF.", "numeric")
  if (identical(code, 9L)) abort(paste("A complete-support IF product would accumulate more than one term with two",
    "nonzero factors; its BLAS-order identity is not established."), "identity")
  if (identical(code, 11L)) cdid_abort("Require finite p-row b and an explicit positive normalizer.",
    "cdid_or_wls_input_error")
  if (identical(code, 12L)) cdid_abort("Nonfinite normalized WLS Gram loading.",
    "cdid_or_wls_numeric_error")
  if (identical(code, 13L)) cdid_abort("Nonfinite compressed WLS factor.",
    "cdid_or_wls_numeric_error")
  abort(paste0("Complete-support native kernel status (", paste(status, collapse = ", "),
    ") is not a route condition."), "integrity")
}

or_baseline_native_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
  length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)

# Fit time (OR and DR alike): the supports of an eligible complete-support fit.
# Each (support, cohort chunk) is that of or_baseline_visit. One kernel call
# checks both posting slots of every chunk row, reads its two endpoints, forms
# the raw right-minus-left difference and checks the declared intercept column
# (the checks and errors of or_baseline_complete_pair_rows and
# or_pooled_support_add), and returns R's sums of w[A], w[obs] and w[obs] * dy
# with the rows [1 * sqrt(w), dy * sqrt(w)] that or_wls_summary compresses.
# The unchanged or_wls_compress then compresses the chunk and merges it into
# the support factor in the original sequential order, so each support's
# summary (factor, counts, merge accounting), rate and moments equal those of
# or_baseline_support. Supports run in specification order.
or_baseline_complete_supports <- function(reader, specs, recipe, ctrl, complete, scratch, z_col,
                                          abort = or_baseline_abort) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  po <- or_baseline_complete_panel(reader, complete, abort)
  if (!identical(recipe$p, 0L) || !isTRUE(recipe$intercept))
    abort("Native complete-support supports require the intercept-only recipe of their facts.", "integrity")
  B <- ctrl$max_chunk_rows
  # One chunk's kernel output (augmented rows, sums, status) within the fit scratch.
  ipw_owner_safe_sum(c(2 * min(B, st$N), 64), scratch, "complete-support native support chunk")
  kernel <- or_baseline_native_kernel("C_gmm_cp08_support_chunk")
  ld <- isTRUE(unname(capabilities("long.double")))
  y <- po$observations$y; post <- po$postings
  lapply(specs, function(spec) {
    acc <- or_pooled_support_start(recipe, spec$transitions)
    s <- acc$summary
    fac <- s$raw$factor; n_raw <- s$raw$n
    merge_rows <- s$accounting$max_chunk_merge_rows; n_comp <- s$accounting$n_compressions
    repeat_mass <- 0; n_repeat <- 0
    first <- st$cohort_start[spec$cohort]; last <- st$cohort_end[spec$cohort]
    dates <- as.integer(c(spec$pair$s_calendar_index, spec$pair$t_calendar_index))
    for (start in seq.int(first, last, by = B)) {
      res <- .Call(kernel, y, post$date_start, post$post_a_index, post$post_obs_index, st$cohort_order,
        st$a_index, st$w, frame$Z_a, as.integer(c(start, min(last, start + B - 1), z_col)), dates, ld,
        ctrl$max_matrix_entries)
      if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, abort)
      repeat_mass <- repeat_mass + res$sums[1L]; n_repeat <- n_repeat + res$n_A
      chunk <- or_wls_compress(res$augmented)
      stacked <- nrow(fac) + nrow(chunk)
      if (stacked > merge_rows) merge_rows <- stacked
      fac <- or_wls_compress(rbind(fac, chunk))
      n <- nrow(res$augmented)
      n_raw <- n_raw + n
      n_comp <- n_comp + 1L + if (stacked > 0L) 1L else 0L
      acc$weight_mass <- acc$weight_mass + res$sums[2L]; acc$n <- acc$n + n
      acc$sumY <- acc$sumY + res$sums[3L]
    }
    s$raw$factor <- fac; s$raw$n <- n_raw; s$n <- n_raw
    s$accounting <- list(input_rows = n_raw, augmented_cols = s$p + 2L, augmented_rows = nrow(fac),
      max_chunk_merge_rows = merge_rows, n_compressions = n_comp)
    acc$summary <- s
    or_baseline_support_rate(acc, repeat_mass, n_repeat, st$N)
  })
}

# Two-pass common-design complete-support state (OR and DR alike):
# - Pass 1 visits each support's cohort chunks via C_gmm_cp08_fit_pass1,
#   checking posting slots, reading endpoints, subtracting raw differences
#   BEFORE summation, verifying the intercept column, and accumulating support
#   sums and (when p > 0) ystar = sum_b (v_b L_b / D) dy_b on the control cohort.
# - When p > 0, ONE common Householder QR summary of (z, ystar, w) is built on
#   the control cohort; when p = 0, no QR is performed.
# - or_pooled_fit_common solves alpha_b and beta before Pass 2, then invokes
#   C_gmm_cp08_fit_pass2 on the G control blocks to accumulate forcing moments
#   and K as non-negative squares and reconstructs the exact pooled fit/summary.
or_baseline_complete_state <- function(reader, specs, pools, links, recipe, ctrl, complete,
                                       scratch, rank_tol, frame_index = seq_along(recipe$columns),
                                       abort = or_baseline_abort) {
  st <- reader$state; N <- st$N; B <- ctrl$max_chunk_rows; p <- recipe$p
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  po <- or_baseline_complete_panel(reader, complete, abort)
  if (!isTRUE(recipe$intercept) || !is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 0L || p != floor(p))
    abort("Native complete-support state requires the intercept recipe of its facts.", "integrity")
  y <- po$observations$y; post <- po$postings
  icpt_col <- as.integer(frame_index[recipe$intercept_index])
  slope_cols <- if (p > 0L) as.integer(frame_index[recipe$slope_index]) else integer(0)
  multi_ctrl <- length(complete$control) > 1L || any(lengths(lapply(pools, `[[`, "support_ids")) != 1L)
  ctrl_first <- as.integer(st$cohort_start[complete$control[1L]])
  ctrl_last <- as.integer(st$cohort_end[complete$control[1L]])
  n_c <- ctrl_last - ctrl_first + 1L
  S <- length(specs); G <- length(pools)
  ipw_owner_safe_sum(c(min(B, N) * (p + 2), if (p > 0L) n_c else 0, 5 * S * (p + 1), 64),
    scratch, "complete-support native fit pass 1")
  sup_mat <- vapply(specs, function(s) as.integer(c(st$cohort_start[s$cohort], st$cohort_end[s$cohort],
    s$pair$s_calendar_index, s$pair$t_calendar_index)), integer(4))
  L_vec <- vapply(pools, `[[`, numeric(1), "transitions")
  ctrl_sids <- if (multi_ctrl) integer(0) else vapply(pools, function(pl) pl$support_ids[1L], integer(1))
  ystar_w <- numeric(S)
  if (!multi_ctrl && p > 0L) {
    D <- sum(L_vec^2)
    ystar_w[ctrl_sids] <- L_vec / D
  }
  ld <- isTRUE(unname(capabilities("long.double")))
  # Rare original-factor validation is operational checking, never another fit.
  # A squared raw augmented entry can overflow even when stable QR is finite;
  # the native pass requests the original check instead of refusing that fact.
  checked_supports <- logical(S)
  checked_cache <- if (multi_ctrl) vector("list", S) else NULL
  validation_work <- 0
  validation_count <- validation_qr <- 0L
  validate_raw_factor <- function(sid) {
    if (checked_supports[sid]) return(invisible(NULL))
    nk <- st$cohort_end[specs[[sid]]$cohort] - st$cohort_start[specs[[sid]]$cohort] + 1
    P2 <- p + 2
    more <- 1024 * nk * P2^2 * (1 + ceiling(log2(st$n_A + 1)))
    validation_work <<- ipw_owner_safe_sum(c(validation_work, more), ctrl$max_work_ops,
      "conditional original-factor validation work")
    cp04_budget(ctrl, complete$fit_retained_reserved,
      c(scratch, 128 * B * P2 + 64 * P2^2),
      c(complete$fit_work_reserved, validation_work), 0, "conditional original-factor validation")
    checked <- if (identical(abort, dr_baseline_abort)) {
      dr_baseline_support(reader, specs[[sid]], recipe, frame_index, ctrl, complete = complete)
    } else or_baseline_support(reader, specs[[sid]], recipe, ctrl, complete = complete)
    validation_count <<- validation_count + 1L
    validation_qr <<- validation_qr + checked$design$summary$accounting$n_compressions
    checked_supports[sid] <<- TRUE
    if (multi_ctrl) checked_cache[[sid]] <<- checked
    invisible(NULL)
  }
  p1 <- .Call(or_baseline_native_kernel("C_gmm_cp08_fit_pass1"), y, post$date_start,
    post$post_a_index, post$post_obs_index, st$cohort_order, st$a_index, st$w, frame$Z_a,
    sup_mat, slope_cols, ystar_w, as.integer(c(icpt_col, B, ctrl_first, ctrl_last)), ld,
    ctrl$max_matrix_entries, validate_raw_factor)
  if (!identical(p1$status[1L], 0L)) or_baseline_native_refuse(p1$status, abort)
  supports <- lapply(seq_len(S), function(s) {
    acc <- list(L = specs[[s]]$transitions, weight_mass = p1$sums[2L, s],
      n = as.double(p1$n_counts[2L, s]), sumZ = if (p > 0L) p1$sumZ[, s] else numeric(0),
      sumY = p1$sums[3L, s], colmax = if (p > 0L) p1$colmax[, s] else numeric(0))
    or_baseline_support_rate(acc, p1$sums[1L, s], as.double(p1$n_counts[1L, s]), N)
  })
  visit_rows <- function(sid, visit) {
    spec <- specs[[sid]]
    first <- st$cohort_start[spec$cohort]; last <- st$cohort_end[spec$cohort]
    for (start in seq.int(first, last, by = B)) {
      rows <- st$cohort_order[seq.int(start, min(last, start + B - 1L))]
      Z_sub <- frame$Z_a[st$a_index[rows], frame_index, drop = FALSE]
      colnames(Z_sub) <- recipe$columns
      visit(Z_sub)
    }
  }
  if (multi_ctrl) {
    ctrl_u <- sort(unique(as.integer(unlist(lapply(pools, `[[`, "support_ids")))))
    for (sid in ctrl_u) {
      supports[[sid]] <- if (!is.null(checked_cache[[sid]])) checked_cache[[sid]] else
        if (identical(abort, dr_baseline_abort)) {
          dr_baseline_support(reader, specs[[sid]], recipe, frame_index, ctrl, complete = complete)
        } else or_baseline_support(reader, specs[[sid]], recipe, ctrl, complete = complete)
    }
    ans <- or_baseline_assemble(recipe, supports, pools, links, rank_tol, N, visit_rows)
    return(c(ans, list(
      common_qr_compressions = if (p > 0L && !is.null(ans$shared$summary))
        ans$shared$summary$accounting$n_compressions else 0L,
      original_factor_validation_work = validation_work,
      original_factor_validations = validation_count,
      original_factor_qr_compressions = validation_qr
    )))
  }
  common_part <- NULL
  if (p > 0L) {
    ord_c <- st$cohort_order[seq.int(ctrl_first, ctrl_last)]
    ystar <- p1$ystar
    common_chunk <- NULL
    for (lo in seq.int(1L, n_c, by = B)) {
      idx <- seq.int(lo, min(n_c, lo + B - 1L))
      r_c <- ord_c[idx]
      z_c <- frame$Z_a[st$a_index[r_c], slope_cols, drop = FALSE]
      colnames(z_c) <- recipe$slope_columns
      next_chunk <- or_pooled_wls_block_summary(z_c, ystar[idx], st$w[r_c])
      common_chunk <- if (is.null(common_chunk)) next_chunk else
        or_pooled_wls_merge_chunks(list(common_chunk, next_chunk))
    }
    common_part <- or_pooled_wls_partial_block(common_chunk)
  }
  pass2_fn <- function(keep, alpha, beta, e_bar) {
    blk_mat <- vapply(pools, function(pl) {
      pr <- specs[[pl$support_ids[1L]]]$pair
      as.integer(c(pr$s_calendar_index, pr$t_calendar_index))
    }, integer(2))
    v_vec <- vapply(ctrl_sids, function(s) supports[[s]]$v, numeric(1))
    blk_par <- rbind(L_vec, v_vec, as.double(alpha), as.double(e_bar))
    kept_cols <- if (length(keep)) slope_cols[keep] else integer(0)
    p2 <- .Call(or_baseline_native_kernel("C_gmm_cp08_fit_pass2"), y, post$date_start,
      post$post_a_index, post$post_obs_index, st$cohort_order, st$a_index, st$w,
      frame$Z_a, blk_mat, blk_par, as.integer(kept_cols), as.double(beta),
      if (length(keep)) p1$ystar else numeric(0),
      as.integer(c(ctrl_first, ctrl_last, B)), ctrl$max_matrix_entries)
    if (!identical(p2$status[1L], 0L)) or_baseline_native_refuse(p2$status, abort)
    p2
  }
  blocks <- lapply(seq_len(G), function(b) list(name = paste0("b", b),
    support_ids = pools[[b]]$support_ids, L = pools[[b]]$transitions))
  fitted <- or_pooled_fit_common(recipe, lapply(supports, function(u)
      list(scale = u$v, weight_mass = u$pair_mass, n = u$n_pair, sumZ = u$design$sumZ,
        sumY = u$design$sumY, colmax = u$design$colmax)),
    blocks, lapply(links, function(l) list(regression = l$regression,
      support_id = l$support_id)), rank_tol, N, visit_rows, common_part, pass2_fn, "baseline")
  p_keep <- length(fitted$keep)
  beta_total <- matrix(0, p_keep, S); energy_total <- numeric(S)
  regressions <- fitted$regressions
  moments <- fitted$pass2$moments; g_slope <- fitted$pass2$g_slope
  for (b in seq_len(G)) {
    sid <- ctrl_sids[b]; rec <- regressions[[b]]; u <- supports[[sid]]
    scale <- if (u$rate > ipw_clip_eps) -1 / (N * u$rate^2) else 0
    g1 <- moments[1L, b]; e2 <- moments[2L, b]
    gs <- if (p_keep > 0L) g_slope[, b] else numeric(0)
    fa <- scale * g1
    fb <- matrix(scale * (gs - rec$mean_d * g1), p_keep, 1L)
    fe <- scale * e2
    cp04_finite(c(fa, fb, fe), "baseline regression rate forcing")
    regressions[[b]]$forcing_alpha <- fa
    regressions[[b]]$forcing_beta <- fb
    regressions[[b]]$forcing_energy <- fe
    beta_total[, sid] <- beta_total[, sid] + fb[, 1L]
    energy_total[sid] <- energy_total[sid] + fe
  }
  for (i in seq_along(links)) {
    l <- links[[i]]; u <- supports[[l$support_id]]; rec <- regressions[[l$regression]]
    m <- or_pooled_link_moments(recipe, fitted$keep, u$pair_mass, u$design$sumZ,
      u$design$sumY, rec$L, rec$alpha, fitted$beta)
    rc <- ipw_clip_floor(u$rate); ra <- u$rate > ipw_clip_eps
    d <- u$pair_mass / (N * rc)
    links[[i]] <- c(l, list(mu = m$mu, c = m$c, d = d, numerator = m$mu * d,
      rate_loading = if (ra) -m$residual / (N * rc^2 * d) else 0,
      log_weight_derivative = if (ra) -1 / rc else 0))
    cp04_finite(c(m$mu, m$c, d, links[[i]]$rate_loading), "baseline link moments")
  }
  beta_named <- fitted$beta; names(beta_named) <- fitted$basis$kept_columns
  list(
    supports = lapply(supports, function(u) list(rate = u$rate, v = u$v, massA = u$massA,
      pair_mass = u$pair_mass, repeat_mass = u$repeat_mass, n_pair = u$n_pair,
      n_repeat = u$n_repeat, sumZ = u$design$sumZ, sumY = u$design$sumY,
      colmax = u$design$colmax, L = u$design$L)),
    regressions = regressions,
    links = links,
    shared = list(p = p_keep, keep = fitted$keep, beta = fitted$beta,
      coefficients = beta_named, fit = fitted$fit, summary = fitted$summary, basis = fitted$basis,
      forcing_beta = beta_total, forcing_energy = energy_total,
      control_support_ids = sort(unique(unlist(lapply(pools, `[[`, "support_ids")))),
      rss = fitted$rss, objective = fitted$objective,
      verification_scan = fitted$verification_scan),
    common_qr_compressions = if (p > 0L) fitted$summary$accounting$n_compressions else 0L,
    original_factor_validation_work = validation_work,
    original_factor_validations = validation_count,
    original_factor_qr_compressions = validation_qr
  )
}

# Replay: compact descriptors of one requested link batch of an eligible
# complete-support fit, in the original pre$sids / pre$gids / link orders,
# derived once per replay from the validated sealed backend and sealed. The
# plan bound plus one native chunk's live entries (B x L output, B x block
# intercept columns, B x gathered dates and their state bytes) are asserted
# inside the unchanged per-chunk concurrent scratch charge, which dominates
# them, before any read.
or_baseline_native_live <- function(B, nl, ng, nd, p = 0L)
  B * (nl + ng + nd + 2 * p) + ceiling(B * nd / 8) + 256 + 4 * p

or_baseline_native_plan_bound <- function(ns, ng, nl, ids, p = 0L, nb = 0L)
  4096 + 48 * (ns + ng + nl) + sum(ceiling(nchar(ids, type = "bytes") / 2)) + 8 * nl +
    if (p > 0L) 2 * p * p + (2 + p) * (nb + ng) + p * (ns + nl) + 8 * (nb + ns + p) + 512 else 0

or_baseline_native_support_arrays <- function(backend, sids, extra_dates = numeric()) {
  sp <- backend$support_specs[sids]; su <- backend$supports[sids]
  pd <- vapply(sp, function(s) as.double(c(s$pair$s_calendar_index, s$pair$t_calendar_index)), numeric(2))
  dates <- sort(unique(c(as.vector(pd), extra_dates)))
  sup <- rbind(vapply(sp, function(s) as.double(s$cohort), numeric(1)), matrix(match(pd, dates), 2L))
  list(dates = dates, sup = sup, sup_par = vapply(su, function(u) c(u$rate, u$massA, u$v), numeric(3)))
}

or_baseline_native_int <- function(x) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x != floor(x)) || any(abs(x) > .Machine$integer.max))
    return(NULL)
  storage.mode(x) <- "integer"
  x
}

or_baseline_native_slope_desc <- function(backend, pre, tail, kept_cols, N, abort = or_baseline_abort) {
  p <- backend$shared$p
  if (p < 1L) return(NULL)
  beta_gids <- pre$beta_gids
  beta_sids <- unlist(lapply(backend$regression_specs[beta_gids], `[[`, "support_ids"), use.names = FALSE)
  beta_reps <- vapply(backend$regression_specs[beta_gids], function(s) length(s$support_ids), integer(1))
  beta_blk <- if (all(beta_reps >= 1L)) as.double(match(beta_sids, pre$sids)) else NA_real_
  beta_par_1 <- vapply(backend$regressions[beta_gids], function(r)
    as.double(c(r$L, r$alpha, r$mean_d)), numeric(2 + p))
  beta_par <- if (all(beta_reps == 1L)) beta_par_1 else
    beta_par_1[, rep(seq_along(beta_gids), beta_reps), drop = FALSE]
  ctrl_sids <- backend$shared$control_support_ids
  ctrl_sup <- match(ctrl_sids, pre$sids)
  forcing_beta <- matrix(as.double(backend$shared$forcing_beta[, ctrl_sids, drop = FALSE]), p, length(ctrl_sids))
  R_mat <- matrix(as.double(backend$shared$fit$R), p, p)
  pivot <- or_baseline_native_int(backend$shared$fit$pivot)
  beta_blk <- or_baseline_native_int(beta_blk)
  ctrl_sup <- or_baseline_native_int(ctrl_sup)
  kept_cols <- or_baseline_native_int(kept_cols)
  if (is.null(pivot) || is.null(beta_blk) || is.null(ctrl_sup) || is.null(kept_cols) ||
      anyNA(beta_par) || anyNA(forcing_beta) || anyNA(R_mat))
    abort("Native slope IF descriptors do not match the fitted shared slope.", "integrity")
  list(list(NULL, kept_cols), as.double(backend$shared$beta), R_mat, pivot, as.double(N),
    beta_blk, beta_par, ctrl_sup, forcing_beta, tail)
}

or_baseline_native_plan <- function(backend, links, B, pre = NULL) {
  validate_or_baseline_backend(backend)
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) || B < 1 ||
      B != floor(B) || B > backend$control$max_chunk_rows)
    or_baseline_abort("Invalid native IF link batch or chunk rows.")
  if (is.null(pre)) pre <- or_baseline_if_budget(backend, B, links, FALSE)
  p <- backend$shared$p
  if (is.null(backend$complete_support) || !is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      p < 0L || p != floor(p) || !isTRUE(backend$outcome_recipe$intercept) || !is.null(backend$reference_ref) ||
      (p == 0L && length(pre$beta_gids)) || (p > 0L && !length(pre$beta_gids)))
    or_baseline_abort("The native IF requires the facts of a complete-support fit.", "integrity")
  sa <- or_baseline_native_support_arrays(backend, pre$sids)
  blk <- vapply(backend$regression_specs[pre$gids], function(s)
    if (length(s$support_ids) == 1L) as.double(match(s$support_ids, pre$sids)) else NA_real_, numeric(1))
  blk_par <- vapply(backend$regressions[pre$gids], function(r) c(r$alpha, r$mass,
    if (length(r$forcing_alpha) == 1L) r$forcing_alpha else NA_real_), numeric(3))
  lk <- backend$links[links]
  lnk <- rbind(match(vapply(lk, `[[`, integer(1), "support_id"), pre$sids),
    match(vapply(lk, `[[`, integer(1), "regression"), pre$gids))
  lnk_par <- vapply(lk, function(l) c(l$mu, l$d, if (length(l$rate_loading) == 1L) l$rate_loading else NA_real_),
    numeric(3))
  ids <- vapply(lk, `[[`, character(1), "link_id")
  dates <- or_baseline_native_int(sa$dates); sup <- or_baseline_native_int(sa$sup)
  blk <- or_baseline_native_int(blk); lnk <- or_baseline_native_int(lnk)
  if (is.null(dates) || is.null(sup) || is.null(blk) || is.null(lnk) || anyNA(blk_par) || anyNA(lnk_par) ||
      anyNA(sa$sup_par))
    or_baseline_abort("Native IF descriptors do not match the fitted supports, blocks and links.", "integrity")
  slope <- NULL
  if (p > 0L) {
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    blk_lm <- vapply(backend$regressions[pre$gids], function(r) as.double(c(r$L, r$mean_d)), numeric(1 + p))
    lnk_c <- matrix(vapply(lk, function(l) as.double(l$c), numeric(p)), p, length(links))
    if (anyNA(blk_lm) || anyNA(lnk_c))
      or_baseline_abort("Native IF descriptors do not match the fitted block and link slopes.", "integrity")
    kept_cols <- backend$outcome_recipe$slope_index[backend$shared$keep]
    slope <- or_baseline_native_slope_desc(backend, pre, list(blk_lm, lnk_c), kept_cols,
      reader$state$N, or_baseline_abort)
  }
  bound <- or_baseline_native_plan_bound(length(pre$sids), length(pre$gids), length(links), ids, p, length(pre$beta_gids))
  ipw_owner_safe_sum(c(bound, or_baseline_native_live(B, length(links), length(pre$gids), length(dates), p)),
    pre$scratch, "complete-support native OR IF reservation")
  plan <- list(kind = "or_baseline_native_plan", backend_ref = backend$ref, links = as.integer(links), B = B,
    sids = pre$sids, gids = pre$gids, dates = dates, sup = sup, sup_par = sa$sup_par, blk = blk,
    blk_par = blk_par, lnk = lnk, lnk_par = lnk_par, dimnames = list(slope, ids), entries = 0)
  plan$entries <- unname(cp04_payload(plan)["entries"])
  if (plan$entries > bound) or_baseline_abort("Native IF plan bound underestimated.", "bound")
  ipw_owner_seal(plan, "or_baseline_native_plan")
}

# A supplied plan must be this backend's sealed plan of exactly this link
# batch and support/block/model sets, for chunks no larger than its rows.
or_baseline_native_check <- function(plan, backend, links, B, pre, kind, abort) {
  if (!is.list(plan) || !identical(plan$kind, kind) || !ipw_owner_seal_ok(plan, kind) ||
      !identical(plan$backend_ref, backend$ref) || !identical(plan$links, as.integer(links)) ||
      !isTRUE(B <= plan$B) || !identical(plan$sids, pre$sids) || !identical(plan$gids, pre$gids) ||
      !identical(plan$tids, pre$tids))
    abort("The native IF plan does not belong to this backend, link batch and chunk.", "integrity")
  invisible(plan)
}

# One eligible non-detailed OR chunk: the kernel fills the B x L IF in the
# requested link order; R keeps the scale, output record and output bound.
or_baseline_native_if_chunk <- function(backend, reader, rows, links, scale, pre, plan = NULL) {
  st <- reader$state; B <- length(rows); p <- backend$shared$p
  if (is.null(plan)) plan <- or_baseline_native_plan(backend, links, B, pre) else
    or_baseline_native_check(plan, backend, links, B, pre, "or_baseline_native_plan", or_baseline_abort)
  ipw_owner_safe_sum(c(plan$entries, or_baseline_native_live(B, length(links), length(plan$gids),
    length(plan$dates), p)), pre$scratch, "complete-support native OR IF chunk")
  po <- or_baseline_complete_panel(reader, backend$complete_support)
  post <- po$postings
  dn <- plan$dimnames
  if (!is.null(dn) && !is.null(dn[[1L]])) {
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    dn[[1L]][[1L]][[1L]] <- frame$Z_a
  }
  res <- .Call(or_baseline_native_kernel("C_gmm_cp08_or_if_chunk"), po$observations$y, post$date_start,
    post$post_a_index, post$post_obs_index, as.integer(rows), st$cohort_index, st$a_index, st$w, plan$dates,
    plan$sup, plan$sup_par, plan$blk, plan$blk_par, plan$lnk, plan$lnk_par, dn,
    backend$control$max_matrix_entries)
  if (!identical(res$status[1L], 0L)) or_baseline_native_refuse(res$status, or_baseline_abort)
  phi <- res$phi
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = length(pre$beta_gids), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_baseline_abort("IF output bound underestimated.", "bound")
  out
}

or_baseline_support_rate <- function(acc, repeat_mass, n_repeat, N) {
  if (!is.finite(repeat_mass) || !is.finite(acc$weight_mass) || repeat_mass <= 0 ||
      acc$weight_mass <= 0 || repeat_mass / N <= 0 || acc$weight_mass / repeat_mass <= 0)
    or_baseline_abort("Positive repeat and actual-pair cohort support required.", "support")
  rate <- acc$weight_mass / repeat_mass
  cp04_finite(c(rate, acc$sumZ, acc$sumY), "baseline compact moments")
  # Raw strictly positive rate (and its IF) kept; only the inverse is floored.
  list(design = acc, rate = rate, v = 1 / ipw_clip_floor(rate), massA = repeat_mass / N,
    pair_mass = acc$weight_mass, repeat_mass = repeat_mass, n_pair = acc$n,
    n_repeat = n_repeat)
}

# q = w/max(r_s, eps) on support s, so t_s = -1{r_s > eps}/r_s on its rows. Per (block, support):
# C_alpha = H(q r_b t), C_beta = H(q (d - dbar_b) r_b t) (uncentered d without
# intercept) and C_energy = H(q r_b^2 t). A single-cohort block's C_alpha is zero
# at the root, its C_beta is not: it rescales the block within the joint fit.
or_baseline_forcing <- function(recipe, fitted, supports, pools, N) {
  p <- length(fitted$keep); S <- length(supports)
  beta_total <- matrix(0, p, S); energy_total <- numeric(S)
  regressions <- fitted$regressions
  for (b in seq_along(pools)) {
    sids <- pools[[b]]$support_ids; rec <- regressions[[b]]
    fa <- fe <- numeric(length(sids)); fb <- matrix(0, p, length(sids))
    for (k in seq_along(sids)) {
      u <- supports[[sids[k]]]
      m <- or_pooled_factor_moments(recipe, u$design$summary, fitted$keep, rec$alpha,
        fitted$beta)
      scale <- if (u$rate > ipw_clip_eps) -1 / (N * u$rate^2) else 0
      if (recipe$intercept) {
        fa[k] <- scale * m$g[1L]
        fb[, k] <- scale * (m$g[-1L] - rec$mean_d * m$g[1L])
      } else fb[, k] <- scale * m$g
      fe[k] <- scale * m$e2
    }
    cp04_finite(c(fa, fb, fe), "baseline regression rate forcing")
    regressions[[b]]$forcing_alpha <- fa
    regressions[[b]]$forcing_beta <- fb
    regressions[[b]]$forcing_energy <- fe
    beta_total[, sids] <- beta_total[, sids] + fb
    energy_total[sids] <- energy_total[sids] + fe
  }
  list(regressions = regressions, beta = beta_total, energy = energy_total)
}

# Charges reserved before any scan. Retained: compact supports, block records
# and forcing, ONE pooled summary and fit, basis map and link moments (support
# factors die with fitting; no full-F masks). Work: support scans, one reserved
# verification rescan, merges, proposal/map, forcing and link moments.
or_baseline_fit_charges <- function(S, G, L, p, N, n_A, B, metadata, header, name_bytes) {
  P2 <- p + 2; log_a <- 1 + ceiling(log2(n_A + 1))
  models <- 16384 + 1024 * (S + G + L) + 64 * (S + G + 4) * P2^2 + 32 * G * S * P2 +
    64 * L * P2 + 64 * (2 * S + 8 * G + 16) * (P2 + name_bytes)
  scratch <- header + metadata + 128 * B * (P2 + 8) + 64 * (2 * S + 2 * G + 8) * P2^2 +
    64 * G * p * (p + 1)
  scan <- 256 * S * N * P2^2 * log_a
  work <- 2 * scan + 256 * S * N * P2 * (p + G + 2) + 256 * G * (S + 1) * P2^3 +
    256 * (G + S) * P2^2 * (G + 2) + 128 * (metadata + header + L * P2)
  list(models = models, scratch = scratch, work = work)
}

# Fit-time state from finished supports. visit_rows(sid, visit) replays bounded
# observed Z_a chunks of support sid, used only for exact alias reconstruction.
or_baseline_assemble <- function(recipe, supports, pools, links, rank_tol, N, visit_rows) {
  blocks <- lapply(seq_along(pools), function(b) list(name = paste0("b", b),
    support_ids = pools[[b]]$support_ids, L = pools[[b]]$transitions))
  fitted <- or_pooled_fit_blocks(recipe, lapply(supports, function(u)
      list(summary = u$design$summary, scale = u$v, colmax = u$design$colmax)),
    blocks, lapply(links, function(l) list(regression = l$regression,
      support_id = l$support_id)), rank_tol, N, visit_rows, "baseline")
  forcing <- or_baseline_forcing(recipe, fitted, supports, pools, N)
  for (i in seq_along(links)) {
    l <- links[[i]]; u <- supports[[l$support_id]]; rec <- forcing$regressions[[l$regression]]
    m <- or_pooled_link_moments(recipe, fitted$keep, u$pair_mass, u$design$sumZ,
      u$design$sumY, rec$L, rec$alpha, fitted$beta)
    rc <- ipw_clip_floor(u$rate); ra <- u$rate > ipw_clip_eps
    d <- u$pair_mass / (N * rc)
    # Retain the ratio's own rate loading (zero at the root), including centering;
    # a floored rate (r <= eps) has no rate path.
    links[[i]] <- c(l, list(mu = m$mu, c = m$c, d = d, numerator = m$mu * d,
      rate_loading = if (ra) -m$residual / (N * rc^2 * d) else 0,
      log_weight_derivative = if (ra) -1 / rc else 0))
    cp04_finite(c(m$mu, m$c, d, links[[i]]$rate_loading), "baseline link moments")
  }
  beta <- fitted$beta; names(beta) <- fitted$basis$kept_columns
  list(supports = lapply(supports, function(u) list(rate = u$rate, v = u$v, massA = u$massA,
      pair_mass = u$pair_mass, repeat_mass = u$repeat_mass, n_pair = u$n_pair,
      n_repeat = u$n_repeat, sumZ = u$design$sumZ, sumY = u$design$sumY,
      colmax = u$design$colmax, L = u$design$L)),
    regressions = forcing$regressions, links = links,
    shared = list(p = length(fitted$keep), keep = fitted$keep, beta = fitted$beta,
      coefficients = beta, fit = fitted$fit, summary = fitted$summary, basis = fitted$basis,
      forcing_beta = forcing$beta, forcing_energy = forcing$energy,
      control_support_ids = sort(unique(unlist(lapply(pools, `[[`, "support_ids")))),
      rss = fitted$rss, objective = fitted$objective,
      verification_scan = fitted$verification_scan))
}

# One block per unique (actual pair, sorted approved pool) within the backend's
# single recipe, whatever its consumer count. Canonical pair-major/pool order,
# independent of link order and target sharing; links keep the supplied order.
# With D60 reference recipes the scientific identity also carries the recipe:
# blocks are (pair, pool, recipe) and supports (pair, cohort, recipe), each
# once with multiplier one; identical recipes are one recipe (one block).
#
# CP08 construction: the pair, pool and support identities are found through
# invocation-local buckets keyed by their calendar/cohort/recipe values instead
# of scanning every accumulated pair, pool and support per record. A key only
# selects candidates: within a bucket the original identical(pair) and
# identical(pool spec) checks decide, in first-appearance order, so attributes,
# ordered controls and reference recipes keep their original identity, and a
# support candidate is accepted only under the original (pair, cohort,
# recipe) predicate (any other case runs the original scan). The pair sort,
# pool radix sort and remapping, support order and supplied link order are the
# original statements. The buckets hold at most L pair keys, L pool keys, S
# support keys and one integer per record; both callers admitted the catalog
# metadata header (4096 + 1024 L (J + 1) entries, 128 L^2 (J + 1) work) before
# this call, which dominates them, and they are dropped on return.
or_baseline_catalog_key <- function(...) {
  paste(vapply(list(...), function(v) paste(as.character(unlist(v)), collapse = ","), character(1)),
    collapse = "|")
}

or_baseline_catalog <- function(records, cohorts, recipe, reference = NULL) {
  dated <- !is.null(reference)
  L <- length(records)
  pairs <- pools <- list()
  links <- vector("list", L)
  pair_bucket <- new.env(hash = TRUE, parent = emptyenv())
  pool_bucket <- new.env(hash = TRUE, parent = emptyenv())
  for (ell in seq_len(L)) {
    x <- records[[ell]]
    kp <- or_baseline_catalog_key(x$pair$s_calendar_index, x$pair$t_calendar_index)
    cand <- get0(kp, envir = pair_bucket, inherits = FALSE)
    pi <- 0L
    for (ci in cand) if (identical(pairs[[ci]], x$pair)) { pi <- ci; break }
    if (!pi) {
      pairs[[length(pairs) + 1L]] <- x$pair; pi <- length(pairs)
      if (pi) assign(kp, c(cand, pi), envir = pair_bucket)
    }
    spec <- if (dated) list(pair_id = pi, cohorts = x$controls, reference = x$reference) else
      list(pair_id = pi, cohorts = x$controls)
    kq <- or_baseline_catalog_key(pi, x$controls, if (dated) x$reference else NULL)
    cand <- get0(kq, envir = pool_bucket, inherits = FALSE)
    ri <- 0L
    for (ci in cand) if (identical(pools[[ci]], spec)) { ri <- ci; break }
    if (!ri) {
      pools[[length(pools) + 1L]] <- spec; ri <- length(pools)
      assign(kq, c(cand, ri), envir = pool_bucket)
    }
    links[[ell]] <- list(regression = as.integer(ri), pair_id = as.integer(pi),
      cohort = x$target, link_id = x$link_id, source_functional_key = x$source_functional_key,
      functional_key = x$functional_key)
    if (dated) links[[ell]]$reference <- x$reference
  }
  pair_bucket <- pool_bucket <- NULL
  transitions <- vapply(pairs, or_pooled_transitions, numeric(1), route = "baseline")
  po <- order(vapply(pairs, function(x) as.double(x$s_calendar_index), numeric(1)),
    vapply(pairs, function(x) as.double(x$t_calendar_index), numeric(1)))
  pairs <- pairs[po]; transitions <- transitions[po]; pair_map <- match(seq_along(po), po)
  for (i in seq_along(pools)) pools[[i]]$pair_id <- as.integer(pair_map[pools[[i]]$pair_id])
  qo <- order(vapply(pools, function(s) paste(sprintf("%012d", c(s$pair_id, s$cohorts, s$reference)),
    collapse = "|"), character(1)), method = "radix")
  pools <- pools[qo]; pool_map <- match(seq_along(qo), qo)
  for (i in seq_along(links)) {
    links[[i]]$pair_id <- as.integer(pair_map[links[[i]]$pair_id])
    links[[i]]$regression <- as.integer(pool_map[links[[i]]$regression])
  }
  # Reverse pair -> pool / link index lists (pool and link order), so each
  # pair's reference and cohort unions read only its own pools and links.
  pools_of_pair <- split(seq_along(pools), factor(vapply(pools, function(s) as.integer(s$pair_id), integer(1)),
    levels = seq_along(pairs)))
  links_of_pair <- split(seq_along(links), factor(vapply(links, function(l) as.integer(l$pair_id), integer(1)),
    levels = seq_along(pairs)))
  specs <- list()
  support_index <- new.env(hash = TRUE, parent = emptyenv())
  support_key <- function(pi, k, r) paste(as.character(c(pi, k, r)), collapse = ":")
  for (pi in seq_along(pairs)) {
    pidx <- pools_of_pair[[pi]]; lidx <- links_of_pair[[pi]]
    refs <- if (dated) sort(unique(c(unlist(lapply(pools[pidx], `[[`, "reference")),
      unlist(lapply(links[lidx], `[[`, "reference"))))) else NA_integer_
    for (r in refs) {
      pr_idx <- if (dated) pidx[vapply(pools[pidx], function(s) s$reference == r, logical(1))] else pidx
      lr_idx <- if (dated) lidx[vapply(links[lidx], function(s) s$reference == r, logical(1))] else lidx
      used <- sort(unique(c(unlist(lapply(pools[pr_idx], `[[`, "cohorts")),
        unlist(lapply(links[lr_idx], `[[`, "cohort")))))
      for (k in used) {
        sp <- list(pair = pairs[[pi]], pair_id = pi, cohort = k, transitions = transitions[pi])
        if (dated) sp$reference <- r
        specs[[length(specs) + 1L]] <- sp
        key <- support_key(pi, k, if (dated) r)
        # A repeated key is never a candidate: the original scan decides it.
        assign(key, if (exists(key, envir = support_index, inherits = FALSE)) NA_integer_ else length(specs),
          envir = support_index)
      }
    }
  }
  find_support <- function(pi, k, r = NULL) which(vapply(specs, function(s)
    s$pair_id == pi && s$cohort == k && (!dated || s$reference == r), logical(1)))
  support_of <- function(pi, k, r = NULL) {
    cand <- get0(support_key(pi, k, r), envir = support_index, inherits = FALSE)
    if (length(cand) == 1L && !is.na(cand)) {
      s <- specs[[cand]]
      if (isTRUE(s$pair_id == pi && s$cohort == k && (!dated || s$reference == r))) return(cand)
    }
    find_support(pi, k, r)
  }
  for (i in seq_along(pools)) {
    s <- pools[[i]]; pr <- pairs[[s$pair_id]]
    pools[[i]]$support_ids <- vapply(s$cohorts, function(k)
      as.integer(support_of(s$pair_id, k, s$reference)), integer(1))
    pools[[i]]$transitions <- transitions[s$pair_id]
    pools[[i]]$block_key <- if (dated) paste0("pair:", pr$s, "-", pr$t, "|pool:",
      paste(cohorts[s$cohorts], collapse = ","), "|reference-design:L_e*z(C^r)|",
      reference$recipes$key[s$reference], "|weights:W_g/r") else paste0("pair:", pr$s, "-", pr$t, "|pool:",
      paste(cohorts[s$cohorts], collapse = ","), "|", recipe$key, "|weights:w/r")
  }
  for (i in seq_along(links)) links[[i]]$support_id <- as.integer(
    support_of(links[[i]]$pair_id, links[[i]]$cohort, links[[i]]$reference))
  support_index <- NULL
  list(specs = specs, pools = pools, links = links)
}

validate_or_baseline_backend <- function(backend) {
  if (!inherits(backend, "cdid_or_baseline") ||
      !ipw_owner_seal_ok(backend, "or_baseline_backend"))
    or_baseline_abort("Expected the immutable private baseline OR backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) &&
    exists("cp04_reader", backend$reader_ref, inherits = FALSE))
      get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) validate_chained_reference_reader(or_baseline_reference(backend))
  if (!is.null(backend$complete_support)) or_baseline_complete_panel(reader, backend$complete_support)
  invisible(backend)
}

# Exact incremental cp04_payload(records) for the catalog loops of
# fit_or_baseline/fit_dr_baseline, which append one record per link: the
# container's slots (length + 8, plus any attributes) and every record's own
# payload, each record walked once when appended instead of re-walking every
# prefix three times. Charges are integer-valued doubles far below 2^53, so the
# sums are exact; the loops compare the total once with the recursive walk.
or_baseline_payload_append <- function(state, records, abort = or_baseline_abort) {
  n <- length(records)
  if (!is.list(records) || n != state$n + 1L)
    abort("Catalog metadata must be charged one appended record at a time.", "integrity")
  state$records <- state$records + cp04_payload(records[[n]])
  state$n <- n
  total <- c(entries = n + 8, char_bytes = 0) + state$records
  if (!is.null(attributes(records))) total <- total + cp04_payload(attributes(records), FALSE)
  state$total <- total
  state
}

or_baseline_payload_verify <- function(state, records, abort = or_baseline_abort) {
  if (!identical(state$n, length(records)) || !identical(state$total, cp04_payload(records)))
    abort("Incremental catalog metadata charge disagrees with its recursive payload.", "integrity")
  invisible(state$total)
}

or_baseline_reference <- function(backend) {
  if (is.null(backend$reference_ref)) return(NULL)
  if (!is.environment(backend$reference_ref) || !exists("reference_reader", backend$reference_ref, inherits = FALSE))
    or_baseline_abort("Invalid sealed reference reader binding.", "integrity")
  get("reference_reader", backend$reference_ref, inherits = FALSE)
}

# CP08 descriptor intake (OR and DR alike). A plain view list returns NULL and
# the view loop runs unchanged. A `cdid_cp04_descriptor_context` (see
# cp04_descriptor_context) must be the sealed object bound to THIS reader and
# its sealed frame/panel, for this estimator, without reference inputs; its
# complete-date facts are checked against the postings the reader reaches now,
# its one cohort descriptor and complete-date support counts are re-established
# from the sealed frame owner (never read from a label), and its recipes must
# be the fit's own checked recipes. It returns the row descriptors the view
# loop would read from M + Q heavy link views: pair descriptor, link ID and
# functional key, with the shared cohorts and support counts on the context.
or_baseline_descriptor_intake <- function(context, reader, frame, recipe, treatment = NULL, reference = NULL,
                                          abort = or_baseline_abort) {
  if (!inherits(context, "cdid_cp04_descriptor_context")) return(NULL)
  if (!is.list(context) || !ipw_owner_seal_ok(context, "cp04_descriptor_context") ||
      !identical(context$route, "descriptor_context"))
    abort("Expected the immutable CP04 descriptor context.", "integrity")
  if (!is.null(reference)) abort("A descriptor context carries no reference recipes.", "integrity")
  if (!identical(context$method, if (is.null(treatment)) "or" else "dr") ||
      !identical(context$reader_ref, reader$ref) || !identical(context$frame_ref, reader$frame_ref) ||
      !identical(context$panel_ref, reader$panel_ref) || !identical(frame$ref, reader$frame_ref) ||
      !ipw_owner_seal_ok(frame, "frame_owner"))
    abort("Descriptor context is bound to other sealed owners or another estimator.", "integrity")
  st <- reader$state
  po <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  post <- if (is.list(po)) po$postings else NULL
  if (!inherits(po, "cdid_ipw_panel_owner") || !ipw_owner_seal_ok(po, "panel_owner") ||
      !identical(as.integer(po$n_a_units), context$n_a) || !identical(st$N, context$N) || st$N != st$n_A ||
      !is.integer(post$post_a_index) || !is.integer(context$dates) || !length(context$dates) ||
      !identical(post$date_start[context$dates], context$date_start) ||
      !identical(post$date_end[context$dates], context$date_end) ||
      !identical(post$date_count[context$dates], context$date_count) ||
      !isTRUE(all(context$date_count == context$n_a)))
    abort("Descriptor context facts do not describe the reader's sealed postings.", "integrity")
  if (!identical(recipe, context$outcome_recipe) || !identical(treatment, context$treatment_recipe))
    abort("Descriptor context recipes differ from the fit's checked recipes.", "integrity")
  cohorts <- ipw_owner_make_cohort_descriptor(frame, context$g, context$K, require_support = TRUE)
  n_risk <- sum(cohorts$group_a_count)
  if (!identical(cohorts, context$cohorts) || !identical(context$support$by_cohort_pair_count, cohorts$group_a_count) ||
      !identical(context$support$n_observed_risk_units, n_risk) || !isTRUE(all(cohorts$group_a_count >= 1L)) ||
      !identical(match(context$g, st$cohorts), context$target) || !identical(match(context$K, st$cohorts), context$control))
    abort("Descriptor context cohort facts disagree with the sealed frame owner.", "integrity")
  desc <- context$descriptors
  if (!is.list(desc) || !length(desc) || length(desc) != context$n_catalog + context$n_auxiliary)
    abort("Descriptor context rows are malformed.", "integrity")
  desc
}

# `reference` (D60) is the sealed recipe-indexed reader of exactly these views
# (post catalog rows, then auxiliary descriptors); NULL keeps the static frame.
# `link_views` may instead be the CP08 descriptor context of the registry's
# eligible route (or_baseline_descriptor_intake): the same records, in the
# same order, enter the same catalog, route gate and fit.
fit_or_baseline <- function(reader, link_views, control, rank_tol, outcome_recipe,
                            retained_external_entries = 0, external_work_ops = 0, reference = NULL) {
  validate_cp04_reader(reader)
  if (missing(rank_tol)) or_baseline_abort("Supply explicit private rank_tol.")
  rank_tol <- or_baseline_rank_tol(rank_tol)
  ctrl <- validate_ipw_owner_control(control)
  if (!gmm_scalar_numeric(retained_external_entries, nonneg = TRUE) ||
      !gmm_scalar_numeric(external_work_ops, nonneg = TRUE))
    or_baseline_abort("Invalid external pipeline reservation.")
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (missing(outcome_recipe))
    or_baseline_abort("Supply the explicit private outcome recipe (intercept flag and Z_a roles).")
  if (!is.null(reference)) {
    validate_chained_reference_reader(reference)
    if (!is.list(link_views) || length(reference$view_recipe) != length(link_views))
      or_baseline_abort("The reference reader must index exactly the supplied link views.", "integrity")
  }
  recipe <- or_pooled_check_recipe(outcome_recipe, if (is.null(reference)) colnames(frame$Z_a) else reference$columns,
    "baseline")
  if (!is.list(link_views) || !length(link_views))
    or_baseline_abort("Supply actual private link views.")
  desc <- or_baseline_descriptor_intake(link_views, reader, frame, recipe, reference = reference)
  L <- if (is.null(desc)) length(link_views) else length(desc); J <- length(st$cohorts); p <- recipe$p
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries +
    retained_external_entries + if (is.null(reference)) 0 else reference$accounting$owned_entries
  header <- 4096 + 1024 * L * (J + 1)
  cp04_budget(ctrl, c(base, header), header, c(external_work_ops, 128 * L * L * (J + 1)), header,
    "baseline OR metadata")
  records <- list(); charge <- list(n = 0L, records = c(entries = 0, char_bytes = 0))
  for (ell in seq_len(L)) {
    if (is.null(desc)) {
      v <- link_views[[ell]]
      validate_ipw_link_view(v, frame_owner = frame)
      source <- v$functional_key
    } else {
      v <- desc[[ell]]
      source <- v$functional_key
    }
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:"))
      or_baseline_abort("Baseline OR requires the supplied comparison-catalog functional key.")
    key_bytes <- nchar(source, type = "bytes")
    cp04_budget(ctrl, c(base, header), header + 16 * key_bytes,
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * key_bytes), header + 16 * key_bytes,
      "baseline catalog key")
    if (is.null(desc)) {
      if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
        or_baseline_abort("Legacy nuisance candidates are unsupported.")
      if (!is.null(v$nuisance_ref)) {
        nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
        if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
          or_baseline_abort("Cross-pair or unverified nuisance sharing is unsupported.")
      }
      g <- v$cohorts$g; K <- v$cohorts$K
    } else {
      g <- link_views$g; K <- link_views$K
    }
    target <- match(g, st$cohorts)
    control_cohorts <- match(sort(K, method = "radix"), st$cohorts)
    if (anyNA(c(target, control_cohorts)))
      or_baseline_abort("Cohorts must belong to the actual common frame.", "support")
    records[[ell]] <- list(pair = v$pair, target = target, controls = control_cohorts,
      link_id = v$link_id, source_functional_key = source,
      functional_key = paste0("cmpcat:OR-baseline-shared-slope-WLS|", substring(source, 8L)))
    if (!is.null(reference)) {
      records[[ell]]$reference <- reference$view_recipe[ell]
      records[[ell]]$functional_key <- paste0(records[[ell]]$functional_key, "|reference:",
        reference$recipes$key[reference$view_recipe[ell]])
    }
    charge <- or_baseline_payload_append(charge, records)
    cp04_budget(ctrl, c(base, header), header + charge$total["entries"],
      c(external_work_ops, 128 * L * L * (J + 1) + 64 * charge$total["char_bytes"]),
      header + charge$total["entries"], "baseline catalog metadata")
  }
  or_baseline_payload_verify(charge, records)
  ids <- vapply(records, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids))
    or_baseline_abort("Distinct nonempty link IDs required.")
  catalog <- or_baseline_catalog(records, st$cohorts, recipe, reference)
  specs <- catalog$specs; pools <- catalog$pools; links <- catalog$links
  metadata <- unname(cp04_payload(list(specs, pools, links))["entries"])
  S <- length(specs); G <- length(pools)
  # CP08 stage 3: the route is chosen here from actual facts, before any scan;
  # neither route catches the other's failures.
  complete <- or_baseline_complete_facts(reader, recipe, records, specs, ctrl, reference,
    retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
    metadata_entries = header + metadata)
  charges <- if (is.null(complete)) or_baseline_fit_charges(S, G, L, p, st$N, st$n_A, min(st$N, ctrl$max_chunk_rows),
    metadata, header, sum(nchar(recipe$columns, type = "bytes"))) else
    or_baseline_complete_fit_charges(st, specs, S, G, L, p, min(st$N, ctrl$max_chunk_rows), metadata, header,
      sum(nchar(recipe$columns, type = "bytes")), complete)
  scratch <- charges$scratch
  # CP08/D78 work ownership. `own_work` is this fit's INCREMENTAL charge (every
  # support scan, pooled fit, factor and validation term of `charges$work`).
  # The returned `accounting$work_ops_bound` stays CUMULATIVE (the caller's
  # `external_work_ops` prefix plus own work), which is how every downstream
  # consumer of this backend reads it; the registry seam uses the declared
  # convention so that the prefix is charged exactly once.
  own_work <- ipw_owner_safe_sum(charges$work, label = "baseline OR own fit work")
  work <- ipw_owner_safe_sum(c(own_work, external_work_ops), ctrl$max_work_ops,
    "baseline OR pipeline fit work")
  bound <- metadata + charges$models
  cp04_budget(ctrl, c(base, bound), scratch, work, bound, "baseline OR fits")
  state <- if (is.null(complete)) {
    supports <- lapply(specs, function(s) or_baseline_support(reader, s, recipe, ctrl, reference, complete))
    ans <- or_baseline_assemble(recipe, supports, pools, links, rank_tol, st$N,
      function(sid, visit) or_baseline_visit(reader, specs[[sid]], ctrl,
        function(rows, S, obs, dy, Z) visit(Z), reference, complete))
    rm(supports)
    ans
  } else {
    complete$fit_work_reserved <- work
    complete$fit_retained_reserved <- c(base, bound)
    or_baseline_complete_state(reader, specs, pools, links, recipe, ctrl, complete, scratch, rank_tol,
      seq_along(recipe$columns), or_baseline_abort)
  }
  obj <- structure(list(reader_ref = reader$ref, control = ctrl, rank_tol = rank_tol,
    outcome_recipe = recipe, support_specs = specs, supports = state$supports,
    regression_specs = pools, regressions = state$regressions, shared = state$shared,
    links = state$links,
    law = paste0("OR:baseline:control-only:shared-one-period-slope:",
      if (recipe$intercept) "block-intercepts" else "no-intercept",
      if (is.null(reference)) ":actual-pair:w/r:common-frame-Z_a" else
        ":actual-pair:W_g/r:recipe-indexed-reference-design"),
    accounting = list(owned_entries = 0, borrowed_entries = base,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      retained_bound = base + bound, fit_scratch_bound = scratch, work_ops_bound = work,
      own_work_ops = own_work, work_ops_convention = "cumulative",
      n_cohort_summaries = S, n_rate_fits = S, n_regression_fits = 1L, n_regression_blocks = G,
      n_local_intercepts = if (recipe$intercept) G else 0L, n_shared_slopes = state$shared$p,
      n_slope_columns = p, n_aliased_slopes = length(state$shared$basis$aliased_columns),
      n_support_scans = S * (1L + state$shared$verification_scan),
      n_selection_fits = 0L, n_treatment_fits = 0L,
      full_F_masks_per_model = 0L, full_F_by_link_IF_retained = FALSE,
      inverse_matrices_retained = FALSE, dense_dummy_regression = FALSE,
      peak_rss_measured = FALSE)),
    class = "cdid_or_baseline")
  if (!is.null(reference)) {
    obj$reference_ref <- reference$ref
    obj$accounting$n_reference_recipes <- nrow(reference$recipes)
  }
  if (!is.null(complete)) {
    complete$fit_work_reserved <- complete$fit_retained_reserved <- NULL
    obj$complete_support <- complete
    obj$accounting$work_ops_bound <- ipw_owner_safe_sum(c(obj$accounting$work_ops_bound,
      state$original_factor_validation_work), ctrl$max_work_ops, "fit plus original-factor validation")
    obj$accounting$own_work_ops <- ipw_owner_safe_sum(c(obj$accounting$own_work_ops,
      state$original_factor_validation_work), ctrl$max_work_ops, "own fit plus original-factor validation")
    obj$accounting$complete_support <- or_baseline_complete_record(charges, complete,
      state$shared$verification_scan, state$shared$basis$rows_checked, state$common_qr_compressions)
    obj$accounting$complete_support$original_factor_validations <- state$original_factor_validations
    obj$accounting$complete_support$original_factor_qr_compressions <- state$original_factor_qr_compressions
    obj$accounting$complete_support$original_factor_validation_work <- state$original_factor_validation_work
  }
  obj["ref"] <- list(NULL)
  actual <- unname(cp04_payload(obj)["entries"])
  if (actual > bound) or_baseline_abort("Retained baseline bound underestimated.", "bound")
  cp04_budget(ctrl, c(base, actual), scratch, work, actual, "retained baseline OR")
  obj$accounting$owned_entries <- actual
  ipw_owner_seal(obj, "or_baseline_backend")
}

# Requested link blocks plus, whenever a shared slope exists (or details), EVERY
# block and control support: xi_beta is global to each link batch. Detail
# matrices use backend support/regression order; slopes follow the kept basis.
# Shares follow reader$state$cohorts. Details stay on the frame scale.
or_baseline_if_budget <- function(backend, B, links, details = FALSE) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- backend$shared$p; G <- length(backend$regression_specs); J <- length(st$cohorts)
  sets <- or_pooled_if_sets(backend, links, details)
  gids <- sets$gids; beta_gids <- sets$beta_gids; sids <- sets$sids
  nl <- length(links); ns <- length(sids)
  # D60 hits also carry their recipe weights (one slot per support row).
  width <- 4 * ns + 4 * (p + 1) * (length(beta_gids) + 1) + 4 * length(gids) +
    4 * nl * (p + 3) + (if (details) 2 * J + 8 else 8) + if (is.null(backend$reference_ref)) 0 else 4 * ns
  metadata <- 8192 + 1024 * (ns + G + nl + J) + B * frame$id_byte_metadata$max +
    8 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (width + 8) + metadata
  complete <- if (is.null(backend$complete_support)) NULL else
    or_baseline_complete_if_charges(backend, B, links, details, sets, J, metadata)
  if (is.null(complete)) {
    scratch <- 64 * B * (width + p + 8) + 128 * (p + 2)^2 + metadata
    work <- 128 * B * ((ns + 1) * (p + 2)^2 * (1 + ceiling(log2(st$n_A + 1))) +
      (length(beta_gids) + length(gids)) * (p + 1) * (p + ns) + J + nl * (p + 2)) +
      64 * metadata
  } else {
    scratch <- complete$scratch; work <- complete$work
  }
  cp04_budget(backend$control, c(backend$accounting$owned_entries,
    backend$accounting$borrowed_entries), scratch, work, output, "baseline OR IF chunk")
  out <- list(gids = gids, beta_gids = beta_gids, sids = sids, output = output, scratch = scratch,
    work = work)
  if (!is.null(complete)) out$complete <- complete$record
  out
}

# CP08 stage 3, eligible IF chunks: live entries per row and work of the
# support stage (rate columns) and the OR-core stage, from the chunk's actual
# supports, their distinct pairs and the requested blocks' own supports (sum of
# |support_ids_b|) instead of dense block-by-support products. `up` (B x ns) is
# the one dense matrix and is charged in full. Complete positions replace the
# endpoint searches; match() lookups over the support list are charged per chunk.
or_baseline_complete_core_charges <- function(backend, links, details, sets) {
  p <- backend$shared$p; rs <- backend$regression_specs
  ns <- length(sets$sids); nl <- length(links); nb <- length(sets$beta_gids); ng <- length(sets$gids)
  pairs <- length(unique(vapply(backend$support_specs[sets$sids], `[[`, integer(1), "pair_id")))
  inc <- function(g) sum(vapply(rs[g], function(s) length(s$support_ids), numeric(1)))
  I_beta <- inc(sets$beta_gids); I_g <- inc(sets$gids)
  # Support stage per row: rows/id/ci/w/A, up and one support's read
  # temporaries. Core per row: hits (obs, dy, z: one support per pair per row),
  # h_beta/xi_beta/rhs/hb, h_alpha and alpha, the resid cache (ix, q, r, d per
  # block-support and link entry), phi and direct, expression temporaries and,
  # with details, block scores and den/num/design.
  list(support_live = 5 + ns + 14, support_units = (ns + 1) * (p + 2)^2,
    core_live = (2 + p) * pairs + 4 * p + 2 * ng + (3 + p) * (I_beta + I_g + nl) + 2 * nl +
      4 * (p + 2) + if (details) p * nb + (2 + p) * nl + 3 + p else 0,
    core_units = (p + 1) * (p * (nb + ng) + I_beta + I_g),
    core_lookups = (I_beta + I_g + ng + 2 * nl + 1) * (ns + ng + 1), pairs = pairs,
    incidence = c(beta = I_beta, local = I_g))
}

# Concurrent payload is twice the chunk's live buffers (expression temporaries
# and copy-on-modify duplicates), plus the generic small-matrix and metadata
# terms; work keeps the generic per-row unit charge (128 per unit).
or_baseline_complete_if_charges <- function(backend, B, links, details, sets, J, metadata) {
  p <- backend$shared$p; nl <- length(links)
  cs <- or_baseline_complete_core_charges(backend, links, details, sets)
  # Plus the named phi copy, its finite check and core scaling; details: tau, eta, lambda.
  live <- cs$support_live + cs$core_live + 3 * nl + if (details) 2 * J + 5 else 0
  units <- cs$support_units + cs$core_units + J + nl * (p + 2)
  list(scratch = 2 * B * live + 128 * (p + 2)^2 + metadata,
    work = 128 * B * units + 64 * cs$core_lookups + 64 * metadata,
    record = list(route = "complete_support", live_entries_per_row = live, row_units = units,
      lookup_elements = cs$core_lookups, distinct_pairs = cs$pairs, block_support_incidence = cs$incidence))
}

evaluate_or_baseline_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE, plan = NULL) {
  validate_or_baseline_backend(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE); st <- reader$state
  valid_index <- function(x, n) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x)) && all(x == floor(x)) &&
    all(x >= 1 & x <= n) && !anyDuplicated(x)
  if (!valid_index(rows, st$N) || length(rows) > backend$control$max_chunk_rows ||
      !valid_index(links, length(backend$links)) || !identical(details, isTRUE(details)) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core"))
    or_baseline_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- or_baseline_if_budget(backend, length(rows), links, details)
  # CP08 stage 4: an eligible non-detailed chunk is one native batch; detailed
  # requests keep the route below.
  if (!details && !is.null(backend$complete_support))
    return(or_baseline_native_if_chunk(backend, reader, rows, links, scale, pre, plan))
  B <- length(rows); frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  kept <- backend$outcome_recipe$slope_index[backend$shared$keep]
  ref <- or_baseline_reference(backend)
  panel <- if (is.null(backend$complete_support)) NULL else
    or_baseline_complete_panel(reader, backend$complete_support)
  # A D60 chunk has no common weight: every support reads its recipe's W_g on
  # its own A rows (NA elsewhere, so any stray common-weight use is refused).
  ci <- st$cohort_index[rows]; w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
  A <- as.double(st$a_index[rows] > 0L)
  up <- matrix(0, B, length(pre$sids))
  hits <- vector("list", length(pre$sids))
  # Sparse hits are temporary and at most one row chunk, never retained by fit.
  for (j in seq_along(pre$sids)) {
    sid <- pre$sids[j]; s <- backend$support_specs[[sid]]; u <- backend$supports[[sid]]
    ix <- which(ci == s$cohort)
    hit <- if (is.null(panel)) cp04_pair_rows(reader, rows[ix], s$pair, TRUE) else
      or_baseline_complete_pair_rows(reader, panel, rows[ix], s$pair, TRUE)
    obs <- ix[hit$S == 1]
    if (is.null(ref)) {
      up[ix, j] <- w[ix] * (hit$S - A[ix] * u$rate) / u$massA
      hits[[j]] <- list(obs = obs, dy = hit$dy[hit$S == 1],
        z = frame$Z_a[st$a_index[rows[obs]], kept, drop = FALSE])
    } else {
      on_a <- A[ix] > 0; ia <- ix[on_a]
      wa <- chained_reference_weight(ref, s$reference, rows[ia])
      up[ia, j] <- wa * (hit$S[on_a] - u$rate) / u$massA
      hits[[j]] <- list(obs = obs, dy = hit$dy[hit$S == 1],
        z = chained_reference_design(ref, s$reference, rows[obs])[, kept, drop = FALSE],
        w = wa[match(obs, ia)])
    }
  }
  core <- or_pooled_if_core(backend, st$N, w, hits, up, pre, links, details)
  phi <- core$phi
  dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1), "link_id"))
  cp04_finite(phi, "complete baseline link IF")
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi, scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_rate_score_evaluations = length(pre$sids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = length(pre$beta_gids), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) {
    # Frame tau/eta summaries need a common fixed weight; a D60 context forms
    # them only where an existing consumer reads them (none here), so NULL.
    tau <- eta <- if (is.null(ref)) matrix(0, B, length(st$cohorts)) else NULL
    if (is.null(ref)) for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(support_ids = pre$sids, rates = up,
      beta_regression_ids = pre$beta_gids, beta = core$beta, beta_score = core$beta_score,
      beta_block_scores = core$beta_block_scores, regression_ids = pre$gids,
      alpha = core$alpha, alpha_score = core$alpha_score, direct = core$direct,
      denominator = core$denominator, numerator = core$numerator, design = core$design,
      objective = core$objective, tau = tau, eta = eta, lambda = A - st$lambda)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output)
    or_baseline_abort("IF output bound underestimated.", "bound")
  out
}

replay_or_baseline_if <- function(backend, chunk_rows, link_batch_size, callback,
                                  scale = "frame", details = FALSE) {
  validate_or_baseline_backend(backend)
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > backend$control$max_chunk_rows || !is.function(callback) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !scale %in% c("frame", "core"))
    or_baseline_abort("Invalid explicit replay controls/callback.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N; L <- length(backend$links); total <- 0
  native <- !details && !is.null(backend$complete_support)
  # Preflight ALL row/link work and each link batch's output before callbacks.
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1))
    pre <- or_baseline_if_budget(backend, min(B, N), ll, details)
    total <- ipw_owner_safe_sum(c(total, ceiling(N / B) * pre$work),
      backend$control$max_work_ops, "total baseline OR replay work")
    if (native) or_baseline_native_plan(backend, ll, min(B, N), pre)
  }
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1))
    plan <- if (native) or_baseline_native_plan(backend, ll, min(B, N)) else NULL
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1))
      callback(evaluate_or_baseline_if_chunk(backend, rr, ll, scale, details, plan))
      batches <- batches + 1L
    }
  }
  invisible(list(n_frame_rows = N, n_links = L, batches = batches, scale = scale,
    normalization_count = if (scale == "core") reader$state$n_A else N,
    work_ops_bound = total, all_F_identities = TRUE))
}
