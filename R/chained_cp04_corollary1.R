# CP04 cohort/actual-pair logit functional, public F/A/B fitting rosters.
# Known pair probabilities remain unsupported; custom callbacks with an explicit
# log-weight gradient run on R/chained_custom_selection.R (O15).
# All factors/partials use N (raw frame count); weights are already in scores.

cp04_abort <- function(message, class = "input", diagnostics = list()) {
  stop(structure(list(message = message, call = NULL, diagnostics = diagnostics),
    class = c(paste0("cdid_cp04_", class, "_error"), "cdid_ipw_fit_error",
              "cdid_error", "error", "condition")))
}

cp04_finite <- function(x, label, positive = FALSE) {
  if (any(!is.finite(x)) || (positive && any(x <= 0)))
    cp04_abort(paste0("Nonfinite", if (positive) " or nonpositive" else "",
                      " ", label, "."), "numeric")
  invisible(x)
}

# Same conservative compact-view convention: list/atomic slots plus eight
# header slots, attributes, and ceil(character bytes/2) per string. Borrowed
# owner arrays are counted once separately. This is a payload charge, not RSS.
#
# CP08 common-boundary performance: the walk runs in the registered native
# kernel `C_gmm_cp08_payload` (src/cdid.c) with the exact semantics of the
# retained recursive definition: environments and NULL charge nothing; every
# other node charges length(x) + 8; lists and pairlists recurse into every
# element (every alias occurrence counts again); character vectors add
# ceil(bytes / 2) entries and the bytes of each non-NA string; attributes(x) is
# charged once as a list whose own names are not scanned, and each attribute
# value receives the normal attribute scan. No R recursion graph is allocated.
# The native symbol address is resolved once per loaded DLL through the same
# verified registration path as the other CP08 kernels.
cp04_payload_native_cache <- new.env(parent = emptyenv())

cp04_payload_native_kernel <- function() {
  syms <- gmm_require_native_dll()
  owner <- environment(sys.function())
  cache <- get0("cp04_payload_native_cache", envir = owner, mode = "environment", inherits = FALSE)
  hit <- if (is.environment(cache)) get0("C_gmm_cp08_payload", envir = cache, inherits = FALSE) else NULL
  if (is.list(hit) && identical(hit$owner, owner) && gmm_sexp_aliases(hit$syms, syms)) return(hit$address)
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  si <- tryCatch(base::getNativeSymbolInfo("C_gmm_cp08_payload", PACKAGE = dlls[["cdid"]]), error = function(e) NULL)
  if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
      !identical(as.integer(si$numParameters), 3L))
    cdid_gmm_abort("Native symbol registration/forcing verification failed for `C_gmm_cp08_payload`.",
      c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
  if (is.environment(cache)) assign("C_gmm_cp08_payload", list(owner = owner, syms = syms, address = si$address), envir = cache)
  si$address
}

cp04_payload <- function(x, scan_attrs = TRUE) {
  .Call(cp04_payload_native_kernel(), x, scan_attrs, environment())
}

# CP08 general-panel performance: the sparse two-date pair-row join of
# `cp04_pair_rows` runs in the registered native kernel
# `C_gmm_cp08_sparse_pair_rows` (src/cdid.c), resolved once per loaded DLL
# through the same verified registration path as `C_gmm_cp08_payload`.
cp04_pair_rows_native_cache <- new.env(parent = emptyenv())

cp04_pair_rows_native_kernel <- function() {
  syms <- gmm_require_native_dll()
  owner <- environment(sys.function())
  cache <- get0("cp04_pair_rows_native_cache", envir = owner, mode = "environment", inherits = FALSE)
  hit <- if (is.environment(cache)) get0("C_gmm_cp08_sparse_pair_rows", envir = cache, inherits = FALSE) else NULL
  if (is.list(hit) && identical(hit$owner, owner) && gmm_sexp_aliases(hit$syms, syms)) return(hit$address)
  dlls <- if (isNamespace(owner)) getNamespaceInfo(owner, "DLLs") else base::getLoadedDLLs()
  si <- tryCatch(base::getNativeSymbolInfo("C_gmm_cp08_sparse_pair_rows", PACKAGE = dlls[["cdid"]]),
    error = function(e) NULL)
  if (is.null(si) || !inherits(si, "NativeSymbolInfo") || is.null(si$address) ||
      !identical(as.integer(si$numParameters), 6L))
    cdid_gmm_abort("Native symbol registration/forcing verification failed for `C_gmm_cp08_sparse_pair_rows`.",
      c("cdid_gmm_native_error", "cdid_gmm_accuracy_error"))
  if (is.environment(cache)) assign("C_gmm_cp08_sparse_pair_rows", list(owner = owner, syms = syms, address = si$address), envir = cache)
  si$address
}

# Structural admission of one pair-row query to the native join: plain integer
# query A-indices without NA, plain integer sealed postings (a long posting
# table carries double observation positions and is not admitted), integer
# non-NA date bounds in order and a scalar logical read flag. Everything else
# keeps the scalar R search below, decided before any native call; a native
# failure is never caught. Only the four bounds and the chunk's A-indices are
# inspected, never the sealed posting vectors' contents.
cp04_pair_rows_native_ok <- function(ai, post, bounds, read_outcomes) {
  is.integer(ai) && !is.object(ai) && length(ai) <= .Machine$integer.max && !anyNA(ai) && all(ai >= 0L) &&
    is.integer(post$post_a_index) && !is.object(post$post_a_index) &&
    is.integer(post$post_obs_index) && !is.object(post$post_obs_index) &&
    length(post$post_a_index) > 0L && length(post$post_a_index) <= .Machine$integer.max &&
    length(post$post_a_index) == length(post$post_obs_index) &&
    is.integer(bounds) && length(bounds) == 4L && !anyNA(bounds) &&
    all(bounds >= 1L & bounds <= length(post$post_a_index)) &&
    bounds[1L] <= bounds[2L] && bounds[3L] <= bounds[4L] &&
    is.logical(read_outcomes) && length(read_outcomes) == 1L && !is.na(read_outcomes)
}

cp04_budget <- function(ctrl, retained, scratch, work, output, label) {
  ipw_owner_safe_sum(retained, ctrl$max_retained_entries, paste(label, "retained payload"))
  ipw_owner_safe_sum(c(retained, scratch), ctrl$max_matrix_entries,
                     paste(label, "concurrent payload"))
  ipw_owner_safe_sum(work, ctrl$max_work_ops, paste(label, "work charge"))
  ipw_owner_safe_sum(output, ctrl$max_output_entries, paste(label, "output payload"))
  invisible(TRUE)
}

validate_cp04_reader <- function(reader) {
  if (!inherits(reader, "cdid_cp04_frame_reader") ||
      !ipw_owner_seal_ok(reader, "cp04_reader"))
    cp04_abort("Expected the immutable CP04 full-F reader.", "integrity")
  frame <- if (is.environment(reader$frame_ref) && exists("frame_owner", reader$frame_ref, inherits = FALSE))
    get("frame_owner", reader$frame_ref, inherits = FALSE) else NULL
  validate_ipw_frame_owner(frame)
  if (!identical(reader$panel_ref, frame$panel_ref))
    cp04_abort("CP04 reader panel reference disagrees with its frame owner.", "integrity")
  invisible(reader)
}

# Narrow borrowed immutable design/status adapter. id is the original supplied
# F order, not A order. frame_g/weights must be explicit, including frame-only
# rows. Existing frame_owner construction validates its own treatment design.
# `reference_carrier` (D60 reference route): after every identity check, the
# mechanical unit weights, the carrier design and every weight-derived cohort
# summary are blanked (NA). Each backend then reads designs, weights and rates
# at their recipe through the sealed reference reader; a stray static read
# fails as a nonfinite value instead of silently using the carrier.
construct_cp04_frame_reader <- function(frame_owner, id, Z, frame_g, weights, control,
                                        U = rep(TRUE, length(id)), rho = rep(1, length(id)),
                                        universe = "F", justification = NULL, allow_unknown_off_A = FALSE,
                                        has_intercept = NULL, reference_carrier = FALSE) {
  hdr <- cp04_reader_headers(frame_owner, id, Z, frame_g, weights, control, U, rho, allow_unknown_off_A, has_intercept)
  panel <- hdr$panel; ctrl <- hdr$ctrl; N <- hdr$N; Z <- hdr$Z
  for (start in seq.int(1L, N, by = min(N, ctrl$max_chunk_rows))) {
    rr <- seq.int(start, min(N, start + ctrl$max_chunk_rows - 1))
    if (anyNA(id[rr]) || (is.numeric(id) && any(!is.finite(id[rr]))) ||
        (is.character(id) && any(!nzchar(id[rr]))) || (!allow_unknown_off_A && anyNA(frame_g[rr])) ||
        any(frame_g[rr] == -Inf, na.rm = TRUE)) cp04_abort("Invalid full-F identity/cohort metadata.")
    cp04_finite(Z[rr[U[rr]], , drop = FALSE], "selection design on U")
    cp04_finite(weights[rr], "fixed full-F weights", TRUE)
  }
  id_key <- if (is.numeric(id)) { id_d <- as.double(id); id_d[id_d == 0.0] <- 0.0; id_d } else as.character(id)
  g_num <- as.double(frame_g); g_num[g_num == 0.0] <- 0.0
  if (anyDuplicated(id_key)) cp04_abort("Duplicate full-F identity.")
  unit <- ipw_owner_match_sorted_id(id_key, panel$unit_id_key, panel$id_type)
  present <- which(!is.na(unit))
  outside <- which(is.na(unit))
  if (length(present) != panel$n_total_units ||
      !identical(g_num[present], panel$units$g[unit[present]]) ||
      !identical(as.double(weights[present]), frame_owner$unit_weights_all[unit[present]]) ||
      !identical(outside, frame_owner$frame_only_units$supplied_pos) ||
      !identical(id[outside], frame_owner$frame_only_units$id) ||
      !identical(g_num[outside], frame_owner$frame_only_units$g) ||
      !identical(as.double(weights[outside]), frame_owner$frame_only_units$weight))
    cp04_abort("Full-F roster/order, cohort or weights disagree with actual owner inputs.", "identity")
  a_index <- integer(N)
  a_index[present] <- panel$unit_to_a[unit[present]]
  cp04_reader_assemble(hdr, frame_owner, id, Z, g_num, weights, U, rho, universe, justification,
    allow_unknown_off_A, reference_carrier, a_index)
}

# Header, declaration, dimension and budget checks of the full-F reader. They
# run before any chunk scan, lookup or allocation and are shared by the
# standalone constructor and the prepared-intake constructor below.
cp04_reader_headers <- function(frame_owner, id, Z, frame_g, weights, control, U, rho,
                                allow_unknown_off_A, has_intercept) {
  validate_ipw_frame_owner(frame_owner)
  panel <- get("panel_owner", frame_owner$panel_ref, inherits = FALSE)
  ctrl <- validate_ipw_owner_control(control)
  N <- length(id)
  if (!identical(allow_unknown_off_A, isTRUE(allow_unknown_off_A)))
    cp04_abort("Invalid unknown-cohort reader declaration.")
  if (!is.null(has_intercept)) {
    if (!is.logical(has_intercept) || is.object(has_intercept) || !is.null(dim(has_intercept)) ||
        length(has_intercept) != 1L || is.na(has_intercept))
      cp04_abort("Invalid selection intercept declaration.")
    if (is.matrix(Z) && !is.object(Z) && is.numeric(Z)) attr(Z, "has_intercept") <- has_intercept
  }
  if (!is.atomic(id) || is.object(id) || !is.null(dim(id)) || N < 1L ||
      !identical(is.numeric(id), identical(panel$id_type, "numeric")) ||
      (!is.numeric(id) && !is.character(id)) ||
      !is.matrix(Z) || is.object(Z) || !is.numeric(Z) || nrow(Z) != N || ncol(Z) < 0L ||
      !is.numeric(frame_g) || is.object(frame_g) || !is.null(dim(frame_g)) || length(frame_g) != N ||
      !is.numeric(weights) || is.object(weights) || !is.null(dim(weights)) || length(weights) != N)
    cp04_abort("Supply full-F id, selection design, cohorts and fixed weights with matching dimensions.")
  if (!is.logical(U) || length(U) != N || anyNA(U) ||
      !is.numeric(rho) || length(rho) != N || !any(U)) cp04_abort("Invalid fitting roster U or known rho.")
  cp04_finite(rho[U], "known normalized rho on U", TRUE)
  d <- ncol(Z)
  if (!is.null(attr(Z, "has_intercept"))) {
    hi <- attr(Z, "has_intercept")
    if (!is.logical(hi) || is.object(hi) || !is.null(dim(hi)) || length(hi) != 1L || is.na(hi) ||
        (isTRUE(hi) && d < 1L))
      cp04_abort("Invalid selection intercept declaration.")
  }
  ipw_owner_check_matrix_dims(N, d, ctrl$max_matrix_entries, "full-F selection design")
  for (item in list(list(x = id, nm = "full-F id"), list(x = Z, nm = "full-F design"),
                    list(x = frame_g, nm = "full-F cohort"), list(x = weights, nm = "full-F weights")))
    ipw_owner_refuse_decorative_attrs(item$x, item$nm)
  borrowed <- panel$accounting$actual_retained_entries +
    frame_owner$accounting$retained_matrix_entries + frame_owner$accounting$retained_vector_entries
  # All long alignment/grouping vectors belong to this one common reader.
  bytes <- frame_owner$id_byte_metadata$total
  extra_z_attrs <- attributes(Z)[setdiff(names(attributes(Z)), c("dim", "dimnames"))]
  supplied_metadata <- cp04_payload(if (length(extra_z_attrs)) list(dimnames(Z), attributes(id), extra_z_attrs) else list(dimnames(Z), attributes(id)))
  reader_bound <- N * (d + 14) + 8 * N + bytes + 4096 +
    unname(supplied_metadata["entries"])
  scratch <- 32 * (N + panel$n_total_units) + 8 * min(N, ctrl$max_chunk_rows) * (d + 1)
  cp04_budget(ctrl, c(borrowed, reader_bound), scratch,
    256 * (N + panel$n_total_units) * (d + 1 + ceiling(log2(N + 1))) +
      64 * (bytes + unname(supplied_metadata["char_bytes"])),
    reader_bound, "full-F reader")
  list(panel = panel, ctrl = ctrl, N = N, Z = Z, d = d, borrowed = borrowed, reader_bound = reader_bound, scratch = scratch)
}

# Prepared-intake reader (CP08 common boundary). The public request layer has
# already validated the complete F = A roster, its cohorts and unit weights on
# exactly these arrays; the sealed prepared record carries the sorted-owner
# position of every frame row. Headers, declarations and the budget are checked
# exactly as in the standalone constructor; the per-chunk identity/cohort/design
# scan, the sorted identity lookup and the owner agreement comparisons are
# replaced by the record's maps (the arrays are the same objects, checked by
# identity, never by a trusted flag). The retained representation is assembled
# and sealed by the same code as the standalone reader.
construct_cp04_frame_reader_prepared <- function(frame_owner, id, Z, frame_g, weights, control, prepared,
                                                 U = rep(TRUE, length(id)), rho = rep(1, length(id)),
                                                 universe = "F", justification = NULL, allow_unknown_off_A = FALSE,
                                                 reference_carrier = FALSE) {
  hdr <- cp04_reader_headers(frame_owner, id, Z, frame_g, weights, control, U, rho, allow_unknown_off_A, NULL)
  panel <- hdr$panel; N <- hdr$N
  if (!inherits(prepared, "cdid_request_prepared_intake") || !is.list(prepared) ||
      !ipw_owner_seal_ok(prepared, "prepared_intake") ||
      !gmm_sexp_aliases(panel$unit_id_key, prepared$unit_key) ||
      !gmm_sexp_aliases(id, prepared$frame_ids) || !gmm_sexp_aliases(weights, prepared$frame_w) ||
      is.null(prepared$cp04) || !gmm_sexp_aliases(Z, prepared$cp04$Z) ||
      !gmm_sexp_aliases(frame_g, prepared$cp04$frame_g) || !gmm_sexp_aliases(U, prepared$cp04$U) ||
      !gmm_sexp_aliases(rho, prepared$cp04$rho) || !identical(prepared$n_frame, N) ||
      !identical(frame_owner$n_frame_only_units, 0L) || !identical(frame_owner$n_a_units, panel$n_total_units) ||
      !isTRUE(reference_carrier == FALSE))
    cp04_abort("Prepared intake record does not describe the supplied full-F reader inputs.", "integrity")
  g_num <- as.double(frame_g); g_num[g_num == 0.0] <- 0.0
  unit <- prepared$frame_unit_sorted
  a_index <- panel$unit_to_a[unit]
  cp04_reader_assemble(hdr, frame_owner, id, Z, g_num, weights, U, rho, universe, justification,
    allow_unknown_off_A, reference_carrier, a_index)
}

# Cohort grouping, state, bound check and seal of the full-F reader (shared).
cp04_reader_assemble <- function(hdr, frame_owner, id, Z, g_num, weights, U, rho, universe, justification,
                                 allow_unknown_off_A, reference_carrier, a_index) {
  panel <- hdr$panel; ctrl <- hdr$ctrl; N <- hdr$N
  borrowed <- hdr$borrowed; reader_bound <- hdr$reader_bound; scratch <- hdr$scratch
  if (any(a_index > 0L & !U)) cp04_abort("U must contain all A structurally.", "support")
  if (anyNA(g_num[a_index > 0L])) cp04_abort("Cohorts must be known on all A rows.", "identity")
  cohorts <- sort(unique(g_num), method = "radix", na.last = if (allow_unknown_off_A) TRUE else NA)
  ci <- match(g_num, cohorts)
  ord <- order(ci, method = "radix")
  counts <- tabulate(ci, length(cohorts))
  ends <- cumsum(as.double(counts)); starts <- ends - counts + 1
  mass <- repeat_mass <- mass_U <- numeric(length(cohorts))
  for (j in seq_along(cohorts)) {
    rr <- ord[seq.int(starts[j], ends[j])]
    mass[j] <- sum(weights[rr]); mass_U[j] <- sum(weights[rr[U[rr]]]); repeat_mass[j] <- sum(weights[rr[a_index[rr] > 0L]])
  }
  cp04_finite(mass, "cohort mass", TRUE)
  cp04_finite(sum(weights), "total weight", TRUE)
  cp04_finite(sum(repeat_mass), "repeat weight", TRUE)
  n_A <- sum(a_index > 0L)
  if (n_A != frame_owner$n_a_units) cp04_abort("A mapping mismatch.", "identity")
  state <- list(id = id, Z = Z, g = g_num, w = as.double(weights),
    U = U, rho = rho, universe = universe, justification = justification, mass_U = mass_U / N,
    a_index = a_index, cohort_index = ci, cohort_order = ord,
    cohort_start = starts, cohort_end = ends, cohorts = cohorts,
    mass_F = mass / N, a = ifelse(mass_U > 0, repeat_mass / mass_U, 0),
    eta = repeat_mass / sum(repeat_mass), tau = mass / sum(mass),
    mean_w = sum(mass) / N, mean_wA = sum(repeat_mass) / N,
    N = N, n_A = n_A, lambda = n_A / N)
  if (isTRUE(reference_carrier)) {
    blank <- function(x) { x[] <- NA_real_; x }
    for (nm in c("w", "Z", "mass_U", "mass_F", "a", "eta", "tau", "mean_w", "mean_wA")) state[[nm]] <- blank(state[[nm]])
    state$reference_carrier <- TRUE
  }
  obj <- structure(list(state = state, frame_ref = frame_owner$ref,
    panel_ref = frame_owner$panel_ref, control = ctrl,
    accounting = list(owned_entries = 0, character_bytes = 0, borrowed_owner_entries = borrowed,
      constructor_scratch_bound = scratch, full_F_designs = 1L,
      full_F_masks_per_model = 0L, peak_rss_measured = FALSE)),
    class = "cdid_cp04_frame_reader")
  obj["ref"] <- list(NULL) # count the existing owner's seal pointer/name slots
  payload <- cp04_payload(obj); actual <- unname(payload["entries"])
  if (actual > reader_bound) cp04_abort("Retained reader bound underestimated.", "bound")
  cp04_budget(ctrl, c(borrowed, actual), scratch, 0, actual, "retained full-F reader")
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  ipw_owner_seal(obj, "cp04_reader")
}

# Uses the owner's sparse calendar postings. Sorted binary searches per chunk
# on date posting bounds avoid per-unit range queries, an N-mask, or a large
# min/max-A query when supplied F order is noncanonical.
#
# CP08 general-panel performance: on structurally admitted inputs
# (`cp04_pair_rows_native_ok`) the two subrange searches of every chunk row run
# in `C_gmm_cp08_sparse_pair_rows`, which returns the observed local rows in
# query order and their two posting observation positions; the R-level
# keep/ms/hs/mt/ht vector chain is not built. S and dy start at zero, the
# outcome endpoints are still read by `ipw_owner_read_y_positions` and
# subtracted as `y_right_obs - y_left_obs` here, and the finite check is
# unchanged. Other inputs take the retained scalar search; the decision is
# made before the native call and no native failure is caught.
cp04_pair_rows <- function(reader, rows, pair, read_outcomes = FALSE) {
  st <- reader$state
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  n <- length(rows)
  out <- list(S = numeric(n))
  if (read_outcomes) out$dy <- numeric(n)
  ai <- st$a_index[rows]
  keep <- which(ai > 0L)
  ds <- pair$s_calendar_index; dt <- pair$t_calendar_index
  if (is.null(ds)) ds <- match(as.double(pair$s), panel$all_times)
  if (is.null(dt)) dt <- match(as.double(pair$t), panel$all_times)
  post <- panel$postings
  if (length(keep) && !is.na(ds) && !is.na(dt) &&
      post$date_count[ds] > 0L && post$date_count[dt] > 0L) {
    bounds <- c(post$date_start[ds], post$date_end[ds], post$date_start[dt], post$date_end[dt])
    if (cp04_pair_rows_native_ok(ai, post, bounds, read_outcomes)) {
      hit <- .Call(cp04_pair_rows_native_kernel(), ai, post$post_a_index, post$post_obs_index,
        bounds, read_outcomes, reader$control$max_matrix_entries)
      obs <- hit$obs
      if (length(obs)) {
        out$S[obs] <- 1
        if (read_outcomes) {
          yy <- ipw_owner_read_y_positions(panel, hit$pos_s, hit$pos_t)
          out$dy[obs] <- yy$y_right_obs - yy$y_left_obs
        }
      }
    } else {
      if (length(keep)) {
        pa <- post$post_a_index
        ms <- ipw_owner_match_sorted_double(ai[keep], pa, post$date_start[ds], post$date_end[ds])
        hs <- which(!is.na(ms))
        if (length(hs)) {
          mt <- ipw_owner_match_sorted_double(ai[keep[hs]], pa, post$date_start[dt], post$date_end[dt])
          ht <- which(!is.na(mt))
          if (length(ht)) {
            obs <- keep[hs[ht]]
            out$S[obs] <- 1
            if (read_outcomes) {
              ps <- post$post_obs_index[ms[hs[ht]]]
              pt <- post$post_obs_index[mt[ht]]
              yy <- ipw_owner_read_y_positions(panel, ps, pt)
              out$dy[obs] <- yy$y_right_obs - yy$y_left_obs
            }
          }
        }
      }
    }
  }
  if (read_outcomes) cp04_finite(out$dy, "observed pair differences")
  out
}

cp04_solve <- function(R, b) backsolve(R, forwardsolve(t(R), b))

# General binary-response adapter to existing stable glm/refusal machinery.
# Selection supplies ALL cohort F rows, including response-zero nonobservers;
# this deliberately does not call the observer-only rate-weighted worker.
cp04_logit <- function(X, y, w, glm_control, N, identity) {
  p <- ncol(X)
  fail <- function(message, type, extra = list()) cp04_abort(message, type,
    c(list(model = identity, n_fit = nrow(X), p = p), extra))
  if (nrow(X) < p) fail("Fitting rows smaller than logit dimension.", "rank")
  if (!is.numeric(y) || length(y) != nrow(X) || anyNA(y) || any(!y %in% c(0, 1)))
    fail("Only binary pair/treatment logit responses are supported.", "unsupported")
  cp04_finite(X, "logit design"); cp04_finite(w, "logit fitting weights", TRUE)
  if (ipw_logit_exact_alias(X)) fail("Logit has deficient model rank.", "rank")
  gf <- ipw_glm_fit_internal(X, y, w, glm_control)
  f <- gf$fit
  diag <- list(converged = f$converged, boundary = f$boundary, rank = f$rank,
    iter = f$iter, warnings = gf$warnings)
  if (is.null(f$rank) || f$rank < p || anyNA(f$coefficients) || is.null(f$qr$qr))
    fail("Logit has deficient model rank.", "rank", diag)
  qr_diag <- abs(diag(f$qr$qr[seq_len(p), seq_len(p), drop = FALSE]))
  if (any(!is.finite(qr_diag) | qr_diag <= 0)) fail("Singular/nonfinite logit QR.", "rank", diag)
  # Give the existing exact finite-root obstruction its specific diagnostic,
  # including observed all-success/all-failure cohorts. Never reinterpret it
  # as a supplied known design probability.
  sep <- ipw_logit_coordinate_separation(X, y, w)
  if (!is.null(sep)) fail("Logit has no finite root: signed-coordinate separation witness.",
    "separation", c(diag, list(separation_witness = sep)))
  if (isTRUE(f$boundary) || any(grepl("fitted probabilities numerically 0 or 1 occurred", gf$warnings, fixed = TRUE)))
    fail("Logit encountered boundary probabilities.", "boundary", diag)
  if (!isTRUE(f$converged) || any(grepl("did not converge", gf$warnings, fixed = TRUE)))
    fail("Logit did not converge under supplied controls.", "convergence", diag)
  beta <- as.double(f$coefficients)
  cp04_finite(beta, "logit coefficients")
  lv <- ipw_logit_values(X, beta)
  cp04_finite(lv$eta, "logit predictors"); cp04_finite(lv$odds, "logit odds", TRUE)
  cp04_finite(lv$prob, "logit probabilities")
  eps <- 10 * .Machine$double.eps # unchanged existing worker boundary rule
  if (any(lv$prob <= eps | lv$prob >= 1 - eps | lv$prob * (1 - lv$prob) <= 0))
    fail("Logit probabilities reached numerical boundary.", "boundary", diag)
  H <- crossprod(X, w * lv$prob * (1 - lv$prob) * X) / N
  cp04_finite(H, "logit information")
  R <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(R) || any(!is.finite(R)) || any(diag(R) <= 0))
    fail("Logit information is singular/nonfinite.", "rank", diag)
  # f, QR/IRLS and all observation-length arrays die with this fitting call.
  list(coefficients = beta, H_chol = R, diagnostics = diag)
}

# `reference` (D60): the cell's recipe selection design, W_g and a_k^{g,U}.
cp04_selection_fit <- function(reader, spec, glm_control, ctrl, reference = NULL) {
  st <- reader$state; k <- spec$cohort
  rows <- st$cohort_order[seq.int(st$cohort_start[k], st$cohort_end[k])]
  rows <- rows[st$U[rows]]
  X <- if (is.null(reference)) st$Z[rows, , drop = FALSE] else chained_reference_s(reference, spec$reference, rows)
  y <- numeric(length(rows))
  for (start in seq.int(1L, length(rows), by = ctrl$max_chunk_rows)) {
    rr <- seq.int(start, min(length(rows), start + ctrl$max_chunk_rows - 1))
    y[rr] <- cp04_pair_rows(reader, rows[rr], spec$pair)$S
  }
  a_k <- if (is.null(reference)) st$a[k] else chained_reference_stat(reference, spec$reference, k, "a")
  if (!is.finite(a_k) || a_k <= 0)
    cp04_abort("A relevant cohort has no repeat-observation mass.", "support", spec)
  w <- if (is.null(reference)) st$w[rows] else chained_reference_weight(reference, spec$reference, rows)
  mass_U <- if (is.null(reference)) st$mass_U[k] else chained_reference_stat(reference, spec$reference, k, "mass_U")
  fit <- cp04_logit(X, y, w, glm_control, st$N, spec)
  fit$pi <- sum(w * y) / (mass_U * st$N)
  fit$r <- fit$pi / a_k # diagnostics only; never another functional argument
  fit
}

cp04_predict_selection <- function(reader, rows, spec, selection, reference = NULL) {
  X <- if (is.null(reference)) reader$state$Z[rows, , drop = FALSE] else chained_reference_s(reference, spec$reference, rows)
  p <- ipw_logit_values(X, selection$coefficients)$prob
  cp04_finite(p, "replayed selection probabilities", TRUE)
  if (any(p >= 1)) cp04_abort("Replayed selection probability at boundary.", "boundary", spec)
  p
}

cp04_treatment_fit <- function(reader, spec, selection_specs, selections, glm_control, ctrl) {
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  p <- frame$n_covariates; d <- ncol(st$Z); m <- length(spec$cohorts)
  # One bounded observed-risk fit buffer. No outcome/model arrays survive return.
  n <- spec$n_observers
  X <- matrix(0, n, p); Z <- matrix(0, n, d)
  y <- w <- pp <- dy <- numeric(n); kk <- integer(n); pos <- 0L
  for (start in seq.int(1L, st$N, by = ctrl$max_chunk_rows)) {
    rows <- seq.int(start, min(st$N, start + ctrl$max_chunk_rows - 1))
    ci <- match(st$cohort_index[rows], spec$cohorts)
    rows <- rows[!is.na(ci)]; ci <- ci[!is.na(ci)]
    if (!length(rows)) next
    hit <- cp04_pair_rows(reader, rows, spec$pair, TRUE)
    keep <- which(hit$S == 1)
    if (!length(keep)) next
    rows <- rows[keep]; ci <- ci[keep]
    ii <- seq.int(pos + 1L, pos + length(rows)); pos <- pos + length(rows)
    X[ii, ] <- frame$Z_a[st$a_index[rows], , drop = FALSE]
    Z[ii, ] <- st$Z[rows, , drop = FALSE]
    y[ii] <- as.double(ci == 1L); kk[ii] <- ci; dy[ii] <- hit$dy[keep]
    for (k in unique(ci)) {
      jj <- which(ci == k); sj <- spec$selection_ids[k]
      pp[ii[jj]] <- cp04_predict_selection(reader, rows[jj], selection_specs[[sj]], selections[[sj]])
    }
    w[ii] <- st$w[rows] * st$a[st$cohort_index[rows]] / (st$rho[rows] * ipw_clip_floor(pp[ii]))
  }
  if (pos != n) cp04_abort("Observed risk support changed from its owner view.", "identity")
  fit <- cp04_logit(X, y, w, glm_control, st$N, spec)
  lv <- ipw_logit_values(X, fit$coefficients)
  wp <- w * y; wm <- w * (1 - y) * lv$odds
  dp <- sum(wp) / st$N; dm <- sum(wm) / st$N
  cp04_finite(c(dp, dm), "treated/control normalizers", TRUE)
  T <- sum(wp * dy) / sum(wp); M <- sum(wm * dy) / sum(wm)
  cp04_finite(c(T, M), "treated/control means")
  # Upstream order is (a_k, alpha_k[1:d]) for each actual model cohort.
  u <- m * (d + 1L)
  C <- matrix(0, p, u)
  U <- matrix(0, 4L, u, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  score <- w * (y - lv$prob) * X
  residual <- cbind(wp * (dy - T) / dp, wm * (dy - M) / dm, wp, wm)
  for (k in seq_len(m)) {
    rr <- which(kk == k); ix <- (k - 1L) * (d + 1L) + seq_len(d + 1L)
    log_derivative <- cbind(1 / st$a[spec$cohorts[k]],
      -((1 - pp[rr]) * ipw_clip_floor_active(pp[rr])) * Z[rr, , drop = FALSE])
    C[, ix] <- crossprod(score[rr, , drop = FALSE], log_derivative) / st$N
    U[, ix] <- crossprod(residual[rr, , drop = FALSE], log_derivative) / st$N
  }
  V <- matrix(0, 4L, p, dimnames = list(c("T", "M", "dp", "dm"), NULL))
  wm_act <- wm * lv$active
  V[2L, ] <- colSums(wm_act * (dy - M) * X) / (st$N * dm)
  V[4L, ] <- colSums(wm_act * X) / st$N
  cp04_finite(c(C, U, V), "analytic cross derivatives")
  fit$cross <- C
  fit$partials <- list(upstream = U, beta = V)
  fit$means <- c(T = T, M = M, dp = dp, dm = dm)
  fit$clipping <- ipw_clip_summary(prob = lv$prob, q = pp, controls = y == 0, weights = w)
  fit
}

validate_cp04_backend <- function(backend) {
  if (!inherits(backend, "cdid_cp04_corollary1") || !ipw_owner_seal_ok(backend, "cp04_backend"))
    cp04_abort("Expected the immutable private CP04 backend.", "integrity")
  reader <- if (is.environment(backend$reader_ref) && exists("cp04_reader", backend$reader_ref, inherits = FALSE))
    get("cp04_reader", backend$reader_ref, inherits = FALSE) else NULL
  validate_cp04_reader(reader)
  if (!is.null(backend$reference_ref)) chained_backend_reference(backend)
  invisible(backend)
}

# A D60 reference request always uses the general selection-model consumer
# (recipe-keyed cells; pooling inside each recipe).
fit_cp04_corollary1 <- function(reader, link_views, control, glm_control,
                                retained_external_entries = 0, external_work_ops = 0,
                                selection_model = "pair_logit", selection_pooling = "separate",
                                reference = NULL) {
  validate_cp04_reader(reader)
  if (!is.null(reference) || !selection_use_legacy_consumer(reader, selection_model, selection_pooling)) {
    return(fit_cp04_selection_models(
      reader = reader, link_views = link_views, control = control, glm_control = glm_control,
      selection_model = selection_model, selection_pooling = selection_pooling,
      retained_external_entries = retained_external_entries, external_work_ops = external_work_ops,
      reference = reference
    ))
  }
  ctrl <- validate_ipw_owner_control(control); gc <- validate_ipw_glm_control(glm_control)
  st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (!is.list(link_views) || !length(link_views)) cp04_abort("Supply actual private link views.")
  # Metadata interned before any numerical fit. All models share THIS actual
  # reader/design/weights/solver/law; no label/count/coef-based reuse is possible.
  L <- length(link_views); d <- ncol(st$Z); p <- frame$n_covariates
  if (d < 1L) cp04_abort("Supply full-F id, selection design, cohorts and fixed weights with matching dimensions.")
  if (!is.numeric(retained_external_entries) || length(retained_external_entries) != 1L ||
      !is.finite(retained_external_entries) || retained_external_entries < 0 ||
      !is.numeric(external_work_ops) || length(external_work_ops) != 1L ||
      !is.finite(external_work_ops) || external_work_ops < 0) cp04_abort("Invalid external pipeline reservation.")
  base <- reader$accounting$owned_entries + reader$accounting$borrowed_owner_entries + retained_external_entries
  header_bound <- 1024 + 512 * L + sum(vapply(link_views, function(v)
    if (is.list(v) && is.list(v$cohorts)) 64 * length(v$cohorts$cohorts) else 0, numeric(1)))
  cp04_budget(ctrl, c(base, header_bound), header_bound, 64 * L * L,
    header_bound, "model metadata")
  ss <- ts <- links <- list()
  max_fit_n <- 0L
  for (ell in seq_along(link_views)) {
    v <- link_views[[ell]]
    validate_ipw_link_view(v, frame_owner = frame)
    source_key <- v$functional_key
    if (!gmm_scalar_character(source_key) || !startsWith(source_key, "cmpcat:"))
      cp04_abort("CP04 requires the supplied comparison-catalog functional key.", "unsupported")
    key_bytes <- nchar(source_key, type = "bytes")
    cp04_budget(ctrl, c(base, header_bound), header_bound + 16 * key_bytes,
      64 * (L * L + key_bytes), header_bound + 16 * key_bytes, "corrected operator key")
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate))
      cp04_abort("Legacy nuisance candidates are unsupported for corrected selection.", "unsupported")
    if (!is.null(v$nuisance_ref)) {
      nv <- get("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair))
        cp04_abort("Cross-pair or unverified legacy nuisance view sharing is unsupported for CP04.", "unsupported")
    }
    ks <- sort(v$cohorts$K, method = "radix")
    cohorts <- match(c(v$cohorts$g, ks), st$cohorts)
    if (anyNA(cohorts) || any(v$support$by_cohort_pair_count < 1L) || v$support$n_observed_risk_units < 1L)
      cp04_abort("Unsupported full-F cohort or pair support.", "support")
    si <- integer(length(cohorts))
    pair <- v$pair
    for (k in seq_along(cohorts)) {
      spec <- list(pair = pair, cohort = cohorts[k])
      same <- which(vapply(ss, identical, logical(1), spec))
      if (length(same)) si[k] <- same[1L] else {
        ss[[length(ss) + 1L]] <- spec; si[k] <- length(ss)
      }
      max_fit_n <- max(max_fit_n, st$cohort_end[cohorts[k]] - st$cohort_start[cohorts[k]] + 1)
    }
    # K order is normalized by actual values; operator labels are provenance.
    spec <- list(pair = pair, cohorts = cohorts, selection_ids = si,
      n_observers = as.integer(v$support$n_observed_risk_units))
    same <- which(vapply(ts, identical, logical(1), spec))
    if (length(same)) ti <- same[1L] else { ts[[length(ts) + 1L]] <- spec; ti <- length(ts) }
    max_fit_n <- max(max_fit_n, spec$n_observers)
    # Bind the corrected law in the existing conventions position, preserving
    # source g-index/pool identity and original/horizon operator provenance.
    links[[ell]] <- list(model = ti, link_id = v$link_id,
      functional_key = paste0("cmpcat:CP04-", st$universe, "-weighted-cohort-pair-logit-a/p-etaQG|",
        substring(source_key, 8L)), source_functional_key = source_key)
  }
  ids <- vapply(links, function(x) if (is.null(x$link_id)) "" else x$link_id, character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) cp04_abort("Distinct nonempty link IDs required.")
  metadata_payload <- cp04_payload(list(ss, ts, links))
  metadata <- unname(metadata_payload["entries"])
  model_bound <- length(ss) * (d * d + d + 2048) + sum(vapply(ts, function(s) {
    u <- length(s$cohorts) * (d + 1)
    p * p + p + p * u + 4 * (u + p) + 2048
  }, numeric(1)))
  # glm.fit's one model buffer, QR/IRLS and derivative temporaries. No live
  # previous model's n_fit arrays are retained by returned compact helpers.
  # The exact-alias scan uses only O(max_fit_n) column/comparison vectors;
  # its O(n_fit * dimension^2) work fits the conservative +2 fit allowance below.
  scratch <- 128 * max_fit_n * (p + d + 16) + 64 * (p + d + 1)^2 + header_bound
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  sel_work <- sum(vapply(ss, function(s) {
    n_s <- st$cohort_end[s$cohort] - st$cohort_start[s$cohort] + 1
    32 * n_s * (d + 1 + log_a) + 32 * (gc$maxit + 2) * n_s * (d + 1)^2 + 32 * (d + 1)^3
  }, numeric(1)))
  trt_work <- sum(vapply(ts, function(s) {
    n_r <- sum(st$cohort_end[s$cohorts] - st$cohort_start[s$cohorts] + 1)
    n_o <- s$n_observers
    u <- length(s$cohorts) * (d + 1)
    16 * st$N + 32 * n_r * (1 + log_a) + 32 * n_o * (p + d + 2) +
      32 * (gc$maxit + 2) * n_o * (p + 1)^2 + 32 * (p + 1)^3 +
      32 * n_o * (p + 4) * (d + 2) + 32 * (p + 4) * (u + p)
  }, numeric(1)))
  work <- sel_work + trt_work + 64 * (L * L + unname(metadata_payload["char_bytes"]))
  retained_bound <- base + metadata + model_bound + 2048
  cp04_budget(ctrl, retained_bound, scratch, c(work, external_work_ops), metadata + model_bound + 2048, "CP04 fits")
  selections <- lapply(ss, function(s) cp04_selection_fit(reader, s, gc, ctrl))
  treatments <- lapply(ts, function(s) cp04_treatment_fit(reader, s, ss, selections, gc, ctrl))
  state <- list(selection_specs = ss, selections = selections, treatment_specs = ts,
    treatments = treatments, links = links, glm_control = gc,
    law = paste0("CP04:", st$universe, ":weighted-cohort-pair-logit:a_U/(rho*p_U):eta=Q(G)"),
    universe = st$universe, transport_justification = st$justification)
  obj <- structure(c(state, list(reader_ref = reader$ref, control = ctrl,
    accounting = list(owned_entries = 0, character_bytes = 0, borrowed_entries = base,
      retained_bound = retained_bound, fit_scratch_bound = scratch, work_ops_bound = work,
      n_selection_fits = length(ss), n_treatment_fits = length(ts),
      full_F_masks_per_model = 0L, full_F_by_link_IF_retained = FALSE,
      inverse_matrices_retained = FALSE, peak_rss_measured = FALSE))),
    class = "cdid_cp04_corollary1")
  obj["ref"] <- list(NULL) # ipw_owner_seal replaces this reserved binding
  payload <- cp04_payload(obj); actual <- unname(payload["entries"])
  if (actual > metadata + model_bound + 2048) cp04_abort("Retained model bound underestimated.", "bound")
  obj$accounting$owned_entries <- actual
  obj$accounting$character_bytes <- unname(payload["char_bytes"])
  ipw_owner_seal(obj, "cp04_backend")
}

# Shared preflight arithmetic for one requested row/link batch. Replay invokes
# it for every link batch BEFORE its first callback, including output limits.
cp04_if_budget <- function(backend, B, links) {
  if (isTRUE(backend$uses_selection_models))
    return(cp04_selection_if_budget(backend, B, links, details = TRUE))
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; ctrl <- backend$control
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  tids <- unique(vapply(backend$links[links], `[[`, integer(1), "model"))
  sids <- unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "selection_ids")))
  d <- ncol(st$Z); p <- frame$n_covariates
  widths <- length(sids) * (d + 2) + sum(vapply(backend$treatment_specs[tids],
    function(s) 2 * p + 10 + length(s$cohorts) * (d + 1), numeric(1))) +
    3 * length(st$cohorts) + length(links) + 16
  max_bytes <- frame$id_byte_metadata$max
  metadata <- 4096 + 1024 * (length(backend$selections) + length(backend$treatments) +
    length(st$cohorts) + length(links)) + B * max_bytes +
    4 * sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  output <- B * (widths + 8) + metadata
  scratch <- 32 * B * (widths + p + d + 4) + 32 * (p + d + 1)^2 + metadata
  retained <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  log_a <- 1 + ceiling(log2(st$n_A + 1))
  per_row_work <- 32 * (length(sids) + length(tids)) * log_a +
    32 * length(sids) * (d + 1)^2 +
    32 * sum(vapply(backend$treatment_specs[tids], function(s)
      (p + 4) * (length(s$cohorts) * (d + 1) + p + 2), numeric(1))) +
    32 * (3 * length(st$cohorts) + length(links) + 4)
  work <- B * per_row_work + 64 * metadata
  cp04_budget(ctrl, retained, scratch, work, output, "CP04 IF chunk")
  list(output = output, scratch = scratch, work = work, tids = tids, sids = sids)
}

# Explicit row and link batches; output is never forced to N x L. Diagnostic
# coordinates permit independent finite-refit checks of the triangular stack.
evaluate_cp04_if_chunk <- function(backend, rows, links, scale = "frame", details = FALSE) {
  validate_cp04_backend(backend)
  if (isTRUE(backend$uses_selection_models))
    return(evaluate_cp04_selection_if_chunk(backend, rows, links, scale = scale, details = details))
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; ctrl <- backend$control
  if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) || anyNA(rows) || any(rows != floor(rows)) ||
      any(rows < 1 | rows > st$N) || anyDuplicated(rows) || length(rows) > ctrl$max_chunk_rows ||
      !is.numeric(links) || is.object(links) || !is.null(dim(links)) || !length(links) || anyNA(links) || any(links != floor(links)) ||
      any(links < 1 | links > length(backend$links)) || anyDuplicated(links) ||
      !identical(details, isTRUE(details)) || !is.character(scale) || length(scale) != 1L || is.na(scale) ||
      !scale %in% c("frame", "core")) cp04_abort("Invalid explicit IF row/link batch, scale or details flag.")
  pre <- cp04_if_budget(backend, length(rows), links)
  tids <- pre$tids; sids <- pre$sids
  B <- length(rows); d <- ncol(st$Z)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE); p <- frame$n_covariates
  # Scores evaluated once for each actual model needed by this requested batch.
  ci <- st$cohort_index[rows]; w <- st$w[rows]; A <- as.double(st$a_index[rows] > 0L)
  a_if <- matrix(0, B, length(st$cohorts))
  for (k in unique(ci)) {
    rr <- which(ci == k & st$U[rows])
    a_if[rr, k] <- w[rr] * (A[rr] - st$a[k]) / st$mass_U[k]
  }
  alpha <- probabilities <- responses <- vector("list", length(backend$selections))
  for (sj in sids) {
    s <- backend$selection_specs[[sj]]; f <- backend$selections[[sj]]
    rr <- which(ci == s$cohort & st$U[rows])
    xi <- matrix(0, B, d); pr <- S <- numeric(B)
    if (length(rr)) {
      pr[rr] <- cp04_predict_selection(reader, rows[rr], s, f)
      S[rr] <- cp04_pair_rows(reader, rows[rr], s$pair)$S
      score <- w[rr] * (S[rr] - pr[rr]) * st$Z[rows[rr], , drop = FALSE]
      xi[rr, ] <- t(cp04_solve(f$H_chol, t(score)))
    }
    alpha[[sj]] <- xi; probabilities[[sj]] <- pr; responses[[sj]] <- S
  }
  models <- vector("list", length(backend$treatments))
  for (tj in tids) {
    spec <- backend$treatment_specs[[tj]]; f <- backend$treatments[[tj]]
    u <- length(spec$cohorts) * (d + 1L)
    upstream <- matrix(0, B, u)
    S <- pp <- numeric(B)
    for (k in seq_along(spec$cohorts)) {
      sj <- spec$selection_ids[k]; rr <- which(ci == spec$cohorts[k])
      ix <- (k - 1L) * (d + 1L) + seq_len(d + 1L)
      upstream[, ix] <- cbind(a_if[, spec$cohorts[k]], alpha[[sj]])
      S[rr] <- responses[[sj]][rr]; pp[rr] <- probabilities[[sj]][rr]
    }
    rr <- which(S == 1)
    beta_score <- matrix(0, B, p); direct <- matrix(0, B, 4L)
    direct[, 3L] <- -f$means["dp"]; direct[, 4L] <- -f$means["dm"]
    if (length(rr)) {
      X <- frame$Z_a[st$a_index[rows[rr]], , drop = FALSE]
      hit <- cp04_pair_rows(reader, rows[rr], spec$pair, TRUE)
      lv <- ipw_logit_values(X, f$coefficients)
      cp04_finite(c(lv$prob, lv$eta, lv$odds), "treatment replay")
      D <- as.double(ci[rr] == spec$cohorts[1L])
      fw <- w[rr] * st$a[ci[rr]] / (st$rho[rows[rr]] * ipw_clip_floor(pp[rr]))
      wp <- fw * D; wm <- fw * (1 - D) * lv$odds
      beta_score[rr, ] <- fw * (D - lv$prob) * X
      direct[rr, 1L] <- wp * (hit$dy - f$means["T"]) / f$means["dp"]
      direct[rr, 2L] <- wm * (hit$dy - f$means["M"]) / f$means["dm"]
      direct[rr, 3L] <- direct[rr, 3L] + wp
      direct[rr, 4L] <- direct[rr, 4L] + wm
    }
    xb <- t(cp04_solve(f$H_chol, t(beta_score + upstream %*% t(f$cross))))
    coord <- direct + upstream %*% t(f$partials$upstream) + xb %*% t(f$partials$beta)
    colnames(coord) <- c("T", "M", "dp", "dm")
    cp04_finite(c(xb, coord), "complete frame IF")
    models[[tj]] <- if (details) list(beta = xb, coordinates = coord, direct = direct,
      beta_score = beta_score, upstream = upstream,
      normalizers_A = coord[, 3:4, drop = FALSE] / st$lambda -
        outer(A - st$lambda, f$means[c("dp", "dm")] / st$lambda^2)) else
      list(phi = coord[, 1L] - coord[, 2L])
  }
  phi <- matrix(0, B, length(links), dimnames = list(NULL,
    vapply(backend$links[links], `[[`, character(1), "link_id")))
  for (j in seq_along(links)) {
    m <- models[[backend$links[[links[j]]]$model]]
    phi[, j] <- if (details) m$coordinates[, 1L] - m$coordinates[, 2L] else m$phi
  }
  if (scale == "core") phi <- st$lambda * phi
  out <- list(rows = rows, id = st$id[rows], links = links, phi = phi,
    scale = scale, normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(n_selection_score_evaluations = length(sids),
      n_treatment_score_evaluations = length(tids), output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch, work_ops_bound = pre$work))
  if (details) {
    # All diagnostic coordinates are on the frame scale, even with core phi.
    tau <- eta <- matrix(0, B, length(st$cohorts))
    for (k in seq_along(st$cohorts)) {
      tau[, k] <- w * (as.double(ci == k) - st$tau[k]) / st$mean_w
      eta[, k] <- w * A * (as.double(ci == k) - st$eta[k]) / st$mean_wA
    }
    out$details <- list(a = a_if, alpha = alpha, models = models,
      tau = tau, eta = eta, lambda = A - st$lambda)
  }
  actual_output <- unname(cp04_payload(out)["entries"])
  if (actual_output > pre$output) cp04_abort("IF output bound underestimated.", "bound")
  out
}

replay_cp04_if <- function(backend, chunk_rows, link_batch_size, callback,
                           scale = "frame", details = FALSE) {
  validate_cp04_backend(backend)
  ctrl <- backend$control
  B <- ipw_owner_int_bound(chunk_rows, "chunk_rows")
  LB <- ipw_owner_int_bound(link_batch_size, "link_batch_size")
  if (B > ctrl$max_chunk_rows || !is.function(callback) ||
      !is.character(scale) || length(scale) != 1L || is.na(scale) || !scale %in% c("frame", "core") ||
      !identical(details, isTRUE(details))) cp04_abort("Invalid replay controls/callback.")
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  N <- reader$state$N; L <- length(backend$links)
  # Preflight cumulative work BEFORE callbacks (not just each individual chunk).
  total_work <- 0
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    pre <- cp04_if_budget(backend, min(B, N), ll)
    total_work <- ipw_owner_safe_sum(c(total_work, ceiling(N / B) * pre$work),
      ctrl$max_work_ops, "total CP04 replay work")
  }
  ipw_owner_safe_sum(total_work, ctrl$max_work_ops, "total CP04 replay work")
  batches <- 0L
  for (ls in seq.int(1L, L, by = LB)) {
    ll <- seq.int(ls, min(L, ls + LB - 1L))
    for (start in seq.int(1L, N, by = B)) {
      rr <- seq.int(start, min(N, start + B - 1L))
      batch <- evaluate_cp04_if_chunk(backend, rr, ll, scale, details)
      callback(batch)
      rm(batch)
      batches <- batches + 1L
    }
  }
  invisible(list(n_frame_rows = N, normalization_count = if (scale == "core") reader$state$n_A else N,
    n_links = L, batches = batches, scale = scale, work_ops_bound = total_work,
    all_F_identities = TRUE))
}
