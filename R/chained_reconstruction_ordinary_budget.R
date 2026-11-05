# Ordinary static generic reconstruction returns only the link phi and compact
# diagnostics. Its cores still construct every rate, OR and treatment channel.
# Each pair below is (constant entries, entries per row). Phase maxima are
# taken COEFFICIENTWISE, so the B=1/2 planner remains an upper affine model.
reconstruction_ordinary_inventory <- function(backend, B, links, pre, selection_pre = NULL) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  augmented <- inherits(backend, "cdid_augmented_backend")
  X <- if (augmented) augmented_backend_X_F(backend) else imputation_backend_X_F(backend, reader)
  add <- ipw_owner_safe_sum; prod <- ipw_owner_safe_product
  A <- function(constant = 0, row = 0) c(constant = constant, row = row)
  plus <- function(...) {
    x <- list(...)
    A(add(vapply(x, `[[`, numeric(1), "constant")), add(vapply(x, `[[`, numeric(1), "row")))
  }
  value <- function(x) add(c(x[["constant"]], prod(B, x[["row"]])))
  nc <- length(backend$cells); ng <- length(backend$groups)
  nl <- length(links); nt <- length(pre$tids); G <- length(st$cohorts)
  pt <- backend$treatment_recipe$p; p <- if (is.null(backend$shared_or)) 0 else backend$shared_or$p
  pg <- vapply(backend$groups, `[[`, integer(1), "p")
  peak_p <- max(c(0, pg)); group_columns <- add(pg)
  cell_columns <- add(vapply(backend$cells, function(cf) pg[cf$group_id], numeric(1)))
  name_bytes <- add(nchar(c(backend$design_columns, backend$outcome_recipe$columns,
    backend$treatment_recipe$columns, colnames(st$Z)), type = "bytes"))
  link_bytes <- add(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  index_attributes <- unname(cp04_payload(attributes(links))["entries"])
  # Existing named designs keep their old cache decline. Fund their complete
  # row-name vector in each possible current/replacement design, without a
  # B-dependent scan or choosing a different cache plan.
  row_names <- unname(cp04_payload(rownames(X))["entries"])
  output <- A(add(c(8192, 256 * nl, 8 * link_bytes, index_attributes)),
    nl + 4 + frame$id_byte_metadata$max)

  # Full-needed-record upper bound, never the selected cache's B-dependent
  # entries. It covers raw/cache/replacement records and their row/column
  # labels. Witnesses, position/seen/membership maps and final/armed views
  # have separate bindings/headers below; source owners remain borrowed.
  cache_constant <- 1024; cache_row <- 0
  for (cid in pre$cids) {
    gf <- backend$groups[[backend$cells[[cid]]$group_id]]
    bytes <- add(nchar(gf$kept_columns, type = "bytes"))
    cache_constant <- add(c(cache_constant, 256, 8 * bytes))
    cache_row <- add(c(cache_row, 4 * (gf$p + 2)))
  }
  cache <- A(cache_constant, cache_row)
  incidence <- add(vapply(backend$groups, function(g) length(g$cell_ids), numeric(1)))
  if (!is.null(backend$shared_or)) incidence <- add(c(incidence,
    vapply(backend$shared_or$blocks, function(b) length(b$control_cell_ids), numeric(1))))
  maps <- A(add(c(8192, 64 * (nc + ng + length(backend$or_blocks) +
    length(backend$treatments) + length(pre$cids) + length(pre$bids) + nt + nl + incidence),
    128 * nc, 256 * ng, 16 * row_names, 32 * name_bytes, index_attributes, 8 * link_bytes)))
  # Dense lambda/nu inputs, solved tilde arrays, inactive zero templates,
  # group/cell subset lists and the current scalar Gram RHS/factor copies.
  forcing <- A(add(c(8192, 64 * nc, 256 * ng, 8 * cell_columns,
    8 * group_columns, 16 * peak_p^2, 16 * name_bytes)))
  rates_design <- A(add(c(2048, 256 * G, 8 * name_bytes)),
    3 * G + ncol(X) + 16)
  # Every selected treatment retains direct/a/gamma/score_direct/score_a.
  # The last loop's five named matrices, X_trt_ix and logit/index/weight
  # locals survive during links, even when a later arm does not replace them.
  treatment <- A(add(c(1024, 256 * (nt + length(backend$treatments)), 40 * nt * name_bytes)),
    5 * nt * pt)
  previous_treatment <- if (nt) A(add(c(2048, 64 * pt, 16 * name_bytes)), 6 * pt + 24) else A()
  n_alpha <- if (is.null(backend$shared_or)) 0 else
    length(unique(vapply(backend$links[links], `[[`, integer(1), "or_block_id")))
  channels <- if (augmented) 4 else 3
  # All three/four returned OR channels, alpha/beta totals, theta, and both
  # returned RHS scores coexist. Alpha RHS arrays exist even without intercept.
  or_return <- if (length(pre$bids)) A(add(c(2048, 256 * (n_alpha + p + channels + 1),
    8 * (channels + 4) * name_bytes)), (channels + 3) * (n_alpha + p)) else A()
  nu <- if (augmented) length(pre$uids) else 0
  nz <- if (augmented) length(pre$zids) else 0
  W <- 0; selection_columns <- 0; selection <- A()
  if (augmented) {
    sel <- backend$selection_backend
    active <- which(vapply(backend$up_map$theta_coords, function(ids) any(ids %in% pre$uids), logical(1)))
    for (m in active) {
      bl <- sel$blocks[[m]]
      W <- add(c(W, if (bl$has_intercept) bl$B_m else 0, bl$d_s,
        length(backend$up_map$theta_coords[[m]])))
    }
    selection_columns <- W
    # Preserve the existing selection envelope in full. up_res constructs ALL
    # alpha/beta/theta block arrays, not merely the requested up coordinates.
    sw <- length(sel$coefficients) + G + 3 * length(pre$q_ids) + length(selection_pre$active_cells) + 32
    sr <- sw + 16 + frame$id_byte_metadata$max +
      32 * (sw + sel$design_spec$d + 8) + frame$id_byte_metadata$max
    selection <- A(selection_pre$output + selection_pre$scratch - B * sr, sr)
  }
  upstream <- A(if (augmented) 1024 + 256 * (length(backend$selection_backend$blocks) + nu) else 0,
    nu + W)
  # Reconstructed Z records and the last cell_Z ev (jhat/Z/obs/q/pi/e),
  # plus replacement prediction/reader/design/clip temporaries.
  zcache <- if (augmented) A(1024 + 256 * (nc + nz), 2 * nz) else A()
  zread <- if (augmented) A(4096 + 64 * name_bytes, 4 * peak_p + 32) else A()
  # Scalar names in the original link loop, including previous arm locals.
  # Full component matrices (3x8 pure, 1x10 augmented) are additional below.
  pure_vectors <- c("w_plus", "j_g_vec", "psi_denplus", "direct_T", "a_T", "j_T_adj", "j_T", "psi_T",
    "direct_Bg", "a_Bg", "or_Bg_dir", "or_Bg_a", "or_Bg_j", "psi_Bg", "psi_numBg", "m_g_vec",
    "direct_C", "a_C", "j_C", "j_C_adj", "trt_C_dir", "trt_C_a", "psi_C", "psi_denminus", "psi_numC",
    "direct_Bo", "a_Bo", "trt_Bo_dir", "trt_Bo_a", "or_Bo_dir", "or_Bo_a", "or_Bo_j", "psi_Bo", "psi_numBo",
    "w_minus", "j_k_vec", "m_k_vec", "a_denminus", "a_numC", "a_C_vec", "j_C_vec", "a_Bo_vec", "xi_ak", "wm_k")
  aug_vectors <- c("w_plus", "Z_g_vec", "m_g_vec", "psi_denplus", "direct_T", "a_T", "j_T_adj", "j_T", "sel_T", "psi_T",
    "direct_Bg", "a_Bg", "or_Bg_dir", "or_Bg_a", "or_Bg_j", "or_Bg_q", "psi_Bg", "w_minus", "Z_k_vec", "m_k_vec",
    "a_denminus", "a_C", "a_Bo", "j_C_fwd", "xi_ak", "direct_C", "j_C_adj", "j_C", "sel_C", "trt_C_dir", "trt_C_a",
    "psi_C", "direct_Bo", "trt_Bo_dir", "trt_Bo_a", "or_Bo_dir", "or_Bo_a", "or_Bo_j", "or_Bo_q", "psi_Bo")
  vectors <- if (augmented) length(aug_vectors) else length(pure_vectors)
  link_state <- A(add(c(8192, 256 * vectors, 64 * name_bytes)),
    vectors + (if (augmented) 10 else 24) + 2 * peak_p + 2 * p + pt + 24)
  # Pure core phi, three psi matrices and scaled phi; augmented phi/scaled phi.
  # Final wrapper/res output copies are charged separately by output below.
  link_matrices <- A(1024 + 32 * link_bytes, (if (augmented) 2 else 5) * nl)
  common <- plus(maps, forcing, cache, rates_design, treatment, previous_treatment,
    or_return, upstream, zcache, link_state, link_matrices)
  # Raw reader ix/hit/obs/rr/dy/Z/d/jhat, centered d, val/out, replacement
  # products and finite checks; scalar RHS/transpose/solve buffers are explicit.
  adjoint <- A(add(c(8192, 256 * (nc + ng), 16 * peak_p^2, 64 * name_bytes)),
    12 * peak_p + 32)
  # RHS channels remain live during nested scalar pullbacks. Pooled solve has
  # input transpose, pivot RHS, forward/back solutions, scale/reorder, returned
  # transpose, then alpha/theta/RHS-sum construction overlaps prior channels.
  or_construct <- if (length(pre$bids)) plus(adjoint, A(add(c(4096, 16 * p^2, 64 * name_bytes)),
    (channels + 4) * (n_alpha + p) + 16 * p + 2 * n_alpha + 24)) else A()
  trt_construct <- if (nt) A(add(c(4096, 16 * pt^2, 4 * pt * G, 64 * name_bytes)),
    12 * pt + 2 * G + 32) else A()
  load_theta <- if (augmented) A(2048 + 32 * nu, 2 * nu + 8) else A()
  phases <- list(selection = selection, adjoint = adjoint, pooled_or = or_construct,
    treatment = trt_construct, reconstructed_cell = zread, theta_loading = load_theta)
  phase <- A(max(vapply(phases, `[[`, numeric(1), "constant")), max(vapply(phases, `[[`, numeric(1), "row")))
  # Coefficientwise phase maximum is affine even if two phase lines cross.
  scratch <- plus(common, phase, output, A(0, 3 * nl))
  coefficients <- list(output = output, scratch = scratch, common = common, phase = phase,
    cache = cache, treatment = treatment, previous_treatment = previous_treatment,
    pooled_or_return = or_return, upstream = upstream, link_matrices = link_matrices)
  list(output = value(output), scratch = value(scratch), coefficients = coefficients, phases = phases,
    cache_full_needed_entries = value(cache), treatment_entries = value(treatment),
    previous_treatment_entries = value(previous_treatment), pooled_or_entries = value(or_return),
    selection_columns = selection_columns, upstream_entries = value(upstream),
    selection_peak = value(selection), rate_entries = value(rates_design),
    link_matrix_entries = value(link_matrices), maps_entries = value(maps), forcing_entries = value(forcing),
    owned_entries = backend$accounting$owned_entries, borrowed_entries = backend$accounting$borrowed_entries)
}

# Retain the legacy bound unless dominance at B=1 and a no-larger slope prove
# the replacement no larger on the evaluator's domain B>=1:
# old(B)-new(B) = old(1)-new(1) + (B-1)*(old.row-new.row) >= 0.
# The decision does not depend on B or optional cache membership. No optional
# optimization adds an admission refusal; the existing work charge is intact.
reconstruction_ordinary_if_envelope <- function(backend, B, links, pre, width, metadata, selection_pre = NULL) {
  x <- reconstruction_ordinary_inventory(backend, B, links, pre, selection_pre)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  augmented <- inherits(backend, "cdid_augmented_backend")
  r <- 64 * (width + length(backend$design_columns) + backend$outcome_recipe$p +
    backend$treatment_recipe$p + 8 + if (augmented) ncol(reader$state$Z) + length(pre$zids) else 0) +
    frame$id_byte_metadata$max
  if (augmented) r <- r + x$phases$selection[["row"]]
  old_output <- c(constant = metadata - B * frame$id_byte_metadata$max,
    row = width + 16 + frame$id_byte_metadata$max)
  old_scratch <- c(constant = pre$scratch - B * r, row = r)
  if (x$coefficients$output[["row"]] <= old_output[["row"]] &&
      sum(x$coefficients$output) <= sum(old_output)) pre$output <- x$output
  if (x$coefficients$scratch[["row"]] <= old_scratch[["row"]] &&
      sum(x$coefficients$scratch) <= sum(old_scratch)) pre$scratch <- x$scratch
  pre
}
