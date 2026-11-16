# Private native fast-kernel dispatch wrappers for selection models,
# Corollary 1 IPW/OR/DR consumers, custom selection, and D60/D61 reference inputs.
# All kernels fall back cleanly to NULL when unavailable or when any downstream
# helper is mocked/traced.

cdid_sel_cache <- new.env(parent = emptyenv())

cdid_sel_fast_kernel <- function(name) {
  owner_env <- environment(cdid_fast_kernel)
  if (cdid_is_default_fn(gmm_require_native_dll, "gmm_require_native_dll")) {
    syms <- gmm_native_dll_cache$syms
    if (is.list(syms) && identical(gmm_native_dll_cache$owner, owner_env)) {
      calls <- attr(syms, "calls")
      fn <- if (is.list(calls)) calls[[name]] else NULL
      if (is.function(fn) && exists(name, envir = environment(fn), inherits = FALSE)) return(fn)
    }
  }
  cdid_fast_kernel(name)
}

cdid_sel_fast_gc <- function(glm_control) {
  if (is.list(glm_control) && length(glm_control) == 3L &&
      identical(names(glm_control), c("epsilon", "maxit", "trace")) &&
      is.double(glm_control$epsilon) && length(glm_control$epsilon) == 1L &&
      is.finite(glm_control$epsilon) && glm_control$epsilon > 0 &&
      is.numeric(glm_control$maxit) && length(glm_control$maxit) == 1L &&
      is.finite(glm_control$maxit) && glm_control$maxit > 0 &&
      is.logical(glm_control$trace) && length(glm_control$trace) == 1L &&
      !is.na(glm_control$trace)) {
    return(glm_control)
  }
  do.call(stats::glm.control, glm_control)
}

cdid_sel_check_validated <- function(obj, class_name, seal_tag) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(FALSE)
  if (!inherits(obj, class_name)) return(FALSE)
  if (!cdid_are_default_fns(c(
    "ipw_owner_seal_ok",
    "validate_ipw_frame_owner",
    "validate_ipw_panel_owner",
    "validate_cp04_reader",
    "validate_chained_reference_reader",
    "validate_ipw_owner_control",
    "validate_ipw_glm_control",
    "validate_selection_models"
  ))) {
    return(FALSE)
  }
  ring <- get0("val_ring", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(ring) || !length(ring)) return(FALSE)
  fk <- cdid_sel_fast_kernel("C_cdid_verify_fast")
  if (is.null(fk)) return(FALSE)
  hit_ok <- FALSE
  for (i in seq_along(ring)) {
    rec <- ring[[i]]
    if (is.list(rec) && identical(rec[[2L]], seal_tag) && isTRUE(fk(4L, obj, rec, NULL, NULL))) {
      hit_ok <- TRUE
      break
    }
  }
  if (!hit_ok) return(FALSE)
  if (!is.null(.subset2(obj, "reference_ref"))) {
    chained_backend_reference(obj)
  }
  sb <- .subset2(obj, "selection_backend")
  if (!is.null(sb) && !identical(seal_tag, "selection_models")) {
    if (!cdid_sel_check_validated(sb, "cdid_selection_models", "selection_models")) return(FALSE)
  }
  TRUE
}

cdid_sel_record_validated <- function(obj, seal_tag, reader, frame, panel, fref) {
  rec <- list(obj, seal_tag, reader, frame, panel, fref)
  ring <- get0("val_ring", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(ring)) ring <- vector("list", 16L)
  oref <- .subset2(obj, "ref")
  slot <- NA_integer_
  for (i in seq_along(ring)) {
    ri <- ring[[i]]
    if (is.null(ri)) {
      if (is.na(slot)) slot <- i
    } else if (is.list(ri) && identical(.subset2(ri[[1L]], "ref"), oref)) {
      slot <- i
      break
    }
  }
  if (is.na(slot)) {
    pos <- get0("val_ring_pos", envir = cdid_sel_cache, inherits = FALSE)
    pos <- if (is.integer(pos) && length(pos) == 1L) (pos %% 16L) + 1L else 1L
    assign("val_ring_pos", pos, envir = cdid_sel_cache)
    slot <- pos
  }
  ring[[slot]] <- rec
  assign("val_ring", ring, envir = cdid_sel_cache)
  invisible(NULL)
}

cdid_sel_get_if_meta <- function(backend, B, links, tag) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  slots <- get0(tag, envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(slots)) return(NULL)
  bref <- .subset2(backend, "ref")
  B_int <- as.integer(B)
  links_int <- as.integer(links)
  for (i in seq_along(slots)) {
    hit <- slots[[i]]
    if (is.list(hit) && identical(hit$backend_ref, bref) &&
        gmm_sexp_aliases(hit$backend, backend) && identical(hit$B, B_int) &&
        identical(hit$links, links_int)) {
      return(hit)
    }
  }
  NULL
}

cdid_sel_set_if_meta <- function(backend, B, links, pre, link_names, tag) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  rec <- list(
    backend = backend,
    backend_ref = .subset2(backend, "ref"),
    B = as.integer(B),
    links = as.integer(links),
    pre = pre,
    link_names = link_names
  )
  slots <- get0(tag, envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(slots) || length(slots) != 2L) {
    slots <- list(rec, NULL)
  } else {
    slots <- list(rec, slots[[1L]])
  }
  assign(tag, slots, envir = cdid_sel_cache)
  rec
}

cdid_sel_valid_idx <- function(x, n) {
  if (is.integer(x)) {
    !is.object(x) && is.null(dim(x)) && length(x) > 0L && !anyNA(x) && (
      (!is.unsorted(x, strictly = TRUE) && x[1L] >= 1L && x[length(x)] <= n) ||
      (min(x) >= 1L && max(x) <= n && !anyDuplicated(x))
    )
  } else {
    is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) > 0L &&
      all(is.finite(x)) && all(x == floor(x)) && all(x >= 1 & x <= n) && !anyDuplicated(x)
  }
}

cdid_sel_slice_id <- function(id_vec, rows, N) {
  if (is.integer(rows) && length(rows) == N && is.null(attributes(rows)) &&
      is.null(attributes(id_vec)) && rows[1L] == 1L && rows[N] == N &&
      !is.unsorted(rows, strictly = TRUE)) {
    id_vec
  } else {
    id_vec[rows]
  }
}

cdid_fast_dr_sep_ord_envelope <- function(backend, B, links) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.null(attributes(links)) || !is.numeric(B) || length(B) != 1L ||
      !is.finite(B) || B < 1 || B > 1e7 ||
      !cdid_are_default_fns(c(
        "dr_conditional_separate_ordinary_if_envelope",
        "dr_conditional_separate_ordinary_inventory",
        "selection_if_budget",
        "chained_backend_reference",
        "cp04_budget",
        "cp04_payload",
        "ipw_owner_safe_sum",
        "ipw_owner_safe_product"
      ))) {
    return(NULL)
  }
  links_int <- as.integer(links)
  bref <- .subset2(backend, "ref")
  rec <- get0("dr_sep_env_static", envir = cdid_sel_cache, inherits = FALSE)
  if (!(is.list(rec) && identical(rec$backend_ref, bref) &&
        gmm_sexp_aliases(rec$backend, backend) && identical(rec$links, links_int))) {
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    G <- length(backend$regression_specs); S <- length(backend$supports)
    p <- backend$outcome_recipe$p; pt <- backend$treatment_recipe$p
    D <- backend$upstream$D; Tn <- length(backend$treatments)
    nl <- length(links_int); J <- length(st$cohorts)
    incidence <- 0
    for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
    for (s in backend$treatment_specs) incidence <- incidence + length(s$support_ids)
    for (u in backend$supports) incidence <- incidence + length(u$channels)
    for (m in backend$treatments) incidence <- incidence + length(m$upstream_ids)
    maps <- 4096 + 32 * (G + S + D + Tn + nl + incidence)
    cp04_budget(
      backend$control,
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      maps,
      c(backend$accounting$work_ops_bound, 64 * maps),
      0,
      "separate DR dependency metadata"
    )
    blinks <- backend$links[links_int]
    link_blocks <- vapply(blinks, `[[`, integer(1L), "regression")
    gids <- sort.int(unique(link_blocks))
    tids <- sort.int(unique(vapply(blinks, `[[`, integer(1L), "model")))
    sids <- sort.int(unique(c(
      vapply(blinks, `[[`, integer(1L), "support_id"),
      unlist(lapply(backend$regression_specs[gids], `[[`, "support_ids"), use.names = FALSE),
      unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"), use.names = FALSE)
    )))
    uids <- sort.int(unique(c(
      unlist(lapply(backend$supports[sids], `[[`, "channels"), use.names = FALSE),
      unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"), use.names = FALSE)
    )))
    ns <- length(sids); nu <- length(uids); ng <- length(gids); nt <- length(tids)
    ps <- vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1L))
    model_u <- vapply(tids, function(t) length(backend$treatments[[t]]$upstream_ids), integer(1L))
    max_u <- max(c(0, model_u))
    ref <- chained_backend_reference(backend)
    p_ref <- if (is.null(ref)) 0L else length(ref$x_columns)
    dr_work_0 <- 128 * maps
    dr_work_B <- 128 * ((ns + 1) * (p + 2)^2 * (1 + ceiling(log2(st$n_A + 1))) +
      sum((ps + 2)^2 * (nu + nl + 2)) + nt * (pt + 2) * (pt + max_u + 8) +
      (ns + nt + 1) * p_ref + J + nl * (p + pt + 8))

    sel <- backend$selection_backend
    sel_q_ids <- unique(backend$q_of_s[sids])
    sp_active_blocks <- unique(unlist(lapply(sel$pair_cohorts[sel_q_ids], `[[`, "block_ids"), use.names = FALSE))
    sp_active_cells <- unique(unlist(lapply(sel$block_specs[sp_active_blocks], `[[`, "cell_ids"), use.names = FALSE))
    n_theta <- length(sel$coefficients)
    n_probes <- length(sel_q_ids)
    sp_widths <- n_theta + J + 3 * n_probes + length(sp_active_cells) + 32
    id_max_bytes <- frame$id_byte_metadata$max
    sp_meta_0 <- 8192 + 1024 * (length(sel$blocks) + length(sel$cells) + J) +
      n_probes * (sel$design_spec$d_s + 1024)
    sp_retained <- sel$accounting$owned_entries + sel$accounting$borrowed_entries
    sp_per_row_work <- 32 * length(sp_active_cells) * (1 + ceiling(log2(st$N + 1))) +
      32 * sum(vapply(sel$blocks[sp_active_blocks], function(bl) bl$B_m * (bl$d_s + 2)^2, numeric(1L))) +
      32 * (n_theta + J + 2 * n_probes + 8)
    sp_target_scan <- 32 * st$N * (sel$design_spec$d + 4)
    sp_d_raw <- sel$design_spec$d
    sp_d_s <- sel$design_spec$d_s

    map <- backend$up_map; M <- length(sel$blocks)
    active <- which(vapply(map$theta_coords, function(ids) any(ids %in% uids), logical(1L)))
    W <- 0; cells <- 0
    for (m in active) {
      bl <- sel$blocks[[m]]
      W <- W + (if (bl$has_intercept) bl$B_m else 0) + bl$d_s + length(map$theta_coords[[m]])
      cells <- cells + length(bl$cell_ids)
    }
    ref_x <- if (is.null(ref)) 0 else length(ref$x_columns)
    ref_s <- if (is.null(ref)) 0 else length(ref$designs[[ref$role_slot["s"]]]$columns)
    d <- max(1L, sp_d_raw)
    or_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
    trt_bytes <- sum(nchar(backend$treatment_recipe$columns, type = "bytes"))
    design_bytes <- sum(nchar(if (is.null(ref)) colnames(st$Z) else
      c(ref$x_columns, ref$designs[[ref$role_slot["s"]]]$columns), type = "bytes"))
    link_bytes <- sum(nchar(vapply(blinks, `[[`, character(1L), "link_id"), type = "bytes"))
    requested_payload <- unname(cp04_payload(blinks)["entries"])
    out_meta_0 <- 8192 + 256 * nl + 8 * link_bytes
    upstream_meta <- 1024 + 64 * (M + cells + length(active) + nu)
    hits_meta <- 1024 + 256 * ns + 8 * ns * or_bytes
    metadata_val <- maps + 8192 + 64 * (ns + nu + ng + nt + nl + M + cells + 1) +
      2 * requested_payload + 8 * link_bytes
    support_state_0 <- 1024 + 16 * (or_bytes + design_bytes)
    common_0 <- metadata_val + 2048 + upstream_meta + hits_meta + support_state_0
    common_B <- (nu + W) + ns * (p + 4) + (p + 16)
    read_0 <- 4096 + 64 * (or_bytes + design_bytes)
    read_B <- 32 * (p + d + ref_x + ref_s + 16)

    supp_ch_len <- vapply(backend$supports, function(u) length(u$channels), integer(1L))
    link_sids <- vapply(blinks, `[[`, integer(1L), "support_id")
    ret0 <- numeric(ng); retB <- numeric(ng)
    con0 <- numeric(ng); conB <- numeric(ng)
    for (idx in seq_len(ng)) {
      b <- gids[idx]
      pp <- ps[idx]
      in_b <- which(link_blocks == b)
      ll <- length(in_b)
      rr <- length(backend$regressions[[b]]$channel_ids)
      r_sids <- backend$regression_specs[[b]]$support_ids
      ss <- length(r_sids)
      h <- max(0L, supp_ch_len[r_sids], supp_ch_len[link_sids[in_b]])
      r0 <- 2048 + 256 * (ll + 1) + 16 * (ll + 4) * or_bytes
      rB <- 4 * ll + ll * pp + 2 * pp + 4
      ret0[idx] <- r0; retB[idx] <- rB
      con0[idx] <- r0 + 2 * pp * rr + 32 * pp * pp + 4096 + 128 * (ss + rr + h + 1) + 16 * (ll + 8) * or_bytes
      conB[idx] <- rB + 20 * pp + 4 * rr + 4 * h + 48
    }
    st0 <- numeric(nt); stB <- numeric(nt)
    mcon0 <- numeric(nt); mconB <- numeric(nt)
    base_bytes_sum <- or_bytes + trt_bytes + design_bytes
    for (idx in seq_len(nt)) {
      t <- tids[idx]
      u <- model_u[idx]
      pp <- length(backend$outcome_groups$groups[[backend$treatment_specs[[t]]$regression]]$keep)
      s0 <- 2048 + 64 * (u + pt + pp + 1) + 16 * base_bytes_sum
      sB <- u + 3 * pt + pp + 32
      st0[idx] <- s0; stB[idx] <- sB
      mcon0[idx] <- s0 + 16 * pt * pt + 4 * pt * u + 4096 + 64 * base_bytes_sum
      mconB[idx] <- sB + 16 * pt + 4 * pp + 4 * u + 2 * ref_x + 2 * ref_s + 4 * d + 64
    }
    rec <- list(
      backend = backend, backend_ref = bref, links = links_int,
      maps = maps, gids = gids, tids = tids, sids = sids, uids = uids,
      ng = ng, nt = nt, nl = nl, id_max_bytes = id_max_bytes,
      dr_work_0 = dr_work_0, dr_work_B = dr_work_B,
      sel_ctrl = sel$control, sp_retained = sp_retained,
      sp_meta_0 = sp_meta_0, sp_widths = sp_widths, sp_d_raw = sp_d_raw, sp_d_s = sp_d_s,
      sp_per_row_work = sp_per_row_work, sp_target_scan = sp_target_scan,
      out_meta_0 = out_meta_0, common_0 = common_0, common_B = common_B,
      read_0 = read_0, read_B = read_B,
      ret0 = ret0, retB = retB, con0 = con0, conB = conB,
      st0 = st0, stB = stB, mcon0 = mcon0, mconB = mconB
    )
    assign("dr_sep_env_static", rec, envir = cdid_sel_cache)
  } else {
    cp04_budget(
      backend$control,
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      rec$maps,
      c(backend$accounting$work_ops_bound, 64 * rec$maps),
      0,
      "separate DR dependency metadata"
    )
  }
  B_dbl <- as.double(B)
  sp_meta <- rec$sp_meta_0 + B_dbl * rec$id_max_bytes
  sp_output <- B_dbl * (rec$sp_widths + 16) + sp_meta
  sp_scratch <- 32 * B_dbl * (rec$sp_widths + rec$sp_d_raw + 8) + 64 * (rec$sp_d_s + 2)^2 + sp_meta
  sp_work <- B_dbl * rec$sp_per_row_work + rec$sp_target_scan + 64 * sp_meta
  cp04_budget(rec$sel_ctrl, rec$sp_retained, sp_scratch, sp_work, sp_output, "selection IF chunk")

  work <- rec$dr_work_0 + B_dbl * rec$dr_work_B + sp_work
  output <- rec$out_meta_0 + B_dbl * rec$id_max_bytes + B_dbl * (rec$nl + 4)
  common <- rec$common_0 + B_dbl * rec$common_B
  read_peak <- rec$read_0 + B_dbl * rec$read_B
  core_return <- if (rec$ng > 0L) max(0, rec$ret0 + B_dbl * rec$retB) else 0
  core_peak <- if (rec$ng > 0L) max(0, rec$con0 + B_dbl * rec$conB) else 0
  model_state <- if (rec$nt > 0L) max(0, rec$st0 + B_dbl * rec$stB) else 0
  model_peak <- if (rec$nt > 0L) max(0, rec$mcon0 + B_dbl * rec$mconB) else 0
  selection_peak <- sp_output + sp_scratch
  phase_peak <- max(selection_peak, read_peak, model_state + core_peak, core_return + model_peak)
  scratch <- common + phase_peak + output + 3 * B_dbl * rec$nl
  if (!is.finite(scratch) || !is.finite(work) || !is.finite(output) ||
      scratch > 1e15 || work > 1e15 || output > 1e15) {
    return(NULL)
  }
  list(
    gids = rec$gids, beta_gids = rec$gids, tids = rec$tids,
    sids = rec$sids, uids = rec$uids,
    output = output, scratch = scratch, work = work
  )
}

cdid_fast_dr_pool_ord_envelope <- function(backend, B, links) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.null(attributes(links)) || !is.numeric(B) || length(B) != 1L ||
      !is.finite(B) || B < 1 || B > 1e7 ||
      !cdid_are_default_fns(c(
        "or_pooled_if_sets",
        "dr_corollary1_ordinary_if_envelope",
        "dr_corollary1_ordinary_inventory",
        "selection_if_budget",
        "chained_backend_reference",
        "cp04_budget",
        "cp04_payload",
        "ipw_owner_safe_sum",
        "ipw_owner_safe_product"
      ))) {
    return(NULL)
  }
  links_int <- as.integer(links)
  bref <- .subset2(backend, "ref")
  rec <- get0("dr_pool_env_static", envir = cdid_sel_cache, inherits = FALSE)
  if (!(is.list(rec) && identical(rec$backend_ref, bref) &&
        gmm_sexp_aliases(rec$backend, backend) && identical(rec$links, links_int))) {
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    p <- backend$shared$p
    pt <- backend$treatment_recipe$p
    d <- max(1L, backend$upstream$d)
    G <- length(backend$regression_specs)
    J <- length(st$cohorts)
    sel <- backend$selection_backend
    up_map <- backend$up_map

    sets <- or_pooled_if_sets(backend, links_int, FALSE)
    gids <- sets$gids
    beta_gids <- sets$beta_gids
    blinks <- backend$links[links_int]
    tids <- sort.int(unique(vapply(blinks, `[[`, integer(1L), "model")))
    trt_sids <- sort.int(unique(unlist(lapply(backend$treatment_specs[tids], `[[`, "support_ids"), use.names = FALSE)))
    sids <- sort.int(unique(c(sets$sids, trt_sids)))
    uids <- sort.int(unique(c(
      unlist(lapply(backend$supports[sids], `[[`, "channels"), use.names = FALSE),
      unlist(lapply(backend$treatments[tids], `[[`, "upstream_ids"), use.names = FALSE)
    )))
    active_blocks <- which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1L)))
    active_cells <- unique(unlist(lapply(sel$block_specs[active_blocks], `[[`, "cell_ids"), use.names = FALSE))

    nl <- length(links_int); ns <- length(sids); nt <- length(tids)
    nu <- length(uids); nc <- length(active_cells)
    M <- length(sel$blocks); Tn <- length(backend$treatments)
    id_max_bytes <- frame$id_byte_metadata$max
    link_bytes <- sum(nchar(vapply(blinks, `[[`, character(1L), "link_id"), type = "bytes"))
    metadata_0 <- 16384 + 1024 * (ns + nu + G + nt + nl + J + M) + 8 * link_bytes
    log_a <- 1 + ceiling(log2(st$n_A + 1))
    dr_work_0 <- 64 * metadata_0
    dr_work_B <- 128 * (
      (ns + 1) * ((p + 2)^2 + (d + 2)^2) * log_a +
        (length(beta_gids) + length(gids)) * (p + 1) * (p + nu) +
        nt * ((pt + 2)^2 * log_a + pt * (nu + 4)) +
        J + nl * (p + pt + d + 4)
    ) + 64 * id_max_bytes

    sel_q_ids <- unique(backend$q_of_s[sids])
    sp_active_blocks <- unique(unlist(lapply(sel$pair_cohorts[sel_q_ids], `[[`, "block_ids"), use.names = FALSE))
    sp_active_cells <- unique(unlist(lapply(sel$block_specs[sp_active_blocks], `[[`, "cell_ids"), use.names = FALSE))
    n_theta <- length(sel$coefficients)
    n_probes <- length(sel_q_ids)
    sp_widths <- n_theta + J + 3 * n_probes + length(sp_active_cells) + 32
    sp_meta_0 <- 8192 + 1024 * (M + length(sel$cells) + J) + n_probes * (sel$design_spec$d_s + 1024)
    sp_retained <- sel$accounting$owned_entries + sel$accounting$borrowed_entries
    sp_per_row_work <- 32 * length(sp_active_cells) * (1 + ceiling(log2(st$N + 1))) +
      32 * sum(vapply(sel$blocks[sp_active_blocks], function(bl) bl$B_m * (bl$d_s + 2)^2, numeric(1L))) +
      32 * (n_theta + J + 2 * n_probes + 8)
    sp_target_scan <- 32 * st$N * (sel$design_spec$d + 4)
    sp_d_raw <- sel$design_spec$d
    sp_d_s <- sel$design_spec$d_s

    W <- 0
    for (m in active_blocks) {
      bl <- sel$blocks[[m]]
      W <- W + (if (bl$has_intercept) bl$B_m else 0) + bl$d_s + length(up_map$theta_coords[[m]])
    }
    I_beta <- 0
    if (p) for (b in beta_gids) I_beta <- I_beta + length(backend$regression_specs[[b]]$support_ids)
    max_u <- 0
    for (t in tids) max_u <- max(max_u, length(backend$treatments[[t]]$upstream_ids))
    ref <- chained_backend_reference(backend)
    ref_x <- if (is.null(ref)) 0 else length(ref$x_columns)
    ref_s <- if (is.null(ref)) 0 else length(ref$designs[[ref$role_slot["s"]]]$columns)
    or_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
    other_bytes <- sum(nchar(c(backend$treatment_recipe$columns,
      if (is.null(ref)) colnames(st$Z) else c(ref$x_columns, ref$designs[[ref$role_slot["s"]]]$columns)), type = "bytes"))

    out_meta_0 <- 8192 + 256 * nl + 8 * link_bytes
    common_0 <- metadata_0 + 256 * (nt + Tn + M + nc + ns + 1)
    common_B <- id_max_bytes + (nu + W) + (if (p) ns * (p + 4) else 0) + 5 * nt + p
    core_0 <- if (p) (8 * p * nu + 16 * p * p + 256 * I_beta + 8 * (ns + I_beta + 4) * or_bytes + 64 * other_bytes) else 0
    core_B <- if (p) (I_beta * (p + 3) + 16 * p + 4 * nu) else 0
    model_0 <- 16 * pt * pt + 8 * pt * max_u + 64 * max_u + 64 * (or_bytes + other_bytes)
    model_B <- 32 * (p + pt + max_u + d + ref_x + ref_s + 32)

    rec <- list(
      backend = backend, backend_ref = bref, links = links_int,
      gids = gids, beta_gids = beta_gids, sids = sids, tids = tids,
      uids = uids, active_blocks = active_blocks, active_cells = active_cells,
      nl = nl, id_max_bytes = id_max_bytes,
      dr_work_0 = dr_work_0, dr_work_B = dr_work_B,
      sel_ctrl = sel$control, sp_retained = sp_retained,
      sp_meta_0 = sp_meta_0, sp_widths = sp_widths, sp_d_raw = sp_d_raw, sp_d_s = sp_d_s,
      sp_per_row_work = sp_per_row_work, sp_target_scan = sp_target_scan,
      out_meta_0 = out_meta_0, common_0 = common_0, common_B = common_B,
      core_0 = core_0, core_B = core_B, model_0 = model_0, model_B = model_B
    )
    assign("dr_pool_env_static", rec, envir = cdid_sel_cache)
  }
  B_dbl <- as.double(B)
  sp_meta <- rec$sp_meta_0 + B_dbl * rec$id_max_bytes
  sp_output <- B_dbl * (rec$sp_widths + 16) + sp_meta
  sp_scratch <- 32 * B_dbl * (rec$sp_widths + rec$sp_d_raw + 8) + 64 * (rec$sp_d_s + 2)^2 + sp_meta
  sp_work <- B_dbl * rec$sp_per_row_work + rec$sp_target_scan + 64 * sp_meta
  cp04_budget(rec$sel_ctrl, rec$sp_retained, sp_scratch, sp_work, sp_output, "selection IF chunk")

  work <- rec$dr_work_0 + B_dbl * rec$dr_work_B + sp_work
  output <- rec$out_meta_0 + B_dbl * rec$id_max_bytes + B_dbl * (rec$nl + 4)
  common <- rec$common_0 + B_dbl * rec$common_B
  selection_peak <- sp_output + sp_scratch
  core_peak <- rec$core_0 + B_dbl * rec$core_B
  model_peak <- rec$model_0 + B_dbl * rec$model_B
  scratch <- common + max(selection_peak, core_peak, model_peak) + output + 3 * B_dbl * rec$nl
  if (!is.finite(scratch) || !is.finite(work) || !is.finite(output) ||
      scratch > 1e15 || work > 1e15 || output > 1e15) {
    return(NULL)
  }
  list(
    gids = rec$gids, beta_gids = rec$beta_gids, sids = rec$sids, tids = rec$tids,
    uids = rec$uids, active_blocks = rec$active_blocks, active_cells = rec$active_cells,
    output = output, scratch = scratch, work = work
  )
}

cdid_fast_sel_ref_lookup <- function(keys, table) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !length(keys) || is.object(keys) || is.object(table) || !is.double(table)) return(NULL)
  fk <- cdid_sel_fast_kernel("C_cdid_sel_ref_lookup")
  if (is.null(fk)) return(NULL)
  fk(keys, table)
}

cdid_fast_sel_ref_slice <- function(mode, rows, base_off, N_val, table, data, cols = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns("chained_reference_lookup")) return(NULL)
  fk <- cdid_sel_fast_kernel("C_cdid_sel_ref_slice")
  if (is.null(fk)) return(NULL)
  fk(as.integer(mode), rows, c(as.double(base_off), as.double(N_val)), table, data, cols)
}

cdid_fast_sel_date_rows <- function(reader, panel, rows, dt, t_val) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
    "ipw_owner_match_sorted_double",
    "ipw_owner_match_sorted_id"
  ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_date_rows")
  if (is.null(fk)) return(NULL)
  st <- reader$state
  post <- panel$postings
  bounds2 <- if (post$date_count[dt] > 0L) as.integer(c(post$date_start[dt], post$date_end[dt])) else integer(0)
  sing <- NULL
  if (isTRUE(panel$n_singleton_units > 0L)) {
    if (!identical(panel$id_type, "numeric")) return(NULL)
    sing <- list(
      st$id,
      panel$unit_id_key,
      panel$units$n_obs,
      panel$units$obs_start,
      panel$observations$calendar_index
    )
  }
  R_vec <- fk(rows, st$a_index, post$post_a_index, bounds2, sing, as.integer(dt))
  if (is.null(R_vec)) return(NULL)
  list(R = R_vec, calendar_index = dt, t = t_val)
}

cdid_fast_sel_cp04_fit <- function(reader, spec, glm_control, reference = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
    "cp04_pair_rows",
    "cp04_logit",
    "ipw_glm_fit_internal",
    "ipw_logit_values",
    "chained_reference_s",
    "chained_reference_stat",
    "chained_reference_weight"
  )) || !is.list(glm_control)) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_fit_block")
  if (is.null(fk)) return(NULL)
  gc <- cdid_sel_fast_gc(glm_control)
  if (isTRUE(gc$trace)) return(NULL)
  params <- c(3, 0, gc$epsilon, as.double(gc$maxit))
  res <- fk(reader, spec, integer(0), NULL, reference, params)
  if (is.null(res)) return(NULL)
  H_chol <- res[[2L]]
  cn <- if (is.null(reference)) colnames(reader$state$Z) else
    reference$designs[[reference$role_slot["s"]]]$columns
  if (!is.null(cn)) dimnames(H_chol) <- list(cn, cn)
  list(
    coefficients = res[[1L]],
    H_chol = H_chol,
    diagnostics = list(
      converged = TRUE,
      boundary = FALSE,
      rank = length(res[[1L]]),
      iter = res[[3L]],
      warnings = character(0)
    ),
    pi = res[[4L]],
    r = res[[5L]]
  )
}

cdid_fast_sel_block_fit <- function(reader, block_spec, cells, design_spec,
                                    glm_control, ctrl, selection_model,
                                    selection_pooling, reference = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
    "selection_owner_block_provider",
    "cp04_selection_fit",
    "cp04_pair_rows",
    "cp04_date_rows",
    "cp04_logit",
    "ipw_glm_fit_internal",
    "ipw_logit_values",
    "fit_selection_block_provider",
    "or_pooled_wls_fit",
    "or_wls_fit",
    "chained_reference_s",
    "chained_reference_stat",
    "chained_reference_weight"
  )) || !is.list(glm_control)) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_fit_block")
  if (is.null(fk)) return(NULL)
  gc <- if (cdid_is_default_fn(validate_ipw_glm_control, "validate_ipw_glm_control") &&
            length(glm_control) == 3L && identical(names(glm_control), c("epsilon", "maxit", "trace")) &&
            is.integer(glm_control$maxit)) {
    cdid_sel_fast_gc(glm_control)
  } else {
    validate_ipw_glm_control(glm_control)
  }
  if (isTRUE(gc$trace)) return(NULL)
  cids <- as.integer(block_spec$cell_ids)
  B_m <- length(cids)
  if (B_m < 1L) return(NULL)
  is_date <- identical(selection_model, "independent_period_logit")
  is_sep_pair <- identical(selection_model, "pair_logit") &&
    identical(selection_pooling, "separate") && B_m == 1L && design_spec$d > 0L

  if (is_sep_pair) {
    cell <- cells[[cids[1L]]]
    sp <- list(pair = cell$pair, cohort = cell$cohort)
    if (!is.null(cell$reference)) sp$reference <- cell$reference
    params <- c(2, 0, gc$epsilon, as.double(gc$maxit))
    res <- fk(reader, sp, cids, design_spec, reference, params)
    if (is.null(res)) return(NULL)
    coef <- res[[1L]]
    H_chol <- res[[2L]]
    cn <- if (is.null(reference)) colnames(reader$state$Z) else
      reference$designs[[reference$role_slot["s"]]]$columns
    if (!is.null(cn)) dimnames(H_chol) <- list(cn, cn)
    has_int <- design_spec$has_intercept
    fit_rec <- list(
      coefficients = coef,
      H_chol = H_chol,
      diagnostics = list(
        converged = TRUE,
        boundary = FALSE,
        rank = length(coef),
        iter = res[[3L]],
        warnings = character(0)
      ),
      pi = res[[4L]],
      r = res[[5L]]
    )
    out <- selection_cp04_block(fit_rec, design_spec, res[[6L]], res[[7L]])
    out$block_id <- block_spec$block_id
    out$cell_ids <- cids
    return(out)
  }

  if (!design_spec$has_intercept && design_spec$d_s == 0L) return(NULL)
  params <- c(1, if (is_date) 1 else 0, gc$epsilon, as.double(gc$maxit))
  res <- fk(reader, cells, cids, design_spec, reference, params)
  if (is.null(res)) return(NULL)

  alpha <- res[[1L]]
  beta <- res[[2L]]
  d_s <- design_spec$d_s
  has_int <- design_spec$has_intercept
  p_m <- (if (has_int) B_m else 0L) + d_s
  slope_fit <- if (d_s > 0L) {
    structure(
      list(
        coefficients = numeric(d_s),
        rss = 0,
        R = res[[5L]],
        pivot = seq_len(d_s),
        n = res[[9L]],
        p = d_s,
        rank = d_s,
        rank_tol = res[[10L]],
        cutoff = res[[11L]]
      ),
      class = "cdid_or_wls_fit"
    )
  } else {
    NULL
  }
  list(
    branch = design_spec$branch,
    has_intercept = has_int,
    B_m = B_m,
    d_s = d_s,
    alpha = alpha,
    beta = beta,
    coefficients = c(alpha, beta),
    H_chol = NULL,
    m_b = res[[3L]],
    z_bar = res[[4L]],
    slope_fit = slope_fit,
    pi = res[[6L]],
    r = res[[7L]],
    diagnostics = list(
      converged = TRUE,
      boundary = FALSE,
      rank = p_m,
      iter = res[[8L]],
      warnings = character(0)
    ),
    block_id = block_spec$block_id,
    cell_ids = cids
  )
}

cdid_fast_sel_upstream <- function(sel_fit, up_map, reader, ref, rows, uids, prepared = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
    "cp04_pair_read_operation",
    "cp04_pair_rows",
    "cp04_date_rows",
    "evaluate_selection_block_if",
    "or_pooled_wls_gram_solve",
    "or_wls_gram_solve",
    "chained_reference_s",
    "chained_reference_stat",
    "chained_reference_weight"
  ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_upstream_chunk")
  if (is.null(fk)) return(NULL)
  res <- fk(sel_fit, up_map, reader, ref, rows, as.integer(uids), prepared)
  if (is.null(res)) return(NULL)
  active_blocks <- if (!is.null(prepared)) prepared$active_blocks else
    which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1)))
  active_cells <- if (!is.null(prepared)) unique(prepared$active_cells) else
    if (length(active_blocks)) unique(unlist(lapply(sel_fit$blocks[active_blocks], `[[`, "cell_ids"), use.names = FALSE)) else integer(0)
  list(
    up = res[[1L]],
    block_ifs = vector("list", length(sel_fit$blocks)),
    active_blocks = active_blocks,
    active_cells = active_cells
  )
}

cdid_fast_sel_trt_fit <- function(mode, reader, spec, sel, up_map = NULL,
                                  glm_control = NULL, reference = NULL,
                                  or_info = NULL, cols_info = NULL) {
  req_fns <- if (mode == 1L) {
    c("cp04_pair_rows", "cp04_predict_selection", "cp04_logit",
      "ipw_glm_fit_internal", "ipw_logit_values", "ipw_clip_summary")
  } else if (mode == 4L) {
    c("cp04_pair_rows", "selection_pair_cohort_log_deriv", "selection_prediction_operation",
      "reference_ipw_endpoint_outcomes", "ipw_owner_read_y_positions",
      "chained_reference_x", "chained_reference_s", "chained_reference_stat", "chained_reference_weight")
  } else {
    c("cp04_pair_rows", "selection_prediction_operation", "cp04_logit",
      "ipw_glm_fit_internal", "ipw_logit_values", "ipw_clip_summary",
      "chained_reference_x", "chained_reference_s", "chained_reference_stat", "chained_reference_weight")
  }
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(req_fns)) return(NULL)
  fk <- cdid_sel_fast_kernel("C_cdid_sel_trt_fit")
  if (is.null(fk)) return(NULL)

  if (mode == 4L) {
    return(fk(4L, reader, spec, sel, up_map, reference, NULL, NULL, NULL))
  }
  if (!is.list(glm_control)) return(NULL)
  gc <- cdid_sel_fast_gc(glm_control)
  if (isTRUE(gc$trace)) return(NULL)
  ctrl_vec <- c(as.double(gc$epsilon), as.double(gc$maxit))

  res <- fk(as.integer(mode), reader, spec, sel, up_map, reference, or_info, cols_info, ctrl_vec)
  if (is.null(res)) return(NULL)

  beta <- res[[1L]]
  H_chol <- res[[2L]]
  U_main <- res[[5L]]
  V_main <- res[[6L]]
  dimnames(U_main) <- list(c("T", "M", "dp", "dm"), NULL)
  trt_col_names <- if (!is.null(or_info) && !is.null(or_info$trt_cols)) or_info$trt_cols else NULL
  dimnames(V_main) <- list(c("T", "M", "dp", "dm"), trt_col_names)
  if (!is.null(trt_col_names)) {
    names(beta) <- trt_col_names
    dimnames(H_chol) <- list(trt_col_names, trt_col_names)
  }

  fit <- list(
    coefficients = beta,
    H_chol = H_chol,
    diagnostics = list(
      converged = TRUE,
      boundary = FALSE,
      rank = length(beta),
      iter = res[[3L]],
      warnings = character(0)
    )
  )
  if (mode >= 2L) fit$upstream_ids <- res[[9L]]
  fit$cross <- res[[4L]]
  fit$partials <- list(upstream = U_main, beta = V_main)

  if (mode != 3L) {
    means <- res[[7L]]
    names(means) <- c("T", "M", "dp", "dm")
    fit$means <- means
  } else {
    U_ipw <- res[[10L]]
    V_ipw <- res[[11L]]
    dimnames(U_ipw) <- list(c("T", "M", "dp", "dm"), NULL)
    dimnames(V_ipw) <- list(c("T", "M", "dp", "dm"), trt_col_names)
    fit$partials_ipw <- list(upstream = U_ipw, beta = V_ipw)
    means <- res[[7L]]
    names(means) <- c("T", "M", "dp", "dm", "T_dr", "M_dr", "mu_dr",
                      "T_ipw", "M_ipw", "mu_ipw", "mu_or", "mu_dr_decomp")
    fit$means <- means
    dbp <- res[[12L]]; dbm <- res[[13L]]; dbd <- res[[14L]]
    if (length(dbp) > 0L && !is.null(or_info$or_cols)) {
      names(dbp) <- names(dbm) <- names(dbd) <- or_info$or_cols
    }
    fit$d_bar_plus <- dbp
    fit$d_bar_minus <- dbm
    fit$d_bar_diff <- dbd
  }
  fit$clipping <- res[[8L]]
  fit
}

cdid_fast_sel_support_collect <- function(mode, reader, spec, sel, up_map = NULL,
                                          q_id = NULL, recipe, ctrl,
                                          reference = NULL, or_frame_index = NULL) {
  req_fns <- if (mode == 1L) {
    c("or_baseline_visit", "cp04_pair_rows", "cp04_pair_read_operation",
      "ipw_owner_read_y_positions", "cp04_predict_selection", "or_pooled_support_add",
      "or_pooled_wls_block_summary", "or_pooled_wls_merge_chunks", "or_wls_summary",
      "or_wls_merge", "or_wls_compress", "or_corollary1_moments_add",
      "ipw_clip_summary", "ipw_clip_summary_merge")
  } else if (mode == 2L) {
    c("cp04_pair_rows", "selection_prediction_operation", "or_pooled_support_add",
      "or_pooled_wls_block_summary", "or_pooled_wls_merge_chunks", "or_wls_summary",
      "or_wls_merge", "or_wls_compress", "or_selection_moments_add",
      "ipw_clip_summary", "ipw_clip_summary_merge",
      "chained_reference_design", "chained_reference_x", "chained_reference_s",
      "chained_reference_stat", "chained_reference_weight")
  } else {
    c("cp04_pair_rows", "selection_prediction_operation", "or_pooled_support_add",
      "or_pooled_wls_block_summary", "or_pooled_wls_merge_chunks", "or_wls_summary",
      "or_wls_merge", "or_wls_compress", "or_selection_moments_add",
      "chained_reference_x", "chained_reference_s", "chained_reference_stat",
      "chained_reference_weight")
  }
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(req_fns)) return(NULL)
  fk <- cdid_sel_fast_kernel("C_cdid_sel_support_collect")
  if (is.null(fk)) return(NULL)

  st <- reader$state
  frame <- get0("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(frame)) return(NULL)
  z_idx <- if (mode == 3L) {
    as.integer(or_frame_index)
  } else if (mode == 2L && !is.null(reference)) {
    if (is.null(reference$or_index)) return(NULL)
    as.integer(reference$or_index)
  } else {
    NULL
  }
  icpt_flag <- isTRUE(recipe$intercept)
  p <- as.integer(recipe$p)
  L_val <- as.double(spec$transitions)
  cols_info <- list(
    Z_a = if (is.null(reference)) frame$Z_a else NULL,
    z_idx = z_idx,
    slope_idx = as.integer(recipe$slope_index),
    icpt_flag = as.integer(icpt_flag),
    icpt_col = if (icpt_flag) as.integer(recipe$intercept_index) else 0L,
    L_val = L_val,
    max_chunk_rows = as.integer(ctrl$max_chunk_rows)
  )
  res <- fk(as.integer(mode), reader, spec, sel, up_map, q_id, reference, cols_info)
  if (is.null(res)) return(NULL)

  a_val <- if (is.null(reference)) st$a[spec$cohort] else
    chained_reference_stat(reference, spec$reference, spec$cohort, "a")
  fac <- res[[1L]]
  sc <- res[[6L]]
  n_obs_int <- as.integer(sc[1L])
  n_obs_dbl <- as.double(sc[1L])
  max_merge <- as.integer(sc[2L])
  n_comp <- as.integer(sc[3L])
  weight_mass <- sc[4L]
  sumY <- sc[5L]
  slope_names <- if (p > 0L) recipe$slope_columns else NULL
  k_aug <- p + 1L + as.integer(icpt_flag)

  raw_sum <- structure(
    list(factor = fac, n = n_obs_int, columns = if (icpt_flag) c("", slope_names) else slope_names),
    class = "cdid_or_wls_summary"
  )
  summary_obj <- if (icpt_flag) {
    structure(
      list(
        raw = raw_sum,
        n = n_obs_int,
        p = p,
        columns = slope_names,
        anchor = res[[2L]],
        constant = res[[3L]],
        accounting = list(
          input_rows = n_obs_int,
          augmented_cols = k_aug,
          augmented_rows = nrow(fac),
          max_chunk_merge_rows = max_merge,
          n_compressions = n_comp
        )
      ),
      class = "cdid_or_pooled_wls_chunk"
    )
  } else {
    raw_sum
  }
  acc <- list(
    L = L_val,
    summary = summary_obj,
    weight_mass = weight_mass,
    n = n_obs_dbl,
    sumZ = res[[4L]],
    sumY = sumY,
    colmax = res[[5L]]
  )
  out <- list(
    design = acc,
    moments = res[[7L]],
    a = a_val,
    pair_mass = weight_mass,
    n_pair = n_obs_dbl
  )
  if (mode <= 2L) out$clipping <- res[[8L]]
  out
}

cdid_fast_sel_if_chunk <- function(mode, backend, reader, frame, ref, rows, links,
                                   scale, pre, link_names = NULL, prep_info = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  req_fns <- switch(
    as.character(mode),
    "1" = c("cp04_pair_read_operation", "cp04_legacy_ipw_native", "cp04_pair_rows",
            "ipw_owner_read_y_positions", "cp04_predict_selection", "cp04_solve", "ipw_logit_values"),
    "2" = c("cp04_pair_consumer_operation", "cp04_pair_read_operation", "selection_prediction_operation",
            "cp04_pair_rows", "cp04_date_rows", "ipw_owner_read_y_positions", "cp04_solve",
            "ipw_logit_values", "chained_reference_x", "chained_reference_s",
            "chained_reference_stat", "chained_reference_weight"),
    "3" = c("cp04_pair_consumer_operation", "cp04_pair_read_operation", "selection_prediction_operation",
            "cp04_pair_rows", "cp04_date_rows", "ipw_owner_read_y_positions",
            "or_outcome_block_fit", "or_selection_separate_if_core", "or_baseline_separate_if_core",
            "or_pooled_wls_gram_solve", "or_wls_gram_solve", "cp04_solve", "ipw_logit_values",
            "chained_reference_x", "chained_reference_s", "chained_reference_stat", "chained_reference_weight"),
    "4" = c("cp04_pair_consumer_operation", "cp04_pair_read_operation", "selection_prediction_operation",
            "cp04_pair_rows", "cp04_date_rows", "ipw_owner_read_y_positions",
            "or_pooled_if_core", "or_pooled_wls_gram_solve", "or_wls_gram_solve",
            "cp04_solve", "ipw_logit_values", "chained_reference_x", "chained_reference_s",
            "chained_reference_stat", "chained_reference_weight"),
    "5" = c("cp04_pair_consumer_operation", "cp04_pair_read_operation", "selection_prediction_operation",
            "cp04_pair_rows", "cp04_date_rows", "ipw_owner_read_y_positions",
            "cp04_predict_selection", "cp04_solve", "or_selection_separate_if_core",
            "or_pooled_if_core", "or_outcome_block_fit", "or_pooled_wls_gram_solve",
            "or_wls_gram_solve", "chained_reference_design", "chained_reference_x",
            "chained_reference_s", "chained_reference_stat", "chained_reference_weight"),
    "6" = c("evaluate_selection_upstream_chunk", "selection_prediction_operation",
            "cp04_pair_rows", "cp04_date_rows", "chained_reference_s",
            "chained_reference_stat", "chained_reference_weight"),
    character(0)
  )
  if (!length(req_fns) || !cdid_are_default_fns(req_fns)) return(NULL)
  fk <- cdid_sel_fast_kernel("C_cdid_sel_if_chunk")
  if (is.null(fk)) return(NULL)

  Za <- if (is.null(ref) && !is.null(frame)) frame$Z_a else NULL
  res <- fk(as.integer(mode), backend, reader, Za, ref, rows, links, pre, prep_info)
  if (is.null(res)) return(NULL)

  phi <- res[[1L]]
  if (is.null(link_names)) {
    link_names <- vapply(backend$links[links], `[[`, character(1), "link_id")
  }
  dimnames(phi) <- list(NULL, link_names)
  st <- reader$state
  if (identical(scale, "core")) phi <- st$lambda * phi
  list(phi = phi, n_act_cells = res[[2L]], n_act_blocks = res[[3L]])
}

cdid_fast_sel_cp04_fit_batch <- function(reader, ss, glm_control, reference = NULL) {
  n_s <- length(ss)
  if (n_s < 1L || !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.list(glm_control) || !cdid_are_default_fns(c(
        "cp04_selection_fit", "cp04_pair_rows", "cp04_logit",
        "ipw_glm_fit_internal", "ipw_logit_values",
        "chained_reference_s", "chained_reference_stat", "chained_reference_weight"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_fit_block")
  if (is.null(fk)) return(NULL)
  gc <- cdid_sel_fast_gc(glm_control)
  if (isTRUE(gc$trace)) return(NULL)
  params <- c(3, 0, gc$epsilon, as.double(gc$maxit))
  cn <- if (is.null(reference)) colnames(reader$state$Z) else
    reference$designs[[reference$role_slot["s"]]]$columns
  dn <- if (!is.null(cn)) list(cn, cn) else NULL
  empty_int <- integer(0L)
  empty_warn <- character(0L)
  out <- vector("list", n_s)
  for (i in seq_len(n_s)) {
    res <- fk(reader, ss[[i]], empty_int, NULL, reference, params)
    if (is.null(res)) return(NULL)
    H_chol <- res[[2L]]
    if (!is.null(dn)) dimnames(H_chol) <- dn
    coef <- res[[1L]]
    out[[i]] <- list(
      coefficients = coef,
      H_chol = H_chol,
      diagnostics = list(
        converged = TRUE,
        boundary = FALSE,
        rank = length(coef),
        iter = res[[3L]],
        warnings = empty_warn
      ),
      pi = res[[4L]],
      r = res[[5L]]
    )
  }
  out
}

cdid_fast_sel_trt_fit_batch_cp04 <- function(reader, ts, selections, glm_control) {
  n_t <- length(ts)
  if (n_t < 1L || !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.list(glm_control) || !cdid_are_default_fns(c(
        "cp04_treatment_fit", "cp04_pair_rows", "cp04_predict_selection",
        "cp04_logit", "ipw_glm_fit_internal", "ipw_logit_values", "ipw_clip_summary"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_trt_fit")
  if (is.null(fk)) return(NULL)
  gc <- cdid_sel_fast_gc(glm_control)
  if (isTRUE(gc$trace)) return(NULL)
  frame <- get0("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(frame)) return(NULL)
  cols_info <- list(Z_a = frame$Z_a)
  ctrl_vec <- c(as.double(gc$epsilon), as.double(gc$maxit))
  rn4 <- c("T", "M", "dp", "dm")
  dn4 <- list(rn4, NULL)
  empty_warn <- character(0L)
  out <- vector("list", n_t)
  for (i in seq_len(n_t)) {
    res <- fk(1L, reader, ts[[i]], selections, NULL, NULL, NULL, cols_info, ctrl_vec)
    if (is.null(res)) return(NULL)
    beta <- res[[1L]]
    U_main <- res[[5L]]
    V_main <- res[[6L]]
    dimnames(U_main) <- dn4
    dimnames(V_main) <- dn4
    means <- res[[7L]]
    names(means) <- rn4
    out[[i]] <- list(
      coefficients = beta,
      H_chol = res[[2L]],
      diagnostics = list(
        converged = TRUE,
        boundary = FALSE,
        rank = length(beta),
        iter = res[[3L]],
        warnings = empty_warn
      ),
      cross = res[[4L]],
      partials = list(upstream = U_main, beta = V_main),
      means = means,
      clipping = res[[8L]]
    )
  }
  out
}

cdid_fast_sel_cor1_support_batch <- function(reader, specs, selections, recipe, ctrl) {
  n_s <- length(specs)
  if (n_s < 1L || !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "or_corollary1_support", "or_baseline_visit", "cp04_pair_rows",
        "cp04_pair_read_operation", "ipw_owner_read_y_positions",
        "cp04_predict_selection", "or_pooled_support_add",
        "or_pooled_wls_block_summary", "or_pooled_wls_merge_chunks",
        "or_wls_summary", "or_wls_merge", "or_wls_compress",
        "or_corollary1_moments_add", "ipw_clip_summary", "ipw_clip_summary_merge"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_support_collect")
  if (is.null(fk)) return(NULL)
  st <- reader$state
  frame <- get0("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(frame)) return(NULL)
  icpt_flag <- isTRUE(recipe$intercept)
  p <- as.integer(recipe$p)
  slope_names <- if (p > 0L) recipe$slope_columns else NULL
  raw_cols <- if (icpt_flag) c("", slope_names) else slope_names
  k_aug <- p + 1L + as.integer(icpt_flag)
  slope_idx <- as.integer(recipe$slope_index)
  icpt_int <- as.integer(icpt_flag)
  icpt_col <- if (icpt_flag) as.integer(recipe$intercept_index) else 0L
  max_chunk <- as.integer(ctrl$max_chunk_rows)
  Za <- frame$Z_a
  a_vec <- st$a
  out <- vector("list", n_s)
  for (j in seq_len(n_s)) {
    spec <- specs[[j]]
    L_val <- as.double(spec$transitions)
    cols_info <- list(
      Z_a = Za, z_idx = NULL, slope_idx = slope_idx,
      icpt_flag = icpt_int, icpt_col = icpt_col,
      L_val = L_val, max_chunk_rows = max_chunk
    )
    res <- fk(1L, reader, spec, selections[[j]], NULL, NULL, NULL, cols_info)
    if (is.null(res)) return(NULL)
    fac <- res[[1L]]
    sc <- res[[6L]]
    n_obs_int <- as.integer(sc[1L])
    n_obs_dbl <- as.double(sc[1L])
    weight_mass <- sc[4L]
    raw_sum <- structure(
      list(factor = fac, n = n_obs_int, columns = raw_cols),
      class = "cdid_or_wls_summary"
    )
    summary_obj <- if (icpt_flag) {
      structure(
        list(
          raw = raw_sum, n = n_obs_int, p = p, columns = slope_names,
          anchor = res[[2L]], constant = res[[3L]],
          accounting = list(
            input_rows = n_obs_int, augmented_cols = k_aug,
            augmented_rows = nrow(fac),
            max_chunk_merge_rows = as.integer(sc[2L]),
            n_compressions = as.integer(sc[3L])
          )
        ),
        class = "cdid_or_pooled_wls_chunk"
      )
    } else {
      raw_sum
    }
    out[[j]] <- list(
      design = list(
        L = L_val, summary = summary_obj, weight_mass = weight_mass,
        n = n_obs_dbl, sumZ = res[[4L]], sumY = sc[5L], colmax = res[[5L]]
      ),
      moments = res[[7L]],
      a = a_vec[spec$cohort],
      pair_mass = weight_mass,
      n_pair = n_obs_dbl,
      clipping = res[[8L]]
    )
  }
  out
}

cdid_sel_fast_or_pooled_fit <- function(recipe, supports, blocks, rank_tol, N, fk_comp = NULL) {
  if (is.null(fk_comp)) {
    if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
        !is.list(recipe) || !isTRUE(recipe$intercept) ||
        !gmm_scalar_numeric(rank_tol) || rank_tol <= 0 || rank_tol >= 1 ||
        !gmm_scalar_numeric(N) || N <= 0 || N != floor(N) ||
        !is.list(supports) || !length(supports) ||
        !is.list(blocks) || !length(blocks) ||
        !cdid_are_default_fns(c(
          "or_pooled_fit_blocks", "or_pooled_wls_fit", "or_pooled_wls_merge_chunks",
          "or_pooled_wls_merge", "or_pooled_wls_partial_block", "or_wls_fit",
          "or_wls_compress", "or_wls_merge", "or_public_basis"
        ))) {
      return(NULL)
    }
    fk_comp <- cdid_sel_fast_kernel("C_cdid_wls_compress")
    if (is.null(fk_comp)) return(NULL)
  }
  p <- as.integer(recipe$p)
  n_b <- length(blocks)
  block_names <- vapply(blocks, `[[`, character(1L), "name")
  if (anyNA(block_names) || any(!nzchar(block_names)) || anyDuplicated(block_names)) return(NULL)
  block_recs <- vector("list", n_b)
  acc_fac <- NULL
  acc_n <- 0L
  max_aug <- 0L
  max_res <- 0L
  max_chunk_merge <- 0L
  max_res_merge <- 0L
  n_compressions <- 0L
  cols <- if (p > 0L) recipe$slope_columns else NULL
  n_sup <- length(supports)
  for (b in seq_len(n_b)) {
    bb <- blocks[[b]]
    sids <- bb$support_ids
    if (!length(sids) || any(sids < 1L | sids > n_sup) || anyDuplicated(sids)) return(NULL)
    sup1 <- supports[[sids[1L]]]
    chunk1 <- if (!is.null(sup1$design)) sup1$design$summary else {
      if (!identical(sup1$scale, 1.0)) return(NULL)
      sup1$summary
    }
    if (!inherits(chunk1, "cdid_or_pooled_wls_chunk")) return(NULL)
    fac <- chunk1$raw$factor
    c_n <- as.integer(chunk1$n)
    c_const <- chunk1$constant
    c_anchor <- chunk1$anchor
    c_merge <- chunk1$accounting$max_chunk_merge_rows
    c_comp <- chunk1$accounting$n_compressions
    if (length(sids) > 1L) {
      for (sid in sids[-1L]) {
        sup_k <- supports[[sid]]
        nxt <- if (!is.null(sup_k$design)) sup_k$design$summary else {
          if (!identical(sup_k$scale, 1.0)) return(NULL)
          sup_k$summary
        }
        if (!inherits(nxt, "cdid_or_pooled_wls_chunk")) return(NULL)
        st_rows <- nrow(fac) + nrow(nxt$raw$factor)
        if (st_rows > c_merge) c_merge <- st_rows
        if (nxt$n > 0L) {
          if (c_n == 0L) {
            c_anchor <- nxt$anchor
            c_const <- nxt$constant
          } else {
            c_const <- c_const & nxt$constant & (c_anchor == nxt$anchor)
          }
        }
        st_mat <- rbind(fac, nxt$raw$factor)
        fac <- if (nrow(st_mat) > 0L) fk_comp(st_mat) else st_mat
        if (is.null(fac)) return(NULL)
        c_n <- c_n + as.integer(nxt$n)
        c_comp <- c_comp + nxt$accounting$n_compressions + if (st_rows > 0L) 1L else 0L
      }
    }
    m_fac <- nrow(fac)
    if (c_n < 1L || m_fac < 1L || ncol(fac) != p + 2L || any(!is.finite(fac))) return(NULL)
    q_int <- qr.default(fac[, 1L, drop = FALSE], LAPACK = TRUE)
    r11 <- qr.R(q_int)[1L, 1L]
    if (!is.finite(r11) || abs(r11) <= 0.0) return(NULL)
    qty <- qr.qty(q_int, fac[, -1L, drop = FALSE])
    if (any(!is.finite(qty))) return(NULL)
    r1d <- if (p > 0L) qty[1L, seq_len(p)] else numeric(0L)
    names(r1d) <- cols
    r1y <- qty[1L, p + 1L]
    res_raw <- qty[-1L, , drop = FALSE]
    if (p > 0L && any(c_const)) {
      res_raw[, which(c_const)] <- 0.0
    }
    res_fac <- if (nrow(res_raw) > 0L) fk_comp(res_raw) else res_raw
    if (is.null(res_fac)) return(NULL)
    m_res <- nrow(res_fac)
    block_recs[[b]] <- list(
      r11 = r11, r1d = r1d, r1y = r1y, n = c_n,
      weight_sum = r11 * r11, mean_d = r1d / r11, mean_y = r1y / r11,
      augmented_rows = m_fac, residual_rows = m_res
    )
    c_comp <- c_comp + if (nrow(res_raw) > 0L) 1L else 0L
    if (b == 1L) {
      acc_fac <- res_fac
      acc_n <- c_n
      max_aug <- m_fac
      max_res <- m_res
      max_chunk_merge <- c_merge
      max_res_merge <- m_res
      n_compressions <- c_comp
    } else {
      stacked <- nrow(acc_fac) + m_res
      if (stacked > max_res_merge) max_res_merge <- stacked
      stacked_mat <- rbind(acc_fac, res_fac)
      acc_fac <- if (nrow(stacked_mat) > 0L) fk_comp(stacked_mat) else stacked_mat
      if (is.null(acc_fac)) return(NULL)
      acc_n <- acc_n + c_n
      if (m_fac > max_aug) max_aug <- m_fac
      if (m_res > max_res) max_res <- m_res
      if (c_merge > max_chunk_merge) max_chunk_merge <- c_merge
      n_compressions <- n_compressions + c_comp + if (stacked > 0L) 1L else 0L
    }
  }
  names(block_recs) <- block_names
  keep <- seq_len(p)
  if (p > 0L) {
    if (nrow(acc_fac) < p) return(NULL)
    X <- acc_fac[, keep, drop = FALSE]
    proposal <- qr.default(X, tol = rank_tol, LAPACK = FALSE)
    if (proposal$rank < p) return(NULL)
    q <- qr.default(X, LAPACK = TRUE)
    R <- qr.R(q)
    diagonal <- abs(diag(R))
    cutoff <- rank_tol * max(diagonal)
    rank <- if (cutoff > 0.0) sum(diagonal > cutoff) else 0L
    if (rank < p) return(NULL)
    qty <- qr.qty(q, acc_fac[, p + 1L])
    coefficients <- numeric(p)
    coefficients[q$pivot] <- backsolve(R, qty[keep])
    names(coefficients) <- cols
    rss <- sum(qty[-keep]^2)
    R_out <- unname(R)
    pivot <- q$pivot
    if (any(!is.finite(c(coefficients, rss, R_out)))) return(NULL)
    slope_fit <- structure(
      list(
        coefficients = coefficients, rss = rss, R = R_out, pivot = pivot,
        n = acc_n, rank = rank, rank_tol = rank_tol, cutoff = cutoff
      ),
      class = "cdid_or_wls_fit"
    )
  } else {
    slope_fit <- NULL
    coefficients <- numeric(0L)
    if (!is.null(cols)) names(coefficients) <- cols
    rss <- sum(acc_fac[, 1L]^2)
    R_out <- matrix(numeric(0L), 0L, 0L)
    pivot <- integer(0L)
    rank <- 0L
    cutoff <- 0.0
    if (!is.finite(rss)) return(NULL)
  }
  intercepts <- numeric(n_b)
  regressions <- vector("list", n_b)
  beta <- unname(coefficients)
  for (b in seq_len(n_b)) {
    rec <- block_recs[[b]]
    rhs <- rec$r1y - sum(rec$r1d * coefficients)
    alpha_b <- as.numeric(backsolve(matrix(rec$r11, 1L, 1L), rhs))
    if (!is.finite(alpha_b)) return(NULL)
    intercepts[b] <- alpha_b
    regressions[[b]] <- list(
      name = block_names[b], L = blocks[[b]]$L, n = rec$n, alpha = alpha_b,
      weight_sum = rec$weight_sum, mass = rec$weight_sum / N,
      mean_d = unname(rec$mean_d), mean_y = rec$mean_y,
      residual_rows = rec$residual_rows
    )
  }
  names(intercepts) <- block_names
  acct <- list(
    n = acc_n, p = p, n_blocks = n_b, residual_dim = dim(acc_fac),
    max_augmented_rows = max_aug, max_residual_rows = max_res,
    max_chunk_merge_rows = max_chunk_merge,
    max_residual_merge_rows = max_res_merge,
    n_compressions = n_compressions
  )
  res_obj <- if (p > 0L) {
    structure(list(factor = acc_fac, n = acc_n, columns = cols), class = "cdid_or_wls_summary")
  } else {
    list(factor = acc_fac, n = acc_n, columns = cols)
  }
  pooled <- structure(
    list(
      residual = res_obj, blocks = block_recs, n = acc_n, p = p,
      columns = cols, accounting = acct
    ),
    class = "cdid_or_pooled_wls_summary"
  )
  fit <- structure(
    list(
      coefficients = coefficients, intercepts = intercepts, rss = rss,
      R = R_out, pivot = pivot, n = acc_n, p = p, n_blocks = n_b,
      rank = rank, rank_tol = rank_tol, cutoff = cutoff,
      blocks = block_recs, slope_fit = slope_fit, accounting = acct
    ),
    class = "cdid_or_pooled_wls_fit"
  )
  map <- matrix(0.0, p, 0L, dimnames = list(cols, character(0L)))
  constants <- matrix(0.0, n_b, 0L, dimnames = list(block_names, character(0L)))
  basis <- list(
    slope_columns = recipe$slope_columns,
    kept_columns = recipe$slope_columns[keep],
    kept_index = keep,
    proposed_columns = character(0L),
    aliased_columns = character(0L),
    nontrivial_candidates_retained = character(0L),
    alias_map = map,
    block_constants = constants,
    reconstruction_control_scaled = numeric(0L),
    reconstruction_treated_scaled = numeric(0L),
    reconstruction_roundoff_threshold = NA_real_,
    rows_checked = 0,
    proposal = "R stats::lm.fit / dqrdc2 pivoted QR at rank_tol on the pooled (partialled) factor",
    rank_check = "unchanged relative-diagonal rule of the accepted kernel on the reduced factor",
    frozen_during_differentiation = TRUE,
    numerical_certificate = FALSE
  )
  list(
    fit = fit, summary = pooled, keep = keep, beta = beta,
    regressions = regressions, basis = basis, rss = rss,
    objective = rss / N, verification_scan = FALSE
  )
}

cdid_sel_fast_or_outcome_fit_blocks <- function(recipe, supports, blocks, links,
                                                rank_tol, N, mode, route = "pooled") {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.list(recipe) || !isTRUE(recipe$intercept) ||
      !gmm_scalar_numeric(rank_tol) || rank_tol <= 0 || rank_tol >= 1 ||
      !gmm_scalar_numeric(N) || N <= 0 || N != floor(N) ||
      !is.character(mode) || length(mode) != 1L || !mode %in% c("pooled", "separate") ||
      !is.list(supports) || !length(supports) ||
      !is.list(blocks) || !length(blocks) || !is.list(links) ||
      !cdid_are_default_fns(c(
        "or_outcome_fit_blocks", "or_outcome_group_layout", "or_pooled_fit_blocks",
        "or_pooled_wls_fit", "or_pooled_wls_merge_chunks", "or_pooled_wls_merge",
        "or_pooled_wls_partial_block", "or_wls_fit", "or_wls_compress",
        "or_wls_merge", "or_public_basis"
      ))) {
    return(NULL)
  }
  fk_comp <- cdid_sel_fast_kernel("C_cdid_wls_compress")
  if (is.null(fk_comp)) return(NULL)
  B <- length(blocks)
  n_sup <- length(supports)
  block_names <- vapply(blocks, `[[`, character(1L), "name")
  if (anyNA(block_names) || any(!nzchar(block_names)) || anyDuplicated(block_names)) return(NULL)
  bids <- vapply(links, function(l) as.integer(l$regression), integer(1L))
  lsup <- vapply(links, function(l) as.integer(l$support_id), integer(1L))
  if (anyNA(bids) || any(bids < 1L | bids > B) || anyNA(lsup) || any(lsup < 1L | lsup > n_sup)) {
    return(NULL)
  }
  if (mode == "pooled") {
    g1 <- cdid_sel_fast_or_pooled_fit(recipe, supports, blocks, rank_tol, N, fk_comp = fk_comp)
    if (is.null(g1)) return(NULL)
    groups <- list(g1)
    group_id <- rep.int(1L, B)
    local_block_id <- seq_len(B)
    lids <- seq_along(links)
    names(lids) <- names(links)
    group_map <- list(list(
      block_ids = seq_len(B),
      support_ids = seq_len(n_sup),
      link_ids = lids,
      slope_dimension_bound = recipe$p,
      coefficient_dimension_bound = recipe$p + B,
      operator_dimension_bound = c(recipe$p, recipe$p),
      summary_dimension_bound = c(recipe$p + 1L, recipe$p + 1L)
    ))
  } else {
    groups <- vector("list", B)
    group_map <- vector("list", B)
    l_names <- names(links)
    for (g in seq_len(B)) {
      lids <- which(bids == g)
      if (!is.null(l_names)) names(lids) <- l_names[lids]
      b_g <- blocks[[g]]
      sids <- sort.int(unique(c(as.integer(b_g$support_ids), lsup[lids])))
      local_b <- b_g
      local_b$support_ids <- match(b_g$support_ids, sids)
      gg <- cdid_sel_fast_or_pooled_fit(recipe, supports[sids], list(local_b), rank_tol, N, fk_comp = fk_comp)
      if (is.null(gg)) return(NULL)
      groups[[g]] <- gg
      group_map[[g]] <- list(
        block_ids = g,
        support_ids = sids,
        link_ids = lids,
        slope_dimension_bound = recipe$p,
        coefficient_dimension_bound = recipe$p + 1L,
        operator_dimension_bound = c(recipe$p, recipe$p),
        summary_dimension_bound = c(recipe$p + 1L, recipe$p + 1L)
      )
    }
    group_id <- seq_len(B)
    local_block_id <- rep.int(1L, B)
  }
  block_map <- structure(
    list(
      block_id = seq_len(B),
      name = block_names,
      group_id = group_id,
      local_block_id = local_block_id
    ),
    class = "data.frame",
    row.names = .set_row_names(B)
  )
  rss <- sum(vapply(groups, `[[`, numeric(1L), "rss"))
  retained <- lapply(groups, function(g) list(
    slope_dimension = length(g$keep),
    operator_dimension = dim(g$fit$R),
    summary_dimension = dim(g$summary$residual$factor),
    alias_map_dimension = dim(g$basis$alias_map),
    block_constants_dimension = dim(g$basis$block_constants),
    verification_scan = g$verification_scan,
    alias_rows_checked = g$basis$rows_checked
  ))
  structure(
    list(
      mode = mode,
      recipe = recipe,
      N = N,
      route = route,
      block_map = block_map,
      group_map = group_map,
      groups = groups,
      rss = rss,
      objective = rss / N,
      retained = retained
    ),
    class = "cdid_or_outcome_groups"
  )
}

cdid_fast_or_sep_envelope <- function(backend, B, links) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.null(attributes(links)) || !is.numeric(B) || length(B) != 1L ||
      !is.finite(B) || B < 1 || B > 1e7 ||
      !cdid_are_default_fns(c(
        "or_selection_separate_if_envelope",
        "selection_if_budget",
        "cp04_budget",
        "ipw_owner_safe_sum",
        "ipw_owner_safe_product"
      ))) {
    return(NULL)
  }
  links_int <- as.integer(links)
  bref <- .subset2(backend, "ref")
  slots <- get0("or_sep_env_static", envir = cdid_sel_cache, inherits = FALSE)
  rec <- NULL
  if (is.list(slots)) {
    for (i in seq_along(slots)) {
      hit <- slots[[i]]
      if (is.list(hit) && identical(hit$backend_ref, bref) &&
          gmm_sexp_aliases(hit$backend, backend) && identical(hit$links, links_int)) {
        rec <- hit
        break
      }
    }
  }
  if (is.null(rec)) {
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    G <- length(backend$regression_specs); S <- length(backend$supports)
    p <- backend$outcome_recipe$p; D <- backend$upstream$D; d <- backend$upstream$d
    J <- length(st$cohorts); nl <- length(links_int); incidence <- 0
    for (r in backend$regression_specs) incidence <- incidence + length(r$support_ids)
    for (u in backend$supports) incidence <- incidence + length(u$channels)
    maps <- 4096 + 32 * (D + S + G + nl + incidence)
    cp04_budget(
      backend$control,
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      maps, c(backend$accounting$work_ops_bound, 64 * maps), 0,
      "separate conditional OR dependency metadata"
    )
    blinks <- backend$links[links_int]
    gids <- sort.int(unique(vapply(blinks, `[[`, integer(1L), "regression")))
    sids <- sort.int(unique(c(
      vapply(blinks, `[[`, integer(1L), "support_id"),
      unlist(lapply(backend$regression_specs[gids], `[[`, "support_ids"), use.names = FALSE)
    )))
    uids <- sort.int(unique(unlist(lapply(backend$supports[sids], `[[`, "channels"), use.names = FALSE)))
    ns <- length(sids); nu <- length(uids)
    ps <- vapply(backend$outcome_groups$groups[gids], function(g) length(g$keep), integer(1L))
    mp <- max(c(0L, ps))
    link_bytes <- sum(nchar(vapply(blinks, `[[`, character(1L), "link_id"), type = "bytes"))
    id_max_bytes <- frame$id_byte_metadata$max
    meta_0 <- 8192 + 16 * link_bytes + 128 * nl
    width <- nl + 8
    scr_0 <- maps + 128 * ((mp + 2)^2 + (d + 1)^2)
    scr_B <- 32 * (ns * (p + 6) + nu + mp * (mp + nu + 6) + d + nl + 16)
    work_0 <- 128 * maps
    work_B <- 128 * (
      (ns + 1) * ((p + 2)^2 + (d + 1)^2) * (1 + ceiling(log2(st$n_A + 1))) +
        sum((ps + 2)^2 * (nu + nl + 2)) + J + nl * (p + nu + 2)
    )
    uses_sel <- isTRUE(backend$uses_selection_models)
    sp_info <- NULL
    if (uses_sel) {
      sel <- backend$selection_backend
      sel_q_ids <- unique(backend$q_of_s[sids])
      sp_active_blocks <- unique(unlist(lapply(sel$pair_cohorts[sel_q_ids], `[[`, "block_ids"), use.names = FALSE))
      sp_active_cells <- unique(unlist(lapply(sel$block_specs[sp_active_blocks], `[[`, "cell_ids"), use.names = FALSE))
      n_theta <- length(sel$coefficients)
      n_probes <- length(sel_q_ids)
      sp_widths <- n_theta + J + 3 * n_probes + length(sp_active_cells) + 32
      sp_meta_0 <- 8192 + 1024 * (length(sel$blocks) + length(sel$cells) + J) + n_probes * (sel$design_spec$d_s + 1024)
      sp_retained <- sel$accounting$owned_entries + sel$accounting$borrowed_entries
      sp_per_row_work <- 32 * length(sp_active_cells) * (1 + ceiling(log2(st$N + 1))) +
        32 * sum(vapply(sel$blocks[sp_active_blocks], function(bl) bl$B_m * (bl$d_s + 2)^2, numeric(1L))) +
        32 * (n_theta + J + 2 * n_probes + 8)
      sp_info <- list(
        ctrl = sel$control, retained = sp_retained,
        meta_0 = sp_meta_0, widths = sp_widths,
        d_raw = sel$design_spec$d, d_s = sel$design_spec$d_s,
        per_row_work = sp_per_row_work,
        target_scan = 32 * st$N * (sel$design_spec$d + 4)
      )
    }
    rec <- list(
      backend = backend, backend_ref = bref, links = links_int,
      gids = gids, sids = sids, uids = uids,
      id_max_bytes = id_max_bytes, meta_0 = meta_0, width = width,
      scr_0 = scr_0, scr_B = scr_B, work_0 = work_0, work_B = work_B,
      uses_sel = uses_sel, sp = sp_info
    )
    slots <- if (!is.list(slots) || length(slots) != 2L) list(rec, NULL) else list(rec, slots[[1L]])
    assign("or_sep_env_static", slots, envir = cdid_sel_cache)
  } else {
    cp04_budget(
      backend$control,
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      rec$work_0 / 128, c(backend$accounting$work_ops_bound, 64 * (rec$work_0 / 128)), 0,
      "separate conditional OR dependency metadata"
    )
  }
  B_dbl <- as.double(B)
  metadata <- rec$meta_0 + B_dbl * rec$id_max_bytes
  output <- B_dbl * rec$width + metadata
  scratch <- rec$scr_0 + output + B_dbl * rec$scr_B
  work <- rec$work_0 + B_dbl * rec$work_B
  if (rec$uses_sel) {
    sp <- rec$sp
    sp_meta <- sp$meta_0 + B_dbl * rec$id_max_bytes
    sp_output <- B_dbl * (sp$widths + 16) + sp_meta
    sp_scratch <- 32 * B_dbl * (sp$widths + sp$d_raw + 8) + 64 * (sp$d_s + 2)^2 + sp_meta
    sp_work <- B_dbl * sp$per_row_work + sp$target_scan + 64 * sp_meta
    cp04_budget(sp$ctrl, sp$retained, sp_scratch, sp_work, sp_output, "selection IF chunk")
    scratch <- scratch + sp_output + sp_scratch
    work <- work + sp_work
  }
  if (!is.finite(scratch) || !is.finite(work) || !is.finite(output) ||
      scratch > 1e15 || work > 1e15 || output > 1e15) {
    return(NULL)
  }
  list(
    gids = rec$gids, beta_gids = rec$gids, sids = rec$sids, uids = rec$uids,
    output = output, scratch = scratch, work = work
  )
}

cdid_fast_eval_or_cor1_if_chunk <- function(backend, rows, links, scale, details) {
  if (!identical(details, FALSE) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !(scale == "frame" || scale == "core") ||
      !inherits(backend, "cdid_or_corollary1") || isTRUE(.subset2(backend, "uses_selection_models")) ||
      !is.null(.subset2(backend, "reference_ref")) || !is.null(.subset2(backend, "layout")) ||
      !is.null(.subset2(backend, "pair_outcomes")) ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "validate_or_corollary1_backend", "ipw_owner_seal_ok",
        "validate_ipw_frame_owner", "validate_ipw_panel_owner",
        "validate_cp04_reader", "validate_chained_reference_reader",
        "validate_ipw_owner_control", "validate_ipw_glm_control",
        "validate_selection_models", "or_corollary1_if_budget",
        "or_corollary1_if_envelope", "or_selection_separate_if_envelope",
        "cp04_budget", "cp04_payload",
        "cp04_pair_consumer_operation", "cp04_pair_read_operation",
        "selection_prediction_operation", "cp04_pair_rows", "cp04_date_rows",
        "ipw_owner_read_y_positions", "cp04_predict_selection", "cp04_solve",
        "or_selection_separate_if_core", "or_pooled_if_core",
        "or_outcome_block_fit", "or_pooled_wls_gram_solve", "or_wls_gram_solve",
        "chained_reference_design", "chained_reference_x",
        "chained_reference_s", "chained_reference_stat", "chained_reference_weight"
      ))) {
    return(NULL)
  }
  ring <- get0("val_ring", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(ring) || !length(ring)) return(NULL)
  fk_vf <- cdid_sel_fast_kernel("C_cdid_verify_fast")
  fk_if <- cdid_sel_fast_kernel("C_cdid_sel_if_chunk")
  if (is.null(fk_vf) || is.null(fk_if)) return(NULL)
  hit_ok <- FALSE
  for (i in seq_along(ring)) {
    rec <- ring[[i]]
    if (is.list(rec) && identical(rec[[2L]], "or_corollary1_backend") &&
        isTRUE(fk_vf(4L, backend, rec, NULL, NULL))) {
      hit_ok <- TRUE
      break
    }
  }
  if (!hit_ok) return(NULL)
  reader <- get("cp04_reader", .subset2(backend, "reader_ref"), inherits = FALSE)
  st <- reader$state
  B <- length(rows)
  if (!cdid_sel_valid_idx(rows, st$N) || B > backend$control$max_chunk_rows ||
      !cdid_sel_valid_idx(links, length(backend$links))) {
    return(NULL)
  }
  meta <- cdid_sel_get_if_meta(backend, B, links, "or_cor1_if_meta")
  if (is.null(meta) || is.null(meta$pre$diag_tmpl)) {
    pre <- or_corollary1_if_budget(backend, B, links, FALSE)
    separate <- identical(backend$outcome_pooling, "separate")
    diag_tmpl <- list(
      n_rate_score_evaluations = sum(pre$uids <= backend$upstream$n_a),
      n_selection_score_evaluations = length(pre$sids),
      n_upstream_coordinates = length(pre$uids),
      n_regression_score_evaluations = length(union(pre$gids, pre$beta_gids)),
      n_global_beta_blocks = if (separate) 0L else length(pre$beta_gids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
    if (separate) {
      diag_tmpl$outcome_pooling <- "separate"
      diag_tmpl$n_local_beta_blocks <- length(pre$beta_gids)
    }
    pre$diag_tmpl <- diag_tmpl
    pre$pay_env <- new.env(parent = emptyenv())
    link_names <- vapply(backend$links[links], `[[`, character(1L), "link_id")
    meta <- cdid_sel_set_if_meta(backend, B, links, pre, link_names, "or_cor1_if_meta")
  } else {
    pre <- meta$pre
    link_names <- meta$link_names
  }
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  res <- fk_if(5L, backend, reader, frame$Z_a, NULL, rows, links, pre, NULL)
  if (is.null(res)) return(NULL)
  phi <- res[[1L]]
  dimnames(phi) <- list(NULL, link_names)
  if (scale == "core") phi <- st$lambda * phi
  out <- list(
    rows = rows,
    id = cdid_sel_slice_id(st$id, rows, st$N),
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = pre$diag_tmpl
  )
  can_skip_pay <- !is.character(st$id) && is.null(attributes(rows)) &&
    is.null(attributes(links)) && is.null(attributes(st$id))
  if (can_skip_pay && isTRUE(get0(scale, envir = pre$pay_env, inherits = FALSE))) {
    return(out)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output) {
    or_corollary1_abort("IF output bound underestimated.", "bound")
  }
  if (can_skip_pay) assign(scale, TRUE, envir = pre$pay_env)
  out
}

cdid_fast_eval_cp04_cor1_if_chunk <- function(backend, rows, links, scale, details) {
  if (!identical(details, FALSE) || !is.character(scale) || length(scale) != 1L ||
      is.na(scale) || !(scale == "frame" || scale == "core") ||
      !inherits(backend, "cdid_cp04_corollary1") || isTRUE(.subset2(backend, "uses_selection_models")) ||
      !is.null(.subset2(backend, "reference_ref")) || !is.null(.subset2(backend, "layout")) ||
      !is.null(.subset2(backend, "pair_outcomes")) ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "validate_cp04_backend", "ipw_owner_seal_ok",
        "validate_ipw_frame_owner", "validate_ipw_panel_owner",
        "validate_cp04_reader", "validate_chained_reference_reader",
        "validate_ipw_owner_control", "validate_ipw_glm_control",
        "validate_selection_models", "cp04_if_budget", "cp04_budget", "cp04_payload",
        "cp04_pair_read_operation", "cp04_legacy_ipw_native", "cp04_pair_rows",
        "ipw_owner_read_y_positions", "cp04_predict_selection", "cp04_solve", "ipw_logit_values"
      ))) {
    return(NULL)
  }
  ring <- get0("val_ring", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(ring) || !length(ring)) return(NULL)
  fk_vf <- cdid_sel_fast_kernel("C_cdid_verify_fast")
  fk_if <- cdid_sel_fast_kernel("C_cdid_sel_if_chunk")
  if (is.null(fk_vf) || is.null(fk_if)) return(NULL)
  hit_ok <- FALSE
  for (i in seq_along(ring)) {
    rec <- ring[[i]]
    if (is.list(rec) && identical(rec[[2L]], "cp04_backend") &&
        isTRUE(fk_vf(4L, backend, rec, NULL, NULL))) {
      hit_ok <- TRUE
      break
    }
  }
  if (!hit_ok) return(NULL)
  reader <- get("cp04_reader", .subset2(backend, "reader_ref"), inherits = FALSE)
  st <- reader$state
  B <- length(rows)
  if (!cdid_sel_valid_idx(rows, st$N) || B > backend$control$max_chunk_rows ||
      !cdid_sel_valid_idx(links, length(backend$links))) {
    return(NULL)
  }
  meta <- cdid_sel_get_if_meta(backend, B, links, "cp04_if_meta")
  if (is.null(meta) || is.null(meta$pre$diag_tmpl)) {
    pre <- cp04_if_budget(backend, B, links)
    pre$diag_tmpl <- list(
      n_selection_score_evaluations = length(pre$sids),
      n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
    pre$pay_env <- new.env(parent = emptyenv())
    link_names <- vapply(backend$links[links], `[[`, character(1L), "link_id")
    meta <- cdid_sel_set_if_meta(backend, B, links, pre, link_names, "cp04_if_meta")
  } else {
    pre <- meta$pre
    link_names <- meta$link_names
  }
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  res <- fk_if(1L, backend, reader, frame$Z_a, NULL, rows, links, pre, NULL)
  if (is.null(res)) return(NULL)
  phi <- res[[1L]]
  dimnames(phi) <- list(NULL, link_names)
  if (scale == "core") phi <- st$lambda * phi
  out <- list(
    rows = rows,
    id = cdid_sel_slice_id(st$id, rows, st$N),
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = pre$diag_tmpl
  )
  can_skip_pay <- !is.character(st$id) && is.null(attributes(rows)) &&
    is.null(attributes(links)) && is.null(attributes(st$id))
  if (can_skip_pay && isTRUE(get0(scale, envir = pre$pay_env, inherits = FALSE))) {
    return(out)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output) {
    cp04_abort("IF output bound underestimated.", "bound")
  }
  if (can_skip_pay) assign(scale, TRUE, envir = pre$pay_env)
  out
}

cdid_fast_or_cor1_catalog <- function(link_views, st, frame, ctrl, base, header,
                                      external_work_ops, recipe, separate,
                                      route = "OR", selection_model = NULL,
                                      selection_pooling = NULL) {
  L <- length(link_views)
  J <- length(st$cohorts)
  if (L < 1L || L > 16384L || J < 1L || J > 16384L ||
      !st$universe %in% c("F", "B", "A") ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "validate_ipw_link_view", "cp04_budget", "cp04_payload",
        "or_baseline_payload_append", "or_baseline_payload_verify",
        "or_pooled_transitions", "or_baseline_catalog",
        "or_baseline_catalog_key", "ipw_owner_safe_sum"
      ))) {
    return(NULL)
  }
  if (!isTRUE(cdid_fast_validate_link_views(link_views, frame))) return(NULL)
  cohorts_vec <- st$cohorts
  fk_prefix <- if (identical(route, "DR")) {
    paste0(
      "cmpcat:DR-Cor1-", st$universe, "-", selection_model, "-", selection_pooling,
      if (separate) "-separate-slopes" else "-shared-slope",
      "-WLS-IPW-a/(rho*p)|"
    )
  } else {
    paste0(
      "cmpcat:OR-Cor1-", st$universe,
      if (separate) "-separate-slopes" else "-shared-slope",
      "-WLS-a/(rho*p)|"
    )
  }
  pairs <- list()
  pair_keys <- integer(0L)
  p_ds <- numeric(0L)
  p_dt <- numeric(0L)
  pools <- list()
  pool_keys <- integer(0L)
  links <- vector("list", L)
  ids <- character(L)
  sum_kb <- 0
  max_kb <- 0
  sum_idb <- 0
  max_clen <- 0L
  for (ell in seq_len(L)) {
    v <- link_views[[ell]]
    source <- v$functional_key
    if (!gmm_scalar_character(source) || !startsWith(source, "cmpcat:")) return(NULL)
    kb <- nchar(source, type = "bytes")
    sum_kb <- sum_kb + kb
    if (kb > max_kb) max_kb <- kb
    if (!is.null(v$candidates) || !is.null(v$loading_beta_candidate)) return(NULL)
    if (!is.null(v$nuisance_ref)) {
      nv <- get0("nuisance_view", v$nuisance_ref, inherits = FALSE)
      if (!is.list(nv) || !isTRUE(v$nuisance_mask_verified) || !identical(nv$defining_pair, v$pair)) {
        return(NULL)
      }
    }
    lid <- v$link_id
    if (!gmm_scalar_character(lid) || !nzchar(lid)) return(NULL)
    ids[ell] <- lid
    sum_idb <- sum_idb + nchar(lid, type = "bytes")
    target <- match(v$cohorts$g, cohorts_vec)
    Kc <- v$cohorts$K
    if (length(Kc) > 1L) Kc <- sort(Kc, method = "radix")
    control_cohorts <- match(Kc, cohorts_vec)
    clen <- length(control_cohorts)
    if (length(target) != 1L || is.na(target) || clen < 1L || anyNA(control_cohorts) ||
        !isTRUE(all(v$support$by_cohort_pair_count >= 1))) {
      return(NULL)
    }
    pr <- v$pair
    ds <- pr$s_calendar_index
    dt <- pr$t_calendar_index
    if (!is.numeric(ds) || !is.numeric(dt) || length(ds) != 1L || length(dt) != 1L ||
        !is.finite(ds) || !is.finite(dt) || ds != floor(ds) || dt != floor(dt) ||
        ds < 1 || dt <= ds || ds > 32767 || dt > 32767) {
      return(NULL)
    }
    pk <- as.integer(ds) * 65536L + as.integer(dt)
    pi <- 0L
    m <- match(pk, pair_keys, 0L)
    if (m > 0L) {
      if (identical(pairs[[m]], pr)) {
        pi <- m
      } else {
        for (ci in which(pair_keys == pk)) {
          if (identical(pairs[[ci]], pr)) { pi <- ci; break }
        }
      }
    }
    if (!pi) {
      pi <- length(pairs) + 1L
      pairs[[pi]] <- pr
      pair_keys[pi] <- pk
      p_ds[pi] <- as.double(ds)
      p_dt[pi] <- as.double(dt)
    }
    spec <- list(pair_id = pi, cohorts = control_cohorts)
    qk <- if (clen == 1L) {
      pi * 65536L + control_cohorts[1L]
    } else {
      -(pi * 65536L + control_cohorts[1L] * 256L + min(clen, 255L))
    }
    ri <- 0L
    mq <- match(qk, pool_keys, 0L)
    if (mq > 0L) {
      if (qk > 0L || identical(pools[[mq]], spec)) {
        ri <- mq
      } else {
        for (ci in which(pool_keys == qk)) {
          if (identical(pools[[ci]], spec)) { ri <- ci; break }
        }
      }
    }
    if (!ri) {
      ri <- length(pools) + 1L
      pools[[ri]] <- spec
      pool_keys[ri] <- qk
      if (clen > max_clen) max_clen <- clen
    }
    links[[ell]] <- list(
      regression = ri,
      pair_id = pi,
      cohort = target,
      link_id = lid,
      source_functional_key = source,
      functional_key = paste0(fk_prefix, substring(source, 8L))
    )
  }
  if (anyDuplicated(ids)) return(NULL)
  est_entries <- (128 + 8 * J) * L + 2 * (3 * sum_kb + sum_idb + 128 * L)
  est_bytes <- 3 * sum_kb + sum_idb + 256 * L
  max_sc <- header + max(16 * max_kb, est_entries)
  max_wk <- external_work_ops + 128 * L * L * (J + 1) + 64 * max(max_kb, est_bytes)
  if (!is.finite(max_sc) || !is.finite(max_wk) ||
      base + header > ctrl$max_retained_entries ||
      base + header + max_sc > ctrl$max_matrix_entries ||
      max_wk > ctrl$max_work_ops ||
      max_sc > ctrl$max_output_entries) {
    return(NULL)
  }
  n_pairs <- length(pairs)
  po <- order(p_ds, p_dt)
  pairs <- pairs[po]
  transitions <- p_dt[po] - p_ds[po]
  pair_map <- integer(n_pairs)
  pair_map[po] <- seq_len(n_pairs)
  n_pools <- length(pools)
  pool_pi <- integer(n_pools)
  for (i in seq_len(n_pools)) {
    npi <- pair_map[pools[[i]]$pair_id]
    pools[[i]]$pair_id <- npi
    pool_pi[i] <- npi
  }
  if (max_clen == 1L) {
    pool_c1 <- vapply(pools, function(s) s$cohorts[1L], integer(1L))
    qo <- order(pool_pi, pool_c1, method = "radix")
  } else {
    cols <- vector("list", 1L + max_clen)
    cols[[1L]] <- pool_pi
    for (j in seq_len(max_clen)) {
      col_j <- integer(n_pools)
      for (i in seq_len(n_pools)) {
        cc <- pools[[i]]$cohorts
        if (length(cc) >= j) col_j[i] <- cc[j]
      }
      cols[[1L + j]] <- col_j
    }
    qo <- do.call(order, c(cols, list(method = "radix")))
  }
  pools <- pools[qo]
  pool_map <- integer(n_pools)
  pool_map[qo] <- seq_len(n_pools)
  used_mat <- matrix(FALSE, J, n_pairs)
  for (i in seq_len(n_pools)) {
    s <- pools[[i]]
    used_mat[s$cohorts, s$pair_id] <- TRUE
  }
  for (ell in seq_len(L)) {
    npi <- pair_map[links[[ell]]$pair_id]
    links[[ell]]$pair_id <- npi
    links[[ell]]$regression <- pool_map[links[[ell]]$regression]
    used_mat[links[[ell]]$cohort, npi] <- TRUE
  }
  n_specs <- sum(used_mat)
  specs <- vector("list", n_specs)
  sup_mat <- matrix(0L, J, n_pairs)
  idx <- 0L
  for (pi in seq_len(n_pairs)) {
    pr <- pairs[[pi]]
    tr <- transitions[pi]
    for (k in which(used_mat[, pi])) {
      idx <- idx + 1L
      specs[[idx]] <- list(pair = pr, pair_id = pi, cohort = k, transitions = tr)
      sup_mat[k, pi] <- idx
    }
  }
  cohorts_chr <- as.character(cohorts_vec)
  bk_suffix <- if (identical(route, "DR") || !is.null(selection_model)) {
    paste0(
      "|", recipe$key, "|weights:w*a_U/(rho_U*p_U)|U:", st$universe,
      "|model:", selection_model, "|pooling:", selection_pooling
    )
  } else {
    paste0("|", recipe$key, "|weights:w*a_U/(rho_U*p_U)|U:", st$universe)
  }
  for (i in seq_len(n_pools)) {
    s <- pools[[i]]
    pi <- s$pair_id
    cc <- s$cohorts
    pr <- pairs[[pi]]
    pool_str <- if (length(cc) == 1L) cohorts_chr[cc] else paste(cohorts_chr[cc], collapse = ",")
    pools[[i]]$support_ids <- as.integer(sup_mat[cc, pi])
    pools[[i]]$transitions <- transitions[pi]
    pools[[i]]$block_key <- paste0("pair:", pr$s, "-", pr$t, "|pool:", pool_str, bk_suffix)
  }
  for (ell in seq_len(L)) {
    l <- links[[ell]]
    links[[ell]]$support_id <- sup_mat[l$cohort, l$pair_id]
  }
  list(specs = specs, pools = pools, links = links)
}

cdid_fast_sel_block_fit_batch <- function(reader, block_specs, cells, design_spec,
                                          gc, ctrl, selection_model,
                                          selection_pooling, reference = NULL) {
  n_m <- length(block_specs)
  if (n_m < 1L || isTRUE(gc$trace) ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "fit_pooled_selection_block",
        "selection_owner_block_provider",
        "cp04_selection_fit",
        "cp04_pair_rows",
        "cp04_date_rows",
        "cp04_logit",
        "ipw_glm_fit_internal",
        "ipw_logit_values",
        "fit_selection_block_provider",
        "or_pooled_wls_fit",
        "or_wls_fit",
        "chained_reference_s",
        "chained_reference_stat",
        "chained_reference_weight"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_fit_block")
  if (is.null(fk)) return(NULL)
  is_date <- identical(selection_model, "independent_period_logit")
  is_pair_sep <- identical(selection_model, "pair_logit") &&
    identical(selection_pooling, "separate") && design_spec$d > 0L
  has_int <- design_spec$has_intercept
  d_s <- design_spec$d_s
  if (!is_pair_sep && !has_int && d_s == 0L) return(NULL)
  cn <- if (is.null(reference)) colnames(reader$state$Z) else
    reference$designs[[reference$role_slot["s"]]]$columns
  dn_h <- if (!is.null(cn)) list(cn, cn) else NULL
  int_col <- design_spec$int_col
  slope_cols <- design_spec$slope_cols
  coord_order <- c(int_col, slope_cols)
  empty_warn <- character(0L)
  empty_num <- numeric(0L)
  piv_ds <- seq_len(d_s)
  zero_ds <- numeric(d_s)
  params_sep <- c(2, 0, gc$epsilon, as.double(gc$maxit))
  params_pool <- c(1, if (is_date) 1 else 0, gc$epsilon, as.double(gc$maxit))
  out <- vector("list", n_m)
  for (m in seq_len(n_m)) {
    bs <- block_specs[[m]]
    cids <- as.integer(bs$cell_ids)
    B_m <- length(cids)
    if (B_m < 1L) return(NULL)
    if (is_pair_sep && B_m == 1L) {
      cell <- cells[[cids[1L]]]
      sp <- if (is.null(cell$reference)) list(pair = cell$pair, cohort = cell$cohort) else
        list(pair = cell$pair, cohort = cell$cohort, reference = cell$reference)
      res <- fk(reader, sp, cids, design_spec, reference, params_sep)
      if (is.null(res)) {
        out[[m]] <- fit_pooled_selection_block(
          reader = reader, block_spec = bs, cells = cells, design_spec = design_spec,
          glm_control = gc, ctrl = ctrl, selection_model = selection_model,
          selection_pooling = selection_pooling, reference = reference
        )
        next
      }
      coef <- res[[1L]]
      H_chol <- res[[2L]]
      if (!is.null(dn_h)) dimnames(H_chol) <- dn_h
      out[[m]] <- list(
        branch = "cp04_separate_pair",
        has_intercept = has_int,
        B_m = 1L,
        d_s = d_s,
        alpha = if (has_int) coef[int_col] else empty_num,
        beta = coef[slope_cols],
        coefficients = coef,
        H_chol = H_chol,
        coordinate_order = coord_order,
        m_b = res[[6L]],
        z_bar = res[[7L]],
        slope_fit = NULL,
        pi = res[[4L]],
        r = res[[5L]],
        diagnostics = list(
          converged = TRUE,
          boundary = FALSE,
          rank = length(coef),
          iter = res[[3L]],
          warnings = empty_warn
        ),
        block_id = bs$block_id,
        cell_ids = cids
      )
    } else {
      res <- fk(reader, cells, cids, design_spec, reference, params_pool)
      if (is.null(res)) {
        out[[m]] <- fit_pooled_selection_block(
          reader = reader, block_spec = bs, cells = cells, design_spec = design_spec,
          glm_control = gc, ctrl = ctrl, selection_model = selection_model,
          selection_pooling = selection_pooling, reference = reference
        )
        next
      }
      alpha <- res[[1L]]
      beta <- res[[2L]]
      p_m <- (if (has_int) B_m else 0L) + d_s
      slope_fit <- if (d_s > 0L) {
        structure(
          list(
            coefficients = zero_ds,
            rss = 0,
            R = res[[5L]],
            pivot = piv_ds,
            n = res[[9L]],
            p = d_s,
            rank = d_s,
            rank_tol = res[[10L]],
            cutoff = res[[11L]]
          ),
          class = "cdid_or_wls_fit"
        )
      } else {
        NULL
      }
      out[[m]] <- list(
        branch = design_spec$branch,
        has_intercept = has_int,
        B_m = B_m,
        d_s = d_s,
        alpha = alpha,
        beta = beta,
        coefficients = c(alpha, beta),
        H_chol = NULL,
        m_b = res[[3L]],
        z_bar = res[[4L]],
        slope_fit = slope_fit,
        pi = res[[6L]],
        r = res[[7L]],
        diagnostics = list(
          converged = TRUE,
          boundary = FALSE,
          rank = p_m,
          iter = res[[8L]],
          warnings = empty_warn
        ),
        block_id = bs$block_id,
        cell_ids = cids
      )
    }
  }
  out
}

cdid_fast_resolve_q_of_s <- function(sel_fit, specs) {
  n_s <- length(specs)
  pcs <- sel_fit$pair_cohorts
  n_pc <- length(pcs)
  if (n_s < 1L || n_pc < 1L || !cdid_is_default_fn(resolve_pair_cohort_id, "resolve_pair_cohort_id")) {
    return(NULL)
  }
  has_ref <- !is.null(sel_fit$reference_ref)
  pc_keys <- numeric(n_pc)
  for (q in seq_len(n_pc)) {
    pc <- pcs[[q]]
    pr <- pc$pair
    r <- if (has_ref) as.integer(pc$reference) else 0L
    pc_keys[q] <- ((as.double(r) * 4096 + as.double(pc$cohort)) * 32768 +
      as.double(pr$s_calendar_index)) * 32768 + as.double(pr$t_calendar_index)
  }
  sp_keys <- numeric(n_s)
  for (j in seq_len(n_s)) {
    sp <- specs[[j]]
    pr <- sp$pair
    r <- if (has_ref) as.integer(sp$reference) else 0L
    sp_keys[j] <- ((as.double(r) * 4096 + as.double(sp$cohort)) * 32768 +
      as.double(pr$s_calendar_index)) * 32768 + as.double(pr$t_calendar_index)
  }
  q_ids <- match(sp_keys, pc_keys)
  if (anyNA(q_ids)) return(NULL)
  for (j in seq_len(n_s)) {
    pc <- pcs[[q_ids[j]]]
    sp <- specs[[j]]
    if (pc$cohort != sp$cohort || pc$pair$s != sp$pair$s || pc$pair$t != sp$pair$t ||
        !identical(as.integer(pc$reference), as.integer(sp$reference))) {
      return(NULL)
    }
  }
  q_ids
}

cdid_fast_sel_dr_support_batch <- function(reader, specs, sel_fit, up_map, q_of_s,
                                           recipe, or_frame_index, ctrl, reference = NULL) {
  n_s <- length(specs)
  if (n_s < 1L || !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "dr_selection_support", "cp04_pair_rows", "selection_prediction_operation",
        "or_pooled_support_add", "or_pooled_wls_block_summary",
        "or_pooled_wls_merge_chunks", "or_wls_summary", "or_wls_merge",
        "or_wls_compress", "or_selection_moments_add",
        "chained_reference_x", "chained_reference_s",
        "chained_reference_stat", "chained_reference_weight"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_support_collect")
  if (is.null(fk)) return(NULL)
  st <- reader$state
  frame <- get0("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(frame)) return(NULL)
  z_idx <- as.integer(or_frame_index)
  slope_idx <- as.integer(recipe$slope_index)
  icpt_flag <- isTRUE(recipe$intercept)
  icpt_int <- as.integer(icpt_flag)
  icpt_col <- if (icpt_flag) as.integer(recipe$intercept_index) else 0L
  max_chunk <- as.integer(ctrl$max_chunk_rows)
  p <- as.integer(recipe$p)
  slope_names <- if (p > 0L) recipe$slope_columns else NULL
  raw_cols <- if (icpt_flag) c("", slope_names) else slope_names
  k_aug <- p + 1L + icpt_int
  Za <- if (is.null(reference)) frame$Z_a else NULL
  a_vec <- st$a
  out <- vector("list", n_s)
  for (j in seq_len(n_s)) {
    spec <- specs[[j]]
    L_val <- as.double(spec$transitions)
    cols_info <- list(
      Z_a = Za, z_idx = z_idx, slope_idx = slope_idx,
      icpt_flag = icpt_int, icpt_col = icpt_col,
      L_val = L_val, max_chunk_rows = max_chunk
    )
    res <- fk(3L, reader, spec, sel_fit, up_map, q_of_s[j], reference, cols_info)
    if (is.null(res)) return(NULL)
    a_val <- if (is.null(reference)) a_vec[spec$cohort] else
      chained_reference_stat(reference, spec$reference, spec$cohort, "a")
    fac <- res[[1L]]
    sc <- res[[6L]]
    n_obs_int <- as.integer(sc[1L])
    n_obs_dbl <- as.double(sc[1L])
    weight_mass <- sc[4L]
    raw_sum <- structure(
      list(factor = fac, n = n_obs_int, columns = raw_cols),
      class = "cdid_or_wls_summary"
    )
    summary_obj <- if (icpt_flag) {
      structure(
        list(
          raw = raw_sum, n = n_obs_int, p = p, columns = slope_names,
          anchor = res[[2L]], constant = res[[3L]],
          accounting = list(
            input_rows = n_obs_int, augmented_cols = k_aug,
            augmented_rows = nrow(fac),
            max_chunk_merge_rows = as.integer(sc[2L]),
            n_compressions = as.integer(sc[3L])
          )
        ),
        class = "cdid_or_pooled_wls_chunk"
      )
    } else {
      raw_sum
    }
    out[[j]] <- list(
      design = list(
        L = L_val, summary = summary_obj, weight_mass = weight_mass,
        n = n_obs_dbl, sumZ = res[[4L]], sumY = sc[5L], colmax = res[[5L]]
      ),
      moments = res[[7L]],
      a = a_val,
      pair_mass = weight_mass,
      n_pair = n_obs_dbl
    )
  }
  out
}

cdid_fast_sel_dr_trt_batch <- function(reader, ts, sel_fit, up_map, state, separate,
                                       recipe, or_frame_index, trt_recipe,
                                       gc, ctrl, reference = NULL) {
  n_t <- length(ts)
  if (n_t < 1L || isTRUE(gc$trace) ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "dr_selection_treatment_fit", "or_outcome_block_fit",
        "cp04_pair_rows", "selection_prediction_operation", "cp04_logit",
        "ipw_glm_fit_internal", "ipw_logit_values", "ipw_clip_summary",
        "chained_reference_x", "chained_reference_s",
        "chained_reference_stat", "chained_reference_weight"
      ))) {
    return(NULL)
  }
  fk <- cdid_sel_fast_kernel("C_cdid_sel_trt_fit")
  if (is.null(fk)) return(NULL)
  frame <- get0("frame_owner", reader$frame_ref, inherits = FALSE)
  if (is.null(frame)) return(NULL)
  Za <- if (is.null(reference)) frame$Z_a else NULL
  trt_idx <- as.integer(trt_recipe$index)
  trt_cols <- trt_recipe$columns
  dn_u <- list(c("T", "M", "dp", "dm"), NULL)
  dn_v <- list(c("T", "M", "dp", "dm"), trt_cols)
  dn_h <- if (!is.null(trt_cols)) list(trt_cols, trt_cols) else NULL
  means_names <- c("T", "M", "dp", "dm", "T_dr", "M_dr", "mu_dr",
                   "T_ipw", "M_ipw", "mu_ipw", "mu_or", "mu_dr_decomp")
  empty_warn <- character(0L)
  ctrl_vec <- c(as.double(gc$epsilon), as.double(gc$maxit))
  has_icpt <- isTRUE(recipe$intercept)
  slope_idx <- recipe$slope_index
  out <- vector("list", n_t)
  for (i in seq_len(n_t)) {
    spec <- ts[[i]]
    reg <- state$regressions[[spec$regression]]
    sh <- if (separate) {
      state$outcome_groups$groups[[state$outcome_groups$block_map$group_id[spec$regression]]]
    } else {
      state$shared
    }
    p_or <- length(sh$keep)
    kept_or_cols <- if (p_or > 0L) as.integer(or_frame_index[slope_idx[sh$keep]]) else integer(0L)
    or_cols <- if (p_or > 0L) sh$basis$kept_columns else NULL
    or_info <- list(
      L = reg$L,
      alpha = if (has_icpt) reg$alpha else 0,
      has_icpt = has_icpt,
      beta = sh$beta,
      trt_cols = trt_cols,
      or_cols = or_cols
    )
    cols_info <- list(Z_a = Za, trt_idx = trt_idx, or_idx = kept_or_cols)
    res <- fk(3L, reader, spec, sel_fit, up_map, reference, or_info, cols_info, ctrl_vec)
    if (is.null(res)) {
      sh_fb <- if (separate) or_outcome_block_fit(state$outcome_groups, spec$regression) else state$shared
      out[[i]] <- dr_selection_treatment_fit(
        reader, spec, sel_fit, up_map, reg, sh_fb,
        recipe, or_frame_index, trt_recipe, gc, ctrl, reference
      )
      next
    }
    beta <- res[[1L]]
    H_chol <- res[[2L]]
    U_main <- res[[5L]]
    V_main <- res[[6L]]
    U_ipw <- res[[10L]]
    V_ipw <- res[[11L]]
    dimnames(U_main) <- dn_u
    dimnames(V_main) <- dn_v
    dimnames(U_ipw) <- dn_u
    dimnames(V_ipw) <- dn_v
    if (!is.null(trt_cols)) {
      names(beta) <- trt_cols
      dimnames(H_chol) <- dn_h
    }
    means <- res[[7L]]
    names(means) <- means_names
    dbp <- res[[12L]]; dbm <- res[[13L]]; dbd <- res[[14L]]
    if (length(dbp) > 0L && !is.null(or_cols)) {
      names(dbp) <- names(dbm) <- names(dbd) <- or_cols
    }
    out[[i]] <- list(
      coefficients = beta,
      H_chol = H_chol,
      diagnostics = list(
        converged = TRUE, boundary = FALSE, rank = length(beta),
        iter = res[[3L]], warnings = empty_warn
      ),
      upstream_ids = res[[9L]],
      cross = res[[4L]],
      partials = list(upstream = U_main, beta = V_main),
      partials_ipw = list(upstream = U_ipw, beta = V_ipw),
      means = means,
      d_bar_plus = dbp,
      d_bar_minus = dbm,
      d_bar_diff = dbd,
      clipping = res[[8L]]
    )
  }
  out
}

cdid_fast_eval_dr_sep_if_chunk <- function(backend, rows, links, scale, details, plan) {
  if (!identical(details, FALSE) || !is.list(backend) ||
      !isTRUE(.subset2(backend, "uses_selection_models")) ||
      !is.null(.subset2(backend, "layout")) ||
      !is.null(.subset2(backend, "pair_outcomes")) ||
      !is.character(scale) || is.object(scale) || length(scale) != 1L ||
      !(scale == "frame" || scale == "core") ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "validate_dr_backend", "validate_dr_corollary1_backend",
        "dr_separate_if_envelope", "dr_conditional_replay_check",
        "dr_conditional_replay_upstream_index", "cp04_budget", "cp04_payload",
        "chained_backend_reference", "cp04_pair_consumer_operation",
        "cp04_pair_read_operation", "selection_prediction_operation",
        "cp04_pair_rows", "cp04_date_rows", "ipw_owner_read_y_positions",
        "or_outcome_block_fit", "or_selection_separate_if_core",
        "or_baseline_separate_if_core", "cp04_solve", "ipw_logit_values",
        "or_corollary1_factor", "chained_reference_x", "chained_reference_s",
        "chained_reference_stat", "chained_reference_weight"
      ))) {
    return(NULL)
  }
  ring <- get0("val_ring", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(ring) || !length(ring)) return(NULL)
  fk_vf <- cdid_sel_fast_kernel("C_cdid_verify_fast")
  fk_if <- cdid_sel_fast_kernel("C_cdid_sel_if_chunk")
  if (is.null(fk_vf) || is.null(fk_if)) return(NULL)
  hit_ok <- FALSE
  for (i in seq_along(ring)) {
    rec <- ring[[i]]
    if (is.list(rec) && identical(rec[[2L]], "dr_corollary1_backend") &&
        isTRUE(fk_vf(4L, backend, rec, NULL, NULL))) {
      hit_ok <- TRUE
      break
    }
  }
  if (!hit_ok) return(NULL)
  reader <- get("cp04_reader", .subset2(backend, "reader_ref"), inherits = FALSE)
  st <- reader$state
  B <- length(rows)
  if (!cdid_sel_valid_idx(rows, st$N) || B > backend$control$max_chunk_rows ||
      !cdid_sel_valid_idx(links, length(backend$links))) {
    return(NULL)
  }
  links_int <- as.integer(links)
  bref <- .subset2(backend, "ref")
  slots <- get0("dr_sep_prep_slots", envir = cdid_sel_cache, inherits = FALSE)
  if (!is.list(slots) || length(slots) != 2L) slots <- list(NULL, NULL)
  hit <- NULL
  for (si in 1:2) {
    sl <- slots[[si]]
    if (is.list(sl) && identical(sl$backend_ref, bref) &&
        gmm_sexp_aliases(sl$backend, backend) &&
        sl$B == B && identical(sl$links, links_int) &&
        ((is.null(plan) && is.null(sl$plan)) ||
         (!is.null(plan) && !is.null(sl$plan) && gmm_sexp_aliases(sl$plan, plan)))) {
      hit <- sl
      break
    }
  }
  if (is.null(hit)) {
    pre <- if (is.null(plan)) dr_separate_if_envelope(backend, B, links, FALSE) else
      dr_conditional_replay_check(plan, backend, links, B)
    cp04_budget(
      backend$control,
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries),
      pre$scratch,
      c(backend$accounting$work_ops_bound, pre$work),
      pre$output,
      "separate DR IF chunk"
    )
    prep_info <- if (!is.null(plan)) list(
      up_idx = dr_conditional_replay_upstream_index(
        plan, backend$selection_backend, backend$up_map, pre$uids, FALSE
      )
    ) else NULL
    link_names <- if (!is.null(plan) && !is.null(plan$link_names)) plan$link_names else
      vapply(backend$links[links], `[[`, character(1L), "link_id")
    ref <- chained_backend_reference(backend)
    hit <- list(
      backend = backend, backend_ref = bref, plan = plan, B = B, links = links_int,
      pre = pre, prep_info = prep_info, link_names = link_names, ref = ref,
      n_rate = sum(pre$uids <= backend$upstream$n_a),
      pay_env = new.env(parent = emptyenv())
    )
    s_idx <- if (!is.null(plan) && isTRUE(B != plan$B)) 2L else 1L
    slots[[s_idx]] <- hit
    assign("dr_sep_prep_slots", slots, envir = cdid_sel_cache)
  }
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  Za <- if (is.null(hit$ref)) frame$Z_a else NULL
  res <- fk_if(3L, backend, reader, Za, hit$ref, rows, links, hit$pre, hit$prep_info)
  if (is.null(res)) return(NULL)
  phi <- res[[1L]]
  dimnames(phi) <- list(NULL, hit$link_names)
  if (scale == "core") phi <- st$lambda * phi
  pre <- hit$pre
  out <- list(
    rows = rows,
    id = cdid_sel_slice_id(st$id, rows, st$N),
    links = links,
    phi = phi,
    scale = scale,
    normalization_count = if (scale == "core") st$n_A else st$N,
    diagnostics = list(
      outcome_pooling = "separate",
      n_rate_score_evaluations = hit$n_rate,
      n_selection_score_evaluations = res[[2L]],
      n_selection_block_evaluations = res[[3L]],
      n_regression_score_evaluations = length(pre$gids),
      n_global_beta_blocks = 0L,
      n_block_beta_blocks = length(pre$gids),
      n_treatment_score_evaluations = length(pre$tids),
      output_entries_bound = pre$output,
      concurrent_scratch_bound = pre$scratch,
      work_ops_bound = pre$work
    )
  )
  can_skip_pay <- !is.character(st$id) && is.null(attributes(rows)) &&
    is.null(attributes(links)) && is.null(attributes(st$id))
  if (can_skip_pay && isTRUE(get0(scale, envir = hit$pay_env, inherits = FALSE))) {
    return(out)
  }
  if (unname(cp04_payload(out)["entries"]) > pre$output) {
    dr_route_abort("corollary1", "IF output bound underestimated.", "bound")
  }
  if (can_skip_pay) assign(scale, TRUE, envir = hit$pay_env)
  out
}

cdid_fast_or_sep_forcing_links <- function(recipe, fitted, supports, pools, links, channels, N) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c(
        "or_pooled_factor_moments", "or_pooled_link_moments", "cp04_finite"
      ))) {
    return(NULL)
  }
  n_b <- length(pools)
  n_l <- length(links)
  has_icpt <- isTRUE(recipe$intercept)
  m <- as.integer(recipe$p) + 2L
  y_col <- as.integer(recipe$p) + 1L + as.integer(has_icpt)
  groups <- fitted$groups
  regressions <- vector("list", n_b)
  for (b in seq_len(n_b)) {
    f <- groups[[b]]
    keep <- f$keep
    p <- length(keep)
    xi <- if (has_icpt) c(1L, 1L + keep) else keep
    rec <- f$regressions[[1L]]
    sids <- pools[[b]]$support_ids
    single_sorted <- length(sids) == 1L && !is.unsorted(channels[[sids[1L]]], strictly = TRUE)
    ids <- if (single_sorted) channels[[sids[1L]]] else sort.int(unique(unlist(channels[sids], use.names = FALSE)))
    n_id <- length(ids)
    fa <- numeric(n_id)
    fe <- numeric(n_id)
    fb <- matrix(0, p, n_id)
    cr <- numeric(m)
    cr[1L] <- -rec$alpha
    if (p > 0L) cr[1L + keep] <- -rec$L * f$beta
    cr[m] <- 1
    c_ab <- if (has_icpt) c(-rec$alpha, -f$beta) else -f$beta
    mean_alpha <- 0
    mean_beta <- numeric(p)
    for (s in sids) {
      u <- supports[[s]]
      ch_s <- channels[[s]]
      if (single_sorted) {
        ka <- 1L
        kl <- if (n_id > 1L) 2L:n_id else integer(0L)
      } else {
        pos <- match(ch_s, ids)
        ka <- pos[1L]
        kl <- pos[-1L]
      }
      Fmat <- if (has_icpt) u$design$summary$raw$factor else u$design$summary$factor
      if (nrow(Fmat) > 0L) {
        Fsub <- Fmat[, xi, drop = FALSE]
        rho <- Fmat[, y_col] + as.vector(Fsub %*% c_ab)
        g_vec <- as.vector(crossprod(Fsub, rho))
        e2 <- sum(rho^2)
      } else {
        g_vec <- numeric(length(xi))
        e2 <- 0
      }
      inv_Na <- 1 / (N * u$a)
      if (has_icpt) {
        ga <- g_vec[1L]
        gb <- g_vec[-1L] - rec$mean_d * ga
        fa[ka] <- fa[ka] + ga * inv_Na
        fb[, ka] <- fb[, ka] + gb * inv_Na
      } else {
        ga <- 0
        gb <- g_vec
        fb[, ka] <- fb[, ka] + gb * inv_Na
      }
      mean_alpha <- mean_alpha + ga / N
      mean_beta <- mean_beta + gb / N
      fe[ka] <- fe[ka] + e2 * inv_Na
      dq_s <- length(kl)
      if (dq_s > 0L) {
        Mr <- matrix(crossprod(matrix(u$moments, nrow = m), cr), m, dq_s)
        if (has_icpt) {
          fa[kl] <- fa[kl] + Mr[1L, ] / N
          if (p > 0L) {
            fb[, kl] <- fb[, kl] + (rec$L * Mr[1L + keep, , drop = FALSE] -
              tcrossprod(rec$mean_d, Mr[1L, ])) / N
          }
        } else if (p > 0L) {
          fb[, kl] <- fb[, kl] + rec$L * Mr[1L + keep, , drop = FALSE] / N
        }
        fe[kl] <- fe[kl] + colSums(cr * Mr) / N
      }
    }
    if (any(!is.finite(fa)) || any(!is.finite(fb)) || any(!is.finite(fe))) {
      cp04_abort("Nonfinite Corollary 1 regression forcing.", "numeric")
    }
    rec$channel_ids <- ids
    rec$forcing_alpha <- fa
    rec$forcing_beta <- fb
    rec$forcing_energy <- fe
    rec$score_mean_alpha <- mean_alpha
    rec$score_mean_beta <- mean_beta
    regressions[[b]] <- rec
  }
  links_out <- links
  for (i in seq_len(n_l)) {
    l <- links[[i]]
    reg_id <- l$regression
    f <- groups[[reg_id]]
    rec <- regressions[[reg_id]]
    u <- supports[[l$support_id]]
    keep <- f$keep
    p <- length(keep)
    dq <- dim(u$moments)[3L]
    mass <- u$pair_mass
    sumZk <- u$design$sumZ[keep]
    sumY <- u$design$sumY
    L_val <- rec$L
    fitted_val <- (if (has_icpt) rec$alpha * mass else 0) + L_val * sum(sumZk * f$beta)
    mu <- (sumY - fitted_val) / mass
    c_vec <- L_val * sumZk / mass
    resid <- sumY - fitted_val - mu * mass
    d_val <- mass / N
    num <- mu * d_val
    if (dq > 0L) {
      cr <- numeric(m)
      cr[1L] <- -rec$alpha
      if (p > 0L) cr[1L + keep] <- -L_val * f$beta
      cr[m] <- 1
      Mr <- matrix(crossprod(matrix(u$moments, nrow = m), cr), m, dq)
      m11 <- u$moments[1L, 1L, ]
      rl <- c(resid / (N * u$a * d_val), (Mr[1L, ] - mu * m11) / (N * d_val))
      wl <- c(d_val / u$a, m11 / N)
      nl <- c(num / u$a, Mr[1L, ] / N)
      dl <- matrix(c(c_vec * d_val / u$a, if (p > 0L) L_val * u$moments[1L + keep, 1L, ] / N else numeric(0L)), p, 1L + dq)
    } else {
      rl <- c(resid / (N * u$a * d_val))
      wl <- c(d_val / u$a)
      nl <- c(num / u$a)
      dl <- matrix(c_vec * d_val / u$a, p, 1L)
    }
    if (!is.finite(mu) || any(!is.finite(c_vec)) || !is.finite(d_val) || !is.finite(num) ||
        any(!is.finite(rl)) || any(!is.finite(wl)) || any(!is.finite(nl)) || any(!is.finite(dl))) {
      cp04_abort("Nonfinite Corollary 1 link moments.", "numeric")
    }
    l$mu <- mu
    l$c <- c_vec
    l$d <- d_val
    l$numerator <- num
    l$rate_loading <- rl
    l$weight_loading <- wl
    l$numerator_loading <- nl
    l$design_loading <- dl
    links_out[[i]] <- l
  }
  list(regressions = regressions, links = links_out)
}
