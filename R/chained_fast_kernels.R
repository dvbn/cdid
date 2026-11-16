# Private native fast-kernel dispatch helpers. All kernels are registered in
# `src/cdid.c` and bound through `gmm_require_native_dll()`'s owner-verified
# `calls` table; when a probe/mock DLL omits optional `C_cdid_*` routines,
# `cdid_fast_kernel()` returns NULL so callers retain exact R fallback parity.

cdid_fast_kernel_cache <- new.env(parent = emptyenv())
cdid_default_interrupt_body <- quote({
  gmm_sexp_aliases(NULL, NULL)
  invisible(NULL)
})

cdid_fast_kernel <- function(name) {
  calls <- cdid_fast_kernel_cache$calls
  if (!is.null(calls) &&
      identical(cdid_fast_kernel_cache$dll, cdid_fast_kernel_cache$ns_info$DLLs[["cdid"]]) &&
      identical(cdid_fast_kernel_cache$syms, gmm_native_dll_cache$syms) &&
      !inherits(gmm_require_native_dll, "functionWithTrace")) {
    return(calls[[name]])
  }
  syms <- tryCatch(gmm_require_native_dll(), error = function(e) NULL)
  if (is.null(syms)) {
    cdid_fast_kernel_cache$calls <- NULL
    return(NULL)
  }
  raw_calls <- attr(syms, "calls")
  owner_env <- environment(cdid_fast_kernel)
  if (is.list(raw_calls) &&
      isNamespace(owner_env) &&
      environmentIsLocked(owner_env) &&
      identical(environment(gmm_require_native_dll), owner_env) &&
      !inherits(gmm_require_native_dll, "functionWithTrace")) {
    ns_info <- get0(".__NAMESPACE__.", envir = owner_env, mode = "environment", inherits = FALSE)
    dll <- if (!is.null(ns_info)) ns_info$DLLs[["cdid"]] else NULL
    if (!is.null(dll) && identical(gmm_native_dll_cache$info, dll[["info"]])) {
      valid_calls <- raw_calls
      for (nm in names(valid_calls)) {
        f <- valid_calls[[nm]]
        if (!is.function(f) || !exists(nm, envir = environment(f), inherits = FALSE)) {
          valid_calls[nm] <- list(NULL)
        }
      }
      cdid_fast_kernel_cache$owner_env <- owner_env
      cdid_fast_kernel_cache$ns_info <- ns_info
      cdid_fast_kernel_cache$dll <- dll
      cdid_fast_kernel_cache$syms <- syms
      cdid_fast_kernel_cache$calls <- valid_calls
      cdid_fast_kernel_cache$vf_fn <- valid_calls[["C_cdid_verify_fast"]]
      return(valid_calls[[name]])
    }
  }
  fn <- if (is.list(raw_calls)) raw_calls[[name]] else NULL
  if (is.function(fn) && exists(name, envir = environment(fn), inherits = FALSE)) fn else NULL
}

cdid_sexp_aliases <- function(a, b) {
  fk <- cdid_fast_kernel("C_gmm_sexp_aliases")
  if (!is.null(fk)) isTRUE(fk(a, b)) else gmm_sexp_aliases(a, b)
}

cdid_interrupt_callback <- function(fn) {
  if (!is.function(fn)) return(fn)
  cached <- get0("default_interrupt", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (identical(fn, cached)) return(NULL)
  owner_env <- environment(cdid_fast_kernel)
  if (identical(environment(fn), owner_env) &&
      identical(body(utils::removeSource(fn)), cdid_default_interrupt_body)) {
    assign("default_interrupt", fn, envir = cdid_fast_kernel_cache)
    return(NULL)
  }
  fn
}

cdid_default_fn_snapshot <- new.env(hash = TRUE, parent = emptyenv())

cdid_snapshot_default_fns <- function(ns = environment(cdid_fast_kernel)) {
  for (nm in ls(ns, all.names = TRUE)) {
    val <- get0(nm, envir = ns, inherits = FALSE)
    if (is.function(val)) assign(nm, val, envir = cdid_default_fn_snapshot)
  }
  invisible(NULL)
}

cdid_is_default_fn <- function(fn, name) {
  if (!is.function(fn) || inherits(fn, "functionWithTrace")) return(FALSE)
  snap <- get0(name, envir = cdid_default_fn_snapshot, inherits = FALSE)
  if (!is.null(snap)) return(identical(fn, snap))
  owner_env <- environment(cdid_fast_kernel)
  if (!identical(environment(fn), owner_env)) return(FALSE)
  orig <- get0(name, envir = owner_env, mode = "function", inherits = FALSE)
  !is.null(orig) && identical(fn, orig)
}

cdid_are_default_fns <- function(names, caller_env = parent.env(parent.frame())) {
  vf_fn <- cdid_fast_kernel_cache$vf_fn
  cached_owner <- cdid_fast_kernel_cache$owner_env
  if (!is.null(vf_fn) && identical(caller_env, cached_owner) &&
      identical(cdid_fast_kernel_cache$syms, gmm_native_dll_cache$syms) &&
      identical(cdid_fast_kernel_cache$dll, cdid_fast_kernel_cache$ns_info$DLLs[["cdid"]]) &&
      length(cdid_default_fn_snapshot) > 0L && is.character(names)) {
    res <- vf_fn(6L, names, cached_owner, cdid_default_fn_snapshot, NULL)
    if (is.logical(res) && length(res) == 1L && !is.na(res)) return(res)
  }
  owner_env <- environment(cdid_fast_kernel)
  has_snap <- length(cdid_default_fn_snapshot) > 0L
  if (identical(caller_env, owner_env) && has_snap && is.character(names)) {
    fk <- cdid_fast_kernel("C_cdid_verify_fast")
    if (is.function(fk)) {
      res <- fk(6L, names, owner_env, cdid_default_fn_snapshot, NULL)
      if (is.logical(res) && length(res) == 1L && !is.na(res)) return(res)
    }
  }
  lookup_env <- if (is.environment(caller_env)) caller_env else owner_env
  use_inh <- !identical(lookup_env, owner_env)
  for (nm in names) {
    fn <- get0(nm, envir = lookup_env, inherits = use_inh)
    if (!is.function(fn) || inherits(fn, "functionWithTrace")) return(FALSE)
    if (has_snap) {
      snap <- get0(nm, envir = cdid_default_fn_snapshot, inherits = FALSE)
      if (is.null(snap) || !identical(fn, snap)) return(FALSE)
    } else if (!identical(environment(fn), owner_env)) {
      return(FALSE)
    }
  }
  TRUE
}

cdid_is_default_glm_fit <- function(fn) {
  cdid_is_default_fn(fn, "ipw_glm_fit_internal")
}

cdid_fast_owner_pair_support <- function(frame_owner, pair_desc, coh_desc, step, require_support = TRUE) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns("ipw_owner_chunked_pair_support")) {
    return(NULL)
  }
  fk <- cdid_fast_kernel("C_cdid_pair_support_scan")
  if (is.null(fk)) return(NULL)
  n_a <- frame_owner$n_a_units
  m <- coh_desc$m_cohorts
  ds <- pair_desc$s_calendar_index
  dt <- pair_desc$t_calendar_index
  post <- frame_owner$postings
  has_pair <- !is.na(ds) && !is.na(dt) && post$date_count[ds] > 0L && post$date_count[dt] > 0L
  bounds4 <- if (has_pair) as.integer(c(post$date_start[ds], post$date_end[ds], post$date_start[dt], post$date_end[dt])) else NULL
  C_own <- length(frame_owner$a_cohorts)
  cmap <- integer(C_own)
  if (length(coh_desc$owner_cohort_index) > 0L) {
    cmap[as.integer(coh_desc$owner_cohort_index)] <- as.integer(coh_desc$view_cohort_index)
  }
  cb <- cdid_interrupt_callback(ipw_owner_check_interrupt)
  fast_res <- fk(post$post_a_index, bounds4, frame_owner$a_cohort_index, frame_owner$w_a,
                 cmap, as.integer(c(m, n_a)), step, cb, environment())
  if (is.null(fast_res)) return(NULL)
  pair_count <- fast_res$pair_count
  pair_mass <- fast_res$pair_mass
  n_obs_pair_all <- fast_res$n_observed
  chunks <- fast_res$chunks
  ops <- fast_res$ops

  n_obs_risk <- sum(pair_count)
  if (isTRUE(require_support)) {
    if (n_obs_pair_all < 1L || n_obs_risk < 1L) {
      cdid_ipw_owner_abort(
        "Requested date pair has zero empirical risk-set observations; requested operator is unavailable.",
        c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
      )
    }
    if (any(pair_count < 1L) || any(!is.finite(pair_mass) | pair_mass <= 0)) {
      cdid_ipw_owner_abort(
        "At least one declared cohort in `{g, K}` has zero empirical pair count or mass on `(s, t)`; zero-support cohorts are not silently discarded.",
        c("cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error", "cdid_ipw_support_error")
      )
    }
  }
  rates <- ifelse(coh_desc$group_a_mass > 0, pair_mass / coh_desc$group_a_mass, NA_real_)
  if (isTRUE(require_support)) {
    rates[pair_count == coh_desc$group_a_count] <- 1
    ops <- ops + 4 * m
  }
  if (isTRUE(require_support) && any(!is.finite(rates) | rates <= 0 | rates > 1)) {
    cdid_ipw_owner_abort(
      "Invalid marginal pair observation rates `r_j`.",
      c("cdid_ipw_owner_pair_support_error", "cdid_ipw_rate_error")
    )
  }
  names(pair_count) <- names(pair_mass) <- names(rates) <- coh_desc$cohort_labels

  list(
    n_observed_pair_all = n_obs_pair_all,
    n_risk_a_units = sum(coh_desc$group_a_count),
    n_outside_risk_a_units = n_a - sum(coh_desc$group_a_count),
    n_observed_risk_units = n_obs_risk,
    n_nonobserved_risk_units = sum(coh_desc$group_a_count) - n_obs_risk,
    by_cohort_pair_count = pair_count,
    by_cohort_pair_mass = pair_mass,
    rates = rates,
    chunks = chunks,
    ops = ops
  )
}

cdid_fast_incomplete_or_supports <- function(reader, specs, recipe, ctrl, qr_buffer = NULL, or_frame_index = NULL) {
  if (!cdid_are_default_fns(c(
    "or_wls_compress",
    "or_pooled_support_add",
    "or_baseline_support",
    "or_baseline_visit",
    "or_baseline_qr_buffered_support",
    "cp04_pair_rows",
    "cp04_pair_read_operation",
    "ipw_owner_read_y_positions"
  ))) {
    return(NULL)
  }
  fk <- cdid_fast_kernel("C_cdid_incomplete_or_supports")
  if (is.null(fk) || !length(specs)) return(NULL)
  st <- reader$state
  frame <- get0("frame_owner", envir = reader$frame_ref, inherits = FALSE)
  panel <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  if (is.null(frame) || is.null(panel) || isTRUE(st$reference_carrier)) return(NULL)
  post <- panel$postings
  y <- panel$observations$y
  Z_a <- frame$Z_a
  if (!is.double(y) || is.object(y) || !is.double(Z_a) || !is.matrix(Z_a) || is.object(Z_a) ||
      !is.integer(post$post_a_index) || is.object(post$post_a_index) ||
      !is.integer(post$post_obs_index) || is.object(post$post_obs_index) ||
      !is.integer(post$date_start) || !is.integer(post$date_end) || !is.integer(post$date_count) ||
      !is.integer(st$cohort_order) || !is.integer(st$a_index) || !is.double(st$w)) {
    return(NULL)
  }

  col_map <- if (is.null(or_frame_index)) seq_len(ncol(Z_a)) else as.integer(or_frame_index)
  icpt_flag <- isTRUE(recipe$intercept)
  icpt_col <- if (icpt_flag) col_map[recipe$intercept_index] else 0L
  slope_cols <- if (recipe$p > 0L) as.integer(col_map[recipe$slope_index]) else integer(0)

  S <- length(specs)
  spec_int <- matrix(0L, nrow = 4L, ncol = S)
  spec_L <- numeric(S)
  for (s in seq_len(S)) {
    sp <- specs[[s]]
    ds <- sp$pair$s_calendar_index
    dt <- sp$pair$t_calendar_index
    if (is.null(ds) || is.na(ds) || is.null(dt) || is.na(dt)) return(NULL)
    spec_int[1L, s] <- as.integer(st$cohort_start[sp$cohort])
    spec_int[2L, s] <- as.integer(st$cohort_end[sp$cohort])
    spec_int[3L, s] <- as.integer(ds)
    spec_int[4L, s] <- as.integer(dt)
    spec_L[s] <- as.double(sp$transitions)
  }

  date_bounds <- rbind(post$date_start, post$date_end, post$date_count)
  storage.mode(date_bounds) <- "integer"
  chunk_rows <- as.integer(min(st$N, ctrl$max_chunk_rows))
  use_qr_buf <- !is.null(qr_buffer)
  qr_buf_rows <- if (use_qr_buf) as.integer(qr_buffer$rows) else chunk_rows
  rparams <- list(
    as.integer(c(as.integer(icpt_flag), as.integer(icpt_col), chunk_rows, qr_buf_rows, as.integer(use_qr_buf))),
    slope_cols,
    as.double(st$N)
  )

  raw_list <- fk(y, post$post_a_index, post$post_obs_index, date_bounds,
                 st$cohort_order, st$a_index, st$w, Z_a, spec_int, spec_L, rparams)
  if (is.null(raw_list) || length(raw_list) != S) return(NULL)

  p <- as.integer(recipe$p)
  slope_names <- if (p > 0L) recipe$slope_columns else NULL
  k_aug <- p + 1L + as.integer(icpt_flag)

  out <- vector("list", S)
  for (s in seq_len(S)) {
    it <- raw_list[[s]]
    fac <- it[[1L]]
    anchor <- it[[2L]]
    constant <- it[[3L]]
    sumZ <- it[[4L]]
    colmax <- it[[5L]]
    sc <- it[[6L]]
    n_obs <- as.integer(sc[1L])
    n_obs_dbl <- as.double(sc[1L])
    max_merge <- as.integer(sc[2L])
    n_comp <- as.integer(sc[3L])
    weight_mass <- sc[4L]
    sumY <- sc[5L]
    repeat_mass <- sc[6L]
    n_repeat_dbl <- as.double(sc[7L])
    rate <- sc[8L]

    raw_sum <- structure(
      list(factor = fac, n = n_obs, columns = if (icpt_flag) c("", slope_names) else slope_names),
      class = "cdid_or_wls_summary"
    )
    summary_obj <- if (icpt_flag) {
      structure(
        list(
          raw = raw_sum,
          n = n_obs,
          p = p,
          columns = slope_names,
          anchor = anchor,
          constant = constant,
          accounting = list(
            input_rows = n_obs,
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
      L = spec_L[s],
      summary = summary_obj,
      weight_mass = weight_mass,
      n = n_obs_dbl,
      sumZ = sumZ,
      sumY = sumY,
      colmax = colmax
    )
    out[[s]] <- list(
      design = acc,
      rate = rate,
      v = 1 / ipw_clip_floor(rate),
      massA = repeat_mass / st$N,
      pair_mass = weight_mass,
      repeat_mass = repeat_mass,
      n_pair = n_obs_dbl,
      n_repeat = n_repeat_dbl
    )
  }
  out
}

cdid_fast_seal <- function(obj, tag) {
  if (cdid_are_default_fns("ipw_owner_seal")) {
    fk <- cdid_fast_kernel("C_cdid_fast_seal")
    if (!is.null(fk) && is.list(obj)) {
      nms <- names(obj)
      if (is.character(nms) && !("ref" %in% nms)) obj["ref"] <- list(NULL)
      res <- fk(obj, tag)
      if (!is.null(res)) return(res)
    }
  }
  ipw_owner_seal(obj, tag)
}

cdid_check_validated_backend <- function(backend, seal_tag) {
  if (!cdid_are_default_fns(c(
    "ipw_owner_seal_ok",
    "validate_ipw_frame_owner",
    "validate_ipw_panel_owner",
    "validate_cp04_reader",
    "validate_ipw_owner_control",
    "validate_ipw_glm_control"
  ))) {
    return(FALSE)
  }
  ring <- get0("val_backends", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(ring) || !length(ring)) return(FALSE)
  fk <- cdid_fast_kernel("C_cdid_verify_fast")
  if (is.null(fk)) return(FALSE)
  for (i in seq_along(ring)) {
    rec <- ring[[i]]
    if (is.list(rec) && identical(rec[[2L]], seal_tag) && isTRUE(fk(4L, backend, rec, NULL, NULL))) {
      return(TRUE)
    }
  }
  FALSE
}

cdid_record_validated_backend <- function(backend, seal_tag, reader, frame, panel, fref) {
  rec <- list(backend, seal_tag, reader, frame, panel, fref)
  ring <- get0("val_backends", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(ring)) ring <- vector("list", 4L)
  bref <- .subset2(backend, "ref")
  slot <- NA_integer_
  for (i in seq_along(ring)) {
    ri <- ring[[i]]
    if (is.null(ri)) {
      if (is.na(slot)) slot <- i
    } else if (is.list(ri) && identical(.subset2(ri[[1L]], "ref"), bref)) {
      slot <- i
      break
    }
  }
  if (is.na(slot)) {
    pos <- get0("val_backends_pos", envir = cdid_fast_kernel_cache, inherits = FALSE)
    pos <- if (is.integer(pos) && length(pos) == 1L) (pos %% 4L) + 1L else 1L
    assign("val_backends_pos", pos, envir = cdid_fast_kernel_cache)
    slot <- pos
  }
  ring[[slot]] <- rec
  assign("val_backends", ring, envir = cdid_fast_kernel_cache)
  invisible(NULL)
}

cdid_fast_validate_link_views <- function(link_views, frame) {
  if (!cdid_are_default_fns(c(
    "validate_ipw_link_view",
    "ipw_owner_seal_ok",
    "validate_ipw_frame_owner",
    "validate_ipw_panel_owner"
  ))) {
    return(FALSE)
  }
  hit <- get0("last_cp04_views", envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (!is.list(hit) || !identical(hit$frame_ref, frame$ref) || !cdid_sexp_aliases(link_views, hit$views)) {
    return(FALSE)
  }
  fk <- cdid_fast_kernel("C_cdid_verify_fast")
  !is.null(fk) && isTRUE(fk(5L, link_views, frame, NULL, NULL))
}

cdid_fast_cp04_public_link_views <- function(frame, specs, aux, control, retained_pipeline_entries) {
  if (!cdid_are_default_fns(c(
    "ipw_owner_link_support_context",
    "cp04_pair_support_context",
    "construct_ipw_link_view",
    "validate_ipw_frame_owner",
    "validate_ipw_panel_owner",
    "validate_ipw_owner_control",
    "preflight_ipw_link_view",
    "ipw_owner_make_pair_descriptor",
    "ipw_owner_make_cohort_descriptor",
    "ipw_owner_context_link_admits",
    "ipw_owner_context_link_lookup",
    "ipw_owner_context_record_link",
    "cp04_pair_support",
    "cp04_pair_support_scan",
    "cp04_pair_support_subset",
    "ipw_owner_chunked_pair_support",
    "ipw_owner_seal"
  ))) {
    return(NULL)
  }
  fk_seal <- cdid_fast_kernel("C_cdid_fast_seal")
  fk_scan <- cdid_fast_kernel("C_cdid_pair_support_scan")
  if (is.null(fk_seal) || is.null(fk_scan)) return(NULL)
  M <- length(specs)
  Q <- length(aux$descriptors)
  if (M + Q < 1L) return(NULL)
  max_cohorts <- 2
  for (s in specs) max_cohorts <- max(max_cohorts, 1 + length(s$K))
  for (s in aux$descriptors) max_cohorts <- max(max_cohorts, 1 + length(s$K))
  ctx <- ipw_owner_link_support_context(frame, M + Q, retained_pipeline_entries, control, max_cohorts = max_cohorts)
  if (is.null(ctx)) return(NULL)
  ctx <- cp04_pair_support_context(ctx, frame, M + Q, retained_pipeline_entries, control, max_cohorts)

  ipw_owner_context_check_frame(ctx, frame)
  ctrl <- validate_ipw_owner_control(control)
  panel_owner <- get("panel_owner", envir = frame$panel_ref, inherits = FALSE)
  n_a <- frame$n_a_units
  p <- frame$n_covariates
  all_times <- panel_owner$all_times
  T_cal <- length(all_times)
  a_cohorts <- frame$a_cohorts
  C_own <- length(a_cohorts)
  post <- frame$postings
  post_dc <- post$date_count
  post_ds <- post$date_start
  post_de <- post$date_end
  post_pa <- post$post_a_index
  aci <- frame$a_cohort_index
  w_a <- frame$w_a
  cb_int <- cdid_interrupt_callback(ipw_owner_check_interrupt)
  ps_env <- ctx$pair_support
  has_ps <- is.environment(ps_env)

  pre_memo <- vector("list", max(32L, as.integer(max_cohorts) + 2L))
  pair_grid <- vector("list", T_cal * T_cal)
  pair_memo <- new.env(hash = TRUE, parent = emptyenv())
  coh_grid1 <- vector("list", C_own * C_own)
  coh_id_grid1 <- integer(C_own * C_own)
  coh_memo <- new.env(hash = TRUE, parent = emptyenv())
  n_coh_cap <- M + Q + 1L
  n_coh_seen <- 0L
  use_grid_cache <- (as.double(n_coh_cap) * as.double(T_cal) * as.double(T_cal)) <= 200000
  link_comp_grid <- if (use_grid_cache) vector("list", n_coh_cap) else NULL
  link_pair_grid <- if (use_grid_cache) vector("list", n_coh_cap * T_cal * T_cal) else NULL
  ps_grid <- if (has_ps && T_cal * T_cal <= 10000L) vector("list", T_cal * T_cal) else NULL
  step0 <- NA_integer_
  fast_links <- new.env(hash = TRUE, parent = emptyenv())
  fast_ps_rec <- if (has_ps) new.env(hash = TRUE, parent = emptyenv()) else NULL

  acct_const <- list(
    per_link_endpoint_copies_retained = FALSE,
    full_n_position_map_retained = FALSE,
    full_n_masks_retained = FALSE,
    endpoints_read = 0L
  )

  build_one <- function(s_val, t_val, g_val, K_val, lid, fkey) {
    if (!is.numeric(s_val) || !is.null(attributes(s_val)) || length(s_val) != 1L || !is.finite(s_val) ||
        !is.numeric(t_val) || !is.null(attributes(t_val)) || length(t_val) != 1L || !is.finite(t_val) ||
        !is.numeric(g_val) || !is.null(attributes(g_val)) || length(g_val) != 1L || !is.finite(g_val) ||
        !is.numeric(K_val) || !is.null(attributes(K_val)) || length(K_val) < 1L || anyNA(K_val) ||
        !gmm_scalar_character(lid) || !gmm_scalar_character(fkey)) {
      return(construct_ipw_link_view(frame, s_val, t_val, g = g_val, K = K_val,
        control = control, link_id = lid, functional_key = fkey, construction_context = ctx))
    }
    m_coh <- max(2L, 1L + length(K_val))
    pre_rec <- if (m_coh <= length(pre_memo)) pre_memo[[m_coh]] else NULL
    if (is.null(pre_rec)) {
      pf <- preflight_ipw_link_view(n_a_units = n_a, p_covariates = p, m_cohorts = m_coh,
        control = ctrl, has_candidates = FALSE)
      adm <- ipw_owner_context_link_admits(ctx, frame, pf, ctrl)
      pre_rec <- list(preflight = pf, admits = adm, step = pf$effective_chunk_rows,
        chunks_out = logit_cert_count_output(if (n_a == 0L) 0L else as.integer(ceiling(as.double(n_a) / as.double(pf$effective_chunk_rows)))))
      if (m_coh <= length(pre_memo)) pre_memo[[m_coh]] <<- pre_rec
      if (is.na(step0)) step0 <<- pf$effective_chunk_rows
    }
    preflight <- pre_rec$preflight
    step <- pre_rec$step
    link_context <- pre_rec$admits

    s_idx <- match(as.double(s_val), all_times)
    t_idx <- match(as.double(t_val), all_times)
    if (!is.na(s_idx) && !is.na(t_idx) && s_idx != t_idx) {
      p_idx <- s_idx + (t_idx - 1L) * T_cal
      pair_desc <- pair_grid[[p_idx]]
      if (is.null(pair_desc)) {
        pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s_val, t_val)
        pair_grid[[p_idx]] <<- pair_desc
      }
    } else {
      p_key <- paste0(sprintf("%a", as.double(s_val)), "|", sprintf("%a", as.double(t_val)))
      pair_desc <- get0(p_key, envir = pair_memo, inherits = FALSE)
      if (is.null(pair_desc)) {
        pair_desc <- ipw_owner_make_pair_descriptor(panel_owner, s_val, t_val)
        assign(p_key, pair_desc, envir = pair_memo)
      }
    }

    cid <- 0L
    g_idx <- match(as.double(g_val), a_cohorts)
    if (length(K_val) == 1L && !is.na(g_idx)) {
      k_idx <- match(as.double(K_val), a_cohorts)
      if (!is.na(k_idx) && g_idx != k_idx) {
        c_idx <- g_idx + (k_idx - 1L) * C_own
        coh_desc <- coh_grid1[[c_idx]]
        if (is.null(coh_desc)) {
          coh_desc <- ipw_owner_make_cohort_descriptor(frame, g_val, K_val, require_support = TRUE)
          coh_grid1[[c_idx]] <<- coh_desc
          n_coh_seen <<- n_coh_seen + 1L
          coh_id_grid1[c_idx] <<- n_coh_seen
        }
        cid <- coh_id_grid1[c_idx]
      } else {
        c_key <- paste0(sprintf("%a", as.double(g_val)), "|", sprintf("%a", as.double(K_val)))
        coh_desc <- get0(c_key, envir = coh_memo, inherits = FALSE)
        if (is.null(coh_desc)) {
          coh_desc <- ipw_owner_make_cohort_descriptor(frame, g_val, K_val, require_support = TRUE)
          assign(c_key, coh_desc, envir = coh_memo)
        }
      }
    } else {
      k_idx <- match(as.double(K_val), a_cohorts)
      if (!is.na(g_idx) && !anyNA(k_idx) && !anyDuplicated(k_idx) && !(g_idx %in% k_idx)) {
        c_key <- paste(c(g_idx, k_idx), collapse = ":")
      } else {
        c_key <- paste0(sprintf("%a", as.double(g_val)), "|", paste(sprintf("%a", as.double(K_val)), collapse = ","))
      }
      coh_desc <- get0(c_key, envir = coh_memo, inherits = FALSE)
      if (is.null(coh_desc)) {
        coh_desc <- ipw_owner_make_cohort_descriptor(frame, g_val, K_val, require_support = TRUE)
        assign(c_key, coh_desc, envir = coh_memo)
      }
    }

    ds <- pair_desc$s_calendar_index
    dt <- pair_desc$t_calendar_index
    valid_ds_dt <- !is.na(ds) && !is.na(dt)
    ds_dt_idx <- if (valid_ds_dt) ds + (dt - 1L) * T_cal else 0L
    counts_full <- valid_ds_dt && post_dc[ds] == n_a && post_dc[dt] == n_a

    sup_info <- NULL
    k_base <- NULL
    k_pair <- NULL
    k_comp <- NULL
    use_fast_slot <- use_grid_cache && cid > 0L && cid <= n_coh_cap && valid_ds_dt && step == step0
    if (link_context) {
      if (use_fast_slot) {
        lp_idx <- cid + (ds_dt_idx - 1L) * n_coh_cap
        rec <- link_pair_grid[[lp_idx]]
        if (is.null(rec)) rec <- link_comp_grid[[cid]]
        if (!is.null(rec) && identical(rec$step, step) && identical(rec$cohorts, coh_desc)) {
          if ((identical(rec$pair, pair_desc) || counts_full) &&
              ipw_owner_pair_masks_equal_fact(frame, rec$pair, pair_desc, step, ctx$dates)) {
            ctx$link_reuses <- ctx$link_reuses + 1
            sup_info <- rec$support
          }
        }
      } else {
        cmap <- coh_desc$a_cohort_map
        k_base <- paste0("ls:", step, ":TRUE:", length(cmap), ":", cmap[1L], ":", cmap[2L], ":", cmap[length(cmap)])
        k_pair <- paste0(k_base, ":pair:", ds, ":", dt)
        k_comp <- paste0(k_base, ":complete")
        for (key in c(k_pair, k_comp)) {
          rec <- get0(key, envir = fast_links, inherits = FALSE)
          if (!is.null(rec) && identical(rec$step, step) && identical(rec$cohorts, coh_desc)) {
            if ((identical(rec$pair, pair_desc) || counts_full) &&
                ipw_owner_pair_masks_equal_fact(frame, rec$pair, pair_desc, step, ctx$dates)) {
              ctx$link_reuses <- ctx$link_reuses + 1
              sup_info <- rec$support
              break
            }
          }
        }
      }
    }
    link_scanned <- is.null(sup_info)
    if (link_scanned && link_context && has_ps && !counts_full && ps_env$n_records <= ps_env$capacity) {
      use_ps_slot <- !is.null(ps_grid) && valid_ds_dt && step == step0
      raw <- if (use_ps_slot) ps_grid[[ds_dt_idx]] else {
        ps_key <- paste0("cp04-pair:", step, ":", ds, ":", dt)
        get0(ps_key, envir = fast_ps_rec, inherits = FALSE)
      }
      if (is.null(raw)) {
        if (ps_env$n_records < ps_env$capacity) {
          has_pair <- valid_ds_dt && post_dc[ds] > 0L && post_dc[dt] > 0L
          bounds4 <- if (has_pair) as.integer(c(post_ds[ds], post_de[ds], post_ds[dt], post_de[dt])) else NULL
          fast_res <- fk_scan(post_pa, bounds4, aci, w_a, NULL, as.integer(c(C_own, n_a)), step, cb_int, environment())
          if (!is.null(fast_res)) {
            raw <- list(pair = pair_desc, step = step, pair_count = fast_res$pair_count,
              pair_mass = fast_res$pair_mass, n_observed = fast_res$n_observed,
              chunks = fast_res$chunks, ops = fast_res$ops)
          } else {
            raw <- cp04_pair_support_scan(frame, pair_desc, step)
          }
          if (use_ps_slot) ps_grid[[ds_dt_idx]] <<- raw else assign(ps_key, raw, envir = fast_ps_rec)
          ps_env$n_records <- ps_env$n_records + 1
          ps_env$scans <- ps_env$scans + 1
          ctx$link_scans <- ctx$link_scans + 1
          last_raw <<- raw; last_coh <<- coh_desc
          sup_info <- cp04_pair_support_subset(raw, frame, coh_desc, TRUE)
          link_scanned <- TRUE
        }
      } else if (identical(raw$pair, pair_desc) && identical(raw$step, step)) {
        ps_env$reuses <- ps_env$reuses + 1
        ctx$link_reuses <- ctx$link_reuses + 1
        last_raw <<- raw; last_coh <<- coh_desc
        sup_info <- cp04_pair_support_subset(raw, frame, coh_desc, TRUE)
        link_scanned <- FALSE
      }
    }
    if (is.null(sup_info)) {
      if (link_context) ctx$link_scans <- ctx$link_scans + 1
      sup_info <- cdid_fast_owner_pair_support(frame, pair_desc, coh_desc, step, TRUE)
      if (is.null(sup_info)) {
        sup_info <- ipw_owner_chunked_pair_support(frame_owner = frame, pair_desc = pair_desc,
          coh_desc = coh_desc, step = step, require_support = TRUE)
      }
    }

    link_view <- structure(
      list(
        link_id = lid,
        functional_key = fkey,
        frame_ref = frame$ref,
        panel_ref = frame$panel_ref,
        nuisance_ref = NULL,
        nuisance_mask_verified = TRUE,
        pair = pair_desc,
        cohorts = coh_desc,
        support = list(
          n_a_units = n_a,
          n_risk_a_units = sup_info$n_risk_a_units,
          n_outside_risk_a_units = sup_info$n_outside_risk_a_units,
          n_observed_pair_all = sup_info$n_observed_pair_all,
          n_observed_risk_units = sup_info$n_observed_risk_units,
          n_nonobserved_risk_units = sup_info$n_nonobserved_risk_units,
          by_cohort_a_count = coh_desc$group_a_count,
          by_cohort_a_mass = coh_desc$group_a_mass,
          by_cohort_pair_count = sup_info$by_cohort_pair_count,
          by_cohort_pair_mass = sup_info$by_cohort_pair_mass,
          q = coh_desc$q,
          rates = sup_info$rates
        ),
        loading_beta_candidate = NULL,
        candidates = NULL,
        control = ctrl,
        preflight = preflight,
        accounting = acct_const,
        work = list(
          support_chunks_processed = if (identical(sup_info$chunks, pre_rec$step)) pre_rec$chunks_out else logit_cert_count_output(sup_info$chunks),
          actual_work_ops = sup_info$ops
        ),
        ref = NULL
      ),
      class = "cdid_ipw_link_view"
    )
    sealed_link <- fk_seal(link_view, "link_view")
    if (is.null(sealed_link)) sealed_link <- ipw_owner_seal(link_view, "link_view")
    if (link_context && link_scanned && ctx$link_records < ctx$link_capacity) {
      complete <- counts_full && ipw_owner_pair_complete(frame, pair_desc, step, ctx$dates)
      if (use_fast_slot) {
        if (complete) {
          if (is.null(link_comp_grid[[cid]])) {
            link_comp_grid[[cid]] <<- list(cohorts = coh_desc, pair = pair_desc, step = step, support = sup_info)
            ctx$link_records <- ctx$link_records + 1
          }
        } else {
          lp_idx <- cid + (ds_dt_idx - 1L) * n_coh_cap
          if (is.null(link_pair_grid[[lp_idx]])) {
            link_pair_grid[[lp_idx]] <<- list(cohorts = coh_desc, pair = pair_desc, step = step, support = sup_info)
            ctx$link_records <- ctx$link_records + 1
          }
        }
      } else {
        rec_key <- if (complete) k_comp else k_pair
        if (!exists(rec_key, envir = fast_links, inherits = FALSE)) {
          assign(rec_key, list(cohorts = coh_desc, pair = pair_desc, step = step, support = sup_info),
            envir = fast_links)
          ctx$link_records <- ctx$link_records + 1
        }
      }
    }
    sealed_link
  }

  last_raw <- NULL
  last_coh <- NULL
  views <- vector("list", M)
  for (j in seq_len(M)) {
    s <- specs[[j]]
    views[[j]] <- build_one(s$s, s$t, s$g, s$K, paste0("cp04:catalog:", j), s$functional_key)
  }
  aux_index <- rep(NA_integer_, Q)
  aux_errors <- vector("list", Q)
  aux_counts <- vector("list", Q)
  err_ps_grid <- if (T_cal * T_cal <= 10000L) vector("list", T_cal * T_cal) else NULL
  for (j in seq_len(Q)) {
    s <- aux$descriptors[[j]]
    last_raw <- NULL
    last_coh <- NULL
    v <- tryCatch(
      build_one(s$s, s$t, s$g, s$K, paste0("cp04:auxiliary:", j),
        paste0("cmpcat:CP04-aux:g", s$g, ":K", paste(s$K, collapse = ","))),
      cdid_ipw_support_error = identity
    )
    if (inherits(v, "error")) {
      aux_errors[[j]] <- v
      if (is.null(last_coh)) {
        last_coh <- tryCatch(
          ipw_owner_make_cohort_descriptor(frame, s$g, s$K, require_support = FALSE),
          error = function(e) NULL
        )
      }
      if (!is.null(last_coh) && is.null(last_raw)) {
        pair_f <- tryCatch(ipw_owner_make_pair_descriptor(panel_owner, s$s, s$t), error = function(e) NULL)
        if (!is.null(pair_f)) {
          ds_f <- pair_f$s_calendar_index
          dt_f <- pair_f$t_calendar_index
          if (!is.na(ds_f) && !is.na(dt_f)) {
            idx_f <- ds_f + (dt_f - 1L) * T_cal
            if (!is.null(ps_grid) && idx_f <= length(ps_grid) && !is.null(ps_grid[[idx_f]])) {
              last_raw <- ps_grid[[idx_f]]
            } else if (!is.null(err_ps_grid) && idx_f <= length(err_ps_grid) && !is.null(err_ps_grid[[idx_f]])) {
              last_raw <- err_ps_grid[[idx_f]]
            } else {
              step_f <- if (!is.na(step0)) step0 else preflight_ipw_link_view(n_a_units = n_a, p_covariates = p,
                m_cohorts = max(2L, 1L + length(s$K)), control = ctrl, has_candidates = FALSE)$effective_chunk_rows
              has_p <- post_dc[ds_f] > 0L && post_dc[dt_f] > 0L
              b4 <- if (has_p) as.integer(c(post_ds[ds_f], post_de[ds_f], post_ds[dt_f], post_de[dt_f])) else NULL
              f_res <- fk_scan(post_pa, b4, aci, w_a, NULL, as.integer(c(C_own, n_a)), step_f, cb_int, environment())
              if (!is.null(f_res)) {
                last_raw <- list(pair = pair_f, step = step_f, pair_count = f_res$pair_count,
                  pair_mass = f_res$pair_mass, n_observed = f_res$n_observed,
                  chunks = f_res$chunks, ops = f_res$ops)
                if (!is.null(err_ps_grid) && idx_f <= length(err_ps_grid)) err_ps_grid[[idx_f]] <- last_raw
              }
            }
          }
        }
      }
      if (!is.null(last_raw) && !is.null(last_coh)) {
        psup <- cp04_pair_support_subset(last_raw, frame, last_coh, FALSE)
        aux_counts[[j]] <- list(
          counts_available = TRUE,
          n_a_units = n_a,
          n_observed_pair_all = psup$n_observed_pair_all,
          n_observed_risk_units = psup$n_observed_risk_units,
          by_cohort_a_count = last_coh$group_a_count,
          by_cohort_pair_count = psup$by_cohort_pair_count
        )
      }
    } else {
      views[[length(views) + 1L]] <- v
      aux_index[j] <- length(views)
    }
  }
  assign("last_cp04_aux", list(frame_ref = frame$ref, aux_errors = aux_errors, counts = aux_counts),
         envir = cdid_fast_kernel_cache)
  work <- list(
    support_scans = ctx$link_scans,
    support_reuses = ctx$link_reuses,
    support_records = ctx$link_records,
    reserved_entries = ctx$link_entries
  )
  if (is.environment(ctx$pair_support)) {
    work$pair_support <- list(
      scans = ctx$pair_support$scans,
      reuses = ctx$pair_support$reuses,
      records = ctx$pair_support$n_records,
      capacity = ctx$pair_support$capacity,
      reserved_entries = ctx$pair_support$entries
    )
  }
  ctx <- NULL
  assign("last_cp04_views", list(views = views, frame_ref = frame$ref), envir = cdid_fast_kernel_cache)
  list(views = views, aux_index = aux_index, aux_errors = aux_errors, construction_work = work)
}

cdid_get_if_chunk_meta <- function(backend, B, links, tag) {
  hit <- get0(tag, envir = cdid_fast_kernel_cache, inherits = FALSE)
  if (is.list(hit) && identical(hit$backend_ref, .subset2(backend, "ref")) &&
      gmm_sexp_aliases(hit$backend, backend) && identical(hit$B, as.integer(B)) &&
      identical(hit$links, as.integer(links))) {
    return(hit)
  }
  NULL
}

cdid_set_if_chunk_meta <- function(backend, B, links, pre, link_names, tag, extra = NULL) {
  rec <- list(
    backend = backend,
    backend_ref = .subset2(backend, "ref"),
    B = as.integer(B),
    links = as.integer(links),
    pre = pre,
    link_names = link_names,
    extra = extra
  )
  assign(tag, rec, envir = cdid_fast_kernel_cache)
  rec
}

cdid_fast_retained_point_cov <- function(carrier, plan, control, mean_phi, n_a, link_ids) {
  if (!cdid_are_default_fns(c(
    "replay_ipw_point_if_retained",
    "ipw_point_if_retained_check",
    "ipw_point_if_retained_plan_verify",
    "custom_ipw_if_plan_verify"
  ))) {
    return(NULL)
  }
  fk_rcov <- cdid_fast_kernel("C_cdid_retained_cov_gram")
  if (is.null(fk_rcov)) return(NULL)
  ipw_point_if_retained_check(carrier, plan, control, expected_state = "filled")
  if (!identical(carrier$state, "filled")) return(NULL)
  ctrl <- plan$control
  pre <- plan$preflight
  N <- carrier$N
  L <- carrier$L
  B <- carrier$B
  ipw_owner_safe_sum(c(pre$total_work_ops_bound, carrier$work), ctrl$max_work_ops, "CP08 retained point IF slice work")
  V <- fk_rcov(carrier$phi, as.double(mean_phi), B, as.double(n_a), 1L)
  if (is.null(V)) return(NULL)
  dimnames(V) <- list(link_ids, link_ids)
  batches <- as.integer(ceiling(as.double(N) / as.double(B)))
  max_b <- min(as.integer(N), as.integer(B))
  carrier$verifications <- carrier$verifications + 2L * batches
  carrier$slice_passes <- carrier$slice_passes + 1L
  carrier$slices_emitted <- carrier$slices_emitted + batches
  replay_info <- list(
    n_a_units = as.integer(N),
    n_links = as.integer(L),
    link_ids = plan$link_ids,
    n_batches = batches,
    max_batch_units = max_b,
    max_batch_entries = as.double(max_b) * L,
    actual_work_ops = 0,
    kernel_calls = 0L,
    normalizations_performed = 0L,
    verifications = 1L + 2L * batches,
    state_verifications = 0L,
    served_from = "retained_point_if",
    preflight = pre
  )
  list(cov_mat = V, replay_info = replay_info)
}

cdid_fast_project_retained_if <- function(if_replay, Lambda, m_R, bounds, materialize_target_if = FALSE) {
  ret <- attr(if_replay$callback, "cdid_retained_if", exact = TRUE)
  if (is.null(ret) || !is.list(ret) ||
      !identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !identical(parent.env(environment(if_replay$callback)), environment(cdid_fast_kernel))) {
    return(NULL)
  }
  J <- ncol(Lambda)
  if (!gmm_if_unit_batch_eligible(if_replay, J)) return(NULL)
  ni <- as.integer(if_replay$n_input_rows)
  no <- as.integer(if_replay$n_output_rows)
  if (is.na(ni) || is.na(no) || ni < 1L || ni != no ||
      !identical(as.integer(if_replay$n_positive_contributing_rows), ni)) {
    return(NULL)
  }
  fk_wb <- cdid_fast_kernel("C_cdid_if_unit_batch_weighted")
  fk_id <- cdid_fast_kernel("C_cdid_num_id_tracker")
  if (is.null(fk_wb) || is.null(fk_id)) return(NULL)

  if (identical(ret$kind, "cp04")) {
    if (!cdid_are_default_fns(c(
      "cp04_replay_links",
      "cp04_retained_if_replay",
      "cp04_retained_if_check",
      cp04_retained_if_default_validators
    ))) {
      return(NULL)
    }
    fk_vf <- cdid_fast_kernel("C_cdid_verify_fast")
    if (is.null(fk_vf)) return(NULL)
    producer <- ret$producer
    if (!inherits(producer, "cdid_cp04_retained_if_bound")) return(NULL)
    store <- producer$store
    backend <- producer$backend
    links <- as.integer(ret$links)
    B <- as.integer(ret$chunk_rows)
    if (!is.environment(store) || is.na(B) || B < 1L || B > bounds$max_batch_rows ||
        !identical(store$state, "filled") || !identical(as.integer(store$N), ni) ||
        !identical(as.integer(store$L), as.integer(m_R))) {
      return(NULL)
    }
    slice_limit <- if (is.null(store$gmm_transfer)) store$max_slice_passes else store$gmm_transfer$max_slice_passes
    if (!is.null(slice_limit) && store$slice_passes + 2L > slice_limit) return(NULL)
    st <- fk_vf(3L, store, backend, links, list(B, "core", 1L))
    if (is.null(st)) return(NULL)
    if (!is.numeric(store$work) || length(store$work) != 1L || !is.finite(store$work) ||
        2 * store$work > backend$control$max_work_ops) {
      pre <- cp04_preflight_links(backend, links, chunk_rows = B)
      ipw_owner_safe_sum(c(pre$work_ops_bound, store$work), backend$control$max_work_ops, "CP04 retained moment IF slice work")
    }
    phi_full <- store$phi
    uid <- st$id
  } else if (identical(ret$kind, "point")) {
    if (!cdid_are_default_fns(c(
      "replay_ipw_point_if_retained",
      "ipw_point_if_retained_check",
      "ipw_point_if_retained_plan_verify",
      "custom_ipw_if_plan_verify",
      "replay_custom_ipw_if_plan",
      "replay_ipw_point_if_general"
    ))) {
      return(NULL)
    }
    carrier <- ret$carrier
    plan <- ret$plan
    control <- ret$control
    if (!is.environment(carrier) || !is.list(plan) || !identical(carrier$state, "filled")) return(NULL)
    B <- as.integer(carrier$B)
    if (is.na(B) || B < 1L || B > bounds$max_batch_rows ||
        !identical(as.integer(carrier$N), ni) || !identical(as.integer(carrier$L), as.integer(m_R))) {
      return(NULL)
    }
    ipw_point_if_retained_check(carrier, plan, control, expected_state = "filled")
    ipw_owner_safe_sum(c(plan$preflight$total_work_ops_bound, carrier$work), plan$control$max_work_ops, "CP08 retained point IF slice work")
    phi_full <- carrier$phi
    uid <- plan$panel_owner$a_id
  } else {
    return(NULL)
  }

  if (!is.null(if_replay$expected_unit_ids) && !identical(if_replay$expected_unit_ids, uid)) return(NULL)
  if (!is.null(if_replay$expected_cluster_ids) && !identical(if_replay$expected_cluster_ids, uid)) return(NULL)

  id_res <- fk_id(5L, NULL, uid, ni, isTRUE(materialize_target_if))
  if (is.null(id_res)) return(NULL)

  res_c <- fk_wb(3L, Lambda, phi_full, isTRUE(if_replay$centered), as.integer(B), NULL, NULL, NULL, NULL)
  if (is.null(res_c)) return(NULL)
  psi <- res_c$psi
  lam_cn <- colnames(Lambda)
  S <- res_c$S
  dimnames(S) <- list(lam_cn, lam_cn)
  sandwich <- if_replay$normalization_factor * S
  if (any(!is.finite(sandwich)) || any(diag(sandwich) < 0)) return(NULL)

  if (isTRUE(materialize_target_if)) {
    dimnames(psi) <- list(id_res, lam_cn)
    materialized <- psi
  } else {
    materialized <- NULL
  }
  psi_mean_out <- res_c$psi_mean
  raw_sum_out <- res_c$raw_sum
  if (!is.null(lam_cn)) {
    names(psi_mean_out) <- lam_cn
    names(raw_sum_out) <- lam_cn
  }

  nb1 <- as.integer(ceiling(as.double(ni) / as.double(B)))
  max_B <- min(ni, B)
  if (identical(ret$kind, "cp04")) {
    store$slice_passes <- store$slice_passes + 2L
  } else {
    carrier$verifications <- carrier$verifications + 1L + 4L * nb1
    carrier$slice_passes <- carrier$slice_passes + 2L
    carrier$slices_emitted <- carrier$slices_emitted + 2L * nb1
  }

  list(
    sandwich_cov = sandwich,
    se_sandwich = sqrt(diag(sandwich)),
    psi_weighted_mean = psi_mean_out,
    psi_raw_sum = raw_sum_out,
    total_weight = res_c$sum_w,
    n_batches = nb1,
    max_batch_rows_seen = max_B,
    n_input_rows = if_replay$n_input_rows,
    n_output_rows = if_replay$n_output_rows,
    n_positive_contributing_rows = ni,
    materialized_target_if = materialized,
    replay_verified_scope = "Exact stored projected rows and complete IDs/weights/counts, across two measured passes. Raw IF/nuisance immutability and sampling independence remain caller assumptions.",
    sampling_provenance_certified = FALSE
  )
}

cdid_fast_df <- function(cols) {
  n <- if (length(cols)) length(cols[[1L]]) else 0L
  structure(cols, class = "data.frame", row.names = .set_row_names(n))
}

cdid_fast_eval_full_if <- function(backend, links, scale, if_plan) {
  if (!is.list(if_plan)) return(NULL)
  links <- as.integer(links)
  kind <- if_plan$kind
  if (inherits(backend, "cdid_or_baseline") && identical(kind, "or_sparse_native_plan")) {
    if (!cdid_are_default_fns(c(
      "or_sparse_native_if", "or_sparse_native_check", "or_baseline_native_check",
      "or_sparse_native_slice", "or_sparse_native_live", "or_baseline_native_kernel",
      "cp04_payload"
    ))) {
      return(NULL)
    }
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    N <- st$N
    pre_b <- cdid_fast_if_budget(backend, if_plan$B, links, "or")
    or_sparse_native_check(if_plan, backend, links, if_plan$B, pre_b)
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
    live_N <- or_sparse_native_live(
      N, length(links), length(if_plan$dates), length(if_plan$sids),
      ncol(if_plan$data$block), length(if_plan$data$beta), max(c(0L, if_plan$data$block[3L, ]))
    )
    if (!is.finite(live_N) || live_N > backend$control$max_matrix_entries) return(NULL)
    fn_if <- or_baseline_native_kernel("C_gmm_or_sparse_if")
    res <- fn_if(
      panel$observations$y, panel$postings$post_a_index, panel$postings$post_obs_index,
      seq_len(N), st$cohort_index, st$a_index, st$w, if_plan$bounds, frame$Z_a,
      if_plan$data, as.double(live_N)
    )
    if (!identical(res$status[1L], 0L)) return(NULL)
    phi <- res$phi
    if (scale == "core" && st$lambda != 1) phi <- st$lambda * phi
    if (!is.finite(sum(phi))) return(NULL)
    dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1L), "link_id"))
    return(phi)
  }

  if (inherits(backend, "cdid_dr_baseline") && identical(kind, "dr_sparse_native_plan")) {
    if (!cdid_are_default_fns(c(
      "dr_sparse_native_if", "dr_sparse_native_check", "or_baseline_native_check",
      "dr_sparse_native_inventory", "or_baseline_native_kernel",
      "or_pooled_wls_gram_solve", "or_wls_gram_solve", "cp04_solve",
      "cp04_payload"
    ))) {
      return(NULL)
    }
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    N <- st$N
    pre_b <- cdid_fast_if_budget(backend, if_plan$B, links, "dr")
    dr_sparse_native_check(if_plan, backend, links, if_plan$B, pre_b)
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
    p <- backend$shared$p
    pt <- backend$treatment_recipe$p
    nt <- length(if_plan$tids)
    inv_N <- dr_sparse_native_inventory(
      N, length(if_plan$sids), length(if_plan$dates), nt, length(links), p, pt, if_plan$metadata
    )
    if (!is.finite(inv_N$prepare) || !is.finite(inv_N$assembler) ||
        inv_N$prepare > backend$control$max_matrix_entries ||
        inv_N$assembler > backend$control$max_matrix_entries) {
      return(NULL)
    }
    fn_prep <- or_baseline_native_kernel("C_gmm_dr_sparse_prepare")
    fn_asm <- or_baseline_native_kernel("C_gmm_dr_sparse_assemble")
    a <- fn_prep(
      panel$observations$y, panel$postings$post_a_index, panel$postings$post_obs_index,
      seq_len(N), st$cohort_index, st$a_index, st$w, if_plan$bounds, frame$Z_a,
      if_plan$data, as.double(inv_N$prepare)
    )
    if (!identical(a$status[1L], 0L)) return(NULL)
    xi_beta <- if (p) {
      sfit <- backend$shared$fit
      if (cdid_are_default_fns(c("or_pooled_wls_gram_solve", "or_wls_gram_solve")) &&
          is.list(sfit) && is.matrix(sfit$R) && nrow(sfit$R) == p && ncol(sfit$R) == p &&
          is.integer(sfit$pivot) && length(sfit$pivot) == p && all(abs(diag(sfit$R)) > 0) &&
          is.finite(sum(a$or_score))) {
        piv <- sfit$pivot
        R_mat <- sfit$R
        Z_sol <- forwardsolve(R_mat, t(a$or_score)[piv, , drop = FALSE], upper.tri = TRUE, transpose = TRUE)
        X_piv <- backsolve(R_mat, Z_sol)
        X_out <- matrix(0.0, N, p)
        X_out[, piv] <- st$N * t(X_piv)
        if (!is.finite(sum(X_out))) {
          X_out <- t(if (backend$outcome_recipe$intercept) {
            or_pooled_wls_gram_solve(sfit, t(a$or_score), normalizer = st$N)
          } else {
            or_wls_gram_solve(sfit, t(a$or_score), normalizer = st$N)
          })
        }
        X_out
      } else {
        t(if (backend$outcome_recipe$intercept) {
          or_pooled_wls_gram_solve(sfit, t(a$or_score), normalizer = st$N)
        } else {
          or_wls_gram_solve(sfit, t(a$or_score), normalizer = st$N)
        })
      }
    } else {
      matrix(0, N, 0L)
    }
    fk_chol <- cdid_fast_kernel("C_cdid_dr_batch_chol_solve")
    xi <- if (!is.null(fk_chol) && pt > 0L && nt > 0L) {
      fk_chol(backend$treatments, as.integer(if_plan$tids), a$rhs, as.integer(c(pt, nt)))
    } else {
      NULL
    }
    if (is.null(xi)) {
      xi <- vector("list", nt)
      for (j in seq_len(nt)) {
        ix <- (j - 1L) * pt + seq_len(pt)
        xi[[j]] <- t(cp04_solve(backend$treatments[[if_plan$tids[j]]]$H_chol, t(a$rhs[, ix, drop = FALSE])))
      }
    }
    b <- fn_asm(a$up, xi_beta, a$direct, xi, if_plan$data, as.integer(nt), as.double(inv_N$assembler))
    if (!identical(b$status[1L], 0L)) return(NULL)
    phi <- b$phi
    if (scale == "core" && st$lambda != 1) phi <- st$lambda * phi
    if (any(!is.finite(phi))) return(NULL)
    dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1L), "link_id"))
    return(phi)
  }

  if (inherits(backend, "cdid_dr_baseline") && identical(kind, "dr_separate_sparse_plan")) {
    if (!cdid_are_default_fns(c(
      "dr_separate_sparse_if", "dr_separate_sparse_check", "or_baseline_native_check",
      "dr_separate_baseline_replay_identity", "dr_separate_sparse_inventory",
      "dr_separate_sparse_live", "or_baseline_native_kernel",
      "cp04_finite", "cp04_payload"
    ))) {
      return(NULL)
    }
    reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
    st <- reader$state
    N <- st$N
    pre_b <- cdid_fast_if_budget(backend, if_plan$B, links, "dr")
    dr_separate_sparse_check(if_plan, backend, links, if_plan$B, pre_b)
    frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
    live_N <- dr_separate_sparse_live(
      N, length(if_plan$sids), length(if_plan$dates), length(if_plan$tids), length(links)
    )
    if (!is.finite(live_N) || live_N > backend$control$max_matrix_entries) return(NULL)
    fn_if <- or_baseline_native_kernel("C_gmm_dr_separate_sparse_if")
    res <- fn_if(
      panel$observations$y, panel$postings$post_a_index, panel$postings$post_obs_index,
      seq_len(N), st$cohort_index, st$a_index, st$w, if_plan$bounds, frame$Z_a,
      if_plan$data, if_plan$v11, as.double(live_N)
    )
    if (!identical(res$status[1L], 0L)) return(NULL)
    phi <- res$phi
    if (scale == "core" && st$lambda != 1) phi <- st$lambda * phi
    if (!is.finite(sum(phi))) return(NULL)
    dimnames(phi) <- list(NULL, vapply(backend$links[links], `[[`, character(1L), "link_id"))
    return(phi)
  }

  if (inherits(backend, "cdid_ipw_raw_backend") && identical(kind, "ipw_raw_if_plan")) {
    if (!cdid_are_default_fns(c(
      "evaluate_ipw_raw_if_chunk", "validate_ipw_raw_backend",
      "ipw_owner_read_y_positions", "ipw_owner_query_pair_mask",
      "cp04_solve", "cp04_finite"
    ))) {
      return(NULL)
    }
    fk_if <- cdid_fast_kernel("C_cdid_ipw_raw_if_chunk")
    if (is.null(fk_if)) return(NULL)
    validate_ipw_raw_backend(backend)
    n_a <- backend$n_a_units
    if (as.double(n_a) * length(links) > backend$control$max_matrix_entries) return(NULL)
    plan_ref <- .subset2(if_plan, "ref")
    if (!is.environment(plan_ref) || !ipw_owner_seal_ok(if_plan, "ipw_raw_if_plan") ||
        !identical(if_plan$backend_ref, backend$ref) || !identical(if_plan$links, links)) {
      return(NULL)
    }
    pre_r <- if_plan$pre
    replay_index <- if_plan$replay_index
    plan_entries <- if_plan$entries - if (is.null(replay_index)) 0 else if_plan$funded_entries
    ipw_owner_safe_sum(
      c(backend$accounting$owned_entries, backend$accounting$borrowed_entries,
        pre_r$output, pre_r$scratch, plan_entries),
      backend$control$max_matrix_entries,
      "raw-IPW sealed plan plus IF output and scratch"
    )
    meta <- if (is.null(replay_index)) ipw_raw_if_metadata(backend, pre_r) else replay_index
    frame <- get("frame_owner", envir = backend$frame_ref, inherits = FALSE)
    panel <- get("panel_owner", envir = backend$panel_ref, inherits = FALSE)
    pack <- cdid_fast_ipw_raw_if_pack(
      backend, pre_r, meta, replay_index, if_plan$model_pos, frame, panel, plan_ref = plan_ref
    )
    if (is.null(pack)) return(NULL)
    fast_res <- fk_if(
      seq_len(n_a),
      panel$postings$post_a_index,
      panel$postings$post_obs_index,
      panel$observations$y,
      frame$a_cohort_index,
      frame$w_a,
      frame$Z_a,
      pack$pair_tab,
      pack$pair_active_coh,
      pack$trt_int,
      pack$trt_dbl_list
    )
    if (is.null(fast_res)) return(NULL)
    phi <- fast_res$phi
    if (!is.finite(sum(phi))) return(NULL)
    dimnames(phi) <- list(NULL, if_plan$link_ids)
    return(phi)
  }

  NULL
}

cdid_fast_dr_sep_env_builder <- function(backend, reader, frame, links_i) {
  st <- reader$state
  ctrl <- backend$control
  own <- backend$accounting$owned_entries + backend$accounting$borrowed_entries
  reg_specs <- backend$regression_specs
  trt_specs <- backend$treatment_specs
  G <- length(reg_specs)
  S <- length(backend$supports)
  Tn <- length(backend$treatments)
  p <- backend$outcome_recipe$p
  pt <- backend$treatment_recipe$p
  J <- length(st$cohorts)
  nl <- length(links_i)
  if (nl < 1L) return(NULL)
  incidence <- 0
  for (i in seq_len(G)) incidence <- incidence + length(reg_specs[[i]]$support_ids)
  for (i in seq_len(Tn)) incidence <- incidence + length(trt_specs[[i]]$support_ids)
  maps <- 4096 + 32 * (G + S + S + Tn + nl + incidence)
  if (own > ctrl$max_retained_entries || own + maps > ctrl$max_matrix_entries ||
      backend$accounting$work_ops_bound + 64 * maps > ctrl$max_work_ops) {
    return(NULL)
  }
  b_links <- backend$links
  link_blocks <- integer(nl)
  link_sids <- integer(nl)
  link_tids <- integer(nl)
  link_bytes <- 0L
  for (i in seq_len(nl)) {
    li <- b_links[[links_i[i]]]
    link_blocks[i] <- li$regression
    link_sids[i] <- li$support_id
    link_tids[i] <- li$model
    link_bytes <- link_bytes + nchar(li$link_id, type = "bytes")
  }
  gids <- sort.int(unique(link_blocks))
  tids <- sort.int(unique(link_tids))
  ng <- length(gids)
  nt <- length(tids)
  dr_sids <- link_sids
  for (i in seq_len(ng)) dr_sids <- c(dr_sids, reg_specs[[gids[i]]]$support_ids)
  for (i in seq_len(nt)) dr_sids <- c(dr_sids, trt_specs[[tids[i]]]$support_ids)
  dr_sids <- sort.int(unique(dr_sids))
  ns <- length(dr_sids)
  nu <- ns
  og_groups <- backend$outcome_groups$groups
  ps <- integer(ng)
  max_p <- 0L
  for (i in seq_len(ng)) {
    pk <- length(og_groups[[gids[i]]]$keep)
    ps[i] <- pk
    if (pk > max_p) max_p <- pk
  }
  model_u <- integer(nt)
  max_u <- 0L
  for (i in seq_len(nt)) {
    uk <- length(trt_specs[[tids[i]]]$support_ids)
    model_u[i] <- uk
    if (uk > max_u) max_u <- uk
  }
  reference <- chained_backend_reference(backend)
  p_ref <- if (is.null(reference)) 0L else length(reference$x_columns)
  id_max <- frame$id_byte_metadata$max
  out_meta_0 <- 8192 + 128 * link_bytes + 128 * nl
  width <- nl + 8
  log_nA <- 1 + ceiling(log2(st$n_A + 1))
  scr_per_B <- 32 * (ns * (p + 6) + nu + max_p * (max_p + max_u + nl + 8) + 4 * pt + max_u + p_ref + nl + 32)
  scr_const <- 128 * ((max_p + 2)^2 + (pt + 2)^2) + 64 * (ns + nu + ng + nt)
  wrk_0 <- 128 * maps
  wrk_slope <- 128 * ((ns + 1) * (p + 2)^2 * log_nA + sum((ps + 2)^2 * (nu + nl + 2)) +
                      nt * (pt + 2) * (pt + max_u + 8) + (ns + nt + 1) * p_ref + J + nl * (p + pt + 8))

  function(B_val, check_dr_budget = FALSE) {
    out <- out_meta_0 + B_val * (id_max + width)
    scr <- maps + out + B_val * scr_per_B + scr_const
    wrk <- wrk_0 + B_val * wrk_slope
    if (own + out + scr > ctrl$max_matrix_entries ||
        backend$accounting$work_ops_bound + wrk > ctrl$max_work_ops ||
        out > ctrl$max_output_entries) {
      return(NULL)
    }
    list(
      gids = gids, beta_gids = gids, tids = tids, sids = dr_sids, uids = dr_sids,
      output = out, scratch = scr, work = wrk
    )
  }
}

cdid_fast_if_budget <- function(backend, B, links, mode = "full") {
  eval_raw <- function() {
    if (identical(mode, "dr_env12")) {
      return(list(one = dr_baseline_if_envelope(backend, 1, links), two = dr_baseline_if_envelope(backend, 2, links)))
    }
    if (identical(mode, "dr")) return(dr_baseline_if_budget(backend, B, links, FALSE))
    if (identical(mode, "or")) return(or_baseline_if_budget(backend, B, links, FALSE))
    chained_full_if_budget(backend, B, links)
  }
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(eval_raw())
  if (!inherits(backend, c("cdid_dr_baseline", "cdid_or_baseline")) ||
      !is.null(backend$reference_ref) || !is.null(backend$selection_backend) ||
      !is.environment(backend$reader_ref) ||
      !cdid_are_default_fns(c(
        "dr_baseline_if_envelope", "dr_baseline_if_budget", "dr_separate_if_envelope",
        "or_baseline_if_budget", "or_baseline_separate_if_budget", "chained_full_if_budget", "cp04_budget"
      ))) {
    return(eval_raw())
  }
  reader <- get0("cp04_reader", envir = backend$reader_ref, inherits = FALSE)
  if (is.null(reader) || !is.environment(reader$frame_ref) || !is.environment(reader$panel_ref)) return(eval_raw())
  frame <- get0("frame_owner", envir = reader$frame_ref, inherits = FALSE)
  panel <- get0("panel_owner", envir = reader$panel_ref, inherits = FALSE)
  if (is.null(frame) || is.null(panel)) return(eval_raw())
  eff_mode <- if (identical(mode, "dr_env12")) "dr_env12" else if (inherits(backend, "cdid_dr_baseline")) "dr" else "or"
  links_i <- as.integer(links)
  B_d <- if (identical(eff_mode, "dr_env12")) 0.0 else as.double(B)
  ring <- cdid_fast_kernel_cache$if_budget_ring
  cached_builder <- NULL
  if (is.list(ring)) {
    for (entry in ring) {
      if (identical(entry$links, links_i) &&
          cdid_sexp_aliases(backend, entry$backend) &&
          cdid_sexp_aliases(reader, entry$reader) &&
          cdid_sexp_aliases(frame, entry$frame) &&
          cdid_sexp_aliases(panel, entry$panel)) {
        if (identical(entry$mode, eff_mode) && identical(entry$B, B_d)) {
          return(entry$val)
        }
        if (is.function(entry$builder)) cached_builder <- entry$builder
      }
    }
  }
  val <- NULL
  if (inherits(backend, "cdid_dr_baseline") && dr_separate(backend) &&
      !isTRUE(backend$uses_selection_models) && is.null(backend$separate_complete)) {
    if (is.null(cached_builder)) {
      cached_builder <- cdid_fast_dr_sep_env_builder(backend, reader, frame, links_i)
    }
    if (is.function(cached_builder)) {
      if (identical(eff_mode, "dr_env12")) {
        one_v <- cached_builder(1.0, FALSE)
        two_v <- if (!is.null(one_v)) cached_builder(2.0, FALSE) else NULL
        if (!is.null(two_v)) val <- list(one = one_v, two = two_v)
      } else if (identical(eff_mode, "dr")) {
        val <- cached_builder(B_d, TRUE)
      }
    }
  }
  if (is.null(val)) val <- eval_raw()
  new_entry <- list(
    mode = eff_mode, B = B_d, links = links_i,
    backend = backend, reader = reader, frame = frame, panel = panel,
    val = val, builder = cached_builder
  )
  cdid_fast_kernel_cache$if_budget_ring <- if (is.list(ring)) c(list(new_entry), utils::head(ring, 5L)) else list(new_entry)
  val
}

cdid_fast_stream_covariance <- function(store, backend, links, plan, st, ctrl) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "cp04_replay_links", "cp04_retained_if_replay", "cp04_retained_if_check",
    "chained_full_if_plan", "chained_full_if_chunk",
    cp04_retained_if_default_validators
  ))) {
    return(NULL)
  }
  fk_vf <- cdid_fast_kernel("C_cdid_verify_fast")
  fk_rcov <- cdid_fast_kernel("C_cdid_retained_cov_gram")
  if (is.null(fk_vf) || is.null(fk_rcov)) return(NULL)
  B <- plan$chunk_rows
  L <- plan$L
  N <- st$N
  if (!is.null(store)) {
    if (!identical(store$state, "reserved")) return(NULL)
    slice_limit <- if (is.null(store$gmm_transfer)) store$max_slice_passes else store$gmm_transfer$max_slice_passes
    if (!is.null(slice_limit) && slice_limit < 1L) return(NULL)
    st_chk <- fk_vf(3L, store, backend, links, list(B, "core", 0L))
    if (is.null(st_chk)) return(NULL)
    pre <- cp04_preflight_links(backend, links, chunk_rows = B)
    ipw_owner_safe_sum(c(pre$work_ops_bound, store$work), ctrl$max_work_ops, "CP04 retained moment IF fill work")
    ipw_owner_check_matrix_dims(N, L, ctrl$max_matrix_entries, "CP04 retained moment IF")
  } else {
    if (as.double(N) * as.double(L) > ctrl$max_matrix_entries) return(NULL)
    if (!cdid_are_default_fns(c("cp04_preflight_links", "cp04_plan_rows", "chained_full_if_budget"))) {
      cp04_preflight_links(backend, links, chunk_rows = B)
    }
  }
  if_plan <- chained_full_if_plan(backend, links, B, phase_store = store)
  on.exit(augmented_adjoint_replay_release(if_plan), add = TRUE)
  phi_full <- cdid_fast_eval_full_if(backend, links, "core", if_plan)
  if (is.null(phi_full)) return(NULL)
  rcov <- fk_rcov(phi_full, NULL, B, as.double(st$n_A), 2L)
  if (is.null(rcov)) return(NULL)
  if (!is.null(store)) {
    raw_ref <- new.env(parent = emptyenv())
    raw_ref$phi <- phi_full
    lockEnvironment(raw_ref, bindings = FALSE)
    store$phi <- phi_full
    store$raw_ref <- raw_ref
    store$state <- "filled"
    store$producer_passes <- store$producer_passes + 1L
    store$batches <- as.integer(ceiling(as.double(N) / as.double(B)))
    store$slice_passes <- store$slice_passes + 1L
  }
  cdid_fast_kernel_cache$last_retained_if <- list(
    backend_ref = backend$ref,
    links = as.integer(links),
    phi_core = phi_full,
    lambda = as.double(st$lambda),
    N = as.integer(N)
  )
  rcov
}

cdid_fast_aux_replay <- function(backend, ids, col_dst_all, col_src_all, signs_all, out, chunk_rows) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(FALSE)
  if (!cdid_are_default_fns(c(
    "cp04_replay_links", "chained_full_if_plan", "chained_full_if_chunk",
    cp04_retained_if_default_validators
  ))) {
    return(FALSE)
  }
  fk_set_aux <- cdid_fast_kernel("C_cdid_matrix_set_rows")
  if (is.null(fk_set_aux)) return(FALSE)
  bound <- cp04_retained_if_unbind(backend)
  if (!is.null(bound$store)) return(FALSE)
  dom <- cp04_backend_domain(backend)
  if (as.double(dom$N) * length(ids) > backend$control$max_matrix_entries) return(FALSE)
  B <- if (!is.null(chunk_rows) && cdid_are_default_fns(c("cp04_preflight_links", "cp04_plan_rows"))) {
    as.integer(chunk_rows)
  } else {
    cp04_preflight_links(backend, ids, chunk_rows = chunk_rows)$chunk_rows
  }

  ids_i <- as.integer(ids)
  last <- cdid_fast_kernel_cache$last_retained_if
  if (!is.null(last) && identical(last$backend_ref, backend$ref) &&
      identical(last$N, as.integer(dom$N)) && is.finite(last$lambda) && last$lambda > 0) {
    m_in_last <- match(ids_i, last$links)
    if (anyNA(m_in_last) &&
        inherits(backend, c("cdid_or_baseline", "cdid_dr_baseline", "cdid_ipw_raw_backend")) &&
        is.null(backend$reference_ref) &&
        cdid_are_default_fns(c("cp04_retained_if_auxiliary_map", "cp04_retained_if_recipe"))) {
      na_pos <- which(is.na(m_in_last))
      tmp_store <- list2env(list(kind = "cp04_retained_if", links = last$links), parent = emptyenv())
      m_in_last[na_pos] <- cp04_retained_if_auxiliary_map(tmp_store, backend, ids_i[na_pos])
    }
    has_cached <- !is.na(m_in_last)
    if (any(has_cached)) {
      ids_rem <- ids_i[!has_cached]
      phi_rem <- NULL
      if (length(ids_rem) > 0L) {
        plan_rem <- chained_full_if_plan(backend, ids_rem, B)
        on.exit(augmented_adjoint_replay_release(plan_rem), add = TRUE)
        phi_rem <- cdid_fast_eval_full_if(backend, ids_rem, "frame", plan_rem)
      }
      if (length(ids_rem) == 0L || !is.null(phi_rem)) {
        in_cached_req <- has_cached[col_src_all]
        if (any(in_cached_req)) {
          idx_c <- which(in_cached_req)
          src_c <- as.integer(m_in_last[col_src_all[idx_c]])
          dst_c <- as.integer(col_dst_all[idx_c])
          sgn_c <- as.double(signs_all[idx_c]) / last$lambda
          if (!isTRUE(fk_set_aux(out, list(NULL, dst_c, src_c, sgn_c), last$phi_core))) return(FALSE)
        }
        if (any(!in_cached_req)) {
          idx_r <- which(!in_cached_req)
          rem_map <- integer(length(ids_i))
          rem_map[!has_cached] <- seq_along(ids_rem)
          src_r <- as.integer(rem_map[col_src_all[idx_r]])
          dst_r <- as.integer(col_dst_all[idx_r])
          sgn_r <- as.double(signs_all[idx_r])
          if (!isTRUE(fk_set_aux(out, list(NULL, dst_r, src_r, sgn_r), phi_rem))) return(FALSE)
        }
        return(TRUE)
      }
    }
  }

  if_plan <- chained_full_if_plan(backend, ids_i, B)
  on.exit(augmented_adjoint_replay_release(if_plan), add = TRUE)
  phi_all <- cdid_fast_eval_full_if(backend, ids_i, "frame", if_plan)
  if (is.null(phi_all)) return(FALSE)
  isTRUE(fk_set_aux(out, list(NULL, as.integer(col_dst_all), as.integer(col_src_all), as.double(signs_all)), phi_all))
}

cdid_fast_panel_assemble <- function(owner_id, ctrl, preflight, cols, n, id_is_num, id_storage, id_byte_metadata,
                                     ids_key, nu, ui, time_ord, y_ord, g_ord, validation_chunks, actual_ops) {
  if (!cdid_are_default_fns(c(
    "ipw_owner_match_sorted_double", "pair_count_safe",
    "ipw_owner_safe_sum", "logit_cert_count_output", "ipw_owner_seal"
  ))) {
    return(NULL)
  }
  fk_vf <- cdid_fast_kernel("C_cdid_verify_fast")
  if (is.null(fk_vf)) return(NULL)
  last_prep <- cdid_fast_kernel_cache$last_panel_prep
  all_times <- if (!is.null(last_prep) && cdid_sexp_aliases(time_ord, last_prep$tt_d)) {
    last_prep$calendar
  } else {
    sort(unique(time_ord), method = "radix")
  }
  n_cal <- length(all_times)
  if (n_cal < 1L || n_cal > .Machine$integer.max) return(NULL)

  c_res <- fk_vf(11L, ui, time_ord, y_ord, list(g_ord, as.integer(nu), all_times))
  if (is.null(c_res)) return(NULL)
  if (!is.finite(c_res$total_pairs) || c_res$total_pairs > 9007199254740991) return(NULL)

  ids <- if (identical(id_storage, "integer")) as.integer(ids_key) else ids_key
  ot <- c_res$ot
  ou <- c_res$ou
  oy <- c_res$oy
  ocal <- c_res$ocal
  na_u <- c_res$na_u
  na_t <- c_res$na_t
  na_cal <- c_res$na_cal
  nr <- c_res$nr
  no <- c_res$no
  n_na_u <- c_res$n_na_u
  ug <- c_res$ug
  first <- c_res$first
  last <- c_res$last
  sy <- c_res$sy
  singleton_time <- c_res$singleton_time
  single <- c_res$single
  unobs <- c_res$unobs
  a_unit_index <- c_res$a_unit_index
  n_a <- length(a_unit_index)
  unit_to_a <- c_res$unit_to_a
  obs_a_index <- c_res$obs_a_index
  q_a <- c_res$q_a
  actual_ops <- actual_ops + 2 * length(ot) + 2 * q_a

  a_id <- ids[a_unit_index]
  a_id_key <- ids_key[a_unit_index]
  a_g <- as.double(ug[a_unit_index])
  a_cohorts <- if (n_a > 0L) sort(unique(a_g), method = "radix") else numeric(0L)
  a_cohort_index <- if (n_a > 0L) ipw_owner_match_sorted_double(a_g, a_cohorts) else integer(0L)
  cohort_a_count <- if (length(a_cohorts) > 0L) tabulate(a_cohort_index, nbins = length(a_cohorts)) else integer(0L)

  units <- cdid_fast_df(list(
    unit_index = seq_len(nu),
    id = ids,
    g = ug,
    n_rows = nr,
    n_obs = no,
    n_na_rows = n_na_u,
    n_pairs = c_res$by_unit,
    is_unobserved = c_res$is_unobserved,
    is_singleton = c_res$is_singleton,
    pair_eligible = c_res$pair_eligible,
    first_obs_time = first,
    last_obs_time = last,
    singleton_time = singleton_time,
    singleton_y = sy,
    obs_start = c_res$obs_start,
    obs_end = c_res$obs_end,
    na_start = c_res$na_start,
    na_end = c_res$na_end
  ))
  singletons <- cdid_fast_df(list(
    unit_index = single,
    id = ids[single],
    g = ug[single],
    time = first[single],
    y = sy[single],
    n_rows = nr[single]
  ))
  unobserved_units <- cdid_fast_df(list(
    unit_index = unobs,
    id = ids[unobs],
    g = ug[unobs],
    n_rows = nr[unobs],
    n_na_rows = n_na_u[unobs]
  ))
  observations <- cdid_fast_df(list(
    time = ot,
    y = oy,
    calendar_index = ocal,
    unit_index = ou,
    a_index = obs_a_index
  ))
  explicit_na_rows <- cdid_fast_df(list(
    unit_index = na_u,
    time = na_t,
    calendar_index = na_cal
  ))
  obs_times <- if (length(na_u) == 0L) all_times else if (length(ot) > 0L) sort(unique(ot), method = "radix") else numeric(0L)
  all_cohorts <- if (n_a == nu) a_cohorts else sort(unique(ug), method = "radix")

  actual_retained_entries <- ipw_owner_safe_sum(
    c(
      20 * nu,
      6 * length(single),
      5 * length(unobs),
      5 * length(ot),
      3 * length(na_u),
      5 * n_a,
      2 * length(a_cohorts),
      n_cal + length(obs_times) + length(all_cohorts),
      3 * n_cal + 3 * q_a,
      256
    ),
    max_bound = ctrl$max_retained_entries,
    label = "actual retained panel owner entries"
  )
  stopifnot(
    actual_retained_entries <= preflight$total_retained_entries,
    actual_ops <= preflight$work_ops_bound
  )

  panel_obj <- structure(
    list(
      owner_id = owner_id,
      id_type = if (id_is_num) "numeric" else "character",
      id_storage = id_storage,
      id_byte_metadata = id_byte_metadata,
      unit_id_key = ids_key,
      n_a_units = n_a,
      n_total_units = nu,
      n_singleton_units = length(single),
      n_unobserved_units = length(unobs),
      n_calendar_dates = n_cal,
      n_observed_dates = length(obs_times),
      n_observed_rows = length(ot),
      n_a_observed_rows = q_a,
      n_explicit_na_rows = length(na_u),
      a_unit_index = a_unit_index,
      unit_to_a = unit_to_a,
      a_id = a_id,
      a_id_key = a_id_key,
      a_g = a_g,
      a_cohorts = a_cohorts,
      a_cohort_index = a_cohort_index,
      cohort_a_count = cohort_a_count,
      units = units,
      singletons = singletons,
      unobserved_units = unobserved_units,
      observations = observations,
      explicit_na_rows = explicit_na_rows,
      times = obs_times,
      all_times = all_times,
      cohorts = all_cohorts,
      postings = list(
        date_start = c_res$date_start,
        date_end = c_res$date_end,
        date_count = c_res$date_count,
        post_a_index = c_res$post_a_index,
        post_unit_index = c_res$post_unit_index,
        post_obs_index = c_res$post_obs_index
      ),
      counts = list(
        n_input_rows = n,
        n_observed_rows = length(ot),
        n_na_rows = length(na_u),
        n_units = nu,
        n_unobserved_units = length(unobs),
        n_singleton_units = length(single),
        n_pair_eligible_units = n_a,
        n_pairs = c_res$total_pairs
      ),
      cols = cols,
      control = ctrl,
      preflight = preflight,
      accounting = list(
        actual_retained_vector_entries = actual_retained_entries,
        actual_retained_matrix_entries = 0,
        actual_retained_entries = actual_retained_entries,
        raw_outcome_validation_performed = TRUE,
        global_pair_subtraction_performed = FALSE,
        endpoints_read_during_ingestion = 0L
      ),
      work = list(
        validation_chunks_processed = logit_cert_count_output(validation_chunks),
        actual_work_ops = actual_ops
      )
    ),
    class = "cdid_ipw_panel_owner"
  )

  ipw_owner_seal(panel_obj, "panel_owner")
}

cdid_fast_panel_owner <- function(data, hdr, preflight, owner_id) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c("ipw_owner_check_interrupt", "ipw_panel_owner_assemble"))) return(NULL)
  if (!is.numeric(hdr$id) || is.logical(hdr$y)) return(NULL)
  n <- hdr$n
  step <- preflight$effective_chunk_rows
  last_prep <- cdid_fast_kernel_cache$last_panel_prep
  if (!is.null(last_prep) && cdid_sexp_aliases(data, last_prep$data_internal) &&
      identical(hdr$cols, last_prep$cols)) {
    ids_key <- last_prep$unit_key
    nu <- length(ids_key)
    ui <- last_prep$pos_unit
    time_ord <- last_prep$tt_d
    y_ord <- last_prep$y_ord
    g_ord <- last_prep$g_ord
  } else {
    fk_vf <- cdid_fast_kernel("C_cdid_verify_fast")
    if (is.null(fk_vf)) return(NULL)
    c_in <- fk_vf(10L, list(hdr$id, hdr$time, hdr$y, hdr$g), NULL, NULL, 1L)
    if (is.null(c_in)) return(NULL)
    ids_key <- c_in$unit_key
    nu <- length(ids_key)
    ui <- c_in$pos_unit
    time_ord <- c_in$tt_d
    y_ord <- c_in$y_ord
    g_ord <- c_in$g_ord
  }
  validation_chunks <- if (n == 0L) 0 else ceiling(as.double(n) / as.double(step))
  actual_ops <- (if (n > 1L) 10 else 8) * as.double(n)
  id_byte_metadata <- list(total = 0, max = 0)
  ipw_panel_owner_assemble(
    owner_id, hdr$ctrl, preflight, hdr$cols, n, TRUE, typeof(hdr$id), id_byte_metadata,
    ids_key, nu, ui, time_ord, y_ord, g_ord, validation_chunks, actual_ops
  )
}

cdid_fast_frame_owner <- function(panel_owner, id, Z, weights, frame_g, allow_frame_only,
                                  frame_id, design_on_A_only, hdr, preflight) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "ipw_owner_check_interrupt", "ipw_owner_match_sorted_id",
    "ipw_owner_match_sorted_double", "ipw_frame_owner_assemble", "ipw_owner_seal"
  ))) {
    return(NULL)
  }
  if (!hdr$id_is_num) return(NULL)
  nsup <- hdr$nsup
  nu <- panel_owner$n_total_units
  if (nsup != nu) return(NULL)
  if (!cdid_sexp_aliases(id, panel_owner$units$id) && !identical(as.double(id), panel_owner$unit_id_key)) {
    return(NULL)
  }
  Z_chk <- if (design_on_A_only && panel_owner$n_a_units < nu) Z[panel_owner$a_unit_index, , drop = FALSE] else Z
  if (length(Z_chk) > 0L && !is.finite(sum(Z_chk))) return(NULL)
  if (hdr$weights_supplied) {
    if (anyNA(weights) || any(!is.finite(weights) | weights <= 0)) return(NULL)
  }
  if (!is.null(frame_g)) {
    if (!cdid_sexp_aliases(frame_g, panel_owner$units$g)) {
      if (anyNA(frame_g) || any(is.nan(frame_g) | frame_g == -Inf)) return(NULL)
      fg_sub <- as.double(frame_g)
      fg_sub[fg_sub == 0.0] <- 0.0
      if (!identical(fg_sub, panel_owner$units$g)) return(NULL)
    }
  }
  step <- preflight$effective_chunk_rows
  validation_chunks <- if (nsup == 0L) 0 else ceiling(as.double(nsup) / as.double(step))
  actual_ops <- as.double(nsup) * (hdr$p + 6)
  id_byte_metadata <- list(total = 0, max = 0)
  unit_pos_of_supplied <- seq_len(nsup)
  absent_from_panel <- integer(0L)
  a_pos_in_supplied <- panel_owner$a_unit_index
  ipw_frame_owner_assemble(
    panel_owner, id, Z, weights, frame_g, hdr$ctrl, preflight, frame_id,
    hdr$n_a, nsup, hdr$p, hdr$weights_supplied, id_byte_metadata,
    unit_pos_of_supplied, absent_from_panel, a_pos_in_supplied,
    validation_chunks, actual_ops
  )
}

cdid_fast_or_baseline_pooled_fit <- function(recipe, supports, blocks, rank_tol, N, fk_comp = NULL) {
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
    fk_comp <- cdid_fast_kernel("C_cdid_wls_compress")
    if (is.null(fk_comp)) return(NULL)
  }
  p <- as.integer(recipe$p)
  n_b <- length(blocks)
  if (n_b == 1L) {
    block_names <- blocks[[1L]]$name
    if (!is.character(block_names) || length(block_names) != 1L || is.na(block_names) || !nzchar(block_names)) return(NULL)
  } else {
    block_names <- vapply(blocks, `[[`, character(1L), "name")
    if (anyNA(block_names) || any(!nzchar(block_names)) || anyDuplicated(block_names)) return(NULL)
  }
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
  block_facs <- vector("list", n_b)
  block_consts <- vector("list", n_b)
  block_ns <- integer(n_b)
  block_merges <- integer(n_b)
  block_comps <- integer(n_b)
  for (b in seq_len(n_b)) {
    bb <- blocks[[b]]
    sids <- bb$support_ids
    if (!length(sids) || any(sids < 1L | sids > n_sup) || anyDuplicated(sids)) return(NULL)
    sup1 <- supports[[sids[1L]]]
    chunk1 <- if (!is.null(sup1$design)) sup1$design$summary else sup1$summary
    sc1 <- if (!is.null(sup1$v)) sup1$v else sup1$scale
    if (!inherits(chunk1, "cdid_or_pooled_wls_chunk") || !gmm_scalar_numeric(sc1, positive = TRUE)) return(NULL)
    fac <- if (identical(sc1, 1.0)) chunk1$raw$factor else chunk1$raw$factor * sqrt(sc1)
    c_n <- as.integer(chunk1$n)
    c_const <- chunk1$constant
    c_anchor <- chunk1$anchor
    c_merge <- chunk1$accounting$max_chunk_merge_rows
    c_comp <- chunk1$accounting$n_compressions
    if (length(sids) > 1L) {
      for (sid in sids[-1L]) {
        sup_k <- supports[[sid]]
        nxt <- if (!is.null(sup_k$design)) sup_k$design$summary else sup_k$summary
        sck <- if (!is.null(sup_k$v)) sup_k$v else sup_k$scale
        if (!inherits(nxt, "cdid_or_pooled_wls_chunk") || !gmm_scalar_numeric(sck, positive = TRUE)) return(NULL)
        nxt_fac <- if (identical(sck, 1.0)) nxt$raw$factor else nxt$raw$factor * sqrt(sck)
        st_rows <- nrow(fac) + nrow(nxt_fac)
        if (st_rows > c_merge) c_merge <- st_rows
        if (nxt$n > 0L) {
          if (c_n == 0L) {
            c_anchor <- nxt$anchor
            c_const <- nxt$constant
          } else {
            c_const <- c_const & nxt$constant & (c_anchor == nxt$anchor)
          }
        }
        st_mat <- rbind(fac, nxt_fac)
        fac <- if (nrow(st_mat) > 0L) fk_comp(st_mat) else st_mat
        if (is.null(fac)) return(NULL)
        c_n <- c_n + as.integer(nxt$n)
        c_comp <- c_comp + nxt$accounting$n_compressions + if (st_rows > 0L) 1L else 0L
      }
    }
    if (c_n < 1L || nrow(fac) < 1L || ncol(fac) != p + 2L) return(NULL)
    block_facs[[b]] <- fac
    block_consts[[b]] <- as.logical(c_const)
    block_ns[b] <- c_n
    block_merges[b] <- as.integer(c_merge)
    block_comps[b] <- as.integer(c_comp)
  }

  keep <- seq_len(p)
  c_pool <- fk_comp(list(block_facs, block_consts, as.integer(p), as.double(rank_tol)))
  if (!is.null(c_pool)) {
    acc_fac <- c_pool[[1L]]
    r11_v <- c_pool[[2L]]
    r1d_m <- c_pool[[3L]]
    r1y_v <- c_pool[[4L]]
    res_rows_v <- c_pool[[5L]]
    acct_v <- c_pool[[6L]]
    coefficients <- c_pool[[7L]]
    R_out <- c_pool[[8L]]
    pivot <- c_pool[[9L]]
    stats_v <- c_pool[[10L]]
    rss <- stats_v[1L]
    cutoff <- stats_v[2L]
    intercepts <- stats_v[seq.int(3L, length.out = n_b)]
    acc_n <- sum(block_ns)
    max_aug <- acct_v[1L]
    max_res <- acct_v[2L]
    max_res_merge <- acct_v[3L]
    max_chunk_merge <- max(block_merges)
    n_compressions <- sum(block_comps) + acct_v[4L]
    if (p > 0L) {
      names(coefficients) <- cols
      rank <- p
      slope_fit <- structure(
        list(
          coefficients = coefficients, rss = rss, R = R_out, pivot = pivot,
          n = acc_n, rank = rank, rank_tol = rank_tol, cutoff = cutoff
        ),
        class = "cdid_or_wls_fit"
      )
    } else {
      if (!is.null(cols)) names(coefficients) <- cols
      rank <- 0L
      slope_fit <- NULL
    }
    regressions <- vector("list", n_b)
    for (b in seq_len(n_b)) {
      r11 <- r11_v[b]
      r1d <- if (p > 0L) r1d_m[, b] else numeric(0L)
      mean_d_un <- r1d / r11
      names(r1d) <- cols
      mean_d_nm <- mean_d_un
      names(mean_d_nm) <- cols
      r1y <- r1y_v[b]
      wsum <- r11 * r11
      mean_y <- r1y / r11
      c_n <- block_ns[b]
      m_res <- res_rows_v[b]
      block_recs[[b]] <- list(
        r11 = r11, r1d = r1d, r1y = r1y, n = c_n,
        weight_sum = wsum, mean_d = mean_d_nm, mean_y = mean_y,
        augmented_rows = nrow(block_facs[[b]]), residual_rows = m_res
      )
      regressions[[b]] <- list(
        name = block_names[b], L = blocks[[b]]$L, n = c_n, alpha = intercepts[b],
        weight_sum = wsum, mass = wsum / N,
        mean_d = mean_d_un, mean_y = mean_y,
        residual_rows = m_res
      )
    }
    names(intercepts) <- block_names
    names(block_recs) <- block_names
    beta <- unname(coefficients)
  } else {
  for (b in seq_len(n_b)) {
    fac <- block_facs[[b]]
    c_n <- block_ns[b]
    c_const <- block_consts[[b]]
    c_merge <- block_merges[b]
    c_comp <- block_comps[b]
    m_fac <- nrow(fac)
    if (any(!is.finite(fac))) return(NULL)
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
    q <- qr.default(X, LAPACK = TRUE)
    R <- q$qr[keep, , drop = FALSE]
    if (p == 2L) {
      R[2L, 1L] <- 0.0
    } else if (p > 2L) {
      R[lower.tri(R)] <- 0.0
    }
    diagonal <- abs(R[1L + (keep - 1L) * (p + 1L)])
    cutoff <- rank_tol * max(diagonal)
    rank <- if (cutoff > 0.0) sum(diagonal > cutoff) else 0L
    if (rank < p) return(NULL)
    if (min(diagonal) <= 16.0 * p * cutoff) {
      proposal <- qr.default(X, tol = rank_tol, LAPACK = FALSE)
      if (proposal$rank < p) return(NULL)
    }
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
  }
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

cdid_fast_or_baseline_outcome_fit <- function(recipe, supports, blocks, links,
                                              rank_tol, N, mode = "separate", route = "baseline") {
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
  fk_comp <- cdid_fast_kernel("C_cdid_wls_compress")
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
    g1 <- cdid_fast_or_baseline_pooled_fit(recipe, supports, blocks, rank_tol, N, fk_comp = fk_comp)
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
      gg <- cdid_fast_or_baseline_pooled_fit(recipe, supports[sids], list(local_b), rank_tol, N, fk_comp = fk_comp)
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

cdid_fast_or_public_basis <- function(Z, active, resource = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !is.matrix(Z) || !is.double(Z) || is.object(Z) ||
      !is.logical(active) || length(active) != nrow(Z) || anyNA(active) ||
      !cdid_are_default_fns(c("or_public_basis", "or_public_rank_tol", "or_pooled_outcome_recipe"))) {
    return(NULL)
  }
  all_act <- all(active)
  n <- if (all_act) nrow(Z) else sum(active)
  p <- ncol(Z)
  cols <- colnames(Z)
  if (n <= 0L || p <= 0L || is.null(cols) || anyNA(cols) || !all(nzchar(cols)) || anyDuplicated(cols)) {
    return(NULL)
  }
  if (!is.null(resource)) {
    w_bytes <- 8 * length(Z) + 64 * n * p + 64 * p * p
    w_ops <- 128 * n * p * p
    if (!is.finite(w_bytes) || w_bytes > resource$numeric_workspace_bytes ||
        !is.finite(w_ops) || w_ops > resource$counted_work_limit) {
      return(NULL)
    }
  }
  X <- if (all_act) Z else Z[active, , drop = FALSE]
  if (anyNA(X) || !all(is.finite(X))) return(NULL)
  if (!is.null(attr(Z, "intercept")) || !is.null(attr(Z, "assign"))) {
    asg <- attr(Z, "assign")
    has_intercept <- if (!is.null(attr(Z, "intercept"))) isTRUE(attr(Z, "intercept")) else (0L %in% as.integer(asg))
    icpt_idx <- if (has_intercept) {
      if (!is.null(asg) && (0L %in% as.integer(asg))) which(as.integer(asg) == 0L)[1L] else 1L
    } else integer(0)
  } else {
    has_intercept <- p >= 1L && (cols[1L] %in% c("(Intercept)", "intercept")) && all(X[, 1L] == 1)
    icpt_idx <- if (has_intercept) 1L else integer(0)
  }
  if (has_intercept) {
    if (length(icpt_idx) != 1L || is.na(icpt_idx) || icpt_idx < 1L || icpt_idx > p || !all(X[, icpt_idx] == 1)) {
      return(NULL)
    }
  }
  if (p == 1L && has_intercept) {
    retained <- 1L
  } else {
    ord <- if (has_intercept && icpt_idx != 1L) c(icpt_idx, setdiff(seq_len(p), icpt_idx)) else seq_len(p)
    q <- qr.default(if (has_intercept && icpt_idx != 1L) X[, ord, drop = FALSE] else X, tol = 1e-7, LAPACK = FALSE)
    if (q$rank < p) return(NULL)
    retained <- seq_len(p)
  }
  proposed <- integer(0L)
  aliased <- integer(0L)
  rejected <- integer(0L)
  map <- diag(p)
  dimnames(map) <- list(cols, cols)
  discrepancy <- numeric(p)
  scaled <- numeric(p)
  threshold <- 64 * .Machine$double.eps * max(1L, p)
  Z_out <- Z
  attr(Z_out, "intercept") <- NULL
  attr(Z_out, "assign") <- NULL
  slope_cols <- if (has_intercept) cols[-icpt_idx] else cols
  slope_index <- if (has_intercept) setdiff(seq_len(p), icpt_idx) else seq_len(p)
  outcome_recipe <- structure(list(
    columns = cols, intercept = has_intercept,
    intercept_index = if (has_intercept) icpt_idx else 0L,
    slope_index = slope_index, slope_columns = slope_cols,
    p = length(slope_index),
    duration = "L_e = t_calendar_index - s_calendar_index on the stored modeled calendar",
    design = "d_e = L_e * z(X) with fixed-unit z",
    key = paste0("fixed-X-summed|intercept:", has_intercept, "|z:", paste(slope_cols, collapse = ","))),
    class = "cdid_or_pooled_recipe")
  list(Z = Z_out, outcome_recipe = outcome_recipe, record = list(
    original_columns = cols, retained_columns = cols,
    aliased_columns = cols[aliased], retained_index = retained,
    intercept = has_intercept,
    intercept_column = if (has_intercept) cols[icpt_idx] else NULL,
    slope_columns = slope_cols,
    outcome_recipe = outcome_recipe,
    alias_map = map, proposed_aliases = cols[proposed],
    nontrivial_candidates_retained = cols[rejected],
    reconstruction_max_abs = stats::setNames(discrepancy, cols),
    reconstruction_max_scaled = stats::setNames(scaled, cols),
    reconstruction_roundoff_threshold = threshold,
    rank_convention = "R stats::lm.fit / dqrdc2 pivoted QR", rank_tol = 1e-7,
    control_fit_rank_check = "private compact LAPACK QR relative-diagonal criterion at the recorded lm.fit tolerance; separate from alias removal",
    preprocessing_work_ops_bound = 128 * n * p * p,
    checked_rows = n, row_universe = "all common repeat-observer rows, including controls and treated prediction rows",
    frozen_during_differentiation = TRUE, numerical_certificate = FALSE))
}

cdid_fast_chr_frame_ids <- function(frame_ids, n) {
  if (is.character(frame_ids)) return(frame_ids)
  if (is.numeric(frame_ids) && length(frame_ids) == n && n >= 1L && n <= 100000L &&
      frame_ids[1L] == 1 && frame_ids[n] == n && all(frame_ids == seq_len(n))) {
    key <- paste0("chr_seq_", n)
    hit <- get0(key, envir = cdid_fast_kernel_cache, inherits = FALSE)
    if (!is.null(hit)) return(hit)
    res <- as.character(seq_len(n))
    assign(key, res, envir = cdid_fast_kernel_cache)
    return(res)
  }
  as.character(frame_ids)
}

cdid_fast_validate_gmm_targets <- function(design, gauge, targets, bounds) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel)) ||
      !cdid_are_default_fns(c("validate_gmm_targets", "gmm_target_spec",
                              "root_path_certificate", "exact_binary64_dyadic_sum",
                              "gmm_safe_product"))) {
    return(NULL)
  }
  if (missing(targets) || is.null(targets) || !is.list(targets) || is.data.frame(targets) ||
      is.object(targets)) {
    return(NULL)
  }
  J <- length(targets)
  if (J < 1L || !is.list(bounds) || !is.finite(bounds$max_targets) || J > bounds$max_targets ||
      !is.finite(bounds$max_matrix_entries)) {
    return(NULL)
  }
  if (length(gauge$unanchored_components) != 0L) return(NULL)
  p <- as.integer(gauge$n_parameters)
  m <- as.integer(gauge$n_links)
  k <- as.integer(gauge$rank_k)
  if (is.na(p) || p < 1L || is.na(m) || m < 1L || is.na(k) || k < 1L) return(NULL)
  raw_entries <- 0
  target_ids <- character(J)
  for (j in seq_len(J)) {
    item <- targets[[j]]
    if (!is.list(item) || is.object(item) ||
        !all(c("target_id", "param_index", "coeff") %in% names(item))) {
      return(NULL)
    }
    tid <- item$target_id
    pi_j <- item$param_index
    co_j <- item$coeff
    if (!is.character(tid) || length(tid) != 1L || is.na(tid) || !nzchar(tid)) return(NULL)
    n_j <- length(pi_j)
    if (!is.numeric(pi_j) || n_j < 1L || anyNA(pi_j) || any(!is.finite(pi_j)) ||
        any(pi_j != floor(pi_j)) || any(pi_j < 1 | pi_j > p) || anyDuplicated(pi_j)) {
      return(NULL)
    }
    if (!is.numeric(co_j) || length(co_j) != n_j || anyNA(co_j) || any(!is.finite(co_j)) ||
        any(co_j == 0) || any(abs(co_j) > 1e9) || any(co_j != floor(co_j))) {
      return(NULL)
    }
    target_ids[j] <- tid
    raw_entries <- raw_entries + n_j
  }
  if (anyDuplicated(target_ids)) return(NULL)
  max_e <- as.double(bounds$max_matrix_entries)
  if (4 * raw_entries > max_e ||
      160 * min(max(raw_entries, p), 1024) > max_e ||
      as.double(p) * J > max_e ||
      as.double(m) * J > max_e ||
      as.double(p) * p > max_e ||
      as.double(k) * J > max_e) {
    return(NULL)
  }
  f <- design$forest
  p_vert <- f$parent_vertex
  p_link <- f$parent_link
  p_sign <- f$parent_sign
  v_verts <- design$parameters$vertex_index
  inc_start <- design$incidence$start
  inc_end <- design$incidence$end
  nv <- length(p_vert)
  if (!is.integer(p_vert) || !is.integer(p_link) || !is.integer(p_sign) ||
      !is.integer(v_verts) || length(v_verts) != p ||
      !is.integer(inc_start) || length(inc_start) != m ||
      !is.integer(inc_end) || length(inc_end) != m) {
    return(NULL)
  }
  c_vectors <- vector("list", J)
  witnesses <- vector("list", J)
  keep_p <- gauge$keep_params
  C_gauge <- matrix(0, nrow = k, ncol = J,
                    dimnames = list(gauge$keep_param_ids, target_ids))
  for (j in seq_len(J)) {
    item <- targets[[j]]
    pi_j <- as.integer(item$param_index)
    co_j <- as.double(item$coeff)
    c_vec <- numeric(p)
    c_vec[pi_j] <- co_j
    c_vectors[[j]] <- c_vec
    C_gauge[, j] <- c_vec[keep_p]
    if (length(pi_j) == 1L) {
      cur <- v_verts[pi_j]
      if (is.na(cur) || cur < 1L || cur > nv) return(NULL)
      k_len <- 0L
      tmp <- cur
      while (p_vert[tmp] != 0L) {
        k_len <- k_len + 1L
        if (k_len > nv) return(NULL)
        tmp <- p_vert[tmp]
        if (is.na(tmp) || tmp < 0L || tmp > nv) return(NULL)
      }
      if (k_len == 0L) return(NULL)
      edges <- integer(k_len)
      weights <- numeric(k_len)
      cur <- v_verts[pi_j]
      for (s_i in seq_len(k_len)) {
        edges[s_i] <- p_link[cur]
        weights[s_i] <- co_j * p_sign[cur]
        cur <- p_vert[cur]
      }
      if (anyNA(edges) || any(edges < 1L | edges > m) || anyDuplicated(edges)) return(NULL)
      acc <- numeric(p)
      acc[pi_j] <- -co_j
      for (s_i in seq_len(k_len)) {
        li <- edges[s_i]
        we <- weights[s_i]
        ie <- inc_end[li]
        is_ <- inc_start[li]
        if (ie > 0L) acc[ie] <- acc[ie] + we
        if (is_ > 0L) acc[is_] <- acc[is_] - we
      }
      if (any(acc != 0)) return(NULL)
      if (k_len > 1L) {
        ord <- order(edges)
        edges <- edges[ord]
        weights <- weights[ord]
      }
      witnesses[[j]] <- cdid_fast_df(list(link_index = edges, weight = weights))
    } else {
      w_full <- numeric(m)
      for (idx in seq_along(pi_j)) {
        cur <- v_verts[pi_j[idx]]
        if (is.na(cur) || cur < 1L || cur > nv) return(NULL)
        cv <- co_j[idx]
        steps <- 0L
        while (p_vert[cur] != 0L) {
          steps <- steps + 1L
          if (steps > nv) return(NULL)
          li <- p_link[cur]
          if (is.na(li) || li < 1L || li > m) return(NULL)
          w_full[li] <- w_full[li] + cv * p_sign[cur]
          cur <- p_vert[cur]
          if (is.na(cur) || cur < 0L || cur > nv) return(NULL)
        }
      }
      nz_l <- which(w_full != 0)
      acc <- -c_vec
      for (li in nz_l) {
        wl <- w_full[li]
        ie <- inc_end[li]
        is_ <- inc_start[li]
        if (ie > 0L) acc[ie] <- acc[ie] + wl
        if (is_ > 0L) acc[is_] <- acc[is_] - wl
      }
      if (any(acc != 0)) return(NULL)
      witnesses[[j]] <- cdid_fast_df(list(link_index = as.integer(nz_l), weight = as.double(w_full[nz_l])))
    }
  }
  list(
    n_targets = J,
    n_estimable = J,
    target_ids = target_ids,
    estimable = rep.int(TRUE, J),
    status = rep.int("estimable", J),
    support_reason = rep.int("exact_zero_sum_on_all_unanchored_components", J),
    c_vectors = c_vectors,
    C_gauge = C_gauge,
    estimable_indices = seq_len(J),
    witnesses = witnesses
  )
}

cdid_fast_extract_catalog_row_specs <- function(metadata_obj, retained_entries = 0,
                                                external_work_ops = 0, control = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "extract_pipeline_catalog_row_specs", "validate_owner_comparison_catalog",
    "pipeline_metadata_row_spec_preflight", "pipeline_metadata_descriptor_source",
    "pipeline_metadata_row_spec_fields"
  ))) {
    return(NULL)
  }
  if (!inherits(metadata_obj, "cdid_pipeline_metadata") || !is.list(metadata_obj) ||
      !ipw_owner_seal_ok(metadata_obj, "pipeline_metadata")) {
    return(NULL)
  }
  rows <- metadata_obj$catalog$rows
  map <- metadata_obj$row_mapping
  pl <- metadata_obj$catalog$pools
  coh <- metadata_obj$catalog$spec$cohorts
  cal <- metadata_obj$catalog$spec$calendar
  n_rows <- nrow(rows)
  if (n_rows < 1L) return(NULL)
  own_ctrl <- metadata_obj$control$owner_control
  eff_ctrl <- if (is.null(control)) own_ctrl else control
  max_k <- max(pl$size[rows$pool_id])
  if (!is.finite(own_ctrl$max_matrix_entries) || !is.finite(own_ctrl$max_output_entries) ||
      !is.finite(own_ctrl$max_work_ops) ||
      128 + 4 * max_k > min(own_ctrl$max_matrix_entries, own_ctrl$max_output_entries) ||
      128 + 16 * max_k > own_ctrl$max_work_ops) {
    return(NULL)
  }
  accounting <- pipeline_metadata_descriptor_preflight(
    metadata_obj, eff_ctrl, retained_entries, external_work_ops, n_rows, include_output = TRUE
  )
  P <- accounting$n_pools
  K_indices_list <- vector("list", P)
  K_vals_list <- vector("list", P)
  pl_entries <- pl$entries
  pl_start <- pl$start
  pl_size <- pl$size
  for (pid in seq_len(P)) {
    ki <- as.integer(pl_entries[seq.int(pl_start[pid], length.out = pl_size[pid])])
    K_indices_list[[pid]] <- ki
    K_vals_list[[pid]] <- as.double(coh[ki])
  }

  r_op <- rows$operator_id
  r_pid <- rows$pool_id
  r_ti <- rows$target_index
  r_si <- rows$s_index
  r_tti <- rows$t_index
  r_refi <- rows$ref_index
  r_srole <- rows$s_role
  r_trole <- rows$t_role
  r_scid <- rows$s_coordinate_id
  r_tcid <- rows$t_coordinate_id
  r_orig <- rows$original
  r_nh <- rows$n_horizons
  r_hfi <- rows$horizon_first_index
  r_hli <- rows$horizon_last_index

  m_fk <- map$functional_key
  m_tocid <- map$target_occupancy_cell_id
  m_tpc <- map$target_pair_count
  m_cppc <- map$control_pool_pair_count
  m_rpc <- map$risk_pair_count
  m_tauc <- map$target_a_unit_count
  m_cpauc <- map$control_pool_a_unit_count
  m_rauc <- map$risk_a_unit_count

  r_g <- coh[r_ti]
  r_s <- cal[r_si]
  r_t <- cal[r_tti]
  r_ref <- cal[r_refi]

  out <- vector("list", n_rows)
  for (ri in seq_len(n_rows)) {
    pid <- r_pid[ri]
    out[[ri]] <- list(
      row_id = ri,
      operator_id = r_op[ri],
      pool_id = pid,
      target_index = r_ti[ri],
      g = r_g[ri],
      K_indices = K_indices_list[[pid]],
      K = K_vals_list[[pid]],
      s_index = r_si[ri],
      t_index = r_tti[ri],
      s = r_s[ri],
      t = r_t[ri],
      ref_index = r_refi[ri],
      ref = r_ref[ri],
      functional_key = m_fk[ri],
      s_role = r_srole[ri],
      t_role = r_trole[ri],
      s_coordinate_id = r_scid[ri],
      t_coordinate_id = r_tcid[ri],
      original = r_orig[ri],
      n_horizons = r_nh[ri],
      horizon_first_index = r_hfi[ri],
      horizon_last_index = r_hli[ri],
      target_occupancy_cell_id = m_tocid[ri],
      target_pair_count = m_tpc[ri],
      control_pool_pair_count = m_cppc[ri],
      risk_pair_count = m_rpc[ri],
      target_a_unit_count = m_tauc[ri],
      control_pool_a_unit_count = m_cpauc[ri],
      risk_a_unit_count = m_rauc[ri]
    )
  }
  list(specs = out, accounting = accounting)
}

cdid_fast_pipeline_metadata_pair <- function(panel_owner, frame_owner, targets, anticipation,
                                             control_group, conventions, requested_post,
                                             control, catalog_bounds, declared_qualified,
                                             ledger, pt_assumption, pt_start) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "build_owner_comparison_catalog", "preflight_owner_comparison_catalog",
    "derive_owner_sparse_occupancy", "validate_owner_sparse_occupancy",
    "pipeline_metadata_bind_catalog_source", "catalog_adapter_preflight",
    "catalog_design_handoff", "build_catalog_handoff_design", "ipw_owner_seal",
    "preflight_owner_sparse_occupancy", "pipeline_occupancy_accumulate",
    "pipeline_metadata_catalog_reservation", "pipeline_metadata_resolve_catalog_mappings",
    "build_comparison_catalog"
  ))) {
    return(NULL)
  }
  fk_seal <- cdid_fast_kernel("C_cdid_fast_seal")
  if (is.null(fk_seal)) return(NULL)

  meta_initial <- NULL
  if (is.null(requested_post) && identical(pt_assumption, "post") && is.null(pt_start) &&
      inherits(panel_owner, "cdid_ipw_panel_owner") && ipw_owner_seal_ok(panel_owner, "panel_owner") &&
      (is.null(frame_owner) || (inherits(frame_owner, "cdid_ipw_frame_owner") &&
       ipw_owner_seal_ok(frame_owner, "frame_owner") && identical(frame_owner$panel_ref, panel_owner$ref))) &&
      is.numeric(targets) && is.null(attributes(targets)) && length(targets) >= 1L &&
      all(is.finite(targets)) && !anyDuplicated(as.double(targets)) &&
      is.numeric(anticipation) && is.null(attributes(anticipation)) && length(anticipation) == 1L &&
      is.finite(anticipation) && anticipation >= 0 && anticipation == floor(anticipation) &&
      is.character(control_group) && is.null(attributes(control_group)) && length(control_group) == 1L &&
      control_group %in% c("nevertreated", "notyettreated") &&
      is.character(conventions) && is.null(attributes(conventions)) && length(conventions) == 1L &&
      !is.na(conventions) && nzchar(conventions)) {
    pipeline_metadata_require_stage(c("owners", "occupancy", "catalog", "design"))
    pm_ctrl <- tryCatch(validate_pipeline_metadata_control(control, catalog_bounds), error = function(e) NULL)
    if (!is.null(pm_ctrl)) {
      own_ctrl <- pm_ctrl$owner_control
      cat_bnd <- pm_ctrl$catalog_bounds
      n_a <- panel_owner$n_a_units
      n_u <- panel_owner$n_total_units
      n_cal <- panel_owner$n_calendar_dates
      n_coh <- length(panel_owner$cohorts)
      n_a_coh <- length(panel_owner$a_cohorts)
      n_pairs <- panel_owner$counts$n_pairs
      q_a <- panel_owner$n_a_observed_rows
      max_ret <- own_ctrl$max_retained_entries
      max_me <- own_ctrl$max_matrix_entries
      max_wo <- own_ctrl$max_work_ops
      max_out <- own_ctrl$max_output_entries
      tg <- sort(match(as.double(targets), panel_owner$cohorts))
      if (!anyNA(tg) && n_cal < 1e6L && n_coh < 1e6L && n_a_coh < 1e6L && n_pairs < 1e9 &&
          is.finite(max_ret) && is.finite(max_me) && is.finite(max_wo) && is.finite(max_out)) {
        max_date_pairs <- if (n_cal >= 2L) as.double(n_cal) * (n_cal - 1) / 2 else 0
        max_grid_cells <- if (n_a_coh > 0L && max_date_pairs > 0) as.double(n_a_coh) * max_date_pairs else 0
        max_occ_cells <- min(as.double(n_pairs), max_grid_cells)
        retained_vec_e <- 9 * max_occ_cells + 128
        header_e <- as.double(n_a_coh + n_coh + n_cal)
        grouping_scr_e <- 20 * as.double(n_a) + 2 * as.double(n_cal) + 256
        accum_fixed_e <- 4 * as.double(n_a) + 11 * as.double(q_a) + 4 * max_occ_cells + 4 * as.double(n_cal) + 256
        t_span <- max(1, n_cal - 1)
        slice_room <- max_me - (retained_vec_e + header_e + accum_fixed_e)
        slice_pairs <- max(t_span, min(65536, max(1, as.double(n_pairs)), if (slice_room > 0) floor(slice_room / 14) else 0))
        accum_scr_e <- accum_fixed_e + 14 * slice_pairs
        assem_scr_e <- 9 * max_occ_cells + 256
        conc_e <- retained_vec_e + header_e + max(grouping_scr_e, accum_scr_e, assem_scr_e)
        sf_na <- 1 + ceiling(log2(n_a + 1))
        sf_qa <- 1 + ceiling(log2(q_a + 1))
        sf_sp <- 1 + ceiling(log2(slice_pairs + 1))
        grp_wo <- 64 * sf_na * n_a + 32 * sf_na * q_a + 256 * n_cal + 256
        row_bnd <- min(as.double(q_a), as.double(n_a_coh) * n_cal)
        acc_wo <- 16 * q_a + 32 * sf_qa * q_a + 32 * sf_sp * n_pairs + 8 * n_pairs + 256 * row_bnd + 16 * max_occ_cells + 256
        asm_wo <- 32 * max_occ_cells + 8 * header_e + 256
        tot_wo <- grp_wo + acc_wo + asm_wo
        nJ0 <- length(targets)
        b_conv <- as.double(nchar(conventions, type = "bytes"))
        H_min <- as.double(n_cal + n_coh + nJ0 + 1)
        H_max <- H_min + max_occ_cells
        min_hdr_scr <- 64 * H_min + 16 * b_conv
        max_hdr_scr <- 64 * H_max + 16 * b_conv
        min_hdr_wrk <- 64 * (1 + ceiling(log2(H_min + 1))) * H_min + 16 * b_conv
        max_hdr_wrk <- 64 * (1 + ceiling(log2(H_max + 1))) * H_max + 16 * b_conv
        if (retained_vec_e <= min(max_ret, max_out) && conc_e <= max_me && tot_wo <= max_wo &&
            min_hdr_wrk <= cat_bnd$max_work && min_hdr_scr <= cat_bnd$max_scratch &&
            retained_vec_e + min_hdr_scr <= max_me && tot_wo + min_hdr_wrk <= max_wo &&
            max(abs(panel_owner$all_times), abs(panel_owner$cohorts[ is.finite(panel_owner$cohorts) ]), anticipation) <= 1e14) {
          occ_pre <- list(
            n_a_units = n_a, n_total_units = n_u, n_calendar_dates = n_cal,
            n_cohorts = n_coh, n_a_cohorts = n_a_coh, n_pairs = n_pairs,
            n_a_observed_rows = q_a, max_date_pairs = max_date_pairs,
            max_occupancy_cells = max_occ_cells,
            accumulation_method = "exact_repeated_history_rowwise_dense_buffer",
            slice_pairs = logit_cert_count_output(slice_pairs),
            max_refinement_positions = logit_cert_count_output(n_cal),
            matrix_dim_limit = as.double(.Machine$integer.max),
            vector_length_limit = 2^52, exact_count_limit = 2^53 - 1,
            retained_matrix_entries = 0, retained_vector_entries = retained_vec_e,
            total_retained_entries = retained_vec_e, output_vector_entries = retained_vec_e,
            header_entries = header_e, grouping_scratch_entries = grouping_scr_e,
            accumulation_scratch_entries = accum_scr_e, assembly_scratch_entries = assem_scr_e,
            concurrent_matrix_entries = 0, concurrent_vector_entries = conc_e,
            concurrent_entries = conc_e, concurrent_bytes_8double = 8 * conc_e,
            grouping_work_ops = grp_wo, accumulation_work_ops = acc_wo,
            assembly_work_ops = asm_wo, work_ops_bound = tot_wo,
            plan = list(
              phase_1 = "exact grouping of A-units by (cohort, observed calendar-index sequence); released before phase 2",
              phase_2 = "distinct-history pairs in slices of at most slice_pairs into one dense length-T buffer per (cohort, s) row, emitted in ascending t",
              distinct_histories_upper = n_a, postings_upper = q_a, pair_visits_upper = n_pairs,
              fill_buffer_cells_upper = max_occ_cells, row_buffer_length = n_cal,
              owner_pair_reader_used = FALSE, unit_id_values_read = FALSE, chunk_rows_control_used = FALSE
            ),
            all_unit_pair_table_materialized = FALSE, endpoints_read = 0L,
            global_pair_subtraction_performed = FALSE, peak_rss_promised = FALSE,
            accounting_scope = "Conservative atomic vector/dictionary slot bounds for the exact grouping/accumulation algorithm before any allocation; excludes R SEXP/list headers and OS process peak RSS."
          )
          a_to_all <- if (n_a_coh > 0L) as.integer(ipw_owner_match_sorted_double(panel_owner$a_cohorts, panel_owner$cohorts)) else integer(0L)
          acc_occ <- tryCatch(pipeline_occupancy_accumulate(panel_owner, occ_pre), error = function(e) NULL)
          if (!is.null(acc_occ) && !anyNA(a_to_all)) {
            acc_diag <- acc_occ$diagnostics
            dict_ki <- a_to_all[acc_occ$cohort_a_index]
            dict_si <- acc_occ$s_index
            dict_ti <- acc_occ$t_index
            dict_count <- acc_occ$count
            n_cells <- length(dict_ki)
            tot_occ_cnt <- sum(dict_count)
            act_ret <- 9 * n_cells + 32
            act_conc <- act_ret + header_e + acc_diag$peak_phase_entries
            act_ops <- 2 * (n_a_coh + n_coh + n_cal) + acc_diag$actual_work_ops
            if (tot_occ_cnt == n_pairs && acc_diag$n_pairs_represented == n_pairs &&
                n_cells <= max_occ_cells && act_ret <= retained_vec_e &&
                act_conc <= conc_e && act_ops <= tot_wo) {
              occ_df <- cdid_fast_df(list(
                cohort = as.double(panel_owner$cohorts[dict_ki]),
                s = as.double(panel_owner$all_times[dict_si]),
                t = as.double(panel_owner$all_times[dict_ti]),
                count = as.double(dict_count)
              ))
              ci_df <- cdid_fast_df(list(
                cell_id = seq_len(n_cells),
                cohort_index = as.integer(dict_ki),
                s_index = as.integer(dict_si),
                t_index = as.integer(dict_ti),
                count = as.double(dict_count)
              ))
              owner_sum <- pipeline_metadata_owner_summary(panel_owner, frame_owner)
              occ_obj_raw <- structure(
                list(
                  occupancy_id = NULL,
                  panel_ref = panel_owner$ref,
                  frame_ref = if (!is.null(frame_owner)) frame_owner$ref else NULL,
                  occupancy = occ_df,
                  cell_indices = ci_df,
                  calendar = panel_owner$all_times,
                  observed_times = panel_owner$times,
                  cohorts = panel_owner$cohorts,
                  a_cohorts = panel_owner$a_cohorts,
                  owner_summary = owner_sum,
                  control = own_ctrl,
                  preflight = occ_pre,
                  diagnostics = list(
                    endpoints_read = 0L, read_outcomes = FALSE,
                    global_pair_subtraction_performed = FALSE, pair_differences_computed = FALSE,
                    all_unit_pair_table_materialized = FALSE, n_by_moment_cache_allocated = FALSE,
                    formatted_float_keys_used = FALSE, dictionary_index_keys_used = TRUE,
                    singleton_metadata_preserved_in_owner = TRUE,
                    unobserved_metadata_preserved_in_owner = TRUE,
                    explicit_na_metadata_preserved_in_owner = TRUE,
                    accumulation_method = occ_pre$accumulation_method,
                    owner_pair_reader_used = FALSE, unit_id_values_read = FALSE,
                    n_distinct_histories = logit_cert_count_output(acc_diag$n_distinct_histories),
                    max_history_multiplicity = logit_cert_count_output(acc_diag$max_history_multiplicity),
                    n_refinement_positions = logit_cert_count_output(acc_diag$n_refinement_positions),
                    refinement_unit_visits = acc_diag$refinement_unit_visits,
                    n_postings = acc_diag$n_postings,
                    n_rows = logit_cert_count_output(acc_diag$n_rows),
                    n_slices = logit_cert_count_output(acc_diag$n_slices),
                    max_slice_pairs = logit_cert_count_output(acc_diag$max_slice_pairs),
                    slice_pairs = occ_pre$slice_pairs,
                    n_pairs_represented = acc_diag$n_pairs_represented,
                    n_history_pair_visits = acc_diag$n_history_pair_visits,
                    pair_visits_avoided = n_pairs - acc_diag$n_history_pair_visits,
                    n_occupancy_cells = logit_cert_count_output(n_cells),
                    total_occupancy_count = logit_cert_count_output(tot_occ_cnt),
                    actual_retained_entries = act_ret,
                    actual_concurrent_entries = act_conc,
                    actual_work_ops = act_ops,
                    temporary_dictionaries = list(
                      cohort_a_to_all_map_length = n_a_coh,
                      row_buffer_length = n_cal,
                      fill_buffer_cells = acc_diag$fill_buffer_cells
                    )
                  )
                ),
                class = "cdid_pipeline_occupancy"
              )
              occ_obj_fast <- fk_seal(occ_obj_raw, "pipeline_occupancy")
              cat_fast <- tryCatch(
                build_comparison_catalog(
                  calendar = panel_owner$all_times, cohorts = panel_owner$cohorts,
                  targets = targets, anticipation = anticipation, occupancy = occ_df,
                  control_group = control_group, conventions = conventions,
                  requested_post = NULL, bounds = cat_bnd
                ),
                error = function(e) NULL
              )
              if (!is.null(occ_obj_fast) && !is.null(cat_fast)) {
                H_res <- H_min + n_cells
                hdr_scr_res <- 64 * H_res + 16 * b_conv
                hdr_wrk_res <- 64 * (1 + ceiling(log2(H_res + 1))) * H_res + 16 * b_conv
                ncand_res <- as.double(cat_fast$resources$candidates)
                ntp_res <- as.double(cat_fast$resources$target_pairs)
                R_res <- min(ncand_res, cat_bnd$max_rows)
                E_res <- min(R_res * n_coh, cat_bnd$max_pool_entries)
                stage_res <- ncand_res + ntp_res + 1
                key_slots_res <- 16 * b_conv * ncand_res
                scr_res <- hdr_scr_res + 9 * n_cells + n_cal + nJ0 + 8 * ntp_res + 256 * stage_res +
                  key_slots_res + 5 * ncand_res + 64 * E_res + 13 * R_res + 4 * nJ0
                cat_wrk_res <- hdr_wrk_res + 128 * (n_coh + 1 + ceiling(log2(stage_res + 1))) * stage_res + key_slots_res
                map_up_res <- retained_vec_e + scr_res + 43 * R_res + 256
                wrk_pre_res <- tot_wo + max_hdr_wrk + 2 * hdr_wrk_res + cat_wrk_res
                n_rows_cat <- nrow(cat_fast$rows)
                n_ops_cat <- nrow(cat_fast$operators)
                n_pools_cat <- cat_fast$pools$n_pools
                n_pe_cat <- length(cat_fast$pools$entries)
                map_ret <- retained_vec_e + cat_fast$resources$scratch_charged + 30 * n_rows_cat + 13 * n_ops_cat + 256
                map_scr <- map_ret + 16 * n_cells + 12 * n_coh + 8 * (n_pools_cat + n_pe_cat + 1) + 256
                rp_vis <- as.double(n_rows_cat) * n_coh
                map_wrk <- wrk_pre_res + 32 * (1 + ceiling(log2(n_cells + 1)) + ceiling(log2(n_coh + 1))) *
                  (n_cells + rp_vis + n_pe_cat + n_ops_cat + n_coh + 1) + 256
                if (map_up_res <= min(max_ret, max_out) && retained_vec_e + scr_res <= max_me &&
                    wrk_pre_res <= max_wo && map_ret <= min(max_ret, max_out) &&
                    map_scr <= max_me && map_wrk <= max_wo && n_rows_cat >= 1L) {
                  coh_a_cnt <- integer(n_coh)
                  coh_a_wt <- if (!is.null(frame_owner)) numeric(n_coh) else rep(NA_real_, n_coh)
                  if (n_a_coh > 0L) {
                    coh_a_cnt[a_to_all] <- panel_owner$cohort_a_count
                    if (!is.null(frame_owner)) coh_a_wt[a_to_all] <- frame_owner$cohort_a_weight
                  }
                  pl <- cat_fast$pools
                  pool_a_cnt <- numeric(n_pools_cat)
                  pool_a_wt <- if (!is.null(frame_owner)) numeric(n_pools_cat) else rep(NA_real_, n_pools_cat)
                  for (p_i in seq_len(n_pools_cat)) {
                    p_idx <- pl$entries[seq.int(pl$start[p_i], length.out = pl$size[p_i])]
                    pool_a_cnt[p_i] <- sum(as.double(coh_a_cnt[p_idx]))
                    if (!is.null(frame_owner)) pool_a_wt[p_i] <- sum(coh_a_wt[p_idx])
                  }
                  ord_pair <- order(dict_si, dict_ti, dict_ki, method = "radix")
                  p_cell_id <- seq_len(n_cells)[ord_pair]
                  p_ki <- dict_ki[ord_pair]
                  p_si <- dict_si[ord_pair]
                  p_ti <- dict_ti[ord_pair]
                  p_cnt <- dict_count[ord_pair]
                  new_p <- if (n_cells > 1L) c(TRUE, p_si[-1L] != p_si[-n_cells] | p_ti[-1L] != p_ti[-n_cells]) else TRUE
                  pair_starts <- which(new_p)
                  pair_ends <- c(pair_starts[-1L] - 1L, n_cells)
                  pair_s_key <- p_si[pair_starts]
                  pair_t_key <- p_ti[pair_starts]
                  rows_cf <- cat_fast$rows
                  fk_map <- cdid_fast_kernel("C_cdid_catalog_row_mappings")
                  fast_map <- if (!is.null(fk_map)) {
                    fk_map(
                      list(as.integer(pair_s_key), as.integer(pair_t_key), as.integer(pair_starts), as.integer(pair_ends), as.integer(p_ki), as.double(p_cnt), as.integer(p_cell_id)),
                      list(as.integer(pl$entries), as.integer(pl$start), as.integer(pl$size)),
                      list(as.integer(rows_cf$target_index), as.integer(rows_cf$s_index), as.integer(rows_cf$t_index), as.integer(rows_cf$pool_id))
                    )
                  } else NULL
                  if (!is.null(fast_map)) {
                    t_a_row <- as.double(coh_a_cnt[rows_cf$target_index])
                    c_a_row <- pool_a_cnt[rows_cf$pool_id]
                    row_map_df <- cdid_fast_df(list(
                      row_id = rows_cf$row_id, operator_id = rows_cf$operator_id,
                      pool_id = rows_cf$pool_id, pool_size = rows_cf$pool_size,
                      target_index = rows_cf$target_index, g = rows_cf$g,
                      s_index = rows_cf$s_index, t_index = rows_cf$t_index,
                      s = rows_cf$s, t = rows_cf$t, ref_index = rows_cf$ref_index,
                      ref = panel_owner$all_times[rows_cf$ref_index],
                      functional_key = cat_fast$operators$functional_key[rows_cf$operator_id],
                      s_role = rows_cf$s_role, t_role = rows_cf$t_role,
                      s_coordinate_id = rows_cf$s_coordinate_id, t_coordinate_id = rows_cf$t_coordinate_id,
                      original = rows_cf$original, n_horizons = rows_cf$n_horizons,
                      horizon_first_index = rows_cf$horizon_first_index, horizon_last_index = rows_cf$horizon_last_index,
                      target_occupancy_cell_id = fast_map$target_occ_cell_id,
                      target_pair_count = fast_map$target_pair_count,
                      control_pool_pair_count = fast_map$control_pool_pair_count,
                      risk_pair_count = fast_map$risk_pair_count,
                      all_cohort_pair_count = fast_map$all_cohort_pair_count,
                      target_a_unit_count = t_a_row, control_pool_a_unit_count = c_a_row,
                      risk_a_unit_count = t_a_row + c_a_row,
                      n_supported_controls = fast_map$n_supported_controls
                    ))
                    ops_cf <- cat_fast$operators
                    op_rc <- if (n_ops_cat > 0L) tabulate(rows_cf$operator_id, nbins = n_ops_cat) else integer(0L)
                    op_ta <- as.double(coh_a_cnt[ops_cf$target_index])
                    op_ca <- if (n_ops_cat > 0L) pool_a_cnt[ops_cf$pool_id] else numeric(0L)
                    op_map_df <- cdid_fast_df(list(
                      operator_id = ops_cf$operator_id, target_index = ops_cf$target_index,
                      g = ops_cf$g, pool_id = ops_cf$pool_id,
                      pool_size = if (n_ops_cat > 0L) pl$size[ops_cf$pool_id] else integer(0L),
                      never_treated_in_pool = if (n_ops_cat > 0L) pl$never_treated[ops_cf$pool_id] else logical(0L),
                      functional_key = ops_cf$functional_key, n_catalog_rows = op_rc,
                      target_a_unit_count = op_ta, control_pool_a_unit_count = op_ca,
                      risk_a_unit_count = op_ta + op_ca,
                      target_a_weight = if (!is.null(frame_owner)) coh_a_wt[ops_cf$target_index] else rep(NA_real_, n_ops_cat),
                      control_pool_a_weight = if (!is.null(frame_owner) && n_ops_cat > 0L) pool_a_wt[ops_cf$pool_id] else rep(NA_real_, n_ops_cat)
                    ))
                    reservation_fast <- list(
                      catalog_scratch_upper = scr_res, catalog_work_upper = cat_wrk_res,
                      catalog_candidate_count = ncand_res, catalog_target_pair_count = ntp_res,
                      mapping_retained_upper = map_up_res, work_prefix = wrk_pre_res
                    )
                    map_acc_fast <- list(retained_entries = map_ret, concurrent_entries = map_scr, work_ops_bound = map_wrk)
                    chp1_fast <- list(
                      n_calendar_dates = n_cal, n_cohorts = n_coh, n_targets = nJ0,
                      n_requested_post = 0L, convention_bytes = b_conv,
                      n_occupancy_lower = 0, n_occupancy_upper = max_occ_cells,
                      min_header_scratch = min_hdr_scr, max_header_scratch = max_hdr_scr,
                      min_header_work = min_hdr_wrk, max_header_work = max_hdr_wrk,
                      catalog_bounds = cat_bnd
                    )
                    pre1_fast <- list(
                      occupancy_supplied = FALSE, occupancy_preflight = occ_pre,
                      catalog_header_preflight = chp1_fast, owner_control = own_ctrl,
                      catalog_bounds = cat_bnd, build_design = FALSE, endpoints_read = 0L,
                      peak_rss_promised = FALSE, catalog_composition_reservation = reservation_fast,
                      catalog_reuse_binding = NULL
                    )
                    diag1_fast <- list(
                      catalog_reuse = list(
                        reused = FALSE, catalog_constructed = TRUE, mappings_constructed = TRUE,
                        borrowed_catalog_mapping_entries = 0, owned_catalog_mapping_entries = map_ret,
                        owned_design_entries = 0, owned_work_ops = map_wrk,
                        borrowed_work_ops = 0, source_metadata_ref = NULL
                      ),
                      endpoints_read = 0L, read_outcomes = FALSE,
                      global_pair_subtraction_performed = FALSE, pair_differences_computed = FALSE,
                      all_unit_pair_table_materialized = FALSE, n_by_moment_cache_allocated = FALSE,
                      formatted_float_keys_used = FALSE, dictionary_index_keys_used = TRUE,
                      all_catalog_rows_retained = TRUE, all_catalog_operators_retained = TRUE,
                      singleton_metadata_preserved_in_owner = TRUE,
                      unobserved_metadata_preserved_in_owner = TRUE,
                      explicit_na_metadata_preserved_in_owner = TRUE,
                      occupancy_diagnostics = occ_obj_fast$diagnostics,
                      catalog_resources = cat_fast$resources,
                      mapping_accounting = map_acc_fast
                    )
                    meta_init_raw <- structure(
                      list(
                        metadata_id = NULL, panel_ref = panel_owner$ref,
                        frame_ref = if (!is.null(frame_owner)) frame_owner$ref else NULL,
                        occupancy_ref = occ_obj_fast$ref, occupancy = occ_df,
                        cell_indices = ci_df, occupancy_object = occ_obj_fast,
                        catalog = cat_fast, row_mapping = row_map_df,
                        operator_mapping = op_map_df, design_links = NULL, design = NULL,
                        owner_summary = owner_sum, control = pm_ctrl,
                        construction = list(
                          targets = targets, anticipation = anticipation,
                          control_group = control_group, conventions = conventions,
                          requested_post = requested_post, pt_assumption = pt_assumption, pt_start = pt_start
                        ),
                        catalog_source_ref = NULL, preflight = pre1_fast, diagnostics = diag1_fast
                      ),
                      class = "cdid_pipeline_metadata"
                    )
                    meta_initial <- fk_seal(meta_init_raw, "pipeline_metadata")
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (is.null(meta_initial)) {
    meta_initial <- build_owner_comparison_catalog(
      panel_owner = panel_owner,
      targets = targets,
      anticipation = anticipation,
      control_group = control_group,
      conventions = conventions,
      requested_post = requested_post,
      control = control,
      catalog_bounds = catalog_bounds,
      frame_owner = frame_owner,
      build_design = FALSE,
      pt_assumption = pt_assumption,
      pt_start = pt_start
    )
  }
  c_rows <- meta_initial$catalog$rows
  if (nrow(c_rows) < 1L) return(NULL)

  if (is.null(declared_qualified)) {
    s_causal <- c_rows$s_role %in% c("anticipation_conditional", "post_conditional")
    t_causal <- c_rows$t_role %in% c("anticipation_conditional", "post_conditional")
    dq_g <- c(c_rows$g[s_causal], c_rows$g[t_causal])
    dq_t <- c(c_rows$s[s_causal], c_rows$t[t_causal])
    if (length(dq_g) > 0L) {
      ord_dq <- order(dq_g, dq_t)
      dq_g <- as.double(dq_g[ord_dq])
      dq_t <- as.double(dq_t[ord_dq])
      n_dq <- length(dq_g)
      keep_dq <- if (n_dq > 1L) c(TRUE, dq_g[-1L] != dq_g[-n_dq] | dq_t[-1L] != dq_t[-n_dq]) else TRUE
      dq_df <- cdid_fast_df(list(g = dq_g[keep_dq], t = dq_t[keep_dq]))
    } else {
      dq_df <- NULL
    }
  } else {
    ipw_owner_refuse_decorative_attrs(declared_qualified, "declared_qualified")
    dq_df <- declared_qualified
  }

  row_cohorts <- unique(as.double(c_rows$g))
  in_design <- ledger$g %in% row_cohorts
  tpt_df <- cdid_fast_df(list(g = ledger$g[in_design], t = ledger$time[in_design]))

  pipeline_metadata_require_stage(c("owners", "occupancy", "catalog", "design"))
  occ_obj <- meta_initial$occupancy_object
  catalog <- meta_initial$catalog
  pm_ctrl <- meta_initial$control
  own_ctrl <- pm_ctrl$owner_control
  cat_bnd <- pm_ctrl$catalog_bounds
  pre1 <- meta_initial$preflight
  chp1 <- pre1$catalog_header_preflight

  n_occ <- nrow(occ_obj$occupancy)
  H_occ <- chp1$n_calendar_dates + chp1$n_cohorts + chp1$n_targets + n_occ + chp1$n_requested_post + 1
  b_conv <- chp1$convention_bytes
  hdr_scr <- 64 * H_occ + 16 * b_conv
  hdr_wrk <- 64 * (1 + ceiling(log2(H_occ + 1))) * H_occ + 16 * b_conv
  occ_ret <- occ_obj$preflight$total_retained_entries
  if (hdr_wrk > cat_bnd$max_work || hdr_scr > cat_bnd$max_scratch ||
      occ_ret > min(own_ctrl$max_retained_entries, own_ctrl$max_output_entries) ||
      occ_ret + hdr_scr > own_ctrl$max_matrix_entries ||
      hdr_wrk > own_ctrl$max_work_ops) {
    return(NULL)
  }

  J_ad <- chp1$n_targets
  lower <- list(
    R = 0, C = as.double(chp1$n_cohorts), T = as.double(chp1$n_calendar_dates),
    P = 0, NP = 0, PR = 0, J = as.double(J_ad),
    D = if (is.data.frame(dq_df)) as.double(nrow(dq_df)) else 0,
    Q = as.double(nrow(tpt_df)),
    O = 0, B = b_conv, Bf = 0
  )
  env_ad <- catalog_design_envelope(lower, design = TRUE)
  if (env_ad$scratch > cat_bnd$max_scratch || env_ad$work > cat_bnd$max_work) {
    return(NULL)
  }

  acc <- meta_initial$diagnostics$mapping_accounting
  ad <- catalog_adapter_preflight(catalog, dq_df, tpt_df, design = TRUE)
  ad_links <- catalog_design_envelope(ad$dimensions, design = FALSE)
  design_scratch <- ipw_owner_safe_sum(c(ad$scratch, ad_links$scratch))
  ipw_owner_safe_sum(c(acc$retained_entries, design_scratch), min(own_ctrl$max_retained_entries, own_ctrl$max_output_entries), "retained catalog/mappings plus design output")
  ipw_owner_safe_sum(c(acc$concurrent_entries, design_scratch), own_ctrl$max_matrix_entries, "catalog/mappings plus design workspace")
  design_work_ops <- ipw_owner_safe_sum(c(acc$work_ops_bound, ad$work, ad_links$work), own_ctrl$max_work_ops, "composed metadata/design work") - acc$work_ops_bound

  handoff <- catalog_design_handoff(catalog, dq_df)
  design_obj <- build_catalog_handoff_design(
    handoff,
    target_post_times = tpt_df,
    pt_assumption = pt_assumption,
    pt_start = pt_start
  )
  cdid_fast_kernel_cache$last_validated_design <- design_obj
  design_links <- design_obj$links
  handoff <- NULL

  chp2 <- chp1
  chp2$n_occupancy_lower <- n_occ
  chp2$n_occupancy_upper <- n_occ
  chp2$min_header_scratch <- hdr_scr
  chp2$max_header_scratch <- hdr_scr
  chp2$min_header_work <- hdr_wrk
  chp2$max_header_work <- hdr_wrk

  cons_names <- c("targets", "anticipation", "control_group", "conventions", "requested_post", "pt_assumption", "pt_start")
  reuse_binding <- list(
    panel_ref = meta_initial$panel_ref,
    frame_ref = meta_initial$frame_ref,
    occupancy_ref = meta_initial$occupancy_ref,
    source_metadata_ref = meta_initial$ref,
    arguments_identical = cons_names,
    control_identical = TRUE
  )
  pre2 <- list(
    occupancy_supplied = TRUE,
    occupancy_preflight = occ_obj$preflight,
    catalog_header_preflight = chp2,
    owner_control = own_ctrl,
    catalog_bounds = cat_bnd,
    build_design = TRUE,
    endpoints_read = 0L,
    peak_rss_promised = FALSE,
    catalog_composition_reservation = pre1$catalog_composition_reservation,
    catalog_reuse_binding = reuse_binding
  )

  diag2 <- meta_initial$diagnostics
  diag2$catalog_reuse <- list(
    reused = TRUE,
    catalog_constructed = FALSE,
    mappings_constructed = FALSE,
    borrowed_catalog_mapping_entries = acc$retained_entries,
    owned_catalog_mapping_entries = 0,
    owned_design_entries = design_scratch,
    owned_work_ops = design_work_ops,
    borrowed_work_ops = acc$work_ops_bound,
    source_metadata_ref = meta_initial$ref
  )

  meta_obj <- structure(
    list(
      metadata_id = NULL,
      panel_ref = meta_initial$panel_ref,
      frame_ref = meta_initial$frame_ref,
      occupancy_ref = meta_initial$occupancy_ref,
      occupancy = meta_initial$occupancy,
      cell_indices = meta_initial$cell_indices,
      occupancy_object = occ_obj,
      catalog = catalog,
      row_mapping = meta_initial$row_mapping,
      operator_mapping = meta_initial$operator_mapping,
      design_links = design_links,
      design = design_obj,
      owner_summary = meta_initial$owner_summary,
      control = pm_ctrl,
      construction = meta_initial$construction,
      catalog_source_ref = meta_initial$ref,
      preflight = pre2,
      diagnostics = diag2
    ),
    class = "cdid_pipeline_metadata"
  )
  sealed_meta <- fk_seal(meta_obj, "pipeline_metadata")
  if (is.null(sealed_meta)) sealed_meta <- ipw_owner_seal(meta_obj, "pipeline_metadata")
  list(meta_initial = meta_initial, metadata_obj = sealed_meta, target_post_times = tpt_df)
}

cdid_fast_dr_separate_sparse_replay_plan <- function(backend, links, B) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "dr_separate_baseline_replay_plan", "dr_separate_sparse_plan",
    "dr_separate_baseline_replay_counts", "dr_separate_sparse_prospective",
    "or_outcome_block_fit", "or_pooled_wls_gram_solve", "or_wls_gram_solve",
    "cp04_solve", "cp04_payload", "cp04_payload_native_kernel", "ipw_owner_seal"
  ))) {
    return(NULL)
  }
  fk_seal <- cdid_fast_kernel("C_cdid_fast_seal")
  if (is.null(fk_seal)) return(NULL)

  validate_dr_baseline_backend(backend)
  if (!dr_separate(backend) || isTRUE(backend$uses_selection_models) ||
      !is.null(backend$reference_ref) || !is.null(backend$separate_complete)) {
    return(NULL)
  }
  if (!or_baseline_native_index(links, length(backend$links)) || !gmm_scalar_numeric(B) ||
      B < 1 || B != floor(B) || B > backend$control$max_chunk_rows) {
    return(NULL)
  }
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state
  N <- st$N
  if (B > N) return(NULL)

  pre <- cdid_fast_if_budget(backend, B, links, "dr")
  if (is.null(pre)) return(NULL)
  if (!is.null(names(pre$sids)) || !is.null(names(pre$uids)) ||
      !is.null(names(pre$gids)) || !is.null(names(pre$tids))) {
    return(NULL)
  }

  counts <- dr_separate_baseline_replay_counts(backend, links, B, pre, N)
  joint <- dr_separate_sparse_prospective(backend, links, B, pre, counts)
  if (is.null(joint)) return(NULL)

  nl <- counts$nl; ns <- counts$ns; ng <- counts$ng; nt <- counts$nt; pt <- counts$pt
  p <- counts$p; control_n <- counts$control_n; model_b <- counts$model_b; model_n <- counts$model_n
  R <- counts$R; C <- counts$C; M <- counts$M; I <- counts$I; per_block <- counts$per_block
  retained11 <- counts$retained; peak11 <- counts$peak; min_B <- counts$min_B; bound11 <- counts$bound

  starts <- function(n) if (length(n) > 0L) cumsum(c(1, n[-length(n)])) else numeric(0L)
  model_order <- order(model_b, seq_along(model_b))
  r_n <- p[model_b]; c_n <- control_n[model_b]
  r_start <- c_start <- numeric(nt)
  r_start[model_order] <- starts(r_n[model_order])
  c_start[model_order] <- starts(c_n[model_order])

  lk_models <- integer(nl)
  b_links <- backend$links
  for (i in seq_len(nl)) lk_models[i] <- b_links[[links[i]]]$model
  requested_models <- match(lk_models, pre$tids)
  if (anyNA(requested_models)) return(NULL)
  link_n <- tabulate(requested_models, nbins = nt)
  link_order <- order(requested_models, seq_along(requested_models))

  block_start <- starts(per_block)
  block_count <- per_block
  control_start <- starts(control_n)
  control_count <- control_n
  control_pos <- numeric(I)
  link_start <- starts(link_n)
  link_count <- link_n
  adj_start <- r_start
  or_adj <- numeric(R)
  or_center <- numeric(nt)
  or_up_start <- c_start
  or_up <- numeric(C)
  eta_mat <- matrix(0, pt, nt)
  support_start <- starts(model_n)
  support_count <- model_n
  support_pos <- numeric(M)
  prop_up <- numeric(M)

  og <- backend$outcome_groups
  og_groups <- og$groups
  og_bm_gid <- og$block_map$group_id
  has_icpt <- isTRUE(backend$outcome_recipe$intercept)
  norm_N <- as.double(og$N)
  reg_specs <- backend$regression_specs
  regs <- backend$regressions
  trts <- backend$treatments
  trt_specs <- backend$treatment_specs
  pre_sids <- pre$sids
  pre_gids <- pre$gids
  pre_tids <- pre$tids

  P_sum <- sum(p)
  start_int <- as.integer(cumsum(c(1, if (ng > 1L) p[-ng] else numeric(0L))))
  block_par <- matrix(0, 2L, ng)
  keep_vec <- integer(P_sum)
  beta_vec <- numeric(P_sum)
  mean_d_vec <- numeric(P_sum)
  or_fidx <- backend$or_frame_index
  slope_idx <- backend$outcome_recipe$slope_index

  for (kb in seq_len(ng)) {
    b <- pre_gids[kb]
    cnt_b <- block_count[kb]
    models <- if (cnt_b > 0L) model_order[seq.int(block_start[kb], length.out = cnt_b)] else integer(0L)
    tids <- pre_tids[models]
    controls <- reg_specs[[b]]$support_ids
    n_ctrl <- control_count[kb]
    if (n_ctrl > 0L) {
      control_pos[seq.int(control_start[kb], length.out = n_ctrl)] <- match(controls, pre_sids)
    }
    g_blk <- og_groups[[og_bm_gid[b]]]
    rec_blk <- g_blk$regressions[[1L]]
    p_b <- length(g_blk$keep)
    block_par[1L, kb] <- rec_blk$L
    block_par[2L, kb] <- rec_blk$alpha
    if (p_b > 0L) {
      ii_b <- seq.int(start_int[kb], length.out = p_b)
      keep_vec[ii_b] <- or_fidx[slope_idx[g_blk$keep]]
      beta_vec[ii_b] <- g_blk$beta
      mean_d_vec[ii_b] <- rec_blk$mean_d
    }
    if (p_b == 0L || length(models) == 0L) next
    n_mod <- length(models)
    loads <- matrix(0.0, nrow = p_b, ncol = n_mod)
    for (j in seq_len(n_mod)) {
      loads[, j] <- as.double(trts[[tids[j]]]$d_bar_diff)
    }
    fit_b <- g_blk$fit
    sfit <- if (has_icpt) fit_b$slope_fit else fit_b
    if (!is.list(sfit) || !is.matrix(sfit$R) || !is.integer(sfit$pivot) ||
        nrow(sfit$R) != p_b || ncol(sfit$R) != p_b || length(sfit$pivot) != p_b ||
        anyNA(sfit$pivot) || any(sfit$pivot < 1L | sfit$pivot > p_b)) {
      return(NULL)
    }
    d_R <- diag(sfit$R)
    if (any(!is.finite(d_R) | d_R == 0)) return(NULL)
    piv_b <- sfit$pivot
    rhs_p <- loads[piv_b, , drop = FALSE]
    if (any(!is.finite(rhs_p))) return(NULL)
    z_sol <- forwardsolve(sfit$R, rhs_p, upper.tri = TRUE, transpose = TRUE)
    u_sol <- backsolve(sfit$R, z_sol)
    adj <- matrix(0.0, nrow = p_b, ncol = n_mod)
    adj[piv_b, ] <- u_sol * norm_N
    if (any(!is.finite(adj))) return(NULL)
    rec <- regs[[b]]
    or_center[models] <- as.vector(crossprod(adj, as.double(rec$score_mean_beta)))
    up_mat <- crossprod(rec$forcing_beta, adj)
    for (j in seq_len(n_mod)) {
      k <- models[j]
      or_adj[seq.int(adj_start[k], length.out = p_b)] <- adj[, j]
      if (n_ctrl > 0L) {
        or_up[seq.int(or_up_start[k], length.out = n_ctrl)] <- up_mat[, j]
      }
    }
  }

  model_block <- integer(nt)
  trt_beta <- matrix(0, pt, nt)
  means_mat <- matrix(0, 4L, nt)
  sel_vec <- numeric(M)

  for (k in seq_len(nt)) {
    tid <- pre_tids[k]
    f <- trts[[tid]]
    spec_t <- trt_specs[[tid]]
    H_c <- f$H_chol
    eta <- as.vector(backsolve(H_c, forwardsolve(H_c, f$partials$beta[2L, ], upper.tri = TRUE, transpose = TRUE)))
    sc_k <- support_count[k]
    ii <- if (sc_k > 0L) seq.int(support_start[k], length.out = sc_k) else integer(0L)
    eta_mat[, k] <- eta
    if (sc_k > 0L) {
      support_pos[ii] <- match(spec_t$support_ids, pre_sids)
      prop_up[ii] <- as.vector(crossprod(f$cross, eta))
      sel_vec[ii] <- f$partials$upstream[1L, ] - f$partials$upstream[2L, ]
    }
    model_block[k] <- match(spec_t$regression, pre_gids)
    trt_beta[, k] <- f$coefficients
    means_mat[, k] <- f$means[c("T_dr", "M_dr", "dp", "dm")]
  }
  if (anyNA(control_pos) || anyNA(support_pos) || anyNA(model_block)) return(NULL)

  fk_pay <- cdid_fast_kernel("C_gmm_cp08_payload")
  v11 <- list(
    kind = "dr_separate_baseline_replay_plan", backend_ref = backend$ref,
    links = as.integer(links), B = B, min_B = min_B, N = N,
    sids = pre_sids, uids = pre$uids, gids = pre_gids, tids = pre_tids,
    p = p, model_order = model_order, block_start = block_start, block_count = block_count,
    control_start = control_start, control_count = control_count, control_pos = control_pos,
    link_order = link_order, link_start = link_start, link_count = link_count,
    adj_start = adj_start, or_adj = or_adj, or_center = or_center,
    or_up_start = or_up_start, or_up = or_up, eta = eta_mat,
    support_start = support_start, support_count = support_count, support_pos = support_pos,
    prop_up = prop_up, retained_bound = retained11, construction_peak = peak11,
    reservation = bound11, entries = 0
  )
  v11_entries <- unname((if (!is.null(fk_pay)) fk_pay(v11, TRUE, environment()) else cp04_payload(v11))["entries"])
  v11$entries <- v11_entries + 64
  if (v11$entries > retained11) return(NULL)
  v11 <- fk_seal(v11, "dr_separate_baseline_replay_plan")
  if (is.null(v11)) return(NULL)

  sup_specs <- backend$support_specs
  sups <- backend$supports
  raw_dates <- integer(2L * ns)
  for (j in seq_len(ns)) {
    pr <- sup_specs[[pre_sids[j]]]$pair
    raw_dates[2L * j - 1L] <- pr$s_calendar_index
    raw_dates[2L * j] <- pr$t_calendar_index
  }
  dates <- sort(unique(raw_dates))
  if (length(dates) > joint$nd) return(NULL)

  support_mat <- matrix(0L, 3L, ns)
  support_par <- matrix(0, 3L, ns)
  for (j in seq_len(ns)) {
    sid <- pre_sids[j]
    s <- sup_specs[[sid]]
    u <- sups[[sid]]
    support_mat[1L, j] <- s$cohort
    support_mat[2L, j] <- match(s$pair$s_calendar_index, dates)
    support_mat[3L, j] <- match(s$pair$t_calendar_index, dates)
    support_par[1L, j] <- u$rate
    support_par[2L, j] <- u$massA
    support_par[3L, j] <- u$v
  }

  block_mat <- rbind(start = start_int, as.integer(p))
  storage.mode(block_mat) <- "integer"
  data <- list(
    intercept = as.integer(has_icpt),
    support = support_mat, support_par = support_par,
    block = block_mat, block_par = block_par,
    keep = keep_vec, beta = beta_vec, mean_d = mean_d_vec, model_block = model_block,
    trt_beta = trt_beta, means = means_mat, selection = sel_vec,
    trt_columns = as.integer(backend$treatment_recipe$index)
  )

  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  post <- panel$postings
  bounds <- rbind(post$date_start[dates], post$date_end[dates])
  bounds[, post$date_count[dates] == 0L] <- 0L
  storage.mode(bounds) <- "integer"

  plan <- list(
    kind = "dr_separate_sparse_plan", backend_ref = backend$ref,
    links = as.integer(links), B = B, min_B = min_B, N = N,
    sids = pre_sids, uids = pre$uids, gids = pre_gids, tids = pre_tids,
    dates = dates, bounds = bounds, data = data, v11 = v11,
    metadata = joint$metadata, retained_bound = joint$retained + retained11,
    construction_peak = joint$peak + peak11 + joint$count_reservation,
    reservation = joint$reservation, entries = 0
  )
  plan_entries <- unname((if (!is.null(fk_pay)) fk_pay(plan, TRUE, environment()) else cp04_payload(plan))["entries"])
  plan$entries <- plan_entries + 64
  if (plan$entries > plan$retained_bound) return(NULL)
  fk_seal(plan, "dr_separate_sparse_plan")
}

cdid_fast_or_dr_baseline_diagnostics <- function(res, control, public_object = NULL) {
  if (!identical(parent.env(parent.frame()), environment(cdid_fast_kernel))) return(NULL)
  if (!cdid_are_default_fns(c(
    "or_public_diagnostics", "dr_public_diagnostics", "cp04_field_diagnostics",
    "cp04_public_diagnostic_retained", "cp04_budget", "outcome_public_blocks",
    "or_outcome_block_fit", "outcome_public_mode", "outcome_public_model",
    "cp04_payload", "cp04_payload_native_kernel", "cp04_public_diagnostic_check"
  ))) {
    return(NULL)
  }
  b <- res$nuisance_registry$cp04
  if (is.null(b) || !is.null(b$reference_ref)) return(NULL)
  is_or <- inherits(b, "cdid_or_baseline")
  is_dr <- inherits(b, "cdid_dr_baseline")
  if (!is_or && !is_dr) return(NULL)
  if (cp04_field_diagnostics(b)) return(NULL)

  label <- if (is_dr) "DR public diagnostic copy" else "OR public diagnostic copy"
  bound <- 4 * b$accounting$owned_entries + 4096
  own_ctrl <- control$owner_control
  if (!is.list(own_ctrl) || own_ctrl$max_retained_entries < 1e8 ||
      own_ctrl$max_matrix_entries < 1e8 || own_ctrl$max_work_ops < 1e8 ||
      own_ctrl$max_output_entries < 1e8 || bound > 5e6 ||
      !is.list(res$frame_owner) || res$frame_owner$n_a_units > 50000L) {
    cp04_budget(
      own_ctrl,
      c(cp04_public_diagnostic_retained(res, public_object), bound),
      4 * bound, 8 * bound, bound, label
    )
  }

  is_sep <- identical(b$outcome_pooling, "separate")
  local_blocks <- NULL
  if (is_sep) {
    og <- if (!is.null(b$outcome_groups)) b$outcome_groups else b$shared_or$outcome_groups
    if (!is.null(og)) {
      bm <- og$block_map
      n_bm <- nrow(bm)
      g_ids <- bm$group_id
      lb_ids <- bm$local_block_id
      grps <- og$groups
      local_blocks <- vector("list", n_bm)
      for (bi in seq_len(n_bm)) {
        gid <- g_ids[bi]
        lbid <- lb_ids[bi]
        grp <- grps[[gid]]
        reg_b <- grp$regressions[[lbid]]
        local_blocks[[bi]] <- list(
          block_id = bi,
          group_id = gid,
          local_block_id = lbid,
          keep = grp$keep,
          columns = grp$basis$kept_columns,
          alpha = reg_b$alpha,
          beta = grp$beta,
          basis = grp$basis,
          fit = grp$fit
        )
      }
      names(local_blocks) <- as.character(seq_len(n_bm))
    }
  }

  sh_names <- c("p", "keep", "beta", "coefficients", "fit", "summary",
                "basis", "control_support_ids", "rss", "objective", "verification_scan")
  shared <- b$shared[sh_names[sh_names %in% names(b$shared)]]

  sels <- b$selections
  sel_specs <- b$selection_specs
  n_sel <- length(sels)
  sel_recs <- vector("list", n_sel)
  for (j in seq_len(n_sel)) {
    sel_recs[[j]] <- c(list(id = j, spec = sel_specs[[j]]), sels[[j]])
  }

  regs <- b$regressions
  reg_specs <- b$regression_specs
  n_reg <- length(regs)
  reg_recs <- vector("list", n_reg)
  block_keys <- character(n_reg)
  local_icpts <- numeric(n_reg)
  reg_L <- numeric(n_reg)
  reg_role <- if (is_sep) "distinct_comparison_outcome_fit" else "local_intercept_block_of_single_joint_fit"
  for (j in seq_len(n_reg)) {
    reg <- regs[[j]]
    spec <- reg_specs[[j]]
    reg_compact <- reg[names(reg) != "qr"]
    reg_recs[[j]] <- c(list(id = j, role = reg_role, spec = spec), reg_compact)
    bk <- spec$block_key
    if (is.null(bk)) bk <- spec$key
    block_keys[j] <- if (is.null(bk)) "" else bk
    local_icpts[j] <- as.double(reg$alpha)
    reg_L[j] <- as.double(reg$L)
  }
  names(local_icpts) <- block_keys

  if (is_dr) {
    trts <- b$treatments
    trt_specs <- b$treatment_specs
    n_trt <- length(trts)
    trt_recs <- vector("list", n_trt)
    trt_coefs <- vector("list", n_trt)
    for (j in seq_len(n_trt)) {
      tr <- trts[[j]]
      trt_recs[[j]] <- c(list(id = j, spec = trt_specs[[j]]), tr)
      trt_coefs[[j]] <- tr$coefficients
    }
    records <- list(selection = sel_recs, regression = reg_recs, treatment = trt_recs)
    coefficients <- list(
      shared_slopes = b$shared$coefficients,
      block_slopes = if (is.null(local_blocks)) NULL else lapply(local_blocks, `[[`, "beta"),
      local_intercepts = local_icpts,
      has_local_intercepts = isTRUE(b$outcome_recipe$intercept),
      treatment = trt_coefs
    )
  } else {
    records <- list(selection = sel_recs, regression = reg_recs)
    coefficients <- list(
      shared_slopes = b$shared$coefficients,
      block_slopes = if (is.null(local_blocks)) NULL else lapply(local_blocks, `[[`, "beta"),
      local_intercepts = local_icpts,
      has_local_intercepts = isTRUE(b$outcome_recipe$intercept)
    )
  }

  sups <- b$supports
  sup_specs <- b$support_specs
  n_sup <- length(sups)
  supports_diag <- vector("list", n_sup)
  sup_s <- numeric(n_sup)
  sup_t <- numeric(n_sup)
  for (j in seq_len(n_sup)) {
    sp <- sup_specs[[j]]
    supports_diag[[j]] <- c(list(id = j, spec = sp), sups[[j]])
    sup_s[j] <- as.double(sp$pair$s)
    sup_t[j] <- as.double(sp$pair$t)
  }

  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  st_cohorts <- reader$state$cohorts
  lks <- b$links
  n_lks <- length(lks)
  link_ids <- character(n_lks)
  lk_coh <- integer(n_lks)
  lk_reg <- integer(n_lks)
  lk_sup <- integer(n_lks)
  lk_mu <- numeric(n_lks)
  if (is_dr) {
    lk_mod <- integer(n_lks)
    lk_mu_dr <- numeric(n_lks)
    lk_mu_or <- numeric(n_lks)
    lk_mu_ipw <- numeric(n_lks)
  }
  for (i in seq_len(n_lks)) {
    l <- lks[[i]]
    link_ids[i] <- l$link_id
    lk_coh[i] <- l$cohort
    lk_reg[i] <- l$regression
    lk_sup[i] <- l$support_id
    lk_mu[i] <- l$mu
    if (is_dr) {
      lk_mod[i] <- l$model
      lk_mu_dr[i] <- l$mu_dr
      lk_mu_or[i] <- l$mu_or
      lk_mu_ipw[i] <- l$mu_ipw
    }
  }
  roles <- rep.int("catalog", n_lks)
  aux_mask <- startsWith(link_ids, "cp04:auxiliary:")
  if (any(aux_mask)) roles[aux_mask] <- "auxiliary"

  fitted_catalog <- if (is_dr) {
    cdid_fast_df(list(
      index = seq_len(n_lks),
      link_id = link_ids,
      role = roles,
      cohort = st_cohorts[lk_coh],
      s = sup_s[lk_sup],
      t = sup_t[lk_sup],
      duration_steps = reg_L[lk_reg],
      regression_block = lk_reg,
      treatment_model = lk_mod,
      support_id = lk_sup,
      mu = lk_mu,
      mu_dr = lk_mu_dr,
      mu_or = lk_mu_or,
      mu_ipw = lk_mu_ipw
    ))
  } else {
    cdid_fast_df(list(
      index = seq_len(n_lks),
      link_id = link_ids,
      role = roles,
      cohort = st_cohorts[lk_coh],
      s = sup_s[lk_sup],
      t = sup_t[lk_sup],
      duration_steps = reg_L[lk_reg],
      regression_block = lk_reg,
      support_id = lk_sup,
      mu = lk_mu
    ))
  }

  diag_copy <- list(
    records = records, shared = shared, outcome_blocks = local_blocks,
    coefficients = coefficients, fitted_catalog = fitted_catalog, supports = supports_diag
  )
  fk_pay <- cdid_fast_kernel("C_gmm_cp08_payload")
  diag_entries <- unname((if (!is.null(fk_pay)) fk_pay(diag_copy, TRUE, environment()) else cp04_payload(diag_copy))["entries"])
  if (diag_entries > bound) {
    return(NULL)
  }

  out <- if (is_dr) {
    list(
      status = "available", method = "dr",
      outcome_model = "shared one-period outcome-change control-mean regression",
      treatment_model = "risk-set treatment propensity logit",
      fit_structure = "single_joint_or_fit_with_block_local_intercepts_plus_deduplicated_treatment_logits",
      duration_convention = "stored-calendar transition count L_e = i_t - i_s",
      selection = "none",
      selection_model = NULL,
      selection_pooling = NULL,
      selection_route = res$selection,
      outcome_recipe = b$outcome_recipe,
      treatment_recipe = b$treatment_recipe,
      fitted_catalog = fitted_catalog,
      coefficients = coefficients,
      basis = list(common_design = res$controls$or_basis, pooled_control = b$shared$basis),
      shared = shared, outcome_pooling = outcome_public_mode(b$outcome_pooling), outcome_blocks = local_blocks,
      records = records,
      supports = supports_diag,
      links = b$links,
      accounting = b$accounting,
      frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
        n_pair_eligible_A_units = res$frame_owner$n_a_units),
      auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
      interpretation = paste("Treatment-DR extension of the shared one-period outcome-change control-mean regression",
        "and residualized Hajek IPW over deduplicated risk-set treatment propensity logits under the declared observation route;",
        "complete F-row influence functions with global OR slope, treatment propensity, and observation/selection paths.",
        "No numerical certificate.")
    )
  } else {
    list(
      status = "available", method = "reg",
      outcome_model = "shared one-period outcome-change control-mean regression",
      fit_structure = "single_joint_fit_with_block_local_intercepts",
      duration_convention = "stored-calendar transition count L_e = i_t - i_s",
      selection = "none",
      selection_model = NULL,
      selection_pooling = NULL,
      selection_route = res$selection,
      outcome_recipe = b$outcome_recipe,
      fitted_catalog = fitted_catalog,
      coefficients = coefficients,
      basis = list(common_design = res$controls$or_basis, pooled_control = b$shared$basis),
      shared = shared, outcome_pooling = outcome_public_mode(b$outcome_pooling), outcome_blocks = local_blocks,
      records = records,
      supports = supports_diag,
      links = b$links,
      accounting = b$accounting,
      frame_counts = list(n_full_frame_units = res$reporting$inference$n_total,
        n_pair_eligible_A_units = res$frame_owner$n_a_units),
      auxiliary_descriptor_status = res$auxiliary_reporting$descriptor_status,
      interpretation = paste("Single joint weighted least-squares fit of shared one-period outcome-change slopes",
        "across all catalog and auxiliary control blocks with local block intercepts (when the outcome formula includes an intercept)",
        "and stored-calendar transition durations L_e = i_t - i_s; block records describe local intercept blocks of that single joint fit,",
        "not independent regression fits. Complete F-row influence functions with global slope dependencies. No treatment fits or numerical certificate.")
    )
  }

  if (!is.null(local_blocks)) {
    out$outcome_model <- outcome_public_model("separate")
    out$fit_structure <- "one outcome fit per unique comparison block"
    out$basis$blocks <- lapply(local_blocks, `[[`, "basis")
    out$interpretation <- paste("Comparison-specific outcome regressions with complete original-unit nuisance influence functions;",
      "local slopes and prediction-preserving bases remain unrestricted across distinct comparison blocks.")
  }
  out
}

