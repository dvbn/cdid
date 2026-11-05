# Ordinary conditional separate DR keeps all upstream rows/hits, one OR block
# result and one current treatment. Only accounting changes; its scalar core,
# full Gamma/direct coordinates, checks and order remain in chained_dr_separate.
dr_conditional_separate_ordinary_inventory <- function(backend, B, links, pre, selection_pre, maps) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  sel <- backend$selection_backend; map <- backend$up_map
  prod <- ipw_owner_safe_product; add <- ipw_owner_safe_sum
  p <- backend$outcome_recipe$p; pt <- backend$treatment_recipe$p
  ns <- length(pre$sids); nu <- length(pre$uids); ng <- length(pre$gids)
  nt <- length(pre$tids); nl <- length(links); M <- length(sel$blocks)
  active <- which(vapply(map$theta_coords, function(ids) any(ids %in% pre$uids), logical(1)))
  W <- cells <- 0
  for (m in active) {
    bl <- sel$blocks[[m]]
    W <- add(c(W, if (bl$has_intercept) bl$B_m else 0, bl$d_s, length(map$theta_coords[[m]])))
    cells <- add(c(cells, length(bl$cell_ids)))
  }
  ref <- chained_backend_reference(backend)
  ref_x <- if (is.null(ref)) 0 else length(ref$x_columns)
  ref_s <- if (is.null(ref)) 0 else length(ref$designs[[ref$role_slot["s"]]]$columns)
  d <- max(1L, sel$design_spec$d)
  or_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
  trt_bytes <- sum(nchar(backend$treatment_recipe$columns, type = "bytes"))
  design_bytes <- sum(nchar(if (is.null(ref)) colnames(st$Z) else
    c(ref$x_columns, ref$designs[[ref$role_slot["s"]]]$columns), type = "bytes"))
  link_bytes <- sum(vapply(backend$links[links], function(l) nchar(l$link_id, type = "bytes"), integer(1)))
  # Count every copied requested-link record and its characters, not a full
  # backend copy. Sealed owners/factors are still in owned/borrowed accounting.
  requested_payload <- unname(cp04_payload(backend$links[links])["entries"])
  attributes_payload <- unname(cp04_payload(attributes(links))["entries"])
  output_meta <- add(c(8192, 256 * nl, 8 * link_bytes,
    B * frame$id_byte_metadata$max, attributes_payload))
  output <- add(c(output_meta, prod(B, nl + 4)))
  # up_res retains ALL alpha/beta/theta columns in each selected block, plus
  # the requested up matrix and empty slots for every global block.
  upstream <- prod(B, add(c(nu, W)))
  upstream_meta <- add(c(1024, 64 * (M + cells + length(active) + nu)))
  # Each support keeps obs/dy/raw-p-column z/v and possibly recipe weights.
  # Matrix/list headers and repeated raw outcome column names are explicit.
  hits <- prod(B * ns, p + 4)
  hits_meta <- add(c(1024, 256 * ns, prod(8 * ns, or_bytes)))
  metadata <- add(c(maps, 8192, 64 * (ns + nu + ng + nt + nl + M + cells + 1),
    2 * requested_payload, 8 * link_bytes))
  # The last support reader's rows/hit/Z/probability locals remain live after
  # the hits loop. Prediction temporaries/full reference designs are separate.
  support_state <- add(c(prod(B, p + 16), 1024, 16 * (or_bytes + design_bytes)))
  # V37's operation-local methods/closures/snapshot bindings hold aliases,
  # not score/probability arrays. Keep their independently inventoried slab
  # dedicated, separate from every map/operator/row term below.
  context_metadata <- 2048
  common <- add(c(metadata, context_metadata, upstream, upstream_meta, hits, hits_meta, support_state))
  read_peak <- add(c(prod(32 * B, p + d + ref_x + ref_s + 16),
    4096, 64 * (or_bytes + design_bytes)))

  link_blocks <- vapply(backend$links[links], `[[`, integer(1), "regression")
  core_return <- core_peak <- 0
  max_links <- max_channels <- max_local_p <- max_supports <- 0
  for (b in pre$gids) {
    pp <- length(backend$outcome_groups$groups[[b]]$keep)
    ll <- sum(link_blocks == b); rr <- length(backend$regressions[[b]]$channel_ids)
    ss <- length(backend$regression_specs[[b]]$support_ids)
    h <- 0
    for (s in backend$regression_specs[[b]]$support_ids) h <- max(h, length(backend$supports[[s]]$channels))
    for (j in which(link_blocks == b)) h <- max(h, length(backend$supports[[backend$links[[links[j]]]$support_id]]$channels))
    # Returned phi/direct/den/num, every local B-by-p design, coefficient
    # beta/beta_score and alpha/alpha_score, and objective stay through models.
    returned <- add(c(prod(B, 4 * ll + ll * pp + 2 * pp + 4),
      2048, 256 * (ll + 1), prod(16 * (ll + 4), or_bytes)))
    # ha/energy/hb/rhs/xb/xa, us/ul, one residual record, design products,
    # sweep/assignment copies and full normalized Gram RHS/solve/reorder.
    construction <- add(c(returned, prod(B, 20 * pp + 4 * rr + 4 * h + 48),
      prod(2 * pp, rr), prod(32 * pp, pp),
      4096, 128 * (ss + rr + h + 1), prod(16 * (ll + 8), or_bytes)))
    core_return <- max(core_return, returned); core_peak <- max(core_peak, construction)
    max_links <- max(max_links, ll); max_channels <- max(max_channels, rr)
    max_local_p <- max(max_local_p, pp); max_supports <- max(max_supports, ss)
  }
  model_state <- model_peak <- 0; max_u <- 0
  for (t in pre$tids) {
    f <- backend$treatments[[t]]; u <- length(f$upstream_ids)
    pp <- length(backend$outcome_groups$groups[[backend$treatment_specs[[t]]$regression]]$keep)
    # Surviving us, score, Gamma, last X/d/residual/logit/weight/row arrays,
    # both four-coordinate direct matrices and pd/po/pp/ps/value vectors.
    state <- add(c(prod(B, u + 3 * pt + pp + 32),
      2048, 64 * (u + pt + pp + 1), 16 * (or_bytes + trt_bytes + design_bytes)))
    # Replacement locals, score+cross RHS, transpose/triangular solves and
    # full Gamma, both direct arrays, current full reference X/S and products.
    construction <- add(c(state, prod(B, 16 * pt + 4 * pp + 4 * u +
      2 * ref_x + 2 * ref_s + 4 * d + 64), prod(16 * pt, pt), prod(4 * pt, u),
      4096, 64 * (or_bytes + trt_bytes + design_bytes)))
    model_state <- max(model_state, state); model_peak <- max(model_peak, construction)
    max_u <- max(max_u, u)
  }
  selection_peak <- add(c(selection_pre$output, selection_pre$scratch))
  # rm(core,bf) does NOT remove the last model's locals. Fund those through
  # the NEXT OR construction, and the current core return through each model.
  phase_peak <- max(selection_peak, read_peak, add(c(model_state, core_peak)),
    add(c(core_return, model_peak)))
  scratch <- add(c(common, phase_peak, output, prod(3 * B, nl)))
  list(output = output, scratch = scratch, common = common, metadata = metadata,
    context_metadata = context_metadata,
    upstream = upstream, upstream_metadata = upstream_meta, selection_columns = W,
    hits = hits, hits_metadata = hits_meta, support_state = support_state,
    selection_peak = selection_peak, read_peak = read_peak, core_return = core_return,
    core_peak = core_peak, previous_model_state = model_state, current_model_peak = model_peak,
    phase_peak = phase_peak, max_block_links = max_links, max_block_channels = max_channels,
    max_block_slopes = max_local_p, max_block_supports = max_supports,
    max_model_channels = max_u, reference_x_columns = ref_x, reference_s_columns = ref_s,
    link_attribute_entries = attributes_payload)
}

dr_conditional_separate_ordinary_if_envelope <- function(backend, B, links, pre, selection_pre, maps) {
  x <- dr_conditional_separate_ordinary_inventory(backend, B, links, pre, selection_pre, maps)
  pre$output <- x$output; pre$scratch <- x$scratch
  pre
}
