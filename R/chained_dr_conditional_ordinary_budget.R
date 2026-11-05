# Ordinary pooled Corollary1 DR returns only phi and row/link diagnostics.
# The unchanged evaluator still constructs full selection/Gamma/coordinate IFs.
dr_corollary1_ordinary_inventory <- function(backend, B, links, pre, selection_pre, metadata) {
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  st <- reader$state; frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  sel <- backend$selection_backend; map <- backend$up_map
  p <- backend$shared$p; pt <- backend$treatment_recipe$p; d <- max(1L, backend$upstream$d)
  ns <- length(pre$sids); nu <- length(pre$uids); nt <- length(pre$tids); nl <- length(links)
  Tn <- length(backend$treatments); M <- length(sel$blocks); nc <- length(pre$active_cells)
  prod <- ipw_owner_safe_product; add <- ipw_owner_safe_sum
  # Each active bif retains all alpha/beta/theta columns, even when uids asks
  # for only part of a pooled selection block. Empty global list slots remain.
  W <- 0
  for (m in pre$active_blocks) {
    bl <- sel$blocks[[m]]
    W <- add(c(W, if (bl$has_intercept) bl$B_m else 0, bl$d_s, length(map$theta_coords[[m]])))
  }
  I_beta <- 0
  if (p) for (b in pre$beta_gids) I_beta <- add(c(I_beta, length(backend$regression_specs[[b]]$support_ids)))
  max_u <- 0
  for (t in pre$tids) max_u <- max(max_u, length(backend$treatments[[t]]$upstream_ids))
  ref <- chained_backend_reference(backend)
  ref_x <- if (is.null(ref)) 0 else length(ref$x_columns)
  ref_s <- if (is.null(ref)) 0 else length(ref$designs[[ref$role_slot["s"]]]$columns)
  link_bytes <- sum(nchar(vapply(backend$links[links], `[[`, character(1), "link_id"), type = "bytes"))
  or_bytes <- sum(nchar(backend$outcome_recipe$columns, type = "bytes"))
  other_bytes <- sum(nchar(c(backend$treatment_recipe$columns,
    if (is.null(ref)) colnames(st$Z) else c(ref$x_columns, ref$designs[[ref$role_slot["s"]]]$columns)), type = "bytes"))
  index_attributes <- unname(cp04_payload(attributes(links))["entries"])
  output_meta <- add(c(8192, 256 * nl, 8 * link_bytes, B * frame$id_byte_metadata$max, index_attributes))
  output <- add(c(output_meta, prod(B, nl + 4)))
  # up_res survives every later phase; hits stores at most one row-width
  # outcome/design/weight record per support. Ordinary models cache FIVE row
  # vectors, not per-model global upstream or coefficient matrices.
  upstream <- prod(B, add(c(nu, W)))
  hits <- if (p) prod(B * ns, p + 4) else 0
  models <- add(c(prod(5 * B, nt), 256 * (nt + Tn + M + nc + ns + 1)))
  core_beta <- prod(B, p)
  common <- add(c(metadata, upstream, hits, models, core_beta))
  # beta_only retains every global (block,support) residual record inside the
  # core, but its cache/closure does not escape the list(beta=xi_beta) return.
  # Count all incidences even for a requested-link subset.
  core <- if (p) add(c(prod(B * I_beta, p + 3), prod(16 * B, p), prod(4 * B, nu),
    prod(8 * p, nu), prod(16 * p, p), 256 * I_beta,
    prod(8 * (ns + I_beta + 4), or_bytes), prod(64, other_bytes))) else 0
  # ONE current model: upstream slice, full score/RHS/transposes/triangular
  # solves/Gamma, direct and four-coordinate matrices/finite concatenation,
  # logits/weights/pair rows and full reference designs before local subsets.
  model <- add(c(prod(32 * B, add(c(p, pt, max_u, d, ref_x, ref_s, 32))),
    prod(16 * pt, pt), prod(8 * pt, max_u), prod(64, max_u), prod(64, or_bytes + other_bytes)))
  selection <- add(c(selection_pre$output, selection_pre$scratch))
  scratch <- add(c(common, max(selection, core, model), output, prod(3 * B, nl)))
  list(output = output, scratch = scratch, common = common, upstream = upstream, hits = hits,
    model_cache = models, core_beta = core_beta, selection_peak = selection, core_peak = core,
    current_model_peak = model, selection_columns = W, beta_support_incidence = I_beta,
    max_model_channels = max_u, reference_x_columns = ref_x, reference_s_columns = ref_s)
}

dr_corollary1_ordinary_if_envelope <- function(backend, B, links, pre, selection_pre, metadata) {
  x <- dr_corollary1_ordinary_inventory(backend, B, links, pre, selection_pre, metadata)
  pre$output <- x$output; pre$scratch <- x$scratch
  pre
}
