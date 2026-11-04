# Private local stage-2 outcome access. Stage-1 group fits and treatment roots
# keep their own identities and dimensions; no fitted-weight deduplication.
reconstruction_separate <- function(sh) !is.null(sh) && identical(sh$outcome_pooling, "separate")

reconstruction_or_block <- function(sh, bid) {
  if (!reconstruction_separate(sh)) return(sh)
  f <- or_outcome_block_fit(sh$outcome_groups, bid)
  list(keep = f$keep, p = f$p, kept_columns = f$basis$kept_columns,
    kept_cols_F = sh$slope_cols_F[f$keep], beta = f$beta, fit = f$fit,
    summary = sh$outcome_groups$groups[[f$group_id]]$summary, basis = f$basis,
    n_blocks = sh$n_blocks)
}

# Reserve all group operators and aliases before stage-2 fitting. The multiplier
# is deliberately conservative; it does not assert measured memory or speed.
reconstruction_bound <- function(x, n_blocks) ipw_owner_safe_sum(
  x * (n_blocks + 2), label = "separate reconstructed outcome envelope")
