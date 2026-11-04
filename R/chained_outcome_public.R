# Chained-only choice. Inactive callers do not evaluate the user's promise.
chained_request_outcome_pooling <- function(value) {
  chained_request_need(chained_request_chr1(value) && value %in% c("separate", "pooled"),
    '`outcome_pooling` must be "separate" (default) or "pooled".',
    "cdid_chained_request_outcome_pooling")
  value
}

outcome_public_mode <- function(x) if (identical(x, "separate")) "separate" else "pooled"
outcome_public_model <- function(x) if (identical(x, "separate"))
  "comparison-specific one-period outcome-change control-mean regressions" else
    "shared one-period outcome-change control-mean regression"

outcome_public_blocks <- function(backend) {
  object <- if (!is.null(backend$outcome_groups)) backend$outcome_groups else backend$shared_or$outcome_groups
  if (is.null(object)) return(NULL)
  out <- lapply(seq_len(nrow(object$block_map)), function(b) {
    f <- or_outcome_block_fit(object, b)
    list(block_id = b, group_id = f$group_id, local_block_id = f$local_block_id,
      keep = f$keep, columns = f$basis$kept_columns, alpha = f$alpha, beta = f$beta,
      basis = f$basis, fit = f$fit)
  })
  names(out) <- as.character(seq_along(out))
  out
}
