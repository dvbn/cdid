# Each caller has already admitted its ORIGINAL envelope. This optional cache
# consumes only an unused compartment of that envelope; it never adds entries
# to a caller/composed budget or changes the issued row partition.
cp04_pair_consumer_operation <- function(backend, reader, rows, links, pre, kind, details) {
  if (details || !is.null(backend$reference_ref) || !is.null(backend$layout) ||
      !is.null(backend$pair_outcomes)) return(NULL)
  B <- length(rows); st <- reader$state
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  if (identical(kind, "ipw")) {
    p <- frame$n_covariates; d <- backend$selection_backend$design_spec$d
    width <- length(links) + length(pre$tids) * (2 * p + length(pre$uids) + 12)
    # Eight of the old 32 row copies. Twenty-four remain for every original
    # score/model/output buffer; selection_pre is an untouched extra reserve.
    fund <- 8 * B * (width + p + d + 8)
  } else if (identical(kind, "or")) {
    if (or_selection_separate(backend)) {
      # Original hits need <=B*ns*(p+4); leave sixteen copies of p+6
      # for hits/replacements and all other original OR terms untouched.
      fund <- 16 * B * length(pre$sids) * (backend$outcome_recipe$p + 6)
    } else {
      p <- backend$shared$p; d <- backend$upstream$d
      width <- 4 * (length(pre$uids) + length(pre$sids)) +
        4 * (p + 1) * (length(pre$beta_gids) + 1) + 4 * length(pre$gids) +
        4 * length(links) * (p + 3) + 8
      fund <- 16 * B * (width + p + 2 * d + 8)
    }
  } else if (identical(kind, "dr_pooled")) {
    # V33 common metadata starts at16384. Keep8192, all1024-per-record
    # slots, ID bytes and link bytes for the unchanged evaluator headers.
    fund <- 8192
  } else if (identical(kind, "dr_separate")) {
    # V38 common metadata:8192+64*records+64*cells, in addition to maps,
    # dedicated selection context2048, hits/upstream headers and full labels.
    # Keep2048+48*records+64*cells. V43's maps reservation is untouched.
    records <- length(pre$sids) + length(pre$uids) + length(pre$gids) +
      length(pre$tids) + length(links) + length(backend$selection_backend$blocks) + 1
    fund <- 6144 + 16 * records
  } else if (identical(kind, "reconstruction")) {
    if (reconstruction_separate(backend$shared_or)) {
      # The old16-fold common-row slab is separate from the 4*nl slice
      # funding score/forward caches. Keep twelve copies, and do not spend
      # either that link slice or any local-OR/forcing term a second time.
      nu <- if (inherits(backend, "cdid_augmented_backend")) length(pre$uids) else 0
      nz <- if (inherits(backend, "cdid_augmented_backend")) length(pre$zids) else 0
      fund <- 4 * B * (3 * length(st$cohorts) + length(backend$design_columns) +
        5 * length(pre$tids) * backend$treatment_recipe$p + nu + nz + 24)
    } else {
      # V44 maps include8192+64 times each of these record counts, plus
      # other global/incidence counts and full row/name/ID bytes. Keep2048
      # and at least48 slots per record; all cache/forward terms stay intact.
      # The unreplaced generic envelope has16384+1024-per-needed/group
      # slots and >=128-per-global-cell, so this smaller slice fits it too.
      records <- length(backend$cells) + length(backend$groups) + length(pre$cids) +
        length(pre$bids) + length(pre$tids) + length(links)
      fund <- 6144 + 16 * records
    }
  } else return(NULL)
  cp04_pair_read_operation(reader, rows, fund)
}
