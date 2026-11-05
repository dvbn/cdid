# Character roles of one raw-registry record bundle, not the whole catalog.
# Every alias occurrence is charged; equal strings imply no storage sharing.
ipw_raw_record_char_inventory <- function(metadata, frame, aux, control,
                                           retained = 0, envelope = control, work = 0) {
  M <- nrow(metadata$catalog$rows)
  descriptors <- if (is.null(aux)) list() else aux$descriptors
  Q <- length(descriptors); D <- ipw_owner_safe_sum(c(M, Q))
  p <- length(colnames(frame$Z_a)); C <- length(frame$a_cohorts)
  max_cohorts <- C
  for (d in descriptors) max_cohorts <- max(max_cohorts, 1 + length(d$K))
  # Only common name vectors and one auxiliary name/key are constructed.
  floor <- ipw_owner_safe_sum(c(4096, 64 * (p + max_cohorts + Q + 1)))
  ipw_owner_safe_sum(c(retained, floor), envelope$max_matrix_entries,
    "raw-IPW character inventory initial scratch")
  peak <- floor; scan_work <- 0
  scan <- function(x) {
    v <- ipw_fitted_view_char_bytes(x, control)
    peak <<- max(peak, floor + v$charged_slots)
    ipw_owner_safe_sum(c(retained, peak), envelope$max_matrix_entries,
      "raw-IPW character inventory source scratch")
    scan_work <<- ipw_owner_safe_sum(c(scan_work, v$work_ops_bound), envelope$max_work_ops,
      "raw-IPW character inventory scan work")
    ipw_owner_safe_sum(c(work, scan_work), envelope$max_work_ops,
      "raw-IPW prefix and character inventory work")
    v
  }
  labels <- c("g", if (max_cohorts > 1) paste0("k_", seq_len(max_cohorts - 1)))
  common <- scan(list(colnames(frame$Z_a), labels, as.character(frame$a_cohorts)))
  catalog <- scan(list(metadata$row_mapping$functional_key))
  op <- if (M) max(metadata$catalog$rows$operator_id) else 0
  generated <- scan(list(c(paste0("link:r", M), paste0("nuisance:op", op, ":r", M),
    paste0("reporting:", Q), paste0("reporting_nuisance:", Q), paste0("op", op, ":mask", D)),
    metadata$row_mapping$link_id, metadata$row_mapping$nuisance_id))
  # Numeric pair/support/model provenance has at most three exact count
  # components (<=16 decimal digits each) plus fixed prefixes.
  max_id <- max(128, generated$max_string_bytes)
  max_key <- catalog$max_string_bytes
  for (d in descriptors) {
    # Same numeric-to-character conversion as the reporting operator factory;
    # count every cohort plus separators before gsub removes spaces.
    cohort <- scan(list(as.character(d$g), as.character(d$K)))
    max_key <- max(max_key, nchar("reporting_operator:", type = "bytes") +
      cohort$source_bytes + max(1, length(d$K)))
    own <- scan(list(d$link_id, d$nuisance_id, d$functional_key))
    max_id <- max(max_id, own$max_string_bytes)
    max_key <- max(max_key, own$max_string_bytes)
  }
  # Source ceilings: nuisance/link records each carry <=8 cohort-name vectors;
  # fitted/point records each <=32 cohort-name and <=16 covariate-name vectors.
  # 128 full common vectors also cover compact support/model/provenance names.
  # Scalar own/defining-pair/view/point/backend/spec/list ID/key roles <=16 each.
  # The existing 20480 plus numeric cardinality terms still fund list/name
  # headers, static disclosures and numeric-index provenance; none are removed.
  name_copies <- 128; id_copies <- 16; key_copies <- 16
  bytes <- ipw_owner_safe_sum(c(ipw_owner_safe_product(name_copies, common$source_bytes),
    ipw_owner_safe_product(id_copies, max_id), ipw_owner_safe_product(key_copies, max_key)),
    label = "raw-IPW local record character bytes")
  # Same byte-to-slot/header convention as ipw_fitted_view_char_bytes, applied
  # once to this local role bundle without allocating repeated strings.
  slots <- ipw_owner_safe_sum(c(ceiling(bytes / 2), 32), label = "raw-IPW local record character slots")
  total <- ipw_owner_safe_product(D, slots, label = "raw-IPW D local character records")
  list(records = D, catalog_records = M, auxiliary_records = Q,
    common_name_bytes = common$source_bytes, max_id_bytes = max_id, max_key_bytes = max_key,
    name_vector_copies = name_copies, id_scalar_copies = id_copies, key_scalar_copies = key_copies,
    per_record_source_bytes = bytes, per_record_slots = slots, total_record_slots = total,
    metadata_entries = 2048, scratch_entries = peak, work_ops_bound = scan_work)
}
