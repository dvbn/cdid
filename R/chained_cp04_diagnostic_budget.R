# Static modeled-selection/reconstruction diagnostics only. Baseline and dated
# reference-layout callers keep their established public diagnostic contract.
cp04_field_diagnostics <- function(b) {
  is.null(b$reference_ref) && inherits(b, c("cdid_imputation_backend", "cdid_augmented_backend",
    "cdid_or_corollary1", "cdid_dr_corollary1"))
}

# Bound the fields actually returned, before creating any diagnostic table.
# E scans one retained field at a time; it never builds an output projection.
# Every alias occurrence is charged again. The factors below are occurrence
# counts, not a discount for physical R sharing. In particular OR/DR return
# complete links, supports, and fit records, while reconstruction returns only
# the listed link scalars and never the link adjoints.
cp04_public_diagnostic_plan <- function(b, selection_route = NULL, common_design = NULL,
                                        auxiliary_entries = 0) {
  if (!cp04_field_diagnostics(b)) {
    output <- 4 * b$accounting$owned_entries + 4096
    return(list(output = output, scratch = 4 * output, peak = 5 * output, work = 8 * output))
  }
  reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
  cp04_registry_field_memory_sources(reader, b)
  E <- function(x) unname(cp04_payload(x)["entries"])
  fields <- function(x, nms) {
    z <- 0
    for (nm in nms) z <- z + E(x[[nm]])
    z
  }
  # Covers the fixed literal text, field names, list/data-frame attributes and
  # row-name headers of both reconstruction and OR/DR schemas. Per-record
  # allowance covers generated scalar columns, role strings, numeric-to-string
  # names, paste separators, and outer/nested list headers. Arbitrarily long
  # source labels/formulas/keys/column names are scanned separately below.
  nc <- length(b$cells); ng <- length(b$groups); nl <- length(b$links)
  nt <- length(b$treatments); nr <- length(b$regressions); ns <- length(b$supports)
  sb <- b$selection_backend; so <- b$shared_or
  nb <- length(so$blocks); nq <- length(sb$blocks); nqc <- length(sb$cells)
  headers <- 65536 + 512 * (nc + ng + nl + nt + nr + ns + nb + nq + nqc)
  auxiliary_entries <- ipw_owner_int_bound(auxiliary_entries, "complete auxiliary diagnostic payload", min_val = 0)
  source <- E(selection_route) + E(common_design) + auxiliary_entries + 2 * fields(b,
    c("outcome_recipe", "treatment_recipe", "accounting", "selection_record", "design_source",
      "design_columns", "selection_columns", "est_method", "imputation_pooling", "outcome_pooling",
      "selection_model", "selection_pooling"))
  reconstruction <- inherits(b, c("cdid_imputation_backend", "cdid_augmented_backend"))
  if (reconstruction) {
    for (x in b$cells) source <- source + 2 * fields(x,
      c("name", "cohort", "pair", "L", "cell_id", "group_id", "q_id", "n_pair", "mass", "alpha", "a_k",
        "M_one", "MW_one", "MD_one", "Y_one", "Yaug_one", "clipping"))
    for (x in b$groups) source <- source + 3 * fields(x,
      c("group_id", "name", "pooling", "cell_ids", "cohorts", "p", "kept_columns", "beta", "rss", "objective")) +
        2 * fields(x$basis, c("aliased_columns", "nontrivial_candidates_retained"))
    # Both the named local intercept vector and coefficient alias are counted.
    source <- source + 3 * fields(so, c("p", "kept_columns", "beta", "intercept", "alpha", "rss", "objective")) +
      2 * fields(so$basis, c("aliased_columns"))
    for (x in so$blocks) source <- source + 3 * fields(x,
      c("block_id", "name", "pair", "L", "control_cohort_labels", "mass", "alpha"))
    for (x in b$treatments) source <- source + 3 * fields(x,
      c("trt_id", "target_cohort", "control_cohorts", "pair", "coefficients", "column_names", "n_risk_F", "clipping"))
    for (x in b$links) source <- source + 2 * fields(x,
      c("link_id", "target_cohort", "control_cohort_labels", "pair", "L", "target_cell_id", "or_block_id",
        "trt_id", "mu", "mu_ipw", "mu_reg", "mu_dr", "T_gl", "C_ol", "B_gl", "B_ol", "d_plus", "d_minus"))
    # Cohort labels may be numeric-to-character values in a group projection.
    reader <- get("cp04_reader", b$reader_ref, inherits = FALSE)
    source <- source + 32 * length(reader$state$cohorts) * (ng + 1)
    if (!is.null(sb)) {
      source <- source + 3 * fields(sb, c("selection_model", "selection_pooling", "coefficients", "design_spec"))
      for (x in sb$cells) source <- source + 2 * fields(x, c("type", "cohort", "label", "block_id"))
      for (x in sb$blocks) source <- source + 2 * fields(x,
        c("branch", "cell_ids", "B_m", "d_s", "has_intercept", "alpha", "beta", "H_chol", "m_b", "z_bar", "pi", "r", "diagnostics"))
    }
  } else {
    source <- source + 3 * fields(b$shared,
      c("p", "keep", "beta", "coefficients", "fit", "summary", "basis", "control_support_ids", "rss", "objective", "verification_scan"))
    for (x in b$regressions) source <- source + fields(x, setdiff(names(x), "qr"))
    source <- source + E(b$regression_specs) + E(b$selections) + E(b$selection_specs) +
      E(b$supports) + E(b$support_specs) + E(b$treatments) + E(b$treatment_specs) + 2 * E(b$links)
    # These vectors also occur outside the full treatment/regression records.
    for (x in b$treatments) source <- source + E(x$coefficients)
    for (x in b$regressions) source <- source + E(x$alpha)
    for (x in b$regression_specs) source <- source + fields(x, c("block_key", "key"))
  }
  # outcome_public_blocks returns the local fit/basis and repeats beta in
  # coefficients; OR/DR additionally repeat every basis in out$basis$blocks.
  object <- if (!is.null(b$outcome_groups)) b$outcome_groups else so$outcome_groups
  if (identical(b$outcome_pooling, "separate") && !is.null(object)) {
    for (j in seq_len(nrow(object$block_map))) {
      g <- object$groups[[object$block_map$group_id[j]]]
      source <- source + 3 * fields(g, c("keep", "beta", "basis", "fit")) +
        E(g$regressions[[object$block_map$local_block_id[j]]]$alpha)
    }
    headers <- headers + 512 * nrow(object$block_map)
  }
  output <- ipw_owner_safe_sum(c(headers, source), label = "source-field public diagnostic output")
  # At most four complete field envelopes coexist during construction: source
  # projections, partially built tables, diag_copy, and final concatenation / a
  # replaced pooled counterfactual record. This deliberately keeps all alias,
  # copy-on-modify, paste/vapply and full-length ifelse buffers in the charge.
  scratch <- ipw_owner_safe_sum(rep(output, 4L))
  list(output = output, scratch = scratch, peak = ipw_owner_safe_sum(c(output, scratch)), work = 8 * output,
    auxiliary_entries = auxiliary_entries)
}

# Complete union-of-auxiliary-descriptors output, including unsupported rows
# absent from b$links. The reporting source is chained_ipw_pipeline.R:
# build_ipw_auxiliary_reporting's desc_status/desc_reason loop and desc_table.
# No table or concatenated reason is built here. The empty seven-column
# schema funds its fixed headers/attributes. Seven vector slots plus at most
# one row-name slot per row (and two compressed-row-name slots) fund growth;
# each non-NA string then adds ceil(actual bytes / 2), exactly as cp04_payload.
cp04_public_auxiliary_diagnostic_bound <- function(registry) {
  aux <- registry$cp04_aux
  if (is.null(aux)) return(0)
  descriptors <- aux$descriptors; Q <- length(descriptors)
  if (length(registry$cp04_aux_index) != Q || length(registry$cp04_aux_errors) != Q)
    cp04_abort("Auxiliary diagnostic bound requires its complete registry descriptor sources.", "integrity")
  empty <- data.frame(index = integer(), g = numeric(), s = numeric(), t = numeric(),
    K = character(), status = character(), reason = character(), stringsAsFactors = FALSE)
  entries <- ipw_owner_safe_sum(c(unname(cp04_payload(empty)["entries"]), 8 * Q, 2))
  for (j in seq_len(Q)) {
    d <- descriptors[[j]]
    # vapply over the unnamed source descriptor list produces plain columns;
    # fmt_design_num is the same %.17g display formatter used by reporting.
    k_bytes <- ipw_owner_safe_sum(c(sum(nchar(fmt_design_num(d$K), type = "bytes")),
      max(0, length(d$K) - 1)))
    K_bytes <- 2 + k_bytes # braces in the returned K column
    status <- "supported"; reason_bytes <- 0
    if (is.na(registry$cp04_aux_index[j])) {
      err <- registry$cp04_aux_errors[[j]]
      if (!inherits(err, "cdid_ipw_support_error"))
        cp04_abort("Auxiliary diagnostic source lacks its captured support error.", "integrity")
      status <- "unsupported_pair"
      # These are every byte appended by reporting, including the complete
      # arbitrary-length condition message and complete formatted control set.
      reason_bytes <- ipw_owner_safe_sum(c(nchar(conditionMessage(err), type = "bytes"),
        nchar(" [g=, s=, t=, K={}]", type = "bytes"),
        nchar(fmt_design_num(d$g), type = "bytes"), nchar(fmt_design_num(d$s), type = "bytes"),
        nchar(fmt_design_num(d$t), type = "bytes"), k_bytes))
    }
    entries <- ipw_owner_safe_sum(c(entries, ceiling(K_bytes / 2),
      ceiling(nchar(status, type = "bytes") / 2), ceiling(reason_bytes / 2)))
  }
  entries
}

cp04_public_diagnostic_admit <- function(res, control, public_object, label) {
  b <- res$nuisance_registry$cp04
  auxiliary_entries <- cp04_public_auxiliary_diagnostic_bound(res$nuisance_registry)
  if (cp04_field_diagnostics(b) &&
      unname(cp04_payload(res$auxiliary_reporting$descriptor_status)["entries"]) > auxiliary_entries)
    cp04_abort("Complete auxiliary diagnostic table exceeds its registry source bound.", "bound")
  plan <- cp04_public_diagnostic_plan(b, res$selection, res$controls$or_basis, auxiliary_entries)
  retained <- cp04_public_diagnostic_retained(res, public_object)
  if (cp04_field_diagnostics(b)) {
    # The recorded envelope is a reservation, not a second copy of res. Take
    # the larger of that envelope and the actual retained result, counting all
    # repeated aliases in res. Its design owner lives behind an environment
    # reference and therefore needs its own explicit payload charge.
    design <- if (is.null(b$design_ref)) 0 else
      get("imputation_design_owner", b$design_ref, inherits = FALSE)$accounting$owned_entries
    actual <- ipw_owner_safe_sum(c(unname(cp04_payload(res)["entries"]), design))
    retained <- ipw_owner_safe_sum(c(max(retained, actual),
      unname(cp04_payload(public_object)["entries"])), label = "retained result and partial public diagnostic copy")
  }
  cp04_budget(control$owner_control, c(retained, plan$output), plan$scratch, plan$work, plan$output, label)
  plan$output
}

cp04_public_diagnostic_check <- function(out, b, bound) {
  if (cp04_field_diagnostics(b) && unname(cp04_payload(out)["entries"]) > bound)
    cp04_abort("Source-field public diagnostic output bound underestimated.", "bound")
  out
}

# Exact source decomposition for static reconstruction and conditional OR/DR.
# Reconstruction's base_entries adds pair payload and removes that SAME term
# from borrowed_entries; reader/design reservations and X remain untouched.
cp04_registry_field_memory_sources <- function(reader, b) {
  validate_cp04_reader(reader)
  if (inherits(b, "cdid_augmented_backend")) validate_augmented_backend(b) else
    if (inherits(b, "cdid_imputation_backend")) validate_imputation_backend(b) else
    if (inherits(b, "cdid_or_corollary1")) validate_or_corollary1_backend(b) else
    if (inherits(b, "cdid_dr_corollary1")) validate_dr_corollary1_backend(b) else
      cp04_abort("Registry field-prefix accounting requires a supported fitted backend.", "integrity")
  if (!identical(b$reader_ref, reader$ref) || !is.null(b$reference_ref) || !is.null(b$layout))
    cp04_abort("Registry field-prefix accounting requires its static reader and backend.", "integrity")
  reconstruction <- inherits(b, c("cdid_imputation_backend", "cdid_augmented_backend"))
  d <- if (is.null(b$design_ref)) NULL else get("imputation_design_owner", b$design_ref, inherits = FALSE)
  if (!is.null(d)) validate_imputation_design_owner(d, reader = reader)
  design_owned <- if (reconstruction) ipw_owner_int_bound(b$accounting$design_owned_entries,
    "backend design owned entries", min_val = 0) else 0
  design_borrowed <- if (reconstruction) ipw_owner_int_bound(b$accounting$design_borrowed_entries,
    "backend design borrowed entries", min_val = 0) else 0
  # X_frame owns the new design; design_owner borrows that owner's exact
  # reservation. A reader-state design has neither. No equal-array test enters.
  if (reconstruction && (if (identical(b$design_source, "X_frame"))
      is.null(d) || design_owned != d$accounting$owned_entries || design_borrowed != 0 else
    if (identical(b$design_source, "design_owner"))
      is.null(d) || design_owned != 0 || design_borrowed != d$accounting$owned_entries else
      !identical(b$design_source, "reader$state$Z") || !is.null(d) || design_owned != 0 || design_borrowed != 0))
    cp04_abort("Registry reconstruction design reservation differs from its sealed source.", "integrity")
  sb <- b$selection_backend
  if (!is.null(sb)) {
    validate_selection_models(sb)
    if (!identical(sb$reader_ref, reader$ref) || !is.null(sb$reference_ref))
      cp04_abort("Registry selection reservation belongs to a different source.", "integrity")
  }
  list(reader_owned_entries = reader$accounting$owned_entries,
    reader_borrowed_owner_entries = reader$accounting$borrowed_owner_entries,
    design_owned_entries = design_owned, design_borrowed_entries = design_borrowed,
    design_ref = b$design_ref, selection_ref = sb$ref)
}

cp04_registry_field_memory_accounting <- function(prefix, reader, b) {
  prefix <- ipw_owner_int_bound(prefix, "registry retained pipeline prefix", min_val = 0)
  src <- cp04_registry_field_memory_sources(reader, b)
  external <- ipw_owner_int_bound(b$accounting$retained_external_entries, "backend retained external prefix", min_val = 0)
  expected <- max(0, prefix - src$reader_borrowed_owner_entries)
  if (external != expected || b$accounting$borrowed_entries != ipw_owner_safe_sum(c(
      src$reader_owned_entries, src$reader_borrowed_owner_entries,
      src$design_owned_entries, src$design_borrowed_entries, external)))
    cp04_abort("Registry field-prefix reservation differs from its source components.", "integrity")
  if (external == 0) return(NULL)
  cp04_budget(b$control, c(b$accounting$owned_entries, b$accounting$borrowed_entries),
    2048, 0, 2048, "registry field-prefix metadata")
  rec <- c(list(kind = "cp04_registry_field_memory_prefix", backend_ref = b$ref,
    reader_ref = reader$ref, pipeline_prefix_entries = prefix, embedded_external_entries = external), src)
  if (unname(cp04_payload(rec)["entries"]) > 2048)
    cp04_abort("Registry field-prefix metadata reservation underestimated.", "bound")
  ipw_owner_seal(rec, "cp04_registry_field_memory_prefix")
}

cp04_registry_field_retained_entries <- function(registry, retained) {
  rec <- registry$cp04_memory_accounting; b <- registry$cp04; reader <- registry$cp04_reader
  if (!ipw_owner_seal_ok(rec, "cp04_registry_field_memory_prefix") ||
      !identical(rec$kind, "cp04_registry_field_memory_prefix") ||
      !identical(b$ref, rec$backend_ref) || !identical(reader$ref, rec$reader_ref))
    cp04_abort("Pipeline field-prefix composition requires its same sealed registry.", "integrity")
  src <- cp04_registry_field_memory_sources(reader, b)
  if (!identical(src, rec[names(src)]) ||
      b$accounting$retained_external_entries != rec$embedded_external_entries ||
      max(0, rec$pipeline_prefix_entries - src$reader_borrowed_owner_entries) != rec$embedded_external_entries ||
      b$accounting$borrowed_entries != ipw_owner_safe_sum(c(src$reader_owned_entries,
        src$reader_borrowed_owner_entries, src$design_owned_entries, src$design_borrowed_entries,
        rec$embedded_external_entries)) || retained < rec$pipeline_prefix_entries)
    cp04_abort("Pipeline field-prefix composition differs from its exact source reservation.", "integrity")
  retained - rec$embedded_external_entries
}
