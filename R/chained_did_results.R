# Private result adapters (CP02): already computed chained IPW-GMM pipeline
# outputs -> did 2.5.1 compatible result subclasses.
#
#   chained_did_MP()        reporting + explicit caller metadata -> cdid_MP
#                           (inherits "MP"; built with did::MP / did::DIDparams)
#   chained_did_AGGTEobj()  aggregation + cdid_MP -> cdid_AGGTEobj
#                           (inherits "AGGTEobj"; built with did::AGGTEobj)
#
# Nothing here estimates, aggregates, bootstraps or re-certifies anything. The
# adapters copy already computed quantities into the legacy fields that did's
# constructors store, validate that the copy is faithful, and keep every D5
# diagnostic (support reasons, certificates, separate band families, basis
# contrasts, aggregation target/exclusion records) in a `chained` slot that the
# narrow methods below summarize. The mapping is documented in
# docs/theory/cp02_did_result_adapters.md. Ordinary did objects are untouched:
# every method is registered only for the two subclasses and defers the
# numerical tabulation/plotting to did's own methods.
#
# Legacy-field conventions (field_map.json):
#   MP$V_analytical   = reporting$inference$Omega   (n * vcov; se = sqrt(diag/n))
#   MP$inffunc        = reporting$if_full with the caller's full-frame IDs as
#                       rownames, NA columns kept for unavailable targets
#   MP$n              = reporting$inference$n_total (full-frame units)
#   AGGTEobj$inf.function = complete aggregation$if_full / overall$if_full
#                       under the type-specific names of upstream compute.aggte

chained_result_abort <- function(message, class = "cdid_result_adapter_error") {
  cdid_abort(message, class)
}

chained_result_need <- function(ok, message, class = "cdid_result_adapter_error") {
  if (!isTRUE(ok)) chained_result_abort(message, class)
  invisible(TRUE)
}

chained_result_pick <- function(x, name, default = NULL) {
  if (is.list(x) && !is.null(x[[name]])) x[[name]] else default
}

chained_result_versions <- function() {
  ver <- function(pkg) tryCatch(as.character(utils::packageVersion(pkg)),
                                error = function(e) NA_character_)
  list(cdid = ver("cdid"), did = ver("did"),
       did_tested = unname(tested_upstream[["did"]]))
}

chained_result_estimator <- function(est_method, selection_est_method = NULL) {
  base <- if (identical(est_method, "reg")) "chained outcome regression/control-mean GMM"
  else if (identical(est_method, "dr")) "chained doubly robust GMM"
  else "chained IPW-GMM"
  # The default probability-weighted label is unchanged.
  if (identical(selection_est_method, "imputation")) paste(base, "(pure actual-change imputation)") else
    if (identical(selection_est_method, "dr")) paste(base, "(augmented observation correction)") else base
}

# tidy()/ggdid() show the estimator label for OR/DR and for every pure-imputation
# or augmented route; the default IPW presentation is unchanged.
chained_result_show_estimator <- function(x) {
  (is.character(x$DIDparams$est_method) && x$DIDparams$est_method %in% c("reg", "dr")) ||
    isTRUE(grepl("(pure actual-change imputation)", x$chained$estimator, fixed = TRUE)) ||
    isTRUE(grepl("(augmented observation correction)", x$chained$estimator, fixed = TRUE))
}

chained_result_references <- function(est_method = "ipw", selection_est_method = NULL) {
  list(
    estimator = paste(chained_result_estimator(est_method, selection_est_method), "difference-in-differences"),
    cdid = paste0("Bellego C., Benatia D. and V. Dortet-Bernadet (2025). ",
                  "\"The Chained Difference-in-Differences.\" Journal of ",
                  "Econometrics 248, 105783. doi:10.1016/j.jeconom.2024.105783."),
    did = paste0("Callaway B. and P. H. C. Sant'Anna (2021). ",
                 "\"Difference-in-Differences with Multiple Time Periods.\" ",
                 "Journal of Econometrics 225(2), 200-230. ",
                 "doi:10.1016/j.jeconom.2020.12.001.")
  )
}

# ------------------------------------------------------------------------------
# Shared validation of already computed inference blocks
# ------------------------------------------------------------------------------

# Effective inference is separate from requested DIDparams. In particular,
# biters remains a parameter for later aggregation, even if no draws ran here.
# Exact-zero and empty selections need no multiplier draws under the existing
# pipeline law (compute_chained_ipw_inference / ipw_pipeline_empty_inference).
#
# Released bstrap = FALSE, cband = TRUE endpoint convention (D6; flag_reference/):
# the only accepted disagreement between requested and effective flags, and only
# when the pipeline recorded that law itself (`inference$flag_law`):
#   endpoint "cells"      (att_gt reported cells, simple aggregation): effective
#                         analytic SE + pointwise normal interval, no draw;
#   endpoint "aggte_bins" (group/dynamic/calendar): effective analytic SE +
#                         multiplier sup-t critical value (draws for the band).
# Any other contradiction between caller metadata and the computed inference
# still refuses. The reported-SE kind follows the pipeline's effective
# `bootstrap` record, never the mere presence of multiplier draws; legacy
# blocks without that record keep the presence rule.
chained_result_released_false_true_law <- function(endpoint = c("cells", "aggte_bins", "simple")) {
  endpoint <- match.arg(endpoint)
  ipw_pipeline_flag_law(if (endpoint %in% c("cells", "simple")) "cells_false_true" else "aggte_bins_false_true")
}

chained_result_inference_flags <- function(inference, certificate, requested = list(),
                                           what = "reported cells", zero_law = FALSE,
                                           endpoint = c("cells", "aggte_bins", "simple")) {
  endpoint <- match.arg(endpoint)
  chained_result_need(is.list(inference),
    paste0("The ", what, " inference block must be the pipeline's inference list."))
  mult <- inference$multiplier
  has_draws <- !is.null(mult)
  law <- inference$flag_law
  released_law <- is.character(law) && length(law) == 1L && !is.na(law) &&
    identical(law, chained_result_released_false_true_law(endpoint))
  bootstrap <- if (is.logical(inference$bootstrap) && length(inference$bootstrap) == 1L && !is.na(inference$bootstrap))
    inference$bootstrap else has_draws
  chained_result_need(bootstrap || !has_draws || released_law,
    paste0("Inconsistent ", what, " inference: multiplier draws are present but the reported standard errors are ",
           "recorded as analytic without the released bstrap = FALSE, cband = TRUE endpoint law."))
  chained_result_need(!bootstrap || has_draws,
    paste0("Inconsistent ", what, " inference: bootstrap IQR standard errors are recorded without a multiplier result."))
  family <- chained_result_pick(certificate, "band_family")
  family_source <- "numerical_certificate"
  if (is.null(family) && !is.null(inference$band_family)) {
    # The pipeline records the effective family itself (bstrap and cband are
    # separate flags: bootstrap SEs may carry a pointwise critical value).
    family_source <- "inference_band_family"
    family <- inference$band_family
  }
  if (is.null(family)) {
    family_source <- "inference_multiplier_presence"
    family <- if (has_draws && (bootstrap || released_law)) "simultaneous_multiplier_sup_t_over_certified_ids" else
      if (identical(inference$status, "empty_supported_selection")) "none" else
      "pointwise_normal_over_certified_ids"
  }
  simultaneous <- startsWith(family, "simultaneous")
  chained_result_need(!simultaneous || has_draws,
    paste0("Inconsistent ", what, " inference: band family `", family,
           "` claims a simultaneous multiplier band without a multiplier bootstrap result."))
  if (!is.null(inference$simultaneous)) {
    chained_result_need(identical(isTRUE(inference$simultaneous), simultaneous),
      paste0("Inconsistent ", what, " inference: the recorded effective `simultaneous` flag disagrees with band family `", family, "`."))
  }
  if (released_law) {
    # The recorded endpoint law must match what was actually computed.
    chained_result_need(!bootstrap && identical(simultaneous, identical(endpoint, "aggte_bins") && has_draws),
      paste0("Inconsistent ", what, " inference: the recorded released bstrap = FALSE, cband = TRUE law does not ",
             "describe the computed inference (expected analytic SEs and ",
             if (identical(endpoint, "aggte_bins")) "a multiplier band whenever draws ran" else "no band", ")."))
  }
  req_b <- requested$bstrap
  req_c <- requested$cband
  released_request <- !is.null(req_b) && !is.null(req_c) && isFALSE(req_b) && isTRUE(req_c)
  if (!is.null(req_b)) {
    chained_result_need(is.logical(req_b) && length(req_b) == 1L && !is.na(req_b),
      "`bstrap` metadata must be TRUE or FALSE when supplied.")
    chained_result_need(zero_law || req_b == bootstrap, if (bootstrap)
      paste0("Caller metadata says bstrap = FALSE but the supplied ", what,
             " inference carries bootstrap IQR standard errors.") else
      paste0("Caller metadata says bstrap = TRUE but the supplied ", what,
             " inference is analytic; bootstrap standard errors cannot be claimed."))
  }
  if (!is.null(req_c)) {
    chained_result_need(is.logical(req_c) && length(req_c) == 1L && !is.na(req_c),
      "`cband` metadata must be TRUE or FALSE when supplied.")
    # Accepted disagreements: requested bstrap = FALSE, cband = TRUE with
    # the pipeline's own record of the released endpoint law, or simple
    # aggregation (scalar/no-band even when bstrap = TRUE, cband = TRUE).
    simple_pointwise <- identical(endpoint, "simple") && !simultaneous && isTRUE(req_b) && bootstrap && !released_law
    chained_result_need(
      zero_law ||
        (req_c == simultaneous && !(identical(endpoint, "simple") && simultaneous)) ||
        (released_request && released_law) ||
        simple_pointwise,
      if (simultaneous)
        paste0("Caller metadata says cband = ", req_c, " but the supplied ", what,
               " critical value is a simultaneous multiplier sup-t quantile.") else
        paste0("Caller metadata says cband = TRUE but the supplied ", what,
               " inference has no simultaneous band; a band cannot be claimed."))
  }
  n_draws <- if (has_draws) {
    nd <- chained_result_pick(chained_result_pick(mult, "metadata"), "n_draws",
                              chained_result_pick(mult, "n_draws", chained_result_pick(inference, "n_draws")))
    if (is.null(nd) || is.na(nd)) NA_integer_ else as.integer(nd)
  } else NA_integer_
  multiplier_summary <- if (has_draws) {
    keep <- setdiff(names(mult), c("raw_draws", "cluster_sums", "cluster_sums_centered"))
    mult[keep]
  } else NULL
  list(bstrap = bootstrap, cband = simultaneous, biters = n_draws,
       kind = if (bootstrap) "multiplier_bootstrap_iqr_se" else "analytic",
       band_family = family, band_family_source = family_source,
       clustered = isTRUE(inference$clustered),
       n_clusters = chained_result_pick(inference, "n_clusters", NA_integer_),
       n_zero_eligible_clusters = chained_result_pick(inference, "n_zero_eligible_clusters", NA_integer_),
       n_total = chained_result_pick(inference, "n_total", NA_integer_),
       rng = inference$rng, multiplier = multiplier_summary,
       requested = list(bstrap = req_b, cband = req_c, biters = requested$biters),
       zero_law = zero_law,
       flag_law = if (is.null(law)) NA_character_ else as.character(law),
       released_false_true = released_law, endpoint = endpoint, draws = has_draws)
}

# One-line description of the actual inference law for summaries.
chained_result_law_text <- function(inf) {
  if (!isTRUE(inf$released_false_true)) return(NULL)
  if (inf$endpoint %in% c("cells", "simple")) {
    return(paste0("Requested bstrap = FALSE, cband = TRUE: released did 2.5.1 att_gt/simple-aggregation convention; ",
                  "analytical standard errors and the normal pointwise interval, no bootstrap draw."))
  }
  if (isTRUE(inf$cband)) {
    paste0("Requested bstrap = FALSE, cband = TRUE: released did 2.5.1 aggte convention; analytical standard errors, ",
           "simultaneous critical value from ", inf$biters, " multiplier draws standardized by their own IQR scale ",
           "(not by the reported SEs).")
  } else {
    paste0("Requested bstrap = FALSE, cband = TRUE: released did 2.5.1 aggte convention, but no multiplier draw was ",
           "needed (exact-zero or empty selection); analytical standard errors, pointwise interval.")
  }
}

# Full-frame unit identity. `frame_ids` is the caller's explicit full-frame ID
# vector in the influence-function ROW order of the pipeline (row i of if_full
# belongs to frame_ids[i]); it is never inferred from data. Optional per-unit
# metadata `units` (id, cohort, weight, cluster, eligible) may arrive in any
# row order and is aligned by identity.
# D61 results (`allow_missing_weight`) keep V = w_t1 by original id with NA
# where a unit has no weight at the first calendar date.
chained_result_frame <- function(frame_ids, units, n, inference, allow_missing_weight = FALSE) {
  chained_result_need(is.atomic(frame_ids) && !is.null(frame_ids) &&
    length(frame_ids) == n && !anyNA(frame_ids) && !anyDuplicated(frame_ids),
    paste0("`frame_ids` must be a complete, unique full-frame ID vector of length n_total = ", n,
           " given in influence-function row order (never inferred from first-period data)."))
  frame <- data.frame(id = frame_ids, stringsAsFactors = FALSE)
  frame$cohort <- NA_real_
  frame$weight <- NA_real_
  frame$cluster <- NA
  frame$eligible <- NA
  if (!is.null(units)) {
    chained_result_need(is.data.frame(units) && "id" %in% names(units),
      "`units` must be a data.frame with an `id` column (any row order).")
    chained_result_need(!anyNA(units$id) && !anyDuplicated(units$id) &&
      length(units$id) == n && all(as.character(units$id) %in% as.character(frame_ids)),
      "`units$id` must contain exactly the full-frame IDs once each.")
    pos <- match(as.character(frame_ids), as.character(units$id))
    if (!is.null(units$cohort)) {
      coh <- as.double(units$cohort[pos])
      chained_result_need(!any(is.nan(coh) | coh == -Inf, na.rm = TRUE),
        "Frame cohorts must be finite dates, Inf (never treated) or NA (unknown frame-only cohort).")
      frame$cohort <- coh
    }
    if (!is.null(units$weight)) {
      w <- as.double(units$weight[pos])
      chained_result_need(if (isTRUE(allow_missing_weight)) all(is.na(w) | (is.finite(w) & w > 0)) else
        all(is.finite(w)) && all(w > 0),
        "Frame weights must be finite and positive for every full-frame unit.")
      frame$weight <- w
    }
    if (!is.null(units$cluster)) {
      cl <- units$cluster[pos]
      chained_result_need(!anyNA(cl), "Frame cluster labels must not be NA.")
      frame$cluster <- cl
    }
    if (!is.null(units$eligible)) {
      el <- units$eligible[pos]
      chained_result_need(is.logical(el) && !anyNA(el),
        "`units$eligible` must be a complete logical vector.")
      frame$eligible <- el
    }
  }
  has_cluster <- !all(is.na(frame$cluster))
  if (isTRUE(inference$clustered)) {
    chained_result_need(has_cluster,
      "The inference is clustered: supply the aligned full-frame cluster labels in `units$cluster`.")
    chained_result_need(length(unique(frame$cluster)) == inference$n_clusters,
      paste0("Supplied cluster labels define ", length(unique(frame$cluster)),
             " clusters but the inference reports n_clusters = ", inference$n_clusters, "."))
  } else {
    chained_result_need(!has_cluster,
      "Cluster labels were supplied but the inference is not clustered; the labels cannot describe it.")
  }
  if (!all(is.na(frame$eligible))) {
    cl <- if (has_cluster) frame$cluster else seq_len(n)
    zero_el <- sum(!unique(cl) %in% cl[frame$eligible])
    chained_result_need(zero_el == inference$n_zero_eligible_clusters,
      paste0("Eligibility flags imply ", zero_el, " zero-eligible cluster(s) but the inference reports ",
             inference$n_zero_eligible_clusters, "."))
  }
  frame
}

# ------------------------------------------------------------------------------
# cdid_MP: reported cells -> MP
# ------------------------------------------------------------------------------

chained_did_DIDparams <- function(params, call, n, frame, flags, versions) {
  chained_result_need(is.list(params), "`params` must be a named list of explicit caller metadata.")
  required <- c("yname", "tname", "idname", "gname", "control_group", "anticipation",
                "alp", "tlist", "glist", "est_method", "base_period")
  missing_p <- required[!required %in% names(params)]
  chained_result_need(!length(missing_p),
    paste0("`params` lacks explicit caller metadata: ", paste(missing_p, collapse = ", "), "."))
  chr1 <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  chained_result_need(chr1(params$yname) && chr1(params$tname) && chr1(params$idname) &&
    chr1(params$gname) && (chr1(params$est_method) || is.function(params$est_method)) && chr1(params$base_period),
    "`yname`, `tname`, `idname`, `gname` and `base_period` must be single strings, and `est_method` must be a single string or function.")
  chained_result_need(chr1(params$control_group) &&
    params$control_group %in% c("nevertreated", "notyettreated"),
    "`control_group` must be \"nevertreated\" or \"notyettreated\".")
  chained_result_need(is.numeric(params$anticipation) && length(params$anticipation) == 1L &&
    is.finite(params$anticipation), "`anticipation` must be a finite number.")
  chained_result_need(is.numeric(params$alp) && length(params$alp) == 1L &&
    is.finite(params$alp) && params$alp > 0 && params$alp < 1, "`alp` must lie in (0, 1).")
  tlist <- as.double(params$tlist)
  glist <- as.double(params$glist)
  chained_result_need(length(tlist) > 0L && all(is.finite(tlist)) && !is.unsorted(tlist, strictly = TRUE),
    "`tlist` must be the strictly increasing calendar of the panel.")
  chained_result_need(length(glist) > 0L && all(is.finite(glist)) && !is.unsorted(glist, strictly = TRUE),
    "`glist` must be the strictly increasing treated cohorts (finite dates).")
  chained_result_need(is.null(params$panel) || isTRUE(params$panel),
    "Chained results are panel results; `panel` must be TRUE when supplied.")
  clustervars <- params$clustervars
  if (flags$clustered) {
    chained_result_need(chr1(clustervars),
      "Clustered inference needs the cluster variable name in `params$clustervars`.")
  } else {
    chained_result_need(is.null(clustervars),
      "`clustervars` was supplied but the inference is not clustered.")
  }
  # Inherited computation/output flags are recorded as honoured (CP08):
  # print_details reported progress, pl/cores admitted multiplier workers and
  # compute_inffunc selected point estimates only. faster_mode = FALSE keeps
  # its meaning: did's data branch was not used (one chained engine).
  flag_or <- function(x, default) if (is.null(x)) default else x
  dp <- did::DIDparams(
    yname = params$yname, tname = params$tname, idname = params$idname, gname = params$gname,
    xformla = params$xformla, data = NULL, control_group = params$control_group,
    anticipation = params$anticipation, weightsname = params$weightsname, fix_weights = params$fix_weights,
    alp = params$alp, bstrap = flags$requested$bstrap, biters = flags$requested$biters, clustervars = clustervars,
    cband = flags$requested$cband, print_details = isTRUE(params$print_details), faster_mode = FALSE,
    pl = isTRUE(params$pl), cores = as.integer(flag_or(params$cores, 1L)),
    est_method = params$est_method, base_period = params$base_period, panel = TRUE,
    true_repeated_cross_sections = FALSE, n = n, nG = length(glist), nT = length(tlist),
    tlist = tlist, glist = glist, call = call)
  # Both count families glance()/nobs() may read, populated consistently.
  dp$id_count <- n
  dp$treated_groups_count <- length(glist)
  dp$time_periods_count <- length(tlist)
  dp$time_periods <- tlist
  dp$treated_groups <- glist
  dp$compute_inffunc <- !identical(params$compute_inffunc, FALSE)
  if (!is.null(params$extra_args)) {
    dp$extra_args <- params$extra_args
  }
  if (flags$clustered) {
    dp$cluster_vector <- frame$cluster
    dp$cluster_vector_var <- clustervars[1L]
  }
  dp$cdid <- list(
    mode = "chained", estimator = chained_result_estimator(params$est_method, params$selection_est_method), versions = versions,
    inference = flags$kind, band_family = flags$band_family,
    requested_bstrap = flags$requested$bstrap, requested_cband = flags$requested$cband,
    requested_biters = flags$requested$biters,
    actual_draws = flags$biters, zero_law = flags$zero_law,
    effective_bstrap = flags$bstrap, effective_cband = flags$cband, flag_law = flags$flag_law,
    data_retained = FALSE,
    note = paste0("DIDparams describes an already computed chained fit: faster_mode = FALSE only ",
                  "means did's faster-mode data branch was not used; no panel data is retained, so ",
                  "unchanged upstream aggregation of this object is not supported; use the ",
                  "reviewed chained aggregation mapped by chained_did_AGGTEobj()."))
  dp
}

chained_did_MP <- function(reporting, frame_ids, call, params, units = NULL,
                           support = NULL, basis = NULL, band_families = NULL,
                           pretest = NULL, covariance_handling = NULL, share_law = NULL) {
  chained_result_need(is.list(reporting) && is.data.frame(reporting$target_meta),
    "`reporting` must be the pipeline reporting list with a `target_meta` data.frame.")
  chained_result_need(is.list(params) && is.numeric(params$alp) && length(params$alp) == 1L &&
    is.finite(params$alp) && params$alp > 0 && params$alp < 1,
    "`params` must be a named list of explicit caller metadata with `alp` in (0, 1).")
  tm <- reporting$target_meta
  need_cols <- c("g", "time", "target_id", "role", "available", "support_status", "support_reason")
  chained_result_need(all(need_cols %in% names(tm)),
    paste0("`target_meta` must carry columns ", paste(need_cols, collapse = ", "), "."))
  J <- nrow(tm)
  ids <- as.character(tm$target_id)
  chained_result_need(J > 0L && !anyNA(ids) && !anyDuplicated(ids),
    "`target_meta` must list at least one uniquely identified reported target.")
  chained_result_need(is.numeric(tm$g) && is.numeric(tm$time) && all(is.finite(tm$g)) &&
    all(is.finite(tm$time)), "Reported cells need finite cohort/date pairs.")
  est <- reporting$estimates
  se <- reporting$se
  chained_result_need(is.numeric(est) && length(est) == J && identical(names(est), ids),
    "`reporting$estimates` must be named exactly by `target_meta$target_id` in the same order.")
  chained_result_need(is.numeric(se) && length(se) == J &&
    (is.null(names(se)) || identical(names(se), ids)),
    "`reporting$se` must align with `target_meta` (same length and order).")
  for (nm in c("z_stat", "p_value", "ci_lower", "ci_upper", "cband_lower", "cband_upper", "analytical_se")) {
    v <- reporting[[nm]]
    chained_result_need(is.numeric(v) && length(v) == J && (is.null(names(v)) || identical(names(v), ids)),
      paste0("`reporting$", nm, "` must align with `target_meta`."))
  }
  avail <- reporting$available
  chained_result_need(is.logical(avail) && length(avail) == J && !anyNA(avail) &&
    identical(unname(avail), unname(as.logical(tm$available))),
    "`reporting$available` must be a complete logical selector equal to `target_meta$available`.")
  inference <- reporting$inference
  chained_result_need(is.list(inference) && is.numeric(inference$n_total) && length(inference$n_total) == 1L,
    "`reporting$inference` must carry `n_total` (full-frame unit count).")
  n <- as.integer(inference$n_total)
  if (identical(params$compute_inffunc, FALSE)) {
    return(chained_did_MP_point_only(reporting, tm, ids, est, se, avail, inference, n, frame_ids, call, params,
      units = units, support = support, basis = basis, covariance_handling = covariance_handling, share_law = share_law))
  }
  IF <- reporting$if_full
  chained_result_need(is.matrix(IF) && is.numeric(IF) && nrow(IF) == n && ncol(IF) == J,
    paste0("`reporting$if_full` must be an n_total x J numeric matrix (", n, " x ", J, ")."))
  chained_result_need(is.null(colnames(IF)) || identical(colnames(IF), ids),
    "`reporting$if_full` columns must follow `target_meta` order.")
  chained_result_need(is.null(rownames(IF)) || identical(rownames(IF), as.character(frame_ids)),
    "`reporting$if_full` rownames disagree with `frame_ids`.")
  Omega <- inference$Omega
  chained_result_need(is.matrix(Omega) && is.numeric(Omega) && all(dim(Omega) == J) &&
    (is.null(dimnames(Omega)) || (identical(rownames(Omega), ids) && identical(colnames(Omega), ids))),
    "`reporting$inference$Omega` must be the J x J asymptotic covariance in `target_meta` order.")
  vc <- inference$vcov
  if (!is.null(vc)) {
    chained_result_need(isTRUE(all.equal(unname(Omega), unname(vc) * n)),
      "`reporting$inference$Omega` must equal n_total * `reporting$inference$vcov`.")
  }
  crit <- reporting$cband_crit_val
  chained_result_need(is.numeric(crit) && length(crit) == 1L,
    "`reporting$cband_crit_val` must be a single critical value (NA when nothing is supported).")
  # Selection-before-arithmetic contract: unavailable targets are NA everywhere,
  # never zero; available targets are finite with finite IF columns.
  if (any(avail)) {
    chained_result_need(all(is.finite(est[avail])) && all(is.finite(IF[, avail, drop = FALSE])),
      "Available targets must have finite estimates and finite influence-function columns.")
  }
  if (any(!avail)) {
    chained_result_need(all(is.na(est[!avail])) && all(is.na(se[!avail])) &&
      all(is.na(IF[, !avail, drop = FALSE])) && all(is.na(Omega[!avail, ])) && all(is.na(Omega[, !avail])),
      "Unavailable targets must stay NA in estimates, SEs, IF columns and covariance (never zero-filled).")
  }
  if (any(avail)) {
    scale_ok <- isTRUE(all.equal(unname(reporting$analytical_se[avail]),
                                 unname(sqrt(pmax(0, diag(Omega)[avail]) / n))))
    chained_result_need(scale_ok,
      "`analytical_se` must equal sqrt(diag(Omega)/n_total): Omega is not on did's V_analytical scale.")
  }
  certificate <- reporting$numerical_certificate
  defaults <- formals(did::DIDparams)
  request <- list(bstrap = chained_result_pick(params, "bstrap", defaults$bstrap),
                  cband = chained_result_pick(params, "cband", defaults$cband),
                  biters = chained_result_pick(params, "biters", defaults$biters))
  zero_law <- !any(avail) || (all(est[avail] == 0) && all(IF[, avail, drop = FALSE] == 0))
  flags <- chained_result_inference_flags(inference, certificate,
    requested = request, what = "reported-cell", zero_law = zero_law, endpoint = "cells")
  if (!flags$cband && any(avail) && !is.na(crit)) {
    deterministic <- all(est[avail] == 0) && all(IF[, avail, drop = FALSE] == 0)
    expected <- if (deterministic) 0 else stats::qnorm(1 - params$alp / 2)
    chained_result_need(isTRUE(all.equal(unname(crit), expected)),
      "The pointwise critical value does not match `params$alp`.")
  }
  frame <- chained_result_frame(frame_ids, units, n, inference,
    allow_missing_weight = identical(share_law, "first_calendar_masses"))

  # Optional separate post-contrast basis: kept apart from ATT(g,t) fields.
  if (!is.null(basis)) {
    chained_result_need(is.list(basis) && is.data.frame(basis$meta) && "basis_id" %in% names(basis$meta) &&
      is.numeric(basis$estimates) && identical(names(basis$estimates), as.character(basis$meta$basis_id)),
      "`basis` must be the pipeline basis result (meta$basis_id aligned with named estimates).")
    chained_result_need(!any(basis$meta$basis_id %in% ids),
      "Basis contrast IDs collide with reported cell IDs; contrasts are never cohort/date cells.")
    cross <- basis$cross_vcov_with_reported_cells
    if (!is.null(cross)) {
      chained_result_need(is.matrix(cross) && nrow(cross) == J && ncol(cross) == nrow(basis$meta) &&
        (is.null(rownames(cross)) || identical(rownames(cross), ids)),
        "`basis$cross_vcov_with_reported_cells` must be J x (number of basis contrasts) in `target_meta` order.")
    }
    bif <- basis$if_full
    if (!is.null(bif)) {
      chained_result_need(is.matrix(bif) && nrow(bif) == n,
        "`basis$if_full` must have one row per full-frame unit.")
      chained_result_need(is.null(rownames(bif)) || identical(rownames(bif), as.character(frame_ids)),
        "`basis$if_full` rownames disagree with `frame_ids`.")
      rownames(bif) <- as.character(frame_ids)
      basis$if_full <- bif
    }
  }
  if (!is.null(band_families)) {
    chained_result_need(is.list(band_families) && !is.null(names(band_families)) &&
      all(vapply(band_families, function(f) is.list(f) && !is.null(f$band_family), logical(1L))),
      "`band_families` must be the pipeline's named list of per-family {certified_target_ids, band_family}.")
  }
  # Optional unconditional pretrend statistic: carried only when supplied
  # explicitly with its description; absence is recorded as such. A record
  # computed by the public D21 binding (`chained_pretest_record`, status
  # "computed") is validated by the same rule and keeps its own status; a
  # record without a status is an externally supplied precomputed statistic.
  W <- NULL
  Wpval <- NULL
  pretest_record <- list(status = "not_computed",
    reason = paste0("no unconditional pretrend Wald statistic was supplied; the public chained att_gt() path binds the ",
                    "full declared nonreference pretrend null"))
  if (!is.null(pretest)) {
    chained_result_need(is.list(pretest) && is.numeric(pretest$W) && length(pretest$W) == 1L &&
      is.finite(pretest$W) && pretest$W >= 0 && is.numeric(pretest$Wpval) && length(pretest$Wpval) == 1L &&
      is.finite(pretest$Wpval) && pretest$Wpval >= 0 && pretest$Wpval <= 1 &&
      is.character(pretest$cells) && length(pretest$cells) > 0L && all(pretest$cells %in% ids) &&
      is.numeric(pretest$df) && length(pretest$df) == 1L && pretest$df >= 1 &&
      is.character(pretest$description) && length(pretest$description) == 1L,
      paste0("`pretest` must supply W, Wpval, the reported cell IDs used, df and a description of the ",
             "covariance/rank rule; nothing is inferred."))
    chained_result_need(is.null(pretest$status) || identical(pretest$status, "computed"),
      "A supplied `pretest` record may only carry status \"computed\"; unavailable tests carry no W/Wpval.")
    chained_result_need(all(tm$role[match(pretest$cells, ids)] == "pre"),
      "`pretest$cells` must be declared pretrend comparisons (role \"pre\"); reference zeros and post cells are never tested restrictions.")
    W <- pretest$W
    Wpval <- pretest$Wpval
    pretest_record <- if (is.null(pretest$status)) c(list(status = "supplied_precomputed"), pretest) else pretest
  }

  versions <- chained_result_versions()
  dp <- chained_did_DIDparams(params, call, n, frame, flags, versions)
  chained_result_need(all(tm$g %in% dp$glist) && all(tm$time %in% dp$tlist),
    "Reported cohorts/dates must belong to the explicit `glist`/`tlist`.")

  dimnames(IF) <- list(as.character(frame_ids), ids)
  dimnames(Omega) <- list(ids, ids)
  mp <- did::MP(group = as.double(tm$g), t = as.double(tm$time), att = unname(as.double(est)),
                V_analytical = Omega, se = unname(as.double(se)), c = unname(crit), inffunc = IF,
                n = n, W = W, Wpval = Wpval, aggte = NULL, alp = params$alp, DIDparams = dp)
  deterministic <- avail & est == 0 & colSums(abs(IF)) == 0
  deterministic[is.na(deterministic)] <- FALSE
  mp$chained <- list(
    estimator = chained_result_estimator(params$est_method, params$selection_est_method),
    references = chained_result_references(params$est_method, params$selection_est_method), versions = versions,
    target_id = ids, target_meta = tm, role = as.character(tm$role), available = unname(avail),
    deterministic_reference = unname(deterministic),
    support_status = as.character(tm$support_status), support_reason = as.character(tm$support_reason),
    support = reporting$support, support_ledger = support,
    z_stat = stats::setNames(unname(as.double(reporting$z_stat)), ids),
    p_value = stats::setNames(unname(as.double(reporting$p_value)), ids),
    analytical_se = stats::setNames(unname(as.double(reporting$analytical_se)), ids),
    ci_lower = stats::setNames(unname(as.double(reporting$ci_lower)), ids),
    ci_upper = stats::setNames(unname(as.double(reporting$ci_upper)), ids),
    cband_lower = stats::setNames(unname(as.double(reporting$cband_lower)), ids),
    cband_upper = stats::setNames(unname(as.double(reporting$cband_upper)), ids),
    base_period = reporting$base_period,
    numerical_certificate = certificate,
    band_family = flags$band_family, band_family_source = flags$band_family_source,
    band_families = band_families,
    basis = basis,
    covariance_handling = covariance_handling,
    inference = flags[c("kind", "bstrap", "cband", "biters", "clustered", "n_clusters",
                        "n_zero_eligible_clusters", "n_total", "rng", "multiplier", "requested", "zero_law",
                        "flag_law", "released_false_true", "endpoint", "draws")],
    frame = frame,
    pretest = pretest_record)
  class(mp) <- c("cdid_MP", "MP")
  mp
}

# compute_inffunc = FALSE (did 2.5.1 att_gt.R convention): the reported cells
# keep their point estimates and support records; inffunc = NULL,
# V_analytical = NULL, se = NA, c = qnorm(1 - alp / 2) and W = Wpval = NULL.
# Nothing is zero-filled, estimated or certified here; aggte() refuses.
chained_did_MP_point_only <- function(reporting, tm, ids, est, se, avail, inference, n, frame_ids, call, params,
                                      units = NULL, support = NULL, basis = NULL, covariance_handling = NULL,
                                      share_law = NULL) {
  J <- length(ids)
  chained_result_need(is.null(reporting$if_full) && is.null(inference$Omega) && is.null(inference$vcov),
    "A compute_inffunc = FALSE result must not carry influence functions or covariance matrices.")
  chained_result_need(all(is.na(se)) && identical(inference$flag_law, "compute_inffunc_false_point_estimates_only"),
    "A compute_inffunc = FALSE result carries no standard errors and records its point-only law.")
  if (any(avail)) chained_result_need(all(is.finite(est[avail])), "Available targets must have finite point estimates.")
  if (any(!avail)) chained_result_need(all(is.na(est[!avail])), "Unavailable targets must stay NA (never zero-filled).")
  crit <- stats::qnorm(1 - params$alp / 2)
  chained_result_need(isTRUE(all.equal(unname(reporting$cband_crit_val), crit)),
    "The point-only critical value must be the pointwise normal constant of `params$alp`.")
  flags <- list(bstrap = FALSE, cband = FALSE, biters = NA_integer_, kind = "not_computed", band_family = "none",
    band_family_source = "compute_inffunc_false", clustered = isTRUE(inference$clustered),
    n_clusters = inference$n_clusters, n_zero_eligible_clusters = inference$n_zero_eligible_clusters,
    n_total = inference$n_total, rng = NULL, multiplier = NULL,
    requested = list(bstrap = FALSE, cband = FALSE, biters = chained_result_pick(params, "biters", formals(did::DIDparams)$biters)),
    zero_law = FALSE, flag_law = inference$flag_law, released_false_true = FALSE, endpoint = "cells", draws = FALSE)
  frame <- chained_result_frame(frame_ids, units, n, inference,
    allow_missing_weight = identical(share_law, "first_calendar_masses"))
  versions <- chained_result_versions()
  dp <- chained_did_DIDparams(params, call, n, frame, flags, versions)
  chained_result_need(all(tm$g %in% dp$glist) && all(tm$time %in% dp$tlist),
    "Reported cohorts/dates must belong to the explicit `glist`/`tlist`.")
  mp <- did::MP(group = as.double(tm$g), t = as.double(tm$time), att = unname(as.double(est)),
                V_analytical = NULL, se = rep(NA_real_, J), c = crit, inffunc = NULL,
                n = n, W = NULL, Wpval = NULL, aggte = NULL, alp = params$alp, DIDparams = dp)
  na_v <- stats::setNames(rep(NA_real_, J), ids)
  mp$chained <- list(
    estimator = chained_result_estimator(params$est_method, params$selection_est_method),
    references = chained_result_references(params$est_method, params$selection_est_method), versions = versions,
    target_id = ids, target_meta = tm, role = as.character(tm$role), available = unname(avail),
    deterministic_reference = unname(avail & tm$role == "reference" & est == 0),
    support_status = as.character(tm$support_status), support_reason = as.character(tm$support_reason),
    support = reporting$support, support_ledger = support,
    z_stat = na_v, p_value = na_v, analytical_se = na_v, ci_lower = na_v, ci_upper = na_v,
    cband_lower = na_v, cband_upper = na_v, base_period = reporting$base_period,
    numerical_certificate = NULL, band_family = "none", band_family_source = flags$band_family_source,
    band_families = NULL, basis = basis, covariance_handling = covariance_handling,
    inference = flags[c("kind", "bstrap", "cband", "biters", "clustered", "n_clusters",
                        "n_zero_eligible_clusters", "n_total", "rng", "multiplier", "requested", "zero_law",
                        "flag_law", "released_false_true", "endpoint", "draws")],
    frame = frame,
    pretest = list(status = "not_computed", reason = "compute_inffunc = FALSE: point estimates only"))
  mp$chained$deterministic_reference[is.na(mp$chained$deterministic_reference)] <- FALSE
  class(mp) <- c("cdid_MP", "MP")
  mp
}

# Presentation grid for did's ggdid.MP, which rebuilds a group-major G x T
# rectangle positionally. Each actual (g, t) row is placed by its own pair;
# cells of the rectangle that were never reported are NA and labelled
# "unrequested", reported-but-unavailable cells "unavailable".
chained_mp_presentation_grid <- function(object) {
  chained_result_need(inherits(object, "cdid_MP"), "Expected a cdid_MP object.")
  key <- paste(object$group, object$t)
  chained_result_need(!anyDuplicated(key), "Reported (g, t) pairs must be unique for plotting.")
  groups <- sort(unique(object$group))
  times <- sort(unique(object$t))
  grid <- data.frame(group = rep(groups, each = length(times)),
                     time = rep(times, times = length(groups)))
  idx <- match(paste(grid$group, grid$time), key)
  ch <- object$chained
  grid$target_id <- ifelse(is.na(idx), NA_character_, ch$target_id[idx])
  grid$att <- object$att[idx]
  grid$se <- object$se[idx]
  grid$status <- ifelse(is.na(idx), "unrequested",
                        ifelse(ch$available[idx], "available", "unavailable"))
  grid$support_status <- ifelse(is.na(idx), NA_character_, ch$support_status[idx])
  grid$support_reason <- ifelse(is.na(idx), NA_character_, ch$support_reason[idx])
  grid$row <- idx
  grid
}

# ------------------------------------------------------------------------------
# cdid_AGGTEobj: reviewed D5 aggregation -> AGGTEobj
# ------------------------------------------------------------------------------

# `requested` (optional): the aggte-level effective request, as resolved by the
# request executor (`bstrap`, `biters`, `cband`, `alp`, `clustervars`, each with
# its `source` = "inherited" or "explicit"). Upstream compute.aggte overwrites
# the MP's DIDparams with these effective values (l. 145-149); the same is done
# here, while the actual computation stays in the separate chained record.
chained_did_AGGTEobj <- function(aggregation, mp, call, min_e = -Inf, max_e = Inf,
                                 balance_e = NULL, requested = NULL) {
  chained_result_need(inherits(aggregation, "cdid_ipw_aggregation_result"),
    "`aggregation` must be a pipeline `cdid_ipw_aggregation_result`.")
  chained_result_need(inherits(mp, "cdid_MP"), "`mp` must be the cdid_MP the aggregation was computed from.")
  type <- aggregation$type
  chained_result_need(is.character(type) && length(type) == 1L &&
    type %in% c("simple", "group", "dynamic", "calendar"), "Unknown aggregation type.")
  # Selectors: the aggregation's own selector record is the only truthful
  # source. A nondefault selector is recorded as applied only when the
  # aggregation carries a record with exactly these values (dynamic/simple/group)
  # or reports them as ignored with the released warning (calendar).
  chained_result_need(is.numeric(min_e) && length(min_e) == 1L && !is.na(min_e) &&
    is.numeric(max_e) && length(max_e) == 1L && !is.na(max_e) &&
    (is.null(balance_e) || (is.numeric(balance_e) && length(balance_e) == 1L && is.finite(balance_e))),
    "`min_e`/`max_e` must be numeric scalars and `balance_e` NULL or a finite number.")
  sel <- aggregation$selectors
  nondefault <- !is.null(balance_e) || is.finite(min_e) || is.finite(max_e)
  if (is.null(sel)) {
    chained_result_need(!nondefault,
      paste0("This aggregation result carries no selector record, so a nondefault `min_e`/`max_e`/`balance_e` ",
             "cannot be recorded as applied."))
    sel <- list(type = type, balance_e = NULL, min_e = -Inf, max_e = Inf, nondefault = FALSE, applied = FALSE, ignored = FALSE,
                warning = NULL, convention = "no selector record: upstream default selectors (no filtering)",
                selected_cell_ids = unlist(aggregation$cell_target_ids))
  } else {
    chained_result_need(is.list(sel) && identical(sel$type, type) &&
      identical(as.double(sel$min_e), as.double(min_e)) && identical(as.double(sel$max_e), as.double(max_e)) &&
      identical(is.null(sel$balance_e), is.null(balance_e)) &&
      (is.null(balance_e) || identical(as.double(sel$balance_e), as.double(balance_e))),
      paste0("The selectors recorded by the aggregation (min_e = ", sel$min_e, ", max_e = ", sel$max_e, ", balance_e = ",
             if (is.null(sel$balance_e)) "NULL" else sel$balance_e, ") differ from the values being bound; nothing is relabelled."))
  }
  n <- mp$n
  frame_ids <- as.character(mp$chained$frame$id)
  labels <- as.character(aggregation$cell_labels)
  keys <- aggregation$cell_keys
  M <- length(labels)
  chained_result_need(M > 0L && length(keys) == M && !anyDuplicated(labels),
    "Aggregation cell keys/labels must be aligned and unique.")
  est <- aggregation$estimates
  se <- aggregation$se
  chained_result_need(is.numeric(est) && length(est) == M && (is.null(names(est)) || identical(names(est), labels)),
    "`aggregation$estimates` must align with `cell_labels`.")
  chained_result_need(is.numeric(se) && length(se) == M && (is.null(names(se)) || identical(names(se), labels)),
    "`aggregation$se` must align with `cell_labels`.")
  IF <- aggregation$if_full
  chained_result_need(is.matrix(IF) && nrow(IF) == n && ncol(IF) == M &&
    (is.null(colnames(IF)) || identical(colnames(IF), labels)),
    paste0("`aggregation$if_full` must be the complete n_total x M aggregate influence function (", n, " x ", M, ")."))
  chained_result_need(is.null(rownames(IF)) || identical(rownames(IF), frame_ids),
    "`aggregation$if_full` rownames disagree with the full-frame IDs of `mp`.")
  cs <- aggregation$cell_status
  chained_result_need(is.data.frame(cs) && nrow(cs) == M && identical(as.character(cs$label), labels) &&
    is.logical(cs$available), "`aggregation$cell_status` must describe each cell in order.")
  cell_avail <- as.logical(cs$available)
  if (any(cell_avail)) {
    chained_result_need(all(is.finite(est[cell_avail])) && all(is.finite(IF[, cell_avail, drop = FALSE])),
      "Available aggregates must have finite estimates and finite complete IF columns.")
  }
  if (any(!cell_avail)) {
    chained_result_need(all(is.na(est[!cell_avail])) && all(is.na(IF[, !cell_avail, drop = FALSE])),
      "Unavailable aggregates must remain NA (never zero-filled).")
  }
  cell_targets <- aggregation$cell_target_ids
  chained_result_need(is.list(cell_targets) && length(cell_targets) == M &&
    all(unlist(cell_targets) %in% mp$chained$target_id),
    "Aggregation intended cells must be reported targets of `mp`.")
  egt <- suppressWarnings(as.double(keys))
  if (type != "simple") {
    chained_result_need(all(is.finite(egt)),
      "Aggregation keys must be numeric cohorts (group), dates (calendar) or event times (dynamic).")
  }
  overall <- aggregation$overall
  chained_result_need(is.list(overall) && is.numeric(overall$estimate) && length(overall$estimate) == 1L &&
    is.numeric(overall$se) && length(overall$se) == 1L && is.numeric(overall$if_full) &&
    length(overall$if_full) == n && is.character(overall$status),
    "`aggregation$overall` must carry estimate, se, status and the complete n_total overall IF.")
  chained_result_need(is.null(names(overall$if_full)) || identical(names(overall$if_full), frame_ids),
    "`aggregation$overall$if_full` names disagree with the full-frame IDs of `mp`.")
  if (identical(overall$status, "available")) {
    chained_result_need(is.finite(overall$estimate) && all(is.finite(overall$if_full)),
      "An available overall summary must have a finite estimate and IF.")
  } else {
    chained_result_need(is.na(overall$estimate) && is.na(overall$se),
      "An unavailable overall summary must be NA.")
  }
  inference <- aggregation$inference
  # Zero/empty law: no draws are needed when nothing stochastic is available.
  zero_law <- !any(cell_avail) || (all(est[cell_avail] == 0) && all(IF[, cell_avail, drop = FALSE] == 0))
  req_flags <- list()
  if (!is.null(requested)) {
    chained_result_need(is.list(requested) && all(c("bstrap", "cband", "biters", "alp") %in% names(requested)),
      "`requested` must list the effective aggte-level `bstrap`, `cband`, `biters` and `alp` (with sources).")
    req_flags <- list(bstrap = requested$bstrap, cband = requested$cband, biters = requested$biters)
  }
  flags <- chained_result_inference_flags(inference, aggregation$numerical_certificate,
    requested = req_flags, what = paste0(type, " aggregation"), zero_law = zero_law,
    endpoint = if (identical(type, "simple")) "simple" else "aggte_bins")
  overall_cert <- overall$numerical_certificate
  overall_family <- chained_result_pick(overall_cert, "band_family",
    if (!identical(overall$status, "available")) "none" else
      if (!isTRUE(overall_cert$certified) || identical(mp$chained$numerical_policy$decision, "D23")) "pointwise_normal_over_supported_ids" else "pointwise_normal_over_certified_ids")
  chained_result_need(!startsWith(overall_family, "simultaneous"),
    "The overall summary interval is pointwise by upstream convention; a simultaneous family cannot be mapped there.")
  crit <- aggregation$cband_crit_val
  chained_result_need(is.numeric(crit) && length(crit) == 1L, "`aggregation$cband_crit_val` must be a scalar.")
  if (!is.null(requested) && identical(overall$status, "available") && !zero_law) {
    ov_boot <- identical(chained_result_pick(overall, "inference_kind", "analytic"), "multiplier_bootstrap_iqr_se")
    chained_result_need(ov_boot == isTRUE(requested$bstrap),
      "The overall summary's SE kind disagrees with the effective `bstrap` request; bootstrap SEs cannot be claimed or hidden.")
  }

  dp <- mp$DIDparams
  # Keep upstream parameter meanings for later aggregation; actual inference
  # belongs to the separate chained record and presentation labels.
  if (!is.null(requested)) {
    chained_result_need(is.numeric(requested$alp) && length(requested$alp) == 1L && requested$alp > 0 && requested$alp < 1,
      "`requested$alp` must lie in (0, 1).")
    dp$bstrap <- requested$bstrap; dp$cband <- requested$cband; dp$biters <- requested$biters; dp$alp <- requested$alp
    if (!is.null(requested$clustervars)) dp$clustervars <- requested$clustervars
    dp$cdid$requested_bstrap <- requested$bstrap; dp$cdid$requested_cband <- requested$cband
    dp$cdid$requested_biters <- requested$biters
    dp$cdid$aggte_request <- requested
  }
  dp$cdid$inference <- flags$kind
  dp$cdid$band_family <- flags$band_family
  dp$cdid$aggregation_type <- type
  dp$cdid$actual_draws <- flags$biters
  dp$cdid$zero_law <- zero_law
  dp$cdid$effective_bstrap <- flags$bstrap; dp$cdid$effective_cband <- flags$cband; dp$cdid$flag_law <- flags$flag_law

  bins <- IF
  dimnames(bins) <- list(frame_ids, labels)
  overall_if <- stats::setNames(as.double(overall$if_full), frame_ids)
  if (type == "simple") {
    inf_fun <- list(simple.att = stats::setNames(as.double(IF[, 1L]), frame_ids))
    overall_att <- unname(as.double(est[1L]))
    overall_se <- unname(as.double(se[1L]))
  } else {
    inf_fun <- switch(type,
      group = list(selective.inf.func.g = bins, selective.inf.func = overall_if),
      dynamic = list(dynamic.inf.func.e = bins, dynamic.inf.func = overall_if),
      calendar = list(calendar.inf.func.t = bins, calendar.inf.func = overall_if))
    overall_att <- unname(as.double(overall$estimate))
    overall_se <- unname(as.double(overall$se))
  }
  obj <- did::AGGTEobj(
    overall.att = overall_att, overall.se = overall_se, type = type,
    egt = if (type == "simple") NULL else egt,
    att.egt = if (type == "simple") NULL else unname(as.double(est)),
    se.egt = if (type == "simple") NULL else unname(as.double(se)),
    crit.val.egt = if (type == "simple") NULL else unname(crit),
    inf.function = inf_fun,
    # Upstream records the selectors on dynamic results only (compute.aggte 564-566).
    min_e = if (type == "dynamic") as.double(min_e) else NULL,
    max_e = if (type == "dynamic") as.double(max_e) else NULL,
    balance_e = if (type == "dynamic") balance_e else NULL,
    call = call, DIDparams = dp)
  col_or <- function(nm, default) if (!is.null(cs[[nm]])) cs[[nm]] else rep(default, length.out = M)
  cell <- data.frame(
    label = labels, key = as.character(keys), available = cell_avail,
    restricted = as.logical(col_or("restricted", FALSE)),
    share_dependent = as.logical(col_or("share_dependent", NA)),
    n_intended = as.integer(col_or("n_intended", lengths(cell_targets))),
    n_used = as.integer(col_or("n_used", NA_integer_)),
    intended_ids = as.character(col_or("intended_ids", vapply(cell_targets, paste, character(1L), collapse = "; "))),
    used_ids = as.character(col_or("used_ids", NA_character_)),
    excluded_ids = as.character(col_or("excluded_ids", NA_character_)),
    reason = as.character(col_or("reason", NA_character_)),
    stringsAsFactors = FALSE)
  overall_record <- list(
    status = overall$status, reason = chained_result_pick(overall, "reason", NA_character_),
    intended_bins = chained_result_pick(overall, "intended_bins", character(0L)),
    used_bins = chained_result_pick(overall, "used_bins", character(0L)),
    excluded_bins = chained_result_pick(overall, "excluded_bins", character(0L)),
    weights = chained_result_pick(overall, "weights"),
    share_dependent = chained_result_pick(overall, "share_dependent", NA),
    share_denominator = chained_result_pick(overall, "share_denominator", NA_real_),
    band_family = overall_family, numerical_certificate = overall_cert,
    inference_kind = chained_result_pick(overall, "inference_kind", "analytic"),
    n_draws = chained_result_pick(overall, "n_draws", NA_integer_),
    draws_reused_from = chained_result_pick(overall, "draws_reused_from"),
    ci_lower = chained_result_pick(overall, "ci_lower", NA_real_),
    ci_upper = chained_result_pick(overall, "ci_upper", NA_real_))
  z <- chained_result_pick(inference, "z_stat", stats::setNames(rep(NA_real_, M), labels))
  p <- chained_result_pick(inference, "p_value", stats::setNames(rep(NA_real_, M), labels))
  obj$chained <- list(
    estimator = mp$chained$estimator, references = mp$chained$references,
    versions = mp$chained$versions, pt = mp$chained$pt, type = type, status = aggregation$status,
    na_rm = isTRUE(aggregation$na_rm), reason = chained_result_pick(aggregation, "reason", NA_character_),
    target_definition = chained_result_pick(aggregation, "target_definition", NA_character_),
    cell = cell, cell_target_ids = cell_targets,
    cell_used_target_ids = chained_result_pick(aggregation, "cell_used_target_ids"),
    overall = overall_record,
    shares_known = isTRUE(aggregation$shares_known), cohort_shares = aggregation$cohort_shares,
    weights_theta = aggregation$weights_theta,
    z_stat = stats::setNames(unname(as.double(z)), labels),
    p_value = stats::setNames(unname(as.double(p)), labels),
    ci_lower = aggregation$ci_lower, ci_upper = aggregation$ci_upper,
    cband_lower = aggregation$cband_lower, cband_upper = aggregation$cband_upper,
    se_fixed_shares_only = aggregation$se_fixed_shares_only,
    if_fixed_full = aggregation$if_fixed_full, if_share_full = aggregation$if_share_full,
    decomposition = aggregation$decomposition,
    numerical_certificate = aggregation$numerical_certificate,
    band_family = flags$band_family, band_family_source = flags$band_family_source,
    inference = flags[c("kind", "bstrap", "cband", "biters", "clustered", "n_clusters",
                        "n_zero_eligible_clusters", "n_total", "rng", "multiplier", "requested", "zero_law",
                        "flag_law", "released_false_true", "endpoint", "draws")],
    requested = requested,
    frame_ids = frame_ids, mp_target_id = mp$chained$target_id,
    selectors = c(sel[setdiff(names(sel), "note")], list(note = chained_result_selector_note(sel))))
  if (!is.null(aggregation$share_law)) {
    # D61: masses H(I_g V) with V = w_t1 by original id (degree-0 share maps).
    obj$chained$share_law <- aggregation$share_law
    obj$chained$cohort_mass_complete <- aggregation$cohort_mass_complete
  }
  class(obj) <- c("cdid_AGGTEobj", "AGGTEobj")
  obj
}

chained_result_selector_note <- function(sel) {
  fmt <- function(x) if (is.null(x)) "NULL" else as.character(x)
  vals <- paste0("min_e = ", fmt(sel$min_e), ", max_e = ", fmt(sel$max_e), ", balance_e = ", fmt(sel$balance_e))
  if (isTRUE(sel$ignored)) return(paste0(vals, " ignored for calendar aggregation (released did 2.5.1 warning)"))
  if (!isTRUE(sel$nondefault)) return(paste0("upstream default selectors (", vals, "): no horizon/cohort filtering"))
  paste0(vals, " applied with the released did 2.5.1 convention (", sel$convention, ")")
}

# ------------------------------------------------------------------------------
# Narrow presentation methods (registered only for the subclasses)
# ------------------------------------------------------------------------------

# The heading records versions and each paper reference once. did's own
# citation is removed below; its call and estimate tables are retained.
# Detailed numerical/resource records remain in the object.
chained_result_banner <- function(object) {
  ch <- object$chained
  v <- ch$versions
  fmt_v <- function(x) if (length(x) != 1L || is.na(x)) "version unavailable" else x
  aggregate <- inherits(object, "cdid_AGGTEobj")
  kind <- if (aggregate) "cdid_AGGTEobj" else "cdid_MP"
  upstream <- if (aggregate) "AGGTEobj" else "MP"
  cat("\ncdid chained difference-in-differences (", ch$estimator, "): ", kind, "\n", sep = "")
  cat("  cdid ", fmt_v(v$cdid), " result built on did ", fmt_v(v$did), " (tested: did ",
      fmt_v(v$did_tested), "); did's ", upstream, " methods tabulate the estimates below.\n", sep = "")
  cat("  Chained estimators: ", ch$references$cdid, "\n", sep = "")
  cat("  Ordinary estimators: ", ch$references$did, "\n", sep = "")
}

# did's own summary of the presentation copy, captured line by line. did
# prints the Callaway and Sant'Anna citation after the call; that line and the
# blank lines after it are dropped, and every other line did prints is kept.
chained_result_upstream_summary <- function(object, class, ...) {
  method <- utils::getS3method("summary", class, envir = asNamespace("did"))
  chained_result_drop_citation(utils::capture.output(invisible(method(object, ...))))
}

chained_result_drop_citation <- function(lines) {
  i <- grep("^Reference: ", lines)
  if (length(i) != 1L) return(lines)
  j <- i
  while (j < length(lines) && !nzchar(trimws(lines[j + 1L]))) j <- j + 1L
  # Keep one blank line between the call and the table.
  if (i > 1L && !nzchar(trimws(lines[i - 1L]))) lines[-(i:j)] else lines[-(i:max(i, j - 1L))]
}

chained_result_cat_lines <- function(lines) {
  if (length(lines)) cat(paste0(lines, "\n"), sep = "")
  invisible(NULL)
}

# Unit table counts: all units, and when stored, those observed at every date
# and those with at least two observed outcomes (the units that form changes).
chained_result_units_text <- function(ch, n) {
  fr <- ch$frame
  parts <- paste0(n, " in the unit table")
  n_obs <- fr$n_obs
  n_dates <- length(ch$calendar)
  if (is.numeric(n_obs) && length(n_obs) == n && !anyNA(n_obs) && n_dates > 0L) {
    parts <- c(parts, paste0(sum(n_obs >= n_dates), " observed at all ", n_dates, " dates"),
               paste0(sum(n_obs >= 2L), " with at least two observed outcomes"))
  } else if (is.logical(fr$eligible) && length(fr$eligible) == n && !anyNA(fr$eligible)) {
    parts <- c(parts, paste0(sum(fr$eligible), " with at least two observed outcomes"))
  }
  cat("  Units: ", paste(parts, collapse = ", "), ".\n", sep = "")
}

# Handling of missing outcomes, read from the stored selection record.
chained_result_selection_text <- function(sel) {
  if (!is.list(sel) || !is.character(sel$method) || length(sel$method) != 1L) return(invisible(NULL))
  q <- function(x) paste0("\"", x, "\"")
  assumption <- if (is.character(sel$requested_selection) && length(sel$requested_selection) == 1L)
    sel$requested_selection else sel$canonical_assumption
  universe <- if (is.character(sel$fitting_universe) && length(sel$fitting_universe) == 1L &&
                  !identical(sel$fitting_universe, "F")) paste0(", selection_universe = ", q(sel$fitting_universe)) else ""
  text <- switch(sel$method,
    none = "none (selection = \"none\")",
    corollary1 = paste0("observation weighting under selection = ", q(assumption), " (selection_model = ",
      q(sel$selection_model), ", selection_pooling = ", q(sel$selection_pooling), universe, ")"),
    imputation = paste0("change imputation under selection = ", q(assumption), " (imputation_pooling = ",
      q(sel$imputation_pooling), "); no observation probability model is fitted"),
    augmented = paste0("augmented correction under selection = ", q(assumption), " (imputation_pooling = ",
      q(sel$imputation_pooling), "; selection_model = ", q(sel$selection_model), ", selection_pooling = ",
      q(sel$selection_pooling), "; selection formula ",
      if (isTRUE(sel$selection_formula_inherited)) "inherited from xformla" else "from selection_xformla", ")"),
    sel$method)
  cat("  Observation correction: ", text, ".\n", sep = "")
}

# GMM weighting actually used and whether covariance regularization was applied.
chained_result_gmm_text <- function(cov_h) {
  if (!is.list(cov_h) || !is.character(cov_h$weighting_canonical)) return(invisible(NULL))
  unweighted <- identical(cov_h$weighting_canonical, "identity")
  weighting <- if (unweighted) "identity" else
    if (identical(cov_h$weighting_canonical, "optimal")) "optimal (two-step)" else cov_h$weighting_canonical
  reason <- chained_result_trigger_text(cov_h$trigger_reason)
  lambda <- cov_h$lambda
  regularization <- if (unweighted) "not applicable" else if (isTRUE(cov_h$handler_invoked)) {
    paste0(if (identical(cov_h$handler_used, "tikhonov")) "Tikhonov ridge" else cov_h$handler_used, " applied",
           if (is.numeric(lambda) && length(lambda) == 1L && is.finite(lambda)) paste0(" with lambda = ", format(lambda, digits = 3)) else "",
           " (", reason, ")")
  } else paste0("not applied (", reason, ")")
  cat("  GMM weighting: ", weighting, "; covariance regularization: ", regularization, ".\n", sep = "")
}

chained_result_trigger_text <- function(code) {
  if (!is.character(code) || length(code) != 1L || is.na(code)) return("reason not recorded")
  switch(code,
    well_conditioned_pd = "moment covariance well conditioned",
    spectral_condition_below_trigger = "moment covariance near singular",
    exact_zero_eigenvalue = "moment covariance singular",
    roundoff_negative_psd = "moment covariance singular up to rounding",
    zero_covariance = "zero moment covariance",
    known_unit_covariance_rank_bound = "moment covariance singular: more moments than contributing units",
    cholesky_failure_on_valid_psd = "Cholesky factorization of the moment covariance failed",
    gsub("_", " ", code, fixed = TRUE))
}

# Inference actually performed (never the request alone).
chained_result_inference_text <- function(inf, what) {
  if (identical(inf$kind, "not_computed")) {
    return(paste0("  Inference: not computed (compute_inffunc = FALSE: point estimates only; no influence functions, ",
                  "standard errors, bands, pre-test or aggte()).\n"))
  }
  draws <- if (isTRUE(inf$draws) && length(inf$biters) == 1L && !is.na(inf$biters)) inf$biters else NA
  se <- if (isTRUE(inf$bstrap)) {
    paste0("multiplier-bootstrap IQR standard errors", if (!is.na(draws)) paste0(" (", draws, " draws)") else "")
  } else "analytic standard errors"
  cluster <- if (isTRUE(inf$clustered)) paste0(", clustered on ", inf$n_clusters, " clusters (",
                                              inf$n_zero_eligible_clusters, " without eligible units)") else ", unit-level"
  interval <- if (isTRUE(inf$cband)) {
    paste0("simultaneous band over the available ", what,
           if (!isTRUE(inf$bstrap) && !is.na(draws)) paste0(" (", draws, " multiplier draws)") else "")
  } else "pointwise intervals"
  paste0("  Inference: ", se, cluster, "; ", interval, ".\n")
}

# Unconditional pretrend Wald test, or the reason it was not computed.
chained_result_pretest_text <- function(pt) {
  if (!is.list(pt) || is.null(pt$status)) return(invisible(NULL))
  if (identical(pt$status, "supplied_precomputed")) {
    cat("  Parallel-trends pre-test: supplied precomputed statistic over ", length(pt$cells),
        " cell(s): ", pt$description, "\n", sep = "")
  } else if (identical(pt$status, "computed")) {
    n_ref <- length(pt$structural_reference_cells)
    cat("  Parallel-trends pre-test: W = ", format(pt$W, digits = 6), ", df = ", pt$df, ", p = ",
        format(pt$Wpval, digits = 4), " (Wald test of ",
        if (identical(as.integer(pt$df), 1L)) "the pretrend cell" else paste0("all ", pt$df, " pretrend cells jointly"),
        if (n_ref) paste0("; ", n_ref, " exact reference zero", if (n_ref == 1L) "" else "s", " excluded") else "",
        ").\n", sep = "")
  } else if (identical(pt$status, "unavailable")) {
    cat("  Parallel-trends pre-test: unavailable; ", pt$reason, ".\n", sep = "")
  } else if (identical(pt$status, "not_applicable")) {
    cat("  Parallel-trends pre-test: not applicable; ", pt$reason, ".\n", sep = "")
  } else {
    cat("  Parallel-trends pre-test: not computed; ", pt$reason, ".\n", sep = "")
  }
  invisible(NULL)
}

# Upstream's table label reads DIDparams$bstrap first ("Pointwise" whenever
# bstrap is FALSE, even when upstream itself bootstrapped the band). The local
# presentation copy therefore labels by the EFFECTIVE interval law: "Simult."
# exactly when the critical value is a multiplier sup-t quantile.
chained_result_presentation_params <- function(object) {
  inf <- object$chained$inference
  dp <- object$DIDparams
  dp$bstrap <- isTRUE(inf$bstrap) || isTRUE(inf$cband)
  dp$cband <- isTRUE(inf$cband)
  dp
}

summary.cdid_MP <- function(object, ...) {
  chained_result_banner(object)
  # Upstream's table labels read DIDparams; use effective flags only on this
  # local presentation copy, preserving the returned object's request fields.
  shown <- object
  shown$DIDparams <- chained_result_presentation_params(object)
  chained_result_cat_lines(chained_result_upstream_summary(shown, "MP", ...))
  ch <- object$chained
  n_av <- sum(ch$available)
  cat("\nChained support: ", n_av, " of ", length(ch$available), " requested cells available",
      if (n_av < length(ch$available)) " (unavailable cells are NA, not zero effects)" else "",
      ".\n", sep = "")
  if (n_av < length(ch$available)) {
    un <- which(!ch$available)
    reason <- ch$support_reason[un]
    reason <- ifelse(!is.na(reason) & nzchar(trimws(reason)), paste0(": ", reason), "")
    cat("  Unavailable: ", paste0(ch$target_id[un], " [", ch$support_status[un], "]", reason, collapse = "; "), "\n", sep = "")
  }
  if (any(ch$deterministic_reference)) {
    cat("  Exact-zero reference cells (no test statistic): ", sum(ch$deterministic_reference), ".\n", sep = "")
  }
  chained_result_units_text(ch, object$n)
  chained_result_selection_text(ch$selection)
  chained_result_pt_text(ch$pt, object$DIDparams$base_period)
  chained_result_gmm_text(ch$covariance_handling)
  inf <- ch$inference
  cat(chained_result_inference_text(inf, "cells"))
  law_text <- chained_result_law_text(inf)
  if (!is.null(law_text)) cat("  ", law_text, "\n", sep = "")
  chained_result_pretest_text(ch$pretest)
  if (!is.null(ch$basis)) {
    cat("  Post-contrast basis: ", nrow(ch$basis$meta), " coordinate(s) retained separately with their own band family;",
        " they are not ATT(g,t) cells and are not tabulated above.\n", sep = "")
  }
  invisible(object)
}

print.cdid_MP <- function(x, ...) {
  summary(x, ...)
}

tidy.cdid_MP <- function(x, ...) {
  out <- NextMethod()
  ch <- x$chained
  # Exact-zero reference cells: the pipeline convention is NA, not 0/0 = NaN.
  out$statistic <- unname(ch$z_stat)
  out$p.value <- unname(ch$p_value)
  cert_ids <- chained_result_pick(ch$numerical_certificate, "certified_target_ids")
  out$target.id <- ch$target_id
  out$role <- ch$role
  out$available <- ch$available
  out$support.status <- ch$support_status
  out$support.reason <- ch$support_reason
  out$certified <- if (is.null(cert_ids)) rep(NA, length(ch$target_id)) else ch$target_id %in% cert_ids
  if (chained_result_show_estimator(x)) out$estimator <- ch$estimator
  out$band.family <- rep(ch$band_family, length(ch$target_id))
  out
}

glance.cdid_MP <- function(x, ...) {
  out <- NextMethod()
  ch <- x$chained
  out$estimator <- ch$estimator
  out$n.requested <- length(ch$available)
  out$n.available <- sum(ch$available)
  out$band.family <- ch$band_family
  out$inference <- ch$inference$kind
  out$clustered <- ch$inference$clustered
  out
}

ggdid.cdid_MP <- function(object, ylim = NULL, xlab = NULL, ylab = NULL, title = "Group",
                          xgap = 1, ncol = 1, legend = TRUE, group = NULL, ref_line = 0,
                          theming = TRUE, grtitle = "Group", ...) {
  grid <- chained_mp_presentation_grid(object)
  plain <- did::MP(group = grid$group, t = grid$time, att = grid$att, V_analytical = NULL,
                   se = grid$se, c = object$c, inffunc = NULL, n = object$n, alp = object$alp,
                   DIDparams = object$DIDparams)
  p <- did::ggdid(plain, ylim = ylim, xlab = xlab, ylab = ylab, title = title, xgap = xgap,
                  ncol = ncol, legend = legend, group = group, ref_line = ref_line,
                  theming = theming, grtitle = grtitle, ...)
  if (chained_result_show_estimator(object))
    p <- p + ggplot2::labs(caption = object$chained$estimator)
  attr(p, "cdid_presentation_grid") <- grid
  p
}

summary.cdid_AGGTEobj <- function(object, ...) {
  chained_result_banner(object)
  shown <- object
  shown$DIDparams <- chained_result_presentation_params(object)
  chained_result_cat_lines(chained_result_upstream_summary(shown, "AGGTEobj", ...))
  ch <- object$chained
  cell <- ch$cell
  cat("\nChained ", ch$type, " aggregation: ", sum(cell$available), " of ", nrow(cell),
      " intended aggregate(s) available; na.rm = ", ch$na_rm, ".\n", sep = "")
  # Restricted and unavailable targets are always disclosed with their reasons.
  for (i in which(!cell$available | cell$restricted)) {
    cat("  ", cell$label[i], if (cell$available[i]) " (restricted target)" else " (unavailable)",
        ": ", cell$reason[i], "\n", sep = "")
  }
  ov <- ch$overall
  if (ch$type != "simple") {
    if (identical(ov$status, "available")) {
      n_int <- length(ov$intended_bins)
      bins <- if (length(ov$excluded_bins)) {
        paste0(length(ov$used_bins), " of its ", n_int, " intended bins (excluded: ",
               paste(ov$excluded_bins, collapse = ", "), ")")
      } else paste0("its ", n_int, " intended bin", if (n_int == 1L) "" else "s")
      cat("  Overall summary: available over ", bins, "; pointwise interval.\n", sep = "")
    } else {
      cat("  Overall summary: ", ov$status, ": ", ov$reason, "\n", sep = "")
    }
  }
  chained_result_pt_text(ch$pt, object$DIDparams$base_period)
  inf <- ch$inference
  cat(chained_result_inference_text(inf, "estimates"))
  law_text <- chained_result_law_text(inf)
  if (!is.null(law_text)) cat("  ", law_text, "\n", sep = "")
  # Shares estimated from the unit table are the default and are not printed.
  if (!isTRUE(ch$shares_known)) {
    cat("  Cohort shares: unknown for some units of the unit table; share-dependent targets are unavailable.\n")
  } else if (identical(ch$share_law, "first_calendar_masses")) {
    cat("  Cohort shares: cohort masses at the first stored date; standard errors include their uncertainty.\n")
  }
  sel <- ch$selectors
  if (isTRUE(sel$nondefault) || isTRUE(sel$ignored)) {
    fmt <- function(x) if (is.null(x)) "NULL" else as.character(x)
    cat("  Selectors: min_e = ", fmt(sel$min_e), ", max_e = ", fmt(sel$max_e), ", balance_e = ", fmt(sel$balance_e),
        if (isTRUE(sel$ignored)) ", ignored for calendar aggregation as in did 2.5.1" else ", applied as in did 2.5.1",
        ".\n", sep = "")
  }
  invisible(object)
}

print.cdid_AGGTEobj <- function(x, ...) {
  summary(x, ...)
}

tidy.cdid_AGGTEobj <- function(x, ...) {
  out <- NextMethod()
  ch <- x$chained
  cell <- ch$cell
  rule <- function(a, s) ifelse(!is.na(s) & s > 0, a / s, NA_real_)
  cert_ids <- chained_result_pick(ch$numerical_certificate, "certified_target_ids")
  bin_cert <- if (is.null(cert_ids)) rep(NA, nrow(cell)) else cell$label %in% cert_ids
  ov_cert <- chained_result_pick(ch$overall$numerical_certificate, "certified")
  ov_cert <- if (is.null(ov_cert)) NA else isTRUE(ov_cert)
  if (ch$type == "group") {
    oz <- rule(x$overall.att, x$overall.se)
    statistic <- c(oz, unname(ch$z_stat))
    p.value <- c(2 * stats::pnorm(-abs(oz)), unname(ch$p_value))
    meta <- data.frame(
      label = c("overall", cell$label),
      available = c(identical(ch$overall$status, "available"), cell$available),
      restricted = c(length(ch$overall$excluded_bins) > 0L, cell$restricted),
      n.intended = c(length(ch$overall$intended_bins), cell$n_intended),
      n.used = c(length(ch$overall$used_bins), cell$n_used),
      excluded.ids = c(paste(ch$overall$excluded_bins, collapse = "; "), cell$excluded_ids),
      reason = c(ch$overall$reason, cell$reason),
      certified = c(ov_cert, bin_cert),
      band.family = c(ch$overall$band_family, rep(ch$band_family, nrow(cell))),
      stringsAsFactors = FALSE)
  } else if (ch$type == "simple") {
    statistic <- unname(ch$z_stat)[1L]
    p.value <- unname(ch$p_value)[1L]
    meta <- data.frame(
      label = cell$label[1L], available = cell$available[1L], restricted = cell$restricted[1L],
      n.intended = cell$n_intended[1L], n.used = cell$n_used[1L], excluded.ids = cell$excluded_ids[1L],
      reason = cell$reason[1L], certified = bin_cert[1L], band.family = ch$band_family,
      stringsAsFactors = FALSE)
  } else {
    statistic <- unname(ch$z_stat)
    p.value <- unname(ch$p_value)
    meta <- data.frame(
      label = cell$label, available = cell$available, restricted = cell$restricted,
      n.intended = cell$n_intended, n.used = cell$n_used, excluded.ids = cell$excluded_ids,
      reason = cell$reason, certified = bin_cert, band.family = rep(ch$band_family, nrow(cell)),
      stringsAsFactors = FALSE)
  }
  chained_result_need(nrow(meta) == nrow(out), "Internal: tidy rows do not match the aggregation record.")
  out$statistic <- statistic
  out$p.value <- p.value
  rownames(meta) <- NULL
  if (chained_result_show_estimator(x)) out$estimator <- ch$estimator
  cbind(out, meta)
}

glance.cdid_AGGTEobj <- function(x, ...) {
  out <- NextMethod()
  ch <- x$chained
  out$estimator <- ch$estimator
  out$status <- ch$status
  out$na.rm <- ch$na_rm
  out$n.bins <- nrow(ch$cell)
  out$n.bins.available <- sum(ch$cell$available)
  out$overall.status <- ch$overall$status
  out$band.family <- ch$band_family
  out$inference <- ch$inference$kind
  out
}

ggdid.cdid_AGGTEobj <- function(object, ylim = NULL, xlab = NULL, ylab = NULL, title = "",
                                xgap = 1, legend = TRUE, ref_line = 0, theming = TRUE, ...) {
  # Upstream refuses type = "simple"; that refusal is kept (no new plot feature).
  plain <- did::AGGTEobj(overall.att = object$overall.att, overall.se = object$overall.se,
                         type = object$type, egt = object$egt, att.egt = object$att.egt,
                         se.egt = object$se.egt, crit.val.egt = object$crit.val.egt,
                         inf.function = NULL, min_e = object$min_e, max_e = object$max_e,
                         balance_e = object$balance_e, call = object$call, DIDparams = object$DIDparams)
  p <- did::ggdid(plain, ylim = ylim, xlab = xlab, ylab = ylab, title = title, xgap = xgap,
                  legend = legend, ref_line = ref_line, theming = theming, ...)
  if (chained_result_show_estimator(object))
    p <- p + ggplot2::labs(caption = object$chained$estimator)
  cell <- object$chained$cell
  attr(p, "cdid_presentation_bins") <- data.frame(
    key = object$egt, label = cell$label, att = object$att.egt, se = object$se.egt,
    status = ifelse(cell$available, ifelse(cell$restricted, "restricted", "available"), "unavailable"),
    reason = cell$reason, stringsAsFactors = FALSE)
  p
}

# Shared by the MP and aggregation summary/print methods: the parallel trends
# assumption and the separate base_period reporting choice.
chained_result_pt_text <- function(pt, base_period = NULL) {
  if (!is.list(pt) || !is.character(pt$pt_assumption)) return(invisible(NULL))
  bp <- if (is.character(base_period) && length(base_period) == 1L) paste0("; base_period = \"", base_period, "\"") else ""
  if (!identical(pt$pt_assumption, "all")) {
    cat("  Parallel trends: PT-post", bp, ".\n", sep = "")
    return(invisible(NULL))
  }
  cat("  Parallel trends: PT-all from common stored date ", fmt_design_num(pt$pt_start),
    ", in addition to every PT-post restriction", bp, ".\n", sep = "")
  if (isTRUE(pt$no_extra_restrictions)) cat("    No additional restrictions (PT-post-equivalent model).\n") else
    cat("    Added pre-period restrictions: ", pt$n_restricted_pre_parameters, ".\n", sep = "")
  by <- pt$by_cohort
  if (is.data.frame(by)) for (j in seq_len(nrow(by))) cat("    Cohort ", fmt_design_num(by$g[j]), ": reference period ",
    if (is.na(by$ref[j])) "unavailable" else fmt_design_num(by$ref[j]),
    "; stored window dates ", by$n_extra_calendar_dates[j], "; represented restrictions ", by$n_restricted_pre_parameters[j], ".\n", sep = "")
  invisible(NULL)
}
