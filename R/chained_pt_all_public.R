# CP03 / CP10 request validation and descriptive assumption records. Estimation
# and restriction algebra remain in the accepted calendar design/GMM helpers.
validate_chained_pt_request <- function(pt_assumption, pt_start = NULL, calendar = NULL) {
  assumption <- tryCatch(match.arg(pt_assumption, c("post", "all", "ctrl")), error = function(e)
    chained_request_abort("`pt_assumption` must be \"post\", \"all\", or \"ctrl\".", "cdid_chained_request_input_error"))
  bad <- function(message) chained_request_abort(message,
    c("cdid_chained_request_pt_start", "cdid_chained_request_input_error"))
  if (identical(assumption, "post")) {
    if (!is.null(pt_start)) bad("`pt_start` is inapplicable under chained `pt_assumption = \"post\"`; use NULL.")
  } else if (identical(assumption, "all")) {
    if (!design_vector(pt_start, 1L))
      bad("PT-all requires an explicit finite numeric scalar `pt_start` in the stored calendar (D33).")
    pt_start <- as.double(pt_start)
    if (!is.null(calendar) && !pt_start %in% calendar)
      bad("`pt_start` must equal a stored calendar date exactly; no rounding, clipping or first-date fallback is used (D33).")
  } else {
    if (!is.null(pt_start)) {
      if (!design_vector(pt_start, 1L))
        bad("Under `pt_assumption = \"ctrl\"`, `pt_start` must be `NULL` (defaulting to min(calendar)) or a finite numeric scalar in the stored calendar.")
      pt_start <- as.double(pt_start)
      if (!is.null(calendar) && !pt_start %in% calendar)
        bad("`pt_start` must equal a stored calendar date exactly under `pt_assumption = \"ctrl\"`.")
    } else if (!is.null(calendar)) {
      pt_start <- as.double(min(calendar))
    }
  }
  list(pt_assumption = assumption, pt_start = pt_start)
}

chained_pt_reporting_tags <- function(meta, pt_assumption, pt_start, anticipation, descriptors, catalog_links) {
  if (!pt_assumption %in% c("all", "ctrl")) return(meta)
  vals_t <- c(meta$time, meta$base_t)
  vals_t <- vals_t[is.finite(vals_t)]
  eff_start <- if (!is.null(pt_start)) pt_start else if (length(vals_t)) min(vals_t) else Inf
  lo <- pmin(meta$time, meta$base_t); hi <- pmax(meta$time, meta$base_t)
  pre <- meta$role == "pre"
  ant <- vapply(seq_len(nrow(meta)), function(j)
    pre[j] && is.finite(hi[j]) && exact_cmp_sum(hi[j], anticipation, meta$g[j]) >= 0L, logical(1L))
  window <- pre & !ant & !is.na(lo) & lo >= eff_start
  meta$pt_window_tag <- ifelse(meta$role == "post", "causal", ifelse(meta$role == "reference", "reference",
    ifelse(is.na(lo), "unavailable_reference", ifelse(ant, "anticipation", ifelse(window, "in_pt_window", "pre_pt_start")))))
  # A calendar tag cannot establish PT for an auxiliary NYT pool. Compare the
  # actual operator's cohort set with the retained catalog's operators; equality
  # uses the same frame/covariate conventions already enforced by the pipeline.
  meta$pt_operator_covered <- rep(NA, nrow(meta))
  for (j in which(pre)) {
    q <- meta$auxiliary_index[j]
    covered <- FALSE
    if (!is.na(q) && q >= 1L && q <= length(descriptors)) {
      d <- descriptors[[q]]
      if (identical(d$K, Inf)) covered <- TRUE
      if (!is.null(catalog_links)) {
        if ("link_role" %in% names(catalog_links)) {
          covered <- covered || (identical(meta$support_status[j], "supported_cross_pair_pretrend") &&
            any(catalog_links$g == d$g & catalog_links$link_role == "control"))
        } else {
          candidates <- which(catalog_links$g == d$g)
          covered <- covered || any(vapply(candidates, function(k) {
            parsed_k <- parse_control_pool_string(catalog_links$control_pool[k])
            !is.null(parsed_k) && identical(sort(as.double(parsed_k)), sort(d$K))
          }, logical(1L)))
        }
      }
    }
    meta$pt_operator_covered[j] <- covered
  }
  meta
}

chained_pt_assumption_record <- function(res, prep, request) {
  if (identical(request$pt_assumption, "ctrl")) {
    cal <- prep$calendar; gs <- prep$target_cohorts
    start <- if (is.null(request$pt_start)) min(cal) else as.double(request$pt_start)
    start_sel <- if (is.null(request$pt_start)) "calendar_minimum_default" else "fixed exogenously"
    ant_val <- if (!is.null(res$controls$anticipation)) res$controls$anticipation else 0
    refs <- vapply(gs, function(g) tryCatch(calendar_reference_date(g, cal, ant_val),
      cdid_error = function(e) NA_real_), numeric(1L))
    by <- data.frame(g = gs, ref = refs, pt_start = start,
      baseline_status = ifelse(is.na(refs), "no_admissible_reference", "stored_reference"),
      n_extra_calendar_dates = vapply(refs, function(b) if (is.na(b)) 0L else sum(cal >= start & cal < b), integer(1L)),
      n_shared_pre_dates = 0L, n_merged_pre_parameters = 0L, stringsAsFactors = FALSE)
    spec <- if (!is.null(res$design_all)) res$design_all$pt_ctrl_spec else NULL
    if (!is.null(spec) && is.data.frame(spec$by_cohort)) {
      hit <- match(gs, spec$by_cohort$g); ok <- !is.na(hit)
      for (nm in intersect(c("n_shared_pre_dates", "n_merged_pre_parameters"), names(spec$by_cohort)))
        by[[nm]][ok] <- spec$by_cohort[[nm]][hit[ok]]
    }
    by$no_extra_calendar_window <- by$n_extra_calendar_dates == 0L
    by$no_extra_restrictions <- by$n_merged_pre_parameters == 0L
    return(list(pt_assumption = "ctrl", pt_start = start, start_selection = start_sel,
      by_cohort = by, n_extra_calendar_dates = sum(by$n_extra_calendar_dates),
      n_shared_pre_dates = sum(by$n_shared_pre_dates),
      n_merged_pre_parameters = sum(by$n_merged_pre_parameters),
      n_restricted_pre_parameters = sum(by$n_merged_pre_parameters),
      no_extra_calendar_window = all(by$no_extra_calendar_window),
      no_extra_restrictions = all(by$no_extra_restrictions),
      restriction_is_identity = all(by$no_extra_restrictions),
      assumptions = paste0("Control-cohort parallel trends (PT-ctrl) from date ", fmt_design_num(start),
        " (", start_sel, "): clean later-treated control cohorts k > g share conditional untreated trends from ",
        fmt_design_num(start), " onward, while target cohort g's own pre-reference trends remain unrestricted."),
      pretest_interpretation = paste0("Target cohort pre-reference trends remain unrestricted under PT-ctrl, so pretrend contrasts test whether target cohort g shares trends with its controls prior to reference.")))
  }
  if (!identical(request$pt_assumption, "all"))
    return(list(pt_assumption = "post", pt_start = NULL))
  cal <- prep$calendar; gs <- prep$target_cohorts; start <- request$pt_start
  refs <- vapply(gs, function(g) tryCatch(calendar_reference_date(g, cal, res$controls$anticipation),
    cdid_error = function(e) NA_real_), numeric(1L))
  by <- data.frame(g = gs, ref = refs, pt_start = start,
    baseline_status = ifelse(is.na(refs), "no_admissible_reference", "stored_reference"),
    n_extra_calendar_dates = vapply(refs, function(b) if (is.na(b)) 0L else sum(cal >= start & cal < b), integer(1L)),
    n_extra_represented_dates = 0L, n_restricted_pre_parameters = 0L,
    n_zero_loading_links = 0L, stringsAsFactors = FALSE)
  spec <- res$design_all$pt_all_spec
  if (!is.null(spec)) {
    hit <- match(gs, spec$by_cohort$g); ok <- !is.na(hit)
    for (nm in c("n_extra_represented_dates", "n_restricted_pre_parameters", "n_zero_loading_links"))
      by[[nm]][ok] <- spec$by_cohort[[nm]][hit[ok]]
  }
  by$no_extra_calendar_window <- by$n_extra_calendar_dates == 0L
  by$no_extra_restrictions <- by$n_restricted_pre_parameters == 0L
  list(pt_assumption = "all", pt_start = start, start_selection = "fixed exogenously",
    by_cohort = by, n_extra_calendar_dates = sum(by$n_extra_calendar_dates),
    n_restricted_pre_parameters = sum(by$n_restricted_pre_parameters),
    n_zero_loading_links = sum(by$n_zero_loading_links),
    no_extra_calendar_window = all(by$no_extra_calendar_window),
    no_extra_restrictions = all(by$no_extra_restrictions),
    restriction_is_identity = all(by$no_extra_restrictions),
    free_post_param_indices = spec$free_post_param_indices,
    restricted_post_param_indices = spec$restricted_post_param_indices,
    assumptions = paste0("Common start fixed exogenously. Baseline PT-post assumptions are retained at every qualified causal date, including dates before the start. ",
      "Added untreated-trend equalities use each cohort's original baseline and the same catalog control operators from the start; sampling, selection, overlap and nuisance premises remain required."),
    pretest_interpretation = paste0("Unrestricted direct pre-comparisons retain the full requested nonreference pretrend null. Pre-start contrasts can reject without contradicting PT-all from the chosen start; anticipation effects need not be zero. ",
      "An in-window tag implies a zero contrast only when the actual control operator and selection/trend premises are covered; an auxiliary NYT pool absent from catalog assumptions is not covered by its date tag. Operator coverage records membership in the declared assumptions, not empirical verification of PT or selection. No subset or overidentification test is introduced."))
}
