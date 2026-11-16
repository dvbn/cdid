# CP08 D60/D61 reference inputs (approved contract and reference TeX, 22 October
# 2026). Request layer: the fix_weights contract, raw provenance of the sampling
# weight and of every outcome/selection-formula variable, the first-calendar
# weight V by original id, and a compact source of genuinely dated values from
# the actual (id, date) rows, including Y = NA rows. Pipeline layer: recipe
# identities, canonical required (date, id) rows of the fixed
# post-plus-both-auxiliary union on each route's required supports, ONE frozen
# transform state per distinct formula on that union, and a sealed
# recipe-indexed reader of designs, anchored weights and recipe cohort rates
# replacing the static Z_a / w / a reads of every backend on the dated path.
#
# Calls the existing static gates accept run the existing path unchanged for
# NULL, "first_period" and "base_period". Unit-attribute X with "first_period"
# and w_t1 on the full roster is the existing static system on w_t1. Any other
# dated recipe runs on the route's full-roster backend: baseline reg/dr/ipw and
# custom callbacks, conditional selection ipw/reg/dr (pair or date logits, every
# pooling), pure imputation ipw/reg/dr and the augmented correction.
#
# Required supports per cohort k of each view (recipe r): designs on observed
# pair rows (baseline, conditional selection), the selection roster U
# (selection designs) or every F row (imputation/augmented); weights on A rows
# (baseline rates), U rows (selection fits and rates a_k^U) or F rows
# (imputation masses and rates). Nothing else is demanded.

chained_reference_lookup <- function(keys, table) {
  # Exact sorted lookup (table strictly increasing); NA when absent.
  fast <- cdid_fast_sel_ref_lookup(keys, table)
  if (!is.null(fast)) return(fast)
  pos <- findInterval(keys, table)
  hit <- pos > 0L
  hit[hit] <- table[pos[hit]] == keys[hit]
  out <- rep(NA_integer_, length(keys))
  out[hit] <- pos[hit]
  out
}

chained_reference_fmt <- function(x) if (is.numeric(x)) format(x, digits = 12) else as.character(x)

# ------------------------------------------------------------------------------
# 1. Request contract
# ------------------------------------------------------------------------------

# NULL, "first_period" and "base_period" are accepted. "varying" (period-specific
# weights) stays refused (D60 defers it); any other value is invalid. Evaluated
# in the scope gate before any data promise is forced.
chained_reference_fix_weights <- function(fix_weights) {
  if (is.null(fix_weights)) return(NULL)
  chained_request_need(chained_request_chr1(fix_weights) &&
      fix_weights %in% c("first_period", "base_period", "varying"),
    '`fix_weights` must be NULL, "first_period" or "base_period" in chained mode.',
    "cdid_chained_request_unsupported_weights")
  chained_request_need(!identical(fix_weights, "varying"),
    paste0('fix_weights = "varying" (period-specific sampling weights) is deferred for chained estimation (D60) and is ',
      'refused rather than ignored; weights that change within unit need "first_period" or "base_period".'),
    "cdid_chained_request_unsupported_weights")
  fix_weights
}

# Public policy for an unavailable required reference-date covariate.
chained_reference_missing_policy <- function(missing_reference) {
  chained_request_need(chained_request_chr1(missing_reference) && missing_reference %in% c("error", "exclude"),
    '`missing_reference` must be "error" or "exclude".',
    c("cdid_chained_request_missing_reference_error", "cdid_chained_request_input_error"))
  missing_reference
}

# ------------------------------------------------------------------------------
# 2. Raw provenance (exact comparisons on supplied values; NA = not supplied)
# ------------------------------------------------------------------------------

# Position of the first data row whose supplied value differs from the first
# supplied value of its unit, or NA. Exact comparison; factors by label.
chained_reference_first_change <- function(x, unit) {
  if (is.factor(x)) x <- as.character(x)
  keep <- which(!is.na(x))
  if (!length(keep)) return(NA_integer_)
  u <- unit[keep]; v <- x[keep]
  bad <- which(v != v[match(u, u)])
  if (length(bad)) keep[bad[1L]] else NA_integer_
}

chained_reference_change_text <- function(x, unit, id, tt, row) {
  if (is.factor(x)) x <- as.character(x)
  keep <- which(!is.na(x) & unit == unit[row])
  first <- keep[1L]
  paste0("id ", chained_reference_fmt(id[row]), ": ", chained_reference_fmt(x[first]), " at date ",
    chained_reference_fmt(tt[first]), " but ", chained_reference_fmt(x[row]), " at date ", chained_reference_fmt(tt[row]))
}

# Sampling weights. The existing static gate and messages are unchanged
# whenever supplied data weights do not change within any unit. Changing data
# weights are returned as a dated source: a frame cannot override them (frame
# supply makes a weight a dateless unit attribute that must agree with data).
chained_reference_weights <- function(data, frame, weightsname, id, tt, calendar, frame_ids,
                                      in_panel, frame_only, first, pos_unit, frame_supplied = TRUE) {
  n_frame <- length(frame_ids)
  frame_w <- rep(1.0, n_frame)
  if (is.null(weightsname)) return(list(kind = "static", frame_w = frame_w, source = "none", name = NULL))
  chained_request_need(chained_request_chr1(weightsname), "`weightsname` must be a single column name.")
  in_data <- weightsname %in% names(data); in_frame <- weightsname %in% names(frame)
  chained_request_need(in_data || in_frame, paste0("`weightsname` column ", weightsname, " not found in `data` or `frame`."),
    "cdid_chained_request_weight_error")
  if (in_data) {
    wd <- data[[weightsname]]
    change <- if (is.numeric(wd) && !is.object(wd)) chained_reference_first_change(wd, pos_unit) else NA_integer_
    if (!is.na(change)) {
      sup <- which(!is.na(wd))
      chained_request_need(all(is.finite(wd[sup])) && all(wd[sup] > 0), "Sampling weights must be finite and positive.",
        "cdid_chained_request_weight_error")
      example <- chained_reference_change_text(wd, pos_unit, id, tt, change)
      chained_request_need(!(frame_supplied && in_frame),
        paste0("`frame` weights disagree with `data` weights for some panel units (aligned by id): data weights change ",
          "within unit (", example, ") and a frame value is a dateless unit attribute that cannot override them."),
        "cdid_chained_request_weight_error")
      pos <- match(id[sup], frame_ids)
      di <- match(as.double(tt[sup]), calendar)
      key <- as.double(di) * (n_frame + 1) + pos
      o <- order(key)
      V <- rep(NA_real_, n_frame)
      r1 <- sup[di == 1L]
      V[match(id[r1], frame_ids)] <- as.double(wd[r1])
      return(list(kind = "varying", frame_w = NULL, source = "data", name = weightsname,
        key = key[o], value = as.double(wd[sup][o]), V = V, example = example))
    }
    chained_request_need(is.numeric(wd) && all(is.finite(wd)) && all(wd > 0), "Sampling weights must be finite and positive.",
      "cdid_chained_request_weight_error")
    w_unit <- wd[first]
    chained_request_need(all(wd == w_unit[pos_unit]),
      "Sampling weights vary within unit across time; supply fix_weights = \"first_period\" or \"base_period\".",
      "cdid_chained_request_weight_error")
    frame_w[!is.na(in_panel)] <- as.double(w_unit[in_panel[!is.na(in_panel)]])
  }
  if (in_frame) {
    wf <- frame[[weightsname]]
    chained_request_need(is.numeric(wf) && all(is.finite(wf)) && all(wf > 0), "Frame weights must be finite and positive.",
      "cdid_chained_request_weight_error")
    if (in_data) chained_request_need(all(as.double(wf[!is.na(in_panel)]) == frame_w[!is.na(in_panel)]),
      "`frame` weights disagree with `data` weights for some panel units (aligned by id).", "cdid_chained_request_weight_error")
    frame_w <- as.double(wf)
  } else {
    chained_request_need(!length(frame_only),
      "Frame-only units need weights: supply the `weightsname` column in `frame` (fixed weights are never imputed).",
      "cdid_chained_request_weight_error")
  }
  list(kind = "static", frame_w = frame_w, source = if (in_frame) "frame" else "data", name = weightsname)
}

# Outcome-formula raw variables on the dated path. Frame-supplied variables
# are dateless unit attributes and must agree with every supplied data value;
# data-only variables are unit attributes unless a supplied value changes
# within some unit (then dated); the date variable takes the extraction date.
# Used only after an existing static covariate gate has failed, and only to
# decide whether a genuinely dated component exists: disagreement or absence is
# recorded here and refused by the caller (the original static refusal stands
# unchanged when nothing is dated).
chained_reference_classify <- function(vars, data, frame, frame_supplied, id, tt, tname, frame_ids, pos_unit, label) {
  # One row per actual variable, including none (e.g. ~0).
  out <- data.frame(variable = vars, source = rep(NA_character_, length(vars)), kind = rep(NA_character_, length(vars)),
    example = rep(NA_character_, length(vars)), stringsAsFactors = FALSE)
  for (j in seq_along(vars)) {
    v <- vars[j]
    if (identical(v, tname)) {
      out$source[j] <- "date"; out$kind[j] <- "date"
      out$example[j] <- "formula occurrence of the date variable"
      next
    }
    if (frame_supplied && v %in% names(frame)) {
      out$source[j] <- "frame"; out$kind[j] <- "unit_attribute"
      if (v %in% names(data)) {
        dv <- data[[v]]; fv <- frame[[v]][match(id, frame_ids)]
        if (is.factor(dv)) dv <- as.character(dv)
        if (is.factor(fv)) fv <- as.character(fv)
        bad <- which(!is.na(dv) & (is.na(fv) | dv != fv))
        if (length(bad)) {
          out$kind[j] <- "disagree"
          out$example[j] <- paste0("id ", chained_reference_fmt(id[bad[1L]]), " at date ", chained_reference_fmt(tt[bad[1L]]))
        }
      }
      next
    }
    if (!v %in% names(data)) { out$kind[j] <- "absent"; next }
    change <- chained_reference_first_change(data[[v]], pos_unit)
    out$source[j] <- "data"
    if (is.na(change)) out$kind[j] <- "unit_attribute" else {
      out$kind[j] <- "dated"
      out$example[j] <- chained_reference_change_text(data[[v]], pos_unit, id, tt, change)
    }
  }
  out
}

# With a genuinely dated component present, a frame value cannot override
# disagreeing data values, and every variable must exist.
chained_reference_check_classes <- function(cls, label) {
  bad <- cls[cls$kind == "disagree", , drop = FALSE]
  chained_request_need(!nrow(bad),
    paste0("Covariate `", bad$variable[1L], "` in `data` disagrees with `frame` (frame supply cannot override failed ",
      "agreement): ", bad$example[1L], " (", label, ")."), "cdid_chained_request_covariate_error")
  absent <- cls$variable[cls$kind == "absent"]
  chained_request_need(!length(absent),
    paste0("Covariate(s) ", paste(absent, collapse = ", "), " of `", label, "` not found in `data` or `frame`."),
    "cdid_chained_request_covariate_error")
  invisible(cls)
}

chained_reference_dated_text <- function(prov) {
  d <- prov[prov$kind %in% c("dated", "date"), , drop = FALSE]
  paste(paste0("`", d$variable, "` (", d$example, ")"), collapse = "; ")
}

# A formula is dated when any of its variables is genuinely dated (or is the
# date variable). Formulas without such a variable keep their static design.
chained_reference_formula_dated <- function(formula, prov) {
  if (is.null(formula) || is.null(prov)) return(FALSE)
  any(all.vars(formula) %in% prov$variable[prov$kind %in% c("dated", "date")])
}

# Full-roster record of a dated pure-imputation or augmented request: the same
# declaration fields as the static constructors (chained_request_imputation_inputs,
# chained_request_augmented_inputs), without evaluating any dated design on the
# first-observed frame. The outcome/treatment columns are resolved by the
# sealed reader on the canonical required rows; nothing reaches X_frame.
chained_request_reference_selection_record <- function(kind, xf, sf, sel_norm, justification, frame_supplied, n_frame,
                                                      imputation_pooling, imputation_est_method = NULL,
                                                      selection_formula_inherited = TRUE,
                                                      augmented_est_method = "dr") {
  population <- if (frame_supplied) "supplied representative frame" else "supplied panel unit roster (including singletons)"
  aug_est <- if (is.null(augmented_est_method)) "dr" else augmented_est_method
  design <- list(route = if (identical(kind, "imputation")) 'pure actual-change imputation (selection_est_method = "imputation")' else
      paste0('augmented observation correction (selection_est_method = "dr", est_method = "', aug_est, '")'),
    formula = xf, columns = NULL, intercept = identical(as.integer(attr(stats::terms(xf), "intercept")), 1L),
    n_rows = n_frame,
    row_universe = paste("every F roster row of each needed cohort, read at each recipe's extraction date on the canonical",
      "required (date, id) rows (D60); columns, factor levels and transform state are frozen on that union"),
    common_A_alias_reduction = FALSE, reference_inputs = TRUE,
    duration_convention = "stored-calendar transition count L_e = i_t - i_s")
  out <- list(Z = matrix(1, n_frame, 1L, dimnames = list(NULL, "reader_intercept")),
    U = rep(TRUE, n_frame), rho = rep(1, n_frame), universe = "F", justification = justification,
    selection_requested = sel_norm$requested, canonical_assumption = "conditional_mtar",
    legacy_alias = sel_norm$legacy_alias, population = population, imputation_pooling = imputation_pooling,
    outcome_recipe = NULL, treatment_columns = NULL, X_frame = NULL)
  if (identical(kind, "imputation")) {
    c(out, list(formula = xf, selection_est_method = "imputation", imputation_est_method = imputation_est_method,
      imputation_design = design))
  } else {
    design$selection_formula <- sf
    design$selection_formula_inherited <- selection_formula_inherited
    c(out, list(formula = sf, selection_model = sel_norm$selection_model, selection_pooling = sel_norm$selection_pooling,
      selection_est_method = "dr", augmented_est_method = aug_est, outcome_formula = xf,
      selection_formula_inherited = selection_formula_inherited, augmented_design = design))
  }
}

# Requirement universes of each route on the dated path (see the file header).
chained_reference_route <- function(kind, selection_universe = "F") {
  switch(kind,
    baseline = list(kind = kind, x_rows = "pair", s_rows = NULL, w_rows = "A"),
    corollary1 = list(kind = kind, x_rows = "pair", s_rows = "U", w_rows = "U", universe = selection_universe),
    imputation = list(kind = kind, x_rows = "F", s_rows = NULL, w_rows = "F"),
    augmented = list(kind = kind, x_rows = "F", s_rows = "U", w_rows = "F", universe = "F"),
    chained_request_abort("Unknown reference-input route.", "cdid_chained_request_input_error"))
}

# Dateless attribute values in frame order: frame column, or the unit's
# common supplied data value (NA when the unit supplies none).
chained_reference_unit_column <- function(v, prov_row, data, frame, frame_supplied, id, frame_ids, pos = NULL) {
  if (identical(prov_row$source, "frame")) return(frame[[v]])
  x <- data[[v]]
  if (is.null(pos)) pos <- match(id, frame_ids)
  keep <- which(!is.na(x))
  rep_row <- keep[match(seq_along(frame_ids), pos[keep])]
  x[rep_row]
}

# Compact request-time source of the dated path. No values are carried forward
# or inferred: dated columns are the actual rows (with their types and factor
# levels), indexed by (calendar date, frame position). `formulas$x` is the
# outcome/treatment/imputation formula and `formulas$s` the selection formula
# (NULL when the route fits no observation probability); each carries its
# formula, whether it is genuinely dated, and for a static formula the existing
# frame-position design of the route's unchanged static gate (with the OR
# alias record when the route has one). `route` gives the required supports.
chained_reference_source <- function(formulas, prov, data, frame, frame_supplied, id, tt, tname, frame_ids, calendar,
                                     anticipation, wres, fix_weights, route, resource = NULL) {
  n_frame <- length(frame_ids)
  for (nm in names(formulas)) formulas[[nm]]$dated <- chained_reference_formula_dated(formulas[[nm]]$formula, prov)
  dated_x <- any(vapply(formulas, `[[`, logical(1), "dated"))
  long <- unit <- NULL; key <- o <- numeric()
  if (dated_x && !is.null(resource)) {
    # Long dated columns, keys, the order and one model frame of the supplied
    # rows, charged before any of them is allocated (no dense design is built).
    nv <- nrow(prov); nr <- nrow(data)
    ipw_owner_safe_sum(c(8 * nr * (2 * nv + 6), 8 * n_frame * (nv + 2)), resource$numeric_workspace_bytes,
      "reference-date source columns and model frame bytes")
    ipw_owner_safe_sum(64 * (nr + n_frame) * (nv + 2) * (1 + ceiling(log2(nr + 2))), resource$counted_work_limit,
      "reference-date source preparation work")
  }
  if (dated_x) {
    pos <- match(id, frame_ids)
    di <- match(as.double(tt), calendar)
    key <- as.double(di) * (n_frame + 1) + pos
    o <- order(key, method = "radix")
    dated_vars <- prov$variable[prov$kind == "dated"]
    unit_vars <- prov$variable[prov$kind == "unit_attribute"]
    if (length(dated_vars)) {
      long <- data[o, dated_vars, drop = FALSE]
      rownames(long) <- NULL
    }
    if (length(unit_vars)) {
      u <- lapply(unit_vars, function(v) chained_reference_unit_column(v, prov[prov$variable == v, , drop = FALSE],
        data, frame, frame_supplied, id, frame_ids, pos = pos))
      names(u) <- unit_vars
      unit <- as.data.frame(u, stringsAsFactors = FALSE, optional = TRUE)
    }
  }
  # Dated transforms must not be evaluated on the supplied long rows or
  # first-observed frame. Their width is resolved on the canonical required
  # union by the reader. Initial public planning uses only the carrier width;
  # the exact reader and backend dimensions receive their own preallocation gates.
  p_bound <- if (dated_x) NULL else max(vapply(formulas, function(f)
    if (is.null(f$static)) 0 else as.double(ncol(f$static)), numeric(1)))
  weight <- if (identical(wres$kind, "varying")) list(kind = "varying", name = wres$name, key = wres$key, value = wres$value) else
    list(kind = "static", name = wres$name, frame_w = wres$frame_w)
  V <- if (identical(wres$kind, "varying")) wres$V else wres$frame_w
  structure(list(formulas = formulas, formula = formulas$x$formula, tname = tname, prov = prov,
    design_kind = if (dated_x) "dated" else "static",
    long = long, long_key = key[o], unit = unit, N = n_frame, frame_ids = frame_ids, calendar = calendar,
    anticipation = as.double(anticipation), fix_weights = fix_weights, weight = weight, V = V,
    route = route, p_bound = p_bound, weight_example = wres$example),
    class = "cdid_reference_source")
}

# One row per requested (date index, frame position) key, raw variables as
# actual values: dated from the exact row (NA if absent), attributes from the
# unit, the date variable equal to the extraction date. `vars` restricts the
# columns to one formula's raw variables.
chained_reference_rows_frame <- function(src, keys, vars = src$prov$variable) {
  N1 <- src$N + 1
  d <- floor(keys / N1); i <- keys - d * N1
  cols <- list()
  for (j in which(src$prov$variable %in% vars)) {
    v <- src$prov$variable[j]
    cols[[v]] <- switch(src$prov$kind[j],
      date = src$calendar[d],
      dated = src$long[[v]][chained_reference_lookup(keys, src$long_key)],
      unit_attribute = src$unit[[v]][i])
  }
  df <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  if (!length(cols)) df <- data.frame(row.names = seq_along(keys))
  df
}

# Rows of `df` (chained_reference_rows_frame) with any unknown raw value.
chained_reference_rows_unknown <- function(df) {
  miss <- logical(nrow(df))
  for (j in seq_along(df)) {
    col <- df[[j]]
    miss <- miss | (if (is.null(dim(col))) is.na(col) else rowSums(is.na(col)) > 0L)
  }
  miss
}

# ------------------------------------------------------------------------------
# Issue #16: cohort-specific reference eligibility (missing_reference = "exclude")
# ------------------------------------------------------------------------------
# An observed-data eligibility rule, never imputation. For each target g with
# an admissible reference date b_g, E_g holds the roster units whose raw values
# of every variable of every formula the sealed reader evaluates at b_g (the
# dated outcome/treatment formula and, when active and dated, the selection
# formula) are supplied: a dated variable at the exact (b_g, id) row (an absent
# row or NA is unknown; a covariate-only row with an NA outcome is known), a
# unit attribute at its established dateless value, the date variable always.
# Nothing is read from another date or inferred outside the data; static
# formulas keep their existing gates and never restrict E_g. Targets sharing
# b_g share one context. Only exclusions are stored (sorted frame positions),
# with compact counts by target and cohort.
chained_reference_eligibility <- function(src, targets, frame_g, resource = NULL,
                                          affected_dates = NULL) {
  if (!inherits(src, "cdid_reference_source") || !identical(src$design_kind, "dated")) return(NULL)
  dated <- Filter(function(f) isTRUE(f$dated), src$formulas)
  vars <- unique(unlist(lapply(dated, function(f) all.vars(f$formula)), use.names = FALSE))
  kinds <- src$prov$kind[match(vars, src$prov$variable)]
  vars <- vars[!is.na(kinds) & kinds %in% c("dated", "unit_attribute")]
  cal <- src$calendar; N <- src$N
  tg <- sort(unique(as.double(targets)))
  b <- vapply(tg, function(g) tryCatch(as.double(calendar_reference_date(g, cal, src$anticipation)),
    cdid_design_reference_error = function(e) NA_real_), numeric(1))
  ref_index <- match(b, cal)
  ctx <- sort(unique(ref_index[!is.na(ref_index)]))
  if (!is.null(affected_dates)) ctx <- intersect(ctx, match(affected_dates, cal))
  if (!is.null(resource) && length(ctx) && length(vars))
    ipw_owner_safe_sum(8 * as.double(N) * (length(vars) + 3), resource$numeric_workspace_bytes,
      "reference eligibility context scan bytes")
  excluded <- vector("list", length(ctx)); n_by_var <- character(length(ctx))
  retained_bytes <- 0
  for (c in seq_along(ctx)) {
    excluded[[c]] <- integer()
    if (!length(vars)) next
    df <- chained_reference_rows_frame(src, ctx[c] * (N + 1) + seq_len(N), vars)
    excluded[[c]] <- which(chained_reference_rows_unknown(df))
    if (!is.null(resource))
      retained_bytes <- ipw_owner_safe_sum(c(retained_bytes, 4 * length(excluded[[c]])),
        resource$numeric_workspace_bytes - 8 * as.double(N) * (length(vars) + 3),
        "reference eligibility retained exclusions")
    n_by_var[c] <- paste0(vars, "=", vapply(vars, function(v) sum(chained_reference_rows_unknown(df[v])), numeric(1)),
      collapse = ", ")
  }
  n_excl <- lengths(excluded)
  ids <- src$frame_ids
  head_ids <- function(pos) if (!length(pos)) "" else paste(c(vapply(utils::head(ids[pos], 5L), chained_reference_fmt,
    character(1)), if (length(pos) > 5L) "..."), collapse = ", ")
  contexts <- data.frame(context = seq_along(ctx), reference_date = cal[ctx], n_excluded = as.integer(n_excl),
    missing_by_variable = n_by_var, excluded_ids_head = vapply(excluded, head_ids, character(1)),
    stringsAsFactors = FALSE)
  tc <- match(ref_index, ctx)
  per_target <- lapply(seq_along(tg), function(i) {
    ex <- if (is.na(tc[i])) integer() else excluded[[tc[i]]]
    own <- which(frame_g == tg[i])
    data.frame(g = tg[i], reference_date = b[i], context = tc[i], n_treated = length(own),
      n_treated_eligible = length(own) - sum(own %in% ex), n_treated_excluded = sum(own %in% ex),
      n_never_treated_excluded = sum(is.infinite(frame_g[ex])),
      n_other_cohort_excluded = sum(is.finite(frame_g[ex]) & frame_g[ex] != tg[i], na.rm = TRUE),
      restricted = length(ex) > 0L, excluded_treated_ids_head = head_ids(ex[ex %in% own]),
      stringsAsFactors = FALSE)
  })
  keep <- n_excl > 0L
  structure(list(policy = "exclude", variables = vars, restricted = any(keep),
    targets = do.call(rbind, per_target), contexts = contexts,
    # Compact sealed-reader form: calendar indices of the restricted contexts
    # and their sorted excluded frame positions (unrestricted contexts absent).
    index = as.integer(ctx[keep]), excluded = excluded[keep],
    target_cohorts = tg, target_ref_index = as.integer(ref_index)),
    class = "cdid_reference_eligibility")
}

# Compact form carried by the reference source and the sealed reader.
chained_reference_eligibility_core <- function(el) {
  if (is.null(el) || !isTRUE(el$restricted)) return(NULL)
  list(index = el$index, excluded = el$excluded, target_cohorts = el$target_cohorts,
    target_ref_index = el$target_ref_index)
}

# Restricted cohort masses for aggregation: for each restricted target cohort g
# the frame positions of its own units outside E_g. Aggregation then uses
# m_g = H(V I_g E_g) with contribution V I_g E_g - m_g (D61 masses otherwise).
chained_reference_share_eligibility <- function(el, frame_g) {
  if (is.null(el) || !length(el$index)) return(NULL)
  ctx <- match(el$target_ref_index, el$index)
  keep <- which(!is.na(ctx))
  if (!length(keep)) return(NULL)
  ex <- lapply(keep, function(i) { e <- el$excluded[[ctx[i]]]; e[which(frame_g[e] == el$target_cohorts[i])] })
  list(cohort = as.double(el$target_cohorts[keep]), excluded = lapply(ex, as.integer))
}

# Compact MP record retained when reference eligibility handling runs.
chained_reference_eligibility_record <- function(prep, res = NULL) {
  el <- prep$reference_eligibility
  aux <- if (is.null(res)) NULL else res$nuisance_registry$cp04_eligibility
  status <- if (isTRUE(el$restricted)) "restricted" else if (NROW(aux$auxiliary_unavailable)) {
    "auxiliary_unavailable"
  } else if (is.null(el)) "not_applicable" else "unrestricted"
  list(policy = "exclude", status = status,
    meaning = switch(status,
      not_applicable = paste("no genuinely dated covariate formula: the existing static gates apply and no unit is",
        "excluded from any context"),
      unrestricted = "No required reference covariate is missing for the fitted comparisons.",
      auxiliary_unavailable = "Some pretrend comparisons need missing reference covariates. The main comparison samples are unchanged.",
      restricted = paste("each target g uses the units with known raw covariates at its reference date b_g (E_g)",
        "for treated and comparison roles in all of its main comparisons; the target is the ATT of its eligible",
        "treated members; units keep their ids, rows and influence-function rows, and N is unchanged")),
    variables = el$variables, targets = el$targets, contexts = el$contexts,
    restricted_targets = if (is.null(el)) numeric() else el$targets$g[el$targets$restricted],
    restricted_treated_cohorts = if (is.null(el)) numeric() else el$targets$g[el$targets$n_treated_excluded > 0L],
    auxiliary_unavailable = aux$auxiliary_unavailable)
}

# Eligibility index (calendar index of b_g, 0 when unrestricted) of target g.
chained_reference_eligibility_of <- function(el, b_index) {
  if (is.null(el) || is.na(b_index) || !(b_index %in% el$index)) 0L else as.integer(b_index)
}

# Frame rows of a recipe's eligibility context (all rows when unrestricted).
chained_reference_keep_eligible <- function(el, e, rows) {
  if (is.null(el) || is.na(e) || e == 0L) return(rows)
  ci <- match(e, el$index)
  if (is.na(ci)) chained_reference_reader_abort("Unknown reference eligibility context.", "integrity")
  rows[!(rows %in% el$excluded[[ci]])]
}

# Per-target occupancy reduction of a restricted fit: for every restricted
# context, the (cohort, s, t) pair counts of its excluded A units, read from
# the sealed panel owner's postings. The catalog subtracts these from the
# owner occupancy for the targets of that context (other targets keep the
# owner occupancy unchanged). Sparse: only excluded units are visited.
chained_reference_eligibility_occupancy <- function(el, panel_owner, frame_id, targets) {
  if (is.null(el) || !length(el$index)) return(NULL)
  post <- panel_owner$postings; cal <- panel_owner$all_times; Tn <- length(cal)
  a_of <- match(frame_id, panel_owner$a_id)
  coh_of_a <- panel_owner$units$g[match(panel_owner$a_id, panel_owner$units$id)]
  cells <- list()
  for (c in seq_along(el$index)) {
    a <- a_of[el$excluded[[c]]]; a <- sort(a[!is.na(a)])
    if (!length(a)) next
    seen <- matrix(FALSE, length(a), Tn)
    for (d in seq_len(Tn)) if (post$date_count[d] > 0L) {
      seen[, d] <- !is.na(ipw_owner_match_sorted_double(as.double(a), post$post_a_index, post$date_start[d], post$date_end[d]))
    }
    for (s in seq_len(Tn - 1L)) for (t in seq.int(s + 1L, Tn)) {
      hit <- seen[, s] & seen[, t]
      if (!any(hit)) next
      tab <- table(coh_of_a[a[hit]])
      cells[[length(cells) + 1L]] <- data.frame(context = c, cohort = as.double(names(tab)), s = cal[s], t = cal[t],
        count = as.integer(tab), stringsAsFactors = FALSE)
    }
  }
  sub <- if (length(cells)) do.call(rbind, cells) else
    data.frame(context = integer(), cohort = numeric(), s = numeric(), t = numeric(), count = integer())
  ti <- match(as.double(targets), el$target_cohorts)
  ctx <- match(el$target_ref_index[ti], el$index)
  ctx[is.na(ctx)] <- 0L
  list(target_context = as.integer(ctx), subtract = sub)
}

# Same column convention as cp04_public_design: backticked single-column terms
# keep their raw names when unambiguous.
chained_reference_model_matrix <- function(tf, df, model_frame = NULL) {
  term_labels <- attr(tf, "term.labels")
  m <- stats::model.matrix(tf, if (is.null(model_frame)) stats::model.frame(tf, df, na.action = stats::na.pass) else model_frame)
  asg <- as.integer(attr(m, "assign"))
  if (ncol(m) > 0L && length(term_labels) > 0L && length(asg) == ncol(m)) {
    cn <- colnames(m)
    for (k in seq_along(term_labels)) {
      idx <- which(asg == k); lbl <- term_labels[k]
      if (length(idx) == 1L && nchar(lbl) >= 3L && startsWith(lbl, "`") && endsWith(lbl, "`")) {
        raw_nm <- substr(lbl, 2L, nchar(lbl) - 1L)
        if (raw_nm %in% names(df) && !(raw_nm %in% cn[-idx])) cn[idx] <- raw_nm
      }
    }
    colnames(m) <- cn
  }
  m
}

# ------------------------------------------------------------------------------
# 3. Recipes, required supports and the sealed reader (every dated route)
# ------------------------------------------------------------------------------

chained_reference_reader_abort <- function(message, class = "input", diagnostics = NULL) {
  stop(structure(list(message = message, call = NULL, diagnostics = diagnostics),
    class = c(paste0("cdid_reference_", class, "_error"), "cdid_reference_error", "cdid_ipw_fit_error",
      "cdid_error", "error", "condition")))
}

# Recipe of each view: post rows of target g extract at b_g; auxiliary
# descriptors at the earlier endpoint of their actual ordered pair (t_{j-1} for
# varying, min(b_g, t_j) for universal reports). The weight anchor is t_1 or the
# target's b_g for every unit of the context, never a control's own b_k. Dates
# and anchors enter the identity only for genuinely dated components; the
# extraction date is common to the outcome and selection formulas (C^r).
chained_reference_recipes <- function(src, views, n_post) {
  cal <- src$calendar
  dated_x <- identical(src$design_kind, "dated")
  dated_w <- identical(src$weight$kind, "varying")
  # Issue #16: a restricted eligibility context E_g is part of the recipe
  # identity of every view of target g (post and auxiliary), so no rate,
  # nuisance, block, layout cohort or cache is shared across contexts. Absent
  # (NULL) on every unrestricted fit, whose recipes are unchanged.
  el <- src$eligibility
  rec <- lapply(seq_along(views), function(j) {
    v <- views[[j]]; g <- as.double(v$cohorts$g)
    b <- calendar_reference_date(g, cal, src$anticipation)
    d <- if (!dated_x) 0L else if (j <= n_post) match(b, cal) else match(as.double(v$pair$s), cal)
    a <- if (!dated_w) 0L else if (identical(src$fix_weights, "first_period")) 1L else match(b, cal)
    if (is.null(el)) c(d, a) else c(d, a, chained_reference_eligibility_of(el, match(b, cal)))
  })
  da <- do.call(rbind, rec)
  key <- if (is.null(el)) paste(da[, 1L], da[, 2L]) else paste(da[, 1L], da[, 2L], da[, 3L])
  u <- !duplicated(key)
  recipes <- data.frame(design_index = as.integer(da[u, 1L]), weight_index = as.integer(da[u, 2L]),
    stringsAsFactors = FALSE)
  if (!is.null(el)) recipes$eligibility_index <- as.integer(da[u, 3L])
  formula_txt <- paste(deparse(src$formulas$x$formula), collapse = " ")
  selection_txt <- if (is.null(src$formulas$s)) "" else
    paste0("|selection_formula:", paste(deparse(src$formulas$s$formula), collapse = " "))
  recipes$key <- paste0("formula:", formula_txt, selection_txt,
    "|design:", ifelse(recipes$design_index == 0L, "unit", paste0("date=", fmt_design_num(cal[pmax(1L, recipes$design_index)]))),
    "|weights:", ifelse(recipes$weight_index == 0L, if (is.null(src$weight$name)) "unit" else "fixed",
      paste0(src$fix_weights, "@", fmt_design_num(cal[pmax(1L, recipes$weight_index)]))),
    if (is.null(el)) "" else ifelse(recipes$eligibility_index == 0L, "",
      paste0("|eligibility:b=", fmt_design_num(cal[pmax(1L, recipes$eligibility_index)]))))
  list(recipes = recipes, view_recipe = match(key, key[u]))
}

# Rows of cohort k in a requirement universe: A (repeat observers), U (the
# declared selection roster) or F (every roster row).
chained_reference_cohort_rows <- function(st, k, universe) {
  rows <- st$cohort_order[seq.int(st$cohort_start[k], st$cohort_end[k])]
  switch(universe, A = rows[st$a_index[rows] > 0L], U = rows[st$U[rows]], F = rows,
    chained_reference_reader_abort("Unknown requirement universe.", "integrity"))
}

# Canonical required (date, id) keys, demanded only on actual supports of the
# route (file header): the dated outcome/treatment design (role x) on observed
# pair rows (S = 1 implies A) or on every F row, the dated selection design
# (role s) on the selection roster, and weights at the context anchor on the
# rate/mass universe of each support cohort. Nothing else is demanded here;
# D61 first-calendar masses are checked at aggregation by original id.
chained_reference_format_labels <- function(codes, ctx) {
  if (!length(codes)) return(character())
  u_codes <- unique(codes)
  J1 <- length(ctx$st$cohorts) + 1L
  k_vec <- u_codes %% J1
  q_vec <- u_codes %/% J1
  p_vec <- (q_vec %% 4L) + 1L
  j_vec <- q_vec %/% 4L
  fmt_u <- character(length(u_codes))
  for (idx in seq_along(u_codes)) {
    jj <- j_vec[idx]; kj <- k_vec[idx]
    v <- ctx$views[[jj]]
    tag <- paste0("support (", fmt_design_num(v$pair$s), "->", fmt_design_num(v$pair$t), ", cohort ",
      fmt_design_num(ctx$st$cohorts[kj]), ") of ", if (jj <= ctx$n_post) "post link " else "auxiliary report ",
      if (jj <= ctx$n_post) jj else jj - ctx$n_post, " (target ", fmt_design_num(v$cohorts$g), ")")
    fmt_u[idx] <- switch(p_vec[idx],
      tag,
      paste0("prediction roster of ", tag),
      paste0("selection fitting roster of ", tag),
      paste0(ctx$w_text, tag)
    )
  }
  fmt_u[match(codes, u_codes)]
}

chained_reference_requirements <- function(src, reader, views, view_recipe, recipes, ctrl, n_post, lazy_labels = FALSE,
                                           keep_raw = FALSE) {
  st <- reader$state; N1 <- src$N + 1; route <- src$route
  # Issue #16: every universe of a restricted recipe is intersected with its
  # eligibility context; cohort and observed-pair caches carry that identity.
  el <- src$eligibility
  x_dated <- isTRUE(src$formulas$x$dated)
  s_dated <- !is.null(src$formulas$s) && isTRUE(src$formulas$s$dated)
  J <- length(st$cohorts); J1 <- J + 1L
  n_cal <- length(src$calendar)
  # Observed-pair rows are cached once per (s, t, cohort): an integer slot on
  # small calendars, the exact string key otherwise. A cross-pair cell view may
  # demand pair rows for some cohorts under any route, so the cache is
  # allocated at the first pair demand instead of from route$x_rows alone.
  use_int_scan <- x_dated && (as.double(n_cal) * n_cal * J <= 65536)
  scan <- NULL
  n_max_slots <- 2L * length(views)
  keys <- list(x = vector("list", n_max_slots), s = vector("list", n_max_slots), w = vector("list", n_max_slots))
  labs <- keys
  cnt <- c(x = 0L, s = 0L, w = 0L)
  add <- function(role, k, code) {
    idx <- cnt[[role]] + 1L
    cnt[[role]] <<- idx
    keys[[role]][[idx]] <<- k
    labs[[role]][[idx]] <<- rep.int(code, length(k))
  }
  w_text <- switch(route$w_rows, A = "rate denominator of ", U = "selection roster and repeat-rate denominator of ",
    F = "roster mass and repeat-rate denominator of ")
  cohort_cache <- list()
  get_cohort_rows <- function(k, u, e = 0L) {
    ck <- if (e == 0L) paste0(u, "_", k) else paste0(u, "_", k, "_e", e)
    res <- cohort_cache[[ck]]
    if (is.null(res)) {
      res <- chained_reference_cohort_rows(st, k, u)
      if (e != 0L) res <- chained_reference_keep_eligible(el, e, res)
      cohort_cache[[ck]] <<- res
    }
    res
  }
  scan_pair <- function(k, e, pair) {
    a_rows <- get_cohort_rows(k, "A", e); obs <- integer()
    if (length(a_rows)) for (s0 in seq.int(1L, length(a_rows), by = max(1L, ctrl$max_chunk_rows))) {
      rr <- a_rows[seq.int(s0, min(length(a_rows), s0 + ctrl$max_chunk_rows - 1L))]
      obs <- c(obs, rr[cp04_pair_rows(reader, rr, pair)$S == 1])
    }
    obs
  }
  scan_e <- list()
  x_is_pair <- identical(route$x_rows, "pair")
  for (j in seq_along(views)) {
    v <- views[[j]]; r <- view_recipe[j]
    e <- if (is.null(el)) 0L else recipes$eligibility_index[r]
    is_cp_cell <- inherits(v, "cdid_cross_pair_cell_view")
    ks <- match(c(v$cohorts$g, v$cohorts$K), st$cohorts)
    if (anyNA(ks)) chained_reference_reader_abort("Recipe support cohort outside the common frame.", "integrity")
    # A cross-pair cell view declares its own per-cohort row universes (aligned
    # with c(g, K)); they need not equal the route defaults. Other views keep
    # the route's universes.
    cp_x <- cp_s <- cp_w <- NULL
    if (is_cp_cell) {
      cp_x <- v$cohort_x_rows; cp_s <- v$cohort_s_rows; cp_w <- v$cohort_w_rows
      for (md in list(cp_x, cp_s, cp_w)) {
        if (!is.null(md) && (!is.character(md) || length(md) != length(ks) || anyNA(md)))
          chained_reference_reader_abort("Cross-pair cell row universes must align with the view's support cohorts.",
            "integrity")
      }
    }
    s_ci <- t_ci <- 0L
    if (x_dated && (x_is_pair || is_cp_cell)) {
      s_ci <- match(as.double(v$pair$s), src$calendar)
      t_ci <- match(as.double(v$pair$t), src$calendar)
    }
    for (idx_k in seq_along(ks)) {
      k <- ks[idx_k]
      x_mode <- if (is.null(cp_x)) route$x_rows else cp_x[idx_k]
      s_mode <- if (is.null(cp_s)) route$s_rows else cp_s[idx_k]
      w_mode <- if (is.null(cp_w)) route$w_rows else cp_w[idx_k]
      if (x_dated && !is.null(x_mode) && !identical(x_mode, "none")) {
        x_pair <- identical(x_mode, "pair")
        xr <- if (x_pair && e != 0L) {
          # Restricted context: its own (pair, cohort, context) cache entry.
          ske <- paste(s_ci, t_ci, k, e)
          obs <- scan_e[[ske]]
          if (is.null(obs)) { obs <- scan_pair(k, e, v$pair); scan_e[[ske]] <- obs }
          obs
        } else if (x_pair) {
          if (is.null(scan)) scan <- if (use_int_scan) vector("list", n_cal * n_cal * J) else list()
          # Integer slot only for in-range calendar indices; otherwise the exact string key.
          sk <- if (use_int_scan && !is.na(s_ci) && !is.na(t_ci)) (((s_ci - 1L) * n_cal + (t_ci - 1L)) * J + k) else
            paste(s_ci, t_ci, k)
          obs <- scan[[sk]]
          if (is.null(obs)) {
            obs <- scan_pair(k, e, v$pair)
            scan[[sk]] <- obs
          }
          obs
        } else get_cohort_rows(k, x_mode, e)
        add("x", recipes$design_index[r] * N1 + xr, ((j * 4L + (if (x_pair) 0L else 1L)) * J1) + k)
      }
      if (s_dated && !is.null(s_mode) && !identical(s_mode, "none"))
        add("s", recipes$design_index[r] * N1 + get_cohort_rows(k, s_mode, e), ((j * 4L + 2L) * J1) + k)
      if (!is.null(w_mode) && !identical(w_mode, "none"))
        add("w", recipes$weight_index[r] * N1 + get_cohort_rows(k, w_mode, e), ((j * 4L + 3L) * J1) + k)
    }
  }
  ctx <- list(views = views, st = st, n_post = n_post, w_text = w_text)
  out <- list()
  for (role in c("x", "s", "w")) {
    nr <- cnt[[role]]
    if (nr == 0L) {
      out[[role]] <- list(keys = numeric(), labels = if (lazy_labels) integer() else character())
      next
    }
    k <- unlist(keys[[role]][seq_len(nr)], use.names = FALSE)
    l <- unlist(labs[[role]][seq_len(nr)], use.names = FALSE)
    if (!length(k)) {
      out[[role]] <- list(keys = numeric(), labels = if (lazy_labels) integer() else character())
      next
    }
    o <- order(k, method = "radix"); k <- k[o]; l <- l[o]; u <- !duplicated(k)
    ku <- k[u]; lu <- l[u]
    out[[role]] <- list(keys = ku, labels = if (lazy_labels) lu else chained_reference_format_labels(lu, ctx))
    # Issue #16 (missing_reference = "exclude"): every (key, demanding view)
    # pair of the design roles, to find the auxiliary reports a missing value
    # makes unavailable.
    if (keep_raw && role != "w") out[[role]]$raw <- list(keys = k, codes = l)
  }
  if (lazy_labels) attr(out, "label_ctx") <- ctx
  out
}

# Issue #16: link views (indices into `views`) that demand an unknown raw value
# of a dated design formula among their required (date, id) rows. Codes are the
# requirement labels ((j * 4 + p) * (J + 1) + k).
chained_reference_unknown_views <- function(src, req, J1) {
  hits <- integer()
  for (role in c("x", "s")) {
    f <- src$formulas[[role]]; raw <- req[[role]]$raw
    if (is.null(f) || !isTRUE(f$dated) || is.null(raw) || !length(raw$keys)) next
    vars <- intersect(all.vars(f$formula), src$prov$variable[src$prov$kind %in% c("dated", "unit_attribute")])
    if (!length(vars)) next
    miss <- chained_reference_rows_unknown(chained_reference_rows_frame(src, req[[role]]$keys, vars))
    if (!any(miss)) next
    codes <- raw$codes[raw$keys %in% req[[role]]$keys[miss]]
    hits <- c(hits, (as.integer(codes) %/% J1) %/% 4L)
  }
  sort(unique(as.integer(hits)))
}

chained_reference_missing <- function(src, keys, labels, variables, what, label_ctx = NULL,
                                      raw_covariate = FALSE) {
  if (is.integer(labels) && !is.null(label_ctx)) {
    labels <- chained_reference_format_labels(labels, label_ctx)
  }
  N1 <- src$N + 1
  d <- floor(keys / N1); i <- keys - d * N1
  date <- rep(NA_real_, length(d)); date[d > 0] <- src$calendar[d[d > 0]]
  tab <- data.frame(id = src$frame_ids[i], date = date, variable = variables, support = labels,
    stringsAsFactors = FALSE)
  show <- utils::head(tab, 5L)
  chained_reference_reader_abort(paste0("Missing required ", what, " (no carry-forward, substitution or deletion): ",
    paste0("id ", vapply(show$id, chained_reference_fmt, character(1)), ", date ",
      ifelse(is.na(show$date), "unit attribute", vapply(show$date, chained_reference_fmt, character(1))),
      ", variable `", show$variable, "`, ", show$support, collapse = "; "),
    if (nrow(tab) > 5L) paste0("; and ", nrow(tab) - 5L, " more") else "", "."),
    c(if (raw_covariate) "raw_covariate_missing", "missing"), diagnostics = list(missing = tab))
}

# One design slot: a static formula keeps the existing frame-position design of
# its unchanged static gate; a dated formula is frozen ONCE on the canonical
# distinct required rows of every role reading it (terms, factor levels,
# contrasts and transform state of that union), its required raw values are
# checked by original id, and its model matrix is built after the width is
# known and charged.
chained_reference_design_slot <- function(src, f, keys, labels, ctrl, retained, work, label, label_ctx = NULL) {
  if (!isTRUE(f$dated)) return(list(kind = "static", static = f$static, columns = colnames(f$static),
    formula = f$formula, n_rows = 0L))
  nd <- length(keys); vars <- intersect(all.vars(f$formula), src$prov$variable); nv <- length(vars)
  cp04_budget(ctrl, c(retained, 8 * nd), 16 * nd * (nv + 2), c(work, 64 * nd * (nv + 2)), 8 * nd,
    paste("reference canonical raw rows of", label))
  tf <- stats::terms(f$formula)
  df <- chained_reference_rows_frame(src, keys, vars)
  if (nv) {
    has_na <- FALSE
    for (j in seq_len(nv)) if (anyNA(df[[j]])) { has_na <- TRUE; break }
    if (has_na) {
      na_mat <- vapply(seq_len(nv), function(j) is.na(df[[j]]), logical(nd))
      if (!is.matrix(na_mat)) na_mat <- matrix(na_mat, nd, nv)
      hit <- which(na_mat, arr.ind = TRUE)
      hit <- hit[order(hit[, 1L], hit[, 2L]), , drop = FALSE]
      chained_reference_missing(src, keys[hit[, 1L]], labels[hit[, 1L]], names(df)[hit[, 2L]],
        "reference covariate values", label_ctx = label_ctx, raw_covariate = TRUE)
    }
  }
  # Freeze stateful transforms once, on this canonical union. A one-row slice
  # preserves matrix-valued transform columns and all factor levels, so its
  # model matrix gives the exact width before the full matrix exists.
  mf <- stats::model.frame(tf, df, na.action = stats::na.pass)
  for (j in seq_along(mf)) {
    if (is.character(mf[[j]])) mf[[j]] <- factor(mf[[j]])
    if (is.logical(mf[[j]])) mf[[j]] <- factor(mf[[j]], levels = c(FALSE, TRUE))
  }
  mf1 <- mf[1L, , drop = FALSE]
  attr(mf1, "terms") <- attr(mf, "terms")
  p0 <- ncol(stats::model.matrix(attr(mf, "terms"), mf1))
  store <- nd * (p0 + 2) + 64 * p0
  scratch <- nd * (8 * (nv + 1) + 4 * p0) + 64 * nd * p0 + 64 * p0 * p0 + unname(cp04_payload(list(df, mf))["entries"])
  cp04_budget(ctrl, c(retained, store), scratch, c(work, 128 * nd * p0 * p0 + 64 * nd * (1 + ceiling(log2(nd + 2)))),
    store, paste("reference frozen design of", label))
  Z <- chained_reference_model_matrix(attr(mf, "terms"), df, model_frame = mf)
  z_bad <- if (length(Z) > 0L) {
    if (anyNA(Z)) TRUE else { rg <- range(Z); !is.finite(rg[1L]) || !is.finite(rg[2L]) }
  } else FALSE
  if (z_bad) {
    bad <- which(rowSums(!is.finite(Z)) > 0L)
    chained_reference_missing(src, keys[bad], labels[bad], rep("(evaluated design column)", length(bad)),
      "finite evaluated design rows", label_ctx = label_ctx)
  }
  asg <- as.integer(attr(Z, "assign"))
  if (typeof(Z) == "double") {
    attributes(Z) <- list(dim = dim(Z), dimnames = list(NULL, colnames(Z)))
  } else {
    Z <- matrix(as.double(Z), nrow(Z), ncol(Z), dimnames = list(NULL, colnames(Z)))
  }
  list(kind = "dated", keys = keys, Z = Z, columns = colnames(Z),
    intercept = identical(as.integer(attr(tf, "intercept")), 1L), assign = asg, formula = f$formula,
    n_rows = nd,
    row_universe = paste("canonical distinct required (reference date, original id) rows of the post plus both",
      "auxiliary catalogs on this formula's required supports"),
    frozen_terms = "terms, factor levels, contrasts and transform state evaluated once on that union")
}

# Builds the sealed recipe-indexed reader of a dated route. All retained and
# scratch entries are charged before the union designs are allocated.
#
# Issue #16: with missing_reference = "exclude" every universe of a restricted
# recipe is read within its eligibility context, and an auxiliary report whose
# own earlier extraction date lacks a required raw value of an eligible unit is
# returned (class "cdid_reference_unavailable_views", nothing sealed) so the
# registry marks it unavailable through the existing auxiliary-support path and
# rebuilds the reader without it. Such values never change any context; a main
# (post-catalog) view can never demand one, and weights stay admission checks.
build_chained_reference_reader <- function(src, reader, views, n_post, ctrl, retained = 0, work = 0,
                                           missing_reference = "error",
                                           recover_missing_reference = FALSE) {
  if (!inherits(src, "cdid_reference_source")) chained_reference_reader_abort("Expected the request reference source.", "integrity")
  report_unknown <- identical(missing_reference, "exclude")
  st <- reader$state
  if (!identical(st$N, src$N) || !identical(as.vector(st$id), as.vector(src$frame_ids)))
    chained_reference_reader_abort("Reference source and full-roster reader disagree on original identities.", "integrity")
  rs <- chained_reference_recipes(src, views, n_post)
  route <- src$route; N1 <- src$N + 1
  xf <- src$formulas$x; sf <- src$formulas$s
  n_roles <- 1 + isTRUE(xf$dated) + (!is.null(sf) && isTRUE(sf$dated))
  # Requirement scan charge: one key and one label slot per (support, row, role).
  # The request's reference source is live during the scan as well: the
  # caller's pipeline envelope (`retained`) reserves the owners, the
  # four-column internal data, views and reports, never this source, so its
  # actual payload is charged here before the scan allocates anything.
  n_req <- 0
  for (v in views) for (k in match(c(v$cohorts$g, v$cohorts$K), st$cohorts))
    n_req <- n_req + n_roles * (st$cohort_end[k] - st$cohort_start[k] + 1)
  source_entries <- unname(cp04_payload(src)["entries"])
  cp04_budget(ctrl, c(retained, source_entries, 8 * n_req), 16 * n_req + 32 * ctrl$max_chunk_rows,
    c(work, 64 * n_req * (1 + ceiling(log2(n_req + 2)))), 8 * n_req, "reference requirement scan")
  if (report_unknown)
    cp04_budget(ctrl, c(retained, source_entries, 10 * n_req), 16 * n_req + 32 * ctrl$max_chunk_rows,
      c(work, 64 * n_req * (2 + ceiling(log2(n_req + 2)))), 10 * n_req, "reference eligibility requirement scan")
  req <- if (report_unknown) {
    chained_reference_requirements(src, reader, views, rs$view_recipe, rs$recipes, ctrl, n_post,
      lazy_labels = TRUE, keep_raw = TRUE)
  } else if (cdid_are_default_fns("chained_reference_requirements") &&
             length(formals(chained_reference_requirements)) >= 8L) {
    chained_reference_requirements(src, reader, views, rs$view_recipe, rs$recipes, ctrl, n_post, lazy_labels = TRUE)
  } else {
    chained_reference_requirements(src, reader, views, rs$view_recipe, rs$recipes, ctrl, n_post)
  }
  if (report_unknown) {
    bad <- chained_reference_unknown_views(src, req, length(st$cohorts) + 1L)
    if (any(bad <= n_post))
      chained_reference_reader_abort(paste("A main comparison demands an unknown reference value inside its",
        "eligibility context; the context and the requirement scan disagree."), "integrity")
    if (length(bad)) return(structure(list(views = bad, n_views = length(views), n_post = n_post),
      class = "cdid_reference_unavailable_views"))
    for (role in c("x", "s")) req[[role]]$raw <- NULL
  }
  label_ctx <- attr(req, "label_ctx")
  nw <- length(req$w$keys)
  base <- retained
  staged <- base + source_entries + 8 * n_req
  # One frozen state per distinct dated formula (identical outcome and
  # selection formulas share it on the union of both roles' rows).
  same <- !is.null(sf) && isTRUE(sf$dated) && isTRUE(xf$dated) &&
    identical(paste(deparse(sf$formula), collapse = " "), paste(deparse(xf$formula), collapse = " "))
  role_slot <- c(x = "x")
  if (!is.null(sf)) role_slot["s"] <- if (same) "x" else "s"
  designs <- list()
  xk <- req$x$keys; xl <- req$x$labels
  if (same) {
    kk <- c(xk, req$s$keys)
    o <- order(kk, method = "radix"); kk <- kk[o]; ll <- c(xl, req$s$labels)[o]
    u <- !duplicated(kk); xk <- kk[u]; xl <- ll[u]
  }
  design_slot <- function(f, keys, labels, label) {
    build <- function() chained_reference_design_slot(src, f, keys, labels, ctrl,
      staged, work, label, label_ctx = label_ctx)
    if (!recover_missing_reference) return(build())
    tryCatch(build(), cdid_reference_raw_covariate_missing_error = function(e) {
      # Reconstruct the demanding views only on a raw omission. The ordinary
      # reader, including valid incomplete panels, keeps its usual scan.
      cp04_budget(ctrl, c(staged, 2 * n_req), 16 * n_req + 32 * ctrl$max_chunk_rows,
        c(work, 64 * n_req * (2 + ceiling(log2(n_req + 2)))), 10 * n_req,
        "missing reference requirement scan")
      req_all <- chained_reference_requirements(src, reader, views, rs$view_recipe,
        rs$recipes, ctrl, n_post, lazy_labels = TRUE, keep_raw = TRUE)
      bad <- chained_reference_unknown_views(src, req_all, length(st$cohorts) + 1L)
      if (!length(bad)) stop(e)
      targets <- vapply(views[bad[bad <= n_post]], function(v) as.double(v$cohorts$g), numeric(1))
      affected <- unique(vapply(targets, calendar_reference_date, numeric(1),
        calendar = src$calendar, anticipation = src$anticipation))
      chained_reference_reader_abort(conditionMessage(e), "exclusion_needed",
        diagnostics = list(missing = e$diagnostics$missing, affected_dates = affected))
    })
  }
  designs$x <- design_slot(xf, xk, xl, "the outcome formula")
  staged <- staged + unname(cp04_payload(designs$x["Z"])["entries"])
  if (!is.null(sf) && !same) {
    designs$s <- design_slot(sf, req$s$keys, req$s$labels, "the selection formula")
    staged <- staged + unname(cp04_payload(designs$s["Z"])["entries"])
  }
  # Common outcome-regression alias map: the existing static record, or one map
  # on the dated outcome design's observed pair rows (fitting and prediction).
  or_index <- NULL; basis <- NULL; outcome_recipe <- NULL
  if (isTRUE(route$or_basis)) {
    if (identical(designs$x$kind, "dated")) {
      idx <- chained_reference_lookup(req$x$keys, designs$x$keys)
      Zb <- designs$x$Z[idx, , drop = FALSE]
      attr(Zb, "intercept") <- designs$x$intercept
      b <- or_public_basis(Zb, rep(TRUE, nrow(Zb)), NULL)
      basis <- b$record
      basis$row_universe <- "canonical distinct required (reference date, original id) observed pair rows of the post plus both auxiliary catalogs"
      basis$frozen_terms <- designs$x$frozen_terms
      or_index <- as.integer(basis$retained_index)
      outcome_recipe <- b$outcome_recipe
    } else {
      basis <- xf$basis
      or_index <- as.integer(basis$retained_index)
    }
  }
  cp04_budget(ctrl, c(staged, 2 * nw + 4096 + 64 * nrow(rs$recipes)), 8 * nw,
    c(work, 64 * nw * (1 + ceiling(log2(nw + 2)))), 2 * nw, "reference anchored weights")
  # Weights at their anchors (static weights keep the existing full-roster vector).
  w_val <- if (identical(src$weight$kind, "varying")) {
    src$weight$value[chained_reference_lookup(req$w$keys, src$weight$key)]
  } else src$weight$frame_w[req$w$keys - floor(req$w$keys / N1) * N1]
  miss_w <- is.na(w_val)
  if (any(miss_w)) chained_reference_missing(src, req$w$keys[miss_w], req$w$labels[miss_w],
    rep(if (is.null(src$weight$name)) "(unit weight)" else src$weight$name, sum(miss_w)),
    "sampling weights at the context anchor", label_ctx = label_ctx)
  # Recipe cohort masses and repeat rates a_k^{g,U} = H(W_g I_k U A) / H(W_g I_k U)
  # on the selection/imputation universes (each recipe's own law P^g).
  stats <- NULL
  if (route$w_rows %in% c("U", "F")) {
    J <- length(st$cohorts); R <- nrow(rs$recipes)
    mass <- rmass <- matrix(NA_real_, R, J)
    for (j in seq_along(views)) {
      r <- rs$view_recipe[j]; v <- views[[j]]
      for (k in match(c(v$cohorts$g, v$cohorts$K), st$cohorts)) {
        if (!is.na(mass[r, k])) next
        rows <- chained_reference_cohort_rows(st, k, route$w_rows)
        if (!is.null(src$eligibility))
          rows <- chained_reference_keep_eligible(src$eligibility, rs$recipes$eligibility_index[r], rows)
        wk <- w_val[chained_reference_lookup(rs$recipes$weight_index[r] * N1 + rows, req$w$keys)]
        mass[r, k] <- sum(wk) / src$N
        rmass[r, k] <- sum(wk[st$a_index[rows] > 0L]) / src$N
      }
    }
    stats <- list(mass_U = mass, repeat_mass = rmass, a = rmass / mass)
  }
  obj <- list(N = src$N, calendar = src$calendar, recipes = rs$recipes, view_recipe = rs$view_recipe,
    route = route, design_kind = src$design_kind, weight_kind = src$weight$kind, fix_weights = src$fix_weights,
    designs = designs, role_slot = role_slot, or_index = or_index, basis = basis, outcome_recipe = outcome_recipe,
    columns = if (is.null(or_index)) designs$x$columns else designs$x$columns[or_index],
    x_columns = designs$x$columns,
    weight_keys = req$w$keys, w = as.double(w_val), stats = stats, provenance = src$prov,
    counts = list(n_recipes = nrow(rs$recipes), n_views = length(views),
      n_required_design_rows = sum(vapply(designs, function(sl) as.double(sl$n_rows), numeric(1))),
      n_required_weight_rows = nw, n_requirement_slots = n_req),
    accounting = list(owned_entries = 0, source_entries = source_entries, scratch_entries_bound = 0,
      dense_unit_by_history_array = FALSE, carry_forward = FALSE))
  # Issue #16: restricted eligibility contexts (absent when unrestricted).
  if (!is.null(src$eligibility)) obj$eligibility <- src$eligibility
  # The request retains its actual-row reference source throughout fitting.
  obj$accounting$owned_entries <- unname(cp04_payload(obj)["entries"]) + source_entries
  cp04_budget(ctrl, c(base, obj$accounting$owned_entries), 0, work,
    obj$accounting$owned_entries, "retained reference reader")
  cdid_fast_seal(structure(obj, class = "cdid_reference_reader"), "reference_reader")
}

validate_chained_reference_reader <- function(ref) {
  if (!inherits(ref, "cdid_reference_reader") || !ipw_owner_seal_ok(ref, "reference_reader"))
    chained_reference_reader_abort("Expected the sealed recipe-indexed reference reader.", "integrity")
  invisible(ref)
}

# The reader bound to a backend (NULL on static paths).
chained_backend_reference <- function(backend) {
  if (is.null(backend) || is.null(backend$reference_ref)) return(NULL)
  if (!is.environment(backend$reference_ref) || !exists("reference_reader", backend$reference_ref, inherits = FALSE))
    chained_reference_reader_abort("Invalid sealed reference reader binding.", "integrity")
  ref <- get("reference_reader", backend$reference_ref, inherits = FALSE)
  validate_chained_reference_reader(ref)
  ref
}

# Sealed readers: rows are original frame positions; any request outside the
# recipe's canonical required support is refused, never filled.
chained_reference_x <- function(ref, recipe, rows, role = "x") {
  nm <- ref$role_slot[role]
  if (is.null(nm) || is.na(nm)) chained_reference_reader_abort(paste0("No `", role, "` design in this reference reader."), "integrity")
  sl <- ref$designs[[nm]]
  if (identical(sl$kind, "static")) return(sl$static[rows, , drop = FALSE])
  d <- ref$recipes$design_index[recipe]
  fast <- cdid_fast_sel_ref_slice(2L, rows, d * (ref$N + 1), ref$N, sl$keys, sl$Z)
  if (!is.null(fast)) return(fast)
  hit <- chained_reference_lookup(d * (ref$N + 1) + rows, sl$keys)
  if (anyNA(hit)) chained_reference_reader_abort("Design requested outside the recipe's required support.", "support")
  sl$Z[hit, , drop = FALSE]
}

chained_reference_s <- function(ref, recipe, rows) chained_reference_x(ref, recipe, rows, "s")

# Issue #16: TRUE for rows inside the recipe's eligibility context (every row
# of an unrestricted recipe or reader). Every same-pair route admitting
# restricted contexts (baseline, reconstruction and observation IPW) reads its
# rosters, supports, risk sets and IF rows through this mask.
chained_reference_eligible <- function(ref, recipe, rows) {
  el <- ref$eligibility
  if (is.null(el)) return(rep.int(TRUE, length(rows)))
  e <- ref$recipes$eligibility_index[recipe]
  if (length(e) != 1L || is.na(e)) chained_reference_reader_abort("Unknown recipe eligibility context.", "integrity")
  if (e == 0L) return(rep.int(TRUE, length(rows)))
  ci <- match(e, el$index)
  if (is.na(ci)) chained_reference_reader_abort("Unknown reference eligibility context.", "integrity")
  !(rows %in% el$excluded[[ci]])
}

chained_reference_restricted <- function(ref) !is.null(ref) && !is.null(ref$eligibility)

# Issue #16: TRUE when any of `recipe` (recipe indices of a sealed reader)
# reads a restricted eligibility context. Native kernels that read the full
# U/F roster of a recipe are bypassed exactly then; FALSE on unrestricted
# readers, recipes and static paths.
chained_reference_recipe_restricted <- function(ref, recipe) {
  if (!chained_reference_restricted(ref) || !length(recipe)) return(FALSE)
  e <- ref$recipes$eligibility_index[as.integer(recipe)]
  if (anyNA(e)) chained_reference_reader_abort("Unknown recipe eligibility context.", "integrity")
  any(e != 0L)
}

# Outcome-regression role: the common alias map's retained columns.
chained_reference_design <- function(ref, recipe, rows) {
  if (is.null(ref$or_index)) chained_reference_reader_abort("No outcome-regression basis in this reference reader.", "integrity")
  if (cdid_are_default_fns("chained_reference_x")) {
    nm <- ref$role_slot["x"]
    if (!is.null(nm) && !is.na(nm)) {
      sl <- ref$designs[[nm]]
      if (identical(sl$kind, "dated")) {
        d <- ref$recipes$design_index[recipe]
        fast <- cdid_fast_sel_ref_slice(3L, rows, d * (ref$N + 1), ref$N, sl$keys, sl$Z, as.integer(ref$or_index))
        if (!is.null(fast)) return(fast)
      }
    }
  }
  chained_reference_x(ref, recipe, rows, "x")[, ref$or_index, drop = FALSE]
}

chained_reference_weight <- function(ref, recipe, rows) {
  a <- ref$recipes$weight_index[recipe]
  fast <- cdid_fast_sel_ref_slice(1L, rows, a * (ref$N + 1), ref$N, ref$weight_keys, ref$w)
  if (!is.null(fast)) return(fast)
  hit <- chained_reference_lookup(a * (ref$N + 1) + rows, ref$weight_keys)
  if (anyNA(hit)) chained_reference_reader_abort("Weight requested outside the recipe's required support.", "support")
  ref$w[hit]
}

# Recipe cohort mass / repeat rate (vectorized over cohort indices k).
chained_reference_stat <- function(ref, recipe, k, what = "a") {
  if (is.null(ref$stats)) chained_reference_reader_abort("No recipe cohort rates in this reference reader.", "integrity")
  x <- ref$stats[[what]][cbind(rep(recipe, length(k)), k)]
  if (!length(x) || anyNA(x)) chained_reference_reader_abort("Recipe cohort rate requested outside the recipe's required support.", "support")
  x
}

# Full-F outcome recipe of a dated imputation/augmented outcome formula: every
# frozen column, the formula intercept and the remaining slope columns (no
# common alias reduction; the backend keeps its group-local preprocessing).
chained_reference_outcome_recipe <- function(ref) {
  sl <- ref$designs[[ref$role_slot["x"]]]
  cols <- sl$columns
  if (identical(sl$kind, "static")) {
    asg <- attr(sl$static, "assign"); intercept <- isTRUE(attr(sl$static, "intercept"))
  } else {
    asg <- sl$assign; intercept <- isTRUE(sl$intercept)
  }
  icpt <- if (intercept) which(as.integer(asg) == 0L) else integer(0)
  list(recipe = or_pooled_outcome_recipe(cols, intercept, cols[setdiff(seq_along(cols), icpt)]),
    columns = cols, intercept = intercept)
}

# Selection design specification of the reader's selection role (intercept,
# slope columns) on its canonical rows, or of the static selection design on U.
chained_reference_selection_spec <- function(ref, U, chunk_rows) {
  nm <- ref$role_slot["s"]
  if (is.null(nm) || is.na(nm)) chained_reference_reader_abort("No selection design in this reference reader.", "integrity")
  sl <- ref$designs[[nm]]
  if (identical(sl$kind, "static")) return(extract_selection_design_spec(sl$static, U, chunk_rows = chunk_rows))
  Z <- sl$Z
  attr(Z, "assign") <- sl$assign
  extract_selection_design_spec(Z, rep(TRUE, nrow(Z)), has_intercept = sl$intercept, chunk_rows = chunk_rows)
}

# Compact retained provenance record (no design rows).
chained_reference_record <- function(ref, src = NULL) {
  if (is.null(ref)) return(NULL)
  cal <- ref$calendar
  rates <- if (is.null(ref$stats)) NULL else {
    hit <- which(!is.na(ref$stats$a), arr.ind = TRUE)
    if (!length(hit)) NULL else data.frame(recipe = hit[, 1L], cohort_index = hit[, 2L],
      repeat_rate = ref$stats$a[hit], mass = ref$stats$mass_U[hit], stringsAsFactors = FALSE)
  }
  recipes <- data.frame(recipe = seq_len(nrow(ref$recipes)), key = ref$recipes$key,
    design_date = ifelse(ref$recipes$design_index == 0L, NA_real_, cal[pmax(1L, ref$recipes$design_index)]),
    weight_anchor = ifelse(ref$recipes$weight_index == 0L, NA_real_, cal[pmax(1L, ref$recipes$weight_index)]),
    stringsAsFactors = FALSE)
  # Issue #16: the reference date of each restricted recipe's eligibility context.
  if (!is.null(ref$eligibility)) recipes$eligibility_reference_date <- ifelse(ref$recipes$eligibility_index == 0L,
    NA_real_, cal[pmax(1L, ref$recipes$eligibility_index)])
  list(design_kind = ref$design_kind, weight_kind = ref$weight_kind, fix_weights = ref$fix_weights,
    route = ref$route$kind,
    recipes = recipes,
    view_recipe = ref$view_recipe, provenance = ref$provenance, columns = ref$columns, basis = ref$basis,
    formulas = lapply(ref$designs, function(sl) list(kind = sl$kind, formula = paste(deparse(sl$formula), collapse = " "),
      columns = sl$columns, n_required_rows = sl$n_rows)),
    role_slot = ref$role_slot, recipe_repeat_rates = rates,
    counts = ref$counts,
    conventions = c(post = "every post-catalog link of target g, including bridges and controls, extracts at b_g",
      varying = "auxiliary varying report (t_{j-1}, t_j) extracts at t_{j-1}",
      universal = "auxiliary universal report extracts at min(b_g, t_j)",
      weights = "W_g = w at t_1 (first_period) or at b_g (base_period) for every unit of target g's context",
      selection = "C^r holds both formulas' raw variables at the same extraction date; a NULL selection_xformla inherits the outcome formula"))
}

# Compact MP record: requested anchor law, the reduction applied, V = w_t1
# completeness by cohort and (reference route) the recipe table, provenance and
# frozen bases. Absent for static calls without an explicit fix_weights.
chained_request_reference_record <- function(prep, res) {
  if (identical(prep$weight_reduction, "static") && is.null(prep$fix_weights)) return(NULL)
  V <- prep$share_weights; g <- prep$frame_g
  coh <- sort(unique(g[!is.na(g)]))
  completeness <- data.frame(cohort = coh,
    n_units = vapply(coh, function(k) sum(g == k, na.rm = TRUE), numeric(1)),
    n_missing_first_period_weight = vapply(coh, function(k) sum(g == k & is.na(V), na.rm = TRUE), numeric(1)),
    missing_ids_head = vapply(coh, function(k) paste(utils::head(prep$frame_ids[which(g == k & is.na(V))], 5L),
      collapse = ", "), character(1)), stringsAsFactors = FALSE)
  b <- res$nuisance_registry$cp04
  recipes <- if (!is.null(b) && !is.null(b$reference_ref)) chained_reference_record(chained_backend_reference(b)) else NULL
  list(fix_weights = prep$fix_weights, weight_source = prep$weight_source, reduction = prep$weight_reduction,
    meaning = switch(prep$weight_reduction,
      static = paste("existing static path: weights and evaluated covariates pass the static gates, so NULL,",
        "first_period and base_period give identical results"),
      first_period_static = paste("unit-level covariates with first_period and w at the first calendar date on every",
        "roster unit: the existing static system evaluated with w_t1"),
      reference_recipes = paste("D60 recipe-indexed route: designs at target b_g (post) or the earlier report",
        "endpoint (auxiliary); W_g = w at t_1 or b_g for every unit of target g's context; every nuisance, rate,",
        "block and link is keyed by its recipe and fitted once")),
    route = prep$reference_route,
    share_law = if (is.null(prep$share_law)) "fixed_weights" else prep$share_law,
    first_period_weight_completeness = completeness, recipes = recipes)
}

# ------------------------------------------------------------------------------
# 4. Recipe-cohort layout of the full-roster imputation and augmented cores
# ------------------------------------------------------------------------------
# The accepted pure-imputation and augmented cores are written per cohort index
# over a reader state (weights, repeat rates, masses) and one full-F design. A
# dated context changes W_g and C^r with its recipe, so each cohort k used under
# recipe r becomes one layout cohort (k, r): its rows are the F rows of k with
# W_g at the recipe anchor and z(C^r) at the recipe extraction date, and its
# rate and mass are that recipe's a_k^{g,F} = H(W_g I_k A) / H(W_g I_k) and
# H(W_g I_k). The empirical normalizer stays the original N, outcome pairs are
# read at each row's original identity, pooling families are recipes (no
# cross-recipe pooling) and the counterfactual OR keeps ONE shared slope over
# every block. Layout rows of one original unit are summed into ONE influence
# row (all channels) before any covariance, cluster sum or GMM step.
chained_reference_layout <- function(reader, reference, views, ctrl, retained = 0, work = 0, selection = FALSE) {
  validate_cp04_reader(reader)
  validate_chained_reference_reader(reference)
  st <- reader$state
  if (length(reference$view_recipe) != length(views))
    chained_reference_reader_abort("The reference reader must index exactly the supplied link views.", "integrity")
  pairs <- lapply(seq_along(views), function(j) {
    v <- views[[j]]
    ks <- match(c(v$cohorts$g, sort(v$cohorts$K, method = "radix")), st$cohorts)
    if (anyNA(ks)) chained_reference_reader_abort("Recipe support cohort outside the common frame.", "integrity")
    cbind(k = ks, r = reference$view_recipe[j])
  })
  kr <- unique(do.call(rbind, pairs))
  kr <- kr[order(kr[, "k"], kr[, "r"], method = "radix"), , drop = FALSE]
  V <- nrow(kr)
  # Count excluded rows before allocating the layout. One tabulation per
  # reference date serves every recipe and cohort using that date.
  n_k <- st$cohort_end[kr[, "k"]] - st$cohort_start[kr[, "k"]] + 1
  el <- reference$eligibility
  e_by_block <- if (!is.null(el)) reference$recipes$eligibility_index[kr[, "r"]]
  for (j in seq_along(el$index)) {
    blocks <- which(e_by_block == el$index[j])
    if (!length(blocks)) next
    n_excluded <- tabulate(st$cohort_index[el$excluded[[j]]], length(st$cohorts))
    n_k[blocks] <- n_k[blocks] - n_excluded[kr[blocks, "k"]]
  }
  if (any(n_k < 1)) chained_reference_reader_abort("Empty recipe-cohort layout block inside an eligibility context.", "support")
  N_v <- sum(n_k)
  p_x <- length(reference$x_columns)
  s_slot <- if (selection) reference$designs[[reference$role_slot["s"]]] else NULL
  p_s <- if (selection) length(s_slot$columns) else 1L
  cp04_budget(ctrl, c(retained, N_v * (p_x + p_s + 10) + 64 * V), 16 * N_v + 8 * max(n_k) * (p_x + p_s),
    c(work, 64 * N_v * (p_x + p_s + 4) * (1 + ceiling(log2(N_v + 2)))), N_v * (p_x + p_s + 10),
    "reference recipe-cohort layout")
  ends <- cumsum(n_k); starts <- ends - n_k + 1
  orig <- integer(N_v); w <- numeric(N_v)
  X <- matrix(0, N_v, p_x, dimnames = list(NULL, reference$x_columns))
  Z <- if (selection) matrix(0, N_v, p_s, dimnames = list(NULL, s_slot$columns)) else
    matrix(1, N_v, 1L, dimnames = list(NULL, "reader_intercept"))
  for (v in seq_len(V)) {
    k <- kr[v, "k"]; r <- kr[v, "r"]
    rows <- st$cohort_order[seq.int(st$cohort_start[k], st$cohort_end[k])]
    if (!is.null(el)) rows <- chained_reference_keep_eligible(el, e_by_block[v], rows)
    ii <- seq.int(starts[v], ends[v])
    orig[ii] <- rows
    w[ii] <- chained_reference_weight(reference, r, rows)
    X[ii, ] <- chained_reference_x(reference, r, rows)
    if (selection) Z[ii, ] <- chained_reference_s(reference, r, rows)
  }
  if (selection) {
    if (identical(s_slot$kind, "dated")) {
      attr(Z, "has_intercept") <- isTRUE(s_slot$intercept)
      attr(Z, "assign") <- s_slot$assign
    } else {
      hi <- attr(s_slot$static, "intercept")
      if (!is.null(hi)) attr(Z, "has_intercept") <- isTRUE(hi)
      if (!is.null(attr(s_slot$static, "assign"))) attr(Z, "assign") <- as.integer(attr(s_slot$static, "assign"))
    }
  }
  a <- vapply(seq_len(V), function(v) chained_reference_stat(reference, kr[v, "r"], kr[v, "k"], "a"), numeric(1))
  mass <- vapply(seq_len(V), function(v) chained_reference_stat(reference, kr[v, "r"], kr[v, "k"], "mass_U"), numeric(1))
  cp04_finite(c(w, X, a, mass), "reference recipe-cohort layout")
  view_cohorts <- lapply(pairs, function(m) match(paste(m[, "k"], m[, "r"]), paste(kr[, "k"], kr[, "r"])))
  state <- list(id = st$id[orig], Z = Z, g = st$g[orig], w = w, U = rep(TRUE, N_v), rho = rep(1, N_v),
    universe = "F", justification = st$justification, mass_U = mass,
    a_index = st$a_index[orig], cohort_index = rep(seq_len(V), n_k), cohort_order = seq_len(N_v),
    cohort_start = starts, cohort_end = ends, cohorts = st$cohorts[kr[, "k"]],
    mass_F = mass, a = a, eta = rep(NA_real_, V), tau = rep(NA_real_, V), mean_w = NA_real_, mean_wA = NA_real_,
    N = st$N, n_A = st$n_A, lambda = st$lambda)
  # Original row -> layout rows (CSR) for summing each unit's rows.
  o <- order(orig, method = "radix")
  cnt <- tabulate(orig, st$N)
  lay <- list(state = state, X = X, view_cohorts = view_cohorts, families = as.integer(kr[, "r"]),
    cohort = as.integer(kr[, "k"]), recipe = as.integer(kr[, "r"]), orig_row = orig,
    by_orig = list(order = o, count = cnt, start = cumsum(cnt) - cnt + 1L), N = st$N, N_layout = N_v,
    frame_id = st$id)
  if (selection) {
    # The observation-probability models of the augmented route are fitted on a
    # sealed reader of this layout (same frame/panel owners and identities).
    rv <- structure(list(state = state, frame_ref = reader$frame_ref, panel_ref = reader$panel_ref,
      control = reader$control, accounting = list(owned_entries = 0, character_bytes = 0,
        borrowed_owner_entries = reader$accounting$borrowed_owner_entries, constructor_scratch_bound = 0,
        full_F_designs = 1L, full_F_masks_per_model = 0L, peak_rss_measured = FALSE,
        recipe_cohort_layout = TRUE)), class = "cdid_cp04_frame_reader")
    rv["ref"] <- list(NULL)
    rv$accounting$owned_entries <- unname(cp04_payload(rv)["entries"])
    lay$reader <- ipw_owner_seal(rv, "cp04_reader")
  }
  lay
}

# Layout rows of a chunk of original frame rows, and the summation map back.
chained_reference_layout_rows <- function(lay, rows) {
  cnt <- lay$by_orig$count[rows]
  idx <- rep(lay$by_orig$start[rows], cnt) + sequence(cnt) - 1L
  list(rows = lay$by_orig$order[idx], owner = rep(seq_along(rows), cnt))
}

# One row per original unit: the sum of its layout rows (zero when absent).
chained_reference_layout_sum <- function(phi_layout, owner, B) {
  out <- matrix(0, B, ncol(phi_layout), dimnames = list(NULL, colnames(phi_layout)))
  if (length(owner)) {
    s <- rowsum(phi_layout, owner, reorder = TRUE)
    out[as.integer(rownames(s)), ] <- s
  }
  out
}

# ------------------------------------------------------------------------------
# 5. Issue #16 registry support: auxiliary reports unavailable in their context
# ------------------------------------------------------------------------------

# The (pair, cohort, context) observed support counts of one link view inside
# its target's eligibility context, or NULL when that context is unrestricted.
chained_reference_view_context_support <- function(src, reader, v, ctrl) {
  el <- src$eligibility
  if (is.null(el)) return(NULL)
  st <- reader$state
  b <- calendar_reference_date(as.double(v$cohorts$g), src$calendar, src$anticipation)
  e <- chained_reference_eligibility_of(el, match(b, src$calendar))
  if (e == 0L) return(NULL)
  ks <- match(c(v$cohorts$g, v$cohorts$K), st$cohorts)
  if (anyNA(ks)) chained_reference_reader_abort("Recipe support cohort outside the common frame.", "integrity")
  vapply(ks, function(k) {
    a_rows <- chained_reference_keep_eligible(el, e, chained_reference_cohort_rows(st, k, "A"))
    n <- 0
    if (length(a_rows)) for (s0 in seq.int(1L, length(a_rows), by = max(1L, ctrl$max_chunk_rows))) {
      rr <- a_rows[seq.int(s0, min(length(a_rows), s0 + ctrl$max_chunk_rows - 1L))]
      n <- n + sum(cp04_pair_rows(reader, rr, v$pair)$S)
    }
    n
  }, numeric(1))
}

chained_reference_aux_condition <- function(message) {
  structure(list(message = message, call = NULL),
    class = c("cdid_reference_eligibility_support_error", "cdid_ipw_owner_pair_support_error", "cdid_ipw_pair_support_error",
      "cdid_ipw_support_error", "cdid_ipw_fit_error", "cdid_error", "error", "condition"))
}

# Drops auxiliary views (indices into `views`, all after the M catalog views)
# and records them as unavailable descriptors through the existing auxiliary
# support-error path; catalog views are never dropped.
chained_reference_drop_aux_views <- function(views, M, aux_index, aux_errors, aux, drop, reason_of) {
  if (!length(drop)) return(list(views = views, aux_index = aux_index, aux_errors = aux_errors, record = NULL))
  if (any(drop <= M)) chained_reference_reader_abort(
    "A main comparison lost its support inside its eligibility context; the restricted catalog and views disagree.",
    "integrity")
  rec <- list()
  for (vi in drop) {
    j <- match(vi, aux_index)
    if (is.na(j)) chained_reference_reader_abort("Auxiliary view without its descriptor.", "integrity")
    d <- aux$descriptors[[j]]; msg <- reason_of(vi)
    aux_errors[[j]] <- chained_reference_aux_condition(msg)
    rec[[length(rec) + 1L]] <- data.frame(descriptor = j, g = as.double(d$g), s = as.double(d$s), t = as.double(d$t),
      K = paste(fmt_design_num(d$K), collapse = ","), reason = msg, stringsAsFactors = FALSE)
  }
  keep <- setdiff(seq_along(views), drop)
  aux_index <- match(aux_index, keep)
  list(views = views[keep], aux_index = aux_index, aux_errors = aux_errors, record = do.call(rbind, rec))
}

# Registry step of a missing_reference = "exclude" fit: auxiliary views without
# eligible support in their target's context, then auxiliary views whose own
# extraction date lacks a required raw value, are marked unavailable; the
# sealed reader is built on the remaining views (the main catalog is intact).
chained_reference_registry_reader <- function(src, reader, views, M, aux_index, aux_errors, aux, ctrl,
                                              retained = 0, work = 0) {
  record <- NULL
  if (!is.null(src$eligibility) && length(views) > M) {
    # Catalog views already use restricted occupancy. Only the separately
    # requested pretrend comparisons can lack support inside the context.
    aux_views <- seq.int(M + 1L, length(views))
    sup <- lapply(views[aux_views], function(v) chained_reference_view_context_support(src, reader, v, ctrl))
    lost <- aux_views[vapply(sup, function(n) !is.null(n) && any(n < 1), logical(1))]
    dr <- chained_reference_drop_aux_views(views, M, aux_index, aux_errors, aux, lost, function(vi) {
      v <- views[[vi]]
      paste0('Unavailable in the restricted eligibility context of target ', fmt_design_num(v$cohorts$g),
        ' (missing_reference = "exclude"): no eligible unit of cohort(s) ',
        paste(fmt_design_num(c(v$cohorts$g, v$cohorts$K)[sup[[vi - M]] < 1]), collapse = ", "),
        " is observed at both dates of (", fmt_design_num(v$pair$s), ", ", fmt_design_num(v$pair$t), ")")
    })
    views <- dr$views; aux_index <- dr$aux_index; aux_errors <- dr$aux_errors; record <- dr$record
  }
  ref <- if (length(views)) build_chained_reference_reader(src, reader, views, M, ctrl, retained = retained, work = work,
    missing_reference = "exclude") else NULL
  if (inherits(ref, "cdid_reference_unavailable_views")) {
    bad <- ref$views
    dr <- chained_reference_drop_aux_views(views, M, aux_index, aux_errors, aux, bad, function(vi) {
      v <- views[[vi]]
      paste0("Unavailable: this auxiliary comparison reads covariates at its own earlier date ", fmt_design_num(v$pair$s),
        ', where a required raw value of an eligible unit is missing (missing_reference = "exclude"; ',
        "eligibility contexts and main comparisons are unchanged)")
    })
    views <- dr$views; aux_index <- dr$aux_index; aux_errors <- dr$aux_errors
    record <- rbind(record, dr$record)
    ref <- if (length(views)) build_chained_reference_reader(src, reader, views, M, ctrl, retained = retained,
      work = work, missing_reference = "exclude") else NULL
    if (inherits(ref, "cdid_reference_unavailable_views"))
      chained_reference_reader_abort("Auxiliary unavailability did not converge after one pass.", "integrity")
  }
  if (is.null(record)) record <- data.frame(descriptor = integer(), g = numeric(), s = numeric(), t = numeric(),
    K = character(), reason = character(), stringsAsFactors = FALSE)
  list(reference = ref, views = views, aux_index = aux_index, aux_errors = aux_errors, record = record)
}
