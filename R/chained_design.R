# Private CP02 augmented incidence algebra. Functional keys are caller-declared
# population-operator identities, NEVER certified from fitted rates or metadata.
# Default storage is O(vertices + links): two endpoints and one parent forest.

fmt_design_num <- function(x) sprintf("%.17g", x) # Display only, never identity.

design_vector <- function(x, n, type = "numeric", finite = TRUE) {
  is.atomic(x) && is.null(dim(x)) && length(x) == n && !is.object(x) &&
    (if (type == "numeric") is.numeric(x) else if (type == "logical") is.logical(x) else is.character(x)) &&
    !anyNA(x) && (type != "numeric" || !finite || all(is.finite(x)))
}

design_columns <- function(data, required, optional = character()) {
  if (!is.data.frame(data) || any(vapply(required, function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L))) ||
      any(vapply(intersect(optional, names(data)), function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L)))) {
    cdid_abort("Missing or ambiguous selected design column.", "cdid_design_input_error")
  }
}

design_bound <- function(x) {
  if (!design_vector(x, 1L) || x < 0 || x != floor(x) || x > 2^53 - 1) {
    cdid_abort("Supply an explicit finite integer entry bound.", "cdid_design_bound_error")
  }
  as.double(x)
}

# Exact real-number sign of (a + b) - c in {-1L, 0L, 1L} for finite doubles a, b
# (with b >= 0) and finite or +Inf double c. Refuses if float cutoff arithmetic
# overflows finite double precision or inputs are non-finite/inconclusive.
exact_cmp_sum <- function(a, b, c) {
  if (!design_vector(a, 1L) || !design_vector(b, 1L) || b < 0 ||
      !design_vector(c, 1L, finite = FALSE) || is.nan(c) || c == -Inf) {
    cdid_abort("Inconclusive or invalid non-finite cutoff comparison.",
      c("cdid_design_arithmetic_error", "cdid_design_qualification_error",
        "cdid_design_reference_error", "cdid_design_horizon_error", "cdid_design_input_error"))
  }
  a <- as.double(a); b <- as.double(b); c <- as.double(c)
  s <- a + b
  if (!is.finite(s) || (is.finite(c) && !is.finite(c - b))) {
    cdid_abort("Cutoff arithmetic overflowed finite double precision; refusing inconclusive cutoff certification.",
      c("cdid_design_arithmetic_error", "cdid_design_qualification_error",
        "cdid_design_reference_error", "cdid_design_horizon_error", "cdid_design_input_error"))
  }
  if (is.infinite(c)) return(-1L)
  if (s < c) return(-1L)
  if (s > c) return(1L)
  if (abs(a) >= abs(b)) {
    bb <- s - a
    err <- (a - (s - bb)) + (b - bb)
  } else {
    aa <- s - b
    err <- (b - (s - aa)) + (a - aa)
  }
  if (!is.finite(err)) {
    cdid_abort("Inconclusive floating-point cutoff comparison.",
      c("cdid_design_arithmetic_error", "cdid_design_qualification_error",
        "cdid_design_reference_error", "cdid_design_horizon_error", "cdid_design_input_error"))
  }
  if (err < 0) -1L else if (err > 0) 1L else 0L
}

validate_calendar_spec <- function(calendar, cohorts, anticipation = 0, pt_assumption = "post", pt_start = NULL) {
  if (!design_vector(calendar, length(calendar)) || length(calendar) < 2L ||
      anyDuplicated(as.double(calendar))) {
    cdid_abort("Supply an explicit full finite numeric `calendar` of at least two unique dates.",
      c("cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (!design_vector(cohorts, length(cohorts), finite = FALSE) || length(cohorts) < 1L ||
      any(is.nan(cohorts) | cohorts == -Inf) || anyDuplicated(as.double(cohorts)) ||
      !any(is.finite(cohorts))) {
    cdid_abort("Supply an explicit non-empty `cohorts` universe of unique finite cohorts and optional `Inf`.",
      c("cdid_design_cohort_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (!design_vector(anticipation, 1L) || anticipation < 0 || anticipation != floor(anticipation)) {
    cdid_abort("`anticipation` must be a finite non-negative whole numeric scalar.",
      c("cdid_design_anticipation_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (!design_vector(pt_assumption, 1L, "character") || !pt_assumption %in% c("post", "all")) {
    cdid_abort("`pt_assumption` must be explicitly \"post\" or \"all\".",
      c("cdid_design_pt_assumption_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  cal <- sort(as.double(calendar))
  coh <- sort(as.double(cohorts))
  if (identical(pt_assumption, "post")) {
    if (!is.null(pt_start)) {
      cdid_abort("`pt_start` must be `NULL` when `pt_assumption = \"post\"`.",
        c("cdid_design_pt_start_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    return(structure(
      list(
        calendar = cal,
        cohorts = coh,
        anticipation = as.double(anticipation),
        max_cohort = max(coh),
        has_never_treated = any(is.infinite(coh))
      ),
      class = "cdid_calendar_spec"
    ))
  }
  if (is.null(pt_start) || !design_vector(pt_start, 1L) || !(as.double(pt_start) %in% cal)) {
    cdid_abort("Under `pt_assumption = \"all\"`, `pt_start` must be an explicit finite scalar date belonging to the supplied full calendar (D33).",
      c("cdid_design_pt_start_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  structure(
    list(
      calendar = cal,
      cohorts = coh,
      anticipation = as.double(anticipation),
      max_cohort = max(coh),
      has_never_treated = any(is.infinite(coh)),
      pt_assumption = "all",
      pt_start = as.double(pt_start)
    ),
    class = "cdid_calendar_spec"
  )
}

calendar_reference_date <- function(g, calendar, anticipation = 0) {
  if (!design_vector(g, length(g)) || length(g) < 1L) {
    cdid_abort("Target cohort `g` must be finite numeric.",
      c("cdid_design_cohort_error", "cdid_design_reference_error", "cdid_design_input_error"))
  }
  if (!design_vector(calendar, length(calendar)) || length(calendar) < 1L ||
      anyDuplicated(as.double(calendar))) {
    cdid_abort("`calendar` must be a non-empty vector of unique finite numeric dates.",
      c("cdid_design_calendar_error", "cdid_design_reference_error", "cdid_design_input_error"))
  }
  if (!design_vector(anticipation, 1L) || anticipation < 0 || anticipation != floor(anticipation)) {
    cdid_abort("`anticipation` must be a finite non-negative whole numeric scalar.",
      c("cdid_design_anticipation_error", "cdid_design_reference_error", "cdid_design_input_error"))
  }
  cal <- sort(as.double(calendar))
  delta <- as.double(anticipation)
  ug <- sort(unique(as.double(g)))
  uref <- vapply(ug, function(gi) {
    pre <- vapply(cal, function(tm) exact_cmp_sum(tm, delta, gi) < 0L, logical(1L))
    if (!any(pre)) {
      cdid_abort("No admissible pre-anticipation reference period t < g - anticipation on the supplied full calendar.",
        c("cdid_design_reference_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    max(cal[pre])
  }, numeric(1L))
  uref[match(as.double(g), ug)]
}

eligible_control_cohorts <- function(g, u, cohorts, anticipation = 0) {
  if (!design_vector(g, 1L) || !design_vector(u, 1L)) {
    cdid_abort("`g` and `u` must be finite numeric scalars.",
      c("cdid_design_cohort_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (!design_vector(cohorts, length(cohorts), finite = FALSE) || length(cohorts) < 1L ||
      any(is.nan(cohorts) | cohorts == -Inf) || anyDuplicated(as.double(cohorts))) {
    cdid_abort("`cohorts` must be a non-empty vector of unique finite cohorts and optional `Inf`.",
      c("cdid_design_cohort_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (!design_vector(anticipation, 1L) || anticipation < 0 || anticipation != floor(anticipation)) {
    cdid_abort("`anticipation` must be a finite non-negative whole numeric scalar.",
      c("cdid_design_anticipation_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  coh <- sort(as.double(cohorts))
  g <- as.double(g); u <- as.double(u); delta <- as.double(anticipation)
  keep <- coh != g & vapply(coh, function(k) exact_cmp_sum(u, delta, k) < 0L, logical(1L))
  coh[keep]
}

parse_control_pool_string <- function(cp_str) {
  if (!is.character(cp_str) || length(cp_str) != 1L || is.na(cp_str) || !grepl("^\\{.*\\}$", cp_str)) {
    return(NULL)
  }
  inner <- substring(cp_str, 2L, nchar(cp_str) - 1L)
  parts <- strsplit(inner, ",", fixed = TRUE)[[1L]]
  vals <- suppressWarnings(as.double(parts))
  if (!length(parts) || grepl(",[[:space:]]*$", inner) || any(!nzchar(parts)) ||
      anyNA(vals) || any(is.nan(vals) | vals == -Inf) || anyDuplicated(vals)) {
    cdid_abort("Invalid `control_pool` cohort set syntax.",
      c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
  }
  vals
}

# Per-call canonical typed cache of serialized `control_pool` sets (CP08 scaling
# action 1). Each distinct string is parsed once, at its first use in row order,
# so every refusal and its order are those of the per-row parse; later rows of
# the same pool reuse its typed vector, and the pool's exact cutoff checks are
# memoized per (pool, date). Nothing persists beyond the calling function.
# `typed` (CP08 catalog handoff only) supplies the interned pools of the same
# invocation that rendered `cp`: `index` is each row's pool, `values` the typed
# cohort set the text was rendered from, `text` its rendering. The caller has
# already bound `text[index]` to the actual `cp` column; only the text parse is
# skipped, every cleanliness/universe/exclusion check runs on the typed values.
design_pool_cache <- function(cp, delta = NULL, typed = NULL) {
  if (is.null(typed)) {
    u <- unique(cp)
    index <- match(cp, u)
    vals <- vector("list", length(u))
    done <- logical(length(u))
  } else {
    u <- typed$text
    index <- typed$index
    vals <- typed$values
    done <- rep(TRUE, length(u))
  }
  clean <- new.env(hash = TRUE, parent = emptyenv())
  values <- function(i) {
    k <- index[i]
    if (!done[k]) {
      vals[k] <<- list(parse_control_pool_string(u[k]))
      done[k] <<- TRUE
    }
    vals[[k]]
  }
  # all(k > u + delta) over the row's pool, exactly as the per-row vapply.
  all_clean <- function(i, date) {
    key <- paste0(index[i], ":", sprintf("%a", as.double(date) + 0))
    hit <- clean[[key]]
    if (is.null(hit)) {
      hit <- all(vapply(values(i), function(k) exact_cmp_sum(date, delta, k) < 0L, logical(1L)))
      assign(key, hit, envir = clean)
    }
    hit
  }
  list(values = values, all_clean = all_clean, index = index, n_pools = length(u))
}

# Truthful metadata validation only. No key generation, equality inference,
# inferred control group, normalized-rate fingerprint, or boolean certification.
validate_functional_metadata <- function(links, check_consistency = TRUE, has_cal = FALSE) {
  n <- nrow(links)
  qual_cols <- c("s_qualified", "t_qualified", "s_coordinate_role", "t_coordinate_role",
                 "max_qualified_horizon", "uses_post_nuisance",
                 if (has_cal) c("uses_anticipation_nuisance", "anticipation"))
  has_qual_any <- has_cal || any(qual_cols %in% names(links))
  # Legacy inputs retain unknown columns as data; select pool metadata only
  # when explicit qualification metadata selects the new branch.
  optional <- c("control_group", "horizon", "profile_id", "composition_corrected",
                "post_shareable", if (has_qual_any) "control_pool", qual_cols)
  design_columns(links, "functional_key", optional)
  fk <- links$functional_key
  if (!design_vector(fk, n, "character") || any(!nzchar(fk))) {
    cdid_abort("Declare nonempty opaque population-functional keys explicitly.", "cdid_design_functional_error")
  }
  if (has_qual_any) {
    if (!all(c("s_qualified", "t_qualified") %in% names(links)) ||
        !design_vector(links$s_qualified, n, "logical") ||
        !design_vector(links$t_qualified, n, "logical")) {
      cdid_abort("Qualified endpoint declarations require non-missing scalar logical `s_qualified` and `t_qualified` columns.",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (any(c("s_coordinate_role", "t_coordinate_role") %in% names(links))) {
      if (!all(c("s_coordinate_role", "t_coordinate_role") %in% names(links)) ||
          !design_vector(links$s_coordinate_role, n, "character") || any(!nzchar(links$s_coordinate_role)) ||
          !design_vector(links$t_coordinate_role, n, "character") || any(!nzchar(links$t_coordinate_role))) {
        cdid_abort("Coordinate role declarations require non-empty scalar character `s_coordinate_role` and `t_coordinate_role` columns.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- links$max_qualified_horizon
      if (!is.atomic(mh) || !is.null(dim(mh)) || is.object(mh) || length(mh) != n ||
          !(is.numeric(mh) || (is.logical(mh) && all(is.na(mh)))) ||
          any(is.nan(mh) | mh == -Inf, na.rm = TRUE)) {
        cdid_abort("Invalid `max_qualified_horizon` metadata.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("uses_post_nuisance" %in% names(links) && !design_vector(links$uses_post_nuisance, n, "logical")) {
      cdid_abort("Invalid `uses_post_nuisance` logical metadata.",
        c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (has_cal && "uses_anticipation_nuisance" %in% names(links) &&
        !design_vector(links$uses_anticipation_nuisance, n, "logical")) {
      cdid_abort("Invalid `uses_anticipation_nuisance` logical metadata.",
        c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (has_cal && "anticipation" %in% names(links) &&
        (!design_vector(links$anticipation, n) || any(links$anticipation < 0))) {
      cdid_abort("Invalid `anticipation` metadata.",
        c("cdid_design_anticipation_error", "cdid_design_calendar_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
  }
  if (has_qual_any && "control_pool" %in% names(links)) {
    cp <- links$control_pool
    if (!design_vector(cp, n, "character") || any(!nzchar(cp) | cp == "{}")) {
      cdid_abort("Invalid `control_pool` metadata.",
        c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
  }
  for (nm in intersect(c("control_group", "profile_id"), names(links))) {
    v <- links[[nm]]
    if (!design_vector(v, n, "character") || any(!nzchar(v))) cdid_abort("Invalid functional metadata.", "cdid_design_functional_error")
  }
  if ("control_group" %in% names(links) && any(!links$control_group %in% c("nevertreated", "notyettreated"))) {
    cdid_abort("Invalid control-group metadata.", "cdid_design_functional_error")
  }
  for (nm in intersect(c("composition_corrected", "post_shareable"), names(links))) {
    if (!design_vector(links[[nm]], n, "logical")) cdid_abort("Invalid logical metadata.", "cdid_design_functional_error")
  }
  if ("horizon" %in% names(links) &&
      (!design_vector(links$horizon, n, finite = FALSE) || any(links$horizon == -Inf))) {
    cdid_abort("Invalid horizon metadata.", "cdid_design_horizon_error")
  }
  rates <- attr(links, "profile_rates", exact = TRUE)
  if (!is.null(rates)) {
    if (is.matrix(rates) && is.numeric(rates) && nrow(rates) == n && ncol(rates) > 0L) {
      good <- all(is.finite(rates) & rates > 0)
    } else if (is.list(rates) && length(rates) == n) {
      good <- all(vapply(rates, function(v) design_vector(v, length(v)) && length(v) > 0L && all(v > 0), logical(1L)))
    } else good <- FALSE
    if (!good) cdid_abort("Invalid rate metadata; rates do not certify operators.", "cdid_design_functional_error")
  }
  # Reject contradictory descriptive metadata; agreement does not prove equality.
  if (check_consistency && "g" %in% names(links) && n > 1L) {
    ord <- order(links$g, fk, method = "radix")
    a <- ord[-n]; b <- ord[-1L]
    same <- links$g[a] == links$g[b] & fk[a] == fk[b]
    for (nm in intersect(c("control_group", "horizon", "profile_id", "composition_corrected",
                           if (has_qual_any) "control_pool"), names(links))) {
      if (any(same & links[[nm]][a] != links[[nm]][b])) {
        cdid_abort("A declared key has contradictory descriptive metadata.",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- as.double(links$max_qualified_horizon)
      same_mh <- (is.na(mh[a]) & is.na(mh[b])) | (!is.na(mh[a]) & !is.na(mh[b]) & mh[a] == mh[b])
      if (any(same & !same_mh)) {
        cdid_abort("A declared key has contradictory `max_qualified_horizon` metadata.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
    # Optional per-key rate descriptions must be internally consistent. This
    # is a conservative metadata contract, not a test of operator equality or
    # inequality. Agreement creates no key; a correction flag cannot bypass
    # conflicting metadata. Raw nuisance rates that are not an operator's
    # profile should not be supplied as its per-key profile description.
    if (!is.null(rates)) {
      for (j in which(same)) {
        ra <- if (is.matrix(rates)) rates[a[j], ] else rates[[a[j]]]
        rb <- if (is.matrix(rates)) rates[b[j], ] else rates[[b[j]]]
        if (!identical(as.double(ra), as.double(rb))) cdid_abort(
          "A declared key has inconsistent rate-profile metadata; no operator equality is certified.",
          c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
  }
  invisible(TRUE)
}

validate_design_links <- function(links, calendar_spec = NULL) {
  design_links_validation(links, calendar_spec = calendar_spec)$links
}

# The one complete semantic validation of declared links (the body of the
# strict `validate_design_links` wrapper above). Returns the validated table
# together with the exact per-link cutoff signs it computed on the calendar
# path, so the graph core (`chained_design_core`) reuses them instead of
# repeating the per-row exact comparisons. `pools` is the CP08 catalog
# handoff's typed pool packet (see `design_pool_cache`); NULL parses the
# `control_pool` text exactly as before.
design_links_validation <- function(links, calendar_spec = NULL, pools = NULL) {
  if (is.null(calendar_spec)) calendar_spec <- attr(links, "calendar_spec", exact = TRUE)
  s_cmp <- t_cmp <- NULL
  has_cal <- !is.null(calendar_spec)
  if (has_cal) {
    if (!inherits(calendar_spec, "cdid_calendar_spec") || !is.list(calendar_spec)) {
      cdid_abort("Invalid `calendar_spec` object.", c("cdid_design_calendar_error", "cdid_design_input_error"))
    }
    cspec <- validate_calendar_spec(
      calendar_spec$calendar,
      calendar_spec$cohorts,
      calendar_spec$anticipation,
      pt_assumption = if (is.null(calendar_spec$pt_assumption)) "post" else calendar_spec$pt_assumption,
      pt_start = calendar_spec$pt_start
    )
    if (!identical(calendar_spec, cspec)) {
      cdid_abort("`calendar_spec` differs from its canonical validated representation.",
        c("cdid_design_calendar_error", "cdid_design_input_error"))
    }
    calendar_spec <- cspec
  }
  if (!is.data.frame(links)) cdid_abort("`links` must be a data.frame.", "cdid_design_input_error")
  if (!("ref" %in% names(links)) && "reference" %in% names(links)) {
    design_columns(links, "reference"); links$ref <- links$reference
  }
  if (!("functional_key" %in% names(links))) cdid_abort(
    "Explicit functional keys are required; rates, labels and booleans cannot certify them.",
    c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
  design_columns(links, c("g", "s", "t", "ref", "functional_key"))
  n <- nrow(links)
  if (n == 0L) cdid_abort("At least one declared link is required.", "cdid_design_input_error")
  for (nm in c("g", "s", "t", "ref")) {
    if (!design_vector(links[[nm]], n)) cdid_abort("Design dates/cohorts must be finite scalar numeric columns.", "cdid_design_input_error")
  }
  validate_functional_metadata(links, check_consistency = FALSE, has_cal = has_cal)
  if (any(links$s >= links$t)) cdid_abort("Each link requires s < t.", "cdid_design_endpoint_error")
  ord <- order(links$g, method = "radix")
  if (any(links$ref >= links$g) || (n > 1L && any(
      links$g[ord[-1L]] == links$g[ord[-n]] & links$ref[ord[-1L]] != links$ref[ord[-n]]))) {
    cdid_abort("Each cohort requires a unique reference preceding treatment.", "cdid_design_reference_error")
  }
  if (!has_cal) {
    has_qual <- "s_qualified" %in% names(links)
    post <- links$t >= links$g | links$s >= links$g
    pre <- (links$s < links$g & links$s != links$ref) | (links$t < links$g & links$t != links$ref)
    if (!has_qual) {
      if ("post_shareable" %in% names(links) && any(post & !links$post_shareable)) cdid_abort(
        "Post sharing requires the declared full CS5/no-anticipation/overlap population qualification.", "cdid_design_post_share_error")
      if ("horizon" %in% names(links) && any(links$t >= links$g & links$horizon < pmax(links$g, links$t))) {
        cdid_abort("Post control horizon must cover max(g,t).", "cdid_design_horizon_error")
      }
      if (any(pre & grepl("^(nyt:)?(h|horizon)=[0-9.+eE-]+$", links$functional_key))) cdid_abort(
        "A bare horizon label does not declare a complete population comparison operator.",
        c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
    } else {
      sq <- links$s_qualified; tq <- links$t_qualified
      if (any((links$s == links$ref & !sq) | (links$t == links$ref & !tq))) {
        cdid_abort("Reference endpoints must be declared qualified.",
          c("cdid_design_qualification_error", "cdid_design_reference_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
      if (any((links$s < links$g & links$s != links$ref & sq) | (links$t < links$g & links$t != links$ref & tq))) {
        cdid_abort("Early non-reference coordinates cannot be shared or qualified under PT-post.",
          c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
      if ("s_coordinate_role" %in% names(links)) {
        check_role <- function(u, g_v, ref_v, q_v, r_v) {
          ok_ref <- u == ref_v & r_v == "reference"
          ok_pre <- u < g_v & u != ref_v & !q_v & r_v %in% c("pre_nuisance", "pre")
          ok_att <- u >= g_v & q_v & r_v %in% c("post_att", "post")
          ok_pnu <- u >= g_v & !q_v & r_v == "post_nuisance"
          all(ok_ref | ok_pre | ok_att | ok_pnu)
        }
        if (!check_role(links$s, links$g, links$ref, sq, links$s_coordinate_role) ||
            !check_role(links$t, links$g, links$ref, tq, links$t_coordinate_role)) {
          cdid_abort("Endpoint coordinate role contradicts its declared qualification or date.",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
      }
      if ("uses_post_nuisance" %in% names(links)) {
        expected_upn <- (links$s >= links$g & !sq) | (links$t >= links$g & !tq)
        if (any(links$uses_post_nuisance != expected_upn)) {
          cdid_abort("`uses_post_nuisance` contradicts endpoint qualification.",
            c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
      }
      if ("post_shareable" %in% names(links)) {
        expected_ps <- (links$s < links$g | sq) & (links$t < links$g | tq)
        if (any(links$post_shareable != expected_ps)) {
          cdid_abort("`post_shareable` must be TRUE iff every post endpoint on the link is qualified.",
            c("cdid_design_post_share_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
      }
      if ("horizon" %in% names(links) && any(links$t >= links$g & links$horizon < pmax(links$g, links$t))) {
        cdid_abort("Post control horizon must cover max(g,t).", "cdid_design_horizon_error")
      }
      if ("max_qualified_horizon" %in% names(links)) {
        mh <- as.double(links$max_qualified_horizon)
        if (any(!is.na(mh) & mh < links$g) ||
            any(links$s >= links$g & sq & (is.na(mh) | mh < links$s)) ||
            any(links$t >= links$g & tq & (is.na(mh) | mh < links$t)) ||
            any(links$s >= links$g & !sq & !is.na(mh) & mh >= links$s) ||
            any(links$t >= links$g & !tq & !is.na(mh) & mh >= links$t)) {
          cdid_abort("`max_qualified_horizon` contradicts endpoint qualification or cohort treatment date.",
            c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
      }
      nuis <- pre | (links$s >= links$g & !sq) | (links$t >= links$g & !tq)
      if (any(nuis & grepl("^(nyt:)?(h|horizon)=[0-9.+eE-]+$", links$functional_key))) cdid_abort(
        "A bare horizon label does not declare a complete population comparison operator.",
        c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
    }
    if (has_qual && "control_pool" %in% names(links)) {
      cp <- links$control_pool
      mh_vec <- if ("max_qualified_horizon" %in% names(links)) as.double(links$max_qualified_horizon) else NULL
      pools <- design_pool_cache(cp)
      braced <- grepl("^\\{.*\\}$", cp)
      for (i in seq_len(n)) {
        if (braced[i]) {
          vals <- pools$values(i)
          if (any(vals <= links$g[i])) {
            cdid_abort("`control_pool` must contain unique control cohorts strictly after target g.",
              c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
          }
          if ("control_group" %in% names(links) &&
              all(is.infinite(vals)) != (links$control_group[i] == "nevertreated")) {
            cdid_abort("`control_pool` contradicts `control_group` metadata.",
              c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
          }
          if (!is.null(mh_vec) && !is.na(mh_vec[i]) && any(is.finite(vals) & vals <= mh_vec[i])) {
            cdid_abort("`control_pool` contradicts `max_qualified_horizon`.",
              c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
          }
          if ((links$s[i] >= links$g[i] && links$s_qualified[i] != all(vals > links$s[i])) ||
              (links$t[i] >= links$g[i] && links$t_qualified[i] != all(vals > links$t[i]))) {
            cdid_abort("Endpoint qualification contradicts declared `control_pool`.",
              c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
          }
        }
      }
    }
  } else {
    has_qual <- TRUE
    delta <- calendar_spec$anticipation
    if (any(!links$g %in% calendar_spec$cohorts)) {
      cdid_abort("Target cohort `g` must belong to the supplied cohort universe.",
        c("cdid_design_cohort_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    if (any(!links$s %in% calendar_spec$calendar) || any(!links$t %in% calendar_spec$calendar)) {
      cdid_abort("Link endpoints `s` and `t` must belong to the supplied full calendar.",
        c("cdid_design_calendar_error", "cdid_design_endpoint_error", "cdid_design_input_error"))
    }
    if (any(!links$ref %in% calendar_spec$calendar)) {
      cdid_abort("Reference date `ref` must belong to the supplied full calendar.",
        c("cdid_design_reference_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    if ("anticipation" %in% names(links) && any(links$anticipation != delta)) {
      cdid_abort("Link `anticipation` column contradicts `calendar_spec$anticipation`.",
        c("cdid_design_anticipation_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    exact_ref <- calendar_reference_date(links$g, calendar_spec$calendar, delta)
    if (any(links$ref != exact_ref)) {
      cdid_abort("Declared reference date does not match max{t in calendar : t < g - anticipation} on the full calendar.",
        c("cdid_design_reference_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
    sq <- links$s_qualified; tq <- links$t_qualified
    s_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$s[i], delta, links$g[i]), integer(1L))
    t_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$t[i], delta, links$g[i]), integer(1L))
    s_is_ref <- links$s == links$ref; t_is_ref <- links$t == links$ref
    s_is_pre <- s_cmp < 0L & !s_is_ref; t_is_pre <- t_cmp < 0L & !t_is_ref
    s_is_ant <- s_cmp >= 0L & links$s < links$g; t_is_ant <- t_cmp >= 0L & links$t < links$g
    s_is_post <- links$s >= links$g; t_is_post <- links$t >= links$g
    s_is_causal <- s_cmp >= 0L; t_is_causal <- t_cmp >= 0L
    if (any((s_is_ref & !sq) | (t_is_ref & !tq))) {
      cdid_abort("Reference endpoints must be declared qualified.",
        c("cdid_design_qualification_error", "cdid_design_reference_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (any((s_is_pre & sq) | (t_is_pre & tq))) {
      cdid_abort("Pre-anticipation non-reference coordinates (u < g - anticipation, u != ref) cannot be shared or qualified.",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    s_bar_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$s[i], delta, calendar_spec$max_cohort), integer(1L))
    t_bar_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$t[i], delta, calendar_spec$max_cohort), integer(1L))
    if (any((s_is_causal & sq & s_bar_cmp >= 0L) | (t_is_causal & tq & t_bar_cmp >= 0L))) {
      cdid_abort("Causal endpoint cannot be qualified outside the parallel-trends cohort boundary u + anticipation < max(cohorts).",
        c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if ("s_coordinate_role" %in% names(links)) {
      is_pt_all_cal <- identical(calendar_spec$pt_assumption, "all")
      check_cal_role <- function(u_v, is_ref, is_pre, is_ant, is_post, q_v, r_v) {
        ok_ref <- is_ref & q_v & r_v == "reference"
        ok_pre <- if (is_pt_all_cal) {
          is_pre & !q_v & (r_v %in% c("pre_nuisance", "pre") |
            (u_v >= calendar_spec$pt_start & r_v == "pre_pt_all_zero"))
        } else {
          is_pre & !q_v & r_v %in% c("pre_nuisance", "pre")
        }
        ok_ant_att <- is_ant & q_v & r_v %in% c("anticipation_att", "anticipation")
        ok_ant_nuis <- is_ant & !q_v & r_v == "anticipation_nuisance"
        ok_post_att <- is_post & q_v & r_v %in% c("post_att", "post")
        ok_post_nuis <- is_post & !q_v & r_v == "post_nuisance"
        all(ok_ref | ok_pre | ok_ant_att | ok_ant_nuis | ok_post_att | ok_post_nuis)
      }
      if (!check_cal_role(links$s, s_is_ref, s_is_pre, s_is_ant, s_is_post, sq, links$s_coordinate_role) ||
          !check_cal_role(links$t, t_is_ref, t_is_pre, t_is_ant, t_is_post, tq, links$t_coordinate_role)) {
        cdid_abort("Endpoint coordinate role contradicts its declared qualification or calendar region.",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("uses_post_nuisance" %in% names(links)) {
      expected_upn <- (s_is_post & !sq) | (t_is_post & !tq)
      if (any(links$uses_post_nuisance != expected_upn)) {
        cdid_abort("`uses_post_nuisance` contradicts endpoint qualification.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("uses_anticipation_nuisance" %in% names(links)) {
      expected_uan <- (s_is_ant & !sq) | (t_is_ant & !tq)
      if (any(links$uses_anticipation_nuisance != expected_uan)) {
        cdid_abort("`uses_anticipation_nuisance` contradicts endpoint qualification.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("post_shareable" %in% names(links)) {
      expected_ps <- (!s_is_post | sq) & (!t_is_post | tq)
      if (any(links$post_shareable != expected_ps)) {
        cdid_abort("`post_shareable` must be TRUE iff every post endpoint on the link is qualified.",
          c("cdid_design_post_share_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("horizon" %in% names(links) && any(t_is_post & links$horizon < pmax(links$g, links$t))) {
      cdid_abort("Post control horizon must cover max(g,t).", "cdid_design_horizon_error")
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- as.double(links$max_qualified_horizon)
      for (i in seq_len(n)) {
        if (!is.na(mh[i])) {
          if (is.finite(mh[i]) && exact_cmp_sum(mh[i], delta, links$g[i]) < 0L) {
            cdid_abort("`max_qualified_horizon` precedes the causal horizon g - anticipation.",
              c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
          }
          if ((is.infinite(mh[i]) && !calendar_spec$has_never_treated) ||
              (is.finite(mh[i]) && exact_cmp_sum(mh[i], delta, calendar_spec$max_cohort) >= 0L)) {
            cdid_abort("`max_qualified_horizon` violates the cohort universe boundary mh + anticipation < max(cohorts).",
              c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
          }
        }
      }
      if (any(s_is_causal & sq & (is.na(mh) | mh < links$s)) ||
          any(t_is_causal & tq & (is.na(mh) | mh < links$t)) ||
          any(s_is_post & !sq & !is.na(mh) & mh >= links$s) ||
          any(t_is_post & !tq & !is.na(mh) & mh >= links$t) ||
          any(s_is_ant & !sq & !is.na(mh) & mh >= links$s & !(t_is_post & tq)) ||
          any(t_is_ant & !tq & !is.na(mh) & mh >= links$t)) {
        cdid_abort("`max_qualified_horizon` contradicts endpoint qualification or causal horizon.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    nuis <- s_is_pre | t_is_pre | (s_is_causal & !sq) | (t_is_causal & !tq)
    if (any(nuis & grepl("^(nyt:)?(h|horizon)=[0-9.+eE-]+$", links$functional_key))) cdid_abort(
      "A bare horizon label does not declare a complete population comparison operator.",
      c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
    if ("control_pool" %in% names(links)) {
      cp <- links$control_pool
      mh_vec <- if ("max_qualified_horizon" %in% names(links)) as.double(links$max_qualified_horizon) else NULL
      # Typed per-call pool cache: one parse per distinct pool and one exact
      # cleanliness evaluation per (pool, date), in the per-row order below.
      # A handoff packet is admitted only when it is bound to this very column:
      # each row's typed pool must be the pool whose rendering is the row's text.
      if (!is.null(pools)) {
        if (!is.list(pools) || !is.list(pools$values) || !is.character(pools$text) ||
            length(pools$text) != length(pools$values) || !is.integer(pools$index) ||
            length(pools$index) != n || anyNA(pools$index) ||
            any(pools$index < 1L | pools$index > length(pools$values)) ||
            !identical(pools$text[pools$index], cp)) {
          cdid_abort("Typed pool handoff is not bound to the actual `control_pool` column; the design is not constructed.",
            c("cdid_design_internal_error", "cdid_design_qualification_error", "cdid_design_input_error"))
        }
        # Bind the actual typed values, not just their claimed text/index. The
        # canonical adapter emits plain doubles and round-trip %.17g strings;
        # validate each unique pool once within the original parsing envelope.
        valid_values <- vapply(pools$values, function(vals) {
          is.double(vals) && is.null(attributes(vals)) && length(vals) > 0L &&
            !anyNA(vals) && !any(vals == -Inf) && !anyDuplicated(vals)
        }, logical(1L))
        rendered <- if (all(valid_values)) vapply(pools$values, function(vals) {
          paste0("{", paste(sprintf("%.17g", vals), collapse = ","), "}")
        }, character(1L)) else NULL
        if (!identical(rendered, pools$text)) {
          cdid_abort("Typed pool values do not match their canonical `control_pool` rendering; the design is not constructed.",
            c("cdid_design_internal_error", "cdid_design_qualification_error", "cdid_design_input_error"))
        }
      }
      pools <- design_pool_cache(cp, delta, typed = pools)
      in_universe <- rep(NA, pools$n_pools)
      for (i in seq_len(n)) {
        vals <- pools$values(i)
        if (is.null(vals)) {
          cdid_abort("On the calendar path, `control_pool` must be an explicit `{...}` cohort set.",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
        k_pool <- pools$index[i]
        if (is.na(in_universe[k_pool])) in_universe[k_pool] <- !any(!vals %in% calendar_spec$cohorts)
        if (!in_universe[k_pool]) {
          cdid_abort("`control_pool` cohorts must belong to the supplied cohort universe.",
            c("cdid_design_cohort_error", "cdid_design_qualification_error", "cdid_design_calendar_error", "cdid_design_input_error"))
        }
        if (any(vals == links$g[i]) || !pools$all_clean(i, links$ref[i])) {
          cdid_abort("`control_pool` must exclude target g and be clean at the reference date (all k > ref + anticipation).",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
        if ("control_group" %in% names(links) &&
            all(is.infinite(vals)) != (links$control_group[i] == "nevertreated")) {
          cdid_abort("`control_pool` contradicts `control_group` metadata.",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
        }
        if (!is.null(mh_vec) && !is.na(mh_vec[i]) && is.finite(mh_vec[i]) && !pools$all_clean(i, mh_vec[i])) {
          cdid_abort("`control_pool` contradicts `max_qualified_horizon` under anticipation.",
            c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
        }
        s_clean <- pools$all_clean(i, links$s[i])
        t_clean <- pools$all_clean(i, links$t[i])
        if ((s_is_causal[i] && sq[i] && !s_clean) || (t_is_causal[i] && tq[i] && !t_clean)) {
          cdid_abort("Qualified causal endpoint requires `control_pool` subset of E_u(anticipation) (all k > u + anticipation).",
            c("cdid_design_qualification_error", "cdid_design_horizon_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
      }
    }
  }
  validate_functional_metadata(links, check_consistency = TRUE, has_cal = has_cal)
  # Injective integer tuple identities; original doubles are matched exactly.
  gi <- match(links$g, sort(unique(links$g)))
  dates <- sort(unique(c(links$s, links$t)))
  fi <- match(links$functional_key, sort(unique(links$functional_key), method = "radix"))
  if (has_qual) {
    ep_gi <- c(gi, gi); ep_fi <- c(fi, fi)
    ep_di <- match(c(links$s, links$t), dates)
    ep_qual <- c(links$s_qualified, links$t_qualified)
    ord_ep <- order(ep_gi, ep_fi, ep_di, method = "radix")
    a_ep <- ord_ep[-length(ord_ep)]; b_ep <- ord_ep[-1L]
    same_ep <- ep_gi[a_ep] == ep_gi[b_ep] & ep_fi[a_ep] == ep_fi[b_ep] & ep_di[a_ep] == ep_di[b_ep]
    if (any(same_ep & ep_qual[a_ep] != ep_qual[b_ep])) {
      cdid_abort("Contradictory endpoint qualification for repeated (g, functional_key, date).",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if ("s_coordinate_role" %in% names(links)) {
      ep_role <- c(links$s_coordinate_role, links$t_coordinate_role)
      if (any(same_ep & ep_role[a_ep] != ep_role[b_ep])) {
        cdid_abort("Contradictory endpoint coordinate role for repeated (g, functional_key, date).",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
  }
  keys <- paste(gi, match(links$s, dates), match(links$t, dates), fi, sep = ":")
  if (anyDuplicated(keys)) cdid_abort("Duplicate declared links.", "cdid_design_duplicate_link_error")
  if (has_cal) attr(links, "calendar_spec") <- calendar_spec
  list(links = links, calendar_spec = if (has_cal) calendar_spec else NULL, s_cmp = s_cmp, t_cmp = t_cmp)
}

adapt_calendar_design_links <- function(links, calendar, cohorts, anticipation = 0,
                                        pt_assumption = "post", pt_start = NULL) {
  adapt_calendar_design_links_validation(links, calendar, cohorts, anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start)$links
}

# Body of the strict `adapt_calendar_design_links` wrapper: the calendar
# adaptation followed by the one complete validation, returned as the
# validation packet the graph core consumes.
adapt_calendar_design_links_validation <- function(links, calendar, cohorts, anticipation = 0,
                                                   pt_assumption = "post", pt_start = NULL) {
  cal_spec <- validate_calendar_spec(calendar, cohorts, anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start)
  if (!is.data.frame(links)) cdid_abort("`links` must be a data.frame.", "cdid_design_input_error")
  if (!("ref" %in% names(links)) && "reference" %in% names(links)) {
    design_columns(links, "reference"); links$ref <- links$reference
  }
  if (!("functional_key" %in% names(links))) cdid_abort(
    "Explicit functional keys are required; rates, labels and booleans cannot certify them.",
    c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
  design_columns(links, c("g", "s", "t", "functional_key"), optional = "ref")
  n <- nrow(links)
  if (n == 0L) cdid_abort("At least one declared link is required.", "cdid_design_input_error")
  for (nm in intersect(c("g", "s", "t", "ref"), names(links))) {
    if (!design_vector(links[[nm]], n)) cdid_abort("Design dates/cohorts must be finite scalar numeric columns.", "cdid_design_input_error")
  }
  if (any(!links$g %in% cal_spec$cohorts)) {
    cdid_abort("Target cohort `g` must belong to the supplied cohort universe.",
      c("cdid_design_cohort_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  if (any(!links$s %in% cal_spec$calendar) || any(!links$t %in% cal_spec$calendar)) {
    cdid_abort("Link endpoints `s` and `t` must belong to the supplied full calendar.",
      c("cdid_design_calendar_error", "cdid_design_endpoint_error", "cdid_design_input_error"))
  }
  exact_ref <- calendar_reference_date(links$g, cal_spec$calendar, cal_spec$anticipation)
  if ("ref" %in% names(links)) {
    if (any(!links$ref %in% cal_spec$calendar) || any(links$ref != exact_ref)) {
      cdid_abort("Declared reference date does not match max{t in calendar : t < g - anticipation} on the full calendar.",
        c("cdid_design_reference_error", "cdid_design_calendar_error", "cdid_design_input_error"))
    }
  } else {
    links$ref <- exact_ref
  }
  has_sq <- "s_qualified" %in% names(links)
  has_tq <- "t_qualified" %in% names(links)
  if (!has_sq && !has_tq) {
    if (all(c("s_coordinate_role", "t_coordinate_role") %in% names(links)) &&
        design_vector(links$s_coordinate_role, n, "character") &&
        design_vector(links$t_coordinate_role, n, "character")) {
      valid_roles <- c("reference", "pre", "pre_nuisance",
                       if (identical(cal_spec$pt_assumption, "all")) "pre_pt_all_zero",
                       "anticipation", "anticipation_att",
                       "anticipation_nuisance", "post", "post_att", "post_nuisance")
      if (any(!links$s_coordinate_role %in% valid_roles) || any(!links$t_coordinate_role %in% valid_roles)) {
        cdid_abort("Invalid endpoint coordinate role name.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
      q_roles <- c("reference", "anticipation", "anticipation_att", "post", "post_att")
      links$s_qualified <- links$s_coordinate_role %in% q_roles
      links$t_qualified <- links$t_coordinate_role %in% q_roles
    } else if ("max_qualified_horizon" %in% names(links)) {
      mh <- links$max_qualified_horizon
      if (!is.atomic(mh) || !is.null(dim(mh)) || is.object(mh) || length(mh) != n ||
          !(is.numeric(mh) || (is.logical(mh) && all(is.na(mh)))) ||
          any(is.nan(mh) | mh == -Inf, na.rm = TRUE)) {
        cdid_abort("Invalid `max_qualified_horizon` metadata.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
      mh_num <- as.double(mh)
      s_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$s[i], cal_spec$anticipation, links$g[i]), integer(1L))
      t_cmp <- vapply(seq_len(n), function(i) exact_cmp_sum(links$t[i], cal_spec$anticipation, links$g[i]), integer(1L))
      links$s_qualified <- links$s == links$ref | (s_cmp >= 0L & !is.na(mh_num) & links$s <= mh_num)
      links$t_qualified <- links$t == links$ref | (t_cmp >= 0L & !is.na(mh_num) & links$t <= mh_num)
    } else if ("control_pool" %in% names(links)) {
      cdid_abort("Clean `control_pool` alone does not declare parallel-trends domain or endpoint qualification; supply `s_qualified`/`t_qualified`, coordinate roles, or `max_qualified_horizon`.",
        c("cdid_design_qualification_error", "cdid_design_horizon_error", "cdid_design_functional_error", "cdid_design_input_error"))
    } else {
      cdid_abort("Declare conditional endpoint qualification explicitly through `s_qualified`/`t_qualified`, coordinate roles, or `max_qualified_horizon`; dates alone do not qualify causal endpoints.",
        c("cdid_design_qualification_error", "cdid_design_horizon_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
  }
  attr(links, "calendar_spec") <- cal_spec
  design_links_validation(links, calendar_spec = cal_spec)
}

normalize_target_post_times <- function(target_post_times, cohorts, refs, calendar = NULL) {
  empty <- data.frame(g = numeric(), time = numeric())
  if (is.null(target_post_times)) return(empty)
  if (design_vector(target_post_times, length(target_post_times))) {
    grid <- expand.grid(time = sort(unique(target_post_times)), g = cohorts,
      KEEP.OUT.ATTRS = FALSE)[, c("g", "time")]
  } else if (is.data.frame(target_post_times)) {
    tc <- if ("time" %in% names(target_post_times)) "time" else "t"
    design_columns(target_post_times, c("g", tc))
    gv <- target_post_times$g; tv <- target_post_times[[tc]]
    if (!design_vector(gv, nrow(target_post_times)) || !design_vector(tv, nrow(target_post_times))) {
      cdid_abort("Targets require finite scalar numeric columns.", "cdid_design_target_error")
    }
    grid <- data.frame(g = as.double(gv), time = as.double(tv))
  } else cdid_abort("Invalid requested post targets.", "cdid_design_target_error")
  if (!is.null(calendar) && nrow(grid) > 0L && any(!grid$time %in% calendar)) {
    cdid_abort("Requested target post date must belong to the supplied full calendar.",
      c("cdid_design_target_error", "cdid_design_calendar_error", "cdid_design_input_error"))
  }
  ci <- match(grid$g, cohorts)
  if (anyNA(ci) || any(grid$time < grid$g | grid$time <= refs[ci])) cdid_abort("Invalid cohort or post target date.", "cdid_design_target_error")
  # Exact tuple deduplication, avoiding data-frame numeric string conversion.
  key <- paste(ci, match(grid$time, sort(unique(grid$time))), sep = ":")
  grid <- grid[!duplicated(key), , drop = FALSE]
  grid <- grid[order(grid$g, grid$time, method = "radix"), , drop = FALSE]
  rownames(grid) <- NULL
  grid
}

# Temporary CSR adjacency is linear size and discarded after BFS. Sorted
# neighbours and canonical link ranks give deterministic parents under row order.
incidence_forest <- function(nv, from, to, anchors, edge_order = seq_along(from)) {
  m <- length(from)
  if (2 * as.double(m) > .Machine$integer.max) cdid_abort("Graph exceeds supported adjacency index range.", "cdid_design_bound_error")
  anchors <- unique(as.integer(anchors))
  u <- c(from, to); v <- c(to, from)
  e <- rep(seq_len(m), 2L); sign <- c(rep(1L, m), rep(-1L, m))
  ord <- order(u, v, rep(edge_order, 2L), method = "radix")
  degree <- tabulate(u, nbins = nv)
  v <- v[ord]; e <- e[ord]; sign <- sign[ord]
  end <- cumsum(as.double(degree)); start <- end - degree + 1
  comp <- integer(nv); parent <- parent_edge <- integer(nv); parent_sign <- integer(nv)
  queue <- integer(nv); roots <- integer(nv); nc <- 0L
  for (root in c(anchors, seq_len(nv))) {
    if (comp[root] != 0L) next
    nc <- nc + 1L; roots[nc] <- root; comp[root] <- nc
    head <- tail <- 1L; queue[1L] <- root
    while (head <= tail) {
      a <- queue[head]; head <- head + 1L
      if (degree[a] == 0L) next
      for (k in seq.int(start[a], end[a])) {
        b <- v[k]
        if (comp[b] == 0L) {
          comp[b] <- nc; parent[b] <- a; parent_edge[b] <- e[k]; parent_sign[b] <- sign[k]
          tail <- tail + 1L; queue[tail] <- b
        }
      }
    }
  }
  if (length(anchors) > 1L && anyDuplicated(comp[anchors]) > 0L) {
    parent <- parent_edge <- parent_sign <- integer(nv)
    vis <- logical(nv)
    # Anchors grouped by component once (stable radix order keeps each group's
    # original anchor order), instead of one length(anchors) scan per component.
    ca <- comp[anchors]
    a_sorted <- anchors[order(ca, method = "radix")]
    a_count <- tabulate(ca, nbins = nc)
    a_end <- cumsum(a_count)
    for (cid in seq_len(nc)) {
      comp_anchors <- if (a_count[cid] > 0L) a_sorted[seq.int(a_end[cid] - a_count[cid] + 1L, a_end[cid])] else integer(0L)
      seeds <- if (length(comp_anchors) > 0L) comp_anchors else roots[cid]
      head <- 1L; tail <- length(seeds)
      queue[seq_len(tail)] <- seeds
      vis[seeds] <- TRUE
      while (head <= tail) {
        a <- queue[head]; head <- head + 1L
        if (degree[a] == 0L) next
        for (k in seq.int(start[a], end[a])) {
          b <- v[k]
          if (!vis[b]) {
            vis[b] <- TRUE
            parent[b] <- a; parent_edge[b] <- e[k]; parent_sign[b] <- sign[k]
            tail <- tail + 1L; queue[tail] <- b
          }
        }
      }
    }
  }
  list(component_id = comp, n_components = nc, roots = roots[seq_len(nc)],
       parent_vertex = parent, parent_link = parent_edge, parent_sign = parent_sign, degree = degree)
}

build_chained_design <- function(links, target_post_times = NULL) {
  cal_spec <- attr(links, "calendar_spec", exact = TRUE)
  chained_design_core(design_links_validation(links, calendar_spec = cal_spec), target_post_times)
}

# Graph core shared by every design constructor (the body of the strict
# `build_chained_design` wrapper above): takes the validation packet of one
# complete `design_links_validation` and builds the unchanged vertices, both
# forests, labels, link order and diagnostics. It validates nothing itself;
# each constructor performs exactly one validation before entering it.
chained_design_core <- function(validated, target_post_times = NULL) {
  links <- validated$links
  cal_spec <- attr(links, "calendar_spec", exact = TRUE)
  has_cal <- !is.null(cal_spec)
  is_pt_all <- has_cal && identical(cal_spec$pt_assumption, "all")
  has_qual <- has_cal || ("s_qualified" %in% names(links))
  g <- as.double(links$g); s <- as.double(links$s); t <- as.double(links$t)
  ref <- as.double(links$ref); fk <- links$functional_key; m <- length(g)
  cohorts <- sort(unique(g)); gi <- match(g, cohorts); refs <- ref[match(cohorts, g)]
  targets <- normalize_target_post_times(target_post_times, cohorts, refs,
    calendar = if (has_cal) cal_spec$calendar else NULL)
  fn <- sort(unique(fk), method = "radix"); fi <- match(fk, fn)
  dates <- if (has_cal) cal_spec$calendar else sort(unique(c(refs, s, t, targets$time)))
  sq <- if (has_qual) links$s_qualified else (s == ref | s >= g)
  tq <- if (has_qual) links$t_qualified else (t == ref | t >= g)
  vertex_key <- function(ci, time, f, qual) {
    # -1 is reference; 0 is shared qualified causal (anticipation or post); positive indices are private nuisance keys.
    kind <- ifelse(time == refs[ci], -1L, ifelse(qual, 0L, f))
    paste(ci, match(time, dates), kind, sep = ":")
  }
  sk <- vertex_key(gi, s, fi, sq); tk <- vertex_key(gi, t, fi, tq)
  ci <- c(seq_along(cohorts), gi, gi, match(targets$g, cohorts))
  time <- c(refs, s, t, targets$time)
  f <- c(rep(0L, length(cohorts)), fi, fi, rep(0L, nrow(targets)))
  qual <- c(rep(TRUE, length(cohorts)), sq, tq, rep(TRUE, nrow(targets)))
  key <- vertex_key(ci, time, f, qual)
  if (!has_cal) {
    kind <- ifelse(time == refs[ci], 1L, ifelse(time < cohorts[ci], 2L, ifelse(!qual, 3L, 4L)))
    fv <- ifelse(kind %in% c(2L, 3L), f, 0L)
    role_names <- c("reference", "pre", "post_nuisance", "post")
    nuis_kinds <- c(2L, 3L)
  } else {
    # Exact per-link cutoff signs of (s + anticipation) - g and (t + anticipation) - g,
    # computed once by the validation on these same columns.
    s_pre_ant <- validated$s_cmp < 0L
    t_pre_ant <- validated$t_cmp < 0L
    pre_ant <- c(rep(TRUE, length(cohorts)), s_pre_ant, t_pre_ant, rep(FALSE, nrow(targets)))
    kind <- ifelse(time == refs[ci], 1L,
      ifelse(pre_ant, 2L,
        ifelse(time < cohorts[ci], ifelse(!qual, 3L, 4L),
          ifelse(!qual, 5L, 6L))))
    fv <- ifelse(kind %in% c(2L, 3L, 5L), f, 0L)
    role_names <- c("reference", "pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post")
    nuis_kinds <- c(2L, 3L, 5L)
  }
  keep <- !duplicated(key)
  ord <- which(keep)[order(ci[keep], kind[keep], time[keep], fv[keep], method = "radix")]
  vertices <- data.frame(g = cohorts[ci[ord]], time = time[ord],
    functional_key = ifelse(kind[ord] %in% nuis_kinds, fn[pmax(fv[ord], 1L)], NA_character_),
    role = role_names[kind[ord]], vertex_index = seq_along(ord),
    stringsAsFactors = FALSE)
  vk <- key[ord]; nv <- nrow(vertices)
  is_anchor_v <- if (is_pt_all) {
    vertices$role == "reference" | (vertices$role == "pre" & vertices$time >= cal_spec$pt_start)
  } else {
    vertices$role == "reference"
  }
  pv <- which(!is_anchor_v); p <- length(pv)
  vertices$param_index <- NA_integer_; vertices$param_index[pv] <- seq_len(p)
  # Unique display labels are not used to construct/match graph vertices.
  label <- ifelse(vertices$role == "reference", paste0("ref(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ")"),
    ifelse(vertices$role %in% c("anticipation", "post"), paste0("ATT(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ")"),
      ifelse(vertices$role == "anticipation_nuisance",
        paste0("beta_anticipation_nuisance(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ",fn=", vertices$functional_key, ")"),
        ifelse(vertices$role == "post_nuisance",
          paste0("beta_post_nuisance(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ",fn=", vertices$functional_key, ")"),
          paste0("beta(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ",fn=", vertices$functional_key, ")")))))
  vertices$vertex_id <- make.unique(label)
  from <- match(sk, vk); to <- match(tk, vk)
  canonical <- order(gi, s, t, fi, method = "radix"); edge_order <- order(canonical)
  anchors <- which(is_anchor_v)
  forest <- incidence_forest(nv, from, to, anchors, edge_order)
  anchored <- unique(forest$component_id[anchors])
  vertices$degree <- forest$degree; vertices$is_isolated <- forest$degree == 0L
  vertices$component_id <- forest$component_id
  vertices$connected_to_ref <- vertices$component_id %in% anchored
  vertices$anchored_value <- ifelse(is_anchor_v, 0, NA_real_)
  parameters <- vertices[pv, , drop = FALSE]; rownames(parameters) <- NULL
  parameters$param_id <- parameters$vertex_id
  parameters$level_estimable <- parameters$connected_to_ref
  # One additional linear-size time graph exposes false collapsed connectivity.
  all_ci <- match(vertices$g, cohorts)
  time_key <- paste(all_ci, match(vertices$time, dates), sep = ":")
  ut <- unique(time_key); ti <- match(time_key, ut)
  naive <- incidence_forest(length(ut), ti[from], ti[to], unique(ti[anchors]), edge_order)
  parameters$time_graph_connected_to_ref <- naive$component_id[ti[pv]] %in% naive$component_id[ti[anchors]]
  links$link_index <- seq_len(m); links$link_id <- paste0("link_", edge_order)
  links$start_vertex_index <- from; links$end_vertex_index <- to
  links$start_param_index <- vertices$param_index[from]; links$end_param_index <- vertices$param_index[to]
  if (!has_cal) {
    links$link_type <- ifelse(t < g, "pre_only", ifelse(s == ref | s >= g, "post_only", "straddling"))
  } else {
    links$link_type <- ifelse(t_pre_ant, "pre_only",
      ifelse(s_pre_ant & s != ref, "straddling",
        ifelse(t < g, "anticipation_only",
          ifelse(s < g & s != ref, "anticipation_to_post", "post_only"))))
  }
  links$component_id <- forest$component_id[from]; rownames(links) <- NULL
  if (has_cal) attr(links, "calendar_spec") <- cal_spec
  post <- parameters[parameters$role == "post", , drop = FALSE]
  post_targets <- post[, c("g", "time", "param_id", "param_index", "degree", "is_isolated", "component_id", "connected_to_ref", "level_estimable", "time_graph_connected_to_ref"), drop = FALSE]
  post_targets$is_supported <- post_targets$degree > 0L
  post_targets$status <- ifelse(post_targets$level_estimable, "estimable",
    ifelse(post_targets$is_isolated, "isolated_unsupported", "disconnected_component"))
  rownames(post_targets) <- NULL
  unanchored <- setdiff(seq_len(forest$n_components), anchored)
  rank <- as.integer(p - length(unanchored))
  # Retain inexpensive per-cohort diagnostics without per-cohort dense graphs.
  nh <- length(cohorts)
  pci <- match(parameters$g, cohorts)
  vcount <- tabulate(all_ci, nh); ccount <- tabulate(all_ci[forest$roots], nh)
  pcount <- tabulate(pci, nh); mcount <- tabulate(gi, nh)
  nullity_coh <- tabulate(all_ci[forest$roots[unanchored]], nh)
  rank_coh <- pcount - nullity_coh
  pre_flag <- parameters$role == "pre"
  ant_nuis_flag <- parameters$role == "anticipation_nuisance"
  ant_flag <- parameters$role == "anticipation"
  post_nuis_flag <- parameters$role == "post_nuisance"
  post_flag <- parameters$role == "post"
  by_cohort <- data.frame(g = cohorts, ref = refs, n_links = mcount,
    n_vertices = vcount, n_parameters = pcount,
    n_pre_parameters = tabulate(pci[pre_flag], nh),
    n_post_parameters = tabulate(pci[post_flag], nh), n_components = ccount,
    rank = rank_coh, nullity = nullity_coh,
    cycle_rank = mcount - rank_coh,
    n_estimable_post = tabulate(pci[post_flag & parameters$level_estimable], nh),
    n_isolated_post = tabulate(pci[post_flag & parameters$is_isolated], nh),
    n_disconnected_post = tabulate(pci[post_flag & !parameters$level_estimable & !parameters$is_isolated], nh))
  if (has_qual) {
    by_cohort$n_post_nuisance_parameters <- tabulate(pci[post_nuis_flag], nh)
  }
  if (has_cal) {
    by_cohort$n_anticipation_nuisance_parameters <- tabulate(pci[ant_nuis_flag], nh)
    by_cohort$n_anticipation_parameters <- tabulate(pci[ant_flag], nh)
    by_cohort$n_estimable_anticipation <- tabulate(pci[ant_flag & parameters$level_estimable], nh)
  }
  if (is_pt_all) {
    pt_start_val <- cal_spec$pt_start
    cal_extra_coh <- vapply(seq_len(nh), function(h) {
      sum(cal_spec$calendar >= pt_start_val & cal_spec$calendar < refs[h])
    }, integer(1L))
    rep_extra_dates_coh <- vapply(seq_len(nh), function(h) {
      length(unique(vertices$time[all_ci == h & vertices$role == "pre" & vertices$time >= pt_start_val]))
    }, integer(1L))
    restr_pre_coh <- tabulate(all_ci[vertices$role == "pre" & is_anchor_v], nh)
    zero_links_coh <- tabulate(gi[is.na(links$start_param_index) & is.na(links$end_param_index)], nh)
    by_cohort$pt_start <- rep(pt_start_val, nh)
    by_cohort$n_extra_calendar_dates <- cal_extra_coh
    by_cohort$n_extra_represented_dates <- rep_extra_dates_coh
    by_cohort$n_restricted_pre_parameters <- restr_pre_coh
    by_cohort$n_zero_loading_links <- zero_links_coh
    by_cohort$no_extra_calendar_window <- cal_extra_coh == 0L
    by_cohort$no_extra_represented_restriction <- restr_pre_coh == 0L
    by_cohort$no_extra_restrictions <- restr_pre_coh == 0L
  }
  rank_diag <- list(n_links = m, n_vertices = nv, n_parameters = p,
    n_cohorts = length(cohorts), n_components = forest$n_components, rank = rank,
    nullity = length(unanchored), cycle_rank = m - rank, has_full_column_rank = length(unanchored) == 0L,
    unanchored_components = unanchored, n_post_targets = nrow(post_targets),
    n_estimable_post_targets = sum(post_targets$level_estimable),
    n_isolated_post_targets = sum(post_targets$is_isolated),
    n_disconnected_post_targets = sum(!post_targets$level_estimable & !post_targets$is_isolated),
    by_cohort = by_cohort)
  if (has_qual) {
    rank_diag$n_pre_parameters <- sum(pre_flag)
    rank_diag$n_post_nuisance_parameters <- sum(post_nuis_flag)
    rank_diag$n_post_parameters <- sum(post_flag)
  }
  if (has_cal) {
    rank_diag$n_anticipation_nuisance_parameters <- sum(ant_nuis_flag)
    rank_diag$n_anticipation_parameters <- sum(ant_flag)
    rank_diag$n_estimable_anticipation <- sum(ant_flag & parameters$level_estimable)
  }
  if (is_pt_all) {
    pv_post <- which(vertices$role != "reference")
    zero_link_idx <- which(is.na(links$start_param_index) & is.na(links$end_param_index))
    rank_diag$n_extra_calendar_dates <- sum(by_cohort$n_extra_calendar_dates)
    rank_diag$n_extra_represented_dates <- sum(by_cohort$n_extra_represented_dates)
    rank_diag$n_restricted_pre_parameters <- sum(by_cohort$n_restricted_pre_parameters)
    rank_diag$n_zero_loading_links <- length(zero_link_idx)
    rank_diag$no_extra_calendar_window <- all(by_cohort$no_extra_calendar_window)
    rank_diag$no_extra_represented_restriction <- all(by_cohort$no_extra_represented_restriction)
    rank_diag$no_extra_restrictions <- all(by_cohort$no_extra_restrictions)
    rank_diag$restriction_is_identity <- sum(by_cohort$n_restricted_pre_parameters) == 0L
  }
  incidence <- list(start = ifelse(is.na(links$start_param_index), 0L, links$start_param_index),
    end = ifelse(is.na(links$end_param_index), 0L, links$end_param_index),
    n_links = m, n_parameters = p, link_ids = links$link_id, param_ids = parameters$param_id)
  out <- list(incidence = incidence, forest = forest, parameters = parameters,
    vertices = vertices, links = links, post_targets = post_targets,
    rank_diagnostics = rank_diag)
  if (has_cal) out$calendar_spec <- cal_spec
  if (is_pt_all) {
    out$pt_all_spec <- list(
      pt_assumption = "all",
      pt_start = cal_spec$pt_start,
      unrestricted_n_parameters = length(pv_post),
      restricted_n_parameters = p,
      free_post_param_indices = match(pv, pv_post),
      restricted_post_param_indices = which(!pv_post %in% pv),
      zero_loading_link_indices = zero_link_idx,
      n_extra_calendar_dates = rank_diag$n_extra_calendar_dates,
      n_extra_represented_dates = rank_diag$n_extra_represented_dates,
      n_restricted_pre_parameters = rank_diag$n_restricted_pre_parameters,
      n_zero_loading_links = rank_diag$n_zero_loading_links,
      no_extra_calendar_window = rank_diag$no_extra_calendar_window,
      no_extra_represented_restriction = rank_diag$no_extra_represented_restriction,
      no_extra_restrictions = rank_diag$no_extra_restrictions,
      restriction_is_identity = rank_diag$restriction_is_identity,
      by_cohort = by_cohort[, c("g", "ref", "pt_start", "n_extra_calendar_dates",
                                "n_extra_represented_dates", "n_restricted_pre_parameters",
                                "n_pre_parameters", "n_zero_loading_links",
                                "no_extra_calendar_window", "no_extra_represented_restriction",
                                "no_extra_restrictions"), drop = FALSE]
    )
  }
  structure(out, class = "cdid_chained_design")
}

construct_chained_design <- build_chained_design

build_calendar_chained_design <- function(links, calendar, cohorts, anticipation = 0,
                                          target_post_times = NULL, pt_assumption = "post",
                                          pt_start = NULL) {
  # The adaptation ends in the complete validation; the already validated
  # packet enters the graph core directly (the former second pass over the
  # same adapted table could neither change it nor refuse).
  adapted <- adapt_calendar_design_links_validation(
    links, calendar = calendar, cohorts = cohorts, anticipation = anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start
  )
  chained_design_core(adapted, target_post_times = target_post_times)
}

construct_calendar_chained_design <- build_calendar_chained_design

# CP08 validated catalog-to-design handoff: the private design-construction
# core for the catalog adapter's normalized packet (`catalog_design_handoff`).
# The packet's links table is exactly the standalone bridge output and already
# carries `ref`, both qualification columns and the rendered pools, so the
# calendar adaptation has nothing to derive; its remaining checks (cohort and
# calendar membership, exact reference) are the calendar branch of the one
# complete validation performed here. The catalog's own canonical calendar,
# cohorts and anticipation form the spec once (with the requested PT
# assumption); the packet's typed pools stand in for re-parsing the pool text
# they were rendered from; the validation's exact cutoff signs feed the graph
# core. Output is identical to `build_calendar_chained_design` on the same links.
build_catalog_handoff_design <- function(handoff, target_post_times = NULL, pt_assumption = "post", pt_start = NULL) {
  if (!is.list(handoff) || is.object(handoff) ||
      !identical(names(handoff), c("links", "calendar", "cohorts", "anticipation", "pools")) ||
      !is.data.frame(handoff$links) || !is.list(handoff$pools)) {
    cdid_abort("Expected the catalog adapter's validated design handoff packet.",
      c("cdid_design_internal_error", "cdid_design_input_error"))
  }
  cal_spec <- validate_calendar_spec(handoff$calendar, handoff$cohorts, handoff$anticipation,
    pt_assumption = pt_assumption, pt_start = pt_start)
  links <- handoff$links
  attr(links, "calendar_spec") <- cal_spec
  validated <- design_links_validation(links, calendar_spec = cal_spec, pools = handoff$pools)
  chained_design_core(validated, target_post_times = target_post_times)
}

restrict_chained_design_pt_all <- function(design, pt_start) {
  if (!inherits(design, "cdid_chained_design") || is.null(design$calendar_spec)) {
    cdid_abort("`restrict_chained_design_pt_all` requires a calendar `cdid_chained_design` object.",
      c("cdid_design_calendar_error", "cdid_design_input_error"))
  }
  cspec <- validate_calendar_spec(
    design$calendar_spec$calendar,
    design$calendar_spec$cohorts,
    design$calendar_spec$anticipation,
    pt_assumption = "all",
    pt_start = pt_start
  )
  links_in <- design$links
  idx_pos <- match("link_index", names(links_in), nomatch = 0L)
  id_pos <- match("link_id", names(links_in), nomatch = 0L)
  drop_link_id <- idx_pos > 0L && id_pos == idx_pos + 1L
  derived_cols <- c("link_index", if (drop_link_id) "link_id",
                    "start_vertex_index", "end_vertex_index",
                    "start_param_index", "end_param_index", "link_type", "component_id")
  keep_cols <- setdiff(names(links_in), derived_cols)
  rates <- attr(links_in, "profile_rates", exact = TRUE)
  links_in <- links_in[, keep_cols, drop = FALSE]
  if ("s_coordinate_role" %in% names(links_in)) {
    s_pre_win <- links_in$s != links_in$ref &
      vapply(seq_len(nrow(links_in)), function(i) exact_cmp_sum(links_in$s[i], cspec$anticipation, links_in$g[i]) < 0L, logical(1L)) &
      !links_in$s_qualified & links_in$s >= cspec$pt_start
    links_in$s_coordinate_role[s_pre_win & links_in$s_coordinate_role %in% c("pre", "pre_nuisance")] <- "pre_pt_all_zero"
  }
  if ("t_coordinate_role" %in% names(links_in)) {
    t_pre_win <- links_in$t != links_in$ref &
      vapply(seq_len(nrow(links_in)), function(i) exact_cmp_sum(links_in$t[i], cspec$anticipation, links_in$g[i]) < 0L, logical(1L)) &
      !links_in$t_qualified & links_in$t >= cspec$pt_start
    links_in$t_coordinate_role[t_pre_win & links_in$t_coordinate_role %in% c("pre", "pre_nuisance")] <- "pre_pt_all_zero"
  }
  if (!is.null(rates)) attr(links_in, "profile_rates") <- rates
  attr(links_in, "calendar_spec") <- cspec
  target_df <- if (is.data.frame(design$post_targets) && nrow(design$post_targets) > 0L) {
    design$post_targets[, c("g", "time"), drop = FALSE]
  } else {
    NULL
  }
  build_chained_design(links_in, target_post_times = target_df)
}

incidence_multiply <- function(design, x, transpose = FALSE) {
  inc <- design$incidence
  n <- if (transpose) inc$n_links else inc$n_parameters
  vector_input <- is.null(dim(x))
  if (!is.numeric(x) || any(!is.finite(x)) ||
      (vector_input && (length(x) != n || is.object(x))) ||
      (!vector_input && (!is.matrix(x) || nrow(x) != n))) cdid_abort("Invalid incidence product input.", "cdid_design_input_error")
  if (vector_input) x <- matrix(x, ncol = 1L)
  if (!transpose) {
    grounded <- rbind(rep(0, ncol(x)), x)
    out <- grounded[inc$end + 1L, , drop = FALSE] - grounded[inc$start + 1L, , drop = FALSE]
    rownames(out) <- inc$link_ids
  } else {
    out <- matrix(0, inc$n_parameters, ncol(x))
    # rowsum aggregates each endpoint without allocating an m-by-p object.
    for (side in c("start", "end")) {
      j <- inc[[side]]; keep <- j > 0L
      if (any(keep)) {
        sums <- rowsum(x[keep, , drop = FALSE], j[keep], reorder = FALSE)
        dest <- as.integer(rownames(sums))
        out[dest, ] <- out[dest, , drop = FALSE] + (if (side == "start") -1 else 1) * sums
      }
    }
    rownames(out) <- inc$param_ids
  }
  if (any(!is.finite(out))) cdid_abort("Incidence product overflow.", "cdid_design_product_error")
  colnames(out) <- colnames(x)
  if (vector_input) drop(out) else out
}

resolve_design_parameter <- function(design, x) {
  if (is.character(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x)) {
    j <- match(x, design$parameters$param_id)
  } else if (design_vector(x, 1L) && x == floor(x) && x >= 1 && x <= nrow(design$parameters)) {
    j <- as.integer(x)
  } else j <- NA_integer_
  if (is.na(j)) cdid_abort("Unknown or noninteger parameter index.", "cdid_design_input_error")
  j
}

# Sparse signed root path, constructed only when requested. Length <= vertices.
root_path_certificate <- function(design, vertex) {
  f <- design$forest; k <- 0L; cur <- vertex
  while (f$parent_vertex[cur] != 0L) { k <- k + 1L; cur <- f$parent_vertex[cur] }
  edges <- signs <- integer(k); cur <- vertex
  for (j in seq_len(k)) {
    edges[j] <- f$parent_link[cur]; signs[j] <- f$parent_sign[cur]; cur <- f$parent_vertex[cur]
  }
  ord <- order(edges)
  data.frame(link_index = edges[ord], sign = signs[ord])
}

check_level_estimable <- function(design, param) {
  j <- resolve_design_parameter(design, param)
  ok <- design$parameters$level_estimable[j]
  list(estimable = ok, witness_links = if (ok) root_path_certificate(design, design$parameters$vertex_index[j]) else NULL)
}

check_contrast_estimable <- function(design, param_a, param_b) {
  if (!inherits(design, "cdid_chained_design")) cdid_abort("Expected chained design.", "cdid_design_input_error")
  a <- resolve_design_parameter(design, param_a); b <- resolve_design_parameter(design, param_b)
  same <- design$parameters$component_id[a] == design$parameters$component_id[b]
  both <- design$parameters$level_estimable[a] && design$parameters$level_estimable[b]
  ok <- same || both; witness <- NULL
  if (ok) {
    aa <- root_path_certificate(design, design$parameters$vertex_index[a])
    bb <- root_path_certificate(design, design$parameters$vertex_index[b])
    edges <- sort(unique(c(aa$link_index, bb$link_index)))
    sign <- integer(length(edges))
    sign[match(aa$link_index, edges)] <- -aa$sign
    sign[match(bb$link_index, edges)] <- sign[match(bb$link_index, edges)] + bb$sign
    witness <- data.frame(link_index = edges[sign != 0L], sign = sign[sign != 0L])
  }
  list(param_a = design$parameters$param_id[a], param_b = design$parameters$param_id[b],
    estimable = ok, same_component = same, both_level_estimable = both, witness_links = witness)
}

nullspace_certificate <- function(design, component) {
  if (!design_vector(component, 1L) || component != floor(component) ||
      !component %in% design$rank_diagnostics$unanchored_components) cdid_abort("Expected unanchored component index.", "cdid_design_input_error")
  list(component_id = component, param_index = which(design$parameters$component_id == component), value = 1L)
}

# Debug/test reference only: check the total dense entries BEFORE any allocation.
materialize_chained_design <- function(design, max_entries) {
  bound <- design_bound(max_entries)
  inc <- design$incidence; m <- as.double(inc$n_links); p <- as.double(inc$n_parameters)
  nv <- as.double(nrow(design$vertices)); q <- as.double(design$rank_diagnostics$nullity)
  terms <- c(m * p, m * p, m * nv, p * q)
  if (any(!is.finite(terms)) || any(terms > bound) || sum(terms) > bound) cdid_abort("Dense reference exceeds explicit entry bound.", "cdid_design_bound_error")
  W <- matrix(0L, m, p, dimnames = list(inc$link_ids, inc$param_ids))
  k <- which(inc$start > 0L); W[cbind(k, inc$start[k])] <- -1L
  k <- which(inc$end > 0L); W[cbind(k, inc$end[k])] <- 1L
  roots <- matrix(0L, m, nv)
  for (v in seq_len(nv)) {
    cert <- root_path_certificate(design, v)
    roots[cert$link_index, v] <- cert$sign
  }
  paths <- matrix(0L, m, p, dimnames = dimnames(W))
  k <- which(design$parameters$level_estimable)
  paths[, k] <- roots[, design$parameters$vertex_index[k], drop = FALSE]
  B <- matrix(0L, p, q, dimnames = list(inc$param_ids, if (q > 0) paste0("comp_", design$rank_diagnostics$unanchored_components) else character()))
  for (j in seq_len(q)) B[nullspace_certificate(design, design$rank_diagnostics$unanchored_components[j])$param_index, j] <- 1L
  out <- design; out$W <- W
  out$rank_diagnostics$root_paths <- roots; out$rank_diagnostics$path_witnesses <- paths
  out$rank_diagnostics$nullspace_basis <- B
  out
}

build_unaugmented_time_design <- function(links, target_post_times = NULL) {
  v <- validate_design_links(links)
  v <- v[, c("g", "s", "t", "ref"), drop = FALSE]
  v$functional_key <- "unaugmented_shared_time"
  dates <- sort(unique(c(v$s, v$t)))
  key <- paste(match(v$g, sort(unique(v$g))), match(v$s, dates), match(v$t, dates), sep = ":")
  out <- build_chained_design(v[!duplicated(key), , drop = FALSE], target_post_times)
  out$is_valid_augmented <- FALSE
  out
}
