# Prototype bounded general fixed-operator bridge-certificate search for CP02 (base R only).
# Generalizes tools/prototypes/chained_bridge_search.R to traverse all calendar edges,
# test recovery at post endpoints v using K intersect E_v, keep K for propagation,
# reinterpret immutable ancestor paths
# with the trimmed final pool, seed newly certified ATT(g, v) vertices with E_v, and apply
# subset dominance across distinct certified seed origins.
# Do NOT modify R/ sources or tools/prototypes/chained_bridge_search.R.

if (!exists("cdid_abort", mode = "function")) {
  cdid_abort <- function(message, class) {
    stop(structure(
      class = c(class, "cdid_error", "error", "condition"),
      list(message = message, call = NULL)
    ))
  }
}

bridge_scalar_num <- function(x, finite = TRUE) {
  is.atomic(x) && is.null(dim(x)) && !is.object(x) && is.numeric(x) &&
    length(x) == 1L && !is.na(x) && (!finite || is.finite(x))
}

bridge_scalar_char <- function(x) {
  is.atomic(x) && is.null(dim(x)) && !is.object(x) && is.character(x) &&
    length(x) == 1L && !is.na(x) && nzchar(x)
}

bridge_bound <- function(x, name = "max_states") {
  if (!bridge_scalar_num(x, finite = TRUE) || x != floor(x) ||
      x < 1 || x > 2^53 - 1) {
    cdid_abort(
      paste0("`", name, "` must be an explicit positive finite integer bound."),
      "cdid_bridge_bound_error"
    )
  }
  as.double(x)
}

pool_canonical_key <- function(mask) {
  idx <- which(mask)
  if (length(idx) == 0L) "" else paste(idx, collapse = ":")
}

pool_is_subset <- function(mask_sub, mask_sup) {
  !any(mask_sub & !mask_sup)
}

pool_is_strict_subset <- function(mask_sub, mask_sup) {
  !any(mask_sub & !mask_sup) && any(mask_sup & !mask_sub)
}

format_cohort_value <- function(x) {
  ifelse(is.infinite(x), "Inf", sprintf("%.17g", x))
}

format_cohort_set <- function(cohorts) {
  if (length(cohorts) == 0L) return("{}")
  paste0("{", paste(format_cohort_value(cohorts), collapse = ","), "}")
}

validate_cohort_vector <- function(cohorts, name = "cohorts") {
  if (!is.atomic(cohorts) || !is.null(dim(cohorts)) || is.object(cohorts) ||
      !is.numeric(cohorts) || length(cohorts) == 0L || anyNA(cohorts) ||
      any(is.nan(cohorts) | cohorts == -Inf) || anyDuplicated(cohorts)) {
    cdid_abort(
      paste0("`", name, "` must be a non-empty unique numeric vector without -Inf."),
      "cdid_bridge_cohort_error"
    )
  }
  sort(as.double(cohorts))
}

resolve_cohort_masses <- function(cohorts, cohort_masses) {
  if (!is.atomic(cohort_masses) || !is.null(dim(cohort_masses)) ||
      is.object(cohort_masses) || !is.numeric(cohort_masses) ||
      length(cohort_masses) == 0L || anyNA(cohort_masses) ||
      any(!is.finite(cohort_masses) | cohort_masses <= 0)) {
    cdid_abort(
      "`cohort_masses` must be a vector of positive finite relative cohort masses.",
      "cdid_bridge_mass_error"
    )
  }
  nms <- names(cohort_masses)
  if (!is.null(nms)) {
    if (any(!nzchar(nms)) || anyDuplicated(nms)) {
      cdid_abort("Named `cohort_masses` must have unique non-empty names.", "cdid_bridge_mass_error")
    }
    lbl <- format_cohort_value(cohorts)
    alt_lbl <- ifelse(is.infinite(cohorts), "Inf", as.character(cohorts))
    idx <- match(lbl, nms)
    miss <- is.na(idx)
    if (any(miss)) idx[miss] <- match(alt_lbl[miss], nms)
    if (anyNA(idx)) {
      cdid_abort("`cohort_masses` is missing entries for one or more cohorts.", "cdid_bridge_mass_error")
    }
    out <- as.double(cohort_masses[idx])
  } else {
    if (length(cohort_masses) != length(cohorts)) {
      cdid_abort(
        "Unnamed `cohort_masses` must have the same length as sorted `cohorts`.",
        "cdid_bridge_mass_error"
      )
    }
    out <- as.double(cohort_masses)
  }
  names(out) <- format_cohort_value(cohorts)
  out
}

# Exact rational ratios of the supplied binary64 masses, without rounded
# normalization, overflowing sums, or a tolerance.
bridge_binary_mass <- function(x) {
  hex <- strsplit(sprintf("%a", x), "p", fixed = TRUE)[[1L]]
  mantissa <- substring(hex[1L], 3L)
  pieces <- strsplit(mantissa, ".", fixed = TRUE)[[1L]]
  fraction <- if (length(pieces) == 2L) pieces[2L] else ""
  digits <- strsplit(paste0(pieces[1L], fraction), "", fixed = TRUE)[[1L]]
  val <- 0
  for (d in digits) val <- val * 16 + match(d, c(0:9, letters[1:6])) - 1L
  exponent <- as.integer(hex[2L]) - 4L * nchar(fraction)
  while (val %% 2 == 0) {
    val <- val / 2
    exponent <- exponent + 1L
  }
  c(significand = val, exponent = exponent)
}

bridge_exact_mass_ratios <- function(masses) {
  if (!length(masses) || any(!is.finite(masses) | masses <= 0)) {
    cdid_abort("Operator masses must be positive and finite.", "cdid_bridge_mass_error")
  }
  parts <- lapply(as.double(masses), bridge_binary_mass)
  base <- parts[[1L]]
  vapply(parts, function(p) {
    a <- p[1L]; b <- base[1L]
    while (b != 0) {
      if (b == 1) { a <- 1; break }
      remainder <- a %% b; a <- b; b <- remainder
    }
    paste0(sprintf("%.0f", p[1L] / a), "/", sprintf("%.0f", base[1L] / a),
           "*2^", p[2L] - base[2L])
  }, character(1L), USE.NAMES = FALSE)
}

bridge_carry_bits <- function(x) {
  carry <- 0
  for (i in seq_along(x)) {
    val <- x[i] + carry
    x[i] <- val %% 2
    carry <- floor(val / 2)
  }
  while (carry > 0) { x <- c(x, carry %% 2); carry <- floor(carry / 2) }
  nz <- which(x != 0)
  if (length(nz)) x[seq_len(max(nz))] else numeric(0L)
}

bridge_dyadic_product <- function(factors) {
  bits <- 1
  exponent <- 0L
  for (factor in factors) {
    p <- bridge_binary_mass(factor)
    val <- p[1L]; next_bits <- numeric(0L)
    while (val > 0) { next_bits <- c(next_bits, val %% 2); val <- floor(val / 2) }
    product <- numeric(length(bits) + length(next_bits) - 1L)
    for (i in which(bits == 1)) {
      at <- i + seq_along(next_bits) - 1L
      product[at] <- product[at] + next_bits
    }
    bits <- bridge_carry_bits(product)
    exponent <- exponent + p[2L]
  }
  list(bits = bits, exponent = exponent)
}

bridge_dyadic_compare_sum <- function(lhs, rhs_terms) {
  start <- min(lhs$exponent, vapply(rhs_terms, function(x) x$exponent, numeric(1L)))
  finish <- max(lhs$exponent + length(lhs$bits),
    vapply(rhs_terms, function(x) x$exponent + length(x$bits), numeric(1L)))
  rhs <- numeric(finish - start)
  for (term in rhs_terms) {
    at <- term$exponent - start + seq_along(term$bits)
    rhs[at] <- rhs[at] + term$bits
  }
  rhs <- bridge_carry_bits(rhs)
  left <- numeric(max(length(rhs), lhs$exponent - start + length(lhs$bits)))
  left[lhs$exponent - start + seq_along(lhs$bits)] <- lhs$bits
  length(rhs) <- length(left)
  rhs[is.na(rhs)] <- 0
  different <- which(left != rhs)
  if (!length(different)) 0L else if (left[max(different)] < rhs[max(different)]) -1L else 1L
}

bridge_exact_odds_leq <- function(target_mass, target_density, control_masses, control_density, bound) {
  lhs <- bridge_dyadic_product(c(target_mass, target_density))
  rhs_terms <- lapply(seq_along(control_masses), function(i)
    bridge_dyadic_product(c(bound, control_masses[i], control_density[i])))
  bridge_dyadic_compare_sum(lhs, rhs_terms) <= 0L
}

bridge_exact_probability_row <- function(row) {
  positive <- row[row > 0]
  if (!length(positive)) return(FALSE)
  bridge_dyadic_compare_sum(bridge_dyadic_product(1),
    lapply(positive, bridge_dyadic_product)) == 0L
}

declare_population_overlap_oracle <- function(cohorts,
                                              cell_density = NULL,
                                              cohort_masses = NULL,
                                              max_odds = Inf,
                                              predicate = NULL,
                                              monotone_declared = FALSE,
                                              from_empirical_counts = FALSE) {
  if (!identical(from_empirical_counts, FALSE)) {
    cdid_abort(
      "Population pooled-X overlap cannot be certified from finite empirical counts.",
      "cdid_bridge_overlap_error"
    )
  }
  input_cohorts <- cohorts
  cohorts <- validate_cohort_vector(cohorts, "cohorts")
  has_density <- !is.null(cell_density)
  has_pred <- !is.null(predicate)
  if (has_density == has_pred) {
    cdid_abort(
      "Supply exactly one of population `cell_density` or explicit `predicate`.",
      "cdid_bridge_overlap_error"
    )
  }
  if (has_density) {
    if (!is.matrix(cell_density) || !is.numeric(cell_density) ||
        nrow(cell_density) != length(cohorts) || ncol(cell_density) < 1L ||
        anyNA(cell_density) || any(!is.finite(cell_density) | cell_density < 0) ||
        any(rowSums(cell_density) != 1) ||
        any(!apply(cell_density, 1L, bridge_exact_probability_row))) {
      cdid_abort(
        "`cell_density` must give finite non-negative population atom probabilities; each row must sum exactly to one (no implicit normalization or tolerance).",
        "cdid_bridge_overlap_error"
      )
    }
    cell_density <- cell_density[match(cohorts, input_cohorts), , drop = FALSE]
    masses <- resolve_cohort_masses(cohorts, cohort_masses)
    if (!bridge_scalar_num(max_odds, finite = FALSE) || max_odds <= 0) {
      cdid_abort("`max_odds` must be a positive numeric scalar (or +Inf).", "cdid_bridge_overlap_error")
    }
    eval_fn <- function(g, control_pool) {
      if (length(control_pool) == 0L) return(FALSE)
      gi <- match(g, cohorts)
      ki <- match(control_pool, cohorts)
      if (is.na(gi) || anyNA(ki) || any(ki == gi)) {
        cdid_abort("Invalid cohort in overlap evaluation.", "cdid_bridge_overlap_error")
      }
      fg <- cell_density[gi, ]
      supp_g <- fg > 0
      if (!any(supp_g)) return(FALSE)
      if (any(colSums(cell_density[ki, supp_g, drop = FALSE] > 0) == 0)) return(FALSE)
      if (is.finite(max_odds)) {
        for (j in which(supp_g)) {
          positive <- ki[cell_density[ki, j] > 0]
          if (!bridge_exact_odds_leq(masses[gi], fg[j], masses[positive],
                                     cell_density[positive, j], max_odds)) return(FALSE)
        }
      }
      TRUE
    }
    kind <- "population_cell_density"
  } else {
    if (!identical(max_odds, Inf)) {
      cdid_abort("`max_odds` applies only to population atom probabilities; a custom predicate must state its own overlap condition.",
                 "cdid_bridge_overlap_error")
    }
    if (!is.function(predicate) || !identical(monotone_declared, TRUE)) {
      cdid_abort(
        "Custom overlap predicate requires a function and `monotone_declared = TRUE`.",
        "cdid_bridge_overlap_error"
      )
    }
    masses <- if (!is.null(cohort_masses)) resolve_cohort_masses(cohorts, cohort_masses) else NULL
    eval_fn <- function(g, control_pool) {
      if (length(control_pool) == 0L) return(FALSE)
      ans <- predicate(g, sort(as.double(control_pool)))
      if (!is.logical(ans) || !is.null(dim(ans)) || length(ans) != 1L || is.na(ans)) {
        cdid_abort(
          "Overlap oracle predicate must return a non-missing logical scalar.",
          "cdid_bridge_overlap_error"
        )
      }
      ans
    }
    kind <- "declared_monotone_predicate"
  }
  structure(
    list(
      cohorts = cohorts,
      kind = kind,
      cohort_masses = masses,
      max_odds = max_odds,
      evaluate = eval_fn,
      population_inputs = if (has_density) list(cell_density = cell_density, cohort_masses = masses) else NULL,
      population_status = "caller_supplied_population_assumption"
    ),
    class = "cdid_population_overlap_oracle"
  )
}

new_overlap_cache <- function() {
  env <- new.env(parent = emptyenv())
  env$masks <- list()
  env$vals <- logical(0L)
  env$by_key <- new.env(parent = emptyenv(), hash = TRUE)
  env$n_evaluations <- 0L
  env
}

evaluate_overlap_mask <- function(oracle, g, control_cohorts, mask, cache) {
  if (!any(mask)) return(FALSE)
  key <- pool_canonical_key(mask)
  if (exists(key, envir = cache$by_key, inherits = FALSE)) {
    return(get(key, envir = cache$by_key, inherits = FALSE))
  }
  val <- oracle$evaluate(g, control_cohorts[mask])
  if (!is.logical(val) || length(val) != 1L || is.na(val)) {
    cdid_abort("Overlap oracle returned an invalid value.", "cdid_bridge_overlap_error")
  }
  cache$n_evaluations <- cache$n_evaluations + 1L
  n_prev <- length(cache$vals)
  if (n_prev > 0L) {
    for (i in seq_len(n_prev)) {
      m_prev <- cache$masks[[i]]
      v_prev <- cache$vals[i]
      if (val && !v_prev && pool_is_subset(mask, m_prev)) {
        cdid_abort(
          "Overlap oracle violated upward closure: subset pool passed while superset failed.",
          c("cdid_bridge_overlap_monotonicity_error", "cdid_bridge_overlap_error")
        )
      }
      if (v_prev && !val && pool_is_subset(m_prev, mask)) {
        cdid_abort(
          "Overlap oracle violated upward closure: subset pool passed while superset failed.",
          c("cdid_bridge_overlap_monotonicity_error", "cdid_bridge_overlap_error")
        )
      }
    }
  }
  idx <- n_prev + 1L
  cache$masks[[idx]] <- mask
  cache$vals[idx] <- val
  assign(key, val, envir = cache$by_key)
  val
}

prime_pool_overlap_checks <- function(oracle, g, control_cohorts, p_mask, cache) {
  if (!any(p_mask)) return(FALSE)
  ok_full <- evaluate_overlap_mask(oracle, g, control_cohorts, p_mask, cache)
  members <- which(p_mask)
  for (j in members) {
    s_mask <- logical(length(control_cohorts))
    s_mask[j] <- TRUE
    evaluate_overlap_mask(oracle, g, control_cohorts, s_mask, cache)
  }
  ok_full
}

induce_pair_support_from_units <- function(data,
                                           idname = "id",
                                           tname = "time",
                                           yname = "y",
                                           gname = "g") {
  if (exists("construct_observed_pairs", mode = "function")) {
    idx <- construct_observed_pairs(data, idname = idname, tname = tname,
                                    yname = yname, gname = gname)
    mat <- materialize_observed_pairs(
      idx,
      max_pairs = max(1, idx$counts$n_pairs),
      chunk_size = max(1, idx$counts$n_pairs)
    )
    cells <- mat$cohort_pair_cells
    calendar <- idx$all_times
    cohorts <- idx$cohorts
    n_units <- idx$counts$n_units
  } else {
    if (!is.data.frame(data)) {
      cdid_abort("`data` must be a data.frame.", "cdid_bridge_input_error")
    }
    req <- c(idname, tname, yname, gname)
    if (!all(vapply(as.list(req), bridge_scalar_char, logical(1L))) ||
        length(unique(req)) != 4L || any(vapply(req, function(nm) sum(names(data) == nm), integer(1L)) != 1L)) {
      cdid_abort("Supply four distinct existing column names.", "cdid_bridge_input_error")
    }
    id <- data[[idname]]; tm <- data[[tname]]
    yv <- data[[yname]]; gv <- data[[gname]]
    scalar_column <- function(x) is.atomic(x) && is.null(dim(x)) && length(x) == nrow(data)
    if (is.factor(id)) id <- as.character(id)
    if (!scalar_column(id) || is.object(id) || !(is.numeric(id) || is.character(id)) ||
        anyNA(id) || (is.numeric(id) && any(!is.finite(id))) ||
        (is.character(id) && any(!nzchar(id)))) cdid_abort("Invalid unit IDs.", "cdid_bridge_input_error")
    if (!scalar_column(tm) || is.object(tm) || !is.numeric(tm) || anyNA(tm) || any(!is.finite(tm)))
      cdid_abort("Invalid observation dates.", "cdid_bridge_input_error")
    if (!scalar_column(gv) || is.object(gv) || !is.numeric(gv) || anyNA(gv) || any(gv == -Inf))
      cdid_abort("Invalid cohort identities.", "cdid_bridge_input_error")
    if (is.logical(yv) && all(is.na(yv))) yv <- as.numeric(yv)
    if (!scalar_column(yv) || is.object(yv) || !is.numeric(yv) || any(is.nan(yv) | is.infinite(yv)))
      cdid_abort("Invalid outcomes.", "cdid_bridge_input_error")
    tm <- as.double(tm); yv <- as.double(yv); gv <- as.double(gv)
    uid <- match(id, sort(unique(id), method = "radix"))
    ord <- order(uid, tm, method = "radix")
    uid <- uid[ord]; tm <- tm[ord]; yv <- yv[ord]; gv <- gv[ord]
    n <- length(uid)
    if (n > 1L) {
      same <- uid[-1L] == uid[-n]
      if (any(same & gv[-1L] != gv[-n])) cdid_abort("Within-unit cohort changes are invalid.", "cdid_bridge_input_error")
      if (any(same & tm[-1L] == tm[-n])) cdid_abort("Duplicate unit/date rows are invalid.", "cdid_bridge_input_error")
    }
    calendar <- sort(unique(tm))
    cohorts <- sort(unique(gv))
    n_units <- length(unique(uid))
    obs <- !is.na(yv)
    uid <- uid[obs]; tm <- tm[obs]; gv <- gv[obs]; yv <- yv[obs]
    sp_t <- split(tm, uid)
    sp_g <- vapply(split(gv, uid), function(v) v[1L], numeric(1L))
    sp_y <- split(yv, uid)
    rows <- list()
    k <- 0L
    for (u in seq_along(sp_t)) {
      tu <- sp_t[[u]]
      m <- length(tu)
      if (m >= 2L) {
        for (a in seq.int(1L, m - 1L)) {
          for (b in seq.int(a + 1L, m)) {
            if (!is.finite(sp_y[[u]][b] - sp_y[[u]][a])) cdid_abort("Outcome pair difference overflow.", "cdid_bridge_input_error")
            k <- k + 1L
            rows[[k]] <- c(sp_g[u], tu[a], tu[b])
          }
        }
      }
    }
    if (k == 0L) {
      cells <- data.frame(g = numeric(), s = numeric(), t = numeric())
    } else {
      raw <- do.call(rbind, rows)
      cells <- unique(data.frame(g = raw[, 1L], s = raw[, 2L], t = raw[, 3L]))
    }
  }
  if (nrow(cells) == 0L) {
    ps <- data.frame(s = numeric(), t = numeric())
    ps$supported_cohorts <- list()
  } else {
    ord <- order(cells$s, cells$t, cells$g, method = "radix")
    cells <- cells[ord, , drop = FALSE]
    pair_key <- paste(match(cells$s, calendar), match(cells$t, calendar), sep = ":")
    ukey <- unique(pair_key)
    first_idx <- match(ukey, pair_key)
    ps <- data.frame(
      s = as.double(cells$s[first_idx]),
      t = as.double(cells$t[first_idx]),
      stringsAsFactors = FALSE
    )
    ps$supported_cohorts <- lapply(ukey, function(k) {
      sort(unique(as.double(cells$g[pair_key == k])))
    })
  }
  list(
    pair_support = ps,
    calendar = as.double(calendar),
    cohorts = as.double(cohorts),
    n_units = as.integer(n_units),
    support_source = "observed_unit_patterns",
    population_support_certified = FALSE
  )
}

normalize_pair_support <- function(pair_support, calendar, cohorts) {
  if (is.list(pair_support) && !is.data.frame(pair_support) &&
      "pair_support" %in% names(pair_support)) {
    pair_support <- pair_support$pair_support
  }
  if (!is.data.frame(pair_support) ||
      !all(c("s", "t", "supported_cohorts") %in% names(pair_support))) {
    cdid_abort(
      "`pair_support` must be a data.frame with columns `s`, `t`, and list-column `supported_cohorts`.",
      "cdid_bridge_input_error"
    )
  }
  m <- nrow(pair_support)
  s <- pair_support$s; t <- pair_support$t
  sc <- pair_support$supported_cohorts
  if (!is.atomic(s) || !is.null(dim(s)) || !is.numeric(s) || anyNA(s) || any(!is.finite(s)) ||
      !is.atomic(t) || !is.null(dim(t)) || !is.numeric(t) || anyNA(t) || any(!is.finite(t)) ||
      !is.list(sc) || length(sc) != m) {
    cdid_abort("Malformed columns in `pair_support`.", "cdid_bridge_input_error")
  }
  s <- as.double(s); t <- as.double(t)
  if (any(!s %in% calendar) || any(!t %in% calendar)) {
    cdid_abort("All `pair_support` dates must belong to `calendar`.", "cdid_bridge_input_error")
  }
  if (any(s >= t)) {
    cdid_abort("Every pair in `pair_support` must satisfy s < t.", "cdid_bridge_input_error")
  }
  if (m > 0L) {
    pkey <- paste(match(s, calendar), match(t, calendar), sep = ":")
    if (anyDuplicated(pkey)) {
      cdid_abort("Duplicate (s, t) rows in `pair_support`.", "cdid_bridge_input_error")
    }
    norm_sc <- vector("list", m)
    for (i in seq_len(m)) {
      v <- sc[[i]]
      if (!is.atomic(v) || !is.null(dim(v)) || !(is.numeric(v) || length(v) == 0L) ||
          anyNA(v) || any(is.nan(v) | v == -Inf)) {
        cdid_abort("Invalid cohort vector in `supported_cohorts`.", "cdid_bridge_input_error")
      }
      v <- sort(unique(as.double(v)))
      if (any(!v %in% cohorts)) {
        cdid_abort("`supported_cohorts` contains a cohort not in `cohorts`.", "cdid_bridge_input_error")
      }
      norm_sc[[i]] <- v
    }
    ord <- order(s, t, method = "radix")
    out <- data.frame(s = s[ord], t = t[ord], stringsAsFactors = FALSE)
    out$supported_cohorts <- norm_sc[ord]
  } else {
    out <- data.frame(s = numeric(), t = numeric(), stringsAsFactors = FALSE)
    out$supported_cohorts <- list()
  }
  rownames(out) <- NULL
  out
}

validate_bridge_search_inputs <- function(g,
                                          calendar,
                                          cohorts,
                                          cohort_masses,
                                          fitting_universe,
                                          pair_support,
                                          target_post_times,
                                          overlap_oracle,
                                          max_states,
                                          reference = NULL,
                                          anticipation = 0) {
  if (!bridge_scalar_num(anticipation, finite = TRUE) || anticipation != 0) {
    cdid_abort(
      "Only zero anticipation (`anticipation = 0`) is validated in this initial certificate class.",
      "cdid_bridge_anticipation_error"
    )
  }
  if (!is.atomic(calendar) || !is.null(dim(calendar)) || is.object(calendar) ||
      !is.numeric(calendar) || length(calendar) < 2L || anyNA(calendar) ||
      any(!is.finite(calendar)) || anyDuplicated(calendar)) {
    cdid_abort(
      "`calendar` must be a unique finite numeric vector with at least two dates.",
      "cdid_bridge_calendar_error"
    )
  }
  calendar <- sort(as.double(calendar))
  if (any(diff(calendar) != 1)) {
    cdid_abort(
      "Initial certificate search validates consecutive unit-spaced calendars only; irregular dates are not claimed.",
      "cdid_bridge_calendar_error"
    )
  }
  cohorts <- validate_cohort_vector(cohorts, "cohorts")
  if (length(cohorts) < 2L) {
    cdid_abort("At least two cohorts are required.", "cdid_bridge_cohort_error")
  }
  finite_cohorts <- cohorts[is.finite(cohorts)]
  if (any(!finite_cohorts %in% calendar)) {
    cdid_abort("Finite cohorts must belong to `calendar`.", "cdid_bridge_cohort_error")
  }
  if (!bridge_scalar_num(g, finite = TRUE) || !as.double(g) %in% finite_cohorts) {
    cdid_abort("Target cohort `g` must be a finite cohort present in `cohorts`.", "cdid_bridge_cohort_error")
  }
  g <- as.double(g)
  unaffected <- calendar[calendar < g]
  if (length(unaffected) == 0L) {
    cdid_abort("Target cohort `g` has no unaffected pre-treatment date in `calendar`.", "cdid_bridge_reference_error")
  }
  expected_ref <- max(unaffected)
  if (is.null(reference)) {
    reference <- expected_ref
  } else {
    if (!bridge_scalar_num(reference, finite = TRUE) || as.double(reference) != expected_ref) {
      cdid_abort(
        "Reference `b` must equal the existing last unaffected date before `g`.",
        "cdid_bridge_reference_error"
      )
    }
    reference <- as.double(reference)
  }
  if (!bridge_scalar_char(fitting_universe)) {
    cdid_abort("`fitting_universe` must be a non-empty character scalar.", "cdid_bridge_mass_error")
  }
  masses <- resolve_cohort_masses(cohorts, cohort_masses)
  if (!inherits(overlap_oracle, "cdid_population_overlap_oracle")) {
    cdid_abort(
      "`overlap_oracle` must be created by `declare_population_overlap_oracle()`.",
      "cdid_bridge_overlap_error"
    )
  }
  if (!identical(as.double(overlap_oracle$cohorts), cohorts)) {
    cdid_abort("`overlap_oracle$cohorts` must match `cohorts`.", "cdid_bridge_overlap_error")
  }
  if (!is.null(overlap_oracle$cohort_masses) &&
      !identical(bridge_exact_mass_ratios(masses), bridge_exact_mass_ratios(overlap_oracle$cohort_masses))) {
    cdid_abort("Overlap and search must use the same declared relative cohort masses.", "cdid_bridge_mass_error")
  }
  max_states <- bridge_bound(max_states, "max_states")
  if (!is.atomic(target_post_times) || !is.null(dim(target_post_times)) ||
      is.object(target_post_times) || !is.numeric(target_post_times) ||
      length(target_post_times) == 0L || anyNA(target_post_times) ||
      any(!is.finite(target_post_times))) {
    cdid_abort(
      "`target_post_times` must be a non-empty finite numeric vector of post dates.",
      "cdid_bridge_target_error"
    )
  }
  target_post_times <- sort(unique(as.double(target_post_times)))
  if (any(!target_post_times %in% calendar) || any(target_post_times < g)) {
    cdid_abort(
      "Every requested target post date must belong to `calendar` and satisfy t >= g.",
      "cdid_bridge_target_error"
    )
  }
  pair_support <- normalize_pair_support(pair_support, calendar, cohorts)
  list(
    g = g,
    reference = reference,
    anticipation = 0,
    calendar = calendar,
    cohorts = cohorts,
    cohort_masses = masses,
    fitting_universe = fitting_universe,
    pair_support = pair_support,
    target_post_times = target_post_times,
    overlap_oracle = overlap_oracle,
    max_states = max_states
  )
}

make_operator_key <- function(g, reference, horizon = NULL, fitting_universe, cohorts, cohort_masses, control_pool) {
  control_pool <- sort(unique(as.double(control_pool)))
  idx <- match(control_pool, cohorts)
  if (!length(idx) || anyNA(idx)) cdid_abort("Invalid operator pool.", "cdid_bridge_cohort_error")
  # Operator L_b depends on (g, reference, universe, pool, control_mass_ratios).
  descriptor <- list(contract = "declared_fixed_mixture_v1", g = as.double(g),
                     reference = as.double(reference), universe = fitting_universe,
                     pool = control_pool,
                     control_mass_ratios = bridge_exact_mass_ratios(cohort_masses[idx]))
  paste0("declared_op:", paste(format(serialize(descriptor, NULL, version = 2L)), collapse = ""))
}

# Qualification helpers for the general fixed-operator class C^gen_g.
# A date u is endpoint-qualified for pool K iff u == reference or (u in H_g and K subset E_u).
endpoint_is_qualified <- function(g, reference, admissible_horizons, control_cohorts, pool_mask, dates) {
  vapply(as.double(dates), function(u) {
    if (u == reference) return(TRUE)
    if (!u %in% admissible_horizons) return(FALSE)
    !any(pool_mask & !(control_cohorts > u))
  }, logical(1L))
}

endpoint_coordinate_role <- function(g, reference, admissible_horizons, control_cohorts, pool_mask, dates) {
  qual <- endpoint_is_qualified(g, reference, admissible_horizons, control_cohorts, pool_mask, dates)
  d <- as.double(dates)
  ifelse(d == reference, "reference",
    ifelse(d < g, "pre_nuisance",
      ifelse(qual, "post_att", "post_nuisance")))
}

max_qualified_post_horizon <- function(admissible_horizons, control_cohorts, pool_mask) {
  if (!any(pool_mask) || length(admissible_horizons) == 0L) return(NA_real_)
  ok <- vapply(admissible_horizons, function(v) !any(pool_mask & !(control_cohorts > v)), logical(1L))
  if (!any(ok)) NA_real_ else max(admissible_horizons[ok])
}

# Reconstruct a single-operator segment from an accepted state ID in the
# append-only state table by reinterpreting its immutable ancestor path with
# the trimmed final pool `final_pool_mask` = K(sid) intersect E_v.
# Never mutates any ancestor record in `state_table`.
reconstruct_general_segment_from_state <- function(sid,
                                                   final_pool_mask,
                                                   state_table,
                                                   g,
                                                   reference,
                                                   fitting_universe,
                                                   cohorts,
                                                   cohort_masses,
                                                   control_cohorts,
                                                   admissible_horizons,
                                                   target_edges) {
  untrimmed_mask <- state_table$pool[[sid]]
  if (!any(final_pool_mask) || !pool_is_subset(final_pool_mask, untrimmed_mask)) {
    cdid_abort("Final trimmed pool must be a non-empty subset of the state's untrimmed pool.",
               "cdid_bridge_input_error")
  }
  control_pool <- control_cohorts[final_pool_mask]
  untrimmed_pool <- control_cohorts[untrimmed_mask]
  pool_key <- pool_canonical_key(final_pool_mask)
  untrimmed_key <- state_table$pool_key[sid]
  trimmed_flag <- any(untrimmed_mask & !final_pool_mask)
  max_qual_h <- max_qualified_post_horizon(admissible_horizons, control_cohorts, final_pool_mask)
  op_key <- make_operator_key(g, reference, max_qual_h, fitting_universe, cohorts, cohort_masses, control_pool)

  # Walk backward through immutable parent pointers to the seed root (parent == 0L).
  chain <- integer(0L)
  cur <- sid
  while (cur != 0L) {
    if (!pool_is_subset(final_pool_mask, state_table$pool[[cur]])) {
      cdid_abort("Ancestor pool containment invariant violated.", "cdid_bridge_input_error")
    }
    chain <- c(cur, chain)
    cur <- state_table$parent[cur]
  }
  k_steps <- length(chain) - 1L
  v_path <- state_table$vertex[chain]
  from_v <- v_path[1L]
  to_v <- v_path[length(v_path)]
  if (!all(endpoint_is_qualified(g, reference, admissible_horizons, control_cohorts,
                                 final_pool_mask, c(from_v, to_v)))) {
    cdid_abort("Segment endpoints must be qualified under the final trimmed pool.",
               "cdid_bridge_input_error")
  }

  v_qual <- endpoint_is_qualified(g, reference, admissible_horizons, control_cohorts, final_pool_mask, v_path)
  v_roles <- endpoint_coordinate_role(g, reference, admissible_horizons, control_cohorts, final_pool_mask, v_path)
  interior_idx <- if (length(v_path) >= 3L) seq.int(2L, length(v_path) - 1L) else integer(0L)
  ancestors <- if (k_steps >= 2L) chain[seq.int(2L, length(chain) - 1L)] else integer(0L)
  n_discarded_anc <- if (length(ancestors) > 0L) sum(!state_table$active[ancestors]) else 0L

  step_ids <- if (k_steps > 0L) chain[-1L] else integer(0L)
  e_idx <- state_table$edge[step_ids]
  signs <- state_table$sign[step_ids]
  s_vec <- if (k_steps > 0L) target_edges$s[e_idx] else numeric(0L)
  t_vec <- if (k_steps > 0L) target_edges$t[e_idx] else numeric(0L)
  from_d <- if (k_steps > 0L) v_path[-length(v_path)] else numeric(0L)
  to_d <- if (k_steps > 0L) v_path[-1L] else numeric(0L)

  s_qual <- endpoint_is_qualified(g, reference, admissible_horizons, control_cohorts, final_pool_mask, s_vec)
  t_qual <- endpoint_is_qualified(g, reference, admissible_horizons, control_cohorts, final_pool_mask, t_vec)
  s_role <- endpoint_coordinate_role(g, reference, admissible_horizons, control_cohorts, final_pool_mask, s_vec)
  t_role <- endpoint_coordinate_role(g, reference, admissible_horizons, control_cohorts, final_pool_mask, t_vec)
  is_early <- (s_vec < g & s_vec != reference) | (t_vec < g & t_vec != reference)
  is_post_nuis <- (s_role == "post_nuisance") | (t_role == "post_nuisance")

  edges_df <- data.frame(
    step = seq_len(k_steps),
    edge_index = e_idx,
    s = s_vec,
    t = t_vec,
    from_date = from_d,
    to_date = to_d,
    sign = signs,
    s_qualified = s_qual,
    t_qualified = t_qual,
    s_coordinate_role = s_role,
    t_coordinate_role = t_role,
    is_early_edge = is_early,
    uses_post_nuisance_endpoint = is_post_nuis,
    stringsAsFactors = FALSE
  )

  list(
    from_vertex = from_v,
    to_vertex = to_v,
    from_discovery_order = state_table$root_disc[sid],
    horizon = max_qual_h,
    max_qualified_horizon = max_qual_h,
    control_pool = control_pool,
    control_pool_mask = final_pool_mask,
    control_pool_key = pool_key,
    untrimmed_pool = untrimmed_pool,
    untrimmed_pool_mask = untrimmed_mask,
    untrimmed_pool_key = untrimmed_key,
    trimmed_at_destination = trimmed_flag,
    operator_key = op_key,
    vertex_path = v_path,
    vertex_qualified = v_qual,
    vertex_roles = v_roles,
    edges = edges_df,
    uses_early_vertex = any(v_roles[interior_idx] == "pre_nuisance"),
    uses_post_nuisance_vertex = any(v_roles[interior_idx] == "post_nuisance"),
    uses_later_date = max(v_path) > max(from_v, to_v),
    max_date_visited = max(v_path),
    state_chain = chain,
    n_discarded_ancestors_used = as.integer(n_discarded_anc)
  )
}

# Multi-source dominance-pruned bridge search over all calendar edges for CP02.
search_general_bridge_certificates <- function(g,
                                               calendar,
                                               cohorts,
                                               cohort_masses,
                                               fitting_universe,
                                               pair_support,
                                               target_post_times,
                                               overlap_oracle,
                                               max_states,
                                               reference = NULL,
                                               anticipation = 0,
                                               retain_state_table = FALSE) {
  if (!is.logical(retain_state_table) || length(retain_state_table) != 1L || is.na(retain_state_table)) {
    cdid_abort("`retain_state_table` must be TRUE or FALSE.", "cdid_bridge_input_error")
  }
  inp <- validate_bridge_search_inputs(
    g = g,
    calendar = calendar,
    cohorts = cohorts,
    cohort_masses = cohort_masses,
    fitting_universe = fitting_universe,
    pair_support = pair_support,
    target_post_times = target_post_times,
    overlap_oracle = overlap_oracle,
    max_states = max_states,
    reference = reference,
    anticipation = anticipation
  )
  g <- inp$g
  reference <- inp$reference
  calendar <- inp$calendar
  cohorts <- inp$cohorts
  cohort_masses <- inp$cohort_masses
  fitting_universe <- inp$fitting_universe
  ps <- inp$pair_support
  targets <- inp$target_post_times
  oracle <- inp$overlap_oracle
  max_states <- inp$max_states

  has_nt <- any(is.infinite(cohorts))
  latest_cohort <- max(cohorts)
  g_bar <- if (has_nt) Inf else latest_cohort
  control_cohorts <- cohorts[cohorts > g]
  n_ctrl <- length(control_cohorts)
  n_cal <- length(calendar)

  # Admissible post horizons v >= g with v < g_bar (where E_v = {k > v} is non-empty).
  admissible_horizons <- calendar[calendar >= g & calendar < g_bar]
  qualified_vertices <- c(reference, admissible_horizons)
  n_qual <- length(qualified_vertices)

  # Filter pair_support to edges supported by target g, in canonical (s, t) order.
  has_target <- vapply(ps$supported_cohorts, function(v) g %in% v, logical(1L))
  target_edges <- ps[has_target, , drop = FALSE]
  rownames(target_edges) <- NULL
  n_edges <- nrow(target_edges)
  edge_ctrl_masks <- lapply(seq_len(n_edges), function(e) {
    control_cohorts %in% target_edges$supported_cohorts[[e]]
  })

  # Full-calendar adjacency list (all dates in `calendar` may be traversed).
  cal_adjacency <- replicate(n_cal, integer(0L), simplify = FALSE)
  if (n_edges > 0L) {
    ends <- c(match(target_edges$s, calendar), match(target_edges$t, calendar))
    groups <- split(rep(seq_len(n_edges), 2L), ends)
    cal_adjacency[as.integer(names(groups))] <- lapply(groups, sort)
  }

  overlap_cache <- new_overlap_cache()

  # Reachability and acyclic discovery dependency tracking over `qualified_vertices`.
  qual_recovered <- logical(n_qual)
  qual_recovered[1L] <- TRUE
  qual_discovery_order <- rep(NA_integer_, n_qual)
  qual_discovery_order[1L] <- 1L
  qual_mechanism <- rep("unrecovered", n_qual)
  qual_mechanism[1L] <- "anchor"
  qual_parent <- rep(NA_real_, n_qual)
  qual_parent_disc <- rep(NA_integer_, n_qual)
  qual_horizon <- rep(NA_real_, n_qual)
  qual_pool_key <- rep(NA_character_, n_qual)
  qual_untrimmed_pool_key <- rep(NA_character_, n_qual)
  qual_trimmed <- rep(FALSE, n_qual)
  qual_cert_state_id <- rep(0L, n_qual)
  qual_seed_state_id <- rep(0L, n_qual)
  qual_segment <- vector("list", n_qual)

  discovery_count <- 1L
  discovery_order_vertices <- numeric(n_qual)
  discovery_order_vertices[1L] <- reference

  # Append-only state table over all calendar dates (never mutates ancestors).
  st_vertex <- numeric(0L)
  st_pool <- list()
  st_pool_key <- character(0L)
  st_parent <- integer(0L)
  st_edge <- integer(0L)
  st_sign <- integer(0L)
  st_root <- numeric(0L)
  st_root_disc <- integer(0L)
  st_depth <- integer(0L)
  st_is_seed <- logical(0L)
  st_active <- logical(0L)
  n_st <- 0L

  active_by_vertex <- replicate(n_cal, integer(0L), simplify = FALSE)
  active_total <- 0L
  peak_active_states <- 0L
  max_active_per_vertex <- 0L
  total_states_created <- 0L

  n_transitions_evaluated <- 0L
  n_rejected_empty <- 0L
  n_rejected_overlap <- 0L
  n_duplicate_skipped <- 0L
  n_dominated_skipped <- 0L
  n_dominated_replaced <- 0L
  n_cross_origin_dominated <- 0L
  n_trim_evaluated <- 0L
  n_trim_rejected_empty <- 0L
  n_trim_rejected_overlap <- 0L
  n_trimmed_certified <- 0L
  n_seeds_submitted <- 0L
  n_seeds_accepted <- 0L
  search_incomplete <- FALSE

  queue <- integer(min(64, max_states))
  q_head <- 1L
  q_tail <- 0L

  memory_samples <- list()
  sample_retained_buffers <- function(phase) {
    frame <- parent.frame()
    collect <- function(nms) {
      present <- nms[vapply(nms, exists, logical(1L), envir = frame, inherits = FALSE)]
      mget(present, envir = frame, inherits = FALSE)
    }
    groups <- list(
      input_and_graph = collect(c("inp", "calendar", "cohorts", "cohort_masses", "ps", "targets",
        "control_cohorts", "admissible_horizons", "qualified_vertices", "target_edges",
        "edge_ctrl_masks", "cal_adjacency")),
      overlap = list(masks = overlap_cache$masks, values = overlap_cache$vals,
        keys = ls(overlap_cache$by_key, all.names = TRUE), lookup = as.list(overlap_cache$by_key, all.names = TRUE)),
      states_queue_and_indexes = collect(c("active_by_vertex", "st_vertex", "st_pool", "st_pool_key",
        "st_parent", "st_edge", "st_sign", "st_root", "st_root_disc", "st_depth", "st_is_seed",
        "st_active", "queue")),
      reconstruction_and_output = collect(c("qual_recovered", "qual_discovery_order", "qual_mechanism",
        "qual_parent", "qual_parent_disc", "qual_horizon", "qual_pool_key", "qual_untrimmed_pool_key",
        "qual_trimmed", "qual_cert_state_id", "qual_seed_state_id", "qual_segment",
        "discovery_order_vertices", "certificates", "segs", "seg", "target_summary",
        "qualified_reachability", "discovery_dependencies", "state_table_out")),
      diagnostics = collect(c("memory_samples"))
    )
    accounted <- unique(unlist(lapply(groups, names), use.names = FALSE))
    other <- collect(setdiff(ls(frame, all.names = TRUE), c(accounted, "overlap_cache")))
    groups$other_retained_data <- other[!vapply(other, function(x) is.function(x) || is.environment(x), logical(1L))]
    sizes <- vapply(groups, function(x) as.double(object.size(x)), numeric(1L))
    sample <- as.data.frame(as.list(c(list(phase = phase), as.list(sizes), list(total_payload_bytes = sum(sizes)))))
    memory_samples[[length(memory_samples) + 1L]] <<- sample
    invisible(NULL)
  }

  # Dominance check and state insertion at calendar date `v_date`.
  # Superset active states at `v_date` dominate subset candidate states even across
  # distinct certified seed origins `root_date`.
  try_insert_state <- function(v_date, p_mask, p_key, parent_id, edge_id, sgn, root_date, root_disc, is_seed = FALSE) {
    vi <- match(v_date, calendar)
    act_v <- active_by_vertex[[vi]]
    if (length(act_v) > 0L) {
      for (aid in act_v) {
        if (st_pool_key[aid] == p_key) {
          n_duplicate_skipped <<- n_duplicate_skipped + 1L
          if (st_root[aid] != root_date) {
            n_cross_origin_dominated <<- n_cross_origin_dominated + 1L
          }
          return(0L)
        }
        if (pool_is_subset(p_mask, st_pool[[aid]])) {
          n_dominated_skipped <<- n_dominated_skipped + 1L
          if (st_root[aid] != root_date) {
            n_cross_origin_dominated <<- n_cross_origin_dominated + 1L
          }
          return(0L)
        }
      }
    }
    if (total_states_created >= max_states) {
      search_incomplete <<- TRUE
      return(0L)
    }
    if (length(act_v) > 0L) {
      keep_act <- logical(length(act_v))
      for (idx_a in seq_along(act_v)) {
        aid <- act_v[idx_a]
        if (pool_is_strict_subset(st_pool[[aid]], p_mask)) {
          st_active[aid] <<- FALSE
          n_dominated_replaced <<- n_dominated_replaced + 1L
          if (st_root[aid] != root_date) {
            n_cross_origin_dominated <<- n_cross_origin_dominated + 1L
          }
          active_total <<- active_total - 1L
          keep_act[idx_a] <- FALSE
        } else {
          keep_act[idx_a] <- TRUE
        }
      }
      act_v <- act_v[keep_act]
    }
    n_st <<- n_st + 1L
    total_states_created <<- total_states_created + 1L
    st_vertex[n_st] <<- v_date
    st_pool[[n_st]] <<- p_mask
    st_pool_key[n_st] <<- p_key
    st_parent[n_st] <<- as.integer(parent_id)
    st_edge[n_st] <<- as.integer(edge_id)
    st_sign[n_st] <<- as.integer(sgn)
    st_root[n_st] <<- root_date
    st_root_disc[n_st] <<- as.integer(root_disc)
    st_depth[n_st] <<- if (parent_id == 0L) 0L else as.integer(st_depth[parent_id] + 1L)
    st_is_seed[n_st] <<- is_seed
    st_active[n_st] <<- TRUE

    act_v <- c(act_v, n_st)
    active_by_vertex[[vi]] <<- act_v
    active_total <<- active_total + 1L
    if (active_total > peak_active_states) peak_active_states <<- active_total
    if (length(act_v) > max_active_per_vertex) max_active_per_vertex <<- length(act_v)

    if (q_tail == length(queue)) length(queue) <<- min(max_states, 2 * length(queue))
    q_tail <<- q_tail + 1L
    queue[q_tail] <<- n_st
    n_st
  }

  submit_seed_for_qualified <- function(q_date, q_disc) {
    seed_mask <- control_cohorts > max(g, q_date)
    if (!any(seed_mask)) return(0L)
    n_seeds_submitted <<- n_seeds_submitted + 1L
    if (!prime_pool_overlap_checks(oracle, g, control_cohorts, seed_mask, overlap_cache)) {
      return(0L)
    }
    seed_key <- pool_canonical_key(seed_mask)
    sid <- try_insert_state(
      v_date = q_date,
      p_mask = seed_mask,
      p_key = seed_key,
      parent_id = 0L,
      edge_id = 0L,
      sgn = 0L,
      root_date = q_date,
      root_disc = q_disc,
      is_seed = TRUE
    )
    if (sid > 0L) {
      n_seeds_accepted <<- n_seeds_accepted + 1L
      qi <- match(q_date, qualified_vertices)
      qual_seed_state_id[qi] <<- sid
    }
    sid
  }

  sample_retained_buffers("initialized")

  # Seed reference anchor b = g - 1 with E_g = {k > g}.
  if (n_ctrl > 0L && !all(qual_recovered)) {
    submit_seed_for_qualified(reference, 1L)
  }

  # Accepted (date, pool) pairs never reappear: replacement preserves an active
  # superset forever. Descendant pools decrease, so ancestor dates cannot recur.
  # Every root was discovered before its descendants can certify a new date.
  # Stop at queue exhaustion, saturation of all qualified ATT levels, or the cap.
  # Saturation proves all requested levels but need not exhaust the state graph.
  while (q_head <= q_tail && !search_incomplete && !all(qual_recovered)) {
    curr_id <- queue[q_head]
    q_head <- q_head + 1L
    if (!st_active[curr_id]) next

    u_date <- st_vertex[curr_id]
    u_mask <- st_pool[[curr_id]]
    u_root <- st_root[curr_id]
    u_root_disc <- st_root_disc[curr_id]
    ui <- match(u_date, calendar)

    for (e in cal_adjacency[[ui]]) {
      se <- target_edges$s[e]
      te <- target_edges$t[e]
      if (se == u_date) {
        w_date <- te; sgn <- 1L
      } else if (te == u_date) {
        w_date <- se; sgn <- -1L
      } else {
        next
      }
      n_transitions_evaluated <- n_transitions_evaluated + 1L
      cand_mask <- u_mask & edge_ctrl_masks[[e]]
      if (!any(cand_mask)) {
        n_rejected_empty <- n_rejected_empty + 1L
        next
      }
      if (!evaluate_overlap_mask(oracle, g, control_cohorts, cand_mask, overlap_cache)) {
        n_rejected_overlap <- n_rejected_overlap + 1L
        next
      }
      cand_key <- pool_canonical_key(cand_mask)

      new_sid <- try_insert_state(
        v_date = w_date,
        p_mask = cand_mask,
        p_key = cand_key,
        parent_id = curr_id,
        edge_id = e,
        sgn = sgn,
        root_date = u_root,
        root_disc = u_root_disc,
        is_seed = FALSE
      )
      if (new_sid == 0L) {
        if (search_incomplete) break
        next
      }

      # At any uncertified post date v in H_g, TEST recovery with K intersect E_v.
      # cand_mask and every stored pool stay untrimmed for further propagation.
      # If non-empty and overlap-valid, reinterpret the immutable ancestor path with
      # that final pool, certify ATT(g, v), and submit its maximal new seed E_v.
      if (w_date %in% admissible_horizons) {
        qi <- match(w_date, qualified_vertices)
        if (!qual_recovered[qi]) {
          n_trim_evaluated <- n_trim_evaluated + 1L
          e_w_mask <- control_cohorts > w_date
          trim_mask <- cand_mask & e_w_mask
          if (!any(trim_mask)) {
            n_trim_rejected_empty <- n_trim_rejected_empty + 1L
          } else if (!evaluate_overlap_mask(oracle, g, control_cohorts, trim_mask, overlap_cache)) {
            n_trim_rejected_overlap <- n_trim_rejected_overlap + 1L
          } else {
            seg <- reconstruct_general_segment_from_state(
              sid = new_sid,
              final_pool_mask = trim_mask,
              state_table = list(
                vertex = st_vertex,
                pool = st_pool,
                pool_key = st_pool_key,
                parent = st_parent,
                edge = st_edge,
                sign = st_sign,
                root = st_root,
                root_disc = st_root_disc,
                active = st_active
              ),
              g = g,
              reference = reference,
              fitting_universe = fitting_universe,
              cohorts = cohorts,
              cohort_masses = cohort_masses,
              control_cohorts = control_cohorts,
              admissible_horizons = admissible_horizons,
              target_edges = target_edges
            )
            discovery_count <- discovery_count + 1L
            discovery_order_vertices[discovery_count] <- w_date
            seg$to_discovery_order <- discovery_count
            if (seg$from_discovery_order >= discovery_count) {
              cdid_abort("Discovery dependency cycle detected.", "cdid_bridge_input_error")
            }

            qual_recovered[qi] <- TRUE
            qual_discovery_order[qi] <- discovery_count
            qual_mechanism[qi] <- if (nrow(seg$edges) == 1L && !seg$uses_early_vertex &&
                                      !seg$uses_post_nuisance_vertex && !seg$trimmed_at_destination) {
              "post_only_edge"
            } else {
              "general_fixed_operator_bridge"
            }
            qual_parent[qi] <- seg$from_vertex
            qual_parent_disc[qi] <- seg$from_discovery_order
            qual_horizon[qi] <- seg$max_qualified_horizon
            qual_pool_key[qi] <- seg$control_pool_key
            qual_untrimmed_pool_key[qi] <- seg$untrimmed_pool_key
            qual_trimmed[qi] <- seg$trimmed_at_destination
            qual_cert_state_id[qi] <- new_sid
            qual_segment[[qi]] <- seg
            if (seg$trimmed_at_destination) {
              n_trimmed_certified <- n_trimmed_certified + 1L
            }

            if (all(qual_recovered)) break
            submit_seed_for_qualified(w_date, discovery_count)
            if (search_incomplete) break
          }
        }
      }
    }
  }

  termination_reason <- if (search_incomplete) "state_cap" else if (all(qual_recovered)) {
    "all_qualified_levels_recovered"
  } else "queue_exhausted"
  state_queue_exhausted <- q_head > q_tail
  sample_retained_buffers("after_search")

  # Refresh discarded-ancestor counts after the full BFS finishes.
  for (qi in seq_len(n_qual)) {
    seg_q <- qual_segment[[qi]]
    if (!is.null(seg_q) && length(seg_q$state_chain) >= 3L) {
      anc_ids <- seg_q$state_chain[seq.int(2L, length(seg_q$state_chain) - 1L)]
      qual_segment[[qi]]$n_discarded_ancestors_used <- as.integer(sum(!st_active[anc_ids]))
    }
  }

  # Assemble full certificate chains from anchor `reference` to each certified target.
  certificates <- list()
  n_targets <- length(targets)
  t_eligible <- targets < g_bar
  t_certified <- logical(n_targets)
  t_status <- character(n_targets)
  t_n_segments <- rep(NA_integer_, n_targets)
  t_n_edges <- rep(NA_integer_, n_targets)
  t_uses_early <- rep(FALSE, n_targets)
  t_uses_post_nuisance <- rep(FALSE, n_targets)
  t_uses_trimmed <- rep(FALSE, n_targets)
  t_uses_later <- rep(FALSE, n_targets)
  t_uses_mixed <- rep(FALSE, n_targets)
  total_discarded_anc_used <- 0L

  for (ti in seq_len(n_targets)) {
    t_val <- targets[ti]
    if (!t_eligible[ti]) {
      t_certified[ti] <- FALSE
      t_status[ti] <- "terminal_ineligible_no_nt"
      next
    }
    qi <- match(t_val, qualified_vertices)
    if (!is.na(qi) && qual_recovered[qi]) {
      segs <- list()
      cur_v <- t_val
      visited_disc <- integer(0L)
      while (cur_v != reference) {
        c_idx <- match(cur_v, qualified_vertices)
        c_disc <- qual_discovery_order[c_idx]
        if (is.na(c_disc) || c_disc %in% visited_disc) {
          cdid_abort("Discovery dependency cycle encountered during certificate assembly.",
                     "cdid_bridge_input_error")
        }
        visited_disc <- c(visited_disc, c_disc)
        seg <- qual_segment[[c_idx]]
        segs <- c(list(seg), segs)
        cur_v <- seg$from_vertex
      }
      for (j in seq_along(segs)) {
        segs[[j]]$segment_index <- j
        total_discarded_anc_used <- total_discarded_anc_used + segs[[j]]$n_discarded_ancestors_used
      }
      n_seg <- length(segs)
      n_ed <- sum(vapply(segs, function(s) nrow(s$edges), integer(1L)))
      u_early <- any(vapply(segs, function(s) s$uses_early_vertex, logical(1L)))
      u_post_nuis <- any(vapply(segs, function(s) s$uses_post_nuisance_vertex, logical(1L)))
      u_trimmed <- any(vapply(segs, function(s) s$trimmed_at_destination, logical(1L)))
      u_later <- any(vapply(segs, function(s) s$max_date_visited > t_val, logical(1L)))
      pool_keys <- vapply(segs, function(s) s$control_pool_key, character(1L))
      u_mixed <- length(unique(pool_keys)) > 1L

      t_certified[ti] <- TRUE
      t_status[ti] <- "certified"
      t_n_segments[ti] <- as.integer(n_seg)
      t_n_edges[ti] <- as.integer(n_ed)
      t_uses_early[ti] <- u_early
      t_uses_post_nuisance[ti] <- u_post_nuis
      t_uses_trimmed[ti] <- u_trimmed
      t_uses_later[ti] <- u_later
      t_uses_mixed[ti] <- u_mixed

      certificates[[format_cohort_value(t_val)]] <- list(
        g = g,
        t = t_val,
        reference = reference,
        segments = segs,
        n_segments = as.integer(n_seg),
        n_total_edges = as.integer(n_ed),
        uses_early_vertex = u_early,
        uses_post_nuisance_vertex = u_post_nuis,
        uses_trimmed_endpoint = u_trimmed,
        uses_later_date = u_later,
        uses_mixed_operators = u_mixed
      )
    } else {
      t_certified[ti] <- FALSE
      t_status[ti] <- if (search_incomplete) "search_incomplete" else "no_certificate_in_C_gen_g"
    }
  }

  target_summary <- data.frame(
    g = rep(g, n_targets),
    t = targets,
    reference = rep(reference, n_targets),
    terminal_eligible = t_eligible,
    certified = t_certified,
    status = t_status,
    n_segments = t_n_segments,
    n_total_edges = t_n_edges,
    uses_early_vertex = t_uses_early,
    uses_post_nuisance_vertex = t_uses_post_nuisance,
    uses_trimmed_endpoint = t_uses_trimmed,
    uses_later_date = t_uses_later,
    uses_mixed_operators = t_uses_mixed,
    stringsAsFactors = FALSE
  )

  qualified_reachability <- data.frame(
    vertex = qualified_vertices,
    role = ifelse(qualified_vertices == reference, "reference", "post"),
    recovered = qual_recovered,
    discovery_order = qual_discovery_order,
    mechanism = qual_mechanism,
    parent_qualified_vertex = qual_parent,
    parent_discovery_order = qual_parent_disc,
    max_qualified_horizon = qual_horizon,
    control_pool_key = qual_pool_key,
    untrimmed_pool_key = qual_untrimmed_pool_key,
    trimmed_at_destination = qual_trimmed,
    certifying_state_id = qual_cert_state_id,
    seed_state_id = qual_seed_state_id,
    stringsAsFactors = FALSE
  )

  disc_vertices <- discovery_order_vertices[seq_len(discovery_count)]
  disc_qi <- match(disc_vertices, qualified_vertices)
  discovery_dependencies <- qualified_reachability[disc_qi, , drop = FALSE]
  rownames(discovery_dependencies) <- NULL

  state_table_out <- if (retain_state_table) {
    list(
      summary = data.frame(
        state_id = seq_len(n_st),
        vertex = st_vertex,
        pool_key = st_pool_key,
        parent_state_id = st_parent,
        edge_index = st_edge,
        sign = st_sign,
        root_vertex = st_root,
        root_discovery_order = st_root_disc,
        depth = st_depth,
        is_seed = st_is_seed,
        active = st_active,
        stringsAsFactors = FALSE
      ),
      pools = st_pool
    )
  } else {
    NULL
  }

  sample_retained_buffers("assembled_certificates")
  structure(
    list(
      operator_spec = list(
        g = g,
        reference = reference,
        anticipation = 0,
        fitting_universe = fitting_universe,
        calendar = calendar,
        cohorts = cohorts,
        control_cohorts = control_cohorts,
        cohort_masses = cohort_masses,
        has_never_treated = has_nt,
        latest_cohort = latest_cohort,
        admissible_horizons = admissible_horizons,
        qualified_vertices = qualified_vertices,
        max_states = max_states,
        certificate_class = "C_gen_g_arbitrary_date_fixed_pool_composed_at_qualified_vertices"
      ),
      target_summary = target_summary,
      certificates = certificates,
      qualified_reachability = qualified_reachability,
      discovery_dependencies = discovery_dependencies,
      state_table = state_table_out,
      scope_conditions = c(
        "Certificates are conditional on supplied population assumptions; observed counts and labels do not certify them.",
        "Fixed common cohort laws and fixed unnormalized positive relative masses independent of K; positive population pair support for every selected cohort/edge.",
        "Joint representative Q sampling, baseline pair MTAR (BBD9), target-g full CS5 through each certified endpoint and its applicable later-horizon population overlap, no anticipation, integrable reconstructed links and coordinates.",
        "Endpoint qualification is sufficient for the ATT identity. Early and unqualified post coordinates are kept private by class convention, even if additional conditions imply equalities.",
        "Completed negative search (no_certificate_in_C_gen_g) concerns anchored ATT levels only in C^gen_g, omitting cross-operator linear mixture relations and cross-cohort systems; it does not prove population nonidentification or absence of estimable unanchored contrasts.",
        "Certificate design output contains witness rows only. Preserving unanchored components and contrasts in a production GMM moment catalog remains the later D5 obligation."
      ),
      search_diagnostics = list(
        search_incomplete = search_incomplete,
        termination_reason = termination_reason,
        state_queue_exhausted = state_queue_exhausted,
        max_states = max_states,
        total_states_created = as.integer(total_states_created),
        peak_active_states = as.integer(peak_active_states),
        max_active_per_vertex = as.integer(max_active_per_vertex),
        n_transitions_evaluated = as.integer(n_transitions_evaluated),
        n_rejected_empty = as.integer(n_rejected_empty),
        n_rejected_overlap = as.integer(n_rejected_overlap),
        n_duplicate_skipped = as.integer(n_duplicate_skipped),
        n_dominated_skipped = as.integer(n_dominated_skipped),
        n_dominated_replaced = as.integer(n_dominated_replaced),
        n_cross_origin_dominated = as.integer(n_cross_origin_dominated),
        n_trim_evaluated = as.integer(n_trim_evaluated),
        n_trim_rejected_empty = as.integer(n_trim_rejected_empty),
        n_trim_rejected_overlap = as.integer(n_trim_rejected_overlap),
        n_trimmed_certified = as.integer(n_trimmed_certified),
        n_seeds_submitted = as.integer(n_seeds_submitted),
        n_seeds_accepted = as.integer(n_seeds_accepted),
        n_discarded_ancestors_used_in_certificates = as.integer(total_discarded_anc_used),
        n_overlap_evaluations = as.integer(overlap_cache$n_evaluations),
        retained_overlap_pools = length(overlap_cache$vals),
        overlap_cache_payload_bytes = as.double(object.size(list(overlap_cache$masks, overlap_cache$vals,
          ls(overlap_cache$by_key, all.names = TRUE), as.list(overlap_cache$by_key, all.names = TRUE)))),
        retained_buffer_samples = do.call(rbind, memory_samples),
        memory_measurement = "R object.size payload samples; duplicate contents counted; hash/allocator/closure/transient overhead excluded. Measure full process peak RSS separately. State cap does not cap overlap-cache entries.",
        custom_predicate_captured_memory_unmeasured = identical(oracle$kind, "declared_monotone_predicate")
      )
    ),
    class = c("cdid_general_bridge_certificate_result", "cdid_bridge_certificate_result")
  )
}

search_bridge_certificates_general <- search_general_bridge_certificates

# Minimal bridge-to-design specification for C^gen_g (without modifying R/chained_design.R).
# Describes required endpoint qualification metadata (`s_qualified`, `t_qualified`,
# `s_coordinate_role`, `t_coordinate_role`, `max_qualified_horizon`) and separates
# shared qualified ATT(g, u) vertices from operator-specific post-treatment nuisance
# coordinates beta_post_nuisance(g, u, fn).
general_bridge_certificates_to_design_links <- function(result) {
  if (!inherits(result, "cdid_general_bridge_certificate_result")) {
    cdid_abort("Expected `cdid_general_bridge_certificate_result`.", "cdid_bridge_input_error")
  }
  spec <- result$operator_spec
  rows <- list()
  k <- 0L
  for (cert in result$certificates) {
    for (seg in cert$segments) {
      if (!length(seg$control_pool) || anyNA(seg$control_pool) ||
          !all(seg$control_pool %in% spec$control_cohorts) || anyDuplicated(seg$control_pool)) {
        cdid_abort("Certificate pool must be a non-empty subset of the target's controls.", "cdid_bridge_input_error")
      }
      expected_key <- make_operator_key(spec$g, spec$reference, seg$max_qualified_horizon,
                                        spec$fitting_universe, spec$cohorts, spec$cohort_masses,
                                        seg$control_pool)
      if (!identical(seg$operator_key, expected_key)) {
        cdid_abort("Certificate key differs from its declared fixed mixture.", "cdid_bridge_input_error")
      }
      ed <- seg$edges
      # Qualification is intrinsic to (g, operator, date), independent of the
      # certificate's destination or discovery order. Reject contradictory copies.
      mask <- spec$control_cohorts %in% seg$control_pool
      intrinsic_horizon <- max_qualified_post_horizon(spec$admissible_horizons, spec$control_cohorts, mask)
      if (!identical(seg$max_qualified_horizon, intrinsic_horizon)) {
        cdid_abort("Qualification horizon contradicts its declared fixed operator.", "cdid_bridge_input_error")
      }
      s_qual <- endpoint_is_qualified(spec$g, spec$reference, spec$admissible_horizons,
                                     spec$control_cohorts, mask, ed$s)
      t_qual <- endpoint_is_qualified(spec$g, spec$reference, spec$admissible_horizons,
                                     spec$control_cohorts, mask, ed$t)
      s_role <- endpoint_coordinate_role(spec$g, spec$reference, spec$admissible_horizons,
                                        spec$control_cohorts, mask, ed$s)
      t_role <- endpoint_coordinate_role(spec$g, spec$reference, spec$admissible_horizons,
                                        spec$control_cohorts, mask, ed$t)
      if (!identical(ed$s_qualified, s_qual) || !identical(ed$t_qualified, t_qual) ||
          !identical(ed$s_coordinate_role, s_role) || !identical(ed$t_coordinate_role, t_role)) {
        cdid_abort("Endpoint metadata contradicts its declared fixed operator.", "cdid_bridge_input_error")
      }
      cg <- if (all(is.infinite(seg$control_pool))) "nevertreated" else "notyettreated"
      pool_str <- format_cohort_set(seg$control_pool)
      for (i in seq_len(nrow(ed))) {
        s_q <- ed$s_qualified[i]
        t_q <- ed$t_qualified[i]
        s_r <- ed$s_coordinate_role[i]
        t_r <- ed$t_coordinate_role[i]
        # Preserve the old private contract: every POST endpoint qualifies;
        # an early nuisance endpoint is allowed on a straddling row.
        post_sh <- isTRUE((ed$s[i] < spec$g || s_q) && (ed$t[i] < spec$g || t_q))
        k <- k + 1L
        rows[[k]] <- data.frame(
          g = spec$g,
          s = ed$s[i],
          t = ed$t[i],
          ref = spec$reference,
          functional_key = seg$operator_key,
          control_group = cg,
          control_pool = pool_str,
          max_qualified_horizon = seg$max_qualified_horizon,
          s_qualified = s_q,
          t_qualified = t_q,
          s_coordinate_role = s_r,
          t_coordinate_role = t_r,
          post_shareable = post_sh,
          uses_post_nuisance = (s_r == "post_nuisance" || t_r == "post_nuisance"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (k == 0L) {
    return(data.frame(
      g = numeric(0L), s = numeric(0L), t = numeric(0L), ref = numeric(0L),
      functional_key = character(0L), control_group = character(0L),
      control_pool = character(0L), max_qualified_horizon = numeric(0L),
      s_qualified = logical(0L), t_qualified = logical(0L),
      s_coordinate_role = character(0L), t_coordinate_role = character(0L),
      post_shareable = logical(0L), uses_post_nuisance = logical(0L),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  dedup_key <- paste(
    match(out$s, spec$calendar),
    match(out$t, spec$calendar),
    out$functional_key,
    sep = ":"
  )
  out <- out[!duplicated(dedup_key), , drop = FALSE]
  ord <- order(out$s, out$t, out$functional_key, method = "radix")
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Build a design specification for certificate WITNESS rows only (links,
# vertices, parameters, sparse incidence). This omits disconnected moments and
# unanchored contrasts; it is not the complete production GMM moment catalog.
general_bridge_certificates_to_design_spec <- function(result) {
  links <- general_bridge_certificates_to_design_links(result)
  spec <- result$operator_spec
  if (nrow(links) == 0L) {
    return(structure(
      list(
        links = links,
        vertices = data.frame(
          vertex_index = integer(0L), g = numeric(0L), time = numeric(0L),
          role = character(0L), qualified = logical(0L),
          functional_key = character(0L), control_pool = character(0L),
          vertex_id = character(0L), param_index = integer(0L),
          stringsAsFactors = FALSE
        ),
        parameters = data.frame(),
        incidence = list(start = integer(0L), end = integer(0L), n_links = 0L, n_parameters = 0L)
      ),
      class = "cdid_general_bridge_design_spec"
    ))
  }
  fn <- sort(unique(links$functional_key), method = "radix")
  fi <- match(links$functional_key, fn)
  dates <- spec$calendar
  endpoint_key <- function(time, role, f_idx) {
    # -1: reference; 0: shared post_att; positive f_idx: operator-specific pre or post nuisance
    kind <- ifelse(role == "reference", -1L, ifelse(role == "post_att", 0L, f_idx))
    paste(match(time, dates), role, kind, sep = ":")
  }
  s_keys <- endpoint_key(links$s, links$s_coordinate_role, fi)
  t_keys <- endpoint_key(links$t, links$t_coordinate_role, fi)

  all_times <- c(spec$reference, links$s, links$t)
  all_roles <- c("reference", links$s_coordinate_role, links$t_coordinate_role)
  all_quals <- c(TRUE, links$s_qualified, links$t_qualified)
  all_f <- c(0L, ifelse(links$s_qualified, 0L, fi), ifelse(links$t_qualified, 0L, fi))
  all_pools <- c(NA_character_, ifelse(links$s_qualified, NA_character_, links$control_pool),
                 ifelse(links$t_qualified, NA_character_, links$control_pool))
  all_keys <- endpoint_key(all_times, all_roles, all_f)

  keep <- !duplicated(all_keys)
  role_rank <- match(all_roles[keep], c("reference", "pre_nuisance", "post_nuisance", "post_att"))
  ord <- which(keep)[order(role_rank, all_times[keep], all_f[keep], method = "radix")]

  v_time <- all_times[ord]
  v_role <- all_roles[ord]
  v_qual <- all_quals[ord]
  v_f <- all_f[ord]
  v_fk <- ifelse(v_f == 0L, NA_character_, fn[pmax(v_f, 1L)])
  v_pool <- all_pools[ord]
  v_keys <- all_keys[ord]

  v_label <- ifelse(v_role == "reference",
    paste0("ref(g=", format_cohort_value(spec$g), ",t=", format_cohort_value(v_time), ")"),
    ifelse(v_role == "post_att",
      paste0("ATT(g=", format_cohort_value(spec$g), ",t=", format_cohort_value(v_time), ")"),
      ifelse(v_role == "pre_nuisance",
        paste0("beta_pre(g=", format_cohort_value(spec$g), ",t=", format_cohort_value(v_time), ",pool=", v_pool, ")"),
        paste0("beta_post_nuisance(g=", format_cohort_value(spec$g), ",t=", format_cohort_value(v_time), ",pool=", v_pool, ")"))))

  nv <- length(ord)
  pv <- which(v_role != "reference")
  p <- length(pv)
  p_idx <- rep(NA_integer_, nv)
  p_idx[pv] <- seq_len(p)

  vertices <- data.frame(
    vertex_index = seq_len(nv),
    g = rep(spec$g, nv),
    time = v_time,
    role = v_role,
    qualified = v_qual,
    functional_key = v_fk,
    control_pool = v_pool,
    vertex_id = make.unique(v_label),
    param_index = p_idx,
    stringsAsFactors = FALSE
  )

  from_v <- match(s_keys, v_keys)
  to_v <- match(t_keys, v_keys)
  links$link_index <- seq_len(nrow(links))
  links$link_id <- paste0("link_", seq_len(nrow(links)))
  links$start_vertex_index <- from_v
  links$end_vertex_index <- to_v
  links$start_param_index <- vertices$param_index[from_v]
  links$end_param_index <- vertices$param_index[to_v]

  parameters <- vertices[pv, , drop = FALSE]
  rownames(parameters) <- NULL
  parameters$param_id <- parameters$vertex_id
  parameters$is_target_att <- parameters$role == "post_att"

  incidence <- list(
    start = ifelse(is.na(links$start_param_index), 0L, links$start_param_index),
    end = ifelse(is.na(links$end_param_index), 0L, links$end_param_index),
    n_links = nrow(links),
    n_parameters = p,
    link_ids = links$link_id,
    param_ids = parameters$param_id
  )

  structure(
    list(
      links = links,
      vertices = vertices,
      parameters = parameters,
      incidence = incidence
    ),
    class = "cdid_general_bridge_design_spec"
  )
}

# Bounded debug/test dense incidence materialization for a general bridge design spec.
materialize_general_design_spec <- function(design_spec, max_entries) {
  if (!inherits(design_spec, "cdid_general_bridge_design_spec")) {
    cdid_abort("Expected `cdid_general_bridge_design_spec`.", "cdid_bridge_input_error")
  }
  bound <- bridge_bound(max_entries, "max_entries")
  inc <- design_spec$incidence
  m <- as.double(inc$n_links)
  p <- as.double(inc$n_parameters)
  if (!is.finite(m * p) || m * p > bound) {
    cdid_abort("Dense design spec materialization exceeds explicit entry bound.",
               "cdid_bridge_bound_error")
  }
  W <- matrix(0L, nrow = as.integer(m), ncol = as.integer(p),
              dimnames = list(inc$link_ids, inc$param_ids))
  s_pos <- which(inc$start > 0L)
  if (length(s_pos) > 0L) W[cbind(s_pos, inc$start[s_pos])] <- -1L
  e_pos <- which(inc$end > 0L)
  if (length(e_pos) > 0L) W[cbind(e_pos, inc$end[e_pos])] <- 1L
  out <- design_spec
  out$W <- W
  out
}
