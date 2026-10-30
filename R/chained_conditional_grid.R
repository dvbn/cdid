# CP08 O14 (D66 choice 2): integration grid of the chained conditional
# pretrend test. One grid row per supplied original unit F (singletons,
# zero-outcome and frame-only units included; duplicates keep their mass), the
# model.matrix(xformla) row of that unit at the cell's approved reference date
# C (the earlier endpoint of the cell's auxiliary pair, the D60 extraction date
# of auxiliary reports). Terms, factor levels, contrasts and transform state are
# those frozen on the route's own reference rows; grid-only rows are evaluated
# through that frozen state and never refit it. A missing required value or an
# unsupported level makes the complete test unavailable (with ids); no unit is
# deleted and no available-only grid is substituted.

# Exact binary64 agreement of two design blocks (attributes ignored; -0 = 0).
cpt_grid_same <- function(a, b) {
  identical(as.integer(dim(a)), as.integer(dim(b))) && all(as.double(a) == as.double(b))
}

cpt_grid_unavailable <- function(reason, missing = NULL, unsupported = NULL) {
  list(status = "unavailable", reason = reason, missing = missing, unsupported = unsupported)
}

cpt_grid_fmt_rows <- function(tab, what) {
  show <- utils::head(tab, 5L)
  txt <- paste0("id ", vapply(show$id, chained_reference_fmt, character(1L)),
    ifelse(is.na(show$date), "", paste0(" at date ", vapply(show$date, function(x)
      if (is.na(x)) "" else chained_reference_fmt(x), character(1L)))),
    ", variable `", show$variable, "`", if (!is.null(show$value)) paste0(" (value ", show$value, ")") else "",
    collapse = "; ")
  paste0(what, ": ", txt, if (nrow(tab) > 5L) paste0("; and ", nrow(tab) - 5L, " more") else "")
}

# Rows on which the route's static xformla design was evaluated (and hence its
# frozen state): every F row for the CP02 IPW gate, baseline reference IPW,
# pure imputation and the augmented correction; the repeat observers A for the
# outcome-regression, DR and conditional_mtar gates (cp04_public_design on A).
cpt_grid_reference_rows <- function(prep) {
  cp <- prep$cp04
  if (is.null(cp) || identical(cp$ipw_method, "baseline") ||
      isTRUE(cp$selection_est_method %in% c("imputation", "dr"))) return(rep(TRUE, prep$n_frame))
  prep$units$n_obs >= 2L
}

# Raw unit-level variables at extraction date C. An explicit frame column is
# a unit attribute. Reference rows retain the route's own frozen frame row;
# data-sourced rows outside that fitting roster require their actual date-C
# value, including an explicitly missing value (no cross-date substitution).
cpt_grid_static_raw <- function(src, vars, date_C, U) {
  N <- length(src$frame_ids)
  raw <- vector("list", length(vars)); names(raw) <- vars
  dated_units <- list()
  first_any <- match(seq_len(N), src$unit)
  at_C <- which(src$tt == date_C)
  rc <- rep(NA_integer_, N); rc[src$unit[at_C]] <- at_C
  for (v in vars) {
    if (src$frame_supplied && v %in% names(src$frame)) {
      raw[[v]] <- src$frame[[v]]
      next
    }
    if (!v %in% names(src$data)) return(list(error = paste0("Covariate `", v, "` of `xformla` is not supplied in `data` or `frame`.")))
    ri <- first_any
    ri[!U] <- rc[!U]
    dated_units[[v]] <- which(!U)
    raw[[v]] <- src$data[[v]][ri]
  }
  df <- as.data.frame(raw, stringsAsFactors = FALSE, optional = TRUE)
  if (!length(vars)) df <- data.frame(row.names = seq_len(N))
  list(df = df, dated_units = dated_units)
}

cpt_grid_missing_table <- function(df, ids, date_of) {
  if (!ncol(df)) return(NULL)
  na <- vapply(seq_along(df), function(j) is.na(df[[j]]), logical(nrow(df)))
  if (!is.matrix(na)) na <- matrix(na, nrow(df), ncol(df))
  if (!any(na)) return(NULL)
  hit <- which(na, arr.ind = TRUE)
  hit <- hit[order(hit[, 1L], hit[, 2L]), , drop = FALSE]
  data.frame(id = ids[hit[, 1L]], date = date_of(hit[, 1L], names(df)[hit[, 2L]]),
    variable = names(df)[hit[, 2L]], stringsAsFactors = FALSE)
}

cpt_grid_level_table <- function(df, xlev, ids) {
  out <- list()
  for (v in intersect(names(xlev), names(df))) {
    val <- df[[v]]
    lab <- if (is.factor(val)) as.character(val) else as.character(val)
    bad <- which(!is.na(lab) & !lab %in% xlev[[v]])
    if (length(bad)) out[[v]] <- data.frame(id = ids[bad], date = NA_real_, variable = v, value = lab[bad],
      stringsAsFactors = FALSE)
  }
  if (!length(out)) NULL else do.call(rbind, out)
}

# Model matrix of `df` through frozen terms/levels/contrasts; NULL plus a
# reason when the frozen state cannot be applied.
cpt_grid_frozen_matrix <- function(tt, df, xlev, contrasts) {
  mf <- tryCatch(stats::model.frame(tt, df, xlev = xlev, na.action = stats::na.pass), error = function(e) e)
  if (inherits(mf, "error")) return(list(error = conditionMessage(mf)))
  Z <- tryCatch(stats::model.matrix(tt, mf, contrasts.arg = contrasts), error = function(e) e)
  if (inherits(Z, "error")) return(list(error = conditionMessage(Z)))
  list(Z = Z)
}

# Admit the full-F matrix and model-frame/transformation temporaries before
# evaluating the grid. Existing cached recipes remain live during construction.
cpt_grid_admit <- function(src, p, n_reference, n_variables) {
  b <- src$grid_budget
  if (is.null(b)) return(invisible(NULL))
  ipw_owner_safe_sum(c(b$retained, b$cached,
    8 * length(src$frame_ids) * (p + n_variables + 4) + 4 * n_reference * (p + n_variables + 1)),
    b$cap, "conditional pretest full-frame grid and simultaneous frozen-design buffers")
  invisible(NULL)
}

# Static recipe: frozen on the route's reference rows U, evaluated on all F.
cpt_grid_static <- function(src, xformla, U, date_C) {
  ids <- src$frame_ids
  vars <- all.vars(xformla)
  tf <- stats::delete.response(stats::terms(xformla))
  rw <- cpt_grid_static_raw(src, vars, date_C, U)
  if (!is.null(rw$error)) return(cpt_grid_unavailable(rw$error))
  df <- rw$df
  date_of <- function(i, v) vapply(seq_along(i), function(k) {
    du <- rw$dated_units[[v[k]]]
    if (!is.null(du) && i[k] %in% du) as.double(date_C) else NA_real_
  }, numeric(1L))
  miss <- cpt_grid_missing_table(df, ids, date_of)
  if (!is.null(miss)) return(cpt_grid_unavailable(cpt_grid_fmt_rows(miss,
    "Missing covariate values required by the conditional test grid (no carry-forward, substitution or deletion)"),
    missing = miss))
  if (!any(U)) return(cpt_grid_unavailable("The route has no reference rows on which its xformla design is frozen."))
  mf_ref <- tryCatch(stats::model.frame(tf, df[U, , drop = FALSE], na.action = stats::na.pass), error = function(e) e)
  if (inherits(mf_ref, "error")) return(cpt_grid_unavailable(paste0("The frozen xformla design could not be ",
    "reconstructed on the route's reference rows: ", conditionMessage(mf_ref))))
  tt <- attr(mf_ref, "terms")
  xlev <- stats::.getXlevels(tt, mf_ref)
  Z_ref <- stats::model.matrix(tt, mf_ref)
  cpt_grid_admit(src, ncol(Z_ref), nrow(Z_ref), length(vars))
  # xlev is keyed by model-frame terms, including factor(x), not raw
  # variable names. Inspect those evaluated columns to retain offending ids.
  mf_grid <- tryCatch(stats::model.frame(tt, df, na.action = stats::na.pass), error = function(e) e)
  if (inherits(mf_grid, "error")) return(cpt_grid_unavailable(paste0("The frozen xformla design cannot be evaluated on ",
    "the grid rows: ", conditionMessage(mf_grid))))
  bad <- cpt_grid_level_table(mf_grid, xlev, ids)
  if (!is.null(bad)) return(cpt_grid_unavailable(cpt_grid_fmt_rows(bad,
    "Grid covariate levels absent from the frozen reference design (unsupported levels; no level is added)"),
    unsupported = bad))
  ev <- cpt_grid_frozen_matrix(tt, df, xlev, attr(Z_ref, "contrasts"))
  if (!is.null(ev$error)) return(cpt_grid_unavailable(paste0("The frozen xformla design cannot be evaluated on ",
    "the grid rows: ", ev$error)))
  Z <- ev$Z
  if (!identical(dim(Z)[2L], dim(Z_ref)[2L]) || any(!is.finite(Z))) {
    nf <- which(rowSums(!is.finite(Z)) > 0L)
    return(cpt_grid_unavailable(paste0("The frozen xformla design is not finite on grid rows (ids ",
      paste(vapply(utils::head(ids[nf], 5L), chained_reference_fmt, character(1L)), collapse = ", "), ").")))
  }
  if (!cpt_grid_same(Z[U, , drop = FALSE], Z_ref))
    return(cpt_grid_unavailable(paste0("The xformla transformation state is not frozen by its terms (a data-dependent ",
      "transformation without a prediction state): grid rows cannot be evaluated without refitting it.")))
  list(status = "available", Z = matrix(as.double(Z), nrow(Z), ncol(Z)), columns = colnames(Z_ref),
    date = as.double(date_C), source = "static", reference_rows = sum(U))
}

# Dated recipe: the reader's canonical union rows are reconstructed from the
# request's reference source, their frozen state is verified against the
# retained design, and every F unit is evaluated at the recipe's date.
cpt_grid_dated <- function(src, ref, slot, design_index) {
  rs <- src$reference_source
  ids <- src$frame_ids; N <- length(ids)
  if (is.null(rs)) return(cpt_grid_unavailable("The dated recipe's reference source is not available."))
  vars <- intersect(all.vars(slot$formula), rs$prov$variable)
  tf <- stats::delete.response(stats::terms(slot$formula))
  cpt_grid_admit(src, ncol(slot$Z), nrow(slot$Z), length(vars))
  df_ref <- chained_reference_rows_frame(rs, slot$keys, vars)
  mf <- stats::model.frame(tf, df_ref, na.action = stats::na.pass)
  for (j in seq_along(mf)) {
    if (is.character(mf[[j]])) mf[[j]] <- factor(mf[[j]])
    if (is.logical(mf[[j]])) mf[[j]] <- factor(mf[[j]], levels = c(FALSE, TRUE))
  }
  tt <- attr(mf, "terms")
  Z_chk <- chained_reference_model_matrix(tt, df_ref, model_frame = mf)
  if (!cpt_grid_same(Z_chk, slot$Z))
    cpt_abort("The reconstructed frozen recipe design differs from the retained reference design.",
      "cdid_conditional_pretest_integrity_error")
  xlev <- lapply(mf, function(x) if (is.factor(x)) levels(x) else NULL)
  xlev <- xlev[!vapply(xlev, is.null, logical(1L))]
  N1 <- N + 1
  keys <- as.double(design_index) * N1 + seq_len(N)
  df <- chained_reference_rows_frame(rs, keys, vars)
  date <- if (design_index > 0L) rs$calendar[design_index] else NA_real_
  date_of <- function(i, v) {
    kind <- rs$prov$kind[match(v, rs$prov$variable)]
    ifelse(kind %in% c("dated", "date"), as.double(date), NA_real_)
  }
  miss <- cpt_grid_missing_table(df, ids, date_of)
  if (!is.null(miss)) return(cpt_grid_unavailable(cpt_grid_fmt_rows(miss,
    "Missing dated covariate values required by the conditional test grid (no carry-forward, substitution or deletion)"),
    missing = miss))
  mf_g <- tryCatch(stats::model.frame(tt, df, na.action = stats::na.pass), error = function(e) e)
  if (inherits(mf_g, "error")) return(cpt_grid_unavailable(paste0("The frozen dated design cannot be evaluated on ",
    "the grid rows: ", conditionMessage(mf_g))))
  bad <- cpt_grid_level_table(mf_g, xlev, ids)
  if (!is.null(bad)) {
    bad$date <- as.double(date)
    return(cpt_grid_unavailable(cpt_grid_fmt_rows(bad,
      "Grid covariate levels absent from the frozen reference design (unsupported levels; no level is added)"),
      unsupported = bad))
  }
  for (v in names(xlev)) mf_g[[v]] <- factor(as.character(mf_g[[v]]), levels = xlev[[v]])
  for (v in names(mf_g)) if (is.logical(mf_g[[v]])) mf_g[[v]] <- factor(mf_g[[v]], levels = c(FALSE, TRUE))
  Z <- chained_reference_model_matrix(tt, df, model_frame = mf_g)
  Z <- unname(matrix(as.double(Z), nrow(Z)))
  if (ncol(Z) != ncol(slot$Z) || any(!is.finite(Z))) {
    nf <- which(rowSums(!is.finite(Z)) > 0L)
    return(cpt_grid_unavailable(paste0("The frozen dated design is not finite on grid rows (ids ",
      paste(vapply(utils::head(ids[nf], 5L), chained_reference_fmt, character(1L)), collapse = ", "), ").")))
  }
  # Rows already on the canonical union must reproduce the retained design.
  hit <- chained_reference_lookup(keys, slot$keys)
  on <- which(!is.na(hit))
  if (length(on) && !cpt_grid_same(Z[on, , drop = FALSE], slot$Z[hit[on], , drop = FALSE]))
    return(cpt_grid_unavailable(paste0("The dated transformation state is not frozen by its terms: grid rows cannot ",
      "be evaluated without refitting it.")))
  list(status = "available", Z = Z, columns = slot$columns, date = as.double(date), source = "dated",
    reference_rows = length(slot$keys))
}

# Grid source of one request (raw columns only; no design is built here).
cpt_grid_source <- function(data, frame, idname, tname, prep) {
  id <- data[[idname]]
  list(data = data, frame = frame, frame_supplied = !is.null(frame), frame_ids = prep$frame_ids,
    unit = match(id, prep$frame_ids), tt = as.double(data[[tname]]), reference_source = prep$reference_source)
}

# Grid of one cell, cached by its recipe key. `ref` is the backend's sealed
# reference reader (NULL on static paths) and `recipe` the cell's recipe.
cpt_grid_for_cell <- function(cache, src, prep, xformla, ref, recipe, date_C) {
  slot <- if (is.null(ref)) NULL else ref$designs[[ref$role_slot["x"]]]
  dated <- !is.null(slot) && identical(slot$kind, "dated")
  key <- if (dated) paste0("dated:", ref$recipes$design_index[recipe]) else paste0("static:", sprintf("%a", as.double(date_C)))
  if (!is.null(cache[[key]])) return(list(cache = cache, grid = cache[[key]], key = key))
  if (!is.null(src$grid_budget)) src$grid_budget$cached <- sum(vapply(cache, function(x) {
    if (!identical(x$status, "available")) 0 else length(x$Z) + length(x$distinct$z) +
      length(x$distinct$first) + length(x$distinct$map) + length(x$distinct$mult)
  }, numeric(1L)))
  g <- if (dated) cpt_grid_dated(src, ref, slot, ref$recipes$design_index[recipe]) else
    cpt_grid_static(src, xformla, cpt_grid_reference_rows(prep), date_C)
  g$key <- key
  if (identical(g$status, "available")) g$distinct <- cpt_distinct_grid(g$Z)
  cache[[key]] <- g
  list(cache = cache, grid = g, key = key)
}
