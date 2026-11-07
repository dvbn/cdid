# Private panel data preparation and transition-cell primitives for the
# binary static `dname` extension.
#
# Public references:
#   - Bellego, C., Benatia, D. and Dortet-Bernadet, V. (2025),
#     "The Chained Difference-in-Differences", Journal of Econometrics 248,
#     105783, doi:10.1016/j.jeconom.2024.105783.
#   - de Chaisemartin, C., D'Haultfoeuille, X., Pasquier, F., Sow, D. and
#     Vazquez-Bare, G. (2026), "Difference-in-Differences Estimators for
#     Treatments Continuously Distributed at Every Period", arXiv:2201.06898v7
#     (static binary specialization).
#
# Scope of this module:
#   1. Long-data normalization (`static_panel_data` / `static_panel_normalize`):
#      explicit `data`, `idname`, `tname`, `yname`, `dname` (no `gname`),
#      numeric finite ordered panel dates, unique original unit IDs (numeric or
#      character, preserving type and deterministic recorded ordering),
#      duplicate `(id, date)` rejection, finite-or-NA outcome `Y`, and binary
#      `0/1` treatment `D` where known. Absent long-data rows and `NA` values in
#      `D` remain explicitly unknown (never assumed untreated `0` and never
#      carried forward or interpolated). Stores a single `N x T` representation
#      with ID and date maps; no `N x T x pairs` allocation.
#   2. Explicit transition-cell construction (`static_transition_cell`,
#      `static_transition_cells`): caller supplies target dates `t` and prior
#      treatment states `d` (or a validated table); verifies `t` and immediate
#      prior panel date `t - 1` lie in the study grid and that `D_{t-1}` and
#      `D_t` are known binary values for all `N` study units (without requiring
#      complete historical `D` at other dates). Produces full-`N` masks
#      `R_c = 1{D_{t-1} = d}`, `G_c = 1{D_{t-1} = d, D_t = 1 - d}`,
#      `C_c = 1{D_{t-1} = d, D_t = d}`, counts, sign `sigma_c = 1 - 2 d`, and
#      support diagnostics without silently dropping empty switcher/stayer cells.
#   3. Reference design extraction at `t - 1` (`static_extract_reference_design`,
#      `static_reference_design`, `static_reference_design_extractor`):
#      evaluates explicit R formulas (separate outcome/treatment and observation/
#      selection formulas allowed) at the reference date `t - 1` using
#      `stats::model.frame` and `stats::model.matrix` with `na.action = na.pass`.
#      Requires finite covariate values only on the requested fixed risk/role
#      population; rows outside the requested population remain `NA_real_`
#      without numerical evaluation. Refuses no-intercept formulas and leaves
#      rank/alias checks to downstream estimator kernels. Supports lazy callable
#      per-cell extraction.
#   4. Lazy observation-pair access (`static_observation_pair`,
#      `static_observation_pairs`, `static_observation_pair_accessor`):
#      for a fixed supplied binary role mask and caller-explicit date pairs
#      `(s, u)` with `s < u` in the study grid, returns `O = is.finite(Y_s) &
#      is.finite(Y_u)`, actual `dy = Y_u - Y_s` on `O == 1` (`NA_real_` when
#      `O == 0`, never invented zeros), counts, and per-period gap length equal
#      to the difference of panel date indices `u_idx - s_idx` (not calendar
#      spacing), preserving full original-unit alignment and returning every
#      requested pair including zero-count pairs. Supports callback/chunked and
#      on-demand access without materializing `N x all_pairs`.

static_panel_abort <- function(message, class = character(), extra = list()) {
  cond <- c(
    list(message = message, call = NULL),
    extra
  )
  stop(structure(
    cond,
    class = unique(c(class, "cdid_static_panel_error", "cdid_error", "error", "condition"))
  ))
}

static_panel_check_scalar_string <- function(x, arg_name, error_class = "cdid_static_panel_input_error") {
  if (missing(x) || is.null(x) || !is.character(x) || !is.null(dim(x)) ||
      is.object(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    static_panel_abort(
      sprintf("`%s` must be a single non-empty character string.", arg_name),
      c(error_class, "cdid_static_input_error")
    )
  }
  x
}

static_panel_check_scalar_numeric <- function(x, arg_name, error_class = "cdid_static_panel_input_error") {
  if (missing(x) || is.null(x) || !is.numeric(x) || !is.null(dim(x)) ||
      is.object(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    static_panel_abort(
      sprintf("`%s` must be a single finite numeric value.", arg_name),
      c(error_class, "cdid_static_input_error")
    )
  }
  as.double(x)
}

static_panel_validate_id_vector <- function(x, arg_name, allow_duplicates = TRUE) {
  if (missing(x) || is.null(x) || !is.atomic(x) || !is.null(dim(x)) ||
      is.object(x) || length(x) < 1L) {
    static_panel_abort(
      sprintf("`%s` must be a non-empty atomic vector of numeric or character unit identifiers (factors and objects are not allowed).", arg_name),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  is_num <- is.numeric(x)
  is_chr <- is.character(x)
  if (!(is_num || is_chr)) {
    static_panel_abort(
      sprintf("`%s` must be numeric or character (got `%s`).", arg_name, typeof(x)),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  if (anyNA(x)) {
    static_panel_abort(
      sprintf("`%s` must not contain missing (`NA` or `NaN`) unit identifiers.", arg_name),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  if (is_num && any(!is.finite(x))) {
    static_panel_abort(
      sprintf("`%s` numeric unit identifiers must all be finite.", arg_name),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  if (is_chr && any(!nzchar(x))) {
    static_panel_abort(
      sprintf("`%s` character unit identifiers must all be non-empty strings.", arg_name),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  if (is_num && is.double(x)) {
    # Normalize -0.0 to +0.0 while preserving double storage mode.
    x <- ifelse(x == 0, 0.0, x)
  }
  match_key <- if (is_num) as.double(x) else x
  chr_key <- as.character(match_key)
  if (!allow_duplicates && (anyDuplicated(match_key) > 0L || anyDuplicated(chr_key) > 0L)) {
    static_panel_abort(
      sprintf("`%s` must not contain duplicate unit identifiers.", arg_name),
      c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
    )
  }
  list(
    values = x,
    match_key = match_key,
    chr_key = chr_key,
    is_num = is_num,
    storage_mode = typeof(x)
  )
}

static_panel_validate_date_vector <- function(x, arg_name, require_ordered_unique = FALSE) {
  if (missing(x) || is.null(x) || !is.atomic(x) || !is.null(dim(x)) ||
      is.object(x) || !is.numeric(x) || length(x) < 1L) {
    static_panel_abort(
      sprintf("`%s` must be a non-empty plain numeric vector of panel dates.", arg_name),
      c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
    )
  }
  if (anyNA(x) || any(!is.finite(x))) {
    static_panel_abort(
      sprintf("`%s` must contain only finite numeric panel dates without `NA`, `NaN`, or `Inf`.", arg_name),
      c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
    )
  }
  if (is.double(x)) {
    x <- ifelse(x == 0, 0.0, x)
  }
  num_key <- as.double(x)
  if (require_ordered_unique) {
    if (length(x) < 2L) {
      static_panel_abort(
        sprintf("`%s` study grid must contain at least 2 panel dates (got %d).", arg_name, length(x)),
        c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
      )
    }
    if (anyDuplicated(num_key) > 0L) {
      static_panel_abort(
        sprintf("`%s` study grid must not contain duplicate panel dates.", arg_name),
        c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
      )
    }
    if (any(diff(num_key) <= 0)) {
      static_panel_abort(
        sprintf("`%s` study grid must be in strictly increasing numeric order.", arg_name),
        c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
      )
    }
  }
  list(
    values = x,
    num_key = num_key,
    chr_key = as.character(num_key),
    storage_mode = typeof(x)
  )
}

# 1. Long-data normalization for the binary static `dname` extension.
#
# Explicit arguments only; no `gname` argument and no default column names.
# Retains all original study units (including units with 0 or 1 observed Y) and
# keeps absent `(unit, date)` rows and `NA` values in `D` explicitly unknown
# (`NA_integer_` in `D`, distinguished by `row_present` and `d_status`).
static_panel_data <- function(data,
                              idname,
                              tname,
                              yname,
                              dname,
                              dates = NULL,
                              unit_ids = NULL) {
  if (missing(data) || is.null(data) || !is.data.frame(data)) {
    static_panel_abort(
      "`data` must be an explicit `data.frame`.",
      "cdid_static_panel_input_error"
    )
  }
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(df) < 1L) {
    static_panel_abort(
      "`data` must have at least one row.",
      "cdid_static_panel_input_error"
    )
  }

  idname <- static_panel_check_scalar_string(idname, "idname")
  tname <- static_panel_check_scalar_string(tname, "tname")
  yname <- static_panel_check_scalar_string(yname, "yname")
  dname <- static_panel_check_scalar_string(dname, "dname")

  req_cols <- c(idname, tname, yname, dname)
  if (anyDuplicated(req_cols) > 0L) {
    static_panel_abort(
      "`idname`, `tname`, `yname`, and `dname` must name four distinct columns.",
      "cdid_static_panel_input_error"
    )
  }
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0L) {
    static_panel_abort(
      sprintf("Required column(s) absent from `data`: %s.", paste(sprintf("`%s`", missing_cols), collapse = ", ")),
      "cdid_static_panel_input_error"
    )
  }

  # Validate unit IDs in `data[[idname]]` and establish canonical original-unit ordering.
  data_id_info <- static_panel_validate_id_vector(df[[idname]], sprintf("data[[%s]]", idname), allow_duplicates = TRUE)
  if (!is.null(unit_ids)) {
    caller_id_info <- static_panel_validate_id_vector(unit_ids, "unit_ids", allow_duplicates = FALSE)
    if (!identical(data_id_info$is_num, caller_id_info$is_num)) {
      static_panel_abort(
        sprintf("`unit_ids` type (`%s`) does not match `data[[%s]]` type (`%s`).",
                caller_id_info$storage_mode, idname, data_id_info$storage_mode),
        c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
      )
    }
    unit_pos <- match(data_id_info$match_key, caller_id_info$match_key)
    if (anyNA(unit_pos)) {
      bad_ids <- unique(data_id_info$values[is.na(unit_pos)])
      static_panel_abort(
        sprintf("`data[[%s]]` contains unit identifier(s) absent from explicit `unit_ids`: %s.",
                idname, paste(utils::head(as.character(bad_ids), 5L), collapse = ", ")),
        c("cdid_static_panel_id_error", "cdid_static_panel_input_error"),
        extra = list(unknown_unit_ids = bad_ids)
      )
    }
    canon_id <- caller_id_info$values
    canon_id_key <- caller_id_info$match_key
    canon_id_chr <- caller_id_info$chr_key
    unit_id_source <- "explicit_caller"
  } else {
    uniq_idx <- which(!duplicated(data_id_info$match_key))
    uniq_vals <- data_id_info$values[uniq_idx]
    uniq_keys <- data_id_info$match_key[uniq_idx]
    ord <- sort.list(uniq_keys, method = "radix")
    canon_id <- uniq_vals[ord]
    canon_id_key <- uniq_keys[ord]
    canon_id_chr <- as.character(canon_id_key)
    if (anyDuplicated(canon_id_chr) > 0L) {
      static_panel_abort(
        sprintf("`data[[%s]]` contains distinct numeric values with indistinguishable character keys.", idname),
        c("cdid_static_panel_id_error", "cdid_static_panel_input_error")
      )
    }
    unit_pos <- match(data_id_info$match_key, canon_id_key)
    unit_id_source <- "sorted_unique_from_data"
  }

  N <- length(canon_id)
  unit_first_row <- match(canon_id_key, data_id_info$match_key)
  names(unit_first_row) <- canon_id_chr
  unit_index_map <- stats::setNames(seq_len(N), canon_id_chr)

  # Validate panel dates in `data[[tname]]` and establish study date grid.
  data_t_info <- static_panel_validate_date_vector(df[[tname]], sprintf("data[[%s]]", tname), require_ordered_unique = FALSE)
  if (!is.null(dates)) {
    grid_t_info <- static_panel_validate_date_vector(dates, "dates", require_ordered_unique = TRUE)
    date_pos <- match(data_t_info$num_key, grid_t_info$num_key)
    if (anyNA(date_pos)) {
      bad_dates <- sort(unique(data_t_info$values[is.na(date_pos)]))
      static_panel_abort(
        sprintf("`data[[%s]]` contains date(s) absent from explicit `dates` study grid: %s.",
                tname, paste(as.character(bad_dates), collapse = ", ")),
        c("cdid_static_panel_date_error", "cdid_static_panel_input_error"),
        extra = list(unknown_dates = bad_dates)
      )
    }
    study_dates <- grid_t_info$values
    study_dates_num <- grid_t_info$num_key
    date_labels <- grid_t_info$chr_key
    date_grid_source <- "explicit_caller"
    explicit_dates <- TRUE
  } else {
    uniq_t_idx <- which(!duplicated(data_t_info$num_key))
    uniq_t_vals <- data_t_info$values[uniq_t_idx]
    uniq_t_keys <- data_t_info$num_key[uniq_t_idx]
    t_ord <- sort.list(uniq_t_keys, method = "radix")
    study_dates <- uniq_t_vals[t_ord]
    study_dates_num <- uniq_t_keys[t_ord]
    if (length(study_dates) < 2L) {
      static_panel_abort(
        sprintf("`data[[%s]]` must contain at least 2 distinct panel dates (found %d).", tname, length(study_dates)),
        c("cdid_static_panel_date_error", "cdid_static_panel_input_error")
      )
    }
    date_labels <- as.character(study_dates_num)
    date_pos <- match(data_t_info$num_key, study_dates_num)
    date_grid_source <- "extracted_unique_sorted"
    explicit_dates <- FALSE
  }

  T_dates <- length(study_dates)
  date_index_map <- stats::setNames(seq_len(T_dates), date_labels)

  # Reject duplicate (unit, date) records before populating N x T matrices.
  cell_linear <- (date_pos - 1L) * N + unit_pos
  dup_mask <- duplicated(cell_linear)
  if (any(dup_mask)) {
    dup_cells <- unique(cell_linear[dup_mask])
    dup_rows <- which(cell_linear %in% dup_cells)
    first_dup_row <- which(dup_mask)[1L]
    dup_id_val <- df[[idname]][first_dup_row]
    dup_t_val <- df[[tname]][first_dup_row]
    static_panel_abort(
      sprintf("Duplicate `(id, date)` record(s) found in `data` (e.g. `%s = %s`, `%s = %s` at row %d).",
              idname, as.character(dup_id_val), tname, as.character(dup_t_val), first_dup_row),
      c("cdid_static_panel_duplicate_error", "cdid_static_panel_input_error"),
      extra = list(
        duplicate_rows = dup_rows,
        duplicate_id = dup_id_val,
        duplicate_date = dup_t_val
      )
    )
  }

  # Validate outcome Y: numeric finite or NA (missing allowed; NaN/Inf/non-numeric rejected).
  raw_y <- df[[yname]]
  y_ok_type <- (is.numeric(raw_y) || (is.logical(raw_y) && all(is.na(raw_y)))) &&
    is.atomic(raw_y) && is.null(dim(raw_y)) && !is.object(raw_y)
  if (!y_ok_type) {
    static_panel_abort(
      sprintf("Outcome column `%s` must be a plain numeric vector (with `NA` allowed for missing outcomes).", yname),
      c("cdid_static_panel_outcome_error", "cdid_static_panel_input_error")
    )
  }
  if (any(is.nan(raw_y)) || any(is.infinite(raw_y))) {
    bad_y_rows <- which(is.nan(raw_y) | is.infinite(raw_y))
    static_panel_abort(
      sprintf("Outcome column `%s` contains `NaN` or `Inf` at row(s): %s. Only finite numbers or `NA` are allowed.",
              yname, paste(utils::head(bad_y_rows, 5L), collapse = ", ")),
      c("cdid_static_panel_outcome_error", "cdid_static_panel_input_error"),
      extra = list(invalid_y_rows = bad_y_rows)
    )
  }
  y_num <- as.double(raw_y)

  # Validate treatment D: binary 0/1 where known, or NA when unknown on a row.
  # Reject logical (except all-NA), character, factor, objects, NaN, Inf, and values not in {0, 1}.
  raw_d <- df[[dname]]
  d_ok_type <- (is.numeric(raw_d) || (is.logical(raw_d) && all(is.na(raw_d)))) &&
    is.atomic(raw_d) && is.null(dim(raw_d)) && !is.object(raw_d)
  if (!d_ok_type) {
    static_panel_abort(
      sprintf("Treatment column `%s` must be a plain numeric vector with binary values in `{0, 1}` (or `NA` where unknown).", dname),
      c("cdid_static_panel_treatment_error", "cdid_static_panel_input_error")
    )
  }
  if (any(is.nan(raw_d)) || any(is.infinite(raw_d))) {
    bad_d_rows <- which(is.nan(raw_d) | is.infinite(raw_d))
    static_panel_abort(
      sprintf("Treatment column `%s` contains `NaN` or `Inf` at row(s): %s.",
              dname, paste(utils::head(bad_d_rows, 5L), collapse = ", ")),
      c("cdid_static_panel_treatment_error", "cdid_static_panel_input_error"),
      extra = list(invalid_d_rows = bad_d_rows)
    )
  }
  non_na_d <- !is.na(raw_d)
  if (any(non_na_d & !(raw_d %in% c(0, 1)))) {
    bad_d_rows <- which(non_na_d & !(raw_d %in% c(0, 1)))
    bad_d_vals <- unique(raw_d[bad_d_rows])
    static_panel_abort(
      sprintf("Treatment column `%s` must be binary in `{0, 1}` where non-missing; found invalid value(s) %s at row(s) %s.",
              dname,
              paste(as.character(utils::head(bad_d_vals, 5L)), collapse = ", "),
              paste(utils::head(bad_d_rows, 5L), collapse = ", ")),
      c("cdid_static_panel_treatment_error", "cdid_static_panel_input_error"),
      extra = list(invalid_d_rows = bad_d_rows, invalid_d_values = bad_d_vals)
    )
  }
  d_int <- as.integer(raw_d)

  # Build single N x T panel representation.
  dimnames_nt <- list(unit_id = canon_id_chr, date = date_labels)
  Y_mat <- matrix(NA_real_, nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  D_mat <- matrix(NA_integer_, nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  row_present <- matrix(FALSE, nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  row_index <- matrix(NA_integer_, nrow = N, ncol = T_dates, dimnames = dimnames_nt)

  Y_mat[cell_linear] <- y_num
  D_mat[cell_linear] <- d_int
  row_present[cell_linear] <- TRUE
  row_index[cell_linear] <- seq_len(nrow(df))

  O_mat <- matrix(as.integer(is.finite(Y_mat)), nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  d_known <- matrix(!is.na(D_mat), nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  d_status <- matrix("absent_row", nrow = N, ncol = T_dates, dimnames = dimnames_nt)
  d_status[row_present & !d_known] <- "missing_d"
  d_status[d_known & D_mat == 0L] <- "untreated_0"
  d_status[d_known & D_mat == 1L] <- "treated_1"

  unit_outcome_counts <- as.integer(rowSums(O_mat))
  names(unit_outcome_counts) <- canon_id_chr
  unit_present_counts <- as.integer(rowSums(row_present))
  names(unit_present_counts) <- canon_id_chr
  unit_known_d_counts <- as.integer(rowSums(d_known))
  names(unit_known_d_counts) <- canon_id_chr

  summary_info <- list(
    n_units = N,
    n_dates = T_dates,
    n_rows_supplied = nrow(df),
    n_cells_grid = as.integer(N * T_dates),
    n_cells_present = as.integer(sum(row_present)),
    n_cells_absent = as.integer(sum(!row_present)),
    n_d_known = as.integer(sum(d_known)),
    n_d_missing_on_present_rows = as.integer(sum(row_present & !d_known)),
    n_d_unknown_total = as.integer(sum(!d_known)),
    n_y_observed = as.integer(sum(O_mat)),
    n_y_missing = as.integer(sum(O_mat == 0L)),
    n_units_zero_outcomes = as.integer(sum(unit_outcome_counts == 0L)),
    n_units_single_outcome = as.integer(sum(unit_outcome_counts == 1L)),
    n_units_two_plus_outcomes = as.integer(sum(unit_outcome_counts >= 2L)),
    unit_id_source = unit_id_source,
    unit_id_storage_mode = data_id_info$storage_mode,
    date_grid_source = date_grid_source,
    explicit_dates = explicit_dates
  )

  structure(
    list(
      n = N,
      T_dates = T_dates,
      unit_id = canon_id,
      unit_key = canon_id_chr,
      unit_match_key = canon_id_key,
      unit_id_is_num = data_id_info$is_num,
      unit_id_storage_mode = data_id_info$storage_mode,
      unit_id_source = unit_id_source,
      unit_first_row = unit_first_row,
      unit_index_map = unit_index_map,
      dates = study_dates,
      dates_num = study_dates_num,
      date_labels = date_labels,
      date_storage_mode = data_t_info$storage_mode,
      date_grid_source = date_grid_source,
      explicit_dates = explicit_dates,
      date_index_map = date_index_map,
      Y = Y_mat,
      O = O_mat,
      D = D_mat,
      row_present = row_present,
      d_known = d_known,
      d_status = d_status,
      row_index = row_index,
      unit_outcome_counts = unit_outcome_counts,
      unit_present_counts = unit_present_counts,
      unit_known_d_counts = unit_known_d_counts,
      columns = list(
        idname = idname,
        tname = tname,
        yname = yname,
        dname = dname
      ),
      data = df,
      summary = summary_info
    ),
    class = "cdid_static_panel_data"
  )
}

static_panel_normalize <- static_panel_data

static_panel_require_object <- function(panel) {
  if (missing(panel) || is.null(panel) || !inherits(panel, "cdid_static_panel_data") || !is.list(panel)) {
    static_panel_abort(
      "Expected a normalized `cdid_static_panel_data` object from `static_panel_data()`.",
      c("cdid_static_panel_input_error", "cdid_static_input_error")
    )
  }
  panel
}

# 2. Explicit transition-cell constructor for a single cell c = (t, d).
#
# Verifies that target date `t` and immediate prior date `t_prev` are in the
# study grid and that `D_{t-1}` and `D_t` are known binary values for all `N`
# original study units. Empty switcher or stayer counts are reported as support
# diagnostics without silently dropping or altering the requested cell.
static_transition_cell <- function(panel,
                                   t,
                                   d,
                                   cell_id = NULL,
                                   t_prev = NULL,
                                   sigma = NULL) {
  panel <- static_panel_require_object(panel)

  t_num <- static_panel_check_scalar_numeric(t, "t", error_class = "cdid_static_cell_date_error")
  t_idx <- match(t_num, panel$dates_num)
  if (is.na(t_idx)) {
    static_panel_abort(
      sprintf("Requested target date `t = %s` is not in the study date grid [%s].",
              as.character(t_num), paste(panel$date_labels, collapse = ", ")),
      c("cdid_static_cell_date_error", "cdid_static_cell_error"),
      extra = list(requested_t = t, study_dates = panel$dates)
    )
  }
  if (t_idx < 2L) {
    static_panel_abort(
      sprintf("Requested target date `t = %s` is the initial study grid date (`t_idx = 1`); no immediate prior panel date exists in the study grid.",
              panel$date_labels[t_idx]),
      c("cdid_static_cell_date_error", "cdid_static_cell_error"),
      extra = list(requested_t = t, t_idx = t_idx, study_dates = panel$dates)
    )
  }
  prev_idx <- t_idx - 1L
  actual_t <- panel$dates[t_idx]
  actual_t_prev <- panel$dates[prev_idx]

  if (!is.null(t_prev)) {
    t_prev_num <- static_panel_check_scalar_numeric(t_prev, "t_prev", error_class = "cdid_static_cell_date_error")
    if (is.na(match(t_prev_num, panel$dates_num))) {
      static_panel_abort(
        sprintf("Supplied prior date `t_prev = %s` is not in the study date grid [%s].",
                as.character(t_prev_num), paste(panel$date_labels, collapse = ", ")),
        c("cdid_static_cell_date_error", "cdid_static_cell_error"),
        extra = list(requested_t_prev = t_prev, study_dates = panel$dates)
      )
    }
    if (t_prev_num != panel$dates_num[prev_idx]) {
      static_panel_abort(
        sprintf("Supplied prior date `t_prev = %s` is not the immediate prior panel date (`%s`) for target date `t = %s`.",
                as.character(t_prev_num), panel$date_labels[prev_idx], panel$date_labels[t_idx]),
        c("cdid_static_cell_date_error", "cdid_static_cell_error"),
        extra = list(
          requested_t = actual_t,
          requested_t_prev = t_prev,
          expected_t_prev = actual_t_prev
        )
      )
    }
  }

  d_num <- static_panel_check_scalar_numeric(d, "d", error_class = "cdid_static_cell_state_error")
  if (!(d_num %in% c(0, 1))) {
    static_panel_abort(
      sprintf("Prior treatment state `d` must be binary in `{0, 1}` (got `%s`).", as.character(d)),
      c("cdid_static_cell_state_error", "cdid_static_cell_error")
    )
  }
  d_val <- as.integer(d_num)
  sigma_val <- as.double(1L - 2L * d_val)

  if (!is.null(sigma)) {
    sigma_num <- static_panel_check_scalar_numeric(sigma, "sigma", error_class = "cdid_static_cell_state_error")
    if (!(sigma_num %in% c(-1, 1)) || sigma_num != sigma_val) {
      static_panel_abort(
        sprintf("Supplied `sigma = %s` is inconsistent with `d = %d` (require `sigma = 1 - 2 * d = %g`).",
                as.character(sigma), d_val, sigma_val),
        c("cdid_static_cell_state_error", "cdid_static_cell_error")
      )
    }
  }

  if (is.null(cell_id)) {
    cell_id <- sprintf("t=%s,d=%d", panel$date_labels[t_idx], d_val)
  } else {
    cell_id <- static_panel_check_scalar_string(cell_id, "cell_id", error_class = "cdid_static_cell_input_error")
  }

  # Verify required treatment statuses at t_prev (prev_idx) and t (t_idx) are known
  # for every unit of the original study population.
  D_prev <- panel$D[, prev_idx]
  D_curr <- panel$D[, t_idx]
  pres_prev <- panel$row_present[, prev_idx]
  pres_curr <- panel$row_present[, t_idx]

  if (anyNA(D_prev) || anyNA(D_curr)) {
    absent_prev_units <- panel$unit_id[!pres_prev]
    na_prev_units <- panel$unit_id[pres_prev & is.na(D_prev)]
    absent_curr_units <- panel$unit_id[!pres_curr]
    na_curr_units <- panel$unit_id[pres_curr & is.na(D_curr)]
    missing_prev_units <- panel$unit_id[is.na(D_prev)]
    missing_curr_units <- panel$unit_id[is.na(D_curr)]
    missing_units <- panel$unit_id[is.na(D_prev) | is.na(D_curr)]

    details <- character(0L)
    if (length(absent_prev_units) > 0L) {
      details <- c(details, sprintf("absent row at t_prev=%s for unit(s) [%s]",
                                    panel$date_labels[prev_idx],
                                    paste(utils::head(as.character(absent_prev_units), 5L), collapse = ", ")))
    }
    if (length(na_prev_units) > 0L) {
      details <- c(details, sprintf("missing D at t_prev=%s for unit(s) [%s]",
                                    panel$date_labels[prev_idx],
                                    paste(utils::head(as.character(na_prev_units), 5L), collapse = ", ")))
    }
    if (length(absent_curr_units) > 0L) {
      details <- c(details, sprintf("absent row at t=%s for unit(s) [%s]",
                                    panel$date_labels[t_idx],
                                    paste(utils::head(as.character(absent_curr_units), 5L), collapse = ", ")))
    }
    if (length(na_curr_units) > 0L) {
      details <- c(details, sprintf("missing D at t=%s for unit(s) [%s]",
                                    panel$date_labels[t_idx],
                                    paste(utils::head(as.character(na_curr_units), 5L), collapse = ", ")))
    }

    static_panel_abort(
      sprintf(
        "Transition cell `%s` (`%s -> %s`, `d = %d`) cannot be constructed: required treatment status `D` is unknown for %d of %d study unit(s) (%s). Absent rows and missing `D` are never assumed untreated or interpolated.",
        cell_id,
        panel$date_labels[prev_idx],
        panel$date_labels[t_idx],
        d_val,
        length(missing_units),
        panel$n,
        paste(details, collapse = "; ")
      ),
      c("cdid_static_missing_treatment_error", "cdid_static_cell_error"),
      extra = list(
        cell_id = cell_id,
        t = actual_t,
        t_prev = actual_t_prev,
        t_idx = t_idx,
        t_prev_idx = prev_idx,
        d = d_val,
        missing_units = missing_units,
        missing_prev_units = missing_prev_units,
        missing_curr_units = missing_curr_units,
        absent_prev_units = absent_prev_units,
        na_prev_units = na_prev_units,
        absent_curr_units = absent_curr_units,
        na_curr_units = na_curr_units
      )
    )
  }

  # Fixed full-N binary masks in original-unit order.
  R_mask <- as.integer(D_prev == d_val)
  G_mask <- as.integer(D_prev == d_val & D_curr == (1L - d_val))
  C_mask <- as.integer(D_prev == d_val & D_curr == d_val)

  n_units <- panel$n
  n_risk <- as.integer(sum(R_mask))
  n_switchers <- as.integer(sum(G_mask))
  n_stayers <- as.integer(sum(C_mask))

  # Informational adjacent outcome availability (never alters R, G, or C).
  O_adj <- as.integer(panel$O[, prev_idx] == 1L & panel$O[, t_idx] == 1L)
  n_adj_obs_all <- as.integer(sum(O_adj))
  n_adj_obs_risk <- as.integer(sum(R_mask == 1L & O_adj == 1L))
  n_adj_obs_switchers <- as.integer(sum(G_mask == 1L & O_adj == 1L))
  n_adj_obs_stayers <- as.integer(sum(C_mask == 1L & O_adj == 1L))

  support_status <- if (n_risk == 0L) {
    "empty_risk"
  } else if (n_switchers == 0L && n_stayers == 0L) {
    "empty_risk"
  } else if (n_switchers == 0L) {
    "empty_switchers"
  } else if (n_stayers == 0L) {
    "empty_stayers"
  } else {
    "supported"
  }

  direction <- if (d_val == 0L) "entry" else "exit"
  label <- sprintf("c=(%s,%d)", panel$date_labels[t_idx], d_val)
  transition_label <- sprintf("%s->%s (d:%d->%d)", panel$date_labels[prev_idx], panel$date_labels[t_idx], d_val, 1L - d_val)

  support <- list(
    n_units = n_units,
    n_risk = n_risk,
    n_switchers = n_switchers,
    n_stayers = n_stayers,
    has_risk = (n_risk > 0L),
    has_switchers = (n_switchers > 0L),
    has_stayers = (n_stayers > 0L),
    is_supported = (n_risk > 0L && n_switchers > 0L && n_stayers > 0L),
    empty_risk = (n_risk == 0L),
    empty_switchers = (n_switchers == 0L),
    empty_stayers = (n_stayers == 0L),
    status = support_status,
    n_adjacent_observed_all = n_adj_obs_all,
    n_adjacent_observed_risk = n_adj_obs_risk,
    n_adjacent_observed_switchers = n_adj_obs_switchers,
    n_adjacent_observed_stayers = n_adj_obs_stayers
  )

  structure(
    list(
      cell_id = cell_id,
      label = label,
      transition_label = transition_label,
      direction = direction,
      t = actual_t,
      t_prev = actual_t_prev,
      t_num = panel$dates_num[t_idx],
      t_prev_num = panel$dates_num[prev_idx],
      t_label = panel$date_labels[t_idx],
      t_prev_label = panel$date_labels[prev_idx],
      t_idx = as.integer(t_idx),
      t_prev_idx = as.integer(prev_idx),
      d = d_val,
      d_prev = d_val,
      d_curr_switcher = as.integer(1L - d_val),
      d_curr_stayer = d_val,
      sigma = sigma_val,
      n = n_units,
      unit_id = panel$unit_id,
      unit_key = panel$unit_key,
      unit_id_is_num = panel$unit_id_is_num,
      R = R_mask,
      G = G_mask,
      C = C_mask,
      n_risk = n_risk,
      n_switchers = n_switchers,
      n_stayers = n_stayers,
      pi_hat = n_switchers / as.double(n_units),
      q0_hat = n_stayers / as.double(n_units),
      risk_share = n_risk / as.double(n_units),
      support = support
    ),
    class = "cdid_static_transition_cell"
  )
}

static_panel_transition_cell <- static_transition_cell

# Construct an ordered collection of explicit requested transition cells.
# Accepts either:
#   - `cells`: a data.frame or list of requested cell specs, OR
#   - explicit parallel vectors `t` and `d` (with optional `cell_id`, `t_prev`).
# Never infers target cells from observed switches and never drops unsupported
# requested cells.
static_transition_cells <- function(panel,
                                    cells = NULL,
                                    t = NULL,
                                    d = NULL,
                                    cell_id = NULL,
                                    t_prev = NULL) {
  panel <- static_panel_require_object(panel)

  has_cells <- !is.null(cells)
  has_td <- !is.null(t) || !is.null(d)
  if (has_cells && has_td) {
    static_panel_abort(
      "Supply either `cells` (table/list) or explicit `t` and `d` vectors, not both.",
      c("cdid_static_cell_input_error", "cdid_static_cell_error")
    )
  }
  if (!has_cells && !has_td) {
    static_panel_abort(
      "Supply explicit requested transition cells via `cells` or (`t`, `d`); cells are never chosen automatically.",
      c("cdid_static_cell_input_error", "cdid_static_cell_error")
    )
  }

  cell_list <- NULL
  if (has_td) {
    if (is.null(t) || is.null(d) || !is.atomic(t) || !is.atomic(d) ||
        !is.null(dim(t)) || !is.null(dim(d)) || is.object(t) || is.object(d) ||
        !is.numeric(t) || !is.numeric(d) || length(t) < 1L || length(t) != length(d)) {
      static_panel_abort(
        "`t` and `d` must be non-empty numeric vectors of equal length.",
        c("cdid_static_cell_input_error", "cdid_static_cell_error")
      )
    }
    n_req <- length(t)
    if (!is.null(cell_id)) {
      if (!is.character(cell_id) || !is.null(dim(cell_id)) || is.object(cell_id) ||
          length(cell_id) != n_req || anyNA(cell_id) || any(!nzchar(cell_id))) {
        static_panel_abort(
          sprintf("`cell_id` must be NULL or a non-empty character vector of length %d.", n_req),
          c("cdid_static_cell_input_error", "cdid_static_cell_error")
        )
      }
    }
    if (!is.null(t_prev)) {
      if (!is.numeric(t_prev) || !is.null(dim(t_prev)) || is.object(t_prev) ||
          length(t_prev) != n_req) {
        static_panel_abort(
          sprintf("`t_prev` must be NULL or a numeric vector of length %d.", n_req),
          c("cdid_static_cell_input_error", "cdid_static_cell_error")
        )
      }
    }
    cell_list <- lapply(seq_len(n_req), function(k) {
      static_transition_cell(
        panel = panel,
        t = t[k],
        d = d[k],
        cell_id = if (!is.null(cell_id)) cell_id[k] else NULL,
        t_prev = if (!is.null(t_prev)) t_prev[k] else NULL
      )
    })
  } else if (is.data.frame(cells)) {
    df_c <- as.data.frame(cells, stringsAsFactors = FALSE)
    if (nrow(df_c) < 1L) {
      static_panel_abort(
        "`cells` data.frame must have at least one row.",
        c("cdid_static_cell_input_error", "cdid_static_cell_error")
      )
    }
    t_col <- if ("t" %in% names(df_c)) "t" else if ("target_date" %in% names(df_c)) "target_date" else NULL
    d_col <- if ("d" %in% names(df_c)) "d" else if ("prior_state" %in% names(df_c)) "prior_state" else if ("d_prev" %in% names(df_c)) "d_prev" else NULL
    if (is.null(t_col) || is.null(d_col)) {
      static_panel_abort(
        "`cells` data.frame must contain target date column `t` (or `target_date`) and prior state column `d` (or `prior_state` / `d_prev`).",
        c("cdid_static_cell_input_error", "cdid_static_cell_error")
      )
    }
    cid_col <- if ("cell_id" %in% names(df_c)) "cell_id" else NULL
    tprev_col <- if ("t_prev" %in% names(df_c)) "t_prev" else if ("ref_date" %in% names(df_c)) "ref_date" else NULL
    sig_col <- if ("sigma" %in% names(df_c)) "sigma" else NULL

    cell_list <- lapply(seq_len(nrow(df_c)), function(k) {
      static_transition_cell(
        panel = panel,
        t = df_c[[t_col]][k],
        d = df_c[[d_col]][k],
        cell_id = if (!is.null(cid_col)) as.character(df_c[[cid_col]][k]) else NULL,
        t_prev = if (!is.null(tprev_col)) df_c[[tprev_col]][k] else NULL,
        sigma = if (!is.null(sig_col)) df_c[[sig_col]][k] else NULL
      )
    })
  } else if (is.list(cells) && !is.object(cells) && length(cells) >= 1L) {
    list_names <- names(cells)
    cell_list <- lapply(seq_along(cells), function(k) {
      el <- cells[[k]]
      fb_id <- if (!is.null(list_names) && !is.na(list_names[k]) && nzchar(list_names[k])) list_names[k] else NULL
      if (inherits(el, "cdid_static_transition_cell")) {
        if (!is.null(fb_id) && !identical(fb_id, el$cell_id)) {
          static_panel_abort(
            sprintf("List name `%s` does not match `cell$cell_id` (`%s`).", fb_id, el$cell_id),
            c("cdid_static_cell_input_error", "cdid_static_cell_error")
          )
        }
        return(static_panel_resolve_cell(panel, el))
      }
      if (!is.list(el) || is.object(el)) {
        static_panel_abort(
          "Every element of `cells` list must be a `cdid_static_transition_cell` or a plain list with `t` and `d`.",
          c("cdid_static_cell_input_error", "cdid_static_cell_error")
        )
      }
      t_val <- if (!is.null(el$t)) el$t else el$target_date
      d_val <- if (!is.null(el$d)) el$d else if (!is.null(el$prior_state)) el$prior_state else el$d_prev
      cid_val <- if (!is.null(el$cell_id)) el$cell_id else fb_id
      tprev_val <- if (!is.null(el$t_prev)) el$t_prev else el$ref_date
      static_transition_cell(
        panel = panel,
        t = t_val,
        d = d_val,
        cell_id = cid_val,
        t_prev = tprev_val,
        sigma = el$sigma
      )
    })
  } else {
    static_panel_abort(
      "`cells` must be a non-empty `data.frame` or list of cell specifications.",
      c("cdid_static_cell_input_error", "cdid_static_cell_error")
    )
  }

  cids <- vapply(cell_list, `[[`, character(1L), "cell_id")
  if (anyDuplicated(cids) > 0L) {
    dup_id <- cids[duplicated(cids)][1L]
    static_panel_abort(
      sprintf("Duplicate `cell_id` (`%s`) in requested transition cells.", dup_id),
      c("cdid_static_cell_input_error", "cdid_static_cell_error")
    )
  }
  td_keys <- vapply(cell_list, function(cl) sprintf("%d:%d", cl$t_idx, cl$d), character(1L))
  if (anyDuplicated(td_keys) > 0L) {
    dup_idx <- which(duplicated(td_keys))[1L]
    dup_cl <- cell_list[[dup_idx]]
    static_panel_abort(
      sprintf("Duplicate transition cell `(t = %s, d = %d)` requested.", dup_cl$t_label, dup_cl$d),
      c("cdid_static_cell_input_error", "cdid_static_cell_error")
    )
  }
  names(cell_list) <- cids

  support_table <- data.frame(
    cell_id = cids,
    label = vapply(cell_list, `[[`, character(1L), "label"),
    direction = vapply(cell_list, `[[`, character(1L), "direction"),
    t_prev = vapply(cell_list, `[[`, numeric(1L), "t_prev_num"),
    t = vapply(cell_list, `[[`, numeric(1L), "t_num"),
    t_prev_idx = vapply(cell_list, `[[`, integer(1L), "t_prev_idx"),
    t_idx = vapply(cell_list, `[[`, integer(1L), "t_idx"),
    d = vapply(cell_list, `[[`, integer(1L), "d"),
    sigma = vapply(cell_list, `[[`, numeric(1L), "sigma"),
    n_units = vapply(cell_list, `[[`, integer(1L), "n"),
    n_risk = vapply(cell_list, `[[`, integer(1L), "n_risk"),
    n_switchers = vapply(cell_list, `[[`, integer(1L), "n_switchers"),
    n_stayers = vapply(cell_list, `[[`, integer(1L), "n_stayers"),
    has_risk = vapply(cell_list, function(cl) cl$support$has_risk, logical(1L)),
    has_switchers = vapply(cell_list, function(cl) cl$support$has_switchers, logical(1L)),
    has_stayers = vapply(cell_list, function(cl) cl$support$has_stayers, logical(1L)),
    is_supported = vapply(cell_list, function(cl) cl$support$is_supported, logical(1L)),
    status = vapply(cell_list, function(cl) cl$support$status, character(1L)),
    n_adjacent_observed_risk = vapply(cell_list, function(cl) cl$support$n_adjacent_observed_risk, integer(1L)),
    n_adjacent_observed_switchers = vapply(cell_list, function(cl) cl$support$n_adjacent_observed_switchers, integer(1L)),
    n_adjacent_observed_stayers = vapply(cell_list, function(cl) cl$support$n_adjacent_observed_stayers, integer(1L)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  structure(
    list(
      cells = cell_list,
      cell_ids = cids,
      n_cells = length(cell_list),
      n = panel$n,
      unit_id = panel$unit_id,
      support_table = support_table,
      supported_cell_ids = cids[support_table$is_supported],
      unsupported_cell_ids = cids[!support_table$is_supported],
      empty_risk_cell_ids = cids[support_table$n_risk == 0L],
      empty_switcher_cell_ids = cids[support_table$n_switchers == 0L],
      empty_stayer_cell_ids = cids[support_table$n_stayers == 0L]
    ),
    class = "cdid_static_transition_cells"
  )
}

static_panel_transition_cells <- static_transition_cells

# Internal helper: validate a formula and build known panel factor levels `xlev`
# without evaluating irrelevant rows.
static_panel_prepare_formula <- function(panel, formula, role = "outcome", xlev = NULL) {
  if (missing(formula) || is.null(formula) || !inherits(formula, "formula")) {
    static_panel_abort(
      sprintf("Supply an explicit R `formula` for `%s` design extraction; no default formula is assumed.", role),
      c("cdid_static_formula_error", "cdid_static_design_error")
    )
  }

  tt_raw <- tryCatch(
    stats::terms(formula, data = panel$data[integer(0L), , drop = FALSE]),
    error = function(e) {
      static_panel_abort(
        sprintf("Invalid `%s` formula `%s`: %s", role, paste(deparse(formula), collapse = " "), conditionMessage(e)),
        c("cdid_static_formula_error", "cdid_static_design_error")
      )
    }
  )
  had_response <- identical(attr(tt_raw, "response"), 1L)
  tt <- if (had_response) stats::delete.response(tt_raw) else tt_raw

  if (!identical(attr(tt, "intercept"), 1L)) {
    static_panel_abort(
      sprintf(
        "Refusing `%s` formula `%s` without an intercept: static switcher/stayer theory requires local cell/pair intercepts and never silently adds one.",
        role, paste(deparse(formula), collapse = " ")
      ),
      c("cdid_static_no_intercept_error", "cdid_static_formula_error", "cdid_static_design_error")
    )
  }

  vars <- all.vars(tt)
  missing_vars <- setdiff(vars, names(panel$data))
  if (length(missing_vars) > 0L) {
    static_panel_abort(
      sprintf("`%s` formula `%s` references variable(s) absent from `data`: %s.",
              role, paste(deparse(formula), collapse = " "),
              paste(sprintf("`%s`", missing_vars), collapse = ", ")),
      c("cdid_static_formula_error", "cdid_static_design_error"),
      extra = list(missing_variables = missing_vars, role = role)
    )
  }

  # Collect deterministic known factor levels from `panel$data` so subsets with
  # fewer observed levels still preserve all known factor levels and column names.
  panel_xlev <- list()
  for (v in vars) {
    col <- panel$data[[v]]
    if (is.factor(col)) {
      panel_xlev[[v]] <- levels(col)
    } else if (is.character(col)) {
      non_na <- col[!is.na(col)]
      panel_xlev[[v]] <- sort.int(unique(non_na), method = "radix")
    }
  }

  # Also handle inline `factor(v)` or `as.factor(v)` expressions in formula terms.
  var_exprs <- attr(tt, "variables")
  if (is.language(var_exprs) && length(var_exprs) > 1L) {
    for (j in seq.int(2L, length(var_exprs))) {
      expr_j <- var_exprs[[j]]
      expr_label <- paste(deparse(expr_j, width.cutoff = 500L), collapse = " ")
      if (is.call(expr_j) && length(expr_j) >= 2L) {
        fn_name <- paste(deparse(expr_j[[1L]]), collapse = "")
        if (fn_name %in% c("factor", "as.factor", "base::factor", "base::as.factor")) {
          arg_vars <- all.vars(expr_j[[2L]])
          if (length(arg_vars) == 1L && arg_vars %in% names(panel$data) && is.null(names(expr_j))) {
            col <- panel$data[[arg_vars]]
            if (is.factor(col)) {
              panel_xlev[[expr_label]] <- levels(col)
            } else {
              non_na <- col[!is.na(col)]
              panel_xlev[[expr_label]] <- as.character(sort.int(unique(non_na), method = "radix"))
            }
          }
        }
      }
    }
  }

  if (!is.null(xlev)) {
    if (!is.list(xlev) || is.object(xlev) || is.null(names(xlev)) || any(!nzchar(names(xlev)))) {
      static_panel_abort(
        "`xlev`, when supplied, must be a named list of character factor levels.",
        c("cdid_static_formula_error", "cdid_static_design_error")
      )
    }
    for (nm in names(xlev)) {
      panel_xlev[[nm]] <- as.character(xlev[[nm]])
    }
  }

  list(
    formula = formula,
    formula_text = paste(deparse(formula, width.cutoff = 500L), collapse = " "),
    terms = tt,
    had_response = had_response,
    variables = vars,
    term_labels = attr(tt, "term.labels"),
    xlev = panel_xlev
  )
}

static_panel_resolve_cell <- function(panel, cell) {
  if (inherits(cell, "cdid_static_transition_cell")) {
    cell_uid <- cell$unit_id
    cell_is_num <- !is.null(cell_uid) && is.atomic(cell_uid) && is.null(dim(cell_uid)) &&
      !is.object(cell_uid) && length(cell_uid) >= 1L && !anyNA(cell_uid) &&
      is.numeric(cell_uid) && all(is.finite(cell_uid))
    cell_is_chr <- !is.null(cell_uid) && is.atomic(cell_uid) && is.null(dim(cell_uid)) &&
      !is.object(cell_uid) && length(cell_uid) >= 1L && !anyNA(cell_uid) &&
      is.character(cell_uid) && all(nzchar(cell_uid))
    cell_match_key <- if (cell_is_num) {
      as.double(ifelse(cell_uid == 0, 0.0, cell_uid))
    } else if (cell_is_chr) {
      cell_uid
    } else {
      NULL
    }
    same_uid <- (cell_is_num || cell_is_chr) &&
      identical(cell_is_num, panel$unit_id_is_num) &&
      identical(length(cell_uid), panel$n) &&
      identical(cell_match_key, panel$unit_match_key) &&
      identical(cell$n, panel$n) &&
      identical(cell$unit_id_is_num, panel$unit_id_is_num) &&
      identical(cell$unit_key, panel$unit_key)
    if (!same_uid) {
      static_panel_abort(
        sprintf("Cell `%s` does not match `panel$unit_id` in exact unit identity and order.",
                as.character(cell$cell_id)[1L]),
        c("cdid_static_cell_input_error", "cdid_static_panel_id_error", "cdid_static_design_error", "cdid_static_cell_error")
      )
    }
    if (is.null(cell$cell_id) || is.null(cell$t) || is.null(cell$t_prev) ||
        is.null(cell$d) || is.null(cell$sigma)) {
      static_panel_abort(
        "Supplied transition cell is missing required transition metadata (`cell_id`, `t`, `t_prev`, `d`, `sigma`).",
        c("cdid_static_cell_input_error", "cdid_static_design_error", "cdid_static_cell_error")
      )
    }
    fresh <- static_transition_cell(
      panel = panel,
      t = cell$t,
      d = cell$d,
      cell_id = cell$cell_id,
      t_prev = cell$t_prev,
      sigma = cell$sigma
    )
    same_transition <- identical(cell$cell_id, fresh$cell_id) &&
      identical(cell$label, fresh$label) &&
      identical(cell$transition_label, fresh$transition_label) &&
      identical(cell$direction, fresh$direction) &&
      identical(as.double(cell$t), as.double(fresh$t)) &&
      identical(as.double(cell$t_prev), as.double(fresh$t_prev)) &&
      identical(cell$t_num, fresh$t_num) &&
      identical(cell$t_prev_num, fresh$t_prev_num) &&
      identical(cell$t_label, fresh$t_label) &&
      identical(cell$t_prev_label, fresh$t_prev_label) &&
      identical(cell$t_idx, fresh$t_idx) &&
      identical(cell$t_prev_idx, fresh$t_prev_idx) &&
      identical(cell$d, fresh$d) &&
      identical(cell$d_prev, fresh$d_prev) &&
      identical(cell$d_curr_switcher, fresh$d_curr_switcher) &&
      identical(cell$d_curr_stayer, fresh$d_curr_stayer) &&
      identical(cell$sigma, fresh$sigma) &&
      identical(cell$R, fresh$R) &&
      identical(cell$G, fresh$G) &&
      identical(cell$C, fresh$C) &&
      identical(cell$n_risk, fresh$n_risk) &&
      identical(cell$n_switchers, fresh$n_switchers) &&
      identical(cell$n_stayers, fresh$n_stayers) &&
      identical(cell$pi_hat, fresh$pi_hat) &&
      identical(cell$q0_hat, fresh$q0_hat) &&
      identical(cell$risk_share, fresh$risk_share) &&
      identical(cell$support, fresh$support)
    if (!same_transition) {
      static_panel_abort(
        sprintf(
          "Cell `%s` has stale or incompatible date indices, reference date, masks, or support relative to `panel`.",
          fresh$cell_id
        ),
        c("cdid_static_cell_input_error", "cdid_static_design_error", "cdid_static_cell_error")
      )
    }
    return(fresh)
  }
  if (is.list(cell) && !is.object(cell) && !is.null(cell$t) && !is.null(cell$d)) {
    return(static_transition_cell(
      panel = panel,
      t = cell$t,
      d = cell$d,
      cell_id = cell$cell_id,
      t_prev = cell$t_prev,
      sigma = cell$sigma
    ))
  }
  static_panel_abort(
    "`cell` must be a `cdid_static_transition_cell` object (or a list with `t` and `d`).",
    c("cdid_static_cell_input_error", "cdid_static_design_error")
  )
}

static_panel_resolve_population_mask <- function(panel, cell, population = "risk", mask = NULL) {
  n <- panel$n
  if (!is.null(mask)) {
    if (!is.atomic(mask) || !is.null(dim(mask)) || is.object(mask) ||
        !(is.numeric(mask) || is.logical(mask)) || length(mask) != n || anyNA(mask)) {
      static_panel_abort(
        sprintf("`mask` must be a non-missing binary vector of length `n = %d`.", n),
        "cdid_static_design_error"
      )
    }
    if (!is.null(names(mask)) && !identical(as.character(names(mask)), panel$unit_key)) {
      static_panel_abort(
        "Names of `mask` do not match `panel$unit_id` in the exact same order.",
        "cdid_static_design_error"
      )
    }
    if (is.numeric(mask) && any(!is.finite(mask) | !(mask %in% c(0, 1)))) {
      static_panel_abort(
        "`mask` numeric entries must lie in `{0, 1}`.",
        "cdid_static_design_error"
      )
    }
    return(list(
      population = "custom",
      mask = as.integer(mask)
    ))
  }

  pop <- static_panel_check_scalar_string(population, "population", error_class = "cdid_static_design_error")
  if (pop %in% c("risk", "R")) {
    list(population = "risk", mask = cell$R)
  } else if (pop %in% c("switcher", "switchers", "G")) {
    list(population = "switcher", mask = cell$G)
  } else if (pop %in% c("stayer", "stayers", "C")) {
    list(population = "stayer", mask = cell$C)
  } else if (pop %in% c("all", "full")) {
    list(population = "all", mask = rep(1L, n))
  } else {
    static_panel_abort(
      sprintf("`population` must be one of `\"risk\"`, `\"switcher\"`, `\"stayer\"`, or `\"all\"` (got `%s`).", pop),
      "cdid_static_design_error"
    )
  }
}

# 3. Reference design extraction for a single formula at cell reference date t - 1.
#
# Evaluates `formula` via `stats::model.frame` and `stats::model.matrix` ONLY on
# the required fixed population rows (`mask == 1L`) at reference date `t - 1`.
# Irrelevant rows (`mask == 0L`) remain `NA_real_` without numerical evaluation.
static_extract_reference_design <- function(panel,
                                            cell,
                                            formula,
                                            role = "outcome",
                                            population = "risk",
                                            mask = NULL,
                                            xlev = NULL,
                                            prepared_formula = NULL) {
  panel <- static_panel_require_object(panel)
  cell <- static_panel_resolve_cell(panel, cell)
  role <- static_panel_check_scalar_string(role, "role", error_class = "cdid_static_design_error")
  pop_info <- static_panel_resolve_population_mask(panel, cell, population = population, mask = mask)
  req_mask <- pop_info$mask
  pop_name <- pop_info$population

  f_info <- if (!is.null(prepared_formula)) {
    prepared_formula
  } else {
    static_panel_prepare_formula(panel, formula = formula, role = role, xlev = xlev)
  }

  n <- panel$n
  ref_idx <- cell$t_prev_idx
  req_units_idx <- which(req_mask == 1L)
  n_req <- length(req_units_idx)

  # Check that every required unit has a present record in `data` at t - 1.
  ref_rows_all <- rep(NA_integer_, n)
  if (n_req > 0L) {
    ref_rows_req <- panel$row_index[req_units_idx, ref_idx]
    if (anyNA(ref_rows_req)) {
      miss_u <- panel$unit_id[req_units_idx[is.na(ref_rows_req)]]
      static_panel_abort(
        sprintf(
          "Cell `%s` `%s` design at reference date `t_prev = %s` is missing long-data row(s) for %d required `%s` unit(s): [%s].",
          cell$cell_id, role, cell$t_prev_label, length(miss_u), pop_name,
          paste(utils::head(as.character(miss_u), 5L), collapse = ", ")
        ),
        c("cdid_static_missing_reference_row_error", "cdid_static_design_error"),
        extra = list(
          cell_id = cell$cell_id,
          role = role,
          population = pop_name,
          t_prev = cell$t_prev,
          missing_units = miss_u
        )
      )
    }
    ref_rows_all[req_units_idx] <- ref_rows_req
  } else {
    ref_rows_req <- integer(0L)
  }

  # Subset ONLY the required reference rows at date t - 1 so irrelevant rows
  # outside `req_mask` are never evaluated numerically.
  sub_df <- panel$data[ref_rows_req, f_info$variables, drop = FALSE]

  # Check raw referenced variables on required rows before formula evaluation,
  # and apply known panel factor levels to character/factor columns.
  for (v in f_info$variables) {
    col_v <- sub_df[[v]]
    if (n_req > 0L) {
      bad_v <- is.na(col_v) | (is.numeric(col_v) & !is.finite(col_v))
      if (any(bad_v)) {
        bad_u <- panel$unit_id[req_units_idx[which(bad_v)]]
        static_panel_abort(
          sprintf(
            "Cell `%s` `%s` design at reference date `t_prev = %s` has missing or non-finite value(s) in variable `%s` for %d required `%s` unit(s): [%s].",
            cell$cell_id, role, cell$t_prev_label, v, length(bad_u), pop_name,
            paste(utils::head(as.character(bad_u), 5L), collapse = ", ")
          ),
          c("cdid_static_missing_covariate_error", "cdid_static_design_error"),
          extra = list(
            cell_id = cell$cell_id,
            role = role,
            population = pop_name,
            t_prev = cell$t_prev,
            variable = v,
            missing_units = bad_u
          )
        )
      }
    }
    if (is.character(col_v) && !is.null(f_info$xlev[[v]])) {
      sub_df[[v]] <- factor(col_v, levels = f_info$xlev[[v]])
    } else if (is.factor(col_v) && !is.null(f_info$xlev[[v]])) {
      sub_df[[v]] <- factor(col_v, levels = f_info$xlev[[v]], ordered = is.ordered(col_v))
    }
  }

  xlev_use <- f_info$xlev[intersect(names(f_info$xlev), c(f_info$variables, rownames(attr(f_info$terms, "factors"))))]
  mf <- tryCatch(
    stats::model.frame(
      f_info$terms,
      data = sub_df,
      na.action = stats::na.pass,
      drop.unused.levels = FALSE,
      xlev = if (length(xlev_use) > 0L) xlev_use else NULL
    ),
    error = function(e) {
      static_panel_abort(
        sprintf("Cell `%s` `%s` `model.frame` failed at `t_prev = %s`: %s",
                cell$cell_id, role, cell$t_prev_label, conditionMessage(e)),
        c("cdid_static_formula_error", "cdid_static_design_error")
      )
    }
  )

  if (nrow(mf) != n_req) {
    static_panel_abort(
      sprintf("Cell `%s` `%s` `model.frame` row count (%d) does not match required population count (%d).",
              cell$cell_id, role, nrow(mf), n_req),
      "cdid_static_design_error"
    )
  }

  # Verify no evaluated term in `mf` produced NA/NaN/Inf on required units.
  if (n_req > 0L && ncol(mf) > 0L) {
    for (j in seq_len(ncol(mf))) {
      col_mf <- mf[[j]]
      bad_mf <- is.na(col_mf) | (is.numeric(col_mf) & !is.finite(col_mf))
      if (any(bad_mf)) {
        bad_u <- panel$unit_id[req_units_idx[which(bad_mf)]]
        term_nm <- names(mf)[j]
        static_panel_abort(
          sprintf(
            "Cell `%s` `%s` formula term `%s` evaluated to missing or non-finite at `t_prev = %s` for required unit(s): [%s].",
            cell$cell_id, role, term_nm, cell$t_prev_label,
            paste(utils::head(as.character(bad_u), 5L), collapse = ", ")
          ),
          c("cdid_static_missing_covariate_error", "cdid_static_design_error"),
          extra = list(
            cell_id = cell$cell_id,
            role = role,
            population = pop_name,
            t_prev = cell$t_prev,
            term = term_nm,
            missing_units = bad_u
          )
        )
      }
    }
  }

  mm <- tryCatch(
    stats::model.matrix(f_info$terms, data = mf),
    error = function(e) {
      static_panel_abort(
        sprintf("Cell `%s` `%s` `model.matrix` failed at `t_prev = %s`: %s",
                cell$cell_id, role, cell$t_prev_label, conditionMessage(e)),
        c("cdid_static_design_error", "cdid_static_formula_error")
      )
    }
  )

  if (nrow(mm) != n_req || ncol(mm) < 1L || !identical(colnames(mm)[1L], "(Intercept)")) {
    static_panel_abort(
      sprintf("Cell `%s` `%s` `model.matrix` must retain all %d required rows and a leading `(Intercept)` column.",
              cell$cell_id, role, n_req),
      c("cdid_static_no_intercept_error", "cdid_static_design_error")
    )
  }
  if (n_req > 0L && (anyNA(mm) || any(!is.finite(mm)))) {
    bad_rows <- which(rowSums(!is.finite(mm)) > 0L)
    bad_u <- panel$unit_id[req_units_idx[bad_rows]]
    static_panel_abort(
      sprintf("Cell `%s` `%s` `model.matrix` contains non-finite values at `t_prev = %s` for required unit(s): [%s].",
              cell$cell_id, role, cell$t_prev_label,
              paste(utils::head(as.character(bad_u), 5L), collapse = ", ")),
      c("cdid_static_missing_covariate_error", "cdid_static_design_error"),
      extra = list(cell_id = cell$cell_id, role = role, missing_units = bad_u)
    )
  }

  mm_assign <- attr(mm, "assign")
  mm_contrasts <- attr(mm, "contrasts")
  mm_colnames <- colnames(mm)
  z_colnames <- mm_colnames[-1L]
  k_cols <- ncol(mm)
  p_cols <- k_cols - 1L

  storage.mode(mm) <- "double"
  X_tilde_req <- mm
  rownames(X_tilde_req) <- panel$unit_key[req_units_idx]
  Z_req <- X_tilde_req[, -1L, drop = FALSE]
  attr(X_tilde_req, "assign") <- mm_assign
  attr(X_tilde_req, "contrasts") <- mm_contrasts

  # Full N-by-columns matrices aligned with original `panel$unit_id`.
  # Rows outside `req_mask` stay `NA_real_` without having been evaluated.
  X_tilde_full <- matrix(
    NA_real_,
    nrow = n,
    ncol = k_cols,
    dimnames = list(panel$unit_key, mm_colnames)
  )
  Z_full <- matrix(
    NA_real_,
    nrow = n,
    ncol = p_cols,
    dimnames = list(panel$unit_key, z_colnames)
  )
  intercept_vec <- rep(NA_real_, n)
  names(intercept_vec) <- panel$unit_key

  if (n_req > 0L) {
    X_tilde_full[req_units_idx, ] <- X_tilde_req
    Z_full[req_units_idx, ] <- Z_req
    intercept_vec[req_units_idx] <- 1.0
  }

  recorded_xlevels <- stats::.getXlevels(f_info$terms, mf)
  if (is.null(recorded_xlevels)) recorded_xlevels <- list()

  diagnostics <- list(
    cell_id = cell$cell_id,
    role = role,
    population = pop_name,
    formula = f_info$formula,
    formula_text = f_info$formula_text,
    had_response = f_info$had_response,
    has_intercept = TRUE,
    variables = f_info$variables,
    term_labels = f_info$term_labels,
    t = cell$t,
    t_prev = cell$t_prev,
    t_idx = cell$t_idx,
    t_prev_idx = ref_idx,
    n_units_total = n,
    n_units_required = n_req,
    n_units_irrelevant = as.integer(n - n_req),
    irrelevant_rows_evaluated = FALSE,
    na_action = "na.pass_checked_on_required",
    rank_preprocessed = FALSE,
    k = k_cols,
    p = p_cols,
    colnames_X_tilde = mm_colnames,
    colnames_Z = z_colnames,
    assign = mm_assign,
    contrasts = mm_contrasts,
    xlevels = recorded_xlevels
  )

  structure(
    list(
      cell_id = cell$cell_id,
      role = role,
      population = pop_name,
      formula = f_info$formula,
      t = cell$t,
      t_prev = cell$t_prev,
      t_idx = cell$t_idx,
      t_prev_idx = ref_idx,
      n = n,
      unit_id = panel$unit_id,
      unit_key = panel$unit_key,
      mask = req_mask,
      required_idx = req_units_idx,
      reference_rows = ref_rows_all,
      n_required = n_req,
      has_intercept = TRUE,
      intercept = intercept_vec,
      X_tilde = X_tilde_full,
      Z = Z_full,
      X_tilde_required = X_tilde_req,
      Z_required = Z_req,
      k = k_cols,
      p = p_cols,
      colnames_X_tilde = mm_colnames,
      colnames_Z = z_colnames,
      z_names = z_colnames,
      assign = mm_assign,
      assign_Z = mm_assign[-1L],
      contrasts = mm_contrasts,
      xlevels = recorded_xlevels,
      terms = f_info$terms,
      diagnostics = diagnostics
    ),
    class = "cdid_static_reference_design"
  )
}

# Build a callable per-cell reference design extractor that validates formulas
# and panel factor levels once and extracts designs on demand per cell without
# storing N-by-columns matrices for all cells up front.
static_reference_design_extractor <- function(panel,
                                              outcome_formula = NULL,
                                              treatment_formula = NULL,
                                              observation_formula = NULL,
                                              selection_formula = NULL,
                                              xformla = NULL,
                                              xlev = NULL) {
  panel <- static_panel_require_object(panel)

  if (!is.null(xformla)) {
    if (!is.null(outcome_formula) && !identical(outcome_formula, xformla)) {
      static_panel_abort(
        "Supply either `outcome_formula` or `xformla` (or identical formulas), not conflicting specifications.",
        c("cdid_static_formula_error", "cdid_static_design_error")
      )
    }
    outcome_formula <- xformla
  }
  if (!is.null(selection_formula) && !is.null(observation_formula) &&
      !identical(selection_formula, observation_formula)) {
    static_panel_abort(
      "`observation_formula` and `selection_formula` are aliases for the observation/selection covariate formula `W` at `t - 1` and must not conflict.",
      c("cdid_static_formula_error", "cdid_static_design_error")
    )
  }
  obs_formula <- if (!is.null(observation_formula)) observation_formula else selection_formula

  if (is.null(outcome_formula) && is.null(treatment_formula) && is.null(obs_formula)) {
    static_panel_abort(
      "Supply at least one explicit covariate formula (`outcome_formula`, `treatment_formula`, or `observation_formula` / `selection_formula`).",
      c("cdid_static_formula_error", "cdid_static_design_error")
    )
  }

  prep_outcome <- if (!is.null(outcome_formula)) {
    static_panel_prepare_formula(panel, outcome_formula, role = "outcome", xlev = xlev)
  } else {
    NULL
  }
  prep_treatment <- if (!is.null(treatment_formula)) {
    static_panel_prepare_formula(panel, treatment_formula, role = "treatment", xlev = xlev)
  } else {
    NULL
  }
  prep_obs <- if (!is.null(obs_formula)) {
    static_panel_prepare_formula(panel, obs_formula, role = "observation", xlev = xlev)
  } else {
    NULL
  }

  extract_fn <- function(cell,
                         role = c("all", "outcome", "treatment", "observation", "selection"),
                         population = "risk",
                         mask = NULL) {
    cell_obj <- static_panel_resolve_cell(panel, cell)
    role_req <- match.arg(role)

    if (identical(role_req, "outcome")) {
      if (is.null(prep_outcome)) {
        static_panel_abort("No `outcome_formula` was configured on this extractor.", "cdid_static_formula_error")
      }
      return(static_extract_reference_design(
        panel = panel, cell = cell_obj, formula = prep_outcome$formula,
        role = "outcome", population = population, mask = mask,
        prepared_formula = prep_outcome
      ))
    }
    if (identical(role_req, "treatment")) {
      prep_use <- if (!is.null(prep_treatment)) prep_treatment else prep_outcome
      if (is.null(prep_use)) {
        static_panel_abort("No `treatment_formula` (or `outcome_formula`) was configured on this extractor.", "cdid_static_formula_error")
      }
      return(static_extract_reference_design(
        panel = panel, cell = cell_obj, formula = prep_use$formula,
        role = "treatment", population = population, mask = mask,
        prepared_formula = prep_use
      ))
    }
    if (role_req %in% c("observation", "selection")) {
      if (is.null(prep_obs)) {
        static_panel_abort("No `observation_formula` / `selection_formula` was configured on this extractor.", "cdid_static_formula_error")
      }
      return(static_extract_reference_design(
        panel = panel, cell = cell_obj, formula = prep_obs$formula,
        role = role_req, population = population, mask = mask,
        prepared_formula = prep_obs
      ))
    }

    # role_req == "all": extract all configured designs for `cell_obj`.
    out_des <- if (!is.null(prep_outcome)) {
      static_extract_reference_design(
        panel = panel, cell = cell_obj, formula = prep_outcome$formula,
        role = "outcome", population = population, mask = mask,
        prepared_formula = prep_outcome
      )
    } else {
      NULL
    }
    trt_des <- if (!is.null(prep_treatment)) {
      if (!is.null(prep_outcome) && identical(prep_treatment$formula, prep_outcome$formula)) {
        des_copy <- out_des
        des_copy$role <- "treatment"
        des_copy$diagnostics$role <- "treatment"
        des_copy
      } else {
        static_extract_reference_design(
          panel = panel, cell = cell_obj, formula = prep_treatment$formula,
          role = "treatment", population = population, mask = mask,
          prepared_formula = prep_treatment
        )
      }
    } else if (!is.null(out_des)) {
      des_copy <- out_des
      des_copy$role <- "treatment"
      des_copy$diagnostics$role <- "treatment"
      des_copy
    } else {
      NULL
    }
    obs_des <- if (!is.null(prep_obs)) {
      static_extract_reference_design(
        panel = panel, cell = cell_obj, formula = prep_obs$formula,
        role = "observation", population = population, mask = mask,
        prepared_formula = prep_obs
      )
    } else {
      NULL
    }

    structure(
      list(
        cell_id = cell_obj$cell_id,
        t = cell_obj$t,
        t_prev = cell_obj$t_prev,
        t_idx = cell_obj$t_idx,
        t_prev_idx = cell_obj$t_prev_idx,
        d = cell_obj$d,
        sigma = cell_obj$sigma,
        n = panel$n,
        unit_id = panel$unit_id,
        population = if (!is.null(out_des)) out_des$population else if (!is.null(trt_des)) trt_des$population else obs_des$population,
        mask = if (!is.null(out_des)) out_des$mask else if (!is.null(trt_des)) trt_des$mask else obs_des$mask,
        outcome = out_des,
        treatment = trt_des,
        observation = obs_des,
        selection = obs_des,
        separate_treatment_formula = !is.null(prep_treatment) && (is.null(prep_outcome) || !identical(prep_treatment$formula, prep_outcome$formula)),
        separate_observation_formula = !is.null(prep_obs) && (is.null(prep_outcome) || !identical(prep_obs$formula, prep_outcome$formula))
      ),
      class = "cdid_static_cell_designs"
    )
  }

  attr(extract_fn, "formulas") <- list(
    outcome_formula = if (!is.null(prep_outcome)) prep_outcome$formula else NULL,
    treatment_formula = if (!is.null(prep_treatment)) prep_treatment$formula else if (!is.null(prep_outcome)) prep_outcome$formula else NULL,
    observation_formula = if (!is.null(prep_obs)) prep_obs$formula else NULL
  )
  attr(extract_fn, "prepared") <- list(
    outcome = prep_outcome,
    treatment = prep_treatment,
    observation = prep_obs
  )
  class(extract_fn) <- c("cdid_static_design_extractor", "function")
  extract_fn
}

# Reference design extraction entry point:
# - If `formula` is supplied alone (with no `outcome_formula` / `treatment_formula` /
#   `observation_formula` / `selection_formula`), extracts a single
#   `cdid_static_reference_design` for `cell` (or returns a single-formula
#   callable if `cell = NULL`).
# - If `outcome_formula`, `treatment_formula`, `observation_formula`, or
#   `selection_formula` is supplied:
#   * when `cell` is supplied, extracts the requested design(s) for `cell`;
#   * when `cell` is `NULL`, returns the lazy callable extractor per cell.
static_reference_design <- function(panel,
                                    cell = NULL,
                                    formula = NULL,
                                    outcome_formula = NULL,
                                    treatment_formula = NULL,
                                    observation_formula = NULL,
                                    selection_formula = NULL,
                                    xformla = NULL,
                                    role = "outcome",
                                    population = "risk",
                                    mask = NULL,
                                    xlev = NULL) {
  panel <- static_panel_require_object(panel)

  has_single <- !is.null(formula)
  has_multi <- !is.null(outcome_formula) || !is.null(treatment_formula) ||
    !is.null(observation_formula) || !is.null(selection_formula) || !is.null(xformla)

  if (has_single && has_multi) {
    static_panel_abort(
      "Supply either `formula` (for single-formula extraction) or explicit `outcome_formula` / `treatment_formula` / `observation_formula` / `selection_formula`, not both.",
      c("cdid_static_formula_error", "cdid_static_design_error")
    )
  }

  if (has_single) {
    prep <- static_panel_prepare_formula(panel, formula = formula, role = role, xlev = xlev)
    if (is.null(cell)) {
      fn <- function(cell, population = "risk", mask = NULL) {
        static_extract_reference_design(
          panel = panel, cell = cell, formula = formula,
          role = role, population = population, mask = mask,
          prepared_formula = prep
        )
      }
      class(fn) <- c("cdid_static_design_extractor", "function")
      return(fn)
    }
    return(static_extract_reference_design(
      panel = panel, cell = cell, formula = formula,
      role = role, population = population, mask = mask,
      prepared_formula = prep
    ))
  }

  extractor <- static_reference_design_extractor(
    panel = panel,
    outcome_formula = outcome_formula,
    treatment_formula = treatment_formula,
    observation_formula = observation_formula,
    selection_formula = selection_formula,
    xformla = xformla,
    xlev = xlev
  )
  if (is.null(cell)) {
    return(extractor)
  }
  extractor(cell = cell, role = "all", population = population, mask = mask)
}

static_panel_reference_design <- static_reference_design

# Internal helper: validate a fixed role mask of length `panel$n`.
static_panel_validate_role_mask <- function(panel, mask, arg_name = "mask") {
  if (missing(mask) || is.null(mask) || !is.atomic(mask) || !is.null(dim(mask)) ||
      is.object(mask) || !(is.numeric(mask) || is.logical(mask)) ||
      length(mask) != panel$n || anyNA(mask)) {
    static_panel_abort(
      sprintf("Role `%s` must be a non-missing binary vector of length `n = %d`.", arg_name, panel$n),
      c("cdid_static_pair_mask_error", "cdid_static_pair_error")
    )
  }
  if (!is.null(names(mask)) && !identical(as.character(names(mask)), panel$unit_key)) {
    static_panel_abort(
      sprintf("Names of `%s` do not match `panel$unit_id` in the exact same order.", arg_name),
      c("cdid_static_pair_mask_error", "cdid_static_pair_error")
    )
  }
  if (is.numeric(mask) && any(!is.finite(mask) | !(mask %in% c(0, 1)))) {
    static_panel_abort(
      sprintf("`%s` must contain only binary values in `{0, 1}`.", arg_name),
      c("cdid_static_pair_mask_error", "cdid_static_pair_error")
    )
  }
  as.integer(mask)
}

# Internal helper: normalize and validate caller-explicit date pairs `(s, u)`.
static_panel_normalize_pairs_spec <- function(panel, pairs = NULL, s = NULL, u = NULL, pair_id = NULL) {
  has_pairs <- !is.null(pairs)
  has_su <- !is.null(s) || !is.null(u)
  if (has_pairs && has_su) {
    static_panel_abort(
      "Supply either `pairs` or explicit `s` and `u` vectors, not both.",
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }
  if (!has_pairs && !has_su) {
    static_panel_abort(
      "Supply explicit date pairs via `pairs` or (`s`, `u`); observation pairs are never chosen automatically.",
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }

  s_vec <- NULL
  u_vec <- NULL
  pid_vec <- pair_id

  if (has_su) {
    if (is.null(s) || is.null(u) || !is.atomic(s) || !is.atomic(u) ||
        !is.null(dim(s)) || !is.null(dim(u)) || is.object(s) || is.object(u) ||
        !is.numeric(s) || !is.numeric(u) || length(s) < 1L || length(s) != length(u)) {
      static_panel_abort(
        "`s` and `u` must be non-empty numeric vectors of equal length.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    s_vec <- s
    u_vec <- u
  } else if (is.matrix(pairs)) {
    if (!is.numeric(pairs) || is.object(pairs) || ncol(pairs) != 2L || nrow(pairs) < 1L) {
      static_panel_abort(
        "`pairs` matrix must be a numeric matrix with 2 columns `(s, u)` and at least one row.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    cn <- colnames(pairs)
    if (!is.null(cn) && setequal(cn, c("s", "u"))) {
      s_vec <- pairs[, "s"]
      u_vec <- pairs[, "u"]
    } else {
      s_vec <- pairs[, 1L]
      u_vec <- pairs[, 2L]
    }
    if (is.null(pid_vec) && !is.null(rownames(pairs)) && all(nzchar(rownames(pairs)))) {
      pid_vec <- rownames(pairs)
    }
  } else if (is.data.frame(pairs)) {
    df_p <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (nrow(df_p) < 1L) {
      static_panel_abort(
        "`pairs` data.frame must have at least one row.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    if (all(c("s", "u") %in% names(df_p))) {
      s_vec <- df_p[["s"]]
      u_vec <- df_p[["u"]]
    } else if (ncol(df_p) >= 2L) {
      s_vec <- df_p[[1L]]
      u_vec <- df_p[[2L]]
    } else {
      static_panel_abort(
        "`pairs` data.frame must have columns `s` and `u` (or at least 2 numeric columns).",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    if (is.null(pid_vec) && "pair_id" %in% names(df_p)) {
      pid_vec <- as.character(df_p[["pair_id"]])
    }
  } else if (is.list(pairs) && !is.object(pairs) && length(pairs) >= 1L) {
    n_p <- length(pairs)
    s_vec <- numeric(n_p)
    u_vec <- numeric(n_p)
    list_names <- names(pairs)
    for (k in seq_len(n_p)) {
      el <- pairs[[k]]
      if (is.list(el) && !is.null(el$s) && !is.null(el$u)) {
        s_vec[k] <- static_panel_check_scalar_numeric(el$s, sprintf("pairs[[%d]]$s", k), "cdid_static_pair_input_error")
        u_vec[k] <- static_panel_check_scalar_numeric(el$u, sprintf("pairs[[%d]]$u", k), "cdid_static_pair_input_error")
      } else if (is.atomic(el) && is.null(dim(el)) && !is.object(el) && is.numeric(el) && length(el) == 2L) {
        s_vec[k] <- static_panel_check_scalar_numeric(el[1L], sprintf("pairs[[%d]][1]", k), "cdid_static_pair_input_error")
        u_vec[k] <- static_panel_check_scalar_numeric(el[2L], sprintf("pairs[[%d]][2]", k), "cdid_static_pair_input_error")
      } else {
        static_panel_abort(
          sprintf("`pairs[[%d]]` must be a length-2 numeric vector `c(s, u)` or a list with `s` and `u`.", k),
          c("cdid_static_pair_input_error", "cdid_static_pair_error")
        )
      }
    }
    if (is.null(pid_vec) && !is.null(list_names) && all(!is.na(list_names)) && all(nzchar(list_names))) {
      pid_vec <- list_names
    }
  } else {
    static_panel_abort(
      "`pairs` must be a non-empty 2-column matrix/data.frame or list of `(s, u)` date pairs.",
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }

  if (!is.numeric(s_vec) || !is.numeric(u_vec) || is.object(s_vec) || is.object(u_vec) ||
      anyNA(s_vec) || anyNA(u_vec) || any(!is.finite(s_vec)) || any(!is.finite(u_vec))) {
    static_panel_abort(
      "Pair endpoints `s` and `u` must be finite numeric values without `NA`, `NaN`, or `Inf`.",
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }

  n_pairs <- length(s_vec)
  s_num <- as.double(s_vec)
  u_num <- as.double(u_vec)

  s_idx <- match(s_num, panel$dates_num)
  u_idx <- match(u_num, panel$dates_num)
  if (anyNA(s_idx) || anyNA(u_idx)) {
    bad_k <- which(is.na(s_idx) | is.na(u_idx))[1L]
    bad_end <- if (is.na(s_idx[bad_k])) sprintf("s = %s", as.character(s_vec[bad_k])) else sprintf("u = %s", as.character(u_vec[bad_k]))
    static_panel_abort(
      sprintf("Requested pair %d (`(%s, %s)`) has endpoint `%s` outside the study date grid [%s].",
              bad_k, as.character(s_vec[bad_k]), as.character(u_vec[bad_k]),
              bad_end, paste(panel$date_labels, collapse = ", ")),
      c("cdid_static_pair_endpoint_error", "cdid_static_pair_error"),
      extra = list(
        pair_index = bad_k,
        s = s_vec[bad_k],
        u = u_vec[bad_k],
        study_dates = panel$dates
      )
    )
  }

  if (any(s_num >= u_num) || any(s_idx >= u_idx)) {
    bad_k <- which(s_num >= u_num | s_idx >= u_idx)[1L]
    static_panel_abort(
      sprintf("Requested pair %d (`(%s, %s)`) violates `s < u`.",
              bad_k, panel$date_labels[s_idx[bad_k]], panel$date_labels[u_idx[bad_k]]),
      c("cdid_static_pair_order_error", "cdid_static_pair_error"),
      extra = list(
        pair_index = bad_k,
        s = panel$dates[s_idx[bad_k]],
        u = panel$dates[u_idx[bad_k]],
        s_idx = s_idx[bad_k],
        u_idx = u_idx[bad_k]
      )
    )
  }

  pair_keys <- sprintf("%d:%d", s_idx, u_idx)
  if (anyDuplicated(pair_keys) > 0L) {
    dup_k <- which(duplicated(pair_keys))[1L]
    static_panel_abort(
      sprintf("Duplicate observation pair `(%s, %s)` requested at index %d.",
              panel$date_labels[s_idx[dup_k]], panel$date_labels[u_idx[dup_k]], dup_k),
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }

  if (is.null(pid_vec)) {
    pid_vec <- sprintf("(%s,%s)", panel$date_labels[s_idx], panel$date_labels[u_idx])
  } else {
    if (!is.character(pid_vec) || !is.null(dim(pid_vec)) || is.object(pid_vec) ||
        length(pid_vec) != n_pairs || anyNA(pid_vec) || any(!nzchar(pid_vec)) ||
        anyDuplicated(pid_vec) > 0L) {
      static_panel_abort(
        sprintf("`pair_id` must be a character vector of %d unique non-empty labels.", n_pairs),
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
  }

  # Per-period gap length is the difference of panel DATE INDICES (u_idx - s_idx),
  # never the calendar numeric difference (u - s).
  gap_vec <- as.integer(u_idx - s_idx)
  calendar_span <- as.double(panel$dates_num[u_idx] - panel$dates_num[s_idx])

  list(
    n_pairs = n_pairs,
    pair_id = pid_vec,
    s = panel$dates[s_idx],
    u = panel$dates[u_idx],
    s_num = panel$dates_num[s_idx],
    u_num = panel$dates_num[u_idx],
    s_label = panel$date_labels[s_idx],
    u_label = panel$date_labels[u_idx],
    s_idx = as.integer(s_idx),
    u_idx = as.integer(u_idx),
    gap = gap_vec,
    calendar_span = calendar_span
  )
}

# Materialize a single observation pair k = (s_idx, u_idx) for a fixed role mask.
static_panel_build_one_pair <- function(panel, mask_int, pair_id, s_idx, u_idx) {
  n <- panel$n
  Ys <- panel$Y[, s_idx]
  Yu <- panel$Y[, u_idx]
  O_s <- panel$O[, s_idx]
  O_u <- panel$O[, u_idx]
  O_pair <- as.integer(O_s == 1L & O_u == 1L)
  A_O <- as.integer(mask_int == 1L & O_pair == 1L)

  # Actual outcome change dy = Yu - Ys on O == 1; NA_real_ on O == 0 (never 0).
  dy <- rep(NA_real_, n)
  obs_idx <- which(O_pair == 1L)
  if (length(obs_idx) > 0L) {
    dy[obs_idx] <- Yu[obs_idx] - Ys[obs_idx]
  }

  n_role <- as.integer(sum(mask_int))
  n_obs_all <- as.integer(length(obs_idx))
  n_obs_role <- as.integer(sum(A_O))
  n_obs_s_all <- as.integer(sum(O_s))
  n_obs_u_all <- as.integer(sum(O_u))
  n_obs_s_role <- as.integer(sum(mask_int == 1L & O_s == 1L))
  n_obs_u_role <- as.integer(sum(mask_int == 1L & O_u == 1L))
  gap <- as.integer(u_idx - s_idx)
  cal_span <- as.double(panel$dates_num[u_idx] - panel$dates_num[s_idx])

  counts <- list(
    n_units = n,
    n_role = n_role,
    n_observed_all = n_obs_all,
    n_observed_role = n_obs_role,
    n_missing_role = as.integer(n_role - n_obs_role),
    n_observed_s_all = n_obs_s_all,
    n_observed_u_all = n_obs_u_all,
    n_observed_s_role = n_obs_s_role,
    n_observed_u_role = n_obs_u_role,
    role_observation_rate = if (n_role > 0L) n_obs_role / as.double(n_role) else NA_real_,
    supported_in_role = (n_role > 0L && n_obs_role > 0L)
  )

  structure(
    list(
      pair_id = pair_id,
      s = panel$dates[s_idx],
      u = panel$dates[u_idx],
      s_num = panel$dates_num[s_idx],
      u_num = panel$dates_num[u_idx],
      s_label = panel$date_labels[s_idx],
      u_label = panel$date_labels[u_idx],
      s_idx = as.integer(s_idx),
      u_idx = as.integer(u_idx),
      gap = gap,
      gap_length = gap,
      index_gap = gap,
      L_k = gap,
      calendar_span = cal_span,
      n = n,
      unit_id = panel$unit_id,
      unit_key = panel$unit_key,
      mask = mask_int,
      O_s = as.integer(O_s),
      O_u = as.integer(O_u),
      O = O_pair,
      A_O = A_O,
      dy = dy,
      n_units = n,
      n_role = n_role,
      n_observed_all = n_obs_all,
      n_observed_role = n_obs_role,
      n_observed_s_role = n_obs_s_role,
      n_observed_u_role = n_obs_u_role,
      supported_in_role = counts$supported_in_role,
      counts = counts
    ),
    class = "cdid_static_observation_pair"
  )
}

# 4a. Single observation-pair extraction for a fixed role mask and explicit (s, u).
static_observation_pair <- function(panel, mask, s, u, pair_id = NULL) {
  panel <- static_panel_require_object(panel)
  mask_int <- static_panel_validate_role_mask(panel, mask, arg_name = "mask")
  spec <- static_panel_normalize_pairs_spec(panel, s = s, u = u, pair_id = pair_id)
  if (spec$n_pairs != 1L) {
    static_panel_abort(
      "`static_observation_pair()` expects a single `(s, u)` pair; use `static_observation_pairs()` for multiple pairs.",
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }
  static_panel_build_one_pair(
    panel = panel,
    mask_int = mask_int,
    pair_id = spec$pair_id[1L],
    s_idx = spec$s_idx[1L],
    u_idx = spec$u_idx[1L]
  )
}

# 4b. Lazy / chunked / callback accessor for observation pairs on a fixed role mask.
# Avoids materializing an N-by-all-pairs array up front.
static_observation_pair_accessor <- function(panel,
                                             mask,
                                             pairs = NULL,
                                             s = NULL,
                                             u = NULL,
                                             pair_id = NULL) {
  panel <- static_panel_require_object(panel)
  mask_int <- static_panel_validate_role_mask(panel, mask, arg_name = "mask")
  spec <- static_panel_normalize_pairs_spec(panel, pairs = pairs, s = s, u = u, pair_id = pair_id)

  n_pairs <- spec$n_pairs
  n_units <- panel$n
  n_role <- as.integer(sum(mask_int))

  # Compute lightweight pair count summary from the N x T binary matrix `panel$O`
  # without allocating N-vectors `dy` or `O` for all pairs simultaneously.
  n_obs_all <- integer(n_pairs)
  n_obs_role <- integer(n_pairs)
  n_obs_s_role <- integer(n_pairs)
  n_obs_u_role <- integer(n_pairs)
  role_idx <- which(mask_int == 1L)

  for (k in seq_len(n_pairs)) {
    si <- spec$s_idx[k]
    ui <- spec$u_idx[k]
    Os <- panel$O[, si]
    Ou <- panel$O[, ui]
    n_obs_all[k] <- as.integer(sum(Os == 1L & Ou == 1L))
    if (n_role > 0L) {
      Os_r <- Os[role_idx]
      Ou_r <- Ou[role_idx]
      n_obs_s_role[k] <- as.integer(sum(Os_r == 1L))
      n_obs_u_role[k] <- as.integer(sum(Ou_r == 1L))
      n_obs_role[k] <- as.integer(sum(Os_r == 1L & Ou_r == 1L))
    }
  }

  summary_table <- data.frame(
    pair_id = spec$pair_id,
    s = spec$s_num,
    u = spec$u_num,
    s_idx = spec$s_idx,
    u_idx = spec$u_idx,
    gap = spec$gap,
    calendar_span = spec$calendar_span,
    n_units = rep(n_units, n_pairs),
    n_role = rep(n_role, n_pairs),
    n_observed_all = n_obs_all,
    n_observed_role = n_obs_role,
    n_missing_role = as.integer(n_role - n_obs_role),
    n_observed_s_role = n_obs_s_role,
    n_observed_u_role = n_obs_u_role,
    supported_in_role = (n_role > 0L & n_obs_role > 0L),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  resolve_pair_index <- function(k) {
    if (is.character(k) && length(k) == 1L && !is.na(k)) {
      idx <- match(k, spec$pair_id)
      if (is.na(idx)) {
        static_panel_abort(
          sprintf("Unknown `pair_id` `%s`.", k),
          c("cdid_static_pair_input_error", "cdid_static_pair_error")
        )
      }
      return(idx)
    }
    if (is.numeric(k) && length(k) == 1L && is.finite(k) && k == floor(k) && k >= 1 && k <= n_pairs) {
      return(as.integer(k))
    }
    static_panel_abort(
      sprintf("Pair index must be an integer in `1..%d` or a declared `pair_id`.", n_pairs),
      c("cdid_static_pair_input_error", "cdid_static_pair_error")
    )
  }

  get_pair <- function(k) {
    idx <- resolve_pair_index(k)
    static_panel_build_one_pair(
      panel = panel,
      mask_int = mask_int,
      pair_id = spec$pair_id[idx],
      s_idx = spec$s_idx[idx],
      u_idx = spec$u_idx[idx]
    )
  }

  get_by_dates <- function(s_val, u_val) {
    s_num <- static_panel_check_scalar_numeric(s_val, "s", "cdid_static_pair_input_error")
    u_num <- static_panel_check_scalar_numeric(u_val, "u", "cdid_static_pair_input_error")
    idx <- which(spec$s_num == s_num & spec$u_num == u_num)
    if (length(idx) != 1L) {
      static_panel_abort(
        sprintf("Date pair `(%s, %s)` is not in the declared pair catalog.", as.character(s_val), as.character(u_val)),
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    get_pair(idx)
  }

  get_chunk <- function(indices) {
    if (missing(indices) || is.null(indices) || length(indices) < 1L) {
      static_panel_abort(
        "`indices` must be a non-empty vector of pair indices or `pair_id`s.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    idx_vec <- vapply(indices, resolve_pair_index, integer(1L))
    out <- lapply(idx_vec, get_pair)
    names(out) <- spec$pair_id[idx_vec]
    out
  }

  iterate <- function(callback) {
    if (missing(callback) || !is.function(callback)) {
      static_panel_abort(
        "`callback` must be a function taking `(pair, index)` or `(pair)`.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    n_formals <- length(formals(callback))
    res <- vector("list", n_pairs)
    names(res) <- spec$pair_id
    for (k in seq_len(n_pairs)) {
      p_k <- get_pair(k)
      res[[k]] <- if (n_formals >= 2L || "..." %in% names(formals(callback))) {
        callback(p_k, k)
      } else {
        callback(p_k)
      }
    }
    res
  }

  iterate_chunks <- function(chunk_size, callback) {
    if (missing(chunk_size) || !is.numeric(chunk_size) || length(chunk_size) != 1L ||
        is.na(chunk_size) || !is.finite(chunk_size) || chunk_size != floor(chunk_size) || chunk_size < 1L) {
      static_panel_abort(
        "`chunk_size` must be a positive integer.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    if (missing(callback) || !is.function(callback)) {
      static_panel_abort(
        "`callback` must be a function taking `(chunk_pairs, chunk_indices)` or `(chunk_pairs)`.",
        c("cdid_static_pair_input_error", "cdid_static_pair_error")
      )
    }
    cs <- as.integer(chunk_size)
    starts <- seq.int(1L, n_pairs, by = cs)
    n_chunks <- length(starts)
    n_formals <- length(formals(callback))
    res <- vector("list", n_chunks)
    for (c_idx in seq_len(n_chunks)) {
      idx_seq <- seq.int(starts[c_idx], min(n_pairs, starts[c_idx] + cs - 1L))
      chunk_k <- get_chunk(idx_seq)
      res[[c_idx]] <- if (n_formals >= 2L || "..." %in% names(formals(callback))) {
        callback(chunk_k, idx_seq)
      } else {
        callback(chunk_k)
      }
    }
    res
  }

  structure(
    list(
      n_units = n_units,
      n_role = n_role,
      n_pairs = n_pairs,
      unit_id = panel$unit_id,
      unit_key = panel$unit_key,
      mask = mask_int,
      pair_ids = spec$pair_id,
      summary_table = summary_table,
      supported_pair_ids = spec$pair_id[summary_table$supported_in_role],
      empty_pair_ids = spec$pair_id[!summary_table$supported_in_role],
      get_pair = get_pair,
      get_by_dates = get_by_dates,
      get_chunk = get_chunk,
      iterate = iterate,
      iterate_chunks = iterate_chunks
    ),
    class = "cdid_static_pair_accessor"
  )
}

# 4c. Observation-pair access for a fixed role mask and caller-explicit date pairs.
# - Returns every requested pair (including zero-count pairs) without graph pruning.
# - When `callback` is supplied, streams pairs (or chunks of size `chunk_size`)
#   to `callback` without materializing all `N x n_pairs` vectors simultaneously.
# - When `materialize = FALSE`, returns the lazy `cdid_static_pair_accessor`.
static_observation_pairs <- function(panel,
                                     mask,
                                     pairs = NULL,
                                     s = NULL,
                                     u = NULL,
                                     pair_id = NULL,
                                     callback = NULL,
                                     chunk_size = NULL,
                                     materialize = TRUE) {
  acc <- static_observation_pair_accessor(
    panel = panel,
    mask = mask,
    pairs = pairs,
    s = s,
    u = u,
    pair_id = pair_id
  )

  if (!is.null(callback)) {
    cb_results <- if (is.null(chunk_size)) {
      acc$iterate(callback)
    } else {
      acc$iterate_chunks(chunk_size = chunk_size, callback = callback)
    }
    return(structure(
      list(
        n_units = acc$n_units,
        n_role = acc$n_role,
        n_pairs = acc$n_pairs,
        unit_id = acc$unit_id,
        mask = acc$mask,
        pair_ids = acc$pair_ids,
        summary_table = acc$summary_table,
        supported_pair_ids = acc$supported_pair_ids,
        empty_pair_ids = acc$empty_pair_ids,
        chunk_size = if (!is.null(chunk_size)) as.integer(chunk_size) else NULL,
        callback_results = cb_results,
        accessor = acc
      ),
      class = "cdid_static_observation_pairs_stream"
    ))
  }

  if (!isTRUE(materialize)) {
    return(acc)
  }

  pair_list <- acc$get_chunk(seq_len(acc$n_pairs))

  structure(
    list(
      pairs = pair_list,
      pair_ids = acc$pair_ids,
      n_units = acc$n_units,
      n_role = acc$n_role,
      n_pairs = acc$n_pairs,
      unit_id = acc$unit_id,
      unit_key = acc$unit_key,
      mask = acc$mask,
      summary_table = acc$summary_table,
      supported_pair_ids = acc$supported_pair_ids,
      empty_pair_ids = acc$empty_pair_ids,
      get_pair = acc$get_pair,
      get_by_dates = acc$get_by_dates,
      get_chunk = acc$get_chunk,
      iterate = acc$iterate,
      iterate_chunks = acc$iterate_chunks
    ),
    class = "cdid_static_observation_pairs"
  )
}

static_panel_observation_pairs <- static_observation_pairs
