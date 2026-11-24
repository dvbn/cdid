# Resource errors at the public boundary.
#
# A fit is admitted stage by stage: every stage compares the numeric arrays it
# plans to allocate with the workspace limit (`control$numeric_workspace_bytes`)
# before allocating them, and the first stage that does not fit stops the call.
# The stage that fails first is often not the one with the largest requirement,
# and several counts depend on the chunk size, which is derived from the
# workspace limit, so the count of the failing check is not what the whole call
# needs, and no value of the limit is stated as sufficient.
#
# The public entry points restate such a refusal once, promptly: nothing is
# rerun, no callback is evaluated again, and neither the caller's time limit
# nor the random number state is touched. The restatement reads the typed
# record that the failing check attached to its condition: its counts
# (`requested` and `bound`, or the GMM payload pair), its label and the control
# its bound comes from (`limit_kind`, see `resource_bound_kind()`), or the
# typed fields of a fixed size limit. It never parses the message, decides
# from a class name, or compares a bound with the limit. A condition without
# such a record passes unchanged. The original diagnostic is kept in the
# message.
# The class vector of the original condition is kept, and
# `cdid_resource_refusal` with one of `cdid_workspace_refusal`,
# `cdid_work_limit_refusal` or `cdid_size_limit_refusal` is added; the numbers
# are in the `resource` field. `sufficient_value` is `NA` and `verified` is
# `FALSE` on every refusal: no route yet checks the complete request at another
# limit.

resource_fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

resource_bytes_text <- function(b) {
  u <- c("bytes", "KiB", "MiB", "GiB", "TiB", "PiB"); k <- 1L; v <- b
  while (v >= 1024 && k < length(u)) { v <- v / 1024; k <- k + 1L }
  paste0(resource_fmt(b), " bytes", if (k > 1L) sprintf(" (%.2f %s)", v, u[k]) else "")
}

# The function that called a count helper (NULL at top level).
resource_caller <- function(n) if (n > 0L) tryCatch(sys.function(n), error = function(e) NULL) else NULL

resource_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

# Governing control of a check, from the expression of its bound at the
# signal site (the shared count helpers record `substitute()` of their bound
# argument when the check fails): "workspace" when the bound is read directly
# from controls derived from the workspace limit (`max_matrix_entries`,
# `max_retained_entries`, `max_output_entries`, `max_scratch_entries`,
# `max_workspace_cells`, or `min()` of them), "work" when it is read directly
# from the work limit (`max_work_ops`, `max_work`), NA otherwise (a formal
# argument, a remainder such as `cap - used`, or a bound derived from the
# dimensions of the panel). A local name is resolved through its assignments
# in the calling function, and `if (cond) NULL else X` has the kind of X.
# Nothing here compares a bound with the limit.
resource_workspace_fields <- c("max_matrix_entries", "max_retained_entries", "max_output_entries", "max_scratch_entries",
  "max_workspace_cells")
resource_work_fields <- c("max_work_ops", "max_work")
resource_bound_kind <- function(expr, caller = NULL) {
  is_assign <- function(y) is.call(y) && length(y) == 3L &&
    (identical(y[[1L]], as.name("<-")) || identical(y[[1L]], as.name("=")))
  # The object of `x$field` must itself be a control: a name such as `ctrl`,
  # `control`, `bounds`, `bnds`, or ending in `_ctrl`, `_control`, `_bounds`,
  # `_bnd`, or a `$` chain whose last component is such a name.
  control_object <- function(o) {
    nm <- if (is.name(o)) as.character(o) else if (is.call(o) && identical(o[[1L]], as.name("$")) && length(o) == 3L)
      as.character(o[[3L]]) else return(FALSE)
    nm %in% c("ctrl", "control", "bounds", "bnds") || grepl("(_ctrl|_control|_bounds|_bnd)$", nm)
  }
  field_kind <- function(x, alias = TRUE) {
    if (is.name(x)) return(if (alias) local_kind(x) else NA_character_)
    if (is.call(x) && identical(x[[1L]], as.name("$")) && length(x) == 3L) {
      if (!control_object(x[[2L]])) return(NA_character_)
      f <- as.character(x[[3L]])
      if (f %in% resource_workspace_fields) return("workspace")
      if (f %in% resource_work_fields) return("work")
      return(NA_character_)
    }
    if (is.call(x) && identical(x[[1L]], as.name("min")) && length(x) > 1L) {
      k <- unique(vapply(as.list(x)[-1L], field_kind, "", alias = alias))
      return(if (length(k) == 1L) k else NA_character_)
    }
    if (is.call(x) && identical(x[[1L]], as.name("(")) && length(x) == 2L) return(field_kind(x[[2L]], alias))
    # `if (cond) NULL else X` or `if (cond) X else NULL`: a NULL bound means the
    # check is not applied, so only X can fail; the kind is that of X.
    if (is.call(x) && identical(x[[1L]], as.name("if")) && length(x) == 4L) {
      if (is.null(x[[3L]])) return(field_kind(x[[4L]], alias))
      if (is.null(x[[4L]])) return(field_kind(x[[3L]], alias))
    }
    NA_character_
  }
  # A local name (for example `max_me <- ctrl$max_matrix_entries` at the top
  # of the calling function) is resolved through its plain assignments in the
  # caller's own body, not in nested function definitions. It is typed only if
  # every such assignment is typed alike, and not at all if the name is a
  # formal argument, is assigned with `<<-` or `assign()` anywhere (nested
  # functions included), or is the target of a subassignment (`x[i] <- v`,
  # `x$a <- v`, `f(x) <- v`).
  local_kind <- function(x) {
    if (!is.name(x) || !is.function(caller)) return(NA_character_)
    nm <- as.character(x); rhs <- list(); unsafe <- FALSE
    target <- function(l) { while (is.call(l) && length(l) >= 2L) l <- l[[2L]]; if (is.name(l)) as.character(l) else "" }
    walk <- function(y, nested) {
      if (!is.call(y)) return(invisible())
      if (identical(y[[1L]], as.name("<<-")) && length(y) == 3L && identical(target(y[[2L]]), nm)) unsafe <<- TRUE
      if (identical(y[[1L]], as.name("assign")) && length(y) >= 2L && identical(y[[2L]], nm)) unsafe <<- TRUE
      if (is_assign(y) && !is.name(y[[2L]]) && identical(target(y[[2L]]), nm)) unsafe <<- TRUE
      if (is_assign(y) && is.name(y[[2L]]) && identical(as.character(y[[2L]]), nm) && !nested) {
        v <- y[[3L]]
        while (is_assign(v)) v <- v[[3L]]
        rhs[[length(rhs) + 1L]] <<- v
      }
      inner <- nested || identical(y[[1L]], as.name("function"))
      for (el in as.list(y)[-1L]) if (!missing(el)) walk(el, inner)
    }
    walk(body(caller), FALSE)
    if (unsafe || nm %in% names(formals(caller)) || !length(rhs)) return(NA_character_)
    k <- unique(vapply(rhs, field_kind, "", alias = FALSE))
    if (length(k) == 1L) k else NA_character_
  }
  tryCatch(field_kind(expr), error = function(e) NA_character_)
}


# Typed record of the failing check: list(kind, requested, allowed, label), or
# NULL when the signal site recorded no governing control.
resource_typed_record <- function(e) {
  d <- e$diagnostics
  if (!is.list(d)) return(NULL)
  lab <- if (is.character(d$label) && length(d$label) == 1L) d$label else NA_character_
  if (resource_scalar(d$exact_count_limit) || resource_scalar(d$matrix_dim_limit))
    return(list(kind = "size", requested = NA_real_, allowed = NA_real_, label = lab))
  kind <- if (is.character(d$limit_kind) && length(d$limit_kind) == 1L) d$limit_kind else NA_character_
  if (!kind %in% c("workspace", "work")) return(NULL)
  for (p in list(c("requested", "bound"), c("total_payload_entries", "max_matrix_entries"),
                 c("requested_entries", "max_matrix_entries")))
    if (resource_scalar(d[[p[1]]]) && resource_scalar(d[[p[2]]]) && d[[p[1]]] > d[[p[2]]])
      return(list(kind = kind, requested = as.double(d[[p[1]]]), allowed = as.double(d[[p[2]]]), label = lab))
  NULL
}

# Limits the call resolved from its control: the workspace in bytes and
# entries, and the work limit (Inf when not set). Used only in the wording.
resource_call_limits <- function(control) {
  ctl <- if (is.list(control)) control else list()
  bytes <- if (resource_scalar(ctl$numeric_workspace_bytes)) as.double(ctl$numeric_workspace_bytes) else
    chained_public_profile_defaults()$numeric_workspace_bytes
  work <- if (resource_scalar(ctl$counted_work_limit)) as.double(ctl$counted_work_limit) else Inf
  list(bytes = bytes, entries = floor(bytes / 8), work = work)
}

resource_condition <- function(e, kind, message, resource) {
  base <- setdiff(class(e), c("cdid_error", "error", "condition"))
  sub <- switch(kind, workspace = "cdid_workspace_refusal", work = "cdid_work_limit_refusal", size = "cdid_size_limit_refusal")
  if (!is.list(e$diagnostics)) e$diagnostics <- list()
  e$diagnostics$internal_message <- conditionMessage(e)
  e$message <- message
  e$resource <- c(resource, list(sufficient_value = NA_real_, verified = FALSE))
  class(e) <- unique(c(base, sub, "cdid_resource_refusal", "cdid_error", "error", "condition"))
  e
}

# Restate a refusal of a call whose limits are `limits` (resource_call_limits()).
# `aggregation = TRUE`: the call is an aggregation of a fit, whose limits were
# fixed when the fit was computed. The original diagnostic is kept in the
# message, as its second sentence ("Check: ...").
resource_restate <- function(e, limits, aggregation = FALSE) {
  if (!inherits(e, "cdid_error") || inherits(e, "cdid_resource_refusal")) return(e)
  rec <- resource_typed_record(e)
  if (is.null(rec)) return(e)
  check <- paste0(" Check: ", conditionMessage(e))
  if (identical(rec$kind, "size")) {
    return(resource_condition(e, "size", paste0("The call exceeds a size limit that no setting changes.", check,
      " Reduce the number of units, dates or comparisons in the call."),
      list(kind = "size_limit", setting = NA_character_, stage = rec$label, check = rec$label,
        units = NA_character_, requested = NA_real_, allowed = NA_real_,
        reason = "a fixed size limit that no setting changes")))
  }
  # A work record when the work limit is off (`Inf`, the default) names no
  # setting the user could change: pass the condition unchanged.
  if (identical(rec$kind, "work") && !is.finite(limits$work)) return(e)
  refit <- if (isTRUE(aggregation)) " The limit comes from the fit, so the fit must be rerun with a larger limit and aggregated again." else ""
  if (identical(rec$kind, "work")) {
    msg <- paste0(if (isTRUE(aggregation)) "The aggregation exceeds the work limit of the fit." else "The call exceeds the work limit.",
      check, " That check counted ", resource_fmt(rec$requested), " estimated work units against ",
      if (identical(rec$allowed, limits$work)) "the work limit `control$counted_work_limit`" else
        paste0("the ", resource_fmt(rec$allowed), " that remain for it under the work limit `control$counted_work_limit`"), ".", refit,
      " The workspace limit does not change it; raise `counted_work_limit` or remove it (`Inf`, the default).")
    return(resource_condition(e, "work", msg, list(kind = "work_limit", setting = "counted_work_limit", stage = rec$label,
      check = rec$label, units = "estimated work units", requested = rec$requested, allowed = rec$allowed, limit = limits$work,
      limit_from = if (isTRUE(aggregation)) "fit" else "call", reason = "the estimated work of a check exceeds the work limit")))
  }
  entries_text <- function(x) paste0(resource_fmt(x), " numeric entries (", sub(" [(]", ", ", sub("[)]$", "", resource_bytes_text(8 * x))), ")")
  msg <- paste0(if (isTRUE(aggregation)) "The aggregation does not fit in the workspace limit of the fit." else
      "The call does not fit in the workspace limit.", check,
    " That check counted ", entries_text(rec$requested), " against ",
    if (rec$allowed == limits$entries) paste0("the workspace limit of ", entries_text(limits$entries)) else
      paste0("the ", entries_text(rec$allowed), " that remain for it under the workspace limit of ", entries_text(limits$entries)),
    ". This is what that check counted, not what the whole call needs; other stages can need more.",
    if (isTRUE(aggregation)) paste0(" The workspace limit comes from the fit, so the fit must be rerun with a larger limit, ",
      "`att_gt(..., control = list(numeric_workspace_bytes = ...))`, and aggregated again.") else
      " The limit is set with `control = list(numeric_workspace_bytes = ...)`.",
    " It counts the numeric arrays that cdid plans to allocate, not the memory of the R process. Increase the workspace limit.")
  resource_condition(e, "workspace", msg, list(kind = "workspace", setting = "numeric_workspace_bytes", stage = rec$label,
    check = rec$label, units = "numeric entries (8 bytes each)", requested = rec$requested, allowed = rec$allowed,
    requested_bytes = 8 * rec$requested, allowed_bytes = 8 * rec$allowed, limit = limits$bytes, limit_entries = limits$entries,
    limit_from = if (isTRUE(aggregation)) "fit" else "call",
    reason = "the numeric arrays counted by a check exceed the workspace limit"))
}

# Evaluate `expr` (a promise: the user's call itself); a refusal is restated
# once and signaled. Nothing is evaluated again.
resource_guard <- function(expr, control) {
  tryCatch(expr, error = function(e) stop(resource_restate(e, resource_call_limits(control))))
}

# Aggregation of a fitted result: the workspace limit is the one stored with
# the fit.
resource_guard_fitted <- function(expr, MP) {
  tryCatch(expr, error = function(e) {
    bytes <- tryCatch(as.double(MP$chained$aggregation_state$public$bytes), error = function(err) NA_real_)
    if (length(bytes) != 1L || !is.finite(bytes)) bytes <- chained_public_profile_defaults()$numeric_workspace_bytes
    work <- tryCatch(as.double(MP$chained$aggregation_state$public$resolved$profile$counted_work_limit), error = function(err) NA_real_)
    if (length(work) != 1L || is.na(work)) work <- Inf
    stop(resource_restate(e, list(bytes = bytes, entries = floor(bytes / 8), work = work), aggregation = TRUE))
  })
}

# An augmented fit with time-varying (dated) inputs evaluates its influence
# functions in passes of at most the automatic chunk size, each with a fixed
# cost; with fewer than 16 units per pass such fits can run for a long time
# especially when several reference dates are in use. Signaled only with an
# automatic chunk size; no workspace value is stated, because the chunk size
# and the other arrays of the changed plan decide the length of the passes.
resource_speed_warning_rows <- function(rows, n_units) {
  passes <- ceiling(n_units / rows)
  msg <- paste0("The influence functions of this fit are evaluated in passes of ", rows, " units (",
    resource_fmt(passes), " passes over the ", resource_fmt(n_units), " units for each stage that needs them); ",
    "on augmented fits with time-varying covariates each pass has a fixed cost, so the fit may take a long time. ",
    "A larger workspace limit, `control = list(numeric_workspace_bytes = ...)`, allows longer passes. ",
    "A different workspace limit changes the chunk size, so results can differ in the last digits and multiplier ",
    "bootstrap draws can differ.")
  warning(structure(class = c("cdid_workspace_speed_warning", "cdid_warning", "warning", "condition"),
    list(message = msg, call = NULL, rows = rows, units = n_units, suggested_value = NA_real_)))
}
