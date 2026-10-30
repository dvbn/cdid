# Formals of did 2.5.1 att_gt(), in order. Everything before `...` is passed
# through unchanged in ordinary mode.
upstream_att_gt_formals <- c(
  "yname", "tname", "idname", "gname", "xformla", "data", "panel",
  "allow_unbalanced_panel", "control_group", "anticipation", "weightsname",
  "fix_weights", "alp", "bstrap", "cband", "biters", "clustervars",
  "est_method", "base_period", "faster_mode", "print_details", "pl", "cores",
  "compute_inffunc"
)

# Chained-only extension formals appended after did's formals and `...`.
# Ordinary mode neither evaluates nor forwards any of them.
cdid_chained_only_formals <- c("mode", "frame", "selection", "selection_model", "selection_pooling",
                               "selection_xformla", "selection_universe", "selection_rho", "selection_justification",
                               "selection_est_method", "imputation_pooling",
                               "gmm_weighting", "gmm_covariance", "pt_assumption", "pt_start", "control")

att_gt <- function(yname,
                   tname,
                   idname = NULL,
                   gname,
                   xformla = NULL,
                   data,
                   panel = TRUE,
                   allow_unbalanced_panel = NULL,
                   control_group = c("nevertreated", "notyettreated"),
                   anticipation = 0,
                   weightsname = NULL,
                   fix_weights = NULL,
                   alp = 0.05,
                   bstrap = TRUE,
                   cband = TRUE,
                   biters = 1000,
                   clustervars = NULL,
                   est_method = "dr",
                   base_period = "varying",
                   faster_mode = TRUE,
                   print_details = FALSE,
                   pl = FALSE,
                   cores = 1,
                   compute_inffunc = TRUE,
                   ...,
                   mode = NULL,
                   frame = NULL,
                   selection = "none",
                   selection_model = "pair_logit",
                   selection_pooling = "separate",
                   selection_xformla = NULL,
                   selection_universe = "F", selection_rho = NULL, selection_justification = NULL,
                   selection_est_method = "ipw",
                   imputation_pooling = "separate",
                   gmm_weighting = c("optimal", "identity"),
                   gmm_covariance = "tikhonov",
                   pt_assumption = c("post", "all"),
                   pt_start = NULL,
                   control = NULL) {
  env <- environment()
  mode <- check_mode(mode)
  if (is.null(mode)) {
    # Only `panel` is needed to resolve the default mode.
    mode <- if (check_panel(panel)) "chained" else "ordinary"
  }

  if (mode == "ordinary") {
    # Chained-only arguments (including `control`) are neither forwarded nor
    # evaluated here. did has no NULL sentinel: an explicit NULL means did's
    # default, FALSE. The replacement stays lazy, so the argument is evaluated
    # by did when did first needs it.
    return(forward_call(
      quote(did::att_gt), upstream_att_gt_formals, env,
      override = list(allow_unbalanced_panel = quote(
        if (is.null(allow_unbalanced_panel)) FALSE else allow_unbalanced_panel
      ))
    ))
  }

  # OR preprocessing records the maximum absolute reconstruction discrepancy
  # and checks the maximum scaled discrepancy at 64 * eps * original columns.
  # Chained mode: est_method = "dr" (default), "ipw", "reg", or a custom
  # did-style panel IPW function. The scope gate refuses every unsupported
  # request (unsupported est_method strings, unsupported selection combinations,
  # time-varying weights, undeclared parallel bootstrap workers, repeated cross
  # sections, allow_unbalanced_panel = FALSE) before the data promise is forced;
  # nothing is substituted silently. compute_inffunc, pl/cores and print_details
  # keep their did meanings (point estimates only, multiplier workers, progress).
  # selection_est_method = "ipw" (default) keeps the
  # probability-weighted correction and never evaluates imputation_pooling;
  # "imputation" never evaluates selection_model/selection_pooling/
  # selection_xformla; "dr" (augmented correction, est_method = "dr" only)
  # evaluates both groups. All stay lazy promises here.
  chained_public_att_gt(
    yname = yname, tname = tname, idname = idname, gname = gname, xformla = xformla, data = data,
    panel = panel, allow_unbalanced_panel = allow_unbalanced_panel, control_group = control_group,
    anticipation = anticipation, weightsname = weightsname, fix_weights = fix_weights, alp = alp,
    bstrap = bstrap, cband = cband, biters = biters, clustervars = clustervars, est_method = est_method,
    base_period = base_period, faster_mode = faster_mode, print_details = print_details, pl = pl,
    cores = cores, compute_inffunc = compute_inffunc, dots = as.list(substitute(list(...)))[-1L],
    dots_eval = function() list(...), mode = mode, frame = frame,
    selection = selection, selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla,
    selection_universe = selection_universe, selection_rho = selection_rho, selection_justification = selection_justification,
    selection_est_method = selection_est_method, imputation_pooling = imputation_pooling, gmm_weighting = gmm_weighting,
    gmm_covariance = gmm_covariance, pt_assumption = pt_assumption, pt_start = pt_start, control = control, call = sys.call())
}

check_mode <- function(mode) {
  if (is.null(mode)) {
    return(NULL)
  }
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !(mode %in% c("chained", "ordinary"))) {
    cdid_abort(
      paste0("`mode` must be NULL, \"chained\" or \"ordinary\". NULL selects ",
             "\"chained\" for panel data and \"ordinary\" for repeated cross ",
             "sections."),
      "cdid_mode_error"
    )
  }
  mode
}

check_panel <- function(panel) {
  if (!is.logical(panel) || length(panel) != 1L || is.na(panel)) {
    cdid_abort(
      paste0("`panel` must be TRUE or FALSE; it determines the default mode ",
             "(chained for panels, ordinary for repeated cross sections)."),
      "cdid_mode_error"
    )
  }
  panel
}

cdid_version <- function() {
  # Source harnesses (no installed namespace) still produce the refusal text.
  tryCatch(as.character(getNamespaceVersion("cdid")), error = function(e) "source tree, no installed namespace")
}
