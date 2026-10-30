# Formals of did 2.5.1 conditional_did_pretest(), in order. Ordinary mode
# forwards exactly these (did's function has no `...`).
upstream_conditional_did_pretest_formals <- c(
  "yname", "tname", "idname", "gname", "xformla", "data", "panel",
  "allow_unbalanced_panel", "control_group", "weightsname", "alp", "bstrap",
  "cband", "biters", "clustervars", "est_method", "print_details", "pl", "cores"
)

# Existing chained arguments appended after did's formals (D66 choice 1).
# Ordinary mode neither evaluates nor forwards any of them.
cdid_conditional_pretest_chained_formals <- c("anticipation", "base_period", "fix_weights", "mode", "frame",
  "selection", "selection_model", "selection_pooling", "selection_xformla", "selection_universe",
  "selection_rho", "selection_justification", "selection_est_method", "imputation_pooling", "control")

conditional_did_pretest <- function(yname,
                                    tname,
                                    idname = NULL,
                                    gname,
                                    xformla = NULL,
                                    data,
                                    panel = TRUE,
                                    allow_unbalanced_panel = NULL,
                                    control_group = c("nevertreated", "notyettreated"),
                                    weightsname = NULL,
                                    alp = 0.05,
                                    bstrap = TRUE,
                                    cband = TRUE,
                                    biters = 1000,
                                    clustervars = NULL,
                                    est_method = "ipw",
                                    print_details = FALSE,
                                    pl = FALSE,
                                    cores = 1,
                                    anticipation = 0,
                                    base_period = "varying",
                                    fix_weights = NULL,
                                    mode = NULL,
                                    frame = NULL,
                                    selection = "none",
                                    selection_model = "pair_logit",
                                    selection_pooling = "separate",
                                    selection_xformla = NULL,
                                    selection_universe = "F",
                                    selection_rho = NULL,
                                    selection_justification = NULL,
                                    selection_est_method = "ipw",
                                    imputation_pooling = "separate",
                                    control = NULL) {
  env <- environment()
  mode <- check_mode(mode)
  if (is.null(mode)) {
    # Only `panel` is needed to resolve the default mode.
    mode <- if (check_panel(panel)) "chained" else "ordinary"
  }

  if (mode == "ordinary") {
    # The appended chained arguments stay lazy and are not forwarded. The NULL
    # unbalanced sentinel maps to did's FALSE without being evaluated here.
    return(forward_call(
      quote(did::conditional_did_pretest), upstream_conditional_did_pretest_formals, env,
      override = list(allow_unbalanced_panel = quote(
        if (is.null(allow_unbalanced_panel)) FALSE else allow_unbalanced_panel
      )),
      dots = FALSE
    ))
  }

  # Chained mode (O14, D66): the chained conditional pretrend CvM test on the
  # accepted union-catalog fit of the requested built-in route; est_method
  # defaults to "ipw" as in did, independently of att_gt()'s "dr" default.
  chained_conditional_did_pretest(
    yname = yname, tname = tname, idname = idname, gname = gname, xformla = xformla, data = data,
    panel = panel, allow_unbalanced_panel = allow_unbalanced_panel, control_group = control_group,
    weightsname = weightsname, alp = alp, bstrap = bstrap, cband = cband, biters = biters,
    clustervars = clustervars, est_method = est_method, print_details = print_details, pl = pl, cores = cores,
    anticipation = anticipation, base_period = base_period, fix_weights = fix_weights, mode = mode, frame = frame,
    selection = selection, selection_model = selection_model, selection_pooling = selection_pooling,
    selection_xformla = selection_xformla, selection_universe = selection_universe, selection_rho = selection_rho,
    selection_justification = selection_justification, selection_est_method = selection_est_method,
    imputation_pooling = imputation_pooling, control = control, call = sys.call())
}
