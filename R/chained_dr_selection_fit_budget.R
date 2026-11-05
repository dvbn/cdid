# Source inventory of ONE conditional DR treatment fit. Previous models remain
# in the unchanged full-D retained reservation; no global D-by-D scratch exists.
dr_selection_treatment_model_scratch <- function(n, p, q, d, B, u, h, ld,
                                                 ref_x = 0, ref_s = 0, name_bytes = 0) {
  prod <- ipw_owner_safe_product
  add <- ipw_owner_safe_sum
  # Six existing observer vectors y/w/pp/dy/kk/obs_rows plus the two raw date
  # endpoint vectors prob_s/prob_t occupy these same eight scalar slots.
  input <- prod(n, add(c(p, q, 8)))
  # The last gather/predict/reference chunk remains bound in the R frame while
  # fitting/derivatives run. Date products retain both endpoint predictions.
  chunk <- prod(32 * B, add(c(p, q, d, h, ref_x, ref_s, 32)))
  # Same QR/IRLS precedent as preflight_ipw_fitted_nuisance_view: eight n*p
  # matrix copies,64 row vectors,eight p*p factors and coefficient work.
  solver <- add(c(prod(8 * n, p), prod(64, n), prod(8 * p, p), prod(16, p)))
  # d_obs/score/residual_dr/residual_ipw and weighted/subset copies; old/new
  # full-cohort ld_mat; C/U/V crossprod,replacement and finite-check copies.
  derivative <- add(c(prod(n, add(c(4 * p, 4 * q, 64))), prod(4, ld),
    prod(8 * p, p), prod(8 * (p + 8), u), prod(64, add(c(p, q, u, 1)))))
  names <- prod(64, name_bytes)
  list(scratch = add(c(8192, input, chunk, max(solver, derivative), names)),
    input = input, chunk = chunk, solver = solver, derivative = derivative, names = names,
    n = n, p = p, q = q, d = d, B = B, channels = u, pair_channels = h, cohort_derivative_entries = ld,
    reference_x_columns = ref_x, reference_s_columns = ref_s)
}

# Count actual model/cohort incidence without allocating expanded channel sets.
# The union used by the fitter is no wider than the sum of its declared sets.
dr_selection_treatment_fit_scratch <- function(st, specs, up_map, p, q, d, B,
                                               reference = NULL, name_bytes = 0) {
  ref_x <- if (is.null(reference)) 0 else length(reference$x_columns)
  ref_s <- if (is.null(reference)) 0 else
    length(reference$designs[[reference$role_slot["s"]]]$columns)
  out <- NULL
  for (s in specs) {
    u <- h <- ld <- 0
    n <- as.double(s$n_observers)
    for (k in seq_along(s$selection_ids)) {
      ck <- length(up_map$pair_channels[[s$selection_ids[k]]])
      nk <- min(n, st$cohort_end[s$cohorts[k]] - st$cohort_start[s$cohorts[k]] + 1)
      u <- ipw_owner_safe_sum(c(u, ck)); h <- max(h, ck)
      ld <- max(ld, ipw_owner_safe_product(nk, ck))
    }
    current <- dr_selection_treatment_model_scratch(n, p, q, d, B, u, h, ld,
      ref_x, ref_s, name_bytes)
    if (is.null(out) || current$scratch > out$scratch) out <- current
  }
  out
}
