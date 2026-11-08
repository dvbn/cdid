# Private static missing-outcome observation kernel: role-pair outcome
# reconstruction and its complete first-order influence function, CONDITIONAL
# on caller-supplied fitted nuisance predictions, parameter Jacobians and
# original-unit coefficient-IF matrices.
#
# Mathematical components implemented here (Bellego, Benatia and Dortet-Bernadet
# 2025; de Chaisemartin, D'Haultfoeuille, Pasquier, Sow and Vazquez-Bare 2026,
# arXiv:2201.06898v7):
#   - Reconstructed normalized mean M_hat(g) and stayer cross-moment covariance C_hat
#   - Final pair-probability lower floor phi = 1 / max(q, eps) and its active-region derivative
#   - Observed-data population ratios and complete original-unit IF of the mean and covariance
#     (with parameter loading tables for observation and treatment nuisances)
#   - Augmented observation correction normalized by the role-group mass E_n[A g],
#     never by the observed inverse-weight mass
#
# One call = one fixed role/target cell and one observed date pair. Routes:
#   baseline    w = A O h,           Y = dy
#   weighting   w = A O h phi,       Y = dy
#   imputation  w = A h,             Y = j
#   augmented   w = A h,             Y = j + O phi (dy - j)
# with phi = 1 / ipw_clip_floor(q). M = sum(w Y) / sum(w).
# Covariance (tilt identically 1): Zbar = sum(w Z) / sum(w),
#   C = sum(w (Z - Zbar) Y) / sum(w), with the route-specific w above (h = 1).
#
# Influence functions are on the original-unit scale over ALL n supplied units
# (n = length(id), units outside the role included):
#   M_hat - M* ~ mean(IF),  IF_i = n w_i (Y_i - M) / S  +  sum_k L_k[i, ] J_k
# where S = sum(w), L_k is the caller's n x p_k coefficient IF for parameter
# channel k (theta_hat_k - theta_k ~ mean of rows of L_k) and J_k is the
# derivative of M_hat in theta_k obtained by contracting the caller's
# prediction-row x p_k Jacobian. No n x n object is ever formed.
#
# Not implemented here (owned elsewhere): nuisance fitting, pair/date-product
# independence assessment, graph selection and flows, treatment adjustment,
# full-role design means for final centering, public interfaces, survey
# weights. Clipping of the treatment tilt h is upstream: h arrives already
# clipped with its Jacobian already carrying the clip indicator.
#
# Depends on `cdid_abort` (R/forward.R) and on `ipw_clip_floor`,
# `ipw_clip_floor_active` and `ipw_clip_eps` (R/chained_ipw_link.R).

static_obs_routes <- c("baseline", "weighting", "imputation", "augmented")

static_obs_abort <- function(message, class) {
  cdid_abort(message, c(class, "cdid_static_obs_error"))
}

static_obs_route <- function(route) {
  if (!is.character(route) || length(route) != 1L || is.na(route) ||
      !(route %in% static_obs_routes)) {
    static_obs_abort(
      paste0("`route` must be exactly one of: ",
             paste(static_obs_routes, collapse = ", "), "."),
      "cdid_static_obs_route_error"
    )
  }
  route
}

static_obs_ids <- function(id) {
  if (missing(id) || is.null(id) || !is.atomic(id) || !is.null(dim(id)) ||
      is.object(id) || length(id) < 1L || anyNA(id) ||
      !(is.numeric(id) || is.character(id))) {
    static_obs_abort(
      "`id` must be a non-empty atomic numeric or character vector of original-unit ids without NA.",
      "cdid_static_obs_input_error"
    )
  }
  ids <- as.character(id)
  if (anyDuplicated(id) || anyDuplicated(ids)) {
    static_obs_abort("`id` must not contain duplicated original units.",
                     "cdid_static_obs_input_error")
  }
  ids
}

static_obs_check_labels <- function(labels, ids, what) {
  if (!is.null(labels) && !identical(as.character(labels), ids)) {
    static_obs_abort(
      paste0("Names of `", what, "` do not equal `id` in the same order; reordered or foreign units are refused."),
      "cdid_static_obs_order_error"
    )
  }
}

static_obs_num_vec <- function(x, n, ids, what) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      !(is.numeric(x) || (is.logical(x) && all(is.na(x)))) || length(x) != n) {
    static_obs_abort(
      paste0("`", what, "` must be a plain numeric vector of length length(id)."),
      "cdid_static_obs_input_error"
    )
  }
  static_obs_check_labels(names(x), ids, what)
  as.double(x)
}

# 0/1 mask; only `rows` are read (NULL = every row). Returns logical with FALSE
# on rows that are not read, so unread values (including NA) never propagate.
static_obs_mask <- function(x, n, ids, what, rows = NULL) {
  if (is.null(x) || !is.atomic(x) || !is.null(dim(x)) || is.object(x) ||
      !(is.numeric(x) || is.logical(x)) || length(x) != n) {
    static_obs_abort(
      paste0("`", what, "` must be a plain 0/1 (or logical) vector of length length(id)."),
      "cdid_static_obs_input_error"
    )
  }
  static_obs_check_labels(names(x), ids, what)
  if (is.null(rows)) rows <- seq_len(n)
  v <- as.double(x[rows])
  if (anyNA(v) || any(v != 0 & v != 1)) {
    static_obs_abort(
      paste0("`", what, "` must be 0/1 without NA on every row it is read."),
      "cdid_static_obs_input_error"
    )
  }
  out <- logical(n)
  out[rows] <- v == 1
  out
}

static_obs_require_finite <- function(v, what, rows_desc) {
  if (length(v) && any(!is.finite(v))) {
    static_obs_abort(
      paste0("`", what, "` must be finite on ", rows_desc, "."),
      "cdid_static_obs_input_error"
    )
  }
}

static_obs_matrix <- function(x, n, ids, what) {
  if (is.atomic(x) && is.null(dim(x)) && !is.object(x) && is.numeric(x) &&
      length(x) == n) {
    nm <- names(x)
    x <- matrix(as.double(x), ncol = 1L)
    static_obs_check_labels(nm, ids, what)
    return(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != n || ncol(x) < 1L) {
    static_obs_abort(
      paste0("`", what, "` must be a numeric matrix with length(id) rows and at least one column."),
      "cdid_static_obs_input_error"
    )
  }
  static_obs_check_labels(rownames(x), ids, what)
  storage.mode(x) <- "double"
  x
}

# Floored inverse pair probability (eq:floor) and its derivative in raw q
# (eq:phidot). The floor is applied once, to the supplied final joint q (for a
# date product the caller passes the product). Strict existing convention:
# active iff q > ipw_clip_eps; at or below the floor the derivative is zero.
static_obs_floor_inverse <- function(q) {
  act <- ipw_clip_floor_active(q)
  qf <- ipw_clip_floor(q)
  list(phi = 1 / qf, act = act, dphi_dq = -act / (qf * qf))
}

# Keyed parameter channels for one nuisance prediction. `NULL` is refused so
# that an estimated nuisance cannot silently lose its IF term; `list()` states
# explicitly that the prediction is known (no estimated parameter channel).
static_obs_channels <- function(jac, n, ids, what) {
  if (is.null(jac) || !is.list(jac) || is.object(jac)) {
    static_obs_abort(
      paste0("Supply `", what, "` explicitly as a keyed list of Jacobian matrices (use list() for a known, non-estimated prediction)."),
      "cdid_static_obs_channel_error"
    )
  }
  if (!length(jac)) return(list())
  keys <- names(jac)
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys))) {
    static_obs_abort(paste0("Every channel of `", what, "` needs an explicit non-empty key."),
                     "cdid_static_obs_channel_error")
  }
  if (anyDuplicated(keys)) {
    static_obs_abort(
      paste0("Repeated channel key in `", what, "`: ",
             paste(unique(keys[duplicated(keys)]), collapse = ", "),
             ". Sum contributions to one parameter before calling; equal arrays are never deduplicated."),
      "cdid_static_obs_channel_error"
    )
  }
  for (k in keys) {
    m <- jac[[k]]
    if (!is.matrix(m) || !is.numeric(m) || nrow(m) != n || ncol(m) < 1L) {
      static_obs_abort(
        paste0("Channel `", k, "` of `", what, "` must be a numeric matrix with length(id) prediction rows and p >= 1 columns."),
        "cdid_static_obs_channel_error"
      )
    }
    static_obs_check_labels(rownames(m), ids, paste0(what, "[[\"", k, "\"]]"))
  }
  jac
}

# Contract one channel: crossprod over the rows actually read. Returns p x m.
static_obs_contract_rows <- function(m, rows, coef, key, what) {
  if (!length(rows)) {
    return(matrix(0, ncol(m), ncol(coef)))
  }
  sub <- m[rows, , drop = FALSE]
  if (any(!is.finite(sub))) {
    static_obs_abort(
      paste0("Channel `", key, "` of `", what, "` has nonfinite Jacobian entries on rows it enters."),
      "cdid_static_obs_channel_error"
    )
  }
  crossprod(sub, coef)
}

static_obs_add_loading <- function(acc, key, value, cols) {
  if (is.null(acc$val[[key]])) {
    acc$val[[key]] <- value
    acc$cols[[key]] <- cols
    return(acc)
  }
  if (!identical(dim(acc$val[[key]]), dim(value)) ||
      (!is.null(cols) && !is.null(acc$cols[[key]]) && !identical(cols, acc$cols[[key]]))) {
    static_obs_abort(
      paste0("Channel key `", key, "` is shared by several nuisances with different parameter dimensions or names."),
      "cdid_static_obs_channel_error"
    )
  }
  acc$val[[key]] <- acc$val[[key]] + value
  if (is.null(acc$cols[[key]])) acc$cols[[key]] <- cols
  acc
}

# Contract keyed loadings with keyed original-unit coefficient IFs:
#   sum_k coef_if[[k]] %*% t(loading_k).
# `loadings[[k]]` is a length-p vector (mean) or m x p matrix (covariance).
# Every loading key must have exactly one coefficient IF, rownames must equal
# `id` in order, and no extra keys are accepted.
static_obs_contract <- function(loadings, coef_if, id) {
  ids <- static_obs_ids(id)
  n <- length(ids)
  if (is.null(coef_if) || !is.list(coef_if) || is.object(coef_if)) {
    static_obs_abort("`coef_if` must be a keyed list of n x p coefficient-IF matrices.",
                     "cdid_static_obs_channel_error")
  }
  lk <- names(loadings)
  if (is.null(lk)) lk <- character()
  ck <- names(coef_if)
  if (length(coef_if) && (is.null(ck) || anyNA(ck) || any(!nzchar(ck)))) {
    static_obs_abort("Every coefficient IF needs an explicit non-empty key.",
                     "cdid_static_obs_channel_error")
  }
  if (is.null(ck)) ck <- character()
  if (anyDuplicated(ck)) {
    static_obs_abort(
      paste0("Repeated coefficient-IF key: ", paste(unique(ck[duplicated(ck)]), collapse = ", "), "."),
      "cdid_static_obs_channel_error"
    )
  }
  miss <- setdiff(lk, ck)
  extra <- setdiff(ck, lk)
  if (length(miss) || length(extra)) {
    static_obs_abort(
      paste0("Coefficient-IF keys must equal the loading keys exactly; missing: [",
             paste(miss, collapse = ", "), "], extra: [", paste(extra, collapse = ", "), "]."),
      "cdid_static_obs_channel_error"
    )
  }
  out <- NULL
  for (k in lk) {
    L <- coef_if[[k]]
    J <- loadings[[k]]
    Jm <- if (is.matrix(J)) J else matrix(J, nrow = 1L, dimnames = list(NULL, names(J)))
    if (!is.matrix(L) || !is.numeric(L) || nrow(L) != n || ncol(L) != ncol(Jm)) {
      static_obs_abort(
        paste0("Coefficient IF `", k, "` must be an n x p numeric matrix with p equal to its Jacobian columns."),
        "cdid_static_obs_channel_error"
      )
    }
    if (is.null(rownames(L)) || !identical(rownames(L), ids)) {
      static_obs_abort(
        paste0("Coefficient IF `", k, "` must carry rownames equal to `id` in the same order; reordered or unlabeled unit rows are refused."),
        "cdid_static_obs_order_error"
      )
    }
    if (!is.null(colnames(L)) && !is.null(colnames(Jm)) &&
        !identical(colnames(L), colnames(Jm))) {
      static_obs_abort(
        paste0("Coefficient IF `", k, "` column names differ from its Jacobian parameter names."),
        "cdid_static_obs_channel_error"
      )
    }
    if (any(!is.finite(L))) {
      static_obs_abort(paste0("Coefficient IF `", k, "` must be finite."),
                       "cdid_static_obs_channel_error")
    }
    term <- L %*% t(Jm)
    out <- if (is.null(out)) term else out + term
  }
  out
}

static_obs_core <- function(kind, id, A, O, dy, h, Z, route, q, j,
                            q_jacobian, j_jacobian, h_jacobian, coef_if) {
  route <- static_obs_route(route)
  ids <- static_obs_ids(id)
  n <- length(ids)
  uses_O <- route != "imputation"
  uses_q <- route %in% c("weighting", "augmented")
  uses_j <- route %in% c("imputation", "augmented")

  # Inputs a route does not read are refused, never silently rerouted.
  unused <- character()
  if (!uses_O && !is.null(O)) unused <- c(unused, "O")
  if (!uses_O && !is.null(dy)) unused <- c(unused, "dy")
  if (!uses_q && !is.null(q)) unused <- c(unused, "q")
  if (!uses_q && !is.null(q_jacobian)) unused <- c(unused, "q_jacobian")
  if (!uses_j && !is.null(j)) unused <- c(unused, "j")
  if (!uses_j && !is.null(j_jacobian)) unused <- c(unused, "j_jacobian")
  if (length(unused)) {
    static_obs_abort(
      paste0("Route `", route, "` does not read: ", paste(unused, collapse = ", "),
             ". Declare the intended route instead of supplying unused inputs."),
      "cdid_static_obs_route_error"
    )
  }

  A <- static_obs_mask(A, n, ids, "A")
  inA <- which(A)
  obs <- logical(n)
  if (uses_O) obs <- static_obs_mask(O, n, ids, "O", rows = inA)
  R_o <- which(A & obs)                       # role x pair-observed rows
  R_w <- if (route %in% c("baseline", "weighting")) R_o else inA
  # Observed support (as:reg-m(e)) is a baseline/weighting requirement only:
  # with A O = 0 throughout the role, the augmented pseudo-change is j and the
  # moment reduces to imputation over the full role (eq:M-hat); the empty-role
  # and positive-denominator checks below still apply.
  if (route %in% c("baseline", "weighting") && !length(R_o)) {
    static_obs_abort(
      paste0("Route `", route, "` has no observed pair in the role population (insufficient observed support)."),
      "cdid_static_obs_support_error"
    )
  }
  if (!length(R_w)) {
    static_obs_abort("The role population is empty.", "cdid_static_obs_support_error")
  }

  # Read exactly what the route uses.
  dy_o <- NULL
  if (uses_O) {
    dy_o <- static_obs_num_vec(dy, n, ids, "dy")[R_o]
    static_obs_require_finite(dy_o, "dy", "role units with an observed pair")
  }
  fl <- NULL
  q_o <- NULL
  if (uses_q) {
    q_o <- static_obs_num_vec(q, n, ids, "q")[R_o]
    static_obs_require_finite(q_o, "q", "role units with an observed pair")
    if (any(q_o < 0 | q_o > 1)) {
      static_obs_abort("Raw pair probability `q` must lie in [0, 1] where read.",
                       "cdid_static_obs_input_error")
    }
    fl <- static_obs_floor_inverse(q_o)
  }
  j_w <- NULL
  if (uses_j) {
    j_w <- static_obs_num_vec(j, n, ids, "j")[R_w]
    static_obs_require_finite(j_w, "j", "role units")
  }
  if (kind == "mean") {
    h_w <- static_obs_num_vec(h, n, ids, "h")[R_w]
    static_obs_require_finite(h_w, "h", "rows entering the weight")
    if (any(h_w < 0)) {
      static_obs_abort("Tilt `h` must be nonnegative where read.", "cdid_static_obs_input_error")
    }
  } else {
    h_w <- rep(1, length(R_w))                # covariance: tilt identically 1
  }

  # Position of observed rows inside R_w.
  o_in_w <- match(R_o, R_w)                   # NA-free when R_o subset of R_w
  nw <- length(R_w)
  phi_w <- numeric(nw)                        # O * phi on R_w (0 when O = 0)
  dphi_w <- numeric(nw)                       # O * dphi/dq on R_w
  act_w <- logical(nw)
  if (uses_q) {
    phi_w[o_in_w] <- fl$phi
    dphi_w[o_in_w] <- fl$dphi_dq
    act_w[o_in_w] <- fl$act > 0
  }
  resid_w <- numeric(nw)                      # O * (dy - j), augmented only
  if (route == "augmented") resid_w[o_in_w] <- dy_o - j_w[o_in_w]

  # Route-specific weight and pseudo-change on R_w (eq:M-hat table).
  if (route == "baseline") {
    w <- h_w
    Y <- dy_o
  } else if (route == "weighting") {
    w <- h_w * phi_w
    Y <- dy_o
  } else if (route == "imputation") {
    w <- h_w
    Y <- j_w
  } else {
    w <- h_w
    Y <- j_w + phi_w * resid_w
  }
  S <- sum(w)
  if (!is.finite(S) || S <= 0) {
    static_obs_abort(
      paste0("Zero or nonfinite denominator sum(A * omega", if (kind == "mean") " * h" else "",
             ") for route `", route, "`."),
      "cdid_static_obs_support_error"
    )
  }
  M <- sum(w * Y) / S

  # Direct term and the factors multiplying weight / pseudo-change derivatives.
  if (kind == "mean") {
    est <- M
    Fw <- matrix(Y - M, ncol = 1L)            # d w  enters with (Y - M)
    Gy <- matrix(1, nw, 1L)                   # d Y  enters with 1
    direct <- numeric(n)
    direct[R_w] <- n * w * (Y - M) / S
    names(direct) <- ids
    out_cols <- NULL
    Zbar <- NULL
  } else {
    Zm <- static_obs_matrix(Z, n, ids, "Z")
    Zw <- Zm[R_w, , drop = FALSE]
    static_obs_require_finite(Zw, "Z", "rows entering the covariance weight")
    Zbar <- colSums(w * Zw) / S               # route-specific weighted design mean
    Zc <- sweep(Zw, 2L, Zbar, "-")
    est <- drop(crossprod(Zc, w * Y)) / S     # eq:M-hat covariance
    Fw <- Zc * (Y - M) - matrix(est, nw, length(est), byrow = TRUE)
    Gy <- Zc
    direct <- matrix(0, n, ncol(Zm), dimnames = list(ids, colnames(Zm)))
    direct[R_w, ] <- n * w * Fw / S
    out_cols <- colnames(Zm)
    names(est) <- out_cols
    names(Zbar) <- out_cols
  }

  # Nuisance channels. Row coefficients (eq:chi-M / eq:chi-C loading table):
  #   h: d w / d h  = (phi for weighting, else 1)           x Fw
  #   q: weighting  d w / d q = h * dphi/dq                  x Fw
  #      augmented  w * d Y / d q = h * (dy - j) * dphi/dq   x Gy
  #   j: w * d Y / d j = h * (1 for imputation, 1 - O phi for augmented) x Gy
  acc <- list(val = list(), cols = list())
  add_channels <- function(acc, jac, rows_local, coef, what) {
    for (k in names(jac)) {
      m <- jac[[k]]
      val <- static_obs_contract_rows(m, R_w[rows_local], coef[rows_local, , drop = FALSE], k, what) / S
      acc <- static_obs_add_loading(acc, k, val, colnames(m))
    }
    acc
  }
  if (kind == "mean") {
    hj <- static_obs_channels(h_jacobian, n, ids, "h_jacobian")
    dw_dh <- if (route == "weighting") phi_w else rep(1, nw)
    acc <- add_channels(acc, hj, seq_len(nw), dw_dh * Fw, "h_jacobian")
  }
  if (uses_q) {
    qj <- static_obs_channels(q_jacobian, n, ids, "q_jacobian")
    rows_q <- which(act_w)                    # floored rows: derivative zero, not read
    coef_q <- if (route == "weighting") h_w * dphi_w * Fw else h_w * resid_w * dphi_w * Gy
    acc <- add_channels(acc, qj, rows_q, coef_q, "q_jacobian")
  }
  if (uses_j) {
    jj <- static_obs_channels(j_jacobian, n, ids, "j_jacobian")
    dY_dj <- if (route == "imputation") rep(1, nw) else 1 - phi_w
    acc <- add_channels(acc, jj, seq_len(nw), h_w * dY_dj * Gy, "j_jacobian")
  }
  loadings <- lapply(names(acc$val), function(k) {
    v <- acc$val[[k]]                         # p x m
    if (kind == "mean") {
      stats::setNames(drop(v[, 1L]), acc$cols[[k]])
    } else {
      matrix(t(v), nrow = ncol(v), dimnames = list(out_cols, acc$cols[[k]]))
    }
  })
  names(loadings) <- names(acc$val)

  full <- NULL
  if (!length(loadings) && (is.null(coef_if) || (is.list(coef_if) && !length(coef_if)))) {
    full <- direct
  } else if (!is.null(coef_if)) {
    contrib <- static_obs_contract(loadings, coef_if, id)
    full <- if (kind == "mean") {
      stats::setNames(direct + drop(contrib), ids)
    } else {
      direct + contrib
    }
  }

  kink <- NULL
  if (uses_q) {
    kink <- list(
      eps = ipw_clip_eps,
      n_q_read = length(q_o),
      n_q_floored = sum(q_o < ipw_clip_eps),
      n_q_at_floor = sum(q_o == ipw_clip_eps),
      n_q_active = sum(q_o > ipw_clip_eps),
      qualification = paste(
        "Floor applied once to the supplied final joint q; derivative zero at or below eps",
        "(strict existing mask), unchanged above. The expansion requires no mass at q* = eps",
        "(Assumption reg-m(f)); exact-boundary rows are reported, not refused."
      )
    )
  }
  components <- NULL
  if (kind == "mean" && route == "augmented") {
    components <- list(
      imputation_part = sum(w * j_w) / S,
      correction_part = sum(w * phi_w * resid_w) / S
    )
  }

  list(
    kind = kind,
    route = route,
    estimate = est,
    denominator_sum = S,
    denominator_mean = S / n,
    pseudo_mean = M,
    weighted_design_mean = Zbar,
    n_units = n,
    n_weight_rows = nw,
    n_observed_rows = length(R_o),
    direct_if = direct,
    loadings = loadings,
    influence = full,
    components = components,
    kink = kink
  )
}

# Normalized reconstructed mean M_hat(h) for one role/target cell and pair
# (eq:M-hat with g = h). `h_jacobian`, and for correction routes `q_jacobian`
# / `j_jacobian`, are required keyed lists (list() = known prediction).
static_obs_mean <- function(id, A, h, route, O = NULL, dy = NULL, q = NULL, j = NULL,
                            h_jacobian = NULL, q_jacobian = NULL, j_jacobian = NULL,
                            coef_if = NULL) {
  static_obs_core("mean", id = id, A = A, O = O, dy = dy, h = h, Z = NULL,
                  route = route, q = q, j = j, q_jacobian = q_jacobian,
                  j_jacobian = j_jacobian, h_jacobian = h_jacobian, coef_if = coef_if)
}

# Reconstructed stayer covariance C_hat between reference Z and the pair change
# (eq:M-hat, tilt identically 1), centered at its route-specific weighted
# design mean, with its IF (eq:chi-C).
static_obs_cov <- function(id, A, Z, route, O = NULL, dy = NULL, q = NULL, j = NULL,
                           q_jacobian = NULL, j_jacobian = NULL, coef_if = NULL) {
  static_obs_core("cov", id = id, A = A, O = O, dy = dy, h = NULL, Z = Z,
                  route = route, q = q, j = j, q_jacobian = q_jacobian,
                  j_jacobian = j_jacobian, h_jacobian = NULL, coef_if = coef_if)
}
