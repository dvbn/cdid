# Bounded, callback-free selection reads. A context is local to one support,
# observer-collection/derivative pass or IF chunk, never stored in a backend or
# used across a caller callback. Its kernels preserve the standalone argument,
# clipping, finite, boundary and reference-support checks. Only repeated source
# validation moves to entry and the exact-source successful-return boundary.
#
# The closure/snapshot holds borrowed immutable source references, no copied
# owner arrays or probabilities. Its fixed binding/function/list headers fit
# the existing fit 8192-slot and replay dependency 4096-slot metadata reserves.
# Prediction/derivative row buffers keep their existing admitted envelopes.
selection_prediction_operation <- function(backend, expected_sources = NULL,
                                           source_reader = NULL, source_reference = NULL, pair_response_reader = NULL) {
  # Resolve caller promises before granting any operation-local permission.
  force(backend); force(expected_sources); force(source_reader); force(source_reference); force(pair_response_reader)
  validate_selection_models(backend)
  reader <- get("cp04_reader", backend$reader_ref, inherits = FALSE)
  frame <- get("frame_owner", reader$frame_ref, inherits = FALSE)
  panel <- get("panel_owner", reader$panel_ref, inherits = FALSE)
  reference <- if (is.null(backend$reference_ref)) NULL else
    get("reference_reader", backend$reference_ref, inherits = FALSE)
  if (!identical(reader$ref, backend$reader_ref) || !identical(frame$ref, reader$frame_ref) ||
      !identical(panel$ref, reader$panel_ref) || !identical(frame$panel_ref, panel$ref) ||
      (!is.null(reference) && !identical(reference$ref, backend$reference_ref)))
    cp04_abort("Selection operation has foreign source bindings.", "integrity")
  if ((!is.null(source_reader) && !gmm_sexp_aliases(source_reader, reader)) ||
      (!is.null(source_reference) && (is.null(reference) || !gmm_sexp_aliases(source_reference, reference))))
    cp04_abort("Selection operation was supplied a foreign reader or reference source.", "integrity")
  st <- reader$state; sel_fit <- backend; ref <- reference
  bindings <- c(backend = "selection_models", reader = "cp04_reader", frame = "frame_owner",
    panel = "panel_owner", reference = "reference_reader")
  check_sources <- function(src) {
    if (!is.environment(src) || !identical(parent.env(src), emptyenv()) ||
        !environmentIsLocked(src) || !identical(sort(ls(src, all.names = TRUE)), sort(names(bindings))))
      cp04_abort("Selection operation lacks its exact source snapshot.", "integrity")
    for (nm in names(bindings)) {
      if (!bindingIsLocked(nm, src) || bindingIsActive(nm, src))
        cp04_abort("Selection operation snapshot has a mutable or active binding.", "integrity")
      original <- get(nm, src, inherits = FALSE)
      if (is.null(original)) next
      ref <- original$ref; binding <- bindings[[nm]]
      if (!is.environment(ref) || !identical(parent.env(ref), emptyenv()) || !environmentIsLocked(ref) ||
          !exists(binding, ref, inherits = FALSE) || !bindingIsLocked(binding, ref) || bindingIsActive(binding, ref) ||
          !gmm_sexp_aliases(original, get(binding, ref, inherits = FALSE)))
        cp04_abort("Selection operation source changed before its successful boundary.", "integrity")
    }
    invisible(NULL)
  }
  # A fresh operation after GLM/warning handlers validates again, then compares
  # the original immutable inputs. A snapshot grants no validation permission.
  if (!is.null(expected_sources)) {
    check_sources(expected_sources)
    for (nm in names(bindings)) {
      current <- get(nm, inherits = FALSE)
      original <- get(nm, expected_sources, inherits = FALSE)
      if (!identical(is.null(current), is.null(original)) ||
          (!is.null(current) && !gmm_sexp_aliases(current, original)))
        cp04_abort("Selection operation restart has different source identities.", "integrity")
    }
  }
  sources <- new.env(parent = emptyenv())
  for (nm in names(bindings)) {
    assign(nm, get(nm, inherits = FALSE), sources)
    lockBinding(nm, sources)
  }
  lockEnvironment(sources, bindings = TRUE)
  check_sources(sources)
  active <- TRUE
  guard <- function() {
    if (!active) cp04_abort("Selection operation is closed.", "integrity")
    invisible(NULL)
  }
  cell_prediction <- function(rows, cell_id) {
    guard()
    if (!is.numeric(rows) || is.object(rows) || !is.null(dim(rows)) || !length(rows) ||
        anyNA(rows) || any(rows != floor(rows)) || any(rows < 1L | rows > length(st$a_index)))
      cp04_abort("Invalid row indices for predict_selection_cell.", "input")
    if (!is.numeric(cell_id) || length(cell_id) != 1L || is.na(cell_id) ||
        cell_id < 1L || cell_id > length(backend$cells))
      cp04_abort("Invalid `cell_id`.", "input")

    cell <- backend$cells[[as.integer(cell_id)]]
    block <- backend$blocks[[cell$block_id]]
    b <- cell$local_b

    if (identical(block$branch, "zero_col")) {
      return(rep.int(0.5, length(rows)))
    }
    if (identical(block$branch, "cp04_separate_pair")) {
      sp <- list(pair = cell$pair, cohort = cell$cohort)
      if (!is.null(cell$reference)) sp$reference <- cell$reference
      return(cp04_predict_selection(reader, rows, sp, block, ref))
    }
    d_s <- block$d_s
    Z_slope <- if (d_s > 0L) {
      if (is.null(ref)) st$Z[rows, backend$design_spec$slope_cols, drop = FALSE] else
        chained_reference_s(ref, cell$reference, rows)[, backend$design_spec$slope_cols, drop = FALSE]
    } else matrix(0, length(rows), 0L)
    eta <- (if (block$has_intercept) block$alpha[b] else 0) +
      (if (d_s > 0L) drop(Z_slope %*% block$beta) else rep.int(0, length(rows)))
    p <- ipw_clip_logit(eta)$prob
    cp04_finite(p, "replayed selection probabilities", TRUE)
    if (any(p >= 1)) cp04_abort("Replayed selection probability at boundary.", "boundary", cell)
    p
  }
  period_prediction <- function(rows, pair_cohort_id) {
    guard()
    if (!identical(backend$selection_model, "independent_period_logit"))
      cp04_abort("`predict_period_pair_selection` requires selection_model == \"independent_period_logit\".", "unsupported")
    q <- resolve_pair_cohort_id(backend, pair_cohort_id)
    pc <- backend$pair_cohorts[[q]]

    ps <- cell_prediction(rows, pc$s_cell_id)
    pt <- cell_prediction(rows, pc$t_cell_id)
    p_pair <- ps * pt
    cp04_finite(p_pair, "replayed product pair selection probabilities", TRUE)
    if (any(p_pair >= 1)) cp04_abort("Replayed product pair selection probability at boundary.", "boundary", pc)
    a_pc <- if (is.null(ref)) st$a[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
    w_rows <- if (is.null(ref)) st$w[rows] else chained_reference_weight(ref, pc$reference, rows)
    # The final raw product is floored once; raw endpoints stay for scores.
    p_weight <- ipw_clip_floor(p_pair)
    w_sel <- w_rows * a_pc / (st$rho[rows] * p_weight)
    cp04_finite(w_sel, "replayed selection weights", TRUE)
    list(prob = p_pair, prob_s = ps, prob_t = pt, prob_weight = p_weight,
         active = ipw_clip_floor_active(p_pair), weight = w_sel, a = a_pc, pair_cohort_id = q)
  }
  pair_prediction <- function(rows, pair_cohort_id) {
    guard()
    if (identical(backend$selection_model, "independent_period_logit")) {
      return(period_prediction(rows, pair_cohort_id))
    }
    q <- resolve_pair_cohort_id(backend, pair_cohort_id)
    pc <- backend$pair_cohorts[[q]]

    p_pair <- cell_prediction(rows, pc$pair_cell_id)
    a_pc <- if (is.null(ref)) st$a[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
    w_rows <- if (is.null(ref)) st$w[rows] else chained_reference_weight(ref, pc$reference, rows)
    p_weight <- ipw_clip_floor(p_pair)
    w_sel <- w_rows * a_pc / (st$rho[rows] * p_weight)
    cp04_finite(w_sel, "replayed selection weights", TRUE)
    list(prob = p_pair, prob_weight = p_weight, active = ipw_clip_floor_active(p_pair),
         weight = w_sel, a = a_pc, pair_cohort_id = q)
  }
  deriv <- function(up_map, q, rows, pred = NULL, Z = NULL, a_vec = NULL) {
    guard()
    pc <- sel_fit$pair_cohorts[[q]]
    ref <- NULL
    if (is.null(Z) || is.null(a_vec)) {
      ref <- reference
      if (is.null(ref)) {
        Z <- st$Z
        a_vec <- st$a
      }
    }
    a_pc <- if (is.null(ref)) a_vec[pc$cohort] else chained_reference_stat(ref, pc$reference, pc$cohort, "a")
    ch <- up_map$pair_channels[[q]]
    dq <- length(ch) - 1L
    n_r <- length(rows)

    if (!n_r) {
      return(list(
        pred = list(prob = numeric(0), weight = numeric(0), a = a_pc),
        T_logit = matrix(0, 0L, dq),
        log_deriv = matrix(0, 0L, length(ch))
      ))
    }

    if (is.null(pred)) {
      pred <- pair_prediction(rows, q)
    }
    act <- ipw_clip_floor_active(pred$prob)
    d_s <- sel_fit$design_spec$d_s
    Z_slope <- if (d_s > 0L) {
      if (is.null(ref)) Z[rows, sel_fit$design_spec$slope_cols, drop = FALSE] else
        chained_reference_s(ref, pc$reference, rows)[, sel_fit$design_spec$slope_cols, drop = FALSE]
    } else matrix(0, n_r, 0L)

    if (identical(sel_fit$selection_model, "pair_logit")) {
      m <- pc$block_id
      bl <- sel_fit$blocks[[m]]
      if (identical(bl$branch, "zero_col")) {
        T_logit <- matrix(0, n_r, 0L)
      } else {
        w_p <- -((1 - pred$prob) * act)
        T_logit <- cbind(
          if (bl$has_intercept) matrix(w_p, n_r, 1L) else matrix(0, n_r, 0L),
          if (d_s > 0L) w_p * Z_slope else matrix(0, n_r, 0L)
        )
      }
    } else {
      ms <- pc$s_block_id; bl_s <- sel_fit$blocks[[ms]]
      mt <- pc$t_block_id; bl_t <- sel_fit$blocks[[mt]]
      if (identical(bl_s$branch, "zero_col")) {
        T_logit <- matrix(0, n_r, 0L)
      } else {
        w_s <- -((1 - pred$prob_s) * act)
        w_t <- -((1 - pred$prob_t) * act)
        if (ms != mt) {
          T_s <- cbind(
            if (bl_s$has_intercept) matrix(w_s, n_r, 1L) else matrix(0, n_r, 0L),
            if (d_s > 0L) w_s * Z_slope else matrix(0, n_r, 0L)
          )
          T_t <- cbind(
            if (bl_t$has_intercept) matrix(w_t, n_r, 1L) else matrix(0, n_r, 0L),
            if (d_s > 0L) w_t * Z_slope else matrix(0, n_r, 0L)
          )
          T_logit <- cbind(T_s, T_t)
        } else {
          T_logit <- cbind(
            if (bl_s$has_intercept) cbind(w_s, w_t) else matrix(0, n_r, 0L),
            if (d_s > 0L) (w_s + w_t) * Z_slope else matrix(0, n_r, 0L)
          )
        }
      }
    }

    log_deriv <- cbind(rep.int(1 / a_pc, n_r), T_logit)
    list(pred = pred, T_logit = T_logit, log_deriv = log_deriv)
  }
  upstream <- function(up_map, rows, uids, details = FALSE, preparation = NULL) {
    guard()
    prepared <- if (is.null(preparation)) NULL else
      dr_conditional_replay_upstream_index(preparation, sel_fit, up_map, uids, details)
    ref <- reference
    B <- length(rows)
    ci <- st$cohort_index[rows]
    # D60: no common weight; each rate coordinate and cell reads its recipe's W_g.
    w <- if (is.null(ref)) st$w[rows] else rep(NA_real_, B)
    U_rows <- st$U[rows]
    A <- as.double(st$a_index[rows] > 0L)
    w_of <- function(r, ix) if (is.null(ref)) w[ix] else chained_reference_weight(ref, r, rows[ix])

    up <- matrix(0, B, length(uids))
    for (j in if (is.null(prepared)) which(uids <= up_map$n_a) else prepared$rate_positions) {
      k <- up_map$cohorts_a[uids[j]]
      r <- up_map$rates_ref[uids[j]]
      ix <- which(ci == k & U_rows)
      if (length(ix)) {
        a_k <- if (is.null(ref)) st$a[k] else chained_reference_stat(ref, r, k, "a")
        m_k <- if (is.null(ref)) st$mass_U[k] else chained_reference_stat(ref, r, k, "mass_U")
        up[ix, j] <- w_of(r, ix) * (A[ix] - a_k) / m_k
      }
    }

    active_blocks <- if (details) {
      seq_along(sel_fit$blocks)
    } else if (is.null(prepared)) {
      which(vapply(up_map$theta_coords, function(tc) any(tc %in% uids), logical(1)))
    } else prepared$active_blocks

    d_s <- sel_fit$design_spec$d_s
    # Selection covariates are defined only on the declared fitting universe.
    # Outside U the nuisance IF is structurally zero; do not read an NA covariate
    # there and then multiply it by a zero score.
    Z_slope <- if (is.null(ref) && d_s > 0L) st$Z[rows[U_rows], sel_fit$design_spec$slope_cols, drop = FALSE] else
      matrix(0, sum(U_rows), d_s)
    block_ifs <- vector("list", length(sel_fit$blocks))
    active_cells <- if (is.null(prepared)) integer(0) else prepared$active_cells
    u_pos <- cumsum(U_rows)

    for (m in active_blocks) {
      bl <- sel_fit$blocks[[m]]
      cids <- bl$cell_ids
      B_m <- bl$B_m
      if (is.null(prepared)) active_cells <- c(active_cells, cids)
      cell_scores <- matrix(0, B, B_m)
      # A D60 block lives inside one recipe: its slope rows are that recipe's
      # selection design on the block cohorts' U rows (zero scores elsewhere).
      Zb <- Z_slope
      if (!is.null(ref) && d_s > 0L) {
        rb <- sel_fit$block_specs[[m]]$reference
        cohort_ids <- if (is.null(prepared)) vapply(sel_fit$cells[cids], `[[`, integer(1), "cohort") else prepared$cohorts[[m]]
        on <- which(U_rows & ci %in% cohort_ids)
        if (length(on)) Zb[u_pos[on], ] <- chained_reference_s(ref, rb, rows[on])[, sel_fit$design_spec$slope_cols, drop = FALSE]
      }
      for (b in seq_len(B_m)) {
        cid <- cids[b]
        cell <- sel_fit$cells[[cid]]
        rr <- which(ci == cell$cohort & U_rows)
        if (!length(rr)) next
        pr <- cell_prediction(rows[rr], cid)
        Y <- if (identical(sel_fit$selection_model, "independent_period_logit")) {
          cp04_date_rows(reader, rows[rr], cell$date_spec)$R
        } else {
          if (is.null(pair_response_reader)) cp04_pair_rows(reader, rows[rr], cell$pair)$S else
            pair_response_reader(rows[rr], cell$pair, FALSE)$S
        }
        cell_scores[rr, b] <- w_of(cell$reference, rr) * (Y - pr)
      }
      bif <- list(
        xi_alpha = matrix(0, B, if (bl$has_intercept) B_m else 0L),
        xi_beta = matrix(0, B, bl$d_s),
        xi_theta = matrix(0, B, length(up_map$theta_coords[[m]]))
      )
      if (any(U_rows)) {
        on_U <- evaluate_selection_block_if(bl, Zb, cell_scores[U_rows, , drop = FALSE], st$N)
        for (name in names(bif)) bif[[name]][U_rows, ] <- on_U[[name]]
      }
      block_ifs[[m]] <- bif
      tc <- up_map$theta_coords[[m]]
      keep_tc <- if (is.null(prepared)) which(tc %in% uids) else prepared$keep[[m]]
      if (length(keep_tc)) {
        destination <- if (is.null(prepared)) match(tc[keep_tc], uids) else prepared$destination[[m]]
        up[, destination] <- bif$xi_theta[, keep_tc, drop = FALSE]
      }
    }

    list(
      up = up,
      block_ifs = block_ifs,
      active_blocks = active_blocks,
      active_cells = unique(active_cells)
    )
  }
  finish <- function() {
    guard()
    active <<- FALSE
    check_sources(sources)
    invisible(NULL)
  }
  abort <- function() { active <<- FALSE; invisible(NULL) }
  list(cell = cell_prediction, period = period_prediction, pair = pair_prediction, deriv = deriv, upstream = upstream,
    finish = finish, abort = abort, sources = sources)
}
