# Private: opt-in fitted-nuisance capture for cross-pair chained DiD
# (cross_pair = TRUE). Owned with R/chained_cross_pair.R and
# R/chained_cross_pair_backends.R; request wiring and the public accessor are
# owned by the request layer.
#
#   res <- run_chained_cross_pair_fit(..., keep_nuisance = TRUE)
#   cap <- chained_cross_pair_nuisance_capture(res, prep, gate, est_label, weightsname)
#
# keep_nuisance = FALSE (default) creates no sink: every hook in the cores is one
# `if (!is.null(nuisance_sink))` branch, the three fit-record retention
# conditions read `compute_if || !is.null(nuisance_sink)` and
# fit_checked_propensity_logit(keep_coef = FALSE) returns its unchanged list.
# With TRUE, a sink environment rides
# along the actual fit. Hooks read the fitted objects of the core frame at the
# point they exist (no refit, no second solve) and store compact records only:
# estimates, diagonal coefficient variances, labels and route diagnostics. The
# sink is reduced to the plain list res$nuisance_capture_input before return.
# Estimation, covariance, RNG use and returned estimates are unchanged.
#
# Covered routes and fits (static designs; dated/recipe reference routes return
# a typed status):
#   baseline (evaluate_cross_pair_moments; A units, normalizer n_A)
#     treatment  checked logit per (g, k, mask): target rows w_i, control rows
#                w_i / r with r the observed share of cohort-k weight (generated)
#     outcome    checked WLS per (k, mask, (s,t)) of the actual change on
#                observed control rows, fixed weights w_i
#   weighted core (selection "none"/"ipw"; frame rows, normalizer N)
#     selection  shared cdid_selection_models fit (post-fit bounded pass through
#                the shared selection tables with the captured upstream map)
#     treatment  cp04_logit per treatment key; weights fw (rates or a/(rho p))
#     outcome    or_outcome_fit_blocks groups of the actual change, weights fw
#   reconstruction core (selection_est_method "imputation"/"dr")
#     imputation actual-change WLS groups on observed pairs, fixed weights
#     treatment  synthetic risk-set logits, weights w_i a_k
#     outcome    OR groups on the reconstructed response Z_c, weights w_i a_k
#     selection  augmented observation logits (dr only), shared tables
#
# Coefficient SEs evaluate the estimator's own reverse-mode IF operator of each
# fit at unit sensitivities: the per-unit score through the fit's own stored
# factor (H_chol of H/N, the WLS Gram solve, the pivoted QR), plus every
# upstream channel the estimator propagates for that fit (observation rates,
# repeat rates, selection log-derivatives against up_if / Xi_theta, and, for Z_c,
# the imputation IF). Only the diagonal is kept:
#   Var(theta_j) = N^{-2} sum_g (S_gj - Sbar_j)^2,
# g = original unit or the requested cluster, centered, no small-sample factor
# (chained_nuisance_moments_*). With compute_if = FALSE the route forms no IF and
# SEs are typed unavailable; coefficients are still reported and checked.

cross_pair_nuisance_tol <- c(wls = 1e-8, logit = 1e-5)

# ------------------------------------------------------------------------------
# 1. Sink and pure helpers
# ------------------------------------------------------------------------------
cross_pair_nuisance_sink <- function(cluster_full, frame_ids, a_to_full, control) {
  s <- new.env(parent = emptyenv())
  s$control <- control
  s$n_frame <- length(frame_ids)
  s$a_to_full <- as.integer(a_to_full)
  s$frame <- NA_character_
  s$clustered <- !is.null(cluster_full)
  s$group_full <- NULL
  s$n_groups <- NA_integer_
  s$cluster_code <- NULL
  if (s$clustered) {
    if (length(cluster_full) != length(frame_ids) || anyNA(cluster_full)) {
      s$cluster_code <- "cluster_alignment"
    } else {
      s$group_full <- match(cluster_full, unique(cluster_full))
      s$n_groups <- max(s$group_full)
      if (s$n_groups < 2L) s$cluster_code <- "insufficient_clusters"
    }
  }
  s$route <- list()
  s$treatment <- list(); s$outcome <- list(); s$imputation <- list()
  s$notes <- list()
  s$selection <- NULL
  s$work <- 0; s$entries <- 0
  s
}

# Rows of `parts` (list(rows, M)) summed per unit into one dense block.
cross_pair_nuisance_dense <- function(parts, p) {
  parts <- Filter(function(x) length(x$rows) > 0L, parts)
  rows <- sort(unique(as.integer(unlist(lapply(parts, `[[`, "rows"), use.names = FALSE))))
  M <- matrix(0, length(rows), p)
  for (pt in parts) {
    m <- pt$M
    if (!is.matrix(m)) m <- matrix(as.double(m), nrow = length(pt$rows), ncol = p)
    ix <- match(as.integer(pt$rows), rows)
    if (anyDuplicated(ix)) {
      m <- rowsum(m, ix)
      ix <- as.integer(rownames(m))
    }
    M[ix, ] <- M[ix, , drop = FALSE] + m
  }
  list(rows = rows, M = M)
}

# Channel part: unit IF of an upstream quantity (F, rows x q) times the
# derivative of the fit's normalized estimating equation (A, p x q).
cross_pair_nuisance_channel <- function(rows, F, A) {
  F <- as.matrix(F); A <- as.matrix(A)
  if (!length(rows) || !ncol(F)) return(NULL)
  list(rows = as.integer(rows), M = F %*% t(A))
}

# Upstream IF columns of an estimated cell rate r = sum_obs w / M_k on the
# unit frame, as the estimator contracts it: (N / M_k)(w_i S_i - r w_i 1{i in k}).
cross_pair_nuisance_rate_if <- function(rows_k, w_k, rows_obs, w_obs, M_k, r_raw, N) {
  inv <- as.double(N) / M_k
  list(rows = c(as.integer(rows_k), as.integer(rows_obs)),
       f = c(-(w_k * inv) * r_raw, w_obs * inv))
}

# Relative estimating-equation residual at the reported coefficients.
cross_pair_nuisance_score_check <- function(S) {
  S <- as.matrix(S)
  if (!length(S)) return(0)
  max(abs(colSums(S))) / max(1, max(colSums(abs(S))))
}

# Sum the block-local intercept scores and the common slope scores without
# filling a dense column of zeros for every other intercept on each row.
cross_pair_nuisance_pooled_check <- function(parts, group, intercept, width, p) {
  total <- absolute <- numeric(width + p)
  slopes <- width + seq_len(p)
  for (j in seq_along(parts)) {
    a <- intercept[[j]]
    if (length(a)) {
      b <- group[j]
      total[b] <- total[b] + sum(a)
      absolute[b] <- absolute[b] + sum(abs(a))
    }
    if (p) {
      total[slopes] <- total[slopes] + colSums(parts[[j]])
      absolute[slopes] <- absolute[slopes] + colSums(abs(parts[[j]]))
    }
  }
  if (length(total)) max(abs(total)) / max(1, absolute) else 0
}

cross_pair_nuisance_group <- function(s, rows) {
  if (is.null(s$group_full)) return(NULL)
  pos <- if (identical(s$frame, "A")) s$a_to_full[rows] else rows
  s$group_full[pos]
}

cross_pair_nuisance_se <- function(var = NULL, clusters = NA_integer_, code = NULL, reason = NULL)
  list(var = var, clusters = as.integer(clusters), code = code, reason = reason)

# Diagonal moments of one fit's coefficient IF (rows unique, xi rows x p).
cross_pair_nuisance_variance <- function(s, rows, xi, N) {
  p <- ncol(xi)
  if (s$clustered && !is.null(s$cluster_code))
    return(cross_pair_nuisance_se(code = s$cluster_code, reason = if (identical(s$cluster_code, "cluster_alignment"))
      "Cluster labels could not be aligned with the frame by position." else
      "At least two clusters are needed for a clustered coefficient standard error."))
  if (!p) return(cross_pair_nuisance_se(var = numeric(0)))
  if (!length(rows)) return(cross_pair_nuisance_se(var = rep(0, p), clusters = if (s$clustered) 0L else NA_integer_))
  grp <- cross_pair_nuisance_group(s, rows)
  ng <- if (is.null(grp)) N else s$n_groups
  v <- tryCatch({
    acc <- chained_nuisance_moments_new(p, if (is.null(grp)) NULL else ng)
    acc <- chained_nuisance_moments_add(acc, xi, grp)
    chained_nuisance_moments_variance(acc, N, ng)
  }, cdid_cp04_integrity_error = function(e) NULL)
  if (is.null(v))
    return(cross_pair_nuisance_se(code = "if_nonfinite", reason = "The coefficient IF columns were not finite."))
  if (is.null(grp)) return(cross_pair_nuisance_se(var = v))
  hit <- rowSums(xi != 0) > 0
  cl <- length(unique(grp[hit]))
  if (cl < 2L)
    return(cross_pair_nuisance_se(clusters = cl, code = "insufficient_clusters",
      reason = "This fit has fewer than two clusters containing nonzero unit coefficient influence contributions."))
  cross_pair_nuisance_se(var = v, clusters = cl)
}

# Admission of one fit's IF evaluation on the core's owner control.
cross_pair_nuisance_admit <- function(s, ctrl, live, scratch, work, label) {
  ok <- tryCatch({
    cp04_budget(ctrl, c(live, s$entries), scratch, c(s$work, work), s$entries, label)
    TRUE
  }, cdid_ipw_owner_bound_error = function(e) conditionMessage(e))
  if (isTRUE(ok)) s$work <- s$work + work
  ok
}

cross_pair_nuisance_no_if <- function()
  cross_pair_nuisance_se(code = "influence_not_computed", reason = paste(
    "This call formed no influence functions (compute_inffunc = FALSE with a weighting other than \"optimal\"),",
    "so the coefficient IF of this fit was not evaluated."))

# Common SE driver: admission, IF assembly (closure), moments.
cross_pair_nuisance_se_run <- function(s, ctrl, live, compute_if, n_rows_bound, p, width, N, label, assemble) {
  if (!isTRUE(compute_if)) return(cross_pair_nuisance_no_if())
  scratch <- 4 * as.double(n_rows_bound) * (p + width + 2) + 8192
  work <- 64 * as.double(n_rows_bound) * (p + width + 2) * (p + 2) + 4096
  ok <- cross_pair_nuisance_admit(s, ctrl, live, scratch, work, label)
  if (!isTRUE(ok))
    return(cross_pair_nuisance_se(code = "budget", reason = paste(
      "The coefficient IF evaluation exceeds the declared workspace or work budget:", ok)))
  d <- assemble()
  cross_pair_nuisance_variance(s, d$rows, d$M, N)
}

cross_pair_nuisance_charge <- function(s, rec) {
  s$entries <- s$entries + 512 + 16 * length(rec$estimate) + 4 * length(rec$links)
  rec
}

cross_pair_nuisance_terms <- function(X, prefix = "x") {
  cn <- colnames(X)
  if (is.null(cn) || anyNA(cn) || !all(nzchar(cn))) paste0(prefix, seq_len(ncol(X))) else cn
}

# ------------------------------------------------------------------------------
# 2. Equations (compact: no multiplication by a literal 1)
# ------------------------------------------------------------------------------
cross_pair_nuisance_rhs <- function(kept, intercept_sym = NULL, L = 1) {
  sl <- chained_nuisance_slope_text(kept)
  one <- is.numeric(L) && length(L) == 1L && isTRUE(L == 1)
  slope <- if (!length(sl)) character(0) else if (one) sl else
    paste0(if (is.numeric(L)) chained_nuisance_num(L) else L, "*(", paste(sl, collapse = " + "), ")")
  terms <- c(intercept_sym, slope)
  if (length(terms)) paste(terms, collapse = " + ") else "0"
}

cross_pair_nuisance_logit_rhs <- function(terms, int_name = "(Intercept)") {
  ct <- vapply(seq_along(terms), function(j)
    if (terms[j] %in% int_name) paste0("c", j - 1L) else paste0("c", j - 1L, "*", terms[j]), "")
  if (length(ct) > 7L) ct <- c(ct[1:6], paste0("... (", length(ct) - 6L, " more terms)"))
  if (length(ct)) paste(ct, collapse = " + ") else "0"
}

cross_pair_nuisance_delta_lhs <- function(k, s, t) {
  if (length(s) == 1L && length(t) == 1L)
    paste0("E[Y_i", chained_nuisance_num(t), " - Y_i", chained_nuisance_num(s), " | X_i, G_i = ",
      chained_nuisance_gtxt(k), ", S_i(", chained_nuisance_num(s), ",", chained_nuisance_num(t), ") = 1]") else
    paste0("E[Y_it - Y_is | X_i, G_i = ", if (length(unique(k)) == 1L) chained_nuisance_gtxt(k[1L]) else "k",
      ", S_i(s,t) = 1]")
}

# ------------------------------------------------------------------------------
# 3. Baseline route hooks (evaluate_cross_pair_moments; A-unit frame, N = n_A)
# ------------------------------------------------------------------------------
cross_pair_nuisance_baseline_live <- function(env) {
  n_a <- as.double(env$n_a); L <- as.double(env$L)
  n_a * L * (1 + !is.null(get0("Omega_scaled", envir = env, inherits = FALSE))) + 1024
}

# Inline: right after the estimator's own solve gamma_piv and residual e_c.
cross_pair_nuisance_baseline_or <- function(s, env) {
  key <- paste0("O|", env$qr_key, "|", env$cs, "|", env$ct)
  ell <- as.integer(env$ell)
  gam <- as.double(env$gamma_piv)
  old <- s$outcome[[key]]
  if (!is.null(old)) {
    old$links <- c(old$links, ell)
    old$reuse_max_diff <- max(old$reuse_max_diff, abs(old$estimate_piv - gam))
    s$outcome[[key]] <- old
    return(invisible(NULL))
  }
  s$frame <- "A"
  qr_rec <- env$qr_rec; piv <- qr_rec$piv; rk <- qr_rec$rk
  cn <- cross_pair_nuisance_terms(env$Z_or)
  X <- qr_rec$X_c[, piv, drop = FALSE]
  w <- as.double(qr_rec$w_c); e <- as.double(env$e_c); n_a <- as.double(env$n_a)
  rows <- as.integer(env$idx_c)
  ord <- order(piv)
  S <- (w * e) * X
  ld <- env$links_df
  se <- cross_pair_nuisance_se_run(s, s$control, cross_pair_nuisance_baseline_live(env), env$compute_if,
    length(rows), rk, 0, n_a, "nuisance capture baseline outcome IF", function() {
      xi <- n_a * t(cp04_solve(qr_rec$R_rk, t(S)))
      list(rows = rows, M = xi[, ord, drop = FALSE])
    })
  int <- cn[piv][ord] == "(Intercept)"
  rec <- list(role = "outcome", key = key, links = ell,
    estimate = gam[ord], estimate_piv = gam, reuse_max_diff = 0,
    term = cn[piv][ord], coef_role = ifelse(int, "intercept", "slope"), cell = rep(NA_integer_, rk),
    kept = cn[piv][ord][!int], aliased = cn[qr_rec$drop_cols], intercept = any(int),
    se = se, check = cross_pair_nuisance_score_check(S), check_kind = "wls",
    meta = list(kind = "baseline_or", cohort = ld$cohort[ell], ci = env$ci, mask = env$mid,
      s = ld$s[ell], t = ld$t[ell], n_rows = length(rows), weight_sum = sum(w), rank = rk, L = 1))
  s$outcome[[key]] <- cross_pair_nuisance_charge(s, rec)
  invisible(NULL)
}

# End of core (and the unconditional early return).
cross_pair_nuisance_baseline_collect <- function(s, env, unconditional = FALSE) {
  s$frame <- "A"
  s$route <- list(core = "baseline", est_method = env$est_method, n_links = as.integer(env$L),
    N_if = as.double(env$n_a), compute_if = isTRUE(env$compute_if), unconditional = isTRUE(unconditional),
    uncond_treatment = isTRUE(env$uncond_trt), uncond_outcome = isTRUE(env$uncond_or),
    n_cell_means = if (isTRUE(unconditional)) as.integer(env$n_u) else NA_integer_)
  if (isTRUE(unconditional)) return(invisible(NULL))
  if (identical(env$est_method, "reg")) return(invisible(NULL))
  ld <- env$links_df; u_g <- env$u_g; n_a <- as.double(env$n_a)
  ctrl_l <- which(ld$link_role != "target")
  key_of <- paste0(match(ld$g, u_g), ":", ld$coh_idx, ":", ld$coh_mask_id)
  live <- cross_pair_nuisance_baseline_live(env)
  n_closed <- 0L
  for (lk in ls(env$logit_cache, all.names = TRUE)) {
    lg <- env$logit_cache[[lk]]
    use <- ctrl_l[key_of[ctrl_l] == lk]
    if (isTRUE(lg$uncond)) { n_closed <- n_closed + 1L; next }
    X <- lg$Z_u; p <- ncol(X); rows <- as.integer(lg$idx_u)
    D <- as.double(lg$D_u); is_c <- D == 0
    beta <- as.double(lg$beta)
    pb <- stats::plogis(as.double(X %*% beta))
    S <- (lg$w_u * (D - lg$p_u)) * X
    se <- cross_pair_nuisance_se_run(s, s$control, live, env$compute_if,
      length(rows) + length(lg$idx_k_all), p, 1, n_a, "nuisance capture baseline treatment IF", function() {
        parts <- list(list(rows = rows, M = S))
        if (!is.null(lg$A_gamma_r)) {
          ri <- cross_pair_nuisance_rate_if(lg$idx_k_all, env$w_a[lg$idx_k_all], rows[is_c], env$w_a[rows[is_c]],
            lg$M_k, lg$r_raw, n_a)
          parts[[2L]] <- cross_pair_nuisance_channel(ri$rows, matrix(ri$f, ncol = 1L), matrix(lg$A_gamma_r, ncol = 1L))
        }
        d <- cross_pair_nuisance_dense(parts, p)
        d$M <- t(cp04_solve(lg$H_chol, t(d$M)))
        d
      })
    gi <- as.integer(strsplit(lk, ":", fixed = TRUE)[[1L]])
    cn <- cross_pair_nuisance_terms(X)
    rec <- list(role = "treatment", key = paste0("T|", lk), links = as.integer(use),
      estimate = beta, term = cn, coef_role = ifelse(cn == "(Intercept)", "intercept", "slope"),
      cell = rep(NA_integer_, p), kept = cn, aliased = character(0), intercept = any(cn == "(Intercept)"),
      se = se, check = cross_pair_nuisance_score_check((lg$w_u * (D - pb)) * X), check_kind = "logit",
      diagnostics = lg$diagnostics,
      meta = list(kind = "baseline_trt", g = u_g[gi[1L]], cohort = env$coh_vals[gi[2L]], mask = gi[3L],
        pairs = unique(paste0("(", chained_nuisance_num(ld$s[use]), ",", chained_nuisance_num(ld$t[use]), ")")),
        n_target = sum(!is_c), n_control = sum(is_c), w_target = sum(lg$w_u[!is_c]),
        w_control = sum(lg$w_u[is_c]), r_raw = lg$r_raw, rate_active = !is.null(lg$A_gamma_r),
        n_rows = length(rows), weight_sum = sum(lg$w_u)))
    s$treatment[[rec$key]] <- cross_pair_nuisance_charge(s, rec)
  }
  s$route$n_closed_form_treatment <- n_closed
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 4. Weighted core hooks (selection "none" / "ipw"; frame rows, N = st$N)
# ------------------------------------------------------------------------------
cross_pair_nuisance_general_live <- function(env, D = 0) {
  N <- as.double(env$N); L <- as.double(env$L)
  N * L + N * D + 1024
}

# Upstream (selection / repeat-rate) channel of rows with log-derivatives `ld`
# on coordinates `ch`: A[, ch] += S' ld / N.
cross_pair_nuisance_ld_add <- function(A, S, ch, ld, N) {
  if (!length(ch) || is.null(ld) || !nrow(S)) return(A)
  A[, ch] <- A[, ch, drop = FALSE] + crossprod(S, as.matrix(ld)) / N
  A
}

cross_pair_nuisance_up_part <- function(A, up_if) {
  if (is.null(A) || is.null(up_if)) return(NULL)
  cols <- which(colSums(A != 0) > 0)
  if (!length(cols)) return(NULL)
  F <- up_if[, cols, drop = FALSE]
  nz <- which(rowSums(F != 0) > 0)
  cross_pair_nuisance_channel(nz, F[nz, , drop = FALSE], A[, cols, drop = FALSE])
}

# Rate channels: one p-vector per active rate id, against its rate IF.
cross_pair_nuisance_rate_parts <- function(Ar, rate_list, N) {
  out <- list()
  for (rid in names(Ar)) {
    rc <- rate_list[[as.integer(rid)]]
    ri <- cross_pair_nuisance_rate_if(rc$rows_A, rc$w_A, rc$rows_obs, rc$w_obs, rc$M_k, rc$r_raw, N)
    out[[length(out) + 1L]] <- cross_pair_nuisance_channel(ri$rows, matrix(ri$f, ncol = 1L), matrix(Ar[[rid]], ncol = 1L))
  }
  out
}

cross_pair_nuisance_rate_add <- function(Ar, rid, a) {
  k <- as.character(rid)
  Ar[[k]] <- if (is.null(Ar[[k]])) a else Ar[[k]] + a
  Ar
}

# Inline: after the treatment fit for `tkey` is obtained (created or cached).
cross_pair_nuisance_weighted_trt <- function(s, env, tkey) {
  key <- paste0("T|", tkey)
  ell <- as.integer(env$ell)
  if (!is.null(s$treatment[[key]])) {
    s$treatment[[key]]$links <- c(s$treatment[[key]]$links, ell)
    return(invisible(NULL))
  }
  s$frame <- "F"
  tf <- env$tfit; d <- env$d; N <- as.double(env$N)
  X <- tf$X_all; p <- ncol(X); rows <- as.integer(tf$rows_all)
  beta <- as.double(tf$coefficients)
  S <- (tf$w_all * (tf$y_all - tf$prob_all)) * X
  pb <- stats::plogis(as.double(X %*% beta))
  sel_mode <- env$sel_mode
  up_if <- get0("up_if", envir = env, inherits = FALSE)
  D_up <- if (is.null(up_if)) 0 else ncol(up_if)
  c_cells <- env$c_cells; starts_k <- env$starts_k; ends_k <- env$ends_k
  rate_g <- env$rate_g
  se <- cross_pair_nuisance_se_run(s, env$ctrl$owner_control, cross_pair_nuisance_general_live(env, D_up),
    env$compute_if, env$N, p, length(c_cells) + 1, N, "nuisance capture treatment IF", function() {
      parts <- list(list(rows = rows, M = S))
      gpos <- tf$g_pos
      if (identical(sel_mode, "ipw")) {
        A <- matrix(0, p, D_up)
        A <- cross_pair_nuisance_ld_add(A, S[gpos, , drop = FALSE], env$ch_g, env$ld_g, N)
        for (kk in seq_along(c_cells)) {
          ix <- tf$c_pos[seq.int(starts_k[kk], ends_k[kk])]
          A <- cross_pair_nuisance_ld_add(A, S[ix, , drop = FALSE], c_cells[[kk]]$ch, c_cells[[kk]]$ld, N)
        }
        parts[[2L]] <- cross_pair_nuisance_up_part(A, up_if)
      } else {
        Ar <- list()
        if (!is.null(rate_g) && isTRUE(rate_g$r_act))
          Ar <- cross_pair_nuisance_rate_add(Ar, rate_g$rate_id, -colSums(S[gpos, , drop = FALSE]) / (N * rate_g$r_c))
        for (kk in seq_along(c_cells)) {
          cc <- c_cells[[kk]]
          if (!isTRUE(cc$r_act)) next
          ix <- tf$c_pos[seq.int(starts_k[kk], ends_k[kk])]
          Ar <- cross_pair_nuisance_rate_add(Ar, cc$rate_id, -colSums(S[ix, , drop = FALSE]) / (N * cc$r_c))
        }
        parts <- c(parts, cross_pair_nuisance_rate_parts(Ar, env$rate_list, N))
      }
      dd <- cross_pair_nuisance_dense(parts, p)
      dd$M <- t(cp04_solve(tf$H_chol, t(dd$M)))
      dd
    })
  cn <- cross_pair_nuisance_terms(X)
  lfit <- get0("lfit", envir = env, inherits = FALSE)
  n_g <- length(tf$g_pos)
  rec <- list(role = "treatment", key = key, links = ell,
    estimate = beta, term = cn, coef_role = ifelse(cn == "(Intercept)", "intercept", "slope"),
    cell = rep(NA_integer_, p), kept = cn, aliased = character(0), intercept = any(cn == "(Intercept)"),
    se = se, check = cross_pair_nuisance_score_check((tf$w_all * (tf$y_all - pb)) * X), check_kind = "logit",
    diagnostics = if (is.list(lfit) && is.list(lfit$diagnostics)) lfit$diagnostics[intersect(c("converged", "iter", "rank"),
      names(lfit$diagnostics))] else NULL,
    meta = list(kind = "weighted_trt", sel_mode = sel_mode, g = d$g, K = d$K, s = d$s, t = d$t,
      pair_tag = env$pair_tag, target_on_pair = isTRUE(d$target_on_pair),
      n_target = n_g, n_control = length(rows) - n_g, w_target = sum(tf$w_all[tf$g_pos]),
      w_control = sum(tf$w_all[tf$c_pos]), n_rows = length(rows), weight_sum = sum(tf$w_all),
      rate_active = identical(sel_mode, "none") && any(vapply(c_cells, function(cc) isTRUE(cc$r_act), TRUE))))
  s$treatment[[key]] <- cross_pair_nuisance_charge(s, rec)
  invisible(NULL)
}

# OR group pieces shared by the weighted and reconstruction cores: one entry
# per (block, support) with rows, weights, residual at the reported
# coefficients, centered design d~ and the alpha score N fw e / W_b.
cross_pair_nuisance_or_piece <- function(bf, rows, fw, y, Zk, intercept, N) {
  p <- length(bf$beta)
  d <- if (p) bf$L * Zk else matrix(0, length(rows), 0L)
  m <- (if (intercept) bf$alpha else 0) + (if (p) as.vector(d %*% bf$beta) else 0)
  e <- as.double(y - m)
  dt <- if (intercept && p) sweep(d, 2L, as.double(bf$mean_d), "-") else d
  list(rows = as.integer(rows), fw = as.double(fw), e = e, dt = dt,
       s_beta = (fw * e) * dt, a_alpha = if (intercept) N * fw * e / bf$weight_sum else NULL)
}

# Assemble beta and alpha IF of one OR group from pieces and channel loads.
#   chans_beta(piece) / chans_alpha(piece, b) are applied by the callers through
#   the A matrices they accumulate; here only the final assembly is done.
cross_pair_nuisance_or_assemble <- function(beta_parts, alpha_parts, Minv, mean_d, p, intercept) {
  xb <- NULL
  if (p) {
    xb <- cross_pair_nuisance_dense(beta_parts, p)
    xb$M <- xb$M %*% Minv
  }
  xa <- list()
  if (intercept) for (b in names(alpha_parts)) {
    parts <- alpha_parts[[b]]
    if (p) parts[[length(parts) + 1L]] <- list(rows = xb$rows, M = -as.vector(xb$M %*% as.double(mean_d[[b]])))
    xa[[b]] <- cross_pair_nuisance_dense(parts, 1L)
  }
  # One dense block over the union rows: alpha columns (block order) then beta.
  all_parts <- list()
  nb <- length(xa)
  W <- nb + p
  for (j in seq_len(nb)) {
    M <- matrix(0, length(xa[[j]]$rows), W); M[, j] <- xa[[j]]$M[, 1L]
    all_parts[[j]] <- list(rows = xa[[j]]$rows, M = M)
  }
  if (p) {
    M <- matrix(0, length(xb$rows), W); M[, nb + seq_len(p)] <- xb$M
    all_parts[[nb + 1L]] <- list(rows = xb$rows, M = M)
  }
  cross_pair_nuisance_dense(all_parts, W)
}

cross_pair_nuisance_weighted_collect <- function(s, env) {
  s$frame <- "F"
  sel_mode <- env$sel_mode
  s$route <- list(core = "weighted", sel_mode = sel_mode, est_method = env$est_method,
    outcome_pooling = env$outcome_pooling, reference = !is.null(env$reference),
    compute_if = isTRUE(env$compute_if), n_links = as.integer(env$L), N_if = as.double(env$N),
    n_rate_cells = length(env$rate_list), n_treatment_fits = as.integer(env$n_trt_fits),
    n_or_blocks = as.integer(env$n_or_blks))
  if (!is.null(env$reference)) return(invisible(NULL))
  if (identical(sel_mode, "ipw")) s$selection <- list(up_map = env$up_map, core = "weighted")
  if (!isTRUE(env$need_or) || !env$n_or_blks) return(invisible(NULL))
  N <- as.double(env$N); rc <- env$recipe_or; icpt <- isTRUE(rc$intercept)
  blk_fits <- env$blk_fits; sup_list <- env$sup_list
  up_if <- get0("up_if", envir = env, inherits = FALSE)
  D_up <- if (is.null(up_if)) 0 else ncol(up_if)
  gids <- vapply(blk_fits, function(bf) as.integer(bf$group_id), integer(1))
  link_blk <- env$link_blk_id
  live <- cross_pair_nuisance_general_live(env, D_up)
  for (gid in sort(unique(gids))) {
    bids <- which(gids == gid)
    bf1 <- blk_fits[[bids[1L]]]
    p <- length(bf1$beta)
    kept_sl <- rc$slope_index[bf1$keep]
    n_rows <- sum(vapply(bids, function(b) sum(vapply(blk_fits[[b]]$support_ids, function(sid)
      length(sup_list[[sid]]$cell$rows_obs), numeric(1))), numeric(1)))
    check <- NA_real_
    n_bound <- n_rows + N * (length(bids) + 1)
    se <- cross_pair_nuisance_se_run(s, env$ctrl$owner_control, live, env$compute_if,
      n_bound, p + length(bids), D_up + 2, N, "nuisance capture outcome-regression IF", function() {
        pieces <- list()
        for (b in bids) for (sid in blk_fits[[b]]$support_ids) {
          sup <- sup_list[[sid]]
          cc <- sup$cell
          pieces[[length(pieces) + 1L]] <- c(list(b = b, cell = cc),
            cross_pair_nuisance_or_piece(blk_fits[[b]], cc$rows_obs, cc$fw, cc$dy,
              sup$Z_obs[, kept_sl, drop = FALSE], icpt, N))
        }
        by_block <- match(vapply(pieces, `[[`, integer(1), "b"), bids)
        check <<- cross_pair_nuisance_pooled_check(lapply(pieces, `[[`, "s_beta"), by_block,
          lapply(pieces, function(pc) if (icpt) pc$fw * pc$e), length(bids) * icpt, p)
        Minv <- if (p) cross_pair_gram_solve(env$or_groups$groups[[gid]]$fit, icpt, diag(p), N) else NULL
        beta_parts <- lapply(pieces, function(pc) list(rows = pc$rows, M = pc$s_beta))
        alpha_parts <- stats::setNames(lapply(bids, function(b)
          lapply(Filter(function(pc) pc$b == b, pieces), function(pc) list(rows = pc$rows, M = pc$a_alpha))),
          as.character(bids))
        if (identical(sel_mode, "ipw")) {
          Ab <- matrix(0, p, D_up); Aa <- lapply(bids, function(b) matrix(0, 1L, D_up)); names(Aa) <- bids
          for (pc in pieces) if (length(pc$cell$ch)) {
            Ab <- cross_pair_nuisance_ld_add(Ab, pc$s_beta, pc$cell$ch, pc$cell$ld, N)
            if (icpt) Aa[[as.character(pc$b)]] <- cross_pair_nuisance_ld_add(Aa[[as.character(pc$b)]],
              matrix(pc$a_alpha, ncol = 1L), pc$cell$ch, pc$cell$ld, N)
          }
          if (p) beta_parts[[length(beta_parts) + 1L]] <- cross_pair_nuisance_up_part(Ab, up_if)
          if (icpt) for (b in names(alpha_parts))
            alpha_parts[[b]][[length(alpha_parts[[b]]) + 1L]] <- cross_pair_nuisance_up_part(Aa[[b]], up_if)
        } else {
          Arb <- list(); Ara <- stats::setNames(vector("list", length(bids)), bids)
          for (pc in pieces) if (isTRUE(pc$cell$r_act)) {
            if (p) Arb <- cross_pair_nuisance_rate_add(Arb, pc$cell$rate_id, -colSums(pc$s_beta) / (N * pc$cell$r_c))
            if (icpt) Ara[[as.character(pc$b)]] <- cross_pair_nuisance_rate_add(Ara[[as.character(pc$b)]],
              pc$cell$rate_id, -sum(pc$a_alpha) / (N * pc$cell$r_c))
          }
          beta_parts <- c(beta_parts, cross_pair_nuisance_rate_parts(Arb, env$rate_list, N))
          if (icpt) for (b in names(alpha_parts))
            alpha_parts[[b]] <- c(alpha_parts[[b]], cross_pair_nuisance_rate_parts(Ara[[b]], env$rate_list, N))
        }
        md <- stats::setNames(lapply(bids, function(b) blk_fits[[b]]$mean_d), bids)
        cross_pair_nuisance_or_assemble(beta_parts, alpha_parts, Minv, md, p, icpt)
      })
    cn_sl <- rc$columns[kept_sl]
    int_name <- if (icpt) rc$columns[rc$intercept_index] else NA_character_
    blk_key <- vapply(bids, function(b) {
      cc <- sup_list[[blk_fits[[b]]$support_ids[1L]]]$cell
      paste0("k=", chained_nuisance_gtxt(env$st$cohorts[cc$ci]), ",",
        chained_nuisance_pair_txt(list(s = env$all_times[cc$cs], t = env$all_times[cc$ct])))
    }, "")
    alpha <- if (icpt) vapply(bids, function(b) as.double(blk_fits[[b]]$alpha), 0) else numeric(0)
    rec <- list(role = "outcome", key = paste0("O|", gid), links = which(link_blk %in% bids),
      estimate = c(alpha, as.double(bf1$beta)),
      term = c(if (icpt) paste0(int_name, "[", blk_key, "]"), cn_sl),
      coef_role = c(rep("local_intercept", length(alpha)), rep(if (length(bids) > 1L) "shared_slope" else "slope", p)),
      cell = c(if (icpt) seq_along(bids), rep(NA_integer_, p)), kept = cn_sl,
      aliased = setdiff(rc$slope_columns, cn_sl), intercept = icpt,
      se = se, check = check, check_kind = "wls",
      meta = list(kind = "weighted_or", sel_mode = sel_mode, blocks = blk_key, block_ids = bids,
        L = vapply(bids, function(b) as.double(blk_fits[[b]]$L), 0),
        cohorts = unique(vapply(bids, function(b) env$st$cohorts[sup_list[[blk_fits[[b]]$support_ids[1L]]]$cell$ci], 0)),
        n_rows = n_rows, weight_sum = sum(vapply(bids, function(b) as.double(blk_fits[[b]]$weight_sum), 0)),
        pooled = length(bids) > 1L))
    s$outcome[[rec$key]] <- cross_pair_nuisance_charge(s, rec)
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 5. Reconstruction core hook (imputation / augmented; frame rows, N = N_orig)
# ------------------------------------------------------------------------------
cross_pair_nuisance_a_if <- function(st, ci, rows) {
  as.double((st$w[rows] / st$mass_F[ci]) * (as.double(st$a_index[rows] > 0L) - st$a[ci]))
}

cross_pair_nuisance_reconstruction_collect <- function(s, env) {
  s$frame <- "F"
  aug <- isTRUE(env$augmented)
  s$route <- list(core = "reconstruction", augmented = aug, est_method = env$est_method,
    outcome_pooling = env$outcome_pooling, imputation_pooling = env$imputation_pooling,
    reference = !is.null(env$reference), compute_if = isTRUE(env$compute_if), n_links = as.integer(env$L),
    N_if = as.double(env$N_orig), n_cells = as.integer(env$n_cells), n_imputation_groups = as.integer(env$n_imp_grps),
    n_or_blocks = as.integer(env$n_or_blks), n_treatment_fits = length(env$trt_list))
  if (!is.null(env$reference) || !is.null(env$lay)) return(invisible(NULL))
  if (aug) s$selection <- list(up_map = env$up_map, core = "reconstruction")
  st <- env$st; N <- as.double(env$N_orig); rc <- env$recipe; icpt <- isTRUE(rc$intercept)
  f_cols <- env$f_cols; X_F <- env$X_F
  cell_data <- env$cell_data; imp_groups <- env$imp_groups
  Xi_theta <- get0("Xi_theta", envir = env, inherits = FALSE)
  D_th <- if (is.null(Xi_theta)) 0 else ncol(Xi_theta)
  live <- cross_pair_nuisance_general_live(env, D_th)
  ctrl <- env$ctrl$owner_control
  pub <- function(ci) st$cohorts[ci]
  cell_key <- function(cf) paste0("g=", chained_nuisance_gtxt(cf$cohort), ",", chained_nuisance_pair_txt(cf$pair))

  # ---- imputation groups (their IF also feeds the Z_c chain of the OR) ------
  imp_xi <- list()
  imp_entries <- 0
  for (gid in seq_along(imp_groups)) {
    gf <- imp_groups[[gid]]; p <- as.integer(gf$p)
    cids <- which(vapply(cell_data, function(cd) cd$gf$group_id == gid, TRUE))
    n_obs <- sum(vapply(cids, function(cid) length(cell_data[[cid]]$obs_local), 0))
    xi_g <- NULL
    check <- NA_real_
    n_bound <- n_obs + N * (length(cids) + 1)
    se <- cross_pair_nuisance_se_run(s, ctrl, c(live, imp_entries), env$compute_if,
      n_bound, p + length(cids), 1, N, "nuisance capture imputation IF", function() {
        Minv <- if (p) cross_pair_gram_solve(gf$fit, icpt, diag(p), N) else NULL
        residual <- lapply(cids, function(cid) {
          cd <- cell_data[[cid]]
          cd$w_raw[cd$obs_local] * cd$e_imp_obs
        })
        beta_parts <- lapply(seq_along(cids), function(j) {
          cd <- cell_data[[cids[j]]]
          list(rows = cd$rows_obs, M = residual[[j]] * cd$d_imp_tilde)
        })
        alpha_parts <- stats::setNames(lapply(seq_along(cids), function(j) {
          cd <- cell_data[[cids[j]]]
          list(list(rows = cd$rows_obs, M = N * residual[[j]] / cd$cf$weight_sum))
        }), as.character(cids))
        check <<- cross_pair_nuisance_pooled_check(lapply(beta_parts, `[[`, "M"),
          seq_along(cids), if (icpt) residual else vector("list", length(cids)),
          length(cids) * icpt, p)
        md <- stats::setNames(lapply(cids, function(cid) cell_data[[cid]]$cf$mean_d), as.character(cids))
        xi_g <<- cross_pair_nuisance_or_assemble(beta_parts, alpha_parts, Minv, md, p, icpt)
        xi_g
      })
    if (!is.null(xi_g)) {
      imp_xi[[gid]] <- list(cids = cids, p = p, rows = xi_g$rows, M = xi_g$M)
      imp_entries <- imp_entries + length(cids) + length(xi_g$rows) + length(xi_g$M)
    }
    kept <- f_cols[gf$kept_cols_F]
    alpha <- if (icpt) vapply(cids, function(cid) as.double(cell_data[[cid]]$cf$alpha), 0) else numeric(0)
    int_name <- if (icpt) rc$columns[rc$intercept_index] else NA_character_
    use <- which(vapply(seq_len(env$L), function(ell)
      any(c(env$link_tgt_cid[ell], env$link_ctrl_cids[[ell]]) %in% cids), TRUE))
    rec <- list(role = "imputation", key = paste0("I|", gid), links = use,
      estimate = c(alpha, as.double(gf$beta)),
      term = c(if (icpt) paste0(int_name, "[", vapply(cids, function(cid) cell_key(cell_data[[cid]]$cf), ""), "]"), kept),
      coef_role = c(rep("local_intercept", length(alpha)), rep(if (length(cids) > 1L) "shared_slope" else "slope", p)),
      cell = c(if (icpt) cids, rep(NA_integer_, p)), kept = kept,
      aliased = setdiff(f_cols[env$slope_cols_F], kept), intercept = icpt,
      se = se, check = check, check_kind = "wls",
      meta = list(kind = "imputation", pooling = env$imputation_pooling, cids = cids,
        cells = lapply(cids, function(cid) {
          cd <- cell_data[[cid]]
          list(cid = cid, cohort = cd$cf$cohort, s = cd$cf$pair$s, t = cd$cf$pair$t, L = cd$cf$L,
            n_cohort_rows = length(cd$rf), n_observed_pairs = length(cd$obs_local), weight_sum = cd$cf$weight_sum)
        }),
        n_rows = n_obs, weight_sum = sum(vapply(cids, function(cid) as.double(cell_data[[cid]]$cf$weight_sum), 0))))
    s$imputation[[rec$key]] <- cross_pair_nuisance_charge(s, rec)
  }

  # ---- synthetic treatment logits (repeat-rate channel only) ----------------
  if (isTRUE(env$need_trt)) {
    keys <- ls(env$trt_env, all.names = TRUE)
    tid_of <- vapply(keys, function(k) as.integer(env$trt_env[[k]]), 0L)
    for (j in seq_along(keys)) {
      tr <- env$trt_list[[tid_of[j]]]
      X <- tr$X_trt; p <- ncol(X); rows <- as.integer(tr$risk_rows); beta <- as.double(tr$coefficients)
      S <- (tr$w_syn * (tr$y_r - tr$prob)) * X
      pb <- stats::plogis(as.double(X %*% beta))
      se <- cross_pair_nuisance_se_run(s, ctrl, live, env$compute_if, length(rows), p, 1, N,
        "nuisance capture treatment IF", function() {
          parts <- list(list(rows = rows, M = S))
          for (kk in seq_along(tr$all_ci)) {
            ci <- tr$all_ci[kk]; rf <- env$rows_by_ci[[ci]]
            parts[[length(parts) + 1L]] <- cross_pair_nuisance_channel(rf,
              matrix(cross_pair_nuisance_a_if(st, ci, rf), ncol = 1L), matrix(tr$C_a[, kk], ncol = 1L))
          }
          dd <- cross_pair_nuisance_dense(parts, p)
          dd$M <- t(cp04_solve(tr$H_chol, t(dd$M)))
          dd
        })
      use <- which(vapply(seq_len(env$L), function(ell) {
        if (!env$descs[[ell]]$has_control) return(FALSE)
        all_ci <- c(env$link_ci_g[ell], sort(unique(env$link_ci_K[[ell]])))
        identical(paste(env$link_ci_g[ell], paste(all_ci[-1L], collapse = ","), sep = "|"), keys[j])
      }, TRUE))
      cn <- if (length(env$trt_recipe$columns) == p) as.character(env$trt_recipe$columns) else cross_pair_nuisance_terms(X)
      int_name <- if (icpt) rc$columns[rc$intercept_index] else "(Intercept)"
      is_t <- tr$y_r == 1
      rec <- list(role = "treatment", key = paste0("T|", keys[j]), links = use,
        estimate = beta, term = cn, coef_role = ifelse(cn %in% int_name, "intercept", "slope"),
        cell = rep(NA_integer_, p), kept = cn, aliased = character(0), intercept = any(cn %in% int_name),
        se = se, check = cross_pair_nuisance_score_check((tr$w_syn * (tr$y_r - pb)) * X), check_kind = "logit",
        meta = list(kind = "synthetic_trt", g = pub(tr$all_ci[1L]), K = pub(tr$all_ci[-1L]),
          repeat_rates = st$a[tr$all_ci], n_target = sum(is_t), n_control = sum(!is_t),
          w_target = sum(tr$w_syn[is_t]), w_control = sum(tr$w_syn[!is_t]),
          n_rows = length(rows), weight_sum = sum(tr$w_syn)))
      s$treatment[[rec$key]] <- cross_pair_nuisance_charge(s, rec)
    }
  }

  # ---- OR on the reconstructed response Z_c (full generated-response chain) --
  if (!isTRUE(env$need_or) || !env$n_or_blks) return(invisible(NULL))
  blk_fits <- env$blk_fits; sup_or <- env$sup_or
  gids <- vapply(blk_fits, function(bf) as.integer(bf$group_id), integer(1))
  slope_F <- env$slope_cols_F
  imp_width <- sum(vapply(imp_xi, function(x) if (is.null(x)) 0 else ncol(x$M), numeric(1)))
  for (gid in sort(unique(gids))) {
    bids <- which(gids == gid); bf1 <- blk_fits[[bids[1L]]]; p <- length(bf1$beta)
    kept_F <- slope_F[bf1$keep]
    n_rows <- sum(vapply(bids, function(b) sum(vapply(blk_fits[[b]]$support_ids, function(sid)
      length(cell_data[[sup_or[[sid]]$cid]]$rf), numeric(1))), numeric(1)))
    imp_ok <- length(imp_xi) == length(imp_groups) && all(!vapply(imp_xi, is.null, TRUE))
    th_ok <- !aug || D_th == 0 || !is.null(Xi_theta)
    check <- NA_real_
    n_bound <- n_rows + N * (length(bids) + 1)
    se <- if (!imp_ok && isTRUE(env$compute_if)) cross_pair_nuisance_se(code = "upstream_if_unavailable",
      reason = "The imputation coefficient IF entering Z_c was not evaluated, so the generated-response chain is incomplete.") else
      cross_pair_nuisance_se_run(s, ctrl, c(live, imp_entries), env$compute_if && th_ok, n_bound,
        p + length(bids), imp_width + D_th + 2, N,
        "nuisance capture counterfactual-OR IF", function() {
        pieces <- list()
        for (b in bids) for (sid in blk_fits[[b]]$support_ids) {
          cid <- sup_or[[sid]]$cid
          cd <- cell_data[[cid]]
          pc <- cross_pair_nuisance_or_piece(blk_fits[[b]], cd$rf, cd$w_syn, cd$Z,
            X_F[cd$rf, kept_F, drop = FALSE], icpt, N)
          pieces[[length(pieces) + 1L]] <- c(list(b = b, cid = cid, cd = cd), pc)
        }
        by_block <- match(vapply(pieces, `[[`, integer(1), "b"), bids)
        check <<- cross_pair_nuisance_pooled_check(lapply(pieces, `[[`, "s_beta"), by_block,
          lapply(pieces, function(pc) if (icpt) pc$fw * pc$e), length(bids) * icpt, p)
        Minv <- if (p) cross_pair_gram_solve(env$or_groups$groups[[gid]]$fit, icpt, diag(p), N) else NULL
        nb <- length(bids); bkey <- as.character(bids)
        beta_parts <- lapply(pieces, function(pc) list(rows = pc$rows, M = pc$s_beta))
        alpha_parts <- stats::setNames(lapply(bids, function(b)
          lapply(Filter(function(pc) pc$b == b, pieces), function(pc) list(rows = pc$rows, M = pc$a_alpha))), bkey)
        # Loads on the repeat rates (by cohort), the imputation IF columns (by
        # imputation group: cell intercepts then slopes) and theta.
        Aa_b <- list(); Aa_a <- stats::setNames(vector("list", nb), bkey)
        Ai_b <- list(); Ai_a <- stats::setNames(lapply(bids, function(b) list()), bkey)
        Ath_b <- matrix(0, p, D_th); Ath_a <- stats::setNames(lapply(bids, function(b) matrix(0, 1L, D_th)), bkey)
        for (pc in pieces) {
          cd <- pc$cd; ci <- cd$ci; bk <- as.character(pc$b); a_k <- st$a[ci]
          if (p) Aa_b <- cross_pair_nuisance_rate_add(Aa_b, ci, colSums(pc$s_beta) / (N * a_k))
          if (icpt) Aa_a[[bk]] <- cross_pair_nuisance_rate_add(Aa_a[[bk]], ci, sum(pc$a_alpha) / (N * a_k))
          # d s / d j: (1 - S pi) fw d~ (beta) and (1 - S pi) N fw / W_b (alpha).
          ig <- cd$gf$group_id; ix <- imp_xi[[ig]]; nci <- if (icpt) length(ix$cids) else 0L
          om <- cd$one_minus_spi
          hb <- (om * pc$fw / N) * pc$dt
          ha <- if (icpt) om * pc$fw / blk_fits[[pc$b]]$weight_sum else NULL
          grow <- function(A, q) if (is.null(A)) matrix(0, q, nci + ix$p) else A
          if (p) {
            A <- grow(Ai_b[[as.character(ig)]], p)
            if (icpt) { j <- match(pc$cid, ix$cids); A[, j] <- A[, j] + colSums(hb) }
            if (ix$p) A[, nci + seq_len(ix$p)] <- A[, nci + seq_len(ix$p), drop = FALSE] + t(crossprod(cd$d_imp, hb))
            Ai_b[[as.character(ig)]] <- A
          }
          if (icpt) {
            A <- grow(Ai_a[[bk]][[as.character(ig)]], 1L)
            j <- match(pc$cid, ix$cids); A[, j] <- A[, j] + sum(ha)
            if (ix$p) A[, nci + seq_len(ix$p)] <- A[, nci + seq_len(ix$p), drop = FALSE] + as.vector(crossprod(cd$d_imp, ha))
            Ai_a[[bk]][[as.character(ig)]] <- A
          }
          # d s / d theta (augmented): pi e_imp T_logit on observed pairs.
          if (aug && D_th && length(cd$theta_ch) && length(cd$obs_local)) {
            ob <- cd$obs_local; f <- cd$pi_obs * cd$e_imp_obs
            if (p) Ath_b[, cd$theta_ch] <- Ath_b[, cd$theta_ch, drop = FALSE] +
              t(crossprod(cd$T_logit, f * (pc$fw[ob] / N) * pc$dt[ob, , drop = FALSE]))
            if (icpt) Ath_a[[bk]][, cd$theta_ch] <- Ath_a[[bk]][, cd$theta_ch, drop = FALSE] +
              as.vector(crossprod(cd$T_logit, f * pc$fw[ob] / blk_fits[[pc$b]]$weight_sum))
          }
        }
        a_parts <- function(Aa) lapply(names(Aa), function(ci) {
          rf <- env$rows_by_ci[[as.integer(ci)]]
          cross_pair_nuisance_channel(rf, matrix(cross_pair_nuisance_a_if(st, as.integer(ci), rf), ncol = 1L),
            matrix(Aa[[ci]], ncol = 1L))
        })
        i_parts <- function(Ai) lapply(names(Ai), function(ig) {
          ix <- imp_xi[[as.integer(ig)]]
          cross_pair_nuisance_channel(ix$rows, ix$M, Ai[[ig]])
        })
        if (p) beta_parts <- c(beta_parts, a_parts(Aa_b), i_parts(Ai_b),
          list(if (D_th) cross_pair_nuisance_up_part(Ath_b, Xi_theta)))
        if (icpt) for (bk in bkey) alpha_parts[[bk]] <- c(alpha_parts[[bk]], a_parts(Aa_a[[bk]]),
          i_parts(Ai_a[[bk]]), list(if (D_th) cross_pair_nuisance_up_part(Ath_a[[bk]], Xi_theta)))
        md <- stats::setNames(lapply(bids, function(b) blk_fits[[b]]$mean_d), bkey)
        cross_pair_nuisance_or_assemble(beta_parts, alpha_parts, Minv, md, p, icpt)
      })
    cn_sl <- f_cols[kept_F]
    int_name <- if (icpt) rc$columns[rc$intercept_index] else NA_character_
    blk_key <- vapply(bids, function(b) {
      cd <- cell_data[[sup_or[[blk_fits[[b]]$support_ids[1L]]]$cid]]
      paste0("k=", chained_nuisance_gtxt(vapply(blk_fits[[b]]$support_ids, function(sid)
        pub(cell_data[[sup_or[[sid]]$cid]]$ci), 0)), ",", chained_nuisance_pair_txt(cd$cf$pair))
    }, "")
    alpha <- if (icpt) vapply(bids, function(b) as.double(blk_fits[[b]]$alpha), 0) else numeric(0)
    rec <- list(role = "outcome", key = paste0("O|", gid), links = which(env$link_or_bid %in% bids),
      estimate = c(alpha, as.double(bf1$beta)),
      term = c(if (icpt) paste0(int_name, "[", blk_key, "]"), cn_sl),
      coef_role = c(rep("local_intercept", length(alpha)), rep(if (length(bids) > 1L) "shared_slope" else "slope", p)),
      cell = c(if (icpt) seq_along(bids), rep(NA_integer_, p)), kept = cn_sl,
      aliased = setdiff(f_cols[slope_F], cn_sl), intercept = icpt,
      se = se, check = check, check_kind = "wls",
      meta = list(kind = "reconstruction_or", augmented = aug, blocks = blk_key, block_ids = bids,
        L = vapply(bids, function(b) as.double(blk_fits[[b]]$L), 0),
        cohorts = unique(unlist(lapply(bids, function(b) vapply(blk_fits[[b]]$support_ids, function(sid)
          pub(cell_data[[sup_or[[sid]]$cid]]$ci), 0)))),
        n_rows = n_rows, weight_sum = sum(vapply(bids, function(b) as.double(blk_fits[[b]]$weight_sum), 0)),
        pooled = length(bids) > 1L))
    s$outcome[[rec$key]] <- cross_pair_nuisance_charge(s, rec)
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# 6. Capture input (plain list) returned on res$nuisance_capture_input
# ------------------------------------------------------------------------------
cross_pair_nuisance_input <- function(s, status = "fitted", reason = NULL) {
  strip <- function(recs) lapply(recs, function(r) { r$estimate_piv <- NULL; r })
  list(schema = "cdid.cross_pair_nuisance_input.v0", status = status, reason = reason,
    frame = s$frame, route = s$route, clustered = s$clustered, n_groups = s$n_groups,
    cluster_code = s$cluster_code, selection = s$selection,
    treatment = strip(s$treatment), outcome = strip(s$outcome), imputation = strip(s$imputation),
    accounting = list(capture_work_ops = s$work, capture_output_entries_bound = s$entries))
}

# ------------------------------------------------------------------------------
# 7. Components and the capture (schema of chained_nuisance_capture_assemble)
# ------------------------------------------------------------------------------
cross_pair_nuisance_link_table <- function(res, inp) {
  dl <- res$design_all$links
  M <- if (is.data.frame(dl)) nrow(dl) else 0L
  kept <- as.integer(res$reduction_summary$kept_indices)
  if (!length(kept)) kept <- seq_len(M)
  list(M = M, kept = kept,
    catalog = function(ell) if (identical(inp$frame, "A")) kept[ell] else as.integer(ell),
    id = function(cat) if (M && all(cat <= M) && !is.null(dl$link_id)) as.character(dl$link_id[cat]) else paste0("cp_link_", cat),
    row = function(cat) dl[cat, intersect(c("g", "s", "t", "link_role", "cohort"), names(dl)), drop = FALSE])
}

cross_pair_nuisance_status_all <- function(st) list(outcome = st, treatment = st, imputation = st)

cross_pair_nuisance_cov_text <- list(
  treatment = list(
    baseline = "the estimated observation rate r of the control cohort's mask in the control weights w_i / r",
    none = "the estimated observation rates r in the weights w_i / r of rows on a sub-sampled pair",
    ipw = "the estimated repeat rates a_k and observation-selection logits in the weights w_i a_k / (rho_i max(p_i, eps))",
    reconstruction = "the estimated repeat rates a_k in its weights w_i a_k"),
  outcome = list(
    baseline = "nothing: its weights w_i are fixed and its response is the observed change",
    none = "the estimated observation rates r in the weights w_i / r",
    ipw = "the estimated repeat rates a_k and observation-selection logits in the weights w_i a_k / (rho_i max(p_i, eps))",
    reconstruction = "the estimated repeat rates a_k in its weights, and the imputation fits j (and, augmented, the observation logits q) that build Z_c"),
  imputation = list(reconstruction = "nothing: its weights are fixed and its response is observed"))

cross_pair_nuisance_covariance <- function(role, route_key, recs, inp, cluster_var, N) {
  has <- vapply(recs, function(r) !is.null(r$se$var), TRUE)
  codes <- unique(vapply(recs[!has], function(r) if (is.null(r$se$code)) "unavailable" else r$se$code, ""))
  reasons <- unique(vapply(recs[!has], function(r) if (is.null(r$se$reason)) "" else r$se$reason, ""))
  clustered <- isTRUE(inp$clustered)
  status <- if (all(has)) "available" else if (any(has)) "partial" else "unavailable"
  list(status = status, code = if (all(has)) "ok" else paste(codes, collapse = ","),
    reason = if (all(has)) NULL else paste(c(reasons, "Coefficients are reported; SEs are not for the listed fits."),
      collapse = " "),
    type = if (clustered) "cluster" else "unit", cluster_var = cluster_var,
    n_groups = if (clustered) inp$n_groups else N, N = N,
    estimator = paste0("Centered sandwich: each original unit's coefficient IF contribution",
      if (clustered) ", summed within its cluster" else "", "; squared centered sums / N^2 over all ",
      if (identical(inp$frame, "A")) "pair-eligible (A) units" else "frame units",
      " (zero-contribution units included), no small-sample factor. Diagonal only."),
    source = paste("The estimator's own reverse-mode IF operator of this fit evaluated at unit sensitivities:",
      "the per-unit score through the fit's stored factor plus every upstream channel it propagates."),
    accounts_for = paste("Accounts for", cross_pair_nuisance_cov_text[[role]][[route_key]]),
    diagonal_only = TRUE)
}

cross_pair_nuisance_route_key <- function(inp) {
  r <- inp$route
  if (identical(r$core, "baseline")) "baseline" else if (identical(r$core, "reconstruction")) "reconstruction" else
    if (identical(r$sel_mode, "ipw")) "ipw" else "none"
}

cross_pair_nuisance_label <- function(role, rec) {
  m <- rec$meta
  switch(m$kind,
    baseline_trt = paste0("risk_set[g=", chained_nuisance_gtxt(m$g), " vs k=", chained_nuisance_gtxt(m$cohort),
      " observed on mask ", m$mask, ": pairs ", paste(m$pairs, collapse = ","), "]"),
    weighted_trt = paste0("risk_set[g=", chained_nuisance_gtxt(m$g), " vs K={", chained_nuisance_gtxt(m$K), "}",
      if (identical(m$pair_tag, "complete_A")) ", complete A rows" else
        paste0(", observed ", chained_nuisance_pair_txt(list(s = m$s, t = m$t))), "]"),
    synthetic_trt = paste0("risk_set[g=", chained_nuisance_gtxt(m$g), " vs K={", chained_nuisance_gtxt(m$K), "}]"),
    baseline_or = paste0("pair[k=", chained_nuisance_gtxt(m$cohort), ",",
      chained_nuisance_pair_txt(list(s = m$s, t = m$t)), "]"),
    weighted_or = , reconstruction_or = if (isTRUE(m$pooled))
      paste0("pooled[", length(m$blocks), " control block(s), common slope]") else paste0("block[", m$blocks, "]"),
    imputation = paste0(m$pooling, "[", paste(vapply(m$cells, function(c) paste0("g=", chained_nuisance_gtxt(c$cohort), ",",
      chained_nuisance_pair_txt(list(s = c$s, t = c$t))), ""), collapse = "; "), "]"))
}

cross_pair_nuisance_equation <- function(role, rec, inp) {
  m <- rec$meta
  if (role == "treatment") {
    # The exact fitted risk set of this fit (target side, control side).
    g <- chained_nuisance_gtxt(m$g)
    K <- chained_nuisance_gtxt(if (!is.null(m$K)) m$K else m$cohort)
    k_set <- if (grepl(",", K, fixed = TRUE)) paste0("G_i in {", K, "}") else paste0("G_i = ", K)
    pr <- function(s, t) paste0("S_i(", chained_nuisance_num(s), ",", chained_nuisance_num(t), ") = 1")
    risk <- switch(m$kind,
      baseline_trt = paste0("G_i = ", g, " or (", k_set, ", ", sub("^\\((.*)\\)$", "S_i(\\1) = 1", m$pairs[1L]), ")"),
      weighted_trt = paste0("(G_i = ", g, if (identical(m$sel_mode, "ipw")) ", U_i = 1)" else ", A_i = 1)",
        " or (", k_set, if (identical(m$pair_tag, "complete_A")) ", A_i = 1)" else paste0(", ", pr(m$s, m$t), ")")),
      paste0("G_i in {", g, ",", K, "}"))
    return(paste0("logit Pr(G_i = ", g, " | X_i, ", risk, ") = ",
      cross_pair_nuisance_logit_rhs(rec$term, rec$term[rec$coef_role == "intercept"])))
  }
  a_sym <- function(single) if (!isTRUE(rec$intercept)) NULL else if (single) "a" else "a[b]"
  if (identical(m$kind, "baseline_or"))
    return(paste(cross_pair_nuisance_delta_lhs(m$cohort, m$s, m$t), "=",
      cross_pair_nuisance_rhs(rec$kept, if (isTRUE(rec$intercept)) "a")))
  if (identical(m$kind, "imputation")) {
    single <- length(m$cells) == 1L
    Ls <- unique(vapply(m$cells, function(c) as.double(c$L), 0))
    lhs <- if (single) cross_pair_nuisance_delta_lhs(m$cells[[1L]]$cohort, m$cells[[1L]]$s, m$cells[[1L]]$t) else
      cross_pair_nuisance_delta_lhs(vapply(m$cells, function(c) c$cohort, 0), NULL, NULL)
    return(paste(lhs, "=", cross_pair_nuisance_rhs(rec$kept, if (isTRUE(rec$intercept)) (if (single) "a" else "a[c]"),
      if (length(Ls) == 1L) Ls else "L(s,t)")))
  }
  single <- length(m$blocks) == 1L
  Ls <- unique(m$L)
  L <- if (length(Ls) == 1L) Ls else "L(s,t)"
  rhs <- cross_pair_nuisance_rhs(rec$kept, a_sym(single), L)
  if (identical(m$kind, "reconstruction_or"))
    return(paste0("E[Z_i(k;s,t) | X_i, G_i = k] = ", rhs, ",  ",
      if (single) paste0("block ", m$blocks) else "k in K_b for each control block b = ((s,t), K_b)"))
  paste0("E[Y_it - Y_is | X_i, G_i = k, S_i(s,t) = 1] = ", rhs, ",  ",
    if (single) paste0("block ", m$blocks) else "for each control block b = (k,(s,t))")
}

cross_pair_nuisance_texts <- function(role, route_key, inp, weights_text) {
  aug <- isTRUE(inp$route$augmented)
  rate_w <- paste0(weights_text, " / r on rows of a sub-sampled pair (r: observed share of the cohort's weight on the pair)")
  ipw_w <- paste0(weights_text, " a_k / (rho_i max(p_i(s,t), eps)) on observed pairs")
  switch(role,
    treatment = list(
      role = "treatment propensity logit of the cross-pair control links",
      template = switch(route_key, reconstruction = "logit Pr(G_i = g | X_i, G_i in {g} U K) = X_i'c[g,K]",
        "logit Pr(G_i = g | X_i, i in R) = X_i'c[R],  R = target cohort g rows + observed control-cohort pair rows"),
      response = "D_i = 1{G_i = g} on the fit's risk set",
      sample = switch(route_key,
        baseline = "target cohort g pair-eligible units and the control cohort's units observed on the link's mask",
        reconstruction = "every frame unit of the target and control cohorts",
        "target cohort g rows (off-pair domain) and the control cohort's rows observed at both pair dates"),
      weights = switch(route_key,
        baseline = list(text = paste0("target ", weights_text, "; controls ", rate_w), estimated = TRUE,
          note = "The control rate r is estimated; its IF is carried in the SE when it is active (r < 1)."),
        none = list(text = paste0(weights_text, ", or ", rate_w), estimated = TRUE,
          note = "Observation rates r on sub-sampled pairs are estimated and carried in the SE."),
        ipw = list(text = paste0("target ", weights_text, " a_g / rho_i; controls ", ipw_w), estimated = TRUE,
          note = "Repeat rates and the selection probabilities are estimated and carried in the SE."),
        reconstruction = list(text = paste0(weights_text, " times a_{G_i}"), estimated = TRUE,
          note = "w_i a_{G_i}: the estimated repeat rate of the unit's own cohort is carried in the SE.")),
      use = "the fitted odds weight the control rows of each listed cross-pair control link (control mean; DR: residual correction)"),
    outcome = list(
      role = if (route_key == "reconstruction") "counterfactual control outcome regression on the reconstructed response Z_c" else
        "control outcome-change regression on the observed actual change",
      template = if (route_key == "reconstruction") "E[Z_i(k;s,t) | X_i, G_i = k] = a[b] + L(s,t)*X_i'b" else
        if (route_key == "baseline") "E[Y_it - Y_is | X_i, G_i = k, S_i(s,t) = 1] = X_i'b[k,mask,(s,t)]" else
        "E[Y_it - Y_is | X_i, G_i = k, S_i(s,t) = 1] = a[b] + L(s,t)*X_i'b",
      response = if (route_key == "reconstruction") (if (aug)
        "Z_c = j_c + S_i(s,t) (Y_it - Y_is - j_c) / max(q_c, eps) on observed pairs, j_c otherwise (generated response)" else
        "Z_c = j_c, the fitted imputation, on every frame row of the control cell (generated response)") else
        "Delta_i(s,t) = Y_it - Y_is on control units observed at both dates (actual change)",
      sample = if (route_key == "reconstruction") "every frame unit of each control cell's cohort; rows are unit x cell" else
        "control-cohort units observed at both pair dates; rows are unit x pair",
      weights = switch(route_key,
        baseline = list(text = weights_text, estimated = FALSE, note = "fixed known weights; no rate in this fit."),
        none = list(text = paste0(weights_text, ", or ", rate_w), estimated = TRUE, note = "Observation rates are carried in the SE."),
        ipw = list(text = ipw_w, estimated = TRUE, note = "Repeat rates and selection probabilities are carried in the SE."),
        reconstruction = list(text = paste0(weights_text, " times a_k"), estimated = TRUE,
          note = "w_i a_k: the estimated repeat rate of the control cohort is carried in the SE.")),
      use = "the fitted m(X) is averaged over target cohort g in each listed control link (DR: minus its odds-weighted control mean)"),
    imputation = list(
      role = "actual-change imputation regression j_c",
      template = "E[Y_it - Y_is | X_i, G_i = g, S_i(s,t) = 1] = a[c] + L(s,t)*X_i'b[r]",
      response = "Delta_i(s,t) = Y_it - Y_is on observed pairs S_i(s,t) = 1",
      sample = "units of the cell's cohort observed at both pair dates; rows are unit x cell",
      weights = list(text = weights_text, estimated = FALSE, note = "fixed known weights on observed pairs only."),
      use = if (aug) "j_c is the predictor in Z_c = j_c + S / max(q_c, eps) (Delta - j_c) of every frame row of cell c" else
        "j_c is the reconstructed change Z_c of every frame row of cell c"))
}

cross_pair_nuisance_role_component <- function(role, recs, inp, lt, weights_text, cluster_var, formula_text) {
  if (!length(recs)) return(NULL)
  route_key <- cross_pair_nuisance_route_key(inp)
  N <- inp$route$N_if
  prefix <- c(treatment = "T", outcome = "O", imputation = "I")[[role]]
  ids <- paste0(prefix, seq_along(recs))
  coefs <- do.call(rbind, c(list(chained_nuisance_coef_frame(character(0), character(0), character(0), integer(0),
    numeric(0), numeric(0))), lapply(seq_along(recs), function(m) {
      r <- recs[[m]]
      v <- if (is.null(r$se$var)) rep(NA_real_, length(r$estimate)) else r$se$var
      chained_nuisance_coef_frame(ids[m], r$term, r$coef_role, r$cell, r$estimate, v)
    })))
  rownames(coefs) <- NULL
  used_links <- lapply(recs, function(r) intersect(lt$catalog(r$links), lt$kept))
  tol <- cross_pair_nuisance_tol
  blocks <- data.frame(block = ids, label = vapply(recs, function(r) cross_pair_nuisance_label(role, r), ""),
    cohorts = vapply(recs, function(r) {
      m <- r$meta
      chained_nuisance_gtxt(sort(unique(c(m$g, m$K, m$cohort, m$cohorts,
        if (!is.null(m$cells)) vapply(m$cells, function(c) c$cohort, 0)))))
    }, ""),
    n_parameters = vapply(recs, function(r) length(r$estimate), integer(1)),
    n_intercepts = vapply(recs, function(r) sum(r$coef_role %in% c("intercept", "local_intercept")), integer(1)),
    n_slopes = vapply(recs, function(r) sum(r$coef_role %in% c("slope", "shared_slope")), integer(1)),
    kept_slopes = vapply(recs, function(r) paste(setdiff(r$kept, r$term[r$coef_role == "intercept"]), collapse = ","), ""),
    aliased_slopes = vapply(recs, function(r) paste(r$aliased, collapse = ","), ""),
    n_fitting_rows = vapply(recs, function(r) as.double(r$meta$n_rows), 0),
    weight_sum = vapply(recs, function(r) as.double(r$meta$weight_sum), 0),
    n_contributing_clusters = vapply(recs, function(r) as.integer(r$se$clusters), integer(1)),
    n_comparisons = lengths(used_links),
    converged = vapply(recs, function(r) if (is.null(r$diagnostics$converged)) TRUE else isTRUE(r$diagnostics$converged), TRUE),
    iterations = vapply(recs, function(r) chained_nuisance_num1(r$diagnostics$iter), 0),
    score_check = vapply(recs, function(r) as.double(r$check), 0),
    verified = vapply(recs, function(r) is.finite(r$check) && r$check <= tol[[r$check_kind]], TRUE),
    se_status = vapply(recs, function(r) if (is.null(r$se$var)) r$se$code else "ok", ""),
    equation = vapply(recs, function(r) cross_pair_nuisance_equation(role, r, inp), ""),
    stringsAsFactors = FALSE)
  cats <- sort(unique(unlist(used_links)))
  fit_of <- lapply(cats, function(cc) ids[vapply(used_links, function(u) cc %in% u, TRUE)])
  lrow <- if (length(cats)) lt$row(cats) else NULL
  comparisons <- data.frame(comparison = if (length(cats)) lt$id(cats) else character(0),
    g = if (length(cats)) chained_nuisance_g(as.double(lrow$g)) else numeric(0),
    s = if (length(cats)) as.double(lrow$s) else numeric(0),
    t = if (length(cats)) as.double(lrow$t) else numeric(0),
    link_role = if (length(cats)) as.character(lrow$link_role) else character(0),
    controls = if (length(cats)) ifelse(lrow$link_role == "target", "", chained_nuisance_num(chained_nuisance_g(as.double(lrow$cohort)))) else character(0),
    fits = vapply(fit_of, function(x) paste(x, collapse = ","), ""),
    stringsAsFactors = FALSE)
  detail <- switch(role,
    imputation = do.call(rbind, c(list(data.frame(cell = integer(0), block = character(0), cohort = numeric(0),
      s = numeric(0), t = numeric(0), L = numeric(0), n_cohort_rows = numeric(0), n_observed_pairs = numeric(0),
      weight_sum = numeric(0))), lapply(seq_along(recs), function(m) do.call(rbind, lapply(recs[[m]]$meta$cells, function(c)
      data.frame(cell = as.integer(c$cid), block = ids[m], cohort = chained_nuisance_g(c$cohort), s = as.double(c$s),
        t = as.double(c$t), L = as.double(c$L), n_cohort_rows = as.double(c$n_cohort_rows),
        n_observed_pairs = as.double(c$n_observed_pairs), weight_sum = as.double(c$weight_sum))))))),
    outcome = data.frame(block = ids, control_blocks = vapply(recs, function(r)
      if (is.null(r$meta$blocks)) cross_pair_nuisance_label(role, r) else paste(r$meta$blocks, collapse = "; "), ""),
      L = vapply(recs, function(r) paste(chained_nuisance_num(unique(r$meta$L)), collapse = ","), ""),
      n_fitting_rows = blocks$n_fitting_rows, weight_sum = blocks$weight_sum,
      reuse_max_diff = vapply(recs, function(r) if (is.null(r$reuse_max_diff)) NA_real_ else as.double(r$reuse_max_diff), 0),
      n_comparisons = blocks$n_comparisons, stringsAsFactors = FALSE),
    treatment = data.frame(block = ids,
      g = vapply(recs, function(r) chained_nuisance_g(as.double(r$meta$g)), 0),
      controls = vapply(recs, function(r) chained_nuisance_gtxt(if (!is.null(r$meta$K)) r$meta$K else r$meta$cohort), ""),
      n_target_rows = vapply(recs, function(r) as.double(r$meta$n_target), 0),
      n_control_rows = vapply(recs, function(r) as.double(r$meta$n_control), 0),
      target_weight_sum = vapply(recs, function(r) as.double(r$meta$w_target), 0),
      control_weight_sum = vapply(recs, function(r) as.double(r$meta$w_control), 0),
      rate_channel = vapply(recs, function(r) isTRUE(r$meta$rate_active), TRUE),
      n_comparisons = blocks$n_comparisons, stringsAsFactors = FALSE))
  if (is.data.frame(detail)) rownames(detail) <- NULL
  txt <- cross_pair_nuisance_texts(role, route_key, inp, weights_text)
  out <- list(status = "available", role = txt$role, route = paste0("cross-pair ", inp$route$core, " route",
      if (!is.null(inp$route$sel_mode)) paste0(" (selection \"", inp$route$sel_mode, "\")") else ""),
    template = txt$template, response = txt$response, use = txt$use, formula = formula_text,
    design_columns = unique(unlist(lapply(recs, function(r) c(r$kept, r$aliased)))),
    intercept_term = NULL, slope_terms = unique(unlist(lapply(recs, `[[`, "kept"))),
    pooling = switch(role, imputation = if (is.null(inp$route$imputation_pooling)) "separate" else inp$route$imputation_pooling,
      outcome = if (identical(inp$route$core, "baseline")) "separate" else inp$route$outcome_pooling, treatment = "risk_set"),
    sample = txt$sample, weights = txt$weights, n_frame = N,
    blocks = blocks, coefficients = coefs, comparisons = comparisons,
    covariance = cross_pair_nuisance_covariance(role, route_key, recs, inp, cluster_var, N),
    accounting = inp$accounting)
  out[[switch(role, imputation = "cells", outcome = "control_blocks", treatment = "models")]] <- detail
  out
}

# Records whose fit is used by at least one retained (kept) cross-pair link.
cross_pair_nuisance_used <- function(recs, lt) Filter(function(r) length(intersect(lt$catalog(r$links), lt$kept)) > 0L, recs)

cross_pair_nuisance_selection <- function(res, inp, prep, weightsname, lt) {
  b <- res$nuisance_registry$cp04
  r <- inp$route
  if (identical(r$core, "baseline") || is.null(b) || !isTRUE(b$uses_selection_models))
    return(chained_nuisance_status("not_applicable", "no_selection_model", paste0(
      "This cross-pair route fits no observation-selection logit (",
      if (identical(r$core, "baseline")) "baseline cross-pair moments" else
        paste0("selection mode \"", if (is.null(r$sel_mode)) (if (isTRUE(r$augmented)) "dr" else "imputation") else r$sel_mode, "\""), ").")))
  if (isTRUE(r$reference) || !is.null(b$reference_ref))
    return(chained_nuisance_status("unavailable", "reference_recipe_route",
      "Dated covariates or recipe-indexed weights fit the selection logits on the recipe layout; outside this capture."))
  sb <- b$selection_backend; up_map <- inp$selection$up_map
  if (!inherits(sb, "cdid_selection_models") || !is.list(up_map) ||
      length(up_map$theta_coords) != length(sb$blocks) || !is.environment(sb$reader_ref))
    return(chained_nuisance_status("unavailable", "unrecognized_backend",
      "The cross-pair backend does not hold the shared selection backend with its captured upstream map."))
  validate_selection_models(sb)
  reader <- get("cp04_reader", sb$reader_ref, inherits = FALSE)
  # Used blocks: those read by a pair-cohort of a retained link.
  kept <- lt$kept[lt$kept <= length(sb$links)]
  used_q <- sort(unique(unlist(lapply(sb$links[kept], function(l) as.integer(l$pair_cohort_ids)))))
  used_b <- sort(unique(unlist(lapply(sb$pair_cohorts[used_q], function(pc) as.integer(pc$block_ids)))))
  cells_k <- sort(unique(unlist(lapply(sb$blocks[used_b], function(bl) as.integer(bl$cell_ids)))))
  re <- function(x, map) match(as.integer(x), map)
  cells <- lapply(sb$cells[cells_k], function(cl) { cl$block_id <- re(cl$block_id, used_b); cl })
  blocks <- lapply(sb$blocks[used_b], function(bl) { bl$cell_ids <- re(bl$cell_ids, cells_k); bl })
  pcs <- lapply(sb$pair_cohorts[used_q], function(pc) {
    for (f in intersect(c("cell_ids", "pair_cell_id", "s_cell_id", "t_cell_id"), names(pc))) pc[[f]] <- re(pc[[f]], cells_k)
    pc$block_ids <- re(pc$block_ids, used_b); pc
  })
  links <- lapply(sb$links[kept], function(l) { l$pair_cohort_ids <- re(l$pair_cohort_ids, used_q); l })
  view <- list(selection_model = sb$selection_model, selection_pooling = sb$selection_pooling,
    universe = sb$universe, design_spec = sb$design_spec, cells = cells, blocks = blocks,
    pair_cohorts = pcs, links = links)
  is_date <- identical(sb$selection_model, "independent_period_logit")
  open <- function(blk) {
    op <- selection_prediction_operation(sb, source_reader = reader)
    list(abort = op$abort, finish = op$finish, fns = list(
      influence = function(rows, m) op$upstream(up_map, rows, blk[[m]]$uids)$up,
      cell = function(rows, cid) {
        old <- cells_k[cid]; cell <- sb$cells[[old]]
        list(p = op$cell(rows, old),
             y = if (is_date) cp04_date_rows(reader, rows, cell$date_spec)$R else cp04_pair_rows(reader, rows, cell$pair)$S)
      },
      pair = function(rows, q) {
        pc <- sb$pair_cohorts[[used_q[q]]]
        op$cell(rows, pc$s_cell_id) * op$cell(rows, pc$t_cell_id)
      }))
  }
  aug <- identical(r$core, "reconstruction")
  src <- list(sel = view, up_map = list(theta_coords = up_map$theta_coords[used_b]), reader = reader,
    ctrl = sb$control, glm_control = sb$glm_control, open = open,
    route = paste0(class(b)[1L], "$selection_backend (shared cdid_selection_models backend; cross-pair ",
      if (aug) "augmented" else "weighted", " route)"),
    role = if (aug) 'observation propensity q_c of the cross-pair augmented observation correction (selection_est_method = "dr")' else
      "observation-selection propensity of the cross-pair weighted route (conditional_mtar)",
    use = if (aug) chained_nuisance_augmented_use else paste("p, this fit's pair probability, enters the observed-pair weight",
      "w a_k / (rho max(p, eps)) of each listed cross-pair link: its cohort pair-change mean, the treatment-logit",
      "fitting weights and the outcome-regression weights"),
    weights_use = paste("The repeat rates a_k and rho enter downstream weights, not this logit fit, so the reported SEs",
      "are not conditional on any estimated weight."))
  comp <- chained_nuisance_selection_run(res, src, prep, weightsname)
  if (identical(comp$status, "available") && length(kept)) {
    lrow <- lt$row(kept)
    comp$comparisons <- data.frame(comparison = lt$id(kept), g = chained_nuisance_g(as.double(lrow$g)),
      s = as.double(lrow$s), t = as.double(lrow$t), link_role = as.character(lrow$link_role),
      controls = ifelse(lrow$link_role == "target", "", chained_nuisance_num(chained_nuisance_g(as.double(lrow$cohort)))),
      pair_cohorts = vapply(links, function(l) paste(l$pair_cohort_ids, collapse = ","), ""), stringsAsFactors = FALSE)
    comp$unused_blocks_after_reduction <- length(sb$blocks) - length(used_b)
  }
  comp
}

cross_pair_nuisance_absent <- function(role, inp) {
  r <- inp$route; em <- r$est_method
  if (role == "imputation")
    return(chained_nuisance_status("not_applicable", "no_imputation", paste0(
      "This cross-pair route (", r$core, ") fits no actual-change imputation.")))
  if (role == "treatment" && identical(em, "reg"))
    return(chained_nuisance_status("not_applicable", "no_treatment_propensity", 'est_method = "reg" fits no treatment propensity.'))
  if (role == "outcome" && identical(em, "ipw"))
    return(chained_nuisance_status("not_applicable", "no_outcome_regression", 'est_method = "ipw" fits no outcome regression.'))
  if (isTRUE(r$unconditional) || (role == "treatment" && isTRUE(r$uncond_treatment)))
    return(structure(chained_nuisance_status("not_applicable", "closed_form_unconditional", paste0(
      "The ", if (role == "treatment") "treatment" else "outcome", " design is a single nonzero constant column, so ",
      "the cross-pair moments use closed-form weighted cell means; no ", if (role == "treatment") "logit" else "regression",
      " was fitted (", r$n_links, " link(s)", if (!is.null(r$n_cell_means) && !is.na(r$n_cell_means))
        paste0(", ", r$n_cell_means, " distinct cohort-mask-pair means") else "",
      if (role == "treatment" && !is.null(r$n_closed_form_treatment)) paste0(", ", r$n_closed_form_treatment,
        " closed-form control weightings") else "", ")."))))
  chained_nuisance_status("not_applicable", "no_fit_used", paste0(
    "No ", role, " fit of this route is read by a retained cross-pair link."))
}

chained_cross_pair_nuisance_capture <- function(res, prep, gate, est_label, weightsname) {
  inp <- res$nuisance_capture_input
  if (!is.list(inp) || !identical(inp$schema, "cdid.cross_pair_nuisance_input.v0")) {
    st <- chained_nuisance_status("not_captured", "not_requested",
      "This cross-pair fit was estimated without keep_nuisance = TRUE; no capture input was retained.")
    cap <- chained_nuisance_capture_assemble(st, gate, est_label, cross_pair_nuisance_status_all(st))
    cap$entries <- unname(cp04_payload(cap)["entries"])
    return(cap)
  }
  if (!identical(inp$status, "fitted")) {
    st <- chained_nuisance_status("not_applicable", inp$status, inp$reason)
    cap <- chained_nuisance_capture_assemble(st, gate, est_label, cross_pair_nuisance_status_all(st))
    cap$cross_pair <- inp$route
    cap$entries <- unname(cp04_payload(cap)["entries"])
    return(cap)
  }
  lt <- cross_pair_nuisance_link_table(res, inp)
  weights_text <- chained_nuisance_aug_weights_text(weightsname, prep)
  formula_text <- if (is.null(prep$cp04$outcome_formula)) NA_character_ else
    paste(deparse(prep$cp04$outcome_formula), collapse = " ")
  # Output admission on the pipeline's recorded live entries (existing rule).
  ctrl <- if (!is.null(res$nuisance_registry$cp04$selection_backend)) res$nuisance_registry$cp04$selection_backend$control else
    res$controls$pipeline_control$owner_control
  ra <- res$resource_accounting
  live <- ceiling(max(c(0, ra$simultaneous_entries, ra$initial_envelope_entries)))
  n_coef <- sum(vapply(c(inp$treatment, inp$outcome, inp$imputation), function(r) length(r$estimate), 0))
  n_fit <- length(inp$treatment) + length(inp$outcome) + length(inp$imputation)
  out_bound <- ceiling(3 * 65536 + 2048 * n_fit + 512 * n_coef + 256 * lt$M)
  if (is.list(ctrl)) {
    prior <- if (is.null(ra$downstream_work_ops)) 0 else ceiling(ra$downstream_work_ops)
    adm <- tryCatch(cp04_budget(ctrl, c(live, out_bound), out_bound, c(prior, 8 * out_bound),
      out_bound, "cross-pair nuisance tables"), cdid_ipw_owner_bound_error = identity)
    if (inherits(adm, "error")) {
      st <- chained_nuisance_status("unavailable", "budget",
        paste("The optional coefficient tables exceed the declared resource limits:", conditionMessage(adm)))
      cap <- chained_nuisance_capture_assemble(st, gate, est_label, cross_pair_nuisance_status_all(st))
      cap$entries <- unname(cp04_payload(cap)["entries"])
      return(cap)
    }
  }
  roles <- list()
  if (isTRUE(inp$route$reference)) {
    roles <- cross_pair_nuisance_status_all(chained_nuisance_status("unavailable", "reference_recipe_route",
      "Dated covariates or recipe-indexed weights fit the cross-pair nuisances on the recipe layout; outside this capture."))
  } else for (role in c("outcome", "treatment", "imputation")) {
    recs <- unname(cross_pair_nuisance_used(inp[[role]], lt))
    roles[[role]] <- if (length(recs))
      cross_pair_nuisance_role_component(role, recs, inp, lt, weights_text, prep$clustervar, formula_text) else
      cross_pair_nuisance_absent(role, inp)
  }
  selection <- cross_pair_nuisance_selection(res, inp, prep, weightsname, lt)
  actual <- unname(cp04_payload(roles)["entries"])
  if (actual > out_bound)
    roles <- cross_pair_nuisance_status_all(chained_nuisance_status("unavailable", "budget",
      "The cross-pair coefficient tables exceeded their reserved output size."))
  cap <- chained_nuisance_capture_assemble(selection, gate, est_label, roles)
  unused <- vapply(c("treatment", "outcome", "imputation"), function(role)
    length(inp[[role]]) - length(cross_pair_nuisance_used(inp[[role]], lt)), 0)
  cap$cross_pair <- c(inp$route[setdiff(names(inp$route), "links")], list(
    n_links_catalog = lt$M, n_links_retained = length(lt$kept),
    fits_unused_after_reduction = as.list(unused),
    note = "Fits read only by links pruned by the verified forest reduction are counted here, not tabulated."))
  cap$entries <- unname(cp04_payload(cap)["entries"])
  cap
}
