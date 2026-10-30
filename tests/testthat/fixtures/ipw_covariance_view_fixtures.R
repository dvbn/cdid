# Self-contained base-R fixtures for CP-02 compact-view IPW covariance transfer tests.
# Uses no external file reads and no runtime glm/Newton fitting.
list(
  build_14unit_3link_view_and_dense = function(
    ns,
    r_ball = 1e-6,
    chunk_rows = 5L,
    operator_block_rows = 4L,
    id_mode = c("character", "numeric")
  ) {
    id_mode <- match.arg(id_mode)
    n <- 14L
    id <- if (identical(id_mode, "character")) {
      sprintf("unit_%02d", 0:13)
    } else {
      as.double(100L + 0:13)
    }
    G <- c(rep(2.0, 4L), rep(5.0, 4L), rep(99.0, 4L), rep(9.0, 2L))
    X <- c(rep(c(0, 1, 1, 0), 3L), 0, 1)
    Z <- cbind(intercept = rep(1.0, n), x = as.double(X))
    S <- as.logical(c(rep(c(1L, 1L, 1L, 0L), 3L), 1L, 0L))
    w <- c(1.0, 2.0, 0.5, 1.5, 2.0, 1.0, 1.5, 1.0, 1.0, 1.5, 2.0, 0.5, 1.0, 1.0)
    y <- matrix(0.0, nrow = n, ncol = 3L)
    for (i in 0:13) {
      for (j in 0:2) {
        y[i + 1L, j + 1L] <- (((i + 1L) * (j + 2L)) %% 17L - 5L) / 3.0 + (i * j) / 11.0
      }
    }

    # Periods 4 and 5 are observed for all 14 units so all 14 units have no >= 2L and belong to A,
    # while periods 1..3 are observed iff S[i] == TRUE.
    pan_rows <- vector("list", n)
    for (i in seq_len(n)) {
      pan_rows[[i]] <- data.frame(
        id = rep(id[i], 5L),
        time = c(1.0, 2.0, 3.0, 4.0, 5.0),
        y = c(if (S[i]) y[i, ] else rep(NA_real_, 3L), 0.0, 0.0),
        g = rep(G[i], 5L),
        stringsAsFactors = FALSE
      )
    }
    pan_df <- do.call(rbind, pan_rows)

    calc_root <- function(cohorts) {
      M <- vapply(cohorts, function(g_val) sum(w[G == g_val]), double(1L))
      P <- vapply(cohorts, function(g_val) sum(w[G == g_val & S]), double(1L))
      r_vec <- P / M
      idx <- which((G %in% cohorts) & S)
      ci <- w[idx] / r_vec[match(G[idx], cohorts)]
      odds <- vapply(c(0, 1), function(xv) {
        sum(ci[X[idx] == xv & G[idx] == 2.0]) / sum(ci[X[idx] == xv & G[idx] != 2.0])
      }, double(1L))
      beta0 <- c(log(odds[1L]), log(odds[2L]) - log(odds[1L]))
      eta <- as.vector(Z[idx, , drop = FALSE] %*% beta0)
      p <- 1.0 / (1.0 + exp(-eta))
      H <- crossprod(Z[idx, , drop = FALSE], Z[idx, , drop = FALSE] * (ci * p * (1.0 - p))) / 14.0
      list(beta0 = beta0, H = H, F_mat = chol(H), ci = ci, idx = idx, odds = odds)
    }
    r12 <- calc_root(c(2.0, 5.0, 99.0))
    r3  <- calc_root(c(2.0, 5.0))

    calc_lambda_beta <- function(rt, s_col, t_col) {
      idx <- rt$idx
      ci <- rt$ci
      delta <- y[idx, t_col] - y[idx, s_col]
      b <- ci * ifelse(X[idx] == 0, rt$odds[1L], rt$odds[2L])
      ctrl <- G[idx] != 2.0
      C_val <- sum(b[ctrl] * delta[ctrl]) / sum(b[ctrl])
      v <- as.vector(crossprod(Z[idx[ctrl], , drop = FALSE], b[ctrl] * (delta[ctrl] - C_val))) / sum(b[ctrl])
      as.vector(solve(rt$H, -v))
    }
    lb_1 <- calc_lambda_beta(r12, 1L, 2L)
    lb_2 <- calc_lambda_beta(r12, 2L, 3L)
    lb_3 <- calc_lambda_beta(r3,  1L, 3L)

    o_ctrl <- ns$ipw_owner_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_retained_entries = 50000L,
      max_matrix_entries = 50000L,
      max_work_ops = 500000L,
      max_output_entries = 50000L
    )
    panel_owner <- ns$construct_ipw_panel_owner(
      data = pan_df, idname = "id", tname = "time", yname = "y", gname = "g", control = o_ctrl
    )
    frame_owner <- ns$construct_ipw_frame_owner(
      panel_owner = panel_owner, id = id, Z = Z, weights = w, control = o_ctrl
    )

    cc <- ns$logit_certificate_control(18L, 40.0, 16L, as.integer(chunk_rows), 5000L, 500000L)
    pc <- ns$ipw_nuisance_problem_control(as.integer(chunk_rows), 10000L, 5000L, 500000L)
    lc <- ns$ipw_link_enclosure_control(as.integer(chunk_rows), 10000L, 5000L, 500000L, 10000L)
    cov_ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_operator_block_rows = as.integer(operator_block_rows),
      max_cluster_rows = 14L,
      max_matrix_entries = 15000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = lc
    )

    nv_12 <- ns$construct_ipw_nuisance_view(
      frame_owner = frame_owner, g = 2.0, K = c(5.0, 99.0), s = 1.0, t = 2.0, control = o_ctrl
    )
    lv_1 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 1.0, t = 2.0, g = 2.0, K = c(5.0, 99.0),
      nuisance_view = nv_12, control = o_ctrl, link_id = "l1"
    )
    lv_2 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 2.0, t = 3.0, g = 2.0, K = c(5.0, 99.0),
      nuisance_view = nv_12, control = o_ctrl, link_id = "l2"
    )
    nv_3 <- ns$construct_ipw_nuisance_view(
      frame_owner = frame_owner, g = 2.0, K = 5.0, s = 1.0, t = 3.0, control = o_ctrl
    )
    lv_3 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 1.0, t = 3.0, g = 2.0, K = 5.0,
      nuisance_view = nv_3, control = o_ctrl, link_id = "l3"
    )

    prob12_view <- ns$construct_ipw_enclosure_problem_view(
      nuisance_view = nv_12, beta0 = r12$beta0, F = r12$F_mat, r = r_ball,
      certificate_control = cc, problem_control = pc, reader_control = o_ctrl
    )
    prob3_view <- ns$construct_ipw_enclosure_problem_view(
      nuisance_view = nv_3, beta0 = r3$beta0, F = r3$F_mat, r = r_ball,
      certificate_control = cc, problem_control = pc, reader_control = o_ctrl
    )

    l1_view <- ns$enclose_ipw_link_view(prob12_view, lv_1, lb_1, lc, o_ctrl)
    l2_view <- ns$enclose_ipw_link_view(prob12_view, lv_2, lb_2, lc, o_ctrl)
    l3_view <- ns$enclose_ipw_link_view(prob3_view,  lv_3, lb_3, lc, o_ctrl)

    prob12_dense <- ns$construct_certified_ipw_nuisance_problem(
      id = id, Z = Z, w = w, G = G, pair_observed = S, cohorts = c(2.0, 5.0, 99.0),
      beta0 = r12$beta0, F = r12$F_mat, r = r_ball, certificate_control = cc, problem_control = pc
    )
    prob3_dense <- ns$construct_certified_ipw_nuisance_problem(
      id = id, Z = Z, w = w, G = G, pair_observed = S, cohorts = c(2.0, 5.0),
      beta0 = r3$beta0, F = r3$F_mat, r = r_ball, certificate_control = cc, problem_control = pc
    )
    l1_dense <- ns$enclose_ipw_link(prob12_dense, id, y[, 1L], y[, 2L], lb_1, lc, link_id = "l1")
    l2_dense <- ns$enclose_ipw_link(prob12_dense, id, y[, 2L], y[, 3L], lb_2, lc, link_id = "l2")
    l3_dense <- ns$enclose_ipw_link(prob3_dense,  id, y[, 1L], y[, 3L], lb_3, lc, link_id = "l3")

    phi_cand <- structure(c(
      -0x1.5d55b01e7f943p+2, 0x1.d80cfac2a3c97p+1, 0x1.c53ccaf4b6bdfp+0, 0x0.0p+0,
      -0x1.f8f4972b29972p-4, 0x1.b9c83d1a5a676p+2, 0x1.277f23965329bp+3, 0x1.c78132a108ecap+1,
      -0x1.05fcdcde5326ep+0, -0x1.a1da954faab28p+2, -0x1.4202fb96db4d2p+3, -0x1.f50e1e1789d11p+0,
      0x0.0p+0, 0x0.0p+0,
      0x1.1941a6b8acb4bp+3, -0x1.ed73e711c2195p+2, -0x1.143d997e5d402p+0, 0x0.0p+0,
      0x1.111061f37112bp+1, -0x1.91ca704c3626bp+1, -0x1.7505e8b409e4bp+2, -0x1.84c9b8c77cd42p+0,
      -0x1.a4da13bd5c0e9p+0, 0x1.20862c6253c16p+2, 0x1.29d0606a4755ap+2, 0x1.abaab1a83c82fp-1,
      0x0.0p+0, 0x0.0p+0,
      0x1.f07c1f07c1f08p-1, -0x1.1111111111111p+1, 0x1.29e4129e4129ep+0, 0x0.0p+0,
      0x0.0p+0, 0x1.04a7904a7904ap+1, -0x1.04a7904a7904ap+1, 0x0.0p+0,
      0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0
    ), dim = c(14L, 3L), dimnames = list(NULL, c("l1", "l2", "l3")))

    C_mat <- structure(c(
      0x1.5555555555555p-2, 0x1.999999999999ap-4,
      -0x1.6666666666666p-1, 0x1.999999999999ap-3,
      0x1.199999999999ap+0, -0x1.3333333333333p-2
    ), dim = c(2L, 3L), dimnames = list(c("m1", "m2"), c("l1", "l2", "l3")))

    V_cand_gmm <- structure(c(
      0x1.004292f3cee41p+4, -0x1.e2cfed42e3ab3p-1,
      -0x1.e2cfed42e3ab3p-1, 0x1.3a4df88b6f8ccp-3
    ), dim = c(2L, 2L), dimnames = list(c("m1", "m2"), c("m1", "m2")))
    ref_exact_V_gmm <- structure(c(
      16.016253425956396, -0.94299260560531761,
      -0.94299260560531761, 0.15346903014682878
    ), dim = c(2L, 2L), dimnames = list(c("m1", "m2"), c("m1", "m2")))

    V_cand_cl <- structure(c(
      0x1.091137650c965p+3, -0x1.0da591370f63ep+0,
      -0x1.0da591370f63ep+0, 0x1.142d7d0b80e55p-3
    ), dim = c(2L, 2L), dimnames = list(c("m1", "m2"), c("m1", "m2")))
    ref_exact_V_cl <- structure(c(
      8.2833516095749449, -1.0533076056884876,
      -1.0533076056884876, 0.13485238734584706
    ), dim = c(2L, 2L), dimnames = list(c("m1", "m2"), c("m1", "m2")))

    labels_cl <- c("b", "d", "a", "c", "b", "d", "a", "c", "b", "d", "a", "c", "zero", "zero")
    order_cl  <- c("c", "zero", "a", "d", "b")
    omega_cl  <- c(0x1.999999999999ap-4, 0x1.6666666666666p+0, 0x1.5555555555555p-2, 0x1.ccccccccccccdp-1, 0x1.b333333333333p+0)
    kappa_14  <- 0x1.2492492492492p-4

    callback <- function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        id = id[idx],
        link_ids = c("l1", "l2", "l3"),
        phi = phi_cand[idx, , drop = FALSE],
        a_index = idx
      )
    }

    op_gmm <- ns$ipw_covariance_operator_spec(
      purpose = "GMM",
      row_level = "unit",
      omega = rep(1.0, n),
      centered = TRUE,
      kappa = kappa_14
    )
    op_cl <- ns$ipw_covariance_operator_spec(
      purpose = "inference",
      row_level = "cluster",
      omega = omega_cl,
      centered = TRUE,
      kappa = kappa_14,
      cluster_id = labels_cl,
      cluster_order = order_cl
    )

    list(
      n = n,
      id = id,
      Z = Z,
      w = w,
      G = G,
      S = S,
      y = y,
      pan_df = pan_df,
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      nv_12 = nv_12,
      nv_3 = nv_3,
      lv_1 = lv_1,
      lv_2 = lv_2,
      lv_3 = lv_3,
      prob12_view = prob12_view,
      prob3_view = prob3_view,
      prob12_dense = prob12_dense,
      prob3_dense = prob3_dense,
      links_view = list(l1 = l1_view, l2 = l2_view, l3 = l3_view),
      links_dense = list(l1 = l1_dense, l2 = l2_dense, l3 = l3_dense),
      frame_dense = ns$ipw_covariance_frame_spec(id = id, normalizer_n = n, embedding = "identity", w = w, G = G),
      phi_cand = phi_cand,
      callback = callback,
      candidate_rows = ns$ipw_candidate_rows_spec(callback = callback, link_ids = c("l1", "l2", "l3")),
      C = C_mat,
      V_cand_gmm = V_cand_gmm,
      ref_exact_V_gmm = ref_exact_V_gmm,
      V_cand_cl = V_cand_cl,
      ref_exact_V_cl = ref_exact_V_cl,
      labels_cl = labels_cl,
      order_cl = order_cl,
      omega_cl = omega_cl,
      kappa = kappa_14,
      op_gmm = op_gmm,
      op_cl = op_cl,
      o_ctrl = o_ctrl,
      cert_ctrl = cc,
      prob_ctrl = pc,
      link_ctrl = lc,
      cov_ctrl = cov_ctrl
    )
  },

  build_singular_observer_q_view_and_dense = function(ns, r_ball = 1e-6, chunk_rows = 5L, operator_block_rows = 3L) {
    n <- 14L
    id <- sprintf("unit_%02d", 0:13)
    G <- c(rep(2.0, 4L), rep(5.0, 4L), rep(99.0, 4L), rep(9.0, 2L))
    X <- c(rep(c(0, 1, 1, 0), 3L), 0, 1)
    Z <- cbind(intercept = rep(1.0, n), x = as.double(X))
    S <- as.logical(c(rep(c(1L, 1L, 1L, 0L), 3L), 1L, 0L))
    w <- c(1.0, 2.0, 0.5, 1.5, 2.0, 1.0, 1.5, 1.0, 1.0, 1.5, 2.0, 0.5, 1.0, 1.0)
    y <- matrix(0.0, nrow = n, ncol = 3L)
    for (i in 0:13) {
      for (j in 0:2) {
        y[i + 1L, j + 1L] <- (((i + 1L) * (j + 2L)) %% 17L - 5L) / 3.0 + (i * j) / 11.0
      }
    }

    pan_rows <- vector("list", n)
    for (i in seq_len(n)) {
      pan_rows[[i]] <- data.frame(
        id = rep(id[i], 5L),
        time = c(1.0, 2.0, 3.0, 4.0, 5.0),
        y = c(if (S[i]) y[i, ] else rep(NA_real_, 3L), 0.0, 0.0),
        g = rep(G[i], 5L),
        stringsAsFactors = FALSE
      )
    }
    pan_df <- do.call(rbind, pan_rows)

    calc_root <- function(cohorts) {
      M <- vapply(cohorts, function(g_val) sum(w[G == g_val]), double(1L))
      P <- vapply(cohorts, function(g_val) sum(w[G == g_val & S]), double(1L))
      r_vec <- P / M
      idx <- which((G %in% cohorts) & S)
      ci <- w[idx] / r_vec[match(G[idx], cohorts)]
      odds <- vapply(c(0, 1), function(xv) {
        sum(ci[X[idx] == xv & G[idx] == 2.0]) / sum(ci[X[idx] == xv & G[idx] != 2.0])
      }, double(1L))
      beta0 <- c(log(odds[1L]), log(odds[2L]) - log(odds[1L]))
      eta <- as.vector(Z[idx, , drop = FALSE] %*% beta0)
      p <- 1.0 / (1.0 + exp(-eta))
      H <- crossprod(Z[idx, , drop = FALSE], Z[idx, , drop = FALSE] * (ci * p * (1.0 - p))) / 14.0
      list(beta0 = beta0, H = H, F_mat = chol(H), ci = ci, idx = idx, odds = odds)
    }
    calc_lambda_beta <- function(rt, s_col, t_col) {
      idx <- rt$idx
      ci <- rt$ci
      delta <- y[idx, t_col] - y[idx, s_col]
      b <- ci * ifelse(X[idx] == 0, rt$odds[1L], rt$odds[2L])
      ctrl <- G[idx] != 2.0
      C_val <- sum(b[ctrl] * delta[ctrl]) / sum(b[ctrl])
      v <- as.vector(crossprod(Z[idx[ctrl], , drop = FALSE], b[ctrl] * (delta[ctrl] - C_val))) / sum(b[ctrl])
      as.vector(solve(rt$H, -v))
    }

    o_ctrl <- ns$ipw_owner_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_retained_entries = 50000L,
      max_matrix_entries = 50000L,
      max_work_ops = 500000L,
      max_output_entries = 50000L
    )
    panel_owner <- ns$construct_ipw_panel_owner(
      data = pan_df, idname = "id", tname = "time", yname = "y", gname = "g", control = o_ctrl
    )
    frame_owner <- ns$construct_ipw_frame_owner(
      panel_owner = panel_owner, id = id, Z = Z, weights = w, control = o_ctrl
    )

    cc <- ns$logit_certificate_control(18L, 40.0, 16L, 14L, 5000L, 500000L)
    pc <- ns$ipw_nuisance_problem_control(14L, 10000L, 5000L, 500000L)
    lc <- ns$ipw_link_enclosure_control(14L, 10000L, 5000L, 500000L, 10000L)
    cov_ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_operator_block_rows = as.integer(operator_block_rows),
      max_cluster_rows = 5L,
      max_matrix_entries = 15000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = lc
    )

    r_a <- calc_root(c(2.0, 5.0, 99.0))
    r_b <- calc_root(c(2.0, 5.0, 9.0))
    lb_a <- calc_lambda_beta(r_a, 1L, 2L)
    lb_b <- calc_lambda_beta(r_b, 1L, 3L)

    nv_a <- ns$construct_ipw_nuisance_view(frame_owner, g = 2.0, K = c(5.0, 99.0), s = 1.0, t = 2.0, control = o_ctrl)
    lv_a <- ns$construct_ipw_link_view(frame_owner, s = 1.0, t = 2.0, g = 2.0, K = c(5.0, 99.0), nuisance_view = nv_a, control = o_ctrl, link_id = "la")
    prob_a_view <- ns$construct_ipw_enclosure_problem_view(nv_a, r_a$beta0, r_a$F_mat, r_ball, cc, pc, o_ctrl)
    la_view <- ns$enclose_ipw_link_view(prob_a_view, lv_a, lb_a, lc, o_ctrl)

    nv_b <- ns$construct_ipw_nuisance_view(frame_owner, g = 2.0, K = c(5.0, 9.0), s = 1.0, t = 3.0, control = o_ctrl)
    lv_b <- ns$construct_ipw_link_view(frame_owner, s = 1.0, t = 3.0, g = 2.0, K = c(5.0, 9.0), nuisance_view = nv_b, control = o_ctrl, link_id = "lb")
    prob_b_view <- ns$construct_ipw_enclosure_problem_view(nv_b, r_b$beta0, r_b$F_mat, r_ball, cc, pc, o_ctrl)
    lb_view <- ns$enclose_ipw_link_view(prob_b_view, lv_b, lb_b, lc, o_ctrl)

    prob_a_dense <- ns$construct_certified_ipw_nuisance_problem(id, Z, w, G, S, c(2.0, 5.0, 99.0), r_a$beta0, r_a$F_mat, r_ball, cc, pc)
    prob_b_dense <- ns$construct_certified_ipw_nuisance_problem(id, Z, w, G, S, c(2.0, 5.0, 9.0),  r_b$beta0, r_b$F_mat, r_ball, cc, pc)
    la_dense <- ns$enclose_ipw_link(prob_a_dense, id, y[, 1L], y[, 2L], lb_a, lc, link_id = "la")
    lb_dense <- ns$enclose_ipw_link(prob_b_dense, id, y[, 1L], y[, 3L], lb_b, lc, link_id = "lb")

    cha <- ns$enclose_ipw_link_if_chunk(la_dense, 1L, n, lc)
    chb <- ns$enclose_ipw_link_if_chunk(lb_dense, 1L, n, lc)
    phi_full <- cbind(
      la = 0.5 * (cha$phi_0[, 1L] + cha$phi_0[, 2L]),
      lb = 0.5 * (chb$phi_0[, 1L] + chb$phi_0[, 2L])
    )
    phi_obs <- phi_full
    phi_obs[!S, ] <- 0.0

    cl_id <- ifelse(S, "c_obs", paste0("c_nonobs_", as.integer(G)))
    cl_ord <- c("c_obs", "c_nonobs_2", "c_nonobs_5", "c_nonobs_99", "c_nonobs_9")
    omega <- rep(1.0, 5L)
    kappa <- 1.0 / 14.0
    C_eye <- diag(1.0, 2L)
    dimnames(C_eye) <- list(c("la", "lb"), c("la", "lb"))

    W_full <- t(vapply(cl_ord, function(cl) colSums(phi_full[cl_id == cl, , drop = FALSE]), double(2L)))
    W_obs  <- t(vapply(cl_ord, function(cl) colSums(phi_obs[cl_id == cl, , drop = FALSE]), double(2L)))
    V_full <- (crossprod(W_full) - tcrossprod(colSums(W_full)) / 5.0) * kappa
    Q_obs  <- (crossprod(W_obs)  - tcrossprod(colSums(W_obs))  / 5.0) * kappa
    V_full <- 0.5 * (V_full + t(V_full))
    Q_obs  <- 0.5 * (Q_obs + t(Q_obs))
    dimnames(V_full) <- dimnames(C_eye)
    dimnames(Q_obs)  <- dimnames(C_eye)

    cand_full <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = c("la", "lb"), phi = phi_full[idx, , drop = FALSE], a_index = idx)
      },
      link_ids = c("la", "lb")
    )
    cand_obs_dropped <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = c("la", "lb"), phi = phi_obs[idx, , drop = FALSE], a_index = idx)
      },
      link_ids = c("la", "lb")
    )

    op_cl <- ns$ipw_covariance_operator_spec(
      purpose = "inference",
      row_level = "cluster",
      omega = omega,
      centered = TRUE,
      kappa = kappa,
      cluster_id = cl_id,
      cluster_order = cl_ord
    )

    list(
      n = n,
      id = id,
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      links_view = list(la = la_view, lb = lb_view),
      links_dense = list(la = la_dense, lb = lb_dense),
      frame_dense = ns$ipw_covariance_frame_spec(id = id, normalizer_n = n, embedding = "identity", w = w, G = G),
      phi_full = phi_full,
      phi_obs = phi_obs,
      C = C_eye,
      V_full = V_full,
      Q_obs = Q_obs,
      candidate_rows_full = cand_full,
      candidate_rows_obs_dropped = cand_obs_dropped,
      op = op_cl,
      o_ctrl = o_ctrl,
      cov_ctrl = cov_ctrl
    )
  },

  build_10unit_enclosure_fixtures_view_and_dense = function(ns, chunk_rows = 4L, operator_block_rows = 4L) {
    # Exact 10-unit rational and nonzero-beta root-ball fixtures matching ipw_enclosure_reference_fixtures.R
    n <- 10L
    id <- as.double(1:10)
    G <- as.double(c(2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L, 10L, 99L))
    po <- c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
    cohorts <- c(2.0, 3.0, 4.0)

    cc <- ns$logit_certificate_control(18L, 40.0, 16L, as.integer(chunk_rows), 2500L, 500000L)
    pc <- ns$ipw_nuisance_problem_control(as.integer(chunk_rows), 5000L, 2500L, 500000L)
    lc <- ns$ipw_link_enclosure_control(10L, 5000L, 2500L, 500000L, 5000L)
    o_ctrl <- ns$ipw_owner_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_retained_entries = 50000L,
      max_matrix_entries = 50000L,
      max_work_ops = 500000L,
      max_output_entries = 50000L
    )
    cov_ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_operator_block_rows = as.integer(operator_block_rows),
      max_cluster_rows = 10L,
      max_matrix_entries = 10000L,
      max_work_ops = 500000L,
      max_output_entries = 5000L,
      link_control = lc
    )
    op_10 <- ns$ipw_covariance_operator_spec(
      purpose = "GMM", row_level = "unit", omega = rep(1.0, n), centered = TRUE, kappa = 0.1
    )

    # --- Fixture 1: Exact rational ---
    Z1 <- structure(c(
      0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0,
      0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0,
      -0x1.0000000000000p-1, 0x1.8000000000000p+0, 0x1.0000000000000p-2, -0x1.0000000000000p+0, 0x1.0000000000000p+0,
      -0x1.0000000000000p-1, 0x1.0000000000000p-1, -0x1.8000000000000p-1, 0x1.0000000000000p+1, -0x1.0000000000000p+1
    ), dim = c(10L, 2L), dimnames = list(NULL, c("z1", "z2")))
    w1 <- c(0x1.8000000000000p+0, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.0000000000000p-2,
            0x1.0000000000000p-1, 0x1.8000000000000p-1, 0x1.8000000000000p-1, 0x1.4000000000000p+0, 0x1.c000000000000p+0)
    F1 <- structure(c(0x1.8000000000000p-2, 0x0.0p+0, 0x1.0000000000000p-5, 0x1.4000000000000p-2), dim = c(2L, 2L))
    yL1_1 <- c(0x1.0000000000000p-1, 0x1.8000000000000p+0, 0x1.f380000000000p+9, 0x1.0000000000000p-2, 0x1.4000000000000p+0,
               -0x1.f380000000000p+9, -0x1.0000000000000p-1, 0x1.8480000000000p+9, 0x1.bc00000000000p+9, -0x1.bc00000000000p+9)
    yR1_1 <- c(0x1.4000000000000p+1, 0x1.4000000000000p+1, -0x1.f380000000000p+9, 0x1.4000000000000p+0, 0x1.a000000000000p+1,
               0x1.f380000000000p+9, 0x1.0000000000000p-1, -0x1.8480000000000p+9, -0x1.bc00000000000p+9, 0x1.bc00000000000p+9)
    yL1_2 <- c(0x0.0p+0, 0x1.0000000000000p+0, 0x1.1580000000000p+9, 0x1.0000000000000p+1, 0x0.0p+0,
               -0x1.1580000000000p+9, 0x1.0000000000000p-2, 0x1.bc00000000000p+8, 0x1.4d00000000000p+8, -0x1.4d00000000000p+8)
    yR1_2 <- c(-0x1.0000000000000p+0, 0x1.0000000000000p+2, -0x1.1580000000000p+9, 0x1.0000000000000p-1, 0x1.0000000000000p-1,
               0x1.1580000000000p+9, 0x1.2000000000000p+1, -0x1.bc00000000000p+8, -0x1.4d00000000000p+8, 0x1.4d00000000000p+8)
    lb1_1 <- c(0x0.0p+0, -0x1.364d9364d9365p+0)
    lb1_2 <- c(0x0.0p+0, -0x1.5d1745d1745d1p+3)
    phi1_1 <- c(0x1.4745d1745d174p+1, -0x1.4745d1745d174p+1, 0x0.0p+0, 0x1.d1745d1745d17p-5, -0x1.8d9364d9364d9p-1,
                -0x1.707c1f07c1f08p-2, 0x1.145d1745d1746p-1, 0x1.145d1745d1746p-1, 0x0.0p+0, 0x0.0p+0)
    phi1_2 <- c(-0x1.5d1745d1745d1p+0, 0x1.5d1745d1745d1p+0, 0x0.0p+0, 0x1.05d1745d1745dp-1, 0x1.c2e8ba2e8ba2fp+0,
                0x1.22e8ba2e8ba2fp+0, -0x1.b45d1745d1746p+0, -0x1.b45d1745d1746p+0, 0x0.0p+0, 0x0.0p+0)
    phi_10 <- cbind(link1 = phi1_1, link2 = phi1_2)
    V_10 <- crossprod(scale(phi_10, center = TRUE, scale = FALSE)) / 10.0
    V_10 <- 0.5 * (V_10 + t(V_10))

    pan1_rows <- vector("list", n)
    for (i in seq_len(n)) {
      pan1_rows[[i]] <- data.frame(
        id = rep(id[i], 6L),
        time = c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
        y = c(if (po[i]) c(yL1_1[i], yR1_1[i], yL1_2[i], yR1_2[i]) else rep(NA_real_, 4L), 0.0, 0.0),
        g = rep(G[i], 6L),
        stringsAsFactors = FALSE
      )
    }
    pan1_df <- do.call(rbind, pan1_rows)
    po1_owner <- ns$construct_ipw_panel_owner(pan1_df, "id", "time", "y", "g", o_ctrl)
    fo1_owner <- ns$construct_ipw_frame_owner(po1_owner, id = id, Z = Z1, weights = w1, control = o_ctrl)
    nv1 <- ns$construct_ipw_nuisance_view(fo1_owner, g = 2.0, K = c(3.0, 4.0), s = 1.0, t = 2.0, control = o_ctrl)
    lv1_1 <- ns$construct_ipw_link_view(fo1_owner, s = 1.0, t = 2.0, g = 2.0, K = c(3.0, 4.0), nuisance_view = nv1, control = o_ctrl, link_id = "link1")
    lv1_2 <- ns$construct_ipw_link_view(fo1_owner, s = 3.0, t = 4.0, g = 2.0, K = c(3.0, 4.0), nuisance_view = nv1, control = o_ctrl, link_id = "link2")
    pr1_v <- ns$construct_ipw_enclosure_problem_view(nv1, c(0.0, 0.0), F1, 1e-5, cc, pc, o_ctrl)
    lk1_1_v <- ns$enclose_ipw_link_view(pr1_v, lv1_1, lb1_1, lc, o_ctrl)
    lk1_2_v <- ns$enclose_ipw_link_view(pr1_v, lv1_2, lb1_2, lc, o_ctrl)

    pr1_d <- ns$construct_certified_ipw_nuisance_problem(id, Z1, w1, G, po, cohorts, c(0.0, 0.0), F1, 1e-5, cc, pc)
    lk1_1_d <- ns$enclose_ipw_link(pr1_d, id, yL1_1, yR1_1, lb1_1, lc, link_id = "link1")
    lk1_2_d <- ns$enclose_ipw_link(pr1_d, id, yL1_2, yR1_2, lb1_2, lc, link_id = "link2")

    cand_10 <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = c("link1", "link2"), phi = phi_10[idx, , drop = FALSE], a_index = idx)
      },
      link_ids = c("link1", "link2")
    )

    # --- Fixture 2: Nonzero beta0 root-ball ---
    Z2 <- structure(c(
      0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0,
      0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0,
      0x1.0000000000000p+0, -0x1.0000000000000p+0, 0x1.0000000000000p-2, -0x1.0000000000000p+0, 0x1.0000000000000p-1,
      -0x1.0000000000000p-1, 0x1.0000000000000p-2, -0x1.8000000000000p-1, 0x1.0000000000000p+1, -0x1.0000000000000p+1
    ), dim = c(10L, 2L), dimnames = list(NULL, c("z1", "z2")))
    w2 <- c(0x1.2000000000000p+1, 0x1.8000000000000p-1, 0x1.0000000000000p+0, 0x1.8000000000000p-1, 0x1.0000000000000p-2,
            0x1.0000000000000p-1, 0x1.8000000000000p-1, 0x1.8000000000000p-1, 0x1.4000000000000p+0, 0x1.c000000000000p+0)
    F2 <- structure(c(0x1.8000000000000p-2, 0x0.0p+0, 0x1.0000000000000p-4, 0x1.4000000000000p-2), dim = c(2L, 2L))
    beta0_2 <- c(0x1.b02c777f03731p-4, 0x1.05282f35e5b34p+0)
    r_ball_2 <- 0x1.47ae147ae147bp-7
    yL2 <- c(0x1.0000000000000p-1, 0x1.8000000000000p+0, 0x1.f380000000000p+9, 0x1.0000000000000p-2, 0x1.4000000000000p+0,
             -0x1.f380000000000p+9, -0x1.0000000000000p-1, 0x1.8480000000000p+9, 0x1.bc00000000000p+9, -0x1.bc00000000000p+9)
    yR2 <- c(0x1.4000000000000p+1, 0x1.c000000000000p+1, -0x1.f380000000000p+9, -0x1.c000000000000p+0, 0x1.d000000000000p+2,
             0x1.f380000000000p+9, 0x1.8000000000000p+0, -0x1.8480000000000p+9, -0x1.bc00000000000p+9, 0x1.bc00000000000p+9)
    lb2 <- c(0x1.aa29634633bb4p+0, -0x1.2d05405441151p+3)
    phi2_0 <- c(-0x1.86b7a23ce44d5p+2, 0x1.f107c8dba4daap+2, 0x1.1b8066fcac238p-1, 0x1.bfd4fda7d106ap+1, -0x1.aa96a4cfcecb3p+2,
                -0x1.95584bf7cc8fcp+0, 0x1.445f8114a9385p+0, 0x1.445f8114a9385p+0, 0x0.0p+0, 0x0.0p+0)
    phi2_star <- c(-0x1.86b5eb224cd17p+2, 0x1.f1134be199de8p+2, 0x1.1ba3aca8cd782p-1, 0x1.bfd65625585d5p+1, -0x1.aa9a6ef59b960p+2,
                   -0x1.955e87c5decebp+0, 0x1.446a1f0000bdfp+0, 0x1.446a1f0000bdfp+0, 0x0.0p+0, 0x0.0p+0)
    phi_f2 <- matrix(phi2_0, ncol = 1L, dimnames = list(NULL, "f2_link"))
    V_f2 <- matrix(sum((phi_f2 - mean(phi_f2))^2) / 10.0, 1L, 1L, dimnames = list("f2_link", "f2_link"))
    V_f2_star <- sum((phi2_star - mean(phi2_star))^2) / 10.0

    pan2_rows <- vector("list", n)
    for (i in seq_len(n)) {
      pan2_rows[[i]] <- data.frame(
        id = rep(id[i], 4L),
        time = c(1.0, 2.0, 3.0, 4.0),
        y = c(if (po[i]) c(yL2[i], yR2[i]) else rep(NA_real_, 2L), 0.0, 0.0),
        g = rep(G[i], 4L),
        stringsAsFactors = FALSE
      )
    }
    pan2_df <- do.call(rbind, pan2_rows)
    po2_owner <- ns$construct_ipw_panel_owner(pan2_df, "id", "time", "y", "g", o_ctrl)
    fo2_owner <- ns$construct_ipw_frame_owner(po2_owner, id = id, Z = Z2, weights = w2, control = o_ctrl)
    nv2 <- ns$construct_ipw_nuisance_view(fo2_owner, g = 2.0, K = c(3.0, 4.0), s = 1.0, t = 2.0, control = o_ctrl)
    lv2 <- ns$construct_ipw_link_view(fo2_owner, s = 1.0, t = 2.0, g = 2.0, K = c(3.0, 4.0), nuisance_view = nv2, control = o_ctrl, link_id = "f2_link")
    pr2_v <- ns$construct_ipw_enclosure_problem_view(nv2, beta0_2, F2, r_ball_2, cc, pc, o_ctrl)
    lk2_v <- ns$enclose_ipw_link_view(pr2_v, lv2, lb2, lc, o_ctrl)

    pr2_d <- ns$construct_certified_ipw_nuisance_problem(id, Z2, w2, G, po, cohorts, beta0_2, F2, r_ball_2, cc, pc)
    lk2_d <- ns$enclose_ipw_link(pr2_d, id, yL2, yR2, lb2, lc, link_id = "f2_link")

    cand_f2 <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = "f2_link", phi = phi_f2[idx, , drop = FALSE], a_index = idx)
      },
      link_ids = "f2_link"
    )

    list(
      n = n,
      id = id,
      o_ctrl = o_ctrl,
      cov_ctrl = cov_ctrl,
      op_10 = op_10,
      fixture1 = list(
        panel_owner = po1_owner,
        frame_owner = fo1_owner,
        frame_dense = ns$ipw_covariance_frame_spec(id = id, normalizer_n = n, embedding = "identity", w = w1, G = G),
        links_view = list(link1 = lk1_1_v, link2 = lk1_2_v),
        links_dense = list(link1 = lk1_1_d, link2 = lk1_2_d),
        candidate_rows = cand_10,
        phi_cand = phi_10,
        V_cand = V_10
      ),
      fixture2 = list(
        panel_owner = po2_owner,
        frame_owner = fo2_owner,
        frame_dense = ns$ipw_covariance_frame_spec(id = id, normalizer_n = n, embedding = "identity", w = w2, G = G),
        links_view = list(f2_link = lk2_v),
        links_dense = list(f2_link = lk2_d),
        candidate_rows = cand_f2,
        phi_cand = phi_f2,
        V_cand = V_f2,
        V_star = V_f2_star
      )
    )
  },

  build_coherent24_sixlink_view_and_dense = function(ns, chunk_rows = 5L, operator_block_rows = 4L) {
    # Self-contained frozen coherent24 24-unit, 5-period, 6-link panel and exact prespecified
    # analytical witness (beta0 = 0, F = 0.5, r = 1e-6; NOT a sealed actual-glm result).
    n <- 24L
    id <- seq_len(n)
    G <- rep(c(3, 3, Inf, Inf), times = 6L)
    signs <- rep(c(-1, 1, -1, 1), times = 6L)
    patterns <- list(
      c(1, 2),
      c(2, 3),
      c(3, 4),
      c(4, 5),
      c(1, 4),
      c(1, 2, 4)
    )
    amps <- list(
      c("1" = 0.0, "2" = 0.25),
      c("2" = 0.0, "3" = 0.50),
      c("3" = 0.0, "4" = 0.75),
      c("4" = 0.0, "5" = 1.00),
      c("1" = 0.0, "4" = 1.25),
      c("1" = 0.0, "2" = 0.75, "4" = -0.50)
    )
    treated_offset <- c(1.5, 0.0, 1.0, 1.0, 1.0)

    pan_rows <- list()
    Y_mat <- matrix(NA_real_, nrow = n, ncol = 5L)
    for (i in seq_len(n)) {
      pat_idx <- ((i - 1L) %/% 4L) + 1L
      obs_times <- patterns[[pat_idx]]
      amp_map <- amps[[pat_idx]]
      for (t_val in 1:5) {
        if (t_val %in% obs_times) {
          y0 <- i / 8.0 + t_val / 2.0 + signs[i] * amp_map[[as.character(t_val)]]
          y_it <- if (is.finite(G[i])) y0 + treated_offset[t_val] else y0
          Y_mat[i, t_val] <- y_it
        }
      }
      pan_rows[[i]] <- data.frame(
        id = rep(as.double(i), 5L),
        time = as.double(1:5),
        y = Y_mat[i, ],
        g = rep(G[i], 5L),
        stringsAsFactors = FALSE
      )
    }
    pan_df <- do.call(rbind, pan_rows)
    Z_a <- matrix(1.0, nrow = n, ncol = 1L, dimnames = list(NULL, "intercept"))
    w_a <- rep(1.0, n)

    o_ctrl <- ns$ipw_owner_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_retained_entries = 50000L,
      max_matrix_entries = 50000L,
      max_work_ops = 500000L,
      max_output_entries = 50000L
    )
    panel_owner <- ns$construct_ipw_panel_owner(
      data = pan_df, idname = "id", tname = "time", yname = "y", gname = "g", control = o_ctrl
    )
    frame_owner <- ns$construct_ipw_frame_owner(
      panel_owner = panel_owner, id = as.double(id), Z = Z_a, weights = w_a, control = o_ctrl
    )

    cert_ctrl <- ns$logit_certificate_control(
      taylor_degree = 18L, max_abs_eta = 40.0, max_squarings = 16L,
      max_chunk_rows = as.integer(chunk_rows), max_matrix_entries = 5000L, max_work_ops = 500000L
    )
    prob_ctrl <- ns$ipw_nuisance_problem_control(
      max_chunk_rows = as.integer(chunk_rows), max_retained_entries = 10000L,
      max_matrix_entries = 10000L, max_work_ops = 500000L
    )
    link_ctrl <- ns$ipw_link_enclosure_control(
      max_chunk_rows = as.integer(chunk_rows), max_retained_entries = 10000L,
      max_matrix_entries = 10000L, max_work_ops = 500000L, max_output_entries = 10000L
    )
    cov_ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_operator_block_rows = as.integer(operator_block_rows),
      max_cluster_rows = 24L,
      max_matrix_entries = 20000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = link_ctrl
    )

    pairs <- rbind(
      c(1, 2),
      c(1, 4),
      c(2, 3),
      c(2, 4),
      c(3, 4),
      c(4, 5)
    )
    L <- nrow(pairs)
    lids <- paste0("pair_", pairs[, 1L], "_", pairs[, 2L])
    links_view <- vector("list", L)
    links_dense <- vector("list", L)

    for (l in seq_len(L)) {
      s_t <- pairs[l, 1L]
      t_t <- pairs[l, 2L]
      nv <- ns$construct_ipw_nuisance_view(
        frame_owner = frame_owner, g = 3, K = Inf, s = s_t, t = t_t, control = o_ctrl
      )
      lv <- ns$construct_ipw_link_view(
        frame_owner = frame_owner, s = s_t, t = t_t, g = 3, K = Inf,
        nuisance_view = nv, control = o_ctrl, link_id = lids[l]
      )
      pr_v <- ns$construct_ipw_enclosure_problem_view(
        nuisance_view = nv, beta0 = 0.0, F = matrix(0.5, 1L, 1L), r = 1e-6,
        certificate_control = cert_ctrl, problem_control = prob_ctrl, reader_control = o_ctrl
      )
      links_view[[l]] <- ns$enclose_ipw_link_view(
        problem = pr_v, link_view = lv, loading_beta_candidate = 0.0,
        control = link_ctrl, reader_control = o_ctrl
      )

      ys <- Y_mat[, s_t]
      yt <- Y_mat[, t_t]
      po <- !is.na(ys) & !is.na(yt)
      pr_d <- ns$construct_ipw_enclosure_problem(
        id = as.double(id), Z = Z_a, w = w_a, G = G, pair_observed = po,
        cohorts = c(3, Inf), beta0 = 0.0, F = matrix(0.5, 1L, 1L), r = 1e-6,
        certificate_control = cert_ctrl, problem_control = prob_ctrl
      )
      links_dense[[l]] <- ns$enclose_ipw_link(
        problem = pr_d, id = as.double(id),
        y_left = ifelse(po, ys, 0.0), y_right = ifelse(po, yt, 0.0),
        loading_beta_candidate = 0.0, control = link_ctrl, link_id = lids[l]
      )
    }

    phi_exact <- matrix(0.0, nrow = n, ncol = L, dimnames = list(NULL, lids))
    phi_exact[1:4, 1L]   <- c(-1.5,  1.5,  1.5, -1.5)
    phi_exact[5:8, 3L]   <- c(-6.0,  6.0,  6.0, -6.0)
    phi_exact[9:12, 5L]  <- c(-9.0,  9.0,  9.0, -9.0)
    phi_exact[13:16, 6L] <- c(-12.0, 12.0, 12.0, -12.0)
    phi_exact[17:20, 2L] <- c(-7.5,  7.5,  7.5, -7.5)
    phi_exact[21:24, 1L] <- c(-4.5,  4.5,  4.5, -4.5)
    phi_exact[21:24, 2L] <- c( 3.0, -3.0, -3.0,  3.0)
    phi_exact[21:24, 4L] <- c(15.0, -15.0, -15.0, 15.0)

    Omega_exact <- matrix(
      c(
         15/4,  -9/4,  0, -45/4,    0,  0,
         -9/4,  87/8,  0,  15/2,    0,  0,
            0,     0,  6,     0,    0,  0,
        -45/4,  15/2,  0,  75/2,    0,  0,
            0,     0,  0,     0, 27/2,  0,
            0,     0,  0,     0,    0, 24
      ),
      nrow = 6L, ncol = 6L, byrow = TRUE,
      dimnames = list(lids, lids)
    )

    C_target <- matrix(
      c(
        -200/1313,  200/1313, 1109/1313, 4/1313, -204/1313, 0,
         -50/101,    50/101,    50/101,  1/101,    50/101,  0,
         -50/101,    50/101,    50/101,  1/101,    50/101,  1
      ),
      nrow = 3L, ncol = 6L, byrow = TRUE,
      dimnames = list(c("att_3_2", "att_3_3", "att_3_4"), lids)
    )
    target_cov_exact <- matrix(
      c(
        1109/5252,  25/202,   25/202,
          25/202,  325/808,  325/808,
          25/202,  325/808, 1133/808
      ),
      nrow = 3L, ncol = 3L, byrow = TRUE,
      dimnames = list(rownames(C_target), rownames(C_target))
    )

    callback <- function(a_start, a_end) {
      idx <- seq.int(a_start, a_end)
      list(
        phi = phi_exact[idx, , drop = FALSE],
        id = as.double(id[idx]),
        link_ids = lids,
        a_index = idx
      )
    }

    list(
      n = n,
      L = L,
      id = as.double(id),
      lids = lids,
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      links_view = links_view,
      links_dense = links_dense,
      frame_dense = ns$ipw_covariance_frame_spec(id = as.double(id), normalizer_n = n, embedding = "identity"),
      phi_exact = phi_exact,
      Omega_exact = Omega_exact,
      C_id = diag(1.0, L),
      C_target = C_target,
      target_cov_exact = target_cov_exact,
      callback = callback,
      o_ctrl = o_ctrl,
      cov_ctrl = cov_ctrl,
      analytical_witness_note = "beta0 = 0 is a prespecified analytical witness, not a sealed actual-glm result."
    )
  },

  build_weighted96_shared_mask_view = function(ns, chunk_rows = 24L, operator_block_rows = 16L) {
    # Frozen 96-unit shared-mask 3-link triangle from helper-cp02-shared-audit.R and
    # cp02_shared_ipw_audit_newton_reference.R (no runtime glm/Newton fitting).
    n <- 96L
    id <- as.double(seq_len(n))
    G <- rep(c(4, 7, Inf, 9), each = 24L)
    x <- rep(seq(-1.4, 1.4, length.out = 24L), 4L)
    Z <- cbind(one = rep(1.0, n), x = x)
    w <- 0.7 + ((seq_len(n) * 13L) %% 29L) / 20.0
    obs <- matrix(FALSE, nrow = n, ncol = 8L)
    for (i in seq_len(n)) {
      k <- ((i - 1L) %% 24L) + 1L
      if (G[i] != 9 && k <= 18L) {
        obs[i, 1:6] <- TRUE
      } else {
        obs[i, 7:8] <- TRUE
        obs[i, 1L + (k %% 6L)] <- TRUE
      }
    }
    Y <- outer(seq_len(n) / 17.0, 1:8, function(a, b) sin(a * b) + a + b^2 / 13.0) +
      outer(x, 1:8, function(a, b) a * b / 3.0)
    pan <- data.frame(
      id = rep(id, each = 8L),
      time = rep(as.double(1:8), n),
      y = as.vector(t(ifelse(obs, Y, NA_real_))),
      g = rep(G, each = 8L),
      stringsAsFactors = FALSE
    )

    beta0 <- c(-0.6796476143628154, -0.098138014861467074)
    F_mat <- matrix(c(0.48787361489041736, 0.0, -0.18791813311566585, 0.31157384415088302), 2L, 2L)
    lb12 <- c(-1.1361384932262233, -2.9496461278979176)
    lb23 <- c(-0.46554179477970431, -1.208640989222411)
    lb13 <- c(-1.6016802880059269, -4.1582871171203282)

    o_ctrl <- ns$ipw_owner_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_retained_entries = 250000L,
      max_matrix_entries = 250000L,
      max_work_ops = 2000000L,
      max_output_entries = 250000L
    )
    panel_owner <- ns$construct_ipw_panel_owner(
      data = pan, idname = "id", tname = "time", yname = "y", gname = "g", control = o_ctrl
    )
    frame_owner <- ns$construct_ipw_frame_owner(
      panel_owner = panel_owner, id = id, Z = Z, weights = w, control = o_ctrl
    )

    nv_sh <- ns$construct_ipw_nuisance_view(
      frame_owner = frame_owner, g = 4, K = c(7, Inf), s = 1, t = 2, control = o_ctrl
    )
    lv_12 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 1, t = 2, g = 4, K = c(7, Inf),
      nuisance_view = nv_sh, control = o_ctrl, link_id = "l12"
    )
    lv_23 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 2, t = 3, g = 4, K = c(7, Inf),
      nuisance_view = nv_sh, control = o_ctrl, link_id = "l23"
    )
    lv_13 <- ns$construct_ipw_link_view(
      frame_owner = frame_owner, s = 1, t = 3, g = 4, K = c(7, Inf),
      nuisance_view = nv_sh, control = o_ctrl, link_id = "l13"
    )

    cert_ctrl <- ns$logit_certificate_control(
      taylor_degree = 18L, max_abs_eta = 40.0, max_squarings = 16L,
      max_chunk_rows = as.integer(chunk_rows), max_matrix_entries = 50000L, max_work_ops = 1000000L
    )
    prob_ctrl <- ns$ipw_nuisance_problem_control(
      max_chunk_rows = as.integer(chunk_rows), max_retained_entries = 50000L,
      max_matrix_entries = 50000L, max_work_ops = 1000000L
    )
    link_ctrl <- ns$ipw_link_enclosure_control(
      max_chunk_rows = as.integer(chunk_rows), max_retained_entries = 50000L,
      max_matrix_entries = 50000L, max_work_ops = 1000000L, max_output_entries = 50000L
    )
    cov_ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = as.integer(chunk_rows),
      max_operator_block_rows = as.integer(operator_block_rows),
      max_cluster_rows = 96L,
      max_matrix_entries = 100000L,
      max_work_ops = 2000000L,
      max_output_entries = 50000L,
      link_control = link_ctrl
    )

    pr_sh <- ns$construct_ipw_enclosure_problem_view(
      nuisance_view = nv_sh, beta0 = beta0, F = F_mat, r = 1e-6,
      certificate_control = cert_ctrl, problem_control = prob_ctrl, reader_control = o_ctrl
    )
    lk_12 <- ns$enclose_ipw_link_view(pr_sh, lv_12, lb12, link_ctrl, o_ctrl)
    lk_23 <- ns$enclose_ipw_link_view(pr_sh, lv_23, lb23, link_ctrl, o_ctrl)
    lk_13 <- ns$enclose_ipw_link_view(pr_sh, lv_13, lb13, link_ctrl, o_ctrl)

    eval_point_if <- function(s, t, lb) {
      K <- c(4, 7, Inf)
      I <- outer(G, K, "==") * 1.0
      S <- as.numeric(obs[, s] & obs[, t] & rowSums(I) > 0)
      D <- I[, 1L]
      C_ind <- rowSums(I) - D
      q <- colSums(I * w) / n
      r_vec <- colSums(I * w * S) / colSums(I * w)
      rG <- drop(I %*% r_vec)
      rG[rG == 0] <- 1.0
      pw <- w * S / rG
      p <- 1.0 / (1.0 + exp(-drop(Z %*% beta0)))
      odds <- exp(drop(Z %*% beta0))
      L_mat <- vapply(
        seq_along(K),
        function(j) -colSums(Z * (w * I[, j] * S * (D - p))) / n / r_vec[j]^2,
        numeric(ncol(Z))
      )
      delta <- Y[, t] - Y[, s]
      a <- w * D * S / r_vec[1L]
      b <- w * C_ind * S * odds / rG
      T_val <- sum(a * delta) / sum(a)
      M_val <- sum(b * delta) / sum(b)
      Gr <- c(0.0, vapply(2:3, function(j) sum(w * I[, j] * S * odds * (delta - M_val)) / r_vec[j]^2 / sum(b), numeric(1L)))
      lr <- Gr + drop(crossprod(L_mat, lb))
      xir <- I * w * (S - rep(r_vec, each = n)) / rep(q, each = n)
      a / mean(a) * (delta - T_val) - b / mean(b) * (delta - M_val) + pw * (D - p) * drop(Z %*% lb) + drop(xir %*% lr)
    }

    phi_mat <- cbind(
      l12 = eval_point_if(1L, 2L, lb12),
      l23 = eval_point_if(2L, 3L, lb23),
      l13 = eval_point_if(1L, 3L, lb13)
    )
    V_3 <- crossprod(phi_mat) / n
    V_3 <- 0.5 * (V_3 + t(V_3))
    dimnames(V_3) <- list(c("l12", "l23", "l13"), c("l12", "l23", "l13"))
    V_2 <- V_3[1:2, 1:2, drop = FALSE]

    list(
      n = n,
      id = id,
      panel_owner = panel_owner,
      frame_owner = frame_owner,
      pr_sh = pr_sh,
      links_3 = list(lk_12, lk_23, lk_13),
      links_2 = list(lk_12, lk_23),
      phi_mat = phi_mat,
      V_3 = V_3,
      V_2 = V_2,
      o_ctrl = o_ctrl,
      cov_ctrl = cov_ctrl
    )
  },

  # Resource-only adversaries. These clone sealed identities while sharing the
  # existing numerical arrays. They are not additional fitted/reference cases.
  resource_links = function(ns, b14, ids, distinct = FALSE) {
    out <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      pr <- b14$prob12_view
      nv <- b14$nv_12
      if (distinct) {
        nv <- ns$ipw_owner_seal(nv, "nuisance_view")
        pr$nuisance_ref <- nv$ref
        pr <- ns$ipw_enclosure_seal(pr, "nuisance_problem")
      }
      lv <- b14$lv_1
      lv$link_id <- ids[i]
      lv$nuisance_ref <- nv$ref
      lv <- ns$ipw_owner_seal(lv, "link_view")
      lk <- b14$links_view$l1
      lk$link_id <- ids[i]
      lk$problem_ref <- pr$ref
      lk$nuisance_ref <- nv$ref
      lk$link_ref <- lv$ref
      out[[i]] <- ns$ipw_enclosure_seal(lk, "link_enclosure")
    }
    names(out) <- ids
    out
  },
  invalid_index_cases = list(
    fractional = function(idx) as.double(idx) + 0.25,
    negative_fraction = function(idx) as.double(idx) - 0.25,
    overflow = function(idx) rep(2^40, length(idx)),
    nonfinite = function(idx) replace(as.double(idx), 1L, Inf),
    shifted = function(idx) idx + 1L
  )
)
