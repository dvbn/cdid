# Self-contained base-R fixtures and 4bd sealed-object builders for CP-02 IPW covariance transfer tests.
list(
  operator_fixtures = list(
    list(
      name = "complete_unit_GMM",
      ids = c("unit_0", "unit_1", "unit_2", "unit_3", "unit_4", "unit_5", "unit_6", "unit_7", "unit_8", "unit_9", "unit_10", "unit_11", "unit_12", "unit_13"),
      labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
      cluster_order = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
      omega = c(0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0),
      kappa = 0x1.2492492492492p-4,
      centered = TRUE,
      purpose = "GMM",
      normalizer = 0x1.c000000000000p+3,
      phi = structure(c(-0x1.5d55b01e7f943p+2, 0x1.d80cfac2a3c97p+1, 0x1.c53ccaf4b6bdfp+0, 0x0.0p+0, -0x1.f8f4972b29972p-4, 0x1.b9c83d1a5a676p+2, 0x1.277f23965329bp+3, 0x1.c78132a108ecap+1, -0x1.05fcdcde5326ep+0, -0x1.a1da954faab28p+2, -0x1.4202fb96db4d2p+3, -0x1.f50e1e1789d11p+0, 0x0.0p+0, 0x0.0p+0, 0x1.1941a6b8acb4bp+3, -0x1.ed73e711c2195p+2, -0x1.143d997e5d402p+0, 0x0.0p+0, 0x1.111061f37112bp+1, -0x1.91ca704c3626bp+1, -0x1.7505e8b409e4bp+2, -0x1.84c9b8c77cd42p+0, -0x1.a4da13bd5c0e9p+0, 0x1.20862c6253c16p+2, 0x1.29d0606a4755ap+2, 0x1.abaab1a83c82fp-1, 0x0.0p+0, 0x0.0p+0, 0x1.f07c1f07c1f08p-1, -0x1.1111111111111p+1, 0x1.29e4129e4129ep+0, 0x0.0p+0, 0x0.0p+0, 0x1.04a7904a7904ap+1, -0x1.04a7904a7904ap+1, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      C = structure(c(0x1.5555555555555p-2, 0x1.999999999999ap-4, -0x1.6666666666666p-1, 0x1.999999999999ap-3, 0x1.199999999999ap+0, -0x1.3333333333333p-2), dim = c(2L, 3L)),
      E_eval = structure(c(0x1.1dd55b01e7f95p-54, 0x1.f7a5fe9051482p-54, 0x1.d5094cdb8d7d5p-55, 0x0.0p+0, 0x1.2e67458804d7fp-58, 0x1.db5f4258e22b8p-55, 0x1.a806483d7dc61p-53, 0x1.50d0b8e415cbcp-59, 0x1.66f16ef3fb668p-56, 0x1.2766e0a0c1532p-54, 0x1.a6d7805f5c430p-52, 0x1.93cf9cd27ad35p-54, 0x0.0p+0, 0x0.0p+0, 0x1.44d7cb28ea697p-54, 0x1.a94f6935e0a20p-52, 0x1.5ea8fff930e9bp-58, 0x0.0p+0, 0x1.8adc237cd58f7p-54, 0x1.2f21f7b22e0b7p-53, 0x1.78eb05cedc0dcp-54, 0x1.786f4fd8e082bp-54, 0x1.6b17fd8c6868ap-55, 0x1.0d4e7ebdb1fdap-53, 0x1.54e533e9849fdp-52, 0x1.9e140b082a296p-55, 0x0.0p+0, 0x0.0p+0, 0x1.f07c1f07c1f08p-56, 0x1.1111111111112p-55, 0x1.04a7904a7904bp-54, 0x0.0p+0, 0x0.0p+0, 0x1.e4129e4129e42p-53, 0x1.e4129e4129e42p-53, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      E_root = structure(c(0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      V_candidate = structure(c(0x1.004292f3cee41p+4, -0x1.e2cfed42e3ab3p-1, -0x1.e2cfed42e3ab3p-1, 0x1.3a4df88b6f8ccp-3), dim = c(2L, 2L)),
      ref_exact_V = structure(c(16.016253425956396, -0.94299260560531761, -0.94299260560531761, 0.15346903014682878), dim = c(2L, 2L)),
      ref_epsilon_acc = 1.0752303651665084e-14,
      ref_nu_eval = 1.9496229819729022e-16,
      ref_nu_root = 2.13657644502059e-05,
      ref_epsilon_V = 0.00017131291479544284
    ),
    list(
      name = "complete_cluster_inference",
      ids = c("unit_0", "unit_1", "unit_2", "unit_3", "unit_4", "unit_5", "unit_6", "unit_7", "unit_8", "unit_9", "unit_10", "unit_11", "unit_12", "unit_13"),
      labels = c("b", "d", "a", "c", "b", "d", "a", "c", "b", "d", "a", "c", "zero", "zero"),
      cluster_order = c("c", "zero", "a", "d", "b"),
      omega = c(0x1.999999999999ap-4, 0x1.6666666666666p+0, 0x1.5555555555555p-2, 0x1.ccccccccccccdp-1, 0x1.b333333333333p+0),
      kappa = 0x1.2492492492492p-4,
      centered = TRUE,
      purpose = "inference",
      normalizer = 0x1.c000000000000p+3,
      phi = structure(c(-0x1.5d55b01e7f943p+2, 0x1.d80cfac2a3c97p+1, 0x1.c53ccaf4b6bdfp+0, 0x0.0p+0, -0x1.f8f4972b29972p-4, 0x1.b9c83d1a5a676p+2, 0x1.277f23965329bp+3, 0x1.c78132a108ecap+1, -0x1.05fcdcde5326ep+0, -0x1.a1da954faab28p+2, -0x1.4202fb96db4d2p+3, -0x1.f50e1e1789d11p+0, 0x0.0p+0, 0x0.0p+0, 0x1.1941a6b8acb4bp+3, -0x1.ed73e711c2195p+2, -0x1.143d997e5d402p+0, 0x0.0p+0, 0x1.111061f37112bp+1, -0x1.91ca704c3626bp+1, -0x1.7505e8b409e4bp+2, -0x1.84c9b8c77cd42p+0, -0x1.a4da13bd5c0e9p+0, 0x1.20862c6253c16p+2, 0x1.29d0606a4755ap+2, 0x1.abaab1a83c82fp-1, 0x0.0p+0, 0x0.0p+0, 0x1.f07c1f07c1f08p-1, -0x1.1111111111111p+1, 0x1.29e4129e4129ep+0, 0x0.0p+0, 0x0.0p+0, 0x1.04a7904a7904ap+1, -0x1.04a7904a7904ap+1, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      C = structure(c(0x1.5555555555555p-2, 0x1.999999999999ap-4, -0x1.6666666666666p-1, 0x1.999999999999ap-3, 0x1.199999999999ap+0, -0x1.3333333333333p-2), dim = c(2L, 3L)),
      E_eval = structure(c(0x1.1dd55b01e7f95p-54, 0x1.f7a5fe9051482p-54, 0x1.d5094cdb8d7d5p-55, 0x0.0p+0, 0x1.2e67458804d7fp-58, 0x1.db5f4258e22b8p-55, 0x1.a806483d7dc61p-53, 0x1.50d0b8e415cbcp-59, 0x1.66f16ef3fb668p-56, 0x1.2766e0a0c1532p-54, 0x1.a6d7805f5c430p-52, 0x1.93cf9cd27ad35p-54, 0x0.0p+0, 0x0.0p+0, 0x1.44d7cb28ea697p-54, 0x1.a94f6935e0a20p-52, 0x1.5ea8fff930e9bp-58, 0x0.0p+0, 0x1.8adc237cd58f7p-54, 0x1.2f21f7b22e0b7p-53, 0x1.78eb05cedc0dcp-54, 0x1.786f4fd8e082bp-54, 0x1.6b17fd8c6868ap-55, 0x1.0d4e7ebdb1fdap-53, 0x1.54e533e9849fdp-52, 0x1.9e140b082a296p-55, 0x0.0p+0, 0x0.0p+0, 0x1.f07c1f07c1f08p-56, 0x1.1111111111112p-55, 0x1.04a7904a7904bp-54, 0x0.0p+0, 0x0.0p+0, 0x1.e4129e4129e42p-53, 0x1.e4129e4129e42p-53, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      E_root = structure(c(0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(14L, 3L)),
      V_candidate = structure(c(0x1.091137650c965p+3, -0x1.0da591370f63ep+0, -0x1.0da591370f63ep+0, 0x1.142d7d0b80e55p-3), dim = c(2L, 2L)),
      ref_exact_V = structure(c(8.2833516095749449, -1.0533076056884876, -1.0533076056884876, 0.13485238734584706), dim = c(2L, 2L)),
      ref_epsilon_acc = 1.443249899377295e-14,
      ref_nu_eval = 2.4741642057128305e-16,
      ref_nu_root = 3.1480092221088686e-05,
      ref_epsilon_V = 0.00018266479912705579
    ),
    list(
      name = "cluster_rounding_cancellation",
      ids = c("unit_0", "unit_1", "unit_2", "unit_3"),
      labels = c("x", "x", "x", "z"),
      cluster_order = c("z", "x"),
      omega = c(0x1.5555555555555p-1, 0x1.5555555555555p-2),
      kappa = 0x1.2492492492492p-3,
      centered = FALSE,
      purpose = "inference",
      normalizer = 0x1.0000000000000p+2,
      phi = structure(c(0x1.0000000000000p+53, 0x1.0000000000000p+0, -0x1.0000000000000p+53, -0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.5555555555555p-2, 0x1.0000000000000p+1, 0x0.0p+0), dim = c(4L, 2L)),
      C = structure(c(0x1.0000000000000p+0, 0x1.999999999999ap-4, 0x1.5555555555555p-2, -0x1.0000000000000p+0), dim = c(2L, 2L)),
      E_eval = structure(c(0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(4L, 2L)),
      E_root = structure(c(0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16), dim = c(4L, 2L)),
      V_candidate = structure(c(0x1.3ad86b27e1e66p-2, -0x1.4317bfb4dedc2p-2, -0x1.4317bfb4dedc2p-2, 0x1.fec0fd20e31f3p-2), dim = c(2L, 2L)),
      ref_exact_V = structure(c(0.30746619635508521, -0.31552028218694883, -0.31552028218694883, 0.49878306878306872), dim = c(2L, 2L)),
      ref_epsilon_acc = 1.7536887283954874,
      ref_nu_eval = 0.0,
      ref_nu_root = 1.9088900357515185e-05,
      ref_epsilon_V = 1.7537489745458417
    ),
    list(
      name = "singular_observer_Q_full_PD",
      ids = c("unit_0", "unit_1", "unit_2", "unit_3"),
      labels = c("0", "1", "2", "3"),
      cluster_order = c("0", "1", "2", "3"),
      omega = c(0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0, 0x1.0000000000000p+0),
      kappa = 0x1.0000000000000p-2,
      centered = TRUE,
      purpose = "inference",
      normalizer = 0x1.0000000000000p+2,
      phi = structure(c(0x1.0000000000000p+0, -0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p+0, -0x1.0000000000000p+0), dim = c(4L, 2L)),
      C = structure(c(0x1.0000000000000p+0, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p+0), dim = c(2L, 2L)),
      E_eval = structure(c(0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0, 0x0.0p+0), dim = c(4L, 2L)),
      E_root = structure(c(0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16, 0x1.0000000000000p-16), dim = c(4L, 2L)),
      V_candidate = structure(c(0x1.0000000000000p-1, 0x0.0p+0, 0x0.0p+0, 0x1.0000000000000p-1), dim = c(2L, 2L)),
      ref_exact_V = structure(c(0.5, 0.0, 0.0, 0.5), dim = c(2L, 2L)),
      ref_epsilon_acc = 0.0,
      ref_nu_eval = 0.0,
      ref_nu_root = 2.1579186437577746e-05,
      ref_epsilon_V = 3.629218670018871e-05
    )
  ),

  build_14unit_3link_4bd = function(ns, r_ball = 1e-6, max_chunk_rows = 5L, max_operator_block_rows = 4L) {
    id <- paste0("unit_", 0:13)
    G <- c(rep(2L, 4L), rep(5L, 4L), rep(99L, 4L), rep(9L, 2L))
    X <- c(rep(c(0, 1, 1, 0), 3L), 0, 1)
    Z <- cbind(intercept = rep(1.0, 14L), x = as.double(X))
    S <- as.logical(c(rep(c(1L, 1L, 1L, 0L), 3L), 1L, 0L))
    w <- c(1.0, 2.0, 0.5, 1.5, 2.0, 1.0, 1.5, 1.0, 1.0, 1.5, 2.0, 0.5, 1.0, 1.0)
    y <- matrix(0.0, nrow = 14L, ncol = 3L)
    for (i in 0:13) {
      for (j in 0:2) {
        y[i + 1L, j + 1L] <- (((i + 1L) * (j + 2L)) %% 17L - 5L) / 3.0 + (i * j) / 11.0
      }
    }

    calc_root <- function(cohorts) {
      M <- vapply(cohorts, function(g) sum(w[G == g]), double(1L))
      P <- vapply(cohorts, function(g) sum(w[G == g & S]), double(1L))
      r <- P / M
      idx <- which((G %in% cohorts) & S)
      ci <- w[idx] / r[match(G[idx], cohorts)]
      odds <- vapply(c(0, 1), function(xv) {
        sum(ci[X[idx] == xv & G[idx] == 2L]) / sum(ci[X[idx] == xv & G[idx] != 2L])
      }, double(1L))
      beta0 <- c(log(odds[1L]), log(odds[2L]) - log(odds[1L]))
      eta <- as.vector(Z[idx, , drop = FALSE] %*% beta0)
      p <- 1 / (1 + exp(-eta))
      H <- crossprod(Z[idx, , drop = FALSE], Z[idx, , drop = FALSE] * (ci * p * (1 - p))) / 14.0
      list(beta0 = beta0, H = H, F_mat = chol(H), ci = ci, idx = idx, odds = odds)
    }
    r12 <- calc_root(c(2L, 5L, 99L))
    r3  <- calc_root(c(2L, 5L))

    calc_lambda_beta <- function(rt, s_col, t_col) {
      idx <- rt$idx
      ci <- rt$ci
      delta <- y[idx, t_col] - y[idx, s_col]
      b <- ci * ifelse(X[idx] == 0, rt$odds[1L], rt$odds[2L])
      ctrl <- G[idx] != 2L
      C_val <- sum(b[ctrl] * delta[ctrl]) / sum(b[ctrl])
      v <- as.vector(crossprod(Z[idx[ctrl], , drop = FALSE], b[ctrl] * (delta[ctrl] - C_val))) / sum(b[ctrl])
      as.vector(solve(rt$H, -v))
    }

    cc <- ns$logit_certificate_control(18L, 40, 16L, 14L, 5000L, 500000L)
    pc <- ns$ipw_nuisance_problem_control(14L, 10000L, 5000L, 500000L)
    lc <- ns$ipw_link_enclosure_control(14L, 10000L, 5000L, 500000L, 10000L)

    prob12 <- ns$construct_certified_ipw_nuisance_problem(
      id = id, Z = Z, w = w, G = G, pair_observed = S, cohorts = c(2L, 5L, 99L),
      beta0 = r12$beta0, F = r12$F_mat, r = r_ball, certificate_control = cc, problem_control = pc
    )
    prob3 <- ns$construct_certified_ipw_nuisance_problem(
      id = id, Z = Z, w = w, G = G, pair_observed = S, cohorts = c(2L, 5L),
      beta0 = r3$beta0, F = r3$F_mat, r = r_ball, certificate_control = cc, problem_control = pc
    )

    l1 <- ns$enclose_ipw_link(prob12, id, y[, 1L], y[, 2L], calc_lambda_beta(r12, 1L, 2L), lc, link_id = "l1")
    l2 <- ns$enclose_ipw_link(prob12, id, y[, 2L], y[, 3L], calc_lambda_beta(r12, 2L, 3L), lc, link_id = "l2")
    l3 <- ns$enclose_ipw_link(prob3,  id, y[, 1L], y[, 3L], calc_lambda_beta(r3,  1L, 3L), lc, link_id = "l3")

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
    ), dim = c(14L, 3L))
    C_mat <- structure(c(
      0x1.5555555555555p-2, 0x1.999999999999ap-4,
      -0x1.6666666666666p-1, 0x1.999999999999ap-3,
      0x1.199999999999ap+0, -0x1.3333333333333p-2
    ), dim = c(2L, 3L))
    V_cand_gmm <- structure(c(
      0x1.004292f3cee41p+4, -0x1.e2cfed42e3ab3p-1,
      -0x1.e2cfed42e3ab3p-1, 0x1.3a4df88b6f8ccp-3
    ), dim = c(2L, 2L))
    V_cand_cl <- structure(c(
      0x1.091137650c965p+3, -0x1.0da591370f63ep+0,
      -0x1.0da591370f63ep+0, 0x1.142d7d0b80e55p-3
    ), dim = c(2L, 2L))
    labels_cl <- c("b", "d", "a", "c", "b", "d", "a", "c", "b", "d", "a", "c", "zero", "zero")
    order_cl  <- c("c", "zero", "a", "d", "b")
    omega_cl  <- c(0x1.999999999999ap-4, 0x1.6666666666666p+0, 0x1.5555555555555p-2, 0x1.ccccccccccccdp-1, 0x1.b333333333333p+0)
    kappa_14  <- 0x1.2492492492492p-4

    cand_spec <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(
          id = id[idx],
          link_ids = c("l1", "l2", "l3"),
          phi = phi_cand[idx, , drop = FALSE]
        )
      },
      link_ids = c("l1", "l2", "l3")
    )

    frame_spec <- ns$ipw_covariance_frame_spec(
      id = id,
      normalizer_n = 14,
      embedding = "identity",
      w = w,
      G = G
    )

    op_gmm <- ns$ipw_covariance_operator_spec(
      purpose = "GMM",
      row_level = "unit",
      omega = rep(1.0, 14L),
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
    ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = max_chunk_rows,
      max_operator_block_rows = max_operator_block_rows,
      max_cluster_rows = 14L,
      max_matrix_entries = 10000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = lc
    )

    list(
      id = id, Z = Z, w = w, G = G, S = S, y = y,
      prob12 = prob12, prob3 = prob3,
      links = list(l1 = l1, l2 = l2, l3 = l3),
      phi_cand = phi_cand,
      C = C_mat,
      V_cand_gmm = V_cand_gmm,
      V_cand_cl = V_cand_cl,
      labels_cl = labels_cl,
      order_cl = order_cl,
      omega_cl = omega_cl,
      kappa = kappa_14,
      candidate_rows = cand_spec,
      frame = frame_spec,
      op_gmm = op_gmm,
      op_cl = op_cl,
      control = ctrl,
      link_control = lc,
      problem_control = pc,
      certificate_control = cc
    )
  },

  build_singular_observer_q_4bd = function(ns, r_ball = 1e-6) {
    id <- paste0("unit_", 0:13)
    G <- c(rep(2L, 4L), rep(5L, 4L), rep(99L, 4L), rep(9L, 2L))
    X <- c(rep(c(0, 1, 1, 0), 3L), 0, 1)
    Z <- cbind(intercept = rep(1.0, 14L), x = as.double(X))
    S <- as.logical(c(rep(c(1L, 1L, 1L, 0L), 3L), 1L, 0L))
    w <- c(1.0, 2.0, 0.5, 1.5, 2.0, 1.0, 1.5, 1.0, 1.0, 1.5, 2.0, 0.5, 1.0, 1.0)
    y <- matrix(0.0, nrow = 14L, ncol = 3L)
    for (i in 0:13) {
      for (j in 0:2) {
        y[i + 1L, j + 1L] <- (((i + 1L) * (j + 2L)) %% 17L - 5L) / 3.0 + (i * j) / 11.0
      }
    }

    calc_root <- function(cohorts) {
      M <- vapply(cohorts, function(g) sum(w[G == g]), double(1L))
      P <- vapply(cohorts, function(g) sum(w[G == g & S]), double(1L))
      r <- P / M
      idx <- which((G %in% cohorts) & S)
      ci <- w[idx] / r[match(G[idx], cohorts)]
      odds <- vapply(c(0, 1), function(xv) {
        sum(ci[X[idx] == xv & G[idx] == 2L]) / sum(ci[X[idx] == xv & G[idx] != 2L])
      }, double(1L))
      beta0 <- c(log(odds[1L]), log(odds[2L]) - log(odds[1L]))
      eta <- as.vector(Z[idx, , drop = FALSE] %*% beta0)
      p <- 1 / (1 + exp(-eta))
      H <- crossprod(Z[idx, , drop = FALSE], Z[idx, , drop = FALSE] * (ci * p * (1 - p))) / 14.0
      list(beta0 = beta0, H = H, F_mat = chol(H), ci = ci, idx = idx, odds = odds)
    }
    calc_lambda_beta <- function(rt, s_col, t_col) {
      idx <- rt$idx
      ci <- rt$ci
      delta <- y[idx, t_col] - y[idx, s_col]
      b <- ci * ifelse(X[idx] == 0, rt$odds[1L], rt$odds[2L])
      ctrl <- G[idx] != 2L
      C_val <- sum(b[ctrl] * delta[ctrl]) / sum(b[ctrl])
      v <- as.vector(crossprod(Z[idx[ctrl], , drop = FALSE], b[ctrl] * (delta[ctrl] - C_val))) / sum(b[ctrl])
      as.vector(solve(rt$H, -v))
    }

    cc <- ns$logit_certificate_control(18L, 40, 16L, 14L, 5000L, 500000L)
    pc <- ns$ipw_nuisance_problem_control(14L, 10000L, 5000L, 500000L)
    lc <- ns$ipw_link_enclosure_control(14L, 10000L, 5000L, 500000L, 10000L)

    r_a <- calc_root(c(2L, 5L, 99L))
    r_b <- calc_root(c(2L, 5L, 9L))
    prob_a <- ns$construct_certified_ipw_nuisance_problem(id, Z, w, G, S, c(2L, 5L, 99L), r_a$beta0, r_a$F_mat, r_ball, cc, pc)
    prob_b <- ns$construct_certified_ipw_nuisance_problem(id, Z, w, G, S, c(2L, 5L, 9L),  r_b$beta0, r_b$F_mat, r_ball, cc, pc)

    la <- ns$enclose_ipw_link(prob_a, id, y[, 1L], y[, 2L], calc_lambda_beta(r_a, 1L, 2L), lc, link_id = "la")
    lb <- ns$enclose_ipw_link(prob_b, id, y[, 1L], y[, 3L], calc_lambda_beta(r_b, 1L, 3L), lc, link_id = "lb")

    cha <- ns$enclose_ipw_link_if_chunk(la, 1L, 14L, lc)
    chb <- ns$enclose_ipw_link_if_chunk(lb, 1L, 14L, lc)
    phi_full <- cbind(
      0.5 * (cha$phi_0[, 1L] + cha$phi_0[, 2L]),
      0.5 * (chb$phi_0[, 1L] + chb$phi_0[, 2L])
    )
    phi_obs <- phi_full
    phi_obs[!S, ] <- 0.0

    cl_id <- ifelse(S, "c_obs", paste0("c_nonobs_", G))
    cl_ord <- c("c_obs", "c_nonobs_2", "c_nonobs_5", "c_nonobs_99", "c_nonobs_9")
    omega <- rep(1.0, 5L)
    kappa <- 1.0 / 14.0
    C_eye <- diag(2L)

    W_full <- t(vapply(cl_ord, function(cl) colSums(phi_full[cl_id == cl, , drop = FALSE]), double(2L)))
    W_obs  <- t(vapply(cl_ord, function(cl) colSums(phi_obs[cl_id == cl, , drop = FALSE]), double(2L)))
    V_full <- (crossprod(W_full) - tcrossprod(colSums(W_full)) / 5.0) * kappa
    Q_obs  <- (crossprod(W_obs)  - tcrossprod(colSums(W_obs))  / 5.0) * kappa
    V_full <- 0.5 * (V_full + t(V_full))
    Q_obs  <- 0.5 * (Q_obs + t(Q_obs))

    cand_full <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = c("la", "lb"), phi = phi_full[idx, , drop = FALSE])
      },
      link_ids = c("la", "lb")
    )
    cand_obs_dropped <- ns$ipw_candidate_rows_spec(
      callback = function(a_start, a_end) {
        idx <- seq.int(a_start, a_end)
        list(id = id[idx], link_ids = c("la", "lb"), phi = phi_obs[idx, , drop = FALSE])
      },
      link_ids = c("la", "lb")
    )

    frame_spec <- ns$ipw_covariance_frame_spec(id = id, normalizer_n = 14, embedding = "identity", w = w, G = G)
    op_cl <- ns$ipw_covariance_operator_spec(
      purpose = "inference",
      row_level = "cluster",
      omega = omega,
      centered = TRUE,
      kappa = kappa,
      cluster_id = cl_id,
      cluster_order = cl_ord
    )
    ctrl <- ns$ipw_covariance_transfer_control(
      max_chunk_rows = 5L,
      max_operator_block_rows = 3L,
      max_cluster_rows = 5L,
      max_matrix_entries = 10000L,
      max_work_ops = 1000000L,
      max_output_entries = 10000L,
      link_control = lc
    )

    list(
      links = list(la = la, lb = lb),
      phi_full = phi_full,
      phi_obs = phi_obs,
      C = C_eye,
      V_full = V_full,
      Q_obs = Q_obs,
      candidate_rows_full = cand_full,
      candidate_rows_obs_dropped = cand_obs_dropped,
      frame = frame_spec,
      op = op_cl,
      control = ctrl
    )
  }
)
