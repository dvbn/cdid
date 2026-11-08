# Package-safe testthat suite for R/static_reconstructed_or.R.
# Resolves production functions from the enclosing test environment / namespace
# (`environment(static_reconstructed_or_fit)`) without sourcing plaintext R/
# files, and loads the independent reference via `testthat::test_path()`.

prod_env <- environment(static_reconstructed_or_fit)
ref_env <- new.env(parent = baseenv())
# Allow base + stats in ref_env without exposing any production symbol
parent.env(ref_env) <- asNamespace("stats")
sys.source(
  testthat::test_path("reference", "static_reconstructed_or_reference.R"),
  envir = ref_env
)

complete_ref_env <- new.env(parent = asNamespace("stats"))
sys.source(
  testthat::test_path("reference", "static_complete_reference.R"),
  envir = complete_ref_env
)

test_that("independent reference calls zero production symbols and LS flows satisfy lem:flows", {
  prod_names <- ls(prod_env, all.names = TRUE)
  ref_names <- ls(ref_env, all.names = TRUE)
  for (rnm in ref_names) {
    fn <- get(rnm, envir = ref_env, inherits = FALSE)
    if (is.function(fn)) {
      used_globals <- codetools::findGlobals(fn, merge = FALSE)$functions
      expect_length(intersect(used_globals, prod_names), 0L)
    }
  }

  # Overidentified diamond + bridge graph for cell c = (t = 3, baseline = 2)
  edges <- rbind(
    c(1, 2),
    c(1, 3),
    c(2, 3),
    c(2, 4),
    c(3, 4)
  )
  fl <- static_reconstructed_or_ls_flow(edges = edges, t = 3, rank_tol = 1e-10, cell_id = "c_flow")
  ref_fl <- ref_env$ref_recon_or_incidence(edges = edges, t = 3, baseline_date = 2)

  expect_equal(fl$D, ref_fl$D, tolerance = 1e-14)
  expect_equal(fl$ell, ref_fl$ell, tolerance = 1e-14)
  expect_equal(fl$varpi, ref_fl$varpi, tolerance = 1e-14)
  expect_equal(drop(crossprod(fl$D, fl$varpi)), fl$ell, tolerance = 1e-14)

  # Incompatible flow check: direct edge (2,3) and 2-step path (1,2),(1,3)
  # LS flow weights are (2/3) on (2,3) and (-1/3, +1/3) on (1,2), (1,3).
  edges_inc <- rbind(c(1, 2), c(1, 3), c(2, 3))
  fl_inc <- static_reconstructed_or_ls_flow(edges = edges_inc, t = 3, rank_tol = 1e-10)
  expect_equal(unname(fl_inc$varpi), c(-1 / 3, 1 / 3, 2 / 3), tolerance = 1e-14)
  psi_inc <- c(0.0, 1.5, 1.0)
  expect_equal(sum(fl_inc$varpi * psi_inc), 7 / 6, tolerance = 1e-14)

  # Flow telescoping identity: for any date vector y, sum_k varpi_k (y_{u_k} - y_{s_k}) == y_t - y_{t-1}
  y_dates <- c("1" = -3.7, "2" = 1.25, "3" = 4.85, "4" = 9.1)
  dy_edges <- y_dates[as.character(edges[, 2L])] - y_dates[as.character(edges[, 1L])]
  expect_equal(sum(fl$varpi * dy_edges), unname(y_dates["3"] - y_dates["2"]), tolerance = 1e-14)
})

test_that("D-only factor blocks do not claim original graph component metadata", {
  # One connected star loses its shared baseline column when grounded, leaving
  # two independent factor blocks. D alone cannot distinguish that star from
  # two separately grounded graph components.
  star_edges <- rbind(c(1, 2), c(1, 3))
  star_graph <- static_reconstructed_or_graph(edges = star_edges, t = 2)
  D_graph <- static_reconstructed_or_graph(D = star_graph$D, ell = star_graph$ell, t = 2)

  expect_identical(unname(star_graph$D), diag(2L))
  expect_identical(unname(star_graph$ell), c(1.0, 0.0))
  expect_identical(star_graph$n_components, 1L)
  expect_identical(unname(star_graph$edge_component_id), c(1L, 1L))
  expect_identical(D_graph$n_components, NA_integer_)
  expect_null(D_graph$component_id)
  expect_null(D_graph$edge_component_id)
  expect_identical(D_graph$factor_col_block_id, c(1L, 2L))
  expect_identical(D_graph$factor_edge_block_id, c(1L, 2L))

  star_flow <- static_reconstructed_or_ls_flow(edges = star_edges, t = 2, rank_tol = 1e-10)
  D_flow <- static_reconstructed_or_ls_flow(D = star_graph$D, ell = star_graph$ell,
                                           t = 2, rank_tol = 1e-10)
  expect_identical(unname(star_flow$varpi), c(1.0, 0.0))
  expect_identical(D_flow$varpi, star_flow$varpi)
  expect_identical(D_flow$n_components, NA_integer_)
  expect_null(D_flow$edge_component_id)
  expect_identical(D_flow$rank_tol, star_flow$rank_tol)
})

test_that("complete-data baseline route reduces to static_complete_core across separate, pooled, intermediate, and p=0 partitions", {
  n <- 30L
  uid <- sprintf("u%02d", seq_len(n))
  idx <- seq_len(n)

  # Four dates t = 1, 2, 3, 4; all outcomes observed (O = 1)
  O_all <- rep(1L, n)
  Z_ab <- cbind(
    z1 = sin(0.37 * idx) + 0.1 * (idx / n),
    z2 = cos(0.53 * idx) - 0.05 * ((idx %% 5L) - 2)
  )
  Z_c <- cbind(
    w1 = 0.4 * sin(0.29 * idx) + 0.3 * cos(0.71 * idx)
  )

  Y1 <- 0.5 + 0.4 * Z_ab[, 1L] - 0.2 * Z_ab[, 2L] + 0.15 * sin(idx)
  Y2 <- 1.1 + 0.9 * Z_ab[, 1L] + 0.3 * Z_ab[, 2L] - 0.20 * cos(idx)
  Y3 <- 2.0 + 1.7 * Z_ab[, 1L] - 0.1 * Z_ab[, 2L] + 0.25 * sin(1.7 * idx)
  Y4 <- 2.8 + 2.1 * Z_ab[, 1L] + 0.5 * Z_ab[, 2L] - 0.18 * cos(1.3 * idx)
  Ymat <- cbind("1" = Y1, "2" = Y2, "3" = Y3, "4" = Y4)

  # Three cells with overlapping units:
  # c1: t = 2, d = 0; c2: t = 3, d = 0; c3: t = 4, d = 1
  G1 <- as.integer(idx %in% c(1:6))
  C1 <- as.integer(idx %in% c(7:22))
  R1 <- G1 + C1

  G2 <- as.integer(idx %in% c(7:12))
  C2 <- as.integer(idx %in% c(13:27))
  R2 <- G2 + C2

  G3 <- as.integer(idx %in% c(3:8))
  C3 <- as.integer(idx %in% c(9:18, 23:28))
  R3 <- G3 + C3

  # Differing, overidentified role graphs for stayers and switchers
  edges_C1 <- rbind(c(1, 2), c(1, 3), c(2, 3), c(1, 4), c(2, 4))
  edges_G1 <- rbind(c(1, 2), c(2, 4), c(1, 4))

  edges_C2 <- rbind(c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(3, 4))
  edges_G2 <- rbind(c(1, 3), c(1, 2), c(2, 3))

  edges_C3 <- rbind(c(2, 3), c(2, 4), c(3, 4), c(1, 3), c(1, 4))
  edges_G3 <- rbind(c(3, 4), c(1, 3), c(1, 4))

  build_recon_cell <- function(cid, t_val, d_val, C_mask, Z_mat, edges_mat, flow = NULL) {
    K <- nrow(edges_mat)
    p <- ncol(Z_mat)
    emeans <- lapply(seq_len(K), function(k) {
      dy_k <- Ymat[, as.character(edges_mat[k, 2L])] - Ymat[, as.character(edges_mat[k, 1L])]
      static_obs_mean(id = uid, A = C_mask, h = rep(1.0, n), route = "baseline",
                      O = O_all, dy = dy_k, h_jacobian = list())
    })
    ecovs <- if (p > 0L) {
      lapply(seq_len(K), function(k) {
        dy_k <- Ymat[, as.character(edges_mat[k, 2L])] - Ymat[, as.character(edges_mat[k, 1L])]
        static_obs_cov(id = uid, A = C_mask, Z = Z_mat, route = "baseline",
                       O = O_all, dy = dy_k)
      })
    } else {
      NULL
    }
    static_reconstructed_or_cell(
      cell_id = cid, unit_id = uid, C = C_mask, Z = Z_mat,
      edges = edges_mat, t = t_val, d = d_val,
      edge_means = emeans, edge_covs = ecovs, flow = flow
    )
  }

  build_complete_cell <- function(cid, t_val, d_val, R_mask, G_mask, C_mask, Z_mat) {
    dy_adj <- Ymat[, as.character(t_val)] - Ymat[, as.character(t_val - 1L)]
    static_complete_cell(
      cell_id = cid, t = t_val, d = d_val, unit_id = uid,
      R = R_mask, G = G_mask, C = C_mask, dy = dy_adj, Z = Z_mat
    )
  }

  # 1. Separate, pooled, and intermediate partitions with covariates
  rc_c1 <- build_recon_cell("c1", 2, 0, C1, Z_ab, edges_C1)
  rc_c2 <- build_recon_cell("c2", 3, 0, C2, Z_ab, edges_C2)
  rc_c3_sameZ <- build_recon_cell("c3", 4, 1, C3, Z_ab, edges_C3)
  rc_c3_diffZ <- build_recon_cell("c3", 4, 1, C3, Z_c, edges_C3)

  cc_c1 <- build_complete_cell("c1", 2, 0, R1, G1, C1, Z_ab)
  cc_c2 <- build_complete_cell("c2", 3, 0, R2, G2, C2, Z_ab)
  cc_c3_sameZ <- build_complete_cell("c3", 4, 1, R3, G3, C3, Z_ab)
  cc_c3_diffZ <- build_complete_cell("c3", 4, 1, R3, G3, C3, Z_c)

  partitions <- list(
    separate = c(c1 = "b1", c2 = "b2", c3 = "b3"),
    pooled = c(c1 = "pool", c2 = "pool", c3 = "pool"),
    intermediate = list(b12 = c("c1", "c2"), b3 = "c3")
  )

  for (pname in names(partitions)) {
    ob <- partitions[[pname]]
    r_cells <- if (pname == "intermediate") list(c1 = rc_c1, c2 = rc_c2, c3 = rc_c3_diffZ) else list(c1 = rc_c1, c2 = rc_c2, c3 = rc_c3_sameZ)
    c_cells <- if (pname == "intermediate") list(c1 = cc_c1, c2 = cc_c2, c3 = cc_c3_diffZ) else list(c1 = cc_c1, c2 = cc_c2, c3 = cc_c3_sameZ)

    fit_rec <- static_reconstructed_or_fit(r_cells, or_blocks = ob, rank_tol = 1e-10, chunk_size = 5L)
    part_comp <- static_complete_validate_partition(c_cells, or_blocks = ob)
    fit_comp <- static_complete_fit_or_partition(c_cells, partition = part_comp, rank_tol = 1e-10, n = n)
    fit_ref <- ref_env$ref_reconstructed_or_eval(r_cells, or_blocks = ob)

    expect_equal(fit_rec$alpha_hat, fit_comp$alpha_hat, tolerance = 1e-13)
    expect_equal(fit_rec$b_hat, fit_comp$slope_hat, tolerance = 1e-13)
    expect_equal(unname(fit_rec$chi_alpha), unname(fit_comp$ell_alpha), tolerance = 1e-12)
    for (bid in names(fit_rec$b_hat)) {
      expect_equal(unname(fit_rec$chi_b[[bid]]), unname(fit_comp$block_ell_slope[[bid]]), tolerance = 1e-12)
    }

    # Agreement with dense full stacked normal system M_full * beta = rhs_full
    expect_equal(fit_rec$alpha_hat, fit_ref$alpha_full, tolerance = 1e-13)
    expect_equal(fit_rec$b_hat, fit_ref$b_full, tolerance = 1e-13)

    # Full cell contrasts (IPW, OR, DR) on complete panel match complete_ref_env$ref_static_complete_eval
    for (m in c("ipw", "or", "dr")) {
      ref_cd <- complete_ref_env$ref_static_complete_eval(c_cells, or_blocks = ob, method = m, polish = TRUE)
      G_edges_list <- list(c1 = edges_G1, c2 = edges_G2, c3 = edges_G3)
      C_edges_list <- list(c1 = edges_C1, c2 = edges_C2, c3 = edges_C3)
      for (cid in names(r_cells)) {
        cl <- c_cells[[cid]]
        Ge <- G_edges_list[[cid]]
        Ce <- C_edges_list[[cid]]
        G_links <- lapply(seq_len(nrow(Ge)), function(k) {
          dy_k <- Ymat[, as.character(Ge[k, 2L])] - Ymat[, as.character(Ge[k, 1L])]
          static_obs_mean(uid, cl$G, rep(1.0, n), "baseline", O = O_all, dy = dy_k, h_jacobian = list())
        })
        if (m %in% c("ipw", "dr")) {
          X_t <- cl$X_tilde
          g_c <- ref_cd$gamma_hat[[cid]]
          u_c <- drop(X_t %*% g_c)
          p_c <- 1.0 / (1.0 + exp(-u_c))
          act_c <- as.double(p_c < (1.0 - 1e-6)) * cl$R
          r_c <- ifelse(act_c == 1, exp(u_c), (1.0 - 1e-6) / 1e-6) * cl$R
          ell_g <- ref_cd$ell_gamma[[cid]]
          rownames(ell_g) <- uid
          h_jac <- list(gamma = X_t * (act_c * r_c))
          rownames(h_jac$gamma) <- uid
          C_links <- lapply(seq_len(nrow(Ce)), function(k) {
            dy_k <- Ymat[, as.character(Ce[k, 2L])] - Ymat[, as.character(Ce[k, 1L])]
            static_obs_mean(uid, cl$C, r_c, "baseline", O = O_all, dy = dy_k,
                            h_jacobian = h_jac, coef_if = list(gamma = ell_g))
          })
        } else {
          C_links <- NULL
        }
        zG_obj <- if (m %in% c("or", "dr")) {
          static_reconstructed_or_role_mean_z(uid, cl$G, cl$Z_or)
        } else {
          NULL
        }
        zCr_obj <- if (m == "dr") {
          static_reconstructed_or_role_mean_z(
            uid, cl$C, cl$Z_or, r = r_c,
            X_tilde = cl$X_tilde, active = act_c,
            ell_gamma = ell_g
          )
        } else {
          NULL
        }
        cent <- static_reconstructed_or_center_cell_links(
          or_fit = fit_rec, cell_id = cid, method = m,
          G_edges = Ge, G_links = G_links,
          C_edges = if (m %in% c("ipw", "dr")) Ce else NULL,
          C_links = C_links, zG = zG_obj, zCr = zCr_obj
        )
        contr <- static_reconstructed_or_cell_contrast(cent, rank_tol = 1e-10)
        expect_equal(contr$theta_hat, ref_cd$theta_hat[[cid]], tolerance = 1e-12)
        expect_equal(unname(contr$phi), unname(ref_cd$phi_theta[, cid]), tolerance = 1e-11)
      }
    }
  }

  # 2. Intercept-only (p = 0) reduction
  Z_empty <- matrix(numeric(0L), nrow = n, ncol = 0L)
  rc0_1 <- build_recon_cell("c1", 2, 0, C1, Z_empty, edges_C1)
  rc0_2 <- build_recon_cell("c2", 3, 0, C2, Z_empty, edges_C2)
  cc0_1 <- build_complete_cell("c1", 2, 0, R1, G1, C1, Z_empty)
  cc0_2 <- build_complete_cell("c2", 3, 0, R2, G2, C2, Z_empty)

  fit_rec0 <- static_reconstructed_or_fit(
    list(c1 = rc0_1, c2 = rc0_2),
    or_blocks = c(c1 = "pool0", c2 = "pool0"),
    rank_tol = 1e-10
  )
  part_comp0 <- static_complete_validate_partition(
    list(c1 = cc0_1, c2 = cc0_2),
    or_blocks = c(c1 = "pool0", c2 = "pool0")
  )
  fit_comp0 <- static_complete_fit_or_partition(
    list(c1 = cc0_1, c2 = cc0_2), partition = part_comp0, rank_tol = 1e-10, n = n
  )
  expect_equal(fit_rec0$alpha_hat, fit_rec0$mC_hat, tolerance = 1e-15)
  expect_equal(fit_rec0$alpha_hat, fit_comp0$alpha_hat, tolerance = 1e-14)
  expect_length(fit_rec0$b_hat[["pool0"]], 0L)
  expect_equal(unname(fit_rec0$chi_alpha), unname(fit_comp0$ell_alpha), tolerance = 1e-13)
})

test_that("early-date and later-date bridges cancel unrestricted non-target trends and level selection (ex:bridge, ex:exit)", {
  # 1. Exact analytical early-date bridge fixture:
  # T = 3, entry cell c = (3, 0). No unit has pair (2, 3) recorded!
  # Types G-A, G-B, C-A, C-B (replicate 5 units per type => n = 20).
  m <- 5L
  n <- 4L * m
  uid <- sprintf("b%02d", seq_len(n))
  G <- c(rep(1L, 2L * m), rep(0L, 2L * m))
  C <- c(rep(0L, 2L * m), rep(1L, 2L * m))
  O_12 <- c(rep(1L, m), rep(0L, m), rep(1L, m), rep(0L, m))
  O_13 <- c(rep(0L, m), rep(1L, m), rep(0L, m), rep(1L, m))

  # Add covariates Z on which the early trend (Y2 - Y1) differs arbitrarily between
  # switchers and stayers, while the target trend (Y3(0) - Y2(0)) shares slope b_true!
  z1 <- rep(seq(-1, 1, length.out = m), 4L)
  Z <- cbind(z1 = z1)
  b_true <- c(z1 = 1.75)
  early_slope_G <- 2.9   # unrestricted switcher slope on (1, 2)
  early_slope_C <- -1.4  # unrestricted stayer slope on (1, 2)

  # Outcome levels depend on type, but pair changes follow:
  # G: dY_12 = 2 + early_slope_G * z1; dY_23 = 5 + b_true * z1 => dY_13 = 7 + (early_slope_G + b_true) * z1
  # C: dY_12 = 3 + early_slope_C * z1; dY_23 = 1 + b_true * z1 => dY_13 = 4 + (early_slope_C + b_true) * z1
  dY_12 <- ifelse(G == 1L, 2.0 + early_slope_G * z1, 3.0 + early_slope_C * z1)
  dY_13 <- ifelse(G == 1L, 7.0 + (early_slope_G + b_true) * z1, 4.0 + (early_slope_C + b_true) * z1)
  # Unobserved pairs carry NA to prove they are never read
  dy_12_obs <- ifelse(O_12 == 1L, dY_12, NA_real_)
  dy_13_obs <- ifelse(O_13 == 1L, dY_13, NA_real_)

  edges_bridge <- rbind(c(1, 2), c(1, 3))
  em_C <- list(
    static_obs_mean(uid, C, rep(1.0, n), "baseline", O = O_12, dy = dy_12_obs, h_jacobian = list()),
    static_obs_mean(uid, C, rep(1.0, n), "baseline", O = O_13, dy = dy_13_obs, h_jacobian = list())
  )
  ec_C <- list(
    static_obs_cov(uid, C, Z, "baseline", O = O_12, dy = dy_12_obs),
    static_obs_cov(uid, C, Z, "baseline", O = O_13, dy = dy_13_obs)
  )
  cell_bridge <- static_reconstructed_or_cell(
    cell_id = "c30", unit_id = uid, C = C, Z = Z,
    edges = edges_bridge, t = 3, d = 0,
    edge_means = em_C, edge_covs = ec_C
  )
  or_fit <- static_reconstructed_or_fit(list(c30 = cell_bridge), or_blocks = c(c30 = "b1"), rank_tol = 1e-10)

  # Flow on {(1,2), (1,3)} for target (2,3) is (-1, +1).
  # Stayer reconstructed mean change is 4 - 3 = 1; slope is b_true = 1.75!
  expect_equal(unname(or_fit$flows$c30$varpi), c(-1.0, 1.0), tolerance = 1e-14)
  expect_equal(or_fit$mC_hat[["c30"]], 1.0, tolerance = 1e-14)
  expect_equal(or_fit$alpha_hat[["c30"]], 1.0, tolerance = 1e-14)
  expect_equal(or_fit$b_hat[["b1"]], b_true, tolerance = 1e-14)

  # Switcher links and centered OR contrast: switcher reconstructed mean is 7 - 2 = 5, theta_hat = 5 - 1 = 4!
  em_G <- list(
    static_obs_mean(uid, G, rep(1.0, n), "baseline", O = O_12, dy = dy_12_obs, h_jacobian = list()),
    static_obs_mean(uid, G, rep(1.0, n), "baseline", O = O_13, dy = dy_13_obs, h_jacobian = list())
  )
  zG <- static_reconstructed_or_role_mean_z(uid, G, Z)
  cent_or <- static_reconstructed_or_center_cell_links(
    or_fit = or_fit, cell_id = "c30", method = "or",
    G_edges = edges_bridge, G_links = em_G, zG = zG
  )
  # Notice edge (1,2) does NOT cross t = 3 (step = 0), while edge (1,3) crosses t = 3 (step = 1)!
  expect_equal(unname(cent_or$step_G), c(0.0, 1.0))
  contr_or <- static_reconstructed_or_cell_contrast(cent_or, rank_tol = 1e-10)
  expect_equal(contr_or$theta_hat, 4.0, tolerance = 1e-14)

  # 2. Exit cell c = (3, 1) with later-date bridge {(2,4), (3,4)}
  edges_exit <- rbind(c(2, 4), c(3, 4))
  fl_exit <- static_reconstructed_or_ls_flow(edges = edges_exit, t = 3, rank_tol = 1e-10)
  expect_equal(unname(fl_exit$varpi), c(1.0, -1.0), tolerance = 1e-14)
})

test_that("missing-outcome cross moments, overlapping cells, cell-specific Z, zero/one-Y units, and bounded chunks match dense reference", {
  # N = 28 units:
  # - units 1..20 appear in various overlapping switcher/stayer roles across cells c1, c2, c3
  # - units 21..24 are stayers in c1 and/or c2 with ONLY ONE recorded outcome (so O_k == 0 on all pairs!)
  # - units 25..28 have ZERO recorded outcomes (O_s == 0 for all dates), with 25..26 in C1/C2 and 27..28 outside all roles!
  n <- 28L
  uid <- sprintf("u%02d", seq_len(n))
  idx <- seq_len(n)

  C1 <- as.integer(idx %in% c(1:14, 21:22, 25))
  G1 <- as.integer(idx %in% c(15:18))
  C2 <- as.integer(idx %in% c(5:18, 23:24, 26))
  G2 <- as.integer(idx %in% c(1:4))
  C3 <- as.integer(idx %in% c(3:16, 21, 25))
  G3 <- as.integer(idx %in% c(17:20))

  # Cell-specific reference covariates Z1 != Z2 (both p = 2, pooled in block "b12") and Z3 (p = 1, block "b3")
  Z1 <- cbind(
    x1 = sin(0.41 * idx) + 0.2 * (idx / n),
    x2 = cos(0.63 * idx) - 0.1 * ((idx %% 4L) - 1.5)
  )
  Z2 <- cbind(
    x1 = cos(0.33 * idx) - 0.15 * (idx / n),
    x2 = sin(0.57 * idx) + 0.25 * ((idx %% 3L) - 1.0)
  )
  Z3 <- cbind(
    w1 = 0.7 * sin(0.49 * idx) - 0.3 * cos(0.27 * idx)
  )
  rownames(Z1) <- rownames(Z2) <- rownames(Z3) <- uid

  # Pair observation indicators: units 21..28 have O_k = 0 on every pair!
  O_pair_A <- as.integer(idx %in% c(1:10, 15:20))
  O_pair_B <- as.integer(idx %in% c(4:18))
  O_pair_C <- as.integer(idx %in% c(1:8, 12:20))

  dy_A <- 1.2 + 0.8 * Z1[, 1L] - 0.4 * Z1[, 2L] + 0.3 * sin(idx)
  dy_B <- -0.5 + 0.5 * Z1[, 1L] + 0.9 * Z1[, 2L] - 0.2 * cos(idx)
  dy_C <- 0.7 + 1.3 * Z1[, 1L] + 0.5 * Z1[, 2L] + 0.15 * sin(2 * idx)

  # Fit real nuisance models (observation logit + imputation OLS) via ref_env helpers
  X_obs_mod <- cbind("(Intercept)" = 1.0, Z1)
  rownames(X_obs_mod) <- uid
  q_fit_c1_e1 <- ref_env$ref_weighted_logit_fit(C1, O_pair_A, X_obs_mod)
  j_fit_c1_e2 <- ref_env$ref_weighted_imputation_fit(C1, O_pair_B, dy_B, X_obs_mod)
  q_fit_c1_e3 <- ref_env$ref_weighted_logit_fit(C1, O_pair_C, X_obs_mod)
  j_fit_c1_e3 <- ref_env$ref_weighted_imputation_fit(C1, O_pair_C, dy_C, X_obs_mod)

  q_fit_c2_e1 <- ref_env$ref_weighted_logit_fit(C2, O_pair_B, cbind("(Intercept)" = 1.0, Z2))
  j_fit_c2_e2 <- ref_env$ref_weighted_imputation_fit(C2, O_pair_C, dy_C, cbind("(Intercept)" = 1.0, Z2))

  coef_if_all <- list(
    eta_c1_e1 = q_fit_c1_e1$coef_if,
    vartheta_c1_e2 = j_fit_c1_e2$coef_if,
    eta_c1_e3 = q_fit_c1_e3$coef_if,
    vartheta_c1_e3 = j_fit_c1_e3$coef_if,
    eta_c2_e1 = q_fit_c2_e1$coef_if,
    vartheta_c2_e2 = j_fit_c2_e2$coef_if
  )

  # Cell 1 uses weighting (edge 1,2), imputation (edge 2,3), and augmented (edge 1,3)
  edges_c1 <- rbind(c(1, 2), c(2, 3), c(1, 3))
  em_c1 <- list(
    static_obs_mean(uid, C1, rep(1, n), "weighting", O = O_pair_A, dy = ifelse(O_pair_A == 1, dy_A, NA_real_),
                    q = q_fit_c1_e1$prob, h_jacobian = list(), q_jacobian = list(eta_c1_e1 = q_fit_c1_e1$prob_jac),
                    coef_if = list(eta_c1_e1 = coef_if_all$eta_c1_e1)),
    static_obs_mean(uid, C1, rep(1, n), "imputation", j = j_fit_c1_e2$j,
                    h_jacobian = list(), j_jacobian = list(vartheta_c1_e2 = j_fit_c1_e2$j_jac),
                    coef_if = list(vartheta_c1_e2 = coef_if_all$vartheta_c1_e2)),
    static_obs_mean(uid, C1, rep(1, n), "augmented", O = O_pair_C, dy = ifelse(O_pair_C == 1, dy_C, NA_real_),
                    q = q_fit_c1_e3$prob, j = j_fit_c1_e3$j, h_jacobian = list(),
                    q_jacobian = list(eta_c1_e3 = q_fit_c1_e3$prob_jac),
                    j_jacobian = list(vartheta_c1_e3 = j_fit_c1_e3$j_jac),
                    coef_if = coef_if_all[c("eta_c1_e3", "vartheta_c1_e3")])
  )
  ec_c1 <- list(
    static_obs_cov(uid, C1, Z1, "weighting", O = O_pair_A, dy = ifelse(O_pair_A == 1, dy_A, NA_real_),
                   q = q_fit_c1_e1$prob, q_jacobian = list(eta_c1_e1 = q_fit_c1_e1$prob_jac),
                   coef_if = list(eta_c1_e1 = coef_if_all$eta_c1_e1)),
    static_obs_cov(uid, C1, Z1, "imputation", j = j_fit_c1_e2$j,
                   j_jacobian = list(vartheta_c1_e2 = j_fit_c1_e2$j_jac),
                   coef_if = list(vartheta_c1_e2 = coef_if_all$vartheta_c1_e2)),
    static_obs_cov(uid, C1, Z1, "augmented", O = O_pair_C, dy = ifelse(O_pair_C == 1, dy_C, NA_real_),
                   q = q_fit_c1_e3$prob, j = j_fit_c1_e3$j,
                   q_jacobian = list(eta_c1_e3 = q_fit_c1_e3$prob_jac),
                   j_jacobian = list(vartheta_c1_e3 = j_fit_c1_e3$j_jac),
                   coef_if = coef_if_all[c("eta_c1_e3", "vartheta_c1_e3")])
  )

  # Cell 2 uses a different role graph {(1,3), (3,4)} for t = 4, baseline = 3
  edges_c2 <- rbind(c(1, 3), c(3, 4), c(1, 4))
  em_c2 <- list(
    static_obs_mean(uid, C2, rep(1, n), "baseline", O = O_pair_A, dy = ifelse(O_pair_A == 1, dy_A, NA_real_), h_jacobian = list()),
    static_obs_mean(uid, C2, rep(1, n), "weighting", O = O_pair_B, dy = ifelse(O_pair_B == 1, dy_B, NA_real_),
                    q = q_fit_c2_e1$prob, h_jacobian = list(), q_jacobian = list(eta_c2_e1 = q_fit_c2_e1$prob_jac),
                    coef_if = list(eta_c2_e1 = coef_if_all$eta_c2_e1)),
    static_obs_mean(uid, C2, rep(1, n), "imputation", j = j_fit_c2_e2$j,
                    h_jacobian = list(), j_jacobian = list(vartheta_c2_e2 = j_fit_c2_e2$j_jac),
                    coef_if = list(vartheta_c2_e2 = coef_if_all$vartheta_c2_e2))
  )
  ec_c2 <- list(
    static_obs_cov(uid, C2, Z2, "baseline", O = O_pair_A, dy = ifelse(O_pair_A == 1, dy_A, NA_real_)),
    static_obs_cov(uid, C2, Z2, "weighting", O = O_pair_B, dy = ifelse(O_pair_B == 1, dy_B, NA_real_),
                   q = q_fit_c2_e1$prob, q_jacobian = list(eta_c2_e1 = q_fit_c2_e1$prob_jac),
                   coef_if = list(eta_c2_e1 = coef_if_all$eta_c2_e1)),
    static_obs_cov(uid, C2, Z2, "imputation", j = j_fit_c2_e2$j,
                   j_jacobian = list(vartheta_c2_e2 = j_fit_c2_e2$j_jac),
                   coef_if = list(vartheta_c2_e2 = coef_if_all$vartheta_c2_e2))
  )

  # Cell 3 in separate block "b3" with p = 1 covariate Z3
  edges_c3 <- rbind(c(2, 3))
  em_c3 <- list(
    static_obs_mean(uid, C3, rep(1, n), "baseline", O = O_pair_C, dy = ifelse(O_pair_C == 1, dy_C, NA_real_), h_jacobian = list())
  )
  ec_c3 <- list(
    static_obs_cov(uid, C3, Z3, "baseline", O = O_pair_C, dy = ifelse(O_pair_C == 1, dy_C, NA_real_))
  )

  cells_list <- list(
    c1 = static_reconstructed_or_cell("c1", uid, C1, Z = Z1, edges = edges_c1, t = 3, d = 0, edge_means = em_c1, edge_covs = ec_c1),
    c2 = static_reconstructed_or_cell("c2", uid, C2, Z = Z2, edges = edges_c2, t = 4, d = 0, edge_means = em_c2, edge_covs = ec_c2),
    c3 = static_reconstructed_or_cell("c3", uid, C3, Z = Z3, edges = edges_c3, t = 3, d = 1, edge_means = em_c3, edge_covs = ec_c3)
  )
  ob_inter <- list(b12 = c("c1", "c2"), b3 = "c3")

  fit_full <- static_reconstructed_or_fit(
    cells_list, or_blocks = ob_inter, rank_tol = 1e-10,
    target_cells = "c1", chunk_size = 4L, materialize_if = TRUE
  )
  ref_out <- ref_env$ref_reconstructed_or_eval(cells_list, or_blocks = ob_inter)

  # Point estimates match both block-wise reference and dense full stacked normal matrix M_full
  expect_equal(fit_full$pi_hat, ref_out$pi_hat, tolerance = 1e-14)
  expect_equal(fit_full$zbar_hat, ref_out$zbar_hat, tolerance = 1e-14)
  expect_equal(fit_full$S_hat, ref_out$S_hat, tolerance = 1e-14)
  expect_equal(fit_full$mC_hat, ref_out$mC_hat, tolerance = 1e-14)
  expect_equal(fit_full$Cov_hat, ref_out$Cov_hat, tolerance = 1e-14)
  expect_equal(fit_full$A_hat, ref_out$A_hat, tolerance = 1e-14)
  expect_equal(fit_full$B_hat, ref_out$B_hat, tolerance = 1e-14)
  expect_equal(fit_full$b_hat, ref_out$b_hat, tolerance = 1e-13)
  expect_equal(fit_full$b_hat, ref_out$b_full, tolerance = 1e-13)
  expect_equal(fit_full$alpha_hat, ref_out$alpha_hat, tolerance = 1e-13)
  expect_equal(fit_full$alpha_hat, ref_out$alpha_full, tolerance = 1e-13)

  # All IF pieces match independent reference across all 28 original units
  if_full_with_A <- fit_full$if_callback(seq_len(n), include_chi_A_tensor = TRUE)
  expect_equal(unname(if_full_with_A$chi_pi), unname(ref_out$chi_pi), tolerance = 1e-14)
  expect_equal(unname(if_full_with_A$chi_mC), unname(ref_out$chi_mC), tolerance = 1e-13)
  expect_equal(unname(if_full_with_A$chi_alpha), unname(ref_out$chi_alpha), tolerance = 1e-13)
  for (cid in names(cells_list)) {
    expect_equal(unname(if_full_with_A$chi_z[[cid]]), unname(ref_out$chi_z[[cid]]), tolerance = 1e-13)
    expect_equal(unname(if_full_with_A$chi_Cov[[cid]]), unname(ref_out$chi_Cov[[cid]]), tolerance = 1e-13)
    expect_equal(unname(if_full_with_A$chi_A_cell[[cid]]), unname(ref_out$chi_A_cell[[cid]]), tolerance = 1e-13)
  }
  for (bid in names(ob_inter)) {
    expect_equal(unname(if_full_with_A$chi_A[[bid]]), unname(ref_out$chi_A[[bid]]), tolerance = 1e-13)
    expect_equal(unname(if_full_with_A$chi_A_b[[bid]]), unname(ref_out$chi_A_b[[bid]]), tolerance = 1e-13)
    expect_equal(unname(if_full_with_A$chi_B[[bid]]), unname(ref_out$chi_B[[bid]]), tolerance = 1e-13)
    expect_equal(unname(if_full_with_A$chi_b[[bid]]), unname(ref_out$chi_b[[bid]]), tolerance = 1e-13)
  }

  # Every IF piece has sample mean zero across all N = 28 original units
  expect_equal(colMeans(if_full_with_A$chi_pi), c(c1 = 0, c2 = 0, c3 = 0), tolerance = 1e-14)
  expect_equal(colMeans(if_full_with_A$chi_mC), c(c1 = 0, c2 = 0, c3 = 0), tolerance = 1e-13)
  expect_equal(colMeans(if_full_with_A$chi_alpha), c(c1 = 0, c2 = 0, c3 = 0), tolerance = 1e-13)
  expect_equal(unname(colMeans(if_full_with_A$chi_b$b12)), c(0, 0), tolerance = 1e-13)
  expect_equal(unname(colMeans(if_full_with_A$chi_b$b3)), 0, tolerance = 1e-13)

  # Units 21..26 (stayers with 0 or 1 recorded outcome, so O_k == 0 on all pairs)
  # have non-trivial known-X and nuisance contributions in chi_b and chi_alpha,
  # while units 27..28 (outside all roles) have non-zero cell-level chi_pi, chi_A_cell,
  # and chi_B_cell that cancel in chi_b by scale homogeneity of b_P in pi_c!
  expect_true(all(apply(abs(if_full_with_A$chi_b$b12[21:26, , drop = FALSE]), 1L, max) > 1e-4))
  expect_true(all(abs(if_full_with_A$chi_alpha[21:26, "c1"]) > 1e-4))
  expect_true(all(abs(if_full_with_A$chi_pi[27:28, "c1"]) > 0.5))
  expect_true(all(apply(abs(if_full_with_A$chi_B_cell$c1[27:28, , drop = FALSE]), 1L, max) > 1e-2))
  expect_equal(unname(if_full_with_A$chi_b$b12[27:28, , drop = FALSE]), matrix(0.0, 2L, 2L), tolerance = 1e-14)

  # Non-target cell c2 in block b12 = {c1, c2} actively influences b_hat$b12 and alpha_hat["c1"]
  cells_pert <- cells_list
  cells_pert$c2$edge_covs[[1L]]$estimate <- cells_pert$c2$edge_covs[[1L]]$estimate + c(0.5, -0.3)
  fit_pert <- static_reconstructed_or_fit(
    cells_pert, or_blocks = ob_inter, rank_tol = 1e-10, target_cells = "c1"
  )
  expect_gt(max(abs(fit_pert$b_hat$b12 - fit_full$b_hat$b12)), 1e-2)
  expect_gt(abs(fit_pert$alpha_hat[["c1"]] - fit_full$alpha_hat[["c1"]]), 1e-3)

  # Bounded chunk sizes and row-subset IF callback without dense IF materialization
  for (cs in c(1L, 2L, 7L, n)) {
    fit_cs <- static_reconstructed_or_fit(
      cells_list, or_blocks = ob_inter, rank_tol = 1e-10,
      chunk_size = cs, materialize_if = FALSE
    )
    expect_null(fit_cs$influence)
    expect_equal(fit_cs$b_hat, fit_full$b_hat, tolerance = 1e-14)
    expect_equal(fit_cs$alpha_hat, fit_full$alpha_hat, tolerance = 1e-14)
    sub_rows <- c(2L, 7L, 19L, 22L, 27L)
    sub_if <- fit_cs$if_callback(sub_rows, include_chi_A_tensor = TRUE)
    expect_equal(sub_if$chi_alpha, if_full_with_A$chi_alpha[sub_rows, , drop = FALSE], tolerance = 1e-14)
    expect_equal(sub_if$chi_b$b12, if_full_with_A$chi_b$b12[sub_rows, , drop = FALSE], tolerance = 1e-14)
    expect_equal(sub_if$chi_A$b12, if_full_with_A$chi_A$b12[sub_rows, , , drop = FALSE], tolerance = 1e-14)
  }
})

test_that("weighted-empirical Gateaux finite differences verify all IF pieces and centered link IFs across all N units", {
  # N = 20 units including stayers, switchers, overlapping roles, 0/1-Y units, and outside units
  n <- 20L
  uid <- sprintf("u%02d", seq_len(n))
  idx <- seq_len(n)

  C1 <- as.integer(idx %in% c(1:11, 17))
  G1 <- as.integer(idx %in% c(12:16))
  R1 <- C1 + G1

  C2 <- as.integer(idx %in% c(4:15, 18))
  G2 <- as.integer(idx %in% c(1:3, 19))
  R2 <- C2 + G2
  # Unit 20 is outside all cells; units 17, 18 have O_k = 0 on all pairs!

  Z1 <- cbind(
    z1 = sin(0.45 * idx) + 0.12 * (idx - 10) / n,
    z2 = cos(0.61 * idx) - 0.18 * ((idx %% 3L) - 1)
  )
  Z2 <- cbind(
    z1 = cos(0.39 * idx) + 0.15 * (idx - 8) / n,
    z2 = sin(0.51 * idx) - 0.10 * ((idx %% 4L) - 1.5)
  )
  rownames(Z1) <- rownames(Z2) <- uid

  O_e1 <- as.integer(idx %in% c(1:8, 12:15))
  O_e2 <- as.integer(idx %in% c(3:11, 13:16))
  O_e3 <- as.integer(idx %in% c(1:6, 9:16))

  dy_e1 <- 0.8 + 0.9 * Z1[, 1L] - 0.4 * Z1[, 2L] + 0.2 * sin(idx)
  dy_e2 <- 1.5 + 1.4 * Z1[, 1L] + 0.3 * Z1[, 2L] - 0.15 * cos(idx)
  dy_e3 <- 0.7 + 0.5 * Z2[, 1L] + 0.7 * Z2[, 2L] + 0.25 * sin(1.5 * idx)

  edges_c1 <- rbind(c(1, 2), c(1, 3))
  edges_c2 <- rbind(c(1, 2), c(2, 3), c(1, 3))

  eval_system <- function(w_u = NULL, use_prod = FALSE) {
    X1_t <- cbind("(Intercept)" = 1.0, Z1)
    X2_t <- cbind("(Intercept)" = 1.0, Z2)
    rownames(X1_t) <- rownames(X2_t) <- uid

    # Fit first-step nuisances under unit weights w_u:
    # - observation logit for c1 edge 1 (weighting) and c1 edge 2 (augmented)
    # - imputation OLS for c1 edge 2 (augmented) and c2 edge 2 (imputation)
    # - treatment logit for c1 on R1
    X_G1 <- cbind("(Intercept)" = 1.0, z1 = Z1[, 1L])
    rownames(X_G1) <- uid
    q_c1_1 <- ref_env$ref_weighted_logit_fit(C1, O_e1, X1_t, unit_weights = w_u)
    q_c1_2 <- ref_env$ref_weighted_logit_fit(C1, O_e2, X1_t, unit_weights = w_u)
    j_c1_2 <- ref_env$ref_weighted_imputation_fit(C1, O_e2, dy_e2, X1_t, unit_weights = w_u)
    j_c2_2 <- ref_env$ref_weighted_imputation_fit(C2, O_e2, dy_e2, X2_t, unit_weights = w_u)
    trt_c1 <- ref_env$ref_weighted_logit_fit(R1, G1, X1_t, unit_weights = w_u)
    q_G1_1 <- ref_env$ref_weighted_logit_fit(G1, O_e1, X_G1, unit_weights = w_u)

    cif <- list(
      q_c1_1 = q_c1_1$coef_if,
      q_c1_2 = q_c1_2$coef_if,
      j_c1_2 = j_c1_2$coef_if,
      j_c2_2 = j_c2_2$coef_if,
      gamma_c1 = trt_c1$coef_if,
      q_G1_1 = q_G1_1$coef_if
    )

    if (use_prod) {
      em_c1 <- list(
        static_obs_mean(uid, C1, rep(1, n), "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                        h_jacobian = list(), q_jacobian = list(q_c1_1 = q_c1_1$prob_jac),
                        coef_if = cif["q_c1_1"]),
        static_obs_mean(uid, C1, rep(1, n), "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob, j = j_c1_2$j,
                        h_jacobian = list(), q_jacobian = list(q_c1_2 = q_c1_2$prob_jac),
                        j_jacobian = list(j_c1_2 = j_c1_2$j_jac),
                        coef_if = cif[c("q_c1_2", "j_c1_2")])
      )
      ec_c1 <- list(
        static_obs_cov(uid, C1, Z1, "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                       q_jacobian = list(q_c1_1 = q_c1_1$prob_jac), coef_if = cif["q_c1_1"]),
        static_obs_cov(uid, C1, Z1, "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob, j = j_c1_2$j,
                       q_jacobian = list(q_c1_2 = q_c1_2$prob_jac),
                       j_jacobian = list(j_c1_2 = j_c1_2$j_jac),
                       coef_if = cif[c("q_c1_2", "j_c1_2")])
      )
      em_c2 <- list(
        static_obs_mean(uid, C2, rep(1, n), "baseline", O = O_e1, dy = dy_e1, h_jacobian = list()),
        static_obs_mean(uid, C2, rep(1, n), "imputation", j = j_c2_2$j,
                        h_jacobian = list(), j_jacobian = list(j_c2_2 = j_c2_2$j_jac),
                        coef_if = cif["j_c2_2"]),
        static_obs_mean(uid, C2, rep(1, n), "baseline", O = O_e3, dy = dy_e3, h_jacobian = list())
      )
      ec_c2 <- list(
        static_obs_cov(uid, C2, Z2, "baseline", O = O_e1, dy = dy_e1),
        static_obs_cov(uid, C2, Z2, "imputation", j = j_c2_2$j,
                       j_jacobian = list(j_c2_2 = j_c2_2$j_jac), coef_if = cif["j_c2_2"]),
        static_obs_cov(uid, C2, Z2, "baseline", O = O_e3, dy = dy_e3)
      )
      cells_p <- list(
        c1 = static_reconstructed_or_cell("c1", uid, C1, Z = Z1, edges = edges_c1, t = 3, d = 0,
                                          edge_means = em_c1, edge_covs = ec_c1),
        c2 = static_reconstructed_or_cell("c2", uid, C2, Z = Z2, edges = edges_c2, t = 3, d = 0,
                                          edge_means = em_c2, edge_covs = ec_c2)
      )
      fit_p <- static_reconstructed_or_fit(cells_p, or_blocks = c(c1 = "pool", c2 = "pool"), rank_tol = 1e-10)
      zG_p <- static_reconstructed_or_role_mean_z(uid, G1, Z1)
      zCr_p <- static_reconstructed_or_role_mean_z(
        uid, C1, Z1, r = trt_c1$odds, X_tilde = X1_t,
        active = trt_c1$active, ell_gamma = trt_c1$coef_if
      )
      G_links_p <- list(
        static_obs_mean(uid, G1, rep(1, n), "weighting", O = O_e1, dy = dy_e1, q = q_G1_1$prob,
                        h_jacobian = list(), q_jacobian = list(q_G1_1 = q_G1_1$prob_jac),
                        coef_if = cif["q_G1_1"]),
        static_obs_mean(uid, G1, rep(1, n), "baseline", O = O_e2, dy = dy_e2, h_jacobian = list())
      )
      C_links_p <- list(
        static_obs_mean(uid, C1, trt_c1$odds, "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                        h_jacobian = list(gamma_c1 = trt_c1$odds_jac),
                        q_jacobian = list(q_c1_1 = q_c1_1$prob_jac),
                        coef_if = cif[c("gamma_c1", "q_c1_1")]),
        static_obs_mean(uid, C1, trt_c1$odds, "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob, j = j_c1_2$j,
                        h_jacobian = list(gamma_c1 = trt_c1$odds_jac),
                        q_jacobian = list(q_c1_2 = q_c1_2$prob_jac),
                        j_jacobian = list(j_c1_2 = j_c1_2$j_jac),
                        coef_if = cif[c("gamma_c1", "q_c1_2", "j_c1_2")])
      )
      cent_or_p <- static_reconstructed_or_center_cell_links(
        fit_p, "c1", "or", G_edges = edges_c1, G_links = G_links_p, zG = zG_p
      )
      cent_dr_p <- static_reconstructed_or_center_cell_links(
        fit_p, "c1", "dr", G_edges = edges_c1, G_links = G_links_p,
        C_edges = edges_c1, C_links = C_links_p, zG = zG_p, zCr = zCr_p
      )
      contr_or_p <- static_reconstructed_or_cell_contrast(cent_or_p, rank_tol = 1e-10)
      contr_dr_p <- static_reconstructed_or_cell_contrast(cent_dr_p, rank_tol = 1e-10)
      return(list(
        fit = fit_p,
        if_full = fit_p$if_callback(seq_len(n), include_chi_A_tensor = TRUE),
        cent_or = cent_or_p,
        cent_dr = cent_dr_p,
        contr_or = contr_or_p,
        contr_dr = contr_dr_p
      ))
    }

    # Weighted reference evaluation for Gateaux finite differences
    em_c1_w <- list(
      ref_env$ref_recon_or_edge_moment("mean", C1, "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                                       h = rep(1, n), unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("mean", C1, "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob,
                                       j = j_c1_2$j, h = rep(1, n), unit_weights = w_u)
    )
    ec_c1_w <- list(
      ref_env$ref_recon_or_edge_moment("cov", C1, "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                                       Z = Z1, unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("cov", C1, "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob,
                                       j = j_c1_2$j, Z = Z1, unit_weights = w_u)
    )
    em_c2_w <- list(
      ref_env$ref_recon_or_edge_moment("mean", C2, "baseline", O = O_e1, dy = dy_e1, h = rep(1, n), unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("mean", C2, "imputation", j = j_c2_2$j, h = rep(1, n), unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("mean", C2, "baseline", O = O_e3, dy = dy_e3, h = rep(1, n), unit_weights = w_u)
    )
    ec_c2_w <- list(
      ref_env$ref_recon_or_edge_moment("cov", C2, "baseline", O = O_e1, dy = dy_e1, Z = Z2, unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("cov", C2, "imputation", j = j_c2_2$j, Z = Z2, unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("cov", C2, "baseline", O = O_e3, dy = dy_e3, Z = Z2, unit_weights = w_u)
    )
    cells_w <- list(
      c1 = list(cell_id = "c1", C = C1, Z = Z1, edges = edges_c1, t = 3, edge_means = em_c1_w, edge_covs = ec_c1_w),
      c2 = list(cell_id = "c2", C = C2, Z = Z2, edges = edges_c2, t = 3, edge_means = em_c2_w, edge_covs = ec_c2_w)
    )
    or_w <- ref_env$ref_reconstructed_or_eval(cells_w, or_blocks = c(c1 = "pool", c2 = "pool"),
                                              unit_weights = w_u, compute_if = FALSE)
    zG_w <- ref_env$ref_role_mean_z(G1, Z1, unit_weights = w_u)
    zCr_w <- ref_env$ref_role_mean_z(C1, Z1, r = trt_c1$odds, unit_weights = w_u)
    G_links_w <- list(
      ref_env$ref_recon_or_edge_moment("mean", G1, "weighting", O = O_e1, dy = dy_e1, q = q_G1_1$prob,
                                       h = rep(1, n), unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("mean", G1, "baseline", O = O_e2, dy = dy_e2,
                                       h = rep(1, n), unit_weights = w_u)
    )
    C_links_w <- list(
      ref_env$ref_recon_or_edge_moment("mean", C1, "weighting", O = O_e1, dy = dy_e1, q = q_c1_1$prob,
                                       h = trt_c1$odds, unit_weights = w_u),
      ref_env$ref_recon_or_edge_moment("mean", C1, "augmented", O = O_e2, dy = dy_e2, q = q_c1_2$prob,
                                       j = j_c1_2$j, h = trt_c1$odds, unit_weights = w_u)
    )
    step_G <- c(0.0, 1.0)
    raw_G <- vapply(G_links_w, `[[`, numeric(1L), "estimate")
    raw_C <- vapply(C_links_w, `[[`, numeric(1L), "estimate")
    center_G <- or_w$alpha_hat[["c1"]] + sum(zG_w$estimate * or_w$b_hat[["pool"]])
    center_C <- or_w$alpha_hat[["c1"]] + sum(zCr_w$estimate * or_w$b_hat[["pool"]])
    adj_G_or <- raw_G - step_G * center_G
    adj_G_dr <- adj_G_or
    adj_C_dr <- raw_C - step_G * center_C
    vp_c1 <- c(-1.0, 1.0)
    list(
      pi = or_w$pi_hat,
      zbar_c1 = or_w$zbar_hat$c1,
      zbar_c2 = or_w$zbar_hat$c2,
      mC = or_w$mC_hat,
      Cov_c1 = or_w$Cov_hat$c1,
      Cov_c2 = or_w$Cov_hat$c2,
      A_pool = as.vector(or_w$A_hat$pool),
      B_pool = or_w$B_hat$pool,
      b_pool = or_w$b_hat$pool,
      alpha = or_w$alpha_hat,
      center_G = center_G,
      center_C = center_C,
      adj_G_or = adj_G_or,
      adj_C_dr = adj_C_dr,
      theta_or = sum(vp_c1 * adj_G_or),
      theta_dr = sum(vp_c1 * adj_G_dr) - sum(vp_c1 * adj_C_dr)
    )
  }

  prod_res <- eval_system(use_prod = TRUE)
  if_full <- prod_res$if_full

  eps <- 1e-5
  fd_pi <- matrix(0.0, n, 2L)
  fd_z_c1 <- matrix(0.0, n, 2L)
  fd_z_c2 <- matrix(0.0, n, 2L)
  fd_mC <- matrix(0.0, n, 2L)
  fd_Cov_c1 <- matrix(0.0, n, 2L)
  fd_Cov_c2 <- matrix(0.0, n, 2L)
  fd_A_pool <- array(0.0, dim = c(n, 2L, 2L))
  fd_B_pool <- matrix(0.0, n, 2L)
  fd_b_pool <- matrix(0.0, n, 2L)
  fd_alpha <- matrix(0.0, n, 2L)
  fd_center_G <- numeric(n)
  fd_center_C <- numeric(n)
  fd_adj_G <- matrix(0.0, n, 2L)
  fd_adj_C <- matrix(0.0, n, 2L)
  fd_theta_or <- numeric(n)
  fd_theta_dr <- numeric(n)

  for (i in seq_len(n)) {
    w_plus <- rep(1.0, n)
    w_minus <- rep(1.0, n)
    w_plus[i] <- 1.0 + eps
    w_minus[i] <- 1.0 - eps
    rp <- eval_system(w_plus, use_prod = FALSE)
    rm <- eval_system(w_minus, use_prod = FALSE)
    scale_fd <- n / (2.0 * eps)

    fd_pi[i, ] <- (rp$pi - rm$pi) * scale_fd
    fd_z_c1[i, ] <- (rp$zbar_c1 - rm$zbar_c1) * scale_fd
    fd_z_c2[i, ] <- (rp$zbar_c2 - rm$zbar_c2) * scale_fd
    fd_mC[i, ] <- (rp$mC - rm$mC) * scale_fd
    fd_Cov_c1[i, ] <- (rp$Cov_c1 - rm$Cov_c1) * scale_fd
    fd_Cov_c2[i, ] <- (rp$Cov_c2 - rm$Cov_c2) * scale_fd
    fd_A_pool[i, , ] <- matrix((rp$A_pool - rm$A_pool) * scale_fd, 2L, 2L)
    fd_B_pool[i, ] <- (rp$B_pool - rm$B_pool) * scale_fd
    fd_b_pool[i, ] <- (rp$b_pool - rm$b_pool) * scale_fd
    fd_alpha[i, ] <- (rp$alpha - rm$alpha) * scale_fd
    fd_center_G[i] <- (rp$center_G - rm$center_G) * scale_fd
    fd_center_C[i] <- (rp$center_C - rm$center_C) * scale_fd
    fd_adj_G[i, ] <- (rp$adj_G_or - rm$adj_G_or) * scale_fd
    fd_adj_C[i, ] <- (rp$adj_C_dr - rm$adj_C_dr) * scale_fd
    fd_theta_or[i] <- (rp$theta_or - rm$theta_or) * scale_fd
    fd_theta_dr[i] <- (rp$theta_dr - rm$theta_dr) * scale_fd
  }

  expect_equal(unname(if_full$chi_pi), fd_pi, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_z$c1), fd_z_c1, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_z$c2), fd_z_c2, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_mC), fd_mC, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_Cov$c1), fd_Cov_c1, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_Cov$c2), fd_Cov_c2, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_A$pool), fd_A_pool, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_B$pool), fd_B_pool, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_b$pool), fd_b_pool, tolerance = 1e-6)
  expect_equal(unname(if_full$chi_alpha), fd_alpha, tolerance = 1e-6)
  expect_equal(unname(prod_res$cent_or$chi_center_G), fd_center_G, tolerance = 1e-6)
  expect_equal(unname(prod_res$cent_dr$chi_center_C), fd_center_C, tolerance = 1e-6)
  expect_equal(unname(prod_res$cent_or$chi_adjusted_G), fd_adj_G, tolerance = 1e-6)
  expect_equal(unname(prod_res$cent_dr$chi_adjusted_C), fd_adj_C, tolerance = 1e-6)
  expect_equal(unname(prod_res$contr_or$phi), fd_theta_or, tolerance = 1e-6)
  expect_equal(unname(prod_res$contr_dr$phi), fd_theta_dr, tolerance = 1e-6)
})

test_that("fixed switch-step centering satisfies prop:invariance, prop:imp-or, and single-correct nuisance algebra", {
  n <- 24L
  uid <- sprintf("u%02d", seq_len(n))
  idx <- seq_len(n)

  G <- as.integer(idx %in% 1:8)
  C <- as.integer(idx %in% 9:22)
  R <- G + C
  Z <- cbind(z1 = sin(0.37 * idx), z2 = cos(0.59 * idx))
  rownames(Z) <- uid
  X_t <- cbind("(Intercept)" = 1.0, Z)
  rownames(X_t) <- uid

  # 1. Non-switch date centering invariance: shifting date-level centering at non-switch
  # dates s != t-1, t alters individual adjusted links touching s, but leaves
  # flow-contracted theta_hat^{or, dr}_c and phi^{or, dr}_c invariant!
  O_all <- rep(1L, n)
  edges <- rbind(c(1, 2), c(2, 3), c(1, 3), c(2, 4), c(3, 4))
  Y_dates <- cbind(
    "1" = 0.3 + 0.5 * Z[, 1L] - 0.2 * Z[, 2L] + 0.1 * sin(idx),
    "2" = 1.0 + 1.1 * Z[, 1L] + 0.4 * Z[, 2L] - 0.1 * cos(idx),
    "3" = 2.5 + 1.8 * Z[, 1L] - 0.3 * Z[, 2L] + 0.2 * sin(2 * idx),
    "4" = 3.2 + 2.2 * Z[, 1L] + 0.6 * Z[, 2L] - 0.15 * cos(2 * idx)
  )
  em_C <- lapply(seq_len(nrow(edges)), function(k) {
    dy_k <- Y_dates[, as.character(edges[k, 2L])] - Y_dates[, as.character(edges[k, 1L])]
    static_obs_mean(uid, C, rep(1, n), "baseline", O = O_all, dy = dy_k, h_jacobian = list())
  })
  ec_C <- lapply(seq_len(nrow(edges)), function(k) {
    dy_k <- Y_dates[, as.character(edges[k, 2L])] - Y_dates[, as.character(edges[k, 1L])]
    static_obs_cov(uid, C, Z, "baseline", O = O_all, dy = dy_k)
  })
  cell <- static_reconstructed_or_cell("c1", uid, C, Z = Z, edges = edges, t = 3, d = 0,
                                       edge_means = em_C, edge_covs = ec_C)
  or_fit <- static_reconstructed_or_fit(list(c1 = cell), or_blocks = c(c1 = "b1"), rank_tol = 1e-10)

  trt <- ref_env$ref_weighted_logit_fit(R, G, X_t)
  em_G <- lapply(seq_len(nrow(edges)), function(k) {
    dy_k <- Y_dates[, as.character(edges[k, 2L])] - Y_dates[, as.character(edges[k, 1L])]
    static_obs_mean(uid, G, rep(1, n), "baseline", O = O_all, dy = dy_k, h_jacobian = list())
  })
  em_Cr <- lapply(seq_len(nrow(edges)), function(k) {
    dy_k <- Y_dates[, as.character(edges[k, 2L])] - Y_dates[, as.character(edges[k, 1L])]
    static_obs_mean(uid, C, trt$odds, "baseline", O = O_all, dy = dy_k,
                    h_jacobian = list(g = trt$odds_jac), coef_if = list(g = trt$coef_if))
  })
  zG <- static_reconstructed_or_role_mean_z(uid, G, Z)
  zCr <- static_reconstructed_or_role_mean_z(uid, C, Z, r = trt$odds, X_tilde = X_t,
                                             active = trt$active, ell_gamma = trt$coef_if)

  # Base switch-step centering vs shifted date-level centering at dates 1 and 4
  shift_if <- cbind(
    "1" = 0.4 * or_fit$chi_alpha[, "c1"] + drop(or_fit$chi_b$b1 %*% c(0.2, -0.3)),
    "2" = rep(0.0, n),
    "3" = rep(0.0, n),
    "4" = -0.7 * or_fit$chi_alpha[, "c1"] + drop(or_fit$chi_b$b1 %*% c(-0.5, 0.1))
  )
  ds_G <- list(dates = c(1, 2, 3, 4), values = c(1.35, 0.0, 0.0, -2.10), if_mat = shift_if)
  ds_C <- list(dates = c(1, 2, 3, 4), values = c(-0.85, 0.0, 0.0, 1.65), if_mat = -0.5 * shift_if)

  for (m in c("or", "dr")) {
    cent_0 <- static_reconstructed_or_center_cell_links(
      or_fit, "c1", m, G_edges = edges, G_links = em_G,
      C_edges = if (m == "dr") edges else NULL,
      C_links = if (m == "dr") em_Cr else NULL,
      zG = zG, zCr = if (m == "dr") zCr else NULL
    )
    cent_s <- static_reconstructed_or_center_cell_links(
      or_fit, "c1", m, G_edges = edges, G_links = em_G,
      C_edges = if (m == "dr") edges else NULL,
      C_links = if (m == "dr") em_Cr else NULL,
      zG = zG, zCr = if (m == "dr") zCr else NULL,
      date_shift_G = ds_G, date_shift_C = if (m == "dr") ds_C else NULL
    )
    # Individual adjusted links at non-adjacent edges differ:
    expect_gt(max(abs(cent_s$adjusted_G - cent_0$adjusted_G)), 1.0)
    # Yet for every valid flow (LS flow and explicit path flow), theta_hat and phi are invariant!
    c0_ls <- static_reconstructed_or_cell_contrast(cent_0, rank_tol = 1e-10)
    cs_ls <- static_reconstructed_or_cell_contrast(cent_s, rank_tol = 1e-10)
    expect_equal(cs_ls$theta_hat, c0_ls$theta_hat, tolerance = 1e-13)
    expect_equal(cs_ls$phi, c0_ls$phi, tolerance = 1e-13)

    flow_tree <- c(-1.0, 0.0, 1.0, 0.0, 0.0) # path (1,3) - (1,2)
    c0_tr <- static_reconstructed_or_cell_contrast(cent_0, rank_tol = 1e-10, flow_G = flow_tree,
                                                    flow_C = if (m == "dr") flow_tree else NULL)
    cs_tr <- static_reconstructed_or_cell_contrast(cent_s, rank_tol = 1e-10, flow_G = flow_tree,
                                                    flow_C = if (m == "dr") flow_tree else NULL)
    expect_equal(cs_tr$theta_hat, c0_tr$theta_hat, tolerance = 1e-13)
    expect_equal(cs_tr$phi, c0_tr$phi, tolerance = 1e-13)
  }

  # 2. Imputation-OR equivalence: when stayer edges use imputation with design (1, Z_c)
  # and P(c) = {c}, DR equals OR whatever the treatment odds r_c are!
  O_miss <- as.integer(idx %in% c(1:6, 9:17))
  edges_2 <- rbind(c(1, 2), c(1, 3))
  j_fits <- lapply(seq_len(nrow(edges_2)), function(k) {
    dy_k <- Y_dates[, as.character(edges_2[k, 2L])] - Y_dates[, as.character(edges_2[k, 1L])]
    ref_env$ref_weighted_imputation_fit(C, O_miss, dy_k, X_t)
  })
  em_imp_C0 <- lapply(seq_len(nrow(edges_2)), function(k) {
    static_obs_mean(uid, C, rep(1, n), "imputation", j = j_fits[[k]]$j,
                    h_jacobian = list(), j_jacobian = list(j = j_fits[[k]]$j_jac),
                    coef_if = list(j = j_fits[[k]]$coef_if))
  })
  ec_imp_C <- lapply(seq_len(nrow(edges_2)), function(k) {
    static_obs_cov(uid, C, Z, "imputation", j = j_fits[[k]]$j,
                   j_jacobian = list(j = j_fits[[k]]$j_jac),
                   coef_if = list(j = j_fits[[k]]$coef_if))
  })
  # Use an arbitrary distorted treatment odds vector r_arb
  r_arb <- exp(0.9 * Z[, 1L] - 1.3 * Z[, 2L])
  em_imp_Cr <- lapply(seq_len(nrow(edges_2)), function(k) {
    static_obs_mean(uid, C, r_arb, "imputation", j = j_fits[[k]]$j,
                    h_jacobian = list(), j_jacobian = list(j = j_fits[[k]]$j_jac),
                    coef_if = list(j = j_fits[[k]]$coef_if))
  })
  em_G2 <- lapply(seq_len(nrow(edges_2)), function(k) {
    dy_k <- Y_dates[, as.character(edges_2[k, 2L])] - Y_dates[, as.character(edges_2[k, 1L])]
    static_obs_mean(uid, G, rep(1, n), "baseline", O = O_miss, dy = dy_k, h_jacobian = list())
  })
  cell_imp <- static_reconstructed_or_cell("c_imp", uid, C, Z = Z, edges = edges_2, t = 3, d = 0,
                                           edge_means = em_imp_C0, edge_covs = ec_imp_C)
  fit_imp <- static_reconstructed_or_fit(list(c_imp = cell_imp), or_blocks = c(c_imp = "b1"), rank_tol = 1e-10)
  zCr_arb <- static_reconstructed_or_role_mean_z(uid, C, Z, r = r_arb)

  cent_imp_or <- static_reconstructed_or_center_cell_links(
    fit_imp, "c_imp", "or", G_edges = edges_2, G_links = em_G2, zG = zG
  )
  cent_imp_dr <- static_reconstructed_or_center_cell_links(
    fit_imp, "c_imp", "dr", G_edges = edges_2, G_links = em_G2,
    C_edges = edges_2, C_links = em_imp_Cr, zG = zG, zCr = zCr_arb
  )
  contr_imp_or <- static_reconstructed_or_cell_contrast(cent_imp_or, rank_tol = 1e-10)
  contr_imp_dr <- static_reconstructed_or_cell_contrast(cent_imp_dr, rank_tol = 1e-10)
  expect_equal(contr_imp_dr$theta_hat, contr_imp_or$theta_hat, tolerance = 1e-13)
  expect_equal(contr_imp_dr$phi, contr_imp_or$phi, tolerance = 1e-12)

  # 3. Single-correct nuisance model / IF generic algebra on augmented route:
  # When actual-change imputation j is exact on observed units (dy == j on O == 1),
  # the observation-probability loading J^eta is identically 0 and exact parameters
  # are recovered for ANY positive observation probability q!
  b_exact <- c(z1 = 1.25, z2 = -0.75)
  alpha_exact <- 2.10
  dy_exact <- alpha_exact + drop(Z %*% b_exact)
  q_wrong <- plogis(-0.4 + 0.8 * Z[, 1L] - 0.5 * Z[, 2L])
  q_jac_dummy <- list(eta_wrong = X_t * (q_wrong * (1 - q_wrong)))
  j_jac_exact <- list(vartheta = X_t)
  coef_if_dummy <- list(
    eta_wrong = matrix(0.1, nrow = n, ncol = 3L, dimnames = list(uid, colnames(X_t))),
    vartheta = matrix(0.0, nrow = n, ncol = 3L, dimnames = list(uid, colnames(X_t)))
  )
  em_aug <- static_obs_mean(uid, C, rep(1, n), "augmented", O = O_miss, dy = dy_exact,
                            q = q_wrong, j = dy_exact, h_jacobian = list(),
                            q_jacobian = q_jac_dummy, j_jacobian = j_jac_exact,
                            coef_if = coef_if_dummy)
  ec_aug <- static_obs_cov(uid, C, Z, "augmented", O = O_miss, dy = dy_exact,
                           q = q_wrong, j = dy_exact,
                           q_jacobian = q_jac_dummy, j_jacobian = j_jac_exact,
                           coef_if = coef_if_dummy)
  expect_equal(unname(em_aug$loadings$eta_wrong), c(0, 0, 0), tolerance = 1e-14)
  expect_equal(unname(ec_aug$loadings$eta_wrong), matrix(0, 2L, 3L), tolerance = 1e-14)

  cell_aug <- static_reconstructed_or_cell("c_aug", uid, C, Z = Z, edges = rbind(c(2, 3)), t = 3,
                                           edge_means = list(em_aug), edge_covs = list(ec_aug))
  fit_aug <- static_reconstructed_or_fit(list(c_aug = cell_aug), or_blocks = c(c_aug = "b1"), rank_tol = 1e-10)
  expect_equal(fit_aug$b_hat$b1, b_exact, tolerance = 1e-13)
  expect_equal(fit_aug$alpha_hat[["c_aug"]], alpha_exact, tolerance = 1e-13)

  # Conversely, under the probability arm of the augmented route (when O_k / q_k == 1
  # on average conditional on Z in the stayer population), the imputation loading J^vartheta
  # is zero and exact stayer OR parameters are recovered for ANY misspecified j_wrong!
  q_exact <- rep(mean(O_miss[C == 1L]), n)
  # Construct balanced O_bal on C == 1 so 1 - O_bal / q_bal is orthogonal to (1, Z) on C == 1
  j_wrong <- -3.5 + 4.2 * Z[, 1L] - 2.8 * Z[, 2L] + 0.9 * sin(idx)
  em_aug_q <- static_obs_mean(uid, C, rep(1, n), "augmented", O = O_all, dy = dy_exact,
                              q = rep(1.0, n), j = j_wrong, h_jacobian = list(),
                              q_jacobian = q_jac_dummy, j_jacobian = j_jac_exact,
                              coef_if = coef_if_dummy)
  ec_aug_q <- static_obs_cov(uid, C, Z, "augmented", O = O_all, dy = dy_exact,
                             q = rep(1.0, n), j = j_wrong,
                             q_jacobian = q_jac_dummy, j_jacobian = j_jac_exact,
                             coef_if = coef_if_dummy)
  expect_equal(unname(em_aug_q$loadings$vartheta), c(0, 0, 0), tolerance = 1e-14)
  expect_equal(unname(ec_aug_q$loadings$vartheta), matrix(0, 2L, 3L), tolerance = 1e-14)
  cell_aug_q <- static_reconstructed_or_cell("c_aug_q", uid, C, Z = Z, edges = rbind(c(2, 3)), t = 3,
                                             edge_means = list(em_aug_q), edge_covs = list(ec_aug_q))
  fit_aug_q <- static_reconstructed_or_fit(list(c_aug_q = cell_aug_q), or_blocks = c(c_aug_q = "b1"), rank_tol = 1e-10)
  expect_equal(fit_aug_q$b_hat$b1, b_exact, tolerance = 1e-13)
  expect_equal(fit_aug_q$alpha_hat[["c_aug_q"]], alpha_exact, tolerance = 1e-13)

  # 4. Direct parameter-channel Jacobians (loadings$B, loadings$b, loadings$alpha)
  # and multi-cell centered catalog `static_reconstructed_or_center_links`:
  eval_theta <- function(theta_vec) {
    q_th <- plogis(drop(X_t %*% theta_vec))
    q_jac_th <- list(shared_theta = X_t * (q_th * (1.0 - q_th)))
    j_th <- drop(X_t %*% (0.5 * theta_vec))
    j_jac_th <- list(shared_theta = 0.5 * X_t)
    cif_th <- list(shared_theta = function(rows) trt$coef_if[rows, , drop = FALSE])

    em1 <- static_obs_mean(uid, C, rep(1, n), "augmented", O = O_miss, dy = dy_exact,
                           q = q_th, j = j_th, h_jacobian = list(),
                           q_jacobian = q_jac_th, j_jacobian = j_jac_th)
    ec1 <- static_obs_cov(uid, C, Z, "augmented", O = O_miss, dy = dy_exact,
                          q = q_th, j = j_th, q_jacobian = q_jac_th, j_jacobian = j_jac_th)
    cl1 <- static_reconstructed_or_cell("c1", uid, C, Z = Z, edges = rbind(c(2, 3)), t = 3,
                                        edge_means = list(em1), edge_covs = list(ec1),
                                        coef_if = cif_th)
    static_reconstructed_or_fit(list(c1 = cl1), or_blocks = c(c1 = "b1"), rank_tol = 1e-10)
  }
  th0 <- c(-0.2, 0.4, -0.3)
  fit_th0 <- eval_theta(th0)
  eps_th <- 1e-6
  jac_B_num <- matrix(0.0, 2L, 3L)
  jac_b_num <- matrix(0.0, 2L, 3L)
  jac_a_num <- numeric(3L)
  for (k in 1:3) {
    d_k <- numeric(3L)
    d_k[k] <- eps_th
    fp <- eval_theta(th0 + d_k)
    fm <- eval_theta(th0 - d_k)
    jac_B_num[, k] <- (fp$B_hat$b1 - fm$B_hat$b1) / (2.0 * eps_th)
    jac_b_num[, k] <- (fp$b_hat$b1 - fm$b_hat$b1) / (2.0 * eps_th)
    jac_a_num[k] <- (fp$alpha_hat[["c1"]] - fm$alpha_hat[["c1"]]) / (2.0 * eps_th)
  }
  expect_equal(unname(fit_th0$loadings$B$b1$shared_theta), jac_B_num, tolerance = 1e-7)
  expect_equal(unname(fit_th0$loadings$b$b1$shared_theta), jac_b_num, tolerance = 1e-7)
  expect_equal(unname(fit_th0$loadings$alpha$c1$shared_theta), jac_a_num, tolerance = 1e-7)

  # Multi-cell catalog wrapper `static_reconstructed_or_center_links`
  cat_dr <- static_reconstructed_or_center_links(
    or_fit = or_fit,
    cell_links = list(
      c1 = list(G_edges = edges, G_links = em_G, C_edges = edges, C_links = em_Cr, zG = zG, zCr = zCr)
    ),
    method = "dr",
    materialize_if = FALSE
  )
  expect_null(cat_dr$influence)
  sub_cat <- cat_dr$if_callback(c(1L, 5L, 12L))
  expect_equal(nrow(sub_cat$chi_stacked), 3L)
  expect_equal(ncol(sub_cat$chi_stacked), 2L * nrow(edges))
})

test_that("deliberate formula mutations fail decisively against the true IF and centering", {
  n <- 20L
  uid <- sprintf("u%02d", seq_len(n))
  idx <- seq_len(n)
  C1 <- as.integer(idx %in% 1:11)
  C2 <- as.integer(idx %in% 6:17)
  G1 <- as.integer(idx %in% 12:16)
  Z1 <- cbind(z1 = sin(0.41 * idx), z2 = cos(0.53 * idx))
  Z2 <- cbind(z1 = cos(0.31 * idx), z2 = sin(0.67 * idx))
  rownames(Z1) <- rownames(Z2) <- uid

  O1 <- as.integer(idx %in% c(1:7, 12:16))
  O2 <- as.integer(idx %in% c(6:14))
  dy1 <- 1.0 + 1.5 * Z1[, 1L] - 0.8 * Z1[, 2L] + 0.3 * sin(idx)
  dy2 <- -0.4 + 0.6 * Z2[, 1L] + 1.1 * Z2[, 2L] - 0.25 * cos(idx)

  c1 <- static_reconstructed_or_cell(
    "c1", uid, C1, Z = Z1, edges = rbind(c(2, 3)), t = 3,
    edge_means = list(static_obs_mean(uid, C1, rep(1, n), "baseline", O = O1, dy = dy1, h_jacobian = list())),
    edge_covs = list(static_obs_cov(uid, C1, Z1, "baseline", O = O1, dy = dy1))
  )
  c2 <- static_reconstructed_or_cell(
    "c2", uid, C2, Z = Z2, edges = rbind(c(2, 3)), t = 3,
    edge_means = list(static_obs_mean(uid, C2, rep(1, n), "baseline", O = O2, dy = dy2, h_jacobian = list())),
    edge_covs = list(static_obs_cov(uid, C2, Z2, "baseline", O = O2, dy = dy2))
  )
  fit <- static_reconstructed_or_fit(list(c1 = c1, c2 = c2), or_blocks = c(c1 = "pool", c2 = "pool"), rank_tol = 1e-10)

  # Mutation 1: Omit random stayer share derivative (C_i - pi_c) Cov_c in chi_B
  chi_B_no_pi <- fit$pi_hat[["c1"]] * fit$chi_Cov$c1 + fit$pi_hat[["c2"]] * fit$chi_Cov$c2
  chi_b_no_pi <- t(solve(fit$A_hat$pool, t(chi_B_no_pi - fit$chi_A_b$pool)))
  expect_gt(max(abs(chi_b_no_pi - fit$chi_b$pool)), 0.05)

  # Mutation 2: Omit known-X Gram derivative chi_A * b_P in chi_b and chi_z' b_P in chi_alpha
  chi_b_no_gram <- t(solve(fit$A_hat$pool, t(fit$chi_B$pool)))
  chi_alpha_no_z <- fit$chi_mC[, "c1"] - drop(fit$chi_b$pool %*% fit$zbar_hat$c1)
  expect_gt(max(abs(chi_b_no_gram - fit$chi_b$pool)), 0.10)
  expect_gt(max(abs(chi_alpha_no_z - fit$chi_alpha[, "c1"])), 0.10)

  # Mutation 3: Use wrong centering boundary 1{u >= t-1} - 1{s >= t-1} or scale by gap length (u - s)
  edges_G <- rbind(c(1, 2), c(1, 3))
  G_links <- list(
    static_obs_mean(uid, G1, rep(1, n), "baseline", O = O1, dy = dy1, h_jacobian = list()),
    static_obs_mean(uid, G1, rep(1, n), "baseline", O = O1, dy = dy1 + 2.0, h_jacobian = list())
  )
  zG <- static_reconstructed_or_role_mean_z(uid, G1, Z1)
  cent_true <- static_reconstructed_or_center_cell_links(
    fit, "c1", "or", G_edges = edges_G, G_links = G_links, zG = zG
  )
  # Wrong boundary at t - 1 = 2 gives step = (1, 1) instead of (0, 1), which cancels along flow (-1, +1)!
  step_wrong_boundary <- as.double(edges_G[, 2L] >= 2) - as.double(edges_G[, 1L] >= 2)
  adj_wrong_boundary <- cent_true$raw_G - step_wrong_boundary * cent_true$center_G
  theta_wrong_boundary <- sum(c(-1, 1) * adj_wrong_boundary)
  contr_true <- static_reconstructed_or_cell_contrast(cent_true, rank_tol = 1e-10)
  expect_gt(abs(theta_wrong_boundary - contr_true$theta_hat), 0.5)

  # Wrong gap-length scaling (u - s) * step_G gives (0, 2) instead of (0, 1) on (1, 3)!
  step_gap_scaled <- (edges_G[, 2L] - edges_G[, 1L]) * cent_true$step_G
  adj_gap_scaled <- cent_true$raw_G - step_gap_scaled * cent_true$center_G
  theta_gap_scaled <- sum(c(-1, 1) * adj_gap_scaled)
  expect_gt(abs(theta_gap_scaled - contr_true$theta_hat), 0.5)

  # Mutation 4: Normalize IF by stayer count sum(C1) instead of original N
  chi_b_wrong_scale <- fit$chi_b$pool * (sum(C1) / n)
  expect_gt(max(abs(chi_b_wrong_scale - fit$chi_b$pool)), 0.10)
})

test_that("deterministic validation refuses invalid graphs, rank deficiency, bad flows, and mismatched inputs", {
  n <- 10L
  uid <- paste0("u", seq_len(n))
  C <- c(rep(1L, 6L), rep(0L, 4L))
  Z <- cbind(z1 = seq_len(n) * 1.0)
  em1 <- list(list(estimate = 1.0, influence = rep(0.0, n)))
  ec1 <- list(list(estimate = c(z1 = 0.5), influence = matrix(0.0, n, 1L, dimnames = list(uid, "z1"))))

  # Missing or invalid rank_tol
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(2, 3)), t = 3),
    class = "cdid_static_reconstructed_or_rank_tol_error"
  )
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(2, 3)), t = 3, rank_tol = 0),
    class = "cdid_static_reconstructed_or_rank_tol_error"
  )

  # Graph missing baseline t-1 or target t
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(1, 2)), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(1, 3)), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )

  # Disconnected required endpoints: {(1,2), (3,4)} for t = 3 (contains 2 and 3 in distinct components -> refused!)
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(1, 2), c(3, 4)), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )
  # Isolated required endpoint (target t = 3 present in dates but has degree 0 -> refused!)
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(1, 2)), dates = c(1, 2, 3), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )
  # Connected required endpoints (t - 1 = 1 and t = 2 in component 1) with extraneous component (3,4)
  # and harmless isolated non-target date 7 are accepted with exact zero flow on (3,4).
  fl_ext <- static_reconstructed_or_ls_flow(
    edges = rbind(c(1, 2), c(3, 4)), dates = c(1, 2, 3, 4, 7), t = 2, rank_tol = 1e-10
  )
  expect_identical(fl_ext$n_components, 3L)
  expect_identical(fl_ext$rank, 2L)
  expect_identical(unname(fl_ext$varpi), c(1.0, 0.0))

  # Invalid edge orientation s >= u or duplicate edge
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(3, 2)), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(2, 3), c(2, 3)), t = 3, rank_tol = 1e-10),
    class = "cdid_static_reconstructed_or_graph_error"
  )

  # Explicit flow violating D' varpi = ell
  expect_error(
    static_reconstructed_or_ls_flow(edges = rbind(c(1, 2), c(1, 3)), t = 3,
                                    rank_tol = 1e-10, flow = c(1.0, 1.0)),
    class = "cdid_static_reconstructed_or_flow_error"
  )

  # Empty stayer mask
  expect_error(
    static_reconstructed_or_cell("c0", uid, rep(0L, n), Z = Z, edges = rbind(c(2, 3)), t = 3,
                                 edge_means = em1, edge_covs = ec1),
    class = "cdid_static_reconstructed_or_support_error"
  )

  # Rank-deficient pooled known-X Gram matrix
  Z_collinear <- cbind(z1 = seq_len(n) * 1.0, z2 = seq_len(n) * 2.0)
  ec_col <- list(list(estimate = c(z1 = 0.5, z2 = 1.0),
                      influence = matrix(0.0, n, 2L, dimnames = list(uid, c("z1", "z2")))))
  cell_col <- static_reconstructed_or_cell("c_col", uid, C, Z = Z_collinear, edges = rbind(c(2, 3)), t = 3,
                                           edge_means = em1, edge_covs = ec_col)
  expect_error(
    static_reconstructed_or_fit(list(c_col = cell_col), or_blocks = c(c_col = "b1"), rank_tol = 1e-10),
    class = "cdid_or_wls_rank_error"
  )

  # Block-constant covariate in c1 identified when pooled with c2
  Z_c1_const <- cbind(z1 = seq_len(n) * 1.0, z2 = rep(5.0, n))
  Z_c2_var <- cbind(z1 = rev(seq_len(n)) * 1.0, z2 = (seq_len(n) %% 2L) * 1.0)
  cell_k1 <- static_reconstructed_or_cell("c1", uid, C, Z = Z_c1_const, edges = rbind(c(2, 3)), t = 3,
                                          edge_means = em1, edge_covs = ec_col)
  cell_k2 <- static_reconstructed_or_cell("c2", uid, C, Z = Z_c2_var, edges = rbind(c(2, 3)), t = 3,
                                          edge_means = em1, edge_covs = ec_col)
  # Separate fit on c1 alone fails rank check:
  expect_error(
    static_reconstructed_or_fit(list(c1 = cell_k1), or_blocks = c(c1 = "b1"), rank_tol = 1e-10),
    class = "cdid_or_wls_rank_error"
  )
  # Pooled fit across {c1, c2} succeeds and matches dense reference!
  fit_pool_k <- static_reconstructed_or_fit(
    list(c1 = cell_k1, c2 = cell_k2), or_blocks = c(c1 = "pool", c2 = "pool"), rank_tol = 1e-10
  )
  ref_pool_k <- ref_env$ref_reconstructed_or_eval(
    list(c1 = cell_k1, c2 = cell_k2), or_blocks = c(c1 = "pool", c2 = "pool")
  )
  expect_equal(fit_pool_k$b_hat$pool, ref_pool_k$b_hat$pool, tolerance = 1e-13)

  # Mislabeled coefficient-IF parameter columns in uncontracted loading contraction are refused,
  # while aligned named and unnamed explicit-position coefficient IFs succeed
  i_seq <- seq_len(n)
  H_sh <- cbind(a = 1.0, b = sin(i_seq), c = cos(0.7 * i_seq))
  rownames(H_sh) <- uid
  L_sh <- cbind(a = i_seq * 1.0, b = cos(i_seq), c = (i_seq %% 3L) * 1.0)
  L_sh <- sweep(L_sh, 2L, colMeans(L_sh))
  rownames(L_sh) <- uid
  Z_sh <- cbind(z1 = i_seq / 3.0, z2 = sin(0.4 * i_seq))
  rownames(Z_sh) <- uid
  j_sh <- drop(H_sh %*% c(a = 0.4, b = 0.8, c = -0.6)) + 0.13 * sin(1.7 * i_seq)
  m_sh <- static_obs_mean(uid, C, rep(1.0, n), "imputation", j = j_sh,
                          h_jacobian = list(), j_jacobian = list(shared = H_sh))
  v_sh <- static_obs_cov(uid, C, Z_sh, "imputation", j = j_sh,
                         j_jacobian = list(shared = H_sh))
  expect_error(
    static_reconstructed_or_cell(
      "c_perm", uid, C, Z = Z_sh, edges = rbind(c(2, 3)), t = 3,
      edge_means = list(m_sh), edge_covs = list(v_sh),
      coef_if = list(shared = L_sh[, c(2L, 3L, 1L), drop = FALSE])
    ),
    class = "cdid_static_reconstructed_or_channel_error"
  )
  cell_sh_named <- static_reconstructed_or_cell(
    "c_sh", uid, C, Z = Z_sh, edges = rbind(c(2, 3)), t = 3,
    edge_means = list(m_sh), edge_covs = list(v_sh),
    coef_if = list(shared = L_sh)
  )
  L_sh_unnamed <- L_sh
  colnames(L_sh_unnamed) <- NULL
  cell_sh_unnamed <- static_reconstructed_or_cell(
    "c_sh", uid, C, Z = Z_sh, edges = rbind(c(2, 3)), t = 3,
    edge_means = list(m_sh), edge_covs = list(v_sh),
    coef_if = list(shared = L_sh_unnamed)
  )
  fit_sh_named <- static_reconstructed_or_fit(list(c_sh = cell_sh_named), or_blocks = c(c_sh = "P"), rank_tol = 1e-10)
  fit_sh_unnamed <- static_reconstructed_or_fit(list(c_sh = cell_sh_unnamed), or_blocks = c(c_sh = "P"), rank_tol = 1e-10)
  expect_equal(fit_sh_unnamed$chi_b$P, fit_sh_named$chi_b$P, tolerance = 1e-15)
  expect_equal(fit_sh_unnamed$chi_alpha, fit_sh_named$chi_alpha, tolerance = 1e-15)

  # Reordered row labels on ell_gamma (and mismatched column labels) in role_mean_z are refused,
  # while aligned named and unlabeled positional ell_gamma succeed
  r_sh <- drop(exp(H_sh %*% c(0.1, -0.2, 0.3)))
  role_good <- static_reconstructed_or_role_mean_z(
    uid, C, Z_sh, r = r_sh, X_tilde = H_sh, active = rep(1L, n), ell_gamma = L_sh
  )
  L_sh_no_rn <- L_sh
  rownames(L_sh_no_rn) <- NULL
  role_unlabeled <- static_reconstructed_or_role_mean_z(
    uid, C, Z_sh, r = r_sh, X_tilde = H_sh, active = rep(1L, n), ell_gamma = L_sh_no_rn
  )
  expect_equal(role_unlabeled$influence, role_good$influence, tolerance = 1e-15)
  expect_error(
    static_reconstructed_or_role_mean_z(
      uid, C, Z_sh, r = r_sh, X_tilde = H_sh, active = rep(1L, n),
      ell_gamma = L_sh[rev(i_seq), , drop = FALSE]
    ),
    class = "cdid_static_reconstructed_or_order_error"
  )
  expect_error(
    static_reconstructed_or_role_mean_z(
      uid, C, Z_sh, r = r_sh, X_tilde = H_sh, active = rep(1L, n),
      ell_gamma = L_sh[, c(2L, 3L, 1L), drop = FALSE]
    ),
    class = "cdid_static_reconstructed_or_channel_error"
  )
})

test_that("multi-component role graphs in reconstructed OR and downstream joint GMM: component gauges, zero non-endpoint flow, row retention, gauge/relabeling invariance, and end-to-end GMM", {
  gmm_ref_env <- new.env(parent = asNamespace("stats"))
  sys.source(
    testthat::test_path("reference", "static_graph_gmm_dense_reference.R"),
    envir = gmm_ref_env
  )
  prod_names <- ls(prod_env, all.names = TRUE)
  for (rnm in ls(gmm_ref_env, all.names = TRUE)) {
    fn <- get(rnm, envir = gmm_ref_env, inherits = FALSE)
    if (is.function(fn)) {
      used_globals <- codetools::findGlobals(fn, merge = FALSE)$functions
      expect_length(intersect(used_globals, prod_names), 0L)
    }
  }

  # 1. Multi-component stayer graph for cell c = (baseline = 17, t = 25) on nonconsecutive dates:
  #    - Component 1 (endpoint triangle): {10, 17, 25} with edges (10,17), (10,25), (17,25)
  #    - Component 2 (disconnected tree): {40, 52, 65} with edges (40,52), (52,65)
  #    - Component 3 (disconnected cycle): {100, 115, 130} with edges (100,115), (115,130), (100,130)
  #    - Component 4 (harmless isolated non-target vertex): {200}
  edges_C_multi <- rbind(
    c(10, 17),
    c(10, 25),
    c(17, 25),
    c(40, 52),
    c(52, 65),
    c(100, 115),
    c(115, 130),
    c(100, 130)
  )
  dates_C_multi <- c(10, 17, 25, 40, 52, 65, 100, 115, 130, 200)

  fl_C <- static_reconstructed_or_ls_flow(
    edges = edges_C_multi,
    t = 25,
    baseline_date = 17,
    dates = dates_C_multi,
    rank_tol = 1e-10,
    cell_id = "c_multi"
  )
  ref_inc_C <- ref_env$ref_recon_or_incidence(
    edges = edges_C_multi,
    t = 25,
    baseline_date = 17,
    dates = dates_C_multi
  )

  # 10 vertices across 4 components -> grounded incidence D has 8 rows and 10 - 4 = 6 columns (rank 6)
  expect_identical(fl_C$n_components, 4L)
  expect_identical(fl_C$isolated_dates, 200.0)
  expect_identical(dim(fl_C$D), c(8L, 6L))
  expect_identical(fl_C$rank, 6L)
  expect_equal(fl_C$D, ref_inc_C$D, tolerance = 1e-14)
  expect_equal(fl_C$ell, ref_inc_C$ell, tolerance = 1e-14)
  expect_equal(fl_C$varpi, ref_inc_C$varpi, tolerance = 1e-14)

  # Exact IEEE 0.0 minimum-norm flow on all edges outside the endpoint component (components 2 and 3)
  expect_identical(unname(fl_C$varpi[4:8]), rep(0.0, 5L))
  expect_equal(unname(fl_C$varpi[1:3]), c(-1 / 3, 1 / 3, 2 / 3), tolerance = 1e-14)

  # 2. Coordinate gauge invariance across all choices of reference date in components 2 and 3
  for (r2 in c(40, 52, 65)) {
    for (r3 in c(100, 115, 130)) {
      fl_alt <- static_reconstructed_or_ls_flow(
        edges = edges_C_multi,
        t = 25,
        baseline_date = 17,
        dates = dates_C_multi,
        reference_dates = c(17, r2, r3),
        rank_tol = 1e-10,
        cell_id = "c_multi"
      )
      ref_alt <- ref_env$ref_recon_or_incidence(
        edges = edges_C_multi,
        t = 25,
        baseline_date = 17,
        dates = dates_C_multi,
        reference_dates = c(17, r2, r3)
      )
      expect_identical(dim(fl_alt$D), c(8L, 6L))
      expect_identical(fl_alt$rank, 6L)
      expect_equal(fl_alt$D, ref_alt$D, tolerance = 1e-14)
      expect_identical(fl_alt$varpi, fl_C$varpi)
      expect_identical(unname(fl_alt$varpi[4:8]), rep(0.0, 5L))
    }
  }

  # 3. Reconstructed OR fit on multi-component stayer graph with direct + nuisance loading channels:
  #    non-endpoint edges remain validated and retained in the cell specification, while their
  #    exact zero minimum-norm flow prevents non-endpoint moments/loadings from polluting (alpha_hat, b_hat)
  n <- 24L
  uid <- sprintf("unit_%02d", seq_len(n))
  i_seq <- seq_len(n)
  C_mask <- as.integer((i_seq %% 3L) != 0L)
  G_mask <- 1L - C_mask
  Z_mat <- cbind(
    z1 = 0.4 * sin(0.37 * i_seq) - 0.1 * cos(0.51 * i_seq),
    z2 = 0.5 * cos(0.29 * i_seq) + 0.2 * ((i_seq %% 2L) - 0.5)
  )
  rownames(Z_mat) <- uid

  # Keyed nuisance coefficient IFs for endpoint component ("nuis_end") and other components ("nuis_other")
  L_end <- cbind(p1 = sin(0.41 * i_seq), p2 = cos(0.63 * i_seq))
  L_end <- sweep(L_end, 2L, colMeans(L_end), "-")
  rownames(L_end) <- uid
  L_other <- cbind(q1 = cos(0.53 * i_seq), q2 = sin(0.79 * i_seq))
  L_other <- sweep(L_other, 2L, colMeans(L_other), "-")
  rownames(L_other) <- uid
  coef_if_list <- list(
    nuis_end = function(rows) L_end[rows, , drop = FALSE],
    nuis_other = function(rows) L_other[rows, , drop = FALSE]
  )

  K_C <- nrow(edges_C_multi)
  edge_means_C <- vector("list", K_C)
  edge_covs_C <- vector("list", K_C)
  for (k in seq_len(K_C)) {
    dir_m <- 0.3 * sin(0.19 * k * i_seq) * C_mask
    dir_m <- dir_m - mean(dir_m)
    names(dir_m) <- uid
    dir_c <- cbind(
      z1 = 0.25 * cos(0.23 * k * i_seq) * C_mask,
      z2 = 0.20 * sin(0.31 * k * i_seq) * C_mask
    )
    dir_c <- sweep(dir_c, 2L, colMeans(dir_c), "-")
    rownames(dir_c) <- uid

    if (k <= 3L) {
      load_m <- list(nuis_end = c(p1 = 0.15 * k, p2 = -0.10 * k))
      load_c <- list(nuis_end = rbind(z1 = c(p1 = 0.08 * k, p2 = -0.05 * k),
                                      z2 = c(p1 = -0.04 * k, p2 = 0.09 * k)))
    } else {
      load_m <- list(nuis_other = c(q1 = 0.50 * k, q2 = -0.35 * k))
      load_c <- list(nuis_other = rbind(z1 = c(q1 = 0.40 * k, q2 = -0.25 * k),
                                        z2 = c(q1 = -0.30 * k, q2 = 0.45 * k)))
    }
    inf_m <- dir_m + drop(if (k <= 3L) L_end %*% load_m$nuis_end else L_other %*% load_m$nuis_other)
    inf_c <- dir_c + (if (k <= 3L) L_end %*% t(load_c$nuis_end) else L_other %*% t(load_c$nuis_other))

    edge_means_C[[k]] <- list(
      estimate = 0.20 * k - 0.05 * (k^2),
      direct_if = dir_m,
      loadings = load_m,
      influence = inf_m
    )
    edge_covs_C[[k]] <- list(
      estimate = c(z1 = 0.12 * k, z2 = -0.08 * k + 0.03),
      direct_if = dir_c,
      loadings = load_c,
      influence = inf_c
    )
  }

  cell_multi <- static_reconstructed_or_cell(
    cell_id = "c_multi",
    unit_id = uid,
    C = C_mask,
    Z = Z_mat,
    edges = edges_C_multi,
    t = 25,
    baseline_date = 17,
    dates = dates_C_multi,
    edge_means = edge_means_C,
    edge_covs = edge_covs_C,
    coef_if = coef_if_list,
    or_block = "P1"
  )
  fit_or_multi <- static_reconstructed_or_fit(
    cells = list(c_multi = cell_multi),
    or_blocks = c(c_multi = "P1"),
    rank_tol = 1e-10,
    chunk_size = 7L,
    materialize_if = TRUE
  )
  ref_or_multi <- ref_env$ref_reconstructed_or_eval(
    cells = list(c_multi = cell_multi),
    or_blocks = c(c_multi = "P1")
  )

  expect_equal(fit_or_multi$alpha_hat, ref_or_multi$alpha_hat, tolerance = 1e-13)
  expect_equal(fit_or_multi$b_hat$P1, ref_or_multi$b_hat$P1, tolerance = 1e-13)
  expect_equal(unname(fit_or_multi$chi_alpha), unname(ref_or_multi$chi_alpha), tolerance = 1e-12)
  expect_equal(unname(fit_or_multi$chi_b$P1), unname(ref_or_multi$chi_b$P1), tolerance = 1e-12)
  # Non-endpoint nuisance loading key "nuis_other" received exact 0 weight and is not accumulated
  expect_identical(names(fit_or_multi$cells$c_multi$loadings_mC), "nuis_end")
  expect_identical(names(fit_or_multi$cells$c_multi$loadings_Cov), "nuis_end")
  expect_identical(names(fit_or_multi$loadings$b$P1), "nuis_end")
  expect_identical(names(fit_or_multi$loadings$alpha$c_multi), "nuis_end")

  # Arbitrary out-of-order positional row subset callback contract preserved
  sub_rows <- c(19L, 4L, 12L, 1L, 23L)
  sub_if <- fit_or_multi$if_callback(sub_rows)
  expect_identical(sub_if$unit_id, uid[sub_rows])
  expect_equal(sub_if$chi_alpha, fit_or_multi$chi_alpha[sub_rows, , drop = FALSE], tolerance = 1e-13)
  expect_equal(sub_if$chi_b$P1, fit_or_multi$chi_b$P1[sub_rows, , drop = FALSE], tolerance = 1e-13)

  # Invariance of reconstructed OR fit under alternative non-endpoint coordinate gauges
  cell_multi_g2 <- static_reconstructed_or_cell(
    cell_id = "c_multi",
    unit_id = uid,
    C = C_mask,
    Z = Z_mat,
    edges = edges_C_multi,
    t = 25,
    baseline_date = 17,
    dates = dates_C_multi,
    reference_dates = c(17, 65, 130),
    edge_means = edge_means_C,
    edge_covs = edge_covs_C,
    coef_if = coef_if_list,
    or_block = "P1"
  )
  fit_or_g2 <- static_reconstructed_or_fit(
    cells = list(c_multi = cell_multi_g2),
    or_blocks = c(c_multi = "P1"),
    rank_tol = 1e-10
  )
  expect_equal(fit_or_g2$alpha_hat, fit_or_multi$alpha_hat, tolerance = 1e-14)
  expect_equal(fit_or_g2$b_hat$P1, fit_or_multi$b_hat$P1, tolerance = 1e-14)
  expect_equal(fit_or_g2$chi_alpha, fit_or_multi$chi_alpha, tolerance = 1e-14)
  expect_equal(fit_or_g2$chi_b$P1, fit_or_multi$chi_b$P1, tolerance = 1e-14)

  # Invariance of reconstructed OR fit under strictly increasing date relabeling s -> 2 * s + 5
  map_d <- function(x) 2.0 * x + 5.0
  cell_multi_relab <- static_reconstructed_or_cell(
    cell_id = "c_multi",
    unit_id = uid,
    C = C_mask,
    Z = Z_mat,
    edges = cbind(map_d(edges_C_multi[, 1L]), map_d(edges_C_multi[, 2L])),
    t = map_d(25),
    baseline_date = map_d(17),
    dates = map_d(dates_C_multi),
    edge_means = edge_means_C,
    edge_covs = edge_covs_C,
    coef_if = coef_if_list,
    or_block = "P1"
  )
  fit_or_relab <- static_reconstructed_or_fit(
    cells = list(c_multi = cell_multi_relab),
    or_blocks = c(c_multi = "P1"),
    rank_tol = 1e-10
  )
  expect_equal(fit_or_relab$alpha_hat, fit_or_multi$alpha_hat, tolerance = 1e-14)
  expect_equal(fit_or_relab$b_hat$P1, fit_or_multi$b_hat$P1, tolerance = 1e-14)
  expect_equal(fit_or_relab$chi_alpha, fit_or_multi$chi_alpha, tolerance = 1e-14)
  expect_equal(fit_or_relab$chi_b$P1, fit_or_multi$chi_b$P1, tolerance = 1e-14)

  # 4. Downstream fixed switch-step link centering and joint GMM across unequal
  #    multi-component role-specific graphs (G vs C), retaining all declared edges
  # G has endpoint edge (17,25) + disconnected cycle {(70,85), (85,95), (70,95)} + isolated date 150
  edges_G_multi <- rbind(
    c(17, 25),
    c(70, 85),
    c(85, 95),
    c(70, 95)
  )
  dates_G_multi <- c(17, 25, 70, 85, 95, 150)
  K_G <- nrow(edges_G_multi)
  G_links <- lapply(seq_len(K_G), function(k) {
    inf_g <- (0.35 * cos(0.27 * k * i_seq) + 0.15 * sin(0.43 * i_seq))
    inf_g <- inf_g - mean(inf_g)
    names(inf_g) <- uid
    list(estimate = 0.90 - 0.15 * k, influence = inf_g)
  })
  C_links <- lapply(seq_len(K_C), function(k) {
    list(estimate = edge_means_C[[k]]$estimate, influence = edge_means_C[[k]]$influence)
  })

  zG_obj <- list(
    estimate = c(z1 = 0.18, z2 = -0.11),
    influence = sweep(
      cbind(z1 = 0.22 * sin(0.31 * i_seq), z2 = 0.19 * cos(0.47 * i_seq)),
      2L,
      colMeans(cbind(z1 = 0.22 * sin(0.31 * i_seq), z2 = 0.19 * cos(0.47 * i_seq))),
      "-"
    )
  )
  rownames(zG_obj$influence) <- uid
  zCr_obj <- list(
    estimate = c(z1 = -0.07, z2 = 0.14),
    influence = sweep(
      cbind(z1 = 0.17 * cos(0.39 * i_seq), z2 = 0.21 * sin(0.53 * i_seq)),
      2L,
      colMeans(cbind(z1 = 0.17 * cos(0.39 * i_seq), z2 = 0.21 * sin(0.53 * i_seq))),
      "-"
    )
  )
  rownames(zCr_obj$influence) <- uid

  centered_dr <- static_reconstructed_or_center_cell_links(
    or_fit = fit_or_multi,
    cell_id = "c_multi",
    method = "dr",
    t = 25,
    G_edges = edges_G_multi,
    G_links = G_links,
    C_edges = edges_C_multi,
    C_links = C_links,
    zG = zG_obj,
    zCr = zCr_obj,
    materialize_if = TRUE
  )
  ref_centered_dr <- ref_env$ref_center_cell_links(
    ref_or = ref_or_multi,
    cell_id = "c_multi",
    block_id = "P1",
    method = "dr",
    t = 25,
    G_edges = edges_G_multi,
    G_links = G_links,
    C_edges = edges_C_multi,
    C_links = C_links,
    zG = zG_obj,
    zCr = zCr_obj
  )

  # All K_G + K_C = 4 + 8 = 12 edge rows retained in original order
  expect_identical(length(centered_dr$adjusted_G), 4L)
  expect_identical(length(centered_dr$adjusted_C), 8L)
  expect_identical(length(centered_dr$stacked_links), 12L)
  expect_identical(dim(centered_dr$chi_stacked), c(n, 12L))
  # Switch-step indicator [1{u >= 25} - 1{s >= 25}] is 0 on all non-endpoint edges (s, u >= 40 > 25)
  expect_identical(unname(centered_dr$step_G), c(1.0, 0.0, 0.0, 0.0))
  expect_identical(unname(centered_dr$step_C), c(0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
  expect_equal(unname(centered_dr$adjusted_G), unname(ref_centered_dr$adjusted_G), tolerance = 1e-13)
  expect_equal(unname(centered_dr$adjusted_C), unname(ref_centered_dr$adjusted_C), tolerance = 1e-13)
  expect_equal(
    unname(centered_dr$chi_stacked),
    unname(cbind(ref_centered_dr$chi_adjusted_G, ref_centered_dr$chi_adjusted_C)),
    tolerance = 1e-12
  )

  # Feed the full 12-moment centered catalog into solve_static_graph_gmm()
  g_link_ids <- paste0("G:", names(centered_dr$adjusted_G))
  c_link_ids <- paste0("C:", names(centered_dr$adjusted_C))
  blk_g_gmm <- static_role_block_spec(
    block_id = "c_multi_G",
    cell_id = "c_multi",
    role = "G",
    baseline_time = 17,
    target_time = 25,
    dates = dates_G_multi,
    edges = data.frame(
      link_id = g_link_ids,
      s = edges_G_multi[, 1L],
      u = edges_G_multi[, 2L],
      stringsAsFactors = FALSE
    )
  )
  blk_c_gmm <- static_role_block_spec(
    block_id = "c_multi_C",
    cell_id = "c_multi",
    role = "C",
    baseline_time = 17,
    target_time = 25,
    dates = dates_C_multi,
    edges = data.frame(
      link_id = c_link_ids,
      s = edges_C_multi[, 1L],
      u = edges_C_multi[, 2L],
      stringsAsFactors = FALSE
    )
  )
  mhat_joint <- centered_dr$stacked_links
  phi_joint <- centered_dr$chi_stacked
  cb_joint <- make_static_if_batch_callback(
    phi_source = function(idx) centered_dr$if_callback(idx)$chi_stacked,
    unit_ids = uid,
    batch_rows = 6L,
    link_ids = names(mhat_joint)
  )
  bounds_j <- gmm_bounds_control(
    max_matrix_entries = 100000L,
    max_batch_rows = 6L,
    max_targets = 20L,
    max_open_clusters = 1L,
    allow_materialize_if = TRUE,
    max_materialized_if_entries = 100000L
  )
  acc_j <- static_graph_gmm_accuracy_control()

  for (wt in c("identity", "two_step")) {
    fit_j <- solve_static_graph_gmm(
      blocks = list(blk_g_gmm, blk_c_gmm),
      method = "dr",
      mhat = mhat_joint,
      if_callback = cb_joint,
      unit_ids = uid,
      cluster_ids = NULL,
      weighting = wt,
      bounds = bounds_j,
      accuracy = acc_j,
      materialize_target_if = TRUE
    )
    ref_j <- gmm_ref_env$dense_reference_static_graph_gmm(
      blocks = list(blk_g_gmm, blk_c_gmm),
      method = "dr",
      mhat = mhat_joint,
      phi = phi_joint,
      weighting = wt
    )
    expect_identical(fit_j$design$links$link_id, names(mhat_joint))
    expect_equal(fit_j$estimates, ref_j$estimates, tolerance = 1e-11)
    expect_equal(fit_j$moment_weights$Lambda, ref_j$Lambda, tolerance = 1e-11)
    expect_equal(fit_j$target_if, ref_j$target_if, tolerance = 1e-11)
    expect_equal(fit_j$vcov, ref_j$vcov_n2, tolerance = 1e-11)
    expect_equal(fit_j$se, ref_j$se_n2, tolerance = 1e-11)
  }
})

