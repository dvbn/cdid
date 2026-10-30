# Portable reliability campaign: single-operation R worker and freeze helper.
# Inert when sourced or parsed; executes only when invoked via Rscript CLI.
say <- function(...) cat(..., "\n", sep = ""); now_utc <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
secs <- function(t0) {
  d <- proc.time() - t0; list(
    elapsed = unname(d[["elapsed"]]), user = unname(d[["user.self"]] + d[["user.child"]]),
    system = unname(d[["sys.self"]] + d[["sys.child"]])
  )
}; write_json_atomic <- function(x, path, overwrite = FALSE) {
  if (!overwrite && file.exists(path)) stop("Refusing to overwrite existing file: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, digits = I(17), null = "null", na = "null", pretty = TRUE)
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path, call. = FALSE); invisible(path)
}; save_rds_atomic <- function(obj, path, overwrite = FALSE) {
  if (!overwrite && file.exists(path)) stop("Refusing to overwrite existing RDS: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE); saveRDS(obj, tmp, version = 3)
  if (!file.rename(tmp, path)) stop("Atomic RDS rename failed for ", path, call. = FALSE); invisible(path)
}; obj_md5 <- function(x) {
  f <- tempfile(); on.exit(unlink(f), add = TRUE); con <- file(f, "wb"); serialize(x, con, xdr = TRUE); close(con); unname(tools::md5sum(f))
}; file_hash <- function(p) {
  if (!file.exists(p)) stop("Missing identity file: ", p, call. = FALSE)
  paste0("md5:", unname(tools::md5sum(p)))
}
# Hash DESCRIPTION excluding build-time Packaged/Built lines so portable installs
# from the same source commit match without pinning an author build timestamp.
description_portable_hash <- function(desc_path) {
  if (!file.exists(desc_path)) return(NA_character_); lines <- readLines(desc_path, warn = FALSE, encoding = "UTF-8")
  keep <- lines[!grepl("^(Packaged|Built):", lines)]; obj_md5(paste(keep, collapse = "\n"))
}; run_captured <- function(f) {
  warns <- list(); err <- NULL; value <- withCallingHandlers(
    tryCatch(f(), error = function(e) { err <<- e; NULL }), warning = function(w) {
      warns[[length(warns) + 1L]] <<- list(class = class(w), message = conditionMessage(w)); invokeRestart("muffleWarning")
    }
  ); list(value = value, error = err, warnings = warns)
}; cond_rec <- function(e) {
  if (is.null(e)) return(NULL); list(
    class = class(e), message = conditionMessage(e),
    call = if (!is.null(conditionCall(e))) paste(utils::head(deparse(conditionCall(e)), 2L), collapse = " ") else NULL
  )
}
# Confirm the installed public interface and source functions against this checkout.
assert_source_identity <- function(repo_path, formals_path, source_sha) {
  fixed <- jsonlite::read_json(formals_path, simplifyVector = FALSE)
  checks <- list(CDID_ATT_GT_FORMALS = cdid::att_gt, DID_ATT_GT_FORMALS = did::att_gt,
                 CDID_AGGTE_FORMALS = cdid::aggte, CDID_PRETEST_FORMALS = cdid::conditional_did_pretest)
  for (key in names(checks)) {
    actual <- setdiff(names(formals(checks[[key]])), "...")
    if (!identical(actual, unlist(fixed[[key]], use.names = FALSE))) stop("Installed formal mismatch: ", key, call. = FALSE)
  }
  desc <- read.dcf(file.path(getNamespaceInfo(asNamespace("cdid"), "path"), "DESCRIPTION"))
  if ("RemoteSha" %in% colnames(desc) && !identical(unname(desc[1L, "RemoteSha"]), source_sha)) {
    stop("Installed cdid RemoteSha differs from frozen source", call. = FALSE)
  }
  # Compare syntax, without sourcing package code or evaluating any estimator.
  syntax <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "\n")
  checked <- character(); ns <- asNamespace("cdid")
  for (path in sort(list.files(file.path(repo_path, "R"), pattern = "\\.[Rr]$", full.names = TRUE))) {
    for (expr in parse(file = path, keep.source = FALSE)) {
      if (is.call(expr) && as.character(expr[[1L]]) %in% c("<-", "=") && length(expr) == 3L &&
          is.symbol(expr[[2L]]) && is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
        name <- as.character(expr[[2L]]); installed <- get0(name, envir = ns, inherits = FALSE)
        expected <- expr[[3L]]
        if (!is.function(installed) || !identical(syntax(formals(installed)), syntax(expected[[2L]])) ||
            !identical(syntax(body(installed)), syntax(expected[[3L]]))) stop("Installed cdid/source mismatch: ", name, call. = FALSE)
        checked <- c(checked, name)
      }
    }
  }
  if (!all(c("att_gt", "aggte", "conditional_did_pretest") %in% checked)) stop("Public source functions were not verified", call. = FALSE)
  list(source_sha = source_sha, n_functions = length(checked), functions_hash = obj_md5(sort(checked)))
}
# Verify R >= 4.1 and exact package guard: cdid 1.1.0.0, did 2.5.1, DRDID 1.3.0
setup_and_assert <- function(lib_paths = character(), frozen_runtime = NULL, repo_path, formals_path, source_sha) {
  if (length(lib_paths) > 0L) {
    norm_libs <- vapply(lib_paths[nzchar(lib_paths)], normalizePath, character(1L), mustWork = TRUE)
    .libPaths(unique(c(norm_libs, .libPaths())))
  }; RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"); r_ver <- getRversion(); r_ok <- r_ver >= "4.1.0"; loadNamespace("DRDID"); loadNamespace("did"); loadNamespace("cdid")
  ver <- function(p) as.character(getNamespaceVersion(p))
  path_of <- function(p) normalizePath(find.package(p, quiet = FALSE, verbose = FALSE), mustWork = TRUE); checks <- list(
    R_ge_4_1 = isTRUE(r_ok), cdid_version = identical(ver("cdid"), "1.1.0.0"), did_version = identical(ver("did"), "2.5.1"),
    DRDID_version = identical(ver("DRDID"), "1.3.0")
  ); pkgs <- c("cdid", "did", "DRDID"); pkg_meta <- lapply(stats::setNames(pkgs, pkgs), function(p) {
    pp <- path_of(p); list(
      version = ver(p), description_portable_hash = description_portable_hash(file.path(pp, "DESCRIPTION")),
      rdb_hash = file_hash(file.path(pp, "R", paste0(p, ".rdb")))
    )
  }); source_identity <- assert_source_identity(repo_path, formals_path, source_sha)
  si <- utils::sessionInfo(); ext <- base::extSoftVersion(); runtime_id <- list(
    source_identity = source_identity, R_version = as.character(r_ver), R_version_string = R.version.string, platform = R.version$platform,
    blas = if (!is.null(si$BLAS)) as.character(si$BLAS) else "", lapack = if (!is.null(si$LAPACK)) as.character(si$LAPACK) else "",
    ext_soft = as.list(ext[intersect(c("BLAS", "LAPACK", "zlib"), names(ext))]), packages = pkg_meta,
    loaded_namespaces = lapply(sort(loadedNamespaces()), function(p) list(package = p, version = ver(p),
      description_hash = description_portable_hash(file.path(path_of(p), "DESCRIPTION")),
      rdb_hash = if (file.exists(file.path(path_of(p), "R", paste0(p, ".rdb")))) file_hash(file.path(path_of(p), "R", paste0(p, ".rdb"))) else NULL)),
    RNG_kind = RNGkind(), session_info = capture.output(print(si))
  ); profile_digest <- obj_md5(runtime_id); if (!is.null(frozen_runtime)) {
    checks$frozen_profile_match <- identical(profile_digest, frozen_runtime$profile_digest)
  }; list(ok = all(unlist(checks)), checks = checks, runtime = runtime_id, profile_digest = profile_digest)
}
# Full IF self-audit following Section 6.3
inference_record <- function(fit, package, full_N, compute_inffunc = TRUE, point_only_rcs = FALSE, original_frame = NULL) {
  if (!isTRUE(compute_inffunc) || isTRUE(point_only_rcs) || is.null(fit[["inffunc"]]) || is.null(fit[["V_analytical"]])) {
    return(list(
      applicable = FALSE, shape = NA, covariance = NA, SE = NA, max_abs_SE_error = NA_real_,
      criterion = "point-only or compute_inffunc=FALSE; IF self-audit not applicable"
    ))
  }; IF <- as.matrix(fit[["inffunc"]]); V <- as.matrix(fit[["V_analytical"]]); J <- length(fit$att); N <- as.integer(fit$n)
  ok <- is.finite(fit$att) & is.finite(fit$se); is_cdid_mp <- inherits(fit, "cdid_MP"); if (is_cdid_mp && !is.null(fit$chained$available)) {
    ok <- ok & as.logical(fit$chained$available)
  }; same <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12, check.attributes = FALSE))
  shape <- identical(dim(IF), c(N, J)) && identical(dim(V), c(J, J)) && (!is_cdid_mp || N == as.integer(full_N)); cov_ok <- se_ok <- FALSE
  max_se <- NA_real_; if (shape && any(ok) && all(is.finite(IF[, ok, drop = FALSE]))) {
    U <- IF[, ok, drop = FALSE]; if (!is_cdid_mp) {
      rec <- crossprod(U) / N
    } else {
      Uc <- sweep(U, 2L, colMeans(U), "-")
      cluster_name <- setdiff(fit$DIDparams$clustervars, fit$DIDparams$idname)
      if (length(cluster_name) == 1L && !is.null(original_frame)) {
        if (nrow(original_frame) != N) stop("Original unit frame/IF shape mismatch")
        labels <- rownames(IF); ix <- if (is.null(labels)) seq_len(N) else match(labels, as.character(original_frame$id))
        if (anyNA(ix)) stop("IF unit labels do not match original unit frame")
        Ucl <- rowsum(Uc, original_frame[[cluster_name]][ix], reorder = FALSE); rec <- crossprod(Ucl) / N
      } else {
        rec <- crossprod(Uc) / N
      }
    }; cov_ok <- same(V[ok, ok, drop = FALSE], rec); se_ana <- sqrt(pmax(0.0, diag(rec)) / N)
    # When bstrap=FALSE, fit$se equals analytical SE; when bstrap=TRUE, compare against V_analytical diag
    se_ref <- if (isTRUE(fit$DIDparams$bstrap)) sqrt(pmax(0.0, diag(V)[ok]) / N) else as.double(fit$se[ok]); se_ok <- same(se_ref, se_ana)
    max_se <- max(abs(se_ref - se_ana))
  }; list(
    applicable = TRUE, shape = shape, n_result_units = N, n_cells = J, positive_cells = sum(ok),
    covariance = cov_ok, SE = se_ok, max_abs_SE_error = max_se,
    criterion = "MP$n normalization; did uncentered crossprod/n; cdid centered (and cluster-summed when clustered) crossprod/n; tol 1e-12"
  )
}; sample_record <- function(fit) {
  if (inherits(fit, "cdid_MP") && !is.null(fit$chained$frame)) {
    f <- fit$chained$frame; return(list(
      n_result_units = as.integer(fit$n), parsed_frame_md5 = obj_md5(f), frame_columns = names(f),
      basis = "chained full population frame; unavailable outcomes retained"
    ))
  }; dp <- fit$DIDparams; d <- if (!is.null(dp)) dp$data else NULL
  if (is.null(d)) return(list(n_result_units = as.integer(fit$n), parsed_data = "not exposed")); d <- as.data.frame(d); list(
    n_result_units = as.integer(fit$n), parsed_data_md5 = obj_md5(d), n_parsed_rows = nrow(d),
    n_parsed_ids = if (!is.null(dp$idname) && dp$idname %in% names(d)) length(unique(d[[dp$idname]])) else NULL,
    columns = names(d), basis = "actual did DIDparams$data after preprocessing"
  )
}; extract_fit <- function(fit, alp = 0.05, point_only_rcs = FALSE) {
  J <- length(fit$att); N <- as.double(fit$n); V <- fit[["V_analytical"]]
  se_ana <- if (!is.null(V) && is.matrix(V) && nrow(V) == J && is.finite(N) && N > 0) {
    sqrt(pmax(0.0, diag(V)) / N)
  } else {
    rep(NA_real_, J)
  }; se_reported <- if (isTRUE(point_only_rcs)) rep(NA_real_, J) else as.double(fit$se)
  c_val <- if (isTRUE(point_only_rcs)) NA_real_ else if (!is.null(fit$c) && length(fit$c) == 1L && is.finite(fit$c)) {
    as.double(fit$c)
  } else {
    stats::qnorm(1.0 - alp / 2.0)
  }; point_c <- stats::qnorm(1.0 - alp / 2.0); cells <- list(
    group = as.double(fit$group), t = as.double(fit$t), att = as.double(fit$att), se = se_reported,
    se_analytical = as.double(se_ana), c = c_val, c_pointwise = stats::qnorm(1.0 - alp / 2.0)
  ); cells$ci_lower <- cells$att - point_c * se_reported; cells$ci_upper <- cells$att + point_c * se_reported
  cells$band_lower <- cells$att - c_val * se_reported; cells$band_upper <- cells$att + c_val * se_reported
  cells <- lapply(cells, I)
  ch <- if (inherits(fit, "cdid_MP")) fit$chained else NULL; if (!is.null(ch)) {
    for (nm in c("available", "support_status", "support_reason", "role")) {
      if (!is.null(ch[[nm]]) && length(ch[[nm]]) == J) cells[[nm]] <- I(as.vector(ch[[nm]]))
    }
  }; num1 <- function(v) if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.double(v) else NULL; list(
    cells = cells, n_cells = J, cband = isTRUE(fit$DIDparams$cband), pretest = list(
      W = num1(fit$W), Wpval = num1(fit$Wpval),
      chained_pretest_pval = if (!is.null(ch) && is.list(ch$pretest)) num1(ch$pretest$Wpval) else NULL,
      status = if (!is.null(ch) && is.list(ch$pretest) && !is.null(ch$pretest$status)) {
        as.character(ch$pretest$status)
      } else if (!is.null(fit$Wpval) && is.finite(fit$Wpval)) {
        "computed"
      } else {
        "unavailable"
      }, reason = if (!is.null(ch) && is.list(ch$pretest)) ch$pretest$reason else NULL
    ), class = class(fit)
  )
}; extract_agg <- function(a, alp = 0.05) {
  c_overall <- stats::qnorm(1.0 - alp / 2.0)
  c_egt <- if (!is.null(a$crit.val.egt) && length(a$crit.val.egt) == 1L && is.finite(a$crit.val.egt)) {
    as.double(a$crit.val.egt)
  } else {
    c_overall
  }; ov <- if (!is.null(a$overall.att)) as.double(a$overall.att) else NA_real_
  ose <- if (!is.null(a$overall.se)) as.double(a$overall.se) else NA_real_
  av <- as.double(a$att.egt); ase <- as.double(a$se.egt)
  list(overall_att = ov, overall_se = ose, overall_c = c_overall,
    overall_ci_lower = ov - c_overall * ose, overall_ci_upper = ov + c_overall * ose,
    egt = I(as.double(a$egt)), att_egt = I(av), se_egt = I(ase), crit_val_egt = c_egt,
    ci_lower_egt = I(av - c_overall * ase), ci_upper_egt = I(av + c_overall * ase),
    band_lower_egt = I(av - c_egt * ase), band_upper_egt = I(av + c_egt * ase),
    cband = isTRUE(a$DIDparams$cband), na_rm = isTRUE(a$DIDparams$na.rm)
  )
}; extract_pretest <- function(pre) {
  ch <- pre$chained; num1 <- function(v) if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.double(v) else NA_real_; list(
    CvM = num1(pre$CvM), CvMcval = num1(pre$CvMcval), CvMpval = num1(pre$CvMpval),
    status = if (!is.null(ch) && !is.null(ch$status)) as.character(ch$status) else if (is.finite(num1(pre$CvMpval))) "computed" else "unavailable",
    reason = if (!is.null(ch) && !is.null(ch$reason)) as.character(ch$reason) else NULL
  )
}; classify_disposition <- function(task, cap, extracted) {
  call_type <- task$call_type; exp_disp <- if (!is.null(task$expected_disposition)) as.character(task$expected_disposition) else ""
  exp_cond <- if (!is.null(task$expected_condition) && nzchar(task$expected_condition)) as.character(task$expected_condition) else "cdid_error"
  if (!is.null(cap$error) || is.null(cap$value)) {
    err_classes <- if (!is.null(cap$error)) class(cap$error) else character()
    if (call_type %in% c("DID_CC", "DID_UNBALANCED", "DID_RCS")) return("unavailable-comparator")
    if (identical(exp_disp, "expected-refusal")) {
      if (exp_cond %in% err_classes || (identical(exp_cond, "cdid_error") && "cdid_error" %in% err_classes)) {
        return("expected-refusal")
      }; return("unexpected-error")
    }; if (grepl("refusal", exp_disp, ignore.case = TRUE) && "cdid_error" %in% err_classes) return("expected-refusal")
    return("unexpected-error")
  }; if (identical(exp_disp, "expected-refusal")) return("unexpected-success"); if (identical(call_type, "PRETEST")) {
    st <- extracted$pretest_result$status; if (identical(st, "not_applicable")) return("not-applicable")
    if (identical(st, "unavailable") || is.na(extracted$pretest_result$CvMpval)) return("partial-support"); return("success")
  }; if (identical(call_type, "AGG")) {
    ag <- extracted$agg_result; if (is.null(ag) || is.na(ag$overall_att) || (!is.null(ag$att_egt) && anyNA(ag$att_egt))) {
      return("partial-support")
    }; return("success")
  }
  # FIT or DID_*
  att <- extracted$fit_result$cells$att; if (is.null(att) || length(att) == 0L) {
    return(if (call_type %in% c("DID_CC", "DID_UNBALANCED", "DID_RCS")) "unavailable-comparator" else "unexpected-error")
  }; prim_coords <- task$primary_coordinates; if (!is.null(prim_coords) && length(prim_coords) > 0L) {
    g_vec <- extracted$fit_result$cells$group; t_vec <- extracted$fit_result$cells$t
    keys <- paste0("gt:", vapply(g_vec, function(x) format(x, trim = TRUE, scientific = FALSE), character(1L)), ":", vapply(t_vec, function(x) format(x, trim = TRUE, scientific = FALSE), character(1L)))
    idx <- match(unlist(prim_coords), keys); if (anyNA(idx) || any(!is.finite(att[idx]))) {
      if (call_type %in% c("DID_CC", "DID_UNBALANCED", "DID_RCS") && all(!is.finite(att[idx[!is.na(idx)]]))) {
        return("unavailable-comparator")
      }; return("partial-support")
    }
  } else if (any(!is.finite(att))) {
    return("partial-support")
  }; "success"
}; run_freeze <- function(cfg_path, script_dir) {
  cfg <- jsonlite::read_json(cfg_path, simplifyVector = FALSE)
  if (is.null(cfg$source_sha) || !nzchar(cfg$source_sha) || identical(cfg$source_sha, "FREEZE_SOURCE_SHA")) {
    stop("Refusing freeze: FREEZE_SOURCE_SHA placeholder must be replaced with the actual source commit SHA.", call. = FALSE)
  }; source(file.path(script_dir, "generate.R"), local = FALSE); idn <- setup_and_assert(unlist(cfg$lib_paths), repo_path = cfg$repo_path, formals_path = cfg$formals_path, source_sha = cfg$source_sha); if (!isTRUE(idn$ok)) {
    write_json_atomic(list(ok = FALSE, checks = idn$checks, runtime = idn$runtime), cfg$runtime_out_json)
    stop("Runtime package/version guard failed during freeze.", call. = FALSE)
  }; states <- expand_rng_streams(n_streams = 59L, n_reps = 100L, n_purposes = 5L, seed = 61012026L)
  save_rds_atomic(states, cfg$seeds_rds_path); out <- list(
    ok = TRUE, frozen_utc = now_utc(), source_sha = cfg$source_sha, seed_master = 61012026L,
    rng_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"), n_streams = 59L, n_reps = 100L, n_purposes = 5L,
    n_states = length(states), seeds_rds_hash = file_hash(cfg$seeds_rds_path), checks = idn$checks,
    runtime = idn$runtime, profile_digest = idn$profile_digest
  ); write_json_atomic(out, cfg$runtime_out_json); say("FREEZE R RUNTIME AND 29500 RNG STATES COMPLETE")
}; run_task <- function(task_path, script_dir) {
  t_startup0 <- proc.time(); options(warn = 1); task <- jsonlite::read_json(task_path, simplifyVector = FALSE)
  source(file.path(script_dir, "generate.R"), local = FALSE); idn <- tryCatch(
    setup_and_assert(unlist(task$lib_paths), frozen_runtime = task$frozen_runtime, repo_path = task$repo_path, formals_path = task$formals_path, source_sha = task$source_sha),
    error = function(e) list(ok = FALSE, error = cond_rec(e))
  ); startup_secs <- secs(t_startup0); if (!isTRUE(idn$ok)) {
    write_json_atomic(list(
      attempt_id = task$attempt_id, task_id = task$task_id, identity_ok = FALSE, checks = idn$checks,
      error = idn$error, startup_seconds = startup_secs$elapsed
    ), task$worker_out_json, overwrite = TRUE); quit(save = "no", status = 3L)
  }; if (!identical(file_hash(task$seeds_rds_path), task$frozen_runtime$seeds_rds_hash)) {
    stop("Frozen seed file hash mismatch", call. = FALSE)
  }
  states <- readRDS(task$seeds_rds_path); call_type <- as.character(task$call_type)
  alp <- if (!is.null(task$alp)) as.double(task$alp) else 0.05
  # ---- AGG child operation: load persisted parent FIT RDS without refitting ----
  if (identical(call_type, "AGG")) {
    t_data0 <- proc.time(); parent_res <- tryCatch(
      list(ok = TRUE, fit = readRDS(task$parent_rds_path)), error = function(e) list(ok = FALSE, error = cond_rec(e))
    ); data_secs <- secs(t_data0); if (!isTRUE(parent_res$ok) || is.null(parent_res$fit)) {
      rec <- list(
        attempt_id = task$attempt_id, task_id = task$task_id, call_type = call_type, identity_ok = TRUE,
        profile_digest = idn$profile_digest, disposition = "unexpected-error", startup_seconds = startup_secs$elapsed,
        data_seconds = data_secs$elapsed, data_generated_now = FALSE, operation_seconds = 0.0,
        error = if (!is.null(parent_res$error)) parent_res$error else list(class = "serialization_error", message = "Failed to read parent RDS"),
        warnings = list(), n_warnings = 0L
      ); write_json_atomic(rec, task$worker_out_json, overwrite = TRUE); return(invisible(NULL))
    }; restore_rng_state(states[[as.integer(task$op_rng_index)]]); call_env <- new.env(parent = globalenv())
    assign("FIT", parent_res$fit, envir = call_env); expr <- parse(text = task$literal_call); t_op0 <- proc.time()
    cap <- run_captured(function() eval(expr, envir = call_env)); op_secs <- secs(t_op0); extracted <- list()
    if (is.null(cap$error) && !is.null(cap$value)) {
      extracted$agg_result <- extract_agg(cap$value, alp = alp)
      extracted$agg_result$na_rm <- grepl("na.rm=TRUE", task$literal_call, fixed = TRUE)
      extracted$agg_result$cband <- isTRUE(parent_res$fit$DIDparams$cband)
      extracted$agg_result$support_changed <- extracted$agg_result$na_rm && any(!is.finite(parent_res$fit$att))
    }; disp <- classify_disposition(task, cap, extracted); rec <- list(
      attempt_id = task$attempt_id, task_id = task$task_id, call_type = call_type, identity_ok = TRUE,
      profile_digest = idn$profile_digest, disposition = disp, startup_seconds = startup_secs$elapsed,
      data_seconds = data_secs$elapsed, data_generated_now = FALSE, operation_seconds = op_secs$elapsed,
      operation_user_seconds = op_secs$user, operation_sys_seconds = op_secs$system,
      error = cond_rec(cap$error), warnings = cap$warnings, n_warnings = length(cap$warnings), agg_result = extracted$agg_result
    ); write_json_atomic(rec, task$worker_out_json, overwrite = TRUE); return(invisible(NULL))
  }
  # ---- FIT / PRETEST / DID_* operation: load or generate cached dataset ----
  t_data0 <- proc.time(); data_rds <- task$data_rds_path; data_json <- task$data_json_path; generated_now <- FALSE
  if (xor(file.exists(data_rds), file.exists(data_json))) stop("Incomplete cached data pair; reconcile the preserved artifacts before resuming.", call. = FALSE)
  if (file.exists(data_rds) && file.exists(data_json)) {
    D <- readRDS(data_rds); dmeta <- jsonlite::read_json(data_json, simplifyVector = FALSE)
    if (!identical(obj_md5(D$PANEL), dmeta$panel_md5) || !identical(obj_md5(D$F), dmeta$frame_md5)) {
      stop("Cached dataset MD5 mismatch for ", data_rds, call. = FALSE)
    }
  } else {
    generated_now <- TRUE; pop_st <- states[[as.integer(task$pop_rng_index)]]; obs_st <- states[[as.integer(task$obs_rng_index)]]
    D <- generate_recipe_draw(task$recipe_id, pop_st, obs_st)
    dmeta <- c(D$meta, list(panel_md5 = obj_md5(D$PANEL), frame_md5 = obj_md5(D$F), created_utc = now_utc()))
    save_rds_atomic(D, data_rds); write_json_atomic(dmeta, data_json)
  }; data_secs <- secs(t_data0); restore_rng_state(states[[as.integer(task$op_rng_index)]])
  seed_before <- get(".Random.seed", envir = .GlobalEnv); call_env <- new.env(parent = globalenv())
  point_only_rcs <- identical(call_type, "DID_RCS") || isTRUE(task$point_only_rcs)
  panel_in <- if (point_only_rcs) D$PANEL[!is.na(D$PANEL$y), , drop = FALSE] else D$PANEL; rownames(panel_in) <- NULL
  assign("PANEL", panel_in, envir = call_env); assign("F", D$F, envir = call_env)
  assign("J_REP", "Within each cohort, membership in the fitting set is independent of the scientific conditioning variables by design.", envir = call_env)
  assign("J_TRANSPORT", "Column holds the known normalized ratio f(X|G,U)/f(X|G) computed from the design law with population cohort means; it is not estimated.", envir = call_env)
  md5_panel_before <- obj_md5(get("PANEL", envir = call_env)); md5_frame_before <- obj_md5(get("F", envir = call_env))
  expr <- parse(text = task$literal_call); t_op0 <- proc.time(); cap <- run_captured(function() eval(expr, envir = call_env))
  op_secs <- secs(t_op0); input_unchanged <- identical(md5_panel_before, obj_md5(get("PANEL", envir = call_env))) &&
    identical(md5_frame_before, obj_md5(get("F", envir = call_env)))
  if (!isTRUE(input_unchanged)) {
    stop("Input PANEL or F mutated during operation ", task$task_id, call. = FALSE)
  }; extracted <- list(); parent_saved <- FALSE; if (is.null(cap$error) && !is.null(cap$value)) {
    if (identical(call_type, "PRETEST")) {
      extracted$pretest_result <- extract_pretest(cap$value)
    } else {
      extracted$fit_result <- extract_fit(cap$value, alp = alp, point_only_rcs = point_only_rcs)
      extracted$sample <- tryCatch(sample_record(cap$value), error = function(e) list(error = cond_rec(e))); extracted$if_audit <- tryCatch(
        inference_record(
          cap$value, task$package, D$meta$n_units, compute_inffunc = !isFALSE(task$compute_inffunc), point_only_rcs = point_only_rcs, original_frame = D$F
        ), error = function(e) list(error = cond_rec(e))
      ); if (isTRUE(task$save_parent_rds) && !is.null(task$parent_rds_path) && nzchar(task$parent_rds_path)) {
        if (file.exists(task$parent_rds_path)) {
          if (!identical(obj_md5(readRDS(task$parent_rds_path)), obj_md5(cap$value))) stop("Preserved parent differs from interrupted replacement", call. = FALSE)
        } else save_rds_atomic(cap$value, task$parent_rds_path)
        parent_saved <- TRUE
      }
    }
  }; disp <- classify_disposition(task, cap, extracted); rec <- list(
    attempt_id = task$attempt_id, task_id = task$task_id, call_type = call_type, package = task$package,
    identity_ok = TRUE, profile_digest = idn$profile_digest, disposition = disp,
    startup_seconds = startup_secs$elapsed, data_seconds = data_secs$elapsed,
    data_generated_now = generated_now, operation_seconds = op_secs$elapsed,
    operation_user_seconds = op_secs$user, operation_sys_seconds = op_secs$system,
    rng_consumed = !identical(seed_before, get(".Random.seed", envir = .GlobalEnv)),
    input_unchanged = input_unchanged, data_panel_md5 = dmeta$panel_md5, data_frame_md5 = dmeta$frame_md5,
    parent_saved = parent_saved, error = cond_rec(cap$error), warnings = cap$warnings,
    n_warnings = length(cap$warnings), fit_result = extracted$fit_result,
    pretest_result = extracted$pretest_result, sample = extracted$sample, if_audit = extracted$if_audit
  ); write_json_atomic(rec, task$worker_out_json, overwrite = TRUE)
}; main <- function() {
  args <- commandArgs(trailingOnly = TRUE); if (length(args) != 2L || !(args[1L] %in% c("--freeze", "--task"))) {
    stop("Usage: worker.R --freeze <freeze.json> | --task <task.json>", call. = FALSE)
  }; if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package is required by worker.R", call. = FALSE)
  }; all_args <- commandArgs(trailingOnly = FALSE); file_arg <- grep("^--file=", all_args, value = TRUE)
  script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(sub("^--file=", "", file_arg[1L]))) else getwd()
  if (identical(args[1L], "--freeze")) {
    run_freeze(args[2L], script_dir)
  } else {
    run_task(args[2L], script_dir)
  }
}; if (sys.nframe() == 0L && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  main()
}
