# Development simulation support for cdid (not package API).
# Requires seed.R (with_seed).
#
# Minimal Monte Carlo harness. One replication per seed:
#   data <- generate(seed)
#   est  <- estimate(data)   # list with `estimate`: named numeric, one per target
#   inf  <- infer(est)       # list with `se`: named numeric; optional `crit`
#                            # (scalar or named numeric; default normal quantile)
# Each replication starts with_seed(seed, ...) in L'Ecuyer-CMRG. Data use
# that stream; estimation and inference use parallel::nextRNGStream() of its
# initial state, even when generate(seed) restores the caller's stream. Callbacks
# must use the supplied seed/with_seed contract or the current stream, without
# reseeding estimation. Replay restores these streams and leaves the caller's
# RNG state unchanged. R version and code still belong to the run provenance.
#
# Every attempted (replication, target) pair is kept. Errors are recorded by
# stage; warnings are recorded (text and count), never dropped. A non-finite
# or missing estimate is an estimation failure. A missing, non-finite or
# non-positive SE is an inference failure; it is not replaced. Nothing is
# redrawn.
#
# mc_summary() keeps failures in the denominators. coverage_unconditional is
# the share of ATTEMPTED replications whose interval was valid and covered the
# truth. coverage_conditional (covered / complete) is reported separately.
# bias, emp_sd and rmse use replications with a valid estimate (n_estimate),
# and mean_se and se_ratio use complete replications (n_complete). These
# conditional quantities are labelled as such. A summary is a development
# record, never a D7 certification.

mc_seeds <- function(n, first = 1L) {
  stopifnot(is.numeric(n), length(n) == 1L, is.finite(n), n >= 1,
            n == round(n), is.numeric(first), length(first) == 1L,
            is.finite(first), first == round(first),
            abs(first) <= .Machine$integer.max,
            first + n - 1 <= .Machine$integer.max)
  as.integer(first + seq_len(n) - 1)
}

mc_run <- function(seeds, generate, estimate, infer, truth, level = 0.95) {
  stopifnot(is.numeric(seeds), length(seeds) >= 1L,
            all(is.finite(seeds)), all(seeds == round(seeds)),
            all(abs(seeds) <= .Machine$integer.max), !anyDuplicated(seeds),
            is.numeric(truth), length(truth) > 0L, all(is.finite(truth)),
            !is.null(names(truth)), !anyNA(names(truth)),
            all(nzchar(names(truth))), !anyDuplicated(names(truth)),
            is.function(generate), is.function(estimate), is.function(infer),
            is.numeric(level), length(level) == 1L, is.finite(level),
            level > 0, level < 1)
  records <- lapply(seq_along(seeds), function(r) {
    mc_replicate(r, seeds[[r]], generate, estimate, infer, truth, level)
  })
  out <- do.call(rbind, records)
  rownames(out) <- NULL
  attr(out, "level") <- level
  out
}

# Re-run the replications for the given seeds exactly as mc_run() did.
mc_replay <- function(records, seeds, generate, estimate, infer) {
  level <- attr(records, "level")
  stopifnot(!is.null(level), all(seeds %in% records$seed))
  targets <- unique(records$target)
  truth <- stats::setNames(records$truth[match(targets, records$target)],
                           targets)
  reps <- unique(records[records$seed %in% seeds, c("rep", "seed")])
  out <- do.call(rbind, lapply(seq_len(nrow(reps)), function(i) {
    mc_replicate(reps$rep[i], reps$seed[i], generate, estimate, infer,
                 truth, level)
  }))
  rownames(out) <- NULL
  out
}

mc_replicate <- function(rep, seed, generate, estimate, infer, truth, level) {
  targets <- names(truth)
  warn <- character()
  record_warning <- function(stage) {
    function(w) {
      warn <<- c(warn, paste0(stage, ": ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  }
  run_stage <- function(stage, expr) {
    tryCatch(
      withCallingHandlers(list(value = expr, error = NA_character_),
                          warning = record_warning(stage)),
      error = function(e) list(value = NULL, error = conditionMessage(e))
    )
  }

  status <- rep("complete", length(targets))
  message <- rep(NA_character_, length(targets))
  est <- se <- rep(NA_real_, length(targets))
  crit <- rep(stats::qnorm(1 - (1 - level) / 2), length(targets))
  names(est) <- names(se) <- names(crit) <- targets

  t0 <- proc.time()[["elapsed"]]
  with_seed(seed, {
    estimation_stream <- parallel::nextRNGStream(get(".Random.seed", globalenv()))
    gen <- run_stage("generate", generate(seed))
    # A self-contained DGP may restore the initial data stream. Never reuse it
    # for bootstrap or other randomized estimation after data generation.
    assign(".Random.seed", estimation_stream, envir = globalenv())
    if (!is.na(gen$error)) {
      status[] <- "generation_failed"
      message[] <- gen$error
    } else {
      fit <- run_stage("estimate", {
        result <- estimate(gen$value)
        if (!is.list(result)) stop("estimate() must return a list")
        list(result = result, estimates = extract_named(result$estimate, targets))
      })
      if (!is.na(fit$error)) {
        status[] <- "estimation_failed"
        message[] <- fit$error
      } else {
        est[] <- fit$value$estimates
        bad_est <- !is.finite(est)
        status[bad_est] <- "estimation_failed"
        message[bad_est] <- "estimate missing or not finite"
        inf <- run_stage("infer", {
          result <- infer(fit$value$result)
          if (!is.list(result)) stop("infer() must return a list")
          result$se <- extract_named(result$se, targets)
          if (!is.null(result$crit)) {
            if (!is.numeric(result$crit)) stop("critical values must be numeric")
            if (length(result$crit) != 1L) {
              result$crit <- extract_named(result$crit, targets)
            }
          }
          result
        })
        if (!is.na(inf$error)) {
          status[!bad_est] <- "inference_failed"
          message[!bad_est] <- inf$error
        } else {
          se[] <- inf$value$se
          if (!is.null(inf$value$crit)) {
            crit[] <- if (length(inf$value$crit) == 1L) {
              inf$value$crit
            } else {
              inf$value$crit
            }
          }
          bad_se <- !bad_est & !(is.finite(se) & se > 0 & is.finite(crit) & crit > 0)
          status[bad_se] <- "inference_failed"
          message[bad_se] <- "standard error or critical value missing, not finite or not positive"
        }
      }
    }
  })
  elapsed <- proc.time()[["elapsed"]] - t0

  complete <- status == "complete"
  lower <- ifelse(complete, est - crit * se, NA_real_)
  upper <- ifelse(complete, est + crit * se, NA_real_)
  bad_bounds <- complete & (!is.finite(lower) | !is.finite(upper))
  status[bad_bounds] <- "inference_failed"
  message[bad_bounds] <- "interval bounds are not finite"
  lower[bad_bounds] <- upper[bad_bounds] <- NA_real_
  complete <- status == "complete"
  covered <- ifelse(complete, lower <= truth & truth <= upper, NA)
  data.frame(
    rep = rep, seed = seed, target = targets, truth = unname(truth),
    estimate = unname(est), se = unname(ifelse(status == "estimation_failed",
                                               NA_real_, se)),
    critical_value = unname(crit), lower = unname(lower),
    upper = unname(upper), covered = unname(covered),
    status = status, message = message, n_warnings = length(warn),
    warnings = if (length(warn)) paste(warn, collapse = " | ") else NA_character_,
    elapsed = elapsed, stringsAsFactors = FALSE
  )
}

extract_named <- function(x, targets) {
  if (is.null(x)) return(rep(NA_real_, length(targets)))
  if (!is.numeric(x) || is.null(names(x)) || anyNA(names(x)) ||
      anyDuplicated(names(x))) {
    stop("returned target values must be uniquely named numeric values")
  }
  unname(x[match(targets, names(x))])
}

mc_summary <- function(records) {
  by_target <- split(records, factor(records$target,
                                     levels = unique(records$target)))
  out <- do.call(rbind, lapply(by_target, function(d) {
    ok_est <- is.finite(d$estimate)
    comp <- d$status == "complete"
    err <- d$estimate[ok_est] - d$truth[ok_est]
    n_att <- nrow(d)
    n_cov <- sum(d$covered %in% TRUE)
    p_unc <- n_cov / n_att
    sd_comp <- if (sum(comp) >= 2L) stats::sd(d$estimate[comp]) else NA_real_
    mean_se <- if (any(comp)) mean(d$se[comp]) else NA_real_
    data.frame(
      target = d$target[1L], truth = d$truth[1L],
      n_attempted = n_att,
      n_generation_failed = sum(d$status == "generation_failed"),
      n_estimation_failed = sum(d$status == "estimation_failed"),
      n_inference_failed = sum(d$status == "inference_failed"),
      n_complete = sum(comp),
      n_with_warnings = sum(d$n_warnings > 0L),
      n_estimate = sum(ok_est),
      bias = if (any(ok_est)) mean(err) else NA_real_,
      emp_sd = if (sum(ok_est) >= 2L) stats::sd(d$estimate[ok_est]) else NA_real_,
      rmse = if (any(ok_est)) sqrt(mean(err^2)) else NA_real_,
      mean_se = mean_se,
      se_ratio = mean_se / sd_comp,
      n_covered = n_cov,
      coverage_unconditional = p_unc,
      mcse_coverage_unconditional = sqrt(p_unc * (1 - p_unc) / n_att),
      coverage_conditional = if (any(comp)) n_cov / sum(comp) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  attr(out, "level") <- attr(records, "level")
  attr(out, "run_type") <- paste0(
    "development run (", length(unique(records$rep)), " replications): ",
    "not a D7 certification")
  out
}
