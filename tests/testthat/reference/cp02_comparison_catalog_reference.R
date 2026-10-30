# Independent brute-force Cartesian reference for the private CP02 approved
# comparison catalog (docs/theory/cp02_approved_comparison_catalog.tex).
# Base R only; never calls the cdid namespace. The exact real comparator
# `cmp(a, b, c)` = sign((a + b) - c) is injected by the caller (the hash-bound
# bit-register comparator of cp02_calendar_design_reference.R). Sets of cohorts
# are compared with identical() on sorted numeric vectors, never via strings.

ref_cat_count <- function(occ, k, s, t) {
  hit <- occ$cohort == k & occ$s == s & occ$t == t
  if (any(hit)) occ$count[hit][1L] else 0
}

ref_cat_reference <- function(g, cal, delta, cmp) {
  pre <- vapply(cal, function(d) cmp(d, delta, g) < 0L, logical(1L))
  if (any(pre)) max(cal[pre]) else NA_real_
}

# Every (target, pair, source) candidate over the full Cartesian product of
# calendar pairs and horizons, then exact grouping of identical pools.
ref_catalog_enumerate <- function(calendar, cohorts, targets, delta, occ, control_group, cmp) {
  cal <- sort(as.double(calendar)); coh <- sort(as.double(cohorts))
  rows <- list(); n_cand <- 0; n_nonempty <- 0
  for (g in sort(as.double(targets))) {
    ref <- ref_cat_reference(g, cal, delta, cmp)
    if (is.na(ref)) next
    for (s in cal) for (t in cal) {
      if (!(s < t) || ref_cat_count(occ, g, s, t) <= 0) next
      Ce <- coh[vapply(coh, function(k) ref_cat_count(occ, k, s, t) > 0, logical(1L))]
      pick <- function(cutoff) Ce[vapply(Ce, function(k) is.infinite(k) ||
        (k > g && cmp(cutoff, delta, k) < 0L), logical(1L))]
      if (control_group == "nevertreated") {
        cands <- list(list(original = TRUE, h = NA_real_, K = Ce[is.infinite(Ce)]))
      } else {
        hs <- cal[cal >= g & cal >= t]
        cands <- c(list(list(original = TRUE, h = NA_real_, K = pick(t))),
                   lapply(hs, function(h) list(original = FALSE, h = h, K = pick(h))))
      }
      n_cand <- n_cand + length(cands)
      found <- list()
      for (cd in cands) {
        if (!length(cd$K)) next
        n_nonempty <- n_nonempty + 1
        j <- which(vapply(found, function(f) identical(f$K, cd$K), logical(1L)))
        if (!length(j)) {
          found[[length(found) + 1L]] <- list(K = cd$K, original = FALSE, horizons = numeric())
          j <- length(found)
        }
        if (cd$original) found[[j]]$original <- TRUE else found[[j]]$horizons <- c(found[[j]]$horizons, cd$h)
      }
      for (f in found) rows[[length(rows) + 1L]] <- list(g = g, s = s, t = t, ref = ref, K = f$K,
        original = f$original, horizons = sort(f$horizons))
    }
  }
  list(rows = rows, n_candidates = n_cand, n_nonempty = n_nonempty)
}

# Conditional endpoint role of date u for target g, reference ref, pool K.
ref_cat_role <- function(u, g, ref, K, delta, cmp) {
  if (u == ref) return("reference")
  if (cmp(u, delta, g) < 0L) return("pre_private")
  clean <- all(vapply(K, function(k) is.infinite(k) || cmp(u, delta, k) < 0L, logical(1L)))
  paste0(if (u < g) "anticipation" else "post", if (clean) "_conditional" else "_private")
}

# Canonical comparison key: exact hex doubles of (g, s, t, K) for sorting only.
ref_cat_key <- function(g, s, t, K) {
  paste(sprintf("%a", g), sprintf("%a", s), sprintf("%a", t), paste(sprintf("%a", K), collapse = ","), sep = "|")
}
