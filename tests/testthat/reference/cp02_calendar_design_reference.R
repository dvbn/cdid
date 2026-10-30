# Self-contained independent reference for the private CP02 calendar/anticipation
# design and six-role endpoint representation (Proposition prop:calendar in
# docs/theory/cp02_moment_catalog.tex and Section 3 of
# docs/theory/cp02_moment_catalog_contract.md).
# Uses base R only; never calls the cdid namespace.

# Exact real-number comparison of (a + b) vs c via hexadecimal float ("0x...p...")
# decomposition into an exact base-2 bit register over 2^[-1080 .. 1030].
ref_add_hex_float <- function(reg, x, mult = 1L) {
  if (x == 0 || mult == 0L) return(reg)
  sgn <- if (x < 0) -as.integer(mult) else as.integer(mult)
  s <- tolower(sprintf("%a", abs(as.double(x))))
  s <- sub("^0x", "", s)
  parts <- strsplit(s, "p", fixed = TRUE)[[1L]]
  mant <- parts[1L]
  exp2 <- as.integer(parts[2L])
  mp <- strsplit(mant, ".", fixed = TRUE)[[1L]]
  lead <- strtoi(mp[1L], base = 16L)
  frac <- if (length(mp) > 1L) mp[2L] else ""
  offset <- 1085L
  if (lead != 0L) {
    idx <- exp2 + offset
    reg[idx] <- reg[idx] + sgn * lead
  }
  if (nzchar(frac)) {
    chars <- strsplit(frac, "", fixed = TRUE)[[1L]]
    for (j in seq_along(chars)) {
      d <- strtoi(chars[j], base = 16L)
      if (d != 0L) {
        idx <- exp2 - 4L * j + offset
        reg[idx] <- reg[idx] + sgn * d
      }
    }
  }
  reg
}

ref_exact_cmp_sum <- function(a, b, c) {
  if (anyNA(c(a, b, c)) || is.nan(a) || is.nan(b) || is.nan(c) ||
      !is.finite(a) || !is.finite(b) || b < 0 || c == -Inf) {
    stop("Invalid non-finite cutoff inputs in reference.", call. = FALSE)
  }
  if (!is.finite(a + b) || (is.finite(c) && !is.finite(c - b))) {
    stop("Overflowing float cutoff in reference.", call. = FALSE)
  }
  if (is.infinite(c)) return(-1L)
  reg <- integer(2120L)
  reg <- ref_add_hex_float(reg, a, 1L)
  reg <- ref_add_hex_float(reg, b, 1L)
  reg <- ref_add_hex_float(reg, c, -1L)
  carry <- 0L
  for (i in seq_along(reg)) {
    val <- reg[i] + carry
    bit <- val %% 2L
    carry <- (val - bit) %/% 2L
    reg[i] <- bit
  }
  if (carry < 0L) return(-1L)
  if (carry > 0L) return(1L)
  nz <- which(reg != 0L)
  if (!length(nz)) 0L else 1L
}

ref_calendar_reference <- function(g, calendar, anticipation = 0) {
  cal <- sort(unique(as.double(calendar)))
  delta <- as.double(anticipation)
  vapply(as.double(g), function(gi) {
    ok <- vapply(cal, function(tm) ref_exact_cmp_sum(tm, delta, gi) < 0L, logical(1L))
    if (!any(ok)) stop("Empty pre-anticipation reference set in reference.", call. = FALSE)
    max(cal[ok])
  }, numeric(1L))
}

ref_eligible_controls <- function(g, u, cohorts, anticipation = 0) {
  coh <- sort(unique(as.double(cohorts)))
  delta <- as.double(anticipation)
  keep <- coh != as.double(g) & vapply(coh, function(k) ref_exact_cmp_sum(u, delta, k) < 0L, logical(1L))
  coh[keep]
}

ref_endpoint_role <- function(u, g, ref, anticipation, qualified) {
  if (u == ref) return("reference")
  cmp <- ref_exact_cmp_sum(u, anticipation, g)
  if (cmp < 0L) return("pre")
  if (u < g) {
    if (qualified) "anticipation" else "anticipation_nuisance"
  } else {
    if (qualified) "post" else "post_nuisance"
  }
}

# Exact integer/rational row-echelon rank and row-space membership for unimodular
# incidence matrices W (entries in {-1, 0, 1}) and integer/rational target vectors.
ref_gcd <- function(a, b) {
  a <- abs(as.numeric(a)); b <- abs(as.numeric(b))
  while (b > 0) {
    r <- a %% b
    a <- b
    b <- r
  }
  if (a == 0) 1 else a
}

ref_rowspace_check <- function(W, c_vec) {
  M <- rbind(as.matrix(W), as.numeric(c_vec))
  nr <- nrow(W); nc <- ncol(W)
  A <- matrix(as.numeric(M), nrow = nr + 1L, ncol = nc)
  piv_row <- 1L
  for (col in seq_len(nc)) {
    if (piv_row > nr) break
    cand <- piv_row:nr
    nz <- cand[which(abs(A[cand, col]) > 0)]
    if (!length(nz)) next
    p <- nz[1L]
    if (p != piv_row) {
      tmp <- A[piv_row, ]; A[piv_row, ] <- A[p, ]; A[p, ] <- tmp
    }
    piv <- A[piv_row, col]
    for (r in setdiff(seq_len(nr + 1L), piv_row)) {
      if (A[r, col] != 0) {
        A[r, ] <- A[r, ] * piv - A[piv_row, ] * A[r, col]
        g <- 0
        for (v in A[r, ]) if (v != 0) g <- if (g == 0) abs(v) else ref_gcd(g, v)
        if (g > 1) A[r, ] <- A[r, ] / g
      }
    }
    piv_row <- piv_row + 1L
  }
  list(rank = as.integer(piv_row - 1L), in_rowspace = all(A[nr + 1L, ] == 0))
}

# Canonical hand-verified reference fixtures for irregular calendar
# {2000, 2002, 2005, 2009}, g = 2009 under delta = 1 and delta = 4.
ref_irregular_2009_fixtures <- function() {
  calendar <- c(2000, 2002, 2005, 2009)
  list(
    calendar = calendar,
    delta1_ref = 2005,
    delta4_ref = 2002,
    delta4_anticipation_dates = 2005,
    # 6-role hand incidence matrix for g = 2009, delta = 4, ref = 2002:
    # Links:
    # 1: (2000, 2002) under op:K2007 (pre -> ref)
    # 2: (2002, 2005) under op:K2007 (ref -> anticipation_nuisance, since 2007 <= 2005 + 4)
    # 3: (2002, 2005) under op:K2010 (ref -> anticipation ATT(2009,2005), since 2010 > 2005 + 4)
    # 4: (2005, 2009) under op:K2010 (anticipation ATT(2009,2005) -> post_nuisance, since 2010 <= 2009 + 4)
    # 5: (2005, 2009) under op:KInf  (anticipation ATT(2009,2005) -> post ATT(2009,2009))
    # Parameters in canonical kind order:
    # 1: beta(g=2009,t=2000,fn=op:K2007)                          [pre]
    # 2: beta_anticipation_nuisance(g=2009,t=2005,fn=op:K2007)    [anticipation_nuisance]
    # 3: ATT(g=2009,t=2005)                                       [anticipation]
    # 4: beta_post_nuisance(g=2009,t=2009,fn=op:K2010)            [post_nuisance]
    # 5: ATT(g=2009,t=2009)                                       [post]
    W_six_role = matrix(
      c(-1L,  0L,  0L, 0L, 0L,
         0L,  1L,  0L, 0L, 0L,
         0L,  0L,  1L, 0L, 0L,
         0L,  0L, -1L, 1L, 0L,
         0L,  0L, -1L, 0L, 1L),
      nrow = 5L, ncol = 5L, byrow = TRUE
    ),
    roles_six_role = c("pre", "anticipation_nuisance", "anticipation", "post_nuisance", "post")
  )
}
