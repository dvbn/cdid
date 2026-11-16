# Test-only roundoff predicate for measured cross-platform floating-point leaves
# Only the listed leaves use numeric equivalence:
#  * bitwise-identical subtrees are accepted without scanning;
#  * every list node must keep identical type, length and attributes (names
#    including duplicates and order, class, dim, ...), and every leaf whose path
#    is not listed in `leaves` must be identical (`num.eq` as at the call site);
#  * a listed leaf must be an unclassed double with identical attributes,
#    identical nonfinite positions/values, identical zero and signed-zero masks,
#    and finite values within the greater of the absolute and relative bounds.
# Paths use the evidence notation: "$" is the compared object, "$/a/b" descends
# by element name, or by position for unnamed elements; duplicate names share a
# path. Returns a single TRUE/FALSE.
roundoff_identical <- function(x, y, leaves = character(), num.eq = TRUE,
                               abs_tolerance = 1e-12) {
  stopifnot(is.character(leaves), !anyNA(leaves), all(startsWith(leaves, "$")),
            is.logical(num.eq), length(num.eq) == 1L, !is.na(num.eq),
            length(abs_tolerance) == 1L, is.finite(abs_tolerance), abs_tolerance >= 0)
  same_attr <- function(a, b) {
    aa <- attributes(a); ab <- attributes(b)
    if (length(aa)) aa <- aa[order(names(aa))]
    if (length(ab)) ab <- ab[order(names(ab))]
    identical(aa, ab, num.eq = num.eq)
  }
  leaf_ok <- function(a, b) {
    if (!is.double(a) || !is.double(b) || is.object(a) || is.object(b) ||
        length(a) != length(b) || !same_attr(a, b)) return(FALSE)
    a <- as.vector(a); b <- as.vector(b)
    fin <- is.finite(a); za <- fin & a == 0
    identical(fin, is.finite(b)) && identical(a[!fin], b[!fin], num.eq = num.eq) &&
      identical(za, fin & b == 0) && identical(1 / a[za], 1 / b[za]) &&
      all(abs(a[fin] - b[fin]) <= pmax(abs_tolerance, 1e-12 * abs(a[fin]), 1e-12 * abs(b[fin])))
  }
  walk <- function(a, b, path) {
    if (identical(a, b, num.eq = FALSE)) return(TRUE)
    if (path %in% leaves) return(leaf_ok(a, b))
    if (!is.list(a) || !is.list(b)) return(identical(a, b, num.eq = num.eq))
    if (!identical(typeof(a), typeof(b)) || length(a) != length(b) || !same_attr(a, b)) return(FALSE)
    nm <- names(a)
    lab <- if (is.null(nm)) seq_along(a) else ifelse(nzchar(nm), nm, seq_along(a))
    kid <- paste0(path, "/", lab)
    for (i in seq_along(a)) if (!walk(.subset2(a, i), .subset2(b, i), kid[i])) return(FALSE)
    TRUE
  }
  walk(x, y, "$")
}
