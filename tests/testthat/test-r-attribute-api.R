# SPDX-License-Identifier: GPL-3.0-only
#
# Focused regression test for R >= 4.6 raw attribute iteration compatibility in
# `C_gmm_cp08_payload` (`src/cdid.c`): version-aware raw attribute lookup, count,
# value iteration (`R_mapAttrib` on R >= 4.6, `ATTRIB` loop on older R) and
# raw-attribute existence predicate (`ANY_ATTRIB` on R >= 4.5, `ATTRIB` on older R).

library(testthat)
if (!requireNamespace("cdid", quietly = TRUE)) stop("Actual installed cdid required; no source fallback.", call. = FALSE)
api <- asNamespace("cdid")

attr_api_ref_payload <- function(x, scan_attrs = TRUE) {
  if (is.environment(x) || is.null(x)) return(c(entries = 0, char_bytes = 0))
  ans <- c(entries = length(x) + 8, char_bytes = 0)
  if (is.list(x)) {
    for (v in x) ans <- ans + attr_api_ref_payload(v)
  } else if (is.character(x)) {
    bytes <- nchar(x, type = "bytes", allowNA = TRUE)
    bytes[is.na(bytes)] <- 0
    ans <- ans + c(entries = sum(ceiling(bytes / 2)), char_bytes = sum(bytes))
  }
  if (scan_attrs && !is.null(attributes(x))) ans <- ans + attr_api_ref_payload(attributes(x), FALSE)
  ans
}

test_that("C_gmm_cp08_payload matches the recursive R reference across plain/attributed nested lists, compact rownames, 1D arrays, tagged pairlists/calls, and Unicode", {
  v <- c(1.25, -3.5, NA_real_)
  utf <- c("ascii", "", NA_character_, "été", "qß", "αβγ", "漢字", "\U0001F600")
  aliased <- list(v, v, list(utf, v), utf)
  nested_attr <- structure(
    list(
      plain = list(1:3, list(), list(NULL, "plain")),
      shared = aliased,
      tagged = structure(list(a = v, b = utf), note = "inner-note", meta = list(k = 1:2, label = "漢字"))
    ),
    class = c("cdid_attr_test_outer", "list"),
    stamp = structure(c("été", "qß"), tag_attr = 10:12)
  )

  df_base <- data.frame(a = 1:5, b = c("u", "v", NA_character_, "ü", "漢字"), stringsAsFactors = FALSE)
  df_compact_auto <- df_base
  attr(df_compact_auto, "row.names") <- c(NA_integer_, -5L)
  df_compact_pos <- df_base
  attr(df_compact_pos, "row.names") <- c(NA_integer_, 5L)
  df_compact_zero <- df_base[0, , drop = FALSE]
  attr(df_compact_zero, "row.names") <- c(NA_integer_, 0L)
  df_char_rn <- df_base
  rownames(df_char_rn) <- c("r1", "r2", "r3", "r4", "r5_é")
  df_int_rn <- df_base
  attr(df_int_rn, "row.names") <- 101:105

  arr_1d_named <- array(1:4, dim = 4L, dimnames = list(c("n1", "n2_é", "n3", "n4")))
  arr_1d_unlisted_dn <- array(1:3, dim = 3L, dimnames = list(NULL))
  arr_1d_plain <- array(1:3, dim = 3L)
  arr_1d_with_raw_names <- structure(1:3, names = c("raw1", "raw2", "raw3"), dim = 3L,
    dimnames = list(c("dn1_é", "dn2", "dn3")))
  mat_2d <- matrix(1:6, nrow = 2L, dimnames = list(c("r1", "r2"), c("c1", "c2", "c3")))

  pl_tagged <- pairlist(a = 1L, 2.5, b = "three", "été" = utf)
  pl_untagged <- pairlist(1L, "two", NULL)
  pl_attributed <- structure(pairlist(a = 1:2, b = "qß"), note = "pairlist-attr", extra = list(1L, "αβγ"))
  call_tagged <- quote(f(x, y = 2, z = "été"))
  call_untagged <- quote(f(1, 2, 3))
  call_attributed <- structure(quote(g(a = 1, 2)), srcref = list("line1", 10L))

  cases <- list(
    NULL,
    new.env(parent = emptyenv()),
    structure(new.env(parent = emptyenv()), note = "ignored-env-attr"),
    v,
    utf,
    aliased,
    nested_attr,
    df_base,
    df_compact_auto,
    df_compact_pos,
    df_compact_zero,
    df_char_rn,
    df_int_rn,
    arr_1d_named,
    arr_1d_unlisted_dn,
    arr_1d_plain,
    arr_1d_with_raw_names,
    mat_2d,
    pl_tagged,
    pl_untagged,
    pl_attributed,
    call_tagged,
    call_untagged,
    call_attributed
  )

  for (i in seq_along(cases)) {
    x <- cases[[i]]
    expect_identical(api$cp04_payload(x, scan_attrs = TRUE), attr_api_ref_payload(x, TRUE), info = paste("scan_attrs=TRUE case", i))
    expect_identical(api$cp04_payload(x, scan_attrs = FALSE), attr_api_ref_payload(x, FALSE), info = paste("scan_attrs=FALSE case", i))
  }

  # Compact row names remain unexpanded in place after payload accounting.
  expect_identical(.row_names_info(df_compact_auto, type = 0L), c(NA_integer_, -5L))
  expect_identical(.row_names_info(df_compact_pos, type = 0L), c(NA_integer_, 5L))
})

test_that("C_gmm_cp08_payload preserves scan_attrs flag coercion and has-attributes short-circuit semantics", {
  plain_vec <- c(1, 2, 3)
  plain_untagged_pl <- pairlist(1, 2)
  attr_vec <- structure(c(1, 2, 3), label = "tagged")
  tagged_pl <- pairlist(a = 1, 2)
  flags <- list(TRUE, FALSE, 1, 0, -2.5, NA, NA_integer_, NULL, character(0), c(TRUE, FALSE), "not_logical")

  for (x in list(plain_vec, plain_untagged_pl, attr_vec, tagged_pl)) {
    for (flag in flags) {
      ref <- tryCatch(suppressWarnings(attr_api_ref_payload(x, flag)), error = identity)
      got <- tryCatch(suppressWarnings(api$cp04_payload(x, flag)), error = identity)
      if (inherits(ref, "error")) {
        expect_s3_class(got, "error")
        expect_identical(conditionMessage(got), conditionMessage(ref))
      } else {
        expect_identical(got, ref)
      }
    }
  }
})

test_that("attribute visitor preserves GC protection across callback R evaluations and allocations", {
  had_method <- exists("length.cdid_attr_gc_probe", envir = globalenv(), inherits = FALSE)
  old_method <- get0("length.cdid_attr_gc_probe", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_method) {
      assign("length.cdid_attr_gc_probe", old_method, envir = globalenv())
    } else if (exists("length.cdid_attr_gc_probe", envir = globalenv(), inherits = FALSE)) {
      rm("length.cdid_attr_gc_probe", envir = globalenv())
    }
  }, add = TRUE)

  calls <- 0L
  assign("length.cdid_attr_gc_probe", function(x) {
    calls <<- calls + 1L
    scratch <- replicate(16L, paste0("alloc_", seq_len(32L)), simplify = FALSE)
    gc(full = FALSE)
    if (length(scratch) < 0L) stop("unreachable")
    3L
  }, envir = globalenv())

  probe_a <- structure(list("alpha", "été"), class = "cdid_attr_gc_probe", inner = "first")
  probe_b <- structure(1:5, class = "cdid_attr_gc_probe", inner = list("beta", "漢字"))
  carrier <- structure(
    list(1:4, structure(list("x"), class = "cdid_attr_gc_probe")),
    first_attr = probe_a,
    second_attr = list(nested = probe_b, tail = c("u", "v", "qß")),
    third_attr = "after-gc-allocations"
  )

  got <- api$cp04_payload(carrier, scan_attrs = TRUE)
  calls_native <- calls
  ref <- attr_api_ref_payload(carrier, scan_attrs = TRUE)
  expect_gt(calls_native, 0L)
  expect_identical(got, ref)
})
