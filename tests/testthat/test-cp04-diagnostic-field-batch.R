# Independent retained R payload/field definitions. No expected count is made
# by either native accounting entry point under test.
dfb_ns <- function(n) getFromNamespace(n, "cdid")
dfb_local_method <- function(name, value, local = parent.frame()) {
  had <- exists(name, globalenv(), inherits = FALSE)
  old <- get0(name, globalenv(), inherits = FALSE)
  withr::defer(if (had) assign(name, old, globalenv()) else rm(list = name, envir = globalenv()), envir = local)
  assign(name, value, globalenv())
}
dfb_ref_payload <- function(x, scan_attrs = TRUE) {
  if (is.environment(x) || is.null(x)) return(c(entries = 0, char_bytes = 0))
  ans <- c(entries = length(x) + 8, char_bytes = 0)
  if (is.list(x)) {
    for (v in x) ans <- ans + dfb_ref_payload(v)
  } else if (is.character(x)) {
    bytes <- nchar(x, type = "bytes", allowNA = TRUE)
    bytes[is.na(bytes)] <- 0
    ans <- ans + c(entries = sum(ceiling(bytes / 2)), char_bytes = sum(bytes))
  }
  if (scan_attrs && !is.null(attributes(x))) ans <- ans + dfb_ref_payload(attributes(x), FALSE)
  ans
}
dfb_ref_fields <- function(x, nms) {
  z <- 0
  for (nm in nms) z <- z + unname(dfb_ref_payload(x[[nm]])["entries"])
  z
}
dfb_native <- function(x, nms, env = asNamespace("cdid"))
  dfb_ns("cp04_diagnostic_field_native")(x, nms, env)

# Exercise admission plus the original scalar loop in an ordinary lexical
# environment. All test overrides live here, never in package bindings.
dfb_field_runner <- function(env = new.env(parent = asNamespace("cdid")), batch = TRUE) {
  f <- if (batch) function(x, nms) {
    ans <- cp04_diagnostic_field_native(x, nms, environment())
    if (!is.null(ans)) return(ans)
    z <- 0
    for (nm in nms) z <- z + unname(cp04_payload(x[[nm]])["entries"])
    z
  } else function(x, nms) {
    z <- 0
    for (nm in nms) z <- z + unname(cp04_payload(x[[nm]])["entries"])
    z
  }
  environment(f) <- env
  f
}

testthat::test_that("diagnostic batches equal independent field sums without deduplicating aliases", {
  shared <- structure(c(1.5, NA_real_, -0), note = list(label = "été", flag = TRUE))
  long <- c(strrep("long-label-", 2000L), NA_character_, "", "漢字")
  records <- list(
    list(a = 1, b = TRUE, c = NULL),
    list(a = shared, b = list(shared, shared, long), c = long),
    list(a = matrix(c(1, 2, 3, 4), 2L, dimnames = list(c("r1", "r2"), c("c1", "c2"))),
      b = structure(list(1, "x"), note = list(inner = shared)), c = raw(3L)),
    list(a = structure(new.env(parent = emptyenv()), class = "ignored_environment"),
      b = c(1 + 2i, NA_complex_), c = list(NULL, list(NULL))))
  for (x in records) {
    for (nms in list(c("a", "b", "c"), c("missing", "a", "a", "b", "missing"))) {
      expected <- dfb_ref_fields(x, nms)
      testthat::expect_identical(dfb_native(x, nms), expected)
      testthat::expect_identical(dfb_field_runner()(x, nms), expected)
      testthat::expect_type(dfb_native(x, nms), "double")
      testthat::expect_null(names(dfb_native(x, nms)))
    }
  }
  # First duplicate record name wins, every requested alias is charged again.
  x <- stats::setNames(list(1, long, shared), c("dup", "dup", "other"))
  testthat::expect_identical(dfb_native(x, c("dup", "dup", "absent")), 18)
  testthat::expect_identical(dfb_native(x, c("other", "dup")), dfb_ref_fields(x, c("other", "dup")))
  testthat::expect_identical(dfb_native(unname(x), "absent"), 0)
  testthat::expect_identical(dfb_native(list(a = NULL), "a"), 0)
  # Arbitrary long payload strings carry their complete bytes, with no new cap.
  testthat::expect_gt(dfb_native(list(a = long), "a"), 10000)
})

testthat::test_that("unsupported metadata declines before custom extraction or length callbacks", {
  events <- character()
  length_method <- function(x) {
    events <<- c(events, paste0("length:", attr(x, "label")))
    length(unclass(x))
  }
  extract_method <- function(x, i, ...) {
    events <<- c(events, paste0("extract:", i))
    .subset2(unclass(x), i)
  }
  dfb_local_method("length.dfb_value", length_method)
  dfb_local_method("[[.dfb_record", extract_method)
  x <- structure(list(a = structure(c(1, 2), class = "dfb_value", label = "a"),
    b = structure(c(3, 4), class = "dfb_value", label = "b")), class = "dfb_record")
  testthat::expect_null(dfb_native(x, c("a", "b")))
  testthat::expect_identical(events, character())
  old <- dfb_field_runner(batch = FALSE)(x, c("a", "b", "a")); old_events <- events
  events <- character()
  now <- dfb_field_runner()(x, c("a", "b", "a"))
  testthat::expect_identical(now, old)
  testthat::expect_identical(events, old_events)
  testthat::expect_identical(events, c("extract:a", "length:a", "extract:b", "length:b", "extract:a", "length:a"))
  # Plain container, unsupported later field: preflight still calls no method.
  events <- character(); plain <- unclass(x)
  testthat::expect_null(dfb_native(plain, c("a", "b")))
  testthat::expect_identical(events, character())
  testthat::expect_identical(dfb_field_runner()(plain, c("a", "b")), dfb_ref_fields(plain, c("a", "b")))
})

testthat::test_that("special attributes and nonplain operands retain the scalar path", {
  deep <- 1
  for (i in seq_len(40L)) deep <- list(deep)
  one_dim <- array(c(1, 2), 2L, dimnames = list(c("a", "b")))
  attr(one_dim, "names") <- c("raw_a", "raw_b")
  rn <- structure(list(1), row.names = c(NA_integer_, -1L))
  unsupported <- list(factor(c("a", "b")), data.frame(a = c(1, 2)), ~ x,
    expression(x), quote(f(x)), pairlist(a = 1), function(x) x, 1:7, deep, rn)
  for (value in unsupported) {
    x <- list(a = value)
    testthat::expect_null(dfb_native(x, "a"))
    testthat::expect_identical(dfb_field_runner()(x, "a"), dfb_ref_fields(x, "a"))
  }
  # A one-dimensional array can expose names through dimnames rather than a
  # raw names attribute. Either admission outcome must preserve that count.
  one <- dfb_native(list(a = one_dim), "a")
  if (!is.null(one)) testthat::expect_identical(one, dfb_ref_fields(list(a = one_dim), "a"))
  testthat::expect_identical(dfb_field_runner()(list(a = one_dim), "a"),
    dfb_ref_fields(list(a = one_dim), "a"))
  for (nms in list(character(), NULL, structure("a", note = "key"), "é", "", NA_character_, rep("a", 129L))) {
    x <- list(a = c(1, 2), "é" = "utf")
    testthat::expect_null(dfb_native(x, nms))
    testthat::expect_identical(dfb_field_runner()(x, nms), dfb_ref_fields(x, nms))
  }
  x <- stats::setNames(list(1), "é")
  testthat::expect_null(dfb_native(x, "absent"))
  testthat::expect_identical(dfb_field_runner()(x, "absent"), 0)
})

testthat::test_that("each record rechecks scalar helpers after a fallback callback", {
  env <- new.env(parent = asNamespace("cdid")); events <- character()
  dfb_local_method("length.dfb_mutate", function(x) {
    events <<- c(events, "mutate")
    assign("cp04_payload", function(x, scan_attrs = TRUE) c(entries = 123, char_bytes = 0), env)
    1L
  })
  run <- dfb_field_runner(env)
  x <- list(a = structure(1, class = "dfb_mutate"), b = 1)
  expected_first <- unname(dfb_ns("cp04_payload")(x$a)["entries"])
  # Reset the private override after measuring only the old first-field charge.
  rm("cp04_payload", envir = env); events <- character()
  testthat::expect_identical(run(x, c("a", "b")), expected_first + 123)
  testthat::expect_identical(events, "mutate")
  testthat::expect_null(dfb_native(list(a = 1), "a", env))
  testthat::expect_identical(run(list(a = 1), "a"), 123)
  # Each helper that changes scalar field behavior independently declines.
  for (nm in c("[[", "+", "[", "unname", "cp04_payload_native_kernel")) {
    e <- new.env(parent = asNamespace("cdid")); assign(nm, function(...) stop("override"), e)
    testthat::expect_null(dfb_native(list(a = 1), "a", e))
  }
})

testthat::test_that("original overflow refusals and large binary64 additions stay on fallback", {
  events <- 0L
  dfb_local_method("length.dfb_huge", function(x) {
    events <<- events + 1L
    attr(x, "declared_length")
  })
  huge <- structure(1, class = "dfb_huge", declared_length = 2^53 - 1)
  x <- list(a = 1, b = huge)
  testthat::expect_null(dfb_native(x, c("a", "b")))
  testthat::expect_identical(events, 0L)
  old <- tryCatch(dfb_field_runner(batch = FALSE)(x, c("a", "b")), error = identity)
  now <- tryCatch(dfb_field_runner()(x, c("a", "b")), error = identity)
  testthat::expect_s3_class(now, "cdid_gmm_bound_error")
  testthat::expect_identical(class(now), class(old))
  testthat::expect_identical(conditionMessage(now), conditionMessage(old))
  # Cheap synthetic length dispatch exercises old ordered additions beyond
  # 2^53 without allocating a giant vector. Native metadata declines it first.
  v <- structure(1, class = "dfb_huge", declared_length = 2^52 - 100)
  x <- list(a = v, b = v, c = v)
  events <- 0L
  testthat::expect_null(dfb_native(x, c("a", "b", "c")))
  testthat::expect_identical(events, 0L)
  testthat::expect_identical(dfb_field_runner()(x, c("a", "b", "c")),
    dfb_field_runner(batch = FALSE)(x, c("a", "b", "c")))
})

testthat::test_that("diagnostic batch uses its registered three-argument symbol", {
  dll <- base::getLoadedDLLs()[["cdid"]]
  reg <- base::getDLLRegisteredRoutines(dll, addNames = TRUE)[[".Call"]]
  testthat::expect_identical(as.integer(reg[["C_gmm_cp04_diagnostic_fields"]]$numParameters), 3L)
  resolved <- base::getNativeSymbolInfo("C_gmm_cp04_diagnostic_fields", PACKAGE = dll)
  kernel <- dfb_ns("cp04_diagnostic_field_native_kernel")()
  # getDLLRegisteredRoutines exposes a RegisteredNativeSymbol wrapper; the
  # resolver returns the NativeSymbol address from getNativeSymbolInfo.
  testthat::expect_identical(kernel, resolved$address)
  testthat::expect_identical(.Call(kernel, list(a = 1), "a", environment()), 9)
  # The common existing probe supplies authoritative arithmetic refusal tests;
  # the new entry calls the same gmm_require_valid_arithmetic gate.
  syms <- dfb_ns("gmm_require_native_dll")()
  testthat::expect_error(.Call(syms$C_gmm_test_arithmetic_probe, -1L, 0, 1, 1 + 2^-51, 0),
    class = "cdid_gmm_native_error")
})
