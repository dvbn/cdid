# Tiny native metadata fixtures: no fitted models or package simulations.
r_api_kernel <- function(name) {
  dll <- base::getLoadedDLLs()[["cdid"]]
  base::getNativeSymbolInfo(name, PACKAGE = dll)$address
}

testthat::test_that("coarse helper admission reads forced bindings without callbacks", {
  create <- r_api_kernel("C_gmm_reconstruction_coarse_create")
  close <- r_api_kernel("C_gmm_reconstruction_coarse_pullback")
  view <- list(cache = list(), position = integer(), seen = logical())
  state <- list(cells = list(), groups = list(), cids = integer(),
    B = 1L, N = 1L, intercept = TRUE)
  expected <- base::identity
  issue <- function(owner) .Call(create, view, state, NULL,
    list(identity = expected, .owner = owner), 1e8)
  release <- function(ptr) .Call(close, ptr, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
  check_value <- function(owner) {
    ptr <- issue(owner)
    testthat::expect_type(ptr, "externalptr")
    if (typeof(ptr) == "externalptr") release(ptr)
  }
  owner <- new.env(parent = baseenv())
  check_value(owner) # Absent local binding still uses the base frame.
  owner$identity <- expected
  check_value(owner)

  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L
  rm("identity", envir = owner)
  delayedAssign("identity", { counter$calls <- counter$calls + 1L; expected },
    eval.env = environment(), assign.env = owner)
  # A present delayed binding shadows base identity; it must not fall through.
  testthat::expect_null(issue(owner))
  testthat::expect_identical(counter$calls, 0L)
  testthat::expect_identical(get("identity", owner, inherits = FALSE), expected)
  check_value(owner) # The same promise is now forced and remains admissible.
  testthat::expect_identical(counter$calls, 1L)

  rm("identity", envir = owner)
  makeActiveBinding("identity", function(value) {
    counter$calls <- counter$calls + 1L
    expected
  }, owner)
  testthat::expect_null(issue(owner))
  testthat::expect_identical(counter$calls, 1L)
  missing <- (function(identity) environment())()
  testthat::expect_null(issue(missing))
})

testthat::test_that("zero-template admission retains raw attribute and closure predicates", {
  kernel <- r_api_kernel("C_gmm_reconstruction_zero_run")
  fit <- function(z) z + 1
  environment(fit) <- baseenv()
  group <- list(cell_ids = 1L, p = 1L,
    fit = structure(fit, note = list(label = "fit")))
  cells <- list(list(value = structure(c(1, 2), note = list(label = "cell"))))
  copy <- function(x) unserialize(serialize(x, NULL))
  template <- list(lambda = 0, nu = 0, group = copy(group), cells = copy(cells), N = 2L)
  run <- function(template, group, cells) {
    .Call(kernel, list(template), list(group), cells, integer(), 2L, 1L)
  }
  good <- run(template, group, cells)
  testthat::expect_identical(good, list(cell_ids = 1L, lambda = 0,
    group_ids = 1L, nu = list(0), next_gid = 2L))
  decorated <- template
  attr(decorated$lambda, "extra") <- TRUE
  testthat::expect_null(run(decorated, group, cells))
  decorated <- group
  attr(decorated, "extra") <- TRUE
  testthat::expect_null(run(template, decorated, cells))
  changed <- copy(cells)
  attr(changed[[1L]]$value, "note") <- list(label = "changed")
  testthat::expect_null(run(template, group, changed))
  lazy_cells <- list(list(value = 1:7))
  lazy <- template
  lazy$cells <- copy(lazy_cells)
  testthat::expect_null(run(lazy, group, lazy_cells))
})
