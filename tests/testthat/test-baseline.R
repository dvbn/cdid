test_that("the loaded upstream versions meet the minimum and the tested versions are recorded", {
  deps <- cdid:::cdid_dependencies()
  expect_identical(deps$package, c("cdid", "did", "DRDID"))
  expect_identical(deps$minimum, c(NA, "2.5.1", "1.3.0"))
  expect_identical(deps$tested, c(NA, "2.5.1", "1.3.0"))
  expect_true(all(deps$ok))
  up <- deps$package != "cdid"
  expect_true(all(package_version(deps$version[up]) >= package_version(deps$minimum[up])))
  expect_identical(deps$matches_tested[up], deps$version[up] == deps$tested[up])
  expect_identical(deps$version[1L],
                   as.character(utils::packageVersion("cdid")))
  expect_true(cdid:::check_upstream_baseline())
})

test_that("versions are compared numerically and the minimum itself is accepted", {
  at_least <- cdid:::upstream_at_least
  expect_identical(at_least(c("2.5.0", "2.5.1", "2.5.2", "2.6", "2.10.0", "3.0.0"), "2.5.1"),
                   c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
  expect_identical(at_least(c("1.2.9", "1.3.0", "1.3.1", "1.10.0"), "1.3.0"),
                   c(FALSE, TRUE, TRUE, TRUE))
  expect_identical(at_least(c("2.4.9", NA), c("2.5.1", "2.5.1")), c(FALSE, NA))
  expect_identical(at_least("2.5.1", NA_character_), NA)
})

test_that("loading refuses only upstream versions below the minimum", {
  versions <- function(did, drdid) {
    cdid:::upstream_version_table(c("cdid", "did", "DRDID"), c("1.0.0.0", did, drdid),
                                  c("/lib/cdid", "/lib/did", "/lib/DRDID"))
  }
  # Equal to the minimum, which is also the tested version.
  equal <- versions("2.5.1", "1.3.0")
  expect_identical(equal$ok, c(TRUE, TRUE, TRUE))
  expect_identical(equal$matches_tested, c(NA, TRUE, TRUE))
  expect_true(cdid:::check_upstream_baseline(equal))
  # Later versions load and are recorded as untested.
  above <- versions("2.10.0", "1.3.1")
  expect_identical(above$ok, c(TRUE, TRUE, TRUE))
  expect_identical(above$matches_tested, c(NA, FALSE, FALSE))
  expect_true(cdid:::check_upstream_baseline(above))
  # Earlier versions are refused, naming the minimum and the loaded version.
  e <- expect_error(cdid:::check_upstream_baseline(versions("2.5.0", "1.3.0")),
                    class = "cdid_untested_upstream")
  expect_s3_class(e, "cdid_error")
  expect_match(conditionMessage(e), "requires did 2.5.1 or later", fixed = TRUE)
  expect_match(conditionMessage(e), "did 2.5.0 (loaded from /lib/did)", fixed = TRUE)
  expect_no_match(conditionMessage(e), "DRDID 1.3.0 or later", fixed = TRUE)
  e <- expect_error(cdid:::check_upstream_baseline(versions("2.6.0", "1.2.9")),
                    class = "cdid_untested_upstream")
  expect_match(conditionMessage(e), "requires DRDID 1.3.0 or later", fixed = TRUE)
  e <- expect_error(cdid:::check_upstream_baseline(versions("2.4.0", "1.2.0")),
                    class = "cdid_untested_upstream")
  expect_match(conditionMessage(e), "requires did 2.5.1 or later and DRDID 1.3.0 or later", fixed = TRUE)
})

test_that("did and DRDID come from the project library when it is declared", {
  lib <- Sys.getenv("CDID_PROJECT_LIB")
  skip_if(!nzchar(lib), "CDID_PROJECT_LIB not set (see tools/CHECKS.md)")
  lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)
  deps <- cdid:::cdid_dependencies()
  for (p in c("did", "DRDID")) {
    path <- deps$path[deps$package == p]
    expect_identical(dirname(path), lib, label = paste(p, "library"))
  }
  # The pinned project library holds exactly the tested versions.
  expect_true(all(deps$matches_tested[deps$package != "cdid"]))
})


test_that("DRDID diagnostics describe the explicitly loaded namespace", {
  expect_true(isNamespaceLoaded("DRDID"))
  deps <- cdid:::cdid_dependencies()
  row <- deps[deps$package == "DRDID", ]
  expect_identical(row$version, as.character(getNamespaceVersion("DRDID")))
  expect_identical(row$path, normalizePath(getNamespaceInfo("DRDID", "path"), winslash = "/", mustWork = FALSE))
})
