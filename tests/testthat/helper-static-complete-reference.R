# Package-test loader for the independent complete-data reference used by
# test-static-complete-core.R (ref_static_complete_eval, ref_static_complete_fd_if).
# testthat evaluates helper files in the shared test environment before the
# test files; this binds the reference there, from the package-safe test path.
# Only the independent reference is read; no production source is read.
sys.source(testthat::test_path("reference", "static_complete_reference.R"), envir = environment())
