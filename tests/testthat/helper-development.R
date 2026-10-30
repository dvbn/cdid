# Load the private simulation helpers shipped in inst/development. They are
# found through the installed package (R CMD check, testthat on an installed
# library). CDID_DEV_DIR may point at a source-tree inst/development for base-R
# runs where cdid itself cannot be installed.
local({
  dev_dir <- Sys.getenv("CDID_DEV_DIR")
  if (!nzchar(dev_dir)) dev_dir <- system.file("development", package = "cdid")
  if (!nzchar(dev_dir) || !file.exists(file.path(dev_dir, "dgp.R"))) {
    stop("cdid development helpers not found; install cdid or set CDID_DEV_DIR.")
  }
  for (f in c("seed.R", "dgp.R", "mc_harness.R")) {
    sys.source(file.path(dev_dir, f), envir = parent.env(environment()))
  }
})
