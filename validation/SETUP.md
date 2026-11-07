# Install cdid

The GitHub repository is public; no special account access is needed to install from it.

## Usual installation

Install `cdid` 1.2.0.1 from [CRAN](https://cran.r-project.org/package=cdid):

```r
install.packages("cdid")
```

The GitHub installation is also available:

```r
# install.packages("remotes")
remotes::install_github("dvbn/cdid")
```

`cdid` requires R 4.1.0 or later, `did` 2.5.1 or later and `DRDID` 1.4.0 or later, and refuses to load with earlier versions of `did` or `DRDID`; if older versions are installed, update them first, for example with `install.packages(c("did", "DRDID"))`. It has been tested with `did` 2.5.1 and `DRDID` 1.4.0 under R 4.5.2. Later versions are accepted, but it is not guaranteed that everything works with them or that results are identical. Source installation compiles C/C++ code, so the compilers described under [Prerequisites](#prerequisites) are needed.

This is all that is needed to use the package. Version 1.2.0.1 was published on CRAN on November 7, 2026. Mirrors and platform binaries may update later; check `packageVersion("cdid")` and, if needed, install the source with `install.packages("cdid", repos = "https://cran.r-project.org", type = "source")`.

## Optional: pinned validation setup

The setup script `validation/setup.R` installs `cdid` together with exactly the tested versions, `did` 2.5.1 and `DRDID` 1.4.0, into a dedicated library and checks the result. Use it to check this revision against its tested baseline; historical guide outputs record their original dependency versions. It is not required to install or load `cdid`.

The repository includes the unchanged CRAN source archives of these two versions in `validation/vendor/`. `remotes::install_github()` or `R CMD INSTALL` alone does not install these exact versions: they resolve `did` and `DRDID` from your libraries or the current CRAN repository, which may provide later versions. The setup script installs the archived versions, the other dependencies and `cdid` itself into the dedicated library.

### Prerequisites

- R **4.1.0** or later. R 4.5.2 is the tested version. Current CRAN dependencies may require a newer R than this minimum; the script then stops with an error instead of changing any version bound.
- Git, to clone the repository. When Git is available, the script checks that the checkout is a clean clone and records its revision.
- C and C++ compilers and `make`, because `DRDID`, `cdid` and some of their dependencies are compiled from source: Rtools on Windows, the command-line developer tools on macOS, and a C/C++ toolchain with the required system development libraries on Linux. The compilers reported by `R CMD config CC` and `R CMD config CXX`, and `make`, must be available to the R process.
- Write access to the checkout and to the target library.
- Network access to a CRAN repository for the ordinary dependencies during `--install`.

### Install with the setup script

```sh
git clone https://github.com/dvbn/cdid.git cdid
cd cdid
Rscript --vanilla validation/setup.R --project . --library validation/.lib --repository https://cloud.r-project.org --check-only
Rscript --vanilla validation/setup.R --project . --library validation/.lib --repository https://cloud.r-project.org --install
Rscript --vanilla validation/setup.R --project . --library validation/.lib --repository https://cloud.r-project.org --check-only
```

On a fresh clone the first `--check-only` normally fails with a list of missing packages; it installs nothing. `--install` installs, and the final `--check-only` confirms the installed versions. The installation needs nothing beyond the clone itself: no previously prepared library or other local files.

### Use the pinned library

Start a fresh R session and put the setup library first:

```r
.libPaths(c(normalizePath("validation/.lib"), .libPaths()))  # from the repository root
packageVersion("cdid")    # 1.2.0.1
?cdid::att_gt
```

Qualify calls as `cdid::att_gt()` so that a later `library(did)` cannot change which function runs. If you chose another `--library`, use that path instead.

### Options

- `--project PATH` is required and selects the root of the clean clone. Relative paths are resolved from the working directory.
- `--library PATH` selects the dedicated target library. The default is `<project>/validation/.lib`; another dedicated directory, such as `../cdid-library`, also works. Base, user and site libraries, the package source directories and Git metadata are refused as targets, and existing package directories linked elsewhere are not replaced.
- `--repository HTTPS_URL` selects one HTTPS CRAN repository (default `https://cloud.r-project.org`).
- `--install` installs; `--check-only` (the default) and its alias `--plan` only check. A check uses no network and installs nothing.

### What the script does

- It checks the MD5 digests and sizes of the two bundled archives before installing anything, and their SHA-256 digests once the `digest` package is available. The archives are the only source of the exact `did` and `DRDID` versions; there is no download fallback.
- Ordinary dependencies come from the `Depends`, `Imports` and `LinkingTo` fields of the two archives and of `cdid`, plus the `testthat`, `digest` and `ggplot2` suggestions. Declared version bounds are kept. Suitable versions already present in your libraries are reused; only missing or inadequate packages are installed from source, and only into the target library. No other library is modified.
- Every `--install` compiles `cdid` from a fresh copy of the checkout, so the checkout itself stays clean. Existing exact `did` and `DRDID` installations in the target are reused.
- A final check verifies the three package versions, that they load from the target library, all dependency bounds, and the package's own upstream version check.
- Each run writes a record to a new folder in `validation/setup-receipts/`: status, errors and warnings, package versions, session information, archive digests, a source file manifest, the Git revision when available and the installation logs. Git ignores these folders and `validation/.lib`.

No profiles, credentials or persistent environment settings are written.

### Tested platforms

The installer has been exercised on macOS arm64 with R 4.5.2, before the version number was set to 1.0.0.0. Other platforms have not been checked. Setup checks installation only; it runs no estimation, tests or simulations.

### Running package checks and the full test suite

Ordinary `R CMD check` executes the bounded regression suite in `tests/testthat/test-cran-regression.R` via `tests/testthat.R`. All 172 comprehensive development, numerical-verification, and Monte Carlo test files remain packaged in `tests/testthat/` and are run unfiltered (`filter = NULL`) as a mandatory CI and validation gate.

To run the complete unfiltered test suite against an installed `cdid` library (such as `check/cdid.Rcheck` after `R CMD check` or `validation/.lib` after `validation/setup.R --install`):

```sh
Rscript validation/run_full_tests.R \
  --lib-loc=check/cdid.Rcheck \
  --test-dir=check/cdid.Rcheck/tests/testthat \
  --output-dir=full-tests
```

This step is for package verification and CI; end users do not need it to use `cdid`.

