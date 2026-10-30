# Source provenance

## Package and authorship

This is version 1.1.0.0 of the existing `cdid` package maintained by David Benatia. It is
an update of the package after version 0.1.1 on CRAN (version 1.0.0.0 was an intermediate
development release not submitted to CRAN). The title and author list in DESCRIPTION are
retained: David Benatia (maintainer), Christophe Bellégo, Joel Cuerrier and Vincent
Dortet-Bernadet. The package is licensed under GPL-3.

The package source was written for this project as an independent implementation. No
source code from the GPL-2 legacy `cdid` estimator is included.

## Use of did and DRDID

`cdid` imports the public interface of `did` 2.5.1. Ordinary mode delegates estimation
to `did::att_gt()`, forwarding supplied upstream arguments lazily. An explicit
`allow_unbalanced_panel = NULL` is translated to `FALSE`; chained-only controls are
not forwarded. Chained estimation and aggregation are computed by `cdid`. Its result
adapters use the exported `did::DIDparams()`, `did::MP()` and `did::AGGTEobj()`
constructors, with `cdid` classes and reporting methods that extend upstream methods,
add diagnostics and adapt presentation metadata where needed. This keeps the familiar
`att_gt()`/`aggte()` interface. No upstream estimation algorithm is copied into
`cdid`'s R implementation.

`DRDID` 1.3.0, which supplies the Sant'Anna and Zhao (2020) estimators used by `did`, is
a dependency of `did`. The tested versions are exactly `did` 2.5.1 and `DRDID` 1.3.0, used
without local changes.

`validation/vendor/` contains the unmodified CRAN source archives of `DRDID` 1.3.0 and
`did` 2.5.1, recorded in `validation/vendor/SOURCES.csv`. They are used only by the
setup script to install those exact versions. Both archives are licensed under GPL-3,
as recorded in `SOURCES.csv`, and retain their source creators' authorship and licensing.

Callaway's development repository of `did`
(https://github.com/bcallaway11/did, commit `74b88eb07f1faa644df1271055a0f77848f6550c`)
was consulted for test and engineering ideas. It is not a runtime dependency, and no code
from it is included.

The simulation helpers use R's standard L'Ecuyer-CMRG random-number streams from the
`parallel` package. They are small synthetic designs written for this project, not a port
of the simulations in the chained DiD paper.

## Citations

Bellégo, C., Benatia, D. and Dortet-Bernadet, V. (2025). The Chained
Difference-in-Differences. *Journal of Econometrics*, 248, 105783.
https://doi.org/10.1016/j.jeconom.2024.105783

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-Differences with Multiple
Time Periods. *Journal of Econometrics*, 225(2), 200–230.
https://doi.org/10.1016/j.jeconom.2020.12.001

Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly Robust Difference-in-Differences
Estimators. *Journal of Econometrics*. https://doi.org/10.1016/j.jeconom.2020.06.003
