# Pinned upstream source dependencies

The repository includes these **unmodified source archives**, not installed packages or binaries:

| Archive | Version | License | Upstream authors |
| --- | --- | --- | --- |
| `DRDID_1.3.0.tar.gz` | 1.3.0 | GPL-3 | Pedro H. C. Sant'Anna; Jun Zhao |
| `did_2.5.1.tar.gz` | 2.5.1 | GPL-3 | Brantly Callaway; Pedro H. C. Sant'Anna |

`SOURCES.csv` is readable with base R `read.csv(..., stringsAsFactors = FALSE)`. It records exact filenames, versions, licenses, byte sizes, MD5/SHA-256 digests, official CRAN current-source and version-archive URLs, and the normalized `Depends`/`Imports`/`LinkingTo` DESCRIPTION fields. The installer uses the bundled files only and does not download them.

The two archives are the original CRAN files, unextracted and not repacked. Their full original DESCRIPTIONs, authorship/copyright information, license declarations, repository/publication metadata, and other source contents remain inside the original archives. Keep these files intact; do not regenerate them from an installed library or replace their metadata. They are already required dependencies, pinned for reproducible installation, and introduce no package feature.

Their digests match `SOURCES.csv`. Keep them unchanged when moving the checkout. Missing or altered files fail setup explicitly.

Setup bootstraps verification with portable `tools::md5sum` before installation, then checks SHA-256 with the existing `digest` validation dependency after ordinary dependencies are available. Pinning these sources does not freeze the rest of CRAN: the setup receipt records the ordinary dependency versions actually used. See [setup instructions](../SETUP.md).
