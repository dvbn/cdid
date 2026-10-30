# Private local pkgdown preview

This builder renders the existing cdid documentation with pkgdown 2.2.0,
Bootstrap 5 and Flatly. It installs nothing and does not publish or deploy the
site. Keep the repository and generated files private. The robots metadata is
an indexing preference; access remains controlled by where you keep the files.

The home page comes from `docs/README.md`. The 20 reading-path pages and
`validation/SETUP.md` become 21 articles. Release
`NEWS.md` becomes pkgdown's `news/index.html`, and the source `man/*.Rd` files
provide the function reference, including the re-exported plotting functions.
The five existing plots and all supplied logos are copied without alteration.
No new vignette is written. Markdown fences remain static; help examples are
shown without running them. Recorded versions and results remain unchanged.

| File | Role |
| --- | --- |
| `build_site.R` | Explicit input paths, staging, pkgdown build and validation. |
| `stage_site.py` | Standard-library Python adapters and source checks. |
| `check_site.py` | Standard-library Python checks of rendered links and assets. |
| `_pkgdown.yml` | Theme, navigation, article index and reference index. |

## Requirements

Use a complete, current cdid source checkout with a release `Version` and
matching ordinary version headings in `NEWS.md` (for example, `# cdid 1.0.0.0`).
Both README files must have identical content after relative links are rebased.
The input documents must be user documentation: replace authoring notes,
development headings and authoring-tool references before building. The builder
reports these problems without rewriting your source or recorded output.

Provide an existing R library containing installed cdid and did. Site tooling
may be in that library or on R's usual library path: pkgdown **2.2.0**,
rmarkdown, knitr, bslib, ragg and jsonlite, plus pandoc and Python **3.9 or
later** on `PATH`. These are optional website tools; the builder makes no
changes to the package's dependencies or installed libraries.

pkgdown 2.2.0 caches third-party build assets and has no supported offline
switch. An offline build requires those assets, including Flatly fonts, to
already be cached. `pkgdown.internet = FALSE` does not disable networking in
this version. Use the build environment's existing network restrictions when
building offline; missing cache entries are a build blocker. This builder
suppresses automatic remote favicon generation and CRAN release-date lookups.
It does not install tools, populate caches or generate favicons.

## Build

From the checkout root, with these files in `site/`:

```sh
Rscript --vanilla site/build_site.R \
  --source . \
  --library validation/.lib \
  --logo path/to/logo.png \
  --output ../cdid-site-preview
```

All four options are required. Quote paths containing spaces. For a separate
existing tooling library, set `R_LIBS` for this command. Nothing is installed.

`--logo` selects the PNG. If a file with the same stem and `.svg` extension is
present beside it, that SVG is also checked and copied. SVGs must be static,
self-contained images. Home-page HTML image links are rebased to the copied
logo. A PNG alone is sufficient.

`--output` must name a new directory whose parent exists, outside both
`--source` and `--library`. Existing outputs, including dangling symbolic
links, are refused. The source checkout is never an output destination.
Staging uses a new directory inside `.site-build/` beside the output. That
folder must be an ordinary directory outside the source and library. Staging
is removed on success and retained on failure for inspection. Choose another
new output directory for a subsequent build. Keep staging and generated files
outside the repository you export.

## Validation and browser review

Before rendering, staging checks the full page mapping, source link graph,
linked headings, navigation and help-topic coverage. Executable R Markdown
chunks, inline R and Rd `Sexpr` expressions are refused. Only the selected
source pages, help files, package metadata, figures and logos are staged;
existing source build hooks and unrelated files are not copied.

After rendering, the checker verifies expected pages, all discovered local
HTML and CSS references, anchors, images, fonts, stylesheets, scripts, logo
hashes, plot hashes and search-index targets and page coverage. It rejects
remote runtime assets, links to Markdown sources, generated adapter source
URLs, machine-specific paths, development headings and authoring-tool
references. External documentation and private repository hyperlinks remain
links; no external requests are made by the checker.

A successful static check does not establish browser behavior. Open the home,
article and reference pages; check Flatly layout, tables, mathematics, plots,
logos, dropdowns and narrow-screen navigation. Search uses pkgdown's local
`search.json` and requires HTTP because browsers restrict `fetch()` on
`file://` pages. When authorized, use a local server and check a search result
from a nested article, including its target heading. The checker does not
start a server or browser. A separate browser review is still required.
