# Independent packaged pooled-OR reference

These machine exports come only from the frozen independent fixed-input
reference. No production backend supplied an expected value. `export_reference.py`
reads the frozen 108-unit CSV, calls only its already-frozen `finite()` and
`analytic()` equations, and verifies exact equality to all frozen point vectors
and complete rows. It disables bytecode writes and never calls the frozen
table-writing or verification entry points. Original independent files remain
unchanged. No RNG, DGP draw, Monte Carlo, fit of the old backend, package build,
namespace injection, dependency installation or compilation occurred.

The original mathematical cohort indices 0,1,2,3 are relabeled 10,20,30,Inf.
This is a bijective fixture-key mapping; it changes no observations, weights,
parameters or influence rows. On calendar 2,5,17,29 the first NYT pool is
20/30/Inf for pair (2,5), and later pools are 30/Inf for pairs (5,29) and
(2,29). All controls are clean at their actual pair endpoints. The fixed
calendar has three modeled one-period transitions, giving static-X multipliers
1,2,3 for these pairs, not numeric time gaps 3,24,27.

`raw.csv`, `calendar.csv`, `cohorts.csv`, `pairs.csv`, `cases.csv`, `blocks.csv`
and `links.csv` define the inputs and catalog. Each of 12 cases has three unique
scientific control blocks. The fourth NYT link is a distinct reporting consumer
of the first block; its link ID and provenance differ, but its block objective
is not repeated. The 12 cases retain exactly the frozen NT/NYT baseline/F/B/A,
NYT baseline/F no-intercept, and NYT baseline/F zero-shared-slope branches.

`selection_F.csv`, `selection_B.csv`, `selection_A.csv` export exact three-cell
logit design Q, U membership and fixed known rho. Q and rho are NA off U;
`transport.csv` retains the clean fixed-law cell values. The loader sets outcome
X to NA off A while retaining required F selection predictors. No NA value is
used through multiplication by zero. Transport is fixed under case derivatives
and is never estimated or sample-normalized.

`points.csv` and `coordinates.csv` specify every coordinate, while each
`*_rows.csv` holds the complete original-unit rows and `*_link_rows.csv` holds
the link columns directly. Gamma order is local block intercepts, when present,
then shared x and x-squared slopes. A no-intercept case has only shared slopes;
a zero-slope case has only local intercepts. Coordinate upstream names retain
the frozen zero-based pair/cohort notation, mapped explicitly by the catalogs.
The complete reference stack also retains unused pair/cohort nuisance
coordinates of its fixed table. A backend need only compare nuisance fits used
by its actual control/target catalog; these extra reference columns do not
require fitting irrelevant models. The `objective` summary is full RSS/N with
N=108; the unnormalized kernel RSS is 108 times this summary. All blocks retain
their original relative weights.

`load_fixture.R` defines `load_pooled_or_reference(directory, case, order)`. It
returns real panel, frame, reader and link views plus raw inputs, recipe/catalog
metadata and independent expectations. It calls only real locked installed
namespace constructors, using the usage patterns in the existing baseline and
Corollary 1 owner tests. It deliberately defines no future backend API and
calls no fit or replay method. The loader keeps one explicit all-NA panel ID;
11 other zero-observer IDs enter only through the linked frame. There are 72
repeats and 24 singletons among the 108 frame identities.

One saved construction pass on R 4.5.2 used an installed cdid namespace and
its dependency libraries. All 12 cases passed
with a fixed frame-order permutation, required NA patterns, full row alignment,
share/count checks and duplicate-consumer catalog checks. See `loader_stdout.txt`
and `loader_construction.csv`. This is owner/fixture readiness, not pooled
backend, scientific-assumption, inference-coverage or package acceptance.

The packaged tests use the loader with an installed backend and compare its
outputs with these independent expectations. Do not edit
fixtures manually; the saved exporter owns every scientific CSV.
