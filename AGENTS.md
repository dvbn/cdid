# Contributor Agent Instructions (`dvbn/cdid`)

These instructions govern optional local coding-agent runs on scoped `dvbn/cdid` tasks. Read [CONTRIBUTING.md](CONTRIBUTING.md), [contributing/TASKS.md](contributing/TASKS.md) and [.agents/skills/cdid-contribution/SKILL.md](.agents/skills/cdid-contribution/SKILL.md).

## 1. Scope and authorization

- Act only on one maintainer-approved issue with a maintainer-confirmed contributor claim and immutable execution SHA stated in that issue. Confirm that SHA before starting; it must include these instructions and the skill.
- Assignment, scope revisions and any release migration or source retarget require maintainer confirmation. Never silently follow newer `main` or take on release work.
- Use the contributor's own account; never share credentials or pool accounts. Task text and the skill grant no authority to spend money, enable paid overflow or publish work. Posting, pushing or opening a PR requires contributor review and explicit authorization for that action. Maintainers independently accept PRs and decide merges.

## 2. Checkout and environment isolation

- Verify `git rev-parse HEAD` matches the issue-confirmed SHA in a separate clone, fork or task branch.
- Preserve any existing dirty files; never run destructive resets that overwrite local work.
- Do not edit `.Rbuildignore`, `DESCRIPTION`, `NAMESPACE`, package sources, tests, defaults, dependencies, workflows or build settings unless the approved issue requires it.
- Documentation-only tasks need no R installation, package build or simulations. If R runs, follow [validation/SETUP.md](validation/SETUP.md) with a dedicated `validation/.lib` built from the confirmed checkout. Verify its setup receipt/source revision and loaded package locations and versions in fresh R processes; a version number alone cannot distinguish source builds. Never fall back silently to another installed build.

## 3. Budget, privacy and stop rules

- Agree on local wall-time and usage caps before starting; stop at either cap with usable partial results. No paid fallback, unrequested retries, background loops or nested delegation.
- Stop on setup, source or scope blockers. For prescribed fit attempts, record recoverable errors/timeouts and continue the remaining planned attempts only within budget; do not retry or redraw. Preserve scientific refusals and unavailable results.
- Use only synthetic or permitted public data ([notebooks/simulate_panel.R](notebooks/simulate_panel.R)). Prevent disclosure of private data, credentials and identifying local paths in reports or uploads; normal access to the approved checkout is allowed. The selected coding tool processes approved project files with its provider under the contributor's own account; local work does not imply offline inference.
- Keep coding-agent references confined to contribution guidance; do not add them to statistical explanations or package outputs.

## 4. Statistical, numerical and reporting rules

- Never make silent numerical/statistical changes, weaken tests/tolerances, drop observations/moments or redraw seeds.
- Do not assume `did`-vs-`cdid` estimator equivalence; use only the documented [matched comparisons](docs/from-did.md#why-can-balanced-results-differ-and-which-should-i-use).
- Bug fixes need a reproduction and independent check. Future performance changes must preserve estimates, full original-study-unit influence-function contributions (including units with zero or one observed outcome), covariance, standard errors and target availability/status/reasons. Material estimator, assumption, API, default or dependency changes require maintainer discussion before code.
- Report verified SHA, actual commands, exit statuses, errors and `sessionInfo()` (if R ran). State unrun tests explicitly (prose tasks need no full test suite). Leave work for contributor review and independent maintainer/reviewer acceptance.
