# Contributing to cdid

Thank you for helping improve `cdid`.

## Ask a question

Use [Issues](https://github.com/dvbn/cdid/issues) for installation help, examples and questions about assumptions or method choices. Include the package version and a small example when useful.

## Report a bug

Open an [issue](https://github.com/dvbn/cdid/issues/new/choose) with the bug report form. Describe the result you expected and what happened. Provide the smallest example that reproduces the problem, together with `sessionInfo()`.

Use synthetic data or an example you are permitted to share. Remove credentials, identifying information and local paths that reveal private information.

## Suggest a change

For a substantial feature, start a discussion so that its scope and assumptions can be considered. For a focused fix, a pull request can include the explanation directly.

Describe the problem and how the change addresses it. Include an example or test that demonstrates the issue, and state which checks you ran. Changes to an estimator, identifying assumption, default, uncertainty calculation or numerical method need an explanation of their statistical implications.

All changes require maintainer review before merging. A suggested fix or feature is not a commitment to include it in a release.

## Scoped contributor tasks

A small queue of self-contained tasks is listed in [contributing/TASKS.md](contributing/TASKS.md), separate from release work and not release blockers:

1. **Issue and claim**: Each task uses a maintainer-approved issue stating its immutable execution SHA, deliverable and checks. Claim one issue and wait for maintainer confirmation of the assignment and SHA before substantial effort. The issue's SHA must include these contributor docs and skill. Maintainers confirm scope revisions and any later release migration or source retarget; never silently follow newer `main`.
2. **Isolated checkout**: Work in a separate clone, fork or task branch at the confirmed SHA. Preserve existing dirty files. Docs-only tasks require no R installation, build or Monte Carlo runs. If R runs, follow [validation/SETUP.md](validation/SETUP.md) with a dedicated library built from that checkout; verify the source receipt and loaded package locations and versions rather than trusting the version alone.
3. **Code and statistical rules**: Material estimator, assumption, API, dependency or default changes require maintainer discussion before code. Bug fixes require a reproduction and independent check. Future performance changes must preserve estimates, full original-study-unit influence-function contributions (including zero- and single-outcome units), covariance, standard errors and target availability/status/reasons. Make no silent numerical/statistical changes, weaker tests/tolerances or seed redraws. Use the [documented matched comparisons](docs/from-did.md#why-can-balanced-results-differ-and-which-should-i-use) rather than assuming estimator equivalence.
4. **Review and evidence**: Review your diff and results before opening a PR. Report `git rev-parse HEAD`, actual commands, statuses, errors and `sessionInfo()`, leaving unrun tests explicit (prose tasks need no full test suite). Independent maintainer or designated-reviewer acceptance is required.

## Optional local coding-agent workflow

Contributors may voluntarily use Codex in their own checkout for a confirmed issue, guided by [AGENTS.md](AGENTS.md) and [.agents/skills/cdid-contribution/SKILL.md](.agents/skills/cdid-contribution/SKILL.md). Keep your own account and credentials; do not share credentials or pool accounts. This workflow accepts completed work, not transferred tokens or unused subscription allowance. The selected coding tool processes approved project files with its provider under your account; working locally does not mean inference is offline. Do not disclose private data, credentials or identifying local paths, or add coding-agent references outside contribution guidance.

Set wall-time and usage caps with paid overflow off and no paid fallback. Stop at a cap or setup/source/scope blocker with usable partial results. Record prescribed recoverable fit errors/timeouts and continue remaining planned attempts within budget, without retries or redraws. Background runs or nested delegation need a separate explicit request within maintainer-approved scope. The skill alone authorizes no public action: review and explicitly authorize any scoped posting, push or PR action. Maintainers independently accept PRs and decide merges.

```text
$cdid-contribution Issue: https://github.com/dvbn/cdid/issues/<NUMBER> | Task: <ANCHOR> in contributing/TASKS.md | SHA: <CONFIRMED_ISSUE_SHA> | Budget: <= 30 min, own usage cap, no paid overflow/fallback; preserve partials and report blockers.
```
