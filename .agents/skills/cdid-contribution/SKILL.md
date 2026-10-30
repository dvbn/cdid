---
name: cdid-contribution
description: Carry out a scoped dvbn/cdid contributor issue, documentation review or profiling report with a maintainer-confirmed claim, issue-pinned execution SHA and verifiable evidence. Not for general cdid package use.
---

# `cdid-contribution` Skill

Use this skill for one confirmed `dvbn/cdid` contributor issue from [../../../contributing/TASKS.md](../../../contributing/TASKS.md), whether matched by its description or invoked with `$cdid-contribution`. Follow [../../../AGENTS.md](../../../AGENTS.md) and [../../../CONTRIBUTING.md](../../../CONTRIBUTING.md).

## Guardrails

- Use the contributor's own account without shared credentials or pooled allowances. The workflow accepts completed work, not transferred tokens. Agree on wall-time/usage caps; no paid fallback, unrequested retries, background runs or nested delegation.
- The skill alone grants no public posting, push, PR, assignment or merge authority. The contributor must review and explicitly authorize any scoped posting/push/PR action; maintainers confirm assignments and scope revisions and independently accept PRs.
- Do not edit `.Rbuildignore`, `DESCRIPTION`, package sources, tests, defaults or release files unless specified by the approved issue. Keep coding-agent references confined to contribution guidance; prevent disclosure of private data, credentials or identifying local paths. Approved project files are processed by the selected tool's provider under the contributor's account.

## Workflow

1. **Verify claim and inputs**: Require issue URL/number, task anchor in [../../../contributing/TASKS.md](../../../contributing/TASKS.md), the immutable execution SHA stated in that issue, wall-time/usage budget, and maintainer confirmation of the contributor's claim. The SHA must include the contributor docs/skill; later source or scope revisions require maintainer confirmation.
2. **Verify isolated checkout**: Check `git rev-parse HEAD` matches that SHA on a separate clone/branch (never retarget silently to newer `main`). Check `git status --porcelain` and preserve existing dirty files. Docs-only tasks need no R/build/MC. If R runs, use [../../../validation/SETUP.md](../../../validation/SETUP.md) with `validation/.lib` built from this checkout; verify the receipt/source revision and actual loaded locations and versions, not just version equality.
3. **Execute bounded task**: Follow only the confirmed task under its caps. Stop on exhausted time/usage or setup/source/scope blockers, retaining partials. Record prescribed fit failures/timeouts and continue remaining planned attempts only within budget, without retries, redraws or weakened tolerances. Preserve the estimates, full original-study-unit uncertainty contributions, covariance, SE and availability invariants in [../../../AGENTS.md](../../../AGENTS.md).
4. **Handoff for review**: Prepare deliverables and evidence (verified SHA, commands, exit statuses, errors, `sessionInfo()` and platform/BLAS if R ran, and unrun checks) for contributor review and independent maintainer/reviewer acceptance. A partial report with failed fits is valid evidence, not certification of statistical properties.
