# Scoped Contributor Tasks (`dvbn/cdid`)

These three report-only tasks are separate from release checks and are not release blockers. Read [../CONTRIBUTING.md](../CONTRIBUTING.md), [../AGENTS.md](../AGENTS.md) and [../.agents/skills/cdid-contribution/SKILL.md](../.agents/skills/cdid-contribution/SKILL.md) first.

| Task | Confirmed issue |
| --- | --- |
| [Offline HTML usability](#task-1-offline-html-usability) | [#1](https://github.com/dvbn/cdid/issues/1) |
| [Docs terminology consistency](#task-2-missingness-terminology-consistency) | [#2](https://github.com/dvbn/cdid/issues/2) |
| [Missing-panel profiling](#task-3-missing-panel-profiling) | [#3](https://github.com/dvbn/cdid/issues/3) |

Each issue owns its exact immutable execution SHA, including these contributor docs and skill. Claim one issue and wait for maintainer confirmation of the assignment and SHA before starting. Do not substitute a parent commit or newer `main`; any later release migration, source retarget or scope revision needs maintainer confirmation.

Do not modify package sources, tests, defaults, dependencies, workflows, build settings or release files. Keep reports in a separate local output directory, not package sources; disclose no private data, credentials or identifying local paths. All reports need the verified SHA, actions/commands, failures, unrun checks and independent maintainer/reviewer acceptance. Partial reports and failed fits are valid evidence; they do not certify statistical properties.

<a id="task-1-offline-html-usability"></a>

## Task 1: Offline usability of standalone practical guide HTML

- **Scope**: At the confirmed issue SHA, disable network access and open existing `notebooks/cdid_practical_guide.html` ([../notebooks/README.md](../notebooks/README.md)) in two available browsers if possible; record a one-browser limitation if only one is available. No R installation, fits, package build or Quarto rerender is needed.
- **Checks**: Across all 16 sections and the function map, inspect code, tables, plots, TeX math, TOC and internal navigation. Separate embedded offline resources (including failed asset requests in browser developer tools) from intentionally online links to the public guide, CRAN and DOIs.
- **Deliverable/acceptance**: Short Markdown report with OS, measured browser version(s), offline pass/fail table for those six features, failed embedded requests, intentionally online links and limitations. No source edits or rerender; no R/build/MC checks required.

<a id="task-2-missingness-terminology-consistency"></a>

## Task 2: Missingness and support terminology consistency

- **Scope**: Read only [../README.md](../README.md), [../notebooks/README.md](../notebooks/README.md), `notebooks/cdid_practical_guide.qmd`, [../docs/help/investigate-missingness.md](../docs/help/investigate-missingness.md), [../docs/help/troubleshooting.md](../docs/help/troubleshooting.md) and [../docs/methods/missing-outcomes.md](../docs/methods/missing-outcomes.md). No R installation, reruns, builds or new statistical claims.
- **Checks**:
  1. Distinguish the full study population `F`/unit table (`frame`, including zero- and single-outcome units), units with at least two observed outcomes, and units observed at both dates of a particular pair. Keep date/pair observation probabilities separate from those realized eligibility counts.
  2. Distinguish conditional date independence used by `selection_model = "independent_period_logit"` on `selection_universe = "F"` from direct pair modeling with `"pair_logit"`, and from sharing slopes across blocks/groups (`selection_pooling` or `imputation_pooling`). Pooling does not establish independence.
  3. Check the public `cdid::tidy(fit)` fields `available`, `support.status`, `support.reason` and unavailable values (`NA`, never zero). Distinguish catalog unavailability under maintained assumptions from population nonidentification; do not audit internal support/pretest objects.
- **Deliverable/acceptance**: Markdown report covering all three distinctions, inspected files and exact file/line locations, current text and proposed minimal wording only. No source edits or new statistical claims; no R/build/MC checks required.

<a id="task-3-missing-panel-profiling"></a>

## Task 3: Starter missing-panel profiling report

- **Scope**: Report-only timing of `cdid::att_gt()` on one fixed synthetic missing panel from [../notebooks/simulate_panel.R](../notebooks/simulate_panel.R). No source/optimizer changes, Monte Carlo/coverage claims or assumed `did` equivalence. No aggregation fits or extra masks.
- **Environment**: Verify `git rev-parse HEAD` against the confirmed issue SHA. Follow [../validation/SETUP.md](../validation/SETUP.md) before fitting: `validation/.lib` must contain `cdid` built from that checkout, with the documented tested dependency versions. Check the setup receipt/source revision and loaded locations/versions, not version equality alone. A stale/missing library is a setup blocker, not permission to benchmark another build.
- **Fixed data and call**: From the repository root, generate once with the unmodified defaults (`n_units = 2000L`, `n_periods = 6L`, `seed = 20261029L`). Save/reload that same `sim` for every worker; never regenerate after failure. Use the documented treatment-comparison calls (`treat-ipw-and-reg`, `treat-dr-and-table`) in `notebooks/cdid_practical_guide.qmd`:

```r
.libPaths(c(normalizePath("validation/.lib", mustWork = TRUE), .libPaths()))
source("notebooks/simulate_panel.R")
sim <- simulate_panel(default_dgp_config())  # once, before worker attempts

# In each fresh worker, select the same library and load the saved sim first.
# method is one of "ipw", "reg", "dr"; change no other option.
fit <- cdid::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x + z, selection_xformla = ~ x,
  data = sim$missing_panel, frame = sim$units,
  panel = TRUE, allow_unbalanced_panel = TRUE,
  control_group = "nevertreated", est_method = method,
  selection = "conditional_mtar", selection_est_method = "ipw",
  selection_model = "independent_period_logit",
  selection_pooling = "separate", selection_universe = "F",
  pt_assumption = "post", base_period = "varying",
  bstrap = FALSE, cband = FALSE
)
```

- **Attempt protocol**:
  1. In order `ipw`, `reg`, `dr`, run one recorded warm-up and then three timed attempts per method: `3 x (1 + 3) = 12` attempts maximum. Each attempt uses a fresh R worker. The first attempt is excluded from timing summaries; it does not warm the later fresh R sessions. Time only `att_gt()`, after library/package loading and input loading; report startup separately if measured.
  2. A supervising process/OS must enforce a hard 300-second wall timeout per worker and a total wall budget of at most 30 minutes (or the contributor's smaller cap), including setup, data generation and workers. Do not rely solely on R's cooperative time limits. Terminate a worker at the earlier per-fit or total deadline; stop at either total time or usage cap.
  3. Write attempt records incrementally. Record recoverable fit errors/timeouts and continue remaining prescribed attempts within budget; no retry, replacement fit, changed estimator/options, seed redraw or tolerance relaxation. Stop on setup/source/scope blockers. Retain completed rows, errors, timeouts and partial output when the budget ends; unstarted attempts are unrun, not successes.
  4. Only after all 12 timing attempts and if budget remains, profile at most one additional attempt with `utils::Rprof()` using the first method with a completed timed fit in the order above. Apply the same call, data and caps; record `role = "profile"` separately. If none completed, skip profiling. Maximum total: 13 attempts.
- **Deliverable/acceptance**: Return only the contributor's driver R script, raw CSV (one row per started attempt, at most 13) and a short Markdown report outside sources. CSV fields: attempt, `est_method`, role (`warmup`, `timed`, `profile`), elapsed seconds, status, error, exit status and successful-fit target counts (`n.requested = nrow(cdid::tidy(fit))`, `n.available` from its `available` field); use `NA` where no fit returned and label timeout elapsed as worker wall time. Report medians of completed timed fits only (or `NA`), counts of completed/error/timeout attempts and unrun attempts per method/role, verified SHA, actual commands, source-matched library evidence, `sessionInfo()`, CPU, OS, available memory and BLAS/LAPACK. If profiled, include a text `summaryRprof()` summary; omit binary/profile dumps. Distinguish `att_gt()` timing from process startup, and completed timing evidence from statistical validation. Retain failures and partials without claiming benchmark completion. Independent review may accept a partial report.
