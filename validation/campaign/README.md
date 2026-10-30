# Portable Reliability Campaign Runbook

This directory (`validation/campaign/`) contains the self-contained preparation and execution toolkit for the 64-recipe, 562-logical-cell synthetic Monte Carlo reliability campaign defined in [`contract/CAMPAIGN_DESIGN.md`](contract/CAMPAIGN_DESIGN.md), [`contract/RECIPES.csv`](contract/RECIPES.csv), [`contract/CELLS.csv`](contract/CELLS.csv), and [`contract/SUBSUMED_MAP.csv`](contract/SUBSUMED_MAP.csv).

## 1. Scope and Scientific Rules

- **Synthetic data only**: All datasets are drawn synthetically from the 64 explicit data-generating laws in [`generate.R`](generate.R) using 59 L'Ecuyer-CMRG RNG streams (`set.seed(61012026)`). No external datasets or network resources are accessed during execution.
- **Descriptive 100-draw scope**: Every setting is scheduled for 100 replicates. Results provide descriptive finite-sample operating characteristics (bias, empirical SD, RMSE, mean SE, SE/SD ratio, interval width, conditional coverage, all-attempt containment, exact 95% Clopper-Pearson intervals, and paired method differences). This campaign is not a formal nominal qualification or PASS/FAIL gate.
- **Common prespecified primary causal truth**: Every estimator (`cdid` and `did` comparators) is evaluated against the same prespecified full-study causal ATT for its requested analysis measure (`ORIGINAL_COHORT_ATT` for unweighted cells, and `WEIGHTED_STATIC_ATT`, `WEIGHTED_FIRST_PERIOD_ATT`, or `WEIGHTED_BASE_PERIOD_ATT` qualified by recipe ID for the five weighted profiles). Selection, complete-case filtering, and method-specific probability limits (`RCS_COMP`, `ATT_CC`, `ATT_R`) are recorded only as separate secondary diagnostics and never redefine primary truth.
- **Prospective status**: No campaign draws or estimator fits are executed during repository preparation. `FREEZE_SOURCE_SHA` remains a placeholder until the target repository commit is supplied at freeze time.

## 2. Prerequisites and Environment Setup

Before freezing or running batches on the target machine, install the required R (`>= 4.1.0`) and Python (`>= 3.9`) environment and pinned package dependencies (`cdid == 1.1.0.0`, `did == 2.5.1`, `DRDID == 1.3.0`, `jsonlite`) with `Rscript --vanilla validation/setup.R` from the clean private checkout; see [`../SETUP.md`](../SETUP.md) for installer arguments. Version 1.1.0.0 identifies the supplied source; no public or CRAN release is implied.

The supervisor requires a POSIX host (Linux or macOS) for process-group termination and advisory locking.

The Python supervisor scripts ([`freeze.py`](freeze.py), [`run_batch.py`](run_batch.py), [`run.py`](run.py), [`summarize.py`](summarize.py)) use only the Python 3 standard library.

## 3. Pre-Freeze Static Verification and Manifest Generation

To verify all 64 recipe dispatches, 562 logical settings, 49 subsumed cell mappings, formal parameter signatures, and exact rational agreement with `validation/campaign/reference/CHECKS.json` (32 adjacent + 80 all-pairs conditional pair-change equalities) and `validation/campaign/reference/POPULATION_TARGETS.csv` (24 exact rational `RCS_COMP` targets), and to write [`manifest.csv`](manifest.csv) and [`targets.json`](targets.json) without running R fits or RNG:

```bash
python3 validation/campaign/freeze.py --check --write-manifest
```

To run the fake-worker supervisor and summary self-tests (exercising lock liveness, cross-host refusal, proven-dead interrupted attempt reconciliation with linked replacement attempts, per-operation timeouts, admission reserve, parent failure handling, and moment-regime summaries without R or RNG):

```bash
python3 validation/campaign/run_batch.py --selftest
python3 validation/campaign/summarize.py --selftest
python3 validation/campaign/check.py --scratch /path/to/check-results
```

## 4. Campaign Freeze

Once the target checkout and R library are prepared per [`../SETUP.md`](../SETUP.md) and the committed source SHA (`<COMMIT_SHA>`) is known, freeze the campaign provenance, runtime package digests, and all 29,500 L'Ecuyer-CMRG RNG stream states (`59 streams x 100 replicates x 5 purposes`):

```bash
python3 validation/campaign/freeze.py \
  --freeze \
  --repo /path/to/checkout \
  --lib /path/to/rlib \
  --out /path/to/results \
  --source-sha COMMIT_SHA
```

- Use an empty results directory. Existing frozen runs cannot be overwritten.
- Passing `--source-sha FREEZE_SOURCE_SHA` (or omitting `--source-sha`) is explicitly refused.
- At freeze time, [`worker.R`](worker.R) verifies `R >= 4.1.0`, `cdid == 1.1.0.0`, `did == 2.5.1`, and `DRDID == 1.3.0`, records the runtime `sessionInfo()`, BLAS/LAPACK paths, and installed source-function identity, dependency `DESCRIPTION`/`.rdb` hashes, the full R/dependency/BLAS/RNG profile into `<RESULTS>/r_runtime.json` and `<RESULTS>/provenance.json`, and saves the 29,500 RNG states to `<RESULTS>/seeds_29500.rds`.

## 5. Running One Bounded Batch

Each invocation of [`run_batch.py`](run_batch.py) (or its alias [`run.py`](run.py)) executes **one bounded batch and exits**:

```bash
python3 validation/campaign/run_batch.py \
  --repo /path/to/checkout \
  --manifest validation/campaign/manifest.csv \
  --out /path/to/results \
  --max-calls 25 \
  --wall-seconds 1200
```

- **Default budget**: At most `25` launched worker operations and `1200` wall-clock seconds per batch, with a `300`-second wall timeout per operation and a `30`-second cleanup reserve. A new operation is admitted only when at least `330` seconds (`300 + 30`) remain in the batch budget.
- **Process isolation**: Every operation (`FIT`, `AGG`, `PRETEST`, `DID_CC`, `DID_UNBALANCED`, `DID_RCS`) runs in a fresh `Rscript --vanilla` process. Parent `FIT` objects with `AGG` children are persisted to `<RESULTS>/parents/` so child `AGG` operations aggregate without refitting the parent.

## 6. Resuming Across Batches and Crash Recovery

To continue pending tasks, re-run the exact same [`run_batch.py`](run_batch.py) command:

```bash
python3 validation/campaign/run_batch.py \
  --repo /path/to/checkout \
  --manifest validation/campaign/manifest.csv \
  --out /path/to/results \
  --max-calls 25 \
  --wall-seconds 1200
```

- **Provenance guard on resume**: Before launching any worker, the supervisor verifies that the frozen contract, harness, manifest, and R package profile match `<RESULTS>/provenance.json`. Any mismatch refuses resumption and requires an explicit new output directory.
- **Immutability and no automatic retries**: Completed, refused, errored, and timed-out (`300s`) operations are final and are never re-run on resume.
- **Lock and interrupted-attempt recovery**: `<RESULTS>/lock.json` records `hostname`, `supervisor_pid`, `worker_pid`, `worker_pgid`, and `owner_id`. Lock age is diagnostic only and never triggers a lock steal.
  - On the **same host**, if `os.kill` proves both the recorded supervisor and worker processes are dead (`ESRCH`), the new batch marks the open attempt as `interrupted` and creates a linked replacement attempt (`attempt_number = k + 1`, `linked_previous_attempt_id`) reusing the identical cached dataset, RNG state index, and parent `.rds`.
  - Across **different hosts**, automatic reclamation is refused because remote process liveness is unknown; once the operator confirms the previous host is stopped, pass `--reconcile-cross-host <OWNER_ID>` to reconcile the interrupted attempt and continue.

## 7. Generating Tidy CSV and Markdown Summaries

During or after execution, build the tidy summary tables from the disk records in `<RESULTS>/records/`:

```bash
python3 validation/campaign/summarize.py \
  --manifest validation/campaign/manifest.csv \
  --out /path/to/results
```

This writes:
- `<RESULTS>/attempts.csv`: Per-attempt ledger including interrupted attempts, separated setup/data/operation timings, warnings, and condition classes.
- `<RESULTS>/cells.csv`: Primary causal ATT error, bias, MCSE (where moments permit), empirical SD, RMSE, mean SE, SE/SD ratio, interval width, conditional coverage, all-attempt containment (`covered / 100`), exact 95% Clopper-Pearson intervals, `FULL_FAMILY` simultaneous containment, `V_analytical` clustered SE diagnostics, and separate secondary comparator metrics.
- `<RESULTS>/aggregates.csv`: Summary of all `aggte` tasks (`AGG_S`, `AGG_6`, `AGG_NA`) against frozen population cohort masses.
- `<RESULTS>/pretests.csv`: Package Wald and conditional pretrend test rejection rates, frozen null classifications, unavailable counts and exact 95% Clopper-Pearson intervals.
- `<RESULTS>/paired.csv`: Paired differences on common completed draws (`cdid` vs `DID_CC`, `DID_UNBALANCED`, `DID_RCS`, plus `equality_ref` pairs).
- `<RESULTS>/SUMMARY.md`: Human-readable Markdown summary report.

Static verification uses only bundled finite references and fixed public formal data under `reference/`. It needs no files outside this toolkit. Freeze and run additionally require the supplied private package checkout and installed library. All 64 recipe dispatches, 562 logical settings, 100 scheduled draws, and five RNG purposes remain fixed. Complete effective comparator keys currently yield 1,494 operations per draw (1,137 cdid operations and 357 did comparisons); the manifest is authoritative.

The checks do not execute R estimators or validate installed package behavior. Before a measured smoke, confirm the reviewed source payload and receiving runtime. Verify the freeze receipt, installed source guard, one persisted-parent aggregation, and package output extraction in that bounded smoke. A self-test alone does not establish these R interfaces.

Intervals are formed from returned package estimates, SEs and critical values. Auxiliary precontrasts, including negative dynamic coordinates, are labelled separately from primary causal ATT. Incomplete runs retain the full scheduled denominator and are identified explicitly. Original unit IF rows and cluster labels determine the audit; explicit RCS is point-only. Raw requests, receipts, logs and failure records are retained, with unknown interrupted costs left undefined.
