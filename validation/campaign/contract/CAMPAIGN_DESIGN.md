# Portable reliability campaign: design contract

This is a prospective design: no fits, draws or results exist, and no table reports any. The campaign is not a work-stop or qualification gate and promises no release, version or platform. `FREEZE_SOURCE_SHA` remains a placeholder until the private source commit is supplied.

`RECIPES.csv` gives one row per data law: law, common causal `truth`, secondary comparator references and RNG stream. `CELLS.csv` holds the targeted cells with literal arguments. Core-sweep cells come from rule 6.1.

## 1. Literal calls

FIT (later entries override earlier ones):

```
cdid::att_gt(yname="y", tname="t", idname="id", gname="G", xformla=<recipe xformla_default>,
  data=PANEL, panel=TRUE, allow_unbalanced_panel=TRUE, control_group="nevertreated",
  anticipation=0, alp=0.05, bstrap=FALSE, cband=FALSE, est_method="dr",
  base_period="varying", compute_inffunc=TRUE, pl=FALSE, cores=1, mode="chained",
  frame=F, selection="none", gmm_weighting="optimal", pt_assumption="post",
  <recipe sweep_fixed_args>, <cell args>)
```

Omitted arguments keep package defaults; `UPSTREAM_ONLY` drops `frame`, `allow_unbalanced_panel` and chained-only arguments.

- `J_REP` = "Within each cohort, membership in the fitting set is independent of the scientific conditioning variables by design." It is true for REPRESENTATIVE_ROSTERS and E_TWO_DATES_A and false for DATE_BASE.
- `J_TRANSPORT` = "Column holds the known normalized ratio f(X|G,U)/f(X|G) computed from the design law with population cohort means; it is not estimated."
- The callback literal is `DRDID::std_ipw_did_panel`: `function(y1, y0, D, covariates, i.weights=NULL, boot=FALSE, boot.type="weighted", nboot=NULL, inffunc=FALSE, trim.level=0.995)`, returning `ATT` and `att.inf.func` when requested. It matches the supplied panel callback calls; the RCS callback instead takes `y, post, D`, so this panel callback is confined to its existing cells. Equality with built-in IPW is not asserted.

PRETEST is `cdid::conditional_did_pretest()` with the cell arguments, minus `pt_assumption`, `pt_start`, `gmm_weighting`, `gmm_covariance`, `compute_inffunc` and `allow_unbalanced_panel`, plus `biters=999`. Children aggregate the FIT object and inherit its inference:

- `AGG_S`: `aggte(type="simple")`.
- `AGG_6`: simple; group; dynamic; dynamic with `balance_e=1`; dynamic with `min_e=0, max_e=2`; calendar.
- `AGG_NA`: simple with `na.rm=FALSE` and with `na.rm=TRUE`.

PANEL has all N×T rows. Unobserved outcomes are NA; dated covariates and weights stay present. F has one row per generated unit (id, G, static covariates and weights, cluster, known transport columns). Dated attributes have no static copy in F. α, U and innovations never enter PANEL or F.

## 2. Population laws

Defaults: N=2000, T=5, dates 1..T, p=5.

- **Covariates.** X1 ~ Bernoulli(1/2), X2 ~ N(0,1), X3 and X5 centered Bernoulli (±.5), X4 ~ Exp(1)−1. For j≥6, (j−6) mod 3 = 0, 1, 2 gives normal, centered binary or Exp−1. All are independent and fixed within unit. X5 is selection-only unless a recipe overrides it. z is everything except X5.
- **Cohorts.** P(G=0,3,4,5 | X1=0) = (.4,.3,.2,.1) and P(·|X1=1) = (.2,.2,.3,.3). Masses are π = (.30,.25,.25,.20); μ_g = P(X1=1|g) = (1/3,.4,.6,.75).
- **Untreated outcomes.** Y0 = α + .1k + k·b'z + ε, with k the stored calendar position, α = ±1 equiprobable, b1 = .3 and b_j = .2/√p otherwise.
- **Effects.** For g>0 and t≥g the effect is θ_g + .25(t−g) + .5X1, with θ_g = 1 + .1(g−3). It is zero otherwise.
- **Errors.** ε = .5hu with h ∈ {1, √(1+.5X1)}.
  - u is N(0,1) in main recipes, or Exp(1)−1 or √(3/5)t5 (finite fourth moments).
  - √(1/3)t3 is a finite-variance sensitivity case; t1.5 has infinite variance and is a failure experiment.
  - AR1: ε_1 = .5hu_1 and ε_t = .7ε_{t−1} + .5h√.51·u_t (stationary).
  - Clusters: ε = .5h(√.2·v_ct + √.8·v_it), with cluster = ((id−1) mod C)+1.
  - Never truncate, clip or redraw.
- **Draw order** per stream and replicate. The population stream draws covariates by coordinate, then U_G, then α, then innovations (N×T, column-major), then cluster shocks. The observation stream draws uniforms in the law's listed order. Discrete draws invert uniforms against cumulative probabilities, with edges in lexicographic order. AR1_MCAR_DATE reuses C_MCAR_DATE's stream, so the two are paired.

## 3. Observation laws

- **COMP.** Each unit is observed at one window e=(s,t) in E (adjacent edges or all pairs).
  - P(e|g,α) = (1+λ_gαh_e)/K, with h_e = (T+1−s−t)/(T−2) and λ = (0,.6,.4,.2). Σh_e = 0 and |λh| ≤ .6, so the probabilities lie in [.4/K, 1.6/K], sum to 1, and P(e|g) = 1/K.
  - α cancels from changes and windows are independent of (X,ε) given G, so `selection="none"` is valid even though levels shift.
  - The observed α mean is m_gd = λ_gH_d, where H_d averages h_e over the edges containing d. For adjacent edges H_d is (1, (T+1−2d)/(T−2) in the interior, −1); for all pairs it is ((T+1)/2−d)/(T−1).
  - P(R_d|G) is the same for every cohort.
- **CC.** A unit has a complete history with probability c_x = (.2,.6); otherwise one pair, uniform or drawn from COMP. q_e = c_x + (1−c_x)/K is saturated in X1.
- **DATES.** Dates are independent with the listed r. Zero- and single-outcome units stay in F, and masks are never conditioned on repeats.
- **ATTRITION.** R_1 = 1, and for t≥2, R_t = 1{U ≤ expit(a_t+.6X1+γ_gα)} with a = (1.5,1,.5,0) and γ = λ.
  - The thresholds decrease, so histories are monotone and S_st = R_t.
  - P(S_st|X,G) is the α-average ρ_t(X1,G), which a per-cell pair logit in X1 represents exactly on F.
  - MTAR holds given X (U and α are independent of ε, and α cancels). `none` is invalid because X1 drives both dropout and changes.
  - The date product is invalid (dependent dates) and is expected to be refused (the first date is always observed).

## 4. Truths and targets

**Common primary causal truth.** Every algorithm is evaluated against the same prespecified full-study causal ATT for its requested analysis measure. Without weights, this is E[Y_t(g)−Y_t(0)|G=g] under the full unweighted population: default θ_g+.25(t−g)+.5μ_g after treatment, zero before unless explicitly overridden. With weights, use the declared weighted full-study measure, not an unweighted substitute. `RECIPES.truth` gives the recipe law and `CELLS.primary_target_override`, when nonempty, specifies the cell’s fixed weighting reference. Observation selection, complete-case membership and method-specific probability limits never redefine primary truth. All truths are frozen before draws; no empirical shares or fitted plug-ins supply them.

**Common primary aggregates.** Apply the requested selectors to the prespecified full-study causal ATT and its declared population masses. Unweighted cells use π_g. W_STATIC uses masses proportional to π_g(1+.5μ_g); DATED_WEIGHTS uses first-calendar masses proportional to π_g(1+.1μ_g), even when its cohort ATT uses baseline weights. Normalize over the groups appropriate to the frozen requested summary. Selected/fitted sample shares never replace these masses. Freeze the full coordinate set and denominators before draws. Partial-support and na.rm=TRUE outputs retain the original full-target reference and report their target/support change separately; never recompute truth only over surviving cells. Full-family containment requires every frozen coordinate.

**Literal did conventions.** The supplied source defaults to `faster_mode=TRUE`; the fast and ordinary branches agree on these rules. Let dates be sorted labels d_1,...,d_T, anticipation a, and b_g = max{d_j: d_j+a<g}. Universal base uses b_g throughout, with a normalized zero at t=b_g. Varying base uses the preceding stored date when t<g and b_g when t≥g; its first-date cell is absent. Anticipation arithmetic uses numeric labels. A cell (g,t) with actual base b has pool C = {0} for NT; NYT is {0} ∪ {k≠g: k > max(t,b)+a}. These are upstream rules, not a change to cdid.

For `aggte`, let q map the sorted union of date and cohort labels to consecutive positions:

- Simple: cells q(g)≤q(t)≤q(g)+max_e, weighted by π_g per cell. Group: the same within-group window, equal date weights within each group, then π_g weights across groups. `min_e` and `balance_e` do not select these types.
- Dynamic: e=t−g uses original numeric labels, including negative event times. `balance_e=B` retains groups with d_T−g≥B and restricts e to B−d_T+d_1≤e≤B; then apply `min_e≤e≤max_e`. Within e use π_g weights; overall average is equal-weight over selected e≥0.
- Calendar: at each date use g≤t with π_g weights, then equal-weight the eligible dates. Event-window/balance arguments are ignored. Simple/group/calendar start at g, without moving the threshold to g−a; anticipation cells may appear in the dynamic coordinates but not its nonnegative overall mean.

For this campaign the full date list determines d_1 and d_T. Upstream can recompute selectors after NA removal; label that support change without changing the primary reference. The source's fitted aggregation masses use first-period unit weights for balanced panels, mean observed-row weights per retained ID for unbalanced panels, and observation weights for explicit RCS; those describe the estimator, not the common population reference.

**Secondary comparator references.** `did_cc_secondary_reference`, `did_unbalanced_secondary_reference` and cell `secondary_target` are diagnostics only. With π_k|x,C = P(k|x)/Σ_{j∈C}P(j|x), the IPW, OR and DR limits of did's RCS construction coincide under the accepted COMPOSITION conditions:

RCS_COMP = ATT(g,t) + (m_gt − m_gb) − Σ_x P(x|g) Σ_{k∈C} π_k|x,C (m_kt − m_kb),

with ATT = 0 in these pre cells. They coincide because observation is independent of X given G with cohort-invariant date rates, and the α shift given X is linear in binary X1. ATT_CC, ATT_R and CC_COMPOSITION's term m_kd(x) = (1−c_x)(2/T)λ_kH_d/r_x remain secondary. DATE_HET, ATTRITION_MTAR and DATED_WEIGHTS have unknown method-specific secondary limits here; that does not block their primary bias or causal-truth interval coverage. A did refusal is recorded as unavailable-comparator. Literal comparator calls are unchanged.

**Field contract.** Every manifest task records primary_target_id, primary_truth, primary_population_masses, primary_coordinate_set, target_match, causal_validity and inference_eligible separately from secondary_reference/status/error. Unweighted cells use ORIGINAL_COHORT_ATT. The five weighting profiles use WEIGHTED_STATIC_ATT, WEIGHTED_FIRST_PERIOD_ATT or WEIGHTED_BASE_PERIOD_ATT, qualified by recipe ID; their primary truths come from primary_target_override, with population masses as specified above. All algorithms linked to one requested comparison use that same prespecified primary reference; a deliberately mismatched procedure remains evaluable against it and is labelled mismatched. Unweighted references for a weighted study are secondary diagnostics only. Unknown secondary limits are blank, never zero; they do not block primary bias or causal-truth interval containment. All primary errors are estimate minus primary_truth, including unexpected successes and invalid-model cases.

## 5. Model classification

The models are p (propensity), m (shared counterfactual change regression), q (observation probability) and j (actual-change imputation). Sufficient branches:

- `none`: representative pair-change means.
- Observation IPW: q.
- Imputation: j.
- Augmented correction: q or j in every needed cell.

In addition, treatment IPW needs p, reg needs m, and DR needs p or m. Every active model must be regular even when another arm is correct. A violation removes the guarantee, not necessarily the accuracy; errors, pseudo-targets and cancellation are all possible.

- **Binary X1** saturates every NT and NYT pool's log odds, and the other covariates have zero coefficients. Rate rescaling changes intercepts for single-cohort pools and coefficients for multi-cohort pools, and the logit stays saturated. Continuous-confounder recipes use cohorts 0 and 3 only, so both pools are {0}.
- **m** has one slope on L_e·z, with L_e counted in stored calendar steps, plus local block intercepts.
- **j.** Separate j slopes absorb crossing-pair effect shifts. Pooled j is correct only when κ=0.
- **Selection-only regressors.** A regressor used only in selection that does not move changes leaves j correct. One that moves changes breaks j, while m stays correct after a valid reconstruction because X5 ⟂ (X,G).
- **Transports.** For r(X1), u_B = 1−(1−r)^T, u_A = u_B − Tr(1−r)^(T−1), ρ_U = u_U(x)/E[u_U|g] and q_U = r²/u_U, using population μ_g; never estimated or normalized. Date products are F-only.
- **Distinct failure types.** Innovation-dependent windows break mean-MTAR despite a correct marginal q, distinct from wrong date independence, wrong nuisance models and invalid PT.

## 6. Matrix

**6.1 Core sweep.** Each `core_sweep=core48` recipe gets the full product, in this order:

- NT, NYT
- PT-post, or PT-all with `pt_start=1`
- varying, universal
- identity, or feasible two-step (literal `"optimal"`; the `"two_step"` alias is checked separately)
- ipw, reg, dr

That is 48 cells, with ID `<recipe>__<NT|NYT>_<post|all>_<var|uni>_<id|opt>_<est>` and child `AGG_S`. PT-all from start 1 adds restrictions for every treated cohort, and they are true in all nine recipes: COMP_adjacent, COMP_arbitrary, CC_BASE, CC_COMPOSITION, BAL_BASE, C_MCAR_DATE, C_MCAR_PAIR, AR1_MCAR_DATE and ATTRITION_MTAR (432 cells).

Sweeps run only at the recipe's own N, T, p and errors. Ladders, stress, observation profiles and refusals are not crossed. A targeted cell whose complete arguments equal a sweep cell on the same data is not created; 49 are subsumed, mapped in SUBSUMED_MAP. Equal post estimates across `base_period` are not separate evidence.

**6.2 Targeted cells.** 130 profiles cover observation, pooling, transport, robustness, PT, anticipation, weights, calendar, inference, aggregation, pretest, alias, mode, formula and refusal behaviour. An `equality_ref` row records its maximum absolute difference from the reference.

**6.3 Comparators.** Per recipe and replicate, deduplicate the complete effective supported upstream argument dictionaries over non-refusal, non-callback cells. Keys include control_group, est_method, base_period, xformla, anticipation, weightsname, fix_weights, alp, bstrap, cband, biters, clustervars, compute_inffunc, faster_mode, pl and cores when active. Preserve matched weighting/inference for matched runtime comparisons; label ignored upstream fields. Chained-only settings are inactive upstream. Each key runs:

- did CC (`allow_unbalanced_panel=FALSE`), unless the law has no complete units;
- did unbalanced (`TRUE`);
- explicit RCS (`panel=FALSE`, observed rows, no idname), point-only, for dr with the varying base.

Comparators inherit the effective matched inference flags, except explicit point-only RCS. Each physical comparator runs once per draw; multiple linked comparison contexts retain their own declared primary reference without being counted as extra replications.

The unbalanced-panel call keeps `panel=TRUE, idname="id"` at its public entry. On genuinely unbalanced observed data, both upstream branches switch internally to RCS estimation, retain original IDs, scale each row IF by M/n_cell (M = retained unique IDs, n_cell = rows in the two-period risk pool), then SUM row contributions per ID. Analytic SE is sqrt(sum(IF_id²))/M; align IF rows by their ID labels, never by row position across branches. Do not apply an extra half-fold or treat the rows as independent units. A balanced input stays on the panel path. Explicit `panel=FALSE` without idname creates one internal ID per row and cannot preserve repeated-unit dependence; it remains POINT-only (`inference_eligible=FALSE`), with no SE or coverage reported. None of these paths changes primary causal truth.

**6.4 Size.** There are 562 logical settings and 56,200 cell-draws. Generate and independently reconcile the actual deduplicated operation/parent-child inventory before freeze. Earlier forecasts of 113,700 chained/pretest/aggregation calls and 36,300 did calls are historical; the narrower old comparator-key expansion gave 34,200 did calls. The effective-key correction above may change physical counts. Do not add duplicate calls, remove logical cells or alter any science to force one forecast. The manifest inventory is authoritative once reviewed.


## 7. Execution and resumption

The freeze writes `manifest.csv` (task IDs, literal calls, the section 4 fields, classifications, parent links), the RNG states and provenance. It runs no fits. A real source commit must replace `FREEZE_SOURCE_SHA` before freeze.

`python3 campaign/run.py --repo CHECKOUT --manifest campaign/manifest.csv --out RESULTS [--max-calls N] [--wall-seconds S]`

- **Batches.** Each invocation runs one small restartable batch: default maximum 25 calls, total wall cap 1200 seconds, then exits. The user may explicitly choose another call maximum or wall budget. No automatic continuation or daemon. Count actual launched operations, including failures/refusals and aggregation/bootstrap/pretest operations, against the call maximum. Resolved not-applicable/parent-unavailable records need no worker.
- **Budget.** Use a monotonic deadline covering lock/setup, workers and cleanup. Reserve 30 seconds for termination, atomic records and lock release. Start an operation only with at least 300+30 seconds remaining; otherwise leave it pending and exit. Each operation has a 300-second wall timeout including bootstrap. Startup/preparation consumes the total budget. If it exhausts the admission reserve, launch no operation. Stop the worker process group at its deadline, reap it and flush records within the reserved cleanup time; no background worker may survive a normal exit.
- **Lock.** An exclusive lock records hostname, supervisor PID, active worker PID/process group when present, and attempt/owner identity. Age is diagnostic only. On the same host, check actual process liveness for the recorded supervisor and any worker; alive or unknown (including denied checks) refuses reclamation. Only proven-dead owners/workers permit atomic reconciliation: preserve the open attempt as interrupted and link its replacement. Another host's liveness is unknown and requires explicit operator reconciliation; never silently steal the lock or duplicate work. Release only the lock owned by this invocation.
- **Order.** Tasks run in frozen order: replicate, stream, recipe, cell, call. Children follow their parent, and the parent object persists until they settle.
- **Calls.** Use the existing validated Python-supervisor/R-worker pattern: Python 3 standard-library `subprocess` and atomic file operations supervise one `Rscript` worker per operation. Load frozen datasets, states and persisted parents; save the result before settling the task. Reuse the project's worker/result conventions during implementation. This needs no `callr`, new R dependency or custom R dependency framework. Setup time and operation time are recorded separately.
- **Dispositions.** success, partial-support, expected-refusal (documented `expected_condition` class, otherwise any `cdid_error`), unexpected-success, unexpected-error, timeout, interrupted, parent-unavailable, unavailable-comparator, not-applicable. Each is written atomically per task, and the index is rebuilt at exit. Records keep condition classes and messages, warnings, primary truth and secondary references separately, counts, timings and small diagnostics. Any usable point/interval returned even with a warning, partial-support or unexpected-success disposition contributes to its eligible primary summaries.
- **Resumption.** A new explicit invocation continues pending IDs. Final dispositions are immutable; errors and timeouts are never retried. An interrupted attempt is retained with a new linked attempt ID under the same task ID. Its replacement reuses identical persisted datasets, parent artifacts and RNG states, so nothing is redrawn, trimmed, clipped or substituted. Refuse changed frozen inputs or runtime profile in the same campaign (section 9).
- **RNG.** `RNGkind("L'Ecuyer-CMRG","Inversion","Rejection")` with `set.seed(61012026)`. The state for stream s, replicate r and purpose p (population, observation, FIT bootstrap, aggregation, pretest) is the n-th successive `nextRNGStream`, with n = ((s−1)·100 + r−1)·5 + p. Every call loads its own state, and stream indices are append-only.
- **Workers.** Package `cores=1` and `pl=FALSE`. This is scheduling, not a performance experiment.

## 8. Summaries

For every returned algorithm, per primary causal coordinate and primary population aggregate, report regardless of assumption validity or whether a secondary limit is known:

- availability;
- `primary_error = estimate − primary_truth`; bias = mean(primary_error) over available points, with MCSE = SD(primary_error)/√n when its moments permit;
- empirical SD, RMSE, mean SE, the SE/SD ratio and mean width.

Primary interval coverage/containment tests whether the interval contains `primary_truth`, including for RCS_COMP, ATT_CC, ATT_R, invalid PT, target mismatch and unknown secondary limits. Label causal validity and target match alongside the metrics; a biased algorithm must not get primary bias zero by centering on its own limit. Conditional coverage uses available eligible intervals. All-attempt containment is covered/100, with missing eligible intervals counted as not covered, and is not called nominal coverage. Explicit point-only workflows have coverage/SE not applicable and retain point metrics; never fabricate intervals. Expected-refusal workflows retain disposition counts; any unexpected returned output is evaluated against the common truth. Rejection rates and simultaneous containment over the frozen full set follow the same rules. Every defined proportion gets an exact 95% `stats::binom.test` interval; there are no Wilson or bootstrap Monte Carlo intervals. Paired differences use common available draws and state that denominator. Report point, interval and refusal denominators explicitly. Secondary-reference error/containment, if reported, must occupy separate labelled fields/tables and cannot suppress primary metrics.

Moment labels: the main laws and t5 have finite fourth moments. t3 has finite variance only, so bias MCSE stands but RMSE and SD carry no uncertainty claim. For t1.5 retain cautious finite-sample descriptive error, SD/RMSE and interval-containment counts, flag infinite variance, and set moment-based MCSE undefined; do not imply finite population variance or reliable moment estimation. Binomial containment intervals remain defined. Hundred-draw cells are descriptive: there is no nominal qualification and no PASS, 2000- or 8000-draw gate. This descriptive interpretation never suppresses the requested primary causal comparisons.

For runtime, report separately:

- the latency of successful calls (mean and median);
- the all-attempt computation cost, including successful operations, refusals, errors, interruptions and timeouts, with supervisor/setup cost separately identified;
- the number of timeouts, right-censored at 300 seconds.

Fast failures never count as estimator speed.

## 9. Provenance and installation

Install the source and dependencies with `validation/setup.R` before freezing. Freeze `FREEZE_SOURCE_SHA` from the supplied source commit, exactly did 2.5.1 and DRDID 1.3.0 (the existing package guard enforces them), recipe/cell/map and harness hashes, and all expanded literal calls. R ≥ 4.1 is the pre-freeze eligibility rule. At the first actual freeze, record the real R version/build, platform, `sessionInfo()`, BLAS/LAPACK, loaded dependency versions and source identities, plus Python and supervisor versions. Build timestamps and installation paths do not determine source identity. Before any resumed worker starts, require the identical frozen source, recipes, harness and recorded R/dependency profile. A changed profile requires an explicit separate new run, never silent continuation or overwrite. The source commit is unresolved until supplied; no placeholder may pass freeze. Python 3 may be required for remote validation tools; this changes no R package API or dependencies. Bitwise equality across platforms is not promised. Harness tools stay outside package dependencies.
