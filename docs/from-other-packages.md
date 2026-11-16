# Estimate familiar DiD models with cdid

The same `cdid::att_gt()` interface can reproduce several widely used difference-in-differences estimators and extend them to panels with incomplete outcomes. Choose the comparison groups and treatment history first. Then decide how to adjust for treatment assignment and whether to correct for selective observation. This page describes the development version. The [CRAN release](https://cran.r-project.org/package=cdid) has its own help and release notes.

## Matched specifications in cdid

For a numerical match, use the same sample, outcome and covariates, untreated reference, and weights on the reported effects. Match the variance convention separately. Here `g` is the first treatment date for an absorbing treatment and `d` is the current binary treatment status for a treatment that can switch on and off.

| To reproduce | Use in `cdid::att_gt()` | Where the match holds |
| --- | --- | --- |
| Callaway–Sant'Anna, R [did](https://bcallaway11.github.io/did/reference/att_gt.html) | `mode = "ordinary"` | Direct delegation, including repeated cross sections. [Matched chained cases](from-did.md#when-chained-and-did-estimates-should-agree) explain when the within-unit version gives the same estimates. |
| Callaway–Sant'Anna, Stata [csdid](https://psantanna.com/csdid/articles/comparing-did-commands.html) | `mode = "ordinary"` for the corresponding R calculation, or a chained fit on the [matched cases](from-did.md#when-chained-and-did-estimates-should-agree) | Align the two-date samples, controls, nuisance fits, and aggregation. In an incomplete panel, `bal(pair)` keeps the units with both outcomes in each comparison. Chaining can also use other supported pairs to estimate an effect. |
| Wooldridge's linear cohort-by-date regression, R [etwfe](https://grantmcdermott.com/etwfe/reference/etwfe.html) and Stata [jwdid](https://github.com/friosavila/stpackages/blob/main/jwdid/readme.md) | `est_method = "reg", control = list(preset = "cohort_time_fe")` | [Matched linear models](methods/cohort-time-fe.md#match-r-etwfe-and-stata-jwdid) with unit fixed effects, never-treated or not-yet-treated controls, and optional time-invariant covariates. The comparison covers complete and incomplete panels. See the [supported specifications](#what-the-linear-etwfe-match-covers). |
| Sun–Abraham, R [fixest::sunab](https://lrberge.github.io/fixest/reference/sunab.html) and Stata [eventstudyinteract](https://github.com/lsun20/EventStudyInteract) | Same cohort-by-date preset and `control_group = "nevertreated"` | The saturated linear cohort-event specification with the same reference period and never-treated controls. We checked the complete-panel post-treatment cells against `sunab()`. To compare the event-study curve, [average the same cohort effects at each horizon](#how-event-study-averages-are-weighted). |
| Gardner's [did2s](https://kylebutts.github.io/did2s/) and Borusyak–Jaravel–Spiess [did_imputation](https://github.com/borusyak/did_imputation) in R and Stata | Same cohort-by-date preset, `control_group = "notyettreated", pt_assumption = "all", pt_start = min(panel$t)` | On a complete, unweighted panel, match the linear untreated regression with unit and date fixed effects, optionally interacting time-invariant covariates with dates. The post-treatment cohort-date effects are unrestricted. [Examples with and without covariates](#two-stage-and-imputation-event-studies). |
| Dube, Girardi, Jordà and Taylor's [LP-DiD](https://github.com/danielegirardi/lpdid) | Fit the [two-date comparison](#two-date-comparison-with-r-did) for the chosen baseline and horizon, or use the linear preset on a complete panel with one adoption cohort and never-treated controls | The unweighted clean-control long difference in that design. Match the risk sets and the horizon weights to compare a pooled staggered-adoption LP-DiD. Its other weighting and baseline options are outside this example. |
| Harmon's [stepwise DiD](https://web.econ.ku.dk/nharmon/docs/harmon2024onlineappendix.pdf) | Estimate each consecutive two-date comparison with `att_gt()` and add the estimates | Stepwise DiD is the consecutive-link case of chained DiD. Each step uses the cohorts still untreated at its endpoint. [Example with a changing set of controls](#harmons-stepwise-did). |
| Static binary switching, R [DIDmultiplegtSTAT](https://cran.r-project.org/package=DIDmultiplegtSTAT) and Stata [did_multiplegt_stat](https://github.com/Credible-Answers/did_multiplegt_stat) | `att_gt(dname = "d", ...)` | [Matched binary switcher effects](#binary-static-switcher-summaries) on complete panels with no covariates. The [static model](methods/static-entry-exit.md) also supports chained estimation with incomplete outcomes under no carryover. |
| Balla-Elliott and Norwich's [Cohort-Chained DiD](https://pdfs.dballaelliott.com/ccdid.pdf) | `control_group = "notyettreated", base_period = "universal", cross_pair = TRUE, pt_assumption = "ctrl", gmm_weighting = "within"` | Reproduces their stacked **cohort** fixed-effects point estimates with common observation dates within each cohort and the same sample and weights. [How the original chaining argument extends to cohorts](#cohort-chained-did). |

## What the linear ETWFE match covers

The preset matches the fully interacted linear panel regression using unit fixed effects, calendar effects and a separate effect for each treated cohort and date. The never-treated option also estimates pre-treatment deviations relative to the date before adoption. The not-yet-treated option fits a common untreated path from the first study date. For covariates, use `xformla = ~ x1 + x2` to allow both the untreated calendar path and the cohort-date effects to vary with time-invariant characteristics. These are point-estimate and covariance matches on the same estimation rows, after aligning the covariance normalization and the covariate values used to evaluate each effect. [R and Stata calls and validation](methods/cohort-time-fe.md#match-r-etwfe-and-stata-jwdid)

The preset requires a never-treated cohort, equal unit weights, zero anticipation and fixed covariates. The linear cohort-fixed-effects regression also gives the same point estimates in the corresponding complete-panel design. When observation dates vary within a cohort, use unit fixed effects in the reference command for this comparison. The upstream packages also offer linear specifications outside this preset, including time-varying regressors, sampling weights, additional fixed effects, and restrictions on treatment effects such as `jwdid`'s `hettype(event)`. That option estimates a common effect at each event time directly. Averaging the separate fitted cohort effects at that event time is a different operation. [R etwfe reference](https://grantmcdermott.com/etwfe/reference/etwfe.html) · [Stata jwdid help](https://github.com/friosavila/stpackages/blob/main/jwdid/jwdid.sthlp) · [Wooldridge, Sections 5 and 10](https://doi.org/10.1007/s00181-025-02807-z)

## How event-study averages are weighted

An event-study point at horizon 1 averages effects measured one period after adoption. For example, if the 2020 cohort has 60 units and an estimated effect of 1 in 2021, and the 2022 cohort has 40 units and an effect of 3 in 2023, the cohort-size-weighted point is `0.6 × 1 + 0.4 × 3 = 1.8`. If only 20 units in the later cohort have an observed outcome in 2023, an average weighted by those observed outcomes would instead be `(60 × 1 + 20 × 3) / 80 = 1.5`, even with the same two estimated cohort effects.

For the comparison, first check the cohort-date effects. Then average the same cohorts at each horizon with the same weights. In the unweighted-panel examples, `aggte(type = "dynamic")` weights each cohort by its study-unit count. If a requested effect is unavailable, it takes an explicit `na.rm = TRUE` to average over the supported effects. A regression post-estimation command can instead average over observed treated rows. On a complete panel, the counts coincide when the estimation samples and supported cohorts are the same. [Aggregate effects](reference/results.md#aggregated-results) · [Match an observed-row average](methods/cohort-time-fe.md#match-the-simple-summary)

## Two-date comparison with R did

Take a complete, unweighted panel with one adoption cohort and never-treated controls. Keep the baseline and one follow-up date, and code them as 1 and 2 with adoption at 2. With `base_period = "universal"` and the same time-invariant covariates, compare the follow-up cell. In this design the cohort-by-date preset gives the corresponding single comparison for `est_method = "reg"`, `"ipw"` or `"dr"`. With no covariates, it also equals the coefficient on the treated indicator in a regression of the long outcome change, the clean-control LP-DiD comparison for these two dates.

With more dates, a cohort-by-date fit combines supported changes under a common conditional untreated calendar path. The standard chained estimator instead follows the [chosen parallel-trends window](methods/parallel-trends.md). In Stata `csdid` and the experimental pair-balancing option of R [fastdid](https://tsailintung.github.io/fastdid/reference/fastdid.html), each two-date comparison has its own sample. `chain_report()` shows which pairs actually contribute to a chained effect. [R did transition guide](from-did.md)

## Two-stage and imputation event studies

On a complete panel, fitting unit and date effects on untreated observations gives the common calendar path for the match below. The treated observations have their own unrestricted cohort-by-date effects. Averaging their residuals by horizon with the same cohort weights gives the same post-treatment point estimates. The R example compares `cdid` with [didimputation](https://github.com/kylebutts/didimputation) on the same data.

```r
set.seed(14603)
panel <- expand.grid(t = 1:6, id = 1:72)
panel$g <- rep(c(0L, 3L, 5L), length.out = 72)[panel$id]
panel$y <- rnorm(72)[panel$id] + sin(panel$t) + rnorm(nrow(panel)) +
  ifelse(panel$g > 0 & panel$t >= panel$g, (panel$t - panel$g + 1) / 2, 0)

fit <- cdid::att_gt(
  yname = "y", tname = "t", idname = "id", gname = "g", data = panel,
  est_method = "reg", control_group = "notyettreated",
  base_period = "universal", pt_assumption = "all", pt_start = min(panel$t),
  bstrap = FALSE, cband = FALSE, control = list(preset = "cohort_time_fe")
)
es <- cdid::aggte(fit, type = "dynamic", bstrap = FALSE, cband = FALSE)
ref <- didimputation::did_imputation(panel, "y", "g", "t", "id", horizon = 0:3)
data.frame(horizon = es$egt, cdid = es$att.egt,
           didimputation = ref$estimate[match(as.character(es$egt), ref$term)])
```

| Horizon | cdid | didimputation |
| ---: | ---: | ---: |
| 0 | 0.607052 | 0.607052 |
| 1 | 0.516187 | 0.516187 |
| 2 | 1.918483 | 1.918483 |
| 3 | 1.807355 | 1.807355 |

Gardner's two-stage estimator also matches the six cohort-date point estimates in this panel. Give each post-treatment cohort-date cell its own second-stage indicator.

```r
panel$treated <- as.integer(panel$g > 0 & panel$t >= panel$g)
panel$cell <- ifelse(panel$treated == 1,
                     paste0("g", panel$g, "t", panel$t), "untreated")
two_stage <- did2s::did2s(
  panel, yname = "y", first_stage = ~ 0 | id + t,
  second_stage = ~ i(cell, ref = "untreated"),
  treatment = "treated", cluster_var = "id"
)
rbind(cdid = fit$att,
      did2s = coef(two_stage)[paste0("cell::g", fit$group, "t", fit$t)])
```

### Add time-invariant covariates

For a baseline covariate `x`, specify `xformla = ~ x` in `cdid`. Match the untreated models in Gardner and Borusyak–Jaravel–Spiess by interacting `x` with date. The unit effects absorb the time-invariant main effect of `x`. The following continues the complete-panel example.

```r
panel$x <- rnorm(72)[panel$id]
panel_x <- transform(panel, y = y + (0.4 * t + 0.2 * treated * (t - g + 1)) * x)

fit_x <- cdid::att_gt(
  yname = "y", tname = "t", idname = "id", gname = "g", data = panel_x,
  xformla = ~ x, est_method = "reg", control_group = "notyettreated",
  base_period = "universal", pt_assumption = "all", pt_start = min(panel_x$t),
  bstrap = FALSE, cband = FALSE, control = list(preset = "cohort_time_fe")
)
es_x <- cdid::aggte(fit_x, type = "dynamic", bstrap = FALSE, cband = FALSE)
ref_x <- didimputation::did_imputation(
  panel_x, "y", "g", "t", "id", first_stage = ~ i(t, x, ref = 1) | id + t,
  horizon = 0:3
)
two_stage_x <- did2s::did2s(
  panel_x, yname = "y", first_stage = ~ i(t, x, ref = 1) | id + t,
  second_stage = ~ i(cell, ref = "untreated"),
  treatment = "treated", cluster_var = "id"
)
c(
  didimputation = max(abs(es_x$att.egt -
    ref_x$estimate[match(as.character(es_x$egt), ref_x$term)])),
  did2s = max(abs(fit_x$att -
    coef(two_stage_x)[paste0("cell::g", fit_x$group, "t", fit_x$t)]))
)
# Largest absolute difference is about 5.8e-8 for each match.
```

These are post-treatment point-estimate matches. Other first-stage specifications, such as time-varying controls, require their own comparison. With incomplete observations, joint fixed-effects regression and untreated-only imputation can use different information to fit the untreated path. The weights and inference can differ as well. [Wooldridge, Section 5.5](https://doi.org/10.1007/s00181-025-02807-z) · [Borusyak, Jaravel and Spiess](https://doi.org/10.1093/restud/rdae007)

## Harmon's stepwise DiD

For cohort 3 and an effect at date 6, the stepwise estimator adds the DiDs for dates 2–3, 3–4, 4–5 and 5–6. Cohort 5 is eligible as a control for the first two steps. The later steps use the never-treated units. The code below imposes exactly these comparisons on the complete example above. All estimates use the same target cohort. Each two-date fit codes that cohort as treated in its second date and the eligible controls as untreated.

```r
g0 <- 3L
steps <- vapply(g0:6L, function(s) {
  keep <- panel$t %in% c(s - 1L, s) &
    (panel$g == g0 | panel$g == 0L | panel$g > s)
  pair <- panel[keep, ]
  pair$g <- ifelse(pair$g == g0, 2L, 0L)
  pair$t <- pair$t - s + 2L
  link <- cdid::att_gt(
    yname = "y", tname = "t", idname = "id", gname = "g", data = pair,
    est_method = "reg", base_period = "universal", bstrap = FALSE, cband = FALSE
  )
  link$att[link$t == 2L]
}, numeric(1))
sum(steps)
# [1] 1.696038
```

This reproduces the point estimate. To reproduce a stepwise result on an incomplete panel, use Harmon's target units eligible for the requested horizon and the same observed controls in each step. To compute uncertainty for the sum, account for covariance across steps because units can contribute to several comparisons. A single standard chained `att_gt()` fit can also combine longer pairs and additional control pools where supported. Its point estimate can then differ from the consecutive-link estimate. [Harmon, online appendix, Equation 12](https://web.econ.ku.dk/nharmon/docs/harmon2024onlineappendix.pdf)

## Cohort-Chained DiD

Cohort-Chained DiD is a simple extension of the original chained DiD argument. Balla-Elliott and Norwich proposed applying it across later-treated control cohorts whose pre-treatment windows are short. For example, one untreated cohort can supply the calendar change from dates 1 to 3 and another can supply the change from 3 to 5. Their stacked regression recovers the change from 1 to 5 through the date effects, even if no single untreated cohort spans it. Their paper cites the original chained DiD argument in footnote 2. [Chained DiD](https://doi.org/10.1016/j.jeconom.2024.105783) · [Cohort-Chained DiD, Section 2](https://pdfs.dballaelliott.com/ccdid.pdf)

In `cdid`, `cross_pair = TRUE` lets the target cohort and the controls supply changes over different pairs, `pt_assumption = "ctrl"` imposes the shared untreated calendar path among eligible controls, and `gmm_weighting = "within"` gives the regression weighting. With common observation dates within each cohort and a common regression sample and weights, these settings recover the paper's cohort fixed-effects point estimates. The cross-pair implementation also supports changes measured within individuals when their dates of observation differ, and allows [observation corrections](methods/missing-outcomes.md#how-to-adjust) for selective observation. [Cross-pair settings and equivalence](methods/cross-pair-chaining.md#options-for-chaining) · [Worked example with short pre-treatment windows](examples/short-pretreatment-windows.md)

## Binary static switcher summaries

For this comparison, supply a complete panel with columns `id`, `time`, `d` and `y`. Use the current binary treatment status `d` in each call. With no covariates, the matched summary averages switches into and out of treatment, expressed in the treated-minus-untreated direction and weighted by the number of switchers.

```r
fit <- cdid::att_gt(
  yname = "y", tname = "time", idname = "id", dname = "d",
  data = panel, est_method = "reg", xformla = ~1,
  switchers = "both", bstrap = FALSE, cband = FALSE
)
fit$aggregate

ref <- DIDmultiplegtSTAT::did_multiplegt_stat(
  df = panel, Y = "y", ID = "id", Time = "time", D = "d",
  estimator = "waoss", estimation_method = "ra", exact_match = FALSE
)
ref$results$table["WAOSS", ]
```

```stata
did_multiplegt_stat y id time d, estimator(was) exact_match graph_off
matrix list e(WAS)
```

The complete-panel example with 384 units and eight dates gives 1.211393 in all three programs. Standard errors depend on the variance normalization. With missing outcomes, `cdid` can also use changes from other supported dates in the static fit. [Static entry and exit](methods/static-entry-exit.md)

## How software uses incomplete panels

The original chained DiD paper formalized a missing-trends assumption. The baseline, called **missing trends completely at random (MTCAR)** here, asks whether units with both outcomes are representative of their cohort with respect to the **mean change** over that pair. This is weaker than missing outcomes completely at random (MCAR). For example, survey response can depend on a persistent difference in levels that cancels in changes. Missingness can still select changes if it depends on characteristics related to outcome growth. [Chained DiD paper, Assumption 2 and Corollary 1](https://doi.org/10.1016/j.jeconom.2024.105783)

The usual options in the matched commands change which observations enter the calculation. A repeated-cross-section DiD compares outcome levels from separate samples at each date, so it needs assumptions about how those samples represent the population. A pair-balanced DiD uses changes from units observed at both dates. To interpret those changes as population means requires mean-change representativeness, as in MTCAR. A linear regression on available rows uses the outcome model to link those rows. Its application to an incomplete sample still relies on the corresponding restrictions on the observed changes or levels. These sample choices do not by themselves fit a model of outcome observation. [R did documentation](https://bcallaway11.github.io/did/reference/att_gt.html) · [Stata csdid sample choices](https://psantanna.com/csdid/articles/comparing-did-commands.html) · [Wooldridge, Section 10](https://doi.org/10.1007/s00181-025-02807-z)

Of the implementations compared here, `cdid` is the one that combines within-unit chaining with explicit correction for covariate-dependent observation. Under **conditional mean trends missing at random**, the observed changes must represent the cohort's mean changes **at the same covariate values**. The corrections average over the study units, including those who are missing one or both outcomes for a given pair. This requires the relevant overlap or prediction support, in addition to the chosen parallel-trends assumption. [Preparing a panel](help/prepare-data.md) · [Handling missing outcomes](methods/missing-outcomes.md)

| How to handle observation | Set in `att_gt()` |
| --- | --- |
| Baseline MTCAR | `selection = "none"` |
| Weight changes by their probability of observation | `selection = "conditional_mtar", selection_est_method = "ipw"` |
| Predict missing actual changes | `selection = "conditional_mtar", selection_est_method = "imputation"` |
| Combine observation weighting and change imputation | `selection = "conditional_mtar", selection_est_method = "dr"` |

The same choices are available with outcome regression (`est_method = "reg"`), treatment weighting (`"ipw"`) and treatment doubly robust adjustment (`"dr"`) in the chained adoption and static switching models, including the cohort-by-date preset. That gives three treatment choices and four ways to handle observation. With both DR options, the treatment step combines an untreated-change model and a treatment propensity. The observation step combines an actual-change model and an observation propensity. The missing treated changes can include treatment effects. This is a different prediction from the untreated counterfactual imputed in the Borusyak–Jaravel–Spiess estimator. [Treatment adjustment](methods/treatment.md) · [Observation models](examples/observation-models.md)

## Measured runtime

These comparisons use matched regression specifications and analytical inference. The cohort-by-date runs use one thread and unit clustering. The table shows the complete and incomplete not-yet-treated examples from the larger benchmark, plus the complete static switcher example.

| Matched specification | Panel | cdid | R reference | Stata reference |
| --- | --- | ---: | ---: | ---: |
| Linear cohort-by-date, 10,000 units × 20 dates | Complete | 0.037 s | 1.245 s, etwfe | 2.739 s, jwdid |
| Linear cohort-by-date, 10,000 units × 20 dates | Incomplete | 0.036 s | 0.842 s, etwfe | 2.385 s, jwdid |
| Binary static regression, 384 units × 8 dates | Complete | 0.033 s | 0.149 s, DIDmultiplegtSTAT | 0.162 s, did_multiplegt_stat |

The cohort-by-date benchmark also covers never-treated controls and reports covariate fits, numerical agreement and peak memory. The examples above use no covariates or observation model. [Benchmark settings and results](methods/cohort-time-fe.md#verified-comparisons)

## Extensions

Dynamic switching with carryover, as in [DIDmultiplegtDYN](https://cran.r-project.org/package=DIDmultiplegtDYN) and Stata [did_multiplegt_dyn](https://github.com/Credible-Answers/did_multiplegt_dyn), is future work. The current binary switching route targets a static effect. Nonlinear cohort-by-date models, such as Poisson in `jwdid`, are also planned extensions. The LP-DiD match above covers the clean-control long difference in the stated design. A dedicated interface for its other horizon and variance-weighting choices would extend that match.

Other methods start from a different design or target. [DIDHAD](https://cran.r-project.org/package=DIDHAD) uses quasi-stayers, [contdid](https://bcallaway11.github.io/contdid/) targets dose responses, and [staggered](https://github.com/jonathandroth/staggered) assumes randomized treatment timing. [PanelMatch](https://cran.r-project.org/package=PanelMatch) matches treatment histories, [fect](https://yiqingxu.org/packages/fect/) also offers factor models, and [synthdid](https://github.com/synth-inference/synthdid) and Stata [sdid](https://github.com/Daniel-Pailanir/sdid) use synthetic weights. [HonestDiD](https://github.com/asheshrambachan/HonestDiD) starts from an estimated event study and its covariance to study sensitivity to parallel trends. Its input can come from a suitably aligned `cdid` event study.

If you would like your method, or another DiD method you use, added to `cdid`, please propose it in [cdid's GitHub issue tracker](https://github.com/dvbn/cdid/issues). A paper or package link and a short example of the target specification help us understand the request.

## References

- Balla-Elliott, D. and I. Norwich (2025). [Cohort-Chained DiD: Long-Run Effects with Limited Pre-Treatment Data](https://pdfs.dballaelliott.com/ccdid.pdf).
- Bellégo, C., D. Benatia and V. Dortet-Bernadet (2025). [The Chained Difference-in-Differences](https://doi.org/10.1016/j.jeconom.2024.105783). *Journal of Econometrics*.
- Borusyak, K., X. Jaravel and J. Spiess (2024). [Revisiting Event-Study Designs: Robust and Efficient Estimation](https://doi.org/10.1093/restud/rdae007). *Review of Economic Studies*.
- Callaway, B. and P. H. C. Sant’Anna (2021). [Difference-in-Differences with Multiple Time Periods](https://doi.org/10.1016/j.jeconom.2020.12.001). *Journal of Econometrics*.
- de Chaisemartin, C., X. D’Haultfœuille, F. Pasquier, D. Sow and G. Vazquez-Bare (2026). [Difference-in-Differences Estimators for Treatments Continuously Distributed at Every Period](https://arxiv.org/abs/2201.06898v7). The binary, static switcher comparison uses this framework.
- de Chaisemartin, C. and X. D’Haultfœuille (2026). [Difference-in-Differences Estimators of Intertemporal Treatment Effects](https://arxiv.org/abs/2007.04267v14). This is the reference for dynamic switching.
- Dube, A., D. Girardi, Ò. Jordà and A. M. Taylor (2025). [A Local Projections Approach to Difference-in-Differences](https://doi.org/10.1002/jae.70000). *Journal of Applied Econometrics*.
- Gardner, J. (2022). [Two-Stage Differences in Differences](https://arxiv.org/abs/2207.05943).
- Harmon, N. A. (2024). [Difference-in-Differences and Efficient Estimation of Treatment Effects](https://web.econ.ku.dk/nharmon/docs/harmon2022difference.pdf).
- Sun, L. and S. Abraham (2021). [Estimating Dynamic Treatment Effects in Event Studies with Heterogeneous Treatment Effects](https://doi.org/10.1016/j.jeconom.2020.09.006). *Journal of Econometrics*.
- Wooldridge, J. M. (2025). [Two-Way Fixed Effects, the Two-Way Mundlak Regression, and Difference-in-Differences Estimators](https://doi.org/10.1007/s00181-025-02807-z). *Empirical Economics*.
