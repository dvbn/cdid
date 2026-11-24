# Practical cdid guide

This editable Quarto notebook works through the main [user guide](https://github.com/dvbn/cdid/blob/main/docs/README.md) topics with synthetic examples: a first balanced-panel fit, missing-outcome corrections, treatment adjustment, cross-pair chains, inference, the static switching model and the cohort-by-date preset.

- [Editable Quarto source](cdid_practical_guide.qmd) for running and modifying the examples.
- [Download standalone HTML](https://github.com/dvbn/cdid/raw/refs/heads/main/notebooks/cdid_practical_guide.html) with the executed code, tables and plots.

## Read it without R

Download the HTML file and open it in a browser. R and Quarto are not needed to read it. On GitHub, use the raw or download action to save the file rather than reading its HTML markup. The portable HTML embeds its resources. Its contents panel and the final function map navigate the 18 sections. Code and results appear together, with links to the public user guide.

## Run it

Open [cdid_practical_guide.qmd](cdid_practical_guide.qmd) in Positron or another Quarto-capable editor.

**Version.** The notebook uses the built-in chained cohort default `gmm_weighting = "within"`. Check that your installed version has that default. Install the current release from CRAN when available or use `remotes::install_github("dvbn/cdid")`. To test a development checkout, run `remotes::install_local(".")` from its root. The optional [pinned validation setup](../validation/SETUP.md) installs exact dependency versions. Start a fresh R session after reinstalling `cdid`, `did` or `DRDID`.

**Packages.** Required packages are `cdid`, `did`, `DRDID`, `knitr` and `ggplot2`. Rendering also requires `rmarkdown` and Quarto. The setup chunk installs nothing. Install the tutorial packages once if needed:

```r
install.packages(c("knitr", "rmarkdown", "ggplot2"))
```

**Order.** Run setup and the data-generation chunk first, then work downward. Chunks run from the repository root or from `notebooks/`. Setup uses `validation/.lib` when present, otherwise the normal installed libraries, and prints exact versions and directory names without full local paths.

## Follow the learning path

| Sections | Stage | What they cover |
| --- | --- | --- |
| 1 to 3 | Set up and inspect | Loaded versions, the editable configuration `cfg`, the draw of study units, the truth benchmarks, and raw missingness and pair-support counts before any fit. |
| 4 | First fit | Balanced-panel ATT(g,t) with `gname`, the four aggregations and the plots. |
| 5 to 8 | Missing outcomes | Why `selection = "none"` lacks its assumption in this design, observation models, treatment IPW, OR and DR, slope pooling, rotating panels, fitted support, change imputation, and the augmented correction with its model reports. |
| 9 | Comparisons and trends | Not-yet-treated controls, `base_period`, PT-all, cross-pair chains on disjoint windows and ordinary `did` delegation. |
| 10 to 12 | Inference and robustness | Multiplier bands, clustering, both pretrend tests, and treatment-DR correctness cases. |
| 13 and 14 | Inputs and weighting | Dated covariates, weight anchors, a comparison of within and two-step GMM on the same incomplete panel, resource controls and point-only fits. |
| 15 and 16 | Tools | Parallel multiplier draws, the custom estimator contract, supported limits, exports and an experiment workflow. |
| 17 and 18 | Other treatment designs | Binary treatment that can stop and restart with `dname`, and the cohort-by-date fixed-effects preset. |

Sections 17 and 18 draw their own data and depend only on the setup chunk. Read them after Section 9 if your treatment can stop or if you want a linear fixed-effects fit.

## Change the experiment

Edit `cfg` in Section 2 to control the main draw (default seed `20261029L`, 2000 study units, six dates). `simulate_panel.R` generates data and truth benchmarks: `population_att` describes the superpopulation distribution, and `study_unit_latent_mean` displays the latent mean effect among the realized study units in `units` (`frame = units`). The helper retains its legacy API field `roster_latent_mean`, and the notebook renames columns only in display copies. Balanced, incomplete and rotating panels share latent outcomes. Run Section 3 before estimation to inspect population denominators, zero- and single-outcome counts, and pair support. Inspect fitted target availability with `cdid::tidy(fit_to_explore)` in Section 7, after the fits exist. Availability does not verify scientific assumptions, and unavailable (`NA`) is never zero.

After changing `cfg`, redraw once and refit downstream examples. When changing only an estimator, refit the same data and rebuild its aggregates. Keep errors and unavailable cells visible.

In Section 8, the DR+DR model-report example shows the fitted equation and coefficient SEs for the observation logit, actual-change imputation, treatment logit and control outcome regression. It requires a build with `nuisance_models()`. Section 9 compares same-pair and cross-pair support on disjoint observation windows and prints the fitted links when the build supports `cross_pair`. Section 13 shows how to retain a known reference covariate on a row with a missing outcome and how to read the warning and counts when the covariate is unknown. Section 17 compares balanced and chained estimates of the static switching effects for a binary treatment that can stop and restart. It requires a build with `att_gt(dname = ...)`. The clustered and continuous-X robustness examples use separate explicitly seeded studies. The multiplier example uses an editable exploratory 199 draws. The missingness, treatment-adjustment and imputation comparisons display six-decimal SEs. A single draw does not establish bias, coverage or efficiency, and nearby SEs do not validate assumptions or inference. Weighted examples have weighted targets and do not reuse the unweighted truth table.

## Render

Render from the repository root:

```bash
quarto render notebooks/cdid_practical_guide.qmd --to html
```

Unexpected errors stop execution. There are 43 R chunks. 39 run automatically and four are disabled (`eval: false`): the opt-in `reg` (OR) and `ipw` comparison on the three continuous-X robustness draws (Section 12), the expected aggregation refusal after a point-only fit (Section 14), the optional parallel multiplier computation (Section 15), and the exports to `guide-output/` (Section 16). Run these individually when wanted. The custom-estimator section documents the supported callable contract. No PDF or TeX compilation is needed.
