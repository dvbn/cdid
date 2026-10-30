# Practical cdid guide

- [Editable Quarto source](cdid_practical_guide.qmd) for running and modifying the examples.
- [Download standalone HTML](https://github.com/dvbn/cdid/raw/refs/heads/main/notebooks/cdid_practical_guide.html) with the executed code, tables and plots.

Download the HTML file and open it in a browser; R and Quarto are not needed to read it. On GitHub, use the raw/download action to save the file rather than reading its HTML markup.

Open [cdid_practical_guide.qmd](cdid_practical_guide.qmd) in Positron or another Quarto-capable editor. The portable HTML embeds its resources. Its contents panel and final function map navigate 16 sections. Code and results appear together, with links to the [public user guide](https://github.com/dvbn/cdid/blob/main/docs/README.md).

Start a fresh R session after reinstalling `cdid`, `did` or `DRDID`. Run setup and the DGP chunk first, then work downward. Chunks run from the repository root or `notebooks/`; setup uses `validation/.lib` when present, otherwise the normal installed libraries, and prints exact versions and directory names without full local paths. It installs nothing. Required packages are `cdid`, `did`, `DRDID`, `knitr` and `ggplot2`; rendering also requires Quarto.

Edit `cfg` in Section 2 to control the main draw (default seed `20261029L`, 2000 study units, six dates). `simulate_panel.R` generates data and truth benchmarks: `population_att` describes the superpopulation distribution; `study_unit_latent_mean` displays the latent mean effect among the realized study units in `units` (`frame = units`). The helper retains its legacy API field `roster_latent_mean`; the notebook renames columns only in display copies. Balanced, incomplete and rotating panels share latent outcomes. Run Section 3 before estimation to inspect population denominators, zero- and single-outcome counts, and pair support. Inspect fitted target availability with `cdid::tidy(fit_to_explore)` in Section 7, after the fits exist. Availability does not verify scientific assumptions; unavailable (`NA`) is never zero.

After changing `cfg`, redraw once and refit downstream examples. When changing only an estimator, refit the same data and rebuild its aggregates. Keep errors and unavailable cells visible.

The clustered and continuous-X robustness examples use separate explicitly seeded studies. The multiplier example uses an editable exploratory 199 draws. The missingness, treatment-adjustment and imputation comparisons display six-decimal SEs. A single draw does not establish bias, coverage or efficiency; nearby SEs do not validate assumptions or inference. Weighted examples have weighted targets and do not reuse the unweighted truth table.

Render from the repository root:

```bash
quarto render notebooks/cdid_practical_guide.qmd --to html
```

Unexpected errors stop execution. There are 33 R chunks: 29 run automatically and four are disabled (`eval: false`): the opt-in `reg` (OR)/`ipw` comparison on the three continuous-X robustness draws (Section 12), expected aggregation refusal after a point-only fit (Section 14), optional parallel multiplier computation (Section 15), and exports to `guide-output/` (Section 16). Run these individually when wanted. The custom-estimator section documents the supported callable contract. No PDF/TeX compilation is needed.
