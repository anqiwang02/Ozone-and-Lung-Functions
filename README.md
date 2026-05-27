# Ozone and Lung Functions

Statistical analysis of DNA methylation at CpG sites associated with lung-function genes under ozone (O3) versus clean air (CA) exposure. The study uses a paired design with **17 individuals** and applies randomization-based (Fisherian) inference alongside classical (Neymanian) confidence intervals.

## Overview

This repository contains R code to analyze methylation changes at CpG sites linked to **primary** and **secondary** lung-function gene lists. For each CpG site, the analysis computes:

1. **Fisher exact p-values** — unadjusted, randomization-based p-values from the full set of treatment assignments
2. **P-value adjustment** — fully randomization-based adjustment for multiple testing and Benjamini–Hochberg FDR
3. **Fisherian confidence intervals** — inversion of randomization tests over a grid of effect sizes
4. **Neymanian confidence intervals** — paired *t*-based 95% intervals for the average causal effect (ACE)

## Study Design

- **Sample:** 17 paired subjects, each observed under ozone (`exp = 1`) and clean air (`exp = -1`)
- **Outcome:** CpG site methylation levels
- **Test statistic:** Absolute value of the mean difference statistic under all possible randomizations of exposure assignments
- **Randomization reference:** `Data/w0true.Rdata` contains `w0.true`, a matrix of all 2^17 = 131,072 possible randomizations of the paired exposure vector

## Statistical Methods

### Fisher exact p-values

For each CpG site, the observed test statistic is compared against the empirical null distribution generated from all randomizations in `w0.true`. The unadjusted Fisher p-value is the proportion of randomization-based statistics at least as extreme as the observed statistic:

```
p = (number of randomization statistics ≥ observed statistic) / 2^17
```

Implemented in **Section 1** of the main analysis notebooks.

### P-value adjustment

Two adjustment procedures are applied to control for multiple comparisons across CpG sites:


| Method                                   | Description                                                                                                                                                                                         |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Fully randomization-based adjustment** | Pool test statistics from the first 1,000 randomizations across all CpG sites, compute minimum p-values per iteration, and calibrate each observed p-value against this null distribution of minima |
| **FDR adjustment**                       | Benjamini–Hochberg FDR applied to unadjusted Fisher p-values via `p.adjust(..., method = "fdr")`                                                                                                    |


Implemented in **Sections 2–3** of the main analysis notebooks.

### Fisherian confidence intervals

Fisherian intervals are constructed by inverting the randomization test. For a grid of candidate effect sizes *a*, potential outcomes are imputed under the sharp null H₀: Yi(1) − Yi(0) = *a*, and a paired *t*-statistic is computed under each randomization. The 95% Fisherian interval consists of all *a* values for which the randomization p-value exceeds 0.025 (two-sided).

Implemented in the Fisherian interval scripts and graphing utilities listed below.

### Neymanian confidence intervals

Neymanian intervals use the standard paired *t* interval for the mean difference (ozone minus clean air):

```
ACE ± t_{0.025, df=16} × (s_d / √17)
```

Applied to the top-ranked CpG sites (top 9 primary, top 5 secondary) by unadjusted Fisher p-value.

Implemented in the Neyman interval scripts listed below.

## Repository Structure

```
Ozone-and-Lung-Functions/
├── Data/
│   └── w0true.Rdata              # All 2^17 randomization matrices
├── primary_genes/
│   ├── Exp_data_primary_CpG.RData
│   ├── primary gene and CPG updated.RData
│   ├── lung_genes_primary_results.xlsx   # Output: p-values and effect sizes
│   ├── top9_95CI_primary.xlsx            # Output: Neymanian CIs (top 9 CpG sites)
│   ├── fisher_int_primary_genes.xlsx     # Output: Fisherian CIs
│   └── code/
│       ├── Analysis for lung genes_primary git.Rmd   # Main analysis (p-values & adjustment)
│       ├── Fisherian interval primary genes.Rmd      # Fisherian CIs (worked examples)
│       ├── Neyman interval primary genes.Rmd         # Neymanian CIs (top 9 sites)
│       ├── fisher_int_prim_graph.R                   # Fisherian CI plots
│       └── fisher_int_prim_graph.Rmd
├── secondary_genes/
│   ├── Exp_data_secondary_CpG.RData
│   ├── secondary gene and CPG.RData
│   ├── lung_genes_secondary_results.xlsx
│   ├── top5_95CI_secondary.xlsx
│   ├── fisher_int_secondary_genes.xlsx
│   └── Code/
│       ├── analysis for lung genes_secondary git.Rmd
│       ├── Fisherian interval sec genes.Rmd
│       ├── Neyman interval secondary genes.Rmd
│       └── fisher_int_sec_graph.Rmd
└── README.md
```

## Requirements

- R (≥ 4.0 recommended)
- **dplyr**, **tidyr**, **readxl**, **writexl**, **ggplot2**, **purrr**, **tidyverse** (for secondary analysis)

Install dependencies:

```r
install.packages(c("dplyr", "tidyr", "readxl", "writexl", "ggplot2", "purrr", "tidyverse"))
```

## Usage

1. Clone the repository and set your working directory to the project root.
2. Update hard-coded file paths in the `.Rmd` and `.R` scripts to match your local directory (paths currently reference `/Users/anqiwang/Documents/GitHub/Ozone-and-Lung-Functions/`).
3. Run the analysis in order:


| Step | Script                                                                                 | Output                                                                                         |
| ---- | -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| 1    | `primary_genes/code/Analysis for lung genes_primary git.Rmd` (or secondary equivalent) | Unadjusted Fisher p-values, randomization-based adjusted p-values, FDR p-values, volcano plots |
| 2    | `primary_genes/code/Fisherian interval primary genes.Rmd`                              | Fisherian 95% CIs for selected CpG sites                                                       |
| 3    | `primary_genes/code/Neyman interval primary genes.Rmd`                                 | Neymanian 95% CIs for top-ranked CpG sites                                                     |
| 4    | `fisher_int_prim_graph.R` / `.Rmd`                                                     | Fisherian interval visualization                                                               |


Open `.Rmd` files in RStudio and knit, or source `.R` scripts from the R console.

## Key Output Files


| File                                                                | Contents                                                                                                                    |
| ------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `lung_genes_primary_results.xlsx`                                   | CpG site, gene name, randomization-adjusted p-value, unadjusted Fisher p-value, FDR-adjusted p-value, mean difference (ACE) |
| `lung_genes_secondary_results.xlsx`                                 | Same structure for secondary gene list (283 CpG sites)                                                                      |
| `top9_95CI_primary.xlsx` / `top5_95CI_secondary.xlsx`               | Neymanian lower and upper bounds                                                                                            |
| `fisher_int_primary_genes.xlsx` / `fisher_int_secondary_genes.xlsx` | Fisherian interval bounds                                                                                                   |


## Notes

- Primary gene list: **232** CpG sites; secondary gene list: **283** CpG sites.
- Randomization-based adjustment uses the first **1,000** rows of `w0.true` for computational efficiency; exact Fisher p-values use all **131,072** randomizations.
- Pre-rendered `.nb.html` notebooks are included alongside source `.Rmd` files for reference.

