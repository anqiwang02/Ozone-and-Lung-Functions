# Ozone and Lung Functions

Statistical analysis of DNA methylation at CpG sites associated with lung-function genes under ozone (O3) versus clean air (CA) exposure. The study uses a paired crossover design with **17 individuals** and combines randomization-based (Fisherian) and classical paired-*t* inference.

## Overview

This repository analyzes CpG methylation changes for **primary** and **secondary** lung-function gene lists. The main analysis notebooks now use **paired absolute t-statistics** as the randomization test statistic and compute:

1. **Unadjusted Fisher exact p-values** from all 2^17 randomizations
2. **Multiple-testing adjustments**:
   - fully randomization-based adjusted p-values
   - Benjamini-Hochberg FDR
3. **Mean differences (ACE)** (ozone minus clean air)
4. **Diagnostic figures** (volcano plots, randomization null histograms, asymptotic |t| overlays, correlation plots)

Separate scripts/notebooks compute Fisherian and Neymanian confidence intervals for selected CpGs.

## Study Design

- **Sample:** 17 paired subjects, each observed once under ozone (`exp = 1`) and once under clean air (`exp = -1`)
- **Outcome:** CpG methylation level
- **Randomization reference:** `Data/w0true.Rdata` contains `w0.true`, the full matrix of all **131,072 (2^17)** paired randomizations
- **Main test statistic:** paired \|t\|:

```text
|T| = | d̄ / (s_d / sqrt(17)) |
```

where `d̄` and `s_d` are the mean and sample SD of pair-level differences.

## Statistical Methods

### 1) Fisher exact p-values (main notebooks)

For each CpG, compare observed paired \|t\| to the randomization null from all rows of `w0.true`:

```text
p_unadjusted = mean(T_rand_abs >= T_obs_abs)
```

Implemented in:
- `primary_genes/code/Analysis for lung genes_primary git.Rmd`
- `secondary_genes/Code/analysis for lung genes_secondary git.Rmd`

### 2) Multiple-testing adjustment

- **Fully randomization-based adjustment:** Uses first 1,000 randomizations per CpG, builds a null distribution of minimum p-values across CpGs, and calibrates each unadjusted p-value against that distribution.
- **FDR adjustment:** `p.adjust(..., method = "fdr")`.

### 3) Fisherian confidence intervals

Constructed by inverting randomization tests over a grid of constant-effect values `a` under the sharp null `Yi(1) - Yi(0) = a`, using paired signed *t* statistics.

Implemented in Fisherian interval `.Rmd` notebooks and graph scripts:
- `primary_genes/code/fisher_int_prim_graph.R`
- `secondary_genes/Code/fisher_int_sec_graph.R`

### 4) Neymanian confidence intervals

Paired-*t* intervals for ACE:

```text
ACE ± t_(0.025,16) * (s_d / sqrt(17))
```

## Repository Structure

```text
Ozone-and-Lung-Functions/
├── Data/
│   └── w0true.Rdata
├── primary_genes/
│   ├── Exp_data_primary_CpG.RData
│   ├── primary gene and CPG updated.RData
│   ├── lung_genes_primary_results.xlsx
│   └── code/
│       ├── Analysis for lung genes_primary git.Rmd
│       ├── Fisherian interval primary genes.Rmd
│       ├── Neyman interval primary genes.Rmd
│       ├── fisher_int_prim_graph.R
│       ├── fisher_int_prim_graph.Rmd
│       ├── run_fisher_int_prim_graph.sh
│       └── graphs/
├── secondary_genes/
│   ├── Exp_data_secondary_CpG.RData
│   ├── secondary gene and CPG.RData
│   ├── lung_genes_secondary_results.xlsx
│   └── Code/
│       ├── analysis for lung genes_secondary git.Rmd
│       ├── Fisherian interval sec genes.Rmd
│       ├── Neyman interval secondary genes.Rmd
│       ├── fisher_int_sec_graph.R
│       ├── fisher_int_sec_graph.Rmd
│       ├── run_fisher_int_sec_graph.sh
│       └── graphs/
├── .gitignore
└── README.md
```

## Requirements

- R (>= 4.0 recommended)
- Core packages used in current analysis/scripts:
  - `dplyr`, `tidyr`, `readxl`, `writexl`, `ggplot2`, `purrr`, `parallel`
  - `tidyverse` (some plotting sections)
  - `corrplot`, `rstatix`, `hrbrthemes`, `viridis`, `esquisse`

Install example:

```r
install.packages(c(
  "dplyr","tidyr","readxl","writexl","ggplot2","purrr","parallel",
  "tidyverse","corrplot","rstatix","hrbrthemes","viridis","esquisse"
))
```

## Usage

1. Clone the repository and set working directory to project root.
   - When running any `.Rmd` notebook, make sure the working directory is `Ozone-and-Lung-Functions` (the repository root).
2. Run main analyses:
   - `primary_genes/code/Analysis for lung genes_primary git.Rmd`
   - `secondary_genes/Code/analysis for lung genes_secondary git.Rmd`
3. Run interval analyses as needed:
   - Fisherian: corresponding Fisherian interval notebooks or graph scripts
   - Neymanian: corresponding Neyman interval notebooks
4. Optional scripted Fisher graph runs (for reproducing graph outputs in batch):
   - `bash primary_genes/code/run_fisher_int_prim_graph.sh`
   - `bash secondary_genes/Code/run_fisher_int_sec_graph.sh`
   - Set core count with `FISHER_N_CORES` (default: 10):
     - `FISHER_N_CORES=12 bash primary_genes/code/run_fisher_int_prim_graph.sh`
     - `FISHER_N_CORES=12 bash secondary_genes/Code/run_fisher_int_sec_graph.sh`

## Key Output Files

- Main analysis tables:
  - `primary_genes/lung_genes_primary_results.xlsx`
  - `secondary_genes/lung_genes_secondary_results.xlsx`
- Figure outputs:
  - `primary_genes/code/graphs/*.pdf`
  - `secondary_genes/Code/graphs/*.pdf`
- Interval summary tables:
  - Neymanian CIs: `primary_genes/top9_95CI_primary.xlsx`, `secondary_genes/top5_95CI_secondary.xlsx`
  - Fisherian interval `.xlsx` files in the repo are manually curated by running the .Rmd files.

## Notes

- Primary list: **232** CpGs; secondary list: **283** CpGs.
- Unadjusted Fisher p-values use all **131,072** randomizations.
- Fully randomization-based adjusted p-values use first **1,000** randomizations for computational efficiency.
- Some scripts still include absolute/local paths; update `setwd()` and file paths for your environment if needed.

