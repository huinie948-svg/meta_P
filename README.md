# Microbial inoculants and soil phosphorus: data and analysis

This repository contains the data, R code, and selected analysis outputs for:

> Context-dependent effects of microbial inoculants on soil phosphorus
> availability and phosphatase activity: a global meta-analysis

The analysis includes 308 studies and 3007 observations. The three primary
soil endpoints are available phosphorus (AP), phosphatase activity, and total
phosphorus (TP).

## Effect-size convention

The primary effect size is the natural-log response ratio:

```text
yi_lnRR = ln(mean_inoculated / mean_control)
```

Its sampling variance is `vi_lnRR`. In the endpoint scripts, `RR` and `Vi` are
retained only as legacy aliases for `yi_lnRR` and `vi_lnRR`. `RR` is therefore
on the lnRR scale and must not be exponentiated before model fitting.

## Repository structure

```text
AP/                         AP data and endpoint analysis
phosphatase/                phosphatase data and endpoint analysis
TP/                         TP data and endpoint analysis
SEM/                        data and exploratory piecewise SEM analysis
map/                        study-location data and map script
robustness_results/         selected sensitivity-analysis outputs
robust_sensitivity_analysis.R
data_dictionary_effect_sizes.csv
```

Main input files:

- `AP/ap.csv`
- `phosphatase/phophatase.csv`
- `TP/TP.csv`
- `AP/ap_effects.csv`
- `phosphatase/ap_effects.csv`
- `TP/ap_effects.csv`
- `SEM/data_collection.csv`

## Analysis scripts

- `AP/ap.R`: AP effect sizes, overall and subgroup meta-analysis,
  dose-response analyses, correlations, and random forest importance.
- `phosphatase/phosphatase.R`: phosphatase analyses.
- `TP/TP.R`: TP analyses.
- `robust_sensitivity_analysis.R`: StudyID random-intercept models, assumed
  within-study sampling correlations, and missing-variance sensitivity
  analyses reported in Supplementary Table S1.
- `SEM/SEM.R`: exploratory observational piecewise SEMs.
- `map/map.R`: geographic distribution map.

The random forest models are used for predictor-importance analysis. The
piecewise SEMs are exploratory observational path analyses and should not be
interpreted as causal mediation tests.

## Running the analyses

Run the endpoint scripts from their respective folders because they use
folder-relative input paths.

Run the robustness analysis from the repository root:

```bash
Rscript robust_sensitivity_analysis.R
```

To run selected endpoints:

```bash
ENDPOINTS=AP Rscript robust_sensitivity_analysis.R
ENDPOINTS=AP,Phosphatase,TP Rscript robust_sensitivity_analysis.R
```

Valid endpoint names are `AP`, `Phosphatase`, and `TP`. Outputs are written to
`robustness_results/`.

The endpoint scripts retain interactive `file.choose()` calls in some
dose-response sections. Select the unit-specific CSV whose filename matches
the object being created in that section. These interactive sections are not
used by `robust_sensitivity_analysis.R`.

## Key outputs

- `robustness_results/overall_sensitivity_results_all.csv`
- `robustness_results/session_info.txt`

## Software

The sensitivity analysis was verified with:

- R 4.5.2
- metafor 4.8-0
- Matrix 1.7-4
- nlme 3.1-168

The complete session information is provided in
`robustness_results/session_info.txt`.

## Data interpretation

Stress and fertilization are broad operational categories that combine
different experimental conditions. The yield SEM contains 27 observations
from three studies and is not used as primary mechanistic evidence.

## Archived version

The fixed archive associated with the manuscript is available at
https://doi.org/10.6084/m9.figshare.31894483.
