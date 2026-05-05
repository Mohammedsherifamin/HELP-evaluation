# HELP Synthetic Control Evaluation Analytic Code

Analytic code for the synthetic control evaluation of the [Health Equity Liverpool Project (HELP)](https://www.lstmed.ac.uk/sites/default/files/FINAL-HELP-report-2024.pdf) childhood vaccination intervention.

This repository corresponds to the following preprint:

> Amin MS, Zhang X, Green MA, Holford D, Hemingway C, Ismail A, Essale N, Doyle V, Taegtmeyer M, Hungerford D. *Evaluating the impact of a community-engagement intervention on the uptake of childhood vaccines in England: A synthetic control analysis.* medRxiv 2026.05.01.26352232. doi: https://doi.org/10.64898/2026.05.01.26352232
>
> Full text: https://www.medrxiv.org/content/10.64898/2026.05.01.26352232v1.full.pdf+html

The study protocol is published separately:

> Amin MS, Zhang X, Green MA, et al. *Impact of a community-led intervention on the uptake of childhood vaccines in Liverpool: a protocol for a synthetic control evaluation.* BMJ Open 2026;16:e111500. doi: 10.1136/bmjopen-2025-111500
>
> Full text: https://bmjopen.bmj.com/content/16/1/e111500

---

## Repository structure

```
.
├── HELP_analysis.Rmd          # Main analysis workflow
├── scripts/
│   ├── functions.R            # Helper functions (plotting, model extraction, etc.)
│   ├── count_0.5km.R          # Sensitivity analysis: 0.5 km buffer
│   ├── count_1.5km.R          # Sensitivity analysis: 1.5 km buffer
│   ├── matching_quality.R     # Pre-intervention balance diagnostics across buffer distances
│   └── weights_map.R          # Figure 5: spatial map of synthetic control weights
├── data/                      # Input data (see Data section below)
└── outputs/                   # Generated figures, tables, and rendered documents
```

### `HELP_analysis.Rmd`

The main analysis file. It runs end-to-end from raw data inputs through to all results, figures, and tables in the paper. It is structured as a sequential pipeline: each section saves an intermediate `.rds` file that the next section reads. Running it top-to-bottom will reproduce the full analysis.

### `scripts/`

- **`functions.R`**: helper functions used throughout the analysis: `extract_syn_results()`, `extract_svy_ratios()`, `plot_svy_ratios()`, `forest_plot_results()`, `run_microsynth()`, `extract_effect()`, and `clean_fingertips_data()`. Sourced by `HELP_analysis.Rmd` and the two standalone scripts below.
- **`count_0.5km.R`**: sensitivity analysis code for the 0.5 km buffer (Supplementary Figure 3). Called from within `HELP_analysis.Rmd`; not intended to be run as a standalone script.
- **`count_1.5km.R`**: sensitivity analysis code for the 1.5 km buffer (Supplementary Figures 4–5). Called from within `HELP_analysis.Rmd`; not intended to be run as a standalone script.
- **`matching_quality.R`** (**run separately after `HELP_analysis.Rmd`**): Loops over buffer distances from 0.5 to 3 km in 0.1 km increments, runs a microsynth model at each, extracts the epsilon statistic, and produces the matching quality plot and table (Supplementary Figure 6 and Supplementary Table 1). Requires `data/analysis_data.rds` and `data/GP_MMR_IMD_ICB_POP_FING_ETHN_dist.rds` to already exist.
- **`weights_map.R`** (**run separately after `HELP_analysis.Rmd`**): Generates Figure 5: the England map showing non-intervention practices coloured by their synthetic control matching weight. Requires `data/analysis_data.rds`, `data/interv_gp/1km_interv.rds`, and `data/GP_MMR_IMD_ICB_POP_FING_ETHN_dist.rds` to already exist.

---

## Data

All data used in this analysis are from publicly available sources. The `data/` folder is organised as follows:

```         
data/
├── uptake/                    # COVER quarterly .ods files (Q1–Q24)
├── IMD/                       # Deprivation data
├── epcn/                      # PCN and ICB lookup files
├── gp_population/             # NHS Digital GP registered population CSVs
├── fingertips_data/           # Fingertips indicator CSVs
├── ethnicity/                 # Census 2021 ethnicity data
├── interv_gp/                 # Spatial files and intervention group definitions
└── flowchart/                 # Saved flow diagram object
```

Most data folders are provided as compressed `.zip` files in this repository. Three large spatial files exceed GitHub's file size limits and must be downloaded separately.

### Step 1: Unzip the included data

After cloning the repository, unzip all `.zip` files in the `data/` folder in place. Each zip expands into the subfolder it corresponds to (e.g. `uptake.zip` → `data/uptake/`).

### Step 2: Download three large spatial files

The following three datasets are too large to include in the repository. Download each one from the links below and place the extracted folder directly inside `data/interv_gp/` before running the analysis.

------------------------------------------------------------------------

#### NSPL Online Latest Centroids

National Statistics Postcode Lookup — postcode point centroids for all UK postcodes (\~500 MB).

**Source:** [ONS Open Geography Portal — NSPL Online Latest Centroids](https://geoportal.statistics.gov.uk/datasets/7a373e28348a4d7aac9892e18fae5706_0/explore?location=55.323015%2C-3.307485%2C6)

Expected location after extraction:

```         
data/interv_gp/NSPL_Online_Latest_Centroids/NSPL_Online_Latest_Centroids.shp
```

(plus the associated `.dbf`, `.prj`, `.shx` sidecar files in the same folder)

------------------------------------------------------------------------

#### Lower Layer Super Output Areas (December 2021) Boundaries EW BFC

LSOA 2021 full boundaries shapefile for England and Wales (\~300 MB).

**Source:** [ONS Open Geography Portal — Lower Layer Super Output Areas (December 2021) Boundaries EW BFC](https://geoportal.statistics.gov.uk/datasets/2bbaef5230694f3abae4f9145a3a9800_0/explore?location=52.837550%2C-2.489483%2C7)

Expected location after extraction:

```         
data/interv_gp/Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10)/
    Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10).shp
```

------------------------------------------------------------------------

#### LSOA (December 2021) Population Weighted Centroids

LSOA 2021 population-weighted centroids for England and Wales.

**Source:** [ONS Open Geography Portal — Lower Layer Super Output Areas (December 2021) Population Weighted Centroids](https://geoportal.statistics.gov.uk/datasets/32729e42d05e4e23bc7e43a36aa4ae8b_0/explore?location=52.234729%2C1.095673%2C6)

Expected location after extraction:

```         
data/interv_gp/LSOA_PopCentroids_EW_2021_V4_-8229881830311519574/
LSOA_PopCentroids_EW_2021_V4.shp
```

------------------------------------------------------------------------

## Reproducing the analysis

1.  Clone this repository.
2.  Unzip all `.zip` files in `data/` as described in Step 1 above.
3.  Download and place the three large spatial files as described in Step 2 above.
4.  Open `HELP_analysis.Rmd` in RStudio (R version 4.4.2 or later recommended).
5.  Install required packages if needed. All `library()` calls are at the top of the Rmd.
6.  Run chunks sequentially from top to bottom, or knit the document.

Once the main Rmd has completed, `matching_quality.R` and `weights_map.R` can be run as standalone scripts — they depend on intermediate `.rds` files saved by the Rmd. `count_0.5km.R` and `count_1.5km.R` are sourced from within the Rmd and do not need to be run separately.

------------------------------------------------------------------------

## Session info

R version 4.4.2. Key packages: `microsynth`, `tidyverse`, `sf`, `patchwork`, `cowplot`, `ggspatial`, `forestploter`, `gt`, `survey`, `janitor`. Full session info is printed at the end of `HELP_analysis.Rmd`.

------------------------------------------------------------------------

## Licence

Code is released under the MIT Licence. Data files sourced from third parties remain subject to their original licences (all are open access; see source links above).
