# Modifiable burden of bacterial antimicrobial resistance associated with human antibiotic consumption

This repository contains the analysis code for the paper:

> **Modifiable burden of bacterial antimicrobial resistance associated with human antibiotic consumption**

## Overview

Human antibiotic use contributes to the evolution and maintenance of antimicrobial resistance (AMR). This analysis estimates this relationship, and thererfore the amount of AMR-attributable mortality that could be prevented by reducing antibiotic consumption, using global surveillance data on resistance prevalence across 17 bacterial pathogens and 7 antibiotic classes.

## Repository structure

```
config.R                          # Central parameters, paths, and scenario definitions
data_processing.r                 # Data preparation: harmonization, covariate joins
regression_models.r               # Model fitting: pathogen-antibiotic and class-level models
avertable_burden.R                # Burden estimation: avertable mortality calculations (writes Outputs CSVs)
plotting.r                        # Figure generation: Figures 1–4 and supplementary figures
generate_supplementary_table3.R   # Supplementary Table 3: variance explained by variable
scripts/
  run_prepare_data.R              # Stage entry: data preparation
  run_fit_models.R                # Stage entry: model fitting
  run_estimate_burden.R           # Stage entry: burden estimation
  run_generate_figures.R          # Stage entry: figure generation
  run_tables.R                    # Stage entry: supplementary tables
  run_wilcoxon_test.R             # Wilcoxon permutation test (Table 2 prerequisite)
  run_pipeline.R                  # End-to-end pipeline runner
  promote_metrics.R               # Promote passing pending metrics to canonical
  validate_outputs.R              # Manifest-based output validation
  validate_stage_data.R           # Merged-data integrity checks
  validate_metrics.R              # Numeric reproducibility checks
  validate_data_manifest.R        # Input file and checksum validation
  generate_supplementary_table2.R # Supplementary Table 2: Wilcoxon permutation results
  plotting_canonical.R            # Canonical figure module
  modules/                        # Modular stage helpers
manifests/
  data_manifest.csv               # Required input files, provenance, and checksums
  figure_manifest.csv             # Canonical figure/table artifact registry
  expected_metrics.csv            # Numeric reproducibility targets
Makefile                          # Orchestration: stage targets and smoke tests
archive/                          # Superseded and exploratory scripts (not part of pipeline)
```

## Required data

The following data files are mentioned in the code. Their paths may need to be edited depending on how the user stores the downloaded data from the sources given in the manuscript.

| File | Description |
|------|-------------|
| `pathogen_abx_analysis_all_variables_(class-specific).csv` | Combined resistance data from WHO GLASS, CAESAR, and ECDC AMR Surveillance atlas |
| `ATLAS_data/ATLAS_data_renamed.csv` | Pfizer ATLAS resistance data |
| `ATLAS_more/ATLAS_more_renamed.csv` | Additional ATLAS extracts |
| `ATLAS_Enterococcus/ATLAS_Enterococcus_renamed.csv` | ATLAS Enterococcus extracts |
| `GASP_N_renamed.csv` | GASP gonococcal surveillance data |
| `antibiotic_consumption_by_ATC3.csv` | Antibiotic consumption by ATC3 class by country and year |
| `Chungman/Chungman_pca_renamed.csv` | GDP and PCA covariate data |
| `IHME_AMR/IHME_AMR_fitted_gammas_v2.csv` | IHME AMR attributable mortality estimates |
| `IHME_AMR/IHME_AMR_PATHOGEN_2019_DATA_COUNTED_AB.CSV` | IHME pathogen-level 2019 burden data |
| `population_by_country_and_year.csv` | Population estimates by country and year |

## Installation

R version 4.2 or later is recommended. Install required packages:

```r
install.packages(c(
  "data.table", "dplyr", "forcats", "forestplot", "ggforce",
  "ggpattern", "ggplot2", "ggtext", "gridExtra", "magrittr",
  "metafor", "patchwork", "scales", "svglite", "tidyr", "tidyverse"
))
```

Confirm `Rscript` is available:

```bash
Rscript --version
```

## Reproducing the manuscript

The pipeline has two execution tracks: the **main paper** (Figures 1–4 plus all core model outputs) and the **supplementary analyses** (sensitivity models, alternative burden scenarios, and supplementary tables). These can be run independently.

### Quick validation (smoke test)

Run a fast check of each stage using reduced bootstrap iterations:

```bash
make smoke SCENARIO=main
```

### Main paper

```bash
make manuscript SCENARIO=main
```

Runs data preparation → model fitting → burden estimation → figure generation (Figures 1–4) → full output validation. Produces all four main-paper figures.

### Stage-by-stage (main paper)

```bash
make data SCENARIO=main      # Prepare and merge input data
make models SCENARIO=main    # Fit resistance-consumption models (Figures 1–2)
make burden SCENARIO=main    # Estimate avertable burden (writes intermediate CSVs)
make figures SCENARIO=main   # Generate Figures 1–4
```

### Supplementary analyses

```bash
make supplementary
```

Runs all supplementary model scenarios, alternative burden scenarios, permutation models, the Wilcoxon permutation test, and supplementary tables. Validates all supplementary outputs on completion.

Individual supplementary targets:

```bash
make models-hic              # HIC-stratified models (Supplementary Figure 1)
make models-lmic             # LMIC-stratified models (Supplementary Figure 1)
make models-nagorsen         # Within-hospital analysis (Supplementary Figure 4)
make burden-region-drug      # Burden by region × antibiotic class
make burden-region-pathogen  # Burden by region × pathogen
make permutation-models      # Permuted bootstrap models, all 7 classes (~7× model runtime)
make wilcoxon-test           # Wilcoxon test comparing bootstrapped vs permuted gradients
make tables-supp             # Supplementary Tables 2 and 3
```

Supplementary tables have a strict dependency order: `permutation-models` must complete before `wilcoxon-test`, and `wilcoxon-test` must complete before `tables-supp`.

### Validate outputs

```bash
make validate-inputs                    # Check required input files and checksums
make validate-stage-data SCENARIO=main  # Check merged data integrity
make validate-all SCENARIO=main         # Main-paper validation: outputs, metrics, manifests
make validate-supp                      # Supplementary table validation
make promote-metrics SCENARIO=main      # Promote passing pending metrics to canonical
```

### All available targets

```bash
make help
```

## Analysis scenarios

The pipeline supports multiple named scenarios. The main manuscript uses `main`.

| Scenario | Description |
|----------|-------------|
| `main` | Main manuscript analysis |
| `hic` | High-income country stratified analysis (Supplementary Figure 1) |
| `lmic` | Low- and middle-income country stratified analysis (Supplementary Figure 1) |
| `hospital_nagorsen` | Within-hospital resistance and consumption (Supplementary Figure 4) |
| `burden_optimistic` | Burden under optimistic consumption-reduction assumptions |
| `burden_pessimistic` | Burden under pessimistic consumption-reduction assumptions |
| `burden_lower_region` | Burden by lower-level IHME regions |
| `burden_upper_region` | Burden by upper-level IHME regions |
| `burden_drug_region` | Burden by drug class and region |
| `burden_pathogen_region` | Burden by pathogen and region |

To run all scenarios:

```bash
make all-scenarios
```

## Expected outputs

### Main paper (`make manuscript SCENARIO=main`)

| File | Description |
|------|-------------|
| `Figure1.pdf` | Figure 1: Pathogen- and class-specific resistance elasticity |
| `Figure2.pdf` | Figure 2: Drug random-effects elasticity by pathogen |
| `Figure3.pdf` | Figure 3: Avertable AMR-attributable mortality by pathogen |
| `Figure4.pdf` | Figure 4: Proportion of avertable burden by region, GDP, and antibiotic use |
| `Outputs/slides/` | Slide-format versions of all figures |

### Supplementary analyses (`make supplementary`)

| File | Description | Prerequisite |
|------|-------------|--------------|
| `Supplementary_Figure1_Slide_narrow.pdf` | Supplementary Figure 1: HIC/LMIC stratified elasticity | `make models-hic models-lmic` |
| `Supplementary_Figure2b_Slide_narrow.pdf` | Supplementary Figure 2: Hospital-based analysis | `make models-nagorsen` |
| `Outputs/wilcoxon_bootstrapped_vs_permuted_gradients.csv` | Wilcoxon test results (Table 2 source) | `make permutation-models wilcoxon-test` |
| `Outputs/supplementary_table_2_wilcoxon.csv` | Supplementary Table 2: permutation test p-values | `make tables-supp` |
| `formatted_variable_explained_table_for_word.txt` | Supplementary Table 3: variance explained by variable | `make tables-supp` |

## Runtime

| Track | Bottleneck | Approximate time |
|-------|-----------|-----------------|
| Main paper (`make manuscript`) | Model fitting: 1,000 bootstrap iterations per pathogen-antibiotic pair | 1-2 hours |
| Supplementary (`make supplementary`) | Permutation models: 7 classes × full model runtime | 7–20 hours additional |

Use `make models-smoke` (reduced iterations) for rapid pipeline validation only. The full `make supplementary` target is designed to run overnight or on a compute cluster.

## Numeric reproducibility

Key model outputs are checked against locked reference values in `manifests/expected_metrics.csv`. Rows with `status=canonical` are enforced by `make validate-metrics`. Rows with `status=pending` record expected values from a prior full run but are not yet enforced.

### Locking metrics after a full run

After running `make manuscript SCENARIO=main` on clean data, promote passing pending metrics to canonical:

```bash
make promote-metrics SCENARIO=main
```

This reads each pending metric row, computes the observed value from the current output file, and promotes rows that are within tolerance to `canonical`. Rows with mismatches are left as pending and reported for manual review. Repeat for supplementary scenarios after running `make supplementary`.

