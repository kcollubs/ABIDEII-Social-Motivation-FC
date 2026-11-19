# Resting-State Functional Connectivity and Social Motivation in Autism

## Project Overview

This project investigates whether individual differences in social motivation among children and adolescents with autism spectrum disorder (ASD) are associated with resting-state functional connectivity in reward-related and salience-related neural circuits.

Rather than comparing ASD and typically developing participants, this study uses a dimensional, within-ASD approach to examine variation in social motivation.

## Research Question

Do nucleus accumbens–anterior cingulate cortex (NAc–ACC) and anterior insula–ACC functional connectivity independently predict SRS Social Motivation T-scores in children and adolescents with ASD?

## Dataset

- **Dataset:** Autism Brain Imaging Data Exchange II (ABIDE II)
- **Modality:** Resting-state functional MRI (rs-fMRI)
- **Population:** Children and adolescents with ASD
- **Imaging sites:** GU, KKI, OHSU, SDSU, UCD
- **Initial imaging sample:** 171 participants
- **Final analysis sample:** 119 participants
- **Preprocessing:** DCAN derivatives
- **Parcellation:** HCP2016 FreeSurfer Subcortical atlas (379 parcels)

## Functional Connectivity Measures

Two primary ROI-to-ROI functional connectivity predictors were examined:

1. **Reward-related connectivity:** NAc–a24 (rostral ACC)
2. **Salience-related connectivity:** AVI–a24pr (anterior insula–ACC)

Functional connectivity was calculated from preprocessed resting-state time series.

## Primary Outcome

Social Responsiveness Scale (SRS) Social Motivation T-score (`SRS_MOTIVATION_T`).

## Statistical Analysis

Multiple linear regression was used to evaluate the associations between functional connectivity and social motivation while adjusting for:

- Age at scan
- Sex
- Mean retained framewise displacement (head motion)
- Imaging site (fixed effects)

The full model included both functional connectivity predictors. A covariate-only baseline model was fitted to evaluate the additional variance explained by functional connectivity.

## Primary Results

The final analysis included 119 participants.

| Predictor | Regression coefficient | p-value |
|---|---:|---:|
| NAc–a24 connectivity | 6.06 | .382 |
| AVI–a24pr connectivity | -1.05 | .864 |

**Model comparison:**

- Covariate-only R² = .1602
- Full model R² = .1666
- ΔR² = .0065
- Nested model comparison: F(2, 109) = 0.424, p = .656

Neither functional connectivity predictor was statistically significant in the primary model. Together, the connectivity predictors did not significantly improve model fit beyond the demographic, motion, and site covariates.

## Sensitivity and Exploratory Analyses

Additional analyses evaluated the robustness of the primary findings, including:

- Alternative minimum usable scan durations (≥4 and ≥5 minutes)
- Restriction to three imaging sites
- Lower-motion participant subsets
- Influential observations using Cook's distance
- Scan-duration interactions
- Alternative anterior insula/ACC parcel combinations
- False discovery rate (FDR) correction for multiple comparisons

An exploratory reward-connectivity association in the three-site, ≥5-minute subset reached nominal statistical significance (p = .032) but did not survive FDR correction (adjusted p = .124).

These findings should be interpreted as exploratory rather than confirmatory.

## Repository Structure

- `scripts/` — Data processing, functional connectivity extraction, regression, sensitivity, and exploratory analysis scripts
- `config/` — ROI specification and methodological configuration
- `results/figures/` — Research figures and regression diagnostics
- `README.md` — Project overview, methodology, and findings

## Project Status

**Completed:**
- Dataset validation and organization
- Motion quality control and analysis sample selection
- ROI definition and functional connectivity extraction
- Primary multiple regression analysis
- Covariate-only versus full-model comparison
- Sensitivity and exploratory analyses
- Results visualization

**Remaining:**
- Final reproducibility verification
- Documentation of software dependencies and execution order
- Final manuscript revisions

## Reproducibility

Analysis scripts and figures are maintained in this repository.

Raw neuroimaging data and participant-level phenotypic or analysis datasets are not included. Reproducing the analyses requires obtaining the appropriate ABIDE II data and generating the intermediate analysis files.

## Data Availability

ABIDE II is a publicly available neuroimaging dataset distributed through the Autism Brain Imaging Data Exchange.

