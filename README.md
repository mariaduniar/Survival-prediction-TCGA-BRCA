# TCGA-BRCA Survival Prediction Using Machine Learning
Machine learning-based prediction of three-year mortality risk in breast cancer using transcriptomic and somatic mutation data.
This proyect is part of a Master's thesis in Bioinformatics.

## Overview

Breast cancer is a highly heterogeneous disease, and molecular information obtained through high-throughput sequencing may provide valuable information for prognostic stratification.

This project investigates whether transcriptomic and somatic mutation data
can be used to predict three-year mortality risk in patients with breast
cancer, and whether the integration of both modalities improves predictive
performance.

Three machine learning approaches were evaluated:

- Linear-kernel Support Vector Machine (SVM)
- Radial-kernel Support Vector Machine (SVM)
- Random Forest

Models were trained using gene expression data, somatic mutation data, and
an integrated multi-omic dataset.

## Dataset

Data for this project were obtained from the TCGA-BRCA cohort through the
Genomic Data Commons using the `TCGAbiolinks` package in R.

The analysis included clinical, transcriptomic, and somatic mutation data.

Raw TCGA data are not included in this repository.

## Workflow

The analysis consisted of the following steps:

1. Data acquisition from TCGA-BRCA
2. Clinical and molecular data preprocessing
3. Definition of three-year mortality risk groups
4. Exploratory analysis of transcriptomic and mutational profiles
5. Train/test split
6. Feature selection
7. Machine learning model training
8. Evaluation and comparison of predictive models

## Repository Structure

```text
R/                  Analysis scripts

data/               Input and processed data generated through TCGAbiolinks

results/
  figures/          Figures generated during the analysis
  tables/           Results tables

docs/               Master's thesis
