# PLOS-one
# Peptide Analysis for PLOS One Publication

This repository contains R code used to analyze peptide data for a publication in PLOS One. The code is designed to process MaxQuant output files, filter peptides based on PEP scores, and prepare data for plogo visualizations and PCA analysis.

## Table of Contents

1. [Introduction](#introduction)
2. [Requirements](#requirements)
3. [Usage](#usage)
4. [Functions](#functions)
5. [Output](#output)


## Introduction

This project involves analyzing peptides from MaxQuant output files. The code filters peptides with a PEP score above 0.95, extracts n-terminal and c-terminal sequences, and prepares data for plogo analysis and PCA.

## Requirements

- R environment with necessary packages installed:
  - `readr`
  - `dplyr`
  - `stringr`
  - `fastDummies`
  - `plsgenomics`
  - `ggplot2`
  - `FactoMineR`
  - `factoextra`
  - `ggfortify`

## Usage

1.Run the R script (`PCA_peptidomics.R`) to perform the analysis.

## Functions

The code is organized into several functions:
- `load_and_filter_peptides()`: Loads and filters peptides based on PEP scores.
- `prepare_nterm_cterm_data()`: Prepares n-terminal and c-terminal data.
- `extract_substrings_and_create_plogo_input()`: Extracts substrings and creates input files for plogo analysis.
- `perform_pca()`: Performs PCA on the prepared data.

## Output

The analysis generates several output files:
- `soaked_nterminal.csv` and `soaked_cterminal.csv`: N-terminal and c-terminal sequences.
- `nterm_soaked.csv` and `cterm_soaked.csv`: Extracted n-terminal and c-terminal windows.
- `n_soaked_new.csv` and `c_soaked_new.csv`: Final data for plogo analysis.
- `nterm.csv` and `cterm.csv`: Data for PCA.
- PCA plots saved as JPEG files.

**Note**: Ensure you have the necessary input files (e.g., `peptides.txt`) in the correct directory before running the script.
