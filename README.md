# HGSVC3_sample_selection

## Overview
- R script for genomic data processing: VCF conversion, PCA, k-means clustering, and centroid distance calculations.
- Designed for analyzing 1KG, HPRC, and HGSVC sample datasets.

## Features
- **Data Import & Preprocessing**
  - Reads sample metadata from 1KG (and optionally HPRC/HGSVC).
  - Converts VCF (`pangenie_merged_bi_nosnvs.vcf.gz`) to GDS format using SNPRelate.
- **PCA Analysis**
  - Performs PCA to capture genomic variation.
  - Merges PCA results with population metadata.
  - Generates PCA plots with superpopulation-specific color coding.
- **Clustering & Centroid Analysis**
  - Applies k-means clustering per superpopulation.
  - Determines optimal cluster numbers (via gap statistic).
  - Calculates sample distances from cluster centroids.
  - Identifies extreme samples (nearest and farthest from centroids).

## Requirements
- **R Packages:**  
  - `dbscan`  
  - `fpc`  
  - `factoextra`  
  - `SNPRelate`  
  - `data.table`  
  - `raster`
- **R Version:** R ≥ 3.5 (or later)

## Input Files
- `1KG_3202_samples.ped`
- `pangenie_merged_bi_nosnvs.vcf.gz`
- Additional sample files:
  - `YEAR_1_LongRead_samples_v2.txt`
  - `YEAR_2_LongRead_samples_v2.txt`
  - `HPRC_combined_120samples.txt`
  - *(Optional)* `all_2504_1KG_samples.txt`

## Usage
- **Preparation:**  
  - Update file paths and sample subset selections as needed.
  - Install required packages.
- **Execution:**  
  - Run the script in R (e.g., `source("your_script.R")`).

## Output
- **Generated Files:**  
  - `merged_1.tab`  
  - `HGSVC_3202_Autosomes_PCA.tab`  
  - `merged_out_distances_from_centroids_21.csv`
- **Plots:**  
  - PCA scatter plots with superpopulation annotations.
  - Cluster visualizations (using `factoextra`).
