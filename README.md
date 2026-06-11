# Fire-Syndrome-Driven Shifts in European Aerofauna

[![DOI](https://img.shields.io/badge/DOI-Pending-blue.svg)](#) 
[![R](https://img.shields.io/badge/R-4.2+-blue.svg)](https://www.r-project.org/)

This repository contains the source code and analytical workflow for the research article: **"Fire-syndrome-driven shifts in community and functional diversity of European vertebrate aerofauna in the early 21st century"**.

## 📖 Project Overview

This study evaluates how different fire regimes (*Fire Syndromes* or *Pyromes*) influence the spatial distribution, species richness, and functional diversity of vertebrate aerofauna (birds and bats) in Europe during the early 21st century.

Using Species Distribution Models (SDMs) powered by **BART (Bayesian Additive Regression Trees)** algorithms, we contrast baseline models (driven solely by climate and land-use) against fire-extended models that incorporate the typology and probability of fire regimes (*Pyromes* and *Macro-Pyromes*). Finally, we project these distributions to compute Functional Diversity indices (FRic, FDis, FEve, etc.) across multidimensional trait spaces.

## 🗂️ Repository Structure

The workflow is divided into 6 sequential stages, separated by taxonomic groups (Birds and Bats). Scripts within each directory are numbered and should be executed in order.

* **`1predictorsEurope/`**: Download, processing, and alignment of baseline environmental variables (WorldClim bioclimatic data, elevation, and ESA WorldCover land-use data).
* **`2fireSyndromes/`**: Processing of historical fire perimeters and integration of the *Pyromes* and *Macro-Pyromes* classifications.
* **`3birdRecordsProcess/` & `3batRecordsProcess/`**: Taxonomic and spatial cleaning of occurrence records (eBird, GBIF, and custom datasets). This includes phenological filtering, presence rasterization, and spatial sampling bias calculation (`sampbias`).
* **`04modelingDataBirds/` & `04modelingDataBats/`**: Integration of all spatial layers (predictors, fire syndromes, sampling bias, and presences) into unified tabular data matrices (CSV) ready for modeling.
* **`05modelEvaluationBirds/` & `05modelEvaluationBats/`**: Spatial cross-validation (`blockCV`), hyperparameter optimization to prevent overfitting, predictive performance evaluation (AUC, Boyce, TSS), and relative Variable Importance calculations for all three model types (Baseline, Pyromes, Macro-Pyromes).
* **`06SDMs_projections_and_functional_diversity/` (Birds & Bats)**: Generation of binary distribution maps using optimal thresholds and parallelized computation of community functional diversity metrics using the `mFD` package alongside morphological and life-history trait databases.

## ⚙️ Prerequisites and Raw Data (Inputs)

To ensure the scripts run correctly, the following raw files must be placed in their respective directories prior to execution (these are not included in this repository due to file size constraints):

1. **Political Boundaries**: Vector file `europaPaises.gpkg`.
2. **Land Use**: Original ESA WorldCover 2021 `.tif` tiles located in `LUC_tilesEuropa/Europa/`.
3. **Fire Scars**: Historical shapefiles downloaded via *firedpy* located in `2fireSyndromes/cicatrices_2024/`.
4. **Raw Occurrences**: 
   * Birds: Massive global text file downloaded from eBird (`ebd_relJun-2025.txt`).
   * Bats: Raw GBIF file (`gbifData27102025.csv`) and custom database (`gracia.csv`).
5. **Trait Databases**: 
   * Birds: `trait_bird_Europe_data.csv`.
   * Bats: `funcTraitEuroBats.csv`.
6. **Phenology (Birds)**: Table containing breeding months and classifications (`tabFeno.csv`).

*Note: Bioclimatic variables and elevation are downloaded automatically in Stage 1 via the `geodata` package.*

## 💻 Dependencies and Software

This code was developed and tested in **R (v4.2+)**. The following main packages are required:

```R
# Spatial and Raster Data Handling
install.packages(c("terra", "sf", "geodata"))

# Modeling and Evaluation
# 'embarcadero' requires installation from GitHub: remotes::install_github("cjcarlson/embarcadero")
install.packages(c("PresenceAbsence", "ecospat", "blockCV"))
remotes::install_github("azizka/sampbias")

# Functional Diversity
install.packages(c("mFD", "ade4"))

# Data Processing and Parallelization
install.packages(c("tidyverse", "data.table", "fastDummies", "doParallel", "foreach", "correlation"))

## 📝 Authors and Citation

**Research Team:**
* Camilo Matus-Olivares
* José Ramón González-Olabarria
* Fulgencio Lisón
* María V. Jiménez-Franco
* Marcelo Miranda-Cavallieri
* Carolina Allendes-Muñoz
* Felipe Ulloa-Fierro
* Jordi Garcia-Gonzalo
* Jaime Carrasco-Barra

**If you use this code or derived data, please cite the original article:**
> *Matus-Olivares, C., González-Olabarria, J. R., Lisón, F., Jiménez-Franco, M. V., Miranda-Cavallieri, M., Allendes-Muñoz, C., Ulloa-Fierro, F., Garcia-Gonzalo, J., & Carrasco-Barra, J. (2026). Fire-syndrome-driven shifts in community and functional diversity of European vertebrate aerofauna in the early 21st century. [Journal Name: pending]. DOI: [pending]*

## 📄 License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
