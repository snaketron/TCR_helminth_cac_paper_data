# Helminth-Driven CD8⁺ T Cell Differentiation Fosters Colorectal Cancer Development

This repository contains the scripts and processed datasets used in the study
"Helminth-Driven CD8⁺ T Cell Differentiation Fosters Colorectal Cancer Development."

## Overview

The project is organized around two main components:

- `R/`: R scripts for TCR repertoire analyses
- `data/`: processed data tables and metadata used in the analyses

## Repository structure

- `R/Clonotyping.R` — clonotype-level analyses
- `R/ClustIRR.R` — clustering and repertoire similarity analyses
- `R/Diversity.R` — diversity metrics and comparisons
- `R/Gene_usage.R` — TCR gene usage analyses
- `data/mixcr_clones.RData` — processed clonotype-level repertoire data
- `data/meta.RData` — sample metadata and study annotations
- `TCR_helminth_cac_paper_data.Rproj` — RStudio project file

## Data

The processed data in `data/` include the clonotype table and sample metadata needed to reproduce the analyses described in the manuscript.

## Reproducing the analyses

1. Open the project in RStudio using `TCR_helminth_cac_paper_data.Rproj`.
2. Ensure the required R packages are installed.
3. Run the scripts in the `R/` directory as needed for your analysis workflow.
4. Use the `data/` objects as inputs to the scripts.

## Notes

This repository is intended to support reproducible analysis of T-cell receptor repertoire data in the context of the study above. For licensing details, see the `LICENSE` file.

## Citation

If you use these scripts or processed data in your work, we encourage you to cite the corresponding paper.
