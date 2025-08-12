# Code for Lysine lactylation & stress responses in canine hemangiosarcoma

This repository contains the analysis and pipelines accompanying our Science Translational Medicine submission on lysine lactylation and ATF4-mediated stress responses in canine hemangiosarcoma (HSA).

## Data availability
- **RNA-seq:** GEO **GSE304507**.
- **CUT&Tag:** GEO **GSE304509**.

> Note: The code in this repo **does not download data automatically**. Please download from GEO yourself and place files in the locations indicated in each script header.

## Pipelines & analyses
- **RNA-seq mapping/quantification:** `pipelines/rnaseq_star_rsem.sh` (STAR → RSEM; outputs counts & TPM matrices).
- **Differential expression (edgeR):** `analysis/edgeR_two_group_DE.R` (two-group QL GLM; outputs full DE table + TMM-CPM).
- **Cross-cohort comparison / PCA / markers:** `analysis/transcriptome_compare.R` (reads local matrices; generates PCA + endothelial marker plot).
- **CUT&Tag end-to-end + heatmaps/profiles:** `pipelines/CUTandTag_pipeline.sh` (produces RPGC bigWigs, gene-body and TSS heatmaps/profiles).

## How to reproduce
1. Download data from GEO (see accessions above).
2. Place files as described in each script header (e.g., `raw/`, `raw_fastq/`, `data/`, `metadata/`).
3. Run the scripts using the example commands in their headers.
4. Results will be written under `results/` (git-ignored).

## Environment
These scripts were run with commonly used versions of STAR (≥2.7), RSEM (≥1.3), samtools (≥1.10), deepTools (≥3.5), MACS2 (≥2.2), and R (≥4.2) with edgeR/limma. Consider pinning with Conda or `renv` if needed.

## License
MIT — see `LICENSE`.

## How to cite
See `CITATION.cff` or cite the archived release: https://doi.org/10.5281/zenodo.16812873

# Code for Lysine lactylation & stress responses in canine hemangiosarcoma  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16812873.svg)](https://doi.org/10.5281/zenodo.16812873) 
