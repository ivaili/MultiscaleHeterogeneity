# MultiscaleHeterogeneity

MATLAB code accompanying the analysis of multiscale functional connectivity heterogeneity in autism. The scripts quantify case–control overlap in the spatial distribution of extreme functional connectivity deviations from a normative model, at three nested scales: individual connections, regions, and canonical functional networks.

This repository contains the code used for the overlap analyses reported in:

> Ilioska, I., *et al.* Multiscale functional connectivity heterogeneity in autism. *Nature Mental Health* (in press).
>
> Preprint: [medRxiv 2024.10.20.24315248](https://doi.org/10.1101/2024.10.20.24315248)

## Overview

Subject-level deviation (Z) scores are obtained from a Gaussian Process Regression normative model fit to functional connectivity edges, and then thresholded at |Z| > 2.3 to identify extreme positive and negative deviations. The three scripts in this repository quantify, at increasing levels of spatial aggregation, how much these extreme-deviation patterns are shared (overlap) across individuals in the autism and neurotypical groups, and test group differences via label-shuffling permutation followed by Benjamini–Hochberg FDR correction.

| Script | Scale | Unit of analysis | Test statistic |
| --- | --- | --- | --- |
| `connection_overlap.m` | Edge | 75,855 unique connections (390 × 390 upper triangle) | Group difference in proportion of subjects with a deviation at each edge |
| `region_overlap.m` | Region | 390 ROIs | Group difference in AUC of the region-level overlap curve across deviation-degree thresholds (1–20) |
| `network_overlap.m` | Network | Canonical functional networks | Group difference in network-aggregated overlap |


## Requirements

- MATLAB (tested on R2021b+)
- [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) 
- MATLAB Bioinformatics Toolbox


If you use this code, please cite the paper above.

## Contact

Iva Ilioska — [@ivaili](https://github.com/ivaili)
email:ii269@cam.ac.uk
