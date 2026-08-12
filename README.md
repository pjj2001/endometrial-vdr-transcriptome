# VDR-related endometrial transcriptome analysis

This repository contains the analysis and figure-generation code accompanying the manuscript **“Perturbation-defined vitamin D receptor programs undergo temporal endometrial remodeling without reproducible disruption in recurrent implantation failure.”**

## Repository structure

- `R/01_perturbation_and_vdr_occupancy.R`: perturbation-defined RNA programs and occupancy integration.
- `R/02_patient_bulk_projection.R`: public-cohort preparation and patient bulk projections.
- `R/03_gse250130_single_cell_analysis.R`: GSE250130 integration, annotation, pseudobulk and donor-level analyses.
- `R/04_gse183837_single_cell_analysis.R`: GSE183837 integration, annotation, pseudobulk and donor-level analyses.
- `R/05_cross_dataset_synthesis.R`: cross-dataset comparison and integrated summaries.
- `figures/`: scripts used to generate the main and supplementary figures.
- `environment/`: recorded R session information from representative analysis stages.
- `run_analysis.R`: sequential entry point for the five analysis modules.

## Public data

Raw data are not redistributed. The analyses use NCBI GEO accessions GSE254251, GSE306127, GSE92324, GSE165004, GSE111974, GSE250130 and GSE183837. Download the source data from GEO and arrange the inputs under the project directory according to the paths declared within the modules.

## Running the analysis

Set `VDR_PROJECT_ROOT` to the project directory, then run the modules individually or execute:

```bash
Rscript run_analysis.R
```

Large single-cell stages require substantial memory and intermediate storage. The modules preserve the staged execution order used in the study. Internal execution records and superseded development scripts are retained in the private project archive rather than published here.

## Reproducibility boundary

The public modules are a path-parameterized derivative of the final executed analysis chain. The statistical specifications and analysis order were not changed. The analysis was not rerun while preparing this repository because the reported statistical results were fixed before release packaging.

## Citation

Please cite the associated article after publication. Repository URL and archived release identifier will be added when the repository is deposited.
