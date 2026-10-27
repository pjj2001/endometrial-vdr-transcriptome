# Inputs and staged data layout

The scripts use public GEO studies and the accepted sample/annotation decisions recorded in this project. They are individual stages: an input to one stage may be an output of an earlier stage. Renaming source datasets is not necessary. Use `RUN_ORDER.csv` to locate the corresponding preparation script.

## Workspace layout

Set `VDR_PROJECT_ROOT` to a writable working folder with:

```text
03_data/
  raw_external_READ_ONLY/
    GSE254251/processed_supplementary/expression/
    GSE92324/series_matrix/
    GSE92324/supplementary/
    GSE165004/series_matrix/
    GSE111974/series_matrix/
    platform_annotations/
  processed/
  metadata_curated/
05_analysis_steps/
  STEP03_VDR_perturbation_regulons/
  STEP04_combined_temporal_patient_single_cell/
```

Preparation scripts create more detailed subfolders. They retain original relative filenames and verify stage-specific metadata. Project paths stored inside historical CSV manifests are relocated by `analysis/analysis_paths.R`; source CSVs are not rewritten. Dependencies outside the project folder are not silently substituted.

## Experimental expression and occupancy

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254251

Processed RPKM matrix:

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254251/suppl/GSE254251_RNAseq-RPKM.txt.gz

Expected matrix location: `03_data/raw_external_READ_ONLY/GSE254251/processed_supplementary/expression/GSE254251_RNAseq-RPKM.txt.gz`.

A01–A04 prepare metadata, the 12-sample column mapping, source inventory and downloaded files. A05 contains the actual factorial limma fit, not only a result-export step. A06/A06b classify candidates and prepare source tables. A07 exports the original 436-gene definition.

Occupancy uses the published author peak/annotation sources from GSE306127:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306127

A07b–A14 prepare and standardize peak annotations; A15–A17 prepare the exact/accepted alias symbol map; A18 performs the RNA–occupancy join. Author peak tables, accepted symbol mapping and overlapping attributes are not replaced with a new gene-annotation release.

## Patient bulk data

GEO studies: GSE92324, GSE165004 and GSE111974. B00–B03 prepare their source inventories, sample designs and platform annotations. B04 contains the expression-to-gene mapping, transformation, program scoring and statistical models. GSE92324 additionally uses the author's normalized probe profile:

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92324/suppl/GSE92324_sample_probe_profile.txt.gz

Platform annotations are GPL10558, GPL16699 and GPL17077 as selected in the preparation code. Do not replace accepted cohort labels, slide covariates or shared-control relationships. Selected design tables supplied under `data/analysis_metadata` can be placed at the relative destinations in `locations.csv`.

## Single-cell inputs and checkpoints

GEO studies: GSE250130 and GSE183837. B01/B02 prepare count archives and bundle maps; B03 supplies donor/sample designs. B04 prepares per-sample pre-integration objects. C01 and D02 contain single-cell count processing and doublet/integration calculations, not precomputed effect-table readers.

C02 was the continuation from the retained GSE250130 post-doublet checkpoint. The original run was staged; it must not be treated as an independently validated, uninterrupted C01-to-C02 run. Its expected checkpoint and status filenames are explicit in `integrate_gse250130.R`. A user resuming at C02 must provide that checkpoint. Large integrated objects are not redistributed in this code archive.

C05 and D05 preserve accepted broad-cell annotations; do not invent new labels from a different clustering run. C06/D06 aggregate raw counts into donor-by-cell-type pseudobulks and calculate scores. C07/D07 perform donor-level tests. C08 requires `pb_raw_counts.rds`, `02_pb_metadata.csv` and `04_model_eligibility.csv` produced by C06, not the final differential-expression tables. E01 requires C07 and D07 model outputs; it keeps the cohorts separate and performs no pooled meta-analysis.

For a working copy of the accepted staged project, the original input and checkpoint folders can be kept in place. For a fresh reconstruction from public data, run the acquisition/preparation stages and inspect their input requirements. This packaging pass did not execute such a reconstruction or certify that every historical checkpoint can be regenerated without an intervening step.

## Raw-read count sensitivity

`data/raw_read_sensitivity/` contains the compact raw gene-count matrix, 12-sample manifest, GENCODE gene mapping and original interaction comparison inputs. `workflows/` records Galaxy preprocessing; no FASTQ/BAM is duplicated here. Use the original sample order and all 12 samples. The standalone sensitivity script does not change the original program membership used in patient analyses.
