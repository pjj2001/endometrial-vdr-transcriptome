# VDR-related endometrial transcriptome analysis

Code accompanying **“Vitamin D receptor related endometrial gene expression across the luteal phase and in recurrent implantation failure.”**

## Analysis scripts

The analyses are supplied as individual scripts, not concatenated modules. `RUN_ORDER.csv` lists each public filename, its original generating-script name, purpose and order. The statistical expressions, thresholds, contrasts and seeds are retained. Filenames inside generated data folders retain their historical names so that the next step can find them.

| Steps | Analysis |
|---|---|
| A01–A07 | GEO metadata, experimental design, processed RPKM factorial model and RNA-program definition |
| A07b–A19 | Published VDR peak annotations, gene-symbol mapping and occupancy integration |
| B00–B05 | Patient inputs, bulk expression and program projections |
| C01–C08 | GSE250130 processing, pseudobulk aggregation, temporal and donor-level models |
| D01–D07 | Independent GSE183837 processing and donor-level models |
| E01–E02 | Cross-cohort comparison and integrated summaries |

Some steps prepare inputs or export results rather than fit models. The order preserves the recorded staged execution, including checkpoint continuation. It is not an automatically validated single-command raw-data workflow. C02 continues from the post-doublet checkpoint; annotation steps use the accepted cell-type assignments. See `DATA_ACCESS.md` for input locations and checkpoint requirements.

## Set paths and select one step

Use a writable working copy, never the only copy of accepted study results. Set `VDR_PROJECT_ROOT` to the data workspace containing `03_data` and `05_analysis_steps`. No username, Desktop location or particular drive is required. Output paths follow that workspace. `VDR_REPORT_ROOT` optionally changes the small-report directory; otherwise reports go to `reports` within the workspace.

Example in R (replace the example path with your working-data location):

```r
Sys.setenv(VDR_PROJECT_ROOT = "/path/to/study-working-copy")
system2(file.path(R.home("bin"), "Rscript"), c("run_analysis.R", "--list"))
```

From a terminal with the environment variable set:

```text
Rscript run_analysis.R --list
Rscript run_analysis.R --check
Rscript run_analysis.R --step A05
```

`--check` parses code and reports directory availability only. Each `--step` starts one script. Choose the required step after its inputs have been prepared; do not run every script blindly on an archive of accepted results. Original administrative checks and checkpoint filenames are retained. The shared path helper also relocates project paths stored inside old CSV manifests in memory, without modifying the original files or their contents.

## Data and software

Public data accessions: GSE254251, GSE306127, GSE92324, GSE165004, GSE111974, GSE250130 and GSE183837. FASTQ, BAM, Seurat checkpoints and large pseudobulk objects are not included. Acquisition/preparation code is included; public raw data and stage-specific inputs must be prepared as described in `DATA_ACCESS.md`. Selected accepted sample designs and identifier/program definitions are included in `data/analysis_metadata/`, with their workspace destinations listed in `data/analysis_metadata/locations.csv`.

The raw-read count sensitivity analysis is separate:

```text
Rscript R/06_raw_read_sensitivity.R data/raw_read_sensitivity results/raw_read_sensitivity
```

Its compact count inputs and Galaxy workflow are included. Alignment used RNA STAR 2.7.11b on hg38/GRCh38; counting used GENCODE v43 comprehensive CHR (GRCh38.p13), reverse strandedness and paired-fragment counting. Recorded R session information and installation guidance are in `environment/`. Package versions are not silently upgraded or inferred from current defaults.

## Figures

The independent plotting scripts and their compact source tables are retained from the preceding package. See `environment/INSTALLATION.md`.

```text
Rscript tests/test_repository.R
python figures/figure1_workflow.py --check
Rscript run_figures.R 2 --check
Rscript run_figures.R 2
Rscript run_figures.R S1
```

All ten figures were generated during the earlier isolated plotting check. That checks plotting, not statistical reproduction. Figure scripts do not require an internal Skill document or the original workstation.

## Scope of this organization

This update restores the omitted computational scripts and separates the staged analysis files. It changes file-location interfaces and naming, not the statistical models. At the authors' request, analyses were not rerun during this organization. Syntax and path checks do not constitute a new numerical reproduction test. Existing manuscript files, accepted figures and statistical results were not changed.

The earlier published code version (v1.0.0) is archived at https://doi.org/10.5281/zenodo.21902833. That DOI identifies the earlier archive, not the updated main-branch code. This main-branch update has not been tagged as a new GitHub release or archived under a new Zenodo DOI. Refer to the relevant commit when using the updated code. Software license: MIT (`LICENSE`).
