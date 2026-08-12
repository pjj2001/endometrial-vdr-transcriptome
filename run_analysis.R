#!/usr/bin/env Rscript
root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
Sys.setenv(VDR_PROJECT_ROOT = root)
modules <- c(
  "01_perturbation_and_vdr_occupancy.R",
  "02_patient_bulk_projection.R",
  "03_gse250130_single_cell_analysis.R",
  "04_gse183837_single_cell_analysis.R",
  "05_cross_dataset_synthesis.R"
)
for (module in modules) {
  message("Running ", module)
  source(file.path(root, "R", module), encoding = "UTF-8", echo = FALSE)
}
