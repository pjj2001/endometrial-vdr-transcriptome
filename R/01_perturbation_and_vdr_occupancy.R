#!/usr/bin/env Rscript
# Public analysis module for the VDR endometrial transcriptome study.
# Set VDR_PROJECT_ROOT to the project directory containing the expected public GEO inputs.
# This module preserves the executed analysis order; raw GEO data are not redistributed.


# ==============================================================================
# Analysis stage 1: 20_step03B1_v10_lock_RNA_defer_final_figures_authorize_ATAC
# ==============================================================================
# ============================================================
# STEP03B1 v10 - Lock RNA results, defer final figure release,
# and authorize STEP03B2 ATAC analysis
# Nutrition_Endometrial_VDR_v2
#
# v10 correction:
# - validates the interaction summary as a class-by-class count table;
# - calculates total candidates by sum(n_genes);
# - calculates attenuation candidates by summing classes beginning with
#   'attenuated_ligand_';
# - reconciles summary counts against the 436-row candidate table.
#
# Candidate-first policy:
# - analytical result lock is separate from final figure release;
# - current figures remain B-class candidates;
# - cosmetic/composite polish is deferred;
# - source package must be complete;
# - ATAC may proceed without C-class figure release.
#
# This script does NOT:
# - refit limma;
# - change samples, thresholds, P values or FDR values;
# - claim the current QC PNGs are final submission figures;
# - generate PDF/TIFF/SVG;
# - assign direct VDR-target status.
#
# Run:
# source(
#   "04_code/R/20_step03B1_v10_lock_RNA_defer_final_figures_authorize_ATAC.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
module_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

run_id <- paste0(
  "STEP03B1_v10_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

freeze_version <- "RNA_PROGRAM_v1.0"
candidate_status <- "B_LOCKED_DEFERRED_POLISH"

dirs <- list(
  control = file.path(module_root, "00_module_control"),
  code = file.path(module_root, "02_code"),
  parameters = file.path(module_root, "03_parameters"),
  environment = file.path(module_root, "04_environment"),
  execution = file.path(module_root, "05_execution"),
  results = file.path(module_root, "06_results"),
  figures = file.path(module_root, "07_figures_tables"),
  qc = file.path(module_root, "08_qc"),
  text = file.path(module_root, "09_manuscript_text"),
  archive_parent = file.path(
    module_root,
    "10_archive_work"
  )
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_root <- file.path(
  dirs$figures,
  "B_submission_candidate",
  "STEP03B1_RNA_candidates_v3"
)

if (!dir.exists(candidate_root)) {
  stop(
    "Candidate package not found:\n",
    candidate_root
  )
}

if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages(
    "digest",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages(
    "jsonlite",
    repos = "https://cloud.r-project.org"
  )
}

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

md5_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

copy_validate <- function(source_path, destination_path) {
  dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)

  if (nchar(destination_path) >= 240L) {
    stop(
      "Destination path exceeds the 240-character hard stop:\n",
      destination_path
    )
  }

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(file.info(source_path)$size),
      as.numeric(file.info(destination_path)$size)
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed:\n",
      source_path
    )
  }

  TRUE
}

log_path <- file.path(
  dirs$execution,
  paste0(run_id, "_execution.log")
)

log_event <- function(...) {
  msg <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )
  cat(msg, "\n")
  cat(msg, "\n", file = log_path, append = TRUE)
}

log_event("Starting", run_id)

# ----------------------------------------------------------------------
# Validate analytical inputs
# ----------------------------------------------------------------------

interaction_path <- file.path(
  dirs$results,
  "STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv"
)

interaction_summary_path <- file.path(
  dirs$results,
  "STEP03B1_v4_interaction_program_summary.csv"
)

all_programs_path <- file.path(
  dirs$results,
  "STEP03B1_program_gene_candidates_NOT_FROZEN.csv"
)

definition_path <- file.path(
  dirs$results,
  "STEP03B1_v4_RNA_program_definition_card.csv"
)

required_analysis_files <- c(
  interaction_path,
  interaction_summary_path,
  all_programs_path,
  definition_path
)

missing_analysis <- required_analysis_files[
  !file.exists(required_analysis_files)
]

if (length(missing_analysis) > 0L) {
  stop(
    "Missing analysis-lock inputs:\n",
    paste(missing_analysis, collapse = "\n")
  )
}

interaction <- read.csv(
  interaction_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

interaction_summary <- read.csv(
  interaction_summary_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_programs <- read.csv(
  all_programs_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

definition <- read.csv(
  definition_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(interaction) != 436L) {
  stop(
    "Expected 436 interaction candidates, found ",
    nrow(interaction)
  )
}

required_summary_columns <- c(
  "program_class",
  "n_genes"
)

missing_summary_columns <- setdiff(
  required_summary_columns,
  names(interaction_summary)
)

if (length(missing_summary_columns) > 0L) {
  stop(
    "Interaction summary is missing columns: ",
    paste(
      missing_summary_columns,
      collapse = ", "
    )
  )
}

interaction_summary$n_genes <- suppressWarnings(
  as.integer(
    as.character(
      interaction_summary$n_genes
    )
  )
)

if (
  any(!is.finite(interaction_summary$n_genes)) ||
    any(interaction_summary$n_genes < 0)
) {
  stop(
    "Interaction summary contains invalid n_genes values."
  )
}

# The summary is a class-by-class table. The reviewed totals therefore
# must be calculated by summing rows, not searched for as literal cells.
summary_total <- sum(
  interaction_summary$n_genes
)

summary_attenuated <- sum(
  interaction_summary$n_genes[
    grepl(
      "^attenuated_ligand_",
      interaction_summary$program_class
    )
  ]
)

candidate_total <- nrow(interaction)

if (!"program_class" %in% names(interaction)) {
  stop(
    "Interaction candidate table is missing program_class."
  )
}

candidate_attenuated <- sum(
  grepl(
    "^attenuated_ligand_",
    interaction$program_class
  )
)

if (
  candidate_total != 436L ||
    candidate_attenuated != 431L
) {
  stop(
    "Candidate-table totals differ from the reviewed lock: total=",
    candidate_total,
    "; attenuated=",
    candidate_attenuated,
    "."
  )
}

if (
  summary_total != candidate_total ||
    summary_attenuated != candidate_attenuated
) {
  stop(
    "Interaction summary does not reconcile with the candidate table: ",
    "summary total=",
    summary_total,
    ", candidate total=",
    candidate_total,
    "; summary attenuated=",
    summary_attenuated,
    ", candidate attenuated=",
    candidate_attenuated,
    "."
  )
}

interaction_count_reconciliation <- data.frame(
  metric = c(
    "candidate_total",
    "summary_total",
    "candidate_attenuated",
    "summary_attenuated"
  ),
  value = c(
    candidate_total,
    summary_total,
    candidate_attenuated,
    summary_attenuated
  ),
  expected = c(
    436L,
    436L,
    431L,
    431L
  ),
  pass = c(
    candidate_total == 436L,
    summary_total == 436L,
    candidate_attenuated == 431L,
    summary_attenuated == 431L
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  interaction_count_reconciliation,
  file.path(
    dirs$results,
    "STEP03B1_v10_INTERACTION_COUNT_RECONCILIATION.csv"
  )
)

# ----------------------------------------------------------------------
# Validate the complete Class B source package
# ----------------------------------------------------------------------

required_candidate_files <- c(
  file.path(
    candidate_root,
    "figures",
    "Figure_STEP03B1_interaction_volcano_candidate_v3_QC.png"
  ),
  file.path(
    candidate_root,
    "figures",
    "Figure_STEP03B1_PCA_candidate_v3_QC.png"
  ),
  file.path(
    candidate_root,
    "figures",
    "Figure_STEP03B1_top_interaction_heatmap_ACCEPTED_v2_QC.png"
  ),
  file.path(
    candidate_root,
    "figures",
    "Figure_STEP03B1_VDR_expression_candidate_v3_QC.png"
  ),
  file.path(
    candidate_root,
    "manifests",
    "STEP03B1_v3_FIGURE_SOURCE_MAP.csv"
  ),
  file.path(
    candidate_root,
    "source_data",
    "README_SOURCE_DATA.md"
  ),
  file.path(
    candidate_root,
    "LANGUAGE_ROUTING_CARD.md"
  ),
  file.path(
    candidate_root,
    "legends",
    "STEP03B1_v3_candidate_legends.md"
  )
)

missing_candidate <- required_candidate_files[
  !file.exists(required_candidate_files)
]

if (length(missing_candidate) > 0L) {
  stop(
    "Class B candidate package is incomplete:\n",
    paste(missing_candidate, collapse = "\n")
  )
}

direct_sources <- list.files(
  file.path(candidate_root, "source_data"),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = FALSE
)

upstream_refs <- list.files(
  file.path(
    candidate_root,
    "source_data",
    "upstream_refs"
  ),
  pattern = "\\.csv$",
  full.names = TRUE
)

generator_code <- list.files(
  file.path(candidate_root, "code"),
  pattern = "\\.R$",
  full.names = TRUE
)

if (
  length(direct_sources) < 6L ||
    length(upstream_refs) < 4L ||
    length(generator_code) < 1L
) {
  stop(
    "The self-contained source package does not meet the Class B minimum."
  )
}

# ----------------------------------------------------------------------
# Freeze RNA statistical membership
# ----------------------------------------------------------------------

freeze_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

add_freeze <- function(x, scope) {
  cbind(
    data.frame(
      program_freeze_version =
        rep(freeze_version, nrow(x)),
      program_freeze_run =
        rep(run_id, nrow(x)),
      program_freeze_date =
        rep(freeze_date, nrow(x)),
      membership_status =
        rep("FROZEN_AT_RNA_STATISTICAL_LAYER", nrow(x)),
      program_scope =
        rep(scope, nrow(x)),
      direct_target_tier =
        rep("PENDING_ATAC_AND_CUTRUN", nrow(x)),
      stringsAsFactors = FALSE
    ),
    x
  )
}

interaction_frozen <- add_freeze(
  interaction,
  paste(
    "Formal siRNA-by-ligand interaction;",
    "FDR < 0.05 and absolute interaction log2FC >= 0.5"
  )
)

all_programs_frozen <- add_freeze(
  all_programs,
  "Five prespecified GSE254251 RNA contrasts"
)

definition$program_freeze_version <- freeze_version
definition$program_freeze_run <- run_id
definition$program_freeze_date <- freeze_date
definition$membership_status <-
  "FROZEN_AT_RNA_STATISTICAL_LAYER"
definition$direct_target_tier <-
  "PENDING_ATAC_AND_CUTRUN"

write_csv_utf8(
  interaction_frozen,
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"
  )
)

write_csv_utf8(
  all_programs_frozen,
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_all_contrast_programs_FROZEN.csv"
  )
)

write_csv_utf8(
  definition,
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_program_definition_FROZEN.csv"
  )
)

# ----------------------------------------------------------------------
# Candidate lock with deferred final polish
# ----------------------------------------------------------------------

figure_paths <- required_candidate_files[1:4]

figure_lock <- data.frame(
  figure_id = c(
    "STEP03B1-RNA-VOLCANO",
    "STEP03B1-RNA-PCA",
    "STEP03B1-RNA-HEATMAP",
    "STEP03B1-RNA-VDR"
  ),
  figure_file = basename(figure_paths),
  current_class = "B_submission_candidate",
  candidate_lock_status = candidate_status,
  scientific_qc = "PASS",
  statistical_qc = "PASS",
  candidate_readability_qc =
    "PASS_FOR_RESULT_VERIFICATION",
  final_visual_qc = "DEFERRED",
  final_release_status = "NOT_STARTED",
  proposed_role = c(
    "MAIN_CANDIDATE",
    "SUPPLEMENTARY_CANDIDATE",
    "MAIN_CANDIDATE",
    "SUPPORTING_MECHANISTIC_CANDIDATE"
  ),
  deferred_defect = c(
    paste(
      "Rebalance final labels and typography",
      "after composite width is known."
    ),
    paste(
      "Reposition centroid labels after",
      "supplementary-panel geometry is known."
    ),
    paste(
      "Recheck font size and color bar",
      "inside the final composite."
    ),
    paste(
      "Refine point/summary geometry",
      "at final physical panel width."
    )
  ),
  md5 = vapply(
    figure_paths,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    figure_paths,
    sha256_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  figure_lock,
  file.path(
    dirs$results,
    "STEP03B1_v10_B_CANDIDATE_LOCK_DECISION.csv"
  )
)

write_text(
  c(
    "# STEP03B1 deferred final figure release",
    "",
    "- Current class: B_submission_candidate",
    paste0("- Candidate lock status: ", candidate_status),
    "- Scientific QC: PASS",
    "- Statistical QC: PASS",
    "- Candidate readability: sufficient to verify the result",
    "- Final publication visual QC: DEFERRED",
    "- PDF/TIFF/SVG release: NOT_STARTED",
    "- Main/supplementary final numbering: DEFERRED",
    "- Composite geometry: DEFERRED",
    "",
    "The current QC PNGs are not final submission figures.",
    "All source tables, upstream references, exact code, legends, QC cards",
    "and manifests are retained for manuscript-wide release without",
    "rerunning the RNA analysis."
  ),
  file.path(
    candidate_root,
    "qc",
    "STEP03B1_v9_DEFERRED_FINAL_RELEASE_CARD.md"
  )
)

write_text(
  c(
    "# STEP03B1 RNA program lock card",
    "",
    paste0("- Version: ", freeze_version),
    paste0("- Run: ", run_id),
    "- Samples: 12/12 retained",
    paste0("- Interaction candidates: ", candidate_total),
    paste0("- Attenuation-pattern candidates: ", candidate_attenuated),
    "- RNA statistical membership: FROZEN",
    "- Direct-target tier: PENDING ATAC/CUT&RUN",
    "- Figure source package: LOCKED",
    paste0("- Candidate visual status: ", candidate_status),
    "- Final submission artwork: DEFERRED",
    "",
    "ATAC and CUT&RUN may add support tiers but may not silently change",
    "RNA membership."
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"
  )
)

gate <- data.frame(
  gate_component = c(
    "RNA scientific/statistical QC",
    "RNA interaction membership",
    "Class B source package",
    "Candidate readability",
    "Final publication visual QC",
    "C-class release",
    "Direct-target tier",
    "STEP03B1 v10 gate"
  ),
  status = c(
    "PASS",
    "FROZEN_RNA_PROGRAM_v1.0",
    "PASS_COMPLETE",
    "PASS_FOR_RESULT_VERIFICATION",
    "DEFERRED_TO_MANUSCRIPT_RELEASE",
    "NOT_STARTED",
    "PENDING_ATAC_AND_CUTRUN",
    "READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED"
  ),
  required_action = c(
    "Retain accepted analysis outputs.",
    "Do not silently change membership.",
    "Use source map for later figure release.",
    "No more standalone candidate polishing now.",
    "Perform once after manuscript figure plan is stable.",
    "Generate PDF/TIFF/PNG/SVG only in RELEASE_FIGURE.",
    "Assign only after cross-omics integration.",
    "Begin STEP03B2 ATAC."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate,
  file.path(
    dirs$results,
    "STEP03B1_v10_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP03B2 ATAC authorization",
    "",
    paste0("- Authorized by: ", run_id),
    "- Authorization: YES",
    "- RNA anchor: RNA_PROGRAM_v1.0",
    "- RNA membership: IMMUTABLE",
    "- Figure C-class release: NOT REQUIRED FOR ATAC START",
    "- ATAC role: chromatin-accessibility support and mechanistic tiering",
    "- Prohibited: count-based differential-accessibility claims from",
    "  normalized BigWig and condition-level merged BED alone",
    "- Required: self-contained source package for every retained ATAC figure"
  ),
  file.path(
    dirs$control,
    "STEP03B2_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP03B1 transition",
    "",
    "RNA scientific/statistical results and program membership are locked.",
    "The four current plots are reproducible B-class candidates with final",
    "visual polish deferred to manuscript-wide figure release.",
    "Proceed to STEP03B2 ATAC without rerunning RNA analysis or continuing",
    "standalone cosmetic iterations."
  ),
  file.path(
    dirs$text,
    "STEP03B1_v10_TRANSITION_TO_STEP03B2.md"
  )
)

# ----------------------------------------------------------------------
# Short-path archive
# ----------------------------------------------------------------------

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S3B1v10_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

for (d in c(
  file.path(archive_root, "B3"),
  file.path(archive_root, "res"),
  file.path(archive_root, "ctl"),
  file.path(archive_root, "txt"),
  file.path(archive_root, "code"),
  file.path(archive_root, "env"),
  file.path(archive_root, "log")
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_files <- list.files(
  candidate_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

copy_rows <- vector("list", length(candidate_files))

for (i in seq_along(candidate_files)) {
  source_path <- candidate_files[i]
  relative_path <- substring(
    normalize_slashes(source_path),
    nchar(normalize_slashes(candidate_root)) + 2L
  )
  destination_path <- file.path(
    archive_root,
    "B3",
    relative_path
  )

  copy_validate(source_path, destination_path)

  copy_rows[[i]] <- data.frame(
    relative_path = relative_path,
    source_md5 = md5_file(source_path),
    archive_md5 = md5_file(destination_path),
    source_sha256 = sha256_file(source_path),
    archive_sha256 = sha256_file(destination_path),
    validation_pass = TRUE,
    stringsAsFactors = FALSE
  )
}

copy_manifest <- do.call(rbind, copy_rows)

write_csv_utf8(
  copy_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B1_v10_B3_COPY_MANIFEST.csv"
  )
)


same_normalized_path <- function(path_a, path_b) {
  if (
    is.null(path_a) ||
      is.null(path_b) ||
      length(path_a) != 1L ||
      length(path_b) != 1L ||
      is.na(path_a) ||
      is.na(path_b)
  ) {
    return(FALSE)
  }

  normalized_a <- tryCatch(
    normalizePath(
      path_a,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
  )
}

safe_file_copy <- function(
    from,
    to,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
) {
  if (same_normalized_path(from, to)) {
    return(TRUE)
  }

  base::file.copy(
    from = from,
    to = to,
    overwrite = overwrite,
    recursive = recursive,
    copy.mode = copy.mode,
    copy.date = copy.date
  )
}

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(e) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  project_code <- file.path(
    project_root,
    "04_code",
    "R",
    "20_step03B1_v10_lock_RNA_defer_final_figures_authorize_ATAC.R"
  )

  module_code <- file.path(
    dirs$code,
    "20_step03B1_v10_lock_RNA_defer_final_figures_authorize_ATAC.R"
  )

  safe_file_copy(
    source_file,
    project_code,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    module_code,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

parameters <- list(
  schema_version = "1.0",
  run_id = run_id,
  policy =
    "candidate-first source-locked release-deferred",
  RNA_program_version = freeze_version,
  RNA_membership_frozen = TRUE,
  figure_class = "B_submission_candidate",
  candidate_lock_status = candidate_status,
  final_figure_release = "DEFERRED",
  ATAC_authorized = TRUE
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP03B1_v10_parameters.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(sessionInfo()),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B1_v10.txt"
  )
)

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP03B1_v10_INTERACTION_COUNT_RECONCILIATION.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_all_contrast_programs_FROZEN.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_program_definition_FROZEN.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_B_CANDIDATE_LOCK_DECISION.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"
  ),
  file.path(
    dirs$results,
    "STEP03B1_v10_GATE_TABLE.csv"
  ),
  file.path(
    dirs$text,
    "STEP03B1_v10_TRANSITION_TO_STEP03B2.md"
  ),
  file.path(
    dirs$control,
    "STEP03B2_AUTHORIZATION.md"
  ),
  file.path(
    dirs$parameters,
    "STEP03B1_v10_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B1_v10.txt"
  ),
  log_path
)

for (p in archive_artifacts[file.exists(archive_artifacts)]) {
  section <- if (grepl("\\.csv$|LOCK_CARD", basename(p))) {
    "res"
  } else if (grepl("TRANSITION|AUTHORIZATION", basename(p))) {
    "txt"
  } else if (grepl("sessionInfo|parameters", basename(p))) {
    "env"
  } else {
    "log"
  }

  copy_validate(
    p,
    file.path(
      archive_root,
      section,
      basename(p)
    )
  )
}

archive_files <- list.files(
  archive_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

archive_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(archive_files),
    nchar(normalize_slashes(archive_root)) + 2L
  ),
  size_bytes = as.numeric(file.info(archive_files)$size),
  md5 = vapply(
    archive_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    archive_files,
    sha256_file,
    character(1)
  ),
  path_length = nchar(archive_files),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B1_v10_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP03B1 v10 archive status",
    "",
    paste0("- Run ID: ", run_id),
    "- RNA result: LOCKED",
    paste0("- Candidate figure status: ", candidate_status),
    "- Final publication visual QC: DEFERRED",
    "- C-class release: NOT_STARTED",
    paste0(
      "- Candidate package files copied: ",
      nrow(copy_manifest)
    ),
    paste0(
      "- Maximum archive path length: ",
      max(archive_manifest$path_length)
    ),
    "- STEP03B2 ATAC: AUTHORIZED"
  ),
  file.path(
    archive_root,
    "ctl",
    "ARCHIVE_STATUS.md"
  )
)

write_text(
  c(
    "# STEP03 Module Status",
    "",
    "- Completed substep: STEP03B1",
    paste0("- Accepted run: ", run_id),
    "- RNA program: FROZEN RNA_PROGRAM_v1.0",
    paste0("- RNA figures: ", candidate_status),
    "- Final figure release: DEFERRED",
    "- Active next substep: STEP03B2 ATAC",
    "- Full STEP03 module: NOT_LOCKED"
  ),
  file.path(
    dirs$control,
    "STATUS.md"
  )
)

log_event(
  "Completed",
  run_id,
  "Gate: READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED"
)

cat("\n============================================================\n")
cat("STEP03B1 v10 candidate-first closeout completed.\n")
cat("RNA program frozen: RNA_PROGRAM_v1.0\n")
cat("Interaction candidates: ", candidate_total, "\n")
cat("Attenuation-pattern candidates: ", candidate_attenuated, "\n")
cat("Current figure class: B_submission_candidate\n")
cat("Candidate lock: B_LOCKED_DEFERRED_POLISH\n")
cat("Final submission figure release: DEFERRED\n")
cat("Gate: READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED\n")
cat("Authoritative archive:\n  ", archive_root, "\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 2: 35_step03B3B_v3_lock_published_VDR_cistrome_source
# ==============================================================================
# ============================================================
# STEP03B3B v3 - Lock published VDR cistrome source
# Nutrition_Endometrial_VDR_v2
#
# This closeout:
# - validates the STEP03B3B v2 chromosome-style repair;
# - confirms exact reconciliation of all four published table checks;
# - classifies the 16 non-coordinate media-6 rows as appendix/blank rows,
#   not failed peaks;
# - locks three primary VDR cistrome coordinate sources;
# - locks two media-2 tables as derived author overlap references only;
# - retains chrM records in the source tables but excludes them from the
#   primary nuclear-genome integration universe;
# - defers local peak calling because the author-processed sources are
#   internally exact and sufficient for primary integration;
# - authorizes STEP03B3C peak-to-gene schema audit and project-specific
#   RNA/CUT&RUN integration preparation.
#
# This script does NOT:
# - change peak membership or coordinates;
# - call peaks;
# - change RNA_PROGRAM_v1.0;
# - assign definitive direct VDR targets;
# - perform RNA/CUT&RUN integration;
# - generate figures.
#
# Run from the project root:
# source(
#   "04_code/R/35_step03B3B_v3_lock_published_VDR_cistrome_source.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

module_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

run_id <- paste0(
  "STEP03B3B_v3_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"
RNA_anchor <- "RNA_PROGRAM_v1.0"
ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

required_packages <- c(
  "digest",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package(s) missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ". This script does not install or update packages."
  )
}

dirs <- list(
  control = file.path(module_root, "00_module_control"),
  code = file.path(module_root, "02_code"),
  parameters = file.path(module_root, "03_parameters"),
  environment = file.path(module_root, "04_environment"),
  execution = file.path(module_root, "05_execution"),
  results = file.path(module_root, "06_results"),
  qc = file.path(module_root, "08_qc"),
  text = file.path(module_root, "09_manuscript_text"),
  archive_parent = file.path(module_root, "10_archive_work")
)

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S3B3Bv3_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

for (directory_path in c(
  file.path(archive_root, "res"),
  file.path(archive_root, "qc"),
  file.path(archive_root, "txt"),
  file.path(archive_root, "code"),
  file.path(archive_root, "env"),
  file.path(archive_root, "log"),
  file.path(archive_root, "ctl")
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  dirs$execution,
  paste0(
    run_id,
    "_execution.log"
  )
)

log_event <- function(...) {
  message_text <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(message_text, "\n")

  cat(
    message_text,
    "\n",
    file = log_path,
    append = TRUE
  )
}

write_text <- function(lines, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
}

write_csv_utf8 <- function(x, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

normalize_slashes <- function(x) {
  gsub(
    "\\\\",
    "/",
    x
  )
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

copy_validate <- function(
    source_path,
    destination_path
) {
  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (nchar(destination_path) >= 240L) {
    stop(
      "Archive destination exceeds 240 characters:\n",
      destination_path
    )
  }

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(
        file.info(source_path)$size
      ),
      as.numeric(
        file.info(destination_path)$size
      )
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed:\n",
      source_path
    )
  }

  TRUE
}

is_nonempty_value <- function(x) {
  !is.na(x) &
    nzchar(
      trimws(
        as.character(x)
      )
    )
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Required v2 outputs
# ----------------------------------------------------------------------

required_paths <- c(
  gate_v2 = file.path(
    dirs$results,
    "STEP03B3B_v2_GATE_TABLE.csv"
  ),
  source_candidate_v2 = file.path(
    dirs$results,
    "STEP03B3B_v2_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"
  ),
  coordinate_summary_v2 = file.path(
    dirs$qc,
    "STEP03B3B_v2_COORDINATE_SCHEMA_SUMMARY.csv"
  ),
  chromosome_style_v2 = file.path(
    dirs$qc,
    "STEP03B3B_v2_CHROMOSOME_STYLE_SUMMARY.csv"
  ),
  invalid_rows_v2 = file.path(
    dirs$qc,
    "STEP03B3B_v2_INVALID_COORDINATE_ROWS.csv"
  ),
  reconciliation_v2 = file.path(
    dirs$qc,
    "STEP03B3B_v2_PUBLISHED_TABLE_RECONCILIATION.csv"
  ),
  RNA_lock_card = file.path(
    dirs$results,
    "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"
  ),
  ATAC_role_policy = file.path(
    dirs$results,
    "STEP03B2A2_ATAC_EVIDENCE_ROLE_POLICY.csv"
  )
)

missing_required_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_required_paths) > 0L) {
  stop(
    "Required source-lock input is missing:\n",
    paste(
      missing_required_paths,
      collapse = "\n"
    )
  )
}

gate_v2 <- read.csv(
  required_paths[["gate_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_candidate_v2 <- read.csv(
  required_paths[["source_candidate_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

coordinate_summary_v2 <- read.csv(
  required_paths[["coordinate_summary_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

chromosome_style_v2 <- read.csv(
  required_paths[["chromosome_style_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

invalid_rows_v2 <- read.csv(
  required_paths[["invalid_rows_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reconciliation_v2 <- read.csv(
  required_paths[["reconciliation_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

RNA_lock_text <- readLines(
  required_paths[["RNA_lock_card"]],
  warn = FALSE,
  encoding = "UTF-8"
)

ATAC_role_policy <- read.csv(
  required_paths[["ATAC_role_policy"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

accepted_gate <- gate_v2$status[
  gate_v2$gate_component ==
    "STEP03B3B v2 gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "PROCEED_TO_MANUAL_PUBLISHED_PEAK_SOURCE_LOCK"
) {
  stop(
    "STEP03B3B v2 does not authorize source lock. Gate: ",
    paste(
      accepted_gate,
      collapse = ", "
    )
  )
}

if (
  !any(
    grepl(
      RNA_anchor,
      RNA_lock_text,
      fixed = TRUE
    )
  ) ||
    !any(
      grepl(
        "RNA statistical membership: FROZEN",
        RNA_lock_text,
        fixed = TRUE
      )
    )
) {
  stop(
    "RNA_PROGRAM_v1.0 is not confirmed as frozen."
  )
}

if (!all(
  ATAC_role_policy$status ==
    ATAC_role
)) {
  stop(
    "ATAC role policy is not uniformly restricted descriptive support."
  )
}

# ----------------------------------------------------------------------
# Reconciliation and source-table identity checks
# ----------------------------------------------------------------------

expected_audits <- c(
  "SHARED_TO_VEHICLE",
  "SHARED_TO_CALCITRIOL",
  "VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME",
  "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME"
)

if (
  nrow(reconciliation_v2) != 4L ||
    !setequal(
      reconciliation_v2$audit_id,
      expected_audits
    ) ||
    !all(
      reconciliation_v2$any_overlap_fraction ==
        1
    ) ||
    !all(
      reconciliation_v2$exact_match_fraction ==
        1
    )
) {
  stop(
    "Published-table reconciliation is not exact for all four locked audits."
  )
}

if (
  nrow(source_candidate_v2) != 5L ||
    !all(
      source_candidate_v2$source_lock_status ==
        "CANDIDATE_FOR_MANUAL_SOURCE_LOCK"
    ) ||
    !all(
      source_candidate_v2$chromosome_style ==
        "canonical_UCSC_chr"
    )
) {
  stop(
    "The five source-lock candidates are incomplete."
  )
}

source_candidate_v2$file_exists <- file.exists(
  source_candidate_v2$clean_table_path
)

source_candidate_v2$observed_md5 <- vapply(
  source_candidate_v2$clean_table_path,
  md5_file,
  character(1)
)

source_candidate_v2$observed_sha256 <- vapply(
  source_candidate_v2$clean_table_path,
  sha256_file,
  character(1)
)

source_candidate_v2$identity_pass <- with(
  source_candidate_v2,
  file_exists &
    clean_table_md5 ==
      observed_md5 &
    clean_table_sha256 ==
      observed_sha256
)

if (!all(
  source_candidate_v2$identity_pass
)) {
  stop(
    "One or more canonical clean source tables failed hash revalidation."
  )
}

# ----------------------------------------------------------------------
# Classify the 16 non-coordinate rows
# ----------------------------------------------------------------------

if (
  nrow(invalid_rows_v2) != 16L ||
    !all(
      invalid_rows_v2$source_id_from_list ==
        "VDR_CALCITRIOL_CISTROME"
    ) ||
    !all(
      !invalid_rows_v2$valid_coordinate
    )
) {
  stop(
    "The invalid-row audit is not the expected 16-row calcitriol appendix."
  )
}

coordinate_fields_empty <- (
  !is_nonempty_value(
    invalid_rows_v2$chromosome_source
  )
) &
  is.na(
    invalid_rows_v2$start0
  ) &
  is.na(
    invalid_rows_v2$end
  )

if (!all(
  coordinate_fields_empty
)) {
  stop(
    "At least one invalid row contains a partial genomic coordinate and requires manual review."
  )
}

has_summary_label <- is_nonempty_value(
  invalid_rows_v2$raw_geneOverlap
)

has_summary_value <- if (
  "raw_promoterOverlap" %in%
    names(invalid_rows_v2)
) {
  is_nonempty_value(
    invalid_rows_v2$raw_promoterOverlap
  )
} else {
  rep(
    FALSE,
    nrow(invalid_rows_v2)
  )
}

invalid_rows_v2$manual_row_class <- ifelse(
  has_summary_label |
    has_summary_value,
  "NON_COORDINATE_APPENDIX_SUMMARY_ROW",
  "BLANK_SEPARATOR_ROW"
)

invalid_rows_v2$peak_membership_status <-
  "NOT_A_PEAK_ROW"

invalid_rows_v2$coordinate_exclusion_reason <- ifelse(
  invalid_rows_v2$manual_row_class ==
    "NON_COORDINATE_APPENDIX_SUMMARY_ROW",
  "AUTHOR_APPENDED_ANNOTATION_SUMMARY",
  "AUTHOR_APPENDED_BLANK_SEPARATOR"
)

invalid_rows_v2$manual_decision <-
  "RETAIN_IN_AUDIT_EXCLUDE_FROM_COORDINATE_SOURCE"

write_csv_utf8(
  invalid_rows_v2,
  file.path(
    dirs$qc,
    "STEP03B3B_v3_NON_COORDINATE_APPENDIX_ROW_DECISIONS.csv"
  )
)

appendix_summary <- as.data.frame(
  table(
    invalid_rows_v2$manual_row_class
  ),
  stringsAsFactors = FALSE
)

names(appendix_summary) <- c(
  "manual_row_class",
  "n_rows"
)

appendix_summary$source_id <-
  "VDR_CALCITRIOL_CISTROME"

appendix_summary$peak_rows_removed <-
  0L

appendix_summary$decision <-
  "NON_PEAK_APPENDIX_ROWS_EXCLUDED_FROM_COORDINATE_SOURCE"

write_csv_utf8(
  appendix_summary,
  file.path(
    dirs$qc,
    "STEP03B3B_v3_NON_COORDINATE_APPENDIX_SUMMARY.csv"
  )
)

if (
  sum(
    invalid_rows_v2$manual_row_class ==
      "NON_COORDINATE_APPENDIX_SUMMARY_ROW"
  ) != 13L ||
    sum(
      invalid_rows_v2$manual_row_class ==
        "BLANK_SEPARATOR_ROW"
    ) != 3L
) {
  stop(
    "The 16-row appendix did not resolve to 13 summary rows and 3 blank separators."
  )
}

# ----------------------------------------------------------------------
# Lock exact source counts and evidence roles
# ----------------------------------------------------------------------

expected_counts <- data.frame(
  source_id = c(
    "VDR_VEHICLE_CISTROME",
    "VDR_CALCITRIOL_CISTROME",
    "VDR_SHARED_CISTROME",
    "VEHICLE_CISTROME_RNA_OVERLAP",
    "CALCITRIOL_CISTROME_RNA_OVERLAP"
  ),
  total_source_rows = c(
    885L,
    6108L,
    519L,
    134L,
    2022L
  ),
  valid_coordinate_rows = c(
    885L,
    6092L,
    519L,
    134L,
    2022L
  ),
  primary_nuclear_rows = c(
    884L,
    6091L,
    518L,
    134L,
    2022L
  ),
  nonprimary_standard_rows = c(
    1L,
    1L,
    1L,
    0L,
    0L
  ),
  stringsAsFactors = FALSE
)

summary_source_map <- c(
  MEDIA6_pLentiVeh =
    "VDR_VEHICLE_CISTROME",
  MEDIA6_pLenti125 =
    "VDR_CALCITRIOL_CISTROME",
  MEDIA6_pLentiVeh_125overlap =
    "VDR_SHARED_CISTROME",
  MEDIA2_OverlapPeaksVeh_ReceptorDEG =
    "VEHICLE_CISTROME_RNA_OVERLAP",
  MEDIA2_OverlapPeaks125_125DEG =
    "CALCITRIOL_CISTROME_RNA_OVERLAP"
)

observed_counts <- coordinate_summary_v2

observed_counts$source_id_locked <-
  unname(
    summary_source_map[
      observed_counts$source_id
    ]
  )

observed_counts$nonprimary_standard_rows <-
  observed_counts$valid_coordinate_rows -
  observed_counts$primary_coordinate_rows

count_check <- merge(
  expected_counts,
  observed_counts[
    c(
      "source_id_locked",
      "total_rows",
      "valid_coordinate_rows",
      "primary_coordinate_rows",
      "nonprimary_standard_rows"
    )
  ],
  by.x = "source_id",
  by.y = "source_id_locked",
  all = TRUE,
  sort = FALSE,
  suffixes = c(
    "_expected",
    "_observed"
  )
)

count_check$count_pass <- with(
  count_check,
  total_source_rows ==
    total_rows &
    valid_coordinate_rows_expected ==
      valid_coordinate_rows_observed &
    primary_nuclear_rows ==
      primary_coordinate_rows &
    nonprimary_standard_rows_expected ==
      nonprimary_standard_rows_observed
)

write_csv_utf8(
  count_check,
  file.path(
    dirs$qc,
    "STEP03B3B_v3_SOURCE_COUNT_LOCK_CHECK.csv"
  )
)

if (!all(
  count_check$count_pass
)) {
  stop(
    "One or more published source counts differ from the reviewed lock."
  )
}

source_lock <- source_candidate_v2

source_lock$source_lock_version <-
  source_lock_version

source_lock$source_class <- ifelse(
  source_lock$source_id %in%
    c(
      "VDR_VEHICLE_CISTROME",
      "VDR_CALCITRIOL_CISTROME",
      "VDR_SHARED_CISTROME"
    ),
  "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES",
  "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_REFERENCE"
)

source_lock$allowed_project_role <- ifelse(
  source_lock$source_class ==
    "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES",
  "PROJECT_SPECIFIC_PEAK_TO_GENE_MAPPING_AND_RNA_OCCUPANCY_SUPPORT",
  "RECONCILIATION_AND_AUTHOR_RESULT_CROSS_CHECK_ONLY"
)

source_lock$prohibited_project_role <- ifelse(
  source_lock$source_class ==
    "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES",
  "CANNOT_ALONE_ESTABLISH_CAUSAL_OR_DEFINITIVE_DIRECT_TARGET_STATUS",
  "CANNOT_REPLACE_RNA_PROGRAM_v1.0_OR_SERVE_AS_INDEPENDENT_VALIDATION"
)

source_lock$nuclear_integration_policy <-
  "USE_chr1_TO_chr22_chrX_chrY;RETAIN_chrM_SEPARATELY"

source_lock$local_peak_calling_status <-
  "DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS"

source_lock$source_lock_status <-
  "LOCKED"

write_csv_utf8(
  source_lock,
  file.path(
    dirs$results,
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"
  )
)

source_role_summary <- data.frame(
  source_class = c(
    "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES",
    "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_REFERENCE"
  ),
  n_sources = c(
    3L,
    2L
  ),
  project_authority = c(
    "PRIMARY_CHROMATIN_OCCUPANCY_SOURCE",
    "DERIVED_CROSS_CHECK_ONLY"
  ),
  independent_validation = c(
    FALSE,
    FALSE
  ),
  changes_RNA_membership = c(
    FALSE,
    FALSE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  source_role_summary,
  file.path(
    dirs$results,
    "STEP03B3B_v3_SOURCE_ROLE_SUMMARY.csv"
  )
)

# ----------------------------------------------------------------------
# Gate and authorization
# ----------------------------------------------------------------------

gate_table <- data.frame(
  gate_component = c(
    "STEP03B3B v2 gate",
    "Canonical source identities",
    "Published-table reconciliation",
    "Chromosome style",
    "Non-coordinate calcitriol rows",
    "Primary VDR cistrome sources",
    "Derived author overlap sources",
    "chrM handling",
    "Local peak calling",
    "RNA membership",
    "ATAC evidence role",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP03B3C authorization",
    "STEP03B3B v3 gate"
  ),
  status = c(
    "PASS",
    "PASS_5_OF_5",
    "PASS_4_OF_4_EXACT",
    "LOCKED_CANONICAL_UCSC_chr",
    "13_SUMMARY_PLUS_3_BLANK_NOT_PEAK_ROWS",
    "LOCKED_3_OF_3",
    "LOCKED_2_OF_2_DERIVED_ONLY",
    "RETAINED_SEPARATELY_EXCLUDED_FROM_PRIMARY_NUCLEAR_JOIN",
    "DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS",
    paste0(
      "IMMUTABLE_",
      RNA_anchor
    ),
    ATAC_role,
    "NOT_PERFORMED",
    "NONE_CREATED",
    "AUTHORIZED",
    "READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP"
  ),
  required_action = c(
    "Retain the v2 repair record.",
    "Use only hash-locked canonical tables.",
    "Preserve the exact-match reconciliation record.",
    "Use chr-prefixed names in all later genomic joins.",
    "Retain these rows in audit only; do not count them as failed peaks.",
    "Use media-6 primary sources for project-specific occupancy support.",
    "Use media-2 tables only as author-derived cross-checks.",
    "Do not silently mix mitochondrial and nuclear integration results.",
    "Do not install or run MACS3 for the primary analysis.",
    "CUT&RUN may add evidence tiers but cannot change RNA membership.",
    "Use ATAC only as optional accessible context.",
    "Use the term direct-target candidate only after project-specific integration.",
    "No Scientific Figure Skill 3 workflow is triggered.",
    "Audit peak annotation fields before creating gene-level evidence tiers.",
    "Proceed to STEP03B3C."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP03B3B_v3_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# Published VDR cistrome source lock",
    "",
    paste0(
      "- Source-lock version: ",
      source_lock_version
    ),
    paste0(
      "- Run ID: ",
      run_id
    ),
    "- Genome build: hg38 as reported by the source article",
    "- Chromosome style: canonical UCSC chr",
    "- Primary source 1: media-6 pLentiVeh, 885 valid coordinate rows",
    "- Primary source 2: media-6 pLenti125, 6,092 valid coordinate rows",
    "- Primary source 3: media-6 pLentiVeh_125overlap, 519 valid coordinate rows",
    "- Nuclear primary rows: 884 vehicle, 6,091 calcitriol, 518 shared",
    "- chrM records: retained separately, one in each primary source",
    "- Derived media-2 overlap sources: 134 vehicle and 2,022 calcitriol rows",
    "- Reconciliation: four of four audits are 100% exact matches",
    "- Non-coordinate appendix: 13 summary rows plus 3 blank separator rows",
    "- Peak rows removed by appendix classification: 0",
    "- Local peak calling: deferred; not required for primary analysis",
    paste0(
      "- RNA anchor: ",
      RNA_anchor,
      " IMMUTABLE"
    ),
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- Definitive direct-target assignment: not performed",
    "- Scientific figures: none",
    "- Next: STEP03B3C peak-to-gene schema audit and integration preparation"
  ),
  file.path(
    dirs$results,
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK_CARD.md"
  )
)

write_text(
  c(
    "# STEP03B3C authorization",
    "",
    paste0(
      "- Authorized by: ",
      run_id
    ),
    paste0(
      "- Published source lock: ",
      source_lock_version
    ),
    "- Authorization: YES",
    "- Authorized task 1: audit author-provided peak annotation fields",
    "- Authorized task 2: build a transparent peak-to-gene mapping table",
    "- Authorized task 3: prepare project-specific RNA/CUT&RUN evidence tiers",
    "- RNA membership changes: PROHIBITED",
    "- Definitive causal or direct-target claims: PROHIBITED",
    "- media-2 author overlap tables: cross-check only",
    "- ATAC: optional accessible context only",
    "- Local MACS peak calling: not required for primary analysis",
    "- Scientific figures: create only after gene-level evidence tables pass QC"
  ),
  file.path(
    dirs$control,
    "STEP03B3C_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP03B3B v3 Methods addendum",
    "",
    "Following chromosome-style harmonization, all four published-table reconciliation tests showed complete exact coordinate matching. The author-processed vehicle, calcitriol and shared VDR cistrome coordinate tables were therefore locked as the primary chromatin-occupancy sources for subsequent project-specific peak-to-gene mapping.",
    "",
    "Sixteen trailing rows in the calcitriol worksheet contained no genomic coordinates. Manual inspection classified 13 as appended annotation-summary rows and three as blank separator rows; none represented failed peak coordinates. These rows were retained in the audit record but excluded from the coordinate source.",
    "",
    "One mitochondrial interval was retained in each primary cistrome table. Mitochondrial intervals remain in the source archive but will be reported separately from the primary nuclear-genome integration."
  ),
  file.path(
    dirs$text,
    "STEP03B3B_v3_Methods_addendum.md"
  )
)

write_text(
  c(
    "# STEP03B3B v3 Results draft",
    "",
    "Canonical chromosome harmonization resolved the apparent failure of the shared-peak reconciliation. All 518 primary-chromosome shared regions were exact matches in both the vehicle and calcitriol cistrome tables. The two author-derived cistrome-transcriptome overlap tables likewise showed complete exact reconciliation with their corresponding primary cistrome sources (134 of 134 vehicle regions and 2,022 of 2,022 calcitriol regions).",
    "",
    "The primary author-processed sources contained 885 valid vehicle regions, 6,092 valid calcitriol regions and 519 valid shared regions. Sixteen additional rows in the calcitriol worksheet were non-coordinate appendix content rather than peak records."
  ),
  file.path(
    dirs$text,
    "STEP03B3B_v3_Results_draft.md"
  )
)

parameters <- list(
  schema_version = "3.0",
  run_id = run_id,
  source_lock_version =
    source_lock_version,
  RNA_anchor = RNA_anchor,
  ATAC_role = ATAC_role,
  primary_sources = 3,
  derived_sources = 2,
  exact_reconciliation_audits = 4,
  non_coordinate_summary_rows = 13,
  blank_separator_rows = 3,
  peak_rows_removed = 0,
  local_peak_calling_required_for_primary_analysis =
    FALSE,
  definitive_direct_target_assignment =
    FALSE,
  scientific_figures_created =
    FALSE,
  STEP03B3C_authorized =
    TRUE,
  gate =
    "READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP"
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP03B3B_v3_parameters.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3B_v3.txt"
  )
)


same_normalized_path <- function(path_a, path_b) {
  if (
    is.null(path_a) ||
      is.null(path_b) ||
      length(path_a) != 1L ||
      length(path_b) != 1L ||
      is.na(path_a) ||
      is.na(path_b)
  ) {
    return(FALSE)
  }

  normalized_a <- tryCatch(
    normalizePath(
      path_a,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
  )
}

safe_file_copy <- function(
    from,
    to,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
) {
  if (same_normalized_path(from, to)) {
    return(TRUE)
  }

  base::file.copy(
    from = from,
    to = to,
    overwrite = overwrite,
    recursive = recursive,
    copy.mode = copy.mode,
    copy.date = copy.date
  )
}

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  safe_file_copy(
    source_file,
    file.path(
      project_root,
      "04_code",
      "R",
      "35_step03B3B_v3_lock_published_VDR_cistrome_source.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    file.path(
      dirs$code,
      "35_step03B3B_v3_lock_published_VDR_cistrome_source.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

# ----------------------------------------------------------------------
# Short-path archive
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3B_v3_SOURCE_ROLE_SUMMARY.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3B_v3_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK_CARD.md"
  ),
  file.path(
    dirs$qc,
    "STEP03B3B_v3_NON_COORDINATE_APPENDIX_ROW_DECISIONS.csv"
  ),
  file.path(
    dirs$qc,
    "STEP03B3B_v3_NON_COORDINATE_APPENDIX_SUMMARY.csv"
  ),
  file.path(
    dirs$qc,
    "STEP03B3B_v3_SOURCE_COUNT_LOCK_CHECK.csv"
  ),
  file.path(
    dirs$control,
    "STEP03B3C_AUTHORIZATION.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3B_v3_Methods_addendum.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3B_v3_Results_draft.md"
  ),
  file.path(
    dirs$parameters,
    "STEP03B3B_v3_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3B_v3.txt"
  ),
  file.path(
    dirs$code,
    "35_step03B3B_v3_lock_published_VDR_cistrome_source.R"
  ),
  log_path
)

for (path_value in archive_artifacts[
  file.exists(
    archive_artifacts
  )
]) {
  section <- if (
    grepl(
      "\\.R$",
      basename(path_value)
    )
  ) {
    "code"
  } else if (
    grepl(
      "sessionInfo|parameters",
      basename(path_value)
    )
  ) {
    "env"
  } else if (
    grepl(
      "\\.log$",
      basename(path_value)
    )
  ) {
    "log"
  } else if (
    grepl(
      "Methods|Results|AUTHORIZATION",
      basename(path_value)
    )
  ) {
    "txt"
  } else if (
    grepl(
      "SUMMARY|DECISION|CHECK|LOCK_CARD",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "qc"
  } else {
    "res"
  }

  copy_validate(
    path_value,
    file.path(
      archive_root,
      section,
      basename(path_value)
    )
  )
}

archive_files <- list.files(
  archive_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

archive_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      archive_files
    ),
    nchar(
      normalize_slashes(
        archive_root
      )
    ) + 2L
  ),
  size_bytes = as.numeric(
    file.info(
      archive_files
    )$size
  ),
  md5 = vapply(
    archive_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    archive_files,
    sha256_file,
    character(1)
  ),
  path_length = nchar(
    archive_files
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B3B_v3_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP03B3B v3 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Source lock: ",
      source_lock_version
    ),
    "- Primary VDR cistrome sources: 3",
    "- Derived author overlap references: 2",
    "- Exact reconciliation audits: 4/4",
    "- Non-coordinate appendix rows: 16",
    "- Peak rows removed: 0",
    "- Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS",
    "- STEP03B3C authorization: YES",
    paste0(
      "- Maximum archive path length: ",
      max(
        archive_manifest$path_length
      )
    )
  ),
  file.path(
    archive_root,
    "ctl",
    "ARCHIVE_STATUS.md"
  )
)

write_text(
  c(
    "# STEP03 Module Status",
    "",
    "- Completed substep: STEP03B3B",
    paste0(
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- Published VDR cistrome source: ",
      source_lock_version,
      " LOCKED"
    ),
    paste0(
      "- RNA anchor: ",
      RNA_anchor,
      " IMMUTABLE"
    ),
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS",
    "- Definitive direct-target assignment: NOT_PERFORMED",
    "- Active next substep: STEP03B3C",
    "- Scientific figures: NONE",
    "- Full STEP03 module: NOT_LOCKED"
  ),
  file.path(
    dirs$control,
    "STATUS.md"
  )
)

log_event(
  "Completed",
  run_id,
  "Gate: READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP"
)

cat("\n============================================================\n")
cat("STEP03B3B v3 published VDR cistrome source lock completed.\n")
cat(
  "Source lock: ",
  source_lock_version,
  "\n"
)
cat("Primary VDR cistrome sources: 3/3 LOCKED\n")
cat("Derived author overlap references: 2/2 LOCKED AS DERIVED\n")
cat("Published reconciliation: 4/4 EXACT\n")
cat("Non-coordinate calcitriol rows: 13 summary + 3 blank\n")
cat("Peak rows removed: 0\n")
cat("Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS\n")
cat(
  "RNA anchor: ",
  RNA_anchor,
  " IMMUTABLE\n"
)
cat(
  "ATAC role: ",
  ATAC_role,
  "\n"
)
cat("Definitive direct-target assignment: NOT PERFORMED\n")
cat("Scientific figures: NONE\n")
cat("STEP03B3C authorization: YES\n")
cat("Gate: READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP\n")
cat("Authoritative archive:\n  ", archive_root, "\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 3: 40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1
# ==============================================================================
# ============================================================
# STEP03B3C2 v3 - Lock PEAK_GENE_MAP_v1.1
# Nutrition_Endometrial_VDR_v2
#
# Scientific purpose:
# - validate the repaired canonical-source peak-gene candidate map;
# - verify exact source-specific relation counts, unique relation IDs,
#   canonical source IDs, occupancy contexts and explicit zero-row
#   exception schemas;
# - copy the accepted candidate artifacts into a locked versioned source;
# - freeze field semantics and mapping membership for downstream RNA
#   integration;
# - authorize STEP03B3C3 RNA gene-ID harmonization and
#   occupancy-evidence integration.
#
# This closeout does NOT:
# - modify PUBLISHED_VDR_CISTROME_SOURCE_v1.0;
# - modify RNA_PROGRAM_v1.0 or its 436 interaction candidates;
# - add, remove or reinterpret peak-gene relations;
# - use media-2 derived tables for discovery;
# - perform RNA/CUT&RUN integration;
# - assign definitive direct VDR targets;
# - run local peak calling;
# - generate scientific figures.
#
# Run from project root:
# source(
#   "04_code/R/40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

module_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

run_id <- paste0(
  "STEP03B3C2_v3_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"
RNA_anchor <- "RNA_PROGRAM_v1.0"
ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"
peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

expected_relation_counts <- c(
  VDR_VEHICLE_CISTROME = 2848L,
  VDR_CALCITRIOL_CISTROME = 19037L,
  VDR_SHARED_CISTROME = 1614L
)

expected_total_relations <- 23499L
expected_unique_gene_symbols <- 10042L
expected_primary_source_n <- 3L
expected_shared_peak_n <- 518L

required_packages <- c(
  "data.table",
  "digest",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package(s) missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ". This script does not install or update packages."
  )
}

dirs <- list(
  control = file.path(module_root, "00_module_control"),
  inputs = file.path(module_root, "01_inputs"),
  code = file.path(module_root, "02_code"),
  parameters = file.path(module_root, "03_parameters"),
  environment = file.path(module_root, "04_environment"),
  execution = file.path(module_root, "05_execution"),
  results = file.path(module_root, "06_results"),
  figures = file.path(module_root, "07_figures_tables"),
  qc = file.path(module_root, "08_qc"),
  text = file.path(module_root, "09_manuscript_text"),
  archive_parent = file.path(module_root, "10_archive_work")
)

required_dir_names <- c(
  "control",
  "inputs",
  "code",
  "parameters",
  "environment",
  "execution",
  "results",
  "figures",
  "qc",
  "text",
  "archive_parent"
)

missing_dir_names <- setdiff(
  required_dir_names,
  names(dirs)
)

invalid_dir_names <- names(dirs)[
  !vapply(
    dirs,
    function(path_value) {
      is.character(path_value) &&
        length(path_value) == 1L &&
        !is.na(path_value) &&
        nzchar(path_value)
    },
    logical(1)
  )
]

if (
  length(missing_dir_names) > 0L ||
    length(invalid_dir_names) > 0L
) {
  stop(
    "Invalid directory registry. Missing: ",
    paste(
      missing_dir_names,
      collapse = ", "
    ),
    "; invalid: ",
    paste(
      invalid_dir_names,
      collapse = ", "
    )
  )
}

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

candidate_root <- file.path(
  dirs$results,
  "STEP03B3C2_v2_peak_gene_map_candidate"
)

locked_root <- file.path(
  dirs$results,
  "STEP03B3C2_locked_peak_gene_map",
  peak_gene_map_lock_version
)

dir.create(
  locked_root,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S3B3C2v3_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

for (directory_path in c(
  file.path(archive_root, "res"),
  file.path(archive_root, "qc"),
  file.path(archive_root, "txt"),
  file.path(archive_root, "code"),
  file.path(archive_root, "env"),
  file.path(archive_root, "log"),
  file.path(archive_root, "ctl")
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  dirs$execution,
  paste0(
    run_id,
    "_execution.log"
  )
)

log_event <- function(...) {
  message_text <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(message_text, "\n")

  cat(
    message_text,
    "\n",
    file = log_path,
    append = TRUE
  )
}

write_text <- function(lines, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
}

write_csv_utf8 <- function(x, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

copy_validate <- function(
    source_path,
    destination_path
) {
  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (nchar(destination_path) >= 240L) {
    stop(
      "Archive destination exceeds 240 characters:\n",
      destination_path
    )
  }

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(
        file.info(source_path)$size
      ),
      as.numeric(
        file.info(destination_path)$size
      )
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed:\n",
      source_path
    )
  }

  TRUE
}

parse_bool <- function(x) {
  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1"
  )
}

read_empty_schema_csv <- function(
    path,
    expected_columns
) {
  table_value <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (
    nrow(table_value) != 0L ||
      !identical(
        names(table_value),
        expected_columns
      )
  ) {
    stop(
      "Zero-row exception table failed schema validation: ",
      path
    )
  }

  table_value
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Required repaired candidate-map outputs
# ----------------------------------------------------------------------

required_paths <- c(
  gate_v2 = file.path(
    dirs$results,
    "STEP03B3C2_v2_GATE_TABLE.csv"
  ),
  transition_lock = file.path(
    dirs$results,
    "STEP03B3C2_v2_SOURCE_ID_TRANSITION_LOCK.csv"
  ),
  field_semantics = file.path(
    dirs$results,
    "STEP03B3C2_FIELD_SEMANTICS_LOCK.csv"
  ),
  relation_reconciliation = file.path(
    dirs$qc,
    "STEP03B3C2_v2_RELATION_COUNT_RECONCILIATION.csv"
  ),
  mapping_qc = file.path(
    dirs$qc,
    "STEP03B3C2_v2_MAPPING_QC_SUMMARY.csv"
  ),
  promoter_exceptions = file.path(
    dirs$qc,
    "STEP03B3C2_v2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"
  ),
  invalid_tokens = file.path(
    dirs$qc,
    "STEP03B3C2_v2_INVALID_GENE_TOKEN_ROWS.csv"
  ),
  candidate_manifest = file.path(
    candidate_root,
    "STEP03B3C2_v2_CANDIDATE_MAP_MANIFEST.csv"
  ),
  candidate_long_map = file.path(
    candidate_root,
    "STEP03B3C2_v2_PEAK_GENE_LONG_CANDIDATE.tsv.gz"
  ),
  candidate_gene_preview = file.path(
    candidate_root,
    "STEP03B3C2_v2_GENE_LEVEL_EVIDENCE_PREVIEW.csv"
  ),
  source_lock = file.path(
    dirs$results,
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"
  ),
  RNA_lock_card = file.path(
    dirs$results,
    "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"
  )
)

missing_required_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_required_paths) > 0L) {
  stop(
    "Required STEP03B3C2 v3 input is missing:\n",
    paste(
      missing_required_paths,
      collapse = "\n"
    )
  )
}

gate_v2 <- read.csv(
  required_paths[["gate_v2"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

transition_lock <- read.csv(
  required_paths[["transition_lock"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

field_semantics <- read.csv(
  required_paths[["field_semantics"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

relation_reconciliation <- read.csv(
  required_paths[["relation_reconciliation"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping_qc <- read.csv(
  required_paths[["mapping_qc"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_manifest <- read.csv(
  required_paths[["candidate_manifest"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_lock <- read.csv(
  required_paths[["source_lock"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

RNA_lock_text <- readLines(
  required_paths[["RNA_lock_card"]],
  warn = FALSE,
  encoding = "UTF-8"
)

accepted_gate <- gate_v2$status[
  gate_v2$gate_component ==
    "STEP03B3C2 v2 gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "PROCEED_TO_MANUAL_PEAK_GENE_MAP_v1_1_LOCK"
) {
  stop(
    "STEP03B3C2 v2 does not authorize map lock. Gate: ",
    paste(
      accepted_gate,
      collapse = ", "
    )
  )
}

if (
  nrow(field_semantics) != 6L ||
    !all(
      field_semantics$status ==
        "LOCKED_FOR_CANDIDATE_MAP"
    )
) {
  stop(
    "The six annotation-field semantics are not fully locked."
  )
}

if (
  nrow(transition_lock) != 3L ||
    !setequal(
      transition_lock$canonical_source_id,
      names(expected_relation_counts)
    ) ||
    any(
      parse_bool(
        transition_lock$relation_membership_changed
      )
    )
) {
  stop(
    "Source-ID transition lock failed."
  )
}

primary_source_lock <- source_lock[
  source_lock$source_class ==
    "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES",
  ,
  drop = FALSE
]

if (
  nrow(primary_source_lock) !=
    expected_primary_source_n ||
    !all(
      primary_source_lock$source_lock_version ==
        source_lock_version
    ) ||
    !all(
      primary_source_lock$source_lock_status ==
        "LOCKED"
    )
) {
  stop(
    "Published VDR cistrome source lock is incomplete."
  )
}

if (
  !any(
    grepl(
      RNA_anchor,
      RNA_lock_text,
      fixed = TRUE
    )
  ) ||
    !any(
      grepl(
        "RNA statistical membership: FROZEN",
        RNA_lock_text,
        fixed = TRUE
      )
    )
) {
  stop(
    "RNA_PROGRAM_v1.0 is not confirmed as frozen."
  )
}

# ----------------------------------------------------------------------
# Validate explicit exception schemas
# ----------------------------------------------------------------------

promoter_exception_columns <- c(
  "source_id",
  "occupancy_context",
  "coordinate_key",
  "promoter_only_gene",
  "raw_gene_overlap_value",
  "raw_promoter_overlap_value"
)

invalid_token_columns <- c(
  "source_id",
  "occupancy_context",
  "coordinate_key",
  "invalid_gene_token"
)

promoter_exceptions <- read_empty_schema_csv(
  required_paths[["promoter_exceptions"]],
  promoter_exception_columns
)

invalid_tokens <- read_empty_schema_csv(
  required_paths[["invalid_tokens"]],
  invalid_token_columns
)

# ----------------------------------------------------------------------
# Validate relation-count reconciliation and mapping QC
# ----------------------------------------------------------------------

if (
  nrow(relation_reconciliation) != 3L ||
    !setequal(
      relation_reconciliation$source_id,
      names(expected_relation_counts)
    ) ||
    !all(
      relation_reconciliation$relation_count_pass
    )
) {
  stop(
    "Source-specific relation reconciliation is incomplete."
  )
}

relation_reconciliation$expected_locked <- unname(
  expected_relation_counts[
    relation_reconciliation$source_id
  ]
)

if (
  !all(
    relation_reconciliation$expected_relation_n ==
      relation_reconciliation$expected_locked
  ) ||
    !all(
      relation_reconciliation$observed_relation_n ==
        relation_reconciliation$expected_locked
    )
) {
  stop(
    "Source-specific relation counts differ from the reviewed lock."
  )
}

qc_value <- setNames(
  mapping_qc$value,
  mapping_qc$metric
)

required_qc_values <- c(
  "Vehicle peak-gene relations" = 2848,
  "Calcitriol peak-gene relations" = 19037,
  "Shared inherited peak-gene relations" = 1614,
  "Total peak-gene relations" = 23499,
  "Unique author gene symbols" = 10042,
  "Promoter-only relation exceptions" = 0,
  "Invalid gene-token rows" = 0,
  "Shared gene-set concordance" = 1,
  "Shared promoter-set concordance" = 1,
  "Shared relative-location concordance" = 1,
  "Shared TSS-distance concordance" = 1,
  "Canonical source IDs represented" = 3
)

missing_qc_metrics <- setdiff(
  names(required_qc_values),
  names(qc_value)
)

if (length(missing_qc_metrics) > 0L) {
  stop(
    "Required QC metrics are missing: ",
    paste(
      missing_qc_metrics,
      collapse = ", "
    )
  )
}

observed_qc_values <- as.numeric(
  qc_value[
    names(required_qc_values)
  ]
)

if (!all(
  observed_qc_values ==
    as.numeric(required_qc_values)
)) {
  stop(
    "Mapping QC values differ from the reviewed lock."
  )
}

if (
  as.numeric(
    qc_value[["Shared nuclear peaks"]]
  ) != expected_shared_peak_n
) {
  stop(
    "Shared nuclear peak count differs from 518."
  )
}

# ----------------------------------------------------------------------
# Validate candidate-manifest identities
# ----------------------------------------------------------------------

if (
  nrow(candidate_manifest) != 2L ||
    !setequal(
      candidate_manifest$artifact_id,
      c(
        "PEAK_GENE_LONG_CANDIDATE",
        "GENE_LEVEL_EVIDENCE_PREVIEW"
      )
    ) ||
    !all(
      candidate_manifest$version ==
        "PEAK_GENE_MAP_CANDIDATE_v1.1"
    ) ||
    !all(
      candidate_manifest$source_id_status ==
        "CANONICAL_LOCKED_SOURCE_IDS"
    ) ||
    any(
      parse_bool(
        candidate_manifest$relation_membership_changed
      )
    )
) {
  stop(
    "Candidate-map manifest metadata failed validation."
  )
}

candidate_manifest$file_exists <- file.exists(
  candidate_manifest$path
)

candidate_manifest$observed_size_bytes <- as.numeric(
  file.info(
    candidate_manifest$path
  )$size
)

candidate_manifest$observed_md5 <- vapply(
  candidate_manifest$path,
  md5_file,
  character(1)
)

candidate_manifest$observed_sha256 <- vapply(
  candidate_manifest$path,
  sha256_file,
  character(1)
)

candidate_manifest$identity_pass <- with(
  candidate_manifest,
  file_exists &
    size_bytes ==
      observed_size_bytes &
    md5 ==
      observed_md5 &
    sha256 ==
      observed_sha256
)

if (!all(
  candidate_manifest$identity_pass
)) {
  stop(
    "One or more candidate-map artifacts failed hash revalidation."
  )
}

write_csv_utf8(
  candidate_manifest,
  file.path(
    dirs$qc,
    "STEP03B3C2_v3_CANDIDATE_ARTIFACT_IDENTITY_RECHECK.csv"
  )
)

# ----------------------------------------------------------------------
# Deep candidate-long-map validation
# ----------------------------------------------------------------------

peak_gene_map <- data.table::fread(
  required_paths[["candidate_long_map"]],
  sep = "\t",
  data.table = FALSE,
  showProgress = FALSE
)

required_map_columns <- c(
  "source_id",
  "source_id_upstream",
  "occupancy_context",
  "coordinate_key",
  "gene_symbol_author",
  "peak_gene_relation_id"
)

missing_map_columns <- setdiff(
  required_map_columns,
  names(peak_gene_map)
)

if (length(missing_map_columns) > 0L) {
  stop(
    "Peak-gene map is missing columns: ",
    paste(
      missing_map_columns,
      collapse = ", "
    )
  )
}

if (nrow(peak_gene_map) != expected_total_relations) {
  stop(
    "Peak-gene map row count differs from 23,499."
  )
}

if (
  any(
    is.na(
      peak_gene_map$source_id
    ) |
      !nzchar(
        trimws(
          peak_gene_map$source_id
        )
      )
  ) ||
    any(
      is.na(
        peak_gene_map$coordinate_key
      ) |
      !nzchar(
        trimws(
          peak_gene_map$coordinate_key
        )
      )
    ) ||
    any(
      is.na(
        peak_gene_map$gene_symbol_author
      ) |
      !nzchar(
        trimws(
          peak_gene_map$gene_symbol_author
        )
      )
    )
) {
  stop(
    "Peak-gene map contains missing required identifiers."
  )
}

if (
  length(
    unique(
      peak_gene_map$gene_symbol_author
    )
  ) != expected_unique_gene_symbols
) {
  stop(
    "Unique author gene-symbol count differs from 10,042."
  )
}

if (
  anyDuplicated(
    peak_gene_map$peak_gene_relation_id
  ) > 0L
) {
  stop(
    "Peak-gene relation IDs are not unique."
  )
}

expected_relation_id <- paste0(
  peak_gene_map$source_id,
  "::",
  peak_gene_map$coordinate_key,
  "::",
  peak_gene_map$gene_symbol_author
)

if (!identical(
  as.character(
    peak_gene_map$peak_gene_relation_id
  ),
  as.character(
    expected_relation_id
  )
)) {
  stop(
    "At least one peak-gene relation ID is inconsistent with its fields."
  )
}

observed_source_counts <- table(
  peak_gene_map$source_id
)

if (
  !setequal(
    names(observed_source_counts),
    names(expected_relation_counts)
  ) ||
    !all(
      as.integer(
        observed_source_counts[
          names(expected_relation_counts)
        ]
      ) ==
        as.integer(
          expected_relation_counts
        )
    )
) {
  stop(
    "Deep source-specific relation counts failed."
  )
}

context_by_source <- c(
  VDR_VEHICLE_CISTROME =
    "vehicle",
  VDR_CALCITRIOL_CISTROME =
    "calcitriol",
  VDR_SHARED_CISTROME =
    "shared_vehicle_and_calcitriol"
)

if (!all(
  peak_gene_map$occupancy_context ==
    unname(
      context_by_source[
        peak_gene_map$source_id
      ]
    )
)) {
  stop(
    "Occupancy contexts are inconsistent with canonical source IDs."
  )
}

upstream_by_source <- c(
  VDR_VEHICLE_CISTROME =
    "MEDIA6_pLentiVeh",
  VDR_CALCITRIOL_CISTROME =
    "MEDIA6_pLenti125",
  VDR_SHARED_CISTROME =
    "VDR_SHARED_CISTROME"
)

if (!all(
  peak_gene_map$source_id_upstream ==
    unname(
      upstream_by_source[
        peak_gene_map$source_id
      ]
    )
)) {
  stop(
    "Upstream source provenance is inconsistent."
  )
}

gene_preview <- read.csv(
  required_paths[["candidate_gene_preview"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_preview_columns <- c(
  "source_id",
  "occupancy_context",
  "gene_symbol_author",
  "supporting_peak_n",
  "gene_overlap_supported",
  "promoter_overlap_supported"
)

missing_preview_columns <- setdiff(
  required_preview_columns,
  names(gene_preview)
)

if (length(missing_preview_columns) > 0L) {
  stop(
    "Gene-level preview is missing columns: ",
    paste(
      missing_preview_columns,
      collapse = ", "
    )
  )
}

preview_key <- paste(
  gene_preview$source_id,
  gene_preview$occupancy_context,
  gene_preview$gene_symbol_author,
  sep = "::"
)

if (
  anyDuplicated(preview_key) > 0L ||
    any(
      gene_preview$supporting_peak_n < 1L
    )
) {
  stop(
    "Gene-level evidence preview failed uniqueness or supporting-peak checks."
  )
}

# ----------------------------------------------------------------------
# Create immutable locked copies
# ----------------------------------------------------------------------

locked_long_map_path <- file.path(
  locked_root,
  "PEAK_GENE_MAP_v1.1_peak_gene_long.tsv.gz"
)

locked_gene_preview_path <- file.path(
  locked_root,
  "PEAK_GENE_MAP_v1.1_gene_level_evidence_preview.csv"
)

copy_validate(
  required_paths[["candidate_long_map"]],
  locked_long_map_path
)

copy_validate(
  required_paths[["candidate_gene_preview"]],
  locked_gene_preview_path
)

locked_manifest <- data.frame(
  artifact_id = c(
    "PEAK_GENE_LONG",
    "GENE_LEVEL_EVIDENCE_PREVIEW"
  ),
  lock_version =
    peak_gene_map_lock_version,
  path = normalize_slashes(
    c(
      locked_long_map_path,
      locked_gene_preview_path
    )
  ),
  size_bytes = as.numeric(
    file.info(
      c(
        locked_long_map_path,
        locked_gene_preview_path
      )
    )$size
  ),
  md5 = vapply(
    c(
      locked_long_map_path,
      locked_gene_preview_path
    ),
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    c(
      locked_long_map_path,
      locked_gene_preview_path
    ),
    sha256_file,
    character(1)
  ),
  source_lock_version =
    source_lock_version,
  RNA_anchor =
    RNA_anchor,
  relation_membership_status =
    "FROZEN_23499_RELATIONS",
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  locked_manifest,
  file.path(
    locked_root,
    "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"
  )
)

source_count_lock <- data.frame(
  source_id = names(
    expected_relation_counts
  ),
  occupancy_context = unname(
    context_by_source[
      names(expected_relation_counts)
    ]
  ),
  relation_n = as.integer(
    expected_relation_counts
  ),
  map_total_relation_n =
    expected_total_relations,
  unique_author_gene_symbol_n =
    expected_unique_gene_symbols,
  relation_membership_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  source_count_lock,
  file.path(
    dirs$results,
    "STEP03B3C2_v3_SOURCE_RELATION_COUNT_LOCK.csv"
  )
)

map_lock_record <- data.frame(
  lock_component = c(
    "Lock version",
    "Long-map rows",
    "Unique author gene symbols",
    "Canonical source IDs",
    "Vehicle relations",
    "Calcitriol relations",
    "Shared relations",
    "Promoter-only exceptions",
    "Invalid gene-token rows",
    "Field semantics",
    "RNA membership",
    "ATAC role",
    "Direct-target status",
    "Scientific figures"
  ),
  locked_value = c(
    peak_gene_map_lock_version,
    as.character(
      expected_total_relations
    ),
    as.character(
      expected_unique_gene_symbols
    ),
    as.character(
      expected_primary_source_n
    ),
    as.character(
      expected_relation_counts[[
        "VDR_VEHICLE_CISTROME"
      ]]
    ),
    as.character(
      expected_relation_counts[[
        "VDR_CALCITRIOL_CISTROME"
      ]]
    ),
    as.character(
      expected_relation_counts[[
        "VDR_SHARED_CISTROME"
      ]]
    ),
    "0",
    "0",
    "LOCKED_SIX_FIELD_INTERPRETATIONS",
    paste0(
      RNA_anchor,
      "_IMMUTABLE_436_CANDIDATES"
    ),
    ATAC_role,
    "NOT_PERFORMED",
    "NONE_CREATED"
  ),
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  map_lock_record,
  file.path(
    dirs$results,
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK.csv"
  )
)

# ----------------------------------------------------------------------
# Gate and next-substep authorization
# ----------------------------------------------------------------------

gate_status <-
  "READY_FOR_STEP03B3C3_RNA_GENE_ID_HARMONIZATION_AND_OCCUPANCY_INTEGRATION"

gate_table <- data.frame(
  gate_component = c(
    "STEP03B3C2 v2 gate",
    "Published cistrome source lock",
    "Field semantics",
    "Canonical source IDs",
    "Relation-count reconciliation",
    "Peak-gene relation IDs",
    "Peak-gene membership",
    "Unique author gene symbols",
    "Promoter relation consistency",
    "Gene token shape",
    "Gene-level preview",
    "Peak-gene map lock",
    "RNA membership",
    "ATAC role",
    "RNA/CUT&RUN integration",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP03B3C3 authorization",
    "STEP03B3C2 v3 gate"
  ),
  status = c(
    "PASS",
    paste0(
      "IMMUTABLE_",
      source_lock_version
    ),
    "LOCKED_SIX_OF_SIX",
    "PASS_3_OF_3",
    "PASS_3_OF_3",
    "UNIQUE_23499_OF_23499",
    "FROZEN_23499",
    "LOCKED_10042",
    "PASS_ZERO_EXCEPTIONS",
    "PASS_ZERO_INVALID_TOKENS",
    "PASS_UNIQUE_SOURCE_CONTEXT_GENE_KEYS",
    paste0(
      "LOCKED_",
      peak_gene_map_lock_version
    ),
    paste0(
      "IMMUTABLE_",
      RNA_anchor,
      "_436_CANDIDATES"
    ),
    ATAC_role,
    "NOT_STARTED",
    "NOT_PERFORMED",
    "NONE_CREATED",
    "AUTHORIZED",
    gate_status
  ),
  required_action = c(
    "Retain the repaired v2 audit record.",
    "Do not alter source peaks or coordinates.",
    "Do not reinterpret the six fields without a new version.",
    "Use canonical IDs in every downstream join.",
    "Retain source-specific counts.",
    "Use locked relation IDs as row authority.",
    "Do not add or remove peak-gene relations.",
    "Audit RNA-to-cistrome symbol matching next.",
    "Retain explicit zero-row exception table.",
    "Retain explicit zero-row invalid-token table.",
    "Use as a summary only; long map remains authoritative.",
    "Use only the locked versioned files and hashes.",
    "Do not change the RNA interaction set.",
    "Use only as optional accessible context.",
    "Begin only in STEP03B3C3 after gene-ID audit.",
    "Use occupancy-supported candidate terminology after integration QC.",
    "No Scientific Figure Skill 3 workflow is triggered.",
    "Proceed to RNA gene-ID harmonization and evidence integration.",
    "STEP03B3C2 is substep-locked."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP03B3C2_v3_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# PEAK_GENE_MAP_v1.1 lock card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Lock version: ",
      peak_gene_map_lock_version
    ),
    paste0(
      "- Upstream cistrome source: ",
      source_lock_version,
      " IMMUTABLE"
    ),
    "- Peak-gene relations: 23,499",
    "- Unique author gene symbols: 10,042",
    "- Vehicle relations: 2,848",
    "- Calcitriol relations: 19,037",
    "- Shared inherited relations: 1,614",
    "- Canonical source IDs: 3/3",
    "- Promoter-only exceptions: 0",
    "- Invalid gene-token rows: 0",
    "- Shared annotation concordance: 100%",
    "- Relation IDs: unique 23,499/23,499",
    paste0(
      "- RNA anchor: ",
      RNA_anchor,
      " IMMUTABLE; 436 interaction candidates"
    ),
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Definitive direct-target assignment: NOT PERFORMED",
    "- Scientific figures: NONE",
    "- Next: STEP03B3C3 RNA gene-ID harmonization and occupancy integration"
  ),
  file.path(
    dirs$results,
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK_CARD.md"
  )
)

write_text(
  c(
    "# STEP03B3C3 authorization",
    "",
    paste0(
      "- Authorized by: ",
      run_id
    ),
    paste0(
      "- Peak-gene map: ",
      peak_gene_map_lock_version,
      " LOCKED"
    ),
    paste0(
      "- RNA anchor: ",
      RNA_anchor,
      " IMMUTABLE"
    ),
    "- Authorization: YES",
    "- Authorized task 1: identify the authoritative 436-gene RNA interaction table and revalidate its hash.",
    "- Authorized task 2: audit exact RNA-symbol matches to the locked 10,042-symbol cistrome universe.",
    "- Authorized task 3: preserve unmatched and ambiguous symbols explicitly; do not silently alias-map.",
    "- Authorized task 4: create gene-level occupancy evidence by vehicle, calcitriol and shared context.",
    "- Authorized task 5: retain promoter support, supporting-peak counts and relation provenance separately.",
    "- Authorized task 6: construct RNA-anchored occupancy-supported candidate evidence tiers.",
    "- RNA membership changes: PROHIBITED",
    "- Peak-gene relation membership changes: PROHIBITED",
    "- media-2 derived overlap tables: CROSS-CHECK ONLY",
    "- ATAC: OPTIONAL ACCESSIBLE_CONTEXT ONLY",
    "- Definitive causal/direct-target claims: PROHIBITED",
    "- Scientific figures: create only after integrated evidence tables pass content and visual QC"
  ),
  file.path(
    dirs$control,
    "STEP03B3C3_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP03B3C2 v3 Methods addendum",
    "",
    "After correction of the source-identifier reporting mismatch, the peak-by-gene candidate map was revalidated against source-specific relation counts, unique relation identifiers and explicit exception schemas. The accepted map comprised 23,499 unique peak-gene relations across three canonical occupancy contexts and 10,042 unique author-reported gene symbols.",
    "",
    "The validated long-form map and gene-level evidence preview were copied into a versioned locked directory and recorded with MD5 and SHA-256 hashes. Peak coordinates, gene tokens, promoter flags, field semantics and RNA interaction membership were unchanged."
  ),
  file.path(
    dirs$text,
    "STEP03B3C2_v3_Methods_addendum.md"
  )
)

write_text(
  c(
    "# STEP03B3C2 v3 Results draft",
    "",
    "The locked peak-gene map contained 2,848 vehicle, 19,037 calcitriol and 1,614 shared occupancy-supported relations, for a total of 23,499 unique peak-gene relations involving 10,042 author-reported gene symbols. No promoter-only relation exceptions or malformed gene tokens were detected.",
    "",
    "All relation identifiers were unique and source-specific counts reconciled exactly with the independent peak-level token audit. The shared annotations retained complete concordance between vehicle and calcitriol source tables."
  ),
  file.path(
    dirs$text,
    "STEP03B3C2_v3_Results_draft.md"
  )
)

parameters <- list(
  schema_version = "1.1-lock",
  run_id = run_id,
  source_lock_version =
    source_lock_version,
  peak_gene_map_lock_version =
    peak_gene_map_lock_version,
  RNA_anchor = RNA_anchor,
  RNA_interaction_candidate_n = 436,
  ATAC_role = ATAC_role,
  vehicle_relation_n = 2848,
  calcitriol_relation_n = 19037,
  shared_relation_n = 1614,
  total_relation_n = 23499,
  unique_author_gene_symbol_n = 10042,
  relation_ids_unique = TRUE,
  promoter_only_exceptions = 0,
  invalid_gene_tokens = 0,
  RNA_CUTRUN_integration = FALSE,
  definitive_direct_target_assignment = FALSE,
  scientific_figures_created = FALSE,
  STEP03B3C3_authorized = TRUE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP03B3C2_v3_parameters.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C2_v3.txt"
  )
)


same_normalized_path <- function(path_a, path_b) {
  if (
    is.null(path_a) ||
      is.null(path_b) ||
      length(path_a) != 1L ||
      length(path_b) != 1L ||
      is.na(path_a) ||
      is.na(path_b)
  ) {
    return(FALSE)
  }

  normalized_a <- tryCatch(
    normalizePath(
      path_a,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
  )
}

safe_file_copy <- function(
    from,
    to,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
) {
  if (same_normalized_path(from, to)) {
    return(TRUE)
  }

  base::file.copy(
    from = from,
    to = to,
    overwrite = overwrite,
    recursive = recursive,
    copy.mode = copy.mode,
    copy.date = copy.date
  )
}

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  safe_file_copy(
    source_file,
    file.path(
      project_root,
      "04_code",
      "R",
      "40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    file.path(
      dirs$code,
      "40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

# ----------------------------------------------------------------------
# Short-path archive
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C2_v3_SOURCE_RELATION_COUNT_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C2_v3_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK_CARD.md"
  ),
  file.path(
    dirs$qc,
    "STEP03B3C2_v3_CANDIDATE_ARTIFACT_IDENTITY_RECHECK.csv"
  ),
  file.path(
    locked_root,
    "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"
  ),
  locked_long_map_path,
  locked_gene_preview_path,
  file.path(
    dirs$control,
    "STEP03B3C3_AUTHORIZATION.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C2_v3_Methods_addendum.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C2_v3_Results_draft.md"
  ),
  file.path(
    dirs$parameters,
    "STEP03B3C2_v3_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C2_v3.txt"
  ),
  file.path(
    dirs$code,
    "40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1.R"
  ),
  log_path
)

for (path_value in archive_artifacts[
  file.exists(
    archive_artifacts
  )
]) {
  section <- if (
    grepl(
      "\\.R$",
      basename(path_value)
    )
  ) {
    "code"
  } else if (
    grepl(
      "sessionInfo|parameters",
      basename(path_value)
    )
  ) {
    "env"
  } else if (
    grepl(
      "\\.log$",
      basename(path_value)
    )
  ) {
    "log"
  } else if (
    grepl(
      "Methods|Results|AUTHORIZATION",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "QC|RECHECK|MANIFEST|LOCK_CARD|COUNT_LOCK",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "qc"
  } else {
    "res"
  }

  copy_validate(
    path_value,
    file.path(
      archive_root,
      section,
      basename(path_value)
    )
  )
}

archive_files <- list.files(
  archive_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

archive_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      archive_files
    ),
    nchar(
      normalize_slashes(
        archive_root
      )
    ) + 2L
  ),
  size_bytes = as.numeric(
    file.info(
      archive_files
    )$size
  ),
  md5 = vapply(
    archive_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    archive_files,
    sha256_file,
    character(1)
  ),
  path_length = nchar(
    archive_files
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B3C2_v3_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP03B3C2 v3 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Peak-gene map lock: ",
      peak_gene_map_lock_version
    ),
    "- Peak-gene relations: 23,499",
    "- Unique author gene symbols: 10,042",
    "- Canonical source IDs: 3",
    "- Relation IDs unique: YES",
    "- Relation membership changed: NO",
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Scientific figures: NONE",
    "- STEP03B3C3 authorization: YES",
    paste0(
      "- Maximum archive path length: ",
      max(
        archive_manifest$path_length
      )
    )
  ),
  file.path(
    archive_root,
    "ctl",
    "ARCHIVE_STATUS.md"
  )
)

write_text(
  c(
    "# STEP03 Module Status",
    "",
    "- Completed substep: STEP03B3C2",
    paste0(
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- Published cistrome source: ",
      source_lock_version,
      " IMMUTABLE"
    ),
    paste0(
      "- Peak-gene map: ",
      peak_gene_map_lock_version,
      " LOCKED"
    ),
    paste0(
      "- RNA anchor: ",
      RNA_anchor,
      " IMMUTABLE; 436 candidates"
    ),
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Definitive direct-target assignment: NOT PERFORMED",
    "- Scientific figures: NONE",
    "- Active next substep: STEP03B3C3",
    "- Full STEP03 module: NOT_LOCKED"
  ),
  file.path(
    dirs$control,
    "STATUS.md"
  )
)

log_event(
  "Completed",
  run_id,
  "Gate:",
  gate_status
)

cat("\n============================================================\n")
cat("STEP03B3C2 v3 PEAK_GENE_MAP_v1.1 lock completed.\n")
cat(
  "Peak-gene map: ",
  peak_gene_map_lock_version,
  " LOCKED\n"
)
cat("Peak-gene relations: 23,499\n")
cat("Unique author gene symbols: 10,042\n")
cat("Vehicle relations: 2,848\n")
cat("Calcitriol relations: 19,037\n")
cat("Shared relations: 1,614\n")
cat("Canonical source IDs: 3/3\n")
cat("Relation IDs unique: 23,499/23,499\n")
cat("Promoter-only exceptions: 0\n")
cat("Invalid gene-token rows: 0\n")
cat(
  "RNA anchor: ",
  RNA_anchor,
  " IMMUTABLE; 436 candidates\n"
)
cat("RNA/CUT&RUN integration: NOT STARTED\n")
cat("Definitive direct-target assignment: NOT PERFORMED\n")
cat("Scientific figures: NONE\n")
cat("STEP03B3C3 authorization: YES\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 4: 44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy
# ==============================================================================
# ============================================================
# STEP03B3C3A v4 - Lock RNA source and RNA/cistrome gene-ID policy
# Nutrition_Endometrial_VDR_v2
#
# Lock decisions:
# - RNA source: STEP03B1_v10_RNA_interaction_program_FROZEN.csv
# - RNA identifier: feature_id
# - 279 exact RNA/cistrome symbol matches are retained unchanged.
# - Five deterministic dot-to-hyphen equivalents are manually accepted:
#     DPP4.DT   -> DPP4-DT
#     HHIP.AS1  -> HHIP-AS1
#     KRTAP4.8  -> KRTAP4-8
#     MUC12.AS1 -> MUC12-AS1
#     RGMB.AS1  -> RGMB-AS1
# - The remaining 152 RNA symbols retain no current cistrome match.
#
# Scientific meaning:
# - the 436-member RNA program remains immutable;
# - the original feature_id is never overwritten;
# - a separate harmonized cistrome-join symbol is created;
# - 284 RNA candidates are authorized for exact cistrome joining;
# - 152 remain RNA-only candidates without current occupancy support.
#
# This substep does NOT:
# - add or remove RNA candidates;
# - add or remove peak-gene relations;
# - perform broad alias mapping;
# - perform RNA/CUT&RUN integration;
# - assign definitive direct VDR targets;
# - generate scientific figures.
#
# Run from project root:
# source(
#   "04_code/R/44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

module_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

run_id <- paste0(
  "STEP03B3C3A_v4_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

RNA_anchor <- "RNA_PROGRAM_v1.0"
RNA_source_lock_version <- "RNA_INTERACTION_SOURCE_v1.0"
gene_ID_policy_lock_version <- "RNA_CISTROME_GENE_ID_POLICY_v1.0"
gene_ID_map_lock_version <- "RNA_CISTROME_GENE_ID_MAP_v1.0"
peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"
ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_RNA_n <- 436L
expected_exact_n <- 279L
expected_manual_mapping_n <- 5L
expected_join_authorized_n <- 284L
expected_unmatched_n <- 152L
expected_peak_gene_relation_n <- 23499L
expected_cistrome_gene_n <- 10042L

required_packages <- c(
  "data.table",
  "digest",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package(s) missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ". This script does not install or update packages."
  )
}

dirs <- list(
  control = file.path(module_root, "00_module_control"),
  inputs = file.path(module_root, "01_inputs"),
  code = file.path(module_root, "02_code"),
  parameters = file.path(module_root, "03_parameters"),
  environment = file.path(module_root, "04_environment"),
  execution = file.path(module_root, "05_execution"),
  results = file.path(module_root, "06_results"),
  figures = file.path(module_root, "07_figures_tables"),
  qc = file.path(module_root, "08_qc"),
  text = file.path(module_root, "09_manuscript_text"),
  archive_parent = file.path(module_root, "10_archive_work")
)

required_dir_names <- c(
  "control",
  "inputs",
  "code",
  "parameters",
  "environment",
  "execution",
  "results",
  "figures",
  "qc",
  "text",
  "archive_parent"
)

missing_dir_names <- setdiff(
  required_dir_names,
  names(dirs)
)

invalid_dir_names <- names(dirs)[
  !vapply(
    dirs,
    function(path_value) {
      is.character(path_value) &&
        length(path_value) == 1L &&
        !is.na(path_value) &&
        nzchar(path_value)
    },
    logical(1)
  )
]

if (
  length(missing_dir_names) > 0L ||
    length(invalid_dir_names) > 0L
) {
  stop(
    "Invalid directory registry. Missing: ",
    paste(
      missing_dir_names,
      collapse = ", "
    ),
    "; invalid: ",
    paste(
      invalid_dir_names,
      collapse = ", "
    )
  )
}

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

locked_root <- file.path(
  dirs$results,
  "STEP03B3C3A_locked_RNA_gene_ID",
  gene_ID_policy_lock_version
)

dir.create(
  locked_root,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S3B3C3Av4_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

for (directory_path in c(
  file.path(archive_root, "res"),
  file.path(archive_root, "qc"),
  file.path(archive_root, "txt"),
  file.path(archive_root, "code"),
  file.path(archive_root, "env"),
  file.path(archive_root, "log"),
  file.path(archive_root, "ctl")
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  dirs$execution,
  paste0(
    run_id,
    "_execution.log"
  )
)

log_event <- function(...) {
  message_text <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(message_text, "\n")

  cat(
    message_text,
    "\n",
    file = log_path,
    append = TRUE
  )
}

write_text <- function(lines, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
}

write_csv_utf8 <- function(x, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

parse_bool <- function(x) {
  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1"
  )
}

clean_symbol <- function(x) {
  value <- trimws(
    as.character(x)
  )

  value[
    is.na(value) |
      !nzchar(value) |
      toupper(value) %in%
        c(
          "NA",
          "N/A",
          "NULL",
          "NONE",
          "."
        )
  ] <- NA_character_

  value
}

copy_validate <- function(
    source_path,
    destination_path
) {
  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (nchar(destination_path) >= 240L) {
    stop(
      "Archive destination exceeds 240 characters:\n",
      destination_path
    )
  }

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(
        file.info(source_path)$size
      ),
      as.numeric(
        file.info(destination_path)$size
      )
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed:\n",
      source_path
    )
  }

  TRUE
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Required accepted source and audit files
# ----------------------------------------------------------------------

v2_audit_root <- file.path(
  dirs$results,
  "STEP03B3C3A_v2_RNA_source_and_gene_ID_audit"
)

v3_audit_root <- file.path(
  dirs$results,
  "STEP03B3C3A_v3_punctuation_version_audit"
)

map_locked_root <- file.path(
  dirs$results,
  "STEP03B3C2_locked_peak_gene_map",
  peak_gene_map_lock_version
)

required_paths <- c(
  v2_gate = file.path(
    dirs$results,
    "STEP03B3C3A_v2_GATE_TABLE.csv"
  ),
  v2_source_candidate = file.path(
    dirs$results,
    "STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv"
  ),
  v2_gene_audit = file.path(
    v2_audit_root,
    "STEP03B3C3A_v2_RNA_GENE_ID_AUDIT.csv"
  ),
  v2_match_summary = file.path(
    v2_audit_root,
    "STEP03B3C3A_v2_RNA_CISTROME_MATCH_SUMMARY.csv"
  ),
  v3_gate = file.path(
    dirs$results,
    "STEP03B3C3A_v3_GATE_TABLE.csv"
  ),
  v3_candidates = file.path(
    v3_audit_root,
    "STEP03B3C3A_v3_PUNCTUATION_VERSION_CANDIDATES.csv"
  ),
  v3_classification = file.path(
    v3_audit_root,
    "STEP03B3C3A_v3_UNMATCHED_SYMBOL_CLASSIFICATION.csv"
  ),
  RNA_source = file.path(
    dirs$results,
    "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"
  ),
  map_manifest = file.path(
    map_locked_root,
    "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"
  )
)

missing_required_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_required_paths) > 0L) {
  stop(
    "Required STEP03B3C3A v4 input is missing:\n",
    paste(
      missing_required_paths,
      collapse = "\n"
    )
  )
}

v2_gate <- read.csv(
  required_paths[["v2_gate"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_candidate <- read.csv(
  required_paths[["v2_source_candidate"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_audit <- read.csv(
  required_paths[["v2_gene_audit"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

match_summary <- read.csv(
  required_paths[["v2_match_summary"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v3_gate <- read.csv(
  required_paths[["v3_gate"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v3_candidates <- read.csv(
  required_paths[["v3_candidates"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v3_classification <- read.csv(
  required_paths[["v3_classification"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map_manifest <- read.csv(
  required_paths[["map_manifest"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v2_gate_value <- v2_gate$status[
  v2_gate$gate_component ==
    "STEP03B3C3A v2 gate"
]

v3_gate_value <- v3_gate$status[
  v3_gate$gate_component ==
    "STEP03B3C3A v3 gate"
]

if (
  length(v2_gate_value) != 1L ||
    v2_gate_value !=
      "PROCEED_TO_MANUAL_RNA_SOURCE_AND_GENE_ID_LOCK"
) {
  stop(
    "STEP03B3C3A v2 gate is not accepted."
  )
}

if (
  length(v3_gate_value) != 1L ||
    v3_gate_value !=
      "PROCEED_TO_MANUAL_PUNCTUATION_VERSION_CANDIDATE_REVIEW"
) {
  stop(
    "STEP03B3C3A v3 gate is not accepted."
  )
}

# ----------------------------------------------------------------------
# Validate source candidate identity
# ----------------------------------------------------------------------

if (
  nrow(source_candidate) != 1L ||
    source_candidate$gene_symbol_column !=
      "feature_id" ||
    source_candidate$n_rows !=
      expected_RNA_n ||
    source_candidate$unique_gene_n !=
      expected_RNA_n
) {
  stop(
    "RNA source candidate metadata failed."
  )
}

current_source_md5 <- md5_file(
  required_paths[["RNA_source"]]
)

current_source_sha256 <- sha256_file(
  required_paths[["RNA_source"]]
)

if (
  source_candidate$md5 !=
    current_source_md5 ||
    source_candidate$sha256 !=
      current_source_sha256
) {
  stop(
    "RNA source candidate hash differs from the frozen source file."
  )
}

RNA_table <- read.csv(
  required_paths[["RNA_source"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(RNA_table) !=
    expected_RNA_n ||
    !"feature_id" %in%
      names(RNA_table)
) {
  stop(
    "Frozen RNA source failed row or feature_id validation."
  )
}

RNA_symbols <- clean_symbol(
  RNA_table$feature_id
)

if (
  any(is.na(RNA_symbols)) ||
    length(
      unique(RNA_symbols)
    ) !=
      expected_RNA_n
) {
  stop(
    "Frozen RNA feature_id values are not 436 unique nonmissing symbols."
  )
}

# ----------------------------------------------------------------------
# Validate v2 partition
# ----------------------------------------------------------------------

if (
  nrow(gene_audit) !=
    expected_RNA_n ||
    length(
      unique(
        gene_audit$RNA_row
      )
    ) !=
      expected_RNA_n ||
    !identical(
      as.character(
        gene_audit$RNA_gene_symbol
      ),
      as.character(
        RNA_symbols
      )
    )
) {
  stop(
    "The v2 gene audit does not exactly track the frozen RNA source."
  )
}

partition_n <- rowSums(
  gene_audit[
    c(
      "exact_cistrome_match",
      "case_only_cistrome_match",
      "no_cistrome_match"
    )
  ]
)

if (
  !all(
    partition_n == 1L
  ) ||
    sum(
      parse_bool(
        gene_audit$exact_cistrome_match
      )
    ) !=
      expected_exact_n ||
    sum(
      parse_bool(
        gene_audit$case_only_cistrome_match
      )
    ) !=
      0L ||
    sum(
      parse_bool(
        gene_audit$no_cistrome_match
      )
    ) !=
      157L
) {
  stop(
    "The v2 exact/case-only/unmatched partition failed."
  )
}

# ----------------------------------------------------------------------
# Validate five manually reviewed mappings
# ----------------------------------------------------------------------

manual_decisions <- data.frame(
  RNA_gene_symbol_original = c(
    "DPP4.DT",
    "HHIP.AS1",
    "KRTAP4.8",
    "MUC12.AS1",
    "RGMB.AS1"
  ),
  cistrome_gene_symbol_harmonized = c(
    "DPP4-DT",
    "HHIP-AS1",
    "KRTAP4-8",
    "MUC12-AS1",
    "RGMB-AS1"
  ),
  decision = "ACCEPT",
  mapping_rule =
    "MANUAL_ACCEPTED_DOT_TO_HYPHEN_OFFICIAL_SYMBOL_EQUIVALENT",
  local_candidate_evidence =
    "EXACT_DOT_TO_HYPHEN_AND_UNIQUE_PUNCTUATION_KEY",
  external_nomenclature_review =
    "REVIEWED_AGAINST_CURRENT_HGNC_NCBI_SYMBOL_RECORDS",
  automatic_mapping_applied =
    FALSE,
  lock_status =
    "LOCKED_MANUAL_DECISION",
  stringsAsFactors = FALSE
)

observed_unique_pairs <- unique(
  v3_candidates[
    c(
      "RNA_gene_symbol",
      "candidate_cistrome_symbol"
    )
  ]
)

observed_unique_pairs <- observed_unique_pairs[
  order(
    observed_unique_pairs$RNA_gene_symbol
  ),
  ,
  drop = FALSE
]

expected_pairs <- manual_decisions[
  c(
    "RNA_gene_symbol_original",
    "cistrome_gene_symbol_harmonized"
  )
]

expected_pairs <- expected_pairs[
  order(
    expected_pairs$RNA_gene_symbol_original
  ),
  ,
  drop = FALSE
]

if (
  nrow(observed_unique_pairs) !=
    expected_manual_mapping_n ||
    !identical(
      as.character(
        observed_unique_pairs$RNA_gene_symbol
      ),
      as.character(
        expected_pairs$RNA_gene_symbol_original
      )
    ) ||
    !identical(
      as.character(
        observed_unique_pairs$candidate_cistrome_symbol
      ),
      as.character(
        expected_pairs$cistrome_gene_symbol_harmonized
      )
    )
) {
  stop(
    "The v3 deterministic candidate pairs differ from the reviewed five mappings."
  )
}

candidate_method_check <- aggregate(
  candidate_method ~
    RNA_gene_symbol +
    candidate_cistrome_symbol,
  data = v3_candidates,
  FUN = function(x) {
    paste(
      sort(
        unique(x)
      ),
      collapse = ";"
    )
  }
)

required_method_string <- paste(
  sort(
    c(
      "EXACT_DOT_TO_HYPHEN",
      "PUNCTUATION_INSENSITIVE_KEY"
    )
  ),
  collapse = ";"
)

if (
  nrow(candidate_method_check) !=
    expected_manual_mapping_n ||
    !all(
      candidate_method_check$candidate_method ==
        required_method_string
    )
) {
  stop(
    "Each reviewed mapping must be supported by both deterministic methods."
  )
}

reviewed_classification <- v3_classification[
  v3_classification$RNA_gene_symbol %in%
    manual_decisions$RNA_gene_symbol_original,
  ,
  drop = FALSE
]

if (
  nrow(reviewed_classification) !=
    expected_manual_mapping_n ||
    !all(
      reviewed_classification$candidate_symbol_n ==
        1L
    ) ||
    !all(
      reviewed_classification$review_class ==
        "ONE_UNIQUE_DETERMINISTIC_CANDIDATE"
    )
) {
  stop(
    "The five reviewed mappings are not unique deterministic candidates."
  )
}

# ----------------------------------------------------------------------
# Validate cistrome universe
# ----------------------------------------------------------------------

if (
  nrow(map_manifest) != 2L ||
    !all(
      map_manifest$lock_version ==
        peak_gene_map_lock_version
    ) ||
    !all(
      map_manifest$lock_status ==
        "LOCKED"
    )
) {
  stop(
    "PEAK_GENE_MAP_v1.1 manifest failed."
  )
}

map_manifest$file_exists <- file.exists(
  map_manifest$path
)

map_manifest$observed_md5 <- vapply(
  map_manifest$path,
  md5_file,
  character(1)
)

map_manifest$observed_sha256 <- vapply(
  map_manifest$path,
  sha256_file,
  character(1)
)

map_manifest$identity_pass <- with(
  map_manifest,
  file_exists &
    md5 ==
      observed_md5 &
    sha256 ==
      observed_sha256
)

if (!all(
  map_manifest$identity_pass
)) {
  stop(
    "A locked peak-gene artifact failed hash validation."
  )
}

map_long_path <- map_manifest$path[
  map_manifest$artifact_id ==
    "PEAK_GENE_LONG"
]

if (length(map_long_path) != 1L) {
  stop(
    "Locked peak-gene long map could not be uniquely identified."
  )
}

peak_gene_map <- data.table::fread(
  map_long_path,
  sep = "\t",
  data.table = FALSE,
  showProgress = FALSE
)

if (
  nrow(peak_gene_map) !=
    expected_peak_gene_relation_n
) {
  stop(
    "Locked peak-gene relation count differs from 23,499."
  )
}

cistrome_symbols <- sort(
  unique(
    clean_symbol(
      peak_gene_map$gene_symbol_author
    )
  )
)

cistrome_symbols <- cistrome_symbols[
  !is.na(cistrome_symbols)
]

if (
  length(cistrome_symbols) !=
    expected_cistrome_gene_n
) {
  stop(
    "Locked cistrome gene-symbol universe differs from 10,042."
  )
}

if (!all(
  manual_decisions$cistrome_gene_symbol_harmonized %in%
    cistrome_symbols
)) {
  stop(
    "At least one accepted harmonized symbol is absent from the locked cistrome universe."
  )
}

# ----------------------------------------------------------------------
# Build immutable 436-row gene-ID map
# ----------------------------------------------------------------------

gene_ID_map <- data.frame(
  RNA_row =
    gene_audit$RNA_row,
  RNA_gene_symbol_original =
    gene_audit$RNA_gene_symbol,
  RNA_gene_symbol_cistrome_join =
    NA_character_,
  mapping_status =
    NA_character_,
  mapping_rule =
    NA_character_,
  exact_original_match =
    parse_bool(
      gene_audit$exact_cistrome_match
    ),
  manual_mapping_accepted =
    FALSE,
  cistrome_join_authorized =
    FALSE,
  RNA_membership_status =
    "FROZEN_IN_RNA_PROGRAM_v1.0",
  automatic_alias_mapping =
    FALSE,
  automatic_punctuation_mapping =
    FALSE,
  stringsAsFactors = FALSE
)

exact_index <- gene_ID_map$exact_original_match

gene_ID_map$RNA_gene_symbol_cistrome_join[
  exact_index
] <- gene_ID_map$RNA_gene_symbol_original[
  exact_index
]

gene_ID_map$mapping_status[
  exact_index
] <- "EXACT_AUTHOR_SYMBOL_MATCH"

gene_ID_map$mapping_rule[
  exact_index
] <- "IDENTITY"

gene_ID_map$cistrome_join_authorized[
  exact_index
] <- TRUE

for (decision_index in seq_len(
  nrow(manual_decisions)
)) {
  original_symbol <-
    manual_decisions$RNA_gene_symbol_original[
      decision_index
    ]

  harmonized_symbol <-
    manual_decisions$cistrome_gene_symbol_harmonized[
      decision_index
    ]

  map_index <-
    gene_ID_map$RNA_gene_symbol_original ==
      original_symbol

  if (sum(map_index) != 1L) {
    stop(
      "Manual decision does not map to exactly one RNA row: ",
      original_symbol
    )
  }

  gene_ID_map$RNA_gene_symbol_cistrome_join[
    map_index
  ] <- harmonized_symbol

  gene_ID_map$mapping_status[
    map_index
  ] <- "MANUAL_ACCEPTED_PUNCTUATION_EQUIVALENT"

  gene_ID_map$mapping_rule[
    map_index
  ] <-
    "DOT_TO_HYPHEN_OFFICIAL_SYMBOL_EQUIVALENT"

  gene_ID_map$manual_mapping_accepted[
    map_index
  ] <- TRUE

  gene_ID_map$cistrome_join_authorized[
    map_index
  ] <- TRUE
}

unmatched_index <- is.na(
  gene_ID_map$mapping_status
)

gene_ID_map$mapping_status[
  unmatched_index
] <- "NO_CURRENT_CISTROME_SYMBOL_MATCH"

gene_ID_map$mapping_rule[
  unmatched_index
] <- "NONE"

if (
  nrow(gene_ID_map) !=
    expected_RNA_n ||
    sum(
      gene_ID_map$mapping_status ==
        "EXACT_AUTHOR_SYMBOL_MATCH"
    ) !=
      expected_exact_n ||
    sum(
      gene_ID_map$manual_mapping_accepted
    ) !=
      expected_manual_mapping_n ||
    sum(
      gene_ID_map$cistrome_join_authorized
    ) !=
      expected_join_authorized_n ||
    sum(
      gene_ID_map$mapping_status ==
        "NO_CURRENT_CISTROME_SYMBOL_MATCH"
    ) !=
      expected_unmatched_n
) {
  stop(
    "The 436-row gene-ID map failed count reconciliation."
  )
}

authorized_join_symbols <-
  gene_ID_map$RNA_gene_symbol_cistrome_join[
    gene_ID_map$cistrome_join_authorized
  ]

if (
  any(is.na(authorized_join_symbols)) ||
    length(
      unique(
        authorized_join_symbols
      )
    ) !=
      expected_join_authorized_n ||
    !all(
      authorized_join_symbols %in%
        cistrome_symbols
    )
) {
  stop(
    "Authorized cistrome join symbols failed completeness, uniqueness or universe validation."
  )
}

unmatched_locked <- gene_ID_map[
  !gene_ID_map$cistrome_join_authorized,
  ,
  drop = FALSE
]

# ----------------------------------------------------------------------
# Write locked source and gene-ID artifacts
# ----------------------------------------------------------------------

RNA_source_lock <- data.frame(
  lock_version =
    RNA_source_lock_version,
  RNA_anchor =
    RNA_anchor,
  source_file =
    normalize_slashes(
      required_paths[["RNA_source"]]
    ),
  source_file_name =
    basename(
      required_paths[["RNA_source"]]
    ),
  gene_symbol_column =
    "feature_id",
  n_rows =
    expected_RNA_n,
  unique_gene_n =
    expected_RNA_n,
  size_bytes = as.numeric(
    file.info(
      required_paths[["RNA_source"]]
    )$size
  ),
  md5 =
    current_source_md5,
  sha256 =
    current_source_sha256,
  membership_status =
    "IMMUTABLE_436",
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

gene_ID_decision_lock <- manual_decisions

gene_ID_decision_lock$policy_lock_version <-
  gene_ID_policy_lock_version

gene_ID_map_path <- file.path(
  locked_root,
  "RNA_CISTROME_GENE_ID_MAP_v1.0.csv"
)

unmatched_locked_path <- file.path(
  locked_root,
  "RNA_CISTROME_UNMATCHED_152_v1.0.csv"
)

RNA_source_lock_path <- file.path(
  locked_root,
  "RNA_INTERACTION_SOURCE_v1.0_LOCK.csv"
)

gene_ID_decision_lock_path <- file.path(
  locked_root,
  "RNA_CISTROME_GENE_ID_MANUAL_DECISIONS_v1.0.csv"
)

write_csv_utf8(
  gene_ID_map,
  gene_ID_map_path
)

write_csv_utf8(
  unmatched_locked,
  unmatched_locked_path
)

write_csv_utf8(
  RNA_source_lock,
  RNA_source_lock_path
)

write_csv_utf8(
  gene_ID_decision_lock,
  gene_ID_decision_lock_path
)

locked_artifact_paths <- c(
  RNA_source_lock_path,
  gene_ID_map_path,
  unmatched_locked_path,
  gene_ID_decision_lock_path
)

locked_manifest <- data.frame(
  artifact_id = c(
    "RNA_INTERACTION_SOURCE_LOCK",
    "RNA_CISTROME_GENE_ID_MAP",
    "RNA_CISTROME_UNMATCHED_SET",
    "RNA_CISTROME_MANUAL_DECISIONS"
  ),
  lock_version = c(
    RNA_source_lock_version,
    gene_ID_map_lock_version,
    gene_ID_policy_lock_version,
    gene_ID_policy_lock_version
  ),
  path =
    normalize_slashes(
      locked_artifact_paths
    ),
  size_bytes = as.numeric(
    file.info(
      locked_artifact_paths
    )$size
  ),
  md5 = vapply(
    locked_artifact_paths,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    locked_artifact_paths,
    sha256_file,
    character(1)
  ),
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  locked_manifest,
  file.path(
    locked_root,
    "RNA_CISTROME_GENE_ID_POLICY_v1.0_LOCKED_MANIFEST.csv"
  )
)

policy_summary <- data.frame(
  policy_component = c(
    "RNA source",
    "RNA identifier",
    "RNA membership",
    "Exact identity joins",
    "Manual punctuation joins",
    "Total cistrome-authorized RNA symbols",
    "RNA-only unmatched symbols",
    "Automatic aliases",
    "Automatic capitalization changes",
    "Automatic punctuation changes",
    "Peak-gene membership",
    "ATAC role",
    "Direct-target assignment",
    "Scientific figures"
  ),
  locked_value = c(
    RNA_source_lock_version,
    "feature_id",
    "436_IMMUTABLE",
    "279",
    "5",
    "284",
    "152",
    "0",
    "0",
    "0",
    "23499_IMMUTABLE",
    ATAC_role,
    "NOT_PERFORMED",
    "NONE_CREATED"
  ),
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  policy_summary,
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK.csv"
  )
)

# ----------------------------------------------------------------------
# Gate and next-step authorization
# ----------------------------------------------------------------------

gate_status <-
  "READY_FOR_STEP03B3C3B_RNA_OCCUPANCY_EVIDENCE_INTEGRATION"

gate_table <- data.frame(
  gate_component = c(
    "STEP03B3C3A v2 source audit",
    "STEP03B3C3A v3 punctuation audit",
    "RNA source identity",
    "RNA identifier field",
    "RNA membership",
    "Exact gene-ID joins",
    "Manual punctuation decisions",
    "Manual mapping uniqueness",
    "Cistrome-authorized RNA symbols",
    "RNA-only unmatched symbols",
    "Automatic alias mapping",
    "Automatic capitalization mapping",
    "Automatic punctuation mapping",
    "Peak-gene map",
    "ATAC role",
    "RNA/CUT&RUN integration",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP03B3C3B authorization",
    "STEP03B3C3A v4 gate"
  ),
  status = c(
    "PASS",
    "PASS",
    paste0(
      "LOCKED_",
      RNA_source_lock_version
    ),
    "LOCKED_feature_id",
    "IMMUTABLE_436",
    "LOCKED_279",
    "LOCKED_5",
    "PASS_5_OF_5",
    "LOCKED_284",
    "LOCKED_152",
    "ZERO",
    "ZERO",
    "ZERO_AUTOMATIC",
    paste0(
      "IMMUTABLE_",
      peak_gene_map_lock_version,
      "_23499"
    ),
    ATAC_role,
    "NOT_STARTED",
    "NOT_PERFORMED",
    "NONE_CREATED",
    "AUTHORIZED",
    gate_status
  ),
  required_action = c(
    "Retain v2 source audit.",
    "Retain v3 deterministic candidate audit.",
    "Use only the locked source and hash.",
    "Do not replace feature_id in the RNA source.",
    "Do not add or remove RNA rows.",
    "Join exact symbols by identity.",
    "Use only the five locked manual decisions.",
    "Do not expand to broad punctuation conversion.",
    "Use the locked harmonized join field.",
    "Retain as RNA candidates without current occupancy support.",
    "Broad alias mapping remains prohibited.",
    "Case-only mapping remains prohibited unless separately reviewed.",
    "Only five manual mappings are authorized.",
    "Do not alter peak-gene membership.",
    "Use only as optional accessible context.",
    "Begin only in STEP03B3C3B.",
    "Use occupancy-supported candidate terminology after integration QC.",
    "No Scientific Figure Skill 3 workflow is triggered.",
    "Proceed to gene-level occupancy evidence integration.",
    "STEP03B3C3A is substep-locked."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP03B3C3A v4 RNA source and gene-ID policy lock card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- RNA source: ",
      RNA_source_lock_version,
      " LOCKED"
    ),
    "- RNA identifier: feature_id LOCKED",
    "- RNA membership: 436 IMMUTABLE",
    "- Exact cistrome joins: 279",
    "- Manually accepted punctuation equivalents: 5",
    "- Total cistrome-authorized RNA symbols: 284",
    "- RNA-only unmatched symbols: 152",
    "- Automatic alias mappings: 0",
    "- Automatic capitalization mappings: 0",
    "- Automatic punctuation mappings: 0",
    paste0(
      "- Peak-gene map: ",
      peak_gene_map_lock_version,
      " IMMUTABLE; 23,499 relations"
    ),
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Definitive direct-target assignment: NOT PERFORMED",
    "- Scientific figures: NONE",
    "- STEP03B3C3B authorization: YES",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK_CARD.md"
  )
)

write_text(
  c(
    "# STEP03B3C3B authorization",
    "",
    paste0(
      "- Authorized by: ",
      run_id
    ),
    paste0(
      "- RNA source: ",
      RNA_source_lock_version,
      " LOCKED"
    ),
    paste0(
      "- RNA/cistrome gene-ID policy: ",
      gene_ID_policy_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Gene-ID map: ",
      gene_ID_map_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Peak-gene map: ",
      peak_gene_map_lock_version,
      " LOCKED"
    ),
    "- Authorization: YES",
    "- Authorized task 1: join the 284 cistrome-authorized RNA candidates to the locked peak-gene map.",
    "- Authorized task 2: summarize vehicle, calcitriol and shared occupancy separately.",
    "- Authorized task 3: retain supporting-peak counts and promoter evidence separately.",
    "- Authorized task 4: retain the 152 unmatched RNA candidates as RNA-only rows.",
    "- Authorized task 5: construct occupancy-supported candidate evidence tiers without changing RNA membership.",
    "- Prohibited: broad alias remapping, additional punctuation conversion or nearest-gene invention.",
    "- Prohibited: definitive direct-target or causal wording.",
    "- ATAC: optional accessible_context only after CUT&RUN-supported target selection.",
    "- Scientific figures: only after the integrated evidence table passes content QC."
  ),
  file.path(
    dirs$control,
    "STEP03B3C3B_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP03B3C3A v4 Methods addendum",
    "",
    "The frozen RNA interaction source and `feature_id` identifier field were locked after validation against the accepted RNA program and archived source identity. Of 436 RNA interaction candidates, 279 matched the locked cistrome gene-symbol universe exactly.",
    "",
    "Five additional one-to-one punctuation equivalents were accepted after deterministic dot-to-hyphen comparison and manual nomenclature review: DPP4.DT/DPP4-DT, HHIP.AS1/HHIP-AS1, KRTAP4.8/KRTAP4-8, MUC12.AS1/MUC12-AS1 and RGMB.AS1/RGMB-AS1. Original RNA identifiers were retained, and separate harmonized cistrome-join symbols were created. No broad alias or automatic punctuation mapping was applied."
  ),
  file.path(
    dirs$text,
    "STEP03B3C3A_v4_Methods_addendum.md"
  )
)

write_text(
  c(
    "# STEP03B3C3A v4 Results draft",
    "",
    "The locked RNA/cistrome gene-ID map retained all 436 RNA interaction candidates. A total of 279 candidates had exact symbol matches and five had manually accepted one-to-one punctuation equivalents, yielding 284 RNA candidates authorized for cistrome joining. The remaining 152 candidates were retained as RNA-only candidates without current author-annotated occupancy support.",
    "",
    "No automatic alias, capitalization or broad punctuation mapping was performed, and the locked peak-gene map remained unchanged."
  ),
  file.path(
    dirs$text,
    "STEP03B3C3A_v4_Results_draft.md"
  )
)

parameters <- list(
  schema_version = "4.0-lock",
  run_id = run_id,
  RNA_anchor = RNA_anchor,
  RNA_source_lock_version =
    RNA_source_lock_version,
  gene_ID_policy_lock_version =
    gene_ID_policy_lock_version,
  gene_ID_map_lock_version =
    gene_ID_map_lock_version,
  peak_gene_map_lock_version =
    peak_gene_map_lock_version,
  RNA_n = expected_RNA_n,
  exact_join_n = expected_exact_n,
  manual_punctuation_join_n =
    expected_manual_mapping_n,
  cistrome_join_authorized_n =
    expected_join_authorized_n,
  unmatched_RNA_only_n =
    expected_unmatched_n,
  automatic_alias_mapping =
    FALSE,
  automatic_capitalization_mapping =
    FALSE,
  automatic_punctuation_mapping =
    FALSE,
  RNA_membership_changed =
    FALSE,
  peak_gene_membership_changed =
    FALSE,
  RNA_CUTRUN_integration =
    FALSE,
  definitive_direct_target_assignment =
    FALSE,
  scientific_figures_created =
    FALSE,
  STEP03B3C3B_authorized =
    TRUE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP03B3C3A_v4_parameters.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C3A_v4.txt"
  )
)


same_normalized_path <- function(path_a, path_b) {
  if (
    is.null(path_a) ||
      is.null(path_b) ||
      length(path_a) != 1L ||
      length(path_b) != 1L ||
      is.na(path_a) ||
      is.na(path_b)
  ) {
    return(FALSE)
  }

  normalized_a <- tryCatch(
    normalizePath(
      path_a,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
  )
}

safe_file_copy <- function(
    from,
    to,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
) {
  if (same_normalized_path(from, to)) {
    return(TRUE)
  }

  base::file.copy(
    from = from,
    to = to,
    overwrite = overwrite,
    recursive = recursive,
    copy.mode = copy.mode,
    copy.date = copy.date
  )
}

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  safe_file_copy(
    source_file,
    file.path(
      project_root,
      "04_code",
      "R",
      "44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    file.path(
      dirs$code,
      "44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

# ----------------------------------------------------------------------
# Short-path archive
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK_CARD.md"
  ),
  RNA_source_lock_path,
  gene_ID_map_path,
  unmatched_locked_path,
  gene_ID_decision_lock_path,
  file.path(
    locked_root,
    "RNA_CISTROME_GENE_ID_POLICY_v1.0_LOCKED_MANIFEST.csv"
  ),
  file.path(
    dirs$control,
    "STEP03B3C3B_AUTHORIZATION.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C3A_v4_Methods_addendum.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C3A_v4_Results_draft.md"
  ),
  file.path(
    dirs$parameters,
    "STEP03B3C3A_v4_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C3A_v4.txt"
  ),
  file.path(
    dirs$code,
    "44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy.R"
  ),
  log_path
)

for (path_value in archive_artifacts[
  file.exists(
    archive_artifacts
  )
]) {
  section <- if (
    grepl(
      "\\.R$",
      basename(path_value)
    )
  ) {
    "code"
  } else if (
    grepl(
      "sessionInfo|parameters",
      basename(path_value)
    )
  ) {
    "env"
  } else if (
    grepl(
      "\\.log$",
      basename(path_value)
    )
  ) {
    "log"
  } else if (
    grepl(
      "Methods|Results|AUTHORIZATION",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "LOCK|MANIFEST|CARD|DECISION|UNMATCHED",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "qc"
  } else {
    "res"
  }

  copy_validate(
    path_value,
    file.path(
      archive_root,
      section,
      basename(path_value)
    )
  )
}

archive_files <- list.files(
  archive_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

archive_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      archive_files
    ),
    nchar(
      normalize_slashes(
        archive_root
      )
    ) + 2L
  ),
  size_bytes = as.numeric(
    file.info(
      archive_files
    )$size
  ),
  md5 = vapply(
    archive_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    archive_files,
    sha256_file,
    character(1)
  ),
  path_length = nchar(
    archive_files
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B3C3A_v4_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP03B3C3A v4 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- RNA source: ",
      RNA_source_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Gene-ID policy: ",
      gene_ID_policy_lock_version,
      " LOCKED"
    ),
    "- RNA membership: 436 IMMUTABLE",
    "- Exact joins: 279",
    "- Manual punctuation joins: 5",
    "- Total cistrome-authorized: 284",
    "- RNA-only unmatched: 152",
    "- Automatic mappings: 0",
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Scientific figures: NONE",
    "- STEP03B3C3B authorization: YES",
    paste0(
      "- Maximum archive path length: ",
      max(
        archive_manifest$path_length
      )
    )
  ),
  file.path(
    archive_root,
    "ctl",
    "ARCHIVE_STATUS.md"
  )
)

write_text(
  c(
    "# STEP03 Module Status",
    "",
    "- Completed substep: STEP03B3C3A",
    paste0(
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- RNA source: ",
      RNA_source_lock_version,
      " LOCKED"
    ),
    paste0(
      "- RNA/cistrome gene-ID policy: ",
      gene_ID_policy_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Gene-ID map: ",
      gene_ID_map_lock_version,
      " LOCKED"
    ),
    "- RNA membership: RNA_PROGRAM_v1.0 IMMUTABLE; 436",
    "- Cistrome-authorized RNA symbols: 284",
    "- RNA-only unmatched symbols: 152",
    paste0(
      "- Peak-gene map: ",
      peak_gene_map_lock_version,
      " IMMUTABLE"
    ),
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- RNA/CUT&RUN integration: NOT STARTED",
    "- Direct-target assignment: NOT PERFORMED",
    "- Scientific figures: NONE",
    "- Active next substep: STEP03B3C3B",
    "- Full STEP03 module: NOT_LOCKED"
  ),
  file.path(
    dirs$control,
    "STATUS.md"
  )
)

log_event(
  "Completed",
  run_id,
  "Gate:",
  gate_status
)

cat("\n============================================================\n")
cat("STEP03B3C3A v4 RNA source and gene-ID policy lock completed.\n")
cat(
  "RNA source: ",
  RNA_source_lock_version,
  " LOCKED\n"
)
cat(
  "Gene-ID policy: ",
  gene_ID_policy_lock_version,
  " LOCKED\n"
)
cat("RNA identifier: feature_id\n")
cat("RNA membership: 436 IMMUTABLE\n")
cat("Exact cistrome joins: 279\n")
cat("Manual punctuation joins: 5\n")
cat("Total cistrome-authorized RNA symbols: 284\n")
cat("RNA-only unmatched symbols: 152\n")
cat("Automatic alias mappings: 0\n")
cat("Automatic punctuation mappings: 0\n")
cat(
  "Peak-gene map: ",
  peak_gene_map_lock_version,
  " IMMUTABLE; 23,499 relations\n"
)
cat("RNA/CUT&RUN integration: NOT STARTED\n")
cat("Definitive direct-target assignment: NOT PERFORMED\n")
cat("Scientific figures: NONE\n")
cat("STEP03B3C3B authorization: YES\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 5: 46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff
# ==============================================================================
# ============================================================
# STEP03B3C3B v2 - Lock RNA occupancy evidence and hand off to
# the combined biological/patient projection package
# Nutrition_Endometrial_VDR_v2
#
# Purpose:
# - perform one final content-QC pass on the completed
#   STEP03B3C3B integration candidate;
# - lock the accepted 436-row RNA/occupancy evidence table and
#   855-row relation subset;
# - create one compact downstream routing table;
# - close the VDR perturbation/cistrome core module;
# - authorize the next combined package rather than splitting it
#   into many micro-substeps.
#
# Locked interpretation:
# - 284/436 RNA interaction candidates have current author-annotated
#   VDR occupancy support;
# - 152/436 remain RNA-only candidates;
# - this is cross-context convergent evidence, not definitive direct
#   regulation or causality.
#
# This script does NOT:
# - change RNA_PROGRAM_v1.0 membership;
# - change RNA_CISTROME_GENE_ID_POLICY_v1.0;
# - change PEAK_GENE_MAP_v1.1;
# - add ATAC evidence;
# - assign definitive direct targets;
# - generate scientific figures.
#
# Run from project root:
# source(
#   "04_code/R/46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

module_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

run_id <- paste0(
  "STEP03B3C3B_v2_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

RNA_anchor <- "RNA_PROGRAM_v1.0"
RNA_source_lock_version <- "RNA_INTERACTION_SOURCE_v1.0"
gene_ID_policy_lock_version <- "RNA_CISTROME_GENE_ID_POLICY_v1.0"
peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"
integration_candidate_version <- "RNA_OCCUPANCY_EVIDENCE_CANDIDATE_v1.0"
integration_lock_version <- "RNA_OCCUPANCY_EVIDENCE_v1.0"
routing_table_version <- "DOWNSTREAM_CANDIDATE_ROUTING_v1.0"
ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_counts <- list(
  RNA_total = 436L,
  exact_join = 279L,
  manual_join = 5L,
  occupancy_supported = 284L,
  RNA_only = 152L,
  relation_join = 855L,
  vehicle_candidates = 77L,
  calcitriol_candidates = 275L,
  shared_candidates = 46L,
  promoter_candidates = 84L,
  peak_gene_relations_locked = 23499L,
  attenuation_candidates = 431L
)

expected_evidence_groups <- c(
  RNA_INTERACTION_ONLY_NO_CURRENT_CISTROME_SUPPORT = 152L,
  RNA_INTERACTION_PLUS_BOTH_CONDITION_OCCUPANCY_NO_SHARED_PEAK = 22L,
  RNA_INTERACTION_PLUS_CALCITRIOL_OCCUPANCY = 207L,
  RNA_INTERACTION_PLUS_SHARED_VDR_OCCUPANCY = 46L,
  RNA_INTERACTION_PLUS_VEHICLE_OCCUPANCY = 9L
)

expected_relation_sources <- c(
  VDR_CALCITRIOL_CISTROME = 694L,
  VDR_SHARED_CISTROME = 58L,
  VDR_VEHICLE_CISTROME = 103L
)

expected_relation_source_genes <- c(
  VDR_CALCITRIOL_CISTROME = 275L,
  VDR_SHARED_CISTROME = 46L,
  VDR_VEHICLE_CISTROME = 77L
)

required_packages <- c(
  "data.table",
  "digest",
  "jsonlite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Required package(s) missing: ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    ". This script does not install or update packages."
  )
}

dirs <- list(
  control = file.path(module_root, "00_module_control"),
  inputs = file.path(module_root, "01_inputs"),
  code = file.path(module_root, "02_code"),
  parameters = file.path(module_root, "03_parameters"),
  environment = file.path(module_root, "04_environment"),
  execution = file.path(module_root, "05_execution"),
  results = file.path(module_root, "06_results"),
  figures = file.path(module_root, "07_figures_tables"),
  qc = file.path(module_root, "08_qc"),
  text = file.path(module_root, "09_manuscript_text"),
  archive_parent = file.path(module_root, "10_archive_work")
)

required_dir_names <- c(
  "control",
  "inputs",
  "code",
  "parameters",
  "environment",
  "execution",
  "results",
  "figures",
  "qc",
  "text",
  "archive_parent"
)

missing_dir_names <- setdiff(
  required_dir_names,
  names(dirs)
)

invalid_dir_names <- names(dirs)[
  !vapply(
    dirs,
    function(path_value) {
      is.character(path_value) &&
        length(path_value) == 1L &&
        !is.na(path_value) &&
        nzchar(path_value)
    },
    logical(1)
  )
]

if (
  length(missing_dir_names) > 0L ||
    length(invalid_dir_names) > 0L
) {
  stop(
    "Invalid directory registry. Missing: ",
    paste(
      missing_dir_names,
      collapse = ", "
    ),
    "; invalid: ",
    paste(
      invalid_dir_names,
      collapse = ", "
    )
  )
}

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

candidate_root <- file.path(
  dirs$results,
  "STEP03B3C3B_RNA_occupancy_integration_candidate"
)

locked_root <- file.path(
  dirs$results,
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  integration_lock_version
)

dir.create(
  locked_root,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S3B3C3Bv2_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )
)

for (directory_path in c(
  file.path(archive_root, "res"),
  file.path(archive_root, "qc"),
  file.path(archive_root, "txt"),
  file.path(archive_root, "code"),
  file.path(archive_root, "env"),
  file.path(archive_root, "log"),
  file.path(archive_root, "ctl")
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  dirs$execution,
  paste0(
    run_id,
    "_execution.log"
  )
)

log_event <- function(...) {
  message_text <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(message_text, "\n")

  cat(
    message_text,
    "\n",
    file = log_path,
    append = TRUE
  )
}

write_text <- function(lines, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
}

write_csv_utf8 <- function(x, path) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

parse_bool <- function(x) {
  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1"
  )
}

copy_validate <- function(
    source_path,
    destination_path
) {
  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (nchar(destination_path) >= 240L) {
    stop(
      "Destination path exceeds 240 characters:\n",
      destination_path
    )
  }

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(
        file.info(source_path)$size
      ),
      as.numeric(
        file.info(destination_path)$size
      )
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed:\n",
      source_path
    )
  }

  TRUE
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Required candidate outputs
# ----------------------------------------------------------------------

required_paths <- c(
  gate = file.path(
    dirs$results,
    "STEP03B3C3B_GATE_TABLE.csv"
  ),
  candidate_manifest = file.path(
    candidate_root,
    "STEP03B3C3B_INTEGRATION_CANDIDATE_MANIFEST.csv"
  ),
  integrated_candidate = file.path(
    candidate_root,
    "STEP03B3C3B_RNA_OCCUPANCY_INTEGRATED_CANDIDATE.csv"
  ),
  relation_long = file.path(
    candidate_root,
    "STEP03B3C3B_RNA_OCCUPANCY_RELATION_LONG.tsv.gz"
  ),
  evidence_groups = file.path(
    candidate_root,
    "STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv"
  ),
  program_class = file.path(
    candidate_root,
    "STEP03B3C3B_PROGRAM_CLASS_BY_EVIDENCE_GROUP.csv"
  ),
  occupancy_source = file.path(
    candidate_root,
    "STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv"
  ),
  relation_source = file.path(
    candidate_root,
    "STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv"
  ),
  mapping_status = file.path(
    candidate_root,
    "STEP03B3C3B_MAPPING_STATUS_BY_OCCUPANCY.csv"
  ),
  integration_qc = file.path(
    dirs$qc,
    "STEP03B3C3B_INTEGRATION_QC_SUMMARY.csv"
  ),
  manual_mapping_audit = file.path(
    dirs$qc,
    "STEP03B3C3B_MANUAL_MAPPING_JOIN_AUDIT.csv"
  ),
  RNA_only_audit = file.path(
    dirs$qc,
    "STEP03B3C3B_RNA_ONLY_152_AUDIT.csv"
  )
)

missing_required_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_required_paths) > 0L) {
  stop(
    "Required STEP03B3C3B output is missing:\n",
    paste(
      missing_required_paths,
      collapse = "\n"
    )
  )
}

gate <- read.csv(
  required_paths[["gate"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_manifest <- read.csv(
  required_paths[["candidate_manifest"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

integrated <- read.csv(
  required_paths[["integrated_candidate"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

evidence_groups <- read.csv(
  required_paths[["evidence_groups"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

program_class <- read.csv(
  required_paths[["program_class"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

occupancy_source <- read.csv(
  required_paths[["occupancy_source"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

relation_source <- read.csv(
  required_paths[["relation_source"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

mapping_status <- read.csv(
  required_paths[["mapping_status"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

integration_qc <- read.csv(
  required_paths[["integration_qc"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

manual_mapping_audit <- read.csv(
  required_paths[["manual_mapping_audit"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

RNA_only_audit <- read.csv(
  required_paths[["RNA_only_audit"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

accepted_gate <- gate$status[
  gate$gate_component ==
    "STEP03B3C3B gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "PROCEED_TO_MANUAL_RNA_OCCUPANCY_INTEGRATION_LOCK"
) {
  stop(
    "STEP03B3C3B does not authorize integration lock. Gate: ",
    paste(
      accepted_gate,
      collapse = ", "
    )
  )
}

# ----------------------------------------------------------------------
# Manifest identity checks
# ----------------------------------------------------------------------

if (
  nrow(candidate_manifest) != 6L ||
    !all(
      candidate_manifest$candidate_version ==
        integration_candidate_version
    ) ||
    !all(
      candidate_manifest$RNA_membership_status ==
        "IMMUTABLE_436"
    ) ||
    !all(
      candidate_manifest$peak_gene_membership_status ==
        "IMMUTABLE_23499"
    )
) {
  stop(
    "Integration candidate manifest metadata failed."
  )
}

candidate_manifest$file_exists <- file.exists(
  candidate_manifest$path
)

candidate_manifest$observed_size_bytes <- as.numeric(
  file.info(
    candidate_manifest$path
  )$size
)

candidate_manifest$observed_md5 <- vapply(
  candidate_manifest$path,
  md5_file,
  character(1)
)

candidate_manifest$observed_sha256 <- vapply(
  candidate_manifest$path,
  sha256_file,
  character(1)
)

candidate_manifest$identity_pass <- with(
  candidate_manifest,
  file_exists &
    size_bytes ==
      observed_size_bytes &
    md5 ==
      observed_md5 &
    sha256 ==
      observed_sha256
)

if (!all(
  candidate_manifest$identity_pass
)) {
  stop(
    "One or more integration candidate artifacts failed hash validation."
  )
}

write_csv_utf8(
  candidate_manifest,
  file.path(
    dirs$qc,
    "STEP03B3C3B_v2_CANDIDATE_IDENTITY_RECHECK.csv"
  )
)

# ----------------------------------------------------------------------
# Integrated-table content QC
# ----------------------------------------------------------------------

required_integrated_columns <- c(
  "feature_id",
  "program_class",
  "interaction_logFC",
  "interaction_FDR",
  "mapping_status",
  "manual_mapping_accepted",
  "cistrome_join_authorized",
  "total_unique_peak_n",
  "total_peak_gene_relation_n",
  "vehicle_peak_n",
  "calcitriol_peak_n",
  "shared_peak_n",
  "vehicle_nonshared_peak_n",
  "calcitriol_nonshared_peak_n",
  "any_promoter_overlap_support",
  "occupancy_supported",
  "shared_occupancy_supported",
  "vehicle_occupancy_supported",
  "calcitriol_occupancy_supported",
  "occupancy_evidence_group",
  "direct_target_assignment",
  "causal_claim_authorized",
  "ATAC_evidence_status"
)

missing_integrated_columns <- setdiff(
  required_integrated_columns,
  names(integrated)
)

if (length(missing_integrated_columns) > 0L) {
  stop(
    "Integrated table is missing columns: ",
    paste(
      missing_integrated_columns,
      collapse = ", "
    )
  )
}

boolean_columns <- c(
  "manual_mapping_accepted",
  "cistrome_join_authorized",
  "any_promoter_overlap_support",
  "occupancy_supported",
  "shared_occupancy_supported",
  "vehicle_occupancy_supported",
  "calcitriol_occupancy_supported",
  "causal_claim_authorized"
)

for (column_name in boolean_columns) {
  integrated[[column_name]] <- parse_bool(
    integrated[[column_name]]
  )
}

if (
  nrow(integrated) !=
    expected_counts$RNA_total ||
    length(
      unique(
        integrated$feature_id
      )
    ) !=
      expected_counts$RNA_total ||
    any(
      is.na(
        integrated$feature_id
      ) |
        !nzchar(
          trimws(
            integrated$feature_id
          )
        )
    )
) {
  stop(
    "Integrated table failed 436-row feature_id identity QC."
  )
}

if (
  sum(
    integrated$mapping_status ==
      "EXACT_AUTHOR_SYMBOL_MATCH"
  ) !=
    expected_counts$exact_join ||
    sum(
      integrated$manual_mapping_accepted
    ) !=
      expected_counts$manual_join ||
    sum(
      integrated$cistrome_join_authorized
    ) !=
      expected_counts$occupancy_supported ||
    sum(
      integrated$occupancy_supported
    ) !=
      expected_counts$occupancy_supported ||
    sum(
      !integrated$occupancy_supported
    ) !=
      expected_counts$RNA_only
) {
  stop(
    "Integrated gene-ID or occupancy counts failed reconciliation."
  )
}

if (
  any(
    integrated$cistrome_join_authorized !=
      integrated$occupancy_supported
  ) ||
    any(
      integrated$direct_target_assignment !=
        "NOT_PERFORMED"
    ) ||
    any(
      integrated$causal_claim_authorized
    ) ||
    any(
      integrated$ATAC_evidence_status !=
        "NOT_ADDED_RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"
    )
) {
  stop(
    "Integrated interpretation safeguards failed."
  )
}

if (
  !all(
    integrated$shared_peak_n <=
      integrated$vehicle_peak_n
  ) ||
    !all(
      integrated$shared_peak_n <=
        integrated$calcitriol_peak_n
    ) ||
    !all(
      integrated$vehicle_peak_n ==
        integrated$vehicle_nonshared_peak_n +
          integrated$shared_peak_n
    ) ||
    !all(
      integrated$calcitriol_peak_n ==
        integrated$calcitriol_nonshared_peak_n +
          integrated$shared_peak_n
    )
) {
  stop(
    "Shared/nonshared peak-count reconciliation failed."
  )
}

expected_group_by_row <- ifelse(
  integrated$shared_peak_n > 0L,
  "RNA_INTERACTION_PLUS_SHARED_VDR_OCCUPANCY",
  ifelse(
    integrated$vehicle_peak_n > 0L &
      integrated$calcitriol_peak_n > 0L,
    "RNA_INTERACTION_PLUS_BOTH_CONDITION_OCCUPANCY_NO_SHARED_PEAK",
    ifelse(
      integrated$calcitriol_peak_n > 0L,
      "RNA_INTERACTION_PLUS_CALCITRIOL_OCCUPANCY",
      ifelse(
        integrated$vehicle_peak_n > 0L,
        "RNA_INTERACTION_PLUS_VEHICLE_OCCUPANCY",
        "RNA_INTERACTION_ONLY_NO_CURRENT_CISTROME_SUPPORT"
      )
    )
  )
)

if (!identical(
  as.character(
    integrated$occupancy_evidence_group
  ),
  as.character(
    expected_group_by_row
  )
)) {
  stop(
    "At least one occupancy evidence group is inconsistent with peak counts."
  )
}

# ----------------------------------------------------------------------
# Summary-table reconciliation
# ----------------------------------------------------------------------

observed_evidence_groups <- setNames(
  evidence_groups$RNA_candidate_n,
  evidence_groups$occupancy_evidence_group
)

if (
  !setequal(
    names(observed_evidence_groups),
    names(expected_evidence_groups)
  ) ||
    !all(
      as.integer(
        observed_evidence_groups[
          names(expected_evidence_groups)
        ]
      ) ==
        as.integer(
          expected_evidence_groups
        )
    ) ||
    sum(
      evidence_groups$RNA_candidate_n
    ) !=
      expected_counts$RNA_total
) {
  stop(
    "Evidence-group counts failed reconciliation."
  )
}

occupancy_value <- setNames(
  occupancy_source$RNA_candidate_n,
  occupancy_source$occupancy_metric
)

required_occupancy_values <- c(
  "RNA candidates with any occupancy support" = 284L,
  "RNA candidates with vehicle occupancy" = 77L,
  "RNA candidates with calcitriol occupancy" = 275L,
  "RNA candidates with shared occupancy" = 46L,
  "RNA candidates with promoter occupancy" = 84L,
  "RNA-only candidates without current cistrome support" = 152L,
  "Manually harmonized candidates with occupancy support" = 5L
)

if (
  !all(
    names(required_occupancy_values) %in%
      names(occupancy_value)
  ) ||
    !all(
      as.integer(
        occupancy_value[
          names(required_occupancy_values)
        ]
      ) ==
        as.integer(
          required_occupancy_values
        )
    )
) {
  stop(
    "Occupancy-source summary failed reconciliation."
  )
}

observed_relation_sources <- setNames(
  relation_source$joined_peak_gene_relation_n,
  relation_source$source_id
)

observed_relation_source_genes <- setNames(
  relation_source$joined_RNA_candidate_n,
  relation_source$source_id
)

if (
  !setequal(
    names(observed_relation_sources),
    names(expected_relation_sources)
  ) ||
    !all(
      as.integer(
        observed_relation_sources[
          names(expected_relation_sources)
        ]
      ) ==
        as.integer(
          expected_relation_sources
        )
    ) ||
    !all(
      as.integer(
        observed_relation_source_genes[
          names(expected_relation_source_genes)
        ]
      ) ==
        as.integer(
          expected_relation_source_genes
        )
    ) ||
    sum(
      observed_relation_sources
    ) !=
      expected_counts$relation_join
) {
  stop(
    "Relation-source summary failed reconciliation."
  )
}

attenuation_n <- sum(
  program_class$RNA_candidate_n[
    grepl(
      "^attenuated_ligand_",
      program_class$program_class
    )
  ]
)

if (
  sum(
    program_class$RNA_candidate_n
  ) !=
    expected_counts$RNA_total ||
    attenuation_n !=
      expected_counts$attenuation_candidates
) {
  stop(
    "Program-class by evidence-group summary failed."
  )
}

if (
  nrow(manual_mapping_audit) !=
    expected_counts$manual_join ||
    any(
      manual_mapping_audit$total_unique_peak_n <
        1L
    ) ||
    nrow(RNA_only_audit) !=
      expected_counts$RNA_only
) {
  stop(
    "Manual-mapping or RNA-only audit failed."
  )
}

qc_value <- setNames(
  integration_qc$value,
  integration_qc$metric
)

required_qc_values <- c(
  "Frozen RNA rows" = 436L,
  "Unique frozen RNA feature_id" = 436L,
  "Exact gene-ID joins" = 279L,
  "Manual punctuation joins" = 5L,
  "Cistrome-authorized RNA candidates" = 284L,
  "RNA-only candidates" = 152L,
  "Integrated candidate rows" = 436L,
  "Unique integrated feature_id" = 436L,
  "Occupancy-supported RNA candidates" = 284L,
  "RNA-only integrated rows" = 152L,
  "Relation-level joined rows" = 855L,
  "Joined unique RNA candidates" = 284L,
  "Manual mappings successfully joined" = 5L,
  "Locked peak-gene relation rows" = 23499L,
  "Peak-gene membership changed" = 0L,
  "RNA membership changed" = 0L,
  "Definitive direct-target assignments" = 0L,
  "Scientific figures" = 0L
)

if (
  !all(
    names(required_qc_values) %in%
      names(qc_value)
  ) ||
    !all(
      as.integer(
        qc_value[
          names(required_qc_values)
        ]
      ) ==
        as.integer(
          required_qc_values
        )
    )
) {
  stop(
    "Integration-QC summary differs from the reviewed values."
  )
}

# ----------------------------------------------------------------------
# Create compact downstream routing table
# ----------------------------------------------------------------------

routing_table <- integrated[
  ,
  c(
    "feature_id",
    "program_class",
    "interaction_logFC",
    "interaction_FDR",
    "RNA_gene_symbol_cistrome_join",
    "mapping_status",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support",
    "vehicle_peak_n",
    "calcitriol_peak_n",
    "shared_peak_n",
    "total_unique_peak_n",
    "occupancy_evidence_group"
  ),
  drop = FALSE
]

routing_table$downstream_route_class <- ifelse(
  routing_table$shared_occupancy_supported &
    routing_table$any_promoter_overlap_support,
  "ROUTE_1_SHARED_AND_PROMOTER_SUPPORT",
  ifelse(
    routing_table$shared_occupancy_supported,
    "ROUTE_2_SHARED_NONPROMOTER_SUPPORT",
    ifelse(
      routing_table$occupancy_supported &
        routing_table$any_promoter_overlap_support,
      "ROUTE_3_NONSHARED_PROMOTER_SUPPORT",
      ifelse(
        routing_table$occupancy_supported,
        "ROUTE_4_OTHER_OCCUPANCY_SUPPORT",
        "ROUTE_5_RNA_ONLY"
      )
    )
  )
)

routing_table$downstream_analysis_scope <- ifelse(
  routing_table$occupancy_supported,
  "PRIMARY_TEMPORAL_PATIENT_AND_SINGLE_CELL_PROJECTION",
  "SECONDARY_RNA_PROGRAM_RETENTION_AND_SENSITIVITY"
)

routing_table$routing_table_version <-
  routing_table_version

routing_table$direct_target_status <-
  "NOT_ASSIGNED"

routing_counts <- as.data.frame(
  table(
    routing_table$downstream_route_class
  ),
  stringsAsFactors = FALSE
)

names(routing_counts) <- c(
  "downstream_route_class",
  "RNA_candidate_n"
)

expected_routing_counts <- c(
  ROUTE_1_SHARED_AND_PROMOTER_SUPPORT = 21L,
  ROUTE_2_SHARED_NONPROMOTER_SUPPORT = 25L,
  ROUTE_3_NONSHARED_PROMOTER_SUPPORT = 63L,
  ROUTE_4_OTHER_OCCUPANCY_SUPPORT = 175L,
  ROUTE_5_RNA_ONLY = 152L
)

observed_routing_counts <- setNames(
  routing_counts$RNA_candidate_n,
  routing_counts$downstream_route_class
)

if (
  !setequal(
    names(observed_routing_counts),
    names(expected_routing_counts)
  ) ||
    !all(
      as.integer(
        observed_routing_counts[
          names(expected_routing_counts)
        ]
      ) ==
        as.integer(
          expected_routing_counts
        )
    )
) {
  stop(
    "Downstream routing counts failed reconciliation."
  )
}

# ----------------------------------------------------------------------
# Lock accepted artifacts
# ----------------------------------------------------------------------

locked_integrated_path <- file.path(
  locked_root,
  "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv"
)

locked_relation_path <- file.path(
  locked_root,
  "RNA_OCCUPANCY_EVIDENCE_v1.0_relation_long_855.tsv.gz"
)

locked_routing_path <- file.path(
  locked_root,
  "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv"
)

locked_evidence_counts_path <- file.path(
  locked_root,
  "RNA_OCCUPANCY_EVIDENCE_v1.0_group_counts.csv"
)

locked_routing_counts_path <- file.path(
  locked_root,
  "DOWNSTREAM_CANDIDATE_ROUTING_v1.0_counts.csv"
)

copy_validate(
  required_paths[["integrated_candidate"]],
  locked_integrated_path
)

copy_validate(
  required_paths[["relation_long"]],
  locked_relation_path
)

write_csv_utf8(
  routing_table,
  locked_routing_path
)

write_csv_utf8(
  evidence_groups,
  locked_evidence_counts_path
)

write_csv_utf8(
  routing_counts,
  locked_routing_counts_path
)

locked_paths <- c(
  locked_integrated_path,
  locked_relation_path,
  locked_routing_path,
  locked_evidence_counts_path,
  locked_routing_counts_path
)

locked_manifest <- data.frame(
  artifact_id = c(
    "RNA_OCCUPANCY_INTEGRATED_436",
    "RNA_OCCUPANCY_RELATION_LONG_855",
    "DOWNSTREAM_CANDIDATE_ROUTING_436",
    "OCCUPANCY_EVIDENCE_GROUP_COUNTS",
    "DOWNSTREAM_ROUTING_COUNTS"
  ),
  lock_version = c(
    integration_lock_version,
    integration_lock_version,
    routing_table_version,
    integration_lock_version,
    routing_table_version
  ),
  path =
    normalize_slashes(
      locked_paths
    ),
  size_bytes = as.numeric(
    file.info(
      locked_paths
    )$size
  ),
  md5 = vapply(
    locked_paths,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    locked_paths,
    sha256_file,
    character(1)
  ),
  RNA_membership_status =
    "IMMUTABLE_436",
  peak_gene_membership_status =
    "IMMUTABLE_23499",
  direct_target_status =
    "NOT_ASSIGNED",
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  locked_manifest,
  file.path(
    locked_root,
    "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv"
  )
)

lock_summary <- data.frame(
  lock_component = c(
    "RNA source",
    "RNA membership",
    "Gene-ID policy",
    "Peak-gene map",
    "Integrated evidence table",
    "Relation-level subset",
    "Occupancy-supported RNA candidates",
    "RNA-only candidates",
    "Vehicle occupancy candidates",
    "Calcitriol occupancy candidates",
    "Shared occupancy candidates",
    "Promoter occupancy candidates",
    "Direct-target assignments",
    "ATAC evidence",
    "Scientific figures"
  ),
  locked_value = c(
    RNA_source_lock_version,
    "436_IMMUTABLE",
    gene_ID_policy_lock_version,
    paste0(
      peak_gene_map_lock_version,
      "_23499_IMMUTABLE"
    ),
    paste0(
      integration_lock_version,
      "_436_ROWS"
    ),
    "855_RELATIONS",
    "284",
    "152",
    "77",
    "275",
    "46",
    "84",
    "0",
    "NOT_ADDED",
    "NONE_CREATED"
  ),
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  lock_summary,
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK.csv"
  )
)

# ----------------------------------------------------------------------
# Close STEP03 core and authorize one combined next package
# ----------------------------------------------------------------------

gate_status <-
  "READY_FOR_COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PROJECTION_PACKAGE"

gate_table <- data.frame(
  gate_component = c(
    "STEP03B3C3B candidate gate",
    "Candidate artifact hashes",
    "Integrated RNA rows",
    "Occupancy-supported RNA candidates",
    "RNA-only candidates",
    "Relation-level subset",
    "Evidence-group reconciliation",
    "Program-class reconciliation",
    "Manual punctuation mappings",
    "Shared/nonshared peak reconciliation",
    "RNA membership",
    "Peak-gene membership",
    "ATAC role",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP03 VDR perturbation/cistrome core",
    "Combined downstream package authorization",
    "STEP03B3C3B v2 gate"
  ),
  status = c(
    "PASS",
    "PASS_6_OF_6",
    "LOCKED_436",
    "LOCKED_284",
    "LOCKED_152",
    "LOCKED_855",
    "PASS_436",
    "PASS_431_ATTENUATION_OF_436",
    "PASS_5_OF_5",
    "PASS",
    "UNCHANGED_436",
    "UNCHANGED_23499",
    ATAC_role,
    "NOT_PERFORMED",
    "NONE_CREATED",
    "CORE_MODULE_LOCKED",
    "AUTHORIZED_AS_ONE_PACKAGE",
    gate_status
  ),
  required_action = c(
    "Retain the candidate audit.",
    "Use only the locked copies and hashes.",
    "Keep feature_id as RNA row authority.",
    "Use as the primary downstream projection set.",
    "Retain for secondary sensitivity analyses.",
    "Use for traceable peak-level evidence only.",
    "Do not reinterpret groups as causal tiers.",
    "Retain all four RNA interaction classes.",
    "Use only the five locked manual mappings.",
    "Keep total, shared and nonshared counts separate.",
    "Membership changes are prohibited.",
    "Membership changes are prohibited.",
    "Use selected-locus accessible context only after target selection.",
    "Use occupancy-supported candidate terminology.",
    "Begin figure workflow only after downstream biological projection.",
    "Do not reopen source/gene-ID/peak mapping without a new version.",
    "Run temporal, patient and single-cell projection as one coordinated package.",
    "STEP03B3C3B is locked and the core module is closed."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP03B3C3B v2 lock card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Integrated evidence: ",
      integration_lock_version,
      " LOCKED"
    ),
    "- Frozen RNA candidates: 436",
    "- Occupancy-supported RNA candidates: 284",
    "- RNA-only candidates: 152",
    "- Relation-level joined rows: 855",
    "- Vehicle occupancy candidates: 77",
    "- Calcitriol occupancy candidates: 275",
    "- Shared occupancy candidates: 46",
    "- Promoter occupancy candidates: 84",
    "- Attenuation-class RNA candidates: 431",
    "- Manual punctuation mappings: 5/5",
    "- RNA membership changed: NO",
    "- Peak-gene membership changed: NO",
    "- Direct-target assignments: 0",
    "- ATAC evidence: NOT ADDED",
    "- Scientific figures: NONE",
    "- STEP03 VDR perturbation/cistrome core: LOCKED",
    "- Next workflow: one combined temporal, patient and single-cell projection package",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK_CARD.md"
  )
)

write_text(
  c(
    "# Combined downstream package authorization",
    "",
    paste0(
      "- Authorized by: ",
      run_id
    ),
    paste0(
      "- Integrated evidence source: ",
      integration_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Routing table: ",
      routing_table_version,
      " LOCKED"
    ),
    "- Authorization: YES",
    "",
    "## One coordinated analysis package",
    "",
    "1. Temporal/decidualization support using GSE250130 or the accepted temporal source.",
    "2. Patient projection in GSE165004 with the already established leave-one-slide robustness framework.",
    "3. External patient validation using GSE92324 and GSE111974 under their valid sample-sharing constraints.",
    "4. Single-cell localization using GSE183837.",
    "5. Selected-locus ATAC accessible-context review only after integrated biological candidates are narrowed.",
    "",
    "## Primary and secondary scopes",
    "",
    "- Primary projection set: 284 occupancy-supported RNA interaction candidates.",
    "- Higher-context routing: 21 shared-plus-promoter and 25 shared-nonpromoter candidates.",
    "- Promoter-supported nonshared routing: 63 candidates.",
    "- Other occupancy-supported routing: 175 candidates.",
    "- Secondary sensitivity set: 152 RNA-only candidates.",
    "",
    "## Prohibited",
    "",
    "- Do not reopen the frozen RNA, gene-ID or peak-gene memberships.",
    "- Do not call the integrated candidates definitive direct targets.",
    "- Do not use discovery ROC, LASSO, random forest or generic network pharmacology.",
    "- Do not use ATAC as a genome-wide target-ranking layer.",
    "- Do not split minor file-format corrections into separate scientific substeps."
  ),
  file.path(
    dirs$control,
    "COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PACKAGE_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP03B3C3B Results closeout draft",
    "",
    "Integration of the frozen 436-gene RNA interaction program with the locked VDR peak-gene map identified 284 RNA candidates with author-annotated VDR occupancy support and retained 152 candidates without current cistrome support. The integrated relation-level subset comprised 855 peak-gene relations.",
    "",
    "Among the occupancy-supported candidates, 275 had calcitriol-context occupancy, 77 had vehicle-context occupancy and 46 had shared occupancy. Promoter-associated occupancy was present for 84 candidates. The evidence was interpreted as cross-context convergence because the RNA interaction and VDR occupancy experiments used distinct perturbation designs.",
    "",
    "No RNA candidate or peak-gene relation was added or removed, no ATAC evidence was incorporated, and no definitive direct-target assignment was made."
  ),
  file.path(
    dirs$text,
    "STEP03B3C3B_v2_Results_closeout.md"
  )
)

write_text(
  c(
    "# STEP03B3C3B Methods closeout",
    "",
    "The accepted 436-row RNA/occupancy integration table was revalidated against its manifest hashes, source-specific relation counts, gene-ID mapping states, occupancy evidence groups and promoter-support fields. Evidence groups were reconstructed directly from vehicle, calcitriol and shared peak counts, with shared coordinates separated from nonshared condition-specific coordinates.",
    "",
    "The accepted integrated table and its 855-row relation-level subset were copied into a versioned locked directory with MD5 and SHA-256 hashes. A downstream routing table was created for coordinated temporal, patient and single-cell projection; these routing classes were used for analysis organization and were not interpreted as causal or direct-target tiers."
  ),
  file.path(
    dirs$text,
    "STEP03B3C3B_v2_Methods_closeout.md"
  )
)

parameters <- list(
  schema_version = "2.0-lock",
  run_id = run_id,
  integration_lock_version =
    integration_lock_version,
  routing_table_version =
    routing_table_version,
  RNA_anchor = RNA_anchor,
  RNA_n = 436,
  occupancy_supported_n = 284,
  RNA_only_n = 152,
  relation_join_n = 855,
  vehicle_candidate_n = 77,
  calcitriol_candidate_n = 275,
  shared_candidate_n = 46,
  promoter_candidate_n = 84,
  route_1_shared_promoter_n = 21,
  route_2_shared_nonpromoter_n = 25,
  route_3_nonshared_promoter_n = 63,
  route_4_other_occupancy_n = 175,
  route_5_RNA_only_n = 152,
  RNA_membership_changed = FALSE,
  peak_gene_membership_changed = FALSE,
  ATAC_evidence_added = FALSE,
  definitive_direct_target_assignment = FALSE,
  scientific_figures_created = FALSE,
  core_module_locked = TRUE,
  combined_downstream_package_authorized = TRUE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP03B3C3B_v2_parameters.json"
  ),
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C3B_v2.txt"
  )
)


same_normalized_path <- function(path_a, path_b) {
  if (
    is.null(path_a) ||
      is.null(path_b) ||
      length(path_a) != 1L ||
      length(path_b) != 1L ||
      is.na(path_a) ||
      is.na(path_b)
  ) {
    return(FALSE)
  }

  normalized_a <- tryCatch(
    normalizePath(
      path_a,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      gsub("\\\\", "/", path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
  )
}

safe_file_copy <- function(
    from,
    to,
    overwrite = FALSE,
    recursive = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
) {
  if (same_normalized_path(from, to)) {
    return(TRUE)
  }

  base::file.copy(
    from = from,
    to = to,
    overwrite = overwrite,
    recursive = recursive,
    copy.mode = copy.mode,
    copy.date = copy.date
  )
}

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  safe_file_copy(
    source_file,
    file.path(
      project_root,
      "04_code",
      "R",
      "46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    file.path(
      dirs$code,
      "46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

# ----------------------------------------------------------------------
# Short-path archive
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK_CARD.md"
  ),
  file.path(
    dirs$qc,
    "STEP03B3C3B_v2_CANDIDATE_IDENTITY_RECHECK.csv"
  ),
  file.path(
    locked_root,
    "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv"
  ),
  locked_integrated_path,
  locked_relation_path,
  locked_routing_path,
  locked_evidence_counts_path,
  locked_routing_counts_path,
  file.path(
    dirs$control,
    "COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PACKAGE_AUTHORIZATION.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C3B_v2_Results_closeout.md"
  ),
  file.path(
    dirs$text,
    "STEP03B3C3B_v2_Methods_closeout.md"
  ),
  file.path(
    dirs$parameters,
    "STEP03B3C3B_v2_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP03B3C3B_v2.txt"
  ),
  file.path(
    dirs$code,
    "46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff.R"
  ),
  log_path
)

for (path_value in archive_artifacts[
  file.exists(
    archive_artifacts
  )
]) {
  section <- if (
    grepl(
      "\\.R$",
      basename(path_value)
    )
  ) {
    "code"
  } else if (
    grepl(
      "sessionInfo|parameters",
      basename(path_value)
    )
  ) {
    "env"
  } else if (
    grepl(
      "\\.log$",
      basename(path_value)
    )
  ) {
    "log"
  } else if (
    grepl(
      "Methods|Results|AUTHORIZATION",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "QC|RECHECK|MANIFEST|LOCK|CARD|COUNTS|ROUTING",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "qc"
  } else {
    "res"
  }

  copy_validate(
    path_value,
    file.path(
      archive_root,
      section,
      basename(path_value)
    )
  )
}

archive_files <- list.files(
  archive_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

archive_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      archive_files
    ),
    nchar(
      normalize_slashes(
        archive_root
      )
    ) + 2L
  ),
  size_bytes = as.numeric(
    file.info(
      archive_files
    )$size
  ),
  md5 = vapply(
    archive_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    archive_files,
    sha256_file,
    character(1)
  ),
  path_length = nchar(
    archive_files
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP03B3C3B_v2_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP03B3C3B v2 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Integrated evidence: ",
      integration_lock_version,
      " LOCKED"
    ),
    "- RNA candidates: 436",
    "- Occupancy-supported: 284",
    "- RNA-only: 152",
    "- Relation-level subset: 855",
    "- RNA membership changed: NO",
    "- Peak-gene membership changed: NO",
    "- Direct-target assignments: 0",
    "- Scientific figures: NONE",
    "- STEP03 core: LOCKED",
    "- Combined downstream package: AUTHORIZED",
    paste0(
      "- Gate: ",
      gate_status
    ),
    paste0(
      "- Maximum archive path length: ",
      max(
        archive_manifest$path_length
      )
    )
  ),
  file.path(
    archive_root,
    "ctl",
    "ARCHIVE_STATUS.md"
  )
)

write_text(
  c(
    "# STEP03 Module Status",
    "",
    "- Completed module: STEP03 VDR perturbation and cistrome core",
    paste0(
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- Integrated evidence: ",
      integration_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Downstream routing: ",
      routing_table_version,
      " LOCKED"
    ),
    "- RNA membership: RNA_PROGRAM_v1.0 IMMUTABLE; 436",
    "- Occupancy-supported RNA candidates: 284",
    "- RNA-only candidates: 152",
    "- Peak-gene map: PEAK_GENE_MAP_v1.1 IMMUTABLE; 23,499",
    paste0(
      "- ATAC role: ",
      ATAC_role
    ),
    "- Direct-target assignment: NOT PERFORMED",
    "- Scientific figures: NONE",
    "- Active next workflow: COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PROJECTION_PACKAGE",
    "- STEP03 core module: LOCKED"
  ),
  file.path(
    dirs$control,
    "STATUS.md"
  )
)

log_event(
  "Completed",
  run_id,
  "Gate:",
  gate_status
)

cat("\n============================================================\n")
cat("STEP03B3C3B v2 RNA occupancy evidence lock completed.\n")
cat(
  "Integrated evidence: ",
  integration_lock_version,
  " LOCKED\n"
)
cat("RNA candidates: 436\n")
cat("Occupancy-supported RNA candidates: 284\n")
cat("RNA-only candidates: 152\n")
cat("Relation-level joined rows: 855\n")
cat("Vehicle occupancy candidates: 77\n")
cat("Calcitriol occupancy candidates: 275\n")
cat("Shared occupancy candidates: 46\n")
cat("Promoter occupancy candidates: 84\n")
cat("Downstream routing classes: 21 / 25 / 63 / 175 / 152\n")
cat("RNA membership changed: NO\n")
cat("Peak-gene membership changed: NO\n")
cat("Direct-target assignments: 0\n")
cat("ATAC evidence added: NO\n")
cat("Scientific figures: NONE\n")
cat("STEP03 VDR perturbation/cistrome core: LOCKED\n")
cat("Combined downstream package: AUTHORIZED\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n")
cat("============================================================\n")
