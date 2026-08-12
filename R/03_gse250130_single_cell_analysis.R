#!/usr/bin/env Rscript
# Public analysis module for the VDR endometrial transcriptome study.
# Set VDR_PROJECT_ROOT to the project directory containing the expected public GEO inputs.
# This module preserves the executed analysis order; raw GEO data are not redistributed.


# ==============================================================================
# Analysis stage 1: 54_step04E1_v2_memory_safe_layered_normalization
# ==============================================================================
# ============================================================
# STEP04E1 v2 - Memory-safe dataset-specific single-cell integration,
# doublet assessment, preliminary broad-cell annotation,
# VDR-program localization, and donor-aware pseudobulk build
# Nutrition_Endometrial_VDR_v2
#
# Biological roles:
# - GSE250130: temporal decidualization support and LH+7 RIF context
# - GSE183837: independent donor-aware RIF localization
#
# Governance:
# - datasets are integrated separately, not pooled into one atlas;
# - sample/donor is the biological replicate;
# - no cell-level causal or diagnostic claims;
# - preliminary broad-cell labels require manual QC;
# - no scientific figures are generated in this step;
# - no definitive direct-target terminology is allowed.
#
# Memory strategy:
# - keep merged samples as separate Assay5 layers during
#   NormalizeData/FindVariableFeatures/ScaleData/RunPCA;
# - join counts/data layers only after PCA/Harmony/UMAP and
#   after scale.data has been removed;
# - persist a post-doublet merged-count checkpoint for resume.
#
# Required CRAN packages:
#   Seurat, SeuratObject, Matrix, harmony, data.table,
#   digest, jsonlite
#
# Strongly recommended Bioconductor packages:
#   scDblFinder, SingleCellExperiment, SummarizedExperiment
#
# If scDblFinder is unavailable, the script still builds the
# integrated objects but returns a HOLD Gate.
#
# Run from project root:
# source(
#   "04_code/R/54_step04E1_v2_memory_safe_layered_normalization.R",
#   encoding = "UTF-8",
#   echo = FALSE
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(3600L, getOption("timeout")))
options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260804)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04E1_",
  format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )
)

integration_version <-
  "SINGLE_CELL_INTEGRATION_v1.1_MEMORY_SAFE"

annotation_version <-
  "PRELIMINARY_BROAD_CELL_ANNOTATION_v1.0"

localization_version <-
  "SINGLE_CELL_VDR_PROGRAM_LOCALIZATION_v1.0"

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "harmony",
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
    ". Install them before rerunning. This script does not install or update packages."
  )
}

doublet_packages <- c(
  "scDblFinder",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

doublet_mode <- if (
  all(
    vapply(
      doublet_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  )
) {
  "SCDBLFINDER"
} else {
  "NOT_RUN_REQUIRED_PACKAGES_MISSING"
}

if (
  requireNamespace(
    "future",
    quietly = TRUE
  )
) {
  future::plan(
    future::sequential
  )
}

dirs <- list(
  control = file.path(
    step04_root,
    "00_module_control"
  ),
  code = file.path(
    step04_root,
    "02_code"
  ),
  parameters = file.path(
    step04_root,
    "03_parameters"
  ),
  environment = file.path(
    step04_root,
    "04_environment"
  ),
  execution = file.path(
    step04_root,
    "05_execution"
  ),
  results = file.path(
    step04_root,
    "06_results"
  ),
  figures = file.path(
    step04_root,
    "07_figures_tables"
  ),
  qc = file.path(
    step04_root,
    "08_qc"
  ),
  text = file.path(
    step04_root,
    "09_manuscript_text"
  ),
  archive_parent = file.path(
    step04_root,
    "10_archive_work"
  )
)

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(directory_path)) {
    stop(
      "Required directory could not be created: ",
      directory_path
    )
  }
}

output_root <- file.path(
  dirs$results,
  "S4E1_sc",
  "SCI_v1_0"
)

object_root <- file.path(
  output_root,
  "objects"
)

table_root <- file.path(
  output_root,
  "tables"
)

pseudobulk_root <- file.path(
  output_root,
  "pseudobulk"
)

checkpoint_root <- file.path(
  output_root,
  "checkpoints"
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4E1_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    )
  )
)

for (directory_path in c(
  output_root,
  object_root,
  table_root,
  pseudobulk_root,
  checkpoint_root,
  file.path(
    archive_root,
    "res"
  ),
  file.path(
    archive_root,
    "qc"
  ),
  file.path(
    archive_root,
    "txt"
  ),
  file.path(
    archive_root,
    "code"
  ),
  file.path(
    archive_root,
    "env"
  ),
  file.path(
    archive_root,
    "log"
  ),
  file.path(
    archive_root,
    "ctl"
  )
)) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(directory_path)) {
    stop(
      "Required output/archive directory could not be created: ",
      directory_path
    )
  }

  if (nchar(directory_path) >= 220L) {
    stop(
      "Directory exceeds the 220-character soft limit: ",
      directory_path
    )
  }
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
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    paste(
      ...,
      collapse = " "
    )
  )

  cat(
    message_text,
    "\n"
  )

  cat(
    message_text,
    "\n",
    file = log_path,
    append = TRUE
  )
}

normalize_slashes <- function(x) {
  gsub(
    "\\\\",
    "/",
    x
  )
}

write_csv_utf8 <- function(
    x,
    path
) {
  if (nchar(path) >= 220L) {
    stop(
      "Output path exceeds the 220-character soft limit: ",
      path
    )
  }

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(
    dirname(path)
  )) {
    stop(
      "Output directory could not be created: ",
      dirname(path)
    )
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

write_text <- function(
    lines,
    path
) {
  if (nchar(path) >= 220L) {
    stop(
      "Output path exceeds the 220-character soft limit: ",
      path
    )
  }

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

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(
      NA_character_
    )
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(
      NA_character_
    )
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

same_normalized_path <- function(
    path_a,
    path_b
) {
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
      normalize_slashes(
        path_a
      )
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      normalize_slashes(
        path_b
      )
    }
  )

  identical(
    tolower(
      normalized_a
    ),
    tolower(
      normalized_b
    )
  )
}

copy_if_different <- function(
    source_path,
    destination_path
) {
  if (same_normalized_path(
    source_path,
    destination_path
  )) {
    return(TRUE)
  }

  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

copy_validate <- function(
    source_path,
    destination_path
) {
  copied <- copy_if_different(
    source_path,
    destination_path
  )

  pass <- isTRUE(copied) &&
    file.exists(
      destination_path
    ) &&
    identical(
      as.numeric(
        file.info(
          source_path
        )$size
      ),
      as.numeric(
        file.info(
          destination_path
        )$size
      )
    ) &&
    identical(
      md5_file(
        source_path
      ),
      md5_file(
        destination_path
      )
    ) &&
    identical(
      sha256_file(
        source_path
      ),
      sha256_file(
        destination_path
      )
    )

  if (!pass) {
    stop(
      "Copy/hash validation failed: ",
      source_path
    )
  }

  TRUE
}

as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  toupper(
    trimws(
      as.character(x)
    )
  ) %in%
    c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
}

get_layer_matrix <- function(
    object,
    layer_name
) {
  matrix_value <- tryCatch(
    SeuratObject::LayerData(
      object,
      assay = "RNA",
      layer = layer_name
    ),
    error = function(error_object) {
      NULL
    }
  )

  if (is.null(matrix_value)) {
    matrix_value <- tryCatch(
      Seurat::GetAssayData(
        object,
        assay = "RNA",
        slot = layer_name
      ),
      error = function(error_object) {
        NULL
      }
    )
  }

  if (is.null(matrix_value)) {
    stop(
      "Could not retrieve RNA ",
      layer_name,
      " matrix."
    )
  }

  matrix_value
}

collapse_duplicate_symbols <- function(
    counts,
    gene_symbol
) {
  gene_symbol <- trimws(
    as.character(
      gene_symbol
    )
  )

  valid <- !is.na(
    gene_symbol
  ) &
    nzchar(
      gene_symbol
    ) &
    gene_symbol != "."

  counts <- counts[
    valid,
    ,
    drop = FALSE
  ]

  gene_symbol <- gene_symbol[
    valid
  ]

  unique_symbols <- unique(
    gene_symbol
  )

  aggregation_matrix <- Matrix::sparseMatrix(
    i =
      seq_along(
        gene_symbol
      ),
    j =
      match(
        gene_symbol,
        unique_symbols
      ),
    x =
      1,
    dims = c(
      length(
        gene_symbol
      ),
      length(
        unique_symbols
      )
    )
  )

  collapsed_counts <- Matrix::t(
    aggregation_matrix
  ) %*%
    counts

  collapsed_counts <- methods::as(
    collapsed_counts,
    "dgCMatrix"
  )

  rownames(
    collapsed_counts
  ) <- unique_symbols

  colnames(
    collapsed_counts
  ) <- colnames(
    counts
  )

  collapsed_counts
}

run_doublet_assessment <- function(
    counts,
    sample_accession
) {
  if (
    doublet_mode !=
      "SCDBLFINDER"
  ) {
    return(
      list(
        keep =
          rep(
            TRUE,
            ncol(
              counts
            )
          ),
        class =
          rep(
            "NOT_RUN",
            ncol(
              counts
            )
          ),
        score =
          rep(
            NA_real_,
            ncol(
              counts
            )
          ),
        method =
          "NOT_RUN_REQUIRED_PACKAGES_MISSING"
      )
    )
  }

  set.seed(
    20260804 +
      sum(
        utf8ToInt(
          sample_accession
        )
      )
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts =
        counts
    )
  )

  sce <- scDblFinder::scDblFinder(
    sce,
    verbose = FALSE
  )

  doublet_metadata <- as.data.frame(
    SummarizedExperiment::colData(
      sce
    )
  )

  class_value <- as.character(
    doublet_metadata$scDblFinder.class
  )

  score_value <- as.numeric(
    doublet_metadata$scDblFinder.score
  )

  list(
    keep =
      class_value ==
        "singlet",
    class =
      class_value,
    score =
      score_value,
    method =
      "SCDBLFINDER"
  )
}

safe_zscore <- function(x) {
  standard_deviation <- stats::sd(
    x,
    na.rm = TRUE
  )

  if (
    !is.finite(
      standard_deviation
    ) ||
      standard_deviation == 0
  ) {
    return(
      rep(
        0,
        length(
          x
        )
      )
    )
  }

  as.numeric(
    scale(
      x
    )
  )
}

calculate_module_score <- function(
    data_matrix,
    genes
) {
  present_genes <- intersect(
    genes,
    rownames(
      data_matrix
    )
  )

  if (length(present_genes) == 0L) {
    return(
      list(
        score =
          rep(
            NA_real_,
            ncol(
              data_matrix
            )
          ),
        present_genes =
          character(0)
      )
    )
  }

  list(
    score =
      as.numeric(
        Matrix::colMeans(
          data_matrix[
            present_genes,
            ,
            drop = FALSE
          ]
        )
      ),
    present_genes =
      present_genes
  )
}

calculate_signed_program_score <- function(
    data_matrix,
    gene_table
) {
  positive_genes <- gene_table$feature_id[
    gene_table$interaction_logFC > 0
  ]

  negative_genes <- gene_table$feature_id[
    gene_table$interaction_logFC < 0
  ]

  positive_present <- intersect(
    positive_genes,
    rownames(
      data_matrix
    )
  )

  negative_present <- intersect(
    negative_genes,
    rownames(
      data_matrix
    )
  )

  positive_score <- if (
    length(
      positive_present
    ) > 0L
  ) {
    as.numeric(
      Matrix::colMeans(
        data_matrix[
          positive_present,
          ,
          drop = FALSE
        ]
      )
    )
  } else {
    rep(
      0,
      ncol(
        data_matrix
      )
    )
  }

  negative_score <- if (
    length(
      negative_present
    ) > 0L
  ) {
    as.numeric(
      Matrix::colMeans(
        data_matrix[
          negative_present,
          ,
          drop = FALSE
        ]
      )
    )
  } else {
    rep(
      0,
      ncol(
        data_matrix
      )
    )
  }

  list(
    score =
      positive_score -
      negative_score,
    positive_present_n =
      length(
        positive_present
      ),
    negative_present_n =
      length(
        negative_present
      ),
    total_present_n =
      length(
        union(
          positive_present,
          negative_present
        )
      )
  )
}

build_sparse_pseudobulk <- function(
    counts,
    metadata,
    minimum_cells = 50L
) {
  group_table <- unique(
    metadata[
      ,
      c(
        "sample_accession",
        "preliminary_broad_cell_type"
      ),
      drop = FALSE
    ]
  )

  group_table$group_id <- paste(
    group_table$sample_accession,
    group_table$preliminary_broad_cell_type,
    sep = "__"
  )

  cell_group_id <- paste(
    metadata$sample_accession,
    metadata$preliminary_broad_cell_type,
    sep = "__"
  )

  group_levels <- unique(
    cell_group_id
  )

  group_factor <- factor(
    cell_group_id,
    levels =
      group_levels
  )

  membership <- Matrix::sparseMatrix(
    i =
      seq_along(
        group_factor
      ),
    j =
      as.integer(
        group_factor
      ),
    x =
      1,
    dims = c(
      length(
        group_factor
      ),
      length(
        group_levels
      )
    ),
    dimnames = list(
      colnames(
        counts
      ),
      group_levels
    )
  )

  cell_n <- as.numeric(
    Matrix::colSums(
      membership
    )
  )

  keep_group <- cell_n >=
    minimum_cells

  pseudobulk_counts <- counts %*%
    membership[
      ,
      keep_group,
      drop = FALSE
    ]

  pseudobulk_counts <- methods::as(
    pseudobulk_counts,
    "dgCMatrix"
  )

  group_metadata <- data.frame(
    group_id =
      group_levels,
    cell_n =
      cell_n,
    stringsAsFactors = FALSE
  )

  group_metadata <- merge(
    group_metadata,
    group_table,
    by =
      "group_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  group_metadata$included_in_pseudobulk <-
    group_metadata$cell_n >=
      minimum_cells

  list(
    counts =
      pseudobulk_counts,
    metadata =
      group_metadata
  )
}

log_event(
  "Starting",
  run_id,
  "Doublet mode:",
  doublet_mode
)

# ----------------------------------------------------------------------
# Validate locked handoff
# ----------------------------------------------------------------------

qc1_root <- file.path(
  dirs$results,
  "S4D_qc_lock",
  "BVP_v1_3"
)

qc1_gate_path <- file.path(
  qc1_root,
  "STEP04D_QC1_GATE.csv"
)

candidate_path <- file.path(
  qc1_root,
  "cross_context_candidates_v1_3.csv"
)

watchlist_path <- file.path(
  qc1_root,
  "sc_sample_watchlist.csv"
)

preQC_root <- file.path(
  dirs$results,
  "S4D_sc_preQC",
  "SCP_v1_0"
)

design_root <- file.path(
  dirs$results,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0"
)

required_inputs <- c(
  qc1_gate_path,
  candidate_path,
  watchlist_path,
  file.path(
    design_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"
  ),
  file.path(
    design_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"
  )
)

missing_inputs <- required_inputs[
  !file.exists(
    required_inputs
  )
]

if (length(missing_inputs) > 0L) {
  stop(
    "Required STEP04E1 input is missing:\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

qc1_gate <- read.csv(
  qc1_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

qc1_status <- qc1_gate$status[
  qc1_gate$gate_component ==
    "STEP04D-QC1 gate"
]

if (
  length(qc1_status) != 1L ||
    qc1_status !=
      "READY_FOR_SINGLE_CELL_INTEGRATION_WITH_GSE183837_SAMPLE_WATCHLIST"
) {
  stop(
    "STEP04D-QC1 did not authorize single-cell integration."
  )
}

candidate_table <- read.csv(
  candidate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_table$occupancy_supported <-
  as_logical_flag(
    candidate_table$occupancy_supported
  )

candidate_table$shared_occupancy_supported <-
  as_logical_flag(
    candidate_table$shared_occupancy_supported
  )

candidate_table$any_promoter_overlap_support <-
  as_logical_flag(
    candidate_table$any_promoter_overlap_support
  )

if (
  nrow(candidate_table) != 436L ||
    length(
      unique(
        candidate_table$feature_id
      )
    ) !=
      436L
) {
  stop(
    "The corrected 436-candidate table failed validation."
  )
}

watchlist <- read.csv(
  watchlist_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

read_locked_design <- function(
    dataset_id
) {
  read.csv(
    file.path(
      design_root,
      paste0(
        "COMBINED_SAMPLE_DESIGN_v1.0_",
        dataset_id,
        ".csv"
      )
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

designs <- list(
  GSE250130 =
    read_locked_design(
      "GSE250130"
    ),
  GSE183837 =
    read_locked_design(
      "GSE183837"
    )
)

expected_sample_n <- c(
  GSE250130 = 28L,
  GSE183837 = 9L
)

for (dataset_id in names(
  designs
)) {
  if (
    nrow(
      designs[[
        dataset_id
      ]]
    ) !=
      expected_sample_n[[
        dataset_id
      ]] ||
      anyDuplicated(
        designs[[
          dataset_id
        ]]$sample_accession
      )
  ) {
    stop(
      dataset_id,
      " locked single-cell design failed sample-count or uniqueness validation."
    )
  }
}

# ----------------------------------------------------------------------
# Marker sets and program definitions
# ----------------------------------------------------------------------

broad_marker_sets <- list(
  Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "CLDN4",
    "MUC1"
  ),
  Stromal_fibroblast = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "COL3A1",
    "LUM",
    "COL6A1",
    "C7"
  ),
  Endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "EMCN",
    "ESAM",
    "RAMP2"
  ),
  Immune = c(
    "PTPRC",
    "TYROBP",
    "LST1",
    "FCER1G",
    "HLA-DRA",
    "CTSS"
  ),
  Pericyte_SMC = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "ACTA2",
    "TAGLN",
    "MYL9",
    "COL4A1"
  ),
  Erythroid = c(
    "HBB",
    "HBA1",
    "HBA2",
    "ALAS2",
    "GYPA"
  )
)

cycling_markers <- c(
  "MKI67",
  "TOP2A",
  "UBE2C",
  "TYMS",
  "PCNA",
  "CENPF"
)

program_definitions <- list(
  RNA_ALL_436 =
    candidate_table,
  OCCUPANCY_SUPPORTED_284 =
    candidate_table[
      candidate_table$occupancy_supported,
      ,
      drop = FALSE
    ],
  SHARED_OCCUPANCY_46 =
    candidate_table[
      candidate_table$shared_occupancy_supported,
      ,
      drop = FALSE
    ],
  PROMOTER_OCCUPANCY_84 =
    candidate_table[
      candidate_table$any_promoter_overlap_support,
      ,
      drop = FALSE
    ]
)

# ----------------------------------------------------------------------
# Dataset processing function
# ----------------------------------------------------------------------

process_dataset <- function(
    dataset_id,
    sample_design
) {
  log_event(
    "Processing dataset",
    dataset_id
  )

  dataset_object_root <- file.path(
    object_root,
    dataset_id
  )

  dataset_table_root <- file.path(
    table_root,
    dataset_id
  )

  dataset_pseudobulk_root <- file.path(
    pseudobulk_root,
    dataset_id
  )

  dataset_checkpoint_root <- file.path(
    checkpoint_root,
    dataset_id
  )

  for (directory_path in c(
    dataset_object_root,
    dataset_table_root,
    dataset_pseudobulk_root,
    dataset_checkpoint_root
  )) {
    dir.create(
      directory_path,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (!dir.exists(directory_path)) {
      stop(
        "Could not create dataset output directory: ",
        directory_path
      )
    }
  }

  integrated_checkpoint_path <- file.path(
    dataset_checkpoint_root,
    paste0(
      dataset_id,
      "_integrated_checkpoint.rds"
    )
  )

  checkpoint_status_path <- file.path(
    dataset_checkpoint_root,
    paste0(
      dataset_id,
      "_checkpoint_status.csv"
    )
  )

  merged_counts_checkpoint_path <- file.path(
    dataset_checkpoint_root,
    paste0(
      dataset_id,
      "_post_doublet_merged_counts_checkpoint.rds"
    )
  )

  merged_counts_status_path <- file.path(
    dataset_checkpoint_root,
    paste0(
      dataset_id,
      "_post_doublet_merged_counts_status.csv"
    )
  )

  checkpoint_reusable <- FALSE

  if (
    file.exists(
      integrated_checkpoint_path
    ) &&
      file.exists(
        checkpoint_status_path
      )
  ) {
    checkpoint_status <- read.csv(
      checkpoint_status_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    checkpoint_reusable <-
      nrow(checkpoint_status) == 1L &&
      checkpoint_status$dataset_id[[1]] ==
        dataset_id &&
      checkpoint_status$doublet_mode[[1]] ==
        doublet_mode &&
      checkpoint_status$sample_n[[1]] ==
        nrow(sample_design) &&
      checkpoint_status$checkpoint_status[[1]] ==
        "COMPLETE"
  }

  if (checkpoint_reusable) {
    log_event(
      "Reusing integrated checkpoint for",
      dataset_id
    )

    object <- readRDS(
      integrated_checkpoint_path
    )

    sample_doublet_qc <- read.csv(
      file.path(
        dataset_table_root,
        paste0(
          dataset_id,
          "_sample_doublet_qc.csv"
        )
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    merged_counts_reusable <- FALSE

    if (
      file.exists(
        merged_counts_checkpoint_path
      ) &&
        file.exists(
          merged_counts_status_path
        )
    ) {
      merged_counts_status <- read.csv(
        merged_counts_status_path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      merged_counts_reusable <-
        nrow(merged_counts_status) == 1L &&
        merged_counts_status$dataset_id[[1]] ==
          dataset_id &&
        merged_counts_status$doublet_mode[[1]] ==
          doublet_mode &&
        merged_counts_status$sample_n[[1]] ==
          nrow(sample_design) &&
        merged_counts_status$checkpoint_status[[1]] ==
          "COMPLETE"
    }

    if (merged_counts_reusable) {
      log_event(
        "Reusing post-doublet merged-count checkpoint for",
        dataset_id
      )

      object <- readRDS(
        merged_counts_checkpoint_path
      )

      sample_doublet_qc <- read.csv(
        file.path(
          dataset_table_root,
          paste0(
            dataset_id,
            "_sample_doublet_qc.csv"
          )
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      sample_objects <- list()
      doublet_qc_rows <- list()

    for (sample_index in seq_len(
      nrow(
        sample_design
      )
    )) {
      sample_accession <-
        sample_design$sample_accession[
          sample_index
        ]

      log_event(
        dataset_id,
        "sample",
        sample_index,
        "/",
        nrow(sample_design),
        sample_accession
      )

      preQC_path <- file.path(
        preQC_root,
        dataset_id,
        paste0(
          sample_accession,
          "_preQC.rds"
        )
      )

      if (!file.exists(preQC_path)) {
        stop(
          "Pre-QC object missing: ",
          preQC_path
        )
      }

      preQC_object <- readRDS(
        preQC_path
      )

      if (
        !all(
          c(
            "counts",
            "feature_table",
            "cell_metadata"
          ) %in%
            names(
              preQC_object
            )
        )
      ) {
        stop(
          dataset_id,
          " ",
          sample_accession,
          " pre-QC object lacks required components."
        )
      }

      counts <- preQC_object$counts

      if (
        nrow(
          preQC_object$feature_table
        ) !=
          nrow(counts)
      ) {
        stop(
          dataset_id,
          " ",
          sample_accession,
          " feature table/count matrix row mismatch."
        )
      }

      counts <- collapse_duplicate_symbols(
        counts =
          counts,
        gene_symbol =
          preQC_object$feature_table$gene_symbol
      )

      doublet_result <- run_doublet_assessment(
        counts =
          counts,
        sample_accession =
          sample_accession
      )

      if (
        length(
          doublet_result$keep
        ) !=
          ncol(counts)
      ) {
        stop(
          dataset_id,
          " ",
          sample_accession,
          " doublet-result length mismatch."
        )
      }

      retained_counts <- counts[
        ,
        doublet_result$keep,
        drop = FALSE
      ]

      retained_metadata <- data.frame(
        row.names =
          colnames(
            retained_counts
          ),
        sample_accession =
          sample_accession,
        dataset_id =
          dataset_id,
        doublet_method =
          doublet_result$method,
        doublet_class =
          doublet_result$class[
            doublet_result$keep
          ],
        doublet_score =
          doublet_result$score[
            doublet_result$keep
          ],
        stringsAsFactors = FALSE
      )

      design_row <- sample_design[
        sample_design$sample_accession ==
          sample_accession,
        ,
        drop = FALSE
      ]

      for (column_name in setdiff(
        names(
          design_row
        ),
        c(
          "sample_accession",
          "dataset_id"
        )
      )) {
        retained_metadata[[
          column_name
        ]] <-
          design_row[[
            column_name
          ]][[1]]
      }

      watch_row <- watchlist[
        watchlist$dataset_id ==
          dataset_id &
          watchlist$sample_accession ==
            sample_accession,
        ,
        drop = FALSE
      ]

      retained_metadata$preQC_watch_tier <-
        if (
          nrow(
            watch_row
          ) == 1L
        ) {
          watch_row$QC_review_tier[[1]]
        } else {
          "PASS"
        }

      sample_object <- SeuratObject::CreateSeuratObject(
        counts =
          retained_counts,
        project =
          dataset_id,
        min.cells =
          0,
        min.features =
          0,
        meta.data =
          retained_metadata
      )

      sample_objects[[
        sample_accession
      ]] <-
        sample_object

      doublet_qc_rows[[
        sample_accession
      ]] <- data.frame(
        dataset_id =
          dataset_id,
        sample_accession =
          sample_accession,
        pre_doublet_cell_n =
          ncol(counts),
        post_doublet_cell_n =
          ncol(retained_counts),
        predicted_doublet_n =
          sum(
            !doublet_result$keep
          ),
        predicted_doublet_fraction =
          mean(
            !doublet_result$keep
          ),
        doublet_method =
          doublet_result$method,
        preQC_watch_tier =
          retained_metadata$preQC_watch_tier[[1]],
        stringsAsFactors = FALSE
      )

      rm(
        preQC_object,
        counts,
        retained_counts,
        sample_object
      )

      gc()
    }

    sample_doublet_qc <- do.call(
      rbind,
      doublet_qc_rows
    )

    write_csv_utf8(
      sample_doublet_qc,
      file.path(
        dataset_table_root,
        paste0(
          dataset_id,
          "_sample_doublet_qc.csv"
        )
      )
    )

    log_event(
      "Merging",
      dataset_id,
      "sample objects"
    )

    first_sample <- names(
      sample_objects
    )[[1]]

    object <- merge(
      x =
        sample_objects[[
          first_sample
        ]],
      y =
        sample_objects[
          setdiff(
            names(
              sample_objects
            ),
            first_sample
          )
        ],
      add.cell.ids =
        NULL,
      merge.data =
        FALSE,
      project =
        dataset_id
    )

      rm(
        sample_objects
      )

      gc()

      saveRDS(
        object,
        merged_counts_checkpoint_path,
        compress =
          FALSE
      )

      write_csv_utf8(
        data.frame(
          dataset_id =
            dataset_id,
          doublet_mode =
            doublet_mode,
          sample_n =
            nrow(sample_design),
          cell_n =
            ncol(object),
          checkpoint_status =
            "COMPLETE",
          checkpoint_path =
            normalize_slashes(
              merged_counts_checkpoint_path
            ),
          stringsAsFactors = FALSE
        ),
        merged_counts_status_path
      )

      log_event(
        "Saved post-doublet merged-count checkpoint for",
        dataset_id
      )
    }

    log_event(
      "Keeping sample-specific RNA layers through normalization for",
      dataset_id
    )

    object <- Seurat::NormalizeData(
      object,
      assay =
        "RNA",
      normalization.method =
        "LogNormalize",
      scale.factor =
        10000,
      block.size =
        2000,
      verbose =
        TRUE
    )

    object <- Seurat::FindVariableFeatures(
      object,
      assay =
        "RNA",
      selection.method =
        "vst",
      nfeatures =
        3000,
      verbose =
        TRUE
    )

    object <- Seurat::ScaleData(
      object,
      assay =
        "RNA",
      features =
        SeuratObject::VariableFeatures(
          object
        ),
      block.size =
        1000,
      verbose =
        TRUE
    )

    object <- Seurat::RunPCA(
      object,
      assay =
        "RNA",
      features =
        SeuratObject::VariableFeatures(
          object
        ),
      npcs =
        50,
      seed.use =
        20260804,
      verbose =
        TRUE
    )

    object <- harmony::RunHarmony(
      object =
        object,
      group.by.vars =
        "sample_accession",
      reduction =
        "pca",
      dims.use =
        1:40,
      reduction.save =
        "harmony",
      plot_convergence =
        FALSE,
      verbose =
        TRUE
    )

    object <- Seurat::FindNeighbors(
      object,
      reduction =
        "harmony",
      dims =
        1:30,
      k.param =
        30,
      verbose =
        TRUE
    )

    object <- Seurat::FindClusters(
      object,
      resolution =
        0.2,
      algorithm =
        1,
      random.seed =
        20260804,
      verbose =
        TRUE
    )

    object <- Seurat::FindClusters(
      object,
      resolution =
        0.6,
      algorithm =
        1,
      random.seed =
        20260804,
      verbose =
        TRUE
    )

    object <- Seurat::FindClusters(
      object,
      resolution =
        0.4,
      algorithm =
        1,
      random.seed =
        20260804,
      verbose =
        TRUE
    )

    object <- Seurat::RunUMAP(
      object,
      reduction =
        "harmony",
      dims =
        1:30,
      reduction.name =
        "umap_harmony",
      reduction.key =
        "UMAPH_",
      seed.use =
        20260804,
      verbose =
        TRUE
    )

    reductions_to_keep <- intersect(
      c(
        "pca",
        "harmony",
        "umap_harmony"
      ),
      SeuratObject::Reductions(
        object
      )
    )

    graphs_to_keep <- SeuratObject::Graphs(
      object
    )

    RNA_layers_before_join <- SeuratObject::Layers(
      object[["RNA"]]
    )

    RNA_layers_to_keep <- grep(
      "^(counts|data)(\\.|$)",
      RNA_layers_before_join,
      value =
        TRUE
    )

    if (length(RNA_layers_to_keep) == 0L) {
      stop(
        dataset_id,
        " has no counts/data RNA layers before memory cleanup."
      )
    }

    object <- Seurat::DietSeurat(
      object,
      assays =
        "RNA",
      layers =
        list(
          RNA =
            RNA_layers_to_keep
        ),
      dimreducs =
        reductions_to_keep,
      graphs =
        graphs_to_keep,
      misc =
        TRUE
    )

    gc()

    if (
      "JoinLayers" %in%
        getNamespaceExports(
          "SeuratObject"
        )
    ) {
      log_event(
        "Joining RNA layers after PCA/Harmony/UMAP for",
        dataset_id
      )

      object <- SeuratObject::JoinLayers(
        object,
        assay =
          "RNA"
      )

      gc()
    }

    object <- Seurat::DietSeurat(
      object,
      assays =
        "RNA",
      layers =
        c(
          "counts",
          "data"
        ),
      dimreducs =
        reductions_to_keep,
      graphs =
        graphs_to_keep,
      misc =
        TRUE
    )

    saveRDS(
      object,
      integrated_checkpoint_path,
      compress =
        FALSE
    )

    write_csv_utf8(
      data.frame(
        dataset_id =
          dataset_id,
        doublet_mode =
          doublet_mode,
        sample_n =
          nrow(sample_design),
        cell_n =
          ncol(object),
        checkpoint_status =
          "COMPLETE",
        checkpoint_path =
          normalize_slashes(
            integrated_checkpoint_path
          ),
        stringsAsFactors = FALSE
      ),
      checkpoint_status_path
    )
  }

  if (
    !all(
      sample_design$sample_accession %in%
        unique(
          object$sample_accession
        )
    )
  ) {
    stop(
      dataset_id,
      " integrated object does not retain all locked samples."
    )
  }

  data_matrix <- get_layer_matrix(
    object,
    "data"
  )

  # ------------------------------------------------------------
  # Preliminary broad-cell annotation
  # ------------------------------------------------------------

  module_score_columns <- character(0)
  module_coverage_rows <- list()

  for (marker_name in names(
    broad_marker_sets
  )) {
    score_value <- calculate_module_score(
      data_matrix =
        data_matrix,
      genes =
        broad_marker_sets[[
          marker_name
        ]]
    )

    column_name <- paste0(
      "lineage_score_",
      marker_name
    )

    object[[
      column_name
    ]] <-
      safe_zscore(
        score_value$score
      )

    module_score_columns <- c(
      module_score_columns,
      column_name
    )

    module_coverage_rows[[
      marker_name
    ]] <- data.frame(
      dataset_id =
        dataset_id,
      module_name =
        marker_name,
      expected_gene_n =
        length(
          broad_marker_sets[[
            marker_name
          ]]
        ),
      present_gene_n =
        length(
          score_value$present_genes
        ),
      present_genes =
        paste(
          score_value$present_genes,
          collapse = ";"
        ),
      stringsAsFactors = FALSE
    )
  }

  cycling_score <- calculate_module_score(
    data_matrix =
      data_matrix,
    genes =
      cycling_markers
  )

  object$cycling_score <-
    safe_zscore(
      cycling_score$score
    )

  cluster_column <- "seurat_clusters"

  metadata_dt <- data.table::as.data.table(
    object@meta.data,
    keep.rownames =
      "cell_barcode"
  )

  cluster_score_dt <- metadata_dt[
    ,
    lapply(
      .SD,
      mean,
      na.rm = TRUE
    ),
    by =
      cluster_column,
    .SDcols =
      module_score_columns
  ]

  score_matrix <- as.matrix(
    cluster_score_dt[
      ,
      module_score_columns,
      with = FALSE
    ]
  )

  top_index <- apply(
    score_matrix,
    1L,
    which.max
  )

  sorted_scores <- t(
    apply(
      score_matrix,
      1L,
      sort,
      decreasing = TRUE
    )
  )

  top_score <- sorted_scores[
    ,
    1L
  ]

  second_score <- if (
    ncol(
      sorted_scores
    ) >= 2L
  ) {
    sorted_scores[
      ,
      2L
    ]
  } else {
    rep(
      NA_real_,
      nrow(
        sorted_scores
      )
    )
  }

  broad_labels <- names(
    broad_marker_sets
  )[
    top_index
  ]

  label_margin <- top_score -
    second_score

  preliminary_label <- ifelse(
    !is.finite(
      top_score
    ) |
      top_score < 0 |
      !is.finite(
        label_margin
      ) |
      label_margin < 0.15,
    "AMBIGUOUS",
    broad_labels
  )

  cluster_annotation <- data.frame(
    dataset_id =
      dataset_id,
    seurat_cluster =
      as.character(
        cluster_score_dt[[
          cluster_column
        ]]
      ),
    preliminary_broad_cell_type =
      preliminary_label,
    top_lineage_score =
      top_score,
    second_lineage_score =
      second_score,
    lineage_score_margin =
      label_margin,
    stringsAsFactors = FALSE
  )

  cluster_annotation <- cbind(
    cluster_annotation,
    as.data.frame(
      score_matrix,
      stringsAsFactors = FALSE
    )
  )

  label_map <- setNames(
    cluster_annotation$preliminary_broad_cell_type,
    cluster_annotation$seurat_cluster
  )

  object$preliminary_broad_cell_type <-
    unname(
      label_map[
        as.character(
          object$seurat_clusters
        )
      ]
    )

  object$cycling_state <- ifelse(
    object$cycling_score >= 1,
    "CYCLING_HIGH",
    "CYCLING_LOW"
  )

  write_csv_utf8(
    do.call(
      rbind,
      module_coverage_rows
    ),
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_broad_marker_coverage.csv"
      )
    )
  )

  write_csv_utf8(
    cluster_annotation,
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_preliminary_cluster_annotation.csv"
      )
    )
  )

  # ------------------------------------------------------------
  # Cluster markers on controlled downsample
  # ------------------------------------------------------------

  SeuratObject::Idents(
    object
  ) <- "seurat_clusters"

  marker_results <- tryCatch(
    Seurat::FindAllMarkers(
      object,
      assay =
        "RNA",
      only.pos =
        TRUE,
      min.pct =
        0.25,
      logfc.threshold =
        0.25,
      test.use =
        "wilcox",
      max.cells.per.ident =
        1500,
      random.seed =
        20260804,
      verbose =
        TRUE
    ),
    error = function(error_object) {
      log_event(
        dataset_id,
        "FindAllMarkers failed:",
        conditionMessage(
          error_object
        )
      )

      data.frame(
        marker_status =
          paste0(
            "FAILED_",
            conditionMessage(
              error_object
            )
          ),
        stringsAsFactors = FALSE
      )
    }
  )

  write_csv_utf8(
    marker_results,
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_cluster_markers.csv"
      )
    )
  )

  if (
    "cluster" %in%
      names(
        marker_results
      )
  ) {
    fold_change_column <- intersect(
      c(
        "avg_log2FC",
        "avg_logFC"
      ),
      names(
        marker_results
      )
    )

    if (length(
      fold_change_column
    ) > 0L) {
      marker_dt <- data.table::as.data.table(
        marker_results
      )

      data.table::setorderv(
        marker_dt,
        cols = c(
          "cluster",
          fold_change_column[[1]]
        ),
        order = c(
          1L,
          -1L
        )
      )

      top_markers <- marker_dt[
        ,
        head(
          .SD,
          25L
        ),
        by =
          cluster
      ]

      write_csv_utf8(
        as.data.frame(
          top_markers
        ),
        file.path(
          dataset_table_root,
          paste0(
            dataset_id,
            "_top25_cluster_markers.csv"
          )
        )
      )
    }
  }

  # ------------------------------------------------------------
  # VDR-program localization
  # ------------------------------------------------------------

  program_coverage_rows <- list()

  for (program_name in names(
    program_definitions
  )) {
    score_value <-
      calculate_signed_program_score(
        data_matrix =
          data_matrix,
        gene_table =
          program_definitions[[
            program_name
          ]]
      )

    score_column <- paste0(
      "VDR_program_",
      program_name
    )

    object[[
      score_column
    ]] <-
      score_value$score

    program_coverage_rows[[
      program_name
    ]] <- data.frame(
      dataset_id =
        dataset_id,
      program_name =
        program_name,
      frozen_gene_n =
        nrow(
          program_definitions[[
            program_name
          ]]
        ),
      positive_present_n =
        score_value$positive_present_n,
      negative_present_n =
        score_value$negative_present_n,
      total_present_n =
        score_value$total_present_n,
      stringsAsFactors = FALSE
    )
  }

  write_csv_utf8(
    do.call(
      rbind,
      program_coverage_rows
    ),
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_VDR_program_coverage.csv"
      )
    )
  )

  # ------------------------------------------------------------
  # Donor, cluster, composition and program audits
  # ------------------------------------------------------------

  metadata_dt <- data.table::as.data.table(
    object@meta.data,
    keep.rownames =
      "cell_barcode"
  )

  cluster_donor_counts <- metadata_dt[
    ,
    .(
      cell_n =
        .N
    ),
    by = .(
      seurat_clusters,
      sample_accession
    )
  ]

  cluster_totals <- cluster_donor_counts[
    ,
    .(
      cluster_cell_n =
        sum(
          cell_n
        ),
      donor_n =
        data.table::uniqueN(
          sample_accession
        ),
      max_donor_fraction =
        max(
          cell_n
        ) /
        sum(
          cell_n
        )
    ),
    by =
      seurat_clusters
  ]

  write_csv_utf8(
    as.data.frame(
      cluster_donor_counts
    ),
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_cluster_by_donor_counts.csv"
      )
    )
  )

  write_csv_utf8(
    as.data.frame(
      cluster_totals
    ),
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_cluster_donor_mixing_audit.csv"
      )
    )
  )

  broad_composition <- metadata_dt[
    ,
    .(
      cell_n =
        .N
    ),
    by = .(
      sample_accession,
      preliminary_broad_cell_type
    )
  ]

  broad_composition[
    ,
    sample_total :=
      sum(
        cell_n
      ),
    by =
      sample_accession
  ]

  broad_composition[
    ,
    cell_fraction :=
      cell_n /
      sample_total
  ]

  write_csv_utf8(
    as.data.frame(
      broad_composition
    ),
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_broad_cell_composition.csv"
      )
    )
  )

  program_columns <- paste0(
    "VDR_program_",
    names(
      program_definitions
    )
  )

  donor_program_summary <- metadata_dt[
    ,
    c(
      list(
        cell_n =
          .N
      ),
      lapply(
        .SD,
        function(x) {
          mean(
            x,
            na.rm = TRUE
          )
        }
      )
    ),
    by = .(
      sample_accession,
      preliminary_broad_cell_type
    ),
    .SDcols =
      program_columns
  ]

  donor_program_summary <- merge(
    as.data.frame(
      donor_program_summary
    ),
    sample_design,
    by =
      "sample_accession",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  write_csv_utf8(
    donor_program_summary,
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_donor_celltype_VDR_program_scores.csv"
      )
    )
  )

  # ------------------------------------------------------------
  # Preliminary donor-aware pseudobulk matrices
  # ------------------------------------------------------------

  counts_matrix <- get_layer_matrix(
    object,
    "counts"
  )

  pseudobulk <- build_sparse_pseudobulk(
    counts =
      counts_matrix,
    metadata =
      object@meta.data,
    minimum_cells =
      50L
  )

  pseudobulk_matrix_path <- file.path(
    dataset_pseudobulk_root,
    paste0(
      dataset_id,
      "_broad_cell_pseudobulk_counts.rds"
    )
  )

  saveRDS(
    pseudobulk$counts,
    pseudobulk_matrix_path,
    compress =
      TRUE
  )

  pseudobulk_metadata <- merge(
    pseudobulk$metadata,
    sample_design,
    by =
      "sample_accession",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  write_csv_utf8(
    pseudobulk_metadata,
    file.path(
      dataset_pseudobulk_root,
      paste0(
        dataset_id,
        "_broad_cell_pseudobulk_metadata.csv"
      )
    )
  )

  # ------------------------------------------------------------
  # Final object and dataset QC
  # ------------------------------------------------------------

  object_path <- file.path(
    dataset_object_root,
    paste0(
      dataset_id,
      "_integrated_preliminary_annotation.rds"
    )
  )

  reductions_to_keep <- intersect(
    c(
      "pca",
      "harmony",
      "umap_harmony"
    ),
    SeuratObject::Reductions(
      object
    )
  )

  graphs_to_keep <- SeuratObject::Graphs(
    object
  )

  object <- Seurat::DietSeurat(
    object,
    assays =
      "RNA",
    layers =
      c(
        "counts",
        "data"
      ),
    dimreducs =
      reductions_to_keep,
    graphs =
      graphs_to_keep,
    misc =
      TRUE
  )

  saveRDS(
    object,
    object_path,
    compress =
      FALSE
  )

  if (
    file.exists(
      integrated_checkpoint_path
    ) &&
      !same_normalized_path(
        integrated_checkpoint_path,
        object_path
      )
  ) {
    unlink(
      integrated_checkpoint_path,
      force =
        TRUE
    )
  }

  ambiguous_fraction <- mean(
    object$preliminary_broad_cell_type ==
      "AMBIGUOUS"
  )

  dataset_qc <- data.frame(
    dataset_id =
      dataset_id,
    locked_sample_n =
      nrow(
        sample_design
      ),
    retained_sample_n =
      data.table::uniqueN(
        object$sample_accession
      ),
    retained_cell_n =
      ncol(
        object
      ),
    doublet_mode =
      doublet_mode,
    predicted_doublet_n =
      sum(
        sample_doublet_qc$predicted_doublet_n
      ),
    predicted_doublet_fraction =
      sum(
        sample_doublet_qc$predicted_doublet_n
      ) /
      sum(
        sample_doublet_qc$pre_doublet_cell_n
      ),
    cluster_n =
      length(
        unique(
          object$seurat_clusters
        )
      ),
    preliminary_broad_cell_type_n =
      length(
        unique(
          object$preliminary_broad_cell_type
        )
      ),
    ambiguous_cell_fraction =
      ambiguous_fraction,
    max_cluster_donor_fraction =
      max(
        cluster_totals$max_donor_fraction
      ),
    pseudobulk_group_n =
      ncol(
        pseudobulk$counts
      ),
    object_path =
      normalize_slashes(
        object_path
      ),
    pseudobulk_matrix_path =
      normalize_slashes(
        pseudobulk_matrix_path
      ),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    dataset_qc,
    file.path(
      dataset_table_root,
      paste0(
        dataset_id,
        "_integration_QC_summary.csv"
      )
    )
  )

  rm(
    object,
    data_matrix,
    counts_matrix
  )

  gc()

  list(
    dataset_qc =
      dataset_qc,
    sample_doublet_qc =
      sample_doublet_qc,
    cluster_annotation =
      cluster_annotation,
    donor_program_summary =
      donor_program_summary
  )
}

# ----------------------------------------------------------------------
# Run both datasets separately
# ----------------------------------------------------------------------

dataset_results <- list()

for (dataset_id in c(
  "GSE250130",
  "GSE183837"
)) {
  dataset_results[[
    dataset_id
  ]] <- process_dataset(
    dataset_id =
      dataset_id,
    sample_design =
      designs[[
        dataset_id
      ]]
  )
}

combined_dataset_qc <- do.call(
  rbind,
  lapply(
    dataset_results,
    function(x) {
      x$dataset_qc
    }
  )
)

combined_doublet_qc <- do.call(
  rbind,
  lapply(
    dataset_results,
    function(x) {
      x$sample_doublet_qc
    }
  )
)

write_csv_utf8(
  combined_dataset_qc,
  file.path(
    table_root,
    "STEP04E1_DATASET_QC_SUMMARY.csv"
  )
)

write_csv_utf8(
  combined_doublet_qc,
  file.path(
    table_root,
    "STEP04E1_SAMPLE_DOUBLET_QC.csv"
  )
)

# ----------------------------------------------------------------------
# Gate
# ----------------------------------------------------------------------

all_samples_retained <- all(
  combined_dataset_qc$locked_sample_n ==
    combined_dataset_qc$retained_sample_n
)

ambiguous_pass <- all(
  combined_dataset_qc$ambiguous_cell_fraction <
    0.20
)

donor_mixing_pass <- all(
  combined_dataset_qc$max_cluster_donor_fraction <
    0.80
)

doublet_pass <-
  doublet_mode ==
    "SCDBLFINDER"

gate_status <- if (
  all_samples_retained &&
    doublet_pass &&
    ambiguous_pass &&
    donor_mixing_pass
) {
  "READY_FOR_MANUAL_CELLTYPE_QC_AND_DONOR_AWARE_PSEUDOBULK_REVIEW"
} else if (
  !doublet_pass
) {
  "HOLD_SCDBLFINDER_PACKAGES_MISSING"
} else {
  "HOLD_STEP04E1_BIOLOGICAL_QC_REVIEW_REQUIRED"
}

gate_table <- data.frame(
  gate_component = c(
    "STEP04D-QC1 authorization",
    "Dataset integration strategy",
    "GSE250130 sample retention",
    "GSE183837 sample retention",
    "Doublet assessment",
    "Preliminary broad-cell annotation",
    "Cluster donor-mixing audit",
    "Donor-aware VDR-program localization",
    "Preliminary broad-cell pseudobulk",
    "GSE183837 watchlist",
    "Automatic donor exclusion",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP04E1 gate"
  ),
  status = c(
    "PASS",
    "DATASETS_PROCESSED_SEPARATELY",
    paste0(
      combined_dataset_qc$retained_sample_n[
        combined_dataset_qc$dataset_id ==
          "GSE250130"
      ],
      "_OF_28"
    ),
    paste0(
      combined_dataset_qc$retained_sample_n[
        combined_dataset_qc$dataset_id ==
          "GSE183837"
      ],
      "_OF_9"
    ),
    doublet_mode,
    ifelse(
      ambiguous_pass,
      "PASS_AMBIGUOUS_CELL_FRACTION_LT_0.20",
      "REVIEW_AMBIGUOUS_CELL_FRACTION"
    ),
    ifelse(
      donor_mixing_pass,
      "PASS_MAX_CLUSTER_DONOR_FRACTION_LT_0.80",
      "REVIEW_DONOR_DOMINATED_CLUSTER"
    ),
    "COMPLETED_SAMPLE_BY_BROAD_CELL_TYPE",
    "BUILT_MINIMUM_50_CELLS_PER_GROUP",
    "ACTIVE_2_REVIEW_REQUIRED_2_CAUTION",
    "NONE",
    "NOT_PERFORMED",
    "NONE_CREATED",
    gate_status
  ),
  required_action = c(
    "Retain corrected v1.3 bulk handoff.",
    "Do not pool the two studies into one integrated embedding.",
    "Review temporal-stage representation by donor and cell type.",
    "Review control/RIF representation and watchlist donors.",
    "Use scDblFinder singlets for downstream work.",
    "Inspect cluster markers and revise preliminary labels before locking.",
    "Inspect any donor-dominated cluster before biological interpretation.",
    "Use donors, not cells, as biological replicates.",
    "Do not run final pseudobulk DE until broad-cell labels are locked.",
    "Inspect GSM5572244 and GSM5572245 most closely; retain GSM5572238 and GSM5572240 as caution donors.",
    "Do not exclude a donor solely from the pre-QC watchlist.",
    "Use occupancy-supported candidate terminology.",
    "Render scientific figures only after cell-type QC and pseudobulk review.",
    "Proceed only according to the Gate."
  ),
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  output_root,
  "STEP04E1_GATE.csv"
)

write_csv_utf8(
  gate_table,
  gate_path
)

result_card <- c(
  "# STEP04E1 single-cell integration handoff card",
  "",
  paste0(
    "- Run ID: ",
    run_id
  ),
  paste0(
    "- Integration version: ",
    integration_version,
    " CANDIDATE"
  ),
  paste0(
    "- Annotation version: ",
    annotation_version,
    " PRELIMINARY"
  ),
  paste0(
    "- Localization version: ",
    localization_version,
    " CANDIDATE"
  ),
  "- Integration strategy: GSE250130 and GSE183837 processed separately",
  paste0(
    "- Doublet mode: ",
    doublet_mode
  ),
  paste0(
    "- GSE250130 retained samples/cells: ",
    combined_dataset_qc$retained_sample_n[
      combined_dataset_qc$dataset_id ==
        "GSE250130"
    ],
    "/",
    combined_dataset_qc$retained_cell_n[
      combined_dataset_qc$dataset_id ==
        "GSE250130"
    ]
  ),
  paste0(
    "- GSE183837 retained samples/cells: ",
    combined_dataset_qc$retained_sample_n[
      combined_dataset_qc$dataset_id ==
        "GSE183837"
    ],
    "/",
    combined_dataset_qc$retained_cell_n[
      combined_dataset_qc$dataset_id ==
        "GSE183837"
    ]
  ),
  paste0(
    "- GSE250130 clusters/ambiguous fraction: ",
    combined_dataset_qc$cluster_n[
      combined_dataset_qc$dataset_id ==
        "GSE250130"
    ],
    "/",
    format(
      combined_dataset_qc$ambiguous_cell_fraction[
        combined_dataset_qc$dataset_id ==
          "GSE250130"
      ],
      digits = 5
    )
  ),
  paste0(
    "- GSE183837 clusters/ambiguous fraction: ",
    combined_dataset_qc$cluster_n[
      combined_dataset_qc$dataset_id ==
        "GSE183837"
    ],
    "/",
    format(
      combined_dataset_qc$ambiguous_cell_fraction[
        combined_dataset_qc$dataset_id ==
          "GSE183837"
      ],
      digits = 5
    )
  ),
  "- Automatic donor exclusions: 0",
  "- Definitive direct-target assignments: 0",
  "- Scientific figures: NONE",
  paste0(
    "- Gate: ",
    gate_status
  )
)

write_text(
  result_card,
  file.path(
    output_root,
    "STEP04E1_CARD.md"
  )
)

write_text(
  c(
    "# STEP04E1 Methods draft",
    "",
    "The two single-cell studies were processed separately because they serve different biological roles: GSE250130 represents a time-series luteal-phase atlas with an LH+7 RIF comparison, whereas GSE183837 provides an independent control-RIF localization cohort. Donor/sample identity was retained throughout.",
    "",
    "Within each sample, duplicated gene symbols were collapsed by summing sparse counts. When available, scDblFinder was applied independently to each sample and predicted doublets were removed before merging. Each dataset was log normalized, highly variable genes were selected, and PCA was calculated. Harmony was applied to the PCA representation using sample accession as the batch variable. Shared-nearest-neighbor clustering and UMAP were calculated from the Harmony embedding.",
    "",
    "Broad cell compartments were assigned preliminarily from cluster-level canonical marker scores for epithelial, stromal fibroblast, endothelial, immune, pericyte/smooth-muscle and erythroid lineages. Cluster markers, donor mixing and ambiguous-label fractions were retained for manual annotation QC.",
    "",
    "Four locked VDR-related transcriptional programs were mapped to single cells using direction-aware signed normalized-expression scores. Scores were summarized at donor-by-broad-cell-type level. Raw counts were aggregated into preliminary donor-by-broad-cell-type pseudobulk matrices for groups containing at least 50 cells. No cell-level significance test was used as evidence of biological replication."
  ),
  file.path(
    dirs$text,
    "STEP04E1_Methods_draft.md"
  )
)

write_text(
  c(
    "# STEP04E1 Results shell",
    "",
    "This Results file remains a shell until preliminary broad-cell labels, cluster markers, donor mixing and watchlist-donor behavior have passed manual QC.",
    "",
    "No scientific figure, causal claim, diagnostic claim or definitive direct-target statement is authorized at this stage."
  ),
  file.path(
    dirs$text,
    "STEP04E1_Results_shell.md"
  )
)

parameters <- list(
  schema_version = "1.0",
  run_id =
    run_id,
  integration_version =
    integration_version,
  annotation_version =
    annotation_version,
  localization_version =
    localization_version,
  doublet_mode =
    doublet_mode,
  variable_feature_n =
    3000,
  PCA_n =
    50,
  harmony_dimension_n =
    40,
  graph_dimension_n =
    30,
  main_cluster_resolution =
    0.4,
  pseudobulk_minimum_cell_n =
    50,
  automatic_donor_exclusion_n =
    0,
  definitive_direct_target_assignment =
    FALSE,
  scientific_figures_created =
    FALSE,
  gate =
    gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP04E1_parameters.json"
  ),
  pretty =
    TRUE,
  auto_unbox =
    TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04E1.txt"
  )
)

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) {
    NULL
  }
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  copy_if_different(
    source_file,
    file.path(
      dirs$code,
      "53_step04E1_single_cell_integration_and_localization.R"
    )
  )
}

# ----------------------------------------------------------------------
# Manifest and short archive
# ----------------------------------------------------------------------

small_artifacts <- c(
  gate_path,
  file.path(
    output_root,
    "STEP04E1_CARD.md"
  ),
  file.path(
    table_root,
    "STEP04E1_DATASET_QC_SUMMARY.csv"
  ),
  file.path(
    table_root,
    "STEP04E1_SAMPLE_DOUBLET_QC.csv"
  ),
  file.path(
    dirs$text,
    "STEP04E1_Methods_draft.md"
  ),
  file.path(
    dirs$text,
    "STEP04E1_Results_shell.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04E1_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04E1.txt"
  ),
  file.path(
    dirs$code,
    "53_step04E1_single_cell_integration_and_localization.R"
  ),
  log_path
)

dataset_table_artifacts <- list.files(
  table_root,
  pattern =
    "\\.(csv|md)$",
  recursive =
    TRUE,
  full.names =
    TRUE
)

pseudobulk_metadata_artifacts <- list.files(
  pseudobulk_root,
  pattern =
    "\\.csv$",
  recursive =
    TRUE,
  full.names =
    TRUE
)

object_artifacts <- list.files(
  object_root,
  pattern =
    "\\.rds$",
  recursive =
    TRUE,
  full.names =
    TRUE
)

pseudobulk_matrix_artifacts <- list.files(
  pseudobulk_root,
  pattern =
    "\\.rds$",
  recursive =
    TRUE,
  full.names =
    TRUE
)

artifact_paths <- unique(
  c(
    small_artifacts,
    dataset_table_artifacts,
    pseudobulk_metadata_artifacts,
    object_artifacts,
    pseudobulk_matrix_artifacts
  )
)

artifact_paths <- artifact_paths[
  file.exists(
    artifact_paths
  )
]

artifact_manifest <- data.frame(
  path =
    normalize_slashes(
      artifact_paths
    ),
  artifact_class =
    ifelse(
      artifact_paths %in%
        object_artifacts,
      "LARGE_INTEGRATED_OBJECT",
      ifelse(
        artifact_paths %in%
          pseudobulk_matrix_artifacts,
        "PSEUDOBULK_MATRIX",
        "SMALL_QC_OR_TEXT_ARTIFACT"
      )
    ),
  size_bytes =
    as.numeric(
      file.info(
        artifact_paths
      )$size
    ),
  md5 =
    vapply(
      artifact_paths,
      md5_file,
      character(1)
    ),
  sha256 =
    vapply(
      artifact_paths,
      sha256_file,
      character(1)
    ),
  artifact_status =
    ifelse(
      gate_status ==
        "READY_FOR_MANUAL_CELLTYPE_QC_AND_DONOR_AWARE_PSEUDOBULK_REVIEW",
      "CANDIDATE_PENDING_MANUAL_CELLTYPE_QC",
      "HOLD_PENDING_QC_OR_PACKAGE_COMPLETION"
    ),
  stringsAsFactors = FALSE
)

manifest_path <- file.path(
  output_root,
  "STEP04E1_ARTIFACT_MANIFEST.csv"
)

write_csv_utf8(
  artifact_manifest,
  manifest_path
)

archive_artifacts <- unique(
  c(
    small_artifacts,
    dataset_table_artifacts,
    pseudobulk_metadata_artifacts,
    manifest_path
  )
)

archive_artifacts <- archive_artifacts[
  file.exists(
    archive_artifacts
  )
]

for (path_value in archive_artifacts) {
  section <- if (
    grepl(
      "\\.R$",
      basename(
        path_value
      )
    )
  ) {
    "code"
  } else if (
    grepl(
      "sessionInfo|parameters",
      basename(
        path_value
      )
    )
  ) {
    "env"
  } else if (
    grepl(
      "\\.log$",
      basename(
        path_value
      )
    )
  ) {
    "log"
  } else if (
    grepl(
      "\\.md$",
      basename(
        path_value
      )
    )
  ) {
    "txt"
  } else if (
    grepl(
      "GATE|CARD|MANIFEST|QC|MARKER|ANNOTATION|MIXING|COMPOSITION",
      basename(
        path_value
      ),
      ignore.case =
        TRUE
    )
  ) {
    "qc"
  } else {
    "res"
  }

  destination_path <- file.path(
    archive_root,
    section,
    basename(
      path_value
    )
  )

  if (
    file.exists(
      destination_path
    )
  ) {
    destination_path <- file.path(
      archive_root,
      section,
      paste0(
        basename(
          dirname(
            path_value
          )
        ),
        "__",
        basename(
          path_value
        )
      )
    )
  }

  copy_validate(
    path_value,
    destination_path
  )
}

archive_files <- list.files(
  archive_root,
  recursive =
    TRUE,
  full.names =
    TRUE,
  include.dirs =
    FALSE
)

archive_manifest <- data.frame(
  relative_path =
    substring(
      normalize_slashes(
        archive_files
      ),
      nchar(
        normalize_slashes(
          archive_root
        )
      ) + 2L
    ),
  size_bytes =
    as.numeric(
      file.info(
        archive_files
      )$size
    ),
  md5 =
    vapply(
      archive_files,
      md5_file,
      character(1)
    ),
  sha256 =
    vapply(
      archive_files,
      sha256_file,
      character(1)
    ),
  stringsAsFactors =
    FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP04E1_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04E1 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Integration version: ",
      integration_version,
      " CANDIDATE"
    ),
    paste0(
      "- Doublet mode: ",
      doublet_mode
    ),
    "- Datasets integrated separately: YES",
    "- Large Seurat objects copied into short archive: NO",
    "- Automatic donor exclusions: 0",
    "- Definitive direct-target assignments: 0",
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
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
    "# STEP04 combined package status",
    "",
    paste0(
      "- Current run: ",
      run_id
    ),
    "- Bulk projection: BULK_VDR_PROGRAM_PROJECTION_v1.3 LOCKED",
    paste0(
      "- Single-cell integration: ",
      integration_version,
      " CANDIDATE"
    ),
    paste0(
      "- Broad-cell annotation: ",
      annotation_version,
      " PRELIMINARY"
    ),
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
    )
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
cat("STEP04E1 single-cell integration completed.\n")
cat(
  "Integration version: ",
  integration_version,
  " CANDIDATE\n"
)
cat(
  "Doublet mode: ",
  doublet_mode,
  "\n"
)
for (dataset_id in combined_dataset_qc$dataset_id) {
  dataset_row <- combined_dataset_qc[
    combined_dataset_qc$dataset_id ==
      dataset_id,
    ,
    drop = FALSE
  ]

  cat(
    dataset_id,
    " samples/cells/clusters: ",
    dataset_row$retained_sample_n,
    " / ",
    dataset_row$retained_cell_n,
    " / ",
    dataset_row$cluster_n,
    "\n",
    sep = ""
  )

  cat(
    dataset_id,
    " ambiguous fraction: ",
    format(
      dataset_row$ambiguous_cell_fraction,
      digits = 6
    ),
    "\n",
    sep = ""
  )
}
cat("Automatic donor exclusions: 0\n")
cat("Definitive direct-target assignments: 0\n")
cat("Scientific figures: NONE\n")
cat(
  "Gate: ",
  gate_status,
  "\n"
)
cat(
  "Authoritative archive:\n  ",
  archive_root,
  "\n\n"
)
cat("Upload:\n")
cat(
  "  06_results/S4E1_sc/SCI_v1_0/STEP04E1_GATE.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/STEP04E1_CARD.md\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/STEP04E1_DATASET_QC_SUMMARY.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/STEP04E1_SAMPLE_DOUBLET_QC.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_preliminary_cluster_annotation.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_preliminary_cluster_annotation.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_top25_cluster_markers.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_top25_cluster_markers.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_donor_celltype_VDR_program_scores.csv\n"
)
cat(
  "  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_donor_celltype_VDR_program_scores.csv\n"
)
cat("============================================================\n")

# ==============================================================================
# Analysis stage 2: 56_step04E1B_GSE250130_integration_resume_v1
# ==============================================================================

# ======================================================================
# 56_step04E1B_GSE250130_integration_resume_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Resume GSE250130 from the locked post-doublet merged-count checkpoint.
# - Do NOT repeat sample-level scDblFinder/doublet processing.
# - Re-run normalization, HVG, scaling, PCA, Harmony, graph, clusters, UMAP.
# - Save a PRE-JOIN checkpoint before any layer joining.
# - Join ONLY normalized data layers; keep counts layers separate.
# - Save a data-joined integrated checkpoint for the next annotation step.
#
# WHY THIS SCRIPT EXISTS
# The prior run reached:
#   "Joining RNA layers after PCA/Harmony/UMAP for GSE250130"
# but no integrated checkpoint was saved. The likely failure point was the
# full JoinLayers operation, which attempted to join both counts and data.
#
# GOVERNANCE
# - New immutable output root: SCI_v1_1_memory_safe
# - Existing SCI_v1_0 is read-only and is never overwritten.
# - This script processes GSE250130 only.
# - No marker testing, annotation, VDR localization, or pseudobulk here.
# - Results remain CANDIDATE until STEP04E1 biological QC is completed.
#
# RUN FROM RSTUDIO CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(7200L, getOption("timeout")))
options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260804)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <- "56_step04E1B_GSE250130_integration_resume_v1.R"
script_version <- "STEP04E1B_GSE250130_INTEGRATION_RESUME_v1.0"
integration_version <- "SINGLE_CELL_INTEGRATION_v1.1_MEMORY_SAFE"
output_version <- "SCI_v1_1_memory_safe"

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1B_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e1) {
      read.csv(
        path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )
}

format_bytes <- function(bytes) {
  if (is.na(bytes)) return("")
  if (bytes < 1024) return(sprintf("%.0f B", bytes))
  if (bytes < 1024^2) return(sprintf("%.1f KB", bytes / 1024))
  if (bytes < 1024^3) return(sprintf("%.1f MB", bytes / 1024^2))
  sprintf("%.2f GB", bytes / 1024^3)
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") return(NA_real_)
  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"$o=Get-CimInstance Win32_OperatingSystem; ',
        '[math]::Round($o.FreePhysicalMemory/1MB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )
  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

get_free_disk_gb <- function(drive_letter = "C") {
  if (.Platform$OS.type != "windows") return(NA_real_)
  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"[math]::Round((Get-PSDrive -Name ',
        drive_letter,
        ').Free/1GB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )
  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

object_size_gb <- function(x) {
  as.numeric(utils::object.size(x)) / 1024^3
}

md5_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

log_event <- NULL
run_status <- "STARTED"
failure_stage <- ""
failure_message <- ""
prejoin_saved <- FALSE
data_join_saved <- FALSE
zip_ok <- FALSE
zip_error <- ""


# ----------------------------------------------------------------------
# 2. Project directories
# ----------------------------------------------------------------------

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

dirs <- list(
  parameters = file.path(step04_root, "03_parameters"),
  environment = file.path(step04_root, "04_environment"),
  execution = file.path(step04_root, "05_execution"),
  results = file.path(step04_root, "06_results"),
  qc = file.path(step04_root, "08_qc"),
  archive = file.path(step04_root, "10_archive_work")
)

for (p in dirs) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(p)) stop("Cannot create required directory: ", p)
}

old_output_root <- file.path(
  dirs$results,
  "S4E1_sc",
  "SCI_v1_0"
)

new_output_root <- file.path(
  dirs$results,
  "S4E1_sc",
  output_version
)

dataset_checkpoint_root <- file.path(
  new_output_root,
  "checkpoints",
  dataset_id
)

dataset_table_root <- file.path(
  new_output_root,
  "tables",
  dataset_id
)

qc_run_root <- file.path(
  dirs$qc,
  paste0(run_id, "_integration_resume")
)

for (p in c(
  new_output_root,
  dataset_checkpoint_root,
  dataset_table_root,
  qc_run_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  dirs$execution,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  dirs$execution,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  dirs$parameters,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  dirs$environment,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)
desktop_candidates <- desktop_candidates[dir.exists(desktop_candidates)]
zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )
  cat(line, "\n")
  cat(line, "\n", file = log_path, append = TRUE)
  flush.console()
}

run_step <- function(label, fun) {
  log_event("START", label)
  t0 <- Sys.time()

  result <- tryCatch(
    fun(),
    error = function(e) {
      failure_stage <<- label
      failure_message <<- conditionMessage(e)
      log_event("ERROR", label, failure_message)
      stop(e)
    }
  )

  elapsed <- round(
    as.numeric(difftime(Sys.time(), t0, units = "mins")),
    2
  )

  log_event("DONE", label, "elapsed_min=", elapsed)
  result
}


# ----------------------------------------------------------------------
# 3. Validate packages, resources, and source checkpoint
# ----------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "harmony",
  "data.table",
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

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

source_checkpoint <- file.path(
  old_output_root,
  "checkpoints",
  dataset_id,
  paste0(
    dataset_id,
    "_post_doublet_merged_counts_checkpoint.rds"
  )
)

source_status_path <- file.path(
  old_output_root,
  "checkpoints",
  dataset_id,
  paste0(
    dataset_id,
    "_post_doublet_merged_counts_status.csv"
  )
)

sample_doublet_qc_path <- file.path(
  old_output_root,
  "tables",
  dataset_id,
  paste0(dataset_id, "_sample_doublet_qc.csv")
)

locked_design_path <- file.path(
  dirs$results,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0",
  "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"
)

required_inputs <- c(
  source_checkpoint,
  source_status_path,
  sample_doublet_qc_path,
  locked_design_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Required input(s) missing: ",
    paste(missing_inputs, collapse = "; ")
  )
}

source_status <- safe_read_csv(source_status_path)
sample_design <- safe_read_csv(locked_design_path)
sample_doublet_qc <- safe_read_csv(sample_doublet_qc_path)

if (is.null(source_status) || nrow(source_status) != 1L) {
  stop("Invalid source checkpoint status table.")
}

if (
  source_status$dataset_id[[1]] != dataset_id ||
    source_status$checkpoint_status[[1]] != "COMPLETE" ||
    source_status$sample_n[[1]] != 28L ||
    source_status$cell_n[[1]] != 236636L
) {
  stop(
    "Source checkpoint status does not match the locked diagnostic facts."
  )
}

if (
  is.null(sample_design) ||
    nrow(sample_design) != 28L ||
    !"sample_accession" %in% names(sample_design)
) {
  stop("Locked GSE250130 sample design is invalid.")
}

if (
  is.null(sample_doublet_qc) ||
    nrow(sample_doublet_qc) != 28L
) {
  stop("GSE250130 sample doublet QC table is invalid.")
}

free_memory_gb <- get_free_memory_gb()
free_disk_gb <- get_free_disk_gb("C")

resource_snapshot <- data.frame(
  metric = c(
    "free_physical_memory_GB_before_load",
    "free_C_drive_GB_before_run",
    "source_checkpoint_size_GB",
    "source_checkpoint_md5"
  ),
  value = c(
    free_memory_gb,
    free_disk_gb,
    round(file.info(source_checkpoint)$size / 1024^3, 3),
    md5_file(source_checkpoint)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resource_snapshot,
  file.path(qc_run_root, "01_RESOURCE_AND_INPUT_SNAPSHOT.csv")
)

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
  stop(
    "Free physical memory is below 18 GB. Restart RStudio and rerun in a clean session."
  )
}

if (is.finite(free_disk_gb) && free_disk_gb < 25) {
  stop(
    "Free C: drive space is below 25 GB. Free disk space before running."
  )
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  integration_version = integration_version,
  output_version = output_version,
  dataset_id = dataset_id,
  source_checkpoint = normalize_slashes(source_checkpoint),
  source_checkpoint_status = "COMPLETE",
  expected_sample_n = 28L,
  expected_cell_n = 236636L,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  normalization_block_size = 2000L,
  variable_feature_method = "vst",
  variable_feature_n = 3000L,
  scaling_block_size = 1000L,
  pca_npcs = 50L,
  harmony_dims = "1:40",
  neighbor_dims = "1:30",
  neighbor_k = 30L,
  cluster_resolutions = c(0.2, 0.6, 0.4),
  umap_dims = "1:30",
  random_seed = 20260804L,
  layer_join_policy = "JOIN_DATA_ONLY_KEEP_COUNTS_SPLIT",
  candidate_status = "NOT_LOCKED"
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

log_event(
  "Starting",
  run_id,
  "source_checkpoint=",
  source_checkpoint,
  "free_memory_GB=",
  free_memory_gb,
  "free_disk_GB=",
  free_disk_gb
)


# ----------------------------------------------------------------------
# 4. Output checkpoint paths
# ----------------------------------------------------------------------

prejoin_checkpoint <- file.path(
  dataset_checkpoint_root,
  paste0(
    dataset_id,
    "_integrated_prejoin_layers_checkpoint.rds"
  )
)

prejoin_status_path <- file.path(
  dataset_checkpoint_root,
  paste0(
    dataset_id,
    "_integrated_prejoin_layers_status.csv"
  )
)

data_join_checkpoint <- file.path(
  dataset_checkpoint_root,
  paste0(
    dataset_id,
    "_integrated_checkpoint.rds"
  )
)

data_join_status_path <- file.path(
  dataset_checkpoint_root,
  paste0(
    dataset_id,
    "_checkpoint_status.csv"
  )
)

# Never silently overwrite accepted/candidate checkpoints from a prior run.
for (p in c(
  prejoin_checkpoint,
  prejoin_status_path,
  data_join_checkpoint,
  data_join_status_path
)) {
  if (file.exists(p)) {
    stop(
      "Target already exists. Refusing silent overwrite: ",
      p
    )
  }
}


# ----------------------------------------------------------------------
# 5. Main pipeline with guaranteed final reporting
# ----------------------------------------------------------------------

pipeline_error <- NULL
object <- NULL

tryCatch(
  {
    run_status <- "RUNNING"

    object <- run_step(
      "LOAD_POST_DOUBLET_CHECKPOINT",
      function() readRDS(source_checkpoint)
    )

    if (!inherits(object, "Seurat")) {
      stop("Source checkpoint is not a Seurat object.")
    }

    if (ncol(object) != 236636L) {
      stop(
        "Unexpected cell count after checkpoint load: ",
        ncol(object)
      )
    }

    if (
      !"sample_accession" %in% colnames(object@meta.data) ||
        data.table::uniqueN(object$sample_accession) != 28L
    ) {
      stop("Checkpoint does not retain all 28 locked samples.")
    }

    SeuratObject::DefaultAssay(object) <- "RNA"

    layer_before <- SeuratObject::Layers(object[["RNA"]])

    write_csv_utf8(
      data.frame(
        phase = "SOURCE_POST_DOUBLET",
        layer = layer_before,
        stringsAsFactors = FALSE
      ),
      file.path(qc_run_root, "02_LAYER_MANIFEST_SOURCE.csv")
    )

    log_event(
      "Loaded source object",
      "cells=", ncol(object),
      "features=", nrow(object),
      "samples=", data.table::uniqueN(object$sample_accession),
      "object_size_GB=", round(object_size_gb(object), 3),
      "RNA_layer_n=", length(layer_before)
    )

    object <- run_step(
      "NORMALIZE_DATA_LAYERED",
      function() {
        Seurat::NormalizeData(
          object,
          assay = "RNA",
          normalization.method = "LogNormalize",
          scale.factor = 10000,
          block.size = 2000,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "FIND_VARIABLE_FEATURES",
      function() {
        Seurat::FindVariableFeatures(
          object,
          assay = "RNA",
          selection.method = "vst",
          nfeatures = 3000,
          verbose = TRUE
        )
      }
    )

    if (length(SeuratObject::VariableFeatures(object)) != 3000L) {
      log_event(
        "WARNING variable feature count=",
        length(SeuratObject::VariableFeatures(object))
      )
    }

    object <- run_step(
      "SCALE_VARIABLE_FEATURES",
      function() {
        Seurat::ScaleData(
          object,
          assay = "RNA",
          features = SeuratObject::VariableFeatures(object),
          block.size = 1000,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "RUN_PCA_50",
      function() {
        Seurat::RunPCA(
          object,
          assay = "RNA",
          features = SeuratObject::VariableFeatures(object),
          npcs = 50,
          seed.use = 20260804,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "RUN_HARMONY_SAMPLE_ACCESSION",
      function() {
        harmony::RunHarmony(
          object = object,
          group.by.vars = "sample_accession",
          reduction = "pca",
          dims.use = 1:40,
          reduction.save = "harmony",
          plot_convergence = FALSE,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "FIND_NEIGHBORS_HARMONY",
      function() {
        Seurat::FindNeighbors(
          object,
          reduction = "harmony",
          dims = 1:30,
          k.param = 30,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "FIND_CLUSTERS_RES_0_2",
      function() {
        Seurat::FindClusters(
          object,
          resolution = 0.2,
          algorithm = 1,
          random.seed = 20260804,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "FIND_CLUSTERS_RES_0_6",
      function() {
        Seurat::FindClusters(
          object,
          resolution = 0.6,
          algorithm = 1,
          random.seed = 20260804,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "FIND_CLUSTERS_RES_0_4",
      function() {
        Seurat::FindClusters(
          object,
          resolution = 0.4,
          algorithm = 1,
          random.seed = 20260804,
          verbose = TRUE
        )
      }
    )

    object <- run_step(
      "RUN_UMAP_HARMONY",
      function() {
        Seurat::RunUMAP(
          object,
          reduction = "harmony",
          dims = 1:30,
          reduction.name = "umap_harmony",
          reduction.key = "UMAPH_",
          seed.use = 20260804,
          verbose = TRUE
        )
      }
    )

    # Remove scale.data before checkpointing/joining.
    reductions_to_keep <- intersect(
      c("pca", "harmony", "umap_harmony"),
      SeuratObject::Reductions(object)
    )

    graphs_to_keep <- SeuratObject::Graphs(object)

    layers_before_diet <- SeuratObject::Layers(object[["RNA"]])

    RNA_layers_to_keep <- grep(
      "^(counts|data)(\\.|$)",
      layers_before_diet,
      value = TRUE
    )

    if (!length(RNA_layers_to_keep)) {
      stop("No counts/data RNA layers found before DietSeurat.")
    }

    object <- run_step(
      "DIET_OBJECT_REMOVE_SCALE_DATA",
      function() {
        Seurat::DietSeurat(
          object,
          assays = "RNA",
          layers = list(RNA = RNA_layers_to_keep),
          dimreducs = reductions_to_keep,
          graphs = graphs_to_keep,
          misc = TRUE
        )
      }
    )

    gc()

    prejoin_layers <- SeuratObject::Layers(object[["RNA"]])

    cluster_columns <- grep(
      "_snn_res\\.",
      colnames(object@meta.data),
      value = TRUE
    )

    umap_embeddings <- Seurat::Embeddings(
      object,
      reduction = "umap_harmony"
    )

    embedding_table <- data.frame(
      cell_barcode = rownames(umap_embeddings),
      UMAPH_1 = umap_embeddings[, 1],
      UMAPH_2 = umap_embeddings[, 2],
      sample_accession = object$sample_accession,
      seurat_clusters = as.character(object$seurat_clusters),
      stringsAsFactors = FALSE
    )

    for (cc in cluster_columns) {
      embedding_table[[cc]] <- object@meta.data[[cc]]
    }

    write_csv_utf8(
      embedding_table,
      file.path(
        dataset_table_root,
        paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv")
      )
    )

    cluster_sample_counts <- data.table::as.data.table(
      object@meta.data
    )[
      ,
      .(cell_n = .N),
      by = .(
        sample_accession,
        seurat_clusters
      )
    ]

    write_csv_utf8(
      as.data.frame(cluster_sample_counts),
      file.path(
        dataset_table_root,
        paste0(dataset_id, "_cluster_by_sample_prejoin.csv")
      )
    )

    write_csv_utf8(
      data.frame(
        phase = "PREJOIN_AFTER_UMAP",
        layer = prejoin_layers,
        stringsAsFactors = FALSE
      ),
      file.path(qc_run_root, "03_LAYER_MANIFEST_PREJOIN.csv")
    )

    log_event(
      "Saving PRE-JOIN checkpoint",
      "object_size_GB=", round(object_size_gb(object), 3),
      "layers=", length(prejoin_layers),
      "reductions=", paste(reductions_to_keep, collapse = ";"),
      "graphs=", paste(graphs_to_keep, collapse = ";")
    )

    run_step(
      "SAVE_PREJOIN_CHECKPOINT",
      function() {
        saveRDS(
          object,
          prejoin_checkpoint,
          compress = FALSE
        )
        TRUE
      }
    )

    write_csv_utf8(
      data.frame(
        dataset_id = dataset_id,
        integration_version = integration_version,
        output_version = output_version,
        checkpoint_stage = "PREJOIN_AFTER_PCA_HARMONY_CLUSTER_UMAP",
        checkpoint_status = "COMPLETE",
        sample_n = data.table::uniqueN(object$sample_accession),
        cell_n = ncol(object),
        feature_n = nrow(object),
        cluster_n = length(unique(object$seurat_clusters)),
        layer_n = length(prejoin_layers),
        checkpoint_size_bytes = file.info(prejoin_checkpoint)$size,
        checkpoint_md5 = md5_file(prejoin_checkpoint),
        checkpoint_path = normalize_slashes(prejoin_checkpoint),
        stringsAsFactors = FALSE
      ),
      prejoin_status_path
    )

    prejoin_saved <- TRUE
    log_event(
      "PRE-JOIN checkpoint safely saved",
      prejoin_checkpoint,
      format_bytes(file.info(prejoin_checkpoint)$size)
    )

    # Join normalized DATA only. Counts remain sample-specific layers.
    data_layers_before <- SeuratObject::Layers(
      object[["RNA"]],
      search = "^data"
    )

    counts_layers_before <- SeuratObject::Layers(
      object[["RNA"]],
      search = "^counts"
    )

    if (!length(data_layers_before)) {
      stop("No normalized data layers found for data-only JoinLayers.")
    }

    log_event(
      "Joining DATA layers only",
      "data_layer_n=", length(data_layers_before),
      "counts_layer_n_kept_split=", length(counts_layers_before)
    )

    object <- run_step(
      "JOIN_NORMALIZED_DATA_LAYERS_ONLY",
      function() {
        SeuratObject::JoinLayers(
          object,
          assay = "RNA",
          layers = "data",
          new = "data"
        )
      }
    )

    gc()

    data_layers_after <- SeuratObject::Layers(
      object[["RNA"]],
      search = "^data"
    )

    counts_layers_after <- SeuratObject::Layers(
      object[["RNA"]],
      search = "^counts"
    )

    write_csv_utf8(
      rbind(
        data.frame(
          phase = "AFTER_DATA_ONLY_JOIN",
          layer_type = "data",
          layer = data_layers_after,
          stringsAsFactors = FALSE
        ),
        data.frame(
          phase = "AFTER_DATA_ONLY_JOIN",
          layer_type = "counts",
          layer = counts_layers_after,
          stringsAsFactors = FALSE
        )
      ),
      file.path(qc_run_root, "04_LAYER_MANIFEST_AFTER_DATA_JOIN.csv")
    )

    if (!identical(data_layers_after, "data")) {
      stop(
        "Data-only JoinLayers did not produce exactly one 'data' layer. Found: ",
        paste(data_layers_after, collapse = ";")
      )
    }

    if (!length(counts_layers_after)) {
      stop("Counts layers were unexpectedly lost after data-only join.")
    }

    log_event(
      "Saving DATA-JOINED integrated checkpoint",
      "object_size_GB=", round(object_size_gb(object), 3),
      "joined_data_layer=", paste(data_layers_after, collapse = ";"),
      "counts_layer_n=", length(counts_layers_after)
    )

    run_step(
      "SAVE_DATA_JOINED_INTEGRATED_CHECKPOINT",
      function() {
        saveRDS(
          object,
          data_join_checkpoint,
          compress = FALSE
        )
        TRUE
      }
    )

    write_csv_utf8(
      data.frame(
        dataset_id = dataset_id,
        integration_version = integration_version,
        output_version = output_version,
        doublet_mode = "SCDBLFINDER_REUSED_FROM_SCI_v1_0",
        sample_n = data.table::uniqueN(object$sample_accession),
        cell_n = ncol(object),
        feature_n = nrow(object),
        cluster_n = length(unique(object$seurat_clusters)),
        checkpoint_status = "COMPLETE",
        join_policy = "DATA_ONLY_COUNTS_REMAIN_SPLIT",
        joined_data_layer_n = length(data_layers_after),
        remaining_counts_layer_n = length(counts_layers_after),
        checkpoint_size_bytes = file.info(data_join_checkpoint)$size,
        checkpoint_md5 = md5_file(data_join_checkpoint),
        checkpoint_path = normalize_slashes(data_join_checkpoint),
        stringsAsFactors = FALSE
      ),
      data_join_status_path
    )

    data_join_saved <- TRUE
    run_status <- "PASS_GSE250130_INTEGRATION_AND_DATA_JOIN_COMPLETE"

    log_event(
      "PASS",
      run_status,
      data_join_checkpoint,
      format_bytes(file.info(data_join_checkpoint)$size)
    )
  },
  error = function(e) {
    pipeline_error <<- e

    if (prejoin_saved && !data_join_saved) {
      run_status <<- "HOLD_PREJOIN_SAVED_DATA_JOIN_FAILED"
    } else if (!prejoin_saved) {
      run_status <<- "FAIL_BEFORE_PREJOIN_CHECKPOINT"
    } else {
      run_status <<- "FAIL_AFTER_DATA_JOIN_CHECKPOINT"
    }

    log_event(
      "PIPELINE_FINAL_STATUS",
      run_status,
      "failure_stage=", failure_stage,
      "message=", failure_message
    )
  }
)


# ----------------------------------------------------------------------
# 6. Final QC, run manifest, sessionInfo, and small upload ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_table <- data.frame(
  package = c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "harmony",
    "data.table",
    "jsonlite",
    "future"
  ),
  version = vapply(
    c(
      "Seurat",
      "SeuratObject",
      "Matrix",
      "harmony",
      "data.table",
      "jsonlite",
      "future"
    ),
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(qc_run_root, "05_PACKAGE_VERSIONS.csv")
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

final_checkpoint_audit <- data.frame(
  artifact = c(
    "source_post_doublet_checkpoint",
    "prejoin_integrated_checkpoint",
    "data_joined_integrated_checkpoint"
  ),
  exists = c(
    file.exists(source_checkpoint),
    file.exists(prejoin_checkpoint),
    file.exists(data_join_checkpoint)
  ),
  size_bytes = c(
    if (file.exists(source_checkpoint)) file.info(source_checkpoint)$size else NA_real_,
    if (file.exists(prejoin_checkpoint)) file.info(prejoin_checkpoint)$size else NA_real_,
    if (file.exists(data_join_checkpoint)) file.info(data_join_checkpoint)$size else NA_real_
  ),
  size_human = c(
    if (file.exists(source_checkpoint)) format_bytes(file.info(source_checkpoint)$size) else "",
    if (file.exists(prejoin_checkpoint)) format_bytes(file.info(prejoin_checkpoint)$size) else "",
    if (file.exists(data_join_checkpoint)) format_bytes(file.info(data_join_checkpoint)$size) else ""
  ),
  md5 = c(
    md5_file(source_checkpoint),
    md5_file(prejoin_checkpoint),
    md5_file(data_join_checkpoint)
  ),
  path = normalize_slashes(c(
    source_checkpoint,
    prejoin_checkpoint,
    data_join_checkpoint
  )),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  final_checkpoint_audit,
  file.path(qc_run_root, "06_FINAL_CHECKPOINT_AUDIT.csv")
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  integration_version = integration_version,
  output_version = output_version,
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_minutes = round(
    as.numeric(difftime(run_end, run_start, units = "mins")),
    2
  ),
  run_status = run_status,
  failure_stage = failure_stage,
  failure_message = failure_message,
  source_checkpoint_reused = TRUE,
  sample_doublet_processing_repeated = FALSE,
  prejoin_checkpoint_saved = prejoin_saved,
  data_join_checkpoint_saved = data_join_saved,
  old_output_root_modified = FALSE,
  new_output_root = normalize_slashes(new_output_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(run_manifest, run_manifest_path)

summary_lines <- c(
  "# STEP04E1B GSE250130 integration resume",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Status: **", run_status, "**"),
  paste0("- Source checkpoint reused: `", normalize_slashes(source_checkpoint), "`"),
  "- Sample-level doublet detection was not repeated.",
  paste0("- Pre-join checkpoint saved: ", prejoin_saved),
  paste0("- Data-joined checkpoint saved: ", data_join_saved),
  paste0(
    "- Failure stage: ",
    ifelse(nzchar(failure_stage), failure_stage, "none")
  ),
  paste0(
    "- Failure message: ",
    ifelse(nzchar(failure_message), failure_message, "none")
  ),
  "",
  "## Interpretation",
  "",
  if (
    identical(run_status, "PASS_GSE250130_INTEGRATION_AND_DATA_JOIN_COMPLETE")
  ) {
    "- GSE250130 normalization, PCA, Harmony, graph, clustering, UMAP, and normalized-data layer joining completed."
  } else if (
    identical(run_status, "HOLD_PREJOIN_SAVED_DATA_JOIN_FAILED")
  ) {
    "- Integration completed and is safely checkpointed before JoinLayers; only data-layer joining remains unresolved."
  } else {
    "- Integration did not reach a safe pre-join checkpoint; inspect the execution log before rerunning."
  },
  "",
  "## Next step",
  "",
  "- Do not start GSE183837 yet.",
  "- Return the small ZIP and review the layer/checkpoint audit.",
  "- After PASS, run the separate annotation/VDR-localization/pseudobulk script.",
  "- Results remain CANDIDATE and are not yet eligible for locked results."
)

write_text(
  summary_lines,
  file.path(qc_run_root, "00_STEP04E1B_RUN_SUMMARY.md")
)

# Copy small evidence only; never include large checkpoints in the ZIP.
small_handoff_root <- file.path(
  qc_run_root,
  "small_handoff"
)
dir.create(small_handoff_root, recursive = TRUE, showWarnings = FALSE)

small_files <- unique(c(
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  list.files(
    qc_run_root,
    full.names = TRUE,
    recursive = FALSE
  ),
  list.files(
    dataset_table_root,
    pattern = "\\.csv$",
    full.names = TRUE
  ),
  prejoin_status_path,
  data_join_status_path,
  normalizePath(
    file.path(project_root, "04_code", "R", script_name),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    file.info(small_files)$size <= 100 * 1024^2
]

copy_rows <- list()

for (i in seq_along(small_files)) {
  src <- small_files[i]
  dst <- file.path(
    small_handoff_root,
    paste0(
      sprintf("%03d", i),
      "_",
      basename(src)
    )
  )
  ok <- file.copy(
    src,
    dst,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  copy_rows[[i]] <- data.frame(
    source = normalize_slashes(src),
    destination = normalize_slashes(dst),
    copied = isTRUE(ok),
    stringsAsFactors = FALSE
  )
}

if (length(copy_rows)) {
  write_csv_utf8(
    do.call(rbind, copy_rows),
    file.path(qc_run_root, "07_SMALL_HANDOFF_COPY_MANIFEST.csv")
  )
}

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = c(
          file.path(qc_run_root, "00_STEP04E1B_RUN_SUMMARY.md"),
          file.path(qc_run_root, "01_RESOURCE_AND_INPUT_SNAPSHOT.csv"),
          file.path(qc_run_root, "02_LAYER_MANIFEST_SOURCE.csv"),
          file.path(qc_run_root, "03_LAYER_MANIFEST_PREJOIN.csv"),
          file.path(qc_run_root, "04_LAYER_MANIFEST_AFTER_DATA_JOIN.csv"),
          file.path(qc_run_root, "05_PACKAGE_VERSIONS.csv"),
          file.path(qc_run_root, "06_FINAL_CHECKPOINT_AUDIT.csv"),
          file.path(qc_run_root, "07_SMALL_HANDOFF_COPY_MANIFEST.csv"),
          small_handoff_root
        ),
        recurse = TRUE,
        include_directories = TRUE
      )
      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1B GSE250130 resume finished\n")
cat("============================================================\n")
cat("Run status: ", run_status, "\n", sep = "")
cat("Execution log: ", log_path, "\n", sep = "")
cat("Run manifest: ", run_manifest_path, "\n", sep = "")
cat("QC folder: ", qc_run_root, "\n", sep = "")
cat("Pre-join checkpoint saved: ", prejoin_saved, "\n", sep = "")
cat("Data-joined checkpoint saved: ", data_join_saved, "\n", sep = "")

if (zip_ok) {
  cat("Upload this small ZIP:\n", handoff_zip, "\n")
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe QC folder is complete; compress it manually if needed.\n",
    sep = ""
  )
}

# Do not throw a final stop after a handled pipeline error, because the
# diagnostic/reporting package is the intended output for troubleshooting.

# ==============================================================================
# Analysis stage 3: 57_step04E1C_GSE250130_integration_QC_v1_1_fixed
# ==============================================================================

# ======================================================================
# 57_step04E1C_GSE250130_integration_QC_v1_1_fixed.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Perform lightweight integration QC for GSE250130.
# - Use CSV/status outputs from STEP04E1B only.
# - Do NOT load the 16.59 GB Seurat checkpoint.
# - Evaluate:
#     * locked cell/sample counts
#     * UMAP coordinate integrity
#     * cluster sizes across resolutions
#     * sample contribution and mixing
#     * sample-dominant and rare clusters
#     * group/timepoint composition when available in locked design
# - Generate QC tables, source data, PNG/PDF figures, Gate and Card.
#
# GOVERNANCE
# - This is QC, not a biological result-lock step.
# - Outputs go to module 08_qc.
# - Figures are QC figures, not manuscript figures.
# - No file is written to 06_locked_results.
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <- "57_step04E1C_GSE250130_integration_QC_v1_1_fixed.R"
script_version <- "STEP04E1C_GSE250130_INTEGRATION_QC_v1.1_FIXED"
input_version <- "SCI_v1_1_memory_safe"
qc_version <- "GSE250130_INTEGRATION_QC_v1.1_FIXED"

set.seed(20260805)

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1C_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

normalized_entropy <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) <= 1L) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p)) / log(length(p))
}

effective_sample_n <- function(counts) {
  counts <- counts[counts > 0]
  if (!length(counts)) return(0)
  p <- counts / sum(counts)
  exp(-sum(p * log(p)))
}

detect_first_column <- function(column_names, candidates) {
  lower <- tolower(column_names)
  candidate_lower <- tolower(candidates)

  for (candidate in candidate_lower) {
    hit <- which(lower == candidate)
    if (length(hit)) return(column_names[hit[1]])
  }

  ""
}

save_plot_pair <- function(plot_object, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")
  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }
  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages and directories
# ----------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them first with:\ninstall.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))"
  )
}

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(step04_root, "05_execution")
results_root <- file.path(step04_root, "06_results")
qc_base_root <- file.path(step04_root, "08_qc")
parameter_root <- file.path(step04_root, "03_parameters")
environment_root <- file.path(step04_root, "04_environment")

input_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version
)

input_table_root <- file.path(
  input_root,
  "tables",
  dataset_id
)

input_checkpoint_root <- file.path(
  input_root,
  "checkpoints",
  dataset_id
)

locked_design_path <- file.path(
  results_root,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0",
  "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"
)

doublet_qc_path <- file.path(
  results_root,
  "S4E1_sc",
  "SCI_v1_0",
  "tables",
  dataset_id,
  paste0(dataset_id, "_sample_doublet_qc.csv")
)

umap_path <- file.path(
  input_table_root,
  paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv")
)

cluster_sample_path <- file.path(
  input_table_root,
  paste0(dataset_id, "_cluster_by_sample_prejoin.csv")
)

checkpoint_status_path <- file.path(
  input_checkpoint_root,
  paste0(dataset_id, "_checkpoint_status.csv")
)

run_root <- file.path(
  qc_base_root,
  paste0(run_id, "_integration_QC")
)

table_root <- file.path(run_root, "tables")
figure_root <- file.path(run_root, "figures")
source_data_root <- file.path(run_root, "source_data")

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  run_root,
  table_root,
  figure_root,
  source_data_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)
desktop_candidates <- desktop_candidates[dir.exists(desktop_candidates)]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)


# ----------------------------------------------------------------------
# 3. Validate inputs
# ----------------------------------------------------------------------

required_inputs <- c(
  umap_path,
  cluster_sample_path,
  checkpoint_status_path,
  locked_design_path,
  doublet_qc_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

log_event("Starting", run_id)

umap <- data.table::fread(
  umap_path,
  data.table = FALSE,
  showProgress = FALSE
)

cluster_sample <- data.table::fread(
  cluster_sample_path,
  data.table = FALSE,
  showProgress = FALSE
)

checkpoint_status <- safe_read_csv(checkpoint_status_path)
locked_design <- safe_read_csv(locked_design_path)
doublet_qc <- safe_read_csv(doublet_qc_path)

required_umap_columns <- c(
  "cell_barcode",
  "UMAPH_1",
  "UMAPH_2",
  "sample_accession",
  "seurat_clusters",
  "RNA_snn_res.0.2",
  "RNA_snn_res.0.4",
  "RNA_snn_res.0.6"
)

missing_umap_columns <- setdiff(required_umap_columns, colnames(umap))

if (length(missing_umap_columns)) {
  stop(
    "UMAP metadata is missing columns: ",
    paste(missing_umap_columns, collapse = ", ")
  )
}

required_cluster_sample_columns <- c(
  "sample_accession",
  "seurat_clusters",
  "cell_n"
)

missing_cluster_sample_columns <- setdiff(
  required_cluster_sample_columns,
  colnames(cluster_sample)
)

if (length(missing_cluster_sample_columns)) {
  stop(
    "Cluster-by-sample table is missing columns: ",
    paste(missing_cluster_sample_columns, collapse = ", ")
  )
}

if (
  is.null(checkpoint_status) ||
    nrow(checkpoint_status) != 1L ||
    checkpoint_status$checkpoint_status[[1]] != "COMPLETE"
) {
  stop("Integrated checkpoint status is absent or not COMPLETE.")
}

if (
  is.null(locked_design) ||
    !"sample_accession" %in% colnames(locked_design)
) {
  stop("Locked GSE250130 design is absent or lacks sample_accession.")
}

if (
  is.null(doublet_qc) ||
    !"sample_accession" %in% colnames(doublet_qc)
) {
  stop("Doublet QC table is absent or lacks sample_accession.")
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  qc_version = qc_version,
  input_version = input_version,
  dataset_id = dataset_id,
  expected_sample_n = 28L,
  expected_cell_n = 236636L,
  expected_feature_n = 36601L,
  selected_resolution = "RNA_snn_res.0.4",
  expected_cluster_n_res_0_4 = 22L,
  rare_cluster_cell_threshold = 200L,
  low_sample_coverage_threshold = 14L,
  moderate_sample_dominance_threshold = 0.30,
  severe_sample_dominance_threshold = 0.50,
  low_normalized_entropy_threshold = 0.65,
  plot_max_cells_per_cluster = 5000L,
  random_seed = 20260805L,
  output_status = "QC_ONLY_NOT_LOCKED"
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)


# ----------------------------------------------------------------------
# 4. Merge locked design and identify design variables
# ----------------------------------------------------------------------

umap$sample_accession <- as.character(umap$sample_accession)
umap$seurat_clusters <- as.character(umap$seurat_clusters)

cluster_sample$sample_accession <- as.character(
  cluster_sample$sample_accession
)
cluster_sample$seurat_clusters <- as.character(
  cluster_sample$seurat_clusters
)

locked_design$sample_accession <- as.character(
  locked_design$sample_accession
)

design_duplicate_n <- sum(duplicated(locked_design$sample_accession))

if (design_duplicate_n > 0L) {
  stop(
    "Locked design contains duplicated sample_accession rows: ",
    design_duplicate_n
  )
}

design_columns <- setdiff(
  colnames(locked_design),
  "sample_accession"
)

group_column <- detect_first_column(
  colnames(locked_design),
  c(
    "group",
    "clinical_group",
    "phenotype",
    "disease_status",
    "condition",
    "fertility_status",
    "patient_group",
    "case_control_status",
    "case_control"
  )
)

time_column <- detect_first_column(
  colnames(locked_design),
  c(
    "timepoint",
    "time_point",
    "temporal_label",
    "cycle_day",
    "lh_day",
    "day",
    "phase",
    "sampling_time"
  )
)

role_column <- detect_first_column(
  colnames(locked_design),
  c(
    "analysis_role",
    "cohort_role",
    "sample_role",
    "role"
  )
)

donor_column <- detect_first_column(
  colnames(locked_design),
  c(
    "donor_id",
    "subject_id",
    "patient_id",
    "individual_id"
  )
)

design_column_inventory <- data.frame(
  column = colnames(locked_design),
  detected_role = ifelse(
    colnames(locked_design) == group_column,
    "GROUP",
    ifelse(
      colnames(locked_design) == time_column,
      "TIME",
      ifelse(
        colnames(locked_design) == role_column,
        "ANALYSIS_ROLE",
        ifelse(
          colnames(locked_design) == donor_column,
          "DONOR",
          "OTHER"
        )
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  design_column_inventory,
  file.path(table_root, "01_LOCKED_DESIGN_COLUMN_INVENTORY.csv")
)

umap_design <- merge(
  umap,
  locked_design,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

unmatched_umap_sample_n <- length(
  setdiff(
    unique(umap$sample_accession),
    locked_design$sample_accession
  )
)

cluster_sample_design <- merge(
  cluster_sample,
  locked_design,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)


# ----------------------------------------------------------------------
# 5. Global QC and resolution summaries
# ----------------------------------------------------------------------

resolution_columns <- c(
  "RNA_snn_res.0.2",
  "RNA_snn_res.0.4",
  "RNA_snn_res.0.6"
)

resolution_summary <- do.call(
  rbind,
  lapply(resolution_columns, function(column_name) {
    tab <- table(umap[[column_name]])

    data.frame(
      resolution_column = column_name,
      cluster_n = length(tab),
      min_cluster_cell_n = min(as.numeric(tab)),
      median_cluster_cell_n = stats::median(as.numeric(tab)),
      max_cluster_cell_n = max(as.numeric(tab)),
      clusters_below_200_cells = sum(as.numeric(tab) < 200L),
      stringsAsFactors = FALSE
    )
  })
)

write_csv_utf8(
  resolution_summary,
  file.path(table_root, "02_RESOLUTION_SUMMARY.csv")
)

finite_umap_fraction <- mean(
  is.finite(umap$UMAPH_1) &
    is.finite(umap$UMAPH_2)
)

global_qc <- data.frame(
  metric = c(
    "checkpoint_status",
    "checkpoint_join_policy",
    "cell_n",
    "sample_n",
    "feature_n",
    "cluster_n_resolution_0_4",
    "finite_UMAP_fraction",
    "locked_design_sample_n",
    "unmatched_UMAP_sample_n",
    "doublet_QC_sample_n",
    "design_duplicate_sample_n",
    "group_column_detected",
    "time_column_detected",
    "role_column_detected",
    "donor_column_detected"
  ),
  observed = c(
    checkpoint_status$checkpoint_status[[1]],
    checkpoint_status$join_policy[[1]],
    nrow(umap),
    length(unique(umap$sample_accession)),
    checkpoint_status$feature_n[[1]],
    length(unique(umap$seurat_clusters)),
    finite_umap_fraction,
    nrow(locked_design),
    unmatched_umap_sample_n,
    nrow(doublet_qc),
    design_duplicate_n,
    group_column,
    time_column,
    role_column,
    donor_column
  ),
  expected_or_interpretation = c(
    "COMPLETE",
    "DATA_ONLY_COUNTS_REMAIN_SPLIT",
    "236636",
    "28",
    "36601",
    "22",
    "1.0",
    "28",
    "0",
    "28",
    "0",
    "Detected when available",
    "Detected when available",
    "Detected when available",
    "Detected when available"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  global_qc,
  file.path(table_root, "03_GLOBAL_INTEGRATION_QC.csv")
)


# ----------------------------------------------------------------------
# 6. Cluster mixing and sample contribution metrics
# ----------------------------------------------------------------------

cluster_sample_dt <- data.table::as.data.table(cluster_sample)

# data.table requires every group to return identical column classes.
# Explicitly coerce grouped numeric summaries to double so that median/min/max
# cannot switch between integer and double across clusters.
cluster_sample_dt[
  ,
  cell_n := as.numeric(cell_n)
]

cluster_mixing <- cluster_sample_dt[
  ,
  {
    counts <- as.numeric(cell_n)
    total <- as.numeric(sum(counts, na.rm = TRUE))
    proportions <- counts / total
    max_index <- which.max(proportions)

    list(
      cell_n = as.numeric(total),
      sample_n_present = as.integer(
        data.table::uniqueN(sample_accession)
      ),
      min_sample_cell_n = as.numeric(
        min(counts, na.rm = TRUE)
      ),
      median_sample_cell_n = as.numeric(
        stats::median(counts, na.rm = TRUE)
      ),
      max_sample_cell_n = as.numeric(
        max(counts, na.rm = TRUE)
      ),
      max_sample_accession = as.character(
        sample_accession[max_index]
      ),
      max_sample_fraction = as.numeric(
        proportions[max_index]
      ),
      normalized_sample_entropy = as.numeric(
        normalized_entropy(counts)
      ),
      effective_sample_n = as.numeric(
        effective_sample_n(counts)
      )
    )
  },
  by = seurat_clusters
]

cluster_mixing$cluster_numeric <- suppressWarnings(
  as.integer(cluster_mixing$seurat_clusters)
)

cluster_mixing <- cluster_mixing[
  order(cluster_mixing$cluster_numeric),
]

cluster_mixing$rare_cluster_flag <- (
  cluster_mixing$cell_n < 200L
)

cluster_mixing$low_sample_coverage_flag <- (
  cluster_mixing$sample_n_present < 14L
)

cluster_mixing$moderate_sample_dominance_flag <- (
  cluster_mixing$max_sample_fraction > 0.30
)

cluster_mixing$severe_sample_dominance_flag <- (
  cluster_mixing$max_sample_fraction > 0.50
)

cluster_mixing$low_entropy_flag <- (
  cluster_mixing$normalized_sample_entropy < 0.65
)

cluster_mixing$review_flag <- (
  cluster_mixing$rare_cluster_flag |
    cluster_mixing$low_sample_coverage_flag |
    cluster_mixing$severe_sample_dominance_flag |
    cluster_mixing$low_entropy_flag
)

cluster_mixing$review_reason <- vapply(
  seq_len(nrow(cluster_mixing)),
  function(i) {
    reasons <- character()

    if (cluster_mixing$rare_cluster_flag[i]) {
      reasons <- c(reasons, "RARE_CLUSTER_LT_200_CELLS")
    }

    if (cluster_mixing$low_sample_coverage_flag[i]) {
      reasons <- c(reasons, "PRESENT_IN_LT_14_SAMPLES")
    }

    if (cluster_mixing$severe_sample_dominance_flag[i]) {
      reasons <- c(reasons, "ONE_SAMPLE_GT_50_PERCENT")
    }

    if (cluster_mixing$low_entropy_flag[i]) {
      reasons <- c(reasons, "NORMALIZED_ENTROPY_LT_0.65")
    }

    collapse_unique(reasons)
  },
  character(1)
)

write_csv_utf8(
  as.data.frame(cluster_mixing),
  file.path(table_root, "04_CLUSTER_SAMPLE_MIXING_METRICS.csv")
)

cluster_total_lookup <- setNames(
  cluster_mixing$cell_n,
  cluster_mixing$seurat_clusters
)

cluster_sample$cluster_fraction <- (
  cluster_sample$cell_n /
    cluster_total_lookup[cluster_sample$seurat_clusters]
)

sample_total_lookup <- tapply(
  cluster_sample$cell_n,
  cluster_sample$sample_accession,
  sum
)

cluster_sample$within_sample_fraction <- (
  cluster_sample$cell_n /
    sample_total_lookup[cluster_sample$sample_accession]
)

write_csv_utf8(
  cluster_sample,
  file.path(source_data_root, "SOURCE_DATA_cluster_by_sample.csv")
)


# ----------------------------------------------------------------------
# 7. Group/time composition summaries when design variables exist
# ----------------------------------------------------------------------

group_time_outputs <- list()

summarize_design_variable <- function(variable_name, output_stem) {
  if (!nzchar(variable_name)) return(NULL)

  dt <- data.table::as.data.table(cluster_sample_design)
  dt[, design_value := as.character(get(variable_name))]

  out <- dt[
    ,
    .(
      cell_n = as.numeric(
        sum(as.numeric(cell_n), na.rm = TRUE)
      ),
      sample_n = as.integer(
        data.table::uniqueN(sample_accession)
      )
    ),
    by = .(
      seurat_clusters,
      design_value
    )
  ]

  out[
    ,
    cluster_fraction := cell_n / sum(cell_n),
    by = seurat_clusters
  ]

  write_csv_utf8(
    as.data.frame(out),
    file.path(
      table_root,
      paste0(output_stem, ".csv")
    )
  )

  out
}

group_composition <- summarize_design_variable(
  group_column,
  "05_CLUSTER_BY_GROUP_COMPOSITION"
)

time_composition <- summarize_design_variable(
  time_column,
  "06_CLUSTER_BY_TIME_COMPOSITION"
)

role_composition <- summarize_design_variable(
  role_column,
  "07_CLUSTER_BY_ANALYSIS_ROLE_COMPOSITION"
)


# ----------------------------------------------------------------------
# 8. Generate QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey95"),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8)
  )

# Deterministic plotting subset: up to 5000 cells per cluster.
umap_dt <- data.table::as.data.table(umap_design)

umap_plot_data <- umap_dt[
  ,
  {
    n_keep <- min(.N, 5000L)
    .SD[sample.int(.N, n_keep)]
  },
  by = seurat_clusters
]

umap_plot_data$seurat_clusters <- factor(
  umap_plot_data$seurat_clusters,
  levels = as.character(
    sort(
      unique(
        suppressWarnings(
          as.integer(umap_plot_data$seurat_clusters)
        )
      )
    )
  )
)

cluster_centroids <- umap_dt[
  ,
  .(
    UMAPH_1 = stats::median(UMAPH_1),
    UMAPH_2 = stats::median(UMAPH_2)
  ),
  by = seurat_clusters
]

plot_umap_cluster <- ggplot2::ggplot(
  umap_plot_data,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = seurat_clusters
  )
) +
  ggplot2::geom_point(
    size = 0.08,
    alpha = 0.55,
    stroke = 0
  ) +
  ggplot2::geom_text(
    data = cluster_centroids,
    ggplot2::aes(
      x = UMAPH_1,
      y = UMAPH_2,
      label = seurat_clusters
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  ggplot2::labs(
    title = "GSE250130 integration QC: Harmony UMAP",
    subtitle = "Resolution 0.4; deterministic plotting subset",
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_qc +
  ggplot2::theme(
    legend.position = "right"
  )

save_plot_pair(
  plot_umap_cluster,
  "QC_Figure_1_UMAP_clusters_res0_4",
  width = 8.5,
  height = 7
)

plot_cluster_size <- ggplot2::ggplot(
  as.data.frame(cluster_mixing),
  ggplot2::aes(
    x = reorder(
      factor(seurat_clusters),
      cluster_numeric
    ),
    y = cell_n,
    fill = review_flag
  )
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Cluster sizes at resolution 0.4",
    x = "Cluster",
    y = "Cell count",
    fill = "Review flag"
  ) +
  theme_qc

save_plot_pair(
  plot_cluster_size,
  "QC_Figure_2_cluster_sizes",
  width = 8,
  height = 5
)

plot_mixing <- ggplot2::ggplot(
  as.data.frame(cluster_mixing),
  ggplot2::aes(
    x = max_sample_fraction,
    y = normalized_sample_entropy,
    size = cell_n,
    label = seurat_clusters
  )
) +
  ggplot2::geom_vline(
    xintercept = c(0.30, 0.50),
    linetype = c("dashed", "solid")
  ) +
  ggplot2::geom_hline(
    yintercept = 0.65,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = review_flag),
    alpha = 0.8
  ) +
  ggplot2::geom_text(
    nudge_y = 0.025,
    size = 3
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_size_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Cluster sample-mixing QC",
    subtitle = "Lower dominance and higher entropy indicate broader sample mixing",
    x = "Largest single-sample fraction",
    y = "Normalized sample entropy",
    size = "Cell count",
    shape = "Review flag"
  ) +
  theme_qc

save_plot_pair(
  plot_mixing,
  "QC_Figure_3_cluster_sample_mixing",
  width = 8,
  height = 6
)

heatmap_data <- cluster_sample
heatmap_data$seurat_clusters <- factor(
  heatmap_data$seurat_clusters,
  levels = as.character(
    sort(
      unique(
        suppressWarnings(
          as.integer(heatmap_data$seurat_clusters)
        )
      )
    )
  )
)

sample_order <- names(
  sort(
    tapply(
      heatmap_data$cell_n,
      heatmap_data$sample_accession,
      sum
    ),
    decreasing = TRUE
  )
)

heatmap_data$sample_accession <- factor(
  heatmap_data$sample_accession,
  levels = rev(sample_order)
)

plot_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = seurat_clusters,
    y = sample_accession,
    fill = within_sample_fraction
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    labels = scales::percent_format(accuracy = 1),
    option = "C"
  ) +
  ggplot2::labs(
    title = "Within-sample cluster composition",
    x = "Cluster",
    y = "Sample",
    fill = "Fraction"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6)
  )

save_plot_pair(
  plot_heatmap,
  "QC_Figure_4_sample_cluster_composition_heatmap",
  width = 9,
  height = 8
)

resolution_long <- do.call(
  rbind,
  lapply(resolution_columns, function(column_name) {
    tab <- as.data.frame(table(umap[[column_name]]))
    colnames(tab) <- c("cluster", "cell_n")
    tab$resolution <- sub(
      "RNA_snn_res.",
      "",
      column_name,
      fixed = TRUE
    )
    tab
  })
)

resolution_long$cluster <- as.character(
  resolution_long$cluster
)

write_csv_utf8(
  resolution_long,
  file.path(source_data_root, "SOURCE_DATA_resolution_cluster_sizes.csv")
)

plot_resolution <- ggplot2::ggplot(
  resolution_long,
  ggplot2::aes(
    x = reorder(cluster, cell_n),
    y = cell_n
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(
    ~ resolution,
    scales = "free_x",
    nrow = 1
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Cluster-size distribution across tested resolutions",
    x = "Cluster",
    y = "Cell count"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    )
  )

save_plot_pair(
  plot_resolution,
  "QC_Figure_5_resolution_cluster_sizes",
  width = 11,
  height = 5
)


# ----------------------------------------------------------------------
# 9. Gate and interpretation
# ----------------------------------------------------------------------

technical_checks <- data.frame(
  check = c(
    "checkpoint_complete",
    "join_policy_correct",
    "cell_count_locked",
    "sample_count_locked",
    "feature_count_locked",
    "cluster_count_res0_4_expected",
    "all_UMAP_coordinates_finite",
    "all_UMAP_samples_in_locked_design",
    "locked_design_has_28_samples",
    "doublet_QC_has_28_samples"
  ),
  pass = c(
    checkpoint_status$checkpoint_status[[1]] == "COMPLETE",
    checkpoint_status$join_policy[[1]] ==
      "DATA_ONLY_COUNTS_REMAIN_SPLIT",
    nrow(umap) == 236636L,
    length(unique(umap$sample_accession)) == 28L,
    checkpoint_status$feature_n[[1]] == 36601L,
    length(unique(umap$seurat_clusters)) == 22L,
    identical(finite_umap_fraction, 1),
    unmatched_umap_sample_n == 0L,
    nrow(locked_design) == 28L,
    nrow(doublet_qc) == 28L
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  technical_checks,
  file.path(table_root, "08_TECHNICAL_GATE_CHECKS.csv")
)

technical_pass <- all(technical_checks$pass)
review_clusters <- cluster_mixing$seurat_clusters[
  cluster_mixing$review_flag
]

gate_status <- if (!technical_pass) {
  "HOLD_TECHNICAL_INTEGRATION_QC_FAILED"
} else if (length(review_clusters)) {
  "PASS_WITH_REVIEW_FLAGS"
} else {
  "PASS"
}

gate_table <- data.frame(
  dataset_id = dataset_id,
  qc_version = qc_version,
  input_version = input_version,
  gate_status = gate_status,
  technical_check_n = nrow(technical_checks),
  technical_pass_n = sum(technical_checks$pass),
  technical_fail_n = sum(!technical_checks$pass),
  cluster_n = nrow(cluster_mixing),
  review_cluster_n = length(review_clusters),
  review_clusters = collapse_unique(review_clusters),
  group_column = group_column,
  time_column = time_column,
  role_column = role_column,
  donor_column = donor_column,
  marker_analysis_authorized = technical_pass,
  biological_annotation_locked = FALSE,
  result_lock_authorized = FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_root,
  "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv"
)

write_csv_utf8(gate_table, gate_path)

review_detail <- cluster_mixing[
  cluster_mixing$review_flag,
  c(
    "seurat_clusters",
    "cell_n",
    "sample_n_present",
    "max_sample_accession",
    "max_sample_fraction",
    "normalized_sample_entropy",
    "review_reason"
  )
]

interpretation_lines <- c(
  "# STEP04E1C GSE250130 integration QC",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  paste0("- Cells: ", format(nrow(umap), big.mark = ",")),
  paste0("- Samples: ", length(unique(umap$sample_accession))),
  paste0("- Resolution 0.4 clusters: ", length(unique(umap$seurat_clusters))),
  paste0("- Finite UMAP fraction: ", finite_umap_fraction),
  paste0(
    "- Review clusters: ",
    ifelse(
      length(review_clusters),
      collapse_unique(review_clusters),
      "none"
    )
  ),
  "",
  "## Interpretation",
  "",
  if (technical_pass) {
    "- The GSE250130 integrated object passes the predefined technical identity and completeness checks."
  } else {
    "- One or more technical identity/completeness checks failed; marker analysis must not proceed."
  },
  if (length(review_clusters)) {
    "- Review flags identify clusters that are rare, present in few samples, dominated by one sample, or show low normalized sample entropy."
  } else {
    "- No cluster crossed the predefined review thresholds."
  },
  "- A review flag does not prove a batch artifact. These clusters must be evaluated using marker genes, sample metadata, and biological plausibility.",
  "- The current Gate authorizes marker analysis when technical checks pass, but it does not authorize final annotation or result locking.",
  "",
  "## Required next analysis",
  "",
  "1. Compute cluster markers using the joined normalized data layer.",
  "2. Apply marker-based broad-cell annotation with explicit marker coverage.",
  "3. Review flagged clusters before accepting labels.",
  "4. Only after annotation QC, calculate VDR program localization and donor-aware pseudobulk."
)

interpretation_path <- file.path(
  run_root,
  "STEP04E1C_GSE250130_QC_INTERPRETATION.md"
)

write_text(interpretation_lines, interpretation_path)

card_lines <- c(
  "# STEP04E1C Result Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Input: ", input_version),
  paste0("- QC version: ", qc_version),
  paste0("- Status: ", gate_status),
  paste0("- Marker analysis authorized: ", technical_pass),
  "- Biological annotation locked: FALSE",
  "- Result lock authorized: FALSE",
  paste0(
    "- Review clusters: ",
    ifelse(
      length(review_clusters),
      collapse_unique(review_clusters),
      "none"
    )
  ),
  "",
  "This card documents integration QC only. It is not a scientific result lock."
)

card_path <- file.path(
  run_root,
  "STEP04E1C_GSE250130_RESULT_CARD.md"
)

write_text(card_lines, card_path)


# ----------------------------------------------------------------------
# 10. Run manifest, sessionInfo, artifact manifest
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_table <- data.frame(
  package = c(
    "data.table",
    "ggplot2",
    "jsonlite",
    "scales"
  ),
  version = vapply(
    c(
      "data.table",
      "ggplot2",
      "jsonlite",
      "scales"
    ),
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(run_root, "PACKAGE_VERSIONS.csv")
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- list.files(
  run_root,
  recursive = TRUE,
  full.names = TRUE
)

artifact_manifest <- data.frame(
  relative_path = normalize_slashes(
    substring(
      artifact_files,
      nchar(run_root) + 2L
    )
  ),
  size_bytes = file.info(artifact_files)$size,
  md5 = unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(run_root, "ARTIFACT_MANIFEST.csv")
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  qc_version = qc_version,
  input_version = input_version,
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_seconds = round(
    as.numeric(difftime(run_end, run_start, units = "secs")),
    2
  ),
  gate_status = gate_status,
  large_Seurat_checkpoint_loaded = FALSE,
  project_analysis_modified = FALSE,
  marker_analysis_authorized = technical_pass,
  result_lock_authorized = FALSE,
  output_root = normalize_slashes(run_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(run_manifest, run_manifest_path)


# ----------------------------------------------------------------------
# 11. Small upload package
# ----------------------------------------------------------------------

small_handoff_root <- file.path(
  run_root,
  "small_handoff"
)

dir.create(
  small_handoff_root,
  recursive = TRUE,
  showWarnings = FALSE
)

small_files <- unique(c(
  list.files(
    run_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <= 100 * 1024^2 &
    !startsWith(
      normalize_slashes(small_files),
      normalize_slashes(small_handoff_root)
    )
]

copy_manifest <- list()

for (i in seq_along(small_files)) {
  source_path <- small_files[i]
  destination_path <- file.path(
    small_handoff_root,
    paste0(
      sprintf("%03d", i),
      "_",
      basename(source_path)
    )
  )

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  copy_manifest[[i]] <- data.frame(
    source = normalize_slashes(source_path),
    destination = normalize_slashes(destination_path),
    copied = isTRUE(copied),
    stringsAsFactors = FALSE
  )
}

if (length(copy_manifest)) {
  write_csv_utf8(
    do.call(rbind, copy_manifest),
    file.path(run_root, "SMALL_HANDOFF_COPY_MANIFEST.csv")
  )
}

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = c(
          gate_path,
          card_path,
          interpretation_path,
          table_root,
          figure_root,
          source_data_root,
          log_path,
          run_manifest_path,
          parameter_path,
          session_path,
          file.path(run_root, "PACKAGE_VERSIONS.csv"),
          file.path(run_root, "ARTIFACT_MANIFEST.csv"),
          file.path(run_root, "SMALL_HANDOFF_COPY_MANIFEST.csv")
        ),
        recurse = TRUE,
        include_directories = TRUE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "review_clusters=",
  collapse_unique(review_clusters)
)

cat("\n============================================================\n")
cat("STEP04E1C GSE250130 integration QC completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat("QC folder: ", run_root, "\n", sep = "")
cat("Run manifest: ", run_manifest_path, "\n", sep = "")
cat(
  "Review clusters: ",
  ifelse(
    length(review_clusters),
    collapse_unique(review_clusters),
    "none"
  ),
  "\n",
  sep = ""
)

if (zip_ok) {
  cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe QC folder is complete; compress it manually if needed.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 4: 58_step04E1D_GSE250130_cluster_markers_v1_1_memory_safe
# ==============================================================================

# ======================================================================
# 58_step04E1D_GSE250130_cluster_markers_v1_1_memory_safe.R
# Nutrition_Endometrial_VDR_v2
#
# SUPERSEDES
# - 58_step04E1D_GSE250130_cluster_markers_v1.R
# - MARKERS_v1_0 (retained as a failed-run record)
#
# FAILURE BEING REPAIRED
# - Seurat::FindMarkers repeatedly attempted an internal rowSums operation
#   on the complete 236,636-cell expression layer.
# - Each cluster failed while requesting an additional 2.0-2.6 GB vector.
#
# PURPOSE
# - Extract annotation-oriented cluster markers without FindMarkers.
# - Explicitly select a small, balanced target and background cell set.
# - Pull only those cells from the joined sparse normalized-data layer.
# - Calculate descriptive marker evidence:
#     * mean log-normalized expression difference
#     * target and background detection fractions
#     * detection-fraction difference
#     * target-sample breadth
#     * composite marker score
# - Process and save one cluster at a time; safe to resume.
#
# SCIENTIFIC SCOPE
# - Marker scores are descriptive evidence for cell-identity annotation.
# - They are NOT donor-level differential-expression statistics.
# - No p-values are reported or interpreted.
# - Clinical/timepoint inference remains reserved for donor-aware pseudobulk.
#
# INPUT
# - SCI_v1_1_memory_safe integrated checkpoint
# - STEP04E1C integration QC Gate
#
# OUTPUT
# - New immutable candidate folder:
#     MARKERS_v1_1_memory_safe
# - No output is written to 06_locked_results.
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)
options(timeout = max(14400L, getOption("timeout")))
options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260805)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <-
  "58_step04E1D_GSE250130_cluster_markers_v1_1_memory_safe.R"
script_version <-
  "STEP04E1D_GSE250130_CLUSTER_MARKERS_v1.1_MEMORY_SAFE"
input_version <- "SCI_v1_1_memory_safe"
marker_version <- "MARKERS_v1_1_memory_safe"
cluster_column <- "RNA_snn_res.0.4"

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1D_",
  dataset_id,
  "_v1_1_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

format_bytes <- function(bytes) {
  if (length(bytes) == 0L || is.na(bytes)) return("")
  if (bytes < 1024) return(sprintf("%.0f B", bytes))
  if (bytes < 1024^2) return(sprintf("%.1f KB", bytes / 1024))
  if (bytes < 1024^3) return(sprintf("%.1f MB", bytes / 1024^2))
  sprintf("%.2f GB", bytes / 1024^3)
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") return(NA_real_)

  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"$o=Get-CimInstance Win32_OperatingSystem; ',
        '[math]::Round($o.FreePhysicalMemory/1MB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )

  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

get_free_disk_gb <- function(drive_letter = "C") {
  if (.Platform$OS.type != "windows") return(NA_real_)

  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"[math]::Round((Get-PSDrive -Name ',
        drive_letter,
        ').Free/1GB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )

  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

numeric_cluster_order <- function(x) {
  nx <- suppressWarnings(as.integer(as.character(x)))
  order(nx, as.character(x), na.last = TRUE)
}

save_plot_pair <- function(plot_object, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

# Balanced sampling retains representation across sample_accession.
balanced_sample_cells <- function(
  cells,
  sample_by_cell,
  max_total
) {
  cells <- unique(as.character(cells))
  cells <- cells[cells %in% names(sample_by_cell)]

  if (!length(cells)) return(character())
  if (length(cells) <= max_total) return(cells)

  sample_values <- as.character(sample_by_cell[cells])
  groups <- split(cells, sample_values, drop = TRUE)

  group_n <- length(groups)
  initial_quota <- max(1L, floor(max_total / group_n))

  chosen <- unlist(
    lapply(
      groups,
      function(one_group) {
        if (length(one_group) <= initial_quota) {
          one_group
        } else {
          sample(
            one_group,
            size = initial_quota,
            replace = FALSE
          )
        }
      }
    ),
    use.names = FALSE
  )

  remaining_n <- max_total - length(chosen)

  if (remaining_n > 0L) {
    remaining_pool <- setdiff(cells, chosen)

    if (length(remaining_pool)) {
      extra <- sample(
        remaining_pool,
        size = min(remaining_n, length(remaining_pool)),
        replace = FALSE
      )
      chosen <- c(chosen, extra)
    }
  }

  if (length(chosen) > max_total) {
    chosen <- sample(
      chosen,
      size = max_total,
      replace = FALSE
    )
  }

  unique(chosen)
}

# Construct target sample-breadth metrics without dense cell-level matrices.
calculate_sample_breadth <- function(
  target_matrix,
  target_samples,
  detection_threshold = 0.10
) {
  target_samples <- factor(as.character(target_samples))
  sample_n <- nlevels(target_samples)

  if (sample_n == 0L) {
    return(
      list(
        breadth_fraction = rep(NA_real_, nrow(target_matrix)),
        sample_n = 0L
      )
    )
  }

  membership <- Matrix::sparseMatrix(
    i = seq_along(target_samples),
    j = as.integer(target_samples),
    x = 1,
    dims = c(
      length(target_samples),
      sample_n
    )
  )

  detected <- target_matrix > 0
  detection_counts <- detected %*% membership

  cells_per_sample <- tabulate(
    as.integer(target_samples),
    nbins = sample_n
  )

  detection_fraction <- as.matrix(detection_counts)

  detection_fraction <- sweep(
    detection_fraction,
    MARGIN = 2,
    STATS = cells_per_sample,
    FUN = "/"
  )

  breadth_fraction <- rowMeans(
    detection_fraction >= detection_threshold
  )

  list(
    breadth_fraction = as.numeric(breadth_fraction),
    sample_n = sample_n
  )
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")

  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages and directories
# ----------------------------------------------------------------------

required_packages <- c(
  "SeuratObject",
  "Matrix",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(step04_root, "05_execution")
results_root <- file.path(step04_root, "06_results")
qc_base_root <- file.path(step04_root, "08_qc")
parameter_root <- file.path(step04_root, "03_parameters")
environment_root <- file.path(step04_root, "04_environment")
manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

checkpoint_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "checkpoints",
  dataset_id
)

checkpoint_path <- file.path(
  checkpoint_root,
  paste0(dataset_id, "_integrated_checkpoint.rds")
)

checkpoint_status_path <- file.path(
  checkpoint_root,
  paste0(dataset_id, "_checkpoint_status.csv")
)

marker_result_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "markers",
  dataset_id,
  marker_version
)

per_cluster_root <- file.path(
  marker_result_root,
  "per_cluster"
)

table_root <- file.path(
  marker_result_root,
  "tables"
)

source_data_root <- file.path(
  marker_result_root,
  "source_data"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(run_id, "_cluster_markers_memory_safe")
)

figure_root <- file.path(
  run_qc_root,
  "figures"
)

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  marker_result_root,
  per_cluster_root,
  table_root,
  source_data_root,
  run_qc_root,
  figure_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)


# ----------------------------------------------------------------------
# 3. Validate checkpoint and integration Gate
# ----------------------------------------------------------------------

if (!file.exists(checkpoint_path)) {
  stop("Integrated checkpoint not found: ", checkpoint_path)
}

checkpoint_status <- safe_read_csv(checkpoint_status_path)

if (
  is.null(checkpoint_status) ||
    nrow(checkpoint_status) != 1L ||
    checkpoint_status$checkpoint_status[[1]] != "COMPLETE" ||
    checkpoint_status$join_policy[[1]] !=
      "DATA_ONLY_COUNTS_REMAIN_SPLIT"
) {
  stop("Integrated checkpoint status is absent or invalid.")
}

qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

qc_dirs <- qc_dirs[
  grepl(
    "^STEP04E1C_GSE250130_[0-9]{8}_[0-9]{6}_integration_QC$",
    basename(qc_dirs),
    perl = TRUE
  )
]

if (!length(qc_dirs)) {
  stop("No STEP04E1C integration-QC folder found.")
}

latest_qc_dir <- qc_dirs[
  which.max(file.info(qc_dirs)$mtime)
]

integration_gate_path <- file.path(
  latest_qc_dir,
  "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv"
)

integration_gate <- safe_read_csv(integration_gate_path)

if (
  is.null(integration_gate) ||
    nrow(integration_gate) != 1L ||
    !integration_gate$gate_status[[1]] %in%
      c("PASS", "PASS_WITH_REVIEW_FLAGS") ||
    !isTRUE(integration_gate$marker_analysis_authorized[[1]])
) {
  stop("STEP04E1C Gate does not authorize marker analysis.")
}

review_clusters <- character()

if (
  "review_clusters" %in% colnames(integration_gate) &&
    !is.na(integration_gate$review_clusters[[1]]) &&
    nzchar(integration_gate$review_clusters[[1]])
) {
  review_clusters <- trimws(
    unlist(
      strsplit(
        integration_gate$review_clusters[[1]],
        ";",
        fixed = TRUE
      )
    )
  )
}

free_memory_gb <- get_free_memory_gb()
free_disk_gb <- get_free_disk_gb("C")

resource_snapshot <- data.frame(
  metric = c(
    "free_physical_memory_GB_before_load",
    "free_C_drive_GB_before_run",
    "checkpoint_size_GB",
    "checkpoint_md5",
    "integration_gate_status",
    "review_clusters",
    "failed_predecessor_marker_version"
  ),
  value = c(
    free_memory_gb,
    free_disk_gb,
    round(file.info(checkpoint_path)$size / 1024^3, 3),
    unname(tools::md5sum(checkpoint_path)),
    integration_gate$gate_status[[1]],
    collapse_unique(review_clusters),
    "MARKERS_v1_0_FINDMARKERS_MEMORY_FAILURE"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resource_snapshot,
  file.path(run_qc_root, "01_RESOURCE_AND_INPUT_SNAPSHOT.csv")
)

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
  stop(
    "Free physical memory is below 18 GB. ",
    "Restart RStudio and rerun in a clean session."
  )
}

if (is.finite(free_disk_gb) && free_disk_gb < 10) {
  stop("Free C: drive space is below 10 GB.")
}

final_marker_path <- file.path(
  table_root,
  paste0(dataset_id, "_all_cluster_markers_descriptive.csv")
)

final_gate_path <- file.path(
  run_qc_root,
  "STEP04E1D_GSE250130_MARKER_GATE.csv"
)

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  marker_version = marker_version,
  dataset_id = dataset_id,
  cluster_column = cluster_column,
  marker_method =
    "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING",
  target_max_cells = 1500L,
  background_max_cells_per_other_cluster = 150L,
  target_min_detection_fraction = 0.10,
  minimum_detection_fraction_difference = 0.05,
  minimum_mean_log_expression_difference = 0.20,
  minimum_target_sample_breadth = 0.25,
  fallback_target_min_detection_fraction = 0.05,
  fallback_detection_fraction_difference = 0.02,
  fallback_mean_log_expression_difference = 0.10,
  primary_minimum_marker_n = 25L,
  sample_breadth_detection_threshold = 0.10,
  canonical_panel_max_cells_per_cluster = 500L,
  p_values_reported = FALSE,
  inferential_role =
    "NONE_ANNOTATION_SUPPORT_ONLY",
  random_seed = 20260805L,
  integration_gate = integration_gate$gate_status[[1]],
  review_clusters = review_clusters,
  output_status =
    "CANDIDATE_MARKERS_FOR_ANNOTATION"
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

log_event(
  "Starting",
  run_id,
  "checkpoint=",
  checkpoint_path,
  "integration_gate=",
  integration_gate$gate_status[[1]],
  "review_clusters=",
  collapse_unique(review_clusters)
)


# ----------------------------------------------------------------------
# 4. Load and validate object
# ----------------------------------------------------------------------

load_start <- Sys.time()
object <- readRDS(checkpoint_path)

log_event(
  "Loaded checkpoint",
  "elapsed_min=",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        load_start,
        units = "mins"
      )
    ),
    2
  )
)

if (!inherits(object, "Seurat")) {
  stop("Checkpoint is not a Seurat object.")
}

if (
  ncol(object) != 236636L ||
    nrow(object) != 36601L
) {
  stop(
    "Unexpected object dimensions: ",
    nrow(object),
    " features x ",
    ncol(object),
    " cells"
  )
}

if (!cluster_column %in% colnames(object@meta.data)) {
  stop("Cluster column missing: ", cluster_column)
}

if (!"sample_accession" %in% colnames(object@meta.data)) {
  stop("sample_accession is missing from object metadata.")
}

SeuratObject::DefaultAssay(object) <- "RNA"

data_layers <- SeuratObject::Layers(
  object[["RNA"]],
  search = "^data$"
)

if (!identical(data_layers, "data")) {
  stop(
    "Expected exactly one normalized data layer named 'data'. Found: ",
    paste(data_layers, collapse = "; ")
  )
}

cell_metadata <- data.frame(
  cell_barcode = rownames(object@meta.data),
  sample_accession = as.character(
    object$sample_accession
  ),
  cluster = as.character(
    object@meta.data[[cluster_column]]
  ),
  stringsAsFactors = FALSE
)

rownames(cell_metadata) <- cell_metadata$cell_barcode

sample_by_cell <- setNames(
  cell_metadata$sample_accession,
  cell_metadata$cell_barcode
)

cells_by_cluster <- split(
  cell_metadata$cell_barcode,
  cell_metadata$cluster
)

cluster_levels <- names(cells_by_cluster)
cluster_levels <- cluster_levels[
  numeric_cluster_order(cluster_levels)
]

cluster_counts <- data.frame(
  cluster = cluster_levels,
  cell_n = vapply(
    cells_by_cluster[cluster_levels],
    length,
    integer(1)
  ),
  sample_n = vapply(
    cells_by_cluster[cluster_levels],
    function(cells) {
      length(unique(sample_by_cell[cells]))
    },
    integer(1)
  ),
  review_cluster = (
    cluster_levels %in% review_clusters
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  cluster_counts,
  file.path(table_root, "01_CLUSTER_CELL_COUNTS.csv")
)

if (
  length(cluster_levels) != 22L ||
    sum(cluster_counts$cell_n) != 236636L
) {
  stop("Cluster identity validation failed.")
}

layer_manifest <- data.frame(
  layer = SeuratObject::Layers(object[["RNA"]]),
  layer_type = ifelse(
    grepl(
      "^counts",
      SeuratObject::Layers(object[["RNA"]])
    ),
    "counts",
    ifelse(
      grepl(
        "^data",
        SeuratObject::Layers(object[["RNA"]])
      ),
      "data",
      "other"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  layer_manifest,
  file.path(run_qc_root, "02_LAYER_MANIFEST.csv")
)


# ----------------------------------------------------------------------
# 5. Per-cluster resumable descriptive marker scoring
# ----------------------------------------------------------------------

marker_status_rows <- list()

for (cluster_id in cluster_levels) {
  cluster_file <- file.path(
    per_cluster_root,
    paste0(
      dataset_id,
      "_cluster_",
      cluster_id,
      "_descriptive_markers.csv"
    )
  )

  status_file <- file.path(
    per_cluster_root,
    paste0(
      dataset_id,
      "_cluster_",
      cluster_id,
      "_marker_status.csv"
    )
  )

  selected_cells_file <- file.path(
    per_cluster_root,
    paste0(
      dataset_id,
      "_cluster_",
      cluster_id,
      "_selected_cells.csv"
    )
  )

  existing_status <- safe_read_csv(status_file)

  if (
    file.exists(cluster_file) &&
      !is.null(existing_status) &&
      nrow(existing_status) == 1L &&
      existing_status$status[[1]] == "COMPLETE"
  ) {
    log_event(
      "SKIP completed cluster",
      cluster_id
    )

    marker_status_rows[[length(marker_status_rows) + 1L]] <-
      existing_status

    next
  }

  cluster_start <- Sys.time()
  cluster_error <- ""

  target_pool <- cells_by_cluster[[cluster_id]]

  target_cells <- balanced_sample_cells(
    cells = target_pool,
    sample_by_cell = sample_by_cell,
    max_total = 1500L
  )

  other_clusters <- setdiff(
    cluster_levels,
    cluster_id
  )

  background_cells <- unlist(
    lapply(
      other_clusters,
      function(other_cluster) {
        balanced_sample_cells(
          cells = cells_by_cluster[[other_cluster]],
          sample_by_cell = sample_by_cell,
          max_total = 150L
        )
      }
    ),
    use.names = FALSE
  )

  background_cells <- unique(background_cells)
  background_cells <- setdiff(
    background_cells,
    target_cells
  )

  selected_cells <- c(
    target_cells,
    background_cells
  )

  selected_cell_manifest <- data.frame(
    cell_barcode = selected_cells,
    role = c(
      rep("TARGET", length(target_cells)),
      rep(
        "CLUSTER_BALANCED_BACKGROUND",
        length(background_cells)
      )
    ),
    sample_accession = sample_by_cell[selected_cells],
    cluster = cell_metadata[
      selected_cells,
      "cluster"
    ],
    target_cluster = cluster_id,
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    selected_cell_manifest,
    selected_cells_file
  )

  log_event(
    "START descriptive marker cluster",
    cluster_id,
    "target_pool_n=",
    length(target_pool),
    "target_selected_n=",
    length(target_cells),
    "background_selected_n=",
    length(background_cells)
  )

  marker_result <- tryCatch(
    {
      selected_matrix <- SeuratObject::LayerData(
        object = object,
        assay = "RNA",
        layer = "data",
        cells = selected_cells
      )

      if (!inherits(selected_matrix, "sparseMatrix")) {
        selected_matrix <- Matrix::Matrix(
          selected_matrix,
          sparse = TRUE
        )
      }

      target_index <- match(
        target_cells,
        colnames(selected_matrix)
      )

      background_index <- match(
        background_cells,
        colnames(selected_matrix)
      )

      if (
        anyNA(target_index) ||
          anyNA(background_index)
      ) {
        stop(
          "Selected cells were not all returned by LayerData."
        )
      }

      target_matrix <- selected_matrix[
        ,
        target_index,
        drop = FALSE
      ]

      background_matrix <- selected_matrix[
        ,
        background_index,
        drop = FALSE
      ]

      mean_target <- Matrix::rowMeans(target_matrix)
      mean_background <- Matrix::rowMeans(
        background_matrix
      )

      pct_target <- Matrix::rowMeans(
        target_matrix > 0
      )

      pct_background <- Matrix::rowMeans(
        background_matrix > 0
      )

      mean_log_expression_difference <- (
        mean_target - mean_background
      )

      detection_fraction_difference <- (
        pct_target - pct_background
      )

      log2_detection_ratio <- log2(
        (pct_target + 0.01) /
          (pct_background + 0.01)
      )

      breadth <- calculate_sample_breadth(
        target_matrix = target_matrix,
        target_samples = sample_by_cell[target_cells],
        detection_threshold = 0.10
      )

      marker_score <- (
        pmax(
          mean_log_expression_difference,
          0
        ) *
          pmax(
            detection_fraction_difference,
            0
          ) *
          (0.5 + breadth$breadth_fraction)
      )

      all_gene_scores <- data.frame(
        cluster = cluster_id,
        gene = rownames(selected_matrix),
        mean_log_expression_target =
          as.numeric(mean_target),
        mean_log_expression_background =
          as.numeric(mean_background),
        mean_log_expression_difference =
          as.numeric(
            mean_log_expression_difference
          ),
        pct_target = as.numeric(pct_target),
        pct_background = as.numeric(
          pct_background
        ),
        detection_fraction_difference =
          as.numeric(
            detection_fraction_difference
          ),
        log2_detection_ratio =
          as.numeric(log2_detection_ratio),
        target_sample_breadth_fraction =
          as.numeric(
            breadth$breadth_fraction
          ),
        target_sample_n =
          as.integer(breadth$sample_n),
        target_pool_cell_n =
          as.integer(length(target_pool)),
        target_selected_cell_n =
          as.integer(length(target_cells)),
        background_selected_cell_n =
          as.integer(length(background_cells)),
        marker_score = as.numeric(marker_score),
        review_cluster = (
          cluster_id %in% review_clusters
        ),
        stringsAsFactors = FALSE
      )

      primary_keep <- (
        all_gene_scores$pct_target >= 0.10 &
          all_gene_scores$
            detection_fraction_difference >= 0.05 &
          all_gene_scores$
            mean_log_expression_difference >= 0.20 &
          all_gene_scores$
            target_sample_breadth_fraction >= 0.25
      )

      relaxed_keep <- (
        all_gene_scores$pct_target >= 0.05 &
          all_gene_scores$
            detection_fraction_difference >= 0.02 &
          all_gene_scores$
            mean_log_expression_difference >= 0.10 &
          all_gene_scores$
            target_sample_breadth_fraction >= 0.20
      )

      all_gene_scores$selection_tier <- ifelse(
        primary_keep,
        "PRIMARY",
        ifelse(
          relaxed_keep,
          "RELAXED",
          "NOT_SELECTED"
        )
      )

      primary_markers <- all_gene_scores[
        primary_keep,
        ,
        drop = FALSE
      ]

      if (nrow(primary_markers) >= 25L) {
        selected_markers <- primary_markers
      } else {
        selected_markers <- all_gene_scores[
          relaxed_keep,
          ,
          drop = FALSE
        ]
      }

      selected_markers <- selected_markers[
        order(
          -selected_markers$marker_score,
          -selected_markers$
            mean_log_expression_difference,
          -selected_markers$
            detection_fraction_difference,
          -selected_markers$pct_target,
          selected_markers$gene
        ),
        ,
        drop = FALSE
      ]

      selected_markers$marker_role <-
        "CLUSTER_ANNOTATION_SUPPORT_ONLY"

      selected_markers$inferential_role <-
        "NONE_NO_P_VALUES_DONOR_INFERENCE_RESERVED_FOR_PSEUDOBULK"

      rm(
        selected_matrix,
        target_matrix,
        background_matrix,
        all_gene_scores,
        primary_markers
      )
      gc()

      selected_markers
    },
    error = function(e) {
      cluster_error <<- conditionMessage(e)
      NULL
    }
  )

  elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        cluster_start,
        units = "mins"
      )
    ),
    2
  )

  if (
    is.null(marker_result) ||
      !nrow(marker_result)
  ) {
    marker_status <- data.frame(
      cluster = cluster_id,
      status = ifelse(
        nzchar(cluster_error),
        "FAILED",
        "NO_MARKERS_RETURNED"
      ),
      marker_n = 0L,
      primary_marker_n = 0L,
      relaxed_marker_n = 0L,
      target_pool_cell_n = length(target_pool),
      target_selected_cell_n = length(target_cells),
      background_selected_cell_n =
        length(background_cells),
      target_sample_n = length(
        unique(sample_by_cell[target_cells])
      ),
      elapsed_minutes = elapsed_minutes,
      error_message = cluster_error,
      output_file =
        normalize_slashes(cluster_file),
      stringsAsFactors = FALSE
    )

    write_csv_utf8(marker_status, status_file)

    marker_status_rows[[length(marker_status_rows) + 1L]] <-
      marker_status

    log_event(
      marker_status$status[[1]],
      "cluster=",
      cluster_id,
      "error=",
      cluster_error
    )

    rm(marker_result)
    gc()
    next
  }

  write_csv_utf8(
    marker_result,
    cluster_file
  )

  marker_status <- data.frame(
    cluster = cluster_id,
    status = "COMPLETE",
    marker_n = nrow(marker_result),
    primary_marker_n = sum(
      marker_result$selection_tier == "PRIMARY"
    ),
    relaxed_marker_n = sum(
      marker_result$selection_tier == "RELAXED"
    ),
    target_pool_cell_n = length(target_pool),
    target_selected_cell_n = length(target_cells),
    background_selected_cell_n =
      length(background_cells),
    target_sample_n = length(
      unique(sample_by_cell[target_cells])
    ),
    elapsed_minutes = elapsed_minutes,
    error_message = "",
    output_file =
      normalize_slashes(cluster_file),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    marker_status,
    status_file
  )

  marker_status_rows[[length(marker_status_rows) + 1L]] <-
    marker_status

  log_event(
    "DONE descriptive marker cluster",
    cluster_id,
    "marker_n=",
    nrow(marker_result),
    "primary_n=",
    sum(marker_result$selection_tier == "PRIMARY"),
    "elapsed_min=",
    elapsed_minutes
  )

  rm(marker_result)
  gc()
}

marker_status_table <- do.call(
  rbind,
  marker_status_rows
)

marker_status_table <- marker_status_table[
  numeric_cluster_order(
    marker_status_table$cluster
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  marker_status_table,
  file.path(
    table_root,
    "02_CLUSTER_MARKER_RUN_STATUS.csv"
  )
)

failed_clusters <- marker_status_table$cluster[
  marker_status_table$status != "COMPLETE"
]


# ----------------------------------------------------------------------
# 6. Combine completed marker tables
# ----------------------------------------------------------------------

completed_cluster_files <- marker_status_table$output_file[
  marker_status_table$status == "COMPLETE"
]

completed_cluster_files <- completed_cluster_files[
  file.exists(completed_cluster_files)
]

marker_tables <- lapply(
  completed_cluster_files,
  safe_read_csv
)

marker_tables <- marker_tables[
  vapply(
    marker_tables,
    function(x) !is.null(x),
    logical(1)
  )
]

all_markers <- if (length(marker_tables)) {
  data.table::rbindlist(
    marker_tables,
    use.names = TRUE,
    fill = TRUE
  )
} else {
  data.table::data.table()
}

if (nrow(all_markers)) {
  all_markers$cluster <- as.character(
    all_markers$cluster
  )

  data.table::setorderv(
    all_markers,
    cols = c(
      "cluster",
      "marker_score",
      "mean_log_expression_difference",
      "detection_fraction_difference"
    ),
    order = c(1L, -1L, -1L, -1L),
    na.last = TRUE
  )

  write_csv_utf8(
    as.data.frame(all_markers),
    final_marker_path
  )

  top25 <- all_markers[
    ,
    head(.SD, 25L),
    by = cluster
  ]

  write_csv_utf8(
    as.data.frame(top25),
    file.path(
      table_root,
      paste0(
        dataset_id,
        "_top25_descriptive_markers_by_cluster.csv"
      )
    )
  )

  top10_summary <- top25[
    ,
    .(
      top_genes = paste(
        head(gene, 10L),
        collapse = "; "
      ),
      median_marker_score = stats::median(
        marker_score,
        na.rm = TRUE
      ),
      minimum_target_sample_breadth =
        min(
          target_sample_breadth_fraction,
          na.rm = TRUE
        )
    ),
    by = cluster
  ]

  write_csv_utf8(
    as.data.frame(top10_summary),
    file.path(
      table_root,
      paste0(
        dataset_id,
        "_top10_marker_summary.csv"
      )
    )
  )
} else {
  top25 <- data.table::data.table()
  top10_summary <- data.table::data.table()
}


# ----------------------------------------------------------------------
# 7. Marker quality audit
# ----------------------------------------------------------------------

if (nrow(top25)) {
  quality_table <- top25[
    ,
    {
      genes <- as.character(gene)
      top20_genes <- head(genes, 20L)

      mitochondrial <- grepl(
        "^MT-",
        top20_genes,
        ignore.case = TRUE
      )

      ribosomal <- grepl(
        "^RP[SL][0-9A-Z]",
        top20_genes,
        ignore.case = TRUE
      )

      stress <- top20_genes %in% c(
        "FOS",
        "FOSB",
        "JUN",
        "JUNB",
        "JUND",
        "ATF3",
        "DDIT3",
        "HSPA1A",
        "HSPA1B",
        "HSP90AA1",
        "NFKBIA",
        "IER2",
        "IER3"
      )

      cycling <- top20_genes %in% c(
        "MKI67",
        "TOP2A",
        "STMN1",
        "TUBB",
        "TYMS",
        "UBE2C",
        "CENPF",
        "CDK1",
        "PCNA"
      )

      list(
        top20_gene_n = length(top20_genes),
        mitochondrial_fraction =
          mean(mitochondrial),
        ribosomal_fraction =
          mean(ribosomal),
        stress_fraction = mean(stress),
        cycling_fraction = mean(cycling),
        minimum_sample_breadth =
          min(
            target_sample_breadth_fraction,
            na.rm = TRUE
          ),
        top20_genes = paste(
          top20_genes,
          collapse = "; "
        )
      )
    },
    by = cluster
  ]

  quality_table$review_cluster <- (
    quality_table$cluster %in%
      review_clusters
  )

  quality_table$technical_marker_flag <- (
    quality_table$mitochondrial_fraction >= 0.30 |
      quality_table$ribosomal_fraction >= 0.50 |
      quality_table$stress_fraction >= 0.30
  )

  quality_table$cell_cycle_dominant_flag <- (
    quality_table$cycling_fraction >= 0.30
  )

  write_csv_utf8(
    as.data.frame(quality_table),
    file.path(
      table_root,
      "03_CLUSTER_MARKER_QUALITY_AUDIT.csv"
    )
  )
} else {
  quality_table <- data.frame()
}


# ----------------------------------------------------------------------
# 8. Canonical broad-lineage marker panel
# ----------------------------------------------------------------------

marker_panels <- list(
  epithelial_general = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7",
    "MUC1",
    "MSLN",
    "WFDC2"
  ),
  ciliated_epithelial = c(
    "FOXJ1",
    "PIFO",
    "TPPP3",
    "CAPS",
    "C20orf85",
    "CDC20B"
  ),
  secretory_epithelial = c(
    "PAEP",
    "GPX3",
    "SCGB2A1",
    "DEFB1",
    "SLPI"
  ),
  stromal_fibroblast = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "COL6A1",
    "PDGFRA",
    "CFD"
  ),
  decidual_stromal = c(
    "PRL",
    "IGFBP1",
    "FOXO1",
    "HAND2",
    "IL15",
    "WNT4",
    "LEFTY2",
    "PROK1"
  ),
  endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "EMCN",
    "ENG",
    "RAMP2",
    "PLVAP"
  ),
  lymphatic_endothelial = c(
    "PROX1",
    "PDPN",
    "LYVE1",
    "CCL21",
    "FLT4"
  ),
  perivascular_smooth_muscle = c(
    "RGS5",
    "MCAM",
    "CSPG4",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "DES",
    "NOTCH3"
  ),
  immune_general = c(
    "PTPRC",
    "TYROBP",
    "LST1",
    "FCER1G"
  ),
  NK_cell = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "PRF1",
    "GZMB",
    "CCL5",
    "XCL1",
    "XCL2"
  ),
  T_cell = c(
    "CD3D",
    "CD3E",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "IL7R",
    "LTB"
  ),
  B_cell = c(
    "MS4A1",
    "CD79A",
    "CD74",
    "CD37",
    "CD22",
    "CD19"
  ),
  plasma_cell = c(
    "MZB1",
    "JCHAIN",
    "SDC1",
    "IGHG1",
    "IGKC",
    "DERL3"
  ),
  myeloid_macrophage = c(
    "LST1",
    "FCER1G",
    "TYROBP",
    "CTSS",
    "LILRB1",
    "CD68",
    "FCGR3A",
    "C1QA",
    "C1QB",
    "C1QC",
    "APOE",
    "CD163",
    "MRC1"
  ),
  dendritic_cell = c(
    "FCER1A",
    "CD1C",
    "CLEC10A",
    "CLEC9A",
    "IRF8",
    "LILRA4"
  ),
  mast_cell = c(
    "TPSAB1",
    "TPSB2",
    "KIT",
    "HDC",
    "CPA3",
    "MS4A2"
  ),
  erythroid = c(
    "HBB",
    "HBA1",
    "HBA2",
    "ALAS2",
    "GYPA"
  ),
  cycling = c(
    "MKI67",
    "TOP2A",
    "STMN1",
    "TYMS",
    "UBE2C",
    "CENPF"
  )
)

marker_panel_manifest <- do.call(
  rbind,
  lapply(
    names(marker_panels),
    function(panel_name) {
      data.frame(
        panel = panel_name,
        gene = marker_panels[[panel_name]],
        panel_role = ifelse(
          panel_name %in% c(
            "immune_general",
            "cycling"
          ),
          "QC_OR_GENERAL_SUPPORT",
          "CANDIDATE_LINEAGE_SUPPORT"
        ),
        evidence_basis =
          "CANONICAL_HUMAN_ENDOMETRIAL_BROAD_LINEAGE_MARKERS",
        stringsAsFactors = FALSE
      )
    }
  )
)

marker_panel_manifest$present_in_object <- (
  marker_panel_manifest$gene %in%
    rownames(object)
)

write_csv_utf8(
  marker_panel_manifest,
  file.path(
    table_root,
    "04_CANONICAL_MARKER_PANEL_MANIFEST.csv"
  )
)

panel_genes_present <- unique(
  marker_panel_manifest$gene[
    marker_panel_manifest$present_in_object
  ]
)

panel_coverage <- aggregate(
  present_in_object ~ panel,
  data = marker_panel_manifest,
  FUN = function(x) {
    c(
      present_n = sum(x),
      total_n = length(x),
      fraction = mean(x)
    )
  }
)

panel_coverage_expanded <- data.frame(
  panel = panel_coverage$panel,
  present_n =
    panel_coverage$present_in_object[, "present_n"],
  total_n =
    panel_coverage$present_in_object[, "total_n"],
  coverage_fraction =
    panel_coverage$present_in_object[, "fraction"],
  stringsAsFactors = FALSE
)

write_csv_utf8(
  panel_coverage_expanded,
  file.path(
    table_root,
    "05_CANONICAL_MARKER_PANEL_COVERAGE.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Sampled canonical marker-panel expression
# ----------------------------------------------------------------------

panel_cells <- unlist(
  lapply(
    cluster_levels,
    function(cluster_id) {
      balanced_sample_cells(
        cells = cells_by_cluster[[cluster_id]],
        sample_by_cell = sample_by_cell,
        max_total = 500L
      )
    }
  ),
  use.names = FALSE
)

panel_matrix <- SeuratObject::LayerData(
  object = object,
  assay = "RNA",
  layer = "data",
  features = panel_genes_present,
  cells = panel_cells
)

if (!inherits(panel_matrix, "sparseMatrix")) {
  panel_matrix <- Matrix::Matrix(
    panel_matrix,
    sparse = TRUE
  )
}

panel_cell_cluster <- cell_metadata[
  colnames(panel_matrix),
  "cluster"
]

panel_rows <- list()

for (cluster_id in cluster_levels) {
  idx <- which(
    panel_cell_cluster == cluster_id
  )

  one_matrix <- panel_matrix[
    ,
    idx,
    drop = FALSE
  ]

  panel_rows[[length(panel_rows) + 1L]] <-
    data.frame(
      cluster = cluster_id,
      gene = rownames(panel_matrix),
      average_log_normalized_expression =
        as.numeric(
          Matrix::rowMeans(one_matrix)
        ),
      detection_fraction = as.numeric(
        Matrix::rowMeans(one_matrix > 0)
      ),
      sampled_cell_n = length(idx),
      stringsAsFactors = FALSE
    )
}

panel_long <- do.call(
  rbind,
  panel_rows
)

panel_long <- merge(
  panel_long,
  marker_panel_manifest[
    marker_panel_manifest$present_in_object,
    c(
      "panel",
      "gene",
      "panel_role"
    )
  ],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

panel_long$cluster_numeric <- suppressWarnings(
  as.integer(panel_long$cluster)
)

panel_long <- panel_long[
  order(
    panel_long$cluster_numeric,
    panel_long$panel,
    panel_long$gene
  ),
  ,
  drop = FALSE
]

scaled_rows <- list()

for (gene_id in unique(panel_long$gene)) {
  one <- panel_long[
    panel_long$gene == gene_id,
    ,
    drop = FALSE
  ]

  values <- one$
    average_log_normalized_expression

  value_sd <- stats::sd(values)

  one$scaled_average_expression <- if (
    !is.finite(value_sd) ||
      value_sd == 0
  ) {
    rep(0, length(values))
  } else {
    as.numeric(scale(values))
  }

  scaled_rows[[length(scaled_rows) + 1L]] <-
    one
}

panel_long <- do.call(
  rbind,
  scaled_rows
)

write_csv_utf8(
  panel_long,
  file.path(
    source_data_root,
    "SOURCE_DATA_canonical_marker_panel_by_cluster.csv"
  )
)

candidate_panels <- setdiff(
  unique(marker_panel_manifest$panel),
  c(
    "immune_general",
    "cycling"
  )
)

panel_support <- aggregate(
  cbind(
    scaled_average_expression,
    detection_fraction
  ) ~ cluster + panel,
  data = panel_long[
    panel_long$panel %in% candidate_panels,
    ,
    drop = FALSE
  ],
  FUN = mean
)

colnames(panel_support)[
  colnames(panel_support) ==
    "scaled_average_expression"
] <- "mean_scaled_expression"

colnames(panel_support)[
  colnames(panel_support) ==
    "detection_fraction"
] <- "mean_detection_fraction"

panel_gene_counts <- aggregate(
  gene ~ cluster + panel,
  data = panel_long[
    panel_long$panel %in% candidate_panels,
    ,
    drop = FALSE
  ],
  FUN = length
)

colnames(panel_gene_counts)[3] <-
  "gene_n_present"

panel_support <- merge(
  panel_support,
  panel_gene_counts,
  by = c("cluster", "panel"),
  all.x = TRUE,
  sort = FALSE
)

panel_support <- panel_support[
  order(
    suppressWarnings(
      as.integer(panel_support$cluster)
    ),
    -panel_support$mean_scaled_expression
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  panel_support,
  file.path(
    table_root,
    "06_BROAD_LINEAGE_PANEL_SUPPORT.csv"
  )
)

top_panel_rows <- lapply(
  split(
    panel_support,
    panel_support$cluster
  ),
  function(x) {
    x <- x[
      order(
        -x$mean_scaled_expression,
        -x$mean_detection_fraction
      ),
      ,
      drop = FALSE
    ]

    first <- x[1, , drop = FALSE]

    second_score <- if (nrow(x) >= 2L) {
      x$mean_scaled_expression[2]
    } else {
      NA_real_
    }

    data.frame(
      cluster = first$cluster,
      top_support_panel = first$panel,
      top_support_score =
        first$mean_scaled_expression,
      second_support_score = second_score,
      support_margin = (
        first$mean_scaled_expression -
          second_score
      ),
      top_panel_gene_n_present =
        first$gene_n_present,
      top_panel_detection_fraction =
        first$mean_detection_fraction,
      annotation_status =
        "ALGORITHMIC_SUPPORT_ONLY_REQUIRES_MANUAL_REVIEW",
      stringsAsFactors = FALSE
    )
  }
)

top_panel_summary <- do.call(
  rbind,
  top_panel_rows
)

top_panel_summary$review_cluster <- (
  top_panel_summary$cluster %in%
    review_clusters
)

top_panel_summary <- top_panel_summary[
  numeric_cluster_order(
    top_panel_summary$cluster
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  top_panel_summary,
  file.path(
    table_root,
    "07_TOP_BROAD_LINEAGE_SUPPORT_BY_CLUSTER.csv"
  )
)


# ----------------------------------------------------------------------
# 10. QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(
  base_size = 9
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    strip.background =
      ggplot2::element_rect(fill = "grey95")
  )

panel_long$cluster <- factor(
  panel_long$cluster,
  levels = rev(cluster_levels)
)

panel_long$gene <- factor(
  panel_long$gene,
  levels = unique(
    marker_panel_manifest$gene[
      marker_panel_manifest$present_in_object
    ]
  )
)

plot_dot <- ggplot2::ggplot(
  panel_long,
  ggplot2::aes(
    x = gene,
    y = cluster,
    size = detection_fraction,
    color = scaled_average_expression
  )
) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(
    . ~ panel,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_size_continuous(
    range = c(0.2, 5),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::scale_color_viridis_c(
    option = "C"
  ) +
  ggplot2::labs(
    title =
      "GSE250130 broad-lineage marker support",
    subtitle = paste(
      "Balanced sampled expression;",
      "not a locked cell-type annotation"
    ),
    x = "Canonical marker",
    y = "Cluster",
    size = "Detected",
    color = "Scaled\naverage"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    ),
    strip.text.x = ggplot2::element_text(
      angle = 90,
      size = 7
    )
  )

save_plot_pair(
  plot_dot,
  "QC_Figure_1_broad_lineage_marker_dotplot",
  width = 20,
  height = 8
)

if (nrow(top25)) {
  top5 <- top25[
    ,
    head(.SD, 5L),
    by = cluster
  ]

  top5$cluster <- factor(
    top5$cluster,
    levels = cluster_levels
  )

  top5$gene_cluster <- paste0(
    top5$gene,
    " [C",
    top5$cluster,
    "]"
  )

  top5$gene_cluster <- factor(
    top5$gene_cluster,
    levels = rev(top5$gene_cluster)
  )

  plot_top_markers <- ggplot2::ggplot(
    as.data.frame(top5),
    ggplot2::aes(
      x = cluster,
      y = gene_cluster,
      fill = marker_score
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      option = "B"
    ) +
    ggplot2::labs(
      title =
        "Top five descriptive markers per cluster",
      subtitle = paste(
        "Balanced sparse marker scoring;",
        "annotation support only"
      ),
      x = "Cluster",
      y = "Marker",
      fill = "Marker\nscore"
    ) +
    theme_qc +
    ggplot2::theme(
      axis.text.y =
        ggplot2::element_text(size = 5)
    )

  save_plot_pair(
    plot_top_markers,
    "QC_Figure_2_top5_markers_per_cluster",
    width = 10,
    height = 18
  )
}

plot_support <- ggplot2::ggplot(
  panel_support,
  ggplot2::aes(
    x = panel,
    y = factor(
      cluster,
      levels = rev(cluster_levels)
    ),
    fill = mean_scaled_expression
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    option = "D"
  ) +
  ggplot2::labs(
    title = "Broad-lineage support matrix",
    subtitle =
      "Panel support requires marker review",
    x = "Marker panel",
    y = "Cluster",
    fill = "Mean scaled\nexpression"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_pair(
  plot_support,
  "QC_Figure_3_broad_lineage_support_matrix",
  width = 11,
  height = 7
)

plot_marker_n <- ggplot2::ggplot(
  marker_status_table,
  ggplot2::aes(
    x = factor(
      cluster,
      levels = cluster_levels
    ),
    y = marker_n,
    fill = status
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title =
      "Descriptive marker yield by cluster",
    x = "Cluster",
    y = "Returned marker count",
    fill = "Status"
  ) +
  theme_qc

save_plot_pair(
  plot_marker_n,
  "QC_Figure_4_marker_yield_by_cluster",
  width = 8,
  height = 5
)


# ----------------------------------------------------------------------
# 11. Gate, Card, and Methods draft
# ----------------------------------------------------------------------

all_cluster_markers_complete <- (
  nrow(marker_status_table) ==
    length(cluster_levels) &&
    all(
      marker_status_table$status ==
        "COMPLETE"
    )
)

panel_coverage_pass <- (
  mean(
    marker_panel_manifest$present_in_object
  ) >= 0.80
)

marker_yield_pass <- (
  all(marker_status_table$marker_n >= 10L)
)

gate_status <- if (
  all_cluster_markers_complete &&
    panel_coverage_pass &&
    marker_yield_pass
) {
  "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED"
} else {
  "HOLD_DESCRIPTIVE_MARKER_OR_PANEL_QC_INCOMPLETE"
}

marker_gate <- data.frame(
  dataset_id = dataset_id,
  input_version = input_version,
  marker_version = marker_version,
  supersedes_marker_version =
    "MARKERS_v1_0_FAILED_FINDMARKERS",
  gate_status = gate_status,
  marker_method =
    "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING",
  p_values_reported = FALSE,
  cluster_n_expected =
    length(cluster_levels),
  cluster_marker_complete_n = sum(
    marker_status_table$status == "COMPLETE"
  ),
  cluster_marker_failed_n = sum(
    marker_status_table$status != "COMPLETE"
  ),
  failed_clusters =
    collapse_unique(failed_clusters),
  canonical_marker_present_fraction = mean(
    marker_panel_manifest$present_in_object
  ),
  minimum_marker_n = if (
    nrow(marker_status_table)
  ) {
    min(marker_status_table$marker_n)
  } else {
    NA_integer_
  },
  review_clusters =
    collapse_unique(review_clusters),
  annotation_review_authorized = (
    gate_status ==
      "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED"
  ),
  biological_annotation_locked = FALSE,
  VDR_localization_authorized = FALSE,
  pseudobulk_authorized = FALSE,
  result_lock_authorized = FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "STEP04E1D_GSE250130_MARKER_GATE.csv"
)

write_csv_utf8(
  marker_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1D Marker Result Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Input: ", input_version),
  paste0("- Marker version: ", marker_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Completed clusters: ",
    sum(
      marker_status_table$status ==
        "COMPLETE"
    ),
    "/",
    length(cluster_levels)
  ),
  paste0(
    "- Review clusters: ",
    ifelse(
      length(review_clusters),
      collapse_unique(review_clusters),
      "none"
    )
  ),
  paste0(
    "- Marker method: ",
    marker_gate$marker_method
  ),
  "- P-values reported: FALSE",
  "- Marker role: cluster annotation support only",
  "- Donor-level inferential role: none",
  paste0(
    "- Annotation review authorized: ",
    marker_gate$annotation_review_authorized
  ),
  "- Biological annotation locked: FALSE",
  "- VDR localization authorized: FALSE",
  "- Pseudobulk authorized: FALSE",
  "- Result lock authorized: FALSE"
)

card_path <- file.path(
  run_qc_root,
  "STEP04E1D_GSE250130_MARKER_RESULT_CARD.md"
)

write_text(
  card_lines,
  card_path
)

interpretation_lines <- c(
  "# STEP04E1D GSE250130 descriptive cluster markers",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  paste0(
    "- Completed clusters: ",
    sum(
      marker_status_table$status ==
        "COMPLETE"
    ),
    "/",
    length(cluster_levels)
  ),
  paste0(
    "- Canonical marker-panel coverage: ",
    round(
      100 *
        mean(
          marker_panel_manifest$
            present_in_object
        ),
      1
    ),
    "%"
  ),
  "",
  "## Scientific interpretation rule",
  "",
  "- The failed MARKERS_v1_0 run is retained as a provenance record.",
  "- MARKERS_v1_1 uses explicitly balanced, sampled sparse matrices and does not call FindMarkers.",
  "- Marker scores are descriptive and used only for cell-identity annotation.",
  "- No cell-level p-value is reported or interpreted.",
  "- Clinical and temporal inference remains reserved for donor-aware pseudobulk.",
  "- Clusters 12, 14, 20, and 21 require explicit review because of prior sample-mixing flags.",
  "",
  "## Next step",
  "",
  "- Review top descriptive markers and canonical marker-panel support.",
  "- Construct an explicit cluster-to-cell-type mapping in STEP04E1E.",
  "- Do not calculate final VDR localization or donor-aware pseudobulk until the annotation Gate passes."
)

interpretation_path <- file.path(
  run_qc_root,
  "STEP04E1D_GSE250130_MARKER_INTERPRETATION.md"
)

write_text(
  interpretation_lines,
  interpretation_path
)

methods_lines <- c(
  "# STEP04E1D Methods draft",
  "",
  paste(
    "Cluster markers were derived from the joined log-normalized RNA data",
    "for the resolution-0.4 clustering solution using an explicitly",
    "memory-bounded descriptive scoring workflow."
  ),
  paste(
    "For each cluster, up to 1,500 target cells were sampled while balancing",
    "representation across contributing samples."
  ),
  paste(
    "The reference set was assembled by sampling up to 150 cells from each",
    "remaining cluster, also balancing cells across samples."
  ),
  paste(
    "Only the selected cells were retrieved from the sparse normalized-data",
    "layer. For every gene, mean log-normalized expression, detection",
    "fractions, detection-fraction difference, target-sample breadth, and a",
    "composite marker score were calculated."
  ),
  paste(
    "Primary marker candidates required a target detection fraction of at",
    "least 10%, a detection-fraction difference of at least 5 percentage",
    "points, a mean log-expression difference of at least 0.20, and detection",
    "in at least 25% of target-contributing samples."
  ),
  paste(
    "A documented relaxed threshold was used only when fewer than 25 primary",
    "markers were available for a cluster."
  ),
  paste(
    "These descriptive marker statistics were used solely for cell-identity",
    "annotation. No cell-level p-values were used for clinical or temporal",
    "inference, which was reserved for subsequent donor-aware pseudobulk",
    "analysis."
  )
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1D_GSE250130_marker_Methods_draft_v1_1.md"
  )
)


# ----------------------------------------------------------------------
# 12. Session, manifests, cleanup, and ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_names <- c(
  "SeuratObject",
  "Matrix",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales",
  "future"
)

package_table <- data.frame(
  package = package_names,
  installed = vapply(
    package_names,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  version = vapply(
    package_names,
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(
    run_qc_root,
    "PACKAGE_VERSIONS.csv"
  )
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- unique(c(
  list.files(
    marker_result_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1D_GSE250130_marker_Methods_draft_v1_1.md"
  )
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(artifact_files),
  size_bytes =
    file.info(artifact_files)$size,
  md5 =
    unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(
    run_qc_root,
    "ARTIFACT_MANIFEST.csv"
  )
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  marker_version = marker_version,
  supersedes_marker_version =
    "MARKERS_v1_0_FINDMARKERS_MEMORY_FAILURE",
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_minutes = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "mins"
      )
    ),
    2
  ),
  gate_status = gate_status,
  cluster_n = length(cluster_levels),
  completed_cluster_n = sum(
    marker_status_table$status == "COMPLETE"
  ),
  failed_cluster_n = sum(
    marker_status_table$status != "COMPLETE"
  ),
  marker_method =
    "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING",
  p_values_reported = FALSE,
  inferential_role =
    "NONE_ANNOTATION_SUPPORT_ONLY",
  checkpoint_modified = FALSE,
  biological_annotation_locked = FALSE,
  result_lock_authorized = FALSE,
  marker_result_root =
    normalize_slashes(marker_result_root),
  qc_root =
    normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "completed_clusters=",
  sum(
    marker_status_table$status ==
      "COMPLETE"
  ),
  "failed_clusters=",
  collapse_unique(failed_clusters)
)

rm(
  panel_matrix,
  object
)
gc()

small_files <- unique(c(
  final_marker_path,
  file.path(
    table_root,
    paste0(
      dataset_id,
      "_top25_descriptive_markers_by_cluster.csv"
    )
  ),
  file.path(
    table_root,
    paste0(
      dataset_id,
      "_top10_marker_summary.csv"
    )
  ),
  list.files(
    table_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    source_data_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <=
      150 * 1024^2
]

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = small_files,
        recurse = FALSE,
        include_directories = FALSE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat(
  "STEP04E1D GSE250130 memory-safe marker analysis completed\n"
)
cat("============================================================\n")
cat(
  "Gate status: ",
  gate_status,
  "\n",
  sep = ""
)
cat(
  "Completed clusters: ",
  sum(
    marker_status_table$status ==
      "COMPLETE"
  ),
  "/",
  length(cluster_levels),
  "\n",
  sep = ""
)
cat(
  "Failed clusters: ",
  ifelse(
    length(failed_clusters),
    collapse_unique(failed_clusters),
    "none"
  ),
  "\n",
  sep = ""
)
cat(
  "Marker result root: ",
  marker_result_root,
  "\n",
  sep = ""
)
cat(
  "QC root: ",
  run_qc_root,
  "\n",
  sep = ""
)
cat(
  "Run manifest: ",
  run_manifest_path,
  "\n",
  sep = ""
)

if (zip_ok) {
  cat(
    "Upload this ZIP:\n",
    handoff_zip,
    "\n"
  )
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe result and QC folders are complete.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 5: 59_step04E1E_GSE250130_broad_annotation_QC_v1_2_pathsafe
# ==============================================================================

# ======================================================================
# 59_step04E1E_GSE250130_broad_annotation_QC_v1_2_pathsafe.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Apply an evidence-reviewed cluster-to-cell-type mapping for GSE250130.
# - Lock broad lineage labels only at the module-QC level.
# - Retain detailed subtype labels as candidates requiring later review.
# - Generate cell-level annotation metadata without loading the 16.59 GB
#   Seurat object.
# - Produce annotation QC figures, source data, Gate, Card, Methods and
#   Results drafts.
#
# BROAD CELL-TYPE NOMENCLATURE
# Aligned to the eight major populations reported for GSE250130:
#   1. Unciliated epithelial cells
#   2. Ciliated epithelial cells
#   3. Stromal cells
#   4. Endothelial cells
#   5. NK/T cells
#   6. Myeloid cells
#   7. B cells
#   8. Mast cells
#
# SCIENTIFIC GUARDRAILS
# - Broad annotations are supported by cluster markers and canonical panels.
# - Fine subtype labels remain candidate labels.
# - Clusters 12, 14, 20 and 21 retain explicit review flags.
# - Sample-level cell composition is descriptive QC only.
# - No composition hypothesis test is performed.
# - No output is written to 06_locked_results.
#
# INPUTS
# - STEP04E1C integration QC Gate
# - STEP04E1D descriptive marker Gate
# - GSE250130 UMAP/cluster metadata
# - MARKERS_v1_1_memory_safe marker tables
# - COMBINED_SAMPLE_DESIGN_v1.0 GSE250130 design
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)

set.seed(20260805)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <- "59_step04E1E_GSE250130_broad_annotation_QC_v1_2_pathsafe.R"
script_version <- "STEP04E1E_GSE250130_BROAD_ANNOTATION_QC_v1.2_PATHSAFE"
input_version <- "SCI_v1_1_memory_safe"
marker_version <- "MARKERS_v1_1_memory_safe"
annotation_version <- "ANN_v1_2"

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1E_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(
  path,
  max_chars = 245L
) {
  normalized <- normalize_slashes(path)
  path_n <- nchar(normalized, type = "chars")

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Output path is too long for safe Windows/R writing (",
      path_n,
      " characters; safe limit ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

split_semicolon <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) return(character())
  trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

detect_first_column <- function(column_names, candidates) {
  lower <- tolower(column_names)

  for (candidate in tolower(candidates)) {
    hit <- which(lower == candidate)
    if (length(hit)) return(column_names[hit[1]])
  }

  ""
}

save_plot_pair <- function(plot_object, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))

  assert_windows_path_length(png_path)
  assert_windows_path_length(pdf_path)

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")

  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages and project directories
# ----------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(step04_root, "05_execution")
results_root <- file.path(step04_root, "06_results")
qc_base_root <- file.path(step04_root, "08_qc")
parameter_root <- file.path(step04_root, "03_parameters")
environment_root <- file.path(step04_root, "04_environment")
manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

input_table_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "tables",
  dataset_id
)

marker_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "markers",
  dataset_id,
  marker_version
)

marker_table_root <- file.path(marker_root, "tables")

annotation_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "annotations",
  dataset_id,
  annotation_version
)

annotation_table_root <- file.path(annotation_root, "tables")
annotation_source_root <- file.path(annotation_root, "source_data")

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1E_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_annQC"
  )
)

figure_root <- file.path(run_qc_root, "figures")

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  annotation_root,
  annotation_table_root,
  annotation_source_root,
  run_qc_root,
  figure_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)

planned_path_audit <- data.frame(
  object = c(
    "annotation_root",
    "annotation_table_root",
    "annotation_source_root",
    "run_qc_root",
    "figure_root",
    "log_path",
    "run_manifest_path",
    "parameter_path",
    "session_path",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    annotation_root,
    annotation_table_root,
    annotation_source_root,
    run_qc_root,
    figure_root,
    log_path,
    run_manifest_path,
    parameter_path,
    session_path,
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_path_audit$character_n <- nchar(
  planned_path_audit$path,
  type = "chars"
)

if (
  .Platform$OS.type == "windows" &&
    any(planned_path_audit$character_n > 245L)
) {
  bad_paths <- planned_path_audit[
    planned_path_audit$character_n > 245L,
    ,
    drop = FALSE
  ]

  stop(
    "One or more planned output paths remain too long:\n",
    paste(
      paste0(
        bad_paths$character_n,
        " chars: ",
        bad_paths$path
      ),
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------------------
# 3. Locate and validate required inputs
# ----------------------------------------------------------------------

umap_path <- file.path(
  input_table_root,
  paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv")
)

marker_gate_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

marker_gate_dirs <- marker_gate_dirs[
  grepl(
    "^STEP04E1D_GSE250130_v1_1_[0-9]{8}_[0-9]{6}_cluster_markers_memory_safe$",
    basename(marker_gate_dirs),
    perl = TRUE
  )
]

if (!length(marker_gate_dirs)) {
  stop("No completed STEP04E1D memory-safe marker QC folder found.")
}

latest_marker_qc <- marker_gate_dirs[
  which.max(file.info(marker_gate_dirs)$mtime)
]

marker_gate_path <- file.path(
  latest_marker_qc,
  "STEP04E1D_GSE250130_MARKER_GATE.csv"
)

integration_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

integration_qc_dirs <- integration_qc_dirs[
  grepl(
    "^STEP04E1C_GSE250130_[0-9]{8}_[0-9]{6}_integration_QC$",
    basename(integration_qc_dirs),
    perl = TRUE
  )
]

if (!length(integration_qc_dirs)) {
  stop("No completed STEP04E1C integration QC folder found.")
}

latest_integration_qc <- integration_qc_dirs[
  which.max(file.info(integration_qc_dirs)$mtime)
]

integration_gate_path <- file.path(
  latest_integration_qc,
  "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv"
)

all_marker_path <- file.path(
  marker_table_root,
  paste0(dataset_id, "_all_cluster_markers_descriptive.csv")
)

top25_path <- file.path(
  marker_table_root,
  paste0(dataset_id, "_top25_descriptive_markers_by_cluster.csv")
)

top10_path <- file.path(
  marker_table_root,
  paste0(dataset_id, "_top10_marker_summary.csv")
)

panel_support_path <- file.path(
  marker_table_root,
  "06_BROAD_LINEAGE_PANEL_SUPPORT.csv"
)

top_panel_path <- file.path(
  marker_table_root,
  "07_TOP_BROAD_LINEAGE_SUPPORT_BY_CLUSTER.csv"
)

quality_path <- file.path(
  marker_table_root,
  "03_CLUSTER_MARKER_QUALITY_AUDIT.csv"
)

locked_design_path <- file.path(
  results_root,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0",
  "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"
)

required_inputs <- c(
  umap_path,
  marker_gate_path,
  integration_gate_path,
  all_marker_path,
  top25_path,
  top10_path,
  panel_support_path,
  top_panel_path,
  quality_path,
  locked_design_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

marker_gate <- safe_read_csv(marker_gate_path)
integration_gate <- safe_read_csv(integration_gate_path)

if (
  is.null(marker_gate) ||
    nrow(marker_gate) != 1L ||
    marker_gate$gate_status[[1]] !=
      "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED" ||
    !isTRUE(marker_gate$annotation_review_authorized[[1]])
) {
  stop("STEP04E1D marker Gate does not authorize annotation review.")
}

if (
  is.null(integration_gate) ||
    nrow(integration_gate) != 1L ||
    !integration_gate$gate_status[[1]] %in%
      c("PASS", "PASS_WITH_REVIEW_FLAGS")
) {
  stop("STEP04E1C integration Gate is not valid.")
}

umap <- data.table::fread(
  umap_path,
  data.table = FALSE,
  showProgress = FALSE
)

all_markers <- data.table::fread(
  all_marker_path,
  data.table = FALSE,
  showProgress = FALSE
)

top25 <- data.table::fread(
  top25_path,
  data.table = FALSE,
  showProgress = FALSE
)

top10 <- safe_read_csv(top10_path)
panel_support <- safe_read_csv(panel_support_path)
top_panel <- safe_read_csv(top_panel_path)
marker_quality <- safe_read_csv(quality_path)
locked_design <- safe_read_csv(locked_design_path)

required_umap_columns <- c(
  "cell_barcode",
  "UMAPH_1",
  "UMAPH_2",
  "sample_accession",
  "seurat_clusters"
)

missing_umap_columns <- setdiff(
  required_umap_columns,
  colnames(umap)
)

if (length(missing_umap_columns)) {
  stop(
    "UMAP table is missing columns: ",
    paste(missing_umap_columns, collapse = ", ")
  )
}

if (
  nrow(umap) != 236636L ||
    length(unique(umap$sample_accession)) != 28L ||
    length(unique(umap$seurat_clusters)) != 22L
) {
  stop("UMAP identity checks failed.")
}

if (
  is.null(locked_design) ||
    !"sample_accession" %in% colnames(locked_design) ||
    nrow(locked_design) != 28L
) {
  stop("Locked GSE250130 sample design is invalid.")
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  marker_version = marker_version,
  annotation_version = annotation_version,
  supersedes_partial_annotation_version =
    "ANNOTATION_v1_1_FIXED_PATH_LENGTH_FAILURE",
  dataset_id = dataset_id,
  broad_nomenclature_source =
    "GSE250130 published eight-major-cell-type framework",
  original_study_doi = "10.1038/s41467-024-55419-z",
  reference_atlas_doi = "10.1038/s41588-024-01873-w",
  older_endometrial_atlas_doi = "10.1038/s41591-020-1040-z",
  broad_annotation_role =
    "ACCEPTED_FOR_BROAD_LEVEL_LOCALIZATION_AND_AGGREGATION_AFTER_GATE",
  subtype_annotation_role =
    "CANDIDATE_REQUIRES_REVIEW",
  review_clusters = c("12", "14", "20", "21"),
  composition_role =
    "DESCRIPTIVE_QC_ONLY_NO_HYPOTHESIS_TEST",
  random_seed = 20260805L,
  result_lock_authorized = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

log_event(
  "Starting",
  run_id,
  "marker_gate=",
  marker_gate$gate_status[[1]],
  "integration_gate=",
  integration_gate$gate_status[[1]]
)


write_csv_utf8(
  planned_path_audit,
  file.path(
    run_qc_root,
    "00_path_length_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 4. Evidence-reviewed cluster annotation map
# ----------------------------------------------------------------------

annotation_map <- data.frame(
  cluster = as.character(0:21),

  broad_cell_type = c(
    "Stromal cells",
    "Unciliated epithelial cells",
    "Stromal cells",
    "NK/T cells",
    "NK/T cells",
    "NK/T cells",
    "Unciliated epithelial cells",
    "NK/T cells",
    "NK/T cells",
    "Myeloid cells",
    "Ciliated epithelial cells",
    "NK/T cells",
    "B cells",
    "Myeloid cells",
    "Mast cells",
    "Endothelial cells",
    "Endothelial cells",
    "NK/T cells",
    "Myeloid cells",
    "Unciliated epithelial cells",
    "Unciliated epithelial cells",
    "Stromal cells"
  ),

  candidate_subtype = c(
    "COL1A1_IGF1_stromal",
    "MT_high_unciliated_epithelial",
    "PTGDS_APOD_stromal",
    "NCAM1_GNLY_NK",
    "conventional_T",
    "KLRB1_TRDC_uNK_like",
    "PAX8_GABRP_unciliated_epithelial",
    "cycling_NK_T",
    "FGFBP2_cytotoxic_NK",
    "HLA_DRA_LYZ_myeloid_APC",
    "ciliated_epithelial",
    "TRDC_GNLY_gamma_delta_like_NK_T",
    "MS4A1_CD79A_B",
    "CSF3R_S100A8_A9_neutrophil",
    "TPSAB1_CPA3_mast",
    "VWF_EMCN_blood_endothelial",
    "CCL21_PROX1_lymphatic_endothelial",
    "IL7R_KLRB1_T",
    "IL3RA_CLEC4C_GZMB_plasmacytoid_DC",
    "KRT19_SLPI_unciliated_epithelial",
    "ESR1_LIFR_LPAR3_luminal_like_epithelial",
    "C11orf96_INHBA_activated_perivascular_stromal"
  ),

  broad_confidence = c(
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "MODERATE",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "MODERATE",
    "MODERATE"
  ),

  subtype_confidence = c(
    "MODERATE",
    "MODERATE",
    "MODERATE",
    "HIGH",
    "HIGH",
    "MODERATE",
    "MODERATE",
    "HIGH",
    "HIGH",
    "MODERATE",
    "HIGH",
    "MODERATE",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "MODERATE",
    "MODERATE_REVIEW",
    "MODERATE_REVIEW"
  ),

  review_flag = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    TRUE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE
  ),

  review_reason = c(
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "NK/T subtype boundary requires conservative naming",
    "B-cell identity is clear; sample-mixing/abundance remains under review",
    "",
    "Mast-cell identity is clear; strong single-sample contribution remains under review",
    "",
    "",
    "",
    "",
    "",
    "Rare early-luteal/luminal-like epithelial cluster; limited to seven samples",
    "Rare activated perivascular/stromal cluster; five-sample and single-sample-dominance review"
  ),

  evidence_genes = c(
    "COL1A1;COL3A1;PDGFRA;DCN;IGF1",
    "KRT19;KRT18;KRT8;SCGB2A1;WFDC2;MT1G",
    "COL1A1;COL3A1;DCN;PTGDS;APOD",
    "GNLY;NCAM1;NKG7;KLRD1;GZMB",
    "CD3D;CD3E;TRAC;CD2;BCL11B",
    "NKG7;KLRB1;TRDC;XCL1;XCL2",
    "EPCAM;KRT8;KRT18;KRT19;PAX8;GABRP",
    "MKI67;TOP2A;STMN1;NKG7;GZMA",
    "NKG7;FGFBP2;GNLY;PRF1;GZMB",
    "LYZ;HLA-DRA;CTSS;TYROBP;AIF1",
    "DNAH12;RSPH1;CAPS;TPPP3;C20orf85",
    "TRDC;GNLY;NKG7;GZMB;CTSW",
    "MS4A1;CD79A;CD74;BANK1;CD37",
    "S100A8;S100A9;CSF3R;FCGR3B;AQP9",
    "TPSAB1;TPSB2;CPA3;KIT;HDC;MS4A2",
    "VWF;EMCN;ADGRL4;PTPRB;ERG",
    "CCL21;PROX1;PDPN;FLT4;CLDN5",
    "IL7R;KLRB1;TRBC1;LTB;ICOS",
    "IL3RA;CLEC4C;GZMB;TCF4;JCHAIN",
    "KRT19;KRT18;KRT7;SLPI;WFDC2;EHF",
    "EPCAM;KRT8;KRT18;PAX8;ERBB4;LIFR;LPAR3;ESR1",
    "COL1A1;COL3A1;PDGFRB;CALD1;PRRX1;TAGLN;NOTCH3;INHBA"
  ),

  annotation_basis = c(
    rep(
      "DESCRIPTIVE_CLUSTER_MARKERS_PLUS_CANONICAL_PANEL_PLUS_PUBLISHED_GSE250130_MAJOR_TYPE_FRAMEWORK",
      22
    )
  ),

  subtype_status = rep(
    "CANDIDATE_NOT_LOCKED",
    22
  ),

  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# 5. Annotation evidence audit
# ----------------------------------------------------------------------

all_markers$cluster <- as.character(all_markers$cluster)
top25$cluster <- as.character(top25$cluster)

evidence_audit_rows <- list()

for (i in seq_len(nrow(annotation_map))) {
  cluster_id <- annotation_map$cluster[i]
  expected <- split_semicolon(
    annotation_map$evidence_genes[i]
  )

  one_all <- all_markers[
    all_markers$cluster == cluster_id,
    ,
    drop = FALSE
  ]

  one_top25 <- top25[
    top25$cluster == cluster_id,
    ,
    drop = FALSE
  ]

  observed_all <- intersect(expected, one_all$gene)
  observed_top25 <- intersect(expected, one_top25$gene)

  evidence_audit_rows[[i]] <- data.frame(
    cluster = cluster_id,
    broad_cell_type = annotation_map$broad_cell_type[i],
    candidate_subtype = annotation_map$candidate_subtype[i],
    expected_evidence_gene_n = length(expected),
    observed_evidence_gene_n_all_markers = length(observed_all),
    observed_evidence_fraction_all_markers =
      length(observed_all) / length(expected),
    observed_evidence_genes_all_markers =
      collapse_unique(observed_all),
    observed_evidence_gene_n_top25 = length(observed_top25),
    observed_evidence_genes_top25 =
      collapse_unique(observed_top25),
    review_flag = annotation_map$review_flag[i],
    stringsAsFactors = FALSE
  )
}

evidence_audit <- do.call(
  rbind,
  evidence_audit_rows
)

annotation_map <- merge(
  annotation_map,
  evidence_audit[
    ,
    c(
      "cluster",
      "observed_evidence_fraction_all_markers",
      "observed_evidence_genes_all_markers",
      "observed_evidence_genes_top25"
    )
  ],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

annotation_map <- annotation_map[
  order(as.integer(annotation_map$cluster)),
  ,
  drop = FALSE
]

annotation_map$broad_annotation_status <- ifelse(
  annotation_map$observed_evidence_fraction_all_markers >= 0.60,
  ifelse(
    annotation_map$review_flag,
    "ACCEPTED_BROAD_WITH_REVIEW",
    "ACCEPTED_BROAD"
  ),
  "HOLD_EVIDENCE_INCOMPLETE"
)

write_csv_utf8(
  annotation_map,
  file.path(
    annotation_table_root,
    "cluster_annotation_map.csv"
  )
)

write_csv_utf8(
  evidence_audit,
  file.path(
    run_qc_root,
    "01_annotation_evidence_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 6. Create cell-level annotations
# ----------------------------------------------------------------------

umap$cluster <- as.character(umap$seurat_clusters)

cell_annotation <- merge(
  umap,
  annotation_map[
    ,
    c(
      "cluster",
      "broad_cell_type",
      "candidate_subtype",
      "broad_confidence",
      "subtype_confidence",
      "review_flag",
      "review_reason",
      "broad_annotation_status",
      "subtype_status"
    )
  ],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(cell_annotation) != 236636L ||
    anyNA(cell_annotation$broad_cell_type)
) {
  stop("Cell-level annotation merge failed.")
}

locked_design$sample_accession <- as.character(
  locked_design$sample_accession
)

cell_annotation$sample_accession <- as.character(
  cell_annotation$sample_accession
)

cell_annotation <- merge(
  cell_annotation,
  locked_design,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(cell_annotation) != 236636L) {
  stop("Locked-design merge changed the cell count.")
}

cell_annotation <- cell_annotation[
  order(match(cell_annotation$cell_barcode, umap$cell_barcode)),
  ,
  drop = FALSE
]

cell_annotation_path <- file.path(
  annotation_table_root,
  "cell_annotation_metadata.csv"
)

write_csv_utf8(
  cell_annotation,
  cell_annotation_path
)


# ----------------------------------------------------------------------
# 7. Broad-type summaries and published-profile contextual QC
# ----------------------------------------------------------------------

major_summary <- aggregate(
  cell_barcode ~ broad_cell_type,
  data = cell_annotation,
  FUN = length
)

colnames(major_summary)[2] <- "cell_n"

major_summary$cell_fraction <- (
  major_summary$cell_n /
    sum(major_summary$cell_n)
)

major_summary <- major_summary[
  order(-major_summary$cell_n),
  ,
  drop = FALSE
]

published_profile <- data.frame(
  broad_cell_type = c(
    "Unciliated epithelial cells",
    "Ciliated epithelial cells",
    "Stromal cells",
    "Endothelial cells",
    "NK/T cells",
    "Myeloid cells",
    "B cells",
    "Mast cells"
  ),
  published_cell_fraction = c(
    0.168,
    0.019,
    0.358,
    0.006,
    0.385,
    0.038,
    0.018,
    0.006
  ),
  source = "Cao_et_al_Nat_Commun_2025_DOI_10.1038_s41467-024-55419-z",
  role = "CONTEXTUAL_QC_NOT_EXACT_REPRODUCTION_TARGET",
  stringsAsFactors = FALSE
)

major_summary <- merge(
  major_summary,
  published_profile,
  by = "broad_cell_type",
  all.x = TRUE,
  sort = FALSE
)

major_summary$absolute_fraction_difference <- abs(
  major_summary$cell_fraction -
    major_summary$published_cell_fraction
)

major_summary <- major_summary[
  order(-major_summary$cell_n),
  ,
  drop = FALSE
]

write_csv_utf8(
  major_summary,
  file.path(
    annotation_table_root,
    "broad_type_summary.csv"
  )
)

write_csv_utf8(
  published_profile,
  file.path(
    annotation_source_root,
    "published_major_profile.csv"
  )
)

cluster_summary <- merge(
  aggregate(
    cell_barcode ~ cluster,
    data = cell_annotation,
    FUN = length
  ),
  annotation_map,
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

colnames(cluster_summary)[
  colnames(cluster_summary) == "cell_barcode"
] <- "cell_n"

cluster_summary <- cluster_summary[
  order(as.integer(cluster_summary$cluster)),
  ,
  drop = FALSE
]

write_csv_utf8(
  cluster_summary,
  file.path(
    annotation_table_root,
    "cluster_annotation_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 8. Sample-level descriptive composition
# ----------------------------------------------------------------------

sample_type_counts <- aggregate(
  cell_barcode ~ sample_accession + broad_cell_type,
  data = cell_annotation,
  FUN = length
)

colnames(sample_type_counts)[3] <- "cell_n"

sample_totals <- aggregate(
  cell_n ~ sample_accession,
  data = sample_type_counts,
  FUN = sum
)

colnames(sample_totals)[2] <- "sample_total_cell_n"

sample_type_counts <- merge(
  sample_type_counts,
  sample_totals,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

sample_type_counts$within_sample_fraction <- (
  sample_type_counts$cell_n /
    sample_type_counts$sample_total_cell_n
)

sample_type_counts <- merge(
  sample_type_counts,
  locked_design,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

write_csv_utf8(
  sample_type_counts,
  file.path(
    annotation_table_root,
    "sample_broad_type_composition.csv"
  )
)

group_column <- detect_first_column(
  colnames(locked_design),
  c(
    "group",
    "clinical_group",
    "phenotype",
    "condition",
    "fertility_status",
    "patient_group",
    "case_control_status",
    "case_control"
  )
)

time_column <- detect_first_column(
  colnames(locked_design),
  c(
    "timepoint",
    "time_point",
    "temporal_label",
    "cycle_day",
    "lh_day",
    "day",
    "phase",
    "sampling_time"
  )
)

role_column <- detect_first_column(
  colnames(locked_design),
  c(
    "analysis_role",
    "cohort_role",
    "sample_role",
    "role"
  )
)

design_detection <- data.frame(
  design_role = c(
    "GROUP",
    "TIME",
    "ANALYSIS_ROLE"
  ),
  detected_column = c(
    group_column,
    time_column,
    role_column
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  design_detection,
  file.path(
    run_qc_root,
    "02_design_columns.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Annotation QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

# Deterministic plotting subset: max 10,000 cells per broad type.
cell_dt <- data.table::as.data.table(cell_annotation)

plot_cells <- cell_dt[
  ,
  {
    keep_n <- min(.N, 10000L)
    .SD[sample.int(.N, keep_n)]
  },
  by = broad_cell_type
]

broad_centroids <- cell_dt[
  ,
  .(
    UMAPH_1 = stats::median(UMAPH_1),
    UMAPH_2 = stats::median(UMAPH_2)
  ),
  by = broad_cell_type
]

plot_umap_broad <- ggplot2::ggplot(
  plot_cells,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = broad_cell_type
  )
) +
  ggplot2::geom_point(
    size = 0.08,
    alpha = 0.60,
    stroke = 0
  ) +
  ggplot2::geom_text(
    data = broad_centroids,
    ggplot2::aes(
      x = UMAPH_1,
      y = UMAPH_2,
      label = broad_cell_type
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold"
  ) +
  ggplot2::labs(
    title = "GSE250130 broad cell-type annotation",
    subtitle = "Eight-major-cell-type framework; deterministic plotting subset",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Broad cell type"
  ) +
  theme_qc +
  ggplot2::theme(
    legend.position = "right"
  )

save_plot_pair(
  plot_umap_broad,
  "QC1_UMAP_broad_types",
  width = 10,
  height = 7
)

# Candidate subtype UMAP, cluster labels retained.
subtype_centroids <- cell_dt[
  ,
  .(
    UMAPH_1 = stats::median(UMAPH_1),
    UMAPH_2 = stats::median(UMAPH_2),
    cluster = unique(cluster)[1]
  ),
  by = candidate_subtype
]

plot_umap_subtype <- ggplot2::ggplot(
  plot_cells,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = candidate_subtype
  )
) +
  ggplot2::geom_point(
    size = 0.06,
    alpha = 0.55,
    stroke = 0
  ) +
  ggplot2::geom_text(
    data = subtype_centroids,
    ggplot2::aes(
      x = UMAPH_1,
      y = UMAPH_2,
      label = paste0("C", cluster)
    ),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "bold"
  ) +
  ggplot2::labs(
    title = "GSE250130 candidate subtype annotation",
    subtitle = "Candidate labels are not locked; C12/C14/C20/C21 retain review flags",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Candidate subtype"
  ) +
  theme_qc +
  ggplot2::theme(
    legend.position = "right",
    legend.text = ggplot2::element_text(size = 6)
  )

save_plot_pair(
  plot_umap_subtype,
  "QC2_UMAP_subtypes",
  width = 13,
  height = 8
)

plot_major_profile <- ggplot2::ggplot(
  reshape(
    major_summary[
      ,
      c(
        "broad_cell_type",
        "cell_fraction",
        "published_cell_fraction"
      )
    ],
    varying = c(
      "cell_fraction",
      "published_cell_fraction"
    ),
    v.names = "fraction",
    timevar = "profile",
    times = c(
      "Current reanalysis",
      "Published study"
    ),
    direction = "long"
  ),
  ggplot2::aes(
    x = reorder(
      broad_cell_type,
      fraction
    ),
    y = fraction,
    fill = profile
  )
) +
  ggplot2::geom_col(
    position = "dodge"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::labs(
    title = "Major cell-type profile: contextual concordance",
    subtitle = "Published proportions are contextual QC, not an exact reproduction target",
    x = "Broad cell type",
    y = "Cell fraction",
    fill = "Profile"
  ) +
  theme_qc

save_plot_pair(
  plot_major_profile,
  "QC3_profile_context",
  width = 9,
  height = 6
)

sample_type_counts$sample_accession <- factor(
  sample_type_counts$sample_accession,
  levels = unique(
    locked_design$sample_accession
  )
)

plot_sample_composition <- ggplot2::ggplot(
  sample_type_counts,
  ggplot2::aes(
    x = sample_accession,
    y = within_sample_fraction,
    fill = broad_cell_type
  )
) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::labs(
    title = "Sample-level broad cell-type composition",
    subtitle = "Descriptive QC only; no composition hypothesis test",
    x = "Sample",
    y = "Within-sample fraction",
    fill = "Broad cell type"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 6
    )
  )

save_plot_pair(
  plot_sample_composition,
  "QC4_sample_composition",
  width = 12,
  height = 7
)

evidence_plot_data <- annotation_map
evidence_plot_data$cluster <- factor(
  evidence_plot_data$cluster,
  levels = as.character(0:21)
)

plot_evidence <- ggplot2::ggplot(
  evidence_plot_data,
  ggplot2::aes(
    x = cluster,
    y = observed_evidence_fraction_all_markers,
    fill = broad_annotation_status
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(
    yintercept = 0.60,
    linetype = "dashed"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::labs(
    title = "Cluster annotation evidence coverage",
    subtitle = "Dashed line marks the predefined 60% broad-evidence threshold",
    x = "Cluster",
    y = "Expected evidence genes observed",
    fill = "Annotation status"
  ) +
  theme_qc

save_plot_pair(
  plot_evidence,
  "QC5_evidence_coverage",
  width = 9,
  height = 5
)


# ----------------------------------------------------------------------
# 10. Annotation Gate
# ----------------------------------------------------------------------

expected_broad_types <- c(
  "Unciliated epithelial cells",
  "Ciliated epithelial cells",
  "Stromal cells",
  "Endothelial cells",
  "NK/T cells",
  "Myeloid cells",
  "B cells",
  "Mast cells"
)

technical_checks <- data.frame(
  check = c(
    "integration_gate_valid",
    "marker_gate_valid",
    "cell_count_locked",
    "sample_count_locked",
    "cluster_count_locked",
    "all_22_clusters_mapped",
    "all_eight_broad_types_present",
    "all_cells_receive_broad_annotation",
    "all_clusters_meet_broad_evidence_threshold",
    "review_clusters_preserved"
  ),
  pass = c(
    integration_gate$gate_status[[1]] %in%
      c("PASS", "PASS_WITH_REVIEW_FLAGS"),
    marker_gate$gate_status[[1]] ==
      "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED",
    nrow(cell_annotation) == 236636L,
    length(unique(cell_annotation$sample_accession)) == 28L,
    length(unique(cell_annotation$cluster)) == 22L,
    setequal(annotation_map$cluster, as.character(0:21)),
    setequal(
      unique(annotation_map$broad_cell_type),
      expected_broad_types
    ),
    !anyNA(cell_annotation$broad_cell_type),
    all(
      annotation_map$observed_evidence_fraction_all_markers >= 0.60
    ),
    setequal(
      annotation_map$cluster[
        annotation_map$review_flag
      ],
      c("12", "14", "20", "21")
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  technical_checks,
  file.path(
    run_qc_root,
    "03_gate_checks.csv"
  )
)

all_technical_pass <- all(technical_checks$pass)

gate_status <- if (all_technical_pass) {
  "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW"
} else {
  "HOLD_BROAD_ANNOTATION_QC_FAILED"
}

annotation_gate <- data.frame(
  dataset_id = dataset_id,
  input_version = input_version,
  marker_version = marker_version,
  annotation_version = annotation_version,
  gate_status = gate_status,
  cell_n = nrow(cell_annotation),
  sample_n = length(
    unique(cell_annotation$sample_accession)
  ),
  cluster_n = length(
    unique(cell_annotation$cluster)
  ),
  broad_cell_type_n = length(
    unique(cell_annotation$broad_cell_type)
  ),
  review_clusters = collapse_unique(
    annotation_map$cluster[
      annotation_map$review_flag
    ]
  ),
  broad_annotation_authorized =
    all_technical_pass,
  broad_VDR_localization_authorized =
    all_technical_pass,
  broad_pseudobulk_aggregation_authorized =
    all_technical_pass,
  broad_composition_descriptive_QC_authorized =
    all_technical_pass,
  composition_inferential_testing_authorized =
    FALSE,
  subtype_inferential_analysis_authorized =
    FALSE,
  subtype_annotation_locked = FALSE,
  result_lock_authorized = FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "ANNOTATION_GATE.csv"
)

write_csv_utf8(
  annotation_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1E Broad Annotation Result Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Input version: ", input_version),
  paste0("- Marker version: ", marker_version),
  paste0("- Annotation version: ", annotation_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Broad cell types: ",
    length(unique(cell_annotation$broad_cell_type))
  ),
  paste0(
    "- Annotated cells: ",
    format(nrow(cell_annotation), big.mark = ",")
  ),
  paste0(
    "- Review clusters: ",
    collapse_unique(
      annotation_map$cluster[
        annotation_map$review_flag
      ]
    )
  ),
  paste0(
    "- Broad VDR localization authorized: ",
    annotation_gate$broad_VDR_localization_authorized
  ),
  paste0(
    "- Broad pseudobulk aggregation authorized: ",
    annotation_gate$broad_pseudobulk_aggregation_authorized
  ),
  "- Composition inferential testing authorized: FALSE",
  "- Subtype inferential analysis authorized: FALSE",
  "- Subtype annotation locked: FALSE",
  "- Result lock authorized: FALSE"
)

card_path <- file.path(
  run_qc_root,
  "ANNOTATION_RESULT_CARD.md"
)

write_text(card_lines, card_path)

interpretation_lines <- c(
  "# STEP04E1E GSE250130 broad annotation interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  paste0(
    "- Broad cell types: ",
    collapse_unique(
      major_summary$broad_cell_type
    )
  ),
  "",
  "## Accepted broad annotation",
  "",
  "- Broad cell types follow the eight-major-cell-type framework reported for GSE250130.",
  "- All 22 clusters met the predefined marker-evidence threshold.",
  "- Broad annotations may be used for VDR program localization and donor-by-cell-type aggregation.",
  "",
  "## Preserved uncertainty",
  "",
  "- Candidate subtype labels are not locked.",
  "- Cluster 12 is accepted as B cells, but its sample-mixing pattern remains under review.",
  "- Cluster 14 is accepted as mast cells, but its abundance is strongly influenced by one sample.",
  "- Cluster 20 is accepted as unciliated epithelial cells and provisionally described as ESR1/LIFR/LPAR3 luminal-like epithelial cells.",
  "- Cluster 21 is accepted as stromal cells and provisionally described as activated perivascular stromal cells.",
  "",
  "## Composition rule",
  "",
  "- Sample-level cell fractions are retained as descriptive QC.",
  "- No group or timepoint composition hypothesis test is authorized at this stage.",
  "",
  "## Next step",
  "",
  "- Use the broad annotation map for VDR RNA-program localization.",
  "- Aggregate raw counts by donor and broad cell type.",
  "- Perform donor-level QC before any clinical or temporal comparison.",
  "- Process GSE183837 only after the GSE250130 annotation and aggregation workflow is stable."
)

interpretation_path <- file.path(
  run_qc_root,
  "ANNOTATION_INTERPRETATION.md"
)

write_text(
  interpretation_lines,
  interpretation_path
)


# ----------------------------------------------------------------------
# 11. Methods and Results drafts
# ----------------------------------------------------------------------

methods_lines <- c(
  "# STEP04E1E Methods draft",
  "",
  paste(
    "Clusters from the resolution-0.4 solution were annotated using an",
    "evidence-reviewed hierarchical procedure."
  ),
  paste(
    "Broad cell-type labels were aligned to the eight major populations",
    "reported in the source GSE250130 study: unciliated epithelial, ciliated",
    "epithelial, stromal, endothelial, NK/T, myeloid, B and mast cells."
  ),
  paste(
    "For each cluster, descriptive positive markers, canonical marker-panel",
    "support and sample-mixing QC were reviewed jointly."
  ),
  paste(
    "Broad labels were accepted only when at least 60% of a predefined",
    "cluster-specific evidence-gene set was recovered among the descriptive",
    "marker candidates."
  ),
  paste(
    "Detailed subtype labels were retained as candidates and were not used",
    "for inferential testing."
  ),
  paste(
    "Clusters previously flagged for sample dominance or restricted donor",
    "coverage retained explicit review flags in the annotation manifest."
  ),
  paste(
    "Cell-type composition summaries were generated for quality-control",
    "purposes only; no cell-fraction hypothesis testing was performed."
  )
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Methods_draft.md"
  )
)

summary_lookup <- setNames(
  as.numeric(major_summary$cell_fraction),
  as.character(major_summary$broad_cell_type)
)

get_fraction <- function(cell_type) {
  value <- unname(summary_lookup[cell_type])

  if (
    length(value) != 1L ||
      is.na(value) ||
      !is.finite(value)
  ) {
    stop(
      "Missing or invalid broad-cell-type fraction for: ",
      cell_type
    )
  }

  as.numeric(value)
}

results_lines <- c(
  "# STEP04E1E Results draft",
  "",
  paste(
    "The 236,636 integrated cells were assigned to eight broad endometrial",
    "cell populations using cluster-level marker evidence."
  ),
  paste0(
    "NK/T cells represented ",
    round(100 * get_fraction("NK/T cells"), 1),
    "% of cells, followed by stromal cells (",
    round(100 * get_fraction("Stromal cells"), 1),
    "%) and unciliated epithelial cells (",
    round(
      100 * get_fraction("Unciliated epithelial cells"),
      1
    ),
    "%)."
  ),
  paste0(
    "Myeloid cells accounted for ",
    round(100 * get_fraction("Myeloid cells"), 1),
    "%, ciliated epithelial cells for ",
    round(
      100 * get_fraction("Ciliated epithelial cells"),
      1
    ),
    "%, B cells for ",
    round(100 * get_fraction("B cells"), 1),
    "%, endothelial cells for ",
    round(100 * get_fraction("Endothelial cells"), 1),
    "% and mast cells for ",
    round(100 * get_fraction("Mast cells"), 1),
    "%."
  ),
  paste(
    "All clusters satisfied the predefined broad-lineage marker-evidence",
    "criterion."
  ),
  paste(
    "Four clusters retained review flags: B-cell cluster 12 and mast-cell",
    "cluster 14 because of sample-mixing characteristics, and rare clusters",
    "20 and 21 because of restricted donor coverage and provisional subtype",
    "interpretation."
  ),
  paste(
    "Broad labels were therefore authorized for subsequent VDR-program",
    "localization and donor-level aggregation, whereas subtype-level",
    "inference remained deferred."
  )
)


write_text(
  results_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Results_draft.md"
  )
)


# ----------------------------------------------------------------------
# 12. Manifests, sessionInfo and ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_names <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

package_table <- data.frame(
  package = package_names,
  version = vapply(
    package_names,
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(
    run_qc_root,
    "PACKAGE_VERSIONS.csv"
  )
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- unique(c(
  list.files(
    annotation_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Results_draft.md"
  )
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(artifact_files),
  size_bytes = file.info(artifact_files)$size,
  md5 = unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(
    run_qc_root,
    "ARTIFACT_MANIFEST.csv"
  )
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  marker_version = marker_version,
  annotation_version = annotation_version,
  supersedes_partial_annotation_version =
    "ANNOTATION_v1_1_FIXED_PATH_LENGTH_FAILURE",
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_seconds = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "secs"
      )
    ),
    2
  ),
  gate_status = gate_status,
  cell_n = nrow(cell_annotation),
  sample_n = length(
    unique(cell_annotation$sample_accession)
  ),
  cluster_n = length(
    unique(cell_annotation$cluster)
  ),
  broad_cell_type_n = length(
    unique(cell_annotation$broad_cell_type)
  ),
  large_Seurat_checkpoint_loaded = FALSE,
  broad_annotation_authorized =
    annotation_gate$broad_annotation_authorized,
  subtype_annotation_locked = FALSE,
  composition_inferential_testing_authorized =
    FALSE,
  result_lock_authorized = FALSE,
  annotation_root =
    normalize_slashes(annotation_root),
  qc_root =
    normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "broad_cell_types=",
  length(
    unique(cell_annotation$broad_cell_type)
  )
)

small_files <- unique(c(
  list.files(
    annotation_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1E_annotation_Results_draft.md"
  ),
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <=
      200 * 1024^2
]

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = small_files,
        recurse = FALSE,
        include_directories = FALSE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1E GSE250130 broad annotation QC completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Annotated cells: ",
  format(nrow(cell_annotation), big.mark = ","),
  "\n",
  sep = ""
)
cat(
  "Broad cell types: ",
  length(unique(cell_annotation$broad_cell_type)),
  "\n",
  sep = ""
)
cat(
  "Review clusters: ",
  collapse_unique(
    annotation_map$cluster[
      annotation_map$review_flag
    ]
  ),
  "\n",
  sep = ""
)
cat(
  "Annotation root: ",
  annotation_root,
  "\n",
  sep = ""
)
cat(
  "QC root: ",
  run_qc_root,
  "\n",
  sep = ""
)

if (zip_ok) {
  cat(
    "Upload this ZIP:\n",
    handoff_zip,
    "\n"
  )
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe annotation and QC folders are complete.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 6: 60_step04E1F_GSE250130_VDR_localization_pseudobulk_v1
# ==============================================================================

# ======================================================================
# 60_step04E1F_GSE250130_VDR_localization_pseudobulk_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# 1. Aggregate raw counts by sample/donor x accepted broad cell type.
# 2. Preserve all raw pseudobulk counts for later donor-level modelling.
# 3. Localize four locked VDR RNA programs across broad cell types using
#    descriptive signed program scores.
# 4. Generate pseudobulk input QC and model-eligibility tables.
#
# LOCKED VDR PROGRAMS
# - RNA_ALL_436
# - OCCUPANCY_SUPPORTED_284
# - SHARED_OCCUPANCY_46
# - PROMOTER_OCCUPANCY_84
#
# SCORE DEFINITION
# - Raw pseudobulk counts are library-size normalized to log2(CPM + 1).
# - Each VDR-program gene is standardized across all donor x broad-type
#   pseudobulk profiles.
# - The standardized expression is multiplied by the sign of the locked
#   VDR interaction logFC.
# - The program score is the mean signed standardized expression.
#
# INTERPRETATION
# - Positive score: expression pattern is aligned with the direction of the
#   locked calcitriol-dependent VDR interaction program.
# - Negative score: expression pattern is relatively opposed to that direction.
# - Scores are descriptive localization evidence, not causal activation and
#   not donor-level differential-expression inference.
#
# MEMORY / RESUME STRATEGY
# - The integrated checkpoint is loaded once.
# - Each counts layer is aggregated separately and saved immediately.
# - Re-running the same script skips completed layer aggregation files.
# - Counts layers remain split in the Seurat object; no JoinLayers(counts).
# - The large Seurat object is removed before score calculation.
#
# INPUT GATE
# - Requires STEP04E1E:
#     PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW
#     broad_VDR_localization_authorized = TRUE
#     broad_pseudobulk_aggregation_authorized = TRUE
#
# OUTPUT STATUS
# - Pseudobulk raw counts: analysis-ready after Gate
# - VDR localization: descriptive candidate evidence
# - No output is written to 06_locked_results.
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)
options(timeout = max(14400L, getOption("timeout")))
options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260805)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <-
  "60_step04E1F_GSE250130_VDR_localization_pseudobulk_v1.R"
script_version <-
  "STEP04E1F_GSE250130_VDR_LOCALIZATION_PSEUDOBULK_v1.0"
input_version <- "SCI_v1_1_memory_safe"
annotation_version <- "ANN_v1_2"
pseudobulk_version <- "PBVDR_v1_0"

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1F_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(
  path,
  max_chars = 245L
) {
  normalized <- normalize_slashes(path)
  path_n <- nchar(normalized, type = "chars")

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Output path is too long for safe Windows/R writing (",
      path_n,
      " characters; safe limit ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

format_bytes <- function(bytes) {
  if (length(bytes) == 0L || is.na(bytes)) return("")
  if (bytes < 1024) return(sprintf("%.0f B", bytes))
  if (bytes < 1024^2) return(sprintf("%.1f KB", bytes / 1024))
  if (bytes < 1024^3) return(sprintf("%.1f MB", bytes / 1024^2))
  sprintf("%.2f GB", bytes / 1024^3)
}

sanitize_file_token <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  sub("^_+|_+$", "", x)
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") return(NA_real_)

  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"$o=Get-CimInstance Win32_OperatingSystem; ',
        '[math]::Round($o.FreePhysicalMemory/1MB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )

  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

get_free_disk_gb <- function(drive_letter = "C") {
  if (.Platform$OS.type != "windows") return(NA_real_)

  lines <- tryCatch(
    system(
      paste0(
        'powershell -NoProfile -Command ',
        '"[math]::Round((Get-PSDrive -Name ',
        drive_letter,
        ').Free/1GB,2)"'
      ),
      intern = TRUE,
      ignore.stderr = TRUE
    ),
    error = function(e) character()
  )

  if (!length(lines)) return(NA_real_)
  suppressWarnings(as.numeric(trimws(lines[1])))
}

row_zscore_dense <- function(x) {
  x <- as.matrix(x)
  row_mean <- rowMeans(x, na.rm = TRUE)
  row_sd <- apply(x, 1L, stats::sd, na.rm = TRUE)

  output <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )

  valid <- is.finite(row_sd) & row_sd > 0

  output[valid, ] <- (
    x[valid, , drop = FALSE] -
      row_mean[valid]
  ) / row_sd[valid]

  output
}

save_plot_pair <- function(plot_object, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))

  assert_windows_path_length(png_path)
  assert_windows_path_length(pdf_path)

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")

  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages and project paths
# ----------------------------------------------------------------------

required_packages <- c(
  "SeuratObject",
  "Matrix",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(step04_root, "05_execution")
results_root <- file.path(step04_root, "06_results")
qc_base_root <- file.path(step04_root, "08_qc")
parameter_root <- file.path(step04_root, "03_parameters")
environment_root <- file.path(step04_root, "04_environment")
manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

checkpoint_path <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "checkpoints",
  dataset_id,
  paste0(dataset_id, "_integrated_checkpoint.rds")
)

checkpoint_status_path <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "checkpoints",
  dataset_id,
  paste0(dataset_id, "_checkpoint_status.csv")
)

annotation_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "annotations",
  dataset_id,
  annotation_version
)

annotation_path <- file.path(
  annotation_root,
  "tables",
  "cell_annotation_metadata.csv"
)

candidate_path <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons",
  "06_results",
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  "RNA_OCCUPANCY_EVIDENCE_v1.0",
  "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv"
)

output_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "pb",
  dataset_id,
  pseudobulk_version
)

per_layer_root <- file.path(output_root, "per_layer")
object_root <- file.path(output_root, "objects")
table_root <- file.path(output_root, "tables")
source_data_root <- file.path(output_root, "source_data")

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1F_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_pbQC"
  )
)

figure_root <- file.path(run_qc_root, "figures")

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  output_root,
  per_layer_root,
  object_root,
  table_root,
  source_data_root,
  run_qc_root,
  figure_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)

planned_paths <- data.frame(
  object = c(
    "output_root",
    "per_layer_root",
    "object_root",
    "table_root",
    "source_data_root",
    "run_qc_root",
    "figure_root",
    "log_path",
    "run_manifest_path",
    "parameter_path",
    "session_path",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    output_root,
    per_layer_root,
    object_root,
    table_root,
    source_data_root,
    run_qc_root,
    figure_root,
    log_path,
    run_manifest_path,
    parameter_path,
    session_path,
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_paths$character_n <- nchar(
  planned_paths$path,
  type = "chars"
)

if (
  .Platform$OS.type == "windows" &&
    any(planned_paths$character_n > 245L)
) {
  bad <- planned_paths[
    planned_paths$character_n > 245L,
    ,
    drop = FALSE
  ]

  stop(
    "One or more planned paths are too long:\n",
    paste(
      paste0(
        bad$character_n,
        " chars: ",
        bad$path
      ),
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------------------
# 3. Validate upstream Gates and inputs
# ----------------------------------------------------------------------

required_inputs <- c(
  checkpoint_path,
  checkpoint_status_path,
  annotation_path,
  candidate_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

annotation_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

annotation_qc_dirs <- annotation_qc_dirs[
  grepl(
    "^STEP04E1E_[0-9]{8}_[0-9]{6}_annQC$",
    basename(annotation_qc_dirs),
    perl = TRUE
  )
]

if (!length(annotation_qc_dirs)) {
  stop("No completed STEP04E1E annotation-QC folder found.")
}

latest_annotation_qc <- annotation_qc_dirs[
  which.max(file.info(annotation_qc_dirs)$mtime)
]

annotation_gate_path <- file.path(
  latest_annotation_qc,
  "ANNOTATION_GATE.csv"
)

annotation_gate <- safe_read_csv(annotation_gate_path)
checkpoint_status <- safe_read_csv(checkpoint_status_path)

if (
  is.null(annotation_gate) ||
    nrow(annotation_gate) != 1L ||
    annotation_gate$gate_status[[1]] !=
      "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW" ||
    !isTRUE(
      annotation_gate$broad_VDR_localization_authorized[[1]]
    ) ||
    !isTRUE(
      annotation_gate$broad_pseudobulk_aggregation_authorized[[1]]
    )
) {
  stop(
    "STEP04E1E Gate does not authorize VDR localization and pseudobulk aggregation."
  )
}

if (
  is.null(checkpoint_status) ||
    nrow(checkpoint_status) != 1L ||
    checkpoint_status$checkpoint_status[[1]] != "COMPLETE" ||
    checkpoint_status$join_policy[[1]] !=
      "DATA_ONLY_COUNTS_REMAIN_SPLIT"
) {
  stop("Integrated checkpoint status is absent or invalid.")
}

annotation <- data.table::fread(
  annotation_path,
  data.table = FALSE,
  showProgress = FALSE
)

candidate_table <- safe_read_csv(candidate_path)

required_annotation_columns <- c(
  "cell_barcode",
  "sample_accession",
  "donor_id",
  "broad_cell_type",
  "clinical_group",
  "LH_day",
  "temporal_fertile_eligible",
  "LH7_RIF_comparison_eligible",
  "LH7_reference_group",
  "LH7_fertile_subgroup"
)

missing_annotation_columns <- setdiff(
  required_annotation_columns,
  colnames(annotation)
)

if (length(missing_annotation_columns)) {
  stop(
    "Annotation table is missing columns: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

required_candidate_columns <- c(
  "feature_id",
  "interaction_logFC",
  "occupancy_supported",
  "shared_occupancy_supported",
  "any_promoter_overlap_support"
)

missing_candidate_columns <- setdiff(
  required_candidate_columns,
  colnames(candidate_table)
)

if (length(missing_candidate_columns)) {
  stop(
    "Locked VDR candidate table is missing columns: ",
    paste(missing_candidate_columns, collapse = ", ")
  )
}

if (
  nrow(annotation) != 236636L ||
    length(unique(annotation$sample_accession)) != 28L ||
    length(unique(annotation$broad_cell_type)) != 8L ||
    anyDuplicated(annotation$cell_barcode) > 0L
) {
  stop("Annotation identity checks failed.")
}

if (
  nrow(candidate_table) != 436L ||
    length(unique(candidate_table$feature_id)) != 436L ||
    sum(candidate_table$occupancy_supported) != 284L ||
    sum(candidate_table$shared_occupancy_supported) != 46L ||
    sum(candidate_table$any_promoter_overlap_support) != 84L
) {
  stop("Locked 436/284/46/84 VDR candidate table failed validation.")
}

free_memory_gb <- get_free_memory_gb()
free_disk_gb <- get_free_disk_gb("C")

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
  stop(
    "Free physical memory is below 18 GB. Restart RStudio and rerun in a clean session."
  )
}

if (is.finite(free_disk_gb) && free_disk_gb < 12) {
  stop("Free C: drive space is below 12 GB.")
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  annotation_version = annotation_version,
  pseudobulk_version = pseudobulk_version,
  dataset_id = dataset_id,
  aggregation_level = "sample_donor_x_broad_cell_type",
  raw_count_layer_policy = "PROCESS_SPLIT_COUNTS_LAYERS_SEPARATELY",
  minimum_cells_for_score = 20L,
  minimum_cells_preferred = 50L,
  localization_normalization = "LOG2_LIBRARY_SIZE_CPM_PLUS_1",
  program_score =
    "MEAN_GENE_ZSCORE_MULTIPLIED_BY_SIGN_INTERACTION_LOGFC",
  score_role =
    "DESCRIPTIVE_LOCALIZATION_NOT_CAUSAL_OR_INFERENTIAL",
  temporal_model_eligibility =
    "AT_LEAST_2_FERTILE_DONORS_AT_EACH_LH3_LH5_LH7_LH9_LH11_AND_TOTAL_GE_12",
  LH7_model_eligibility =
    "AT_LEAST_3_FERTILE_LH7_AND_3_RIF_LH7_PSEUDOBULKS",
  composition_testing_authorized = FALSE,
  subtype_analysis_authorized = FALSE,
  random_seed = 20260805L,
  result_lock_authorized = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_csv_utf8(
  planned_paths,
  file.path(run_qc_root, "00_path_audit.csv")
)

resource_snapshot <- data.frame(
  metric = c(
    "free_physical_memory_GB_before_load",
    "free_C_drive_GB_before_run",
    "checkpoint_size_GB",
    "checkpoint_md5",
    "annotation_cell_n",
    "annotation_sample_n",
    "annotation_broad_type_n",
    "candidate_gene_n",
    "occupancy_supported_gene_n",
    "shared_occupancy_gene_n",
    "promoter_occupancy_gene_n",
    "annotation_gate"
  ),
  value = c(
    free_memory_gb,
    free_disk_gb,
    round(file.info(checkpoint_path)$size / 1024^3, 3),
    unname(tools::md5sum(checkpoint_path)),
    nrow(annotation),
    length(unique(annotation$sample_accession)),
    length(unique(annotation$broad_cell_type)),
    nrow(candidate_table),
    sum(candidate_table$occupancy_supported),
    sum(candidate_table$shared_occupancy_supported),
    sum(candidate_table$any_promoter_overlap_support),
    annotation_gate$gate_status[[1]]
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resource_snapshot,
  file.path(run_qc_root, "01_input_resource_snapshot.csv")
)

log_event(
  "Starting",
  run_id,
  "annotation_gate=",
  annotation_gate$gate_status[[1]],
  "free_memory_GB=",
  free_memory_gb,
  "free_disk_GB=",
  free_disk_gb
)


# ----------------------------------------------------------------------
# 4. Load integrated object and validate cell identity
# ----------------------------------------------------------------------

load_start <- Sys.time()
object <- readRDS(checkpoint_path)

log_event(
  "Loaded checkpoint",
  "elapsed_min=",
  round(
    as.numeric(
      difftime(Sys.time(), load_start, units = "mins")
    ),
    2
  )
)

if (!inherits(object, "Seurat")) {
  stop("Checkpoint is not a Seurat object.")
}

if (
  ncol(object) != 236636L ||
    nrow(object) != 36601L
) {
  stop(
    "Unexpected checkpoint dimensions: ",
    nrow(object),
    " genes x ",
    ncol(object),
    " cells"
  )
}

object_cells <- colnames(object)
annotation_cells <- as.character(annotation$cell_barcode)

if (!setequal(object_cells, annotation_cells)) {
  stop("Checkpoint cells and annotation cells are not identical.")
}

annotation <- annotation[
  match(object_cells, annotation$cell_barcode),
  ,
  drop = FALSE
]

if (!identical(as.character(annotation$cell_barcode), object_cells)) {
  stop("Failed to align annotation rows to object cells.")
}

annotation_lookup <- annotation
rownames(annotation_lookup) <- annotation_lookup$cell_barcode

SeuratObject::DefaultAssay(object) <- "RNA"

count_layers <- SeuratObject::Layers(
  object[["RNA"]],
  search = "^counts"
)

if (length(count_layers) != 28L) {
  stop(
    "Expected 28 split counts layers; found ",
    length(count_layers),
    "."
  )
}

write_csv_utf8(
  data.frame(
    layer_index = seq_along(count_layers),
    layer = count_layers,
    stringsAsFactors = FALSE
  ),
  file.path(run_qc_root, "02_counts_layer_manifest.csv")
)


# ----------------------------------------------------------------------
# 5. Resumable per-layer raw-count aggregation
# ----------------------------------------------------------------------

layer_status_rows <- list()
processed_cell_lists <- list()
raw_count_total_from_layers <- 0

for (layer_index in seq_along(count_layers)) {
  layer_name <- count_layers[layer_index]
  layer_token <- sanitize_file_token(layer_name)

  layer_rds <- file.path(
    per_layer_root,
    paste0(
      sprintf("%02d", layer_index),
      "_",
      layer_token,
      ".rds"
    )
  )

  layer_status_path <- file.path(
    per_layer_root,
    paste0(
      sprintf("%02d", layer_index),
      "_",
      layer_token,
      "_status.csv"
    )
  )

  existing_status <- safe_read_csv(layer_status_path)

  if (
    file.exists(layer_rds) &&
      !is.null(existing_status) &&
      nrow(existing_status) == 1L &&
      existing_status$status[[1]] == "COMPLETE"
  ) {
    log_event(
      "SKIP completed counts layer",
      layer_index,
      "/",
      length(count_layers),
      layer_name
    )

    saved <- readRDS(layer_rds)

    layer_status_rows[[length(layer_status_rows) + 1L]] <-
      existing_status

    processed_cell_lists[[length(processed_cell_lists) + 1L]] <-
      saved$cell_barcodes

    raw_count_total_from_layers <-
      raw_count_total_from_layers +
      as.numeric(saved$raw_count_total)

    rm(saved)
    gc()
    next
  }

  log_event(
    "START counts layer",
    layer_index,
    "/",
    length(count_layers),
    layer_name
  )

  layer_start <- Sys.time()
  layer_error <- ""

  result <- tryCatch(
    {
      count_matrix <- SeuratObject::LayerData(
        object = object,
        assay = "RNA",
        layer = layer_name
      )

      if (!inherits(count_matrix, "sparseMatrix")) {
        count_matrix <- Matrix::Matrix(
          count_matrix,
          sparse = TRUE
        )
      }

      layer_cells <- colnames(count_matrix)

      if (
        is.null(layer_cells) ||
          !length(layer_cells)
      ) {
        stop("Counts layer has no cell names.")
      }

      layer_annotation <- annotation_lookup[
        layer_cells,
        ,
        drop = FALSE
      ]

      if (
        nrow(layer_annotation) != length(layer_cells) ||
          anyNA(layer_annotation$sample_accession) ||
          anyNA(layer_annotation$broad_cell_type)
      ) {
        stop("Counts-layer cells could not be fully annotated.")
      }

      group_id <- paste(
        layer_annotation$sample_accession,
        layer_annotation$broad_cell_type,
        sep = "||"
      )

      group_factor <- factor(
        group_id,
        levels = unique(group_id)
      )

      membership <- Matrix::sparseMatrix(
        i = seq_along(layer_cells),
        j = as.integer(group_factor),
        x = 1,
        dims = c(
          length(layer_cells),
          nlevels(group_factor)
        ),
        dimnames = list(
          layer_cells,
          levels(group_factor)
        )
      )

      pb_counts <- count_matrix %*% membership

      colnames(pb_counts) <- levels(group_factor)

      cell_count_table <- as.data.frame(
        table(group_id),
        stringsAsFactors = FALSE
      )

      colnames(cell_count_table) <- c(
        "pb_id",
        "cell_n"
      )

      cell_count_table$cell_n <- as.integer(
        cell_count_table$cell_n
      )

      group_rows <- match(
        cell_count_table$pb_id,
        group_id
      )

      pb_metadata <- layer_annotation[
        group_rows,
        c(
          "sample_accession",
          "donor_id",
          "broad_cell_type",
          "clinical_group",
          "LH_day",
          "temporal_fertile_eligible",
          "LH7_RIF_comparison_eligible",
          "LH7_reference_group",
          "LH7_fertile_subgroup",
          "analysis_role"
        ),
        drop = FALSE
      ]

      pb_metadata$pb_id <- cell_count_table$pb_id
      pb_metadata$cell_n <- cell_count_table$cell_n
      pb_metadata$source_counts_layer <- layer_name

      pb_metadata <- pb_metadata[
        match(colnames(pb_counts), pb_metadata$pb_id),
        ,
        drop = FALSE
      ]

      if (
        anyNA(pb_metadata$pb_id) ||
          !identical(
            as.character(pb_metadata$pb_id),
            colnames(pb_counts)
          )
      ) {
        stop("Pseudobulk metadata did not align with count columns.")
      }

      raw_count_total <- as.numeric(
        sum(count_matrix)
      )

      aggregated_count_total <- as.numeric(
        sum(pb_counts)
      )

      if (!isTRUE(all.equal(
        raw_count_total,
        aggregated_count_total,
        tolerance = 0
      ))) {
        stop(
          "Raw-count total changed during aggregation: ",
          raw_count_total,
          " vs ",
          aggregated_count_total
        )
      }

      saved_object <- list(
        layer_index = layer_index,
        layer_name = layer_name,
        counts = pb_counts,
        metadata = pb_metadata,
        cell_barcodes = layer_cells,
        raw_count_total = raw_count_total,
        gene_n = nrow(pb_counts),
        pseudobulk_n = ncol(pb_counts)
      )

      saveRDS(
        saved_object,
        layer_rds,
        compress = FALSE
      )

      list(
        saved_object = saved_object,
        raw_count_total = raw_count_total
      )
    },
    error = function(e) {
      layer_error <<- conditionMessage(e)
      NULL
    }
  )

  elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        layer_start,
        units = "mins"
      )
    ),
    2
  )

  if (is.null(result)) {
    status <- data.frame(
      layer_index = layer_index,
      layer_name = layer_name,
      status = "FAILED",
      cell_n = NA_integer_,
      sample_n = NA_integer_,
      broad_type_n = NA_integer_,
      pseudobulk_n = NA_integer_,
      raw_count_total = NA_real_,
      elapsed_minutes = elapsed_minutes,
      error_message = layer_error,
      output_rds = normalize_slashes(layer_rds),
      stringsAsFactors = FALSE
    )

    write_csv_utf8(status, layer_status_path)

    layer_status_rows[[length(layer_status_rows) + 1L]] <-
      status

    log_event(
      "FAILED counts layer",
      layer_name,
      "error=",
      layer_error
    )

    stop(
      "Counts-layer aggregation failed. Completed layers are safely saved. ",
      "Re-run the same script after reviewing the error."
    )
  }

  saved_object <- result$saved_object

  status <- data.frame(
    layer_index = layer_index,
    layer_name = layer_name,
    status = "COMPLETE",
    cell_n = length(saved_object$cell_barcodes),
    sample_n = length(
      unique(saved_object$metadata$sample_accession)
    ),
    broad_type_n = length(
      unique(saved_object$metadata$broad_cell_type)
    ),
    pseudobulk_n = ncol(saved_object$counts),
    raw_count_total = saved_object$raw_count_total,
    elapsed_minutes = elapsed_minutes,
    error_message = "",
    output_rds = normalize_slashes(layer_rds),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(status, layer_status_path)

  layer_status_rows[[length(layer_status_rows) + 1L]] <-
    status

  processed_cell_lists[[length(processed_cell_lists) + 1L]] <-
    saved_object$cell_barcodes

  raw_count_total_from_layers <-
    raw_count_total_from_layers +
    saved_object$raw_count_total

  log_event(
    "DONE counts layer",
    layer_index,
    "/",
    length(count_layers),
    layer_name,
    "cells=",
    length(saved_object$cell_barcodes),
    "pseudobulks=",
    ncol(saved_object$counts),
    "elapsed_min=",
    elapsed_minutes
  )

  rm(
    result,
    saved_object,
    count_matrix,
    membership,
    pb_counts
  )
  gc()
}

layer_status <- do.call(
  rbind,
  layer_status_rows
)

layer_status <- layer_status[
  order(layer_status$layer_index),
  ,
  drop = FALSE
]

write_csv_utf8(
  layer_status,
  file.path(table_root, "01_layer_aggregation_status.csv")
)

if (
  nrow(layer_status) != 28L ||
    any(layer_status$status != "COMPLETE")
) {
  stop("Not all 28 counts layers completed aggregation.")
}

processed_cells <- unlist(
  processed_cell_lists,
  use.names = FALSE
)

cell_identity_audit <- data.frame(
  metric = c(
    "processed_cell_n",
    "unique_processed_cell_n",
    "annotation_cell_n",
    "duplicated_processed_cell_n",
    "missing_annotation_cell_n",
    "unexpected_processed_cell_n"
  ),
  value = c(
    length(processed_cells),
    length(unique(processed_cells)),
    nrow(annotation),
    sum(duplicated(processed_cells)),
    length(
      setdiff(
        annotation$cell_barcode,
        processed_cells
      )
    ),
    length(
      setdiff(
        processed_cells,
        annotation$cell_barcode
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  cell_identity_audit,
  file.path(run_qc_root, "03_cell_identity_audit.csv")
)

if (
  length(processed_cells) != 236636L ||
    length(unique(processed_cells)) != 236636L ||
    anyDuplicated(processed_cells) > 0L ||
    !setequal(processed_cells, annotation$cell_barcode)
) {
  stop("Counts-layer cell identity audit failed.")
}


# ----------------------------------------------------------------------
# 6. Remove large object and combine saved pseudobulks
# ----------------------------------------------------------------------

rm(
  object,
  annotation_lookup,
  processed_cell_lists,
  processed_cells
)
gc()

layer_rds_files <- layer_status$output_rds

saved_layers <- lapply(
  layer_rds_files,
  readRDS
)

reference_genes <- rownames(
  saved_layers[[1]]$counts
)

gene_identity_pass <- vapply(
  saved_layers,
  function(x) {
    identical(
      rownames(x$counts),
      reference_genes
    )
  },
  logical(1)
)

if (!all(gene_identity_pass)) {
  stop("Gene order differs across saved counts-layer pseudobulks.")
}

pb_matrices <- lapply(
  saved_layers,
  function(x) x$counts
)

pb_metadata_list <- lapply(
  saved_layers,
  function(x) x$metadata
)

pb_counts <- do.call(
  cbind,
  pb_matrices
)

pb_metadata <- do.call(
  rbind,
  pb_metadata_list
)

rownames(pb_metadata) <- NULL

if (
  anyDuplicated(colnames(pb_counts)) > 0L ||
    anyDuplicated(pb_metadata$pb_id) > 0L
) {
  stop("Duplicated pseudobulk IDs were created.")
}

pb_metadata <- pb_metadata[
  match(colnames(pb_counts), pb_metadata$pb_id),
  ,
  drop = FALSE
]

if (
  anyNA(pb_metadata$pb_id) ||
    !identical(
      colnames(pb_counts),
      as.character(pb_metadata$pb_id)
    )
) {
  stop("Combined pseudobulk metadata alignment failed.")
}

library_size <- as.numeric(
  Matrix::colSums(pb_counts)
)

detected_gene_n <- as.numeric(
  Matrix::colSums(pb_counts > 0)
)

pb_metadata$library_size <- library_size
pb_metadata$detected_gene_n <- detected_gene_n
pb_metadata$minimum_cell_threshold_pass <- (
  pb_metadata$cell_n >= 20L
)
pb_metadata$preferred_cell_threshold_pass <- (
  pb_metadata$cell_n >= 50L
)

combined_raw_count_total <- as.numeric(
  sum(pb_counts)
)

if (!isTRUE(all.equal(
  raw_count_total_from_layers,
  combined_raw_count_total,
  tolerance = 0
))) {
  stop(
    "Combined raw-count total mismatch: ",
    raw_count_total_from_layers,
    " vs ",
    combined_raw_count_total
  )
}

pb_counts_path <- file.path(
  object_root,
  "pb_raw_counts.rds"
)

pb_metadata_path <- file.path(
  table_root,
  "02_pb_metadata.csv"
)

saveRDS(
  pb_counts,
  pb_counts_path,
  compress = FALSE
)

write_csv_utf8(
  pb_metadata,
  pb_metadata_path
)

combined_object_status <- data.frame(
  dataset_id = dataset_id,
  pseudobulk_version = pseudobulk_version,
  status = "COMPLETE",
  gene_n = nrow(pb_counts),
  pseudobulk_n = ncol(pb_counts),
  sample_n = length(
    unique(pb_metadata$sample_accession)
  ),
  donor_n = length(
    unique(pb_metadata$donor_id)
  ),
  broad_type_n = length(
    unique(pb_metadata$broad_cell_type)
  ),
  aggregated_cell_n = sum(pb_metadata$cell_n),
  raw_count_total = combined_raw_count_total,
  counts_file_size_bytes = file.info(pb_counts_path)$size,
  counts_file_md5 = unname(tools::md5sum(pb_counts_path)),
  counts_path = normalize_slashes(pb_counts_path),
  metadata_path = normalize_slashes(pb_metadata_path),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  combined_object_status,
  file.path(object_root, "pb_status.csv")
)

rm(
  saved_layers,
  pb_matrices,
  pb_metadata_list
)
gc()


# ----------------------------------------------------------------------
# 7. Pseudobulk QC and model eligibility
# ----------------------------------------------------------------------

broad_types <- sort(
  unique(pb_metadata$broad_cell_type)
)

pb_qc_by_type <- do.call(
  rbind,
  lapply(
    broad_types,
    function(cell_type) {
      one <- pb_metadata[
        pb_metadata$broad_cell_type == cell_type,
        ,
        drop = FALSE
      ]

      data.frame(
        broad_cell_type = cell_type,
        pseudobulk_n = nrow(one),
        sample_n = length(unique(one$sample_accession)),
        donor_n = length(unique(one$donor_id)),
        median_cell_n = stats::median(one$cell_n),
        min_cell_n = min(one$cell_n),
        median_library_size = stats::median(one$library_size),
        min_library_size = min(one$library_size),
        median_detected_gene_n =
          stats::median(one$detected_gene_n),
        pseudobulk_n_cells_ge_20 =
          sum(one$cell_n >= 20L),
        pseudobulk_n_cells_ge_50 =
          sum(one$cell_n >= 50L),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_csv_utf8(
  pb_qc_by_type,
  file.path(table_root, "03_pb_QC_by_type.csv")
)

eligibility_rows <- list()

for (cell_type in broad_types) {
  one <- pb_metadata[
    pb_metadata$broad_cell_type == cell_type &
      pb_metadata$cell_n >= 20L,
    ,
    drop = FALSE
  ]

  fertile_temporal <- one[
    one$clinical_group == "FERTILE" &
      one$temporal_fertile_eligible,
    ,
    drop = FALSE
  ]

  day_counts <- table(
    factor(
      fertile_temporal$LH_day,
      levels = c(3, 5, 7, 9, 11)
    )
  )

  temporal_eligible <- (
    nrow(fertile_temporal) >= 12L &&
      all(day_counts >= 2L)
  )

  lh7 <- one[
    one$LH7_RIF_comparison_eligible,
    ,
    drop = FALSE
  ]

  fertile_lh7_n <- sum(
    lh7$LH7_reference_group == "FERTILE_LH7",
    na.rm = TRUE
  )

  rif_lh7_n <- sum(
    lh7$LH7_reference_group == "RIF",
    na.rm = TRUE
  )

  lh7_eligible <- (
    fertile_lh7_n >= 3L &&
      rif_lh7_n >= 3L
  )

  localization_eligible <- (
    nrow(one) >= 5L &&
      length(unique(one$donor_id)) >= 5L
  )

  eligibility_rows[[length(eligibility_rows) + 1L]] <-
    data.frame(
      broad_cell_type = cell_type,
      adequate_pseudobulk_n = nrow(one),
      adequate_donor_n = length(unique(one$donor_id)),
      fertile_LH3_n = as.integer(day_counts["3"]),
      fertile_LH5_n = as.integer(day_counts["5"]),
      fertile_LH7_n = as.integer(day_counts["7"]),
      fertile_LH9_n = as.integer(day_counts["9"]),
      fertile_LH11_n = as.integer(day_counts["11"]),
      RIF_LH7_n = rif_lh7_n,
      localization_eligible = localization_eligible,
      fertile_temporal_model_eligible = temporal_eligible,
      LH7_RIF_vs_fertile_model_eligible = lh7_eligible,
      stringsAsFactors = FALSE
    )
}

model_eligibility <- do.call(
  rbind,
  eligibility_rows
)

write_csv_utf8(
  model_eligibility,
  file.path(table_root, "04_model_eligibility.csv")
)


# ----------------------------------------------------------------------
# 8. Locked VDR program definitions and descriptive scores
# ----------------------------------------------------------------------

program_definitions <- list(
  RNA_ALL_436 =
    candidate_table$feature_id,

  OCCUPANCY_SUPPORTED_284 =
    candidate_table$feature_id[
      candidate_table$occupancy_supported
    ],

  SHARED_OCCUPANCY_46 =
    candidate_table$feature_id[
      candidate_table$shared_occupancy_supported
    ],

  PROMOTER_OCCUPANCY_84 =
    candidate_table$feature_id[
      candidate_table$any_promoter_overlap_support
    ]
)

all_candidate_genes <- unique(
  candidate_table$feature_id
)

present_candidate_genes <- intersect(
  all_candidate_genes,
  rownames(pb_counts)
)

candidate_counts <- as.matrix(
  pb_counts[
    present_candidate_genes,
    ,
    drop = FALSE
  ]
)

safe_library_size <- library_size
safe_library_size[safe_library_size <= 0] <- NA_real_

candidate_log_cpm <- log2(
  sweep(
    candidate_counts,
    MARGIN = 2,
    STATS = safe_library_size / 1e6,
    FUN = "/"
  ) + 1
)

candidate_z <- row_zscore_dense(
  candidate_log_cpm
)

score_rows <- list()
coverage_rows <- list()

for (program_name in names(program_definitions)) {
  frozen_genes <- unique(
    program_definitions[[program_name]]
  )

  measurable_genes <- intersect(
    frozen_genes,
    rownames(candidate_z)
  )

  candidate_subset <- candidate_table[
    match(
      measurable_genes,
      candidate_table$feature_id
    ),
    ,
    drop = FALSE
  ]

  weights <- sign(
    candidate_subset$interaction_logFC
  )

  valid <- is.finite(weights) & weights != 0

  measurable_genes <- measurable_genes[valid]
  weights <- weights[valid]

  if (length(measurable_genes) < 10L) {
    stop(
      "Program ",
      program_name,
      " has fewer than 10 measurable genes."
    )
  }

  signed_matrix <- sweep(
    candidate_z[
      measurable_genes,
      ,
      drop = FALSE
    ],
    MARGIN = 1,
    STATS = weights,
    FUN = "*"
  )

  score_value <- colMeans(
    signed_matrix,
    na.rm = TRUE
  )

  score_rows[[program_name]] <- data.frame(
    pb_id = colnames(pb_counts),
    program_name = program_name,
    score = as.numeric(score_value),
    frozen_gene_n = length(frozen_genes),
    measurable_gene_n = length(measurable_genes),
    coverage_fraction =
      length(measurable_genes) / length(frozen_genes),
    score_definition =
      "MEAN_SIGNED_GENE_ZSCORE_FROM_LOG2_CPM_PLUS_1",
    interpretation =
      "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION",
    stringsAsFactors = FALSE
  )

  coverage_rows[[program_name]] <- data.frame(
    program_name = program_name,
    frozen_gene_n = length(frozen_genes),
    measurable_gene_n = length(measurable_genes),
    coverage_fraction =
      length(measurable_genes) / length(frozen_genes),
    missing_genes = collapse_unique(
      setdiff(frozen_genes, measurable_genes)
    ),
    stringsAsFactors = FALSE
  )
}

program_scores <- do.call(
  rbind,
  score_rows
)

program_coverage <- do.call(
  rbind,
  coverage_rows
)

program_scores <- merge(
  program_scores,
  pb_metadata,
  by = "pb_id",
  all.x = TRUE,
  sort = FALSE
)

if (
  anyNA(program_scores$broad_cell_type) ||
    nrow(program_scores) !=
      nrow(pb_metadata) * length(program_definitions)
) {
  stop("Program-score metadata merge failed.")
}

program_scores$score_eligible <- (
  program_scores$cell_n >= 20L
)

write_csv_utf8(
  program_scores,
  file.path(table_root, "05_VDR_program_scores.csv")
)

write_csv_utf8(
  program_coverage,
  file.path(table_root, "06_VDR_program_coverage.csv")
)

localization_summary <- do.call(
  rbind,
  lapply(
    split(
      program_scores[
        program_scores$score_eligible,
        ,
        drop = FALSE
      ],
      list(
        program_scores$program_name[
          program_scores$score_eligible
        ],
        program_scores$broad_cell_type[
          program_scores$score_eligible
        ]
      ),
      drop = TRUE
    ),
    function(one) {
      data.frame(
        program_name = one$program_name[1],
        broad_cell_type = one$broad_cell_type[1],
        pseudobulk_n = nrow(one),
        donor_n = length(unique(one$donor_id)),
        mean_score = mean(one$score, na.rm = TRUE),
        median_score = stats::median(one$score, na.rm = TRUE),
        score_sd = stats::sd(one$score, na.rm = TRUE),
        score_IQR = stats::IQR(one$score, na.rm = TRUE),
        min_score = min(one$score, na.rm = TRUE),
        max_score = max(one$score, na.rm = TRUE),
        median_cell_n = stats::median(one$cell_n),
        stringsAsFactors = FALSE
      )
    }
  )
)

localization_summary <- localization_summary[
  order(
    localization_summary$program_name,
    -localization_summary$median_score
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  localization_summary,
  file.path(table_root, "07_VDR_localization_summary.csv")
)

# VDR gene expression itself, when present.
if ("VDR" %in% rownames(pb_counts)) {
  VDR_counts <- as.numeric(
    pb_counts["VDR", , drop = TRUE]
  )

  VDR_log_cpm <- log2(
    VDR_counts / (safe_library_size / 1e6) + 1
  )

  VDR_expression <- pb_metadata
  VDR_expression$VDR_raw_count <- VDR_counts
  VDR_expression$VDR_log2_CPM_plus_1 <- VDR_log_cpm

  write_csv_utf8(
    VDR_expression,
    file.path(table_root, "08_VDR_expression.csv")
  )
} else {
  VDR_expression <- data.frame()
  write_text(
    "VDR was not present in the pseudobulk count matrix row names.",
    file.path(run_qc_root, "VDR_gene_absent_note.txt")
  )
}


# ----------------------------------------------------------------------
# 9. QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank()
  )

plot_score_data <- program_scores[
  program_scores$score_eligible,
  ,
  drop = FALSE
]

plot_score_data$broad_cell_type <- factor(
  plot_score_data$broad_cell_type,
  levels = broad_types
)

plot_scores <- ggplot2::ggplot(
  plot_score_data,
  ggplot2::aes(
    x = broad_cell_type,
    y = score
  )
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(
      shape = clinical_group
    ),
    width = 0.15,
    height = 0,
    size = 1.4,
    alpha = 0.65
  ) +
  ggplot2::facet_wrap(
    ~ program_name,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::labs(
    title = "VDR-program localization across broad cell types",
    subtitle = paste(
      "Signed standardized pseudobulk scores;",
      "descriptive localization only"
    ),
    x = "Broad cell type",
    y = "Signed VDR-program score",
    shape = "Clinical group"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_pair(
  plot_scores,
  "QC1_VDR_program_scores",
  width = 12,
  height = 9
)

heatmap_data <- localization_summary
heatmap_data$broad_cell_type <- factor(
  heatmap_data$broad_cell_type,
  levels = rev(broad_types)
)

plot_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = program_name,
    y = broad_cell_type,
    fill = median_score
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.2f", median_score)
    ),
    size = 3
  ) +
  ggplot2::scale_fill_gradient2(
    midpoint = 0
  ) +
  ggplot2::labs(
    title = "Median VDR-program localization score",
    subtitle = "Positive values align with the locked VDR interaction direction",
    x = "VDR program",
    y = "Broad cell type",
    fill = "Median\nscore"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_pair(
  plot_heatmap,
  "QC2_VDR_localization_heatmap",
  width = 9,
  height = 6
)

plot_pb_cells <- ggplot2::ggplot(
  pb_metadata,
  ggplot2::aes(
    x = broad_cell_type,
    y = cell_n,
    shape = clinical_group
  )
) +
  ggplot2::geom_hline(
    yintercept = c(20, 50),
    linetype = c("solid", "dashed")
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    height = 0,
    size = 1.8,
    alpha = 0.7
  ) +
  ggplot2::scale_y_log10(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Pseudobulk cell-count QC",
    subtitle = "Solid line: 20-cell minimum; dashed line: 50-cell preferred threshold",
    x = "Broad cell type",
    y = "Aggregated cell count (log10 scale)",
    shape = "Clinical group"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_pair(
  plot_pb_cells,
  "QC3_pseudobulk_cell_counts",
  width = 10,
  height = 6
)

plot_eligibility_data <- reshape(
  model_eligibility[
    ,
    c(
      "broad_cell_type",
      "localization_eligible",
      "fertile_temporal_model_eligible",
      "LH7_RIF_vs_fertile_model_eligible"
    )
  ],
  varying = c(
    "localization_eligible",
    "fertile_temporal_model_eligible",
    "LH7_RIF_vs_fertile_model_eligible"
  ),
  v.names = "eligible",
  timevar = "analysis",
  times = c(
    "Localization",
    "Fertile temporal",
    "LH7 RIF vs fertile"
  ),
  direction = "long"
)

plot_eligibility <- ggplot2::ggplot(
  plot_eligibility_data,
  ggplot2::aes(
    x = analysis,
    y = broad_cell_type,
    fill = eligible
  )
) +
  ggplot2::geom_tile() +
  ggplot2::labs(
    title = "Broad-cell-type model eligibility",
    x = "Planned analysis",
    y = "Broad cell type",
    fill = "Eligible"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1
    )
  )

save_plot_pair(
  plot_eligibility,
  "QC4_model_eligibility",
  width = 8,
  height = 6
)

if (nrow(VDR_expression)) {
  plot_VDR <- ggplot2::ggplot(
    VDR_expression[
      VDR_expression$cell_n >= 20L,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      x = broad_cell_type,
      y = VDR_log2_CPM_plus_1,
      shape = clinical_group
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      height = 0,
      size = 1.5,
      alpha = 0.65
    ) +
    ggplot2::labs(
      title = "VDR pseudobulk expression across broad cell types",
      subtitle = "Descriptive expression localization",
      x = "Broad cell type",
      y = "VDR log2(CPM + 1)",
      shape = "Clinical group"
    ) +
    theme_qc +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )

  save_plot_pair(
    plot_VDR,
    "QC5_VDR_expression",
    width = 10,
    height = 6
  )
}


# ----------------------------------------------------------------------
# 10. Gate, Card and text drafts
# ----------------------------------------------------------------------

technical_checks <- data.frame(
  check = c(
    "annotation_gate_authorizes_broad_analysis",
    "checkpoint_complete",
    "counts_layer_n_is_28",
    "all_counts_layers_aggregated",
    "all_236636_cells_aggregated_once",
    "raw_count_total_preserved",
    "sample_n_is_28",
    "donor_n_is_28",
    "broad_type_n_is_8",
    "four_VDR_programs_scored",
    "all_program_coverage_ge_0_90",
    "at_least_one_type_eligible_for_temporal_model",
    "at_least_one_type_eligible_for_LH7_model"
  ),
  pass = c(
    annotation_gate$gate_status[[1]] ==
      "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW",
    checkpoint_status$checkpoint_status[[1]] == "COMPLETE",
    length(count_layers) == 28L,
    nrow(layer_status) == 28L &&
      all(layer_status$status == "COMPLETE"),
    sum(pb_metadata$cell_n) == 236636L,
    isTRUE(all.equal(
      raw_count_total_from_layers,
      combined_raw_count_total,
      tolerance = 0
    )),
    length(unique(pb_metadata$sample_accession)) == 28L,
    length(unique(pb_metadata$donor_id)) == 28L,
    length(unique(pb_metadata$broad_cell_type)) == 8L,
    length(unique(program_scores$program_name)) == 4L,
    all(program_coverage$coverage_fraction >= 0.90),
    any(model_eligibility$fertile_temporal_model_eligible),
    any(model_eligibility$LH7_RIF_vs_fertile_model_eligible)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  technical_checks,
  file.path(run_qc_root, "04_gate_checks.csv")
)

technical_pass <- all(technical_checks$pass)

gate_status <- if (technical_pass) {
  "PASS_PSEUDOBULK_AND_VDR_LOCALIZATION_READY_FOR_DONOR_MODELS"
} else {
  "HOLD_PSEUDOBULK_OR_LOCALIZATION_QC_FAILED"
}

gate_table <- data.frame(
  dataset_id = dataset_id,
  input_version = input_version,
  annotation_version = annotation_version,
  pseudobulk_version = pseudobulk_version,
  gate_status = gate_status,
  gene_n = nrow(pb_counts),
  pseudobulk_n = ncol(pb_counts),
  sample_n = length(unique(pb_metadata$sample_accession)),
  donor_n = length(unique(pb_metadata$donor_id)),
  broad_type_n = length(unique(pb_metadata$broad_cell_type)),
  aggregated_cell_n = sum(pb_metadata$cell_n),
  low_cell_pseudobulk_n = sum(pb_metadata$cell_n < 20L),
  preferred_cell_threshold_fail_n =
    sum(pb_metadata$cell_n < 50L),
  program_n = length(unique(program_scores$program_name)),
  minimum_program_coverage =
    min(program_coverage$coverage_fraction),
  temporal_model_eligible_type_n =
    sum(model_eligibility$fertile_temporal_model_eligible),
  LH7_model_eligible_type_n =
    sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible),
  donor_level_modeling_authorized = technical_pass,
  cell_composition_inference_authorized = FALSE,
  subtype_modeling_authorized = FALSE,
  result_lock_authorized = FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "PBVDR_GATE.csv"
)

write_csv_utf8(gate_table, gate_path)

card_lines <- c(
  "# STEP04E1F Pseudobulk and VDR Localization Result Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Input version: ", input_version),
  paste0("- Annotation version: ", annotation_version),
  paste0("- Pseudobulk version: ", pseudobulk_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Aggregated cells: ",
    format(sum(pb_metadata$cell_n), big.mark = ",")
  ),
  paste0(
    "- Pseudobulks: ",
    ncol(pb_counts)
  ),
  paste0(
    "- Broad cell types: ",
    length(unique(pb_metadata$broad_cell_type))
  ),
  paste0(
    "- Minimum VDR-program coverage: ",
    round(
      100 * min(program_coverage$coverage_fraction),
      1
    ),
    "%"
  ),
  paste0(
    "- Temporal-model eligible broad types: ",
    sum(model_eligibility$fertile_temporal_model_eligible)
  ),
  paste0(
    "- LH7-model eligible broad types: ",
    sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible)
  ),
  paste0(
    "- Donor-level modelling authorized: ",
    gate_table$donor_level_modeling_authorized
  ),
  "- Cell-composition inference authorized: FALSE",
  "- Subtype modelling authorized: FALSE",
  "- Result lock authorized: FALSE"
)

card_path <- file.path(
  run_qc_root,
  "PBVDR_RESULT_CARD.md"
)

write_text(card_lines, card_path)

interpretation_lines <- c(
  "# STEP04E1F interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Pseudobulk aggregation",
  "",
  paste0(
    "- All ",
    format(sum(pb_metadata$cell_n), big.mark = ","),
    " cells were aggregated exactly once across 28 samples and eight broad cell types."
  ),
  "- Raw count totals were preserved during layer-wise aggregation.",
  "- Raw pseudobulk counts were saved for subsequent donor-level modelling.",
  "",
  "## VDR localization",
  "",
  "- Four locked VDR RNA programs were scored using signed standardized pseudobulk expression.",
  "- Scores provide relative cell-type localization of expression aligned with the locked VDR interaction direction.",
  "- Scores do not establish direct VDR binding, causal activation or clinical differential expression.",
  "",
  "## Inferential boundary",
  "",
  "- Only broad cell types passing donor/cell-count eligibility may enter temporal or LH+7 clinical models.",
  "- Cell-composition hypothesis testing remains unauthorized.",
  "- Candidate subtype labels remain excluded from inferential analysis.",
  "",
  "## Next step",
  "",
  "- Fit donor-level fertile temporal models for eligible broad cell types.",
  "- Fit LH+7 RIF-versus-fertile models for eligible broad cell types.",
  "- Control multiple testing across broad cell types and programs.",
  "- Preserve VDR localization and donor-level clinical effects as distinct evidence layers."
)

interpretation_path <- file.path(
  run_qc_root,
  "PBVDR_INTERPRETATION.md"
)

write_text(
  interpretation_lines,
  interpretation_path
)

methods_lines <- c(
  "# STEP04E1F Methods draft",
  "",
  paste(
    "Raw RNA counts were aggregated separately within each split counts",
    "layer using the accepted sample-level broad cell-type annotations."
  ),
  paste(
    "For each sample and broad cell type, a sparse cell-membership matrix",
    "was multiplied by the raw gene-by-cell count matrix to obtain donor-level",
    "pseudobulk counts."
  ),
  paste(
    "Aggregation was performed layer by layer, and each intermediate",
    "pseudobulk object was saved immediately to support reproducible restart."
  ),
  paste(
    "Pseudobulk profiles containing fewer than 20 cells were excluded from",
    "program-score summaries and were retained with explicit QC flags."
  ),
  paste(
    "For descriptive VDR-program localization, pseudobulk counts were",
    "converted to log2 counts per million plus one."
  ),
  paste(
    "Each program gene was standardized across all donor-by-cell-type",
    "pseudobulk profiles, multiplied by the sign of its locked VDR",
    "interaction log fold-change, and averaged within each prespecified",
    "program."
  ),
  paste(
    "The four programs comprised the complete 436-gene RNA interaction",
    "program, the 284-gene occupancy-supported subset, the 46-gene",
    "shared-occupancy subset and the 84-gene promoter-occupancy subset."
  ),
  paste(
    "These scores were used as descriptive localization measures and not",
    "as evidence of causal VDR activation or donor-level clinical effects."
  )
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Methods_draft.md"
  )
)

results_lines <- c(
  "# STEP04E1F Results shell",
  "",
  paste0(
    "Raw counts from ",
    format(sum(pb_metadata$cell_n), big.mark = ","),
    " cells were aggregated across ",
    length(unique(pb_metadata$sample_accession)),
    " samples and ",
    length(unique(pb_metadata$broad_cell_type)),
    " broad cell populations."
  ),
  paste0(
    "The resulting pseudobulk matrix contained ",
    format(nrow(pb_counts), big.mark = ","),
    " genes and ",
    ncol(pb_counts),
    " donor-by-cell-type profiles."
  ),
  paste0(
    "All four locked VDR programs were measurable, with a minimum gene",
    " coverage of ",
    round(
      100 * min(program_coverage$coverage_fraction),
      1
    ),
    "%."
  ),
  paste0(
    sum(model_eligibility$fertile_temporal_model_eligible),
    " broad cell types met the predefined fertile temporal modelling",
    " criteria, and ",
    sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible),
    " met the LH+7 RIF-versus-fertile criteria."
  ),
  paste(
    "Cell-type localization scores were retained as a descriptive evidence",
    "layer, whereas clinical and temporal comparisons were deferred to",
    "donor-level models."
  )
)

write_text(
  results_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Results_shell.md"
  )
)


# ----------------------------------------------------------------------
# 11. Manifests, cleanup and small upload ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_names <- c(
  "SeuratObject",
  "Matrix",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales",
  "future"
)

package_table <- data.frame(
  package = package_names,
  installed = vapply(
    package_names,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ),
  version = vapply(
    package_names,
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(run_qc_root, "PACKAGE_VERSIONS.csv")
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- unique(c(
  list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Results_shell.md"
  )
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(artifact_files),
  size_bytes = file.info(artifact_files)$size,
  md5 = unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(run_qc_root, "ARTIFACT_MANIFEST.csv")
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  annotation_version = annotation_version,
  pseudobulk_version = pseudobulk_version,
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_minutes = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "mins"
      )
    ),
    2
  ),
  gate_status = gate_status,
  gene_n = nrow(pb_counts),
  pseudobulk_n = ncol(pb_counts),
  sample_n = length(unique(pb_metadata$sample_accession)),
  donor_n = length(unique(pb_metadata$donor_id)),
  broad_type_n = length(unique(pb_metadata$broad_cell_type)),
  aggregated_cell_n = sum(pb_metadata$cell_n),
  counts_layer_n = length(count_layers),
  program_n = length(unique(program_scores$program_name)),
  donor_level_modeling_authorized =
    gate_table$donor_level_modeling_authorized,
  composition_inference_authorized = FALSE,
  subtype_modeling_authorized = FALSE,
  result_lock_authorized = FALSE,
  output_root = normalize_slashes(output_root),
  qc_root = normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "pseudobulk_n=",
  ncol(pb_counts),
  "temporal_eligible_types=",
  sum(model_eligibility$fertile_temporal_model_eligible),
  "LH7_eligible_types=",
  sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible)
)

# Large count objects are not included in the upload ZIP.
small_files <- unique(c(
  list.files(
    table_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    source_data_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  file.path(object_root, "pb_status.csv"),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1F_PBVDR_Results_shell.md"
  ),
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <=
      200 * 1024^2
]

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = small_files,
        recurse = FALSE,
        include_directories = FALSE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

rm(
  candidate_counts,
  candidate_log_cpm,
  candidate_z,
  pb_counts
)
gc()

cat("\n============================================================\n")
cat("STEP04E1F GSE250130 VDR localization/pseudobulk completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Aggregated cells: ",
  format(sum(pb_metadata$cell_n), big.mark = ","),
  "\n",
  sep = ""
)
cat(
  "Pseudobulks: ",
  nrow(pb_metadata),
  "\n",
  sep = ""
)
cat(
  "Broad cell types: ",
  length(unique(pb_metadata$broad_cell_type)),
  "\n",
  sep = ""
)
cat(
  "Temporal-model eligible types: ",
  sum(model_eligibility$fertile_temporal_model_eligible),
  "\n",
  sep = ""
)
cat(
  "LH7-model eligible types: ",
  sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible),
  "\n",
  sep = ""
)
cat("Output root: ", output_root, "\n", sep = "")
cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
  cat(
    "Upload this ZIP:\n",
    handoff_zip,
    "\n"
  )
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe output and QC folders are complete.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 7: 61_step04E1G_GSE250130_donor_VDR_models_v1_1_unbalanced_fixed
# ==============================================================================

# ======================================================================
# 61_step04E1G_GSE250130_donor_VDR_models_v1_1_unbalanced_fixed.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Fit donor-level models for the four locked VDR RNA programs.
# - Fit parallel donor-level models for VDR expression.
# - Analyse two predefined questions:
#
#   A. Fertile temporal trajectory:
#      LH+3, LH+5, LH+7, LH+9 and LH+11
#
#   B. LH+7 RIF versus fertile comparison:
#      primary: RIF versus all fertile LH+7 donors
#      sensitivity: RIF versus fertile LH7_STANDARD
#      sensitivity: RIF versus fertile LH7_PGT
#
# PRIMARY INFERENCE
# - Temporal omnibus: permutation one-way ANOVA.
# - LH+7 group difference: exact permutation test for mean difference.
# - Classical ANOVA/Welch tests are retained as secondary diagnostics.
# - Benjamini-Hochberg correction is applied within predefined model families.
#
# SCIENTIFIC SCOPE
# - One pseudobulk profile represents one biological donor x broad cell type.
# - No cell is treated as an independent replicate.
# - Program scores represent relative alignment with the locked VDR
#   interaction direction; they do not prove pathway activation or causality.
# - Cell-composition and subtype inference remain unauthorized.
#
# INPUT
# - PBVDR_v1_0 Gate:
#   PASS_PSEUDOBULK_AND_VDR_LOCALIZATION_READY_FOR_DONOR_MODELS
#
# OUTPUT STATUS
# - DONOR_LEVEL_VDR_MODEL_CANDIDATES
# - No output is written to 06_locked_results.
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)

set.seed(20260805)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <- "61_step04E1G_GSE250130_donor_VDR_models_v1_1_unbalanced_fixed.R"
script_version <- "STEP04E1G_GSE250130_DONOR_VDR_MODELS_v1.1_UNBALANCED_FIXED"
input_version <- "SCI_v1_1_memory_safe"
pseudobulk_version <- "PBVDR_v1_0"
model_version <- "VDRMODELS_v1_1"

run_start <- Sys.time()
run_id <- paste0(
  "STEP04E1G_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(
  path,
  max_chars = 245L
) {
  normalized <- normalize_slashes(path)
  path_n <- nchar(normalized, type = "chars")

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Output path is too long for safe Windows/R writing (",
      path_n,
      " characters; safe limit ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) "" else paste(x, collapse = sep)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x < 0.001,
      format(x, scientific = TRUE, digits = 2),
      sprintf("%.3f", x)
    )
  )
}

sign_label <- function(effect) {
  ifelse(
    is.na(effect),
    "NA",
    ifelse(
      effect > 0,
      "POSITIVE",
      ifelse(effect < 0, "NEGATIVE", "ZERO")
    )
  )
}

save_plot_pair <- function(plot_object, stem, width, height) {
  png_path <- file.path(figure_root, paste0(stem, ".png"))
  pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))

  assert_windows_path_length(png_path)
  assert_windows_path_length(pdf_path)

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

oneway_F <- function(y, group) {
  keep <- is.finite(y) & !is.na(group)
  y <- y[keep]
  group <- droplevels(factor(group[keep]))

  n <- length(y)
  k <- nlevels(group)

  if (n <= k || k < 2L) return(NA_real_)

  grand_mean <- mean(y)
  group_means <- tapply(y, group, mean)
  group_n <- table(group)

  ss_between <- sum(
    as.numeric(group_n) *
      (as.numeric(group_means) - grand_mean)^2
  )

  fitted_group_mean <- group_means[as.character(group)]
  ss_within <- sum((y - fitted_group_mean)^2)

  df_between <- k - 1L
  df_within <- n - k

  if (df_within <= 0L || ss_within <= 0) {
    return(NA_real_)
  }

  (ss_between / df_between) /
    (ss_within / df_within)
}

permutation_oneway <- function(
  y,
  group,
  n_perm = 5000L,
  seed = 20260805L
) {
  keep <- is.finite(y) & !is.na(group)
  y <- as.numeric(y[keep])
  group <- droplevels(factor(group[keep]))

  observed_F <- oneway_F(y, group)

  if (!is.finite(observed_F)) {
    return(
      list(
        F = observed_F,
        p = NA_real_,
        permutation_n = 0L
      )
    )
  }

  set.seed(seed)

  perm_F <- replicate(
    n_perm,
    oneway_F(
      sample(y, length(y), replace = FALSE),
      group
    )
  )

  perm_F <- perm_F[is.finite(perm_F)]

  p_value <- (
    1 + sum(perm_F >= observed_F)
  ) / (
    1 + length(perm_F)
  )

  list(
    F = observed_F,
    p = p_value,
    permutation_n = length(perm_F)
  )
}

exact_two_group_test <- function(
  y,
  group,
  reference_level,
  comparison_level,
  maximum_exact_combinations = 50000L,
  Monte_Carlo_n = 20000L,
  seed = 20260805L
) {
  keep <- (
    is.finite(y) &
      !is.na(group) &
      group %in% c(
        reference_level,
        comparison_level
      )
  )

  y <- as.numeric(y[keep])
  group <- as.character(group[keep])

  reference_values <- y[group == reference_level]
  comparison_values <- y[group == comparison_level]

  n_reference <- length(reference_values)
  n_comparison <- length(comparison_values)
  n_total <- length(y)

  if (
    n_reference < 2L ||
      n_comparison < 2L
  ) {
    return(
      list(
        estimate = NA_real_,
        p = NA_real_,
        test_mode = "INSUFFICIENT_REPLICATION",
        permutation_n = 0L
      )
    )
  }

  observed <- (
    mean(comparison_values) -
      mean(reference_values)
  )

  total_sum <- sum(y)
  combination_n <- choose(n_total, n_comparison)

  if (
    is.finite(combination_n) &&
      combination_n <= maximum_exact_combinations
  ) {
    combinations <- utils::combn(
      n_total,
      n_comparison
    )

    comparison_sums <- colSums(
      matrix(
        y[combinations],
        nrow = n_comparison
      )
    )

    permuted_differences <- (
      comparison_sums / n_comparison
    ) - (
      (total_sum - comparison_sums) /
        n_reference
    )

    p_value <- mean(
      abs(permuted_differences) >=
        abs(observed) - 1e-12
    )

    test_mode <- "EXACT"
    permutation_n <- length(permuted_differences)
  } else {
    set.seed(seed)

    permuted_differences <- replicate(
      Monte_Carlo_n,
      {
        comparison_index <- sample.int(
          n_total,
          n_comparison,
          replace = FALSE
        )

        comparison_sum <- sum(y[comparison_index])

        (
          comparison_sum / n_comparison
        ) - (
          (total_sum - comparison_sum) /
            n_reference
        )
      }
    )

    p_value <- (
      1 +
        sum(
          abs(permuted_differences) >=
            abs(observed) - 1e-12
        )
    ) / (
      1 + length(permuted_differences)
    )

    test_mode <- "MONTE_CARLO"
    permutation_n <- length(permuted_differences)
  }

  list(
    estimate = observed,
    p = p_value,
    test_mode = test_mode,
    permutation_n = permutation_n
  )
}

welch_summary <- function(
  y,
  group,
  reference_level,
  comparison_level
) {
  keep <- (
    is.finite(y) &
      !is.na(group) &
      group %in% c(
        reference_level,
        comparison_level
      )
  )

  y <- as.numeric(y[keep])
  group <- as.character(group[keep])

  reference_values <- y[group == reference_level]
  comparison_values <- y[group == comparison_level]

  if (
    length(reference_values) < 2L ||
      length(comparison_values) < 2L
  ) {
    return(
      list(
        estimate = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p = NA_real_
      )
    )
  }

  test <- stats::t.test(
    comparison_values,
    reference_values,
    var.equal = FALSE,
    conf.level = 0.95
  )

  list(
    estimate = (
      mean(comparison_values) -
        mean(reference_values)
    ),
    ci_low = unname(test$conf.int[1]),
    ci_high = unname(test$conf.int[2]),
    p = unname(test$p.value)
  )
}

hedges_g <- function(
  y,
  group,
  reference_level,
  comparison_level
) {
  keep <- (
    is.finite(y) &
      !is.na(group) &
      group %in% c(
        reference_level,
        comparison_level
      )
  )

  y <- as.numeric(y[keep])
  group <- as.character(group[keep])

  x_ref <- y[group == reference_level]
  x_cmp <- y[group == comparison_level]

  n_ref <- length(x_ref)
  n_cmp <- length(x_cmp)

  if (
    n_ref < 2L ||
      n_cmp < 2L
  ) {
    return(NA_real_)
  }

  pooled_variance <- (
    (n_ref - 1) * stats::var(x_ref) +
      (n_cmp - 1) * stats::var(x_cmp)
  ) / (
    n_ref + n_cmp - 2
  )

  if (
    !is.finite(pooled_variance) ||
      pooled_variance <= 0
  ) {
    return(NA_real_)
  }

  cohen_d <- (
    mean(x_cmp) - mean(x_ref)
  ) / sqrt(pooled_variance)

  correction <- (
    1 -
      3 / (
        4 * (n_ref + n_cmp) - 9
      )
  )

  correction * cohen_d
}

leave_one_out_difference <- function(
  y,
  group,
  reference_level,
  comparison_level
) {
  keep <- (
    is.finite(y) &
      !is.na(group) &
      group %in% c(
        reference_level,
        comparison_level
      )
  )

  y <- as.numeric(y[keep])
  group <- as.character(group[keep])

  if (
    sum(group == reference_level) < 3L ||
      sum(group == comparison_level) < 3L
  ) {
    return(
      c(
        loo_min = NA_real_,
        loo_max = NA_real_,
        loo_sign_consistency = NA_real_
      )
    )
  }

  full_effect <- (
    mean(y[group == comparison_level]) -
      mean(y[group == reference_level])
  )

  loo_effects <- vapply(
    seq_along(y),
    function(i) {
      y_i <- y[-i]
      group_i <- group[-i]

      (
        mean(y_i[group_i == comparison_level]) -
          mean(y_i[group_i == reference_level])
      )
    },
    numeric(1)
  )

  full_sign <- sign(full_effect)

  c(
    loo_min = min(loo_effects),
    loo_max = max(loo_effects),
    loo_sign_consistency = mean(
      sign(loo_effects) == full_sign
    )
  )
}

expected_temporal_counts <- function(
  eligibility_table,
  cell_type
) {
  one <- eligibility_table[
    eligibility_table$broad_cell_type == cell_type,
    ,
    drop = FALSE
  ]

  if (nrow(one) != 1L) {
    stop(
      "Eligibility table does not contain exactly one row for: ",
      cell_type
    )
  }

  counts <- c(
    `3` = safe_numeric(one$fertile_LH3_n),
    `5` = safe_numeric(one$fertile_LH5_n),
    `7` = safe_numeric(one$fertile_LH7_n),
    `9` = safe_numeric(one$fertile_LH9_n),
    `11` = safe_numeric(one$fertile_LH11_n)
  )

  if (
    any(!is.finite(counts)) ||
      any(counts < 2L)
  ) {
    stop(
      "Temporal eligibility counts are invalid for ",
      cell_type,
      ": ",
      paste(
        names(counts),
        counts,
        sep = "=",
        collapse = ";"
      )
    )
  }

  as.integer(counts)
}


group_summary <- function(
  data,
  value_column,
  group_columns
) {
  dt <- data.table::as.data.table(data)

  dt[
    ,
    {
      values <- safe_numeric(
        get(value_column)
      )

      values <- values[is.finite(values)]

      n <- length(values)
      mean_value <- if (n) mean(values) else NA_real_
      sd_value <- if (n >= 2L) stats::sd(values) else NA_real_
      se_value <- if (n >= 2L) sd_value / sqrt(n) else NA_real_
      critical <- if (n >= 2L) {
        stats::qt(0.975, df = n - 1L)
      } else {
        NA_real_
      }

      list(
        n = n,
        mean = mean_value,
        sd = sd_value,
        se = se_value,
        ci_low = mean_value - critical * se_value,
        ci_high = mean_value + critical * se_value,
        median = if (n) stats::median(values) else NA_real_,
        IQR = if (n) stats::IQR(values) else NA_real_
      )
    },
    by = group_columns
  ]
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")

  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages and paths
# ----------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(step04_root, "05_execution")
results_root <- file.path(step04_root, "06_results")
qc_base_root <- file.path(step04_root, "08_qc")
parameter_root <- file.path(step04_root, "03_parameters")
environment_root <- file.path(step04_root, "04_environment")
manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

pb_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "pb",
  dataset_id,
  pseudobulk_version
)

pb_table_root <- file.path(pb_root, "tables")
pb_object_root <- file.path(pb_root, "objects")

program_score_path <- file.path(
  pb_table_root,
  "05_VDR_program_scores.csv"
)

VDR_expression_path <- file.path(
  pb_table_root,
  "08_VDR_expression.csv"
)

eligibility_path <- file.path(
  pb_table_root,
  "04_model_eligibility.csv"
)

pb_metadata_path <- file.path(
  pb_table_root,
  "02_pb_metadata.csv"
)

pb_status_path <- file.path(
  pb_object_root,
  "pb_status.csv"
)

output_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "models",
  dataset_id,
  model_version
)

table_root <- file.path(output_root, "tables")
source_data_root <- file.path(output_root, "source_data")

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1G_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_modelQC"
  )
)

figure_root <- file.path(run_qc_root, "figures")

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  output_root,
  table_root,
  source_data_root,
  run_qc_root,
  figure_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)

planned_paths <- data.frame(
  object = c(
    "output_root",
    "table_root",
    "source_data_root",
    "run_qc_root",
    "figure_root",
    "log_path",
    "run_manifest_path",
    "parameter_path",
    "session_path",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    output_root,
    table_root,
    source_data_root,
    run_qc_root,
    figure_root,
    log_path,
    run_manifest_path,
    parameter_path,
    session_path,
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_paths$character_n <- nchar(
  planned_paths$path,
  type = "chars"
)

if (
  .Platform$OS.type == "windows" &&
    any(planned_paths$character_n > 245L)
) {
  bad <- planned_paths[
    planned_paths$character_n > 245L,
    ,
    drop = FALSE
  ]

  stop(
    "One or more planned paths are too long:\n",
    paste(
      paste0(
        bad$character_n,
        " chars: ",
        bad$path
      ),
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------------------
# 3. Validate upstream Gate and inputs
# ----------------------------------------------------------------------

required_inputs <- c(
  program_score_path,
  VDR_expression_path,
  eligibility_path,
  pb_metadata_path,
  pb_status_path
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

pb_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

pb_qc_dirs <- pb_qc_dirs[
  grepl(
    "^STEP04E1F_[0-9]{8}_[0-9]{6}_pbQC$",
    basename(pb_qc_dirs),
    perl = TRUE
  )
]

if (!length(pb_qc_dirs)) {
  stop("No completed STEP04E1F pseudobulk-QC folder found.")
}

latest_pb_qc <- pb_qc_dirs[
  which.max(file.info(pb_qc_dirs)$mtime)
]

pb_gate_path <- file.path(
  latest_pb_qc,
  "PBVDR_GATE.csv"
)

pb_gate <- safe_read_csv(pb_gate_path)
pb_status <- safe_read_csv(pb_status_path)

if (
  is.null(pb_gate) ||
    nrow(pb_gate) != 1L ||
    pb_gate$gate_status[[1]] !=
      "PASS_PSEUDOBULK_AND_VDR_LOCALIZATION_READY_FOR_DONOR_MODELS" ||
    !isTRUE(
      pb_gate$donor_level_modeling_authorized[[1]]
    )
) {
  stop("STEP04E1F Gate does not authorize donor-level models.")
}

if (
  is.null(pb_status) ||
    nrow(pb_status) != 1L ||
    pb_status$status[[1]] != "COMPLETE" ||
    pb_status$pseudobulk_n[[1]] != 222L ||
    pb_status$aggregated_cell_n[[1]] != 236636L
) {
  stop("PBVDR pseudobulk status failed validation.")
}

scores <- data.table::fread(
  program_score_path,
  data.table = FALSE,
  showProgress = FALSE
)

VDR_expression <- data.table::fread(
  VDR_expression_path,
  data.table = FALSE,
  showProgress = FALSE
)

eligibility <- safe_read_csv(eligibility_path)
pb_metadata <- safe_read_csv(pb_metadata_path)

required_score_columns <- c(
  "pb_id",
  "program_name",
  "score",
  "sample_accession",
  "donor_id",
  "broad_cell_type",
  "clinical_group",
  "LH_day",
  "temporal_fertile_eligible",
  "LH7_RIF_comparison_eligible",
  "LH7_reference_group",
  "LH7_fertile_subgroup",
  "minimum_cell_threshold_pass",
  "score_eligible"
)

missing_score_columns <- setdiff(
  required_score_columns,
  colnames(scores)
)

if (length(missing_score_columns)) {
  stop(
    "Program-score table is missing columns: ",
    paste(missing_score_columns, collapse = ", ")
  )
}

required_VDR_columns <- c(
  "pb_id",
  "sample_accession",
  "donor_id",
  "broad_cell_type",
  "clinical_group",
  "LH_day",
  "temporal_fertile_eligible",
  "LH7_RIF_comparison_eligible",
  "LH7_reference_group",
  "LH7_fertile_subgroup",
  "minimum_cell_threshold_pass",
  "VDR_log2_CPM_plus_1"
)

missing_VDR_columns <- setdiff(
  required_VDR_columns,
  colnames(VDR_expression)
)

if (length(missing_VDR_columns)) {
  stop(
    "VDR-expression table is missing columns: ",
    paste(missing_VDR_columns, collapse = ", ")
  )
}

if (
  nrow(scores) != 888L ||
    length(unique(scores$pb_id)) != 222L ||
    length(unique(scores$program_name)) != 4L ||
    nrow(VDR_expression) != 222L ||
    nrow(eligibility) != 8L ||
    nrow(pb_metadata) != 222L
) {
  stop("PBVDR input identity checks failed.")
}

temporal_types <- as.character(
  eligibility$broad_cell_type[
    eligibility$fertile_temporal_model_eligible
  ]
)

LH7_types <- as.character(
  eligibility$broad_cell_type[
    eligibility$LH7_RIF_vs_fertile_model_eligible
  ]
)

program_names <- sort(
  unique(as.character(scores$program_name))
)

expected_programs <- c(
  "OCCUPANCY_SUPPORTED_284",
  "PROMOTER_OCCUPANCY_84",
  "RNA_ALL_436",
  "SHARED_OCCUPANCY_46"
)

if (!setequal(program_names, expected_programs)) {
  stop("The four locked VDR programs were not recovered exactly.")
}

if (
  length(temporal_types) != 5L ||
    length(LH7_types) != 6L
) {
  stop("Model-eligible cell-type counts do not match the Gate.")
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  model_version = model_version,
  supersedes_failed_model_version =
    "VDRMODELS_v1_0_HARD_CODED_SAMPLE_COUNT_QC_FAILURE",
  dataset_id = dataset_id,
  temporal_question =
    "FERTILE_LH3_LH5_LH7_LH9_LH11",
  temporal_replication_policy =
    "USE_ACTUAL_GATE_ELIGIBLE_COUNTS_MINIMUM_2_DONORS_PER_TIMEPOINT",
  temporal_design_balance =
    "UNBALANCED_ALLOWED_AND_EXPECTED_FOR_CILIATED_EPITHELIAL_LH11_N2",
  supersedes_failed_model_version =
    "VDRMODELS_v1_0_HARD_CODED_SAMPLE_COUNT_QC_FAILURE",
  temporal_primary_test =
    "5000_LABEL_PERMUTATION_ONE_WAY_ANOVA",
  temporal_secondary_test =
    "CLASSICAL_ONE_WAY_ANOVA",
  temporal_planned_pairwise = c(
    "LH5_vs_LH3",
    "LH7_vs_LH5",
    "LH9_vs_LH7",
    "LH11_vs_LH9",
    "LH7_vs_LH3",
    "LH11_vs_LH3",
    "LH11_vs_LH7"
  ),
  LH7_primary_question =
    "RIF_VS_ALL_FERTILE_LH7",
  LH7_primary_test =
    "EXACT_PERMUTATION_MEAN_DIFFERENCE",
  LH7_secondary_test =
    "WELCH_T_TEST",
  LH7_sensitivity_questions = c(
    "RIF_VS_FERTILE_LH7_STANDARD",
    "RIF_VS_FERTILE_LH7_PGT"
  ),
  multiple_testing =
    "BENJAMINI_HOCHBERG_WITHIN_PREDEFINED_MODEL_FAMILY",
  minimum_cell_policy =
    "minimum_cell_threshold_pass_TRUE",
  program_score_role =
    "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION",
  composition_inference_authorized = FALSE,
  subtype_inference_authorized = FALSE,
  random_seed = 20260805L,
  result_lock_authorized = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_csv_utf8(
  planned_paths,
  file.path(run_qc_root, "00_path_audit.csv")
)

log_event(
  "Starting",
  run_id,
  "temporal_types=",
  collapse_unique(temporal_types),
  "LH7_types=",
  collapse_unique(LH7_types)
)


# ----------------------------------------------------------------------
# 4. Input sample and model-family audit
# ----------------------------------------------------------------------

scores$LH_day <- safe_numeric(scores$LH_day)
VDR_expression$LH_day <- safe_numeric(VDR_expression$LH_day)

scores$score <- safe_numeric(scores$score)
VDR_expression$VDR_log2_CPM_plus_1 <-
  safe_numeric(VDR_expression$VDR_log2_CPM_plus_1)

scores$minimum_cell_threshold_pass <- as.logical(
  scores$minimum_cell_threshold_pass
)
scores$score_eligible <- as.logical(scores$score_eligible)

VDR_expression$minimum_cell_threshold_pass <- as.logical(
  VDR_expression$minimum_cell_threshold_pass
)

temporal_score_data <- scores[
  scores$clinical_group == "FERTILE" &
    scores$temporal_fertile_eligible &
    scores$minimum_cell_threshold_pass &
    scores$score_eligible &
    scores$broad_cell_type %in% temporal_types,
  ,
  drop = FALSE
]

LH7_score_data <- scores[
  scores$LH7_RIF_comparison_eligible &
    scores$minimum_cell_threshold_pass &
    scores$score_eligible &
    scores$broad_cell_type %in% LH7_types,
  ,
  drop = FALSE
]

temporal_VDR_data <- VDR_expression[
  VDR_expression$clinical_group == "FERTILE" &
    VDR_expression$temporal_fertile_eligible &
    VDR_expression$minimum_cell_threshold_pass &
    VDR_expression$broad_cell_type %in% temporal_types,
  ,
  drop = FALSE
]

LH7_VDR_data <- VDR_expression[
  VDR_expression$LH7_RIF_comparison_eligible &
    VDR_expression$minimum_cell_threshold_pass &
    VDR_expression$broad_cell_type %in% LH7_types,
  ,
  drop = FALSE
]

input_audit_rows <- list()

for (cell_type in temporal_types) {
  one <- temporal_VDR_data[
    temporal_VDR_data$broad_cell_type == cell_type,
    ,
    drop = FALSE
  ]

  observed <- table(factor(
    one$LH_day,
    levels = c(3, 5, 7, 9, 11)
  ))

  expected <- expected_temporal_counts(
    eligibility_table = eligibility,
    cell_type = cell_type
  )

  input_audit_rows[[length(input_audit_rows) + 1L]] <-
    data.frame(
      model_family = "FERTILE_TEMPORAL",
      broad_cell_type = cell_type,
      expected_groups = paste0(
        "LH",
        names(expected),
        "=",
        as.integer(expected),
        collapse = ";"
      ),
      observed_groups = paste0(
        "LH",
        names(observed),
        "=",
        as.integer(observed),
        collapse = ";"
      ),
      donor_n = length(unique(one$donor_id)),
      minimum_group_n = min(as.integer(observed)),
      unbalanced_design = length(
        unique(as.integer(observed))
      ) > 1L,
      identity_pass = identical(
        as.integer(observed),
        as.integer(expected)
      ) &&
        all(as.integer(observed) >= 2L),
      stringsAsFactors = FALSE
    )
}

for (cell_type in LH7_types) {
  one <- LH7_VDR_data[
    LH7_VDR_data$broad_cell_type == cell_type,
    ,
    drop = FALSE
  ]

  observed <- table(factor(
    one$clinical_group,
    levels = c("FERTILE", "RIF")
  ))

  input_audit_rows[[length(input_audit_rows) + 1L]] <-
    data.frame(
      model_family = "LH7_RIF_VS_FERTILE",
      broad_cell_type = cell_type,
      expected_groups = "FERTILE>=3;RIF>=3",
      observed_groups = paste0(
        names(observed),
        "=",
        as.integer(observed),
        collapse = ";"
      ),
      donor_n = length(unique(one$donor_id)),
      minimum_group_n = min(as.integer(observed)),
      unbalanced_design = length(
        unique(as.integer(observed))
      ) > 1L,
      identity_pass = all(as.integer(observed) >= 3L),
      stringsAsFactors = FALSE
    )
}

input_audit <- do.call(
  rbind,
  input_audit_rows
)

write_csv_utf8(
  input_audit,
  file.path(run_qc_root, "01_model_input_audit.csv")
)

if (!all(input_audit$identity_pass)) {
  stop("One or more donor-level model inputs failed sample-count QC.")
}


# ----------------------------------------------------------------------
# 5. Temporal program-score models
# ----------------------------------------------------------------------

temporal_program_rows <- list()
temporal_pairwise_rows <- list()
temporal_summary_rows <- list()

planned_pairs <- data.frame(
  comparison = c(
    "LH5_vs_LH3",
    "LH7_vs_LH5",
    "LH9_vs_LH7",
    "LH11_vs_LH9",
    "LH7_vs_LH3",
    "LH11_vs_LH3",
    "LH11_vs_LH7"
  ),
  reference_day = c(3, 5, 7, 9, 3, 3, 7),
  comparison_day = c(5, 7, 9, 11, 7, 11, 11),
  stringsAsFactors = FALSE
)

model_counter <- 0L

for (cell_type in temporal_types) {
  for (program_name in program_names) {
    model_counter <- model_counter + 1L

    one <- temporal_score_data[
      temporal_score_data$broad_cell_type == cell_type &
        temporal_score_data$program_name == program_name,
      ,
      drop = FALSE
    ]

    one$LH_factor <- factor(
      one$LH_day,
      levels = c(3, 5, 7, 9, 11)
    )

    expected_counts <- expected_temporal_counts(
      eligibility_table = eligibility,
      cell_type = cell_type
    )

    observed_counts <- as.integer(
      table(one$LH_factor)
    )

    if (
      nrow(one) != sum(expected_counts) ||
        !identical(
          observed_counts,
          as.integer(expected_counts)
        ) ||
        any(observed_counts < 2L)
    ) {
      stop(
        "Temporal score input mismatch for ",
        cell_type,
        " / ",
        program_name,
        ". Expected ",
        paste(
          names(expected_counts),
          expected_counts,
          sep = "=",
          collapse = ";"
        ),
        "; observed ",
        paste(
          levels(one$LH_factor),
          observed_counts,
          sep = "=",
          collapse = ";"
        )
      )
    }

    permutation <- permutation_oneway(
      y = one$score,
      group = one$LH_factor,
      n_perm = 5000L,
      seed = 20260805L + model_counter
    )

    fit <- stats::lm(
      score ~ LH_factor,
      data = one
    )

    anova_table <- stats::anova(fit)

    classical_p <- anova_table$`Pr(>F)`[1]
    classical_F <- anova_table$`F value`[1]

    group_means <- tapply(
      one$score,
      one$LH_factor,
      mean
    )

    total_ss <- sum(
      (one$score - mean(one$score))^2
    )

    fitted_ss <- sum(
      (
        fitted(fit) -
          mean(one$score)
      )^2
    )

    eta_squared <- if (
      is.finite(total_ss) &&
        total_ss > 0
    ) {
      fitted_ss / total_ss
    } else {
      NA_real_
    }

    temporal_program_rows[[
      length(temporal_program_rows) + 1L
    ]] <- data.frame(
      broad_cell_type = cell_type,
      program_name = program_name,
      donor_n = nrow(one),
      LH3_n = sum(one$LH_day == 3),
      LH5_n = sum(one$LH_day == 5),
      LH7_n = sum(one$LH_day == 7),
      LH9_n = sum(one$LH_day == 9),
      LH11_n = sum(one$LH_day == 11),
      LH3_mean = unname(group_means["3"]),
      LH5_mean = unname(group_means["5"]),
      LH7_mean = unname(group_means["7"]),
      LH9_mean = unname(group_means["9"]),
      LH11_mean = unname(group_means["11"]),
      permutation_F = permutation$F,
      permutation_p = permutation$p,
      permutation_n = permutation$permutation_n,
      classical_F = classical_F,
      classical_p = classical_p,
      eta_squared = eta_squared,
      stringsAsFactors = FALSE
    )

    summary_one <- group_summary(
      data = one,
      value_column = "score",
      group_columns = c(
        "broad_cell_type",
        "program_name",
        "LH_day"
      )
    )

    temporal_summary_rows[[
      length(temporal_summary_rows) + 1L
    ]] <- as.data.frame(summary_one)

    for (pair_index in seq_len(nrow(planned_pairs))) {
      pair <- planned_pairs[pair_index, ]

      pair_data <- one[
        one$LH_day %in% c(
          pair$reference_day,
          pair$comparison_day
        ),
        ,
        drop = FALSE
      ]

      pair_group <- paste0("LH", pair_data$LH_day)

      reference_level <- paste0(
        "LH",
        pair$reference_day
      )

      comparison_level <- paste0(
        "LH",
        pair$comparison_day
      )

      exact_test <- exact_two_group_test(
        y = pair_data$score,
        group = pair_group,
        reference_level = reference_level,
        comparison_level = comparison_level,
        seed = 20260805L +
          model_counter * 100L +
          pair_index
      )

      welch <- welch_summary(
        y = pair_data$score,
        group = pair_group,
        reference_level = reference_level,
        comparison_level = comparison_level
      )

      temporal_pairwise_rows[[
        length(temporal_pairwise_rows) + 1L
      ]] <- data.frame(
        broad_cell_type = cell_type,
        program_name = program_name,
        comparison = pair$comparison,
        reference_day = pair$reference_day,
        comparison_day = pair$comparison_day,
        reference_n = sum(
          pair_group == reference_level
        ),
        comparison_n = sum(
          pair_group == comparison_level
        ),
        mean_difference = exact_test$estimate,
        exact_permutation_p = exact_test$p,
        exact_test_mode = exact_test$test_mode,
        exact_permutation_n =
          exact_test$permutation_n,
        Welch_ci_low = welch$ci_low,
        Welch_ci_high = welch$ci_high,
        Welch_p = welch$p,
        Hedges_g = hedges_g(
          y = pair_data$score,
          group = pair_group,
          reference_level = reference_level,
          comparison_level = comparison_level
        ),
        stringsAsFactors = FALSE
      )
    }

    log_event(
      "DONE temporal program model",
      model_counter,
      "/",
      length(temporal_types) *
        length(program_names),
      cell_type,
      program_name,
      "permutation_p=",
      format_p(permutation$p)
    )
  }
}

temporal_program <- do.call(
  rbind,
  temporal_program_rows
)

temporal_program$FDR_BH <- stats::p.adjust(
  temporal_program$permutation_p,
  method = "BH"
)

temporal_program$trajectory_direction <- ifelse(
  temporal_program$LH11_mean >
    temporal_program$LH3_mean,
  "LH11_ABOVE_LH3",
  ifelse(
    temporal_program$LH11_mean <
      temporal_program$LH3_mean,
    "LH11_BELOW_LH3",
    "LH11_EQUALS_LH3"
  )
)

temporal_pairwise <- do.call(
  rbind,
  temporal_pairwise_rows
)

temporal_pairwise$FDR_BH_within_comparison <- ave(
  temporal_pairwise$exact_permutation_p,
  temporal_pairwise$comparison,
  FUN = function(p) {
    stats::p.adjust(p, method = "BH")
  }
)

temporal_pairwise$effect_direction <- sign_label(
  temporal_pairwise$mean_difference
)

temporal_program_summary <- do.call(
  rbind,
  temporal_summary_rows
)

write_csv_utf8(
  temporal_program,
  file.path(
    table_root,
    "01_temporal_program_omnibus.csv"
  )
)

write_csv_utf8(
  temporal_pairwise,
  file.path(
    table_root,
    "02_temporal_program_pairwise.csv"
  )
)

write_csv_utf8(
  temporal_program_summary,
  file.path(
    source_data_root,
    "source_temporal_program_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 6. LH+7 program-score models
# ----------------------------------------------------------------------

LH7_primary_rows <- list()
LH7_sensitivity_rows <- list()

model_counter <- 0L

for (cell_type in LH7_types) {
  for (program_name in program_names) {
    model_counter <- model_counter + 1L

    one <- LH7_score_data[
      LH7_score_data$broad_cell_type == cell_type &
        LH7_score_data$program_name == program_name,
      ,
      drop = FALSE
    ]

    one$clinical_group <- factor(
      one$clinical_group,
      levels = c("FERTILE", "RIF")
    )

    exact_test <- exact_two_group_test(
      y = one$score,
      group = one$clinical_group,
      reference_level = "FERTILE",
      comparison_level = "RIF",
      seed = 20260805L + model_counter
    )

    welch <- welch_summary(
      y = one$score,
      group = one$clinical_group,
      reference_level = "FERTILE",
      comparison_level = "RIF"
    )

    loo <- leave_one_out_difference(
      y = one$score,
      group = one$clinical_group,
      reference_level = "FERTILE",
      comparison_level = "RIF"
    )

    LH7_primary_rows[[
      length(LH7_primary_rows) + 1L
    ]] <- data.frame(
      broad_cell_type = cell_type,
      program_name = program_name,
      fertile_n = sum(
        one$clinical_group == "FERTILE"
      ),
      RIF_n = sum(
        one$clinical_group == "RIF"
      ),
      fertile_mean = mean(
        one$score[
          one$clinical_group == "FERTILE"
        ]
      ),
      RIF_mean = mean(
        one$score[
          one$clinical_group == "RIF"
        ]
      ),
      RIF_minus_fertile = exact_test$estimate,
      exact_permutation_p = exact_test$p,
      exact_test_mode = exact_test$test_mode,
      exact_permutation_n =
        exact_test$permutation_n,
      Welch_ci_low = welch$ci_low,
      Welch_ci_high = welch$ci_high,
      Welch_p = welch$p,
      Hedges_g = hedges_g(
        y = one$score,
        group = one$clinical_group,
        reference_level = "FERTILE",
        comparison_level = "RIF"
      ),
      LOO_min = unname(loo["loo_min"]),
      LOO_max = unname(loo["loo_max"]),
      LOO_sign_consistency = unname(
        loo["loo_sign_consistency"]
      ),
      stringsAsFactors = FALSE
    )

    for (
      fertile_subgroup in c(
        "LH7_STANDARD",
        "LH7_PGT"
      )
    ) {
      sensitivity <- one[
        one$clinical_group == "RIF" |
          (
            one$clinical_group == "FERTILE" &
              one$LH7_fertile_subgroup ==
                fertile_subgroup
          ),
        ,
        drop = FALSE
      ]

      sensitivity_group <- ifelse(
        sensitivity$clinical_group == "RIF",
        "RIF",
        fertile_subgroup
      )

      sensitivity_test <- exact_two_group_test(
        y = sensitivity$score,
        group = sensitivity_group,
        reference_level = fertile_subgroup,
        comparison_level = "RIF",
        seed = 20260805L +
          model_counter * 100L +
          match(
            fertile_subgroup,
            c("LH7_STANDARD", "LH7_PGT")
          )
      )

      sensitivity_welch <- welch_summary(
        y = sensitivity$score,
        group = sensitivity_group,
        reference_level = fertile_subgroup,
        comparison_level = "RIF"
      )

      LH7_sensitivity_rows[[
        length(LH7_sensitivity_rows) + 1L
      ]] <- data.frame(
        broad_cell_type = cell_type,
        program_name = program_name,
        sensitivity_reference =
          fertile_subgroup,
        fertile_subgroup_n = sum(
          sensitivity_group ==
            fertile_subgroup
        ),
        RIF_n = sum(
          sensitivity_group == "RIF"
        ),
        RIF_minus_reference =
          sensitivity_test$estimate,
        exact_permutation_p =
          sensitivity_test$p,
        exact_test_mode =
          sensitivity_test$test_mode,
        exact_permutation_n =
          sensitivity_test$permutation_n,
        Welch_ci_low =
          sensitivity_welch$ci_low,
        Welch_ci_high =
          sensitivity_welch$ci_high,
        Welch_p =
          sensitivity_welch$p,
        Hedges_g = hedges_g(
          y = sensitivity$score,
          group = sensitivity_group,
          reference_level = fertile_subgroup,
          comparison_level = "RIF"
        ),
        stringsAsFactors = FALSE
      )
    }

    log_event(
      "DONE LH7 program model",
      model_counter,
      "/",
      length(LH7_types) *
        length(program_names),
      cell_type,
      program_name,
      "effect=",
      round(exact_test$estimate, 4),
      "p=",
      format_p(exact_test$p)
    )
  }
}

LH7_primary <- do.call(
  rbind,
  LH7_primary_rows
)

LH7_primary$FDR_BH <- stats::p.adjust(
  LH7_primary$exact_permutation_p,
  method = "BH"
)

LH7_primary$effect_direction <- sign_label(
  LH7_primary$RIF_minus_fertile
)

LH7_primary$LOO_stable_sign <- (
  LH7_primary$LOO_sign_consistency >= 0.80
)

LH7_sensitivity <- do.call(
  rbind,
  LH7_sensitivity_rows
)

LH7_sensitivity$FDR_BH_within_reference <- ave(
  LH7_sensitivity$exact_permutation_p,
  LH7_sensitivity$sensitivity_reference,
  FUN = function(p) {
    stats::p.adjust(p, method = "BH")
  }
)

primary_lookup <- LH7_primary[
  ,
  c(
    "broad_cell_type",
    "program_name",
    "RIF_minus_fertile"
  )
]

colnames(primary_lookup)[3] <-
  "primary_RIF_minus_fertile"

LH7_sensitivity <- merge(
  LH7_sensitivity,
  primary_lookup,
  by = c(
    "broad_cell_type",
    "program_name"
  ),
  all.x = TRUE,
  sort = FALSE
)

LH7_sensitivity$sign_concordant_with_primary <- (
  sign(
    LH7_sensitivity$RIF_minus_reference
  ) ==
    sign(
      LH7_sensitivity$
        primary_RIF_minus_fertile
    )
)

write_csv_utf8(
  LH7_primary,
  file.path(
    table_root,
    "03_LH7_program_primary.csv"
  )
)

write_csv_utf8(
  LH7_sensitivity,
  file.path(
    table_root,
    "04_LH7_program_sensitivity.csv"
  )
)


# ----------------------------------------------------------------------
# 7. VDR-expression models
# ----------------------------------------------------------------------

temporal_VDR_rows <- list()
temporal_VDR_pair_rows <- list()
temporal_VDR_summary_rows <- list()

for (
  type_index in seq_along(temporal_types)
) {
  cell_type <- temporal_types[type_index]

  one <- temporal_VDR_data[
    temporal_VDR_data$broad_cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]

  one$LH_factor <- factor(
    one$LH_day,
    levels = c(3, 5, 7, 9, 11)
  )

  expected_counts <- expected_temporal_counts(
    eligibility_table = eligibility,
    cell_type = cell_type
  )

  observed_counts <- as.integer(
    table(one$LH_factor)
  )

  if (
    nrow(one) != sum(expected_counts) ||
      !identical(
        observed_counts,
        as.integer(expected_counts)
      ) ||
      any(observed_counts < 2L)
  ) {
    stop(
      "Temporal VDR-expression input mismatch for ",
      cell_type
    )
  }

  permutation <- permutation_oneway(
    y = one$VDR_log2_CPM_plus_1,
    group = one$LH_factor,
    n_perm = 5000L,
    seed = 20260805L + 10000L + type_index
  )

  fit <- stats::lm(
    VDR_log2_CPM_plus_1 ~ LH_factor,
    data = one
  )

  anova_table <- stats::anova(fit)

  group_means <- tapply(
    one$VDR_log2_CPM_plus_1,
    one$LH_factor,
    mean
  )

  total_ss <- sum(
    (
      one$VDR_log2_CPM_plus_1 -
        mean(one$VDR_log2_CPM_plus_1)
    )^2
  )

  fitted_ss <- sum(
    (
      fitted(fit) -
        mean(one$VDR_log2_CPM_plus_1)
    )^2
  )

  temporal_VDR_rows[[
    length(temporal_VDR_rows) + 1L
  ]] <- data.frame(
    broad_cell_type = cell_type,
    donor_n = nrow(one),
    LH3_mean = unname(group_means["3"]),
    LH5_mean = unname(group_means["5"]),
    LH7_mean = unname(group_means["7"]),
    LH9_mean = unname(group_means["9"]),
    LH11_mean = unname(group_means["11"]),
    permutation_F = permutation$F,
    permutation_p = permutation$p,
    permutation_n = permutation$permutation_n,
    classical_F =
      anova_table$`F value`[1],
    classical_p =
      anova_table$`Pr(>F)`[1],
    eta_squared = if (
      is.finite(total_ss) &&
        total_ss > 0
    ) {
      fitted_ss / total_ss
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  summary_one <- group_summary(
    data = one,
    value_column = "VDR_log2_CPM_plus_1",
    group_columns = c(
      "broad_cell_type",
      "LH_day"
    )
  )

  temporal_VDR_summary_rows[[
    length(temporal_VDR_summary_rows) + 1L
  ]] <- as.data.frame(summary_one)

  for (pair_index in seq_len(nrow(planned_pairs))) {
    pair <- planned_pairs[pair_index, ]

    pair_data <- one[
      one$LH_day %in% c(
        pair$reference_day,
        pair$comparison_day
      ),
      ,
      drop = FALSE
    ]

    pair_group <- paste0("LH", pair_data$LH_day)
    reference_level <- paste0("LH", pair$reference_day)
    comparison_level <- paste0(
      "LH",
      pair$comparison_day
    )

    exact_test <- exact_two_group_test(
      y = pair_data$VDR_log2_CPM_plus_1,
      group = pair_group,
      reference_level = reference_level,
      comparison_level = comparison_level,
      seed = 20260805L +
        11000L +
        type_index * 100L +
        pair_index
    )

    welch <- welch_summary(
      y = pair_data$VDR_log2_CPM_plus_1,
      group = pair_group,
      reference_level = reference_level,
      comparison_level = comparison_level
    )

    temporal_VDR_pair_rows[[
      length(temporal_VDR_pair_rows) + 1L
    ]] <- data.frame(
      broad_cell_type = cell_type,
      comparison = pair$comparison,
      reference_day = pair$reference_day,
      comparison_day = pair$comparison_day,
      mean_difference = exact_test$estimate,
      exact_permutation_p = exact_test$p,
      exact_test_mode = exact_test$test_mode,
      exact_permutation_n =
        exact_test$permutation_n,
      Welch_ci_low = welch$ci_low,
      Welch_ci_high = welch$ci_high,
      Welch_p = welch$p,
      Hedges_g = hedges_g(
        y = pair_data$VDR_log2_CPM_plus_1,
        group = pair_group,
        reference_level = reference_level,
        comparison_level = comparison_level
      ),
      stringsAsFactors = FALSE
    )
  }
}

temporal_VDR <- do.call(
  rbind,
  temporal_VDR_rows
)

temporal_VDR$FDR_BH <- stats::p.adjust(
  temporal_VDR$permutation_p,
  method = "BH"
)

temporal_VDR_pairwise <- do.call(
  rbind,
  temporal_VDR_pair_rows
)

temporal_VDR_pairwise$FDR_BH_within_comparison <- ave(
  temporal_VDR_pairwise$exact_permutation_p,
  temporal_VDR_pairwise$comparison,
  FUN = function(p) {
    stats::p.adjust(p, method = "BH")
  }
)

temporal_VDR_summary <- do.call(
  rbind,
  temporal_VDR_summary_rows
)

write_csv_utf8(
  temporal_VDR,
  file.path(
    table_root,
    "05_temporal_VDR_omnibus.csv"
  )
)

write_csv_utf8(
  temporal_VDR_pairwise,
  file.path(
    table_root,
    "06_temporal_VDR_pairwise.csv"
  )
)

write_csv_utf8(
  temporal_VDR_summary,
  file.path(
    source_data_root,
    "source_temporal_VDR_summary.csv"
  )
)

LH7_VDR_primary_rows <- list()
LH7_VDR_sensitivity_rows <- list()

for (type_index in seq_along(LH7_types)) {
  cell_type <- LH7_types[type_index]

  one <- LH7_VDR_data[
    LH7_VDR_data$broad_cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]

  one$clinical_group <- factor(
    one$clinical_group,
    levels = c("FERTILE", "RIF")
  )

  exact_test <- exact_two_group_test(
    y = one$VDR_log2_CPM_plus_1,
    group = one$clinical_group,
    reference_level = "FERTILE",
    comparison_level = "RIF",
    seed = 20260805L + 12000L + type_index
  )

  welch <- welch_summary(
    y = one$VDR_log2_CPM_plus_1,
    group = one$clinical_group,
    reference_level = "FERTILE",
    comparison_level = "RIF"
  )

  loo <- leave_one_out_difference(
    y = one$VDR_log2_CPM_plus_1,
    group = one$clinical_group,
    reference_level = "FERTILE",
    comparison_level = "RIF"
  )

  LH7_VDR_primary_rows[[
    length(LH7_VDR_primary_rows) + 1L
  ]] <- data.frame(
    broad_cell_type = cell_type,
    fertile_n = sum(
      one$clinical_group == "FERTILE"
    ),
    RIF_n = sum(
      one$clinical_group == "RIF"
    ),
    fertile_mean = mean(
      one$VDR_log2_CPM_plus_1[
        one$clinical_group == "FERTILE"
      ]
    ),
    RIF_mean = mean(
      one$VDR_log2_CPM_plus_1[
        one$clinical_group == "RIF"
      ]
    ),
    RIF_minus_fertile = exact_test$estimate,
    exact_permutation_p = exact_test$p,
    exact_test_mode = exact_test$test_mode,
    exact_permutation_n =
      exact_test$permutation_n,
    Welch_ci_low = welch$ci_low,
    Welch_ci_high = welch$ci_high,
    Welch_p = welch$p,
    Hedges_g = hedges_g(
      y = one$VDR_log2_CPM_plus_1,
      group = one$clinical_group,
      reference_level = "FERTILE",
      comparison_level = "RIF"
    ),
    LOO_min = unname(loo["loo_min"]),
    LOO_max = unname(loo["loo_max"]),
    LOO_sign_consistency = unname(
      loo["loo_sign_consistency"]
    ),
    stringsAsFactors = FALSE
  )

  for (
    fertile_subgroup in c(
      "LH7_STANDARD",
      "LH7_PGT"
    )
  ) {
    sensitivity <- one[
      one$clinical_group == "RIF" |
        (
          one$clinical_group == "FERTILE" &
            one$LH7_fertile_subgroup ==
              fertile_subgroup
        ),
      ,
      drop = FALSE
    ]

    sensitivity_group <- ifelse(
      sensitivity$clinical_group == "RIF",
      "RIF",
      fertile_subgroup
    )

    exact_sensitivity <- exact_two_group_test(
      y = sensitivity$VDR_log2_CPM_plus_1,
      group = sensitivity_group,
      reference_level = fertile_subgroup,
      comparison_level = "RIF",
      seed = 20260805L +
        13000L +
        type_index * 100L +
        match(
          fertile_subgroup,
          c("LH7_STANDARD", "LH7_PGT")
        )
    )

    welch_sensitivity <- welch_summary(
      y = sensitivity$VDR_log2_CPM_plus_1,
      group = sensitivity_group,
      reference_level = fertile_subgroup,
      comparison_level = "RIF"
    )

    LH7_VDR_sensitivity_rows[[
      length(LH7_VDR_sensitivity_rows) + 1L
    ]] <- data.frame(
      broad_cell_type = cell_type,
      sensitivity_reference =
        fertile_subgroup,
      fertile_subgroup_n = sum(
        sensitivity_group ==
          fertile_subgroup
      ),
      RIF_n = sum(
        sensitivity_group == "RIF"
      ),
      RIF_minus_reference =
        exact_sensitivity$estimate,
      exact_permutation_p =
        exact_sensitivity$p,
      exact_test_mode =
        exact_sensitivity$test_mode,
      exact_permutation_n =
        exact_sensitivity$permutation_n,
      Welch_ci_low =
        welch_sensitivity$ci_low,
      Welch_ci_high =
        welch_sensitivity$ci_high,
      Welch_p =
        welch_sensitivity$p,
      Hedges_g = hedges_g(
        y = sensitivity$VDR_log2_CPM_plus_1,
        group = sensitivity_group,
        reference_level = fertile_subgroup,
        comparison_level = "RIF"
      ),
      stringsAsFactors = FALSE
    )
  }
}

LH7_VDR_primary <- do.call(
  rbind,
  LH7_VDR_primary_rows
)

LH7_VDR_primary$FDR_BH <- stats::p.adjust(
  LH7_VDR_primary$exact_permutation_p,
  method = "BH"
)

LH7_VDR_primary$effect_direction <- sign_label(
  LH7_VDR_primary$RIF_minus_fertile
)

LH7_VDR_primary$LOO_stable_sign <- (
  LH7_VDR_primary$LOO_sign_consistency >= 0.80
)

LH7_VDR_sensitivity <- do.call(
  rbind,
  LH7_VDR_sensitivity_rows
)

LH7_VDR_sensitivity$FDR_BH_within_reference <- ave(
  LH7_VDR_sensitivity$exact_permutation_p,
  LH7_VDR_sensitivity$sensitivity_reference,
  FUN = function(p) {
    stats::p.adjust(p, method = "BH")
  }
)

write_csv_utf8(
  LH7_VDR_primary,
  file.path(
    table_root,
    "07_LH7_VDR_primary.csv"
  )
)

write_csv_utf8(
  LH7_VDR_sensitivity,
  file.path(
    table_root,
    "08_LH7_VDR_sensitivity.csv"
  )
)


# ----------------------------------------------------------------------
# 8. Cross-model summary
# ----------------------------------------------------------------------

significant_temporal_program <- temporal_program[
  temporal_program$FDR_BH < 0.05,
  ,
  drop = FALSE
]

significant_LH7_program <- LH7_primary[
  LH7_primary$FDR_BH < 0.05,
  ,
  drop = FALSE
]

significant_temporal_VDR <- temporal_VDR[
  temporal_VDR$FDR_BH < 0.05,
  ,
  drop = FALSE
]

significant_LH7_VDR <- LH7_VDR_primary[
  LH7_VDR_primary$FDR_BH < 0.05,
  ,
  drop = FALSE
]

cross_model_summary <- data.frame(
  model_family = c(
    "TEMPORAL_VDR_PROGRAM",
    "LH7_RIF_VS_FERTILE_VDR_PROGRAM",
    "TEMPORAL_VDR_EXPRESSION",
    "LH7_RIF_VS_FERTILE_VDR_EXPRESSION"
  ),
  tested_hypothesis_n = c(
    nrow(temporal_program),
    nrow(LH7_primary),
    nrow(temporal_VDR),
    nrow(LH7_VDR_primary)
  ),
  FDR_lt_0_05_n = c(
    nrow(significant_temporal_program),
    nrow(significant_LH7_program),
    nrow(significant_temporal_VDR),
    nrow(significant_LH7_VDR)
  ),
  minimum_raw_p = c(
    min(
      temporal_program$permutation_p,
      na.rm = TRUE
    ),
    min(
      LH7_primary$exact_permutation_p,
      na.rm = TRUE
    ),
    min(
      temporal_VDR$permutation_p,
      na.rm = TRUE
    ),
    min(
      LH7_VDR_primary$exact_permutation_p,
      na.rm = TRUE
    )
  ),
  minimum_FDR = c(
    min(temporal_program$FDR_BH, na.rm = TRUE),
    min(LH7_primary$FDR_BH, na.rm = TRUE),
    min(temporal_VDR$FDR_BH, na.rm = TRUE),
    min(LH7_VDR_primary$FDR_BH, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  cross_model_summary,
  file.path(
    table_root,
    "09_cross_model_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 9. QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(
      fill = "grey95"
    )
  )

temporal_program_summary$LH_day <- factor(
  temporal_program_summary$LH_day,
  levels = c(3, 5, 7, 9, 11),
  labels = c(
    "LH+3",
    "LH+5",
    "LH+7",
    "LH+9",
    "LH+11"
  )
)

plot_temporal_program <- ggplot2::ggplot(
  temporal_program_summary,
  ggplot2::aes(
    x = LH_day,
    y = mean,
    group = 1
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = ci_low,
      ymax = ci_high,
      group = 1
    ),
    alpha = 0.15
  ) +
  ggplot2::geom_line(linewidth = 0.55) +
  ggplot2::geom_point(size = 1.3) +
  ggplot2::facet_grid(
    program_name ~ broad_cell_type,
    scales = "free_y"
  ) +
  ggplot2::labs(
    title = "Fertile temporal trajectories of locked VDR programs",
    subtitle = "Donor-level pseudobulk scores; ribbons show 95% mean confidence intervals",
    x = "Biopsy timing",
    y = "Signed VDR-program score"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text = ggplot2::element_text(
      size = 6
    )
  )

save_plot_pair(
  plot_temporal_program,
  "QC1_temporal_programs",
  width = 16,
  height = 11
)

LH7_primary$FDR_label <- ifelse(
  LH7_primary$FDR_BH < 0.05,
  "FDR < 0.05",
  "FDR >= 0.05"
)

plot_LH7_program <- ggplot2::ggplot(
  LH7_primary,
  ggplot2::aes(
    x = RIF_minus_fertile,
    y = broad_cell_type
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = Welch_ci_low,
      xend = Welch_ci_high,
      y = broad_cell_type,
      yend = broad_cell_type
    ),
    linewidth = 0.5
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      shape = FDR_label,
      size = abs(Hedges_g)
    )
  ) +
  ggplot2::facet_wrap(
    ~ program_name,
    scales = "free_x",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "LH+7 RIF versus fertile VDR-program effects",
    subtitle = "Positive values indicate greater alignment in RIF; intervals are Welch 95% CIs",
    x = "RIF minus fertile score",
    y = "Broad cell type",
    shape = "Permutation FDR",
    size = "|Hedges g|"
  ) +
  theme_qc

save_plot_pair(
  plot_LH7_program,
  "QC2_LH7_program_forest",
  width = 11,
  height = 9
)

temporal_VDR_summary$LH_day <- factor(
  temporal_VDR_summary$LH_day,
  levels = c(3, 5, 7, 9, 11),
  labels = c(
    "LH+3",
    "LH+5",
    "LH+7",
    "LH+9",
    "LH+11"
  )
)

plot_temporal_VDR <- ggplot2::ggplot(
  temporal_VDR_summary,
  ggplot2::aes(
    x = LH_day,
    y = mean,
    group = 1
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = ci_low,
      ymax = ci_high,
      group = 1
    ),
    alpha = 0.15
  ) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::facet_wrap(
    ~ broad_cell_type,
    scales = "free_y"
  ) +
  ggplot2::labs(
    title = "Fertile temporal trajectory of VDR expression",
    subtitle = "Donor-level pseudobulk log2(CPM + 1)",
    x = "Biopsy timing",
    y = "VDR log2(CPM + 1)"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot_pair(
  plot_temporal_VDR,
  "QC3_temporal_VDR",
  width = 11,
  height = 7
)

LH7_VDR_primary$FDR_label <- ifelse(
  LH7_VDR_primary$FDR_BH < 0.05,
  "FDR < 0.05",
  "FDR >= 0.05"
)

plot_LH7_VDR <- ggplot2::ggplot(
  LH7_VDR_primary,
  ggplot2::aes(
    x = RIF_minus_fertile,
    y = broad_cell_type
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = Welch_ci_low,
      xend = Welch_ci_high,
      y = broad_cell_type,
      yend = broad_cell_type
    ),
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      shape = FDR_label,
      size = abs(Hedges_g)
    )
  ) +
  ggplot2::labs(
    title = "LH+7 RIF versus fertile VDR expression",
    subtitle = "Positive values indicate higher VDR expression in RIF",
    x = "RIF minus fertile log2(CPM + 1)",
    y = "Broad cell type",
    shape = "Permutation FDR",
    size = "|Hedges g|"
  ) +
  theme_qc

save_plot_pair(
  plot_LH7_VDR,
  "QC4_LH7_VDR_forest",
  width = 9,
  height = 6
)

sensitivity_plot_data <- LH7_sensitivity
sensitivity_plot_data$comparison_label <- ifelse(
  sensitivity_plot_data$sensitivity_reference ==
    "LH7_STANDARD",
  "RIF vs standard fertile",
  "RIF vs PGT fertile"
)

plot_sensitivity <- ggplot2::ggplot(
  sensitivity_plot_data,
  ggplot2::aes(
    x = primary_RIF_minus_fertile,
    y = RIF_minus_reference,
    shape = comparison_label
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(
    ~ program_name,
    scales = "free"
  ) +
  ggplot2::labs(
    title = "LH+7 fertile-subgroup sensitivity",
    subtitle = "Concordance of primary and subgroup-specific RIF effects",
    x = "Primary: RIF minus all fertile",
    y = "Sensitivity: RIF minus fertile subgroup",
    shape = "Sensitivity comparison"
  ) +
  theme_qc

save_plot_pair(
  plot_sensitivity,
  "QC5_LH7_sensitivity",
  width = 10,
  height = 8
)


# ----------------------------------------------------------------------
# 10. Gate, Card, interpretation, Methods and Results shell
# ----------------------------------------------------------------------

gate_checks <- data.frame(
  check = c(
    "upstream_gate_authorizes_donor_models",
    "pseudobulk_status_complete",
    "temporal_type_n_is_5",
    "LH7_type_n_is_6",
    "four_locked_programs_tested",
    "temporal_program_model_n_is_20",
    "LH7_program_model_n_is_24",
    "temporal_VDR_model_n_is_5",
    "LH7_VDR_model_n_is_6",
    "all_temporal_program_models_have_permutation_p",
    "all_LH7_program_models_have_exact_p",
    "all_primary_LH7_models_have_leave_one_out_QC",
    "all_FDR_values_generated",
    "composition_inference_remains_disabled",
    "subtype_inference_remains_disabled"
  ),
  pass = c(
    pb_gate$donor_level_modeling_authorized[[1]],
    pb_status$status[[1]] == "COMPLETE",
    length(temporal_types) == 5L,
    length(LH7_types) == 6L,
    setequal(program_names, expected_programs),
    nrow(temporal_program) == 20L,
    nrow(LH7_primary) == 24L,
    nrow(temporal_VDR) == 5L,
    nrow(LH7_VDR_primary) == 6L,
    all(is.finite(
      temporal_program$permutation_p
    )),
    all(is.finite(
      LH7_primary$exact_permutation_p
    )),
    all(is.finite(
      LH7_primary$LOO_sign_consistency
    )),
    all(is.finite(c(
      temporal_program$FDR_BH,
      LH7_primary$FDR_BH,
      temporal_VDR$FDR_BH,
      LH7_VDR_primary$FDR_BH
    ))),
    !pb_gate$cell_composition_inference_authorized[[1]],
    !pb_gate$subtype_modeling_authorized[[1]]
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_checks,
  file.path(run_qc_root, "02_gate_checks.csv")
)

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
  "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE"
} else {
  "HOLD_DONOR_LEVEL_VDR_MODEL_QC_FAILED"
}

model_gate <- data.frame(
  dataset_id = dataset_id,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  model_version = model_version,
  gate_status = gate_status,
  temporal_cell_type_n =
    length(temporal_types),
  LH7_cell_type_n = length(LH7_types),
  program_n = length(program_names),
  temporal_program_model_n =
    nrow(temporal_program),
  temporal_program_FDR_lt_0_05_n =
    nrow(significant_temporal_program),
  LH7_program_model_n =
    nrow(LH7_primary),
  LH7_program_FDR_lt_0_05_n =
    nrow(significant_LH7_program),
  temporal_VDR_model_n =
    nrow(temporal_VDR),
  temporal_VDR_FDR_lt_0_05_n =
    nrow(significant_temporal_VDR),
  LH7_VDR_model_n =
    nrow(LH7_VDR_primary),
  LH7_VDR_FDR_lt_0_05_n =
    nrow(significant_LH7_VDR),
  gene_level_pseudobulk_DE_authorized =
    all_gate_checks_pass,
  composition_inference_authorized =
    FALSE,
  subtype_inference_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "VDR_MODEL_GATE.csv"
)

write_csv_utf8(
  model_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1G Donor-level VDR Model Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Pseudobulk version: ", pseudobulk_version),
  paste0("- Model version: ", model_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Temporal program tests: ",
    nrow(temporal_program)
  ),
  paste0(
    "- Temporal program FDR < 0.05: ",
    nrow(significant_temporal_program)
  ),
  paste0(
    "- LH+7 program tests: ",
    nrow(LH7_primary)
  ),
  paste0(
    "- LH+7 program FDR < 0.05: ",
    nrow(significant_LH7_program)
  ),
  paste0(
    "- Temporal VDR-expression FDR < 0.05: ",
    nrow(significant_temporal_VDR)
  ),
  paste0(
    "- LH+7 VDR-expression FDR < 0.05: ",
    nrow(significant_LH7_VDR)
  ),
  paste0(
    "- Gene-level pseudobulk DE authorized: ",
    model_gate$gene_level_pseudobulk_DE_authorized
  ),
  "- Cell-composition inference authorized: FALSE",
  "- Subtype inference authorized: FALSE",
  "- Result lock authorized: FALSE"
)

card_path <- file.path(
  run_qc_root,
  "VDR_MODEL_RESULT_CARD.md"
)

write_text(card_lines, card_path)

interpretation_lines <- c(
  "# STEP04E1G donor-level VDR model interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Model scope",
  "",
  "- Each observation represents one donor-by-broad-cell-type pseudobulk.",
  "- Fertile temporal models compare cross-sectional LH+3, LH+5, LH+7, LH+9 and LH+11 donor groups.",
  "- The primary LH+7 comparison contrasts RIF with all fertile LH+7 donors.",
  "- Separate standard-fertile and PGT-fertile comparisons are sensitivity analyses.",
  "- Program scores indicate relative alignment with the locked VDR interaction direction and are not causal activation scores.",
  "",
  "## Statistical scope",
  "",
  "- Permutation tests provide the primary p-values.",
  "- Benjamini-Hochberg correction is applied within predefined model families.",
  "- Welch confidence intervals and Hedges g are retained as secondary effect-size diagnostics.",
  "- Leave-one-out sign stability is reported for primary LH+7 comparisons.",
  "",
  "## Next step",
  "",
  "- Review temporal and LH+7 VDR-program effects.",
  "- Proceed to gene-level pseudobulk differential expression only after this Gate passes.",
  "- Integrate gene-level effects with the locked 436/284/46/84 candidate hierarchy.",
  "- Do not perform cell-composition or subtype-level inference."
)

interpretation_path <- file.path(
  run_qc_root,
  "VDR_MODEL_INTERPRETATION.md"
)

write_text(
  interpretation_lines,
  interpretation_path
)

methods_lines <- c(
  "# STEP04E1G Methods draft",
  "",
  paste(
    "Donor-level analyses were performed using one pseudobulk profile per",
    "donor and broad cell type."
  ),
  paste(
    "Four frozen VDR RNA programs were examined: the full 436-gene",
    "interaction program, the 284-gene occupancy-supported program, the",
    "46-gene shared-occupancy subset and the 84-gene promoter-occupancy",
    "subset."
  ),
  paste(
    "Fertile temporal trajectories were evaluated across LH+3, LH+5, LH+7,",
    "LH+9 and LH+11 using the actual pseudobulk profiles that passed the",
    "predefined cell-count threshold."
  ),
  paste(
    "The temporal design was allowed to be unbalanced provided that each",
    "time point retained at least two independent donors; specifically,",
    "ciliated epithelial cells included two eligible LH+11 donors whereas",
    "the remaining eligible temporal cell types retained three LH+11 donors."
  ),
  paste(
    "Temporal omnibus evidence was assessed using a 5,000-label permutation",
    "one-way analysis of variance, with classical one-way ANOVA retained as",
    "a secondary diagnostic."
  ),
  paste(
    "At LH+7, RIF and fertile donors were compared using an exact",
    "two-group permutation test of the mean difference whenever the number",
    "of possible allocations was computationally tractable."
  ),
  paste(
    "Welch confidence intervals, Hedges standardized effect sizes and",
    "leave-one-out sign stability were reported as complementary",
    "diagnostics."
  ),
  paste(
    "Sensitivity analyses compared RIF donors separately with standard",
    "fertile LH+7 donors and PGT-associated fertile LH+7 donors."
  ),
  paste(
    "Benjamini-Hochberg adjustment was applied within each predefined",
    "model family. Parallel analyses were performed for VDR expression",
    "measured as log2 counts per million plus one."
  )
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Methods_draft.md"
  )
)

results_shell <- c(
  "# STEP04E1G Results shell",
  "",
  paste0(
    "The donor-level VDR-program analysis included ",
    length(temporal_types),
    " broad cell types in the fertile temporal series and ",
    length(LH7_types),
    " broad cell types in the LH+7 RIF comparison."
  ),
  paste0(
    "Across ",
    nrow(temporal_program),
    " temporal program tests, ",
    nrow(significant_temporal_program),
    " met FDR < 0.05."
  ),
  paste0(
    "Across ",
    nrow(LH7_primary),
    " primary LH+7 program comparisons, ",
    nrow(significant_LH7_program),
    " met FDR < 0.05."
  ),
  paste0(
    "For VDR expression, ",
    nrow(significant_temporal_VDR),
    " temporal cell-type model(s) and ",
    nrow(significant_LH7_VDR),
    " LH+7 RIF-versus-fertile model(s) met FDR < 0.05."
  ),
  "",
  "Detailed biological interpretation will be drafted after review of effect directions, sensitivity analyses and gene-level pseudobulk results."
)

write_text(
  results_shell,
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Results_shell.md"
  )
)


# ----------------------------------------------------------------------
# 11. Manifests, sessionInfo and ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_names <- c(
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

package_table <- data.frame(
  package = package_names,
  version = vapply(
    package_names,
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(run_qc_root, "PACKAGE_VERSIONS.csv")
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- unique(c(
  list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Results_shell.md"
  )
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(artifact_files),
  size_bytes = file.info(artifact_files)$size,
  md5 = unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(run_qc_root, "ARTIFACT_MANIFEST.csv")
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  model_version = model_version,
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_seconds = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "secs"
      )
    ),
    2
  ),
  gate_status = gate_status,
  temporal_program_model_n =
    nrow(temporal_program),
  LH7_program_model_n =
    nrow(LH7_primary),
  temporal_VDR_model_n =
    nrow(temporal_VDR),
  LH7_VDR_model_n =
    nrow(LH7_VDR_primary),
  donor_level_modeling_completed =
    all_gate_checks_pass,
  gene_level_pseudobulk_DE_authorized =
    model_gate$gene_level_pseudobulk_DE_authorized,
  composition_inference_authorized =
    FALSE,
  subtype_inference_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  output_root =
    normalize_slashes(output_root),
  qc_root =
    normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "temporal_program_FDR_lt_0_05_n=",
  nrow(significant_temporal_program),
  "LH7_program_FDR_lt_0_05_n=",
  nrow(significant_LH7_program)
)

small_files <- unique(c(
  list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1G_VDR_models_Results_shell.md"
  ),
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <=
      200 * 1024^2
]

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = small_files,
        recurse = FALSE,
        include_directories = FALSE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1G GSE250130 donor-level VDR models completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Temporal program models: ",
  nrow(temporal_program),
  "\n",
  sep = ""
)
cat(
  "Temporal program FDR < 0.05: ",
  nrow(significant_temporal_program),
  "\n",
  sep = ""
)
cat(
  "LH7 program models: ",
  nrow(LH7_primary),
  "\n",
  sep = ""
)
cat(
  "LH7 program FDR < 0.05: ",
  nrow(significant_LH7_program),
  "\n",
  sep = ""
)
cat(
  "Temporal VDR-expression FDR < 0.05: ",
  nrow(significant_temporal_VDR),
  "\n",
  sep = ""
)
cat(
  "LH7 VDR-expression FDR < 0.05: ",
  nrow(significant_LH7_VDR),
  "\n",
  sep = ""
)
cat("Output root: ", output_root, "\n", sep = "")
cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
  cat(
    "Upload this ZIP:\n",
    handoff_zip,
    "\n"
  )
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe model and QC folders are complete.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 8: 62_step04E1H_GSE250130_gene_pseudobulk_DE_v1_2_pathsafe_resume
# ==============================================================================

# ======================================================================
# 62_step04E1H_GSE250130_gene_pseudobulk_DE_v1_2_pathsafe_resume.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Perform gene-level donor-aware pseudobulk differential-expression
#   analysis for GSE250130.
#
# PREDEFINED MODEL FAMILIES
#
# A. Fertile temporal omnibus:
#    LH+3, LH+5, LH+7, LH+9 and LH+11
#    Eligible broad cell types: 5
#
# B. Fertile temporal planned pairwise contrasts:
#    LH5_vs_LH3
#    LH7_vs_LH5
#    LH9_vs_LH7
#    LH11_vs_LH9
#    LH7_vs_LH3
#    LH11_vs_LH3
#    LH11_vs_LH7
#
# C. LH+7 RIF versus all fertile:
#    Eligible broad cell types: 6
#
# D. LH+7 sensitivity:
#    RIF versus fertile LH7_STANDARD
#    RIF versus fertile LH7_PGT
#
# STATISTICAL WORKFLOW
# - One pseudobulk profile = one biological donor x broad cell type.
# - edgeR filterByExpr with the cell-type-specific design.
# - TMM normalization.
# - Robust negative-binomial quasi-likelihood GLM.
# - glmQLFTest for omnibus and one-degree-of-freedom contrasts.
# - Genome-wide BH FDR is the primary gene-level multiplicity control.
# - Separate focused BH values are reported for the prospectively locked
#   436/284/46/84 VDR candidate hierarchies; these are clearly labelled
#   candidate-restricted and do not replace genome-wide FDR.
# - Directionally weighted limma fry tests assess whether each locked VDR
#   program changes concordantly with the calcitriol-dependent interaction
#   direction for one-degree-of-freedom contrasts.
#
# SCIENTIFIC GUARDRAILS
# - Cell composition and subtype-level inference remain unauthorized.
# - The temporal design is cross-sectional, not repeated-measures.
# - No disease claim is made solely from candidate-restricted FDR.
# - No result is written to 06_locked_results.
#
# INPUT GATE
# - Requires STEP04E1G:
#     PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE
#     gene_level_pseudobulk_DE_authorized = TRUE
#
# OUTPUT STATUS
# - CANDIDATE_GENE_LEVEL_PSEUDOBULK_RESULTS
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)
options(timeout = max(14400L, getOption("timeout")))

set.seed(20260805)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE250130"

script_name <- "62_step04E1H_GSE250130_gene_pseudobulk_DE_v1_2_pathsafe_resume.R"
script_version <- "STEP04E1H_GSE250130_GENE_PSEUDOBULK_DE_v1.2_PATHSAFE_RESUME"
input_version <- "SCI_v1_1_memory_safe"
pseudobulk_version <- "PBVDR_v1_0"
upstream_model_version <- "VDRMODELS_v1_1"
DE_version <- "DE_v1_2"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1H_",
  dataset_id,
  "_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(
  path,
  max_chars = 245L
) {
  normalized <- normalize_slashes(path)
  path_n <- nchar(normalized, type = "chars")

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Output path is too long for safe Windows/R writing (",
      path_n,
      " characters; safe limit ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_windows_path_length(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(dirname(path))) {
    stop("Could not create output directory:\n", dirname(path))
  }

  writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  out <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

safe_package_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]

  if (!length(x)) {
    ""
  } else {
    paste(x, collapse = sep)
  }
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

  toupper(trimws(as.character(x))) %in%
    c("TRUE", "T", "1", "YES", "Y")
}

slugify <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  substr(out, 1L, 60L)
}

cell_type_code <- function(cell_type) {
  code_map <- c(
    "Ciliated epithelial cells" = "CE",
    "Endothelial cells" = "ENDO",
    "Myeloid cells" = "MY",
    "NK/T cells" = "NKT",
    "Stromal cells" = "STR",
    "Unciliated epithelial cells" = "UE"
  )

  code <- unname(code_map[as.character(cell_type)])

  if (length(code) != 1L || is.na(code) || !nzchar(code)) {
    stop("No short cell-type code defined for: ", cell_type)
  }

  code
}

contrast_code <- function(contrast_name) {
  code_map <- c(
    "LH5_vs_LH3" = "5v3",
    "LH7_vs_LH5" = "7v5",
    "LH9_vs_LH7" = "9v7",
    "LH11_vs_LH9" = "11v9",
    "LH7_vs_LH3" = "7v3",
    "LH11_vs_LH3" = "11v3",
    "LH11_vs_LH7" = "11v7",
    "RIF_vs_LH7_STANDARD" = "STD",
    "RIF_vs_LH7_PGT" = "PGT"
  )

  code <- unname(code_map[as.character(contrast_name)])

  if (length(code) != 1L || is.na(code) || !nzchar(code)) {
    stop("No short contrast code defined for: ", contrast_name)
  }

  code
}

sign_label <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(
      x > 0,
      "UP",
      ifelse(x < 0, "DOWN", "ZERO")
    )
  )
}

save_plot_pair <- function(
  plot_object,
  stem,
  width,
  height
) {
  png_path <- file.path(
    figure_root,
    paste0(stem, ".png")
  )

  pdf_path <- file.path(
    figure_root,
    paste0(stem, ".pdf")
  )

  assert_windows_path_length(png_path)
  assert_windows_path_length(pdf_path)

  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  c(png_path, pdf_path)
}

log_path <- NULL

log_event <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste(..., collapse = " ")
  )

  cat(line, "\n")

  if (!is.null(log_path)) {
    cat(line, "\n", file = log_path, append = TRUE)
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Package preflight
# ----------------------------------------------------------------------

required_CRAN_packages <- c(
  "Matrix",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

required_Bioc_packages <- c(
  "edgeR",
  "limma"
)

missing_CRAN <- required_CRAN_packages[
  !vapply(
    required_CRAN_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

missing_Bioc <- required_Bioc_packages[
  !vapply(
    required_Bioc_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_CRAN) || length(missing_Bioc)) {
  message_lines <- c(
    "Required package(s) are missing."
  )

  if (length(missing_CRAN)) {
    message_lines <- c(
      message_lines,
      paste0(
        "Missing CRAN: ",
        paste(missing_CRAN, collapse = ", ")
      ),
      paste0(
        "Run: install.packages(c(",
        paste(
          paste0('"', missing_CRAN, '"'),
          collapse = ", "
        ),
        "))"
      )
    )
  }

  if (length(missing_Bioc)) {
    message_lines <- c(
      message_lines,
      paste0(
        "Missing Bioconductor: ",
        paste(missing_Bioc, collapse = ", ")
      ),
      "Run:",
      'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")',
      paste0(
        "BiocManager::install(c(",
        paste(
          paste0('"', missing_Bioc, '"'),
          collapse = ", "
        ),
        "), ask = FALSE, update = FALSE)"
      )
    )
  }

  stop(paste(message_lines, collapse = "\n"))
}


if (
  !exists(
    "fry",
    envir = asNamespace("limma"),
    inherits = FALSE
  ) ||
    !is.function(
      get(
        "fry",
        envir = asNamespace("limma"),
        inherits = FALSE
      )
    )
) {
  stop(
    "limma::fry is unavailable in the installed limma namespace. ",
    "Reinstall limma with BiocManager, restart R, and rerun this script."
  )
}


# ----------------------------------------------------------------------
# 3. Project paths
# ----------------------------------------------------------------------

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

execution_root <- file.path(
  step04_root,
  "05_execution"
)

results_root <- file.path(
  step04_root,
  "06_results"
)

qc_base_root <- file.path(
  step04_root,
  "08_qc"
)

parameter_root <- file.path(
  step04_root,
  "03_parameters"
)

environment_root <- file.path(
  step04_root,
  "04_environment"
)

manuscript_text_root <- file.path(
  step04_root,
  "09_manuscript_text"
)

pb_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "pb",
  dataset_id,
  pseudobulk_version
)

pb_counts_path <- file.path(
  pb_root,
  "objects",
  "pb_raw_counts.rds"
)

pb_status_path <- file.path(
  pb_root,
  "objects",
  "pb_status.csv"
)

pb_metadata_path <- file.path(
  pb_root,
  "tables",
  "02_pb_metadata.csv"
)

eligibility_path <- file.path(
  pb_root,
  "tables",
  "04_model_eligibility.csv"
)

candidate_path <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons",
  "06_results",
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  "RNA_OCCUPANCY_EVIDENCE_v1.0",
  "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv"
)

output_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "de",
  dataset_id,
  DE_version
)

previous_partial_root <- file.path(
  results_root,
  "S4E1_sc",
  input_version,
  "de",
  dataset_id,
  "GENEDE_v1_1"
)

previous_per_model_root <- file.path(
  previous_partial_root,
  "per_model"
)

per_model_root <- file.path(
  output_root,
  "m"
)

full_result_root <- file.path(
  output_root,
  "full"
)

table_root <- file.path(
  output_root,
  "tab"
)

candidate_root <- file.path(
  output_root,
  "cand"
)

source_data_root <- file.path(
  output_root,
  "src"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1H_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_deQC"
  )
)

figure_root <- file.path(
  run_qc_root,
  "figures"
)

for (p in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  output_root,
  per_model_root,
  full_result_root,
  table_root,
  candidate_root,
  source_data_root,
  run_qc_root,
  figure_root
)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
)

session_path <- file.path(
  environment_root,
  paste0("sessionInfo_", run_id, ".txt")
)

desktop_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
  file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

zip_parent <- if (length(desktop_candidates)) {
  desktop_candidates[1]
} else {
  dirname(project_root)
}

handoff_zip <- file.path(
  zip_parent,
  paste0(run_id, "_SMALL_REPORT.zip")
)

planned_paths <- data.frame(
  object = c(
    "output_root",
    "per_model_root",
    "full_result_root",
    "table_root",
    "candidate_root",
    "source_data_root",
    "run_qc_root",
    "figure_root",
    "log_path",
    "run_manifest_path",
    "parameter_path",
    "session_path",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    output_root,
    per_model_root,
    full_result_root,
    table_root,
    candidate_root,
    source_data_root,
    run_qc_root,
    figure_root,
    log_path,
    run_manifest_path,
    parameter_path,
    session_path,
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_paths$character_n <- nchar(
  planned_paths$path,
  type = "chars"
)

if (
  .Platform$OS.type == "windows" &&
    any(planned_paths$character_n > 245L)
) {
  bad_paths <- planned_paths[
    planned_paths$character_n > 245L,
    ,
    drop = FALSE
  ]

  stop(
    "One or more planned output paths are too long:\n",
    paste(
      paste0(
        bad_paths$character_n,
        " chars: ",
        bad_paths$path
      ),
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream Gate and inputs
# ----------------------------------------------------------------------

required_inputs <- c(
  pb_counts_path,
  pb_status_path,
  pb_metadata_path,
  eligibility_path,
  candidate_path
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_inputs, collapse = "\n")
  )
}

upstream_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

upstream_qc_dirs <- upstream_qc_dirs[
  grepl(
    "^STEP04E1G_[0-9]{8}_[0-9]{6}_modelQC$",
    basename(upstream_qc_dirs),
    perl = TRUE
  )
]

if (!length(upstream_qc_dirs)) {
  stop("No completed STEP04E1G model-QC folder found.")
}

latest_upstream_qc <- upstream_qc_dirs[
  which.max(file.info(upstream_qc_dirs)$mtime)
]

upstream_gate_path <- file.path(
  latest_upstream_qc,
  "VDR_MODEL_GATE.csv"
)

upstream_gate <- safe_read_csv(upstream_gate_path)
pb_status <- safe_read_csv(pb_status_path)
pb_metadata <- safe_read_csv(pb_metadata_path)
eligibility <- safe_read_csv(eligibility_path)
candidate_table <- safe_read_csv(candidate_path)

if (
  is.null(upstream_gate) ||
    nrow(upstream_gate) != 1L ||
    upstream_gate$gate_status[[1]] !=
      "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE" ||
    !isTRUE(
      upstream_gate$gene_level_pseudobulk_DE_authorized[[1]]
    )
) {
  stop("STEP04E1G Gate does not authorize gene-level pseudobulk DE.")
}

if (
  is.null(pb_status) ||
    nrow(pb_status) != 1L ||
    pb_status$status[[1]] != "COMPLETE" ||
    pb_status$gene_n[[1]] != 36601L ||
    pb_status$pseudobulk_n[[1]] != 222L ||
    pb_status$aggregated_cell_n[[1]] != 236636L
) {
  stop("Pseudobulk status failed identity validation.")
}

required_metadata_columns <- c(
  "pb_id",
  "sample_accession",
  "donor_id",
  "broad_cell_type",
  "clinical_group",
  "LH_day",
  "temporal_fertile_eligible",
  "LH7_RIF_comparison_eligible",
  "LH7_fertile_subgroup",
  "minimum_cell_threshold_pass",
  "cell_n",
  "library_size"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(pb_metadata)
)

if (length(missing_metadata_columns)) {
  stop(
    "Pseudobulk metadata is missing columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

required_candidate_columns <- c(
  "feature_id",
  "interaction_logFC",
  "occupancy_supported",
  "shared_occupancy_supported",
  "any_promoter_overlap_support"
)

missing_candidate_columns <- setdiff(
  required_candidate_columns,
  colnames(candidate_table)
)

if (length(missing_candidate_columns)) {
  stop(
    "Candidate table is missing columns: ",
    paste(missing_candidate_columns, collapse = ", ")
  )
}

if (
  nrow(pb_metadata) != 222L ||
    nrow(eligibility) != 8L ||
    nrow(candidate_table) != 436L
) {
  stop("Input table row-count identity checks failed.")
}

candidate_table$feature_id <- as.character(
  candidate_table$feature_id
)

candidate_table$interaction_logFC <- safe_numeric(
  candidate_table$interaction_logFC
)

candidate_table$occupancy_supported <- as_logical_safe(
  candidate_table$occupancy_supported
)

candidate_table$shared_occupancy_supported <- as_logical_safe(
  candidate_table$shared_occupancy_supported
)

candidate_table$any_promoter_overlap_support <- as_logical_safe(
  candidate_table$any_promoter_overlap_support
)

if (
  anyDuplicated(candidate_table$feature_id) > 0L ||
    sum(candidate_table$occupancy_supported) != 284L ||
    sum(candidate_table$shared_occupancy_supported) != 46L ||
    sum(candidate_table$any_promoter_overlap_support) != 84L ||
    any(
      !is.finite(candidate_table$interaction_logFC)
    ) ||
    any(candidate_table$interaction_logFC == 0)
) {
  stop("Locked VDR candidate hierarchy failed validation.")
}

pb_metadata$LH_day <- safe_numeric(pb_metadata$LH_day)

logical_columns <- c(
  "temporal_fertile_eligible",
  "LH7_RIF_comparison_eligible",
  "minimum_cell_threshold_pass"
)

for (column_name in logical_columns) {
  pb_metadata[[column_name]] <- as_logical_safe(
    pb_metadata[[column_name]]
  )
}

temporal_types <- as.character(
  eligibility$broad_cell_type[
    as_logical_safe(
      eligibility$fertile_temporal_model_eligible
    )
  ]
)

LH7_types <- as.character(
  eligibility$broad_cell_type[
    as_logical_safe(
      eligibility$LH7_RIF_vs_fertile_model_eligible
    )
  ]
)

if (
  length(temporal_types) != 5L ||
    length(LH7_types) != 6L
) {
  stop("Eligible cell-type counts do not match the upstream Gate.")
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  upstream_model_version = upstream_model_version,
  DE_version = DE_version,
  supersedes_partial_DE_versions = c(
    "GENEDE_v1_0_EDGE_R_NAMESPACE_FRY_FAILURE",
    "GENEDE_v1_1_DYNAMIC_PATH_LENGTH_FAILURE"
  ),
  previous_partial_root =
    normalize_slashes(previous_partial_root),
  resume_policy =
    "IMPORT_COMPLETE_PER_MODEL_RDS_FROM_GENEDE_v1_1",
  dataset_id = dataset_id,
  primary_gene_level_FDR =
    "GENOME_WIDE_BH_WITHIN_EACH_CELL_TYPE_AND_CONTRAST",
  focused_candidate_FDR =
    "SEPARATE_BH_WITHIN_PROSPECTIVELY_LOCKED_436_284_46_84_HIERARCHIES",
  reporting_effect_threshold_pairwise =
    "ABS_LOG2_FOLD_CHANGE_GE_0.5",
  normalization = "TMM",
  expression_filter = "edgeR_filterByExpr_DEFAULT_WITH_MODEL_DESIGN",
  model = "ROBUST_EDGE_R_QUASI_LIKELIHOOD_GLM",
  temporal_design =
    "CROSS_SECTIONAL_FACTOR_LH3_LH5_LH7_LH9_LH11_UNBALANCED_ALLOWED",
  temporal_omnibus =
    "MULTI_DF_LH5_LH7_LH9_LH11_VERSUS_LH3",
  LH7_primary =
    "RIF_VERSUS_ALL_FERTILE_LH7",
  LH7_sensitivity = c(
    "RIF_VERSUS_LH7_STANDARD",
    "RIF_VERSUS_LH7_PGT"
  ),
  locked_program_test =
    "DIRECTIONALLY_WEIGHTED_EDGE_R_FRY_FOR_ONE_DF_CONTRASTS",
  temporal_cell_types = temporal_types,
  LH7_cell_types = LH7_types,
  composition_inference_authorized = FALSE,
  subtype_inference_authorized = FALSE,
  random_seed = 20260805L,
  result_lock_authorized = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_csv_utf8(
  planned_paths,
  file.path(run_qc_root, "00_path_audit.csv")
)

log_event(
  "Starting",
  run_id,
  "temporal_types=",
  collapse_unique(temporal_types),
  "LH7_types=",
  collapse_unique(LH7_types)
)


# ----------------------------------------------------------------------
# 5. Load and align pseudobulk count matrix
# ----------------------------------------------------------------------

pb_counts <- readRDS(pb_counts_path)

if (
  !inherits(pb_counts, "Matrix") &&
    !is.matrix(pb_counts)
) {
  stop("pb_raw_counts.rds is not a matrix-like object.")
}

if (
  nrow(pb_counts) != 36601L ||
    ncol(pb_counts) != 222L ||
    anyDuplicated(rownames(pb_counts)) > 0L ||
    anyDuplicated(colnames(pb_counts)) > 0L
) {
  stop("Pseudobulk count matrix identity checks failed.")
}

pb_metadata <- pb_metadata[
  match(
    colnames(pb_counts),
    pb_metadata$pb_id
  ),
  ,
  drop = FALSE
]

if (
  anyNA(pb_metadata$pb_id) ||
    !identical(
      as.character(pb_metadata$pb_id),
      colnames(pb_counts)
    )
) {
  stop("Pseudobulk metadata did not align with count columns.")
}

candidate_coverage <- data.frame(
  candidate_hierarchy = c(
    "RNA_ALL_436",
    "OCCUPANCY_SUPPORTED_284",
    "SHARED_OCCUPANCY_46",
    "PROMOTER_OCCUPANCY_84"
  ),
  locked_gene_n = c(436L, 284L, 46L, 84L),
  present_in_counts_n = c(
    sum(
      candidate_table$feature_id %in%
        rownames(pb_counts)
    ),
    sum(
      candidate_table$feature_id[
        candidate_table$occupancy_supported
      ] %in% rownames(pb_counts)
    ),
    sum(
      candidate_table$feature_id[
        candidate_table$shared_occupancy_supported
      ] %in% rownames(pb_counts)
    ),
    sum(
      candidate_table$feature_id[
        candidate_table$any_promoter_overlap_support
      ] %in% rownames(pb_counts)
    )
  ),
  stringsAsFactors = FALSE
)

candidate_coverage$coverage_fraction <- (
  candidate_coverage$present_in_counts_n /
    candidate_coverage$locked_gene_n
)

write_csv_utf8(
  candidate_coverage,
  file.path(
    run_qc_root,
    "01_candidate_coverage.csv"
  )
)

if (any(candidate_coverage$coverage_fraction < 0.90)) {
  stop("Locked candidate coverage is below 90% in pseudobulk counts.")
}


# ----------------------------------------------------------------------
# 6. Result annotation and model helpers
# ----------------------------------------------------------------------

candidate_annotation <- candidate_table[
  ,
  c(
    "feature_id",
    "interaction_logFC",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support"
  ),
  drop = FALSE
]

colnames(candidate_annotation)[1] <- "gene"

candidate_annotation$RNA_ALL_436 <- TRUE

add_candidate_annotation <- function(result_table) {
  result_table$gene <- as.character(result_table$gene)

  out <- merge(
    result_table,
    candidate_annotation,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )

  out <- out[
    match(result_table$gene, out$gene),
    ,
    drop = FALSE
  ]

  logical_flags <- c(
    "RNA_ALL_436",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support"
  )

  for (flag in logical_flags) {
    out[[flag]][is.na(out[[flag]])] <- FALSE
  }

  out
}

add_candidate_restricted_FDR <- function(result_table) {
  result_table$FDR_locked436 <- NA_real_
  result_table$FDR_occupancy284 <- NA_real_
  result_table$FDR_shared46 <- NA_real_
  result_table$FDR_promoter84 <- NA_real_

  hierarchy_map <- list(
    FDR_locked436 = result_table$RNA_ALL_436,
    FDR_occupancy284 =
      result_table$occupancy_supported,
    FDR_shared46 =
      result_table$shared_occupancy_supported,
    FDR_promoter84 =
      result_table$any_promoter_overlap_support
  )

  for (column_name in names(hierarchy_map)) {
    index <- which(
      hierarchy_map[[column_name]] &
        is.finite(result_table$PValue)
    )

    if (length(index)) {
      result_table[[column_name]][index] <-
        stats::p.adjust(
          result_table$PValue[index],
          method = "BH"
        )
    }
  }

  result_table
}

finalize_one_df_result <- function(
  result_table,
  cell_type,
  model_family,
  contrast_name
) {
  result_table <- add_candidate_annotation(result_table)
  result_table <- add_candidate_restricted_FDR(result_table)

  result_table$broad_cell_type <- cell_type
  result_table$model_family <- model_family
  result_table$contrast <- contrast_name

  result_table$effect_direction <- sign_label(
    result_table$logFC
  )

  result_table$whole_genome_DE_status <- ifelse(
    result_table$FDR < 0.05 &
      abs(result_table$logFC) >= 0.5,
    ifelse(
      result_table$logFC > 0,
      "UP_FDR_LT_0.05_ABS_LOGFC_GE_0.5",
      "DOWN_FDR_LT_0.05_ABS_LOGFC_GE_0.5"
    ),
    "NOT_WHOLE_GENOME_DE"
  )

  result_table$interaction_direction_concordant <- ifelse(
    result_table$RNA_ALL_436 &
      is.finite(result_table$interaction_logFC),
    sign(result_table$logFC) ==
      sign(result_table$interaction_logFC),
    NA
  )

  first_columns <- c(
    "gene",
    "broad_cell_type",
    "model_family",
    "contrast",
    "logFC",
    "logCPM",
    "F",
    "PValue",
    "FDR",
    "whole_genome_DE_status",
    "effect_direction",
    "RNA_ALL_436",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support",
    "interaction_logFC",
    "interaction_direction_concordant",
    "FDR_locked436",
    "FDR_occupancy284",
    "FDR_shared46",
    "FDR_promoter84"
  )

  result_table[
    ,
    c(
      first_columns,
      setdiff(
        colnames(result_table),
        first_columns
      )
    ),
    drop = FALSE
  ]
}

finalize_omnibus_result <- function(
  result_table,
  cell_type
) {
  result_table <- add_candidate_annotation(result_table)
  result_table <- add_candidate_restricted_FDR(result_table)

  result_table$broad_cell_type <- cell_type
  result_table$model_family <- "FERTILE_TEMPORAL_OMNIBUS"
  result_table$contrast <- "ANY_LH_TIME_EFFECT"

  result_table$whole_genome_DE_status <- ifelse(
    result_table$FDR < 0.05,
    "TEMPORAL_OMNIBUS_FDR_LT_0.05",
    "NOT_WHOLE_GENOME_TEMPORAL_DE"
  )

  first_columns <- c(
    "gene",
    "broad_cell_type",
    "model_family",
    "contrast",
    "logCPM",
    "F",
    "PValue",
    "FDR",
    "whole_genome_DE_status",
    "RNA_ALL_436",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support",
    "interaction_logFC",
    "FDR_locked436",
    "FDR_occupancy284",
    "FDR_shared46",
    "FDR_promoter84"
  )

  result_table[
    ,
    c(
      first_columns[
        first_columns %in%
          colnames(result_table)
      ],
      setdiff(
        colnames(result_table),
        first_columns
      )
    ),
    drop = FALSE
  ]
}

top_table_from_test <- function(test_object) {
  table <- edgeR::topTags(
    test_object,
    n = Inf,
    sort.by = "PValue"
  )$table

  table$gene <- rownames(table)
  rownames(table) <- NULL

  if (!"FDR" %in% colnames(table)) {
    table$FDR <- stats::p.adjust(
      table$PValue,
      method = "BH"
    )
  }

  table
}

make_directional_sets <- function(filtered_genes) {
  set_definitions <- list(
    RNA_ALL_436 = rep(
      TRUE,
      nrow(candidate_table)
    ),
    OCCUPANCY_SUPPORTED_284 =
      candidate_table$occupancy_supported,
    SHARED_OCCUPANCY_46 =
      candidate_table$shared_occupancy_supported,
    PROMOTER_OCCUPANCY_84 =
      candidate_table$any_promoter_overlap_support
  )

  output <- list()

  for (set_name in names(set_definitions)) {
    one <- candidate_table[
      set_definitions[[set_name]],
      c(
        "feature_id",
        "interaction_logFC"
      ),
      drop = FALSE
    ]

    one$index <- match(
      one$feature_id,
      filtered_genes
    )

    one <- one[
      !is.na(one$index),
      ,
      drop = FALSE
    ]

    one$weight <- sign(one$interaction_logFC)

    output[[set_name]] <- data.frame(
      index = as.integer(one$index),
      weight = as.numeric(one$weight),
      stringsAsFactors = FALSE
    )
  }

  output
}

run_directional_fry <- function(
  y,
  design,
  contrast,
  cell_type,
  model_family,
  contrast_name
) {
  sets <- make_directional_sets(
    rownames(y$counts)
  )

  fry_result <- limma::fry(
    y = y,
    index = sets,
    design = design,
    contrast = contrast,
    sort = "none"
  )

  fry_result$program_name <- rownames(fry_result)
  rownames(fry_result) <- NULL

  fry_result$broad_cell_type <- cell_type
  fry_result$model_family <- model_family
  fry_result$contrast <- contrast_name

  first_columns <- c(
    "broad_cell_type",
    "model_family",
    "contrast",
    "program_name"
  )

  fry_result[
    ,
    c(
      first_columns,
      setdiff(
        colnames(fry_result),
        first_columns
      )
    ),
    drop = FALSE
  ]
}

fit_edgeR_QL <- function(
  count_matrix,
  metadata,
  group_factor
) {
  counts_dense <- as.matrix(count_matrix)

  y <- edgeR::DGEList(
    counts = counts_dense,
    samples = metadata,
    group = group_factor
  )

  design <- stats::model.matrix(
    ~ 0 + group_factor
  )

  colnames(design) <- levels(group_factor)
  rownames(design) <- colnames(y$counts)

  keep <- edgeR::filterByExpr(
    y,
    design = design
  )

  filtering_QC <- data.frame(
    input_gene_n = nrow(y),
    retained_gene_n = sum(keep),
    filtered_gene_n = sum(!keep),
    retained_fraction = mean(keep),
    stringsAsFactors = FALSE
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )

  y <- edgeR::estimateDisp(
    y,
    design = design,
    robust = TRUE
  )

  fit <- edgeR::glmQLFit(
    y,
    design = design,
    robust = TRUE
  )

  normalization_QC <- data.frame(
    pb_id = colnames(y$counts),
    group = as.character(group_factor),
    library_size = y$samples$lib.size,
    norm_factor = y$samples$norm.factors,
    effective_library_size = (
      y$samples$lib.size *
        y$samples$norm.factors
    ),
    stringsAsFactors = FALSE
  )

  log_CPM <- edgeR::cpm(
    y,
    log = TRUE,
    prior.count = 2
  )

  PCA_fit <- stats::prcomp(
    t(log_CPM),
    center = TRUE,
    scale. = FALSE
  )

  variance_fraction <- (
    PCA_fit$sdev^2 /
      sum(PCA_fit$sdev^2)
  )

  PCA_table <- data.frame(
    pb_id = rownames(PCA_fit$x),
    PC1 = PCA_fit$x[, 1],
    PC2 = PCA_fit$x[, 2],
    PC1_variance = variance_fraction[1],
    PC2_variance = variance_fraction[2],
    group = as.character(group_factor),
    stringsAsFactors = FALSE
  )

  list(
    y = y,
    design = design,
    fit = fit,
    filtering_QC = filtering_QC,
    normalization_QC = normalization_QC,
    PCA_table = PCA_table
  )
}

expected_temporal_counts <- function(
  eligibility_table,
  cell_type
) {
  one <- eligibility_table[
    eligibility_table$broad_cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]

  if (nrow(one) != 1L) {
    stop(
      "Eligibility row was not unique for ",
      cell_type
    )
  }

  counts <- c(
    LH3 = safe_numeric(one$fertile_LH3_n),
    LH5 = safe_numeric(one$fertile_LH5_n),
    LH7 = safe_numeric(one$fertile_LH7_n),
    LH9 = safe_numeric(one$fertile_LH9_n),
    LH11 = safe_numeric(one$fertile_LH11_n)
  )

  if (
    any(!is.finite(counts)) ||
      any(counts < 2L)
  ) {
    stop(
      "Invalid temporal eligibility counts for ",
      cell_type
    )
  }

  as.integer(counts)
}


# ----------------------------------------------------------------------
# 7. Contrast definitions
# ----------------------------------------------------------------------

temporal_levels <- c(
  "LH3",
  "LH5",
  "LH7",
  "LH9",
  "LH11"
)

planned_temporal_contrasts <- list(
  LH5_vs_LH3 = c(
    LH3 = -1,
    LH5 = 1,
    LH7 = 0,
    LH9 = 0,
    LH11 = 0
  ),
  LH7_vs_LH5 = c(
    LH3 = 0,
    LH5 = -1,
    LH7 = 1,
    LH9 = 0,
    LH11 = 0
  ),
  LH9_vs_LH7 = c(
    LH3 = 0,
    LH5 = 0,
    LH7 = -1,
    LH9 = 1,
    LH11 = 0
  ),
  LH11_vs_LH9 = c(
    LH3 = 0,
    LH5 = 0,
    LH7 = 0,
    LH9 = -1,
    LH11 = 1
  ),
  LH7_vs_LH3 = c(
    LH3 = -1,
    LH5 = 0,
    LH7 = 1,
    LH9 = 0,
    LH11 = 0
  ),
  LH11_vs_LH3 = c(
    LH3 = -1,
    LH5 = 0,
    LH7 = 0,
    LH9 = 0,
    LH11 = 1
  ),
  LH11_vs_LH7 = c(
    LH3 = 0,
    LH5 = 0,
    LH7 = -1,
    LH9 = 0,
    LH11 = 1
  )
)

temporal_omnibus_contrast <- cbind(
  LH5_vs_LH3 =
    planned_temporal_contrasts$LH5_vs_LH3,
  LH7_vs_LH3 =
    planned_temporal_contrasts$LH7_vs_LH3,
  LH9_vs_LH3 = c(
    LH3 = -1,
    LH5 = 0,
    LH7 = 0,
    LH9 = 1,
    LH11 = 0
  ),
  LH11_vs_LH3 =
    planned_temporal_contrasts$LH11_vs_LH3
)


dynamic_paths <- character()

for (cell_type in temporal_types) {
  code <- cell_type_code(cell_type)
  dynamic_paths <- c(
    dynamic_paths,
    file.path(per_model_root, paste0("T_", code, ".rds")),
    file.path(full_result_root, paste0("T_", code, "_OMNI.csv"))
  )

  for (contrast_name in names(planned_temporal_contrasts)) {
    dynamic_paths <- c(
      dynamic_paths,
      file.path(
        full_result_root,
        paste0("T_", code, "_", contrast_code(contrast_name), ".csv")
      )
    )
  }
}

for (cell_type in LH7_types) {
  code <- cell_type_code(cell_type)
  dynamic_paths <- c(
    dynamic_paths,
    file.path(per_model_root, paste0("L7_", code, ".rds")),
    file.path(full_result_root, paste0("L7_", code, "_PRI.csv")),
    file.path(full_result_root, paste0("L7_", code, "_STD.csv")),
    file.path(full_result_root, paste0("L7_", code, "_PGT.csv"))
  )
}

dynamic_path_audit <- data.frame(
  path = normalize_slashes(dynamic_paths),
  character_n = nchar(normalize_slashes(dynamic_paths), type = "chars"),
  stringsAsFactors = FALSE
)

if (.Platform$OS.type == "windows" && any(dynamic_path_audit$character_n > 245L)) {
  bad_dynamic_paths <- dynamic_path_audit[
    dynamic_path_audit$character_n > 245L,
    ,
    drop = FALSE
  ]

  stop(
    "One or more dynamic output paths remain too long:\n",
    paste(
      paste0(bad_dynamic_paths$character_n, " chars: ", bad_dynamic_paths$path),
      collapse = "\n"
    )
  )
}

write_csv_utf8(
  dynamic_path_audit,
  file.path(run_qc_root, "01_dynamic_path_audit.csv")
)

materialize_temporal_outputs <- function(model_object, cell_type) {
  code <- cell_type_code(cell_type)

  write_csv_utf8(
    model_object$omnibus,
    file.path(full_result_root, paste0("T_", code, "_OMNI.csv"))
  )

  for (contrast_name in names(model_object$pairwise)) {
    write_csv_utf8(
      model_object$pairwise[[contrast_name]],
      file.path(
        full_result_root,
        paste0("T_", code, "_", contrast_code(contrast_name), ".csv")
      )
    )
  }

  invisible(TRUE)
}

materialize_LH7_outputs <- function(model_object, cell_type) {
  code <- cell_type_code(cell_type)

  write_csv_utf8(
    model_object$primary,
    file.path(full_result_root, paste0("L7_", code, "_PRI.csv"))
  )

  for (subgroup in names(model_object$sensitivity)) {
    contrast_name <- paste0("RIF_vs_", subgroup)

    write_csv_utf8(
      model_object$sensitivity[[subgroup]],
      file.path(
        full_result_root,
        paste0("L7_", code, "_", contrast_code(contrast_name), ".csv")
      )
    )
  }

  invisible(TRUE)
}


# ----------------------------------------------------------------------
# 8. Temporal models: five broad cell types
# ----------------------------------------------------------------------

temporal_model_objects <- list()
temporal_status_rows <- list()

for (
  model_index in seq_along(temporal_types)
) {
  cell_type <- temporal_types[model_index]
  cell_code <- cell_type_code(cell_type)
  old_cell_slug <- slugify(cell_type)

  model_rds <- file.path(
    per_model_root,
    paste0("T_", cell_code, ".rds")
  )

  old_model_rds <- file.path(
    previous_per_model_root,
    paste0("temporal_", old_cell_slug, ".rds")
  )

  if (!file.exists(model_rds) && file.exists(old_model_rds)) {
    old_object <- tryCatch(readRDS(old_model_rds), error = function(e) NULL)

    if (!is.null(old_object) && identical(old_object$status, "COMPLETE")) {
      saveRDS(old_object, model_rds, compress = "gzip")
      log_event("IMPORTED completed temporal model", cell_type, "from GENEDE_v1_1")
    }

    rm(old_object)
    gc()
  }

  if (file.exists(model_rds)) {
    existing <- tryCatch(
      readRDS(model_rds),
      error = function(e) NULL
    )

    if (
      !is.null(existing) &&
        identical(existing$status, "COMPLETE")
    ) {
      materialize_temporal_outputs(existing, cell_type)

      temporal_model_objects[[cell_type]] <- existing

      temporal_status_rows[[length(temporal_status_rows) + 1L]] <-
        existing$status_table

      log_event(
        "SKIP completed temporal model",
        cell_type
      )

      next
    }
  }

  model_start <- Sys.time()

  metadata_one <- pb_metadata[
    pb_metadata$broad_cell_type == cell_type &
      pb_metadata$clinical_group == "FERTILE" &
      pb_metadata$temporal_fertile_eligible &
      pb_metadata$minimum_cell_threshold_pass,
    ,
    drop = FALSE
  ]

  metadata_one$time_group <- factor(
    paste0("LH", metadata_one$LH_day),
    levels = temporal_levels
  )

  observed_counts <- as.integer(
    table(metadata_one$time_group)
  )

  expected_counts <- expected_temporal_counts(
    eligibility,
    cell_type
  )

  if (
    !identical(
      observed_counts,
      as.integer(expected_counts)
    ) ||
      any(observed_counts < 2L) ||
      anyDuplicated(metadata_one$donor_id) > 0L
  ) {
    stop(
      "Temporal donor identity failed for ",
      cell_type,
      ". Observed ",
      paste(
        temporal_levels,
        observed_counts,
        sep = "=",
        collapse = ";"
      )
    )
  }

  column_index <- match(
    metadata_one$pb_id,
    colnames(pb_counts)
  )

  if (anyNA(column_index)) {
    stop(
      "Temporal count columns were not found for ",
      cell_type
    )
  }

  fit_result <- fit_edgeR_QL(
    count_matrix = pb_counts[
      ,
      column_index,
      drop = FALSE
    ],
    metadata = metadata_one,
    group_factor = metadata_one$time_group
  )

  omnibus_test <- edgeR::glmQLFTest(
    fit_result$fit,
    contrast = temporal_omnibus_contrast
  )

  omnibus_table <- top_table_from_test(
    omnibus_test
  )

  omnibus_table <- finalize_omnibus_result(
    omnibus_table,
    cell_type
  )

  omnibus_path <- file.path(
    full_result_root,
    paste0("T_", cell_code, "_OMNI.csv")
  )

  write_csv_utf8(
    omnibus_table,
    omnibus_path
  )

  pairwise_tables <- list()
  fry_tables <- list()

  for (
    contrast_name in names(
      planned_temporal_contrasts
    )
  ) {
    contrast_vector <-
      planned_temporal_contrasts[[contrast_name]]

    pairwise_test <- edgeR::glmQLFTest(
      fit_result$fit,
      contrast = contrast_vector
    )

    pairwise_table <- top_table_from_test(
      pairwise_test
    )

    pairwise_table <- finalize_one_df_result(
      pairwise_table,
      cell_type = cell_type,
      model_family =
        "FERTILE_TEMPORAL_PAIRWISE",
      contrast_name = contrast_name
    )

    pairwise_path <- file.path(
      full_result_root,
      paste0("T_", cell_code, "_", contrast_code(contrast_name), ".csv")
    )

    write_csv_utf8(
      pairwise_table,
      pairwise_path
    )

    pairwise_tables[[contrast_name]] <-
      pairwise_table

    fry_tables[[contrast_name]] <-
      run_directional_fry(
        y = fit_result$y,
        design = fit_result$design,
        contrast = contrast_vector,
        cell_type = cell_type,
        model_family =
          "FERTILE_TEMPORAL_PAIRWISE",
        contrast_name = contrast_name
      )
  }

  filtering_QC <- fit_result$filtering_QC
  filtering_QC$broad_cell_type <- cell_type
  filtering_QC$model_family <-
    "FERTILE_TEMPORAL"

  normalization_QC <- fit_result$normalization_QC
  normalization_QC$broad_cell_type <- cell_type
  normalization_QC$model_family <-
    "FERTILE_TEMPORAL"

  PCA_table <- fit_result$PCA_table
  PCA_table$broad_cell_type <- cell_type
  PCA_table$model_family <-
    "FERTILE_TEMPORAL"
  PCA_table$model_label <- paste0(
    "Temporal | ",
    cell_type
  )

  elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        model_start,
        units = "mins"
      )
    ),
    2
  )

  status_table <- data.frame(
    model_family = "FERTILE_TEMPORAL",
    broad_cell_type = cell_type,
    status = "COMPLETE",
    donor_n = nrow(metadata_one),
    minimum_group_n = min(observed_counts),
    filtered_gene_n =
      filtering_QC$retained_gene_n,
    omnibus_genome_FDR_lt_0_05_n = sum(
      omnibus_table$FDR < 0.05
    ),
    elapsed_minutes = elapsed_minutes,
    stringsAsFactors = FALSE
  )

  output_object <- list(
    status = "COMPLETE",
    cell_type = cell_type,
    metadata = metadata_one,
    filtering_QC = filtering_QC,
    normalization_QC = normalization_QC,
    PCA_table = PCA_table,
    omnibus = omnibus_table,
    pairwise = pairwise_tables,
    fry = fry_tables,
    status_table = status_table
  )

  saveRDS(
    output_object,
    model_rds,
    compress = "gzip"
  )

  temporal_model_objects[[cell_type]] <-
    output_object

  temporal_status_rows[[length(temporal_status_rows) + 1L]] <-
    status_table

  log_event(
    "DONE temporal gene model",
    model_index,
    "/",
    length(temporal_types),
    cell_type,
    "retained_genes=",
    filtering_QC$retained_gene_n,
    "omnibus_FDR_lt_0.05=",
    status_table$omnibus_genome_FDR_lt_0_05_n,
    "elapsed_min=",
    elapsed_minutes
  )

  rm(
    fit_result,
    omnibus_test,
    omnibus_table,
    pairwise_tables,
    fry_tables,
    output_object
  )
  gc()
}


# ----------------------------------------------------------------------
# 9. LH+7 models: primary and subgroup sensitivities
# ----------------------------------------------------------------------

LH7_model_objects <- list()
LH7_status_rows <- list()

for (
  model_index in seq_along(LH7_types)
) {
  cell_type <- LH7_types[model_index]
  cell_code <- cell_type_code(cell_type)
  old_cell_slug <- slugify(cell_type)

  model_rds <- file.path(
    per_model_root,
    paste0("L7_", cell_code, ".rds")
  )

  old_model_rds <- file.path(
    previous_per_model_root,
    paste0("LH7_", old_cell_slug, ".rds")
  )

  if (!file.exists(model_rds) && file.exists(old_model_rds)) {
    old_object <- tryCatch(readRDS(old_model_rds), error = function(e) NULL)

    if (!is.null(old_object) && identical(old_object$status, "COMPLETE")) {
      saveRDS(old_object, model_rds, compress = "gzip")
      log_event("IMPORTED completed LH7 model", cell_type, "from GENEDE_v1_1")
    }

    rm(old_object)
    gc()
  }

  if (file.exists(model_rds)) {
    existing <- tryCatch(
      readRDS(model_rds),
      error = function(e) NULL
    )

    if (
      !is.null(existing) &&
        identical(existing$status, "COMPLETE")
    ) {
      materialize_LH7_outputs(existing, cell_type)

      LH7_model_objects[[cell_type]] <- existing

      LH7_status_rows[[length(LH7_status_rows) + 1L]] <-
        existing$status_table

      log_event(
        "SKIP completed LH7 model",
        cell_type
      )

      next
    }
  }

  model_start <- Sys.time()

  metadata_primary <- pb_metadata[
    pb_metadata$broad_cell_type == cell_type &
      pb_metadata$LH7_RIF_comparison_eligible &
      pb_metadata$minimum_cell_threshold_pass,
    ,
    drop = FALSE
  ]

  metadata_primary$group <- factor(
    metadata_primary$clinical_group,
    levels = c("FERTILE", "RIF")
  )

  primary_group_counts <- as.integer(
    table(metadata_primary$group)
  )

  if (
    any(primary_group_counts < 3L) ||
      anyDuplicated(metadata_primary$donor_id) > 0L
  ) {
    stop(
      "LH+7 donor identity failed for ",
      cell_type
    )
  }

  primary_column_index <- match(
    metadata_primary$pb_id,
    colnames(pb_counts)
  )

  if (anyNA(primary_column_index)) {
    stop(
      "LH+7 count columns were not found for ",
      cell_type
    )
  }

  primary_fit <- fit_edgeR_QL(
    count_matrix = pb_counts[
      ,
      primary_column_index,
      drop = FALSE
    ],
    metadata = metadata_primary,
    group_factor = metadata_primary$group
  )

  primary_contrast <- c(
    FERTILE = -1,
    RIF = 1
  )

  primary_test <- edgeR::glmQLFTest(
    primary_fit$fit,
    contrast = primary_contrast
  )

  primary_table <- top_table_from_test(
    primary_test
  )

  primary_table <- finalize_one_df_result(
    primary_table,
    cell_type = cell_type,
    model_family =
      "LH7_RIF_VS_ALL_FERTILE",
    contrast_name =
      "RIF_vs_all_fertile_LH7"
  )

  primary_path <- file.path(
    full_result_root,
    paste0("L7_", cell_code, "_PRI.csv")
  )

  write_csv_utf8(
    primary_table,
    primary_path
  )

  primary_fry <- run_directional_fry(
    y = primary_fit$y,
    design = primary_fit$design,
    contrast = primary_contrast,
    cell_type = cell_type,
    model_family =
      "LH7_RIF_VS_ALL_FERTILE",
    contrast_name =
      "RIF_vs_all_fertile_LH7"
  )

  sensitivity_tables <- list()
  sensitivity_fry <- list()
  sensitivity_filtering <- list()
  sensitivity_normalization <- list()
  sensitivity_PCA <- list()

  for (
    subgroup in c(
      "LH7_STANDARD",
      "LH7_PGT"
    )
  ) {
    metadata_sensitivity <- metadata_primary[
      metadata_primary$clinical_group == "RIF" |
        (
          metadata_primary$clinical_group ==
            "FERTILE" &
            metadata_primary$LH7_fertile_subgroup ==
              subgroup
        ),
      ,
      drop = FALSE
    ]

    metadata_sensitivity$group <- factor(
      ifelse(
        metadata_sensitivity$clinical_group ==
          "RIF",
        "RIF",
        subgroup
      ),
      levels = c(
        subgroup,
        "RIF"
      )
    )

    sensitivity_counts <- as.integer(
      table(metadata_sensitivity$group)
    )

    if (
      any(sensitivity_counts < 2L) ||
        anyDuplicated(
          metadata_sensitivity$donor_id
        ) > 0L
    ) {
      stop(
        "LH+7 sensitivity donor identity failed for ",
        cell_type,
        " / ",
        subgroup
      )
    }

    sensitivity_columns <- match(
      metadata_sensitivity$pb_id,
      colnames(pb_counts)
    )

    sensitivity_fit <- fit_edgeR_QL(
      count_matrix = pb_counts[
        ,
        sensitivity_columns,
        drop = FALSE
      ],
      metadata = metadata_sensitivity,
      group_factor =
        metadata_sensitivity$group
    )

    sensitivity_contrast <- c(-1, 1)
    names(sensitivity_contrast) <-
      levels(metadata_sensitivity$group)

    sensitivity_test <- edgeR::glmQLFTest(
      sensitivity_fit$fit,
      contrast = sensitivity_contrast
    )

    sensitivity_table <- top_table_from_test(
      sensitivity_test
    )

    sensitivity_table <- finalize_one_df_result(
      sensitivity_table,
      cell_type = cell_type,
      model_family =
        "LH7_RIF_SENSITIVITY",
      contrast_name = paste0(
        "RIF_vs_",
        subgroup
      )
    )

    sensitivity_path <- file.path(
      full_result_root,
      paste0(
        "L7_",
        cell_code,
        "_",
        contrast_code(paste0("RIF_vs_", subgroup)),
        ".csv"
      )
    )

    write_csv_utf8(
      sensitivity_table,
      sensitivity_path
    )

    sensitivity_tables[[subgroup]] <-
      sensitivity_table

    sensitivity_fry[[subgroup]] <-
      run_directional_fry(
        y = sensitivity_fit$y,
        design = sensitivity_fit$design,
        contrast = sensitivity_contrast,
        cell_type = cell_type,
        model_family =
          "LH7_RIF_SENSITIVITY",
        contrast_name = paste0(
          "RIF_vs_",
          subgroup
        )
      )

    one_filtering <-
      sensitivity_fit$filtering_QC

    one_filtering$broad_cell_type <-
      cell_type

    one_filtering$model_family <-
      "LH7_RIF_SENSITIVITY"

    one_filtering$contrast <- paste0(
      "RIF_vs_",
      subgroup
    )

    sensitivity_filtering[[subgroup]] <-
      one_filtering

    one_normalization <-
      sensitivity_fit$normalization_QC

    one_normalization$broad_cell_type <-
      cell_type

    one_normalization$model_family <-
      "LH7_RIF_SENSITIVITY"

    one_normalization$contrast <- paste0(
      "RIF_vs_",
      subgroup
    )

    sensitivity_normalization[[subgroup]] <-
      one_normalization

    one_PCA <- sensitivity_fit$PCA_table

    one_PCA$broad_cell_type <- cell_type
    one_PCA$model_family <-
      "LH7_RIF_SENSITIVITY"
    one_PCA$model_label <- paste0(
      "LH7 ",
      subgroup,
      " | ",
      cell_type
    )

    sensitivity_PCA[[subgroup]] <- one_PCA

    rm(
      sensitivity_fit,
      sensitivity_test,
      sensitivity_table
    )
    gc()
  }

  primary_filtering <- primary_fit$filtering_QC
  primary_filtering$broad_cell_type <- cell_type
  primary_filtering$model_family <-
    "LH7_RIF_VS_ALL_FERTILE"
  primary_filtering$contrast <-
    "RIF_vs_all_fertile_LH7"

  primary_normalization <-
    primary_fit$normalization_QC
  primary_normalization$broad_cell_type <-
    cell_type
  primary_normalization$model_family <-
    "LH7_RIF_VS_ALL_FERTILE"
  primary_normalization$contrast <-
    "RIF_vs_all_fertile_LH7"

  primary_PCA <- primary_fit$PCA_table
  primary_PCA$broad_cell_type <- cell_type
  primary_PCA$model_family <-
    "LH7_RIF_VS_ALL_FERTILE"
  primary_PCA$model_label <- paste0(
    "LH7 primary | ",
    cell_type
  )

  elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        model_start,
        units = "mins"
      )
    ),
    2
  )

  status_table <- data.frame(
    model_family =
      "LH7_RIF_VS_ALL_FERTILE",
    broad_cell_type = cell_type,
    status = "COMPLETE",
    fertile_n = primary_group_counts[1],
    RIF_n = primary_group_counts[2],
    filtered_gene_n =
      primary_filtering$retained_gene_n,
    primary_genome_DE_n = sum(
      primary_table$FDR < 0.05 &
        abs(primary_table$logFC) >= 0.5
    ),
    elapsed_minutes = elapsed_minutes,
    stringsAsFactors = FALSE
  )

  output_object <- list(
    status = "COMPLETE",
    cell_type = cell_type,
    metadata_primary = metadata_primary,
    primary_filtering = primary_filtering,
    primary_normalization =
      primary_normalization,
    primary_PCA = primary_PCA,
    primary = primary_table,
    primary_fry = primary_fry,
    sensitivity = sensitivity_tables,
    sensitivity_fry = sensitivity_fry,
    sensitivity_filtering =
      sensitivity_filtering,
    sensitivity_normalization =
      sensitivity_normalization,
    sensitivity_PCA = sensitivity_PCA,
    status_table = status_table
  )

  saveRDS(
    output_object,
    model_rds,
    compress = "gzip"
  )

  LH7_model_objects[[cell_type]] <-
    output_object

  LH7_status_rows[[length(LH7_status_rows) + 1L]] <-
    status_table

  log_event(
    "DONE LH7 gene model",
    model_index,
    "/",
    length(LH7_types),
    cell_type,
    "retained_genes=",
    primary_filtering$retained_gene_n,
    "primary_DE_n=",
    status_table$primary_genome_DE_n,
    "elapsed_min=",
    elapsed_minutes
  )

  rm(
    primary_fit,
    primary_test,
    primary_table,
    primary_fry,
    sensitivity_tables,
    sensitivity_fry,
    output_object
  )
  gc()
}


# ----------------------------------------------------------------------
# 10. Combine temporal outputs
# ----------------------------------------------------------------------

temporal_status <- do.call(
  rbind,
  temporal_status_rows
)

LH7_status <- do.call(
  rbind,
  LH7_status_rows
)

write_csv_utf8(
  temporal_status,
  file.path(
    table_root,
    "01_temporal_model_status.csv"
  )
)

write_csv_utf8(
  LH7_status,
  file.path(
    table_root,
    "02_LH7_model_status.csv"
  )
)

temporal_omnibus_all <- data.table::rbindlist(
  lapply(
    temporal_model_objects,
    function(x) x$omnibus
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_pairwise_all <- data.table::rbindlist(
  unlist(
    lapply(
      temporal_model_objects,
      function(x) x$pairwise
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_fry_all <- data.table::rbindlist(
  unlist(
    lapply(
      temporal_model_objects,
      function(x) x$fry
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_filtering_all <- data.table::rbindlist(
  lapply(
    temporal_model_objects,
    function(x) x$filtering_QC
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_normalization_all <- data.table::rbindlist(
  lapply(
    temporal_model_objects,
    function(x) x$normalization_QC
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_PCA_all <- data.table::rbindlist(
  lapply(
    temporal_model_objects,
    function(x) x$PCA_table
  ),
  use.names = TRUE,
  fill = TRUE
)

temporal_omnibus_significant <-
  temporal_omnibus_all[
    temporal_omnibus_all$FDR < 0.05,
    ,
    drop = FALSE
  ]

temporal_pairwise_significant <-
  temporal_pairwise_all[
    temporal_pairwise_all$FDR < 0.05 &
      abs(temporal_pairwise_all$logFC) >= 0.5,
    ,
    drop = FALSE
  ]

write_csv_utf8(
  as.data.frame(
    temporal_omnibus_significant
  ),
  file.path(
    table_root,
    "03_temporal_omnibus_significant.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    temporal_pairwise_significant
  ),
  file.path(
    table_root,
    "04_temporal_pairwise_significant.csv"
  )
)


# ----------------------------------------------------------------------
# 11. Combine LH+7 outputs
# ----------------------------------------------------------------------

LH7_primary_all <- data.table::rbindlist(
  lapply(
    LH7_model_objects,
    function(x) x$primary
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_sensitivity_all <- data.table::rbindlist(
  unlist(
    lapply(
      LH7_model_objects,
      function(x) x$sensitivity
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_fry_primary_all <- data.table::rbindlist(
  lapply(
    LH7_model_objects,
    function(x) x$primary_fry
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_fry_sensitivity_all <- data.table::rbindlist(
  unlist(
    lapply(
      LH7_model_objects,
      function(x) x$sensitivity_fry
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_filtering_primary_all <- data.table::rbindlist(
  lapply(
    LH7_model_objects,
    function(x) x$primary_filtering
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_normalization_primary_all <- data.table::rbindlist(
  lapply(
    LH7_model_objects,
    function(x) x$primary_normalization
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_PCA_primary_all <- data.table::rbindlist(
  lapply(
    LH7_model_objects,
    function(x) x$primary_PCA
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_filtering_sensitivity_all <- data.table::rbindlist(
  unlist(
    lapply(
      LH7_model_objects,
      function(x) x$sensitivity_filtering
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_normalization_sensitivity_all <- data.table::rbindlist(
  unlist(
    lapply(
      LH7_model_objects,
      function(x) x$sensitivity_normalization
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_PCA_sensitivity_all <- data.table::rbindlist(
  unlist(
    lapply(
      LH7_model_objects,
      function(x) x$sensitivity_PCA
    ),
    recursive = FALSE
  ),
  use.names = TRUE,
  fill = TRUE
)

LH7_primary_significant <- LH7_primary_all[
  LH7_primary_all$FDR < 0.05 &
    abs(LH7_primary_all$logFC) >= 0.5,
  ,
  drop = FALSE
]

LH7_sensitivity_significant <-
  LH7_sensitivity_all[
    LH7_sensitivity_all$FDR < 0.05 &
      abs(LH7_sensitivity_all$logFC) >= 0.5,
    ,
    drop = FALSE
  ]

write_csv_utf8(
  as.data.frame(LH7_primary_significant),
  file.path(
    table_root,
    "05_LH7_primary_significant.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    LH7_sensitivity_significant
  ),
  file.path(
    table_root,
    "06_LH7_sensitivity_significant.csv"
  )
)


# ----------------------------------------------------------------------
# 12. Candidate-focused and VDR-specific result tables
# ----------------------------------------------------------------------

candidate_temporal_omnibus <-
  temporal_omnibus_all[
    temporal_omnibus_all$RNA_ALL_436,
    ,
    drop = FALSE
  ]

candidate_temporal_pairwise <-
  temporal_pairwise_all[
    temporal_pairwise_all$RNA_ALL_436,
    ,
    drop = FALSE
  ]

candidate_LH7_primary <- LH7_primary_all[
  LH7_primary_all$RNA_ALL_436,
  ,
  drop = FALSE
]

candidate_LH7_sensitivity <-
  LH7_sensitivity_all[
    LH7_sensitivity_all$RNA_ALL_436,
    ,
    drop = FALSE
  ]

write_csv_utf8(
  as.data.frame(
    candidate_temporal_omnibus
  ),
  file.path(
    candidate_root,
    "01_locked436_temporal_omnibus.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    candidate_temporal_pairwise
  ),
  file.path(
    candidate_root,
    "02_locked436_temporal_pairwise.csv"
  )
)

write_csv_utf8(
  as.data.frame(candidate_LH7_primary),
  file.path(
    candidate_root,
    "03_locked436_LH7_primary.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    candidate_LH7_sensitivity
  ),
  file.path(
    candidate_root,
    "04_locked436_LH7_sensitivity.csv"
  )
)

all_result_tables <- list(
  temporal_omnibus =
    as.data.frame(temporal_omnibus_all),
  temporal_pairwise =
    as.data.frame(temporal_pairwise_all),
  LH7_primary =
    as.data.frame(LH7_primary_all),
  LH7_sensitivity =
    as.data.frame(LH7_sensitivity_all)
)

VDR_rows <- data.table::rbindlist(
  lapply(
    names(all_result_tables),
    function(table_name) {
      one <- all_result_tables[[table_name]]

      one <- one[
        one$gene == "VDR",
        ,
        drop = FALSE
      ]

      if (nrow(one)) {
        one$result_family <- table_name
      }

      one
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

write_csv_utf8(
  as.data.frame(VDR_rows),
  file.path(
    candidate_root,
    "05_VDR_gene_all_models.csv"
  )
)

fry_all <- data.table::rbindlist(
  list(
    temporal_fry_all,
    LH7_fry_primary_all,
    LH7_fry_sensitivity_all
  ),
  use.names = TRUE,
  fill = TRUE
)

write_csv_utf8(
  as.data.frame(fry_all),
  file.path(
    candidate_root,
    "06_directional_VDR_program_fry.csv"
  )
)


# ----------------------------------------------------------------------
# 13. Summaries and top-result handoff tables
# ----------------------------------------------------------------------

temporal_omnibus_summary <- data.table::as.data.table(
  temporal_omnibus_all
)[
  ,
  .(
    retained_gene_n = .N,
    whole_genome_FDR_lt_0_05_n =
      sum(FDR < 0.05),
    locked436_FDR_lt_0_05_n =
      sum(
        RNA_ALL_436 &
          FDR_locked436 < 0.05,
        na.rm = TRUE
      ),
    occupancy284_FDR_lt_0_05_n =
      sum(
        occupancy_supported &
          FDR_occupancy284 < 0.05,
        na.rm = TRUE
      ),
    shared46_FDR_lt_0_05_n =
      sum(
        shared_occupancy_supported &
          FDR_shared46 < 0.05,
        na.rm = TRUE
      ),
    promoter84_FDR_lt_0_05_n =
      sum(
        any_promoter_overlap_support &
          FDR_promoter84 < 0.05,
        na.rm = TRUE
      )
  ),
  by = broad_cell_type
]

temporal_pairwise_summary <- data.table::as.data.table(
  temporal_pairwise_all
)[
  ,
  .(
    whole_genome_up_n = sum(
      FDR < 0.05 &
        logFC >= 0.5
    ),
    whole_genome_down_n = sum(
      FDR < 0.05 &
        logFC <= -0.5
    ),
    locked436_focused_FDR_lt_0_05_n =
      sum(
        RNA_ALL_436 &
          FDR_locked436 < 0.05,
        na.rm = TRUE
      ),
    locked436_concordant_focused_n =
      sum(
        RNA_ALL_436 &
          FDR_locked436 < 0.05 &
          interaction_direction_concordant,
        na.rm = TRUE
      )
  ),
  by = .(
    broad_cell_type,
    contrast
  )
]

LH7_primary_summary <- data.table::as.data.table(
  LH7_primary_all
)[
  ,
  .(
    whole_genome_up_n = sum(
      FDR < 0.05 &
        logFC >= 0.5
    ),
    whole_genome_down_n = sum(
      FDR < 0.05 &
        logFC <= -0.5
    ),
    locked436_focused_FDR_lt_0_05_n =
      sum(
        RNA_ALL_436 &
          FDR_locked436 < 0.05,
        na.rm = TRUE
      ),
    locked436_concordant_focused_n =
      sum(
        RNA_ALL_436 &
          FDR_locked436 < 0.05 &
          interaction_direction_concordant,
        na.rm = TRUE
      )
  ),
  by = broad_cell_type
]

write_csv_utf8(
  as.data.frame(
    temporal_omnibus_summary
  ),
  file.path(
    table_root,
    "07_temporal_omnibus_summary.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    temporal_pairwise_summary
  ),
  file.path(
    table_root,
    "08_temporal_pairwise_summary.csv"
  )
)

write_csv_utf8(
  as.data.frame(LH7_primary_summary),
  file.path(
    table_root,
    "09_LH7_primary_summary.csv"
  )
)

top_temporal_omnibus <- data.table::as.data.table(
  temporal_omnibus_all
)[
  order(FDR, -F),
  head(.SD, 100L),
  by = broad_cell_type
]

top_temporal_pairwise <- data.table::as.data.table(
  temporal_pairwise_all
)[
  order(FDR, -abs(logFC)),
  head(.SD, 50L),
  by = .(
    broad_cell_type,
    contrast
  )
]

top_LH7_primary <- data.table::as.data.table(
  LH7_primary_all
)[
  order(FDR, -abs(logFC)),
  head(.SD, 100L),
  by = broad_cell_type
]

write_csv_utf8(
  as.data.frame(top_temporal_omnibus),
  file.path(
    table_root,
    "10_top_temporal_omnibus.csv"
  )
)

write_csv_utf8(
  as.data.frame(top_temporal_pairwise),
  file.path(
    table_root,
    "11_top_temporal_pairwise.csv"
  )
)

write_csv_utf8(
  as.data.frame(top_LH7_primary),
  file.path(
    table_root,
    "12_top_LH7_primary.csv"
  )
)


# ----------------------------------------------------------------------
# 14. Model QC tables
# ----------------------------------------------------------------------

filtering_all <- data.table::rbindlist(
  list(
    temporal_filtering_all,
    LH7_filtering_primary_all,
    LH7_filtering_sensitivity_all
  ),
  use.names = TRUE,
  fill = TRUE
)

normalization_all <- data.table::rbindlist(
  list(
    temporal_normalization_all,
    LH7_normalization_primary_all,
    LH7_normalization_sensitivity_all
  ),
  use.names = TRUE,
  fill = TRUE
)

PCA_all <- data.table::rbindlist(
  list(
    temporal_PCA_all,
    LH7_PCA_primary_all,
    LH7_PCA_sensitivity_all
  ),
  use.names = TRUE,
  fill = TRUE
)

write_csv_utf8(
  as.data.frame(filtering_all),
  file.path(
    run_qc_root,
    "02_filtering_QC.csv"
  )
)

write_csv_utf8(
  as.data.frame(normalization_all),
  file.path(
    run_qc_root,
    "03_normalization_QC.csv"
  )
)

write_csv_utf8(
  as.data.frame(PCA_all),
  file.path(
    source_data_root,
    "source_model_PCA.csv"
  )
)


# ----------------------------------------------------------------------
# 15. QC figures
# ----------------------------------------------------------------------

theme_qc <- ggplot2::theme_bw(
  base_size = 9
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank(),
    strip.background =
      ggplot2::element_rect(fill = "grey95")
  )

PCA_plot_data <- as.data.frame(PCA_all)

plot_PCA <- ggplot2::ggplot(
  PCA_plot_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    shape = group
  )
) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::facet_wrap(
    ~ model_label,
    scales = "free",
    ncol = 3
  ) +
  ggplot2::labs(
    title = "Gene-level pseudobulk model PCA",
    subtitle = "PCA calculated after cell-type-specific expression filtering and TMM normalization",
    x = "PC1",
    y = "PC2",
    shape = "Model group"
  ) +
  theme_qc +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      size = 6
    )
  )

save_plot_pair(
  plot_PCA,
  "QC1_model_PCA",
  width = 13,
  height = 12
)

temporal_count_plot <- as.data.frame(
  temporal_omnibus_summary
)

temporal_count_plot <- reshape(
  temporal_count_plot[
    ,
    c(
      "broad_cell_type",
      "whole_genome_FDR_lt_0_05_n",
      "locked436_FDR_lt_0_05_n"
    )
  ],
  varying = c(
    "whole_genome_FDR_lt_0_05_n",
    "locked436_FDR_lt_0_05_n"
  ),
  v.names = "gene_n",
  timevar = "result_scope",
  times = c(
    "Genome-wide FDR",
    "Locked 436 focused FDR"
  ),
  direction = "long"
)

plot_temporal_counts <- ggplot2::ggplot(
  temporal_count_plot,
  ggplot2::aes(
    x = broad_cell_type,
    y = gene_n,
    fill = result_scope
  )
) +
  ggplot2::geom_col(
    position = "dodge"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Temporal omnibus gene-level findings",
    subtitle = "Candidate-restricted FDR is secondary and does not replace genome-wide FDR",
    x = "Broad cell type",
    y = "Gene count",
    fill = "Multiplicity scope"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  )

save_plot_pair(
  plot_temporal_counts,
  "QC2_temporal_DE_counts",
  width = 10,
  height = 6
)

fry_plot_data <- as.data.frame(fry_all)

direction_column <- if (
  "Direction" %in% colnames(fry_plot_data)
) {
  "Direction"
} else {
  ""
}

p_column <- if (
  "FDR" %in% colnames(fry_plot_data)
) {
  "FDR"
} else if (
  "PValue" %in% colnames(fry_plot_data)
) {
  "PValue"
} else {
  ""
}

if (nzchar(p_column)) {
  fry_plot_data$plot_FDR <- safe_numeric(
    fry_plot_data[[p_column]]
  )

  fry_plot_data$minus_log10_FDR <- -log10(
    pmax(
      fry_plot_data$plot_FDR,
      .Machine$double.xmin
    )
  )

  fry_plot_data$direction_for_plot <- if (
    nzchar(direction_column)
  ) {
    as.character(
      fry_plot_data[[direction_column]]
    )
  } else {
    "NA"
  }

  plot_fry <- ggplot2::ggplot(
    fry_plot_data,
    ggplot2::aes(
      x = contrast,
      y = broad_cell_type,
      size = minus_log10_FDR,
      shape = direction_for_plot
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(
      ~ program_name,
      scales = "free_x",
      ncol = 2
    ) +
    ggplot2::geom_hline(
      yintercept = NA_real_
    ) +
    ggplot2::labs(
      title = "Directionally weighted locked VDR-program tests",
      subtitle = "limma fry; direction is defined relative to the locked calcitriol-dependent interaction signs",
      x = "Contrast",
      y = "Broad cell type",
      size = "-log10(FDR)",
      shape = "Direction"
    ) +
    theme_qc +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 5
      ),
      strip.text = ggplot2::element_text(
        size = 7
      )
    )

  save_plot_pair(
    plot_fry,
    "QC3_directional_VDR_fry",
    width = 16,
    height = 10
  )
}

volcano_data <- as.data.frame(
  LH7_primary_all
)

volcano_data$minus_log10_FDR <- -log10(
  pmax(
    volcano_data$FDR,
    .Machine$double.xmin
  )
)

volcano_data$plot_class <- ifelse(
  volcano_data$FDR < 0.05 &
    abs(volcano_data$logFC) >= 0.5,
  "Whole-genome DE",
  ifelse(
    volcano_data$RNA_ALL_436,
    "Locked VDR candidate",
    "Other gene"
  )
)

plot_LH7_volcano <- ggplot2::ggplot(
  volcano_data,
  ggplot2::aes(
    x = logFC,
    y = minus_log10_FDR,
    shape = plot_class
  )
) +
  ggplot2::geom_point(
    size = 0.45,
    alpha = 0.55
  ) +
  ggplot2::geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  ggplot2::facet_wrap(
    ~ broad_cell_type,
    scales = "free",
    ncol = 3
  ) +
  ggplot2::labs(
    title = "LH+7 RIF versus fertile gene-level effects",
    subtitle = "Dashed lines show FDR 0.05 and absolute log2 fold-change 0.5 reporting thresholds",
    x = "RIF minus fertile log2 fold-change",
    y = "-log10 genome-wide FDR",
    shape = "Gene class"
  ) +
  theme_qc

save_plot_pair(
  plot_LH7_volcano,
  "QC4_LH7_volcano",
  width = 12,
  height = 8
)

top_candidate_LH7 <- data.table::as.data.table(
  candidate_LH7_primary
)[
  order(PValue, -abs(logFC)),
  head(.SD, 8L),
  by = broad_cell_type
]

top_candidate_LH7$gene_label <- paste0(
  top_candidate_LH7$gene,
  " | ",
  top_candidate_LH7$broad_cell_type
)

top_candidate_LH7$gene_label <- factor(
  top_candidate_LH7$gene_label,
  levels = rev(unique(
    top_candidate_LH7$gene_label
  ))
)

plot_candidate_LH7 <- ggplot2::ggplot(
  as.data.frame(top_candidate_LH7),
  ggplot2::aes(
    x = logFC,
    y = gene_label,
    size = -log10(
      pmax(PValue, .Machine$double.xmin)
    ),
    shape = interaction_direction_concordant
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_point() +
  ggplot2::labs(
    title = "Top locked VDR candidates in the LH+7 comparison",
    subtitle = "Ranking uses nominal gene-level evidence; formal claims require genome-wide or explicitly labelled focused FDR",
    x = "RIF minus fertile log2 fold-change",
    y = "Candidate gene and broad cell type",
    size = "-log10(P)",
    shape = "Concordant with\ninteraction direction"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(
      size = 5
    )
  )

save_plot_pair(
  plot_candidate_LH7,
  "QC5_LH7_top_locked_candidates",
  width = 11,
  height = 13
)

LH7_count_plot <- as.data.frame(
  LH7_primary_summary
)

LH7_count_long <- data.table::melt(
  data.table::as.data.table(LH7_count_plot),
  id.vars = "broad_cell_type",
  measure.vars = c(
    "whole_genome_up_n",
    "whole_genome_down_n",
    "locked436_focused_FDR_lt_0_05_n"
  ),
  variable.name = "result_class",
  value.name = "gene_n"
)

plot_LH7_counts <- ggplot2::ggplot(
  as.data.frame(LH7_count_long),
  ggplot2::aes(
    x = broad_cell_type,
    y = gene_n,
    fill = result_class
  )
) +
  ggplot2::geom_col(
    position = "dodge"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "LH+7 RIF versus fertile gene-level summary",
    subtitle = "Locked-436 focused FDR is displayed separately from genome-wide DE",
    x = "Broad cell type",
    y = "Gene count",
    fill = "Result class"
  ) +
  theme_qc +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    )
  )

save_plot_pair(
  plot_LH7_counts,
  "QC6_LH7_DE_counts",
  width = 10,
  height = 6
)


# ----------------------------------------------------------------------
# 16. Gate, Card and interpretation
# ----------------------------------------------------------------------

expected_fry_n <- (
  length(temporal_types) *
    length(planned_temporal_contrasts) *
    4L
) + (
  length(LH7_types) *
    4L
) + (
  length(LH7_types) *
    2L *
    4L
)

gate_checks <- data.frame(
  check = c(
    "upstream_gate_authorizes_gene_DE",
    "pseudobulk_matrix_dimensions_locked",
    "pseudobulk_metadata_aligned",
    "locked_candidate_hierarchy_valid",
    "candidate_coverage_above_90_percent",
    "five_temporal_models_complete",
    "six_LH7_primary_models_complete",
    "twelve_LH7_sensitivity_models_present",
    "all_models_retain_at_least_1000_genes",
    "all_temporal_omnibus_FDR_values_finite",
    "all_LH7_primary_FDR_values_finite",
    "directional_fry_result_count_complete",
    "composition_inference_remains_disabled",
    "subtype_inference_remains_disabled"
  ),
  pass = c(
    upstream_gate$
      gene_level_pseudobulk_DE_authorized[[1]],
    nrow(pb_counts) == 36601L &&
      ncol(pb_counts) == 222L,
    identical(
      colnames(pb_counts),
      as.character(pb_metadata$pb_id)
    ),
    nrow(candidate_table) == 436L &&
      sum(candidate_table$occupancy_supported) ==
        284L &&
      sum(
        candidate_table$shared_occupancy_supported
      ) == 46L &&
      sum(
        candidate_table$any_promoter_overlap_support
      ) == 84L,
    all(
      candidate_coverage$coverage_fraction >= 0.90
    ),
    nrow(temporal_status) == 5L &&
      all(temporal_status$status == "COMPLETE"),
    nrow(LH7_status) == 6L &&
      all(LH7_status$status == "COMPLETE"),
    nrow(
      unique(
        LH7_sensitivity_all[
          ,
          c(
            "broad_cell_type",
            "contrast"
          )
        ]
      )
    ) == 12L,
    all(
      filtering_all$retained_gene_n >= 1000L
    ),
    all(
      is.finite(
        temporal_omnibus_all$FDR
      )
    ),
    all(
      is.finite(
        LH7_primary_all$FDR
      )
    ),
    nrow(fry_all) == expected_fry_n,
    !upstream_gate$
      composition_inference_authorized[[1]],
    !upstream_gate$
      subtype_inference_authorized[[1]]
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_checks,
  file.path(
    run_qc_root,
    "04_gate_checks.csv"
  )
)

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
  "PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE"
} else {
  "HOLD_GENE_LEVEL_PSEUDOBULK_DE_QC_FAILED"
}

DE_gate <- data.frame(
  dataset_id = dataset_id,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  upstream_model_version =
    upstream_model_version,
  DE_version = DE_version,
  gate_status = gate_status,
  temporal_model_n = nrow(temporal_status),
  LH7_primary_model_n = nrow(LH7_status),
  LH7_sensitivity_model_n = nrow(
    unique(
      LH7_sensitivity_all[
        ,
        c(
          "broad_cell_type",
          "contrast"
        )
      ]
    )
  ),
  temporal_omnibus_genome_FDR_lt_0_05_n =
    nrow(temporal_omnibus_significant),
  temporal_pairwise_genome_DE_n =
    nrow(temporal_pairwise_significant),
  LH7_primary_genome_DE_n =
    nrow(LH7_primary_significant),
  LH7_sensitivity_genome_DE_n =
    nrow(LH7_sensitivity_significant),
  locked_program_fry_test_n =
    nrow(fry_all),
  integrated_interpretation_authorized =
    all_gate_checks_pass,
  GSE183837_independent_localization_authorized =
    all_gate_checks_pass,
  composition_inference_authorized =
    FALSE,
  subtype_inference_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "GENE_DE_GATE.csv"
)

write_csv_utf8(
  DE_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1H Gene-level Pseudobulk DE Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0("- Pseudobulk version: ", pseudobulk_version),
  paste0("- DE version: ", DE_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Temporal cell-type models: ",
    nrow(temporal_status)
  ),
  paste0(
    "- LH+7 primary cell-type models: ",
    nrow(LH7_status)
  ),
  paste0(
    "- Temporal omnibus genes with genome-wide FDR < 0.05: ",
    nrow(temporal_omnibus_significant)
  ),
  paste0(
    "- Temporal planned-pairwise genes meeting genome-wide FDR and effect threshold: ",
    nrow(temporal_pairwise_significant)
  ),
  paste0(
    "- LH+7 primary genes meeting genome-wide FDR and effect threshold: ",
    nrow(LH7_primary_significant)
  ),
  paste0(
    "- Directional locked-program fry tests: ",
    nrow(fry_all)
  ),
  paste0(
    "- Integrated interpretation authorized: ",
    DE_gate$integrated_interpretation_authorized
  ),
  paste0(
    "- GSE183837 independent localization authorized: ",
    DE_gate$GSE183837_independent_localization_authorized
  ),
  "- Composition inference authorized: FALSE",
  "- Subtype inference authorized: FALSE",
  "- Result lock authorized: FALSE"
)

card_path <- file.path(
  run_qc_root,
  "GENE_DE_RESULT_CARD.md"
)

write_text(
  card_lines,
  card_path
)

interpretation_lines <- c(
  "# STEP04E1H gene-level pseudobulk interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Primary inference",
  "",
  "- Genome-wide BH FDR within each cell type and contrast is the primary gene-level multiplicity control.",
  "- Pairwise DE reporting additionally uses an absolute log2 fold-change threshold of 0.5.",
  "- The temporal omnibus test asks whether any of the five fertile LH time groups differ.",
  "- The LH+7 primary contrast compares RIF with all eligible fertile LH+7 donors.",
  "",
  "## VDR-focused inference",
  "",
  "- The 436/284/46/84 VDR hierarchies were locked before the present gene-level analysis.",
  "- Candidate-restricted FDR values are therefore reported as focused secondary evidence.",
  "- Candidate-restricted FDR does not replace genome-wide FDR and must be labelled explicitly.",
  "- Directionally weighted fry tests evaluate whether a locked gene set changes in the same direction as the calcitriol-dependent VDR interaction program.",
  "",
  "## Preserved limitations",
  "",
  "- The fertile temporal series is cross-sectional.",
  "- Cell-composition inference remains outside the current Gate.",
  "- Candidate subtype labels are not used for differential-expression testing.",
  "- Results remain candidates until content, statistical and visual QC are reviewed."
)

interpretation_path <- file.path(
  run_qc_root,
  "GENE_DE_INTERPRETATION.md"
)

write_text(
  interpretation_lines,
  interpretation_path
)


# ----------------------------------------------------------------------
# 17. Methods and Results shell
# ----------------------------------------------------------------------

methods_lines <- c(
  "# STEP04E1H Methods draft",
  "",
  paste(
    "Gene-level differential-expression analyses were performed separately",
    "for each eligible broad cell type using donor-level raw-count",
    "pseudobulk profiles."
  ),
  paste(
    "Genes were filtered using edgeR filterByExpr with the corresponding",
    "cell-type-specific design matrix, followed by trimmed mean of M-values",
    "normalization."
  ),
  paste(
    "Negative-binomial dispersions were estimated with robust empirical",
    "Bayes procedures, and robust quasi-likelihood generalized linear",
    "models were fitted."
  ),
  paste(
    "Fertile temporal effects were tested using a multi-degree-of-freedom",
    "contrast comparing LH+5, LH+7, LH+9 and LH+11 with LH+3."
  ),
  paste(
    "Seven planned temporal pairwise contrasts were evaluated separately."
  ),
  paste(
    "At LH+7, the primary comparison contrasted RIF with all fertile donors;",
    "sensitivity models contrasted RIF separately with standard fertile and",
    "PGT-associated fertile donors."
  ),
  paste(
    "Benjamini-Hochberg correction was applied genome-wide within each cell",
    "type and contrast. Pairwise genes were summarized as differentially",
    "expressed when FDR was below 0.05 and the absolute log2 fold-change was",
    "at least 0.5."
  ),
  paste(
    "Because the VDR candidate hierarchy was frozen before this analysis,",
    "separate focused BH values were also calculated within the locked",
    "436-, 284-, 46- and 84-gene sets; these values were interpreted as",
    "secondary candidate-focused evidence."
  ),
  paste(
    "Directionally weighted limma fry tests were used for one-degree-of-",
    "freedom contrasts, with gene weights defined by the sign of the locked",
    "calcitriol-dependent VDR interaction log-fold change."
  )
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Methods_draft.md"
  )
)

results_shell <- c(
  "# STEP04E1H Results shell",
  "",
  paste0(
    "Gene-level pseudobulk models were completed for ",
    nrow(temporal_status),
    " fertile temporal cell types and ",
    nrow(LH7_status),
    " LH+7 RIF-versus-fertile cell types."
  ),
  paste0(
    "Across the five temporal omnibus models, ",
    nrow(temporal_omnibus_significant),
    " cell-type-specific gene results met genome-wide FDR < 0.05."
  ),
  paste0(
    "Across the seven planned temporal contrasts, ",
    nrow(temporal_pairwise_significant),
    " gene results met genome-wide FDR < 0.05 and absolute log2 fold-change >= 0.5."
  ),
  paste0(
    "In the primary LH+7 RIF-versus-fertile analyses, ",
    nrow(LH7_primary_significant),
    " gene results met the corresponding genome-wide FDR and effect-size criteria."
  ),
  "",
  "Biological interpretation will be finalized after reviewing cell-type specificity, locked VDR-candidate concordance, directional gene-set tests and sensitivity analyses."
)

write_text(
  results_shell,
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Results_shell.md"
  )
)


# ----------------------------------------------------------------------
# 18. Manifests, sessionInfo and handoff ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_names <- c(
  "Matrix",
  "edgeR",
  "limma",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scales"
)

package_table <- data.frame(
  package = package_names,
  version = vapply(
    package_names,
    safe_package_version,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_table,
  file.path(
    run_qc_root,
    "PACKAGE_VERSIONS.csv"
  )
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

artifact_files <- unique(c(
  list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Results_shell.md"
  )
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(artifact_files),
  size_bytes = file.info(artifact_files)$size,
  md5 = unname(
    tools::md5sum(artifact_files)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(
    run_qc_root,
    "ARTIFACT_MANIFEST.csv"
  )
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  input_version = input_version,
  pseudobulk_version = pseudobulk_version,
  upstream_model_version =
    upstream_model_version,
  DE_version = DE_version,
  supersedes_partial_DE_versions =
    "GENEDE_v1_0;GENEDE_v1_1",
  resume_source_root =
    normalize_slashes(previous_partial_root),
  dataset_id = dataset_id,
  started_at = as.character(run_start),
  completed_at = as.character(run_end),
  elapsed_minutes = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "mins"
      )
    ),
    2
  ),
  gate_status = gate_status,
  temporal_model_n = nrow(temporal_status),
  LH7_primary_model_n = nrow(LH7_status),
  temporal_omnibus_significant_n =
    nrow(temporal_omnibus_significant),
  temporal_pairwise_significant_n =
    nrow(temporal_pairwise_significant),
  LH7_primary_significant_n =
    nrow(LH7_primary_significant),
  locked_program_fry_test_n =
    nrow(fry_all),
  large_Seurat_checkpoint_loaded =
    FALSE,
  composition_inference_authorized =
    FALSE,
  subtype_inference_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  output_root =
    normalize_slashes(output_root),
  qc_root =
    normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "temporal_omnibus_significant_n=",
  nrow(temporal_omnibus_significant),
  "temporal_pairwise_significant_n=",
  nrow(temporal_pairwise_significant),
  "LH7_primary_significant_n=",
  nrow(LH7_primary_significant)
)

small_handoff_files <- unique(c(
  list.files(
    table_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    candidate_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    source_data_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    run_qc_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  log_path,
  run_manifest_path,
  parameter_path,
  session_path,
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Methods_draft.md"
  ),
  file.path(
    manuscript_text_root,
    "STEP04E1H_gene_DE_Results_shell.md"
  ),
  normalizePath(
    file.path(
      project_root,
      "04_code",
      "R",
      script_name
    ),
    winslash = "/",
    mustWork = FALSE
  )
))

small_handoff_files <- small_handoff_files[
  file.exists(small_handoff_files) &
    !dir.exists(small_handoff_files) &
    file.info(small_handoff_files)$size <=
      150 * 1024^2
]

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos = "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

zip_ok <- FALSE
zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
  zip_ok <- tryCatch(
    {
      zip::zipr(
        zipfile = handoff_zip,
        files = small_handoff_files,
        recurse = FALSE,
        include_directories = FALSE
      )

      file.exists(handoff_zip)
    },
    error = function(e) {
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1H GSE250130 gene-level pseudobulk DE completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Temporal models: ",
  nrow(temporal_status),
  "\n",
  sep = ""
)
cat(
  "Temporal omnibus FDR < 0.05 results: ",
  nrow(temporal_omnibus_significant),
  "\n",
  sep = ""
)
cat(
  "Temporal pairwise DE results: ",
  nrow(temporal_pairwise_significant),
  "\n",
  sep = ""
)
cat(
  "LH7 primary DE results: ",
  nrow(LH7_primary_significant),
  "\n",
  sep = ""
)
cat(
  "Directional VDR-program fry tests: ",
  nrow(fry_all),
  "\n",
  sep = ""
)
cat("Output root: ", output_root, "\n", sep = "")
cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
  cat(
    "Upload this ZIP:\n",
    handoff_zip,
    "\n"
  )
} else {
  cat(
    "ZIP creation failed: ",
    zip_error,
    "\nThe results and QC folders are complete.\n",
    sep = ""
  )
}
