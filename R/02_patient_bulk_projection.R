#!/usr/bin/env Rscript
# Public analysis module for the VDR endometrial transcriptome study.
# Set VDR_PROJECT_ROOT to the project directory containing the expected public GEO inputs.
# This module preserves the executed analysis order; raw GEO data are not redistributed.


# ==============================================================================
# Analysis stage 1: 48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate
# ==============================================================================
# ============================================================
# STEP04A v2 - Acquire missing GSE250130 source, correct the
# temporal dataset role, narrow GSE165004 slide-QC provenance,
# and build the five-dataset source-lock candidate
# Nutrition_Endometrial_VDR_v2
#
# Scientific correction:
# - GSE250130 is a 28-donor time-series single-cell endometrial
#   atlas across LH+3, +5, +7, +9 and +11, with RIF samples at LH+7.
# - It is therefore locked as
#   TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT,
#   not as a generic bulk temporal matrix.
#
# Technical correction:
# - STEP04A v1 counted every file under STEP00 as a possible
#   GSE165004 slide-QC artifact, producing 407 rows.
# - v2 retains only four accepted canonical statistical outputs.
#
# This script:
# - downloads the GSE250130 series matrix and 3.4-GB RAW outer TAR
#   with resume/retry support;
# - inventories the outer TAR without extracting all nested bundles;
# - creates a corrected five-dataset role lock;
# - creates a source-lock candidate for all five datasets;
# - writes a narrow four-file GSE165004 slide-QC source lock;
# - does not run temporal, patient or single-cell analyses.
#
# Run from project root:
# source(
#   "04_code/R/48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(3600L, getOption("timeout")))

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04A_v2_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

download_large_GSE250130_RAW <- TRUE

GSE250130_series_matrix_url <-
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250130/matrix/GSE250130_series_matrix.txt.gz"

GSE250130_RAW_url <-
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250130/suppl/GSE250130_RAW.tar"

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
  control = file.path(step04_root, "00_module_control"),
  inputs = file.path(step04_root, "01_inputs"),
  code = file.path(step04_root, "02_code"),
  parameters = file.path(step04_root, "03_parameters"),
  environment = file.path(step04_root, "04_environment"),
  execution = file.path(step04_root, "05_execution"),
  results = file.path(step04_root, "06_results"),
  figures = file.path(step04_root, "07_figures_tables"),
  qc = file.path(step04_root, "08_qc"),
  text = file.path(step04_root, "09_manuscript_text"),
  archive_parent = file.path(step04_root, "10_archive_work")
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

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4Av2_",
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

download_with_resume <- function(
    url,
    destination,
    large_file = FALSE
) {
  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  curl_candidates <- c(
    Sys.which("curl.exe"),
    Sys.which("curl")
  )

  curl_candidates <- curl_candidates[
    nzchar(curl_candidates)
  ]

  if (length(curl_candidates) > 0L) {
    curl_binary <- curl_candidates[[1]]

    args <- c(
      "--fail",
      "--location",
      "--retry",
      "20",
      "--retry-all-errors",
      "--retry-delay",
      "10",
      "--connect-timeout",
      "60",
      "--continue-at",
      "-",
      "--output",
      shQuote(destination),
      shQuote(url)
    )

    log_event(
      "Downloading with resume support:",
      basename(destination)
    )

    status <- system2(
      command = curl_binary,
      args = args,
      stdout = "",
      stderr = ""
    )

    if (
      status != 0L ||
        !file.exists(destination)
    ) {
      stop(
        "curl download failed for: ",
        url
      )
    }

    return(TRUE)
  }

  if (
    large_file &&
      file.exists(destination)
  ) {
    stop(
      "A partial large file exists but curl is unavailable for resume: ",
      destination
    )
  }

  log_event(
    "Downloading with R download.file:",
    basename(destination)
  )

  download.file(
    url = url,
    destfile = destination,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )

  if (!file.exists(destination)) {
    stop(
      "download.file failed for: ",
      url
    )
  }

  TRUE
}

parse_GEO_series_metadata <- function(path) {
  connection <- gzfile(
    path,
    open = "rt"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  all_lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "UTF-8"
  )

  sample_lines <- all_lines[
    grepl(
      "^!Sample_",
      all_lines
    )
  ]

  if (length(sample_lines) == 0L) {
    return(
      data.frame()
    )
  }

  parsed_rows <- lapply(
    sample_lines,
    function(line_value) {
      fields <- strsplit(
        line_value,
        "\t",
        fixed = TRUE
      )[[1]]

      field_name <- sub(
        "^!Sample_",
        "",
        fields[[1]]
      )

      values <- fields[-1]

      values <- gsub(
        '^"|"$',
        "",
        values
      )

      data.frame(
        field_name =
          field_name,
        sample_index =
          seq_along(values),
        value =
          values,
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    parsed_rows
  )
}

build_source_row <- function(
    dataset_id,
    source_role,
    path,
    source_status =
      "CANDIDATE_PENDING_MANUAL_LOCK"
) {
  data.frame(
    dataset_id =
      dataset_id,
    source_role =
      source_role,
    path =
      normalize_slashes(
        path
      ),
    file_name =
      basename(path),
    file_exists =
      file.exists(path),
    size_bytes = if (
      file.exists(path)
    ) {
      as.numeric(
        file.info(path)$size
      )
    } else {
      NA_real_
    },
    md5 = if (
      file.exists(path)
    ) {
      md5_file(path)
    } else {
      NA_character_
    },
    sha256 = if (
      file.exists(path)
    ) {
      sha256_file(path)
    } else {
      NA_character_
    },
    source_status =
      source_status,
    stringsAsFactors = FALSE
  )
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Validate STEP04A v1 hold
# ----------------------------------------------------------------------

v1_gate_path <- file.path(
  dirs$results,
  "STEP04A_GATE_TABLE.csv"
)

v1_readiness_path <- file.path(
  dirs$results,
  "STEP04A_INPUT_READINESS_MATRIX.csv"
)

if (
  !file.exists(v1_gate_path) ||
    !file.exists(v1_readiness_path)
) {
  stop(
    "STEP04A v1 gate or readiness matrix is missing."
  )
}

v1_gate <- read.csv(
  v1_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v1_readiness <- read.csv(
  v1_readiness_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v1_gate_value <- v1_gate$status[
  v1_gate$gate_component ==
    "STEP04A gate"
]

if (
  length(v1_gate_value) != 1L ||
    v1_gate_value !=
      "HOLD_MISSING_REQUIRED_DOWNSTREAM_DATASET_FILES"
) {
  stop(
    "Unexpected STEP04A v1 gate: ",
    paste(
      v1_gate_value,
      collapse = ", "
    )
  )
}

if (
  v1_readiness$preliminary_readiness[
    v1_readiness$dataset_id ==
      "GSE250130"
  ] !=
    "MISSING"
) {
  stop(
    "STEP04A v2 correction is not applicable because GSE250130 was not marked missing."
  )
}

# ----------------------------------------------------------------------
# Corrected dataset-role lock
# ----------------------------------------------------------------------

dataset_role_lock_v2 <- data.frame(
  dataset_id = c(
    "GSE250130",
    "GSE165004",
    "GSE92324",
    "GSE111974",
    "GSE183837"
  ),
  analysis_role = c(
    "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT",
    "PRIMARY_PATIENT_PROJECTION",
    "EXTERNAL_PATIENT_SUPPORT",
    "EXTERNAL_PATIENT_SUPPORT_WITH_SHARED_CONTROL_CONSTRAINT",
    "INDEPENDENT_SINGLE_CELL_LOCALIZATION"
  ),
  priority = c(
    1L,
    1L,
    2L,
    2L,
    1L
  ),
  required_for_main_manuscript =
    TRUE,
  primary_input_expected = c(
    "28 sample-level single-cell count bundles plus LH-timepoint, fertility/RIF and donor metadata",
    "processed endometrial expression matrix plus phenotype and slide metadata",
    "processed endometrial expression matrix plus phenotype metadata",
    "processed endometrial expression matrix plus phenotype metadata and explicit shared-control audit",
    "9 sample-level single-cell count triplets plus Ctrl/RIF and donor metadata"
  ),
  locked_constraints = c(
    "Use fertile LH+3/+5/+7/+9/+11 samples for temporal support; RIF is available at LH+7. Use donor-aware cell-type or pseudobulk analyses and do not pool cells as independent replicates.",
    "Use the accepted slide-aware and leave-one-slide robustness framework.",
    "Use as external supportive validation; do not train a diagnostic classifier.",
    "Account explicitly for controls shared with GSE165004.",
    "Use for independent cell-type/state localization; do not infer causality from localization."
  ),
  version =
    "STEP04A_DATASET_ROLE_LOCK_v2.0",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  dataset_role_lock_v2,
  file.path(
    dirs$results,
    "STEP04A_v2_DATASET_ROLE_LOCK.csv"
  )
)

# ----------------------------------------------------------------------
# Acquire GSE250130
# ----------------------------------------------------------------------

GSE250130_root <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE250130"
)

GSE250130_matrix_dir <- file.path(
  GSE250130_root,
  "series_matrix"
)

GSE250130_raw_dir <- file.path(
  GSE250130_root,
  "raw_counts"
)

dir.create(
  GSE250130_matrix_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  GSE250130_raw_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

GSE250130_series_matrix_path <- file.path(
  GSE250130_matrix_dir,
  "GSE250130_series_matrix.txt.gz"
)

GSE250130_RAW_path <- file.path(
  GSE250130_raw_dir,
  "GSE250130_RAW.tar"
)

if (!file.exists(
  GSE250130_series_matrix_path
)) {
  download_with_resume(
    GSE250130_series_matrix_url,
    GSE250130_series_matrix_path,
    large_file = FALSE
  )
}

if (
  download_large_GSE250130_RAW &&
    !file.exists(
      GSE250130_RAW_path
    )
) {
  download_with_resume(
    GSE250130_RAW_url,
    GSE250130_RAW_path,
    large_file = TRUE
  )
}

if (
  download_large_GSE250130_RAW &&
    file.exists(
      GSE250130_RAW_path
    ) &&
    as.numeric(
      file.info(
        GSE250130_RAW_path
      )$size
    ) < 2500000000
) {
  stop(
    "GSE250130_RAW.tar is smaller than 2.5 GB and is likely incomplete. Re-run the script to resume."
  )
}

# ----------------------------------------------------------------------
# Parse GSE250130 series metadata
# ----------------------------------------------------------------------

GSE250130_series_metadata_long <-
  parse_GEO_series_metadata(
    GSE250130_series_matrix_path
  )

if (
  nrow(
    GSE250130_series_metadata_long
  ) == 0L
) {
  stop(
    "GSE250130 series metadata could not be parsed."
  )
}

write_csv_utf8(
  GSE250130_series_metadata_long,
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_SERIES_METADATA_LONG.csv"
  )
)

GSE250130_sample_accessions <-
  GSE250130_series_metadata_long$value[
    GSE250130_series_metadata_long$field_name ==
      "geo_accession"
  ]

GSE250130_sample_accessions <- unique(
  GSE250130_sample_accessions[
    nzchar(
      GSE250130_sample_accessions
    )
  ]
)

if (
  length(
    GSE250130_sample_accessions
  ) != 28L
) {
  stop(
    "Expected 28 GSE250130 samples but parsed ",
    length(
      GSE250130_sample_accessions
    ),
    "."
  )
}

# ----------------------------------------------------------------------
# Inventory the outer TAR without extracting nested bundles
# ----------------------------------------------------------------------

tar_member_list <- character(0)
tar_inventory_status <- "RAW_TAR_NOT_PRESENT"

if (file.exists(
  GSE250130_RAW_path
)) {
  tar_member_list <- tryCatch(
    utils::untar(
      GSE250130_RAW_path,
      list = TRUE
    ),
    error = function(error_object) {
      character(0)
    }
  )

  if (length(tar_member_list) > 0L) {
    tar_inventory_status <-
      "RAW_TAR_LIST_OK"
  } else {
    tar_inventory_status <-
      "RAW_TAR_PRESENT_BUT_LIST_FAILED"
  }
}

write_text(
  tar_member_list,
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt"
  )
)

tar_summary <- data.frame(
  metric = c(
    "RAW tar exists",
    "RAW tar size bytes",
    "RAW tar member-list status",
    "RAW tar member count",
    "Parsed GEO sample count",
    "Series matrix exists"
  ),
  value = c(
    as.character(
      file.exists(
        GSE250130_RAW_path
      )
    ),
    if (
      file.exists(
        GSE250130_RAW_path
      )
    ) {
      as.character(
        as.numeric(
          file.info(
            GSE250130_RAW_path
          )$size
        )
      )
    } else {
      NA_character_
    },
    tar_inventory_status,
    as.character(
      length(
        tar_member_list
      )
    ),
    as.character(
      length(
        GSE250130_sample_accessions
      )
    ),
    as.character(
      file.exists(
        GSE250130_series_matrix_path
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  tar_summary,
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_ACQUISITION_QC.csv"
  )
)

# ----------------------------------------------------------------------
# Narrow accepted GSE165004 slide-QC provenance
# ----------------------------------------------------------------------

step00_locked_stat_root <- file.path(
  project_root,
  "06_locked_results",
  "statistical_outputs",
  "STEP00"
)

slide_QC_paths <- file.path(
  step00_locked_stat_root,
  c(
    "GSE165004_leave_one_slide_out_metrics.csv",
    "GSE165004_slide_aggregate_top1000_per_contrast.csv",
    "GSE165004_slide_by_group_counts.csv",
    "GSE165004_slide_confounding_summary.csv"
  )
)

missing_slide_QC <- slide_QC_paths[
  !file.exists(slide_QC_paths)
]

if (length(missing_slide_QC) > 0L) {
  stop(
    "Accepted GSE165004 slide-QC file is missing:\n",
    paste(
      missing_slide_QC,
      collapse = "\n"
    )
  )
}

slide_QC_lock <- do.call(
  rbind,
  lapply(
    slide_QC_paths,
    function(path_value) {
      build_source_row(
        dataset_id = "GSE165004",
        source_role =
          "LOCKED_SLIDE_ROBUSTNESS_STATISTICAL_OUTPUT",
        path = path_value,
        source_status =
          "LOCKED_ACCEPTED_STEP00_SOURCE"
      )
    }
  )
)

write_csv_utf8(
  slide_QC_lock,
  file.path(
    dirs$results,
    "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv"
  )
)

# ----------------------------------------------------------------------
# Build five-dataset source-lock candidate
# ----------------------------------------------------------------------

source_rows <- list()

add_source <- function(
    dataset_id,
    source_role,
    path,
    status =
      "CANDIDATE_PENDING_MANUAL_LOCK"
) {
  source_rows[[
    length(source_rows) + 1L
  ]] <<- build_source_row(
    dataset_id =
      dataset_id,
    source_role =
      source_role,
    path =
      path,
    source_status =
      status
  )
}

# GSE250130
add_source(
  "GSE250130",
  "SERIES_MATRIX_METADATA",
  GSE250130_series_matrix_path
)

add_source(
  "GSE250130",
  "RAW_OUTER_TAR_OF_SAMPLE_BUNDLES",
  GSE250130_RAW_path
)

# GSE165004
add_source(
  "GSE165004",
  "SERIES_MATRIX_EXPRESSION",
  file.path(
    project_root,
    "03_data",
    "raw_external_READ_ONLY",
    "GSE165004",
    "series_matrix",
    "GSE165004_series_matrix.txt.gz"
  )
)

add_source(
  "GSE165004",
  "LOCKED_SERIES_METADATA",
  file.path(
    step00_locked_stat_root,
    "GSE165004_series_metadata_parsed.csv"
  )
)

add_source(
  "GSE165004",
  "LOCKED_GROUP_METADATA",
  file.path(
    step00_locked_stat_root,
    "GSE165004_title_group_audit.csv"
  )
)

for (path_value in slide_QC_paths) {
  add_source(
    "GSE165004",
    "LOCKED_SLIDE_QC",
    path_value,
    status =
      "LOCKED_ACCEPTED_STEP00_SOURCE"
  )
}

# GSE92324
add_source(
  "GSE92324",
  "SERIES_MATRIX_EXPRESSION_AND_METADATA",
  file.path(
    project_root,
    "03_data",
    "raw_external_READ_ONLY",
    "GSE92324",
    "series_matrix",
    "GSE92324_series_matrix.txt.gz"
  )
)

add_source(
  "GSE92324",
  "PLATFORM_ANNOTATION",
  file.path(
    project_root,
    "03_data",
    "raw_external_READ_ONLY",
    "GSE92324",
    "platform_annotation",
    "GPL10558.soft.gz"
  )
)

# GSE111974
add_source(
  "GSE111974",
  "SERIES_MATRIX_EXPRESSION",
  file.path(
    project_root,
    "03_data",
    "raw_external_READ_ONLY",
    "GSE111974",
    "series_matrix",
    "GSE111974_series_matrix.txt.gz"
  )
)

add_source(
  "GSE111974",
  "LOCKED_SERIES_METADATA",
  file.path(
    step00_locked_stat_root,
    "GSE111974_series_metadata_parsed.csv"
  )
)

add_source(
  "GSE111974",
  "LOCKED_GROUP_METADATA",
  file.path(
    step00_locked_stat_root,
    "GSE111974_title_group_audit.csv"
  )
)

# GSE183837
GSE183837_root <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE183837"
)

GSE183837_triplet_files <- list.files(
  file.path(
    GSE183837_root,
    "raw_counts",
    "matrix_market_by_donor"
  ),
  pattern =
    "\\.(barcodes|genes|matrix)\\.(tsv|mtx)\\.gz$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (
  length(
    GSE183837_triplet_files
  ) != 27L
) {
  stop(
    "Expected 27 GSE183837 triplet files but found ",
    length(
      GSE183837_triplet_files
    ),
    "."
  )
}

for (path_value in sort(
  GSE183837_triplet_files
)) {
  add_source(
    "GSE183837",
    "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER",
    path_value
  )
}

add_source(
  "GSE183837",
  "SERIES_MATRIX_METADATA",
  file.path(
    GSE183837_root,
    "series_matrix",
    "GSE183837_series_matrix.txt.gz"
  )
)

source_lock_candidate <- do.call(
  rbind,
  source_rows
)

missing_source_rows <- source_lock_candidate[
  !source_lock_candidate$file_exists,
  ,
  drop = FALSE
]

write_csv_utf8(
  source_lock_candidate,
  file.path(
    dirs$results,
    "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv"
  )
)

write_csv_utf8(
  missing_source_rows,
  file.path(
    dirs$qc,
    "STEP04A_v2_MISSING_SOURCE_ROWS.csv"
  )
)

source_summary <- aggregate(
  file_exists ~
    dataset_id,
  data = source_lock_candidate,
  FUN = function(x) {
    paste0(
      sum(x),
      "/",
      length(x)
    )
  }
)

names(source_summary)[
  names(source_summary) ==
    "file_exists"
] <- "files_present"

source_summary$all_required_files_present <-
  vapply(
    source_summary$dataset_id,
    function(dataset_id) {
      all(
        source_lock_candidate$file_exists[
          source_lock_candidate$dataset_id ==
            dataset_id
        ]
      )
    },
    logical(1)
  )

write_csv_utf8(
  source_summary,
  file.path(
    dirs$results,
    "STEP04A_v2_DATASET_SOURCE_READINESS.csv"
  )
)

# ----------------------------------------------------------------------
# Gate
# ----------------------------------------------------------------------

all_required_present <- all(
  source_lock_candidate$file_exists
)

GSE250130_RAW_ready <-
  file.exists(
    GSE250130_RAW_path
  ) &&
  as.numeric(
    file.info(
      GSE250130_RAW_path
    )$size
  ) >=
    2500000000 &&
  tar_inventory_status ==
    "RAW_TAR_LIST_OK"

gate_status <- if (
  !GSE250130_RAW_ready
) {
  "HOLD_GSE250130_RAW_DOWNLOAD_OR_TAR_INVENTORY_INCOMPLETE"
} else if (
  !all_required_present
) {
  "HOLD_OTHER_REQUIRED_SOURCE_FILES_MISSING"
} else {
  "PROCEED_TO_MANUAL_STEP04_INPUT_SOURCE_LOCK"
}

gate_table <- data.frame(
  gate_component = c(
    "STEP04A v1 hold",
    "GSE250130 scientific role",
    "GSE250130 series metadata",
    "GSE250130 sample count",
    "GSE250130 RAW tar",
    "GSE250130 tar inventory",
    "GSE165004 source candidate",
    "GSE165004 slide-QC source lock",
    "GSE92324 source candidate",
    "GSE111974 source candidate",
    "GSE183837 triplet source candidate",
    "All five dataset source rows",
    "Scientific figures",
    "STEP04A v2 gate"
  ),
  status = c(
    "CONFIRMED_ONLY_GSE250130_MISSING",
    "LOCKED_TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT",
    ifelse(
      file.exists(
        GSE250130_series_matrix_path
      ),
      "PASS",
      "MISSING"
    ),
    paste0(
      "PARSED_",
      length(
        GSE250130_sample_accessions
      )
    ),
    ifelse(
      file.exists(
        GSE250130_RAW_path
      ),
      paste0(
        "PRESENT_BYTES_",
        as.numeric(
          file.info(
            GSE250130_RAW_path
          )$size
        )
      ),
      "MISSING"
    ),
    tar_inventory_status,
    ifelse(
      all(
        source_lock_candidate$file_exists[
          source_lock_candidate$dataset_id ==
            "GSE165004"
        ]
      ),
      "PASS",
      "MISSING_SOURCE"
    ),
    "LOCKED_4_CANONICAL_FILES",
    ifelse(
      all(
        source_lock_candidate$file_exists[
          source_lock_candidate$dataset_id ==
            "GSE92324"
        ]
      ),
      "PASS",
      "MISSING_SOURCE"
    ),
    ifelse(
      all(
        source_lock_candidate$file_exists[
          source_lock_candidate$dataset_id ==
            "GSE111974"
        ]
      ),
      "PASS",
      "MISSING_SOURCE"
    ),
    ifelse(
      all(
        source_lock_candidate$file_exists[
          source_lock_candidate$dataset_id ==
            "GSE183837"
        ]
      ),
      "PASS_27_TRIPLET_FILES_PLUS_METADATA",
      "MISSING_SOURCE"
    ),
    ifelse(
      all_required_present,
      "PASS",
      paste0(
        "MISSING_ROWS_",
        nrow(
          missing_source_rows
        )
      )
    ),
    "NONE_CREATED",
    gate_status
  ),
  required_action = c(
    "Retain v1 as the discovery audit.",
    "Use donor-aware temporal single-cell analysis.",
    "Use for sample metadata only.",
    "Require exactly 28 GEO samples.",
    "Resume download until the full outer TAR is present.",
    "Require a readable TAR member list before source lock.",
    "Use series matrix plus locked STEP00 metadata and slide evidence.",
    "Ignore the overinclusive 407-row v1 inventory.",
    "Use series matrix and GPL10558 annotation.",
    "Use series matrix and locked metadata with shared-control audit.",
    "Use nine complete donor triplets.",
    "Lock only after manual review of this candidate.",
    "No figure workflow is triggered.",
    "Upload the source candidate, readiness, GSE250130 QC, narrow slide lock and gate."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP04A_v2_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP04A v2 acquisition and source-lock candidate card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    "- GSE250130 role: TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT",
    "- GSE250130 design: 28 biopsies across LH+3/+5/+7/+9/+11; RIF at LH+7",
    paste0(
      "- GSE250130 GEO samples parsed: ",
      length(
        GSE250130_sample_accessions
      )
    ),
    paste0(
      "- GSE250130 RAW tar present: ",
      file.exists(
        GSE250130_RAW_path
      )
    ),
    paste0(
      "- GSE250130 RAW tar bytes: ",
      if (
        file.exists(
          GSE250130_RAW_path
        )
      ) {
        as.numeric(
          file.info(
            GSE250130_RAW_path
          )$size
        )
      } else {
        NA
      }
    ),
    paste0(
      "- GSE250130 tar member count: ",
      length(
        tar_member_list
      )
    ),
    "- GSE165004 slide-QC provenance: 4 canonical locked files",
    paste0(
      "- Combined source rows: ",
      nrow(
        source_lock_candidate
      )
    ),
    paste0(
      "- Missing source rows: ",
      nrow(
        missing_source_rows
      )
    ),
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_ACQUISITION_SOURCE_CANDIDATE_CARD.md"
  )
)

write_text(
  c(
    "# STEP04A v2 correction note",
    "",
    "STEP04A v1 correctly identified GSE250130 as the only missing mandatory dataset. However, its role was described too generically. GSE250130 is a time-series single-cell endometrial atlas and must be analyzed with donor-aware cell-type or pseudobulk methods.",
    "",
    "The v1 GSE165004 slide-QC inventory was also overinclusive because every file under the STEP00 path was collected. The accepted provenance is narrowed to four canonical statistical outputs in the locked STEP00 directory.",
    "",
    "No downstream biological analysis or scientific figure is generated in this acquisition step."
  ),
  file.path(
    dirs$text,
    "STEP04A_v2_Correction_note.md"
  )
)

parameters <- list(
  schema_version = "2.0",
  run_id = run_id,
  GSE250130_role =
    "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT",
  GSE250130_expected_sample_n =
    28,
  GSE250130_parsed_sample_n =
    length(
      GSE250130_sample_accessions
    ),
  GSE250130_RAW_download_enabled =
    download_large_GSE250130_RAW,
  GSE250130_RAW_present =
    file.exists(
      GSE250130_RAW_path
    ),
  GSE250130_RAW_size_bytes =
    if (
      file.exists(
        GSE250130_RAW_path
      )
    ) {
      as.numeric(
        file.info(
          GSE250130_RAW_path
        )$size
      )
    } else {
      NA
    },
  GSE250130_tar_member_n =
    length(
      tar_member_list
    ),
  GSE165004_slide_QC_locked_file_n =
    nrow(
      slide_QC_lock
    ),
  combined_source_row_n =
    nrow(
      source_lock_candidate
    ),
  missing_source_row_n =
    nrow(
      missing_source_rows
    ),
  biological_analysis_started =
    FALSE,
  scientific_figures_created =
    FALSE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP04A_v2_parameters.json"
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
    "sessionInfo_POST_RUN_STEP04A_v2.txt"
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
      "48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  safe_file_copy(
    source_file,
    file.path(
      dirs$code,
      "48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R"
    ),
    overwrite = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
}

# ----------------------------------------------------------------------
# Short archive
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP04A_v2_DATASET_ROLE_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv"
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_DATASET_SOURCE_READINESS.csv"
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP04A_v2_ACQUISITION_SOURCE_CANDIDATE_CARD.md"
  ),
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_SERIES_METADATA_LONG.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt"
  ),
  file.path(
    dirs$qc,
    "STEP04A_v2_GSE250130_ACQUISITION_QC.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04A_v2_MISSING_SOURCE_ROWS.csv"
  ),
  file.path(
    dirs$text,
    "STEP04A_v2_Correction_note.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04A_v2_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04A_v2.txt"
  ),
  file.path(
    dirs$code,
    "48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R"
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
      "Correction",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "QC|READINESS|LOCK|CANDIDATE|CARD|MEMBER|MISSING",
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
    "STEP04A_v2_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04A v2 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    "- GSE250130 role corrected: YES",
    paste0(
      "- GSE250130 RAW present: ",
      file.exists(
        GSE250130_RAW_path
      )
    ),
    paste0(
      "- GSE250130 GEO samples: ",
      length(
        GSE250130_sample_accessions
      )
    ),
    "- GSE165004 slide-QC canonical files: 4",
    paste0(
      "- Missing source rows: ",
      nrow(
        missing_source_rows
      )
    ),
    "- Biological analysis: NOT STARTED",
    "- Scientific figures: NONE",
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
    "# STEP04 combined package status",
    "",
    paste0(
      "- Active run: ",
      run_id
    ),
    "- Dataset role lock: v2 CANDIDATE",
    "- GSE250130: ACQUISITION_AND_TAR_INVENTORY",
    "- GSE165004 slide-QC provenance: NARROWED_TO_4_CANONICAL_FILES",
    "- Five-dataset source lock: CANDIDATE_PENDING_MANUAL_REVIEW",
    "- Biological analysis: NOT STARTED",
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
cat("STEP04A v2 GSE250130 acquisition and source-candidate build completed.\n")
cat("GSE250130 role: TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT\n")
cat(
  "GSE250130 GEO samples parsed: ",
  length(
    GSE250130_sample_accessions
  ),
  "\n"
)
cat(
  "GSE250130 RAW tar present: ",
  file.exists(
    GSE250130_RAW_path
  ),
  "\n"
)
cat(
  "GSE250130 RAW tar bytes: ",
  if (
    file.exists(
      GSE250130_RAW_path
    )
  ) {
    as.numeric(
      file.info(
        GSE250130_RAW_path
      )$size
    )
  } else {
    NA
  },
  "\n"
)
cat(
  "GSE250130 tar member count: ",
  length(
    tar_member_list
  ),
  "\n"
)
cat("GSE165004 slide-QC canonical files: 4\n")
cat(
  "Combined source rows: ",
  nrow(
    source_lock_candidate
  ),
  "\n"
)
cat(
  "Missing source rows: ",
  nrow(
    missing_source_rows
  ),
  "\n"
)
cat("Biological analysis: NOT STARTED\n")
cat("Scientific figures: NONE\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n\n")
cat("Upload:\n")
cat("  06_results/STEP04A_v2_GATE_TABLE.csv\n")
cat("  06_results/STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv\n")
cat("  06_results/STEP04A_v2_DATASET_SOURCE_READINESS.csv\n")
cat("  06_results/STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv\n")
cat("  08_qc/STEP04A_v2_GSE250130_ACQUISITION_QC.csv\n")
cat("  08_qc/STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 2: 49_step04B_lock_inputs_extract_sc_and_audit_designs
# ==============================================================================
# ============================================================
# STEP04B - Lock five-dataset inputs, extract single-cell bundles,
# and audit sample designs for the combined downstream package
# Nutrition_Endometrial_VDR_v2
#
# Upstream accepted state:
# - STEP03 core locked
# - STEP04A v2 acquired GSE250130
# - 42/42 combined source rows present
# - GSE250130: 28 GEO samples and 28 outer-TAR members
#
# This script:
# - validates and locks the 42 source references without duplicating
#   the 3.6-GB GSE250130 TAR;
# - extracts the 28 GSE250130 sample bundles once into a short path;
# - validates 28 complete single-cell matrix triplets;
# - validates the 9 GSE183837 matrix triplets;
# - parses series-matrix metadata for all five datasets;
# - audits sample counts, metadata fields and shared-control overlap;
# - creates sample-design candidates for manual review;
# - does not run biological association analyses or create figures.
#
# Run from project root:
# source(
#   "04_code/R/49_step04B_lock_inputs_extract_sc_and_audit_designs.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(3600L, getOption("timeout")))

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04B_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

source_lock_version <- "COMBINED_INPUT_SOURCE_v1.0"
sample_design_candidate_version <- "COMBINED_SAMPLE_DESIGN_CANDIDATE_v1.0"

expected_source_rows <- c(
  GSE250130 = 2L,
  GSE165004 = 7L,
  GSE92324 = 2L,
  GSE111974 = 3L,
  GSE183837 = 28L
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
  control = file.path(step04_root, "00_module_control"),
  inputs = file.path(step04_root, "01_inputs"),
  code = file.path(step04_root, "02_code"),
  parameters = file.path(step04_root, "03_parameters"),
  environment = file.path(step04_root, "04_environment"),
  execution = file.path(step04_root, "05_execution"),
  results = file.path(step04_root, "06_results"),
  figures = file.path(step04_root, "07_figures_tables"),
  qc = file.path(step04_root, "08_qc"),
  text = file.path(step04_root, "09_manuscript_text"),
  archive_parent = file.path(step04_root, "10_archive_work")
)

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

locked_root <- file.path(
  dirs$results,
  "STEP04B_locked_combined_inputs",
  source_lock_version
)

design_root <- file.path(
  dirs$results,
  "STEP04B_sample_design_candidates"
)

dir.create(
  locked_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  design_root,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4B_",
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
      normalize_slashes(path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      normalize_slashes(path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
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

parse_GEO_series_matrix <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Series matrix is missing: ",
      path
    )
  }

  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  all_lines <- readLines(
    connection,
    warn = FALSE,
    encoding = "UTF-8"
  )

  sample_line_indices <- grep(
    "^!Sample_",
    all_lines
  )

  metadata_rows <- list()
  metadata_index <- 1L
  field_occurrence_counter <- list()

  for (line_index in sample_line_indices) {
    fields <- strsplit(
      all_lines[[line_index]],
      "\t",
      fixed = TRUE
    )[[1]]

    field_name <- sub(
      "^!Sample_",
      "",
      fields[[1]]
    )

    occurrence <- if (
      is.null(
        field_occurrence_counter[[field_name]]
      )
    ) {
      1L
    } else {
      field_occurrence_counter[[field_name]] + 1L
    }

    field_occurrence_counter[[field_name]] <-
      occurrence

    values <- fields[-1]

    values <- gsub(
      '^"|"$',
      "",
      values
    )

    metadata_rows[[metadata_index]] <- data.frame(
      field_name =
        field_name,
      field_occurrence =
        occurrence,
      sample_index =
        seq_along(values),
      value =
        values,
      stringsAsFactors = FALSE
    )

    metadata_index <- metadata_index + 1L
  }

  metadata_long <- if (
    length(metadata_rows) > 0L
  ) {
    do.call(
      rbind,
      metadata_rows
    )
  } else {
    data.frame(
      field_name = character(0),
      field_occurrence = integer(0),
      sample_index = integer(0),
      value = character(0),
      stringsAsFactors = FALSE
    )
  }

  table_begin <- grep(
    "^!series_matrix_table_begin",
    all_lines
  )

  table_end <- grep(
    "^!series_matrix_table_end",
    all_lines
  )

  expression_feature_n <- 0L
  expression_sample_headers <- character(0)

  if (
    length(table_begin) == 1L &&
      length(table_end) == 1L &&
      table_end > table_begin + 1L
  ) {
    expression_feature_n <-
      table_end - table_begin - 2L

    header_line <- all_lines[
      table_begin + 1L
    ]

    header_fields <- strsplit(
      header_line,
      "\t",
      fixed = TRUE
    )[[1]]

    header_fields <- gsub(
      '^"|"$',
      "",
      header_fields
    )

    expression_sample_headers <-
      header_fields[-1]
  }

  geo_accessions <- metadata_long$value[
    metadata_long$field_name ==
      "geo_accession"
  ]

  geo_accessions <- unique(
    geo_accessions[
      nzchar(geo_accessions)
    ]
  )

  list(
    metadata_long =
      metadata_long,
    sample_accessions =
      geo_accessions,
    sample_n =
      length(geo_accessions),
    expression_feature_n =
      expression_feature_n,
    expression_sample_headers =
      expression_sample_headers
  )
}

metadata_to_sample_sheet <- function(
    dataset_id,
    parsed_series
) {
  metadata_long <- parsed_series$metadata_long

  if (
    nrow(metadata_long) == 0L ||
      parsed_series$sample_n == 0L
  ) {
    return(
      data.frame()
    )
  }

  geo_rows <- metadata_long[
    metadata_long$field_name ==
      "geo_accession",
    ,
    drop = FALSE
  ]

  if (nrow(geo_rows) == 0L) {
    stop(
      "No GEO accession row was found for ",
      dataset_id
    )
  }

  sample_lookup <- data.frame(
    sample_index =
      geo_rows$sample_index,
    sample_accession =
      geo_rows$value,
    stringsAsFactors = FALSE
  )

  selected_fields <- metadata_long[
    metadata_long$field_name %in%
      c(
        "title",
        "source_name_ch1",
        "organism_ch1",
        "characteristics_ch1",
        "description",
        "relation"
      ),
    ,
    drop = FALSE
  ]

  selected_fields$field_key <- paste0(
    selected_fields$field_name,
    "__",
    selected_fields$field_occurrence
  )

  sample_sheet <- sample_lookup

  for (field_key in unique(
    selected_fields$field_key
  )) {
    field_rows <- selected_fields[
      selected_fields$field_key ==
        field_key,
      c(
        "sample_index",
        "value"
      ),
      drop = FALSE
    ]

    names(field_rows)[
      names(field_rows) ==
        "value"
    ] <- field_key

    sample_sheet <- merge(
      sample_sheet,
      field_rows,
      by = "sample_index",
      all.x = TRUE,
      sort = FALSE
    )
  }

  sample_sheet <- sample_sheet[
    order(
      sample_sheet$sample_index
    ),
    ,
    drop = FALSE
  ]

  sample_sheet$dataset_id <-
    dataset_id

  sample_sheet
}

characteristics_to_key_value <- function(
    dataset_id,
    sample_sheet
) {
  if (
    nrow(sample_sheet) == 0L
  ) {
    return(
      data.frame(
        dataset_id = character(0),
        sample_accession = character(0),
        characteristic_field = character(0),
        raw_value = character(0),
        parsed_key = character(0),
        parsed_value = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  characteristic_columns <- grep(
    "^characteristics_ch1__",
    names(sample_sheet),
    value = TRUE
  )

  rows <- list()
  row_index <- 1L

  for (column_name in characteristic_columns) {
    for (sample_index in seq_len(
      nrow(sample_sheet)
    )) {
      raw_value <- as.character(
        sample_sheet[
          sample_index,
          column_name
        ]
      )

      if (
        is.na(raw_value) ||
          !nzchar(
            trimws(raw_value)
          )
      ) {
        next
      }

      colon_position <- regexpr(
        ":",
        raw_value,
        fixed = TRUE
      )[[1]]

      if (colon_position > 0L) {
        parsed_key <- trimws(
          substr(
            raw_value,
            1L,
            colon_position - 1L
          )
        )

        parsed_value <- trimws(
          substr(
            raw_value,
            colon_position + 1L,
            nchar(raw_value)
          )
        )
      } else {
        parsed_key <- NA_character_
        parsed_value <- raw_value
      }

      rows[[row_index]] <- data.frame(
        dataset_id =
          dataset_id,
        sample_accession =
          sample_sheet$sample_accession[
            sample_index
          ],
        characteristic_field =
          column_name,
        raw_value =
          raw_value,
        parsed_key =
          parsed_key,
        parsed_value =
          parsed_value,
        stringsAsFactors = FALSE
      )

      row_index <- row_index + 1L
    }
  }

  if (length(rows) == 0L) {
    return(
      data.frame(
        dataset_id = character(0),
        sample_accession = character(0),
        characteristic_field = character(0),
        raw_value = character(0),
        parsed_key = character(0),
        parsed_value = character(0),
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(
    rbind,
    rows
  )
}

identify_triplet_component <- function(path) {
  file_name <- tolower(
    basename(path)
  )

  if (grepl(
    "matrix.*\\.mtx(\\.gz)?$|counts.*\\.mtx(\\.gz)?$",
    file_name,
    perl = TRUE
  )) {
    return("matrix")
  }

  if (grepl(
    "barcodes.*\\.(tsv|txt)(\\.gz)?$",
    file_name,
    perl = TRUE
  )) {
    return("barcodes")
  }

  if (grepl(
    "(features|genes).*\\.(tsv|txt)(\\.gz)?$",
    file_name,
    perl = TRUE
  )) {
    return("features_or_genes")
  }

  "other"
}

extract_GSM <- function(path) {
  match_value <- regmatches(
    basename(path),
    regexpr(
      "GSM[0-9]+",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  )

  if (
    length(match_value) == 0L ||
      !nzchar(match_value)
  ) {
    return(NA_character_)
  }

  toupper(match_value)
}

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Validate STEP04A v2 outputs
# ----------------------------------------------------------------------

step04A_gate_path <- file.path(
  dirs$results,
  "STEP04A_v2_GATE_TABLE.csv"
)

source_candidate_path <- file.path(
  dirs$results,
  "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv"
)

role_lock_path <- file.path(
  dirs$results,
  "STEP04A_v2_DATASET_ROLE_LOCK.csv"
)

slide_lock_path <- file.path(
  dirs$results,
  "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv"
)

required_upstream <- c(
  step04A_gate_path,
  source_candidate_path,
  role_lock_path,
  slide_lock_path
)

missing_upstream <- required_upstream[
  !file.exists(required_upstream)
]

if (length(missing_upstream) > 0L) {
  stop(
    "Required STEP04A v2 output is missing:\n",
    paste(
      missing_upstream,
      collapse = "\n"
    )
  )
}

step04A_gate <- read.csv(
  step04A_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_candidate <- read.csv(
  source_candidate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

role_lock <- read.csv(
  role_lock_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

slide_lock <- read.csv(
  slide_lock_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

accepted_gate <- step04A_gate$status[
  step04A_gate$gate_component ==
    "STEP04A v2 gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "PROCEED_TO_MANUAL_STEP04_INPUT_SOURCE_LOCK"
) {
  stop(
    "STEP04A v2 does not authorize source lock."
  )
}

if (
  nrow(source_candidate) !=
    sum(expected_source_rows) ||
    any(
      !source_candidate$file_exists
    ) ||
    any(
      is.na(
        source_candidate$md5
      ) |
        is.na(
          source_candidate$sha256
        )
    )
) {
  stop(
    "The 42-row source candidate is incomplete."
  )
}

observed_source_rows <- table(
  source_candidate$dataset_id
)

if (
  !setequal(
    names(observed_source_rows),
    names(expected_source_rows)
  ) ||
    !all(
      as.integer(
        observed_source_rows[
          names(expected_source_rows)
        ]
      ) ==
        as.integer(
          expected_source_rows
        )
    )
) {
  stop(
    "Dataset-specific source-row counts differ from 2/7/2/3/28."
  )
}

if (
  nrow(role_lock) != 5L ||
    role_lock$analysis_role[
      role_lock$dataset_id ==
        "GSE250130"
    ] !=
      "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT"
) {
  stop(
    "The corrected five-dataset role lock failed."
  )
}

if (
  nrow(slide_lock) != 4L ||
    any(
      !slide_lock$file_exists
    )
) {
  stop(
    "The four-file GSE165004 slide-QC lock failed."
  )
}

# ----------------------------------------------------------------------
# Revalidate source references without rehashing the 3.6-GB TAR
# ----------------------------------------------------------------------

source_candidate$current_exists <- file.exists(
  source_candidate$path
)

source_candidate$current_size_bytes <- as.numeric(
  file.info(
    source_candidate$path
  )$size
)

source_candidate$size_identity_pass <-
  source_candidate$current_exists &
  source_candidate$current_size_bytes ==
    source_candidate$size_bytes

source_candidate$hash_recheck_mode <- ifelse(
  source_candidate$size_bytes <=
    500 * 1024^2,
  "MD5_SHA256_RECHECKED",
  "LARGE_FILE_SIZE_AND_ARCHIVE_STRUCTURE_RECHECKED"
)

source_candidate$observed_md5 <- NA_character_
source_candidate$observed_sha256 <- NA_character_

small_file_index <-
  source_candidate$size_bytes <=
    500 * 1024^2

source_candidate$observed_md5[
  small_file_index
] <- vapply(
  source_candidate$path[
    small_file_index
  ],
  md5_file,
  character(1)
)

source_candidate$observed_sha256[
  small_file_index
] <- vapply(
  source_candidate$path[
    small_file_index
  ],
  sha256_file,
  character(1)
)

source_candidate$hash_identity_pass <-
  TRUE

source_candidate$hash_identity_pass[
  small_file_index
] <- with(
  source_candidate[
    small_file_index,
    ,
    drop = FALSE
  ],
  md5 ==
    observed_md5 &
  sha256 ==
    observed_sha256
)

if (
  !all(
    source_candidate$size_identity_pass
  ) ||
    !all(
      source_candidate$hash_identity_pass
    )
) {
  stop(
    "At least one combined source reference failed identity validation."
  )
}

write_csv_utf8(
  source_candidate,
  file.path(
    dirs$qc,
    "STEP04B_SOURCE_IDENTITY_RECHECK.csv"
  )
)

# ----------------------------------------------------------------------
# Extract and validate GSE250130 sample bundles
# ----------------------------------------------------------------------

GSE250130_RAW_path <- source_candidate$path[
  source_candidate$dataset_id ==
    "GSE250130" &
    source_candidate$source_role ==
      "RAW_OUTER_TAR_OF_SAMPLE_BUNDLES"
]

if (length(GSE250130_RAW_path) != 1L) {
  stop(
    "GSE250130 outer TAR could not be uniquely identified."
  )
}

GSE250130_extract_root <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE250130",
  "raw_counts",
  "outer_extracted"
)

GSE250130_sample_root <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE250130",
  "raw_counts",
  "mm_by_sample"
)

dir.create(
  GSE250130_extract_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  GSE250130_sample_root,
  recursive = TRUE,
  showWarnings = FALSE
)

outer_marker <- file.path(
  GSE250130_extract_root,
  "_OUTER_EXTRACTION_COMPLETE.txt"
)

if (!file.exists(outer_marker)) {
  log_event(
    "Extracting GSE250130 outer TAR."
  )

  extraction_status <- tryCatch(
    {
      utils::untar(
        GSE250130_RAW_path,
        exdir =
          GSE250130_extract_root
      )

      TRUE
    },
    error = function(error_object) {
      stop(
        "GSE250130 outer TAR extraction failed: ",
        conditionMessage(error_object)
      )
    }
  )

  if (isTRUE(extraction_status)) {
    write_text(
      c(
        paste0(
          "Run ID: ",
          run_id
        ),
        paste0(
          "Source: ",
          normalize_slashes(
            GSE250130_RAW_path
          )
        ),
        paste0(
          "Completed: ",
          format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S"
          )
        )
      ),
      outer_marker
    )
  }
}

nested_bundles <- list.files(
  GSE250130_extract_root,
  pattern =
    "\\.(tar|tar\\.gz|tgz)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(nested_bundles) != 28L) {
  stop(
    "Expected 28 GSE250130 nested sample bundles but found ",
    length(nested_bundles),
    "."
  )
}

GSE250130_bundle_rows <- list()
bundle_row_index <- 1L

for (bundle_path in sort(
  nested_bundles
)) {
  sample_id <- extract_GSM(
    bundle_path
  )

  if (is.na(sample_id)) {
    stop(
      "No GSM accession was found in GSE250130 bundle name: ",
      basename(bundle_path)
    )
  }

  sample_dir <- file.path(
    GSE250130_sample_root,
    sample_id
  )

  marker_path <- file.path(
    sample_dir,
    "_EXTRACTION_COMPLETE.txt"
  )

  dir.create(
    sample_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!file.exists(marker_path)) {
    log_event(
      "Extracting",
      sample_id
    )

    tryCatch(
      {
        utils::untar(
          bundle_path,
          exdir = sample_dir
        )
      },
      error = function(error_object) {
        stop(
          "Nested extraction failed for ",
          sample_id,
          ": ",
          conditionMessage(error_object)
        )
      }
    )

    write_text(
      c(
        paste0(
          "Run ID: ",
          run_id
        ),
        paste0(
          "Bundle: ",
          normalize_slashes(
            bundle_path
          )
        )
      ),
      marker_path
    )
  }

  extracted_files <- list.files(
    sample_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )

  components <- vapply(
    extracted_files,
    identify_triplet_component,
    character(1)
  )

  matrix_files <- extracted_files[
    components == "matrix"
  ]

  barcode_files <- extracted_files[
    components == "barcodes"
  ]

  feature_files <- extracted_files[
    components ==
      "features_or_genes"
  ]

  GSE250130_bundle_rows[[bundle_row_index]] <- data.frame(
    dataset_id =
      "GSE250130",
    sample_accession =
      sample_id,
    nested_bundle_path =
      normalize_slashes(
        bundle_path
      ),
    matrix_file = if (
      length(matrix_files) == 1L
    ) {
      normalize_slashes(
        matrix_files
      )
    } else {
      paste(
        normalize_slashes(
          matrix_files
        ),
        collapse = ";"
      )
    },
    barcode_file = if (
      length(barcode_files) == 1L
    ) {
      normalize_slashes(
        barcode_files
      )
    } else {
      paste(
        normalize_slashes(
          barcode_files
        ),
        collapse = ";"
      )
    },
    feature_file = if (
      length(feature_files) == 1L
    ) {
      normalize_slashes(
        feature_files
      )
    } else {
      paste(
        normalize_slashes(
          feature_files
        ),
        collapse = ";"
      )
    },
    matrix_file_n =
      length(matrix_files),
    barcode_file_n =
      length(barcode_files),
    feature_file_n =
      length(feature_files),
    complete_triplet =
      length(matrix_files) == 1L &&
      length(barcode_files) == 1L &&
      length(feature_files) == 1L,
    stringsAsFactors = FALSE
  )

  bundle_row_index <- bundle_row_index + 1L
}

GSE250130_bundle_map <- do.call(
  rbind,
  GSE250130_bundle_rows
)

if (
  nrow(GSE250130_bundle_map) != 28L ||
    anyDuplicated(
      GSE250130_bundle_map$sample_accession
    ) > 0L ||
    !all(
      GSE250130_bundle_map$complete_triplet
    )
) {
  stop(
    "GSE250130 did not yield 28 unique complete matrix triplets."
  )
}

write_csv_utf8(
  GSE250130_bundle_map,
  file.path(
    dirs$results,
    "STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv"
  )
)

# ----------------------------------------------------------------------
# Validate GSE183837 nine donor triplets
# ----------------------------------------------------------------------

GSE183837_triplet_paths <- source_candidate$path[
  source_candidate$dataset_id ==
    "GSE183837" &
    source_candidate$source_role ==
      "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER"
]

if (length(GSE183837_triplet_paths) != 27L) {
  stop(
    "Expected 27 locked GSE183837 triplet members."
  )
}

GSE183837_triplet_inventory <- data.frame(
  dataset_id =
    "GSE183837",
  sample_accession = vapply(
    GSE183837_triplet_paths,
    extract_GSM,
    character(1)
  ),
  component = vapply(
    GSE183837_triplet_paths,
    identify_triplet_component,
    character(1)
  ),
  path =
    normalize_slashes(
      GSE183837_triplet_paths
    ),
  stringsAsFactors = FALSE
)

if (
  any(is.na(
    GSE183837_triplet_inventory$sample_accession
  )) ||
    any(
      GSE183837_triplet_inventory$component ==
        "other"
    )
) {
  stop(
    "GSE183837 triplet filename classification failed."
  )
}

GSE183837_bundle_summary <- aggregate(
  component ~
    sample_accession,
  data =
    GSE183837_triplet_inventory,
  FUN = function(x) {
    paste(
      sort(x),
      collapse = ";"
    )
  }
)

GSE183837_component_count <- aggregate(
  component ~
    sample_accession,
  data =
    GSE183837_triplet_inventory,
  FUN = length
)

names(
  GSE183837_component_count
)[
  names(
    GSE183837_component_count
  ) == "component"
] <- "component_n"

GSE183837_bundle_summary <- merge(
  GSE183837_bundle_summary,
  GSE183837_component_count,
  by = "sample_accession",
  all = TRUE,
  sort = FALSE
)

GSE183837_bundle_summary$complete_triplet <-
  GSE183837_bundle_summary$component_n ==
    3L &
  GSE183837_bundle_summary$component ==
    paste(
      sort(
        c(
          "barcodes",
          "features_or_genes",
          "matrix"
        )
      ),
      collapse = ";"
    )

if (
  nrow(GSE183837_bundle_summary) != 9L ||
    !all(
      GSE183837_bundle_summary$complete_triplet
    )
) {
  stop(
    "GSE183837 did not yield nine complete sample triplets."
  )
}

write_csv_utf8(
  GSE183837_triplet_inventory,
  file.path(
    dirs$qc,
    "STEP04B_GSE183837_TRIPLET_INVENTORY.csv"
  )
)

write_csv_utf8(
  GSE183837_bundle_summary,
  file.path(
    dirs$results,
    "STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv"
  )
)

# ----------------------------------------------------------------------
# Parse five series matrices and build sample-design candidates
# ----------------------------------------------------------------------

series_matrix_rows <- source_candidate[
  grepl(
    "SERIES_MATRIX",
    source_candidate$source_role
  ),
  ,
  drop = FALSE
]

dataset_ids <- c(
  "GSE250130",
  "GSE165004",
  "GSE92324",
  "GSE111974",
  "GSE183837"
)

series_audit_rows <- list()
sample_sheet_list <- list()
characteristics_list <- list()

for (dataset_id in dataset_ids) {
  dataset_series_rows <- series_matrix_rows[
    series_matrix_rows$dataset_id ==
      dataset_id,
    ,
    drop = FALSE
  ]

  if (nrow(dataset_series_rows) != 1L) {
    stop(
      "Expected exactly one series matrix for ",
      dataset_id,
      " but found ",
      nrow(dataset_series_rows),
      "."
    )
  }

  series_path <- dataset_series_rows$path[[1]]

  parsed_series <- parse_GEO_series_matrix(
    series_path
  )

  sample_sheet <- metadata_to_sample_sheet(
    dataset_id,
    parsed_series
  )

  characteristic_long <-
    characteristics_to_key_value(
      dataset_id,
      sample_sheet
    )

  sample_sheet_list[[dataset_id]] <-
    sample_sheet

  characteristics_list[[dataset_id]] <-
    characteristic_long

  series_audit_rows[[dataset_id]] <- data.frame(
    dataset_id =
      dataset_id,
    series_matrix_path =
      normalize_slashes(
        series_path
      ),
    parsed_sample_n =
      parsed_series$sample_n,
    expression_feature_n =
      parsed_series$expression_feature_n,
    expression_header_sample_n =
      length(
        parsed_series$expression_sample_headers
      ),
    metadata_field_n =
      length(
        unique(
          parsed_series$metadata_long$field_name
        )
      ),
    characteristic_key_n =
      length(
        unique(
          characteristic_long$parsed_key[
            !is.na(
              characteristic_long$parsed_key
            )
          ]
        )
      ),
    sample_design_status =
      "CANDIDATE_PENDING_MANUAL_LOCK",
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    parsed_series$metadata_long,
    file.path(
      dirs$qc,
      paste0(
        "STEP04B_",
        dataset_id,
        "_SERIES_METADATA_LONG.csv"
      )
    )
  )

  write_csv_utf8(
    sample_sheet,
    file.path(
      design_root,
      paste0(
        "STEP04B_",
        dataset_id,
        "_SAMPLE_SHEET_CANDIDATE.csv"
      )
    )
  )

  write_csv_utf8(
    characteristic_long,
    file.path(
      design_root,
      paste0(
        "STEP04B_",
        dataset_id,
        "_CHARACTERISTICS_KEY_VALUE.csv"
      )
    )
  )
}

series_audit <- do.call(
  rbind,
  series_audit_rows
)

write_csv_utf8(
  series_audit,
  file.path(
    dirs$results,
    "STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv"
  )
)

# ----------------------------------------------------------------------
# Cross-check single-cell bundle accessions against GEO metadata
# ----------------------------------------------------------------------

GSE250130_geo_samples <-
  sample_sheet_list[[
    "GSE250130"
  ]]$sample_accession

GSE183837_geo_samples <-
  sample_sheet_list[[
    "GSE183837"
  ]]$sample_accession

GSE250130_accession_pass <- setequal(
  GSE250130_geo_samples,
  GSE250130_bundle_map$sample_accession
)

GSE183837_accession_pass <- setequal(
  GSE183837_geo_samples,
  GSE183837_bundle_summary$sample_accession
)

if (
  !GSE250130_accession_pass ||
    !GSE183837_accession_pass
) {
  stop(
    "Single-cell sample-bundle accessions do not reconcile with GEO metadata."
  )
}

# ----------------------------------------------------------------------
# Shared-control and title-overlap audit for GSE165004/GSE111974
# ----------------------------------------------------------------------

sheet_165004 <- sample_sheet_list[[
  "GSE165004"
]]

sheet_111974 <- sample_sheet_list[[
  "GSE111974"
]]

title_column_165004 <- grep(
  "^title__",
  names(sheet_165004),
  value = TRUE
)

title_column_111974 <- grep(
  "^title__",
  names(sheet_111974),
  value = TRUE
)

normalize_title <- function(x) {
  tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      as.character(x),
      perl = TRUE
    )
  )
}

shared_GSM <- intersect(
  sheet_165004$sample_accession,
  sheet_111974$sample_accession
)

shared_title_rows <- data.frame(
  GSE165004_sample = character(0),
  GSE111974_sample = character(0),
  normalized_title = character(0),
  stringsAsFactors = FALSE
)

if (
  length(title_column_165004) >= 1L &&
    length(title_column_111974) >= 1L
) {
  titles_165004 <- data.frame(
    sample_accession =
      sheet_165004$sample_accession,
    normalized_title =
      normalize_title(
        sheet_165004[[
          title_column_165004[[1]]
        ]]
      ),
    stringsAsFactors = FALSE
  )

  titles_111974 <- data.frame(
    sample_accession =
      sheet_111974$sample_accession,
    normalized_title =
      normalize_title(
        sheet_111974[[
          title_column_111974[[1]]
        ]]
      ),
    stringsAsFactors = FALSE
  )

  shared_title_rows <- merge(
    titles_165004,
    titles_111974,
    by = "normalized_title",
    suffixes = c(
      "_GSE165004",
      "_GSE111974"
    ),
    all = FALSE
  )

  names(shared_title_rows) <- c(
    "normalized_title",
    "GSE165004_sample",
    "GSE111974_sample"
  )

  shared_title_rows <- shared_title_rows[
    nzchar(
      shared_title_rows$normalized_title
    ),
    ,
    drop = FALSE
  ]
}

shared_control_audit <- data.frame(
  metric = c(
    "Exact shared GSM accessions",
    "Normalized shared sample titles",
    "GSE165004 parsed samples",
    "GSE111974 parsed samples",
    "Shared-control design status"
  ),
  value = c(
    length(shared_GSM),
    nrow(shared_title_rows),
    nrow(sheet_165004),
    nrow(sheet_111974),
    "PENDING_MANUAL_CONTROL_SHARING_LOCK"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  shared_control_audit,
  file.path(
    dirs$results,
    "STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv"
  )
)

write_csv_utf8(
  shared_title_rows,
  file.path(
    dirs$qc,
    "STEP04B_GSE165004_GSE111974_SHARED_TITLE_ROWS.csv"
  )
)

write_text(
  shared_GSM,
  file.path(
    dirs$qc,
    "STEP04B_GSE165004_GSE111974_SHARED_GSM.txt"
  )
)

# ----------------------------------------------------------------------
# Lock source references and derived bundle maps
# ----------------------------------------------------------------------

locked_source_reference_path <- file.path(
  locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_REFERENCE_LOCK.csv"
)

write_csv_utf8(
  source_candidate,
  locked_source_reference_path
)

locked_role_path <- file.path(
  locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_DATASET_ROLE_LOCK.csv"
)

write_csv_utf8(
  role_lock,
  locked_role_path
)

locked_GSE250130_bundle_path <- file.path(
  locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_GSE250130_BUNDLE_MAP.csv"
)

write_csv_utf8(
  GSE250130_bundle_map,
  locked_GSE250130_bundle_path
)

locked_GSE183837_bundle_path <- file.path(
  locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_GSE183837_BUNDLE_MAP.csv"
)

write_csv_utf8(
  GSE183837_bundle_summary,
  locked_GSE183837_bundle_path
)

locked_series_audit_path <- file.path(
  locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_SERIES_AUDIT.csv"
)

write_csv_utf8(
  series_audit,
  locked_series_audit_path
)

locked_paths <- c(
  locked_source_reference_path,
  locked_role_path,
  locked_GSE250130_bundle_path,
  locked_GSE183837_bundle_path,
  locked_series_audit_path
)

locked_manifest <- data.frame(
  artifact_id = c(
    "COMBINED_SOURCE_REFERENCE_LOCK",
    "DATASET_ROLE_LOCK",
    "GSE250130_SAMPLE_BUNDLE_MAP",
    "GSE183837_SAMPLE_BUNDLE_MAP",
    "SERIES_MATRIX_SAMPLE_AUDIT"
  ),
  lock_version =
    source_lock_version,
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
  lock_status =
    "LOCKED",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  locked_manifest,
  file.path(
    locked_root,
    "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv"
  )
)

# ----------------------------------------------------------------------
# Metadata-field inventory and Gate
# ----------------------------------------------------------------------

metadata_field_inventory <- do.call(
  rbind,
  lapply(
    dataset_ids,
    function(dataset_id) {
      characteristic_long <-
        characteristics_list[[
          dataset_id
        ]]

      if (
        nrow(characteristic_long) == 0L
      ) {
        return(
          data.frame(
            dataset_id =
              dataset_id,
            parsed_key =
              NA_character_,
            nonmissing_value_n =
              0L,
            unique_value_n =
              0L,
            example_values =
              NA_character_,
            stringsAsFactors = FALSE
          )
        )
      }

      split_rows <- split(
        characteristic_long,
        characteristic_long$parsed_key,
        drop = TRUE
      )

      do.call(
        rbind,
        lapply(
          names(split_rows),
          function(key_name) {
            key_rows <- split_rows[[
              key_name
            ]]

            values <- unique(
              key_rows$parsed_value[
                !is.na(
                  key_rows$parsed_value
                ) &
                  nzchar(
                    trimws(
                      key_rows$parsed_value
                    )
                  )
              ]
            )

            data.frame(
              dataset_id =
                dataset_id,
              parsed_key =
                key_name,
              nonmissing_value_n =
                sum(
                  !is.na(
                    key_rows$parsed_value
                  ) &
                    nzchar(
                      trimws(
                        key_rows$parsed_value
                      )
                    )
                ),
              unique_value_n =
                length(values),
              example_values =
                paste(
                  head(values, 8L),
                  collapse = ";"
                ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

write_csv_utf8(
  metadata_field_inventory,
  file.path(
    dirs$qc,
    "STEP04B_METADATA_FIELD_INVENTORY.csv"
  )
)

gate_status <-
  "PROCEED_TO_MANUAL_SAMPLE_DESIGN_LOCK_AND_COMBINED_ANALYSIS_BUILD"

gate_table <- data.frame(
  gate_component = c(
    "STEP04A v2 gate",
    "Combined source rows",
    "Source identity recheck",
    "GSE250130 outer extraction",
    "GSE250130 complete triplets",
    "GSE250130 GEO/bundle accession reconciliation",
    "GSE183837 complete triplets",
    "GSE183837 GEO/bundle accession reconciliation",
    "Five series matrices parsed",
    "GSE165004/GSE111974 shared-control audit",
    "Combined source lock",
    "Sample-design candidates",
    "Biological analysis",
    "Scientific figures",
    "STEP04B gate"
  ),
  status = c(
    "PASS",
    "LOCKED_42",
    "PASS",
    "PASS_28_NESTED_BUNDLES",
    "PASS_28_OF_28",
    "PASS",
    "PASS_9_OF_9",
    "PASS",
    "PASS_5_OF_5",
    paste0(
      "EXACT_GSM_",
      length(shared_GSM),
      "_NORMALIZED_TITLE_",
      nrow(shared_title_rows)
    ),
    paste0(
      "LOCKED_",
      source_lock_version
    ),
    paste0(
      "CREATED_",
      sample_design_candidate_version
    ),
    "NOT_STARTED",
    "NONE_CREATED",
    gate_status
  ),
  required_action = c(
    "Retain the acquisition audit.",
    "Do not change source paths or hashes without a new version.",
    "Use the locked source-reference table.",
    "Do not re-extract unless the extraction markers are removed intentionally.",
    "Use donor-aware single-cell analysis.",
    "Use the 28-sample GEO design.",
    "Use nine donor-level sample bundles.",
    "Use the nine-sample GEO design.",
    "Review group, timepoint, fertility/RIF and tissue fields.",
    "Lock shared-control handling before external validation.",
    "Do not duplicate the 3.6-GB TAR inside results or archives.",
    "Manually approve sample/group assignments before analysis.",
    "Build one coordinated biological-analysis package.",
    "Render only after biological QC.",
    "Upload the sample audits, metadata inventory, shared-control audit and Gate."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP04B_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP04B source-lock and sample-design audit card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Combined input source: ",
      source_lock_version,
      " LOCKED"
    ),
    "- Source references: 42",
    "- GSE250130 complete sample triplets: 28/28",
    "- GSE183837 complete sample triplets: 9/9",
    paste0(
      "- GSE165004/GSE111974 exact shared GSM: ",
      length(shared_GSM)
    ),
    paste0(
      "- GSE165004/GSE111974 normalized shared titles: ",
      nrow(shared_title_rows)
    ),
    "- Five series matrices parsed: YES",
    "- Sample-design candidates: CREATED",
    "- Biological analysis: NOT STARTED",
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "STEP04B_SOURCE_LOCK_SAMPLE_DESIGN_AUDIT_CARD.md"
  )
)

write_text(
  c(
    "# STEP04B Methods note",
    "",
    "The five-dataset source set was locked by reference rather than by duplicating large raw files. All 42 source rows were revalidated by path and file size, with MD5 and SHA-256 recomputed for files no larger than 500 MB. The 3.6-GB GSE250130 outer TAR retained its acquisition hashes and was additionally validated by archive structure and complete extraction.",
    "",
    "The 28 GSE250130 nested bundles and nine GSE183837 sample triplets were reconciled with GEO sample accessions. GEO sample metadata, titles and characteristics were parsed into dataset-specific sample-design candidates. No biological association analysis or scientific figure was generated."
  ),
  file.path(
    dirs$text,
    "STEP04B_Methods_note.md"
  )
)

parameters <- list(
  schema_version = "1.0",
  run_id = run_id,
  source_lock_version =
    source_lock_version,
  source_row_n =
    nrow(source_candidate),
  GSE250130_complete_triplet_n =
    nrow(GSE250130_bundle_map),
  GSE183837_complete_triplet_n =
    nrow(GSE183837_bundle_summary),
  parsed_series_matrix_n =
    nrow(series_audit),
  exact_shared_GSM_n =
    length(shared_GSM),
  normalized_shared_title_n =
    nrow(shared_title_rows),
  biological_analysis_started =
    FALSE,
  scientific_figures_created =
    FALSE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP04B_parameters.json"
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
    "sessionInfo_POST_RUN_STEP04B.txt"
  )
)

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  copy_if_different(
    source_file,
    file.path(
      dirs$code,
      "49_step04B_lock_inputs_extract_sc_and_audit_designs.R"
    )
  )
}

# ----------------------------------------------------------------------
# Short archive: derived metadata and manifests only; no raw data copy
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP04B_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP04B_SOURCE_LOCK_SAMPLE_DESIGN_AUDIT_CARD.md"
  ),
  file.path(
    dirs$results,
    "STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv"
  ),
  file.path(
    dirs$results,
    "STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv"
  ),
  file.path(
    dirs$results,
    "STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv"
  ),
  file.path(
    dirs$results,
    "STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04B_SOURCE_IDENTITY_RECHECK.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04B_METADATA_FIELD_INVENTORY.csv"
  ),
  file.path(
    locked_root,
    "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv"
  ),
  locked_source_reference_path,
  locked_role_path,
  locked_GSE250130_bundle_path,
  locked_GSE183837_bundle_path,
  locked_series_audit_path,
  file.path(
    dirs$text,
    "STEP04B_Methods_note.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04B_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04B.txt"
  ),
  file.path(
    dirs$code,
    "49_step04B_lock_inputs_extract_sc_and_audit_designs.R"
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
      "Methods",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "QC|AUDIT|CARD|MANIFEST|LOCK|MAP|SUMMARY|INVENTORY",
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
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP04B_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04B archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Source lock: ",
      source_lock_version
    ),
    "- Source rows: 42",
    "- GSE250130 complete triplets: 28",
    "- GSE183837 complete triplets: 9",
    "- Raw files duplicated into archive: NO",
    "- Biological analysis: NOT STARTED",
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
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- Combined input source: ",
      source_lock_version,
      " LOCKED"
    ),
    "- GSE250130: 28 complete sample bundles",
    "- GSE183837: 9 complete sample bundles",
    "- Five series matrices: PARSED",
    "- Sample designs: CANDIDATE_PENDING_MANUAL_LOCK",
    "- Biological analysis: NOT STARTED",
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
cat("STEP04B combined input lock and sample-design audit completed.\n")
cat(
  "Combined input source: ",
  source_lock_version,
  " LOCKED\n"
)
cat("Source rows: 42\n")
cat("GSE250130 complete sample triplets: 28/28\n")
cat("GSE183837 complete sample triplets: 9/9\n")
cat(
  "GSE165004/GSE111974 exact shared GSM: ",
  length(shared_GSM),
  "\n"
)
cat(
  "GSE165004/GSE111974 normalized shared titles: ",
  nrow(shared_title_rows),
  "\n"
)
cat("Five series matrices parsed: YES\n")
cat("Sample-design candidates: CREATED\n")
cat("Biological analysis: NOT STARTED\n")
cat("Scientific figures: NONE\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n\n")
cat("Upload:\n")
cat("  06_results/STEP04B_GATE_TABLE.csv\n")
cat("  06_results/STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv\n")
cat("  06_results/STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv\n")
cat("  06_results/STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv\n")
cat("  06_results/STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv\n")
cat("  08_qc/STEP04B_METADATA_FIELD_INVENTORY.csv\n")
cat("  06_results/STEP04B_sample_design_candidates/\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 3: 50_step04C_v3_lock_sample_designs_and_platform_preflight
# ==============================================================================
# ============================================================
# STEP04C - Lock five-dataset sample designs, verify shared
# controls, derive GSE165004 slide IDs, and acquire platform
# annotations for the combined biological-analysis build
# Nutrition_Endometrial_VDR_v2
#
# Key scientific decisions:
# - GSE250130 temporal analysis uses 18 fertile samples across
#   LH+3/+5/+7/+9/+11; LH+7 contains six fertile samples
#   (three LH7 and three LH7_PGT). The 10 RIF samples are compared
#   only with the six fertile LH+7 samples.
# - GSE165004 is the primary patient cohort:
#   Control 24, RPL 24, UIF 24, with slide-aware sensitivity.
# - GSE92324 is an independent stimulated-cycle validation cohort:
#   control 8, implantation failure 10; HCG+6/+7 is retained.
# - GSE111974 contains RIF 24 plus the same 24 control arrays used
#   in GSE165004. It is a shared-control disease-context extension,
#   not an independent external cohort.
# - GSE183837 provides independent single-cell localization:
#   control 3, RIF 6.
#
# This script:
# - validates STEP04B source lock and sample candidates;
# - locks all five sample-design tables;
# - derives GSE165004 slide IDs from raw Agilent filenames;
# - confirms GSE111974/GSE165004 control-array identity by title,
#   array identifier and decompressed file hash;
# - downloads GPL16699 and GPL17077 annotation files if absent;
# - audits the GSE92324 group/day design matrix;
# - authorizes one combined bulk and single-cell analysis build.
#
# This script does NOT:
# - run patient association models;
# - run single-cell QC/clustering;
# - select final genes;
# - create scientific figures.
#
# Run from project root:
# source(
#   "04_code/R/50_step04C_v3_lock_sample_designs_and_platform_preflight.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(3600L, getOption("timeout")))

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04C_v3_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

sample_design_lock_version <- "COMBINED_SAMPLE_DESIGN_v1.0"
platform_source_version <- "BULK_PLATFORM_ANNOTATION_SOURCE_v1.0"

GPL16699_self_url <-
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16699&targ=self&view=full&form=text"

GPL17077_self_url <-
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17077&targ=self&view=full&form=text"

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
  control = file.path(step04_root, "00_module_control"),
  inputs = file.path(step04_root, "01_inputs"),
  code = file.path(step04_root, "02_code"),
  parameters = file.path(step04_root, "03_parameters"),
  environment = file.path(step04_root, "04_environment"),
  execution = file.path(step04_root, "05_execution"),
  results = file.path(step04_root, "06_results"),
  figures = file.path(step04_root, "07_figures_tables"),
  qc = file.path(step04_root, "08_qc"),
  text = file.path(step04_root, "09_manuscript_text"),
  archive_parent = file.path(step04_root, "10_archive_work")
)

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

locked_root <- file.path(
  dirs$results,
  "STEP04C_locked_sample_designs",
  sample_design_lock_version
)

platform_root <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "platform_annotations"
)

dir.create(
  locked_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  platform_root,
  recursive = TRUE,
  showWarnings = FALSE
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4C_",
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
      normalize_slashes(path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      normalize_slashes(path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
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

try_download_with_retry <- function(
    url,
    destination
) {
  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (
    file.exists(destination) &&
      as.numeric(
        file.info(destination)$size
      ) > 10000
  ) {
    return(TRUE)
  }

  if (file.exists(destination)) {
    unlink(
      destination,
      force = TRUE
    )
  }

  curl_candidates <- c(
    Sys.which("curl.exe"),
    Sys.which("curl")
  )

  curl_candidates <- curl_candidates[
    nzchar(curl_candidates)
  ]

  if (length(curl_candidates) > 0L) {
    status <- suppressWarnings(
      system2(
        curl_candidates[[1]],
        args = c(
          "--fail",
          "--location",
          "--retry",
          "10",
          "--retry-all-errors",
          "--retry-delay",
          "5",
          "--connect-timeout",
          "60",
          "--output",
          shQuote(destination),
          shQuote(url)
        )
      )
    )

    return(
      identical(
        status,
        0L
      ) &&
        file.exists(destination) &&
        as.numeric(
          file.info(destination)$size
        ) > 10000
    )
  }

  status <- tryCatch(
    {
      download.file(
        url,
        destination,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )

      TRUE
    },
    error = function(error_object) {
      FALSE
    }
  )

  isTRUE(status) &&
    file.exists(destination) &&
    as.numeric(
      file.info(destination)$size
    ) > 10000
}

read_platform_soft_header <- function(path) {
  if (!file.exists(path)) {
    return(
      data.frame(
        read_status = "MISSING",
        n_columns = NA_integer_,
        first_columns = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  table_begin_found <- FALSE
  header_line <- NA_character_
  line_counter <- 0L

  repeat {
    line_value <- readLines(
      connection,
      n = 1L,
      warn = FALSE,
      encoding = "UTF-8"
    )

    if (length(line_value) == 0L) {
      break
    }

    line_counter <- line_counter + 1L

    if (
      identical(
        trimws(line_value),
        "!platform_table_begin"
      )
    ) {
      table_begin_found <- TRUE

      header_value <- readLines(
        connection,
        n = 1L,
        warn = FALSE,
        encoding = "UTF-8"
      )

      if (length(header_value) == 1L) {
        header_line <- header_value
      }

      break
    }

    if (line_counter > 500000L) {
      break
    }
  }

  if (
    !table_begin_found ||
      is.na(header_line) ||
      !nzchar(header_line)
  ) {
    return(
      data.frame(
        read_status = "PLATFORM_TABLE_NOT_FOUND",
        n_columns = NA_integer_,
        first_columns = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }

  header_fields <- strsplit(
    header_line,
    "\t",
    fixed = TRUE
  )[[1]]

  header_fields <- gsub(
    '^"|"$',
    "",
    header_fields
  )

  data.frame(
    read_status = "READ_OK",
    n_columns = length(header_fields),
    first_columns = paste(
      header_fields,
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

acquire_platform_self <- function(
    platform_id,
    source_url,
    destination_root
) {
  destination_path <- file.path(
    destination_root,
    paste0(
      platform_id,
      "_self_full.soft"
    )
  )

  download_ok <- try_download_with_retry(
    source_url,
    destination_path
  )

  if (!download_ok) {
    stop(
      "Official GEO self-only platform download failed for ",
      platform_id,
      "."
    )
  }

  header_audit <- read_platform_soft_header(
    destination_path
  )

  if (
    header_audit$read_status !=
      "READ_OK"
  ) {
    stop(
      "Downloaded GEO self-only platform file failed table-header validation for ",
      platform_id,
      "."
    )
  }

  data.frame(
    platform_id =
      platform_id,
    source_scope =
      "SELF_ONLY_PLATFORM_RECORD",
    source_url =
      source_url,
    path =
      normalize_slashes(
        destination_path
      ),
    download_ok =
      TRUE,
    read_status =
      header_audit$read_status,
    n_columns =
      header_audit$n_columns,
    first_columns =
      header_audit$first_columns,
    stringsAsFactors = FALSE
  )
}

strip_prefix <- function(x, prefix) {
  sub(
    paste0(
      "^",
      prefix
    ),
    "",
    as.character(x)
  )
}

normalize_title <- function(x) {
  tolower(
    gsub(
      "[^a-z0-9]+",
      "",
      as.character(x),
      perl = TRUE
    )
  )
}

extract_array_suffix <- function(path) {
  sub(
    "^GSM[0-9]+_",
    "",
    basename(path),
    ignore.case = TRUE
  )
}

extract_slide_id <- function(path) {
  file_name <- basename(path)

  match_value <- regmatches(
    file_name,
    regexpr(
      "SG[0-9]+_[0-9]+",
      file_name,
      ignore.case = TRUE,
      perl = TRUE
    )
  )

  if (
    length(match_value) == 0L ||
      !nzchar(match_value)
  ) {
    return(NA_character_)
  }

  toupper(match_value)
}

extract_GSM <- function(path) {
  match_value <- regmatches(
    basename(path),
    regexpr(
      "GSM[0-9]+",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  )

  if (
    length(match_value) == 0L ||
      !nzchar(match_value)
  ) {
    return(NA_character_)
  }

  toupper(match_value)
}

digest_decompressed_gzip <- function(path) {
  connection <- gzfile(
    path,
    open = "rb"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  chunks <- list()
  chunk_index <- 1L

  repeat {
    raw_chunk <- readBin(
      connection,
      what = "raw",
      n = 4L * 1024L * 1024L
    )

    if (length(raw_chunk) == 0L) {
      break
    }

    chunks[[chunk_index]] <- raw_chunk
    chunk_index <- chunk_index + 1L
  }

  raw_value <- do.call(
    c,
    chunks
  )

  digest::digest(
    raw_value,
    algo = "sha256",
    serialize = FALSE
  )
}

cramers_v <- function(table_value) {
  chi_value <- suppressWarnings(
    chisq.test(
      table_value,
      correct = FALSE
    )$statistic
  )

  n_value <- sum(table_value)
  min_dimension <- min(
    nrow(table_value) - 1L,
    ncol(table_value) - 1L
  )

  if (
    n_value <= 0L ||
      min_dimension <= 0L
  ) {
    return(NA_real_)
  }

  sqrt(
    as.numeric(chi_value) /
      (
        n_value *
          min_dimension
      )
  )
}


log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Validate STEP04B accepted state
# ----------------------------------------------------------------------

step04B_gate_path <- file.path(
  dirs$results,
  "STEP04B_GATE_TABLE.csv"
)

step04B_locked_root <- file.path(
  dirs$results,
  "STEP04B_locked_combined_inputs",
  "COMBINED_INPUT_SOURCE_v1.0"
)

step04B_manifest_path <- file.path(
  step04B_locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv"
)

design_candidate_root <- file.path(
  dirs$results,
  "STEP04B_sample_design_candidates"
)

required_upstream <- c(
  step04B_gate_path,
  step04B_manifest_path,
  design_candidate_root
)

missing_upstream <- required_upstream[
  !file.exists(required_upstream) &
    !dir.exists(required_upstream)
]

if (length(missing_upstream) > 0L) {
  stop(
    "Required STEP04B input is missing:\n",
    paste(
      missing_upstream,
      collapse = "\n"
    )
  )
}

step04B_gate <- read.csv(
  step04B_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

step04B_manifest <- read.csv(
  step04B_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

accepted_gate <- step04B_gate$status[
  step04B_gate$gate_component ==
    "STEP04B gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "PROCEED_TO_MANUAL_SAMPLE_DESIGN_LOCK_AND_COMBINED_ANALYSIS_BUILD"
) {
  stop(
    "STEP04B does not authorize sample-design lock."
  )
}

if (
  nrow(step04B_manifest) != 5L ||
    !all(
      step04B_manifest$lock_status ==
        "LOCKED"
    )
) {
  stop(
    "STEP04B locked manifest is incomplete."
  )
}

step04B_manifest$file_exists <- file.exists(
  step04B_manifest$path
)

step04B_manifest$observed_md5 <- vapply(
  step04B_manifest$path,
  md5_file,
  character(1)
)

step04B_manifest$observed_sha256 <- vapply(
  step04B_manifest$path,
  sha256_file,
  character(1)
)

step04B_manifest$identity_pass <- with(
  step04B_manifest,
  file_exists &
    md5 ==
      observed_md5 &
    sha256 ==
      observed_sha256
)

if (!all(
  step04B_manifest$identity_pass
)) {
  stop(
    "A STEP04B locked artifact failed hash validation."
  )
}

read_design <- function(dataset_id) {
  path <- file.path(
    design_candidate_root,
    paste0(
      "STEP04B_",
      dataset_id,
      "_SAMPLE_SHEET_CANDIDATE.csv"
    )
  )

  if (!file.exists(path)) {
    stop(
      "Sample-design candidate is missing: ",
      path
    )
  }

  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

sheet_250130 <- read_design("GSE250130")
sheet_165004 <- read_design("GSE165004")
sheet_92324 <- read_design("GSE92324")
sheet_111974 <- read_design("GSE111974")
sheet_183837 <- read_design("GSE183837")

# ----------------------------------------------------------------------
# GSE250130 locked design
# ----------------------------------------------------------------------

design_250130 <- data.frame(
  dataset_id = "GSE250130",
  sample_accession =
    sheet_250130$sample_accession,
  sample_title =
    sheet_250130$title__1,
  raw_group = strip_prefix(
    sheet_250130$characteristics_ch1__3,
    "group: "
  ),
  donor_id =
    sheet_250130$sample_accession,
  stringsAsFactors = FALSE
)

design_250130$clinical_group <- ifelse(
  design_250130$raw_group == "RIF",
  "RIF",
  "FERTILE"
)

design_250130$LH_day <- ifelse(
  design_250130$raw_group == "LH3",
  3L,
  ifelse(
    design_250130$raw_group == "LH5",
    5L,
    ifelse(
      design_250130$raw_group %in%
        c(
          "LH7",
          "LH7_PGT",
          "RIF"
        ),
      7L,
      ifelse(
        design_250130$raw_group == "LH9",
        9L,
        ifelse(
          design_250130$raw_group == "LH11",
          11L,
          NA_integer_
        )
      )
    )
  )
)

design_250130$temporal_fertile_eligible <-
  design_250130$clinical_group ==
    "FERTILE"

design_250130$LH7_RIF_comparison_eligible <-
  design_250130$LH_day == 7L

design_250130$LH7_reference_group <- ifelse(
  design_250130$LH_day != 7L,
  NA_character_,
  ifelse(
    design_250130$clinical_group ==
      "RIF",
    "RIF",
    "FERTILE_LH7"
  )
)

design_250130$LH7_fertile_subgroup <- ifelse(
  design_250130$raw_group ==
    "LH7_PGT",
  "LH7_PGT",
  ifelse(
    design_250130$raw_group ==
      "LH7",
    "LH7_STANDARD",
    NA_character_
  )
)

design_250130$analysis_role <-
  "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT"

expected_250130_counts <- c(
  LH3 = 3L,
  LH5 = 3L,
  LH7 = 3L,
  LH7_PGT = 3L,
  LH9 = 3L,
  LH11 = 3L,
  RIF = 10L
)

observed_250130_counts <- table(
  design_250130$raw_group
)

if (
  !setequal(
    names(observed_250130_counts),
    names(expected_250130_counts)
  ) ||
    !all(
      as.integer(
        observed_250130_counts[
          names(expected_250130_counts)
        ]
      ) ==
        as.integer(
          expected_250130_counts
        )
    ) ||
    sum(
      design_250130$temporal_fertile_eligible
    ) !=
      18L ||
    sum(
      design_250130$LH7_reference_group ==
        "FERTILE_LH7",
      na.rm = TRUE
    ) !=
      6L ||
    sum(
      design_250130$LH7_reference_group ==
        "RIF",
      na.rm = TRUE
    ) !=
      10L
) {
  stop(
    "GSE250130 design failed 18-fertile/6-LH7/10-RIF reconciliation."
  )
}

# ----------------------------------------------------------------------
# GSE165004 locked design and slide IDs
# ----------------------------------------------------------------------

design_165004 <- data.frame(
  dataset_id = "GSE165004",
  sample_accession =
    sheet_165004$sample_accession,
  sample_title =
    sheet_165004$title__1,
  group = strip_prefix(
    sheet_165004$characteristics_ch1__1,
    "subject status/group: "
  ),
  stringsAsFactors = FALSE
)

design_165004$group[
  design_165004$group ==
    "patient with RPL"
] <- "RPL"

design_165004$group[
  design_165004$group ==
    "patient with UIF"
] <- "UIF"

design_165004$group[
  design_165004$group ==
    "Control"
] <- "CONTROL"

design_165004$analysis_role <-
  "PRIMARY_PATIENT_PROJECTION"

step00_raw_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP00_dataset_identity_and_technical_audit",
  "01_inputs",
  "raw_extracted"
)

raw_165004 <- list.files(
  file.path(
    step00_raw_root,
    "GSE165004"
  ),
  pattern = "\\.txt\\.gz$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (length(raw_165004) != 72L) {
  stop(
    "Expected 72 GSE165004 raw TXT.GZ files but found ",
    length(raw_165004),
    "."
  )
}

raw_165004_map <- data.frame(
  sample_accession = vapply(
    raw_165004,
    extract_GSM,
    character(1)
  ),
  raw_file =
    normalize_slashes(
      raw_165004
    ),
  array_suffix = vapply(
    raw_165004,
    extract_array_suffix,
    character(1)
  ),
  slide_id = vapply(
    raw_165004,
    extract_slide_id,
    character(1)
  ),
  stringsAsFactors = FALSE
)

design_165004 <- merge(
  design_165004,
  raw_165004_map,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

design_165004 <- design_165004[
  match(
    sheet_165004$sample_accession,
    design_165004$sample_accession
  ),
  ,
  drop = FALSE
]

if (
  any(is.na(
    design_165004$slide_id
  )) ||
    length(
      unique(
        design_165004$slide_id
      )
    ) !=
      9L
) {
  stop(
    "GSE165004 slide-ID derivation failed."
  )
}

slide_group_table <- table(
  design_165004$slide_id,
  design_165004$group
)

slide_group_nonzero_n <- rowSums(
  slide_group_table > 0L
)

slide_purity_pass <- all(
  slide_group_nonzero_n == 1L
)

slide_cramers_v <- cramers_v(
  slide_group_table
)

if (
  !slide_purity_pass ||
    abs(
      slide_cramers_v - 1
    ) > 1e-12 ||
    any(
      table(
        design_165004$group
      ) != 24L
    )
) {
  stop(
    "GSE165004 slide/group structure failed the locked phenotype-pure design."
  )
}

slide_group_counts <- as.data.frame(
  slide_group_table,
  stringsAsFactors = FALSE
)

names(slide_group_counts) <- c(
  "slide_id",
  "group",
  "sample_n"
)

# ----------------------------------------------------------------------
# GSE92324 locked design
# ----------------------------------------------------------------------

design_92324 <- data.frame(
  dataset_id = "GSE92324",
  sample_accession =
    sheet_92324$sample_accession,
  sample_title =
    sheet_92324$title__1,
  sample_identifier = strip_prefix(
    sheet_92324$characteristics_ch1__2,
    "sample group and identifier: "
  ),
  age = as.numeric(
    strip_prefix(
      sheet_92324$characteristics_ch1__3,
      "age: "
    )
  ),
  collection_day = strip_prefix(
    sheet_92324$characteristics_ch1__4,
    "sample collection day: "
  ),
  stringsAsFactors = FALSE
)

design_92324$group <- ifelse(
  grepl(
    "^control_",
    design_92324$sample_identifier
  ),
  "CONTROL",
  ifelse(
    grepl(
      "^case_",
      design_92324$sample_identifier
    ),
    "IMPLANTATION_FAILURE",
    NA_character_
  )
)

design_92324$analysis_role <-
  "INDEPENDENT_STIMULATED_CYCLE_VALIDATION"

if (
  sum(
    design_92324$group ==
      "CONTROL"
  ) !=
    8L ||
    sum(
      design_92324$group ==
        "IMPLANTATION_FAILURE"
    ) !=
      10L ||
    any(is.na(
      design_92324$group
    ))
) {
  stop(
    "GSE92324 group assignment failed 8-control/10-case reconciliation."
  )
}

design_matrix_92324 <- model.matrix(
  ~ group + collection_day,
  data = design_92324
)

design_rank_92324 <- qr(
  design_matrix_92324
)$rank

design_full_rank_92324 <-
  design_rank_92324 ==
  ncol(
    design_matrix_92324
  )

if (!design_full_rank_92324) {
  stop(
    "GSE92324 group plus collection-day design matrix is not full rank."
  )
}

group_day_92324 <- as.data.frame(
  table(
    design_92324$group,
    design_92324$collection_day
  ),
  stringsAsFactors = FALSE
)

names(group_day_92324) <- c(
  "group",
  "collection_day",
  "sample_n"
)

# ----------------------------------------------------------------------
# GSE111974 locked design and shared-control identity
# ----------------------------------------------------------------------

design_111974 <- data.frame(
  dataset_id = "GSE111974",
  sample_accession =
    sheet_111974$sample_accession,
  sample_title =
    sheet_111974$title__1,
  stringsAsFactors = FALSE
)

design_111974$group <- ifelse(
  grepl(
    "_Control_",
    design_111974$sample_title,
    fixed = TRUE
  ),
  "CONTROL",
  ifelse(
    grepl(
      "_RIF_",
      design_111974$sample_title,
      fixed = TRUE
    ),
    "RIF",
    NA_character_
  )
)

design_111974$analysis_role <-
  "SHARED_CONTROL_RIF_DISEASE_CONTEXT_EXTENSION"

if (
  sum(
    design_111974$group ==
      "CONTROL"
  ) !=
    24L ||
    sum(
      design_111974$group ==
        "RIF"
    ) !=
      24L ||
    any(is.na(
      design_111974$group
    ))
) {
  stop(
    "GSE111974 group assignment failed 24-control/24-RIF reconciliation."
  )
}

raw_111974 <- list.files(
  file.path(
    step00_raw_root,
    "GSE111974"
  ),
  pattern = "\\.txt\\.gz$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)

if (length(raw_111974) != 48L) {
  stop(
    "Expected 48 GSE111974 raw TXT.GZ files but found ",
    length(raw_111974),
    "."
  )
}

raw_111974_map <- data.frame(
  sample_accession = vapply(
    raw_111974,
    extract_GSM,
    character(1)
  ),
  raw_file =
    normalize_slashes(
      raw_111974
    ),
  array_suffix = vapply(
    raw_111974,
    extract_array_suffix,
    character(1)
  ),
  slide_id = vapply(
    raw_111974,
    extract_slide_id,
    character(1)
  ),
  stringsAsFactors = FALSE
)

design_111974 <- merge(
  design_111974,
  raw_111974_map,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

design_111974 <- design_111974[
  match(
    sheet_111974$sample_accession,
    design_111974$sample_accession
  ),
  ,
  drop = FALSE
]

control_165004 <- design_165004[
  design_165004$group ==
    "CONTROL",
  c(
    "sample_accession",
    "sample_title",
    "raw_file",
    "array_suffix",
    "slide_id"
  ),
  drop = FALSE
]

control_111974 <- design_111974[
  design_111974$group ==
    "CONTROL",
  c(
    "sample_accession",
    "sample_title",
    "raw_file",
    "array_suffix",
    "slide_id"
  ),
  drop = FALSE
]

control_165004$normalized_title <-
  normalize_title(
    control_165004$sample_title
  )

control_111974$normalized_title <-
  normalize_title(
    control_111974$sample_title
  )

shared_control_pairs <- merge(
  control_165004,
  control_111974,
  by = "normalized_title",
  suffixes = c(
    "_GSE165004",
    "_GSE111974"
  ),
  all = FALSE,
  sort = TRUE
)

if (nrow(shared_control_pairs) != 24L) {
  stop(
    "Expected 24 normalized-title shared control pairs but found ",
    nrow(shared_control_pairs),
    "."
  )
}

shared_control_pairs$array_suffix_match <-
  shared_control_pairs$array_suffix_GSE165004 ==
  shared_control_pairs$array_suffix_GSE111974

shared_control_pairs$slide_id_match <-
  shared_control_pairs$slide_id_GSE165004 ==
  shared_control_pairs$slide_id_GSE111974

log_event(
  "Computing decompressed hashes for 24 shared-control array pairs."
)

shared_control_pairs$decompressed_sha256_GSE165004 <-
  vapply(
    shared_control_pairs$raw_file_GSE165004,
    digest_decompressed_gzip,
    character(1)
  )

shared_control_pairs$decompressed_sha256_GSE111974 <-
  vapply(
    shared_control_pairs$raw_file_GSE111974,
    digest_decompressed_gzip,
    character(1)
  )

shared_control_pairs$decompressed_content_match <-
  shared_control_pairs$decompressed_sha256_GSE165004 ==
  shared_control_pairs$decompressed_sha256_GSE111974

shared_control_pairs$shared_biological_control_lock <-
  shared_control_pairs$array_suffix_match &
  shared_control_pairs$slide_id_match &
  shared_control_pairs$decompressed_content_match

if (!all(
  shared_control_pairs$shared_biological_control_lock
)) {
  stop(
    "At least one GSE165004/GSE111974 control pair failed exact array-content identity."
  )
}

design_111974$shared_control_with_GSE165004 <-
  FALSE

design_111974$paired_GSE165004_sample <-
  NA_character_

for (row_index in seq_len(
  nrow(shared_control_pairs)
)) {
  target_index <-
    design_111974$sample_accession ==
    shared_control_pairs$sample_accession_GSE111974[
      row_index
    ]

  design_111974$shared_control_with_GSE165004[
    target_index
  ] <- TRUE

  design_111974$paired_GSE165004_sample[
    target_index
  ] <-
    shared_control_pairs$sample_accession_GSE165004[
      row_index
    ]
}

design_111974$independent_external_validation <-
  FALSE

design_111974$validation_scope <- ifelse(
  design_111974$group ==
    "RIF",
  "RIF_DISEASE_CONTEXT_EXTENSION_WITH_SHARED_CONTROLS",
  "SHARED_CONTROL_REUSE_NOT_INDEPENDENT"
)

# ----------------------------------------------------------------------
# GSE183837 locked design
# ----------------------------------------------------------------------

design_183837 <- data.frame(
  dataset_id = "GSE183837",
  sample_accession =
    sheet_183837$sample_accession,
  sample_title =
    sheet_183837$title__1,
  group = toupper(
    strip_prefix(
      sheet_183837$characteristics_ch1__3,
      "group: "
    )
  ),
  stage = strip_prefix(
    sheet_183837$characteristics_ch1__2,
    "Stage: "
  ),
  donor_id =
    sheet_183837$sample_accession,
  stringsAsFactors = FALSE
)

design_183837$group[
  design_183837$group ==
    "CONTROL"
] <- "CONTROL"

design_183837$analysis_role <-
  "INDEPENDENT_SINGLE_CELL_LOCALIZATION"

if (
  sum(
    design_183837$group ==
      "CONTROL"
  ) !=
    3L ||
    sum(
      design_183837$group ==
        "RIF"
    ) !=
      6L ||
    length(
      unique(
        design_183837$stage
      )
    ) !=
      1L
) {
  stop(
    "GSE183837 design failed 3-control/6-RIF stage-matched reconciliation."
  )
}

# ----------------------------------------------------------------------
# Audit stale family-SOFT downloads from the failed v2 attempt.
stale_family_paths <- file.path(
  platform_root,
  c(
    "GPL16699_family.soft.gz",
    "GPL17077_family.soft.gz"
  )
)

stale_family_audit <- data.frame(
  path =
    normalize_slashes(
      stale_family_paths
    ),
  exists =
    file.exists(
      stale_family_paths
    ),
  size_bytes = as.numeric(
    file.info(
      stale_family_paths
    )$size
  ),
  analysis_use =
    "NOT_USED_DELETE_AFTER_MANUAL_CONFIRMATION",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  stale_family_audit,
  file.path(
    dirs$qc,
    "STEP04C_v3_STALE_FAMILY_SOFT_AUDIT.csv"
  )
)

# Acquire bulk platform annotations
# ----------------------------------------------------------------------
#
# Important: use GEO `targ=self` only. Do not download the platform
# family SOFT files, because those may bundle large numbers of linked
# Sample and Series records and can be gigabytes in size.

GPL16699_self <- acquire_platform_self(
  platform_id = "GPL16699",
  source_url =
    GPL16699_self_url,
  destination_root =
    platform_root
)

GPL17077_self <- acquire_platform_self(
  platform_id = "GPL17077",
  source_url =
    GPL17077_self_url,
  destination_root =
    platform_root
)

GPL16699_path <-
  GPL16699_self$path[[1]]

GPL17077_path <-
  GPL17077_self$path[[1]]

GPL16699_url <-
  GPL16699_self$source_url[[1]]

GPL17077_url <-
  GPL17077_self$source_url[[1]]

platform_download_attempts <- rbind(
  GPL16699_self,
  GPL17077_self
)

write_csv_utf8(
  platform_download_attempts,
  file.path(
    dirs$qc,
    "STEP04C_v3_PLATFORM_SELF_DOWNLOAD_AUDIT.csv"
  )
)

GPL10558_path <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE92324",
  "platform_annotation",
  "GPL10558.soft.gz"
)

if (!file.exists(GPL10558_path)) {
  stop(
    "Locked GSE92324 GPL10558 annotation is missing."
  )
}

platform_rows <- list(
  GPL16699 = cbind(
    data.frame(
      platform_id = "GPL16699",
      source_url = GPL16699_url,
      path = normalize_slashes(
        GPL16699_path
      ),
      stringsAsFactors = FALSE
    ),
    read_platform_soft_header(
      GPL16699_path
    )
  ),
  GPL17077 = cbind(
    data.frame(
      platform_id = "GPL17077",
      source_url = GPL17077_url,
      path = normalize_slashes(
        GPL17077_path
      ),
      stringsAsFactors = FALSE
    ),
    read_platform_soft_header(
      GPL17077_path
    )
  ),
  GPL10558 = cbind(
    data.frame(
      platform_id = "GPL10558",
      source_url = "LOCKED_LOCAL_GEO_SOFT_SOURCE",
      path = normalize_slashes(
        GPL10558_path
      ),
      stringsAsFactors = FALSE
    ),
    read_platform_soft_header(
      GPL10558_path
    )
  )
)

platform_audit <- do.call(
  rbind,
  platform_rows
)

platform_audit$size_bytes <- as.numeric(
  file.info(
    platform_audit$path
  )$size
)

platform_audit$md5 <- vapply(
  platform_audit$path,
  md5_file,
  character(1)
)

platform_audit$sha256 <- vapply(
  platform_audit$path,
  sha256_file,
  character(1)
)

platform_audit$source_version <-
  platform_source_version

platform_audit$lock_status <-
  "LOCKED_SOURCE_REFERENCE"

if (
  any(
    platform_audit$read_status !=
      "READ_OK"
  ) ||
    any(
      platform_audit$size_bytes <
        10000
    )
) {
  stop(
    "At least one platform annotation failed acquisition or header audit."
  )
}

# ----------------------------------------------------------------------
# Write locked designs
# ----------------------------------------------------------------------

design_paths <- c(
  GSE250130 = file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"
  ),
  GSE165004 = file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE165004.csv"
  ),
  GSE92324 = file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE92324.csv"
  ),
  GSE111974 = file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE111974.csv"
  ),
  GSE183837 = file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"
  )
)

write_csv_utf8(
  design_250130,
  design_paths[["GSE250130"]]
)

write_csv_utf8(
  design_165004,
  design_paths[["GSE165004"]]
)

write_csv_utf8(
  design_92324,
  design_paths[["GSE92324"]]
)

write_csv_utf8(
  design_111974,
  design_paths[["GSE111974"]]
)

write_csv_utf8(
  design_183837,
  design_paths[["GSE183837"]]
)

write_csv_utf8(
  slide_group_counts,
  file.path(
    dirs$results,
    "STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv"
  )
)

write_csv_utf8(
  group_day_92324,
  file.path(
    dirs$results,
    "STEP04C_GSE92324_GROUP_DAY_COUNTS.csv"
  )
)

write_csv_utf8(
  shared_control_pairs,
  file.path(
    dirs$results,
    "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"
  )
)

write_csv_utf8(
  platform_audit,
  file.path(
    dirs$results,
    "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"
  )
)

sample_design_summary <- data.frame(
  dataset_id = c(
    "GSE250130",
    "GSE165004",
    "GSE92324",
    "GSE111974",
    "GSE183837"
  ),
  total_sample_n = c(
    nrow(design_250130),
    nrow(design_165004),
    nrow(design_92324),
    nrow(design_111974),
    nrow(design_183837)
  ),
  locked_design = c(
    "18 fertile temporal; 6 fertile LH7 vs 10 RIF",
    "24 control; 24 RPL; 24 UIF; 9 phenotype-pure slides",
    "8 control; 10 implantation failure; HCG day retained",
    "24 shared controls; 24 RIF; not independent external validation",
    "3 control; 6 RIF; same post-ovulation stage"
  ),
  analysis_role = c(
    "TEMPORAL_AND_LH7_RIF_SINGLE_CELL_SUPPORT",
    "PRIMARY_PATIENT_PROJECTION",
    "INDEPENDENT_STIMULATED_CYCLE_VALIDATION",
    "SHARED_CONTROL_RIF_CONTEXT_EXTENSION",
    "INDEPENDENT_SINGLE_CELL_LOCALIZATION"
  ),
  lock_version =
    sample_design_lock_version,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  sample_design_summary,
  file.path(
    dirs$results,
    "STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv"
  )
)

locked_artifact_paths <- c(
  unname(design_paths),
  file.path(
    dirs$results,
    "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"
  )
)

locked_manifest <- data.frame(
  artifact_id = c(
    "GSE250130_SAMPLE_DESIGN",
    "GSE165004_SAMPLE_DESIGN",
    "GSE92324_SAMPLE_DESIGN",
    "GSE111974_SAMPLE_DESIGN",
    "GSE183837_SAMPLE_DESIGN",
    "SHARED_CONTROL_IDENTITY_LOCK",
    "BULK_PLATFORM_ANNOTATION_SOURCE_LOCK"
  ),
  lock_version = c(
    rep(
      sample_design_lock_version,
      6L
    ),
    platform_source_version
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
    "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv"
  )
)

# ----------------------------------------------------------------------
# Gate and authorization
# ----------------------------------------------------------------------

gate_status <-
  "READY_FOR_COMBINED_BULK_AND_SINGLE_CELL_ANALYSIS_BUILD"

gate_table <- data.frame(
  gate_component = c(
    "STEP04B gate",
    "GSE250130 temporal fertile design",
    "GSE250130 LH7 RIF comparison",
    "GSE165004 groups",
    "GSE165004 slide derivation",
    "GSE165004 slide/group confounding",
    "GSE92324 groups",
    "GSE92324 group plus HCG-day model rank",
    "GSE111974 groups",
    "GSE111974/GSE165004 shared controls",
    "GSE111974 external-independence status",
    "GSE183837 groups and stage",
    "Bulk platform annotations",
    "Sample-design lock",
    "Biological analysis",
    "Scientific figures",
    "STEP04C gate"
  ),
  status = c(
    "PASS",
    "LOCKED_18_FERTILE_ACROSS_5_TIMEPOINTS",
    "LOCKED_6_FERTILE_LH7_VS_10_RIF",
    "LOCKED_24_CONTROL_24_RPL_24_UIF",
    "PASS_72_SAMPLES_9_SLIDES",
    paste0(
      "PHENOTYPE_PURE_CRAMERS_V_",
      format(
        slide_cramers_v,
        digits = 6
      )
    ),
    "LOCKED_8_CONTROL_10_IF",
    paste0(
      "PASS_FULL_RANK_",
      design_rank_92324,
      "_OF_",
      ncol(
        design_matrix_92324
      )
    ),
    "LOCKED_24_CONTROL_24_RIF",
    "LOCKED_24_IDENTICAL_CONTROL_ARRAYS",
    "NOT_INDEPENDENT_SHARED_CONTROL_EXTENSION",
    "LOCKED_3_CONTROL_6_RIF_STAGE_MATCHED",
    "LOCKED_GPL16699_GPL17077_GPL10558",
    paste0(
      "LOCKED_",
      sample_design_lock_version
    ),
    "NOT_STARTED",
    "NONE_CREATED",
    gate_status
  ),
  required_action = c(
    "Retain the STEP04B source lock.",
    "Use donor-aware temporal single-cell summaries.",
    "Compare RIF only with fertile LH+7 donors.",
    "Use separate Control-vs-RPL and Control-vs-UIF projections.",
    "Use the derived slide IDs for sensitivity analyses.",
    "Do not include slide as a standard covariate because it is phenotype-pure; use leave-one-slide sensitivity.",
    "Retain HCG+6/+7 metadata.",
    "Use group plus collection day when gene-level modeling is required.",
    "Use as RIF disease-context support.",
    "Do not double-count controls across datasets.",
    "Do not call GSE111974 an independent external validation.",
    "Use donor-aware single-cell summaries.",
    "Build gene-level matrices only from locked platform sources.",
    "Do not change assignments without a new version.",
    "Proceed as one combined analysis package.",
    "Render only after biological QC.",
    "Build bulk and single-cell analysis scripts."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP04C_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP04C sample-design and platform lock card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Sample-design lock: ",
      sample_design_lock_version,
      " LOCKED"
    ),
    "- GSE250130: 18 fertile temporal samples; 6 fertile LH+7 versus 10 RIF",
    "- GSE165004: 24 control, 24 RPL, 24 UIF; 9 phenotype-pure slides",
    "- GSE92324: 8 control, 10 implantation failure; HCG+6/+7 retained",
    "- GSE111974: 24 RIF plus 24 controls identical to GSE165004 controls",
    "- GSE111974 independent external validation: NO",
    "- GSE183837: 3 control, 6 RIF; stage matched",
    "- Shared-control decompressed content matches: 24/24",
    "- Platform sources: GPL16699, GPL17077 and GPL10558 LOCKED",
    "- Biological analysis: NOT STARTED",
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "STEP04C_SAMPLE_DESIGN_PLATFORM_LOCK_CARD.md"
  )
)

write_text(
  c(
    "# Combined biological analysis authorization",
    "",
    paste0(
      "- Authorized by: ",
      run_id
    ),
    paste0(
      "- Sample-design source: ",
      sample_design_lock_version,
      " LOCKED"
    ),
    "- Authorization: YES",
    "",
    "## Bulk patient analyses",
    "",
    "- GSE165004: primary projection, Control vs UIF and Control vs RPL separately, with leave-one-slide sensitivity.",
    "- GSE92324: independent stimulated-cycle validation, retaining HCG collection day.",
    "- GSE111974: RIF disease-context extension sharing all 24 control arrays with GSE165004; never count it as an independent cohort.",
    "",
    "## Single-cell analyses",
    "",
    "- GSE250130: fertile time series across LH+3/+5/+7/+9/+11 and RIF-vs-fertile comparison restricted to LH+7.",
    "- GSE183837: independent Ctrl/RIF localization at five days after ovulation.",
    "- Use donor-aware cell-type pseudobulk or sample-level scores; cells are not independent replicates.",
    "",
    "## Scope limits",
    "",
    "- Use the locked 284 occupancy-supported RNA candidates as the primary projection set.",
    "- Retain the 152 RNA-only candidates for secondary sensitivity analysis.",
    "- Do not use discovery ROC, LASSO, random forest or generic network pharmacology.",
    "- Do not define definitive direct targets.",
    "- Do not generate scientific figures until analysis QC passes."
  ),
  file.path(
    dirs$control,
    "COMBINED_BULK_SINGLE_CELL_ANALYSIS_AUTHORIZATION.md"
  )
)

write_text(
  c(
    "# STEP04C Methods note",
    "",
    "Sample assignments were locked from GEO metadata and raw array identifiers. For GSE250130, fertile samples at LH+3, LH+5, LH+7, LH+9 and LH+11 were retained for temporal analysis, whereas the RIF comparison was restricted to the six fertile and ten RIF samples collected at LH+7.",
    "",
    "For GSE165004, raw Agilent filenames were used to reconstruct nine slide identifiers. Slides were phenotype-pure, so slide was treated through leave-one-slide sensitivity rather than as a conventional adjustment covariate. The 24 control arrays in GSE111974 were confirmed to be the same arrays used in GSE165004 by normalized title, array identifier and decompressed raw-file content.",
    "",
    "GSE92324 retained collection day as a design variable, and GSE183837 retained the donor-level 3-control/6-RIF design. Platform annotations were locked before gene-level expression mapping."
  ),
  file.path(
    dirs$text,
    "STEP04C_Methods_note.md"
  )
)

parameters <- list(
  schema_version = "1.0-lock",
  run_id = run_id,
  sample_design_lock_version =
    sample_design_lock_version,
  platform_source_version =
    platform_source_version,
  GSE250130_temporal_fertile_n = 18,
  GSE250130_fertile_LH7_n = 6,
  GSE250130_RIF_LH7_n = 10,
  GSE165004_control_n = 24,
  GSE165004_RPL_n = 24,
  GSE165004_UIF_n = 24,
  GSE165004_slide_n = 9,
  GSE165004_slide_cramers_v =
    slide_cramers_v,
  GSE92324_control_n = 8,
  GSE92324_IF_n = 10,
  GSE92324_model_full_rank =
    design_full_rank_92324,
  GSE111974_control_n = 24,
  GSE111974_RIF_n = 24,
  shared_control_identity_n =
    sum(
      shared_control_pairs$shared_biological_control_lock
    ),
  GSE111974_independent_external_validation =
    FALSE,
  GSE183837_control_n = 3,
  GSE183837_RIF_n = 6,
  biological_analysis_started =
    FALSE,
  scientific_figures_created =
    FALSE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP04C_parameters.json"
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
    "sessionInfo_POST_RUN_STEP04C.txt"
  )
)

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  copy_if_different(
    source_file,
    file.path(
      dirs$code,
      "50_step04C_v3_lock_sample_designs_and_platform_preflight.R"
    )
  )
}

# ----------------------------------------------------------------------
# Short archive: no raw arrays or large single-cell matrices
# ----------------------------------------------------------------------

archive_artifacts <- c(
  file.path(
    dirs$results,
    "STEP04C_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_SAMPLE_DESIGN_PLATFORM_LOCK_CARD.md"
  ),
  file.path(
    dirs$results,
    "STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_GSE92324_GROUP_DAY_COUNTS.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"
  ),
  file.path(
    dirs$results,
    "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04C_v3_PLATFORM_SELF_DOWNLOAD_AUDIT.csv"
  ),
  file.path(
    dirs$qc,
    "STEP04C_v3_STALE_FAMILY_SOFT_AUDIT.csv"
  ),
  file.path(
    locked_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv"
  ),
  unname(design_paths),
  file.path(
    dirs$control,
    "COMBINED_BULK_SINGLE_CELL_ANALYSIS_AUTHORIZATION.md"
  ),
  file.path(
    dirs$text,
    "STEP04C_Methods_note.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04C_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04C.txt"
  ),
  file.path(
    dirs$code,
    "50_step04C_v3_lock_sample_designs_and_platform_preflight.R"
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
      "Methods|AUTHORIZATION",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "LOCK|CARD|COUNTS|MANIFEST|SUMMARY|AUDIT|GATE",
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
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP04C_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04C archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Sample-design lock: ",
      sample_design_lock_version
    ),
    "- Shared control identity: 24/24",
    "- Platform annotations: 3 locked sources",
    "- Raw data copied into archive: NO",
    "- Biological analysis: NOT STARTED",
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
      "- Accepted run: ",
      run_id
    ),
    paste0(
      "- Sample designs: ",
      sample_design_lock_version,
      " LOCKED"
    ),
    paste0(
      "- Platform sources: ",
      platform_source_version,
      " LOCKED"
    ),
    "- GSE111974 role: SHARED_CONTROL_RIF_CONTEXT_EXTENSION",
    "- Biological analysis: AUTHORIZED_NOT_STARTED",
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
cat("STEP04C v3 sample-design and platform preflight completed.\n")
cat(
  "Sample-design lock: ",
  sample_design_lock_version,
  " LOCKED\n"
)
cat("GSE250130 temporal fertile samples: 18\n")
cat("GSE250130 LH7 comparison: 6 fertile vs 10 RIF\n")
cat("GSE165004: 24 control / 24 RPL / 24 UIF\n")
cat(
  "GSE165004 phenotype-pure slides: 9; Cramer's V = ",
  format(
    slide_cramers_v,
    digits = 6
  ),
  "\n"
)
cat("GSE92324: 8 control / 10 implantation failure\n")
cat("GSE92324 group + HCG-day model: FULL RANK\n")
cat("GSE111974: 24 shared controls / 24 RIF\n")
cat("Shared control decompressed identity: 24/24\n")
cat("GSE111974 independent external validation: NO\n")
cat("GSE183837: 3 control / 6 RIF; stage matched\n")
cat("Platform annotations locked: GPL16699 / GPL17077 / GPL10558\n")
cat("Biological analysis: AUTHORIZED_NOT_STARTED\n")
cat("Scientific figures: NONE\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n\n")
cat("Upload:\n")
cat("  06_results/STEP04C_GATE_TABLE.csv\n")
cat("  06_results/STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv\n")
cat("  06_results/STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv\n")
cat("  06_results/STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv\n")
cat("  06_results/STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv\n")
cat("  06_results/STEP04C_GSE92324_GROUP_DAY_COUNTS.csv\n")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 4: 51_step04D_v6_bulk_projection_and_sc_preQC
# ==============================================================================
# ============================================================
# STEP04D - Build bulk gene-level matrices, run patient
# projection analyses, and prepare donor-level single-cell QC
# objects for the combined downstream package
# Nutrition_Endometrial_VDR_v2
#
# Inputs:
# - RNA_OCCUPANCY_EVIDENCE_v1.0 (436 candidates)
# - COMBINED_INPUT_SOURCE_v1.0
# - COMBINED_SAMPLE_DESIGN_v1.0
# - BULK_PLATFORM_ANNOTATION_SOURCE_v1.0
#
# Bulk analyses:
# - GSE165004: Control vs RPL and Control vs UIF
# - GSE92324: implantation failure vs Control, adjusted for HCG day
# - GSE111974: RIF vs shared Control; context extension only
# - four prespecified VDR-program scores
# - GSE165004 leave-one-slide sensitivity for genes and scores
#
# Single-cell preparation:
# - GSE250130: 28 donor/sample matrices
# - GSE183837: 9 donor/sample matrices
# - per-sample sparse-matrix validation and conservative QC
# - filtered donor-level RDS objects for later cell-type integration
#
# Interpretation safeguards:
# - no diagnostic classifier or discovery ROC
# - no definitive direct-target assignment
# - GSE111974 is not counted as independent external validation
# - single cells are not treated as independent biological replicates
# - no scientific figures are generated in this step
#
# Run from project root:
# source(
#   "04_code/R/51_step04D_v6_bulk_projection_and_sc_preQC.R",
#   encoding = "UTF-8"
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)
options(timeout = max(3600L, getOption("timeout")))

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04D_v6_",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

bulk_matrix_version <- "BULK_GENE_MATRIX_v1.2"
bulk_projection_version <- "BULK_VDR_PROGRAM_PROJECTION_v1.2"
sc_preQC_version <- "SINGLE_CELL_PREQC_OBJECTS_v1.0"

GSE92324_profile_url <-
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92324/suppl/GSE92324_sample_probe_profile.txt.gz"

GSE92324_profile_path <- file.path(
  project_root,
  "03_data",
  "raw_external_READ_ONLY",
  "GSE92324",
  "supplementary",
  "GSE92324_sample_probe_profile.txt.gz"
)

required_packages <- c(
  "data.table",
  "digest",
  "jsonlite",
  "limma",
  "Matrix"
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
  control = file.path(step04_root, "00_module_control"),
  inputs = file.path(step04_root, "01_inputs"),
  code = file.path(step04_root, "02_code"),
  parameters = file.path(step04_root, "03_parameters"),
  environment = file.path(step04_root, "04_environment"),
  execution = file.path(step04_root, "05_execution"),
  results = file.path(step04_root, "06_results"),
  figures = file.path(step04_root, "07_figures_tables"),
  qc = file.path(step04_root, "08_qc"),
  text = file.path(step04_root, "09_manuscript_text"),
  archive_parent = file.path(step04_root, "10_archive_work")
)

for (directory_path in dirs) {
  dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

bulk_matrix_root <- file.path(
  dirs$results,
  "S4D_bulk_matrix",
  "BGM_v1_2"
)

bulk_result_root <- file.path(
  dirs$results,
  "S4D_bulk",
  "BVP_v1_2"
)

sc_object_root <- file.path(
  dirs$results,
  "S4D_sc_preQC",
  "SCP_v1_0"
)

previous_bulk_matrix_root <- file.path(
  dirs$results,
  "STEP04D_bulk_gene_matrices",
  "BULK_GENE_MATRIX_v1.1"
)

for (directory_path in c(
  bulk_matrix_root,
  bulk_result_root,
  sc_object_root
)) {
  directory_created <- dir.create(
    directory_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (
    !dir.exists(directory_path)
  ) {
    stop(
      "Required STEP04D output directory could not be created: ",
      directory_path
    )
  }

  if (nchar(directory_path) >= 220L) {
    stop(
      "STEP04D output directory exceeds the 220-character soft limit: ",
      directory_path
    )
  }
}

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4D_",
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
      "Output directory does not exist and could not be created: ",
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
      normalize_slashes(path_a)
    }
  )

  normalized_b <- tryCatch(
    normalizePath(
      path_b,
      winslash = "/",
      mustWork = FALSE
    ),
    error = function(error_object) {
      normalize_slashes(path_b)
    }
  )

  identical(
    tolower(normalized_a),
    tolower(normalized_b)
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

download_with_resume <- function(
    url,
    destination
) {
  dir.create(
    dirname(destination),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (
    file.exists(destination) &&
      as.numeric(
        file.info(destination)$size
      ) > 10000
  ) {
    return(TRUE)
  }

  curl_candidates <- c(
    Sys.which("curl.exe"),
    Sys.which("curl")
  )

  curl_candidates <- curl_candidates[
    nzchar(curl_candidates)
  ]

  if (length(curl_candidates) > 0L) {
    status <- system2(
      curl_candidates[[1]],
      args = c(
        "--fail",
        "--location",
        "--retry",
        "15",
        "--retry-all-errors",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "60",
        "--continue-at",
        "-",
        "--output",
        shQuote(destination),
        shQuote(url)
      )
    )

    if (
      status != 0L ||
        !file.exists(destination) ||
        as.numeric(
          file.info(destination)$size
        ) <= 10000
    ) {
      stop(
        "Download failed: ",
        url
      )
    }

    return(TRUE)
  }

  download.file(
    url = url,
    destfile = destination,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )

  if (
    !file.exists(destination) ||
      as.numeric(
        file.info(destination)$size
      ) <= 10000
  ) {
    stop(
      "Download failed: ",
      url
    )
  }

  TRUE
}

read_tabular_maybe_gz <- function(path) {
  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  read.delim(
    connection,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

extract_soft_table <- function(
    path,
    begin_marker,
    end_marker
) {
  if (!file.exists(path)) {
    stop(
      "SOFT/matrix source is missing: ",
      path
    )
  }

  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  temp_path <- tempfile(
    fileext = ".tsv"
  )

  output_connection <- file(
    temp_path,
    open = "wt"
  )

  output_connection_open <- TRUE

  on.exit(
    {
      if (isTRUE(output_connection_open)) {
        try(
          close(output_connection),
          silent = TRUE
        )
      }

      if (file.exists(temp_path)) {
        unlink(
          temp_path,
          force = TRUE
        )
      }
    },
    add = TRUE
  )

  in_table <- FALSE
  found_begin <- FALSE
  found_end <- FALSE

  repeat {
    lines <- readLines(
      connection,
      n = 5000L,
      warn = FALSE,
      encoding = "UTF-8"
    )

    if (length(lines) == 0L) {
      break
    }

    if (!in_table) {
      begin_index <- which(
        trimws(lines) ==
          begin_marker
      )

      if (length(begin_index) > 0L) {
        found_begin <- TRUE
        in_table <- TRUE

        start_index <-
          begin_index[[1]] + 1L

        if (start_index <= length(lines)) {
          lines <- lines[
            start_index:length(lines)
          ]
        } else {
          lines <- character(0)
        }
      } else {
        next
      }
    }

    if (
      in_table &&
        length(lines) > 0L
    ) {
      end_index <- which(
        trimws(lines) ==
          end_marker
      )

      if (length(end_index) > 0L) {
        end_position <-
          end_index[[1]] - 1L

        if (end_position >= 1L) {
          writeLines(
            lines[
              seq_len(
                end_position
              )
            ],
            output_connection,
            useBytes = TRUE
          )
        }

        found_end <- TRUE
        break
      }

      writeLines(
        lines,
        output_connection,
        useBytes = TRUE
      )
    }
  }

  close(output_connection)
  output_connection_open <- FALSE

  if (
    !found_begin ||
      !found_end
  ) {
    stop(
      "Could not locate complete table markers in: ",
      path
    )
  }

  table_value <- data.table::fread(
    temp_path,
    data.table = FALSE,
    showProgress = FALSE,
    check.names = FALSE
  )

  table_value
}

read_series_expression <- function(path) {
  table_value <- extract_soft_table(
    path = path,
    begin_marker =
      "!series_matrix_table_begin",
    end_marker =
      "!series_matrix_table_end"
  )

  if (ncol(table_value) < 3L) {
    stop(
      "Series matrix has fewer than three columns: ",
      path
    )
  }

  probe_id <- trimws(
    gsub(
      '^"|"$',
      "",
      as.character(
        table_value[[1]]
      )
    )
  )

  expression_value <- as.matrix(
    table_value[
      ,
      -1,
      drop = FALSE
    ]
  )

  storage.mode(
    expression_value
  ) <- "numeric"

  rownames(
    expression_value
  ) <- probe_id

  sample_names <- gsub(
    '^"|"$',
    "",
    colnames(
      expression_value
    )
  )

  colnames(
    expression_value
  ) <- sample_names

  list(
    expression =
      expression_value,
    probe_id =
      probe_id
  )
}

read_platform_table <- function(path) {
  extract_soft_table(
    path = path,
    begin_marker =
      "!platform_table_begin",
    end_marker =
      "!platform_table_end"
  )
}

choose_column <- function(
    table_value,
    candidate_names
) {
  normalized_names <- toupper(
    gsub(
      "[^A-Z0-9]+",
      "_",
      names(table_value)
    )
  )

  normalized_candidates <- toupper(
    gsub(
      "[^A-Z0-9]+",
      "_",
      candidate_names
    )
  )

  match_index <- match(
    normalized_candidates,
    normalized_names
  )

  match_index <- match_index[
    !is.na(match_index)
  ]

  if (length(match_index) == 0L) {
    return(NA_character_)
  }

  names(table_value)[
    match_index[[1]]
  ]
}

clean_single_gene_symbol <- function(x) {
  values <- as.character(x)

  vapply(
    values,
    function(value) {
      if (
        is.na(value) ||
          !nzchar(
            trimws(value)
          ) ||
          trimws(value) %in%
            c(
              "---",
              "NA",
              "N/A",
              "."
            )
      ) {
        return(NA_character_)
      }

      tokens <- unlist(
        strsplit(
          value,
          "\\s*(///|//|;|\\||,)\\s*",
          perl = TRUE
        )
      )

      tokens <- unique(
        trimws(tokens)
      )

      tokens <- tokens[
        nzchar(tokens) &
          !tokens %in%
            c(
              "---",
              "NA",
              "N/A",
              "."
            ) &
          grepl(
            "^[A-Za-z0-9][A-Za-z0-9._-]*$",
            tokens
          )
      ]

      if (length(tokens) != 1L) {
        return(NA_character_)
      }

      tokens[[1]]
    },
    character(1)
  )
}

evaluate_annotation_table <- function(
    annotation_table,
    expression_probe_ids,
    symbol_candidates,
    source_label
) {
  if (
    nrow(annotation_table) == 0L ||
      ncol(annotation_table) < 2L
  ) {
    stop(
      "Annotation table is empty or has fewer than two columns: ",
      source_label
    )
  }

  expression_probe_ids <- trimws(
    as.character(
      expression_probe_ids
    )
  )

  normalized_names <- toupper(
    gsub(
      "[^A-Z0-9]+",
      "_",
      names(annotation_table)
    )
  )

  id_preferred <- normalized_names %in%
    c(
      "ID",
      "ID_REF",
      "PROBE_ID",
      "PROBEID",
      "ARRAY_ADDRESS_ID",
      "ARRAYADDRESSID"
    )

  id_overlap_n <- vapply(
    seq_len(
      ncol(annotation_table)
    ),
    function(column_index) {
      values <- trimws(
        gsub(
          '^"|"$',
          "",
          as.character(
            annotation_table[[
              column_index
            ]]
          )
        )
      )

      length(
        intersect(
          unique(values),
          expression_probe_ids
        )
      )
    },
    integer(1)
  )

  if (any(id_preferred)) {
    preferred_indices <- which(
      id_preferred
    )

    preferred_best <- preferred_indices[
      which.max(
        id_overlap_n[
          preferred_indices
        ]
      )
    ]

    if (
      id_overlap_n[
        preferred_best
      ] >=
        max(
          1000L,
          floor(
            0.02 *
              length(
                expression_probe_ids
              )
          )
        )
    ) {
      id_column_index <-
        preferred_best
    } else {
      id_column_index <-
        which.max(
          id_overlap_n
        )
    }
  } else {
    id_column_index <- which.max(
      id_overlap_n
    )
  }

  id_overlap_best <- id_overlap_n[
    id_column_index
  ]

  if (
    id_overlap_best <
      max(
        1000L,
        floor(
          0.02 *
            length(
              expression_probe_ids
            )
        )
      )
  ) {
    stop(
      "No annotation ID column had adequate overlap with expression probes for ",
      source_label,
      ". Best overlap: ",
      id_overlap_best
    )
  }

  id_column <- names(annotation_table)[
    id_column_index
  ]

  annotation_probe_ids <- trimws(
    gsub(
      '^"|"$',
      "",
      as.character(
        annotation_table[[
          id_column_index
        ]]
      )
    )
  )

  normalized_candidates <- toupper(
    gsub(
      "[^A-Z0-9]+",
      "_",
      symbol_candidates
    )
  )

  symbol_candidate_index <- which(
    normalized_names %in%
      normalized_candidates |
      grepl(
        "SYMBOL",
        normalized_names,
        fixed = TRUE
      ) |
      normalized_names %in%
        c(
          "ILMN_GENE",
          "GENE",
          "GENE_NAME"
        )
  )

  symbol_candidate_index <- setdiff(
    symbol_candidate_index,
    id_column_index
  )

  if (length(symbol_candidate_index) == 0L) {
    stop(
      "No plausible gene-symbol column was found for ",
      source_label,
      "."
    )
  }

  expression_match_index <- match(
    expression_probe_ids,
    annotation_probe_ids
  )

  symbol_audit_rows <- list()

  for (column_index in symbol_candidate_index) {
    raw_symbol <- as.character(
      annotation_table[[
        column_index
      ]]
    )

    mapped_raw_symbol <- raw_symbol[
      expression_match_index
    ]

    cleaned_symbol <- clean_single_gene_symbol(
      mapped_raw_symbol
    )

    symbol_audit_rows[[
      as.character(
        column_index
      )
    ]] <- data.frame(
      annotation_source =
        source_label,
      id_column =
        id_column,
      symbol_column =
        names(annotation_table)[
          column_index
        ],
      probe_overlap_n =
        id_overlap_best,
      accepted_probe_n =
        sum(
          !is.na(
            cleaned_symbol
          )
        ),
      unique_gene_n =
        length(
          unique(
            cleaned_symbol[
              !is.na(
                cleaned_symbol
              )
            ]
          )
        ),
      stringsAsFactors = FALSE
    )
  }

  symbol_audit <- do.call(
    rbind,
    symbol_audit_rows
  )

  best_row <- symbol_audit[
    order(
      -symbol_audit$unique_gene_n,
      -symbol_audit$accepted_probe_n
    ),
    ,
    drop = FALSE
  ][
    1L,
    ,
    drop = FALSE
  ]

  best_symbol_column <-
    best_row$symbol_column[[1]]

  platform_map <- data.frame(
    probe_id =
      annotation_probe_ids,
    raw_gene_symbol =
      as.character(
        annotation_table[[
          best_symbol_column
        ]]
      ),
    stringsAsFactors = FALSE
  )

  platform_map$gene_symbol <-
    clean_single_gene_symbol(
      platform_map$raw_gene_symbol
    )

  platform_map <- platform_map[
    !is.na(
      platform_map$probe_id
    ) &
      nzchar(
        platform_map$probe_id
      ),
    ,
    drop = FALSE
  ]

  platform_map <- platform_map[
    !duplicated(
      platform_map$probe_id
    ),
    ,
    drop = FALSE
  ]

  list(
    platform_map =
      platform_map,
    audit =
      symbol_audit,
    selected_id_column =
      id_column,
    selected_symbol_column =
      best_symbol_column,
    selected_unique_gene_n =
      best_row$unique_gene_n[[1]],
    selected_accepted_probe_n =
      best_row$accepted_probe_n[[1]],
    source_label =
      source_label,
    annotation_row_n =
      nrow(annotation_table)
  )
}

transform_expression_if_needed <- function(
    expression_value
) {
  finite_values <- as.numeric(
    expression_value[
      is.finite(
        expression_value
      )
    ]
  )

  if (length(finite_values) == 0L) {
    stop(
      "Expression matrix has no finite values."
    )
  }

  quantiles <- quantile(
    finite_values,
    probs = c(
      0.01,
      0.5,
      0.99
    ),
    na.rm = TRUE,
    names = FALSE
  )

  log2_required <-
    quantiles[[3]] > 100 ||
    (
      quantiles[[3]] -
        quantiles[[1]]
    ) > 50

  transformed <- expression_value

  if (log2_required) {
    transformed[
      transformed < 0
    ] <- NA_real_

    transformed <- log2(
      transformed + 1
    )
  }

  list(
    expression =
      transformed,
    log2_transform_applied =
      log2_required,
    q01 =
      quantiles[[1]],
    median =
      quantiles[[2]],
    q99 =
      quantiles[[3]]
  )
}

build_gene_matrix <- function(
    dataset_id,
    series_path,
    platform_path,
    sample_design,
    symbol_candidates,
    supplementary_annotation_path =
      NA_character_
) {
  log_event(
    "Building bulk gene matrix for",
    dataset_id
  )

  series_value <- read_series_expression(
    series_path
  )

  transformed_value <-
    transform_expression_if_needed(
      series_value$expression
    )

  expression_value <-
    transformed_value$expression

  if (!setequal(
    colnames(expression_value),
    sample_design$sample_accession
  )) {
    stop(
      dataset_id,
      " expression sample headers do not match the locked sample design."
    )
  }

  expression_value <- expression_value[
    ,
    sample_design$sample_accession,
    drop = FALSE
  ]

  platform_value <- read_platform_table(
    platform_path
  )

  annotation_candidates <- list(
    GEO_PLATFORM_TABLE =
      evaluate_annotation_table(
        annotation_table =
          platform_value,
        expression_probe_ids =
          rownames(
            expression_value
          ),
        symbol_candidates =
          symbol_candidates,
        source_label =
          "GEO_PLATFORM_TABLE"
      )
  )

  if (
    !is.na(
      supplementary_annotation_path
    ) &&
      nzchar(
        supplementary_annotation_path
      ) &&
      file.exists(
        supplementary_annotation_path
      )
  ) {
    supplementary_value <-
      read_tabular_maybe_gz(
        supplementary_annotation_path
      )

    annotation_candidates[[
      "AUTHOR_SUPPLEMENTARY_PROFILE"
    ]] <- evaluate_annotation_table(
      annotation_table =
        supplementary_value,
      expression_probe_ids =
        rownames(
          expression_value
        ),
      symbol_candidates =
        c(
          symbol_candidates,
          "Gene Symbol",
          "Gene_Symbol",
          "SYMBOL",
          "ILMN_Gene",
          "Gene"
        ),
      source_label =
        "AUTHOR_SUPPLEMENTARY_PROFILE"
    )
  }

  annotation_candidate_summary <- do.call(
    rbind,
    lapply(
      annotation_candidates,
      function(candidate_value) {
        data.frame(
          dataset_id =
            dataset_id,
          annotation_source =
            candidate_value$source_label,
          annotation_row_n =
            candidate_value$annotation_row_n,
          selected_id_column =
            candidate_value$selected_id_column,
          selected_symbol_column =
            candidate_value$selected_symbol_column,
          accepted_probe_n =
            candidate_value$selected_accepted_probe_n,
          unique_gene_n =
            candidate_value$selected_unique_gene_n,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  selected_source <- annotation_candidate_summary$annotation_source[
    order(
      -annotation_candidate_summary$unique_gene_n,
      -annotation_candidate_summary$accepted_probe_n
    )
  ][[1]]

  selected_candidate <- annotation_candidates[[
    selected_source
  ]]

  platform_map <-
    selected_candidate$platform_map

  probe_map <- platform_map[
    match(
      rownames(
        expression_value
      ),
      platform_map$probe_id
    ),
    ,
    drop = FALSE
  ]

  unique_symbol_pass <-
    !is.na(
      probe_map$gene_symbol
    ) &
    nzchar(
      probe_map$gene_symbol
    )

  mapped_expression_all <- expression_value[
    unique_symbol_pass,
    ,
    drop = FALSE
  ]

  mapped_symbols_all <- probe_map$gene_symbol[
    unique_symbol_pass
  ]

  complete_probe_pass <- rowSums(
    is.finite(
      mapped_expression_all
    )
  ) ==
    ncol(
      mapped_expression_all
    )

  mapped_expression <- mapped_expression_all[
    complete_probe_pass,
    ,
    drop = FALSE
  ]

  mapped_symbols <- mapped_symbols_all[
    complete_probe_pass
  ]

  mapping_audit_path <- file.path(
    dirs$qc,
    paste0(
      "STEP04D_v3_",
      dataset_id,
      "_ANNOTATION_SOURCE_COMPARISON.csv"
    )
  )

  write_csv_utf8(
    annotation_candidate_summary,
    mapping_audit_path
  )

  symbol_column_audit <- do.call(
    rbind,
    lapply(
      annotation_candidates,
      function(candidate_value) {
        candidate_value$audit
      }
    )
  )

  write_csv_utf8(
    symbol_column_audit,
    file.path(
      dirs$qc,
      paste0(
        "STEP04D_v3_",
        dataset_id,
        "_SYMBOL_COLUMN_AUDIT.csv"
      )
    )
  )

  if (
    nrow(mapped_expression) < 5000L ||
      length(
        unique(
          mapped_symbols
        )
      ) < 5000L
  ) {
    stop(
      dataset_id,
      " yielded fewer than 5,000 uniquely symbol-mapped probes or genes after evaluating all available annotation sources. Review: ",
      mapping_audit_path
    )
  }

  gene_expression_before_model_filter <- limma::avereps(
    mapped_expression,
    ID = mapped_symbols
  )

  gene_complete_pass <- rowSums(
    is.finite(
      gene_expression_before_model_filter
    )
  ) ==
    ncol(
      gene_expression_before_model_filter
    )

  gene_variance <- apply(
    gene_expression_before_model_filter,
    1L,
    stats::var
  )

  gene_variable_pass <- is.finite(
    gene_variance
  ) &
    gene_variance > 0

  gene_model_pass <- gene_complete_pass &
    gene_variable_pass

  gene_expression <- gene_expression_before_model_filter[
    gene_model_pass,
    ,
    drop = FALSE
  ]

  gene_expression <- gene_expression[
    order(
      rownames(
        gene_expression
      )
    ),
    ,
    drop = FALSE
  ]

  qc <- data.frame(
    dataset_id =
      dataset_id,
    input_probe_n =
      nrow(expression_value),
    selected_annotation_source =
      selected_source,
    selected_annotation_row_n =
      selected_candidate$annotation_row_n,
    matched_unique_symbol_probe_n_before_complete_filter =
      nrow(
        mapped_expression_all
      ),
    complete_mapped_probe_n =
      nrow(
        mapped_expression
      ),
    incomplete_mapped_probe_removed_n =
      sum(
        !complete_probe_pass
      ),
    gene_n_before_model_filter =
      nrow(
        gene_expression_before_model_filter
      ),
    incomplete_gene_removed_n =
      sum(
        !gene_complete_pass
      ),
    invariant_gene_removed_n =
      sum(
        gene_complete_pass &
          !gene_variable_pass
      ),
    gene_n =
      nrow(gene_expression),
    sample_n =
      ncol(gene_expression),
    log2_transform_applied =
      transformed_value$log2_transform_applied,
    input_q01 =
      transformed_value$q01,
    input_median =
      transformed_value$median,
    input_q99 =
      transformed_value$q99,
    platform_id_column =
      selected_candidate$selected_id_column,
    platform_symbol_column =
      selected_candidate$selected_symbol_column,
    stringsAsFactors = FALSE
  )

  complete_expression_accepted_full <- rep(
    FALSE,
    length(
      unique_symbol_pass
    )
  )

  complete_expression_accepted_full[
    unique_symbol_pass
  ] <- complete_probe_pass

  list(
    gene_expression =
      gene_expression,
    qc =
      qc,
    probe_map =
      data.frame(
        probe_id =
          rownames(
            expression_value
          ),
        gene_symbol =
          probe_map$gene_symbol,
        unique_symbol_accepted =
          unique_symbol_pass,
        complete_expression_accepted =
          complete_expression_accepted_full,
        selected_annotation_source =
          selected_source,
        stringsAsFactors = FALSE
      ),
    annotation_candidate_summary =
      annotation_candidate_summary,
    symbol_column_audit =
      symbol_column_audit
  )
}

fit_limma_contrasts <- function(
    expression_value,
    design_matrix,
    contrast_matrix,
    dataset_id
) {
  if (
    ncol(expression_value) !=
      nrow(design_matrix)
  ) {
    stop(
      dataset_id,
      " expression/design dimensions do not match."
    )
  }

  complete_row_pass <- rowSums(
    is.finite(
      expression_value
    )
  ) ==
    ncol(
      expression_value
    )

  row_variance <- apply(
    expression_value,
    1L,
    stats::var
  )

  variable_row_pass <- is.finite(
    row_variance
  ) &
    row_variance > 0

  model_row_pass <- complete_row_pass &
    variable_row_pass

  model_expression <- expression_value[
    model_row_pass,
    ,
    drop = FALSE
  ]

  if (nrow(model_expression) == 0L) {
    stop(
      dataset_id,
      " has no complete variable rows for limma modeling."
    )
  }

  log_event(
    "limma input",
    dataset_id,
    "rows:",
    nrow(expression_value),
    "retained:",
    nrow(model_expression),
    "removed_nonfinite:",
    sum(
      !complete_row_pass
    ),
    "removed_invariant:",
    sum(
      complete_row_pass &
        !variable_row_pass
    )
  )

  fit <- limma::lmFit(
    model_expression,
    design_matrix
  )

  fit <- limma::contrasts.fit(
    fit,
    contrast_matrix
  )

  finite_sigma_pass <- is.finite(
    fit$sigma
  ) &
    fit$sigma > 0 &
    is.finite(
      fit$df.residual
    ) &
    fit$df.residual > 0

  if (!all(
    finite_sigma_pass
  )) {
    log_event(
      "limma post-fit filtering",
      dataset_id,
      "removed_invalid_sigma_or_df:",
      sum(
        !finite_sigma_pass
      )
    )

    fit <- fit[
      finite_sigma_pass,
      ]
  }

  if (nrow(fit$coefficients) == 0L) {
    stop(
      dataset_id,
      " has no valid rows after limma residual-variance QC."
    )
  }

  fit <- limma::eBayes(
    fit,
    trend = TRUE,
    robust = TRUE
  )

  result_rows <- list()

  for (contrast_name in colnames(
    contrast_matrix
  )) {
    table_value <- limma::topTable(
      fit,
      coef = contrast_name,
      number = Inf,
      sort.by = "none"
    )

    table_value$gene_symbol <-
      rownames(table_value)

    table_value$dataset_id <-
      dataset_id

    table_value$contrast <-
      contrast_name

    result_rows[[
      contrast_name
    ]] <- table_value
  }

  do.call(
    rbind,
    result_rows
  )
}

row_zscore <- function(expression_value) {
  row_mean <- rowMeans(
    expression_value,
    na.rm = TRUE
  )

  row_sd <- apply(
    expression_value,
    1L,
    sd,
    na.rm = TRUE
  )

  valid <- is.finite(row_sd) &
    row_sd > 0

  output <- matrix(
    NA_real_,
    nrow = nrow(expression_value),
    ncol = ncol(expression_value),
    dimnames = dimnames(expression_value)
  )

  output[
    valid,
    ] <- (
      expression_value[
        valid,
        ,
        drop = FALSE
      ] -
        row_mean[
          valid
        ]
    ) /
      row_sd[
        valid
      ]

  output
}

calculate_program_scores <- function(
    gene_expression,
    candidate_table,
    program_definitions,
    dataset_id
) {
  z_expression <- row_zscore(
    gene_expression
  )

  score_rows <- list()
  coverage_rows <- list()

  for (program_name in names(
    program_definitions
  )) {
    program_genes <- program_definitions[[
      program_name
    ]]

    program_genes <- unique(
      program_genes[
        !is.na(program_genes) &
          nzchar(program_genes)
      ]
    )

    present_genes <- intersect(
      program_genes,
      rownames(z_expression)
    )

    candidate_subset <- candidate_table[
      match(
        present_genes,
        candidate_table$feature_id
      ),
      ,
      drop = FALSE
    ]

    weights <- sign(
      candidate_subset$interaction_logFC
    )

    valid_weight <- is.finite(weights) &
      weights != 0

    present_genes <- present_genes[
      valid_weight
    ]

    weights <- weights[
      valid_weight
    ]

    if (length(present_genes) < 10L) {
      stop(
        dataset_id,
        " program ",
        program_name,
        " has fewer than 10 measurable genes."
      )
    }

    score_value <- colMeans(
      sweep(
        z_expression[
          present_genes,
          ,
          drop = FALSE
        ],
        1L,
        weights,
        "*"
      ),
      na.rm = TRUE
    )

    score_rows[[
      program_name
    ]] <- data.frame(
      dataset_id =
        dataset_id,
      sample_accession =
        colnames(gene_expression),
      program_name =
        program_name,
      score =
        as.numeric(score_value),
      measurable_gene_n =
        length(present_genes),
      stringsAsFactors = FALSE
    )

    coverage_rows[[
      program_name
    ]] <- data.frame(
      dataset_id =
        dataset_id,
      program_name =
        program_name,
      frozen_gene_n =
        length(program_genes),
      measurable_gene_n =
        length(present_genes),
      coverage_fraction =
        length(present_genes) /
        length(program_genes),
      stringsAsFactors = FALSE
    )
  }

  list(
    scores = do.call(
      rbind,
      score_rows
    ),
    coverage = do.call(
      rbind,
      coverage_rows
    )
  )
}

fit_score_model <- function(
    score_table,
    sample_design,
    design_matrix,
    contrast_matrix,
    dataset_id
) {
  programs <- unique(
    score_table$program_name
  )

  result_rows <- list()

  for (program_name in programs) {
    program_scores <- score_table[
      score_table$program_name ==
        program_name,
      ,
      drop = FALSE
    ]

    program_scores <- program_scores[
      match(
        sample_design$sample_accession,
        program_scores$sample_accession
      ),
      ,
      drop = FALSE
    ]

    if (any(is.na(
      program_scores$score
    ))) {
      stop(
        dataset_id,
        " score alignment failed for ",
        program_name
      )
    }

    score_matrix <- matrix(
      program_scores$score,
      nrow = 1L,
      dimnames = list(
        program_name,
        sample_design$sample_accession
      )
    )

    fitted <- fit_limma_contrasts(
      score_matrix,
      design_matrix,
      contrast_matrix,
      dataset_id
    )

    fitted$program_name <-
      fitted$gene_symbol

    fitted$gene_symbol <- NULL

    result_rows[[
      program_name
    ]] <- fitted
  }

  do.call(
    rbind,
    result_rows
  )
}

read_text_table_maybe_gz <- function(
    path,
    header = FALSE
) {
  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rt"
    )
  } else {
    file(
      path,
      open = "rt"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  read.delim(
    connection,
    header = header,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

read_mtx_maybe_gz <- function(path) {
  connection <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      open = "rb"
    )
  } else {
    file(
      path,
      open = "rb"
    )
  }

  on.exit(
    close(connection),
    add = TRUE
  )

  Matrix::readMM(
    connection
  )
}

prepare_single_cell_sample <- function(
    dataset_id,
    sample_accession,
    matrix_path,
    barcode_path,
    feature_path,
    output_root
) {
  log_event(
    "Single-cell pre-QC:",
    dataset_id,
    sample_accession
  )

  matrix_value <- read_mtx_maybe_gz(
    matrix_path
  )

  feature_table <- read_text_table_maybe_gz(
    feature_path,
    header = FALSE
  )

  barcode_table <- read_text_table_maybe_gz(
    barcode_path,
    header = FALSE
  )

  feature_id <- as.character(
    feature_table[[1]]
  )

  gene_symbol <- if (
    ncol(feature_table) >= 2L
  ) {
    as.character(
      feature_table[[2]]
    )
  } else {
    feature_id
  }

  gene_symbol[
    is.na(gene_symbol) |
      !nzchar(
        trimws(gene_symbol)
      )
  ] <- feature_id[
    is.na(gene_symbol) |
      !nzchar(
        trimws(gene_symbol)
      )
  ]

  barcodes <- as.character(
    barcode_table[[1]]
  )

  if (
    nrow(matrix_value) ==
      length(feature_id) &&
      ncol(matrix_value) ==
        length(barcodes)
  ) {
    orientation <-
      "FEATURES_BY_CELLS"
  } else if (
    ncol(matrix_value) ==
      length(feature_id) &&
      nrow(matrix_value) ==
        length(barcodes)
  ) {
    matrix_value <- Matrix::t(
      matrix_value
    )

    orientation <-
      "TRANSPOSED_TO_FEATURES_BY_CELLS"
  } else {
    stop(
      dataset_id,
      " ",
      sample_accession,
      " matrix dimensions do not match features and barcodes."
    )
  }

  matrix_value <- methods::as(
    matrix_value,
    "dgCMatrix"
  )

  unique_feature_names <- make.unique(
    gene_symbol,
    sep = "__"
  )

  rownames(
    matrix_value
  ) <- unique_feature_names

  prefixed_barcodes <- paste0(
    sample_accession,
    "__",
    barcodes
  )

  colnames(
    matrix_value
  ) <- prefixed_barcodes

  total_counts <- Matrix::colSums(
    matrix_value
  )

  detected_features <- diff(
    matrix_value@p
  )

  mito_index <- grepl(
    "^MT-",
    gene_symbol,
    ignore.case = TRUE
  )

  mito_counts <- if (any(mito_index)) {
    Matrix::colSums(
      matrix_value[
        mito_index,
        ,
        drop = FALSE
      ]
    )
  } else {
    rep(
      0,
      ncol(matrix_value)
    )
  }

  percent_mito <- ifelse(
    total_counts > 0,
    100 *
      mito_counts /
      total_counts,
    NA_real_
  )

  upper_feature_limit <- min(
    7500,
    max(
      2500,
      median(
        detected_features,
        na.rm = TRUE
      ) +
        4 *
        stats::mad(
          detected_features,
          na.rm = TRUE
        )
    )
  )

  keep_cell <-
    total_counts >= 500 &
    detected_features >= 200 &
    detected_features <=
      upper_feature_limit &
    percent_mito <= 25

  filtered_matrix <- matrix_value[
    ,
    keep_cell,
    drop = FALSE
  ]

  cell_metadata <- data.frame(
    cell_barcode =
      prefixed_barcodes,
    dataset_id =
      dataset_id,
    sample_accession =
      sample_accession,
    total_counts =
      as.numeric(
        total_counts
      ),
    detected_features =
      as.numeric(
        detected_features
      ),
    percent_mito =
      as.numeric(
        percent_mito
      ),
    preQC_keep =
      keep_cell,
    stringsAsFactors = FALSE
  )

  sample_object <- list(
    counts =
      filtered_matrix,
    feature_table = data.frame(
      feature_id =
        feature_id,
      gene_symbol =
        gene_symbol,
      matrix_row_name =
        unique_feature_names,
      stringsAsFactors = FALSE
    ),
    cell_metadata =
      cell_metadata[
        keep_cell,
        ,
        drop = FALSE
      ],
    sample_metadata = data.frame(
      dataset_id =
        dataset_id,
      sample_accession =
        sample_accession,
      stringsAsFactors = FALSE
    ),
    preQC_parameters = list(
      minimum_counts = 500,
      minimum_features = 200,
      maximum_features =
        upper_feature_limit,
      maximum_percent_mito = 25
    )
  )

  dataset_output_root <- file.path(
    output_root,
    dataset_id
  )

  dir.create(
    dataset_output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_path <- file.path(
    dataset_output_root,
    paste0(
      sample_accession,
      "_preQC.rds"
    )
  )

  saveRDS(
    sample_object,
    output_path,
    compress = TRUE
  )

  data.frame(
    dataset_id =
      dataset_id,
    sample_accession =
      sample_accession,
    orientation =
      orientation,
    feature_n =
      nrow(matrix_value),
    raw_cell_n =
      ncol(matrix_value),
    retained_cell_n =
      ncol(filtered_matrix),
    retained_fraction =
      ncol(filtered_matrix) /
      ncol(matrix_value),
    median_counts =
      median(
        total_counts,
        na.rm = TRUE
      ),
    median_features =
      median(
        detected_features,
        na.rm = TRUE
      ),
    median_percent_mito =
      median(
        percent_mito,
        na.rm = TRUE
      ),
    adaptive_max_features =
      upper_feature_limit,
    output_path =
      normalize_slashes(
        output_path
      ),
    stringsAsFactors = FALSE
  )
}

log_event(
  "Starting",
  run_id
)

planned_path_examples <- c(
  file.path(
    bulk_result_root,
    "loso_gene_results.csv"
  ),
  file.path(
    bulk_result_root,
    "cross_context_candidates.csv"
  ),
  file.path(
    sc_object_root,
    "GSE250130",
    "GSM0000000_preQC.rds"
  )
)

path_length_audit <- data.frame(
  path =
    normalize_slashes(
      planned_path_examples
    ),
  path_length =
    nchar(
      planned_path_examples
    ),
  under_soft_limit =
    nchar(
      planned_path_examples
    ) < 220L,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  path_length_audit,
  file.path(
    dirs$qc,
    "S4D_path_length_audit.csv"
  )
)

if (!all(
  path_length_audit$under_soft_limit
)) {
  stop(
    "At least one planned STEP04D path exceeds the 220-character soft limit."
  )
}

# ----------------------------------------------------------------------
# Validate locked upstream artifacts
# ----------------------------------------------------------------------

step04C_gate_path <- file.path(
  dirs$results,
  "STEP04C_GATE_TABLE.csv"
)

step04C_locked_root <- file.path(
  dirs$results,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0"
)

step04C_manifest_path <- file.path(
  step04C_locked_root,
  "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv"
)

step04B_locked_root <- file.path(
  dirs$results,
  "STEP04B_locked_combined_inputs",
  "COMBINED_INPUT_SOURCE_v1.0"
)

source_reference_path <- file.path(
  step04B_locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_REFERENCE_LOCK.csv"
)

GSE250130_bundle_path <- file.path(
  step04B_locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_GSE250130_BUNDLE_MAP.csv"
)

GSE183837_bundle_path <- file.path(
  step04B_locked_root,
  "COMBINED_INPUT_SOURCE_v1.0_GSE183837_BUNDLE_MAP.csv"
)

platform_lock_path <- file.path(
  dirs$results,
  "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"
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

routing_path <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons",
  "06_results",
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  "RNA_OCCUPANCY_EVIDENCE_v1.0",
  "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv"
)

required_upstream <- c(
  step04C_gate_path,
  step04C_manifest_path,
  source_reference_path,
  GSE250130_bundle_path,
  GSE183837_bundle_path,
  platform_lock_path,
  candidate_path,
  routing_path
)

missing_upstream <- required_upstream[
  !file.exists(required_upstream)
]

if (length(missing_upstream) > 0L) {
  stop(
    "Required locked upstream artifact is missing:\n",
    paste(
      missing_upstream,
      collapse = "\n"
    )
  )
}

step04C_gate <- read.csv(
  step04C_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

accepted_gate <- step04C_gate$status[
  step04C_gate$gate_component ==
    "STEP04C gate"
]

if (
  length(accepted_gate) != 1L ||
    accepted_gate !=
      "READY_FOR_COMBINED_BULK_AND_SINGLE_CELL_ANALYSIS_BUILD"
) {
  stop(
    "STEP04C does not authorize the combined analysis build."
  )
}

step04C_manifest <- read.csv(
  step04C_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

step04C_manifest$file_exists <- file.exists(
  step04C_manifest$path
)

step04C_manifest$observed_md5 <- vapply(
  step04C_manifest$path,
  md5_file,
  character(1)
)

step04C_manifest$observed_sha256 <- vapply(
  step04C_manifest$path,
  sha256_file,
  character(1)
)

step04C_manifest$identity_pass <- with(
  step04C_manifest,
  file_exists &
    md5 ==
      observed_md5 &
    sha256 ==
      observed_sha256
)

if (!all(
  step04C_manifest$identity_pass
)) {
  stop(
    "At least one STEP04C locked artifact failed hash validation."
  )
}

read_locked_design <- function(dataset_id) {
  read.csv(
    file.path(
      step04C_locked_root,
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

design_250130 <- read_locked_design(
  "GSE250130"
)

design_165004 <- read_locked_design(
  "GSE165004"
)

design_92324 <- read_locked_design(
  "GSE92324"
)

design_111974 <- read_locked_design(
  "GSE111974"
)

design_183837 <- read_locked_design(
  "GSE183837"
)

source_reference <- read.csv(
  source_reference_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

platform_lock <- read.csv(
  platform_lock_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_table <- read.csv(
  candidate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

routing_table <- read.csv(
  routing_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(candidate_table) != 436L ||
    length(
      unique(
        candidate_table$feature_id
      )
    ) !=
      436L ||
    sum(
      candidate_table$occupancy_supported
    ) !=
      284L
) {
  stop(
    "The locked 436/284 RNA occupancy candidate table failed."
  )
}

# ----------------------------------------------------------------------
# Build three bulk gene-level matrices
# ----------------------------------------------------------------------

get_series_path <- function(dataset_id) {
  rows <- source_reference[
    source_reference$dataset_id ==
      dataset_id &
      grepl(
        "SERIES_MATRIX",
        source_reference$source_role
      ),
    ,
    drop = FALSE
  ]

  if (nrow(rows) != 1L) {
    stop(
      "Expected one series matrix for ",
      dataset_id
    )
  }

  rows$path[[1]]
}

get_platform_path <- function(platform_id) {
  rows <- platform_lock[
    platform_lock$platform_id ==
      platform_id,
    ,
    drop = FALSE
  ]

  if (nrow(rows) != 1L) {
    stop(
      "Expected one platform source for ",
      platform_id
    )
  }

  rows$path[[1]]
}

download_with_resume(
  url =
    GSE92324_profile_url,
  destination =
    GSE92324_profile_path
)

GSE92324_profile_source_audit <- data.frame(
  dataset_id =
    "GSE92324",
  source_role =
    "AUTHOR_SUPPLEMENTARY_NORMALIZED_PROFILE_WITH_GENE_ANNOTATION",
  source_url =
    GSE92324_profile_url,
  path =
    normalize_slashes(
      GSE92324_profile_path
    ),
  size_bytes = as.numeric(
    file.info(
      GSE92324_profile_path
    )$size
  ),
  md5 =
    md5_file(
      GSE92324_profile_path
    ),
  sha256 =
    sha256_file(
      GSE92324_profile_path
    ),
  source_status =
    "CANDIDATE_FOR_PROBE_GENE_MAPPING",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  GSE92324_profile_source_audit,
  file.path(
    dirs$results,
    "STEP04D_v3_GSE92324_SUPPLEMENTARY_PROFILE_SOURCE.csv"
  )
)

bulk_specs <- list(
  GSE165004 = list(
    series_path =
      get_series_path(
        "GSE165004"
      ),
    platform_path =
      get_platform_path(
        "GPL16699"
      ),
    design =
      design_165004,
    symbol_candidates =
      c(
        "GENE_SYMBOL",
        "Gene Symbol",
        "Symbol"
      ),
    supplementary_annotation_path =
      NA_character_
  ),
  GSE92324 = list(
    series_path =
      get_series_path(
        "GSE92324"
      ),
    platform_path =
      get_platform_path(
        "GPL10558"
      ),
    design =
      design_92324,
    symbol_candidates =
      c(
        "Symbol",
        "ILMN_Gene",
        "GENE_SYMBOL"
      ),
    supplementary_annotation_path =
      GSE92324_profile_path
  ),
  GSE111974 = list(
    series_path =
      get_series_path(
        "GSE111974"
      ),
    platform_path =
      get_platform_path(
        "GPL17077"
      ),
    design =
      design_111974,
    symbol_candidates =
      c(
        "GENE_SYMBOL",
        "Gene Symbol",
        "Symbol"
      ),
    supplementary_annotation_path =
      NA_character_
  )
)

bulk_objects <- list()
bulk_qc_rows <- list()

previous_bulk_qc_path <- file.path(
  dirs$results,
  "STEP04D_BULK_GENE_MATRIX_QC.csv"
)

previous_matrix_paths <- file.path(
  previous_bulk_matrix_root,
  paste0(
    "BULK_GENE_MATRIX_v1.1_",
    names(
      bulk_specs
    ),
    ".rds"
  )
)

reuse_previous_bulk_matrices <-
  all(
    file.exists(
      previous_matrix_paths
    )
  ) &&
  file.exists(
    previous_bulk_qc_path
  )

if (reuse_previous_bulk_matrices) {
  log_event(
    "Reusing three complete v1.1 bulk matrices from the interrupted v5 run."
  )

  previous_bulk_qc <- read.csv(
    previous_bulk_qc_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!setequal(
    previous_bulk_qc$dataset_id,
    names(
      bulk_specs
    )
  )) {
    stop(
      "The previous bulk-matrix QC table does not contain all three datasets."
    )
  }

  for (dataset_id in names(
    bulk_specs
  )) {
    previous_path <- file.path(
      previous_bulk_matrix_root,
      paste0(
        "BULK_GENE_MATRIX_v1.1_",
        dataset_id,
        ".rds"
      )
    )

    bulk_objects[[dataset_id]] <-
      readRDS(
        previous_path
      )

    if (
      !is.matrix(
        bulk_objects[[dataset_id]]
      ) ||
        nrow(
          bulk_objects[[dataset_id]]
        ) < 5000L
    ) {
      stop(
        "Reusable bulk matrix failed validation: ",
        dataset_id
      )
    }

    new_matrix_path <- file.path(
      bulk_matrix_root,
      paste0(
        "BGM_v1_2_",
        dataset_id,
        ".rds"
      )
    )

    saveRDS(
      bulk_objects[[dataset_id]],
      new_matrix_path,
      compress = TRUE
    )

    bulk_qc_rows[[dataset_id]] <-
      previous_bulk_qc[
        previous_bulk_qc$dataset_id ==
          dataset_id,
        ,
        drop = FALSE
      ]

    bulk_qc_rows[[dataset_id]]$matrix_source <-
      "REUSED_FROM_INTERRUPTED_V5_BULK_GENE_MATRIX_v1.1"
  }
} else {
  log_event(
    "Reusable v1.1 bulk matrices were not all available; rebuilding all three matrices."
  )

  for (dataset_id in names(
    bulk_specs
  )) {
    spec <- bulk_specs[[
      dataset_id
    ]]

    built <- build_gene_matrix(
      dataset_id =
        dataset_id,
      series_path =
        spec$series_path,
      platform_path =
        spec$platform_path,
      sample_design =
        spec$design,
      symbol_candidates =
        spec$symbol_candidates,
      supplementary_annotation_path =
        spec$supplementary_annotation_path
    )

    bulk_objects[[dataset_id]] <-
      built$gene_expression

    built$qc$matrix_source <-
      "REBUILT_IN_STEP04D_V6"

    bulk_qc_rows[[dataset_id]] <-
      built$qc

    saveRDS(
      built$gene_expression,
      file.path(
        bulk_matrix_root,
        paste0(
          "BGM_v1_2_",
          dataset_id,
          ".rds"
        )
      ),
      compress = TRUE
    )

    write_csv_utf8(
      built$probe_map,
      file.path(
        dirs$qc,
        paste0(
          "S4D_",
          dataset_id,
          "_probe_gene_map.csv"
        )
      )
    )
  }
}

bulk_matrix_qc <- do.call(
  rbind,
  bulk_qc_rows
)

write_csv_utf8(
  bulk_matrix_qc,
  file.path(
    dirs$results,
    "S4D_bulk_matrix_qc.csv"
  )
)

# ----------------------------------------------------------------------
# Prespecified program scores
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

score_tables <- list()
coverage_tables <- list()

for (dataset_id in names(
  bulk_objects
)) {
  score_value <- calculate_program_scores(
    gene_expression =
      bulk_objects[[
        dataset_id
      ]],
    candidate_table =
      candidate_table,
    program_definitions =
      program_definitions,
    dataset_id =
      dataset_id
  )

  score_tables[[dataset_id]] <-
    score_value$scores

  coverage_tables[[dataset_id]] <-
    score_value$coverage
}

all_scores <- do.call(
  rbind,
  score_tables
)

all_coverage <- do.call(
  rbind,
  coverage_tables
)

write_csv_utf8(
  all_scores,
  file.path(
    bulk_result_root,
    "sample_scores.csv"
  )
)

write_csv_utf8(
  all_coverage,
  file.path(
    bulk_result_root,
    "gene_coverage.csv"
  )
)

# ----------------------------------------------------------------------
# Gene-level and score-level models
# ----------------------------------------------------------------------

# GSE165004
group_165004 <- factor(
  design_165004$group,
  levels = c(
    "CONTROL",
    "RPL",
    "UIF"
  )
)

design_matrix_165004 <- model.matrix(
  ~ 0 + group_165004
)

colnames(
  design_matrix_165004
) <- levels(
  group_165004
)

contrast_matrix_165004 <- limma::makeContrasts(
  RPL_vs_CONTROL =
    RPL - CONTROL,
  UIF_vs_CONTROL =
    UIF - CONTROL,
  levels =
    design_matrix_165004
)

gene_results_165004 <- fit_limma_contrasts(
  bulk_objects[[
    "GSE165004"
  ]],
  design_matrix_165004,
  contrast_matrix_165004,
  "GSE165004"
)

score_results_165004 <- fit_score_model(
  score_tables[[
    "GSE165004"
  ]],
  design_165004,
  design_matrix_165004,
  contrast_matrix_165004,
  "GSE165004"
)

# GSE92324
group_92324 <- factor(
  design_92324$group,
  levels = c(
    "CONTROL",
    "IMPLANTATION_FAILURE"
  )
)

day_92324 <- factor(
  design_92324$collection_day
)

design_matrix_92324 <- model.matrix(
  ~ 0 + group_92324 + day_92324
)

colnames(
  design_matrix_92324
) <- make.names(
  colnames(
    design_matrix_92324
  )
)

contrast_matrix_92324 <- limma::makeContrasts(
  IF_vs_CONTROL =
    group_92324IMPLANTATION_FAILURE -
    group_92324CONTROL,
  levels =
    design_matrix_92324
)

gene_results_92324 <- fit_limma_contrasts(
  bulk_objects[[
    "GSE92324"
  ]],
  design_matrix_92324,
  contrast_matrix_92324,
  "GSE92324"
)

score_results_92324 <- fit_score_model(
  score_tables[[
    "GSE92324"
  ]],
  design_92324,
  design_matrix_92324,
  contrast_matrix_92324,
  "GSE92324"
)

# GSE111974
group_111974 <- factor(
  design_111974$group,
  levels = c(
    "CONTROL",
    "RIF"
  )
)

design_matrix_111974 <- model.matrix(
  ~ 0 + group_111974
)

colnames(
  design_matrix_111974
) <- levels(
  group_111974
)

contrast_matrix_111974 <- limma::makeContrasts(
  RIF_vs_CONTROL =
    RIF - CONTROL,
  levels =
    design_matrix_111974
)

gene_results_111974 <- fit_limma_contrasts(
  bulk_objects[[
    "GSE111974"
  ]],
  design_matrix_111974,
  contrast_matrix_111974,
  "GSE111974"
)

score_results_111974 <- fit_score_model(
  score_tables[[
    "GSE111974"
  ]],
  design_111974,
  design_matrix_111974,
  contrast_matrix_111974,
  "GSE111974"
)

all_gene_results <- rbind(
  gene_results_165004,
  gene_results_92324,
  gene_results_111974
)

all_score_results <- rbind(
  score_results_165004,
  score_results_92324,
  score_results_111974
)

write_csv_utf8(
  all_gene_results,
  file.path(
    bulk_result_root,
    "full_gene_results.csv"
  )
)

write_csv_utf8(
  all_score_results,
  file.path(
    bulk_result_root,
    "score_models.csv"
  )
)

candidate_gene_results <- merge(
  candidate_table,
  all_gene_results,
  by.x = "feature_id",
  by.y = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)

write_csv_utf8(
  candidate_gene_results,
  file.path(
    bulk_result_root,
    "candidate_gene_results.csv"
  )
)

# ----------------------------------------------------------------------
# GSE165004 leave-one-slide sensitivity
# ----------------------------------------------------------------------

full_165004_results <- gene_results_165004[
  gene_results_165004$gene_symbol %in%
    candidate_table$feature_id,
  c(
    "gene_symbol",
    "contrast",
    "logFC"
  ),
  drop = FALSE
]

slide_ids <- unique(
  design_165004$slide_id
)

slide_gene_rows <- list()
slide_score_rows <- list()

for (slide_id in slide_ids) {
  keep_samples <-
    design_165004$slide_id !=
    slide_id

  subset_design <- design_165004[
    keep_samples,
    ,
    drop = FALSE
  ]

  subset_group <- factor(
    subset_design$group,
    levels = c(
      "CONTROL",
      "RPL",
      "UIF"
    )
  )

  subset_design_matrix <- model.matrix(
    ~ 0 + subset_group
  )

  colnames(
    subset_design_matrix
  ) <- levels(
    subset_group
  )

  subset_contrast_matrix <-
    limma::makeContrasts(
      RPL_vs_CONTROL =
        RPL - CONTROL,
      UIF_vs_CONTROL =
        UIF - CONTROL,
      levels =
        subset_design_matrix
    )

  subset_gene_result <- fit_limma_contrasts(
    bulk_objects[[
      "GSE165004"
    ]][
      ,
      keep_samples,
      drop = FALSE
    ],
    subset_design_matrix,
    subset_contrast_matrix,
    "GSE165004"
  )

  subset_gene_result <- subset_gene_result[
    subset_gene_result$gene_symbol %in%
      candidate_table$feature_id,
    c(
      "gene_symbol",
      "contrast",
      "logFC",
      "P.Value",
      "adj.P.Val"
    ),
    drop = FALSE
  ]

  subset_gene_result$left_out_slide <-
    slide_id

  slide_gene_rows[[slide_id]] <-
    subset_gene_result

  subset_scores <- score_tables[[
    "GSE165004"
  ]]

  subset_scores <- subset_scores[
    subset_scores$sample_accession %in%
      subset_design$sample_accession,
    ,
    drop = FALSE
  ]

  subset_score_result <- fit_score_model(
    subset_scores,
    subset_design,
    subset_design_matrix,
    subset_contrast_matrix,
    "GSE165004"
  )

  subset_score_result$left_out_slide <-
    slide_id

  slide_score_rows[[slide_id]] <-
    subset_score_result
}

slide_gene_results <- do.call(
  rbind,
  slide_gene_rows
)

slide_score_results <- do.call(
  rbind,
  slide_score_rows
)

write_csv_utf8(
  slide_gene_results,
  file.path(
    bulk_result_root,
    "loso_gene_results.csv"
  )
)

write_csv_utf8(
  slide_score_results,
  file.path(
    bulk_result_root,
    "loso_score_results.csv"
  )
)

slide_robustness_rows <- list()

for (contrast_name in unique(
  full_165004_results$contrast
)) {
  full_subset <- full_165004_results[
    full_165004_results$contrast ==
      contrast_name,
    ,
    drop = FALSE
  ]

  for (slide_id in slide_ids) {
    left_subset <- slide_gene_results[
      slide_gene_results$contrast ==
        contrast_name &
        slide_gene_results$left_out_slide ==
          slide_id,
      ,
      drop = FALSE
    ]

    joined <- merge(
      full_subset,
      left_subset,
      by = c(
        "gene_symbol",
        "contrast"
      ),
      suffixes = c(
        "_full",
        "_leaveout"
      ),
      all = FALSE
    )

    slide_robustness_rows[[
      paste(
        contrast_name,
        slide_id,
        sep = "__"
      )
    ]] <- data.frame(
      contrast =
        contrast_name,
      left_out_slide =
        slide_id,
      candidate_gene_n =
        nrow(joined),
      spearman_rho =
        suppressWarnings(
          cor(
            joined$logFC_full,
            joined$logFC_leaveout,
            method = "spearman",
            use = "pairwise.complete.obs"
          )
        ),
      direction_concordance =
        mean(
          sign(
            joined$logFC_full
          ) ==
            sign(
              joined$logFC_leaveout
            ),
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )
  }
}

slide_robustness <- do.call(
  rbind,
  slide_robustness_rows
)

write_csv_utf8(
  slide_robustness,
  file.path(
    bulk_result_root,
    "loso_robustness.csv"
  )
)

# ----------------------------------------------------------------------
# Cross-context candidate convergence table
# ----------------------------------------------------------------------

compact_candidate_results <- candidate_gene_results[
  ,
  c(
    "feature_id",
    "program_class",
    "interaction_logFC",
    "interaction_FDR",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support",
    "dataset_id",
    "contrast",
    "logFC",
    "P.Value",
    "adj.P.Val"
  ),
  drop = FALSE
]

compact_candidate_results$context_id <- paste(
  compact_candidate_results$dataset_id,
  compact_candidate_results$contrast,
  sep = "__"
)

context_ids <- unique(
  compact_candidate_results$context_id[
    !is.na(
      compact_candidate_results$context_id
    )
  ]
)

candidate_convergence <- unique(
  compact_candidate_results[
    ,
    c(
      "feature_id",
      "program_class",
      "interaction_logFC",
      "interaction_FDR",
      "occupancy_supported",
      "shared_occupancy_supported",
      "any_promoter_overlap_support"
    ),
    drop = FALSE
  ]
)

for (context_id in context_ids) {
  context_rows <- compact_candidate_results[
    compact_candidate_results$context_id ==
      context_id,
    c(
      "feature_id",
      "logFC",
      "P.Value",
      "adj.P.Val"
    ),
    drop = FALSE
  ]

  names(context_rows)[
    -1
  ] <- paste0(
    context_id,
    "__",
    c(
      "logFC",
      "P.Value",
      "adj.P.Val"
    )
  )

  candidate_convergence <- merge(
    candidate_convergence,
    context_rows,
    by = "feature_id",
    all.x = TRUE,
    sort = FALSE
  )
}

logFC_columns <- grep(
  "__logFC$",
  names(candidate_convergence),
  value = TRUE
)

candidate_convergence$measured_context_n <-
  rowSums(
    is.finite(
      as.matrix(
        candidate_convergence[
          ,
          logFC_columns,
          drop = FALSE
        ]
      )
    )
  )

candidate_convergence$positive_context_n <-
  rowSums(
    as.matrix(
      candidate_convergence[
        ,
        logFC_columns,
        drop = FALSE
      ]
    ) > 0,
    na.rm = TRUE
  )

candidate_convergence$negative_context_n <-
  rowSums(
    as.matrix(
      candidate_convergence[
        ,
        logFC_columns,
        drop = FALSE
      ]
    ) < 0,
    na.rm = TRUE
  )

candidate_convergence$cross_context_direction_label <-
  ifelse(
    candidate_convergence$positive_context_n ==
      candidate_convergence$measured_context_n,
    "POSITIVE_ALL_MEASURED_CONTEXTS",
    ifelse(
      candidate_convergence$negative_context_n ==
        candidate_convergence$measured_context_n,
      "NEGATIVE_ALL_MEASURED_CONTEXTS",
      "CONTEXT_DEPENDENT_OR_MIXED_DIRECTION"
    )
  )

write_csv_utf8(
  candidate_convergence,
  file.path(
    bulk_result_root,
    "cross_context_candidates.csv"
  )
)

# ----------------------------------------------------------------------
# Single-cell donor-level pre-QC
# ----------------------------------------------------------------------

GSE250130_bundle_map <- read.csv(
  GSE250130_bundle_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

GSE183837_bundle_map <- read.csv(
  GSE183837_bundle_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

source_183837 <- source_reference[
  source_reference$dataset_id ==
    "GSE183837" &
    source_reference$source_role ==
      "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER",
  ,
  drop = FALSE
]

source_183837$sample_accession <-
  vapply(
    source_183837$path,
    function(path_value) {
      match_value <- regmatches(
        basename(path_value),
        regexpr(
          "GSM[0-9]+",
          basename(path_value),
          ignore.case = TRUE,
          perl = TRUE
        )
      )

      toupper(match_value)
    },
    character(1)
  )

source_183837$component <- vapply(
  source_183837$file_name,
  function(file_name) {
    file_lower <- tolower(file_name)

    if (grepl(
      "matrix.*\\.mtx(\\.gz)?$",
      file_lower,
      perl = TRUE
    )) {
      "matrix"
    } else if (grepl(
      "barcodes.*\\.(tsv|txt)(\\.gz)?$",
      file_lower,
      perl = TRUE
    )) {
      "barcodes"
    } else if (grepl(
      "(genes|features).*\\.(tsv|txt)(\\.gz)?$",
      file_lower,
      perl = TRUE
    )) {
      "features"
    } else {
      "other"
    }
  },
  character(1)
)

sc_qc_rows <- list()
sc_qc_index <- 1L

for (row_index in seq_len(
  nrow(GSE250130_bundle_map)
)) {
  map_row <- GSE250130_bundle_map[
    row_index,
    ,
    drop = FALSE
  ]

  sc_qc_rows[[sc_qc_index]] <-
    prepare_single_cell_sample(
      dataset_id =
        "GSE250130",
      sample_accession =
        map_row$sample_accession,
      matrix_path =
        map_row$matrix_file,
      barcode_path =
        map_row$barcode_file,
      feature_path =
        map_row$feature_file,
      output_root =
        sc_object_root
    )

  sc_qc_index <- sc_qc_index + 1L

  gc()
}

for (sample_accession in unique(
  source_183837$sample_accession
)) {
  sample_rows <- source_183837[
    source_183837$sample_accession ==
      sample_accession,
    ,
    drop = FALSE
  ]

  matrix_path <- sample_rows$path[
    sample_rows$component ==
      "matrix"
  ]

  barcode_path <- sample_rows$path[
    sample_rows$component ==
      "barcodes"
  ]

  feature_path <- sample_rows$path[
    sample_rows$component ==
      "features"
  ]

  if (
    length(matrix_path) != 1L ||
      length(barcode_path) != 1L ||
      length(feature_path) != 1L
  ) {
    stop(
      "GSE183837 triplet path reconstruction failed for ",
      sample_accession
    )
  }

  sc_qc_rows[[sc_qc_index]] <-
    prepare_single_cell_sample(
      dataset_id =
        "GSE183837",
      sample_accession =
        sample_accession,
      matrix_path =
        matrix_path,
      barcode_path =
        barcode_path,
      feature_path =
        feature_path,
      output_root =
        sc_object_root
    )

  sc_qc_index <- sc_qc_index + 1L

  gc()
}

sc_qc_summary <- do.call(
  rbind,
  sc_qc_rows
)

write_csv_utf8(
  sc_qc_summary,
  file.path(
    dirs$results,
    "S4D_sc_sample_qc.csv"
  )
)

sc_dataset_summary <- aggregate(
  cbind(
    raw_cell_n,
    retained_cell_n
  ) ~ dataset_id,
  data = sc_qc_summary,
  FUN = sum
)

sc_dataset_summary$sample_n <- vapply(
  sc_dataset_summary$dataset_id,
  function(dataset_id) {
    sum(
      sc_qc_summary$dataset_id ==
        dataset_id
    )
  },
  integer(1)
)

sc_dataset_summary$retained_fraction <-
  sc_dataset_summary$retained_cell_n /
  sc_dataset_summary$raw_cell_n

write_csv_utf8(
  sc_dataset_summary,
  file.path(
    dirs$results,
    "S4D_sc_dataset_qc.csv"
  )
)

# ----------------------------------------------------------------------
# Final QC and Gate
# ----------------------------------------------------------------------

minimum_candidate_coverage <-
  min(
    all_coverage$measurable_gene_n[
      all_coverage$program_name ==
        "OCCUPANCY_SUPPORTED_284"
    ]
  )

minimum_slide_spearman <- min(
  slide_robustness$spearman_rho,
  na.rm = TRUE
)

minimum_slide_direction <- min(
  slide_robustness$direction_concordance,
  na.rm = TRUE
)

expected_sc_samples <- c(
  GSE250130 = 28L,
  GSE183837 = 9L
)

observed_sc_samples <- table(
  sc_qc_summary$dataset_id
)

sc_sample_pass <-
  all(
    as.integer(
      observed_sc_samples[
        names(expected_sc_samples)
      ]
    ) ==
      as.integer(
        expected_sc_samples
      )
  )

sc_retention_pass <-
  all(
    sc_qc_summary$retained_cell_n > 0L
  )

bulk_matrix_pass <-
  all(
    bulk_matrix_qc$gene_n >=
      5000L
  )

candidate_coverage_pass <-
  minimum_candidate_coverage >=
    200L

gate_status <- if (
  bulk_matrix_pass &&
    candidate_coverage_pass &&
    sc_sample_pass &&
    sc_retention_pass
) {
  "READY_FOR_SINGLE_CELL_INTEGRATION_AND_FINAL_MULTILAYER_CANDIDATE_PRIORITIZATION"
} else {
  "HOLD_STEP04D_QC_REVIEW_REQUIRED"
}

gate_table <- data.frame(
  gate_component = c(
    "STEP04C gate",
    "Bulk gene matrices",
    "Occupancy-supported candidate coverage",
    "Bulk gene-level contrasts",
    "Program-score models",
    "GSE165004 leave-one-slide candidate robustness",
    "GSE250130 donor/sample pre-QC",
    "GSE183837 donor/sample pre-QC",
    "Single-cell retained cells",
    "GSE111974 independence interpretation",
    "Definitive direct-target assignment",
    "Scientific figures",
    "STEP04D gate"
  ),
  status = c(
    "PASS",
    ifelse(
      bulk_matrix_pass,
      "PASS_3_OF_3",
      "FAIL"
    ),
    paste0(
      "MIN_MEASURABLE_",
      minimum_candidate_coverage,
      "_OF_284"
    ),
    "COMPLETED_4_CONTRASTS",
    "COMPLETED_4_PROGRAMS_ACROSS_4_CONTRASTS",
    paste0(
      "MIN_SPEARMAN_",
      format(
        minimum_slide_spearman,
        digits = 6
      ),
      "_MIN_DIRECTION_",
      format(
        minimum_slide_direction,
        digits = 6
      )
    ),
    paste0(
      "PASS_",
      sum(
        sc_qc_summary$dataset_id ==
          "GSE250130"
      ),
      "_OF_28"
    ),
    paste0(
      "PASS_",
      sum(
        sc_qc_summary$dataset_id ==
          "GSE183837"
      ),
      "_OF_9"
    ),
    ifelse(
      sc_retention_pass,
      "PASS_ALL_SAMPLES_RETAIN_CELLS",
      "FAIL"
    ),
    "SHARED_CONTROL_CONTEXT_EXTENSION_ONLY",
    "NOT_PERFORMED",
    "NONE_CREATED",
    gate_status
  ),
  required_action = c(
    "Retain locked designs and source references.",
    "Use only locked gene-level matrices.",
    "Review coverage before interpreting absent genes.",
    "Interpret cross-context direction as supportive or mixed, not causal.",
    "Do not interpret scores as diagnostic classifiers.",
    "Retain slide sensitivity in the primary patient result.",
    "Use donor-aware integration and pseudobulk.",
    "Use donor-aware integration and pseudobulk.",
    "Review sample-level QC before clustering.",
    "Do not count as independent external validation.",
    "Use occupancy-supported candidate terminology.",
    "Render figures only after single-cell integration and candidate prioritization.",
    "Proceed only if the Gate is READY."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    dirs$results,
    "STEP04D_GATE_TABLE.csv"
  )
)

write_text(
  c(
    "# STEP04D bulk projection and single-cell pre-QC card",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Bulk matrix version: ",
      bulk_matrix_version
    ),
    paste0(
      "- Bulk projection version: ",
      bulk_projection_version
    ),
    paste0(
      "- Single-cell pre-QC version: ",
      sc_preQC_version
    ),
    paste0(
      "- Bulk genes per dataset: ",
      paste(
        paste0(
          bulk_matrix_qc$dataset_id,
          "=",
          bulk_matrix_qc$gene_n
        ),
        collapse = "; "
      )
    ),
    paste0(
      "- Minimum measurable occupancy-supported candidates: ",
      minimum_candidate_coverage,
      "/284"
    ),
    paste0(
      "- Minimum GSE165004 leave-one-slide Spearman: ",
      format(
        minimum_slide_spearman,
        digits = 6
      )
    ),
    paste0(
      "- Minimum leave-one-slide direction concordance: ",
      format(
        minimum_slide_direction,
        digits = 6
      )
    ),
    paste0(
      "- GSE250130 raw/retained cells: ",
      sc_dataset_summary$raw_cell_n[
        sc_dataset_summary$dataset_id ==
          "GSE250130"
      ],
      "/",
      sc_dataset_summary$retained_cell_n[
        sc_dataset_summary$dataset_id ==
          "GSE250130"
      ]
    ),
    paste0(
      "- GSE183837 raw/retained cells: ",
      sc_dataset_summary$raw_cell_n[
        sc_dataset_summary$dataset_id ==
          "GSE183837"
      ],
      "/",
      sc_dataset_summary$retained_cell_n[
        sc_dataset_summary$dataset_id ==
          "GSE183837"
      ]
    ),
    "- GSE111974 independent external validation: NO",
    "- Definitive direct-target assignments: 0",
    "- Scientific figures: NONE",
    paste0(
      "- Gate: ",
      gate_status
    )
  ),
  file.path(
    dirs$results,
    "S4D_result_card.md"
  )
)

write_text(
  c(
    "# STEP04D Methods draft",
    "",
    "Processed bulk-expression matrices were mapped to gene symbols using locked GEO platform tables. Only probes with one unambiguous gene symbol were retained, and duplicate probes were averaged with limma. A log2 transformation was applied only when matrix-scale diagnostics indicated that the submitted values were not already log transformed.",
    "",
    "The prespecified 436-gene RNA interaction program and its 284-gene occupancy-supported subset were projected into the three bulk cohorts using signed standardized expression scores based on the direction of the locked RNA interaction effect. GSE165004 was analyzed separately for RPL versus control and UIF versus control, with leave-one-slide sensitivity. GSE92324 was adjusted for HCG collection day. GSE111974 was treated as a shared-control RIF context extension.",
    "",
    "For GSE250130 and GSE183837, each donor/sample matrix was read independently. Cells with fewer than 500 counts, fewer than 200 detected features, excessive detected features under a sample-adaptive upper threshold, or more than 25% mitochondrial counts were removed. Filtered sparse matrices were retained as donor-level preprocessing objects; cells were not treated as independent biological replicates."
  ),
  file.path(
    dirs$text,
    "STEP04D_Methods_draft.md"
  )
)

write_text(
  c(
    "# STEP04D Results shell",
    "",
    "This file is intentionally a result shell. Numeric bulk projection, leave-one-slide and single-cell preprocessing results must be inserted only after the STEP04D Gate and QC tables are reviewed.",
    "",
    "No diagnostic, causal or definitive direct-target claim is authorized."
  ),
  file.path(
    dirs$text,
    "STEP04D_Results_shell.md"
  )
)

parameters <- list(
  schema_version = "1.0",
  run_id = run_id,
  bulk_matrix_version =
    bulk_matrix_version,
  bulk_projection_version =
    bulk_projection_version,
  sc_preQC_version =
    sc_preQC_version,
  frozen_candidate_n = 436,
  occupancy_supported_n = 284,
  minimum_measurable_occupancy_supported_n =
    minimum_candidate_coverage,
  minimum_leave_one_slide_spearman =
    minimum_slide_spearman,
  minimum_leave_one_slide_direction_concordance =
    minimum_slide_direction,
  GSE250130_sample_n =
    sum(
      sc_qc_summary$dataset_id ==
        "GSE250130"
    ),
  GSE183837_sample_n =
    sum(
      sc_qc_summary$dataset_id ==
        "GSE183837"
    ),
  definitive_direct_target_assignment =
    FALSE,
  scientific_figures_created =
    FALSE,
  gate = gate_status
)

jsonlite::write_json(
  parameters,
  file.path(
    dirs$parameters,
    "STEP04D_parameters.json"
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
    "sessionInfo_POST_RUN_STEP04D.txt"
  )
)

source_file <- tryCatch(
  sys.frame(1)$ofile,
  error = function(error_object) NULL
)

if (
  !is.null(source_file) &&
    file.exists(source_file)
) {
  copy_if_different(
    source_file,
    file.path(
      dirs$code,
      "51_step04D_v6_bulk_projection_and_sc_preQC.R"
    )
  )
}

# ----------------------------------------------------------------------
# Archive derived summaries and manifests only
# ----------------------------------------------------------------------

matrix_paths <- list.files(
  bulk_matrix_root,
  pattern = "\\.rds$",
  full.names = TRUE
)

result_paths <- c(
  file.path(
    dirs$results,
    "STEP04D_GATE_TABLE.csv"
  ),
  file.path(
    dirs$results,
    "S4D_result_card.md"
  ),
  file.path(
    dirs$results,
    "S4D_bulk_matrix_qc.csv"
  ),
  file.path(
    dirs$results,
    "STEP04D_v3_GSE92324_SUPPLEMENTARY_PROFILE_SOURCE.csv"
  ),
  file.path(
    dirs$results,
    "S4D_sc_sample_qc.csv"
  ),
  file.path(
    dirs$results,
    "S4D_sc_dataset_qc.csv"
  ),
  file.path(
    bulk_result_root,
    "score_models.csv"
  ),
  file.path(
    bulk_result_root,
    "gene_coverage.csv"
  ),
  file.path(
    bulk_result_root,
    "loso_robustness.csv"
  ),
  file.path(
    bulk_result_root,
    "cross_context_candidates.csv"
  ),
  file.path(
    dirs$text,
    "STEP04D_Methods_draft.md"
  ),
  file.path(
    dirs$text,
    "STEP04D_Results_shell.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04D_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04D.txt"
  ),
  file.path(
    dirs$code,
    "51_step04D_v6_bulk_projection_and_sc_preQC.R"
  ),
  log_path
)

artifact_manifest_paths <- c(
  matrix_paths,
  result_paths[
    file.exists(
      result_paths
    )
  ]
)

artifact_manifest <- data.frame(
  path =
    normalize_slashes(
      artifact_manifest_paths
    ),
  size_bytes = as.numeric(
    file.info(
      artifact_manifest_paths
    )$size
  ),
  md5 = vapply(
    artifact_manifest_paths,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    artifact_manifest_paths,
    sha256_file,
    character(1)
  ),
  artifact_status =
    "CANDIDATE_PENDING_MANUAL_QC",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  artifact_manifest,
  file.path(
    dirs$results,
    "S4D_artifact_manifest.csv"
  )
)

archive_artifacts <- c(
  result_paths,
  file.path(
    dirs$results,
    "S4D_artifact_manifest.csv"
  )
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
      "Methods|Results",
      basename(path_value),
      ignore.case = TRUE
    )
  ) {
    "txt"
  } else if (
    grepl(
      "QC|CARD|MANIFEST|ROBUSTNESS|COVERAGE|GATE|SUMMARY",
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
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP04D_SHORT_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04D archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    "- Bulk gene matrices: 3",
    "- Bulk contrasts: 4",
    "- Program-score definitions: 4",
    "- Single-cell sample objects: 37",
    "- Large single-cell RDS objects copied into archive: NO",
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
      "- Accepted run candidate: ",
      run_id
    ),
    paste0(
      "- Bulk matrices: ",
      bulk_matrix_version,
      " CANDIDATE"
    ),
    paste0(
      "- Bulk projection: ",
      bulk_projection_version,
      " CANDIDATE"
    ),
    paste0(
      "- Single-cell pre-QC: ",
      sc_preQC_version,
      " CANDIDATE"
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
cat("STEP04D v6 bulk projection and single-cell pre-QC completed.\n")
cat(
  "Bulk matrix version: ",
  bulk_matrix_version,
  "\n"
)
cat(
  "Bulk projection version: ",
  bulk_projection_version,
  "\n"
)
cat(
  "Single-cell pre-QC version: ",
  sc_preQC_version,
  "\n"
)
cat(
  "Bulk gene counts: ",
  paste(
    paste0(
      bulk_matrix_qc$dataset_id,
      "=",
      bulk_matrix_qc$gene_n
    ),
    collapse = "; "
  ),
  "\n"
)
cat(
  "Minimum measurable occupancy-supported genes: ",
  minimum_candidate_coverage,
  "/284\n"
)
cat(
  "Minimum GSE165004 leave-one-slide Spearman: ",
  format(
    minimum_slide_spearman,
    digits = 6
  ),
  "\n"
)
cat(
  "Minimum leave-one-slide direction concordance: ",
  format(
    minimum_slide_direction,
    digits = 6
  ),
  "\n"
)
cat(
  "GSE250130 samples processed: ",
  sum(
    sc_qc_summary$dataset_id ==
      "GSE250130"
  ),
  "/28\n"
)
cat(
  "GSE183837 samples processed: ",
  sum(
    sc_qc_summary$dataset_id ==
      "GSE183837"
  ),
  "/9\n"
)
cat("Definitive direct-target assignments: 0\n")
cat("Scientific figures: NONE\n")
cat("Gate: ", gate_status, "\n")
cat("Authoritative archive:\n  ", archive_root, "\n\n")
cat("Upload:\n")
cat("  06_results/STEP04D_GATE_TABLE.csv\n")
cat("  06_results/S4D_bulk_matrix_qc.csv\n")
cat("  06_results/S4D_sc_sample_qc.csv\n")
cat("  06_results/S4D_sc_dataset_qc.csv\n")
cat("  06_results/S4D_bulk/BVP_v1_2/score_models.csv
")
cat("  06_results/S4D_bulk/BVP_v1_2/loso_robustness.csv
")
cat("  06_results/S4D_bulk/BVP_v1_2/cross_context_candidates.csv
")
cat("============================================================\n")

# ==============================================================================
# Analysis stage 5: 52_step04D_QC1_v2_repair_and_lock_outputs
# ==============================================================================
# ============================================================
# STEP04D-QC1 - Repair and lock STEP04D bulk projection outputs
# Nutrition_Endometrial_VDR_v2
#
# Purpose:
# 1. recompute the four program-score models with preserved names;
# 2. apply BH correction across the four programs per contrast;
# 3. rebuild the cross-context candidate table without NA__NA;
# 4. assign an explicit label to candidates not measured anywhere;
# 5. create a donor/sample single-cell QC watchlist;
# 6. lock the corrected STEP04D handoff without rerunning bulk
#    matrices or the 37 single-cell pre-QC objects.
#
# No scientific figures are generated.
#
# Run from the project root:
# source(
#   "04_code/R/52_step04D_QC1_v2_repair_and_lock_outputs.R",
#   encoding = "UTF-8",
#   echo = FALSE
# )
# ============================================================

options(stringsAsFactors = FALSE)
options(warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

step04_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP04_combined_temporal_patient_single_cell"
)

run_id <- paste0(
  "STEP04D_QC1_",
  format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )
)

repair_version <-
  "STEP04D_QC_REPAIR_v1.1"

projection_version <-
  "BULK_VDR_PROGRAM_PROJECTION_v1.3"

required_packages <- c(
  "digest",
  "jsonlite",
  "limma"
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
  "S4D_qc_lock",
  "BVP_v1_3"
)

archive_root <- file.path(
  dirs$archive_parent,
  paste0(
    "S4DQC_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    )
  )
)

for (directory_path in c(
  output_root,
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
    dirname(
      destination_path
    ),
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

log_event(
  "Starting",
  run_id
)

# ----------------------------------------------------------------------
# Locate and validate STEP04D v6 inputs
# ----------------------------------------------------------------------

v6_projection_root <- file.path(
  dirs$results,
  "S4D_bulk",
  "BVP_v1_2"
)

v6_gate_path <- file.path(
  dirs$results,
  "STEP04D_GATE_TABLE.csv"
)

v6_manifest_path <- file.path(
  dirs$results,
  "S4D_artifact_manifest.csv"
)

sample_scores_path <- file.path(
  v6_projection_root,
  "sample_scores.csv"
)

candidate_gene_results_path <- file.path(
  v6_projection_root,
  "candidate_gene_results.csv"
)

sc_sample_qc_path <- file.path(
  dirs$results,
  "S4D_sc_sample_qc.csv"
)

sc_dataset_qc_path <- file.path(
  dirs$results,
  "S4D_sc_dataset_qc.csv"
)

locked_design_root <- file.path(
  dirs$results,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0"
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

required_paths <- c(
  v6_gate_path,
  v6_manifest_path,
  sample_scores_path,
  candidate_gene_results_path,
  sc_sample_qc_path,
  sc_dataset_qc_path,
  candidate_path,
  file.path(
    locked_design_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE165004.csv"
  ),
  file.path(
    locked_design_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE92324.csv"
  ),
  file.path(
    locked_design_root,
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE111974.csv"
  )
)

missing_paths <- required_paths[
  !file.exists(
    required_paths
  )
]

if (length(missing_paths) > 0L) {
  stop(
    "Required STEP04D-QC1 input is missing:\n",
    paste(
      missing_paths,
      collapse = "\n"
    )
  )
}

v6_gate <- read.csv(
  v6_gate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

v6_gate_status <- v6_gate$status[
  v6_gate$gate_component ==
    "STEP04D gate"
]

if (
  length(v6_gate_status) != 1L ||
    v6_gate_status !=
      "READY_FOR_SINGLE_CELL_INTEGRATION_AND_FINAL_MULTILAYER_CANDIDATE_PRIORITIZATION"
) {
  stop(
    "STEP04D v6 did not finish with the expected READY Gate."
  )
}

v6_manifest <- read.csv(
  v6_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

repair_input_paths <- c(
  sample_scores_path,
  candidate_gene_results_path,
  sc_sample_qc_path,
  sc_dataset_qc_path
)

repair_input_names <- c(
  "sample_scores",
  "candidate_gene_results",
  "sc_sample_qc",
  "sc_dataset_qc"
)

repair_input_normalized <- normalize_slashes(
  repair_input_paths
)

manifest_path_normalized <- normalize_slashes(
  v6_manifest$path
)

manifest_match_index <- match(
  repair_input_normalized,
  manifest_path_normalized
)

input_identity_audit <- data.frame(
  input_name =
    repair_input_names,
  path =
    repair_input_normalized,
  file_exists =
    file.exists(
      repair_input_paths
    ),
  size_bytes =
    as.numeric(
      file.info(
        repair_input_paths
      )$size
    ),
  listed_in_v6_manifest =
    !is.na(
      manifest_match_index
    ),
  manifest_md5 =
    NA_character_,
  observed_md5 =
    vapply(
      repair_input_paths,
      md5_file,
      character(1)
    ),
  manifest_sha256 =
    NA_character_,
  observed_sha256 =
    vapply(
      repair_input_paths,
      sha256_file,
      character(1)
    ),
  identity_status =
    NA_character_,
  stringsAsFactors = FALSE
)

listed_rows <- which(
  input_identity_audit$listed_in_v6_manifest
)

if (length(listed_rows) > 0L) {
  matched_manifest_rows <- manifest_match_index[
    listed_rows
  ]

  input_identity_audit$manifest_md5[
    listed_rows
  ] <- as.character(
    v6_manifest$md5[
      matched_manifest_rows
    ]
  )

  input_identity_audit$manifest_sha256[
    listed_rows
  ] <- as.character(
    v6_manifest$sha256[
      matched_manifest_rows
    ]
  )
}

input_identity_audit$identity_status <- ifelse(
  !input_identity_audit$file_exists,
  "FAIL_FILE_MISSING",
  ifelse(
    input_identity_audit$size_bytes <= 0 |
      is.na(
        input_identity_audit$size_bytes
      ),
    "FAIL_EMPTY_FILE",
    ifelse(
      input_identity_audit$listed_in_v6_manifest &
        input_identity_audit$manifest_md5 ==
          input_identity_audit$observed_md5 &
        input_identity_audit$manifest_sha256 ==
          input_identity_audit$observed_sha256,
      "PASS_V6_MANIFEST_HASH_MATCH",
      ifelse(
        input_identity_audit$listed_in_v6_manifest,
        "FAIL_V6_MANIFEST_HASH_MISMATCH",
        "PASS_DIRECT_HASH_CAPTURE_V6_MANIFEST_OMISSION"
      )
    )
  )
)

write_csv_utf8(
  input_identity_audit,
  file.path(
    output_root,
    "STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv"
  )
)

if (any(
  grepl(
    "^FAIL_",
    input_identity_audit$identity_status
  )
)) {
  stop(
    "At least one STEP04D-QC1 repair input failed identity validation. Review STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv."
  )
}

manifest_omission_n <- sum(
  !input_identity_audit$listed_in_v6_manifest
)

log_event(
  "Repair input identity audit passed.",
  "Inputs listed in v6 manifest:",
  sum(
    input_identity_audit$listed_in_v6_manifest
  ),
  "/4;",
  "direct hash captures for manifest omissions:",
  manifest_omission_n
)

# ----------------------------------------------------------------------
# Read locked designs and source tables
# ----------------------------------------------------------------------

read_locked_design <- function(
    dataset_id
) {
  read.csv(
    file.path(
      locked_design_root,
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

design_165004 <- read_locked_design(
  "GSE165004"
)

design_92324 <- read_locked_design(
  "GSE92324"
)

design_111974 <- read_locked_design(
  "GSE111974"
)

sample_scores <- read.csv(
  sample_scores_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_gene_results <- read.csv(
  candidate_gene_results_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

candidate_table <- read.csv(
  candidate_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sc_sample_qc <- read.csv(
  sc_sample_qc_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sc_dataset_qc <- read.csv(
  sc_dataset_qc_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_programs <- c(
  "RNA_ALL_436",
  "OCCUPANCY_SUPPORTED_284",
  "SHARED_OCCUPANCY_46",
  "PROMOTER_OCCUPANCY_84"
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
    "The locked 436-candidate table failed validation."
  )
}

if (!setequal(
  unique(
    sample_scores$program_name
  ),
  expected_programs
)) {
  stop(
    "The sample-score file does not contain the four expected program names."
  )
}

if (
  nrow(sc_sample_qc) != 37L ||
    sum(
      sc_sample_qc$dataset_id ==
        "GSE250130"
    ) !=
      28L ||
    sum(
      sc_sample_qc$dataset_id ==
        "GSE183837"
    ) !=
      9L
) {
  stop(
    "The 37-sample single-cell QC table failed validation."
  )
}

# ----------------------------------------------------------------------
# Recompute score models with preserved row names
# ----------------------------------------------------------------------

build_score_matrix <- function(
    dataset_id,
    sample_design
) {
  score_subset <- sample_scores[
    sample_scores$dataset_id ==
      dataset_id,
    ,
    drop = FALSE
  ]

  score_matrix <- matrix(
    NA_real_,
    nrow =
      length(
        expected_programs
      ),
    ncol =
      nrow(
        sample_design
      ),
    dimnames = list(
      expected_programs,
      sample_design$sample_accession
    )
  )

  for (program_name in expected_programs) {
    program_rows <- score_subset[
      score_subset$program_name ==
        program_name,
      ,
      drop = FALSE
    ]

    program_rows <- program_rows[
      match(
        sample_design$sample_accession,
        program_rows$sample_accession
      ),
      ,
      drop = FALSE
    ]

    if (
      nrow(program_rows) !=
        nrow(sample_design) ||
        any(
          is.na(
            program_rows$score
          )
        )
    ) {
      stop(
        dataset_id,
        " sample-score alignment failed for ",
        program_name
      )
    }

    score_matrix[
      program_name,
    ] <- program_rows$score
  }

  score_matrix
}

fit_score_matrix <- function(
    score_matrix,
    design_matrix,
    contrast_matrix,
    dataset_id
) {
  if (
    ncol(score_matrix) !=
      nrow(design_matrix)
  ) {
    stop(
      dataset_id,
      " score/design dimensions do not match."
    )
  }

  fit <- limma::lmFit(
    score_matrix,
    design_matrix
  )

  fit <- limma::contrasts.fit(
    fit,
    contrast_matrix
  )

  fit <- limma::eBayes(
    fit,
    trend = FALSE,
    robust = FALSE
  )

  result_rows <- list()

  for (contrast_name in colnames(
    contrast_matrix
  )) {
    table_value <- limma::topTable(
      fit,
      coef = contrast_name,
      number = Inf,
      sort.by = "none",
      adjust.method = "BH"
    )

    table_value$program_name <-
      rownames(
        table_value
      )

    table_value$dataset_id <-
      dataset_id

    table_value$contrast <-
      contrast_name

    table_value$score_model <-
      "FOUR_PRESPECIFIED_PROGRAMS_MODELED_TOGETHER"

    result_rows[[
      contrast_name
    ]] <- table_value
  }

  do.call(
    rbind,
    result_rows
  )
}

score_matrix_165004 <- build_score_matrix(
  "GSE165004",
  design_165004
)

group_165004 <- factor(
  design_165004$group,
  levels = c(
    "CONTROL",
    "RPL",
    "UIF"
  )
)

design_matrix_165004 <- model.matrix(
  ~ 0 + group_165004
)

colnames(
  design_matrix_165004
) <- levels(
  group_165004
)

contrast_matrix_165004 <- limma::makeContrasts(
  RPL_vs_CONTROL =
    RPL - CONTROL,
  UIF_vs_CONTROL =
    UIF - CONTROL,
  levels =
    design_matrix_165004
)

score_result_165004 <- fit_score_matrix(
  score_matrix_165004,
  design_matrix_165004,
  contrast_matrix_165004,
  "GSE165004"
)

score_matrix_92324 <- build_score_matrix(
  "GSE92324",
  design_92324
)

group_92324 <- factor(
  design_92324$group,
  levels = c(
    "CONTROL",
    "IMPLANTATION_FAILURE"
  )
)

day_92324 <- factor(
  design_92324$collection_day
)

design_matrix_92324 <- model.matrix(
  ~ 0 + group_92324 + day_92324
)

colnames(
  design_matrix_92324
) <- make.names(
  colnames(
    design_matrix_92324
  )
)

contrast_matrix_92324 <- limma::makeContrasts(
  IF_vs_CONTROL =
    group_92324IMPLANTATION_FAILURE -
    group_92324CONTROL,
  levels =
    design_matrix_92324
)

score_result_92324 <- fit_score_matrix(
  score_matrix_92324,
  design_matrix_92324,
  contrast_matrix_92324,
  "GSE92324"
)

score_matrix_111974 <- build_score_matrix(
  "GSE111974",
  design_111974
)

group_111974 <- factor(
  design_111974$group,
  levels = c(
    "CONTROL",
    "RIF"
  )
)

design_matrix_111974 <- model.matrix(
  ~ 0 + group_111974
)

colnames(
  design_matrix_111974
) <- levels(
  group_111974
)

contrast_matrix_111974 <- limma::makeContrasts(
  RIF_vs_CONTROL =
    RIF - CONTROL,
  levels =
    design_matrix_111974
)

score_result_111974 <- fit_score_matrix(
  score_matrix_111974,
  design_matrix_111974,
  contrast_matrix_111974,
  "GSE111974"
)

score_models_corrected <- rbind(
  score_result_165004,
  score_result_92324,
  score_result_111974
)

score_models_corrected <- score_models_corrected[
  ,
  c(
    "dataset_id",
    "contrast",
    "program_name",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "score_model"
  ),
  drop = FALSE
]

if (
  nrow(score_models_corrected) != 16L ||
    length(
      unique(
        score_models_corrected$program_name
      )
    ) !=
      4L ||
    any(
      score_models_corrected$program_name ==
        "1"
    )
) {
  stop(
    "Corrected program-score models failed the 16-row/4-program validation."
  )
}

write_csv_utf8(
  score_models_corrected,
  file.path(
    output_root,
    "score_models_v1_3.csv"
  )
)

# ----------------------------------------------------------------------
# Rebuild clean cross-context candidate table
# ----------------------------------------------------------------------

candidate_base_columns <- c(
  "feature_id",
  "program_class",
  "interaction_logFC",
  "interaction_FDR",
  "occupancy_supported",
  "shared_occupancy_supported",
  "any_promoter_overlap_support"
)

missing_candidate_columns <- setdiff(
  candidate_base_columns,
  names(
    candidate_table
  )
)

if (length(missing_candidate_columns) > 0L) {
  stop(
    "Locked candidate table lacks required column(s): ",
    paste(
      missing_candidate_columns,
      collapse = ", "
    )
  )
}

candidate_convergence <- candidate_table[
  ,
  candidate_base_columns,
  drop = FALSE
]

expected_contexts <- data.frame(
  dataset_id = c(
    "GSE165004",
    "GSE165004",
    "GSE92324",
    "GSE111974"
  ),
  contrast = c(
    "UIF_vs_CONTROL",
    "RPL_vs_CONTROL",
    "IF_vs_CONTROL",
    "RIF_vs_CONTROL"
  ),
  context_id = c(
    "GSE165004__UIF_vs_CONTROL",
    "GSE165004__RPL_vs_CONTROL",
    "GSE92324__IF_vs_CONTROL",
    "GSE111974__RIF_vs_CONTROL"
  ),
  stringsAsFactors = FALSE
)

for (context_index in seq_len(
  nrow(
    expected_contexts
  )
)) {
  dataset_id <- expected_contexts$dataset_id[
    context_index
  ]

  contrast_name <- expected_contexts$contrast[
    context_index
  ]

  context_id <- expected_contexts$context_id[
    context_index
  ]

  context_rows <- candidate_gene_results[
    !is.na(
      candidate_gene_results$dataset_id
    ) &
      !is.na(
        candidate_gene_results$contrast
      ) &
      candidate_gene_results$dataset_id ==
        dataset_id &
      candidate_gene_results$contrast ==
        contrast_name,
    c(
      "feature_id",
      "logFC",
      "P.Value",
      "adj.P.Val"
    ),
    drop = FALSE
  ]

  if (anyDuplicated(
    context_rows$feature_id
  )) {
    stop(
      "Duplicate candidate/context rows detected for ",
      context_id
    )
  }

  names(context_rows)[
    -1
  ] <- paste0(
    context_id,
    "__",
    c(
      "logFC",
      "P.Value",
      "adj.P.Val"
    )
  )

  candidate_convergence <- merge(
    candidate_convergence,
    context_rows,
    by = "feature_id",
    all.x = TRUE,
    sort = FALSE
  )
}

logFC_columns <- paste0(
  expected_contexts$context_id,
  "__logFC"
)

adj_columns <- paste0(
  expected_contexts$context_id,
  "__adj.P.Val"
)

logFC_matrix <- as.matrix(
  candidate_convergence[
    ,
    logFC_columns,
    drop = FALSE
  ]
)

adj_matrix <- as.matrix(
  candidate_convergence[
    ,
    adj_columns,
    drop = FALSE
  ]
)

candidate_convergence$measured_context_n <-
  rowSums(
    is.finite(
      logFC_matrix
    )
  )

candidate_convergence$positive_context_n <-
  rowSums(
    logFC_matrix > 0,
    na.rm = TRUE
  )

candidate_convergence$negative_context_n <-
  rowSums(
    logFC_matrix < 0,
    na.rm = TRUE
  )

candidate_convergence$FDR05_context_n <-
  rowSums(
    adj_matrix < 0.05,
    na.rm = TRUE
  )

candidate_convergence$cross_context_direction_label <-
  ifelse(
    candidate_convergence$measured_context_n ==
      0L,
    "NOT_MEASURED_IN_ANY_BULK_CONTEXT",
    ifelse(
      candidate_convergence$positive_context_n ==
        candidate_convergence$measured_context_n,
      "POSITIVE_ALL_MEASURED_CONTEXTS",
      ifelse(
        candidate_convergence$negative_context_n ==
          candidate_convergence$measured_context_n,
        "NEGATIVE_ALL_MEASURED_CONTEXTS",
        "CONTEXT_DEPENDENT_OR_MIXED_DIRECTION"
      )
    )
  )

if (
  nrow(candidate_convergence) != 436L ||
    any(
      grepl(
        "^NA__NA",
        names(
          candidate_convergence
        )
      )
    ) ||
    any(
      candidate_convergence$measured_context_n ==
        0L &
        candidate_convergence$cross_context_direction_label !=
          "NOT_MEASURED_IN_ANY_BULK_CONTEXT"
    )
) {
  stop(
    "Corrected cross-context candidate table failed semantic validation."
  )
}

write_csv_utf8(
  candidate_convergence,
  file.path(
    output_root,
    "cross_context_candidates_v1_3.csv"
  )
)

context_coverage <- data.frame(
  measured_context_n =
    0:4,
  candidate_n =
    as.integer(
      table(
        factor(
          candidate_convergence$measured_context_n,
          levels = 0:4
        )
      )
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  context_coverage,
  file.path(
    output_root,
    "candidate_context_coverage.csv"
  )
)

# ----------------------------------------------------------------------
# Create single-cell sample watchlist
# ----------------------------------------------------------------------

sc_sample_qc$QC_review_tier <- ifelse(
  sc_sample_qc$retained_fraction <
    0.65 |
    sc_sample_qc$median_percent_mito >
      20,
  "REVIEW_REQUIRED",
  ifelse(
    sc_sample_qc$retained_fraction <
      0.75 |
      sc_sample_qc$median_percent_mito >
        15,
    "CAUTION",
    "PASS"
  )
)

sc_sample_qc$integration_inclusion <-
  ifelse(
    sc_sample_qc$QC_review_tier ==
      "PASS",
    "INCLUDE",
    "INCLUDE_WITH_DONOR_LEVEL_QC_REVIEW"
  )

sc_watchlist <- sc_sample_qc[
  sc_sample_qc$QC_review_tier !=
    "PASS",
  ,
  drop = FALSE
]

write_csv_utf8(
  sc_sample_qc,
  file.path(
    output_root,
    "sc_sample_qc_reviewed.csv"
  )
)

write_csv_utf8(
  sc_watchlist,
  file.path(
    output_root,
    "sc_sample_watchlist.csv"
  )
)

watchlist_summary <- aggregate(
  sample_accession ~
    dataset_id +
    QC_review_tier,
  data =
    sc_sample_qc,
  FUN =
    length
)

names(
  watchlist_summary
)[
  names(
    watchlist_summary
  ) ==
    "sample_accession"
] <- "sample_n"

write_csv_utf8(
  watchlist_summary,
  file.path(
    output_root,
    "sc_watchlist_summary.csv"
  )
)

# ----------------------------------------------------------------------
# QC summary and lock decision
# ----------------------------------------------------------------------

score_significant_summary <- aggregate(
  adj.P.Val ~
    dataset_id +
    contrast,
  data =
    transform(
      score_models_corrected,
      adj.P.Val =
        adj.P.Val < 0.05
    ),
  FUN =
    sum
)

names(
  score_significant_summary
)[
  names(
    score_significant_summary
  ) ==
    "adj.P.Val"
] <- "FDR05_program_n"

direction_counts <- as.data.frame(
  table(
    candidate_convergence$cross_context_direction_label
  ),
  stringsAsFactors = FALSE
)

names(
  direction_counts
) <- c(
  "direction_label",
  "candidate_n"
)

write_csv_utf8(
  score_significant_summary,
  file.path(
    output_root,
    "score_significance_summary.csv"
  )
)

write_csv_utf8(
  direction_counts,
  file.path(
    output_root,
    "candidate_direction_summary.csv"
  )
)

review_required_n <- sum(
  sc_sample_qc$QC_review_tier ==
    "REVIEW_REQUIRED"
)

caution_n <- sum(
  sc_sample_qc$QC_review_tier ==
    "CAUTION"
)

gate_status <-
  "READY_FOR_SINGLE_CELL_INTEGRATION_WITH_GSE183837_SAMPLE_WATCHLIST"

gate_table <- data.frame(
  gate_component = c(
    "STEP04D v6 computational completion",
    "STEP04D v6 repair-input identity audit",
    "Program-score name integrity",
    "Program-score multiplicity correction",
    "Cross-context table NA__NA removal",
    "Unmeasured-candidate direction semantics",
    "Candidate row integrity",
    "Single-cell sample count",
    "Single-cell review-required samples",
    "Single-cell caution samples",
    "Automatic sample exclusion",
    "Definitive direct-target assignment",
    "Scientific figures",
    "Corrected bulk projection version",
    "STEP04D-QC1 lock",
    "STEP04D-QC1 gate"
  ),
  status = c(
    "PASS",
    paste0(
      "PASS_",
      4L - manifest_omission_n,
      "_MANIFEST_HASH_MATCH_",
      manifest_omission_n,
      "_DIRECT_HASH_CAPTURE"
    ),
    "PASS_4_NAMED_PROGRAMS",
    "BH_WITHIN_EACH_CONTRAST_ACROSS_4_PROGRAMS",
    "PASS_NO_NA_NA_COLUMNS",
    "PASS_EXPLICIT_NOT_MEASURED_LABEL",
    "PASS_436_UNIQUE_CANDIDATES",
    "PASS_37_OF_37",
    paste0(
      "WATCHLIST_",
      review_required_n
    ),
    paste0(
      "CAUTION_",
      caution_n
    ),
    "NONE",
    "NOT_PERFORMED",
    "NONE_CREATED",
    paste0(
      projection_version,
      "_LOCKED"
    ),
    paste0(
      repair_version,
      "_LOCKED"
    ),
    gate_status
  ),
  required_action = c(
    "Retain the v6 bulk matrices and single-cell pre-QC objects.",
    "Retain STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv; direct hashes repair the v6 manifest omission without altering source files.",
    "Use only the corrected v1.3 score table.",
    "Interpret FDR within each prespecified four-program contrast family.",
    "Use only the corrected v1.3 candidate table.",
    "Do not classify unmeasured candidates as positive or negative.",
    "Preserve the locked 436-candidate universe.",
    "Proceed donor-aware.",
    "Inspect these donors after normalization, doublet scoring and cell-type annotation.",
    "Retain and inspect these donors during integration.",
    "Do not exclude donors solely from this pre-QC watchlist.",
    "Continue occupancy-supported candidate terminology.",
    "Generate figures only after integrated biological QC.",
    "Use v1.3 for downstream prioritization.",
    "Retain repair manifest and archive.",
    "Proceed to single-cell integration with the watchlist active."
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_table,
  file.path(
    output_root,
    "STEP04D_QC1_GATE.csv"
  )
)

result_card <- c(
  "# STEP04D-QC1 corrected handoff card",
  "",
  paste0(
    "- Run ID: ",
    run_id
  ),
  paste0(
    "- Repair version: ",
    repair_version,
    " LOCKED"
  ),
  paste0(
    "- Corrected bulk projection: ",
    projection_version,
    " LOCKED"
  ),
  "- Bulk matrices retained: BULK_GENE_MATRIX_v1.2",
  "- Single-cell pre-QC objects retained: SINGLE_CELL_PREQC_OBJECTS_v1.0",
  paste0(
    "- Repair inputs listed in v6 manifest/direct-hash captured: ",
    4L - manifest_omission_n,
    "/",
    manifest_omission_n
  ),
  "- Program-score rows: 16",
  "- Named programs: 4",
  "- Multiple-testing rule: BH across four prespecified programs within each contrast",
  "- Cross-context candidates: 436",
  paste0(
    "- Candidates measured in 0/1/2/3/4 contexts: ",
    paste(
      context_coverage$candidate_n,
      collapse = "/"
    )
  ),
  paste0(
    "- Single-cell review-required samples: ",
    review_required_n
  ),
  paste0(
    "- Single-cell caution samples: ",
    caution_n
  ),
  "- Automatic donor exclusions: 0",
  "- GSE111974 independent external validation: NO",
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
    "STEP04D_QC1_CARD.md"
  )
)

write_text(
  c(
    "# STEP04D-QC1 interpretation note",
    "",
    "The STEP04D v6 computational run completed successfully, but manual output review identified two table-construction defects: program names were lost in the score-model table, and the cross-context candidate table contained an artificial NA__NA context and misclassified candidates with zero measurable bulk contexts.",
    "",
    "The corrected v1.3 score models preserve all four prespecified program names and control the false-discovery rate across the four programs within each dataset contrast. The corrected candidate table contains only the four authorized biological contexts and explicitly labels candidates not measured in any bulk cohort.",
    "",
    "Single-cell preprocessing is technically adequate for integration. Samples on the review watchlist remain included; they require donor-level review after normalization, doublet assessment and cell-type annotation. No donor is excluded at this stage."
  ),
  file.path(
    dirs$text,
    "STEP04D_QC1_interpretation_note.md"
  )
)

parameters <- list(
  schema_version = "1.0",
  run_id = run_id,
  repair_version =
    repair_version,
  projection_version =
    projection_version,
  candidate_n =
    nrow(
      candidate_convergence
    ),
  score_model_row_n =
    nrow(
      score_models_corrected
    ),
  score_program_n =
    length(
      unique(
        score_models_corrected$program_name
      )
    ),
  sc_sample_n =
    nrow(
      sc_sample_qc
    ),
  sc_review_required_n =
    review_required_n,
  sc_caution_n =
    caution_n,
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
    "STEP04D_QC1_parameters.json"
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
    "sessionInfo_POST_RUN_STEP04D_QC1.txt"
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
      "52_step04D_QC1_v2_repair_and_lock_outputs.R"
    )
  )
}

# ----------------------------------------------------------------------
# Manifest and short archive
# ----------------------------------------------------------------------

artifact_paths <- c(
  file.path(
    output_root,
    "STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv"
  ),
  file.path(
    output_root,
    "score_models_v1_3.csv"
  ),
  file.path(
    output_root,
    "cross_context_candidates_v1_3.csv"
  ),
  file.path(
    output_root,
    "candidate_context_coverage.csv"
  ),
  file.path(
    output_root,
    "sc_sample_qc_reviewed.csv"
  ),
  file.path(
    output_root,
    "sc_sample_watchlist.csv"
  ),
  file.path(
    output_root,
    "sc_watchlist_summary.csv"
  ),
  file.path(
    output_root,
    "score_significance_summary.csv"
  ),
  file.path(
    output_root,
    "candidate_direction_summary.csv"
  ),
  file.path(
    output_root,
    "STEP04D_QC1_GATE.csv"
  ),
  file.path(
    output_root,
    "STEP04D_QC1_CARD.md"
  ),
  file.path(
    dirs$text,
    "STEP04D_QC1_interpretation_note.md"
  ),
  file.path(
    dirs$parameters,
    "STEP04D_QC1_parameters.json"
  ),
  file.path(
    dirs$environment,
    "sessionInfo_POST_RUN_STEP04D_QC1.txt"
  ),
  file.path(
    dirs$code,
    "52_step04D_QC1_v2_repair_and_lock_outputs.R"
  ),
  log_path
)

missing_artifacts <- artifact_paths[
  !file.exists(
    artifact_paths
  )
]

if (length(missing_artifacts) > 0L) {
  stop(
    "Expected repair artifact is missing:\n",
    paste(
      missing_artifacts,
      collapse = "\n"
    )
  )
}

artifact_manifest <- data.frame(
  path =
    normalize_slashes(
      artifact_paths
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
    "LOCKED",
  stringsAsFactors = FALSE
)

manifest_path <- file.path(
  output_root,
  "STEP04D_QC1_LOCKED_MANIFEST.csv"
)

write_csv_utf8(
  artifact_manifest,
  manifest_path
)

archive_artifacts <- c(
  artifact_paths,
  manifest_path
)

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
      "GATE|CARD|MANIFEST|SUMMARY|WATCHLIST|COVERAGE",
      basename(
        path_value
      ),
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
      basename(
        path_value
      )
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
  stringsAsFactors = FALSE
)

write_csv_utf8(
  archive_manifest,
  file.path(
    archive_root,
    "ctl",
    "STEP04D_QC1_ARCHIVE_MANIFEST.csv"
  )
)

write_text(
  c(
    "# STEP04D-QC1 archive status",
    "",
    paste0(
      "- Run ID: ",
      run_id
    ),
    paste0(
      "- Repair version: ",
      repair_version,
      " LOCKED"
    ),
    paste0(
      "- Projection version: ",
      projection_version,
      " LOCKED"
    ),
    "- Bulk matrices rerun: NO",
    "- Single-cell pre-QC rerun: NO",
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
      "- Accepted correction run: ",
      run_id
    ),
    "- Bulk matrices: BULK_GENE_MATRIX_v1.2 RETAINED",
    paste0(
      "- Bulk projection: ",
      projection_version,
      " LOCKED"
    ),
    "- Single-cell pre-QC: SINGLE_CELL_PREQC_OBJECTS_v1.0 RETAINED_WITH_WATCHLIST",
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
cat("STEP04D-QC1 v2 repair and lock completed.\n")
cat(
  "Repair version: ",
  repair_version,
  " LOCKED\n"
)
cat(
  "Corrected projection: ",
  projection_version,
  " LOCKED\n"
)
cat(
  "Program-score rows/names: ",
  nrow(
    score_models_corrected
  ),
  " / ",
  length(
    unique(
      score_models_corrected$program_name
    )
  ),
  "\n",
  sep = ""
)
cat(
  "Candidate rows: ",
  nrow(
    candidate_convergence
  ),
  "\n"
)
cat(
  "Candidate context coverage 0/1/2/3/4: ",
  paste(
    context_coverage$candidate_n,
    collapse = "/"
  ),
  "\n"
)
cat(
  "Single-cell REVIEW_REQUIRED / CAUTION: ",
  review_required_n,
  " / ",
  caution_n,
  "\n"
)
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
  "  06_results/S4D_qc_lock/BVP_v1_3/STEP04D_QC1_GATE.csv\n"
)
cat(
  "  06_results/S4D_qc_lock/BVP_v1_3/score_models_v1_3.csv\n"
)
cat(
  "  06_results/S4D_qc_lock/BVP_v1_3/cross_context_candidates_v1_3.csv\n"
)
cat(
  "  06_results/S4D_qc_lock/BVP_v1_3/sc_sample_watchlist.csv\n"
)
cat(
  "  06_results/S4D_qc_lock/BVP_v1_3/STEP04D_QC1_CARD.md\n"
)
cat("============================================================\n")
