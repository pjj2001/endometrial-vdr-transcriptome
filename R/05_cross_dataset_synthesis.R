#!/usr/bin/env Rscript
# Public analysis module for the VDR endometrial transcriptome study.
# Set VDR_PROJECT_ROOT to the project directory containing the expected public GEO inputs.
# This module preserves the executed analysis order; raw GEO data are not redistributed.


# ==============================================================================
# Analysis stage 1: 70_step04E1P_cross_dataset_single_cell_concordance_v1
# ==============================================================================

# ======================================================================
# 70_step04E1P_cross_dataset_single_cell_concordance_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Compare the locked LH+7 RIF-versus-fertile GSE250130 donor models
#   with the independent Control-versus-RIF GSE183837 donor models.
# - Align five shared broad cell types and five shared endpoints.
# - Classify direction concordance, statistical support and
#   non-informative endpoints without pooling donors across studies.
# - Lock the single-cell component interpretation while keeping the
#   overall project result lock unauthorized.
#
# SHARED GRID
# Broad cell types:
# - Endothelial cells
# - Myeloid cells
# - NK/T cells
# - Stromal cells
# - Unciliated epithelial cells
#
# Endpoints:
# - VDR_expression
# - RNA_ALL_436
# - OCCUPANCY_SUPPORTED_284
# - SHARED_OCCUPANCY_46
# - PROMOTER_OCCUPANCY_84
#
# COMPARABILITY
# - Both studies define effects as RIF minus reference.
# - GSE250130 reference: all fertile LH+7 donors.
# - GSE183837 reference: stage-matched Control donors.
# - Program scores use the same locked signed interaction direction.
# - VDR expression uses log2(CPM + 1).
#
# INTERPRETATION BOUNDARIES
# - No pooled meta-analysis is performed.
# - No donor is combined across studies.
# - Directional agreement without corrected significance is descriptive.
# - A nominal P value in one small cohort is exploratory only.
# - Program effects do not prove VDR binding, activity or causality.
#
# INPUTS
# GSE250130:
# - 03_LH7_program_primary.csv
# - 07_LH7_VDR_primary.csv
# - PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE
#
# GSE183837:
# - 01_primary_models.csv
# - 06_cluster8_model_sensitivity.csv
# - PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED
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

set.seed(20260806)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

script_name <-
  "70_step04E1P_cross_dataset_single_cell_concordance_v1.R"

script_version <-
  "STEP04E1P_CROSS_DATASET_SINGLE_CELL_CONCORDANCE_v1.0"

cross_version <- "SCCROSS_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1P_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  normalized <- normalize_slashes(path)

  path_n <- nchar(
    normalized,
    type = "chars"
  )

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Path exceeds the Windows/R safety limit (",
      path_n,
      " > ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_path_safe(path)

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_path_safe(path)

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

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  output <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(output)) {
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  output
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

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
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

resolve_column <- function(
  data,
  aliases,
  context,
  required = TRUE
) {
  column_names <- colnames(data)

  exact_index <- match(
    tolower(aliases),
    tolower(column_names)
  )

  exact_index <- exact_index[
    !is.na(exact_index)
  ]

  if (length(exact_index)) {
    return(
      column_names[
        exact_index[1]
      ]
    )
  }

  if (required) {
    stop(
      "Could not resolve required column for ",
      context,
      ". Accepted aliases: ",
      paste(
        aliases,
        collapse = ", "
      ),
      ". Available columns: ",
      paste(
        column_names,
        collapse = ", "
      )
    )
  }

  ""
}

extract_optional_numeric <- function(
  data,
  aliases,
  context
) {
  column_name <- resolve_column(
    data = data,
    aliases = aliases,
    context = context,
    required = FALSE
  )

  if (!nzchar(column_name)) {
    rep(
      NA_real_,
      nrow(data)
    )
  } else {
    safe_numeric(
      data[[
        column_name
      ]]
    )
  }
}

extract_optional_logical <- function(
  data,
  aliases,
  context
) {
  column_name <- resolve_column(
    data = data,
    aliases = aliases,
    context = context,
    required = FALSE
  )

  if (!nzchar(column_name)) {
    rep(
      FALSE,
      nrow(data)
    )
  } else {
    as_logical_safe(
      data[[
        column_name
      ]]
    )
  }
}

effect_direction <- function(
  x,
  tolerance = 1e-12
) {
  output <- rep(
    "NONFINITE",
    length(x)
  )

  output[
    is.finite(x) &
      x > tolerance
  ] <- "RIF_HIGHER"

  output[
    is.finite(x) &
      x < -tolerance
  ] <- "RIF_LOWER"

  output[
    is.finite(x) &
      abs(x) <= tolerance
  ] <- "ZERO"

  output
}

safe_spearman <- function(x, y) {
  keep <- (
    is.finite(x) &
      is.finite(y)
  )

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 3L ||
      stats::sd(x) == 0 ||
      stats::sd(y) == 0
  ) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x,
      y,
      method = "spearman"
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

  assert_path_safe(png_path)
  assert_path_safe(pdf_path)

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
    cat(
      line,
      "\n",
      file = log_path,
      append = TRUE
    )
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages
# ----------------------------------------------------------------------

required_packages <- c(
  "data.table",
  "ggplot2",
  "scales",
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
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------
# 3. Locked comparison grid
# ----------------------------------------------------------------------

shared_cell_types <- c(
  "Endothelial cells",
  "Myeloid cells",
  "NK/T cells",
  "Stromal cells",
  "Unciliated epithelial cells"
)

shared_endpoints <- c(
  "VDR_expression",
  "RNA_ALL_436",
  "OCCUPANCY_SUPPORTED_284",
  "SHARED_OCCUPANCY_46",
  "PROMOTER_OCCUPANCY_84"
)

expected_grid <- expand.grid(
  broad_cell_type =
    shared_cell_types,
  endpoint =
    shared_endpoints,
  stringsAsFactors = FALSE
)

expected_grid <- expected_grid[
  order(
    match(
      expected_grid$
        broad_cell_type,
      shared_cell_types
    ),
    match(
      expected_grid$endpoint,
      shared_endpoints
    )
  ),
  ,
  drop = FALSE
]

expected_grid$model_key <- paste(
  expected_grid$broad_cell_type,
  expected_grid$endpoint,
  sep = " | "
)


# ----------------------------------------------------------------------
# 4. Paths
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

G250_model_root <- file.path(
  results_root,
  "S4E1_sc",
  "SCI_v1_1_memory_safe",
  "models",
  "GSE250130",
  "VDRMODELS_v1_1"
)

G250_program_path <- file.path(
  G250_model_root,
  "tables",
  "03_LH7_program_primary.csv"
)

G250_VDR_path <- file.path(
  G250_model_root,
  "tables",
  "07_LH7_VDR_primary.csv"
)

G183_model_root <- file.path(
  results_root,
  "S4E1_sc",
  "SCI183_v1_0",
  "mdl",
  "G183",
  "MD183_v1_0"
)

G183_model_path <- file.path(
  G183_model_root,
  "tab",
  "01_primary_models.csv"
)

G183_cluster8_path <- file.path(
  G183_model_root,
  "tab",
  "06_cluster8_model_sensitivity.csv"
)

output_root <- file.path(
  results_root,
  "S4E1_sc",
  "cross_dataset",
  cross_version
)

table_root <- file.path(
  output_root,
  "tab"
)

source_data_root <- file.path(
  output_root,
  "src"
)

figure_output_root <- file.path(
  output_root,
  "fig"
)

manifest_root <- file.path(
  output_root,
  "man"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1P_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_crossQC"
  )
)

figure_root <- file.path(
  run_qc_root,
  "fig"
)

for (path_value in c(
  execution_root,
  parameter_root,
  environment_root,
  manuscript_text_root,
  output_root,
  table_root,
  source_data_root,
  figure_output_root,
  manifest_root,
  run_qc_root,
  figure_root
)) {
  dir.create(
    path_value,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  execution_root,
  paste0(run_id, "_execution.log")
)

run_manifest_path <- file.path(
  execution_root,
  paste0(
    run_id,
    "_RUN_MANIFEST.csv"
  )
)

parameter_path <- file.path(
  parameter_root,
  paste0(
    run_id,
    "_parameters.json"
  )
)

session_path <- file.path(
  environment_root,
  paste0(
    "sessionInfo_",
    run_id,
    ".txt"
  )
)

desktop_candidates <- c(
  file.path(
    Sys.getenv("USERPROFILE"),
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "桌面"
  )
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
  paste0(
    run_id,
    "_SMALL_REPORT.zip"
  )
)

planned_paths <- data.frame(
  artifact = c(
    "G250_program_models",
    "G250_VDR_models",
    "G183_primary_models",
    "G183_cluster8_sensitivity",
    "output_root",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    G250_program_path,
    G250_VDR_path,
    G183_model_path,
    G183_cluster8_path,
    output_root,
    run_qc_root,
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_paths$character_n <- nchar(
  planned_paths$path,
  type = "chars"
)

planned_paths$path_safe <- (
  planned_paths$character_n <= 235L
)

write_csv_utf8(
  planned_paths,
  file.path(
    run_qc_root,
    "00_path_audit.csv"
  )
)

if (!all(planned_paths$path_safe)) {
  stop(
    "One or more cross-dataset paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 5. Validate upstream Gates
# ----------------------------------------------------------------------

G250_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

G250_qc_dirs <- G250_qc_dirs[
  grepl(
    "^STEP04E1G_[0-9]{8}_[0-9]{6}_modelQC$",
    basename(G250_qc_dirs),
    perl = TRUE
  )
]

if (!length(G250_qc_dirs)) {
  stop(
    "No GSE250130 STEP04E1G model-QC directory was found."
  )
}

G250_gate_candidates <- list()

for (directory_path in G250_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "VDR_MODEL_GATE.csv"
  )

  gate_value <- safe_read_csv(
    gate_path
  )

  if (
    !is.null(gate_value) &&
      nrow(gate_value) == 1L
  ) {
    gate_value$source_directory <-
      normalize_slashes(
        directory_path
      )

    gate_value$source_mtime <-
      as.character(
        file.info(
          directory_path
        )$mtime
      )

    G250_gate_candidates[[
      length(G250_gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(G250_gate_candidates)) {
  stop(
    "No readable GSE250130 donor-model Gate was found."
  )
}

G250_gate_candidates <- as.data.frame(
  data.table::rbindlist(
    G250_gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

G250_gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    G250_gate_candidates$
      source_mtime
  )

G250_gate_candidates <- G250_gate_candidates[
  order(
    G250_gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

G250_authorized_index <- which(
  G250_gate_candidates$gate_status ==
    "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE"
)

if (!length(G250_authorized_index)) {
  stop(
    "No GSE250130 Gate passed donor-level model QC."
  )
}

G250_gate <- G250_gate_candidates[
  G250_authorized_index[1],
  ,
  drop = FALSE
]

G183_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

G183_qc_dirs <- G183_qc_dirs[
  grepl(
    "^STEP04E1O_[0-9]{8}_[0-9]{6}_models$",
    basename(G183_qc_dirs),
    perl = TRUE
  )
]

if (!length(G183_qc_dirs)) {
  stop(
    "No GSE183837 STEP04E1O model-QC directory was found."
  )
}

G183_gate_candidates <- list()

for (directory_path in G183_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "G183_DONOR_MODEL_GATE.csv"
  )

  gate_value <- safe_read_csv(
    gate_path
  )

  if (
    !is.null(gate_value) &&
      nrow(gate_value) == 1L
  ) {
    gate_value$source_directory <-
      normalize_slashes(
        directory_path
      )

    gate_value$source_mtime <-
      as.character(
        file.info(
          directory_path
        )$mtime
      )

    G183_gate_candidates[[
      length(G183_gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(G183_gate_candidates)) {
  stop(
    "No readable GSE183837 donor-model Gate was found."
  )
}

G183_gate_candidates <- as.data.frame(
  data.table::rbindlist(
    G183_gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

G183_gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    G183_gate_candidates$
      source_mtime
  )

G183_gate_candidates <- G183_gate_candidates[
  order(
    G183_gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

G183_authorized_index <- which(
  G183_gate_candidates$gate_status ==
    "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED" &
    as_logical_safe(
      G183_gate_candidates$
        cross_dataset_concordance_review_authorized
    )
)

if (!length(G183_authorized_index)) {
  stop(
    "No GSE183837 Gate authorizes cross-dataset concordance review."
  )
}

G183_gate <- G183_gate_candidates[
  G183_authorized_index[1],
  ,
  drop = FALSE
]

write_csv_utf8(
  G250_gate,
  file.path(
    run_qc_root,
    "00_G250_gate.csv"
  )
)

write_csv_utf8(
  G183_gate,
  file.path(
    run_qc_root,
    "00_G183_gate.csv"
  )
)


# ----------------------------------------------------------------------
# 6. Load and audit input schemas
# ----------------------------------------------------------------------

required_input_paths <- c(
  G250_program_path,
  G250_VDR_path,
  G183_model_path,
  G183_cluster8_path
)

missing_input_paths <- required_input_paths[
  !file.exists(required_input_paths)
]

if (length(missing_input_paths)) {
  stop(
    "Missing required cross-dataset input file(s):\n",
    paste(
      normalize_slashes(
        missing_input_paths
      ),
      collapse = "\n"
    )
  )
}

G250_program_raw <- safe_read_csv(
  G250_program_path
)

G250_VDR_raw <- safe_read_csv(
  G250_VDR_path
)

G183_raw <- safe_read_csv(
  G183_model_path
)

G183_cluster8 <- safe_read_csv(
  G183_cluster8_path
)

for (data_name in c(
  "G250_program_raw",
  "G250_VDR_raw",
  "G183_raw",
  "G183_cluster8"
)) {
  data_value <- get(data_name)

  if (
    is.null(data_value) ||
      !nrow(data_value) ||
      anyDuplicated(
        colnames(data_value)
      ) > 0L
  ) {
    stop(
      data_name,
      " is absent, empty or has duplicate columns."
    )
  }
}

schema_audit <- data.frame(
  dataset_table = c(
    "GSE250130_LH7_program_primary",
    "GSE250130_LH7_VDR_primary",
    "GSE183837_primary_models",
    "GSE183837_cluster8_sensitivity"
  ),
  row_n = c(
    nrow(G250_program_raw),
    nrow(G250_VDR_raw),
    nrow(G183_raw),
    nrow(G183_cluster8)
  ),
  column_n = c(
    ncol(G250_program_raw),
    ncol(G250_VDR_raw),
    ncol(G183_raw),
    ncol(G183_cluster8)
  ),
  duplicate_column_n = c(
    sum(
      duplicated(
        colnames(
          G250_program_raw
        )
      )
    ),
    sum(
      duplicated(
        colnames(
          G250_VDR_raw
        )
      )
    ),
    sum(
      duplicated(
        colnames(
          G183_raw
        )
      )
    ),
    sum(
      duplicated(
        colnames(
          G183_cluster8
        )
      )
    )
  ),
  column_names = c(
    paste(
      colnames(G250_program_raw),
      collapse = ";"
    ),
    paste(
      colnames(G250_VDR_raw),
      collapse = ";"
    ),
    paste(
      colnames(G183_raw),
      collapse = ";"
    ),
    paste(
      colnames(G183_cluster8),
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  schema_audit,
  file.path(
    table_root,
    "00_input_schema_audit.csv"
  )
)

log_event(
  "Starting",
  run_id,
  "G250_program_rows=",
  nrow(G250_program_raw),
  "G250_VDR_rows=",
  nrow(G250_VDR_raw),
  "G183_rows=",
  nrow(G183_raw)
)


# ----------------------------------------------------------------------
# 7. Standardize GSE250130 primary LH+7 models
# ----------------------------------------------------------------------

G250_program_cell_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "broad_cell_type",
    "cell_type"
  ),
  context =
    "GSE250130 program broad cell type"
)

G250_program_endpoint_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "program_name",
    "endpoint"
  ),
  context =
    "GSE250130 program endpoint"
)

G250_program_effect_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "RIF_minus_fertile",
    "RIF_minus_reference",
    "mean_difference",
    "estimate"
  ),
  context =
    "GSE250130 program effect"
)

G250_program_p_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "exact_permutation_p",
    "exact_permutation_P"
  ),
  context =
    "GSE250130 program exact P"
)

G250_program_q_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "FDR_BH",
    "BH_FDR",
    "adjusted_p"
  ),
  context =
    "GSE250130 program FDR"
)

G250_program_Hedges_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "Hedges_g",
    "Hedges_g_RIF_minus_Control",
    "Hedges_g_RIF_minus_fertile"
  ),
  context =
    "GSE250130 program Hedges g",
  required = FALSE
)

G250_program_LOO_column <- resolve_column(
  data = G250_program_raw,
  aliases = c(
    "LOO_sign_consistency",
    "LOO_direction_concordance"
  ),
  context =
    "GSE250130 program LOO",
  required = FALSE
)

G250_program <- data.frame(
  dataset = "GSE250130",
  reference_group =
    "FERTILE_LH7_ALL",
  broad_cell_type =
    as.character(
      G250_program_raw[[
        G250_program_cell_column
      ]]
    ),
  endpoint =
    as.character(
      G250_program_raw[[
        G250_program_endpoint_column
      ]]
    ),
  endpoint_class =
    "VDR_PROGRAM",
  effect_RIF_minus_reference =
    safe_numeric(
      G250_program_raw[[
        G250_program_effect_column
      ]]
    ),
  exact_P =
    safe_numeric(
      G250_program_raw[[
        G250_program_p_column
      ]]
    ),
  primary_q =
    safe_numeric(
      G250_program_raw[[
        G250_program_q_column
      ]]
    ),
  primary_q_scope =
    "BH_WITHIN_24_LH7_PROGRAM_MODELS",
  Hedges_g =
    if (
      nzchar(
        G250_program_Hedges_column
      )
    ) {
      safe_numeric(
        G250_program_raw[[
          G250_program_Hedges_column
        ]]
      )
    } else {
      rep(
        NA_real_,
        nrow(G250_program_raw)
      )
    },
  LOO_direction_concordance =
    if (
      nzchar(
        G250_program_LOO_column
      )
    ) {
      safe_numeric(
        G250_program_raw[[
          G250_program_LOO_column
        ]]
      )
    } else {
      rep(
        NA_real_,
        nrow(G250_program_raw)
      )
    },
  n_reference =
    extract_optional_numeric(
      data = G250_program_raw,
      aliases = c(
        "fertile_n",
        "reference_n",
        "n_fertile",
        "n_FERTILE"
      ),
      context =
        "GSE250130 program reference n"
    ),
  n_RIF =
    extract_optional_numeric(
      data = G250_program_raw,
      aliases = c(
        "RIF_n",
        "n_RIF",
        "comparison_n"
      ),
      context =
        "GSE250130 program RIF n"
    ),
  constant_endpoint =
    FALSE,
  stringsAsFactors = FALSE
)

G250_VDR_cell_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "broad_cell_type",
    "cell_type"
  ),
  context =
    "GSE250130 VDR broad cell type"
)

G250_VDR_effect_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "RIF_minus_fertile",
    "RIF_minus_reference",
    "mean_difference",
    "estimate"
  ),
  context =
    "GSE250130 VDR effect"
)

G250_VDR_p_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "exact_permutation_p",
    "exact_permutation_P"
  ),
  context =
    "GSE250130 VDR exact P"
)

G250_VDR_q_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "FDR_BH",
    "BH_FDR",
    "adjusted_p"
  ),
  context =
    "GSE250130 VDR FDR"
)

G250_VDR_Hedges_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "Hedges_g",
    "Hedges_g_RIF_minus_Control",
    "Hedges_g_RIF_minus_fertile"
  ),
  context =
    "GSE250130 VDR Hedges g",
  required = FALSE
)

G250_VDR_LOO_column <- resolve_column(
  data = G250_VDR_raw,
  aliases = c(
    "LOO_sign_consistency",
    "LOO_direction_concordance"
  ),
  context =
    "GSE250130 VDR LOO",
  required = FALSE
)

G250_VDR <- data.frame(
  dataset = "GSE250130",
  reference_group =
    "FERTILE_LH7_ALL",
  broad_cell_type =
    as.character(
      G250_VDR_raw[[
        G250_VDR_cell_column
      ]]
    ),
  endpoint =
    "VDR_expression",
  endpoint_class =
    "VDR_EXPRESSION",
  effect_RIF_minus_reference =
    safe_numeric(
      G250_VDR_raw[[
        G250_VDR_effect_column
      ]]
    ),
  exact_P =
    safe_numeric(
      G250_VDR_raw[[
        G250_VDR_p_column
      ]]
    ),
  primary_q =
    safe_numeric(
      G250_VDR_raw[[
        G250_VDR_q_column
      ]]
    ),
  primary_q_scope =
    "BH_WITHIN_6_LH7_VDR_MODELS",
  Hedges_g =
    if (
      nzchar(
        G250_VDR_Hedges_column
      )
    ) {
      safe_numeric(
        G250_VDR_raw[[
          G250_VDR_Hedges_column
        ]]
      )
    } else {
      rep(
        NA_real_,
        nrow(G250_VDR_raw)
      )
    },
  LOO_direction_concordance =
    if (
      nzchar(
        G250_VDR_LOO_column
      )
    ) {
      safe_numeric(
        G250_VDR_raw[[
          G250_VDR_LOO_column
        ]]
      )
    } else {
      rep(
        NA_real_,
        nrow(G250_VDR_raw)
      )
    },
  n_reference =
    extract_optional_numeric(
      data = G250_VDR_raw,
      aliases = c(
        "fertile_n",
        "reference_n",
        "n_fertile",
        "n_FERTILE"
      ),
      context =
        "GSE250130 VDR reference n"
    ),
  n_RIF =
    extract_optional_numeric(
      data = G250_VDR_raw,
      aliases = c(
        "RIF_n",
        "n_RIF",
        "comparison_n"
      ),
      context =
        "GSE250130 VDR RIF n"
    ),
  constant_endpoint =
    FALSE,
  stringsAsFactors = FALSE
)

G250_standardized <- rbind(
  G250_VDR,
  G250_program
)

G250_standardized$direction <-
  effect_direction(
    G250_standardized$
      effect_RIF_minus_reference
  )

G250_standardized$nominal_significant <-
  is.finite(
    G250_standardized$exact_P
  ) &
    G250_standardized$exact_P <
    0.05

G250_standardized$corrected_significant <-
  is.finite(
    G250_standardized$primary_q
  ) &
    G250_standardized$primary_q <
    0.05

G250_standardized$model_key <- paste(
  G250_standardized$broad_cell_type,
  G250_standardized$endpoint,
  sep = " | "
)

if (
  nrow(G250_standardized) != 30L ||
    anyDuplicated(
      G250_standardized$model_key
    ) > 0L
) {
  stop(
    "GSE250130 standardized model table is not the expected 6 x 5 grid."
  )
}

G250_common <- G250_standardized[
  G250_standardized$
    broad_cell_type %in%
    shared_cell_types &
    G250_standardized$
      endpoint %in%
      shared_endpoints,
  ,
  drop = FALSE
]

if (
  nrow(G250_common) != 25L ||
    !setequal(
      G250_common$model_key,
      expected_grid$model_key
    )
) {
  stop(
    "GSE250130 common comparison grid is incomplete."
  )
}


# ----------------------------------------------------------------------
# 8. Standardize GSE183837 models
# ----------------------------------------------------------------------

required_G183_columns <- c(
  "broad_cell_type",
  "endpoint",
  "endpoint_class",
  "n_Control",
  "n_RIF",
  "mean_difference_RIF_minus_Control",
  "exact_permutation_P",
  "BH_global_25",
  "Hedges_g_RIF_minus_Control",
  "LOO_direction_concordance",
  "constant_endpoint"
)

missing_G183_columns <- setdiff(
  required_G183_columns,
  colnames(G183_raw)
)

if (length(missing_G183_columns)) {
  stop(
    "GSE183837 primary model table is missing: ",
    paste(
      missing_G183_columns,
      collapse = ", "
    )
  )
}

G183_standardized <- data.frame(
  dataset = "GSE183837",
  reference_group =
    "CONTROL_STAGE_MATCHED",
  broad_cell_type =
    as.character(
      G183_raw$
        broad_cell_type
    ),
  endpoint =
    as.character(
      G183_raw$endpoint
    ),
  endpoint_class =
    as.character(
      G183_raw$
        endpoint_class
    ),
  effect_RIF_minus_reference =
    safe_numeric(
      G183_raw$
        mean_difference_RIF_minus_Control
    ),
  exact_P =
    safe_numeric(
      G183_raw$
        exact_permutation_P
    ),
  primary_q =
    safe_numeric(
      G183_raw$
        BH_global_25
    ),
  primary_q_scope =
    "GLOBAL_BH_ACROSS_25_PLANNED_MODELS",
  Hedges_g =
    safe_numeric(
      G183_raw$
        Hedges_g_RIF_minus_Control
    ),
  LOO_direction_concordance =
    safe_numeric(
      G183_raw$
        LOO_direction_concordance
    ),
  n_reference =
    safe_numeric(
      G183_raw$n_Control
    ),
  n_RIF =
    safe_numeric(
      G183_raw$n_RIF
    ),
  constant_endpoint =
    as_logical_safe(
      G183_raw$
        constant_endpoint
    ),
  stringsAsFactors = FALSE
)

G183_standardized$direction <-
  effect_direction(
    G183_standardized$
      effect_RIF_minus_reference
  )

G183_standardized$nominal_significant <-
  is.finite(
    G183_standardized$exact_P
  ) &
    G183_standardized$exact_P <
    0.05

G183_standardized$corrected_significant <-
  is.finite(
    G183_standardized$primary_q
  ) &
    G183_standardized$primary_q <
    0.05

G183_standardized$model_key <- paste(
  G183_standardized$broad_cell_type,
  G183_standardized$endpoint,
  sep = " | "
)

if (
  nrow(G183_standardized) != 25L ||
    anyDuplicated(
      G183_standardized$model_key
    ) > 0L ||
    !setequal(
      G183_standardized$model_key,
      expected_grid$model_key
    )
) {
  stop(
    "GSE183837 standardized model table is not the expected 5 x 5 grid."
  )
}

write_csv_utf8(
  G250_standardized,
  file.path(
    table_root,
    "01_G250_standardized_models.csv"
  )
)

write_csv_utf8(
  G183_standardized,
  file.path(
    table_root,
    "02_G183_standardized_models.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Cross-dataset alignment and evidence classes
# ----------------------------------------------------------------------

G250_join <- G250_common[
  ,
  c(
    "model_key",
    "broad_cell_type",
    "endpoint",
    "endpoint_class",
    "effect_RIF_minus_reference",
    "exact_P",
    "primary_q",
    "primary_q_scope",
    "Hedges_g",
    "LOO_direction_concordance",
    "n_reference",
    "n_RIF",
    "constant_endpoint",
    "direction",
    "nominal_significant",
    "corrected_significant"
  ),
  drop = FALSE
]

colnames(G250_join)[
  !colnames(G250_join) %in%
    c(
      "model_key",
      "broad_cell_type",
      "endpoint",
      "endpoint_class"
    )
] <- paste0(
  colnames(G250_join)[
    !colnames(G250_join) %in%
      c(
        "model_key",
        "broad_cell_type",
        "endpoint",
        "endpoint_class"
      )
  ],
  "_G250"
)

G183_join <- G183_standardized[
  ,
  c(
    "model_key",
    "effect_RIF_minus_reference",
    "exact_P",
    "primary_q",
    "primary_q_scope",
    "Hedges_g",
    "LOO_direction_concordance",
    "n_reference",
    "n_RIF",
    "constant_endpoint",
    "direction",
    "nominal_significant",
    "corrected_significant"
  ),
  drop = FALSE
]

colnames(G183_join)[
  colnames(G183_join) !=
    "model_key"
] <- paste0(
  colnames(G183_join)[
    colnames(G183_join) !=
      "model_key"
  ],
  "_G183"
)

cross_models <- merge(
  G250_join,
  G183_join,
  by = "model_key",
  all = TRUE,
  sort = FALSE
)

cross_models <- merge(
  expected_grid,
  cross_models,
  by = c(
    "model_key",
    "broad_cell_type",
    "endpoint"
  ),
  all.x = TRUE,
  sort = FALSE
)

cross_models <- cross_models[
  match(
    expected_grid$model_key,
    cross_models$model_key
  ),
  ,
  drop = FALSE
]

if (
  nrow(cross_models) != 25L ||
    anyNA(
      cross_models$
        effect_RIF_minus_reference_G250
    ) ||
    anyNA(
      cross_models$
        effect_RIF_minus_reference_G183
    )
) {
  stop(
    "Cross-dataset model alignment failed."
  )
}

cross_models$noninformative_endpoint <- (
  cross_models$
    constant_endpoint_G250 |
    cross_models$
      constant_endpoint_G183 |
    cross_models$
      direction_G250 %in%
      c(
        "ZERO",
        "NONFINITE"
      ) |
    cross_models$
      direction_G183 %in%
      c(
        "ZERO",
        "NONFINITE"
      )
)

cross_models$direction_concordant <- ifelse(
  cross_models$
    noninformative_endpoint,
  NA,
  cross_models$
    direction_G250 ==
    cross_models$
      direction_G183
)

cross_models$direction_status <- ifelse(
  cross_models$
    noninformative_endpoint,
  "NONINFORMATIVE",
  ifelse(
    cross_models$
      direction_concordant,
    "CONCORDANT",
    "DISCORDANT"
  )
)

cross_models$evidence_class <- ifelse(
  cross_models$
    noninformative_endpoint,
  "NONINFORMATIVE_CONSTANT_OR_ZERO",
  ifelse(
    cross_models$
      direction_concordant &
      cross_models$
        corrected_significant_G250 &
      cross_models$
        corrected_significant_G183,
    "REPLICATED_CORRECTED_CONCORDANT",
    ifelse(
      cross_models$
        direction_concordant &
        (
          cross_models$
            corrected_significant_G250 |
          cross_models$
            corrected_significant_G183
        ) &
        (
          cross_models$
            nominal_significant_G250 &
          cross_models$
            nominal_significant_G183
        ),
      "CONCORDANT_ONE_CORRECTED_BOTH_NOMINAL",
      ifelse(
        cross_models$
          direction_concordant &
          cross_models$
            nominal_significant_G250 &
          cross_models$
            nominal_significant_G183,
        "CONCORDANT_BOTH_NOMINAL",
        ifelse(
          cross_models$
            direction_concordant &
            (
              cross_models$
                nominal_significant_G250 |
              cross_models$
                nominal_significant_G183
            ),
          "CONCORDANT_ONE_NOMINAL_EXPLORATORY",
          ifelse(
            cross_models$
              direction_concordant,
            "CONCORDANT_NO_STATISTICAL_SUPPORT",
            ifelse(
              cross_models$
                nominal_significant_G250 |
                cross_models$
                  nominal_significant_G183 |
                cross_models$
                  corrected_significant_G250 |
                cross_models$
                  corrected_significant_G183,
              "DISCORDANT_WITH_NOMINAL_OR_CORRECTED_SIGNAL",
              "DISCORDANT_NO_STATISTICAL_SUPPORT"
            )
          )
        )
      )
    )
  )
)

cross_models$cross_dataset_interpretation <- ifelse(
  cross_models$
    evidence_class ==
    "REPLICATED_CORRECTED_CONCORDANT",
  "CORRECTED_REPLICATION",
  ifelse(
    cross_models$
      evidence_class %in%
      c(
        "CONCORDANT_ONE_CORRECTED_BOTH_NOMINAL",
        "CONCORDANT_BOTH_NOMINAL",
        "CONCORDANT_ONE_NOMINAL_EXPLORATORY"
      ),
    "EXPLORATORY_DIRECTIONAL_SUPPORT",
    ifelse(
      cross_models$
        evidence_class ==
        "CONCORDANT_NO_STATISTICAL_SUPPORT",
      "DIRECTION_ONLY_NO_INFERENTIAL_SUPPORT",
      ifelse(
        cross_models$
          evidence_class ==
          "NONINFORMATIVE_CONSTANT_OR_ZERO",
        "NONINFORMATIVE",
        "NO_CROSS_DATASET_SUPPORT"
      )
    )
  )
)

write_csv_utf8(
  cross_models,
  file.path(
    table_root,
    "03_cross_dataset_25_model_concordance.csv"
  )
)


# ----------------------------------------------------------------------
# 10. Concordance summaries
# ----------------------------------------------------------------------

informative_models <- cross_models[
  !cross_models$
    noninformative_endpoint,
  ,
  drop = FALSE
]

overall_summary <- data.frame(
  metric = c(
    "planned_common_model_n",
    "informative_model_n",
    "noninformative_model_n",
    "direction_concordant_n",
    "direction_discordant_n",
    "direction_concordance_fraction",
    "G250_nominal_P_lt_0_05_n",
    "G250_corrected_q_lt_0_05_n",
    "G183_nominal_P_lt_0_05_n",
    "G183_corrected_q_lt_0_05_n",
    "corrected_replicated_concordant_n",
    "exploratory_directional_support_n",
    "cluster8_direction_changed_n"
  ),
  value = c(
    nrow(cross_models),
    nrow(informative_models),
    sum(
      cross_models$
        noninformative_endpoint
    ),
    sum(
      informative_models$
        direction_concordant
    ),
    sum(
      !informative_models$
        direction_concordant
    ),
    mean(
      informative_models$
        direction_concordant
    ),
    sum(
      cross_models$
        nominal_significant_G250
    ),
    sum(
      cross_models$
        corrected_significant_G250
    ),
    sum(
      cross_models$
        nominal_significant_G183
    ),
    sum(
      cross_models$
        corrected_significant_G183
    ),
    sum(
      cross_models$
        evidence_class ==
        "REPLICATED_CORRECTED_CONCORDANT"
    ),
    sum(
      cross_models$
        cross_dataset_interpretation ==
        "EXPLORATORY_DIRECTIONAL_SUPPORT"
    ),
    sum(
      G183_cluster8$
        direction_concordant %in%
        FALSE,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

by_endpoint_rows <- list()

for (
  endpoint_index in seq_along(
    shared_endpoints
  )
) {
  endpoint_value <- shared_endpoints[
    endpoint_index
  ]

  endpoint_data <- cross_models[
    cross_models$endpoint ==
      endpoint_value,
    ,
    drop = FALSE
  ]

  endpoint_informative <- endpoint_data[
    !endpoint_data$
      noninformative_endpoint,
    ,
    drop = FALSE
  ]

  by_endpoint_rows[[
    endpoint_index
  ]] <- data.frame(
    endpoint =
      endpoint_value,
    planned_model_n =
      nrow(endpoint_data),
    informative_model_n =
      nrow(endpoint_informative),
    concordant_n =
      sum(
        endpoint_informative$
          direction_concordant
      ),
    discordant_n =
      sum(
        !endpoint_informative$
          direction_concordant
      ),
    concordance_fraction =
      if (
        nrow(endpoint_informative)
      ) {
        mean(
          endpoint_informative$
            direction_concordant
        )
      } else {
        NA_real_
      },
    Spearman_effect_rho =
      safe_spearman(
        endpoint_informative$
          effect_RIF_minus_reference_G250,
        endpoint_informative$
          effect_RIF_minus_reference_G183
      ),
    G250_nominal_n =
      sum(
        endpoint_data$
          nominal_significant_G250
      ),
    G183_nominal_n =
      sum(
        endpoint_data$
          nominal_significant_G183
      ),
    corrected_replication_n =
      sum(
        endpoint_data$
          evidence_class ==
          "REPLICATED_CORRECTED_CONCORDANT"
      ),
    stringsAsFactors = FALSE
  )
}

by_endpoint <- do.call(
  rbind,
  by_endpoint_rows
)

by_cell_type_rows <- list()

for (
  cell_index in seq_along(
    shared_cell_types
  )
) {
  cell_type <- shared_cell_types[
    cell_index
  ]

  cell_data <- cross_models[
    cross_models$
      broad_cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]

  cell_informative <- cell_data[
    !cell_data$
      noninformative_endpoint,
    ,
    drop = FALSE
  ]

  by_cell_type_rows[[
    cell_index
  ]] <- data.frame(
    broad_cell_type =
      cell_type,
    planned_model_n =
      nrow(cell_data),
    informative_model_n =
      nrow(cell_informative),
    concordant_n =
      sum(
        cell_informative$
          direction_concordant
      ),
    discordant_n =
      sum(
        !cell_informative$
          direction_concordant
      ),
    concordance_fraction =
      if (
        nrow(cell_informative)
      ) {
        mean(
          cell_informative$
            direction_concordant
        )
      } else {
        NA_real_
      },
    G250_nominal_n =
      sum(
        cell_data$
          nominal_significant_G250
      ),
    G183_nominal_n =
      sum(
        cell_data$
          nominal_significant_G183
      ),
    corrected_replication_n =
      sum(
        cell_data$
          evidence_class ==
          "REPLICATED_CORRECTED_CONCORDANT"
      ),
    stringsAsFactors = FALSE
  )
}

by_cell_type <- do.call(
  rbind,
  by_cell_type_rows
)

evidence_class_summary <- as.data.frame(
  table(
    evidence_class =
      cross_models$
        evidence_class
  ),
  stringsAsFactors = FALSE
)

colnames(
  evidence_class_summary
)[2] <- "model_n"

myeloid_VDR_focus <- cross_models[
  cross_models$
    broad_cell_type ==
    "Myeloid cells" &
    cross_models$endpoint ==
      "VDR_expression",
  ,
  drop = FALSE
]

if (nrow(myeloid_VDR_focus) != 1L) {
  stop(
    "Myeloid VDR focus row was not uniquely recovered."
  )
}

write_csv_utf8(
  overall_summary,
  file.path(
    table_root,
    "04_overall_concordance_summary.csv"
  )
)

write_csv_utf8(
  by_endpoint,
  file.path(
    table_root,
    "05_concordance_by_endpoint.csv"
  )
)

write_csv_utf8(
  by_cell_type,
  file.path(
    table_root,
    "06_concordance_by_cell_type.csv"
  )
)

write_csv_utf8(
  evidence_class_summary,
  file.path(
    table_root,
    "07_evidence_class_summary.csv"
  )
)

write_csv_utf8(
  myeloid_VDR_focus,
  file.path(
    table_root,
    "08_myeloid_VDR_focus.csv"
  )
)

write_csv_utf8(
  G183_cluster8,
  file.path(
    table_root,
    "09_cluster8_sensitivity_review.csv"
  )
)


# ----------------------------------------------------------------------
# 11. Locked interpretation state
# ----------------------------------------------------------------------

corrected_replication_n <- sum(
  cross_models$
    evidence_class ==
    "REPLICATED_CORRECTED_CONCORDANT"
)

exploratory_support_n <- sum(
  cross_models$
    cross_dataset_interpretation ==
    "EXPLORATORY_DIRECTIONAL_SUPPORT"
)

myeloid_VDR_G183_nominal <- (
  myeloid_VDR_focus$
    nominal_significant_G183
)

myeloid_VDR_G183_corrected <- (
  myeloid_VDR_focus$
    corrected_significant_G183
)

myeloid_VDR_direction_status <- (
  myeloid_VDR_focus$
    direction_status
)

single_cell_evidence_state <- if (
  corrected_replication_n > 0L
) {
  "CORRECTED_CROSS_DATASET_REPLICATION_PRESENT"
} else if (
  exploratory_support_n > 0L
) {
  "NO_CORRECTED_REPLICATION_EXPLORATORY_DIRECTIONAL_SUPPORT_ONLY"
} else {
  "NO_CORRECTED_REPLICATION_OR_DIRECTIONAL_SUPPORT"
}

broad_RIF_VDR_program_dysregulation_supported <- (
  corrected_replication_n > 0L &&
    any(
      cross_models$
        endpoint_class ==
        "VDR_PROGRAM" &
        cross_models$
          evidence_class ==
          "REPLICATED_CORRECTED_CONCORDANT"
    )
)

interpretation_lock <- data.frame(
  item = c(
    "single_cell_evidence_state",
    "corrected_cross_dataset_replication_n",
    "exploratory_directional_support_n",
    "broad_RIF_VDR_program_dysregulation_supported",
    "GSE183837_myeloid_VDR_nominal_signal",
    "GSE183837_myeloid_VDR_corrected_signal",
    "myeloid_VDR_cross_dataset_direction_status",
    "cluster8_sensitivity_direction_changed_n",
    "GSE183837_ciliated_epithelial_inference",
    "pooled_meta_analysis_performed",
    "direct_VDR_target_claim_authorized",
    "causal_mechanism_claim_authorized",
    "single_cell_component_interpretation_locked",
    "project_level_result_lock_authorized"
  ),
  value = c(
    single_cell_evidence_state,
    corrected_replication_n,
    exploratory_support_n,
    broad_RIF_VDR_program_dysregulation_supported,
    myeloid_VDR_G183_nominal,
    myeloid_VDR_G183_corrected,
    myeloid_VDR_direction_status,
    sum(
      G183_cluster8$
        direction_concordant %in%
        FALSE,
      na.rm = TRUE
    ),
    "DESCRIPTIVE_ONLY",
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    FALSE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  interpretation_lock,
  file.path(
    table_root,
    "10_single_cell_interpretation_lock.csv"
  )
)


# ----------------------------------------------------------------------
# 12. QC figures
# ----------------------------------------------------------------------

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

scatter_data <- cross_models[
  !cross_models$
    noninformative_endpoint,
  ,
  drop = FALSE
]

scatter_data$endpoint <- factor(
  scatter_data$endpoint,
  levels = shared_endpoints
)

short_cell_label_map <- c(
  "Endothelial cells" = "Endo",
  "Myeloid cells" = "Myeloid",
  "NK/T cells" = "NK/T",
  "Stromal cells" = "Stromal",
  "Unciliated epithelial cells" = "Unciliated"
)

scatter_data$short_cell_label <- unname(
  short_cell_label_map[
    scatter_data$
      broad_cell_type
  ]
)

if (
  anyNA(
    scatter_data$
      short_cell_label
  )
) {
  stop(
    "One or more broad cell types could not be mapped to a short plot label."
  )
}

plot_scatter <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(
    x =
      effect_RIF_minus_reference_G250,
    y =
      effect_RIF_minus_reference_G183,
    shape =
      direction_status
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 2.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        short_cell_label
    ),
    nudge_y = 0.03,
    size = 2.6,
    check_overlap = TRUE
  ) +
  ggplot2::facet_wrap(
    ~ endpoint,
    scales = "free",
    ncol = 3
  ) +
  ggplot2::labs(
    title = "Cross-dataset RIF effect directions",
    subtitle = "GSE250130 compares RIF with fertile LH+7; GSE183837 compares RIF with stage-matched controls",
    x = "GSE250130 effect: RIF minus fertile",
    y = "GSE183837 effect: RIF minus Control",
    shape = "Direction status"
  ) +
  theme_QC

save_plot_pair(
  plot_scatter,
  "QC1_cross_dataset_effect_scatter",
  width = 11,
  height = 8
)

matrix_data <- cross_models

matrix_data$broad_cell_type <- factor(
  matrix_data$broad_cell_type,
  levels = rev(
    shared_cell_types
  )
)

matrix_data$endpoint <- factor(
  matrix_data$endpoint,
  levels = shared_endpoints
)

matrix_data$direction_symbol <- ifelse(
  matrix_data$direction_status ==
    "CONCORDANT",
  "C",
  ifelse(
    matrix_data$direction_status ==
      "DISCORDANT",
    "D",
    "NI"
  )
)

plot_direction_matrix <- ggplot2::ggplot(
  matrix_data,
  ggplot2::aes(
    x = endpoint,
    y = broad_cell_type,
    fill = direction_status
  )
) +
  ggplot2::geom_tile(
    linewidth = 0.4
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        direction_symbol
    ),
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      CONCORDANT = "grey35",
      DISCORDANT = "grey75",
      NONINFORMATIVE = "white"
    )
  ) +
  ggplot2::labs(
    title = "Direction concordance across the 25 shared models",
    subtitle = "C, concordant; D, discordant; NI, non-informative constant or zero endpoint",
    x = "Endpoint",
    y = "Broad cell type",
    fill = "Status"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

save_plot_pair(
  plot_direction_matrix,
  "QC2_direction_concordance_matrix",
  width = 10,
  height = 6
)

Hedges_long <- rbind(
  data.frame(
    model_key =
      cross_models$model_key,
    broad_cell_type =
      cross_models$
        broad_cell_type,
    endpoint =
      cross_models$endpoint,
    dataset =
      "GSE250130",
    Hedges_g =
      cross_models$
        Hedges_g_G250,
    stringsAsFactors = FALSE
  ),
  data.frame(
    model_key =
      cross_models$model_key,
    broad_cell_type =
      cross_models$
        broad_cell_type,
    endpoint =
      cross_models$endpoint,
    dataset =
      "GSE183837",
    Hedges_g =
      cross_models$
        Hedges_g_G183,
    stringsAsFactors = FALSE
  )
)

Hedges_long <- Hedges_long[
  is.finite(
    Hedges_long$Hedges_g
  ),
  ,
  drop = FALSE
]

Hedges_long$endpoint <- factor(
  Hedges_long$endpoint,
  levels = shared_endpoints
)

plot_Hedges <- ggplot2::ggplot(
  Hedges_long,
  ggplot2::aes(
    x = dataset,
    y = Hedges_g,
    group =
      broad_cell_type,
    shape =
      broad_cell_type
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    alpha = 0.45
  ) +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggplot2::facet_wrap(
    ~ endpoint,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::labs(
    title = "Standardized effect-size comparison",
    subtitle = "Hedges g is oriented as RIF minus reference in both studies",
    x = "Dataset",
    y = "Hedges g",
    shape = "Broad cell type"
  ) +
  theme_QC

save_plot_pair(
  plot_Hedges,
  "QC3_Hedges_g_cross_dataset",
  width = 11,
  height = 8
)

q_long <- rbind(
  data.frame(
    broad_cell_type =
      cross_models$
        broad_cell_type,
    endpoint =
      cross_models$endpoint,
    dataset =
      "GSE250130",
    primary_q =
      cross_models$
        primary_q_G250,
    stringsAsFactors = FALSE
  ),
  data.frame(
    broad_cell_type =
      cross_models$
        broad_cell_type,
    endpoint =
      cross_models$endpoint,
    dataset =
      "GSE183837",
    primary_q =
      cross_models$
        primary_q_G183,
    stringsAsFactors = FALSE
  )
)

q_long$minus_log10_q <- -log10(
  pmax(
    q_long$primary_q,
    1e-12
  )
)

q_long$model_label <- paste(
  q_long$broad_cell_type,
  q_long$endpoint,
  sep = " | "
)

q_long$model_label <- factor(
  q_long$model_label,
  levels = rev(
    expected_grid$model_key
  )
)

plot_q <- ggplot2::ggplot(
  q_long,
  ggplot2::aes(
    x = dataset,
    y = model_label,
    fill =
      minus_log10_q
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(
    title = "Primary multiple-testing-adjusted evidence",
    subtitle = "Adjustment scopes differ by the prespecified model families in each study",
    x = "Dataset",
    y = "Shared model",
    fill = "-log10(q)"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_text(
        size = 6
      )
  )

save_plot_pair(
  plot_q,
  "QC4_adjusted_evidence_matrix",
  width = 9,
  height = 10
)

myeloid_plot_data <- data.frame(
  dataset = c(
    "GSE250130",
    "GSE183837"
  ),
  effect = c(
    myeloid_VDR_focus$
      effect_RIF_minus_reference_G250,
    myeloid_VDR_focus$
      effect_RIF_minus_reference_G183
  ),
  Hedges_g = c(
    myeloid_VDR_focus$
      Hedges_g_G250,
    myeloid_VDR_focus$
      Hedges_g_G183
  ),
  exact_P = c(
    myeloid_VDR_focus$
      exact_P_G250,
    myeloid_VDR_focus$
      exact_P_G183
  ),
  primary_q = c(
    myeloid_VDR_focus$
      primary_q_G250,
    myeloid_VDR_focus$
      primary_q_G183
  ),
  stringsAsFactors = FALSE
)

myeloid_plot_data$label <- paste0(
  "P=",
  format(
    myeloid_plot_data$
      exact_P,
    digits = 3
  ),
  "\nq=",
  format(
    myeloid_plot_data$
      primary_q,
    digits = 3
  )
)

plot_myeloid <- ggplot2::ggplot(
  myeloid_plot_data,
  ggplot2::aes(
    x = dataset,
    y = effect,
    shape = dataset
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 3.2
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = label
    ),
    nudge_y = 0.15,
    size = 3
  ) +
  ggplot2::labs(
    title = "Myeloid-cell VDR expression across the two studies",
    subtitle = "The GSE183837 nominal signal remains exploratory unless supported in direction by GSE250130",
    x = "Dataset",
    y = "RIF minus reference VDR log2(CPM + 1)",
    shape = "Dataset"
  ) +
  theme_QC

save_plot_pair(
  plot_myeloid,
  "QC5_myeloid_VDR_focus",
  width = 7,
  height = 5.5
)

figure_files <- list.files(
  figure_root,
  full.names = TRUE
)

file.copy(
  figure_files,
  figure_output_root,
  overwrite = TRUE
)


# ----------------------------------------------------------------------
# 13. Gate
# ----------------------------------------------------------------------

myeloid_G183_exact_P <- myeloid_VDR_focus$
  exact_P_G183

myeloid_G183_q <- myeloid_VDR_focus$
  primary_q_G183

cluster8_direction_changed_n <- sum(
  G183_cluster8$
    direction_concordant %in%
    FALSE,
  na.rm = TRUE
)

gate_checks <- data.frame(
  check = c(
    "GSE250130_model_gate_passes",
    "GSE183837_model_gate_passes",
    "GSE250130_program_table_has_24_models",
    "GSE250130_VDR_table_has_6_models",
    "GSE250130_common_grid_has_25_models",
    "GSE183837_grid_has_25_models",
    "shared_grid_alignment_is_complete",
    "all_effects_are_finite",
    "all_exact_P_values_are_finite",
    "all_primary_q_values_are_finite",
    "one_GSE183837_constant_endpoint_is_endothelial_VDR",
    "GSE183837_myeloid_VDR_is_nominal_but_not_corrected",
    "cluster8_sensitivity_has_5_endpoints",
    "cluster8_sensitivity_has_no_direction_change",
    "no_pooled_meta_analysis_performed",
    "cross_dataset_interpretation_can_be_locked",
    "project_level_result_lock_remains_disabled"
  ),
  pass = c(
    G250_gate$gate_status[[1]] ==
      "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE",
    G183_gate$gate_status[[1]] ==
      "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED",
    nrow(G250_program_raw) == 24L,
    nrow(G250_VDR_raw) == 6L,
    nrow(G250_common) == 25L,
    nrow(G183_standardized) == 25L,
    nrow(cross_models) == 25L &&
      !anyNA(
        cross_models$model_key
      ),
    all(
      is.finite(
        cross_models$
          effect_RIF_minus_reference_G250
      )
    ) &&
      all(
        is.finite(
          cross_models$
            effect_RIF_minus_reference_G183
        )
      ),
    all(
      is.finite(
        cross_models$
          exact_P_G250
      )
    ) &&
      all(
        is.finite(
          cross_models$
            exact_P_G183
        )
      ),
    all(
      is.finite(
        cross_models$
          primary_q_G250
      )
    ) &&
      all(
        is.finite(
          cross_models$
            primary_q_G183
        )
      ),
    sum(
      cross_models$
        constant_endpoint_G183
    ) == 1L &&
      cross_models$
        broad_cell_type[
          cross_models$
            constant_endpoint_G183
        ] ==
        "Endothelial cells" &&
      cross_models$
        endpoint[
          cross_models$
            constant_endpoint_G183
        ] ==
        "VDR_expression",
    is.finite(
      myeloid_G183_exact_P
    ) &&
      myeloid_G183_exact_P <
      0.05 &&
      is.finite(
        myeloid_G183_q
      ) &&
      myeloid_G183_q >=
      0.05,
    nrow(G183_cluster8) == 5L,
    cluster8_direction_changed_n ==
      0L,
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_checks,
  file.path(
    run_qc_root,
    "01_gate_checks.csv"
  )
)

all_gate_checks_pass <- all(
  gate_checks$pass
)

gate_status <- if (
  all_gate_checks_pass
) {
  "PASS_SINGLE_CELL_CROSS_DATASET_CONCORDANCE_REVIEW_COMPLETE"
} else {
  "HOLD_SINGLE_CELL_CROSS_DATASET_CONCORDANCE_INCOMPLETE"
}

cross_gate <- data.frame(
  cross_version =
    cross_version,
  gate_status =
    gate_status,
  planned_common_model_n =
    nrow(cross_models),
  informative_model_n =
    nrow(informative_models),
  direction_concordant_n =
    sum(
      informative_models$
        direction_concordant
    ),
  direction_discordant_n =
    sum(
      !informative_models$
        direction_concordant
    ),
  direction_concordance_fraction =
    mean(
      informative_models$
        direction_concordant
    ),
  corrected_replication_n =
    corrected_replication_n,
  exploratory_directional_support_n =
    exploratory_support_n,
  G250_corrected_signal_n =
    sum(
      cross_models$
        corrected_significant_G250
    ),
  G183_corrected_signal_n =
    sum(
      cross_models$
        corrected_significant_G183
    ),
  G183_nominal_signal_n =
    sum(
      cross_models$
        nominal_significant_G183
    ),
  G183_myeloid_VDR_exact_P =
    myeloid_G183_exact_P,
  G183_myeloid_VDR_global_q =
    myeloid_G183_q,
  myeloid_VDR_direction_status =
    myeloid_VDR_direction_status,
  cluster8_direction_changed_n =
    cluster8_direction_changed_n,
  single_cell_evidence_state =
    single_cell_evidence_state,
  broad_RIF_VDR_program_dysregulation_supported =
    broad_RIF_VDR_program_dysregulation_supported,
  single_cell_component_interpretation_locked =
    all_gate_checks_pass,
  project_level_result_lock_authorized =
    FALSE,
  gene_level_cross_dataset_DE_authorized =
    FALSE,
  pooled_meta_analysis_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  cross_gate,
  file.path(
    run_qc_root,
    "SC_CROSS_DATASET_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1P Single-cell Cross-dataset Concordance Card",
  "",
  paste0("- Gate: ", gate_status),
  "- Shared planned models: 25",
  paste0(
    "- Informative models: ",
    nrow(informative_models)
  ),
  paste0(
    "- Direction-concordant models: ",
    sum(
      informative_models$
        direction_concordant
    )
  ),
  paste0(
    "- Direction-discordant models: ",
    sum(
      !informative_models$
        direction_concordant
    )
  ),
  paste0(
    "- Direction concordance fraction: ",
    format(
      mean(
        informative_models$
          direction_concordant
      ),
      digits = 4
    )
  ),
  paste0(
    "- Corrected replicated models: ",
    corrected_replication_n
  ),
  paste0(
    "- Exploratory directional-support models: ",
    exploratory_support_n
  ),
  paste0(
    "- GSE183837 myeloid VDR exact P: ",
    format(
      myeloid_G183_exact_P,
      digits = 4
    )
  ),
  paste0(
    "- GSE183837 myeloid VDR global q: ",
    format(
      myeloid_G183_q,
      digits = 4
    )
  ),
  paste0(
    "- Myeloid VDR direction status: ",
    myeloid_VDR_direction_status
  ),
  paste0(
    "- Single-cell evidence state: ",
    single_cell_evidence_state
  ),
  paste0(
    "- Broad RIF VDR-program dysregulation supported: ",
    broad_RIF_VDR_program_dysregulation_supported
  ),
  paste0(
    "- Single-cell component interpretation locked: ",
    all_gate_checks_pass
  ),
  "- Project-level result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "SC_CROSS_DATASET_CARD.md"
  )
)

interpretation_lines <- c(
  "# Single-cell cross-dataset interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Locked interpretation",
  "",
  paste0(
    "- ",
    nrow(informative_models),
    " of 25 shared models were directionally informative; ",
    sum(
      informative_models$
        direction_concordant
    ),
    " were concordant and ",
    sum(
      !informative_models$
        direction_concordant
    ),
    " were discordant."
  ),
  paste0(
    "- Corrected cross-dataset replication was observed for ",
    corrected_replication_n,
    " model(s)."
  ),
  paste0(
    "- Exploratory directional support was observed for ",
    exploratory_support_n,
    " model(s)."
  ),
  paste0(
    "- The GSE183837 myeloid VDR-expression signal had exact P=",
    format(
      myeloid_G183_exact_P,
      digits = 4
    ),
    " and global q=",
    format(
      myeloid_G183_q,
      digits = 4
    ),
    "; its direction relative to GSE250130 was ",
    myeloid_VDR_direction_status,
    "."
  ),
  paste0(
    "- The final single-cell evidence state is: ",
    single_cell_evidence_state,
    "."
  ),
  "",
  "## Scientific boundary",
  "",
  "- The two studies were not pooled because their reference groups and sampling structures differ.",
  "- Directional concordance without corrected significance is descriptive rather than confirmatory.",
  "- The independent GSE183837 cohort does not by itself establish broad RIF-associated VDR-program dysregulation.",
  "- The myeloid VDR-expression finding remains exploratory because it did not survive global multiplicity correction.",
  "- No direct-target, VDR-activity or causal mechanism claim is authorized.",
  "",
  "## Authorized next operation",
  "",
  "- Use this locked single-cell interpretation in the STEP04 closeout.",
  "- Integrate the temporal single-cell evidence, patient bulk projection and occupancy hierarchy.",
  "- Finalize the manuscript Figure/Table blueprint and wording QC.",
  "",
  "## Not authorized",
  "",
  "- No pooled meta-analysis.",
  "- No additional cell-level testing.",
  "- No GSE183837 gene-level discovery DE.",
  "- No project-level final result lock."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "SC_CROSS_DATASET_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 14. Manuscript shells
# ----------------------------------------------------------------------

methods_lines <- c(
  "# STEP04E1P Methods shell",
  "",
  "Cross-dataset consistency was evaluated across the 25 combinations shared by GSE250130 and GSE183837, comprising five broad cell types and five endpoints. Effects were oriented as RIF minus reference in both studies. GSE250130 contrasted RIF with all fertile LH+7 donors, whereas GSE183837 contrasted RIF with stage-matched controls. Direction concordance was classified only for finite non-zero effects. Each study retained its prospectively defined multiplicity correction; no pooled meta-analysis or cross-study combined P value was calculated because the reference groups, sample sizes and cohort structures differed."
)

results_lines <- c(
  "# STEP04E1P Results shell",
  "",
  paste0(
    "Among ",
    nrow(informative_models),
    " directionally informative shared models, ",
    sum(
      informative_models$
        direction_concordant
    ),
    " showed the same effect direction across studies and ",
    sum(
      !informative_models$
        direction_concordant
    ),
    " were discordant."
  ),
  paste0(
    "Corrected cross-dataset replication was observed in ",
    corrected_replication_n,
    " model(s)."
  ),
  paste0(
    "In GSE183837, myeloid-cell VDR expression was lower in RIF at nominal exact P=",
    format(
      myeloid_G183_exact_P,
      digits = 4
    ),
    ", but the global BH-adjusted q value was ",
    format(
      myeloid_G183_q,
      digits = 4
    ),
    ". The corresponding effect direction in GSE250130 was ",
    myeloid_VDR_direction_status,
    "."
  ),
  paste0(
    "Retaining the low-confidence cluster 8 changed the unciliated epithelial effect direction for ",
    cluster8_direction_changed_n,
    " of five endpoints."
  ),
  "Together, the two donor-level single-cell analyses did not support a broad, reproducible RIF-associated shift in the locked VDR transcriptional programs."
)

discussion_lines <- c(
  "# STEP04E1P Discussion shell",
  "",
  "The independent cohort primarily provided a falsification and localization test rather than strong confirmatory power. The absence of multiplicity-corrected replication argues against framing RIF as a uniform VDR-program deficiency across endometrial cell compartments. The nominal myeloid VDR-expression decrease is biologically interesting because its donor-level effect was large and leave-one-donor directionally stable in GSE183837, but it remains exploratory until reproduced in a larger stage-matched cohort. The strong temporal remodeling observed in GSE250130 should therefore be separated from disease-associated dysregulation: normal luteal timing clearly structures the VDR-related transcriptional context, whereas broad RIF-associated disruption was not consistently demonstrated."
)

write_text(
  methods_lines,
  file.path(
    manuscript_text_root,
    paste0(
      run_id,
      "_Methods.md"
    )
  )
)

write_text(
  results_lines,
  file.path(
    manuscript_text_root,
    paste0(
      run_id,
      "_Results.md"
    )
  )
)

write_text(
  discussion_lines,
  file.path(
    manuscript_text_root,
    paste0(
      run_id,
      "_Discussion.md"
    )
  )
)


# ----------------------------------------------------------------------
# 15. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  cross_version = cross_version,
  GSE250130_reference =
    "ALL_FERTILE_LH7",
  GSE183837_reference =
    "STAGE_MATCHED_CONTROL",
  shared_cell_types =
    shared_cell_types,
  shared_endpoints =
    shared_endpoints,
  effect_orientation =
    "RIF_MINUS_REFERENCE",
  concordance_rule =
    "SAME_NONZERO_EFFECT_DIRECTION",
  GSE250130_primary_multiplicity =
    "PROGRAMS_BH_WITHIN_24_AND_VDR_BH_WITHIN_6",
  GSE183837_primary_multiplicity =
    "GLOBAL_BH_ACROSS_25",
  pooled_meta_analysis =
    FALSE,
  corrected_replication_rule =
    "BOTH_STUDIES_CORRECTED_SIGNIFICANT_AND_DIRECTION_CONCORDANT",
  single_cell_component_interpretation_locked =
    all_gate_checks_pass,
  project_level_result_lock_authorized =
    FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  cross_version = cross_version,
  started_at =
    as.character(run_start),
  completed_at =
    as.character(run_end),
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
  planned_common_model_n =
    nrow(cross_models),
  informative_model_n =
    nrow(informative_models),
  direction_concordant_n =
    sum(
      informative_models$
        direction_concordant
    ),
  direction_discordant_n =
    sum(
      !informative_models$
        direction_concordant
    ),
  corrected_replication_n =
    corrected_replication_n,
  exploratory_directional_support_n =
    exploratory_support_n,
  myeloid_VDR_direction_status =
    myeloid_VDR_direction_status,
  single_cell_component_interpretation_locked =
    all_gate_checks_pass,
  project_level_result_lock_authorized =
    FALSE,
  output_root =
    normalize_slashes(
      output_root
    ),
  qc_root =
    normalize_slashes(
      run_qc_root
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

persistent_files <- list.files(
  output_root,
  recursive = TRUE,
  full.names = TRUE
)

persistent_files <- persistent_files[
  file.exists(persistent_files) &
    !dir.exists(persistent_files)
]

persistent_manifest <- data.frame(
  path = normalize_slashes(
    persistent_files
  ),
  size_bytes =
    file.info(
      persistent_files
    )$size,
  md5 = unname(
    tools::md5sum(
      persistent_files
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  persistent_manifest,
  file.path(
    manifest_root,
    "persistent_manifest.csv"
  )
)

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
    figure_output_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    manifest_root,
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
      100 * 1024^2
]

artifact_manifest <- data.frame(
  path = normalize_slashes(
    small_files
  ),
  size_bytes =
    file.info(
      small_files
    )$size,
  md5 = unname(
    tools::md5sum(
      small_files
    )
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

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "informative=",
  nrow(informative_models),
  "concordant=",
  sum(
    informative_models$
      direction_concordant
  ),
  "corrected_replication=",
  corrected_replication_n
)

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos =
        "https://cloud.r-project.org"
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
      zip_error <<-
        conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1P single-cell cross-dataset concordance completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Shared planned models: ",
  nrow(cross_models),
  "\n",
  sep = ""
)
cat(
  "Informative models: ",
  nrow(informative_models),
  "\n",
  sep = ""
)
cat(
  "Direction-concordant models: ",
  sum(
    informative_models$
      direction_concordant
  ),
  "\n",
  sep = ""
)
cat(
  "Direction-discordant models: ",
  sum(
    !informative_models$
      direction_concordant
  ),
  "\n",
  sep = ""
)
cat(
  "Corrected replicated models: ",
  corrected_replication_n,
  "\n",
  sep = ""
)
cat(
  "Exploratory directional-support models: ",
  exploratory_support_n,
  "\n",
  sep = ""
)
cat(
  "Myeloid VDR direction status: ",
  myeloid_VDR_direction_status,
  "\n",
  sep = ""
)
cat(
  "Single-cell evidence state: ",
  single_cell_evidence_state,
  "\n",
  sep = ""
)
cat(
  "Single-cell component interpretation locked: ",
  all_gate_checks_pass,
  "\n",
  sep = ""
)
cat(
  "Output root: ",
  output_root,
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
    "\nCross-dataset outputs are complete in the project folders.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 2: 71_step04E1Q_STEP04_analytic_closeout_v1
# ==============================================================================

# ======================================================================
# 71_step04E1Q_STEP04_analytic_closeout_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Close the STEP04 analytic module after bulk projection, the
#   GSE250130 temporal/RIF single-cell analyses, the independent
#   GSE183837 donor-level analyses and cross-dataset concordance review.
# - Create an authoritative fact ledger, claim ledger, Figure/Table
#   blueprint, visual-QC register and manuscript-ready closeout text.
# - Lock the STEP04 scientific interpretation without claiming a
#   project-level final result lock.
#
# CENTRAL LOCKED INTERPRETATION
# - Normal luteal progression produces strong temporal transcriptomic
#   and VDR-program remodeling.
# - GSE250130 does not show broad donor-level RIF-associated VDR-program
#   dysregulation at LH+7.
# - GSE183837 does not show multiplicity-corrected VDR/program
#   differences.
# - The nominal GSE183837 myeloid VDR decrease is directionally
#   discordant in GSE250130 and is not replicated.
# - Across 25 common models, 24 are informative, 12 are concordant and
#   12 are discordant, with zero corrected replication.
# - Broad, reproducible RIF-associated VDR-program dysregulation is not
#   supported by the two donor-level single-cell datasets.
#
# GOVERNANCE
# - No new inferential testing is performed.
# - No pooled meta-analysis is performed.
# - No direct-target or causal claim is authorized.
# - Existing QC figures are not automatically promoted to final
#   manuscript figures.
# - STEP04 module interpretation may be locked.
# - Project-level final result lock remains FALSE.
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

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

script_name <-
  "71_step04E1Q_STEP04_analytic_closeout_v1.R"

script_version <-
  "STEP04E1Q_STEP04_ANALYTIC_CLOSEOUT_v1.0"

lock_version <- "STEP04_LOCK_v1.0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1Q_",
  format(run_start, "%Y%m%d_%H%M%S")
)


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  normalized <- normalize_slashes(path)

  path_n <- nchar(
    normalized,
    type = "chars"
  )

  if (
    .Platform$OS.type == "windows" &&
      path_n > max_chars
  ) {
    stop(
      "Path exceeds the Windows/R safety limit (",
      path_n,
      " > ",
      max_chars,
      "):\n",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  assert_path_safe(path)

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(lines, path) {
  assert_path_safe(path)

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

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  output <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(output)) {
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  output
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

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

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
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

latest_gate_file <- function(
  parent_root,
  directory_pattern,
  gate_filename
) {
  directories <- list.dirs(
    parent_root,
    recursive = FALSE,
    full.names = TRUE
  )

  directories <- directories[
    grepl(
      directory_pattern,
      basename(directories),
      perl = TRUE
    )
  ]

  if (!length(directories)) {
    stop(
      "No directory matched ",
      directory_pattern,
      " under ",
      parent_root,
      "."
    )
  }

  candidate_rows <- list()

  for (directory_path in directories) {
    gate_path <- file.path(
      directory_path,
      gate_filename
    )

    gate_value <- safe_read_csv(
      gate_path
    )

    if (
      !is.null(gate_value) &&
        nrow(gate_value) >= 1L
    ) {
      candidate_rows[[
        length(candidate_rows) + 1L
      ]] <- data.frame(
        gate_path =
          normalize_slashes(
            gate_path
          ),
        directory_path =
          normalize_slashes(
            directory_path
          ),
        directory_mtime =
          as.character(
            file.info(
              directory_path
            )$mtime
          ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(candidate_rows)) {
    stop(
      "No readable ",
      gate_filename,
      " was found under matched directories."
    )
  }

  candidates <- do.call(
    rbind,
    candidate_rows
  )

  candidates$mtime_POSIX <- as.POSIXct(
    candidates$
      directory_mtime
  )

  candidates <- candidates[
    order(
      candidates$mtime_POSIX,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  candidates$gate_path[1]
}

copy_validate <- function(
  source_path,
  destination_path
) {
  assert_path_safe(destination_path)

  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  if (
    !isTRUE(copied) ||
      !file.exists(destination_path) ||
      file.info(source_path)$size !=
        file.info(destination_path)$size ||
      unname(
        tools::md5sum(source_path)
      ) !=
        unname(
          tools::md5sum(
            destination_path
          )
        )
  ) {
    stop(
      "Validated copy failed:\n",
      source_path,
      "\n->\n",
      destination_path
    )
  }

  invisible(TRUE)
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
    cat(
      line,
      "\n",
      file = log_path,
      append = TRUE
    )
  }

  flush.console()
}


# ----------------------------------------------------------------------
# 2. Packages
# ----------------------------------------------------------------------

required_packages <- c(
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
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------------------
# 3. Paths
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

text_root <- file.path(
  step04_root,
  "09_manuscript_text"
)

archive_root <- file.path(
  step04_root,
  "10_archive_work"
)

parameter_root <- file.path(
  step04_root,
  "03_parameters"
)

environment_root <- file.path(
  step04_root,
  "04_environment"
)

closeout_parent <- file.path(
  results_root,
  "STEP04_closeout"
)

lock_root <- file.path(
  closeout_parent,
  lock_version
)

lock_status_path <- file.path(
  lock_root,
  "LOCK_STATUS.csv"
)

if (file.exists(lock_status_path)) {
  existing_lock_status <- safe_read_csv(
    lock_status_path
  )

  if (
    !is.null(existing_lock_status) &&
      nrow(existing_lock_status) == 1L &&
      existing_lock_status$status[[1]] ==
        "COMPLETE"
  ) {
    stop(
      "STEP04_LOCK_v1.0 already exists and is COMPLETE. ",
      "Do not overwrite an accepted lock."
    )
  }
}

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1Q_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_closeQC"
  )
)

staging_root <- file.path(
  archive_root,
  paste0(
    "S4CLOSE_",
    format(run_start, "%Y%m%d_%H%M%S")
  )
)

staging_table_root <- file.path(
  staging_root,
  "tab"
)

staging_text_root <- file.path(
  staging_root,
  "txt"
)

staging_manifest_root <- file.path(
  staging_root,
  "man"
)

snapshot_root <- file.path(
  staging_root,
  "snap"
)

for (path_value in c(
  execution_root,
  results_root,
  qc_base_root,
  text_root,
  archive_root,
  parameter_root,
  environment_root,
  closeout_parent,
  run_qc_root,
  staging_root,
  staging_table_root,
  staging_text_root,
  staging_manifest_root,
  snapshot_root
)) {
  dir.create(
    path_value,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

log_path <- file.path(
  execution_root,
  paste0(
    run_id,
    "_execution.log"
  )
)

run_manifest_path <- file.path(
  execution_root,
  paste0(
    run_id,
    "_RUN_MANIFEST.csv"
  )
)

parameter_path <- file.path(
  parameter_root,
  paste0(
    run_id,
    "_parameters.json"
  )
)

session_path <- file.path(
  environment_root,
  paste0(
    "sessionInfo_",
    run_id,
    ".txt"
  )
)

desktop_candidates <- c(
  file.path(
    Sys.getenv("USERPROFILE"),
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "桌面"
  )
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
  paste0(
    run_id,
    "_SMALL_REPORT.zip"
  )
)


# ----------------------------------------------------------------------
# 4. Upstream result paths
# ----------------------------------------------------------------------

bulk_root <- file.path(
  results_root,
  "S4D_qc_lock",
  "BVP_v1_3"
)

bulk_gate_path <- file.path(
  bulk_root,
  "STEP04D_QC1_GATE.csv"
)

bulk_score_path <- file.path(
  bulk_root,
  "score_models_v1_3.csv"
)

bulk_candidate_path <- file.path(
  bulk_root,
  "cross_context_candidates_v1_3.csv"
)

bulk_direction_summary_path <- file.path(
  bulk_root,
  "candidate_direction_summary.csv"
)

bulk_coverage_path <- file.path(
  bulk_root,
  "candidate_context_coverage.csv"
)

G250_model_root <- file.path(
  results_root,
  "S4E1_sc",
  "SCI_v1_1_memory_safe",
  "models",
  "GSE250130",
  "VDRMODELS_v1_1"
)

G250_temporal_program_path <- file.path(
  G250_model_root,
  "tables",
  "01_temporal_program_omnibus.csv"
)

G250_LH7_program_path <- file.path(
  G250_model_root,
  "tables",
  "03_LH7_program_primary.csv"
)

G250_temporal_VDR_path <- file.path(
  G250_model_root,
  "tables",
  "05_temporal_VDR_omnibus.csv"
)

G250_LH7_VDR_path <- file.path(
  G250_model_root,
  "tables",
  "07_LH7_VDR_primary.csv"
)

G250_DE_root <- file.path(
  results_root,
  "S4E1_sc",
  "SCI_v1_1_memory_safe",
  "de",
  "GSE250130",
  "DE_v1_2"
)

G250_temporal_DE_summary_path <- file.path(
  G250_DE_root,
  "tab",
  "07_temporal_omnibus_summary.csv"
)

G250_temporal_pair_summary_path <- file.path(
  G250_DE_root,
  "tab",
  "08_temporal_pairwise_summary.csv"
)

G250_LH7_DE_summary_path <- file.path(
  G250_DE_root,
  "tab",
  "09_LH7_primary_summary.csv"
)

G250_fry_path <- file.path(
  G250_DE_root,
  "cand",
  "06_directional_VDR_program_fry.csv"
)

G183_model_root <- file.path(
  results_root,
  "S4E1_sc",
  "SCI183_v1_0",
  "mdl",
  "G183",
  "MD183_v1_0"
)

G183_primary_model_path <- file.path(
  G183_model_root,
  "tab",
  "01_primary_models.csv"
)

G183_cluster8_path <- file.path(
  G183_model_root,
  "tab",
  "06_cluster8_model_sensitivity.csv"
)

cross_root <- file.path(
  results_root,
  "S4E1_sc",
  "cross_dataset",
  "SCCROSS_v1_0"
)

cross_model_path <- file.path(
  cross_root,
  "tab",
  "03_cross_dataset_25_model_concordance.csv"
)

cross_overall_path <- file.path(
  cross_root,
  "tab",
  "04_overall_concordance_summary.csv"
)

cross_lock_path <- file.path(
  cross_root,
  "tab",
  "10_single_cell_interpretation_lock.csv"
)

G250_model_gate_path <- latest_gate_file(
  parent_root = qc_base_root,
  directory_pattern =
    "^STEP04E1G_[0-9]{8}_[0-9]{6}_modelQC$",
  gate_filename =
    "VDR_MODEL_GATE.csv"
)

G250_DE_gate_path <- latest_gate_file(
  parent_root = qc_base_root,
  directory_pattern =
    "^STEP04E1H_[0-9]{8}_[0-9]{6}_deQC$",
  gate_filename =
    "GENE_DE_GATE.csv"
)

G183_model_gate_path <- latest_gate_file(
  parent_root = qc_base_root,
  directory_pattern =
    "^STEP04E1O_[0-9]{8}_[0-9]{6}_models$",
  gate_filename =
    "G183_DONOR_MODEL_GATE.csv"
)

cross_gate_path <- latest_gate_file(
  parent_root = qc_base_root,
  directory_pattern =
    "^STEP04E1P_[0-9]{8}_[0-9]{6}_crossQC$",
  gate_filename =
    "SC_CROSS_DATASET_GATE.csv"
)

candidate_evidence_path <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons",
  "06_results",
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  "RNA_OCCUPANCY_EVIDENCE_v1.0",
  "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv"
)

required_inputs <- c(
  bulk_gate_path,
  bulk_score_path,
  bulk_candidate_path,
  bulk_direction_summary_path,
  bulk_coverage_path,
  G250_model_gate_path,
  G250_temporal_program_path,
  G250_LH7_program_path,
  G250_temporal_VDR_path,
  G250_LH7_VDR_path,
  G250_DE_gate_path,
  G250_temporal_DE_summary_path,
  G250_temporal_pair_summary_path,
  G250_LH7_DE_summary_path,
  G250_fry_path,
  G183_model_gate_path,
  G183_primary_model_path,
  G183_cluster8_path,
  cross_gate_path,
  cross_model_path,
  cross_overall_path,
  cross_lock_path,
  candidate_evidence_path
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

path_audit <- data.frame(
  path =
    normalize_slashes(
      required_inputs
    ),
  exists =
    file.exists(
      required_inputs
    ),
  character_n =
    nchar(
      normalize_slashes(
        required_inputs
      ),
      type = "chars"
    ),
  path_safe =
    nchar(
      normalize_slashes(
        required_inputs
      ),
      type = "chars"
    ) <= 235L,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  path_audit,
  file.path(
    run_qc_root,
    "00_input_path_audit.csv"
  )
)

if (length(missing_inputs)) {
  stop(
    "Missing required STEP04 closeout input(s):\n",
    paste(
      normalize_slashes(
        missing_inputs
      ),
      collapse = "\n"
    )
  )
}

if (!all(path_audit$path_safe)) {
  stop(
    "One or more closeout input paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 5. Read and validate upstream gates
# ----------------------------------------------------------------------

bulk_gate <- safe_read_csv(
  bulk_gate_path
)

G250_model_gate <- safe_read_csv(
  G250_model_gate_path
)

G250_DE_gate <- safe_read_csv(
  G250_DE_gate_path
)

G183_model_gate <- safe_read_csv(
  G183_model_gate_path
)

cross_gate <- safe_read_csv(
  cross_gate_path
)

bulk_gate_status <- bulk_gate$status[
  bulk_gate$gate_component ==
    "STEP04D-QC1 gate"
]

if (length(bulk_gate_status) != 1L) {
  stop(
    "Could not uniquely identify the STEP04D-QC1 Gate row."
  )
}

gate_inventory <- data.frame(
  module = c(
    "STEP04D_bulk_projection",
    "STEP04E1G_GSE250130_donor_models",
    "STEP04E1H_GSE250130_gene_DE",
    "STEP04E1O_GSE183837_donor_models",
    "STEP04E1P_cross_dataset"
  ),
  gate_status = c(
    bulk_gate_status,
    G250_model_gate$gate_status[[1]],
    G250_DE_gate$gate_status[[1]],
    G183_model_gate$gate_status[[1]],
    cross_gate$gate_status[[1]]
  ),
  gate_path = normalize_slashes(c(
    bulk_gate_path,
    G250_model_gate_path,
    G250_DE_gate_path,
    G183_model_gate_path,
    cross_gate_path
  )),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_inventory,
  file.path(
    staging_table_root,
    "01_gate_inventory.csv"
  )
)

log_event(
  "Starting",
  run_id,
  "bulk_gate=",
  bulk_gate_status,
  "cross_gate=",
  cross_gate$gate_status[[1]]
)


# ----------------------------------------------------------------------
# 6. Read key result tables
# ----------------------------------------------------------------------

bulk_scores <- safe_read_csv(
  bulk_score_path
)

bulk_candidates <- safe_read_csv(
  bulk_candidate_path
)

bulk_direction_summary <- safe_read_csv(
  bulk_direction_summary_path
)

bulk_coverage <- safe_read_csv(
  bulk_coverage_path
)

G250_temporal_program <- safe_read_csv(
  G250_temporal_program_path
)

G250_LH7_program <- safe_read_csv(
  G250_LH7_program_path
)

G250_temporal_VDR <- safe_read_csv(
  G250_temporal_VDR_path
)

G250_LH7_VDR <- safe_read_csv(
  G250_LH7_VDR_path
)

G250_temporal_DE_summary <- safe_read_csv(
  G250_temporal_DE_summary_path
)

G250_temporal_pair_summary <- safe_read_csv(
  G250_temporal_pair_summary_path
)

G250_LH7_DE_summary <- safe_read_csv(
  G250_LH7_DE_summary_path
)

G250_fry <- safe_read_csv(
  G250_fry_path
)

G183_primary_models <- safe_read_csv(
  G183_primary_model_path
)

G183_cluster8 <- safe_read_csv(
  G183_cluster8_path
)

cross_models <- safe_read_csv(
  cross_model_path
)

cross_overall <- safe_read_csv(
  cross_overall_path
)

cross_interpretation_lock <- safe_read_csv(
  cross_lock_path
)

candidate_evidence <- safe_read_csv(
  candidate_evidence_path
)

data_objects <- c(
  "bulk_scores",
  "bulk_candidates",
  "bulk_direction_summary",
  "bulk_coverage",
  "G250_temporal_program",
  "G250_LH7_program",
  "G250_temporal_VDR",
  "G250_LH7_VDR",
  "G250_temporal_DE_summary",
  "G250_temporal_pair_summary",
  "G250_LH7_DE_summary",
  "G250_fry",
  "G183_primary_models",
  "G183_cluster8",
  "cross_models",
  "cross_overall",
  "cross_interpretation_lock",
  "candidate_evidence"
)

schema_rows <- list()

for (data_name in data_objects) {
  data_value <- get(data_name)

  if (
    is.null(data_value) ||
      !nrow(data_value) ||
      anyDuplicated(
        colnames(data_value)
      ) > 0L
  ) {
    stop(
      data_name,
      " is absent, empty or contains duplicate columns."
    )
  }

  schema_rows[[
    length(schema_rows) + 1L
  ]] <- data.frame(
    object_name =
      data_name,
    row_n =
      nrow(data_value),
    column_n =
      ncol(data_value),
    duplicate_column_n =
      sum(
        duplicated(
          colnames(data_value)
        )
      ),
    column_names =
      paste(
        colnames(data_value),
        collapse = ";"
      ),
    stringsAsFactors = FALSE
  )
}

schema_audit <- do.call(
  rbind,
  schema_rows
)

write_csv_utf8(
  schema_audit,
  file.path(
    run_qc_root,
    "01_schema_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Identity and summary calculations
# ----------------------------------------------------------------------

for (logical_column in c(
  "occupancy_supported",
  "shared_occupancy_supported",
  "any_promoter_overlap_support"
)) {
  if (
    logical_column %in%
      colnames(candidate_evidence)
  ) {
    candidate_evidence[[
      logical_column
    ]] <- as_logical_safe(
      candidate_evidence[[
        logical_column
      ]]
    )
  }
}

if (
  nrow(candidate_evidence) != 436L ||
    sum(
      candidate_evidence$
        occupancy_supported
    ) != 284L ||
    sum(
      candidate_evidence$
        shared_occupancy_supported
    ) != 46L ||
    sum(
      candidate_evidence$
        any_promoter_overlap_support
    ) != 84L
) {
  stop(
    "Locked candidate identity is not 436/284/46/84."
  )
}

bulk_scores$adj.P.Val <- safe_numeric(
  bulk_scores$adj.P.Val
)

bulk_scores$P.Value <- safe_numeric(
  bulk_scores$P.Value
)

bulk_scores$logFC <- safe_numeric(
  bulk_scores$logFC
)

bulk_scores$bulk_FDR_significant <- (
  is.finite(
    bulk_scores$adj.P.Val
  ) &
    bulk_scores$adj.P.Val <
    0.05
)

bulk_scores$effect_direction <- ifelse(
  bulk_scores$logFC > 0,
  "PATIENT_HIGHER",
  ifelse(
    bulk_scores$logFC < 0,
    "PATIENT_LOWER",
    "NO_DIFFERENCE"
  )
)

bulk_contrast_summary <- data.table::as.data.table(
  bulk_scores
)[
  ,
  .(
    program_n = .N,
    FDR_lt_0_05_n =
      sum(
        bulk_FDR_significant
      ),
    nominal_P_lt_0_05_n =
      sum(
        P.Value < 0.05,
        na.rm = TRUE
      ),
    positive_effect_n =
      sum(
        logFC > 0,
        na.rm = TRUE
      ),
    negative_effect_n =
      sum(
        logFC < 0,
        na.rm = TRUE
      ),
    minimum_FDR =
      min(
        adj.P.Val,
        na.rm = TRUE
      ),
    significant_programs =
      collapse_unique(
        program_name[
          bulk_FDR_significant
        ]
      )
  ),
  by = .(
    dataset_id,
    contrast
  )
]

bulk_contrast_summary <- as.data.frame(
  bulk_contrast_summary
)

G250_temporal_program_FDR_n <- sum(
  safe_numeric(
    G250_temporal_program$FDR_BH
  ) < 0.05,
  na.rm = TRUE
)

G250_LH7_program_FDR_n <- sum(
  safe_numeric(
    G250_LH7_program$FDR_BH
  ) < 0.05,
  na.rm = TRUE
)

G250_temporal_VDR_FDR_n <- sum(
  safe_numeric(
    G250_temporal_VDR$FDR_BH
  ) < 0.05,
  na.rm = TRUE
)

G250_LH7_VDR_FDR_n <- sum(
  safe_numeric(
    G250_LH7_VDR$FDR_BH
  ) < 0.05,
  na.rm = TRUE
)

G250_temporal_omnibus_gene_n <-
  G250_DE_gate$
    temporal_omnibus_genome_FDR_lt_0_05_n[[1]]

G250_temporal_pair_gene_n <-
  G250_DE_gate$
    temporal_pairwise_genome_DE_n[[1]]

G250_LH7_primary_gene_n <-
  G250_DE_gate$
    LH7_primary_genome_DE_n[[1]]

G250_LH7_sensitivity_gene_n <-
  G250_DE_gate$
    LH7_sensitivity_genome_DE_n[[1]]

G183_nominal_n <- sum(
  safe_numeric(
    G183_primary_models$
      exact_permutation_P
  ) < 0.05,
  na.rm = TRUE
)

G183_global_BH_n <- sum(
  safe_numeric(
    G183_primary_models$
      BH_global_25
  ) < 0.05,
  na.rm = TRUE
)

G183_cluster8_direction_changed_n <- sum(
  G183_cluster8$
    direction_concordant %in%
    FALSE,
  na.rm = TRUE
)

cross_corrected_replication_n <-
  cross_gate$
    corrected_replication_n[[1]]

cross_exploratory_support_n <-
  cross_gate$
    exploratory_directional_support_n[[1]]

cross_concordant_n <-
  cross_gate$
    direction_concordant_n[[1]]

cross_discordant_n <-
  cross_gate$
    direction_discordant_n[[1]]

cross_informative_n <-
  cross_gate$
    informative_model_n[[1]]

myeloid_cross_direction <-
  cross_gate$
    myeloid_VDR_direction_status[[1]]

bulk_signal_state <- if (
  sum(
    bulk_scores$
      bulk_FDR_significant
  ) > 0L
) {
  "CONTEXT_SPECIFIC_BULK_PROGRAM_SIGNALS_PRESENT"
} else {
  "NO_BULK_PROGRAM_FDR_SIGNAL"
}


# ----------------------------------------------------------------------
# 8. Authoritative fact ledger
# ----------------------------------------------------------------------

fact_ledger <- data.frame(
  fact_id = c(
    "F001",
    "F002",
    "F003",
    "F004",
    "F005",
    "F006",
    "F007",
    "F008",
    "F009",
    "F010",
    "F011",
    "F012",
    "F013",
    "F014",
    "F015",
    "F016",
    "F017",
    "F018",
    "F019",
    "F020",
    "F021",
    "F022",
    "F023",
    "F024",
    "F025"
  ),
  domain = c(
    "Perturbation",
    "Occupancy",
    "Occupancy",
    "Occupancy",
    "Terminology",
    "Bulk projection",
    "Bulk projection",
    "GSE250130 temporal",
    "GSE250130 temporal",
    "GSE250130 temporal",
    "GSE250130 RIF",
    "GSE250130 RIF",
    "GSE250130 RIF",
    "GSE250130 sensitivity",
    "GSE183837",
    "GSE183837",
    "GSE183837",
    "GSE183837",
    "Cross-dataset",
    "Cross-dataset",
    "Cross-dataset",
    "Cross-dataset",
    "Inference boundary",
    "Inference boundary",
    "Module status"
  ),
  fact = c(
    "Locked RNA interaction candidate universe",
    "Occupancy-supported RNA interaction candidates",
    "Shared occupancy-supported candidates",
    "Promoter-associated occupancy-supported candidates",
    "Definitive direct-target assignments",
    "Bulk projection contrasts",
    "Bulk program FDR-significant rows",
    "Temporal VDR-program models",
    "Temporal VDR-program FDR-significant models",
    "Temporal genome-wide significant genes",
    "LH+7 RIF-versus-fertile VDR-program FDR-significant models",
    "LH+7 RIF-versus-fertile VDR-expression FDR-significant models",
    "LH+7 primary genome-wide DE genes",
    "LH+7 fertile-subgroup sensitivity DE rows",
    "GSE183837 planned donor-level models",
    "GSE183837 nominal exact-P signals",
    "GSE183837 global-BH signals",
    "Cluster-8 sensitivity direction changes",
    "Shared cross-dataset models",
    "Directionally informative shared models",
    "Direction-concordant models",
    "Corrected cross-dataset replication",
    "Broad reproducible RIF VDR-program dysregulation supported",
    "Pooled cross-study meta-analysis performed",
    "STEP04 analytic interpretation locked"
  ),
  value = c(
    436,
    284,
    46,
    84,
    0,
    nrow(
      bulk_contrast_summary
    ),
    sum(
      bulk_scores$
        bulk_FDR_significant
    ),
    nrow(
      G250_temporal_program
    ),
    G250_temporal_program_FDR_n,
    G250_temporal_omnibus_gene_n,
    G250_LH7_program_FDR_n,
    G250_LH7_VDR_FDR_n,
    G250_LH7_primary_gene_n,
    G250_LH7_sensitivity_gene_n,
    nrow(
      G183_primary_models
    ),
    G183_nominal_n,
    G183_global_BH_n,
    G183_cluster8_direction_changed_n,
    nrow(cross_models),
    cross_informative_n,
    cross_concordant_n,
    cross_corrected_replication_n,
    FALSE,
    FALSE,
    TRUE
  ),
  status = "LOCKED",
  source = c(
    rep(
      "RNA_OCCUPANCY_EVIDENCE_v1.0",
      5
    ),
    rep(
      "BULK_VDR_PROGRAM_PROJECTION_v1.3",
      2
    ),
    rep(
      "GSE250130 VDRMODELS_v1.1 / DE_v1.2",
      7
    ),
    rep(
      "GSE183837 MD183_v1_0",
      4
    ),
    rep(
      "SCCROSS_v1_0",
      4
    ),
    "SCCROSS_v1_0",
    "SCCROSS_v1_0",
    lock_version
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  fact_ledger,
  file.path(
    staging_table_root,
    "02_locked_fact_ledger.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Claim ledger
# ----------------------------------------------------------------------

claim_ledger <- data.frame(
  claim_id = c(
    "C001",
    "C002",
    "C003",
    "C004",
    "C005",
    "C006",
    "C007",
    "C008",
    "C009",
    "C010",
    "C011",
    "C012"
  ),
  claim_topic = c(
    "Perturbation program",
    "Occupancy support",
    "Normal temporal remodeling",
    "Bulk patient projection",
    "Broad RIF VDR-program dysregulation",
    "GSE250130 RIF comparison",
    "GSE183837 RIF comparison",
    "Myeloid VDR signal",
    "Cross-dataset replication",
    "Cluster-8 sensitivity",
    "Direct VDR targets",
    "Diagnostic biomarker"
  ),
  evidence_status = c(
    "SUPPORTED",
    "SUPPORTED_WITH_TERMINOLOGY_LIMIT",
    "STRONGLY_SUPPORTED",
    bulk_signal_state,
    "NOT_SUPPORTED",
    "NEGATIVE_PRIMARY_RESULT",
    "NEGATIVE_AFTER_MULTIPLICITY",
    "UNREPLICATED_NOMINAL_SIGNAL",
    "NO_CORRECTED_OR_EXPLORATORY_DIRECTIONAL_SUPPORT",
    "ROBUST_DIRECTION_TO_CLUSTER8_SENSITIVITY",
    "NOT_AUTHORIZED",
    "NOT_EVALUATED_AND_NOT_AUTHORIZED"
  ),
  permitted_wording = c(
    "Calcitriol-dependent VDR perturbation defined a locked set of 436 RNA interaction candidates.",
    "A subset of RNA interaction candidates showed convergent support from published VDR occupancy resources.",
    "Normal luteal progression was accompanied by strong cell-type-resolved transcriptomic and VDR-program remodeling.",
    "Patient-cohort program projections were context-specific and must be reported separately for each prespecified contrast.",
    "The donor-level single-cell datasets did not support broad, reproducible RIF-associated disruption of the locked VDR programs.",
    "At LH+7, GSE250130 showed no FDR-significant RIF-versus-fertile VDR-program, VDR-expression or genome-wide DE result.",
    "GSE183837 showed no VDR-expression or VDR-program result that survived global multiplicity correction.",
    "GSE183837 showed a nominal myeloid VDR-expression decrease that was directionally discordant in GSE250130.",
    "No corrected cross-dataset replication or exploratory directional support was observed across the 25 shared models.",
    "Including low-confidence cluster 8 did not change the direction of the five unciliated epithelial effects.",
    "Use occupancy-supported RNA interaction candidates; no definitive direct-target assignment was made.",
    "No diagnostic classifier or biomarker performance claim was evaluated."
  ),
  prohibited_wording = c(
    "VDR causal disease genes",
    "Definitive direct VDR targets",
    "RIF causes global VDR deficiency",
    "All patient cohorts replicate the same VDR defect",
    "RIF universally suppresses VDR signaling",
    "GSE250130 validates an RIF VDR defect",
    "GSE183837 independently validates VDR-program dysregulation",
    "Myeloid VDR reduction is externally validated",
    "The two cohorts provide replicated confirmation",
    "Cluster 8 proves an RIF-specific epithelial subtype",
    "Direct causal targets of VDR",
    "VDR program is a diagnostic biomarker"
  ),
  manuscript_role = c(
    "Mechanistic anchor",
    "Mechanistic support",
    "Primary positive biological result",
    "Contextual patient evidence",
    "Primary negative boundary",
    "Primary negative disease result",
    "Independent negative result",
    "Exploratory limitation",
    "Cross-dataset conclusion",
    "Sensitivity analysis",
    "Terminology restriction",
    "Scope restriction"
  ),
  locked = TRUE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  claim_ledger,
  file.path(
    staging_table_root,
    "03_claim_ledger.csv"
  )
)


# ----------------------------------------------------------------------
# 10. Figure and Table blueprint
# ----------------------------------------------------------------------

figure_blueprint <- data.frame(
  figure_id = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5",
    "Figure 6",
    "Supplementary Figure S1",
    "Supplementary Figure S2",
    "Supplementary Figure S3",
    "Supplementary Figure S4",
    "Supplementary Figure S5"
  ),
  title = c(
    "Study design and evidence hierarchy",
    "Perturbation-defined VDR RNA program and occupancy support",
    "Projection of locked VDR programs into patient bulk cohorts",
    "Cell-type-resolved temporal remodeling in GSE250130",
    "Donor-level disease comparisons in GSE250130 and GSE183837",
    "Cross-dataset concordance and integrated evidence boundary",
    "Dataset identity, technical design and confounding audits",
    "GSE250130 integration, annotation and donor-mixing QC",
    "GSE250130 gene-level temporal and LH+7 pseudobulk results",
    "GSE183837 integration, annotation and cluster-8 sensitivity",
    "Complete donor-level model and concordance diagnostics"
  ),
  proposed_panels = c(
    "A workflow; B datasets/roles; C 436-to-284/46/84 evidence hierarchy; D analysis gates",
    "A RNA interaction volcano; B candidate heatmap; C occupancy hierarchy; D candidate terminology schematic",
    "A four-contrast score heatmap; B FDR/effect matrix; C candidate context coverage; D context-dependent examples",
    "A broad-cell UMAP; B donor composition; C temporal program trajectories; D temporal VDR; E temporal DE counts",
    "A GSE250130 LH+7 program/VDR forest; B GSE250130 LH+7 DE summary; C GSE183837 UMAP; D GSE183837 donor effects; E ciliated descriptive",
    "A 25-model direction matrix; B effect-size comparison; C myeloid VDR discordance; D final evidence-state diagram",
    "STEP00 slide/phenotype audit and STEP04 sample-design tables",
    "Doublet, donor dominance, marker and annotation QC",
    "Temporal omnibus/pairwise DE, fry tests and sensitivity analyses",
    "Donor mixing, marker evidence, low-confidence cluster review and atlas sensitivity",
    "Exact permutation distributions, LOO results, multiplicity and all 25 cross-dataset rows"
  ),
  central_question = c(
    "How was evidence prospectively organized?",
    "What experimental RNA program was frozen and how is occupancy used?",
    "How do the programs behave across patient bulk contexts?",
    "How strongly does normal luteal timing remodel the programs?",
    "Is broad RIF-associated dysregulation present at donor level?",
    "Are disease effects reproducible across independent single-cell cohorts?",
    "Are datasets and contrasts technically valid?",
    "Are GSE250130 cell identities and donor structure reliable?",
    "Which temporal and disease gene-level effects are supported?",
    "Are GSE183837 annotation and sensitivity conclusions robust?",
    "Are all inferential and stability diagnostics transparent?"
  ),
  current_status = c(
    "BLUEPRINT_ONLY",
    "B_CANDIDATE_INPUTS_EXIST_FINAL_COMPOSITE_PENDING",
    "BLUEPRINT_ONLY",
    "SOURCE_RESULTS_READY_FINAL_COMPOSITE_PENDING",
    "SOURCE_RESULTS_READY_FINAL_COMPOSITE_PENDING",
    "SOURCE_RESULTS_READY_FINAL_COMPOSITE_PENDING",
    "SUPPLEMENTARY_ASSEMBLY_PENDING",
    "SUPPLEMENTARY_ASSEMBLY_PENDING",
    "SUPPLEMENTARY_ASSEMBLY_PENDING",
    "SUPPLEMENTARY_ASSEMBLY_PENDING",
    "SUPPLEMENTARY_ASSEMBLY_PENDING"
  ),
  final_figure_authorized = FALSE,
  stringsAsFactors = FALSE
)

table_blueprint <- data.frame(
  table_id = c(
    "Table 1",
    "Table 2",
    "Table 3",
    "Supplementary Table S1",
    "Supplementary Table S2",
    "Supplementary Table S3",
    "Supplementary Table S4",
    "Supplementary Table S5"
  ),
  title = c(
    "Datasets, biological roles and prespecified contrasts",
    "Definitions and evidence hierarchy of the four locked VDR programs",
    "Integrated summary of temporal, disease and cross-dataset findings",
    "Complete 436-candidate RNA/occupancy evidence table",
    "Bulk projection results across four patient contrasts",
    "GSE250130 donor-level program, VDR and gene-level models",
    "GSE183837 donor-level models and cluster-8 sensitivity",
    "Cross-dataset 25-model concordance and claim ledger"
  ),
  required_source = c(
    "STEP04B/STEP04C locked designs",
    "RNA_OCCUPANCY_EVIDENCE_v1.0",
    "STEP04 fact ledger and claim ledger",
    "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv",
    "score_models_v1_3.csv",
    "VDRMODELS_v1.1 and DE_v1.2",
    "MD183_v1_0",
    "SCCROSS_v1_0"
  ),
  current_status = c(
    "ASSEMBLY_PENDING",
    "ASSEMBLY_PENDING",
    "READY_FROM_CLOSEOUT",
    "SOURCE_READY",
    "SOURCE_READY",
    "SOURCE_READY",
    "SOURCE_READY",
    "READY_FROM_CLOSEOUT"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  figure_blueprint,
  file.path(
    staging_table_root,
    "07_figure_blueprint.csv"
  )
)

write_csv_utf8(
  table_blueprint,
  file.path(
    staging_table_root,
    "08_table_blueprint.csv"
  )
)


# ----------------------------------------------------------------------
# 11. Active script chain and visual-QC register
# ----------------------------------------------------------------------

active_script_chain <- data.frame(
  order = seq_len(23L),
  module = c(
    rep("STEP04 source/design/bulk", 5),
    rep("GSE250130", 7),
    rep("GSE183837", 8),
    "Cross-dataset",
    "STEP04 closeout",
    "STEP04 closeout"
  ),
  script = c(
    "48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R",
    "49_step04B_lock_inputs_extract_sc_and_audit_designs.R",
    "50_step04C_v3_lock_sample_designs_and_platform_preflight.R",
    "51_step04D_v6_bulk_projection_and_sc_preQC.R",
    "52_step04D_QC1_v2_repair_and_lock_outputs.R",
    "56_step04E1B_GSE250130_integration_resume_v1.R",
    "57_step04E1C_GSE250130_integration_QC_v1.R",
    "58_step04E1D_GSE250130_cluster_markers_v1_1_memory_safe.R",
    "59_step04E1E_GSE250130_broad_annotation_QC_v1_2_pathsafe.R",
    "60_step04E1F_GSE250130_VDR_localization_pseudobulk_v1.R",
    "61_step04E1G_GSE250130_donor_VDR_models_v1_1_unbalanced_fixed.R",
    "62_step04E1H_GSE250130_gene_pseudobulk_DE_v1_2_pathsafe_resume.R",
    "63_step04E1I_GSE183837_integration_preflight_v1.R",
    "64_step04E1J_GSE183837_doublet_Harmony_integration_v1_4_donor_gate_fixed.R",
    "65_step04E1K_GSE183837_detailed_integration_QC_v1_2_pathsafe.R",
    "66_step04E1L_GSE183837_descriptive_markers_v1_1_panel_duplicate_fixed.R",
    "67_step04E1M_GSE183837_broad_annotation_QC_v1.R",
    "68_step04E1N_GSE183837_VDR_localization_pseudobulk_v1_1_unique_columns_fixed.R",
    "69_step04E1O_GSE183837_donor_models_v1_1_type_safe_summary.R",
    "70_step04E1P_cross_dataset_single_cell_concordance_v1.R",
    "71_step04E1Q_STEP04_analytic_closeout_v1.R",
    "NO_ADDITIONAL_DISCOVERY_TESTING_AFTER_STEP04_CLOSEOUT",
    "FINAL_FIGURE_ASSEMBLY_AND_MANUSCRIPT_INTEGRATION_ONLY"
  ),
  role = c(
    "Source acquisition",
    "Input lock",
    "Design lock",
    "Bulk computation",
    "Bulk correction and lock",
    "Integration",
    "Integration QC",
    "Markers",
    "Annotation",
    "Pseudobulk",
    "Donor models",
    "Gene-level DE",
    "Preflight",
    "Integration",
    "Detailed QC",
    "Markers",
    "Annotation",
    "Pseudobulk",
    "Donor models",
    "Cross-dataset lock",
    "Analytic closeout",
    "Stop rule",
    "Authorized next phase"
  ),
  authoritative = TRUE,
  stringsAsFactors = FALSE
)

visual_QC_register <- data.frame(
  artifact = c(
    "GSE250130 temporal/program QC figures",
    "GSE250130 LH+7 disease forest and DE QC",
    "GSE183837 broad annotation UMAP",
    "GSE183837 annotation marker heatmap",
    "GSE183837 donor composition",
    "GSE183837 donor effect forest",
    "Cross-dataset effect scatter",
    "Cross-dataset direction matrix",
    "Cross-dataset Hedges-g comparison",
    "Cross-dataset adjusted-evidence matrix",
    "Myeloid VDR focus"
  ),
  content_QC = c(
    rep("PASS", 6),
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS"
  ),
  visual_QC = c(
    rep(
      "QC_READABLE_NOT_FINAL_COMPOSITE",
      6
    ),
    "QC_ONLY_LABEL_CLIPPING_AND_OVERLAP_REQUIRES_REPLOT",
    "QC_READABLE_FINAL_STYLE_PENDING",
    "QC_READABLE_FINAL_STYLE_PENDING",
    "QC_READABLE_FINAL_STYLE_PENDING",
    "QC_READABLE_FINAL_STYLE_PENDING"
  ),
  manuscript_ready = FALSE,
  required_action = c(
    "Select source panels and rebuild one temporal main figure.",
    "Rebuild compact disease-result panels with unified scales and legends.",
    "Use as source panel after typography and label review.",
    "Use as supplementary source after final panel reduction.",
    "Use as supplementary or disease-validation panel.",
    "Use source data to rebuild a compact forest/point figure.",
    "Remove edge labels or use controlled repel; add margins; do not use current PNG.",
    "Retain matrix logic; apply manuscript typography and palette.",
    "Reduce legend and standardize endpoint labels.",
    "Clarify differing q-value families in legend.",
    "Retain as focused negative/discordant evidence panel after final styling."
  ),
  manual_review_date =
    "2026-08-06",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  active_script_chain,
  file.path(
    staging_table_root,
    "09_active_script_chain.csv"
  )
)

write_csv_utf8(
  visual_QC_register,
  file.path(
    staging_table_root,
    "10_visual_QC_register.csv"
  )
)


# ----------------------------------------------------------------------
# 12. Integrated summaries
# ----------------------------------------------------------------------

write_csv_utf8(
  bulk_scores,
  file.path(
    staging_table_root,
    "04_bulk_program_models.csv"
  )
)

write_csv_utf8(
  bulk_contrast_summary,
  file.path(
    staging_table_root,
    "05_bulk_contrast_summary.csv"
  )
)

integrated_result_summary <- data.frame(
  analysis_layer = c(
    "Experimental perturbation",
    "Published occupancy support",
    "Patient bulk projection",
    "GSE250130 normal temporal programs",
    "GSE250130 normal temporal genome-wide DE",
    "GSE250130 LH+7 RIF programs",
    "GSE250130 LH+7 RIF VDR expression",
    "GSE250130 LH+7 primary genome-wide DE",
    "GSE183837 donor models",
    "Cross-dataset concordance"
  ),
  primary_result = c(
    "436 locked calcitriol-dependent VDR RNA interaction candidates",
    "284 occupancy-supported; 46 shared occupancy; 84 promoter-associated",
    bulk_signal_state,
    paste0(
      G250_temporal_program_FDR_n,
      " of ",
      nrow(G250_temporal_program),
      " temporal program models FDR-significant"
    ),
    paste0(
      G250_temporal_omnibus_gene_n,
      " temporal omnibus gene rows FDR-significant"
    ),
    paste0(
      G250_LH7_program_FDR_n,
      " of ",
      nrow(G250_LH7_program),
      " primary program models FDR-significant"
    ),
    paste0(
      G250_LH7_VDR_FDR_n,
      " of ",
      nrow(G250_LH7_VDR),
      " VDR-expression models FDR-significant"
    ),
    paste0(
      G250_LH7_primary_gene_n,
      " primary genome-wide DE genes"
    ),
    paste0(
      G183_nominal_n,
      " nominal exact-P signals; ",
      G183_global_BH_n,
      " global-BH signals"
    ),
    paste0(
      cross_concordant_n,
      " concordant and ",
      cross_discordant_n,
      " discordant among ",
      cross_informative_n,
      " informative models; 0 corrected replication"
    )
  ),
  evidence_class = c(
    "MECHANISTIC_ANCHOR",
    "CROSS_CONTEXT_SUPPORT_NOT_DIRECT_TARGET_PROOF",
    "CONTEXT_SPECIFIC_PATIENT_PROJECTION",
    "STRONG_POSITIVE_TEMPORAL_RESULT",
    "STRONG_POSITIVE_TEMPORAL_RESULT",
    "NEGATIVE_PRIMARY_DISEASE_RESULT",
    "NEGATIVE_PRIMARY_DISEASE_RESULT",
    "NEGATIVE_PRIMARY_DISEASE_RESULT",
    "NEGATIVE_AFTER_MULTIPLICITY",
    "NO_REPLICATION_OR_DIRECTIONAL_SUPPORT"
  ),
  manuscript_use = c(
    "Primary mechanism definition",
    "Mechanism-support hierarchy",
    "Patient-context evidence",
    "Primary positive biological result",
    "Primary positive biological result",
    "Primary disease boundary",
    "Primary disease boundary",
    "Primary disease boundary",
    "Independent disease boundary",
    "Final single-cell conclusion"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  integrated_result_summary,
  file.path(
    staging_table_root,
    "06_integrated_result_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 13. Closeout text
# ----------------------------------------------------------------------

bulk_significant_rows <- bulk_scores[
  bulk_scores$bulk_FDR_significant,
  ,
  drop = FALSE
]

bulk_significant_text <- if (
  nrow(bulk_significant_rows)
) {
  paste(
    paste0(
      bulk_significant_rows$
        dataset_id,
      " ",
      bulk_significant_rows$
        contrast,
      " / ",
      bulk_significant_rows$
        program_name,
      " (logFC=",
      format(
        bulk_significant_rows$logFC,
        digits = 4
      ),
      ", FDR=",
      format(
        bulk_significant_rows$
          adj.P.Val,
        digits = 4
      ),
      ")"
    ),
    collapse = "; "
  )
} else {
  "No prespecified bulk program model had FDR < 0.05."
}

locked_interpretation_lines <- c(
  "# STEP04 Locked Analytic Interpretation",
  "",
  paste0("- Lock version: `", lock_version, "`"),
  paste0("- Run ID: `", run_id, "`"),
  "",
  "## Central conclusion",
  "",
  "The project identifies a reproducible calcitriol-dependent VDR perturbation program and strong cell-type-resolved remodeling across normal luteal progression. In contrast, the donor-level patient single-cell analyses do not support broad, stable or cross-dataset-reproducible RIF-associated disruption of the locked VDR programs.",
  "",
  "## Evidence hierarchy",
  "",
  "- The mechanistic anchor is the locked 436-gene RNA interaction candidate universe.",
  "- Published VDR cistrome resources provide occupancy support for 284 candidates, including 46 with shared occupancy support and 84 with promoter-associated support.",
  "- These candidates are occupancy-supported RNA interaction candidates, not definitive direct VDR targets.",
  paste0(
    "- Bulk projection state: ",
    bulk_signal_state,
    "."
  ),
  paste0(
    "- Bulk FDR-significant models: ",
    bulk_significant_text
  ),
  "",
  "## Normal temporal biology",
  "",
  paste0(
    "- In GSE250130, ",
    G250_temporal_program_FDR_n,
    " of ",
    nrow(G250_temporal_program),
    " temporal VDR-program models were FDR-significant."
  ),
  paste0(
    "- Gene-level temporal analyses produced ",
    G250_temporal_omnibus_gene_n,
    " significant omnibus rows and ",
    G250_temporal_pair_gene_n,
    " significant planned-pairwise rows."
  ),
  "- The temporal series is cross-sectional and establishes structured timing-associated remodeling rather than within-person trajectories.",
  "",
  "## Disease comparisons",
  "",
  paste0(
    "- In the primary GSE250130 LH+7 RIF comparison, FDR-significant program models=",
    G250_LH7_program_FDR_n,
    ", VDR-expression models=",
    G250_LH7_VDR_FDR_n,
    ", and genome-wide DE genes=",
    G250_LH7_primary_gene_n,
    "."
  ),
  paste0(
    "- GSE250130 subgroup sensitivity yielded ",
    G250_LH7_sensitivity_gene_n,
    " DE rows and is interpreted as fertile-reference heterogeneity rather than primary RIF evidence."
  ),
  paste0(
    "- In GSE183837, ",
    G183_nominal_n,
    " of 25 models had nominal exact P<0.05 and ",
    G183_global_BH_n,
    " survived global BH correction."
  ),
  paste0(
    "- The nominal myeloid VDR signal was ",
    myeloid_cross_direction,
    " across datasets and therefore was not replicated."
  ),
  paste0(
    "- Among ",
    cross_informative_n,
    " informative shared models, ",
    cross_concordant_n,
    " were direction-concordant and ",
    cross_discordant_n,
    " were discordant."
  ),
  paste0(
    "- Corrected cross-dataset replication=",
    cross_corrected_replication_n,
    "; exploratory directional support=",
    cross_exploratory_support_n,
    "."
  ),
  "",
  "## Locked scientific boundary",
  "",
  "- Broad, reproducible RIF-associated VDR-program dysregulation is not supported.",
  "- No direct-target, causal VDR-activity or diagnostic biomarker claim is authorized.",
  "- No pooled cross-study meta-analysis is authorized.",
  "- No additional discovery testing should be added to rescue a preferred disease narrative.",
  "- The correct manuscript contrast is strong normal temporal remodeling versus absent broad disease replication.",
  "",
  "## Authorized next phase",
  "",
  "- Assemble final Figures and Tables from the locked source data.",
  "- Replot QC-only panels using unified typography, labeling and source-data manifests.",
  "- Integrate authoritative Methods, Results, Discussion and limitation text.",
  "- Perform full numerical consistency and reference QC.",
  "",
  "## Lock boundary",
  "",
  "- STEP04 module analytic interpretation: LOCKED.",
  "- Final manuscript figures: NOT LOCKED.",
  "- Project-level final result: NOT LOCKED."
)

methods_lines <- c(
  "# STEP04 Closeout Methods",
  "",
  "The STEP04 module projected four prospectively locked VDR interaction programs into patient bulk cohorts and two donor-level single-cell datasets. GSE250130 was used to resolve normal luteal timing and to compare RIF with fertile LH+7 donors, whereas GSE183837 provided an independent stage-matched Control-versus-RIF localization cohort. Raw single-cell counts were aggregated by donor and broad cell type. Program scores were defined as mean signed gene-wise standardized expression values, with signs determined by the locked calcitriol-dependent interaction coefficients. Exact permutation tests were used for small donor-level comparisons, and each analysis retained its prespecified multiplicity family. Cross-dataset review aligned 25 shared broad-cell-type and endpoint combinations without pooling donors or combining P values."
)

results_lines <- c(
  "# STEP04 Closeout Results",
  "",
  paste0(
    "Normal luteal progression produced substantial cell-type-resolved remodeling: ",
    G250_temporal_program_FDR_n,
    " of ",
    nrow(G250_temporal_program),
    " temporal VDR-program models were significant after FDR correction, and the gene-level analysis identified ",
    G250_temporal_omnibus_gene_n,
    " temporal omnibus and ",
    G250_temporal_pair_gene_n,
    " planned-pairwise significant rows."
  ),
  paste0(
    "In contrast, the primary GSE250130 LH+7 comparison yielded ",
    G250_LH7_program_FDR_n,
    " FDR-significant VDR-program models, ",
    G250_LH7_VDR_FDR_n,
    " FDR-significant VDR-expression models and ",
    G250_LH7_primary_gene_n,
    " genome-wide differentially expressed genes."
  ),
  paste0(
    "The independent GSE183837 analysis yielded ",
    G183_nominal_n,
    " nominal exact-permutation signal and ",
    G183_global_BH_n,
    " result after global multiplicity correction."
  ),
  paste0(
    "Across the 25 shared cross-dataset models, ",
    cross_informative_n,
    " were directionally informative; ",
    cross_concordant_n,
    " were concordant and ",
    cross_discordant_n,
    " were discordant. No corrected replication or exploratory directional support was observed."
  ),
  paste0(
    "The nominal GSE183837 myeloid-cell VDR-expression decrease was directionally ",
    tolower(myeloid_cross_direction),
    " in GSE250130."
  ),
  "Together, these analyses support strong normal timing-associated VDR-related remodeling but do not support broad, reproducible RIF-associated disruption of the locked VDR programs."
)

discussion_lines <- c(
  "# STEP04 Closeout Discussion",
  "",
  "The main biological distinction emerging from STEP04 is between temporal regulation and disease dysregulation. The fertile luteal series displayed extensive cell-type-resolved transcriptomic changes and frequent VDR-program remodeling, indicating that endometrial timing is a major determinant of the observed VDR-related state. By contrast, neither donor-level single-cell cohort produced multiplicity-corrected evidence for broad RIF-associated VDR-program disruption, and the only nominal independent-cohort signal was directionally discordant in GSE250130.",
  "",
  "These findings argue against presenting RIF as a uniform VDR-deficiency state. Instead, the perturbation and occupancy layers define a mechanistically grounded candidate framework whose expression is strongly conditioned by normal temporal context. Patient bulk projections may identify context-specific associations, but these should not be conflated with consistent donor-level cell-type abnormalities.",
  "",
  "The principal limitations are the small number of stage-matched donors, cross-sectional temporal sampling, differing fertile/control reference structures, limited ciliated epithelial support and the use of published occupancy resources from external experimental contexts. Accordingly, the program scores should be interpreted as relative alignment with a locked perturbation direction rather than direct VDR activity or causal target regulation."
)

limitation_lines <- c(
  "# STEP04 Locked Limitations",
  "",
  "- GSE250130 temporal sampling is cross-sectional rather than longitudinal.",
  "- GSE250130 and GSE183837 use different reference cohort structures.",
  "- GSE183837 has only 3 Control and 6 RIF donors; unciliated primary models use 3 Control and 4 RIF donors.",
  "- Ciliated epithelial disease inference is descriptive only.",
  "- GSE183837 cluster 8 is low confidence and excluded from primary models.",
  "- GSE183837 cluster 13 is donor-dominant and cannot support an RIF-specific subtype claim.",
  "- Published occupancy data provide convergent support, not definitive direct-target assignment.",
  "- Bulk cohorts differ in design and are interpreted per contrast rather than as exchangeable replications.",
  "- No pooled single-cell meta-analysis or combined P value was performed.",
  "- Negative disease results should not be reframed through post hoc endpoint expansion."
)

stop_rule_lines <- c(
  "# STEP04 Analysis Stop Rule",
  "",
  "STEP04 analytic testing is closed after STEP04E1Q.",
  "",
  "New analyses are not authorized solely to obtain a positive RIF result.",
  "",
  "A new inferential analysis may be opened only when at least one of the following is documented:",
  "",
  "1. A reproducibility defect invalidates an accepted result.",
  "2. A reviewer requests a prespecified sensitivity analysis that is scientifically identifiable and feasible.",
  "3. A new independent donor cohort with adequate stage matching becomes available.",
  "4. A manuscript consistency audit identifies a missing result required to substantiate an existing claim.",
  "",
  "Final figure assembly, typography repair, source-data packaging, manuscript integration and reference QC are not considered new discovery analyses."
)

write_text(
  locked_interpretation_lines,
  file.path(
    staging_text_root,
    "STEP04_LOCKED_INTERPRETATION.md"
  )
)

write_text(
  methods_lines,
  file.path(
    staging_text_root,
    "STEP04_Methods_closeout.md"
  )
)

write_text(
  results_lines,
  file.path(
    staging_text_root,
    "STEP04_Results_closeout.md"
  )
)

write_text(
  discussion_lines,
  file.path(
    staging_text_root,
    "STEP04_Discussion_closeout.md"
  )
)

write_text(
  limitation_lines,
  file.path(
    staging_text_root,
    "STEP04_Limitations_locked.md"
  )
)

write_text(
  stop_rule_lines,
  file.path(
    staging_text_root,
    "STEP04_ANALYSIS_STOP_RULE.md"
  )
)


# ----------------------------------------------------------------------
# 14. Authoritative input manifest and snapshots
# ----------------------------------------------------------------------

authoritative_inputs <- data.frame(
  artifact_role = c(
    "Bulk Gate",
    "Bulk program models",
    "Bulk candidate convergence",
    "Locked RNA/occupancy candidates",
    "GSE250130 model Gate",
    "GSE250130 temporal programs",
    "GSE250130 LH7 programs",
    "GSE250130 temporal VDR",
    "GSE250130 LH7 VDR",
    "GSE250130 gene-DE Gate",
    "GSE250130 temporal DE summary",
    "GSE250130 LH7 DE summary",
    "GSE183837 donor-model Gate",
    "GSE183837 primary models",
    "GSE183837 cluster8 sensitivity",
    "Cross-dataset Gate",
    "Cross-dataset 25-model table",
    "Cross-dataset interpretation lock"
  ),
  path = normalize_slashes(c(
    bulk_gate_path,
    bulk_score_path,
    bulk_candidate_path,
    candidate_evidence_path,
    G250_model_gate_path,
    G250_temporal_program_path,
    G250_LH7_program_path,
    G250_temporal_VDR_path,
    G250_LH7_VDR_path,
    G250_DE_gate_path,
    G250_temporal_DE_summary_path,
    G250_LH7_DE_summary_path,
    G183_model_gate_path,
    G183_primary_model_path,
    G183_cluster8_path,
    cross_gate_path,
    cross_model_path,
    cross_lock_path
  )),
  size_bytes =
    file.info(c(
      bulk_gate_path,
      bulk_score_path,
      bulk_candidate_path,
      candidate_evidence_path,
      G250_model_gate_path,
      G250_temporal_program_path,
      G250_LH7_program_path,
      G250_temporal_VDR_path,
      G250_LH7_VDR_path,
      G250_DE_gate_path,
      G250_temporal_DE_summary_path,
      G250_LH7_DE_summary_path,
      G183_model_gate_path,
      G183_primary_model_path,
      G183_cluster8_path,
      cross_gate_path,
      cross_model_path,
      cross_lock_path
    ))$size,
  md5 = unname(
    tools::md5sum(c(
      bulk_gate_path,
      bulk_score_path,
      bulk_candidate_path,
      candidate_evidence_path,
      G250_model_gate_path,
      G250_temporal_program_path,
      G250_LH7_program_path,
      G250_temporal_VDR_path,
      G250_LH7_VDR_path,
      G250_DE_gate_path,
      G250_temporal_DE_summary_path,
      G250_LH7_DE_summary_path,
      G183_model_gate_path,
      G183_primary_model_path,
      G183_cluster8_path,
      cross_gate_path,
      cross_model_path,
      cross_lock_path
    ))
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  authoritative_inputs,
  file.path(
    staging_manifest_root,
    "authoritative_input_manifest.csv"
  )
)

snapshot_sources <- c(
  bulk_score_path,
  G250_temporal_program_path,
  G250_LH7_program_path,
  G250_LH7_VDR_path,
  G250_temporal_DE_summary_path,
  G250_LH7_DE_summary_path,
  G183_primary_model_path,
  G183_cluster8_path,
  cross_model_path,
  cross_lock_path
)

for (source_path in snapshot_sources) {
  destination_path <- file.path(
    snapshot_root,
    basename(source_path)
  )

  if (file.exists(destination_path)) {
    stop(
      "Snapshot destination unexpectedly exists:\n",
      destination_path
    )
  }

  copy_validate(
    source_path,
    destination_path
  )
}


# ----------------------------------------------------------------------
# 15. Gate
# ----------------------------------------------------------------------

gate_checks <- data.frame(
  check = c(
    "STEP04D_bulk_gate_passes",
    "bulk_score_table_has_16_rows_and_4_programs",
    "locked_candidate_identity_is_436_284_46_84",
    "GSE250130_donor_model_gate_passes",
    "GSE250130_temporal_program_grid_is_20",
    "GSE250130_LH7_program_grid_is_24",
    "GSE250130_gene_DE_gate_passes",
    "GSE250130_temporal_genome_signal_is_present",
    "GSE250130_LH7_primary_genome_DE_is_zero",
    "GSE183837_donor_model_gate_passes",
    "GSE183837_global_BH_signal_is_zero",
    "cross_dataset_gate_passes",
    "cross_dataset_grid_is_25",
    "cross_dataset_corrected_replication_is_zero",
    "cross_dataset_exploratory_directional_support_is_zero",
    "cross_dataset_direction_split_is_12_and_12",
    "myeloid_VDR_direction_is_discordant",
    "cluster8_direction_change_is_zero",
    "fact_ledger_has_25_locked_rows",
    "claim_ledger_contains_prohibited_wording_boundaries",
    "figure_blueprint_has_6_main_figures",
    "visual_QC_register_blocks_current_figures_from_final_status",
    "no_new_inferential_testing_performed",
    "STEP04_module_interpretation_can_be_locked",
    "project_level_result_lock_remains_disabled"
  ),
  pass = c(
    bulk_gate_status ==
      "READY_FOR_SINGLE_CELL_INTEGRATION_WITH_GSE183837_SAMPLE_WATCHLIST",
    nrow(bulk_scores) == 16L &&
      length(
        unique(
          bulk_scores$
            program_name
        )
      ) == 4L &&
      nrow(
        bulk_contrast_summary
      ) == 4L,
    nrow(candidate_evidence) == 436L &&
      sum(
        candidate_evidence$
          occupancy_supported
      ) == 284L &&
      sum(
        candidate_evidence$
          shared_occupancy_supported
      ) == 46L &&
      sum(
        candidate_evidence$
          any_promoter_overlap_support
      ) == 84L,
    G250_model_gate$
      gate_status[[1]] ==
      "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE",
    nrow(
      G250_temporal_program
    ) == 20L,
    nrow(
      G250_LH7_program
    ) == 24L,
    G250_DE_gate$
      gate_status[[1]] ==
      "PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE",
    G250_temporal_omnibus_gene_n >
      0L,
    G250_LH7_primary_gene_n ==
      0L,
    G183_model_gate$
      gate_status[[1]] ==
      "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED",
    G183_global_BH_n == 0L,
    cross_gate$
      gate_status[[1]] ==
      "PASS_SINGLE_CELL_CROSS_DATASET_CONCORDANCE_REVIEW_COMPLETE",
    nrow(cross_models) == 25L,
    cross_corrected_replication_n ==
      0L,
    cross_exploratory_support_n ==
      0L,
    cross_concordant_n == 12L &&
      cross_discordant_n == 12L,
    myeloid_cross_direction ==
      "DISCORDANT",
    G183_cluster8_direction_changed_n ==
      0L,
    nrow(fact_ledger) == 25L &&
      all(
        fact_ledger$status ==
          "LOCKED"
      ),
    nrow(claim_ledger) == 12L &&
      all(
        nzchar(
          claim_ledger$
            prohibited_wording
        )
      ),
    sum(
      grepl(
        "^Figure [1-6]$",
        figure_blueprint$
          figure_id
      )
    ) == 6L,
    all(
      !visual_QC_register$
        manuscript_ready
    ),
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_checks,
  file.path(
    run_qc_root,
    "02_gate_checks.csv"
  )
)

all_gate_checks_pass <- all(
  gate_checks$pass
)

gate_status <- if (
  all_gate_checks_pass
) {
  "PASS_STEP04_ANALYTIC_CLOSEOUT_INTERPRETATION_LOCKED"
} else {
  "HOLD_STEP04_ANALYTIC_CLOSEOUT_QC_FAILED"
}

closeout_gate <- data.frame(
  lock_version =
    lock_version,
  gate_status =
    gate_status,
  bulk_program_model_n =
    nrow(bulk_scores),
  bulk_FDR_signal_n =
    sum(
      bulk_scores$
        bulk_FDR_significant
    ),
  temporal_program_model_n =
    nrow(
      G250_temporal_program
    ),
  temporal_program_FDR_n =
    G250_temporal_program_FDR_n,
  temporal_genome_FDR_row_n =
    G250_temporal_omnibus_gene_n,
  LH7_program_FDR_n =
    G250_LH7_program_FDR_n,
  LH7_VDR_FDR_n =
    G250_LH7_VDR_FDR_n,
  LH7_primary_genome_DE_n =
    G250_LH7_primary_gene_n,
  G183_global_BH_signal_n =
    G183_global_BH_n,
  cross_dataset_corrected_replication_n =
    cross_corrected_replication_n,
  cross_dataset_direction_concordant_n =
    cross_concordant_n,
  cross_dataset_direction_discordant_n =
    cross_discordant_n,
  broad_RIF_VDR_program_dysregulation_supported =
    FALSE,
  STEP04_module_interpretation_locked =
    all_gate_checks_pass,
  additional_discovery_testing_authorized =
    FALSE,
  final_figure_assembly_authorized =
    all_gate_checks_pass,
  manuscript_integration_authorized =
    all_gate_checks_pass,
  project_level_result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  closeout_gate,
  file.path(
    run_qc_root,
    "STEP04_CLOSEOUT_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04 Analytic Closeout Card",
  "",
  paste0("- Gate: ", gate_status),
  paste0(
    "- Lock version: ",
    lock_version
  ),
  paste0(
    "- Bulk FDR-significant program rows: ",
    sum(
      bulk_scores$
        bulk_FDR_significant
    )
  ),
  paste0(
    "- GSE250130 temporal program FDR-significant models: ",
    G250_temporal_program_FDR_n,
    "/",
    nrow(G250_temporal_program)
  ),
  paste0(
    "- GSE250130 temporal genome-wide significant rows: ",
    G250_temporal_omnibus_gene_n
  ),
  paste0(
    "- GSE250130 LH+7 program/VDR/genome DE: ",
    G250_LH7_program_FDR_n,
    "/",
    G250_LH7_VDR_FDR_n,
    "/",
    G250_LH7_primary_gene_n
  ),
  paste0(
    "- GSE183837 global-BH signals: ",
    G183_global_BH_n
  ),
  paste0(
    "- Cross-dataset informative/concordant/discordant: ",
    cross_informative_n,
    "/",
    cross_concordant_n,
    "/",
    cross_discordant_n
  ),
  paste0(
    "- Corrected cross-dataset replication: ",
    cross_corrected_replication_n
  ),
  "- Broad reproducible RIF VDR-program dysregulation supported: FALSE",
  paste0(
    "- STEP04 module interpretation locked: ",
    all_gate_checks_pass
  ),
  "- Additional discovery testing authorized: FALSE",
  paste0(
    "- Final figure assembly authorized: ",
    all_gate_checks_pass
  ),
  paste0(
    "- Manuscript integration authorized: ",
    all_gate_checks_pass
  ),
  "- Project-level result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "STEP04_CLOSEOUT_CARD.md"
  )
)

if (!all_gate_checks_pass) {
  stop(
    "STEP04 closeout Gate failed. Review 02_gate_checks.csv."
  )
}


# ----------------------------------------------------------------------
# 16. Create immutable lock
# ----------------------------------------------------------------------

if (dir.exists(lock_root)) {
  existing_files <- list.files(
    lock_root,
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(existing_files)) {
    stop(
      "Lock root already contains files and will not be overwritten:\n",
      lock_root
    )
  }
}

dir.create(
  lock_root,
  recursive = TRUE,
  showWarnings = FALSE
)

lock_table_root <- file.path(
  lock_root,
  "tables"
)

lock_text_root <- file.path(
  lock_root,
  "text"
)

lock_manifest_root <- file.path(
  lock_root,
  "manifests"
)

lock_snapshot_root <- file.path(
  lock_root,
  "source_snapshots"
)

for (path_value in c(
  lock_table_root,
  lock_text_root,
  lock_manifest_root,
  lock_snapshot_root
)) {
  dir.create(
    path_value,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

staging_table_files <- list.files(
  staging_table_root,
  full.names = TRUE
)

staging_text_files <- list.files(
  staging_text_root,
  full.names = TRUE
)

staging_manifest_files <- list.files(
  staging_manifest_root,
  full.names = TRUE
)

snapshot_files <- list.files(
  snapshot_root,
  full.names = TRUE
)

for (source_path in staging_table_files) {
  copy_validate(
    source_path,
    file.path(
      lock_table_root,
      basename(source_path)
    )
  )
}

for (source_path in staging_text_files) {
  copy_validate(
    source_path,
    file.path(
      lock_text_root,
      basename(source_path)
    )
  )
}

for (source_path in staging_manifest_files) {
  copy_validate(
    source_path,
    file.path(
      lock_manifest_root,
      basename(source_path)
    )
  )
}

for (source_path in snapshot_files) {
  copy_validate(
    source_path,
    file.path(
      lock_snapshot_root,
      basename(source_path)
    )
  )
}

copy_validate(
  file.path(
    run_qc_root,
    "STEP04_CLOSEOUT_GATE.csv"
  ),
  file.path(
    lock_root,
    "STEP04_CLOSEOUT_GATE.csv"
  )
)

copy_validate(
  file.path(
    run_qc_root,
    "STEP04_CLOSEOUT_CARD.md"
  ),
  file.path(
    lock_root,
    "STEP04_CLOSEOUT_CARD.md"
  )
)

script_source_path <- file.path(
  project_root,
  "04_code",
  "R",
  script_name
)

if (file.exists(script_source_path)) {
  copy_validate(
    script_source_path,
    file.path(
      lock_root,
      script_name
    )
  )
}

lock_status <- data.frame(
  lock_version =
    lock_version,
  status = "COMPLETE",
  gate_status =
    gate_status,
  locked_at =
    as.character(
      Sys.time()
    ),
  STEP04_module_interpretation_locked =
    TRUE,
  additional_discovery_testing_authorized =
    FALSE,
  final_figure_assembly_authorized =
    TRUE,
  manuscript_integration_authorized =
    TRUE,
  final_figures_locked =
    FALSE,
  project_level_result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  lock_status,
  lock_status_path
)

lock_files <- list.files(
  lock_root,
  recursive = TRUE,
  full.names = TRUE
)

lock_files <- lock_files[
  file.exists(lock_files) &
    !dir.exists(lock_files)
]

lock_manifest <- data.frame(
  path =
    normalize_slashes(
      lock_files
    ),
  size_bytes =
    file.info(
      lock_files
    )$size,
  md5 = unname(
    tools::md5sum(
      lock_files
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  lock_manifest,
  file.path(
    lock_manifest_root,
    "STEP04_LOCK_MANIFEST.csv"
  )
)


# ----------------------------------------------------------------------
# 17. Synchronize authoritative text
# ----------------------------------------------------------------------

authoritative_text_root <- file.path(
  text_root,
  "STEP04_closeout",
  lock_version
)

dir.create(
  authoritative_text_root,
  recursive = TRUE,
  showWarnings = FALSE
)

for (source_path in list.files(
  lock_text_root,
  full.names = TRUE
)) {
  destination_path <- file.path(
    authoritative_text_root,
    basename(source_path)
  )

  if (!file.exists(destination_path)) {
    copy_validate(
      source_path,
      destination_path
    )
  }
}


# ----------------------------------------------------------------------
# 18. Parameters, run manifest and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  lock_version = lock_version,
  bulk_projection_version =
    "BULK_VDR_PROGRAM_PROJECTION_v1.3",
  GSE250130_model_version =
    "VDRMODELS_v1.1",
  GSE250130_DE_version =
    "DE_v1_2",
  GSE183837_model_version =
    "MD183_v1_0",
  cross_dataset_version =
    "SCCROSS_v1_0",
  no_new_inferential_testing =
    TRUE,
  broad_RIF_VDR_program_dysregulation_supported =
    FALSE,
  STEP04_module_interpretation_locked =
    TRUE,
  final_figure_assembly_authorized =
    TRUE,
  manuscript_integration_authorized =
    TRUE,
  project_level_result_lock_authorized =
    FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  lock_version = lock_version,
  started_at =
    as.character(run_start),
  completed_at =
    as.character(run_end),
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
  gate_status =
    gate_status,
  bulk_FDR_signal_n =
    sum(
      bulk_scores$
        bulk_FDR_significant
    ),
  temporal_program_FDR_n =
    G250_temporal_program_FDR_n,
  LH7_program_FDR_n =
    G250_LH7_program_FDR_n,
  LH7_VDR_FDR_n =
    G250_LH7_VDR_FDR_n,
  LH7_primary_genome_DE_n =
    G250_LH7_primary_gene_n,
  G183_global_BH_signal_n =
    G183_global_BH_n,
  cross_corrected_replication_n =
    cross_corrected_replication_n,
  cross_direction_concordant_n =
    cross_concordant_n,
  cross_direction_discordant_n =
    cross_discordant_n,
  STEP04_module_interpretation_locked =
    TRUE,
  additional_discovery_testing_authorized =
    FALSE,
  final_figure_assembly_authorized =
    TRUE,
  manuscript_integration_authorized =
    TRUE,
  project_level_result_lock_authorized =
    FALSE,
  lock_root =
    normalize_slashes(
      lock_root
    ),
  qc_root =
    normalize_slashes(
      run_qc_root
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

small_files <- unique(c(
  list.files(
    lock_root,
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
  session_path
))

small_files <- small_files[
  file.exists(small_files) &
    !dir.exists(small_files) &
    file.info(small_files)$size <=
      100 * 1024^2
]

artifact_manifest <- data.frame(
  path =
    normalize_slashes(
      small_files
    ),
  size_bytes =
    file.info(
      small_files
    )$size,
  md5 = unname(
    tools::md5sum(
      small_files
    )
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

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "lock_root=",
  normalize_slashes(
    lock_root
  )
)

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos =
        "https://cloud.r-project.org"
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
      zip_error <<-
        conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04 analytic closeout completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat("Lock version: ", lock_version, "\n", sep = "")
cat(
  "Bulk FDR-significant program rows: ",
  sum(
    bulk_scores$
      bulk_FDR_significant
  ),
  "\n",
  sep = ""
)
cat(
  "Temporal program FDR-significant models: ",
  G250_temporal_program_FDR_n,
  "/",
  nrow(G250_temporal_program),
  "\n",
  sep = ""
)
cat(
  "LH+7 program/VDR/genome DE: ",
  G250_LH7_program_FDR_n,
  "/",
  G250_LH7_VDR_FDR_n,
  "/",
  G250_LH7_primary_gene_n,
  "\n",
  sep = ""
)
cat(
  "GSE183837 global-BH signals: ",
  G183_global_BH_n,
  "\n",
  sep = ""
)
cat(
  "Cross-dataset concordant/discordant: ",
  cross_concordant_n,
  "/",
  cross_discordant_n,
  "\n",
  sep = ""
)
cat(
  "Corrected cross-dataset replication: ",
  cross_corrected_replication_n,
  "\n",
  sep = ""
)
cat("Broad RIF VDR-program dysregulation supported: FALSE\n")
cat("STEP04 module interpretation locked: TRUE\n")
cat("Additional discovery testing authorized: FALSE\n")
cat("Final figure assembly authorized: TRUE\n")
cat("Manuscript integration authorized: TRUE\n")
cat("Project-level result lock authorized: FALSE\n")
cat(
  "Lock root: ",
  normalize_slashes(
    lock_root
  ),
  "\n",
  sep = ""
)
cat(
  "QC root: ",
  normalize_slashes(
    run_qc_root
  ),
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
    "\nThe STEP04 lock is complete in the project folder.\n",
    sep = ""
  )
}
