#!/usr/bin/env Rscript
# Public analysis module for the VDR endometrial transcriptome study.
# Set VDR_PROJECT_ROOT to the project directory containing the expected public GEO inputs.
# This module preserves the executed analysis order; raw GEO data are not redistributed.


# ==============================================================================
# Analysis stage 1: 63_step04E1I_GSE183837_integration_preflight_v1
# ==============================================================================

# ======================================================================
# 63_step04E1I_GSE183837_integration_preflight_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Verify that GSE183837 is ready for independent single-cell integration.
# - Validate the locked 3-Control / 6-RIF donor design.
# - Validate all nine donor-level pre-QC RDS objects.
# - Audit package availability, existing checkpoints, resource headroom,
#   output paths and the planned integration parameters.
# - Do not modify any upstream object and do not perform integration.
#
# UPSTREAM AUTHORIZATION
# - Requires STEP04E1H:
#   PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE
#   GSE183837_independent_localization_authorized = TRUE
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
dataset_id <- "GSE183837"

script_name <- "63_step04E1I_GSE183837_integration_preflight_v1.R"
script_version <- "STEP04E1I_GSE183837_INTEGRATION_PREFLIGHT_v1.0"
planned_integration_version <- "SCI183_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1I_",
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

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  path_n <- nchar(
    normalize_slashes(path),
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
      normalize_slashes(path)
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
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
  assert_path_safe(path)
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
  }

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
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
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(pkg)
  )
}

detect_first_column <- function(
  column_names,
  candidates
) {
  lower_names <- tolower(column_names)

  for (candidate in tolower(candidates)) {
    hit <- which(lower_names == candidate)

    if (length(hit)) {
      return(column_names[hit[1]])
    }
  }

  ""
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

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^2, 2)
  }
}

get_free_disk_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_LogicalDisk",
    "-Filter \"DeviceID='C:'\").FreeSpace"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^3, 2)
  }
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
# 2. Project paths
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

design_path <- file.path(
  results_root,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0",
  "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"
)

preQC_root <- file.path(
  results_root,
  "S4D_sc_preQC",
  "SCP_v1_0",
  dataset_id
)

planned_output_root <- file.path(
  results_root,
  "S4E1_sc",
  planned_integration_version
)

planned_checkpoint_root <- file.path(
  planned_output_root,
  "checkpoints",
  dataset_id
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1I_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_preflight"
  )
)

for (path_value in c(
  execution_root,
  parameter_root,
  environment_root,
  run_qc_root
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
# 3. Validate STEP04E1H authorization
# ----------------------------------------------------------------------

DE_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

DE_qc_dirs <- DE_qc_dirs[
  grepl(
    "^STEP04E1H_[0-9]{8}_[0-9]{6}_deQC$",
    basename(DE_qc_dirs),
    perl = TRUE
  )
]

if (!length(DE_qc_dirs)) {
  stop(
    "No completed STEP04E1H gene-level DE QC directory found."
  )
}

latest_DE_qc <- DE_qc_dirs[
  which.max(file.info(DE_qc_dirs)$mtime)
]

DE_gate_path <- file.path(
  latest_DE_qc,
  "GENE_DE_GATE.csv"
)

DE_gate <- safe_read_csv(DE_gate_path)

if (
  is.null(DE_gate) ||
    nrow(DE_gate) != 1L ||
    DE_gate$gate_status[[1]] !=
      "PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE" ||
    !isTRUE(
      DE_gate$
        GSE183837_independent_localization_authorized[[1]]
    )
) {
  stop(
    "STEP04E1H Gate does not authorize GSE183837 independent localization."
  )
}


# ----------------------------------------------------------------------
# 4. Validate locked design
# ----------------------------------------------------------------------

if (!file.exists(design_path)) {
  stop(
    "Locked GSE183837 sample design is missing:\n",
    design_path
  )
}

design <- safe_read_csv(design_path)

if (
  is.null(design) ||
    !"sample_accession" %in% colnames(design) ||
    nrow(design) != 9L ||
    anyDuplicated(design$sample_accession) > 0L
) {
  stop(
    "Locked GSE183837 sample design failed identity validation."
  )
}

group_column <- detect_first_column(
  colnames(design),
  c(
    "clinical_group",
    "group",
    "condition",
    "phenotype",
    "case_control_status",
    "case_control"
  )
)

donor_column <- detect_first_column(
  colnames(design),
  c(
    "donor_id",
    "donor",
    "subject_id",
    "patient_id"
  )
)

if (!nzchar(group_column)) {
  stop(
    "Could not identify the Control/RIF group column in the locked design."
  )
}

group_value <- toupper(
  trimws(as.character(design[[group_column]]))
)

group_value[
  group_value %in% c(
    "CTRL",
    "CONTROL",
    "FERTILE",
    "NORMAL"
  )
] <- "CONTROL"

group_value[
  group_value %in% c(
    "RIF",
    "RECURRENT_IMPLANTATION_FAILURE"
  )
] <- "RIF"

group_counts <- table(
  factor(
    group_value,
    levels = c("CONTROL", "RIF")
  )
)

design_identity_pass <- identical(
  as.integer(group_counts),
  c(3L, 6L)
)

design_audit <- data.frame(
  metric = c(
    "locked_sample_n",
    "unique_sample_n",
    "group_column",
    "donor_column",
    "CONTROL_n",
    "RIF_n",
    "design_identity_pass"
  ),
  value = c(
    nrow(design),
    length(unique(design$sample_accession)),
    group_column,
    donor_column,
    as.integer(group_counts["CONTROL"]),
    as.integer(group_counts["RIF"]),
    design_identity_pass
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  design_audit,
  file.path(run_qc_root, "01_design_audit.csv")
)

if (!design_identity_pass) {
  stop(
    "Locked GSE183837 design is not exactly 3 Control and 6 RIF donors."
  )
}


# ----------------------------------------------------------------------
# 5. Validate nine donor-level pre-QC RDS objects
# ----------------------------------------------------------------------

if (!dir.exists(preQC_root)) {
  stop(
    "GSE183837 pre-QC directory is missing:\n",
    preQC_root
  )
}

preQC_rows <- list()

for (sample_index in seq_len(nrow(design))) {
  sample_accession <- as.character(
    design$sample_accession[sample_index]
  )

  object_path <- file.path(
    preQC_root,
    paste0(sample_accession, "_preQC.rds")
  )

  if (!file.exists(object_path)) {
    preQC_rows[[sample_index]] <- data.frame(
      sample_accession = sample_accession,
      object_exists = FALSE,
      object_size_MB = NA_real_,
      feature_n = NA_integer_,
      retained_cell_n = NA_integer_,
      feature_table_row_n = NA_integer_,
      cell_metadata_row_n = NA_integer_,
      required_components_pass = FALSE,
      matrix_dimnames_pass = FALSE,
      feature_alignment_pass = FALSE,
      cell_alignment_pass = FALSE,
      sample_identity_pass = FALSE,
      md5 = NA_character_,
      object_path = normalize_slashes(object_path),
      stringsAsFactors = FALSE
    )

    next
  }

  object_value <- tryCatch(
    readRDS(object_path),
    error = function(e) NULL
  )

  required_components <- c(
    "counts",
    "feature_table",
    "cell_metadata"
  )

  components_pass <- (
    !is.null(object_value) &&
      all(
        required_components %in%
          names(object_value)
      )
  )

  if (!components_pass) {
    preQC_rows[[sample_index]] <- data.frame(
      sample_accession = sample_accession,
      object_exists = TRUE,
      object_size_MB = round(
        file.info(object_path)$size / 1024^2,
        2
      ),
      feature_n = NA_integer_,
      retained_cell_n = NA_integer_,
      feature_table_row_n = NA_integer_,
      cell_metadata_row_n = NA_integer_,
      required_components_pass = FALSE,
      matrix_dimnames_pass = FALSE,
      feature_alignment_pass = FALSE,
      cell_alignment_pass = FALSE,
      sample_identity_pass = FALSE,
      md5 = unname(tools::md5sum(object_path)),
      object_path = normalize_slashes(object_path),
      stringsAsFactors = FALSE
    )

    next
  }

  counts <- object_value$counts
  feature_table <- object_value$feature_table
  cell_metadata <- object_value$cell_metadata

  matrix_like <- (
    inherits(counts, "Matrix") ||
      is.matrix(counts)
  )

  matrix_dimnames_pass <- (
    matrix_like &&
      !is.null(rownames(counts)) &&
      !is.null(colnames(counts)) &&
      !anyDuplicated(colnames(counts))
  )

  feature_alignment_pass <- (
    is.data.frame(feature_table) &&
      nrow(feature_table) == nrow(counts)
  )

  cell_alignment_pass <- (
    is.data.frame(cell_metadata) &&
      nrow(cell_metadata) == ncol(counts)
  )

  metadata_sample_pass <- TRUE

  if (
    "sample_accession" %in%
      colnames(cell_metadata)
  ) {
    metadata_sample_pass <- (
      length(unique(
        as.character(
          cell_metadata$sample_accession
        )
      )) == 1L &&
        unique(
          as.character(
            cell_metadata$sample_accession
          )
        ) == sample_accession
    )
  }

  preQC_rows[[sample_index]] <- data.frame(
    sample_accession = sample_accession,
    object_exists = TRUE,
    object_size_MB = round(
      file.info(object_path)$size / 1024^2,
      2
    ),
    feature_n = nrow(counts),
    retained_cell_n = ncol(counts),
    feature_table_row_n = nrow(feature_table),
    cell_metadata_row_n = nrow(cell_metadata),
    required_components_pass = components_pass,
    matrix_dimnames_pass = matrix_dimnames_pass,
    feature_alignment_pass = feature_alignment_pass,
    cell_alignment_pass = cell_alignment_pass,
    sample_identity_pass = metadata_sample_pass,
    md5 = unname(tools::md5sum(object_path)),
    object_path = normalize_slashes(object_path),
    stringsAsFactors = FALSE
  )

  rm(
    object_value,
    counts,
    feature_table,
    cell_metadata
  )
  gc()
}

preQC_audit <- do.call(
  rbind,
  preQC_rows
)

write_csv_utf8(
  preQC_audit,
  file.path(run_qc_root, "02_preQC_object_audit.csv")
)

preQC_total_cells <- sum(
  preQC_audit$retained_cell_n,
  na.rm = TRUE
)

preQC_all_pass <- (
  nrow(preQC_audit) == 9L &&
    all(preQC_audit$object_exists) &&
    all(preQC_audit$required_components_pass) &&
    all(preQC_audit$matrix_dimnames_pass) &&
    all(preQC_audit$feature_alignment_pass) &&
    all(preQC_audit$cell_alignment_pass) &&
    all(preQC_audit$sample_identity_pass) &&
    preQC_total_cells == 50997L
)


# ----------------------------------------------------------------------
# 6. Package, checkpoint and resource audit
# ----------------------------------------------------------------------

package_names <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "harmony",
  "data.table",
  "jsonlite",
  "future",
  "scDblFinder",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "zip"
)

package_audit <- data.frame(
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
  required_for_integration = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  package_audit,
  file.path(run_qc_root, "03_package_audit.csv")
)

required_packages_pass <- all(
  package_audit$installed[
    package_audit$required_for_integration
  ]
)

existing_GSE183837_files <- list.files(
  file.path(
    results_root,
    "S4E1_sc"
  ),
  pattern = "GSE183837|183837",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

checkpoint_inventory <- data.frame(
  path = normalize_slashes(
    existing_GSE183837_files
  ),
  size_MB = if (
    length(existing_GSE183837_files)
  ) {
    round(
      file.info(existing_GSE183837_files)$size /
        1024^2,
      2
    )
  } else {
    numeric()
  },
  modified_time = if (
    length(existing_GSE183837_files)
  ) {
    as.character(
      file.info(existing_GSE183837_files)$mtime
    )
  } else {
    character()
  },
  stringsAsFactors = FALSE
)

write_csv_utf8(
  checkpoint_inventory,
  file.path(run_qc_root, "04_existing_output_inventory.csv")
)

planned_paths <- data.frame(
  artifact = c(
    "planned_output_root",
    "planned_checkpoint_root",
    "planned_integrated_checkpoint",
    "planned_checkpoint_status",
    "planned_UMAP_metadata",
    "planned_doublet_QC"
  ),
  path = normalize_slashes(c(
    planned_output_root,
    planned_checkpoint_root,
    file.path(
      planned_checkpoint_root,
      "integrated.rds"
    ),
    file.path(
      planned_checkpoint_root,
      "status.csv"
    ),
    file.path(
      planned_output_root,
      "tables",
      "UMAP_metadata.csv"
    ),
    file.path(
      planned_output_root,
      "tables",
      "doublet_QC.csv"
    )
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
  file.path(run_qc_root, "05_path_audit.csv")
)

resource_audit <- data.frame(
  metric = c(
    "free_physical_memory_GB",
    "free_C_drive_GB",
    "preQC_object_total_size_GB",
    "preQC_total_cells",
    "preQC_sample_n"
  ),
  value = c(
    get_free_memory_gb(),
    get_free_disk_gb(),
    round(
      sum(preQC_audit$object_size_MB) /
        1024,
      3
    ),
    preQC_total_cells,
    nrow(preQC_audit)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resource_audit,
  file.path(run_qc_root, "06_resource_audit.csv")
)

integration_plan <- data.frame(
  parameter = c(
    "dataset_processing",
    "biological_replicate",
    "doublet_method",
    "doublet_scope",
    "normalization",
    "variable_features",
    "PCA_components",
    "Harmony_batch_variable",
    "Harmony_dimensions",
    "UMAP_dimensions",
    "cluster_resolutions",
    "primary_resolution",
    "counts_layer_policy",
    "data_layer_policy",
    "scale_data_policy",
    "annotation_policy",
    "disease_testing_policy"
  ),
  value = c(
    "GSE183837 processed independently from GSE250130",
    "donor/sample",
    "scDblFinder",
    "run independently within each donor",
    "LogNormalize scale.factor=10000",
    "3000",
    "50",
    "sample_accession",
    "1:30",
    "1:30",
    "0.2;0.4;0.6",
    "0.4",
    "retain 9 split raw-count layers",
    "join normalized data into one layer after integration",
    "remove scale.data before final checkpoint",
    "cluster markers plus canonical panels and GSE250130 broad reference",
    "no cell-level test; donor-level pseudobulk only after annotation QC"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  integration_plan,
  file.path(run_qc_root, "07_integration_plan.csv")
)


# ----------------------------------------------------------------------
# 7. Preflight Gate
# ----------------------------------------------------------------------

memory_value <- suppressWarnings(
  as.numeric(
    resource_audit$value[
      resource_audit$metric ==
        "free_physical_memory_GB"
    ]
  )
)

disk_value <- suppressWarnings(
  as.numeric(
    resource_audit$value[
      resource_audit$metric ==
        "free_C_drive_GB"
    ]
  )
)

memory_pass <- (
  is.na(memory_value) ||
    memory_value >= 12
)

disk_pass <- (
  is.na(disk_value) ||
    disk_value >= 15
)

gate_checks <- data.frame(
  check = c(
    "STEP04E1H_authorizes_GSE183837",
    "locked_design_has_9_donors",
    "locked_design_is_3_Control_6_RIF",
    "nine_preQC_objects_present",
    "all_preQC_object_structures_valid",
    "preQC_total_cells_is_50997",
    "required_integration_packages_installed",
    "planned_paths_safe",
    "free_memory_at_least_12GB_or_unknown",
    "free_disk_at_least_15GB_or_unknown"
  ),
  pass = c(
    DE_gate$
      GSE183837_independent_localization_authorized[[1]],
    nrow(design) == 9L,
    design_identity_pass,
    nrow(preQC_audit) == 9L &&
      all(preQC_audit$object_exists),
    all(
      preQC_audit$required_components_pass &
        preQC_audit$matrix_dimnames_pass &
        preQC_audit$feature_alignment_pass &
        preQC_audit$cell_alignment_pass &
        preQC_audit$sample_identity_pass
    ),
    preQC_total_cells == 50997L,
    required_packages_pass,
    all(planned_paths$path_safe),
    memory_pass,
    disk_pass
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  gate_checks,
  file.path(run_qc_root, "08_gate_checks.csv")
)

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
  "PASS_GSE183837_INTEGRATION_PREFLIGHT"
} else {
  "HOLD_GSE183837_INTEGRATION_PREFLIGHT_FAILED"
}

preflight_gate <- data.frame(
  dataset_id = dataset_id,
  script_version = script_version,
  planned_integration_version =
    planned_integration_version,
  gate_status = gate_status,
  locked_sample_n = nrow(design),
  Control_n = as.integer(
    group_counts["CONTROL"]
  ),
  RIF_n = as.integer(
    group_counts["RIF"]
  ),
  preQC_cell_n = preQC_total_cells,
  preQC_object_n = nrow(preQC_audit),
  required_packages_pass =
    required_packages_pass,
  integration_authorized =
    all_gate_checks_pass,
  automatic_annotation_authorized =
    FALSE,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "GSE183837_PREFLIGHT_GATE.csv"
)

write_csv_utf8(
  preflight_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1I GSE183837 Integration Preflight Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0(
    "- Planned integration version: ",
    planned_integration_version
  ),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Locked donors: ",
    nrow(design),
    " (Control ",
    as.integer(group_counts["CONTROL"]),
    "; RIF ",
    as.integer(group_counts["RIF"]),
    ")"
  ),
  paste0(
    "- Donor-level pre-QC objects: ",
    nrow(preQC_audit)
  ),
  paste0(
    "- Pre-QC retained cells: ",
    format(
      preQC_total_cells,
      big.mark = ","
    )
  ),
  paste0(
    "- Required packages complete: ",
    required_packages_pass
  ),
  paste0(
    "- Integration authorized: ",
    all_gate_checks_pass
  ),
  "- Automatic annotation authorized: FALSE",
  "- Disease pseudobulk testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "GSE183837_PREFLIGHT_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 integration preflight interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Authorized next operation",
  "",
  "- Run sample-wise scDblFinder on the nine pre-QC donor objects.",
  "- Integrate GSE183837 independently using sample accession as the Harmony batch variable.",
  "- Save staged checkpoints before marker extraction or biological annotation.",
  "",
  "## Not authorized yet",
  "",
  "- No automatic broad-cell annotation is locked at this stage.",
  "- No cell-level Control-versus-RIF test is permitted.",
  "- No donor-level pseudobulk disease test is permitted until broad-cell annotation QC passes.",
  "- No GSE250130 and GSE183837 cell pooling is permitted."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "GSE183837_PREFLIGHT_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 8. Parameters, session, manifest and ZIP
# ----------------------------------------------------------------------

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop(
    "jsonlite is required to write the parameter manifest."
  )
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  dataset_id = dataset_id,
  upstream_gate =
    DE_gate$gate_status[[1]],
  planned_integration_version =
    planned_integration_version,
  expected_sample_n = 9L,
  expected_Control_n = 3L,
  expected_RIF_n = 6L,
  expected_preQC_cell_n = 50997L,
  minimum_free_memory_GB = 12,
  minimum_free_disk_GB = 15,
  integration_plan = as.list(
    setNames(
      integration_plan$value,
      integration_plan$parameter
    )
  ),
  result_lock_authorized = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
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
  locked_sample_n = nrow(design),
  preQC_object_n = nrow(preQC_audit),
  preQC_cell_n = preQC_total_cells,
  required_packages_pass =
    required_packages_pass,
  integration_authorized =
    all_gate_checks_pass,
  planned_output_root =
    normalize_slashes(planned_output_root),
  qc_root =
    normalize_slashes(run_qc_root),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)

artifact_files <- unique(c(
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

log_event(
  "Completed",
  run_id,
  "gate_status=",
  gate_status,
  "preQC_objects=",
  nrow(preQC_audit),
  "preQC_cells=",
  preQC_total_cells
)

small_files <- unique(c(
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
    !dir.exists(small_files)
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
cat("STEP04E1I GSE183837 integration preflight completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Locked donors: ",
  nrow(design),
  " (Control ",
  as.integer(group_counts["CONTROL"]),
  "; RIF ",
  as.integer(group_counts["RIF"]),
  ")\n",
  sep = ""
)
cat(
  "Pre-QC objects: ",
  nrow(preQC_audit),
  "\n",
  sep = ""
)
cat(
  "Pre-QC retained cells: ",
  format(preQC_total_cells, big.mark = ","),
  "\n",
  sep = ""
)
cat(
  "Integration authorized: ",
  all_gate_checks_pass,
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
    "\nThe preflight report is complete in the QC folder.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 2: 64_step04E1J_GSE183837_doublet_Harmony_integration_v1_4_donor_gate_fixed
# ==============================================================================

# ======================================================================
# 64_step04E1J_GSE183837_doublet_Harmony_integration_v1_4_donor_gate_fixed.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Run scDblFinder independently for each of the nine GSE183837 donors.
# - Save resumable donor-level singlet Seurat checkpoints.
# - Merge the nine singlet objects while retaining split raw-count layers.
# - Perform LogNormalize, variable-feature selection, PCA and Harmony.
# - Generate Harmony UMAP and clusters at resolutions 0.2, 0.4 and 0.6.
# - Use resolution 0.4 as the primary cluster solution.
# - Join normalized data layers only; preserve nine raw-count layers.
# - Save a staged merged-count checkpoint and a final integrated checkpoint.
# - Produce execution-level technical QC and authorize a separate detailed
#   integration-QC step. Marker analysis is not authorized by this script.
#
# UPSTREAM GATE
# - PASS_GSE183837_INTEGRATION_PREFLIGHT
#
# SCIENTIFIC GUARDRAILS
# - GSE183837 is processed independently from GSE250130.
# - Doublet detection is performed separately within each donor/capture.
# - No cell-level Control-versus-RIF hypothesis test is performed.
# - No biological annotation is locked.
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

set.seed(20260806)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE183837"

script_name <-
  "64_step04E1J_GSE183837_doublet_Harmony_integration_v1_4_donor_gate_fixed.R"

script_version <-
  "STEP04E1J_GSE183837_DOUBLET_HARMONY_INTEGRATION_v1.4_DONOR_GATE_FIXED"

integration_version <- "SCI183_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1J_",
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

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  path_n <- nchar(
    normalize_slashes(path),
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
      normalize_slashes(path)
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
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
  assert_path_safe(path)

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
  }

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

safe_package_version <- function(package_name) {
  if (
    !requireNamespace(
      package_name,
      quietly = TRUE
    )
  ) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(package_name)
  )
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

  toupper(trimws(as.character(x))) %in%
    c("TRUE", "T", "1", "YES", "Y")
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

detect_first_column <- function(
  column_names,
  candidates
) {
  lower_names <- tolower(column_names)

  for (candidate in tolower(candidates)) {
    hit <- which(lower_names == candidate)

    if (length(hit)) {
      return(column_names[hit[1]])
    }
  }

  ""
}

safe_factor_group <- function(x) {
  values <- toupper(
    trimws(as.character(x))
  )

  values[
    values %in% c(
      "CTRL",
      "CONTROL",
      "FERTILE",
      "NORMAL"
    )
  ] <- "CONTROL"

  values[
    values %in% c(
      "RIF",
      "RECURRENT_IMPLANTATION_FAILURE"
    )
  ] <- "RIF"

  factor(
    values,
    levels = c("CONTROL", "RIF")
  )
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^2, 2)
  }
}

get_free_disk_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_LogicalDisk",
    "-Filter \"DeviceID='C:'\").FreeSpace"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^3, 2)
  }
}

prepare_unique_UMAP_metadata <- function(
  metadata,
  UMAP_embedding,
  audit_path
) {
  if (
    is.null(rownames(metadata)) ||
      is.null(rownames(UMAP_embedding))
  ) {
    stop("Metadata or UMAP embedding lacks row names.")
  }

  if (
    !all(
      rownames(UMAP_embedding) %in%
        rownames(metadata)
    )
  ) {
    stop("UMAP cells are not fully represented in Seurat metadata.")
  }

  metadata_aligned <- metadata[
    rownames(UMAP_embedding),
    ,
    drop = FALSE
  ]

  original_column_names <- colnames(
    metadata_aligned
  )

  duplicate_name_table <- data.frame(
    column_name = unique(
      original_column_names[
        duplicated(original_column_names) |
          duplicated(
            original_column_names,
            fromLast = TRUE
          )
      ]
    ),
    occurrence_n = vapply(
      unique(
        original_column_names[
          duplicated(original_column_names) |
            duplicated(
              original_column_names,
              fromLast = TRUE
            )
        ]
      ),
      function(column_name) {
        sum(
          original_column_names ==
            column_name
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  reserved_names <- c(
    "cell_barcode",
    "UMAPH_1",
    "UMAPH_2"
  )

  reserved_collision <- data.frame(
    reserved_name = reserved_names,
    present_in_original_metadata =
      reserved_names %in%
        original_column_names,
    occurrence_n = vapply(
      reserved_names,
      function(column_name) {
        sum(
          original_column_names ==
            column_name
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  # Keep the first occurrence of each original metadata column.
  metadata_aligned <- metadata_aligned[
    ,
    !duplicated(
      colnames(metadata_aligned)
    ),
    drop = FALSE
  ]

  # New export columns take precedence over any existing names.
  metadata_aligned <- metadata_aligned[
    ,
    !colnames(metadata_aligned) %in%
      reserved_names,
    drop = FALSE
  ]

  output <- data.frame(
    cell_barcode =
      rownames(UMAP_embedding),
    UMAPH_1 =
      UMAP_embedding[, 1],
    UMAPH_2 =
      UMAP_embedding[, 2],
    metadata_aligned,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(output) <- make.unique(
    colnames(output),
    sep = "__duplicate_"
  )

  if (anyDuplicated(colnames(output)) > 0L) {
    stop(
      "QC UMAP metadata still contains duplicate column names after repair."
    )
  }

  required_columns <- c(
    "cell_barcode",
    "UMAPH_1",
    "UMAPH_2",
    "sample_accession",
    "donor_id",
    "analysis_group",
    "seurat_clusters"
  )

  missing_required <- setdiff(
    required_columns,
    colnames(output)
  )

  if (length(missing_required)) {
    stop(
      "QC UMAP metadata is missing required columns: ",
      paste(
        missing_required,
        collapse = ", "
      )
    )
  }

  audit_summary <- data.frame(
    metric = c(
      "original_metadata_column_n",
      "original_duplicate_name_n",
      "reserved_collision_n",
      "export_column_n",
      "export_duplicate_name_n",
      "required_columns_complete"
    ),
    value = c(
      length(original_column_names),
      nrow(duplicate_name_table),
      sum(
        reserved_collision$
          present_in_original_metadata
      ),
      ncol(output),
      anyDuplicated(colnames(output)),
      !length(missing_required)
    ),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    audit_summary,
    audit_path
  )

  write_csv_utf8(
    duplicate_name_table,
    file.path(
      dirname(audit_path),
      "00_metadata_duplicate_names.csv"
    )
  )

  write_csv_utf8(
    reserved_collision,
    file.path(
      dirname(audit_path),
      "00_metadata_reserved_collisions.csv"
    )
  )

  output
}


join_normalized_data_layers_only <- function(
  object,
  assay_name = "RNA"
) {
  layers_before <- SeuratObject::Layers(
    object,
    assay = assay_name,
    search = NA
  )

  count_layers_before <- grep(
    "^counts",
    layers_before,
    value = TRUE
  )

  data_layers_before <- grep(
    "^data",
    layers_before,
    value = TRUE
  )

  if (!length(data_layers_before)) {
    stop(
      "No normalized data layers were found before JoinLayers."
    )
  }

  if (length(data_layers_before) > 1L) {
    # JoinLayers expects one layer search category for one new layer name.
    # layers = "data" finds and stitches all data.* layers.
    object <- SeuratObject::JoinLayers(
      object = object,
      assay = assay_name,
      layers = "data",
      new = "data"
    )
  }

  layers_after <- SeuratObject::Layers(
    object,
    assay = assay_name,
    search = NA
  )

  count_layers_after <- grep(
    "^counts",
    layers_after,
    value = TRUE
  )

  data_layers_after <- grep(
    "^data",
    layers_after,
    value = TRUE
  )

  if (
    !identical(
      count_layers_before,
      count_layers_after
    )
  ) {
    stop(
      "Raw count-layer identities changed during normalized-data joining."
    )
  }

  if (
    length(data_layers_after) != 1L ||
      data_layers_after[[1]] != "data"
  ) {
    stop(
      "Normalized data layers did not resolve to exactly one 'data' layer."
    )
  }

  attr(
    object,
    "STEP04E1J_join_audit"
  ) <- list(
    count_layers_before = count_layers_before,
    count_layers_after = count_layers_after,
    data_layers_before = data_layers_before,
    data_layers_after = data_layers_after
  )

  object
}


run_harmony_Seurat_compatible <- function(
  object,
  group_variable,
  dimensions = 1:30,
  reduction_name = "pca",
  output_reduction = "harmony",
  verbose = FALSE
) {
  if (
    !requireNamespace(
      "harmony",
      quietly = TRUE
    )
  ) {
    stop("The harmony package is not installed.")
  }

  harmony_method <- tryCatch(
    utils::getS3method(
      "RunHarmony",
      "Seurat",
      envir = asNamespace("harmony")
    ),
    error = function(e) NULL
  )

  method_formals <- if (
    is.null(harmony_method)
  ) {
    character()
  } else {
    names(formals(harmony_method))
  }

  call_arguments <- list(
    object = object,
    group.by.vars = group_variable,
    reduction = reduction_name,
    dims.use = dimensions,
    reduction.save = output_reduction,
    plot_convergence = FALSE,
    verbose = verbose
  )

  # Do not pass assay.use or project.dim. Some current Harmony builds
  # reject assay.use through check_legacy_args(), even though older
  # documentation and releases accepted it.
  output <- do.call(
    harmony::RunHarmony,
    call_arguments
  )

  attr(
    output,
    "STEP04E1J_Harmony_method_formals"
  ) <- method_formals

  output
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
# 2. Package preflight
# ----------------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "harmony",
  "data.table",
  "ggplot2",
  "jsonlite",
  "scDblFinder",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "BiocParallel"
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
    paste(missing_packages, collapse = ", ")
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

design_path <- file.path(
  results_root,
  "STEP04C_locked_sample_designs",
  "COMBINED_SAMPLE_DESIGN_v1.0",
  "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"
)

preQC_root <- file.path(
  results_root,
  "S4D_sc_preQC",
  "SCP_v1_0",
  dataset_id
)

output_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

sample_checkpoint_root <- file.path(
  output_root,
  "sample_ckpt"
)

checkpoint_root <- file.path(
  output_root,
  "checkpoints",
  dataset_id
)

table_root <- file.path(
  output_root,
  "tables"
)

manifest_root <- file.path(
  output_root,
  "manifests"
)

merged_counts_checkpoint <- file.path(
  checkpoint_root,
  "merged_singlets.rds"
)

prejoin_checkpoint_path <- file.path(
  checkpoint_root,
  "harmony_prejoin.rds"
)

final_checkpoint_path <- file.path(
  checkpoint_root,
  "integrated.rds"
)

final_status_path <- file.path(
  checkpoint_root,
  "status.csv"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1J_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_integration"
  )
)

figure_root <- file.path(
  run_qc_root,
  "figures"
)

for (path_value in c(
  execution_root,
  parameter_root,
  environment_root,
  output_root,
  sample_checkpoint_root,
  checkpoint_root,
  table_root,
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
  artifact = c(
    "output_root",
    "sample_checkpoint_root",
    "checkpoint_root",
    "merged_counts_checkpoint",
    "prejoin_checkpoint",
    "final_checkpoint",
    "final_status",
    "table_root",
    "manifest_root",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    output_root,
    sample_checkpoint_root,
    checkpoint_root,
    merged_counts_checkpoint,
    prejoin_checkpoint_path,
    final_checkpoint_path,
    final_status_path,
    table_root,
    manifest_root,
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

if (!all(planned_paths$path_safe)) {
  stop(
    "One or more planned paths exceed the safety limit:\n",
    paste(
      planned_paths$path[
        !planned_paths$path_safe
      ],
      collapse = "\n"
    )
  )
}

write_csv_utf8(
  planned_paths,
  file.path(run_qc_root, "00_path_audit.csv")
)


# ----------------------------------------------------------------------
# 4. Validate upstream preflight Gate and locked design
# ----------------------------------------------------------------------

preflight_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

preflight_qc_dirs <- preflight_qc_dirs[
  grepl(
    "^STEP04E1I_[0-9]{8}_[0-9]{6}_preflight$",
    basename(preflight_qc_dirs),
    perl = TRUE
  )
]

if (!length(preflight_qc_dirs)) {
  stop(
    "No completed STEP04E1I preflight-QC directory found."
  )
}

latest_preflight_qc <- preflight_qc_dirs[
  which.max(
    file.info(preflight_qc_dirs)$mtime
  )
]

preflight_gate_path <- file.path(
  latest_preflight_qc,
  "GSE183837_PREFLIGHT_GATE.csv"
)

preflight_gate <- safe_read_csv(
  preflight_gate_path
)

if (
  is.null(preflight_gate) ||
    nrow(preflight_gate) != 1L ||
    preflight_gate$gate_status[[1]] !=
      "PASS_GSE183837_INTEGRATION_PREFLIGHT" ||
    !isTRUE(
      preflight_gate$integration_authorized[[1]]
    )
) {
  stop(
    "STEP04E1I Gate does not authorize GSE183837 integration."
  )
}

if (!file.exists(design_path)) {
  stop(
    "Locked sample design is missing:\n",
    design_path
  )
}

design <- safe_read_csv(design_path)

if (
  is.null(design) ||
    nrow(design) != 9L ||
    !"sample_accession" %in% colnames(design) ||
    anyDuplicated(design$sample_accession) > 0L
) {
  stop(
    "Locked GSE183837 design failed identity validation."
  )
}

group_column <- detect_first_column(
  colnames(design),
  c(
    "group",
    "clinical_group",
    "condition",
    "phenotype",
    "case_control_status",
    "case_control"
  )
)

donor_column <- detect_first_column(
  colnames(design),
  c(
    "donor_id",
    "donor",
    "subject_id",
    "patient_id"
  )
)

if (
  !nzchar(group_column) ||
    !nzchar(donor_column)
) {
  stop(
    "Could not identify group and donor columns in the locked design."
  )
}

design$analysis_group <- safe_factor_group(
  design[[group_column]]
)

design$analysis_donor_id <- as.character(
  design[[donor_column]]
)

if (
  anyNA(design$analysis_group) ||
    anyNA(design$analysis_donor_id) ||
    anyDuplicated(design$analysis_donor_id) > 0L ||
    !identical(
      as.integer(table(design$analysis_group)),
      c(3L, 6L)
    )
) {
  stop(
    "Locked GSE183837 donor/group design is not 3 Control and 6 RIF."
  )
}

free_memory_before_GB <- get_free_memory_gb()
free_disk_before_GB <- get_free_disk_gb()

if (
  is.finite(free_memory_before_GB) &&
    free_memory_before_GB < 12
) {
  stop(
    "Free physical memory is below 12 GB: ",
    free_memory_before_GB
  )
}

if (
  is.finite(free_disk_before_GB) &&
    free_disk_before_GB < 15
) {
  stop(
    "Free C-drive space is below 15 GB: ",
    free_disk_before_GB
  )
}

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  integration_version = integration_version,
  dataset_id = dataset_id,
  upstream_gate =
    preflight_gate$gate_status[[1]],
  expected_sample_n = 9L,
  expected_Control_n = 3L,
  expected_RIF_n = 6L,
  expected_preQC_cell_n = 50997L,
  doublet_method = "scDblFinder",
  doublet_scope =
    "INDEPENDENT_WITHIN_EACH_DONOR_CAPTURE",
  scDblFinder_clusters = TRUE,
  scDblFinder_parallel =
    "BiocParallel_SerialParam",
  normalization =
    "LogNormalize_scale_factor_10000",
  variable_features = 3000L,
  PCA_components = 50L,
  Harmony_batch_variable =
    "sample_accession",
  Harmony_dimensions = "1:30",
  Harmony_API_policy =
    "MINIMAL_CURRENT_API_NO_ASSAY_USE_NO_PROJECT_DIM",
  JoinLayers_policy =
    "LAYERS_DATA_SEARCH_CATEGORY_TO_SINGLE_DATA_LAYER_COUNTS_UNCHANGED",
  QC_metadata_policy =
    "DROP_DUPLICATE_OR_RESERVED_METADATA_NAMES_BEFORE_UMAP_EXPORT",
  donor_gate_policy =
    "COUNT_UNIQUE_DONOR_GROUP_PAIRS_NOT_CELLS",
  supersedes_failed_executions = c(
    "V1_0_ASSAY_USE_REJECTED_BY_HARMONY_CHECK_LEGACY_ARGS",
    "V1_1_EXACT_DATA_LAYER_VECTOR_WITH_SINGLE_NEW_NAME",
    "V1_2_GGPLOT_DUPLICATE_QC_METADATA_COLUMNS",
    "V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS"
  ),
  UMAP_dimensions = "1:30",
  cluster_resolutions = c(0.2, 0.4, 0.6),
  primary_resolution = 0.4,
  counts_layer_policy =
    "RETAIN_NINE_SPLIT_RAW_COUNT_LAYERS",
  data_layer_policy =
    "JOIN_NORMALIZED_DATA_LAYERS_ONLY",
  scale_data_policy =
    "REMOVE_BEFORE_FINAL_CHECKPOINT_WHEN_SUPPORTED",
  marker_analysis_authorized = FALSE,
  disease_testing_authorized = FALSE,
  result_lock_authorized = FALSE,
  random_seed = 20260806L
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

harmony_method <- tryCatch(
  utils::getS3method(
    "RunHarmony",
    "Seurat",
    envir = asNamespace("harmony")
  ),
  error = function(e) NULL
)

harmony_API_audit <- data.frame(
  metric = c(
    "harmony_version",
    "Seurat_method_found",
    "Seurat_method_formals",
    "assay_use_passed_by_script",
    "project_dim_passed_by_script",
    "reduction",
    "dims_use",
    "reduction_save"
  ),
  value = c(
    safe_package_version("harmony"),
    !is.null(harmony_method),
    if (
      is.null(harmony_method)
    ) {
      ""
    } else {
      paste(
        names(formals(harmony_method)),
        collapse = ";"
      )
    },
    FALSE,
    FALSE,
    "pca",
    "1:30",
    "harmony"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  harmony_API_audit,
  file.path(
    run_qc_root,
    "00_Harmony_API_audit.csv"
  )
)

log_event(
  "Starting",
  run_id,
  "free_memory_GB=",
  free_memory_before_GB,
  "free_disk_GB=",
  free_disk_before_GB
)


# ----------------------------------------------------------------------
# 5. Validate input objects and run sample-wise scDblFinder
# ----------------------------------------------------------------------

reference_features <- NULL
sample_status_rows <- list()
doublet_metadata_paths <- character()

for (sample_index in seq_len(nrow(design))) {
  sample_accession <- as.character(
    design$sample_accession[sample_index]
  )

  donor_id <- as.character(
    design$analysis_donor_id[sample_index]
  )

  analysis_group <- as.character(
    design$analysis_group[sample_index]
  )

  sample_dir <- file.path(
    sample_checkpoint_root,
    sample_accession
  )

  dir.create(
    sample_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  sample_checkpoint_path <- file.path(
    sample_dir,
    "singlet.rds"
  )

  sample_status_path <- file.path(
    sample_dir,
    "status.csv"
  )

  sample_doublet_path <- file.path(
    sample_dir,
    "doublet_metadata.csv"
  )

  existing_status <- safe_read_csv(
    sample_status_path
  )

  if (
    file.exists(sample_checkpoint_path) &&
      !is.null(existing_status) &&
      nrow(existing_status) == 1L &&
      existing_status$status[[1]] == "COMPLETE"
  ) {
    checkpoint_object <- tryCatch(
      readRDS(sample_checkpoint_path),
      error = function(e) NULL
    )

    checkpoint_valid <- (
      inherits(checkpoint_object, "Seurat") &&
        ncol(checkpoint_object) ==
          existing_status$singlet_cell_n[[1]] &&
        length(
          unique(
            as.character(
              checkpoint_object$sample_accession
            )
          )
        ) == 1L &&
        unique(
          as.character(
            checkpoint_object$sample_accession
          )
        ) == sample_accession
    )

    if (checkpoint_valid) {
      sample_status_rows[[sample_index]] <-
        existing_status

      doublet_metadata_paths <- c(
        doublet_metadata_paths,
        sample_doublet_path
      )

      if (is.null(reference_features)) {
        reference_features <- rownames(
          checkpoint_object
        )
      } else if (
        !identical(
          reference_features,
          rownames(checkpoint_object)
        )
      ) {
        stop(
          "Feature order differs in existing sample checkpoint: ",
          sample_accession
        )
      }

      log_event(
        "SKIP completed donor",
        sample_index,
        "/",
        nrow(design),
        sample_accession,
        "singlets=",
        ncol(checkpoint_object)
      )

      rm(checkpoint_object)
      gc()

      next
    }
  }

  sample_start <- Sys.time()

  input_path <- file.path(
    preQC_root,
    paste0(
      sample_accession,
      "_preQC.rds"
    )
  )

  if (!file.exists(input_path)) {
    stop(
      "Pre-QC object is missing:\n",
      input_path
    )
  }

  input_object <- readRDS(input_path)

  required_components <- c(
    "counts",
    "feature_table",
    "cell_metadata"
  )

  if (
    !all(
      required_components %in%
        names(input_object)
    )
  ) {
    stop(
      "Pre-QC object lacks required components: ",
      sample_accession
    )
  }

  counts <- input_object$counts
  cell_metadata <- input_object$cell_metadata

  if (
    (!inherits(counts, "Matrix") &&
      !is.matrix(counts)) ||
      nrow(counts) != 27984L ||
      ncol(counts) < 1000L ||
      is.null(rownames(counts)) ||
      is.null(colnames(counts)) ||
      anyDuplicated(colnames(counts)) > 0L
  ) {
    stop(
      "Count-matrix identity failed for ",
      sample_accession
    )
  }

  if (is.null(reference_features)) {
    reference_features <- rownames(counts)
  } else if (
    identical(
      reference_features,
      rownames(counts)
    )
  ) {
    invisible(TRUE)
  } else if (
    setequal(
      reference_features,
      rownames(counts)
    )
  ) {
    counts <- counts[
      reference_features,
      ,
      drop = FALSE
    ]
  } else {
    stop(
      "Feature identities differ across donors at ",
      sample_accession
    )
  }

  if (
    !is.data.frame(cell_metadata) ||
      nrow(cell_metadata) != ncol(counts)
  ) {
    stop(
      "Cell metadata does not align with counts for ",
      sample_accession
    )
  }

  if (
    !is.null(rownames(cell_metadata)) &&
      setequal(
        rownames(cell_metadata),
        colnames(counts)
      )
  ) {
    cell_metadata <- cell_metadata[
      colnames(counts),
      ,
      drop = FALSE
    ]
  } else {
    rownames(cell_metadata) <- colnames(counts)
  }

  original_cell_barcode <- colnames(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  SummarizedExperiment::colData(sce)$sample_accession <-
    sample_accession

  set.seed(
    20260806L + sample_index
  )

  sce <- scDblFinder::scDblFinder(
    sce = sce,
    clusters = TRUE,
    samples = NULL,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(
      progressbar = FALSE
    )
  )

  scDblFinder_metadata <- as.data.frame(
    SummarizedExperiment::colData(sce)
  )

  if (
    !"scDblFinder.class" %in%
      colnames(scDblFinder_metadata) ||
      !"scDblFinder.score" %in%
        colnames(scDblFinder_metadata)
  ) {
    stop(
      "scDblFinder output columns are missing for ",
      sample_accession
    )
  }

  predicted_class <- tolower(
    as.character(
      scDblFinder_metadata$
        scDblFinder.class
    )
  )

  singlet_index <- which(
    predicted_class == "singlet"
  )

  doublet_index <- which(
    predicted_class == "doublet"
  )

  input_cell_n <- ncol(counts)
  singlet_cell_n <- length(singlet_index)
  doublet_cell_n <- length(doublet_index)

  if (
    singlet_cell_n < 1000L ||
      singlet_cell_n +
        doublet_cell_n != input_cell_n
  ) {
    stop(
      "Invalid singlet/doublet partition for ",
      sample_accession
    )
  }

  doublet_fraction <- (
    doublet_cell_n /
      input_cell_n
  )

  doublet_metadata <- data.frame(
    original_cell_barcode =
      original_cell_barcode,
    integrated_cell_barcode =
      paste0(
        sample_accession,
        "__",
        original_cell_barcode
      ),
    sample_accession =
      sample_accession,
    donor_id = donor_id,
    analysis_group =
      analysis_group,
    scDblFinder_class =
      predicted_class,
    scDblFinder_score =
      as.numeric(
        scDblFinder_metadata$
          scDblFinder.score
      ),
    stringsAsFactors = FALSE
  )

  optional_scDblFinder_columns <- grep(
    "^scDblFinder\\.",
    colnames(scDblFinder_metadata),
    value = TRUE
  )

  optional_scDblFinder_columns <- setdiff(
    optional_scDblFinder_columns,
    c(
      "scDblFinder.class",
      "scDblFinder.score"
    )
  )

  for (
    column_name in optional_scDblFinder_columns
  ) {
    doublet_metadata[[column_name]] <-
      scDblFinder_metadata[[column_name]]
  }

  write_csv_utf8(
    doublet_metadata,
    sample_doublet_path
  )

  singlet_counts <- counts[
    ,
    singlet_index,
    drop = FALSE
  ]

  singlet_metadata <- cell_metadata[
    singlet_index,
    ,
    drop = FALSE
  ]

  new_cell_names <- paste0(
    sample_accession,
    "__",
    original_cell_barcode[
      singlet_index
    ]
  )

  colnames(singlet_counts) <- new_cell_names
  rownames(singlet_metadata) <-
    new_cell_names

  singlet_metadata$
    original_cell_barcode <-
    original_cell_barcode[
      singlet_index
    ]

  singlet_metadata$sample_accession <-
    sample_accession

  singlet_metadata$donor_id <-
    donor_id

  singlet_metadata$analysis_group <-
    analysis_group

  singlet_metadata$scDblFinder_class <-
    predicted_class[
      singlet_index
    ]

  singlet_metadata$scDblFinder_score <-
    as.numeric(
      scDblFinder_metadata$
        scDblFinder.score[
          singlet_index
        ]
    )

  design_row <- design[
    design$sample_accession ==
      sample_accession,
    ,
    drop = FALSE
  ]

  for (
    design_column in colnames(design_row)
  ) {
    if (
      !design_column %in%
        colnames(singlet_metadata)
    ) {
      singlet_metadata[[design_column]] <-
        rep(
          as.character(
            design_row[[design_column]][1]
          ),
          singlet_cell_n
        )
    }
  }

  singlet_object <- Seurat::CreateSeuratObject(
    counts = singlet_counts,
    project = sample_accession,
    assay = "RNA",
    min.cells = 0,
    min.features = 0,
    meta.data = singlet_metadata
  )

  SeuratObject::Project(
    singlet_object
  ) <- sample_accession

  if (
    any(
      grepl(
        "^MT-",
        rownames(singlet_object),
        ignore.case = FALSE
      )
    )
  ) {
    singlet_object[["percent.mt"]] <-
      Seurat::PercentageFeatureSet(
      singlet_object,
      pattern = "^MT-"
    )
  }

  singlet_object@misc$
    GSE183837_sample_checkpoint <- list(
      sample_accession =
        sample_accession,
      donor_id = donor_id,
      analysis_group =
        analysis_group,
      input_cell_n =
        input_cell_n,
      singlet_cell_n =
        singlet_cell_n,
      doublet_cell_n =
        doublet_cell_n,
      doublet_fraction =
        doublet_fraction,
      scDblFinder_version =
        safe_package_version(
          "scDblFinder"
        ),
      random_seed =
        20260806L + sample_index
    )

  saveRDS(
    singlet_object,
    sample_checkpoint_path,
    compress = "gzip"
  )

  checkpoint_md5 <- unname(
    tools::md5sum(
      sample_checkpoint_path
    )
  )

  elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        sample_start,
        units = "mins"
      )
    ),
    2
  )

  sample_status <- data.frame(
    sample_accession =
      sample_accession,
    donor_id = donor_id,
    analysis_group =
      analysis_group,
    status = "COMPLETE",
    input_cell_n =
      input_cell_n,
    singlet_cell_n =
      singlet_cell_n,
    doublet_cell_n =
      doublet_cell_n,
    doublet_fraction =
      doublet_fraction,
    high_doublet_review_flag =
      doublet_fraction > 0.20,
    singlet_checkpoint_size_MB =
      round(
        file.info(
          sample_checkpoint_path
        )$size / 1024^2,
        2
      ),
    singlet_checkpoint_md5 =
      checkpoint_md5,
    elapsed_minutes =
      elapsed_minutes,
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    sample_status,
    sample_status_path
  )

  sample_status_rows[[sample_index]] <-
    sample_status

  doublet_metadata_paths <- c(
    doublet_metadata_paths,
    sample_doublet_path
  )

  log_event(
    "DONE donor",
    sample_index,
    "/",
    nrow(design),
    sample_accession,
    "input=",
    input_cell_n,
    "doublets=",
    doublet_cell_n,
    "singlets=",
    singlet_cell_n,
    "doublet_fraction=",
    round(doublet_fraction, 4),
    "elapsed_min=",
    elapsed_minutes
  )

  rm(
    input_object,
    counts,
    cell_metadata,
    sce,
    scDblFinder_metadata,
    doublet_metadata,
    singlet_counts,
    singlet_metadata,
    singlet_object
  )

  gc()
}

sample_status_all <- do.call(
  rbind,
  sample_status_rows
)

sample_status_all <- sample_status_all[
  match(
    design$sample_accession,
    sample_status_all$sample_accession
  ),
  ,
  drop = FALSE
]

if (
  nrow(sample_status_all) != 9L ||
    anyNA(
      sample_status_all$sample_accession
    ) ||
    !all(
      sample_status_all$status ==
        "COMPLETE"
    ) ||
    sum(
      sample_status_all$input_cell_n
    ) != 50997L
) {
  stop(
    "Sample-wise scDblFinder completion audit failed."
  )
}

write_csv_utf8(
  sample_status_all,
  file.path(
    table_root,
    "01_doublet_summary.csv"
  )
)

doublet_metadata_all <- data.table::rbindlist(
  lapply(
    doublet_metadata_paths,
    function(path_value) {
      data.table::fread(
        path_value,
        data.table = FALSE,
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

write_csv_utf8(
  as.data.frame(doublet_metadata_all),
  file.path(
    table_root,
    "02_doublet_metadata_all.csv"
  )
)

total_singlets <- sum(
  sample_status_all$singlet_cell_n
)

total_doublets <- sum(
  sample_status_all$doublet_cell_n
)


# ----------------------------------------------------------------------
# 6. Merge the nine singlet checkpoints
# ----------------------------------------------------------------------

if (file.exists(final_status_path)) {
  final_status_existing <- safe_read_csv(
    final_status_path
  )
} else {
  final_status_existing <- NULL
}

final_checkpoint_complete <- (
  file.exists(final_checkpoint_path) &&
    !is.null(final_status_existing) &&
    nrow(final_status_existing) == 1L &&
    final_status_existing$status[[1]] ==
      "COMPLETE"
)

if (final_checkpoint_complete) {
  integrated_object <- readRDS(
    final_checkpoint_path
  )

  if (
    !inherits(integrated_object, "Seurat") ||
      ncol(integrated_object) !=
        total_singlets ||
      length(
        unique(
          integrated_object$
            sample_accession
        )
      ) != 9L
  ) {
    stop(
      "Existing final checkpoint failed validation."
    )
  }

  log_event(
    "SKIP integration; valid final checkpoint found",
    "cells=",
    ncol(integrated_object)
  )
} else {
  if (file.exists(merged_counts_checkpoint)) {
    merged_object <- tryCatch(
      readRDS(
        merged_counts_checkpoint
      ),
      error = function(e) NULL
    )

    merged_checkpoint_valid <- (
      inherits(
        merged_object,
        "Seurat"
      ) &&
        ncol(merged_object) ==
          total_singlets &&
        length(
          unique(
            merged_object$
              sample_accession
          )
        ) == 9L
    )

    if (!merged_checkpoint_valid) {
      rm(merged_object)
      merged_object <- NULL
    }
  } else {
    merged_object <- NULL
  }

  if (is.null(merged_object)) {
    log_event(
      "START merging nine singlet checkpoints"
    )

    object_list <- vector(
      "list",
      nrow(design)
    )

    names(object_list) <-
      as.character(
        design$sample_accession
      )

    for (
      sample_index in seq_len(
        nrow(design)
      )
    ) {
      sample_accession <- as.character(
        design$sample_accession[
          sample_index
        ]
      )

      checkpoint_path <- file.path(
        sample_checkpoint_root,
        sample_accession,
        "singlet.rds"
      )

      object_list[[sample_index]] <-
        readRDS(checkpoint_path)
    }

    merged_object <- merge(
      x = object_list[[1]],
      y = object_list[-1],
      merge.data = FALSE,
      merge.dr = FALSE,
      collapse = FALSE,
      project = dataset_id
    )

    if (
      ncol(merged_object) !=
        total_singlets ||
        length(
          unique(
            merged_object$
              sample_accession
          )
        ) != 9L
    ) {
      stop(
        "Merged singlet checkpoint failed cell/sample identity."
      )
    }

    saveRDS(
      merged_object,
      merged_counts_checkpoint,
      compress = "gzip"
    )

    log_event(
      "DONE merging singlets",
      "cells=",
      ncol(merged_object),
      "features=",
      nrow(merged_object),
      "checkpoint_MB=",
      round(
        file.info(
          merged_counts_checkpoint
        )$size / 1024^2,
        2
      )
    )

    rm(object_list)
    gc()
  } else {
    log_event(
      "RESUME from merged singlet checkpoint",
      "cells=",
      ncol(merged_object)
    )
  }


# ----------------------------------------------------------------------
# 7. Normalize, PCA, Harmony, UMAP and clustering
# ----------------------------------------------------------------------

  first_preQC_path <- file.path(
    preQC_root,
    paste0(
      as.character(
        design$sample_accession[1]
      ),
      "_preQC.rds"
    )
  )

  first_preQC_object <- readRDS(
    first_preQC_path
  )

  original_features <- rownames(
    first_preQC_object$counts
  )

  expected_Seurat_features <- gsub(
    "_",
    "-",
    original_features,
    fixed = TRUE
  )

  expected_Seurat_features <- make.unique(
    expected_Seurat_features
  )

  current_features <- rownames(
    merged_object
  )

  feature_name_audit <- data.frame(
    metric = c(
      "original_feature_n",
      "merged_feature_n",
      "original_feature_names_with_underscore_n",
      "expected_transformed_feature_unique_n",
      "merged_feature_unique_n",
      "expected_transformation_matches_merged"
    ),
    value = c(
      length(original_features),
      length(current_features),
      sum(
        grepl(
          "_",
          original_features,
          fixed = TRUE
        )
      ),
      length(
        unique(expected_Seurat_features)
      ),
      length(
        unique(current_features)
      ),
      identical(
        expected_Seurat_features,
        current_features
      )
    ),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    feature_name_audit,
    file.path(
      table_root,
      "00_feature_name_audit.csv"
    )
  )

  if (
    length(current_features) !=
      length(original_features) ||
      anyDuplicated(current_features) > 0L
  ) {
    stop(
      "Feature-name transformation produced a dimension or uniqueness problem."
    )
  }

  rm(
    first_preQC_object,
    original_features,
    expected_Seurat_features,
    current_features
  )
  gc()

  prejoin_checkpoint_valid <- FALSE
  prejoin_object <- NULL

  if (file.exists(prejoin_checkpoint_path)) {
    prejoin_object <- tryCatch(
      readRDS(prejoin_checkpoint_path),
      error = function(e) NULL
    )

    if (
      inherits(prejoin_object, "Seurat") &&
        ncol(prejoin_object) == total_singlets &&
        length(
          unique(
            prejoin_object$sample_accession
          )
        ) == 9L &&
        all(
          c(
            "pca",
            "harmony",
            "umap.harmony"
          ) %in%
            names(prejoin_object@reductions)
        ) &&
        "harmony_snn_res.0.4" %in%
          colnames(prejoin_object[[]])
    ) {
      prejoin_checkpoint_valid <- TRUE
    }
  }

  integration_start <- Sys.time()

  if (prejoin_checkpoint_valid) {
    merged_object <- prejoin_object
    primary_cluster_column <-
      "harmony_snn_res.0.4"

    merged_object$seurat_clusters <-
      as.factor(
        merged_object[[primary_cluster_column]][, 1]
      )

    SeuratObject::Idents(
      merged_object
    ) <- merged_object$
      seurat_clusters

    log_event(
      "RESUME from Harmony prejoin checkpoint",
      "cells=",
      ncol(merged_object),
      "clusters_res0.4=",
      length(
        unique(
          merged_object$seurat_clusters
        )
      )
    )
  } else {
    rm(prejoin_object)
    gc()

  SeuratObject::DefaultAssay(
    merged_object
  ) <- "RNA"

  merged_object <- Seurat::NormalizeData(
    object = merged_object,
    assay = "RNA",
    normalization.method =
      "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  merged_object <- Seurat::FindVariableFeatures(
    object = merged_object,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  if (
    length(
      SeuratObject::VariableFeatures(
        merged_object
      )
    ) < 2000L
  ) {
    stop(
      "Too few variable features were recovered."
    )
  }

  merged_object <- Seurat::ScaleData(
    object = merged_object,
    assay = "RNA",
    features = SeuratObject::VariableFeatures(
      merged_object
    ),
    verbose = FALSE
  )

  merged_object <- Seurat::RunPCA(
    object = merged_object,
    assay = "RNA",
    features = SeuratObject::VariableFeatures(
      merged_object
    ),
    npcs = 50,
    seed.use = 20260806,
    verbose = FALSE
  )

  merged_object <- run_harmony_Seurat_compatible(
    object = merged_object,
    group_variable = "sample_accession",
    dimensions = 1:30,
    reduction_name = "pca",
    output_reduction = "harmony",
    verbose = FALSE
  )

  merged_object <- Seurat::RunUMAP(
    object = merged_object,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony",
    reduction.key = "UMAPH_",
    seed.use = 20260806,
    verbose = FALSE
  )

  merged_object <- Seurat::FindNeighbors(
    object = merged_object,
    reduction = "harmony",
    dims = 1:30,
    graph.name = c(
      "harmony_nn",
      "harmony_snn"
    ),
    verbose = FALSE
  )

  merged_object <- Seurat::FindClusters(
    object = merged_object,
    graph.name = "harmony_snn",
    resolution = c(
      0.2,
      0.4,
      0.6
    ),
    algorithm = 1,
    random.seed = 20260806,
    verbose = FALSE
  )

  primary_cluster_column <-
    "harmony_snn_res.0.4"

  if (
    !primary_cluster_column %in%
      colnames(merged_object[[]])
  ) {
    candidate_primary_columns <- grep(
      "res\\.0\\.4$",
      colnames(merged_object[[]]),
      value = TRUE
    )

    if (length(candidate_primary_columns) != 1L) {
      stop(
        "Could not resolve the primary 0.4 cluster column."
      )
    }

    primary_cluster_column <-
      candidate_primary_columns[1]
  }

  merged_object$seurat_clusters <-
    as.factor(
      merged_object[[primary_cluster_column]][, 1]
    )

  SeuratObject::Idents(
    merged_object
  ) <- merged_object$
    seurat_clusters

  merged_object@misc$
    GSE183837_prejoin_checkpoint <- list(
      script_version = script_version,
      cell_n = ncol(merged_object),
      sample_n = length(
        unique(
          merged_object$sample_accession
        )
      ),
      primary_cluster_column =
        primary_cluster_column,
      primary_cluster_n = length(
        unique(
          merged_object$seurat_clusters
        )
      ),
      reductions =
        names(merged_object@reductions),
      checkpoint_role =
        "POST_HARMONY_UMAP_CLUSTERING_PRE_JOINLAYERS"
    )

  saveRDS(
    merged_object,
    prejoin_checkpoint_path,
    compress = "gzip"
  )

  log_event(
    "SAVED Harmony prejoin checkpoint",
    "checkpoint_MB=",
    round(
      file.info(
        prejoin_checkpoint_path
      )$size / 1024^2,
      2
    )
  )
  }

  data_layers_before_join <- SeuratObject::Layers(
    merged_object,
    assay = "RNA",
    search = "^data"
  )

  count_layers_before_join <- SeuratObject::Layers(
    merged_object,
    assay = "RNA",
    search = "^counts"
  )

  merged_object <- join_normalized_data_layers_only(
    object = merged_object,
    assay_name = "RNA"
  )

  data_layers_after_join <- SeuratObject::Layers(
    merged_object,
    assay = "RNA",
    search = "^data"
  )

  count_layers_after_join <- SeuratObject::Layers(
    merged_object,
    assay = "RNA",
    search = "^counts"
  )

  join_audit <- data.frame(
    metric = c(
      "data_layers_before_n",
      "data_layers_before",
      "data_layers_after_n",
      "data_layers_after",
      "count_layers_before_n",
      "count_layers_before",
      "count_layers_after_n",
      "count_layers_after",
      "count_layer_identity_preserved"
    ),
    value = c(
      length(data_layers_before_join),
      collapse_unique(data_layers_before_join),
      length(data_layers_after_join),
      collapse_unique(data_layers_after_join),
      length(count_layers_before_join),
      collapse_unique(count_layers_before_join),
      length(count_layers_after_join),
      collapse_unique(count_layers_after_join),
      identical(
        count_layers_before_join,
        count_layers_after_join
      )
    ),
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    join_audit,
    file.path(
      table_root,
      "00_JoinLayers_audit.csv"
    )
  )

  scale_layers_before <- SeuratObject::Layers(
    merged_object,
    assay = "RNA",
    search = "^scale\\.data"
  )

  scale_data_removed <- FALSE

  if (
    "scale.data" %in%
      scale_layers_before
  ) {
    assay_object <- merged_object[["RNA"]]

    scale_data_removed <- tryCatch(
      {
        assay_object$scale.data <- NULL
        merged_object[["RNA"]] <-
          assay_object
        TRUE
      },
      error = function(e) FALSE
    )
  }

  integrated_object <- merged_object

  integrated_object@misc$
    GSE183837_integration <- list(
      integration_version =
        integration_version,
      script_version =
        script_version,
      Harmony_API_policy =
        "MINIMAL_CURRENT_API_NO_ASSAY_USE_NO_PROJECT_DIM",
      JoinLayers_policy =
        "LAYERS_DATA_TO_SINGLE_DATA_COUNTS_UNCHANGED",
      QC_metadata_policy =
        "UNIQUE_COLUMN_NAMES_FOR_EXPORT_AND_GGPLOT",
      donor_gate_policy =
        "UNIQUE_DONOR_GROUP_PAIRS",
      supersedes_failed_executions = c(
        "V1_0_HARMONY_ASSAY_USE_FAILURE",
        "V1_1_JOINLAYERS_LENGTH_MISMATCH",
        "V1_2_GGPLOT_DUPLICATE_METADATA_COLUMNS",
        "V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS"
      ),
      upstream_preflight_gate =
        preflight_gate$
          gate_status[[1]],
      sample_n = 9L,
      preQC_cell_n = 50997L,
      singlet_cell_n =
        total_singlets,
      predicted_doublet_n =
        total_doublets,
      Harmony_batch_variable =
        "sample_accession",
      Harmony_dimensions =
        1:30,
      cluster_resolutions =
        c(0.2, 0.4, 0.6),
      primary_resolution =
        0.4,
      primary_cluster_column =
        primary_cluster_column,
      scale_data_removed =
        scale_data_removed,
      marker_analysis_authorized =
        FALSE,
      result_lock_authorized =
        FALSE
    )

  gc()

  saveRDS(
    integrated_object,
    final_checkpoint_path,
    compress = "gzip"
  )

  final_checkpoint_md5 <- unname(
    tools::md5sum(
      final_checkpoint_path
    )
  )

  integration_elapsed_minutes <- round(
    as.numeric(
      difftime(
        Sys.time(),
        integration_start,
        units = "mins"
      )
    ),
    2
  )

  final_status <- data.frame(
    dataset_id = dataset_id,
    integration_version =
      integration_version,
    status = "COMPLETE",
    feature_n = nrow(integrated_object),
    singlet_cell_n =
      ncol(integrated_object),
    sample_n = length(
      unique(
        integrated_object$
          sample_accession
      )
    ),
    donor_n = length(
      unique(
        integrated_object$
          donor_id
      )
    ),
    group_n = length(
      unique(
        integrated_object$
          analysis_group
      )
    ),
    primary_cluster_column =
      primary_cluster_column,
    primary_cluster_n = length(
      unique(
        integrated_object$
          seurat_clusters
      )
    ),
    final_checkpoint_size_GB =
      round(
        file.info(
          final_checkpoint_path
        )$size / 1024^3,
        3
      ),
    final_checkpoint_md5 =
      final_checkpoint_md5,
    integration_elapsed_minutes =
      integration_elapsed_minutes,
    stringsAsFactors = FALSE
  )

  write_csv_utf8(
    final_status,
    final_status_path
  )

  log_event(
    "DONE Harmony integration",
    "cells=",
    ncol(integrated_object),
    "clusters_res0.4=",
    length(
      unique(
        integrated_object$
          seurat_clusters
      )
    ),
    "elapsed_min=",
    integration_elapsed_minutes
  )

  rm(merged_object)
  gc()
}


# ----------------------------------------------------------------------
# 8. Execution-level technical QC tables
# ----------------------------------------------------------------------

if (!exists("primary_cluster_column")) {
  primary_cluster_column <-
    integrated_object@misc$
      GSE183837_integration$
        primary_cluster_column
}

RNA_layers <- SeuratObject::Layers(
  integrated_object,
  assay = "RNA",
  search = NA
)

count_layers <- grep(
  "^counts",
  RNA_layers,
  value = TRUE
)

data_layers_final <- grep(
  "^data",
  RNA_layers,
  value = TRUE
)

scale_layers_final <- grep(
  "^scale\\.data",
  RNA_layers,
  value = TRUE
)

layer_audit <- data.frame(
  layer = RNA_layers,
  layer_type = ifelse(
    grepl("^counts", RNA_layers),
    "counts",
    ifelse(
      grepl("^data", RNA_layers),
      "data",
      ifelse(
        grepl(
          "^scale\\.data",
          RNA_layers
        ),
        "scale.data",
        "other"
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  layer_audit,
  file.path(
    table_root,
    "03_layer_audit.csv"
  )
)

metadata <- integrated_object[[]]

UMAP_embedding <- SeuratObject::Embeddings(
  integrated_object,
  reduction = "umap.harmony"
)

UMAP_metadata <- prepare_unique_UMAP_metadata(
  metadata = metadata,
  UMAP_embedding = UMAP_embedding,
  audit_path = file.path(
    table_root,
    "00_metadata_name_audit.csv"
  )
)

write_csv_utf8(
  UMAP_metadata,
  file.path(
    table_root,
    "04_UMAP_metadata.csv"
  )
)

sample_cell_summary <- aggregate(
  cell_barcode ~
    sample_accession +
    donor_id +
    analysis_group,
  data = UMAP_metadata,
  FUN = length
)

colnames(
  sample_cell_summary
)[
  colnames(sample_cell_summary) ==
    "cell_barcode"
] <- "singlet_cell_n"

sample_cell_summary <- merge(
  sample_cell_summary,
  sample_status_all[
    ,
    c(
      "sample_accession",
      "input_cell_n",
      "doublet_cell_n",
      "doublet_fraction"
    )
  ],
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

sample_cell_summary <- sample_cell_summary[
  match(
    design$sample_accession,
    sample_cell_summary$
      sample_accession
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  sample_cell_summary,
  file.path(
    table_root,
    "05_sample_cell_summary.csv"
  )
)

cluster_sample_counts <- aggregate(
  cell_barcode ~
    seurat_clusters +
    sample_accession +
    analysis_group,
  data = UMAP_metadata,
  FUN = length
)

colnames(
  cluster_sample_counts
)[
  colnames(cluster_sample_counts) ==
    "cell_barcode"
] <- "cell_n"

write_csv_utf8(
  cluster_sample_counts,
  file.path(
    table_root,
    "06_cluster_sample_counts.csv"
  )
)

cluster_summary <- aggregate(
  cell_barcode ~
    seurat_clusters,
  data = UMAP_metadata,
  FUN = length
)

colnames(
  cluster_summary
)[
  colnames(cluster_summary) ==
    "cell_barcode"
] <- "cell_n"

cluster_summary$cell_fraction <- (
  cluster_summary$cell_n /
    sum(cluster_summary$cell_n)
)

cluster_summary <- cluster_summary[
  order(
    as.integer(
      as.character(
        cluster_summary$
          seurat_clusters
      )
    )
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  cluster_summary,
  file.path(
    table_root,
    "07_cluster_summary.csv"
  )
)

resolution_columns <- grep(
  "harmony_snn_res\\.",
  colnames(metadata),
  value = TRUE
)

resolution_summary_rows <- list()

for (
  resolution_column in
    resolution_columns
) {
  resolution_summary_rows[[
    length(resolution_summary_rows) + 1L
  ]] <- data.frame(
    resolution_column =
      resolution_column,
    cluster_n = length(
      unique(
        metadata[[resolution_column]]
      )
    ),
    minimum_cluster_n = min(
      table(
        metadata[[resolution_column]]
      )
    ),
    maximum_cluster_n = max(
      table(
        metadata[[resolution_column]]
      )
    ),
    stringsAsFactors = FALSE
  )
}

resolution_summary <- do.call(
  rbind,
  resolution_summary_rows
)

write_csv_utf8(
  resolution_summary,
  file.path(
    table_root,
    "08_resolution_summary.csv"
  )
)

PCA_embedding <- SeuratObject::Embeddings(
  integrated_object,
  reduction = "pca"
)

Harmony_embedding <- SeuratObject::Embeddings(
  integrated_object,
  reduction = "harmony"
)

reduction_audit <- data.frame(
  reduction = c(
    "pca",
    "harmony",
    "umap.harmony"
  ),
  cell_n = c(
    nrow(PCA_embedding),
    nrow(Harmony_embedding),
    nrow(UMAP_embedding)
  ),
  dimension_n = c(
    ncol(PCA_embedding),
    ncol(Harmony_embedding),
    ncol(UMAP_embedding)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  reduction_audit,
  file.path(
    table_root,
    "09_reduction_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Execution-level QC figures
# ----------------------------------------------------------------------

if (anyDuplicated(colnames(UMAP_metadata)) > 0L) {
  stop("UMAP plotting metadata contains duplicate column names.")
}

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

sample_status_plot <- sample_status_all[
  ,
  !duplicated(
    colnames(sample_status_all)
  ),
  drop = FALSE
]

if (anyDuplicated(colnames(sample_status_plot)) > 0L) {
  stop("Sample-status plotting data contains duplicate column names.")
}

sample_status_plot$sample_accession <-
  factor(
    sample_status_plot$
      sample_accession,
    levels = design$sample_accession
  )

plot_doublet <- ggplot2::ggplot(
  sample_status_plot,
  ggplot2::aes(
    x = sample_accession,
    y = doublet_fraction,
    shape = analysis_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0.20,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  ggplot2::labs(
    title = "GSE183837 donor-wise predicted doublet fractions",
    subtitle = "scDblFinder was run independently within each donor; dashed line is a review threshold",
    x = "Sample",
    y = "Predicted doublet fraction",
    shape = "Group"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  )

save_plot_pair(
  plot_doublet,
  "QC1_doublet_fraction",
  width = 9,
  height = 5
)

plot_singlets <- ggplot2::ggplot(
  sample_status_plot,
  ggplot2::aes(
    x = sample_accession,
    y = singlet_cell_n,
    shape = analysis_group
  )
) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "GSE183837 retained singlet cells",
    x = "Sample",
    y = "Retained singlet cells",
    shape = "Group"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  )

save_plot_pair(
  plot_singlets,
  "QC2_singlet_counts",
  width = 9,
  height = 5
)

plot_UMAP_sample <- ggplot2::ggplot(
  UMAP_metadata,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = sample_accession
  )
) +
  ggplot2::geom_point(
    size = 0.12,
    alpha = 0.65,
    stroke = 0
  ) +
  ggplot2::labs(
    title = "GSE183837 Harmony UMAP by sample",
    x = "Harmony UMAP 1",
    y = "Harmony UMAP 2",
    color = "Sample"
  ) +
  theme_QC

save_plot_pair(
  plot_UMAP_sample,
  "QC3_UMAP_sample",
  width = 10,
  height = 7
)

plot_UMAP_group <- ggplot2::ggplot(
  UMAP_metadata,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = analysis_group
  )
) +
  ggplot2::geom_point(
    size = 0.12,
    alpha = 0.65,
    stroke = 0
  ) +
  ggplot2::labs(
    title = "GSE183837 Harmony UMAP by clinical group",
    subtitle = "Visualization only; no cell-level group inference",
    x = "Harmony UMAP 1",
    y = "Harmony UMAP 2",
    color = "Group"
  ) +
  theme_QC

save_plot_pair(
  plot_UMAP_group,
  "QC4_UMAP_group",
  width = 9,
  height = 7
)

cluster_centroids <- aggregate(
  cbind(
    UMAPH_1,
    UMAPH_2
  ) ~ seurat_clusters,
  data = UMAP_metadata,
  FUN = median
)

plot_UMAP_cluster <- ggplot2::ggplot(
  UMAP_metadata,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    color = seurat_clusters
  )
) +
  ggplot2::geom_point(
    size = 0.12,
    alpha = 0.65,
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
    title = "GSE183837 Harmony UMAP by primary clusters",
    subtitle = "Primary resolution = 0.4; biological identities are not assigned",
    x = "Harmony UMAP 1",
    y = "Harmony UMAP 2",
    color = "Cluster"
  ) +
  theme_QC

save_plot_pair(
  plot_UMAP_cluster,
  "QC5_UMAP_cluster",
  width = 10,
  height = 7
)


# ----------------------------------------------------------------------
# 10. Integration execution Gate
# ----------------------------------------------------------------------

free_memory_after_GB <- get_free_memory_gb()
free_disk_after_GB <- get_free_disk_gb()

group_cell_counts_final <- table(
  factor(
    as.character(
      integrated_object$
        analysis_group
    ),
    levels = c(
      "CONTROL",
      "RIF"
    )
  )
)

donor_group_table <- unique(
  data.frame(
    donor_id = as.character(
      integrated_object$donor_id
    ),
    analysis_group = as.character(
      integrated_object$analysis_group
    ),
    stringsAsFactors = FALSE
  )
)

donor_group_table <- donor_group_table[
  order(donor_group_table$donor_id),
  ,
  drop = FALSE
]

if (
  nrow(donor_group_table) != 9L ||
    anyDuplicated(
      donor_group_table$donor_id
    ) > 0L ||
    anyNA(donor_group_table$donor_id) ||
    anyNA(donor_group_table$analysis_group)
) {
  stop(
    "Final donor/group identity table is not one row per nine unique donors."
  )
}

group_donor_counts_final <- table(
  factor(
    donor_group_table$analysis_group,
    levels = c(
      "CONTROL",
      "RIF"
    )
  )
)

group_count_audit <- data.frame(
  metric = c(
    "unique_donor_n",
    "Control_donor_n",
    "RIF_donor_n",
    "Control_cell_n",
    "RIF_cell_n",
    "donor_group_design_is_3_Control_6_RIF"
  ),
  value = c(
    nrow(donor_group_table),
    as.integer(
      group_donor_counts_final["CONTROL"]
    ),
    as.integer(
      group_donor_counts_final["RIF"]
    ),
    as.integer(
      group_cell_counts_final["CONTROL"]
    ),
    as.integer(
      group_cell_counts_final["RIF"]
    ),
    identical(
      as.integer(
        group_donor_counts_final
      ),
      c(3L, 6L)
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  donor_group_table,
  file.path(
    table_root,
    "00_donor_group_identity.csv"
  )
)

write_csv_utf8(
  group_count_audit,
  file.path(
    table_root,
    "00_group_count_audit.csv"
  )
)

review_flags <- c(
  if (
    any(
      sample_status_all$
        high_doublet_review_flag
    )
  ) {
    "SAMPLE_DOUBLET_FRACTION_ABOVE_20_PERCENT"
  } else {
    character()
  },
  if (
    any(
      cluster_summary$cell_n < 100L
    )
  ) {
    "PRIMARY_CLUSTER_BELOW_100_CELLS"
  } else {
    character()
  },
  if (
    length(scale_layers_final) > 0L
  ) {
    "SCALE_DATA_LAYER_RETAINED"
  } else {
    character()
  }
)

review_flags <- unique(
  review_flags[
    nzchar(review_flags)
  ]
)

gate_checks <- data.frame(
  check = c(
    "upstream_preflight_authorizes_integration",
    "nine_donor_doublet_jobs_complete",
    "preQC_cells_partitioned_exactly",
    "all_donors_retain_at_least_1000_singlets",
    "final_cells_equal_singlet_sum",
    "final_sample_n_is_9",
    "final_donor_n_is_9",
    "final_group_counts_are_3_Control_6_RIF",
    "PCA_has_50_dimensions",
    "Harmony_has_at_least_30_dimensions",
    "UMAP_has_2_dimensions",
    "cluster_resolutions_0.2_0.4_0.6_present",
    "primary_resolution_0.4_present",
    "nine_split_count_layers_retained",
    "one_joined_data_layer_present",
    "QC_UMAP_metadata_columns_unique",
    "final_checkpoint_exists",
    "final_checkpoint_md5_available",
    "free_disk_after_run_at_least_10GB_or_unknown"
  ),
  pass = c(
    preflight_gate$
      integration_authorized[[1]],
    nrow(sample_status_all) == 9L &&
      all(
        sample_status_all$status ==
          "COMPLETE"
      ),
    sum(
      sample_status_all$
        input_cell_n
    ) == 50997L &&
      sum(
        sample_status_all$
          singlet_cell_n +
          sample_status_all$
            doublet_cell_n
      ) == 50997L,
    all(
      sample_status_all$
        singlet_cell_n >= 1000L
    ),
    ncol(integrated_object) ==
      total_singlets,
    length(
      unique(
        integrated_object$
          sample_accession
      )
    ) == 9L,
    length(
      unique(
        integrated_object$
          donor_id
      )
    ) == 9L,
    identical(
      as.integer(
        group_donor_counts_final
      ),
      c(3L, 6L)
    ),
    ncol(PCA_embedding) == 50L,
    ncol(Harmony_embedding) >= 30L,
    ncol(UMAP_embedding) == 2L,
    all(
      c(
        "harmony_snn_res.0.2",
        "harmony_snn_res.0.4",
        "harmony_snn_res.0.6"
      ) %in%
        colnames(metadata)
    ),
    primary_cluster_column %in%
      colnames(metadata),
    length(count_layers) == 9L,
    length(data_layers_final) == 1L,
    anyDuplicated(
      colnames(UMAP_metadata)
    ) == 0L,
    file.exists(final_checkpoint_path),
    nzchar(
      unname(
        tools::md5sum(
          final_checkpoint_path
        )
      )
    ),
    is.na(free_disk_after_GB) ||
      free_disk_after_GB >= 10
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
  all_gate_checks_pass &&
    !length(review_flags)
) {
  "PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED"
} else if (
  all_gate_checks_pass
) {
  "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS"
} else {
  "HOLD_GSE183837_INTEGRATION_EXECUTION_FAILED"
}

integration_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  gate_status = gate_status,
  preQC_cell_n = 50997L,
  predicted_doublet_n =
    total_doublets,
  retained_singlet_n =
    total_singlets,
  predicted_doublet_fraction =
    total_doublets / 50997L,
  sample_n = 9L,
  Control_donor_n =
    as.integer(
      group_donor_counts_final["CONTROL"]
    ),
  RIF_donor_n =
    as.integer(
      group_donor_counts_final["RIF"]
    ),
  Control_cell_n =
    as.integer(
      group_cell_counts_final["CONTROL"]
    ),
  RIF_cell_n =
    as.integer(
      group_cell_counts_final["RIF"]
    ),
  feature_n =
    nrow(integrated_object),
  primary_cluster_n =
    length(
      unique(
        integrated_object$
          seurat_clusters
      )
    ),
  count_layer_n =
    length(count_layers),
  data_layer_n =
    length(data_layers_final),
  review_flags =
    collapse_unique(review_flags),
  detailed_integration_QC_authorized =
    all_gate_checks_pass,
  marker_analysis_authorized =
    FALSE,
  biological_annotation_locked =
    FALSE,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

gate_path <- file.path(
  run_qc_root,
  "GSE183837_INTEGRATION_GATE.csv"
)

write_csv_utf8(
  integration_gate,
  gate_path
)

card_lines <- c(
  "# STEP04E1J GSE183837 Integration Execution Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0(
    "- Integration version: ",
    integration_version
  ),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Pre-QC cells: ",
    format(50997L, big.mark = ",")
  ),
  paste0(
    "- Predicted doublets: ",
    format(total_doublets, big.mark = ","),
    " (",
    round(
      100 * total_doublets / 50997L,
      2
    ),
    "%)"
  ),
  paste0(
    "- Retained singlets: ",
    format(total_singlets, big.mark = ",")
  ),
  paste0(
    "- Primary resolution-0.4 clusters: ",
    length(
      unique(
        integrated_object$
          seurat_clusters
      )
    )
  ),
  paste0(
    "- Raw count layers retained: ",
    length(count_layers)
  ),
  paste0(
    "- Joined normalized data layers: ",
    length(data_layers_final)
  ),
  paste0(
    "- Review flags: ",
    if (
      length(review_flags)
    ) {
      collapse_unique(review_flags)
    } else {
      "NONE"
    }
  ),
  paste0(
    "- Detailed integration QC authorized: ",
    all_gate_checks_pass
  ),
  "- Marker analysis authorized: FALSE",
  "- Biological annotation locked: FALSE",
  "- Disease pseudobulk testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "GSE183837_INTEGRATION_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 integration execution interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Completed operations",
  "",
  "- scDblFinder was run independently within each donor/capture.",
  "- Predicted doublets were excluded before integration.",
  "- Nine donor-level singlet checkpoints were retained for reproducibility and resume support.",
  "- PCA and Harmony were performed on the merged singlet data.",
  "- UMAP and clustering were generated from Harmony dimensions 1-30.",
  "- Resolution 0.4 is the primary technical cluster solution.",
  "- Nine split raw-count layers were retained and normalized data layers were joined.",
  "",
  "## Not authorized yet",
  "",
  "- Marker extraction is deferred until detailed integration QC passes.",
  "- Broad cell-type identities are not assigned or locked.",
  "- No cell-level Control-versus-RIF inference is permitted.",
  "- Donor-level disease testing is deferred until annotation and pseudobulk QC pass.",
  "",
  "## Next step",
  "",
  "- Quantify sample mixing, cluster dominance, rare-cluster composition and group/sample confounding.",
  "- Review UMAP and cluster-resolution stability.",
  "- Authorize memory-safe marker extraction only after the detailed QC Gate."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "GSE183837_INTEGRATION_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 11. Manifests, sessionInfo and ZIP
# ----------------------------------------------------------------------

run_end <- Sys.time()

package_table <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
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

resource_summary <- data.frame(
  metric = c(
    "free_memory_before_GB",
    "free_memory_after_GB",
    "free_C_drive_before_GB",
    "free_C_drive_after_GB",
    "merged_checkpoint_size_GB",
    "prejoin_checkpoint_size_GB",
    "final_checkpoint_size_GB"
  ),
  value = c(
    free_memory_before_GB,
    free_memory_after_GB,
    free_disk_before_GB,
    free_disk_after_GB,
    if (
      file.exists(
        merged_counts_checkpoint
      )
    ) {
      round(
        file.info(
          merged_counts_checkpoint
        )$size / 1024^3,
        3
      )
    } else {
      NA_real_
    },
    if (
      file.exists(
        prejoin_checkpoint_path
      )
    ) {
      round(
        file.info(
          prejoin_checkpoint_path
        )$size / 1024^3,
        3
      )
    } else {
      NA_real_
    },
    round(
      file.info(
        final_checkpoint_path
      )$size / 1024^3,
      3
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resource_summary,
  file.path(
    run_qc_root,
    "02_resource_summary.csv"
  )
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

checkpoint_manifest <- data.frame(
  artifact = c(
    "merged_singlet_checkpoint",
    "Harmony_prejoin_checkpoint",
    "final_integrated_checkpoint"
  ),
  path = normalize_slashes(c(
    merged_counts_checkpoint,
    prejoin_checkpoint_path,
    final_checkpoint_path
  )),
  size_bytes = file.info(c(
    merged_counts_checkpoint,
    prejoin_checkpoint_path,
    final_checkpoint_path
  ))$size,
  md5 = unname(
    tools::md5sum(c(
      merged_counts_checkpoint,
      prejoin_checkpoint_path,
      final_checkpoint_path
    ))
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  checkpoint_manifest,
  file.path(
    manifest_root,
    "checkpoint_manifest.csv"
  )
)

artifact_files <- unique(c(
  list.files(
    table_root,
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
  session_path
))

artifact_files <- artifact_files[
  file.exists(artifact_files) &
    !dir.exists(artifact_files)
]

artifact_manifest <- data.frame(
  path = normalize_slashes(
    artifact_files
  ),
  size_bytes =
    file.info(artifact_files)$size,
  md5 = unname(
    tools::md5sum(
      artifact_files
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

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  supersedes_failed_executions =
    "V1_0_HARMONY_ASSAY_USE_FAILURE;V1_1_JOINLAYERS_LENGTH_MISMATCH;V1_2_GGPLOT_DUPLICATE_METADATA_COLUMNS;V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS",
  integration_version =
    integration_version,
  dataset_id = dataset_id,
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
  preQC_cell_n = 50997L,
  predicted_doublet_n =
    total_doublets,
  retained_singlet_n =
    total_singlets,
  sample_n = 9L,
  primary_cluster_n =
    length(
      unique(
        integrated_object$
          seurat_clusters
      )
    ),
  count_layer_n =
    length(count_layers),
  data_layer_n =
    length(data_layers_final),
  final_checkpoint_size_GB =
    round(
      file.info(
        final_checkpoint_path
      )$size / 1024^3,
      3
    ),
  final_checkpoint_md5 =
    unname(
      tools::md5sum(
        final_checkpoint_path
      )
    ),
  detailed_integration_QC_authorized =
    all_gate_checks_pass,
  marker_analysis_authorized =
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
  "preQC_cells=",
  50997L,
  "predicted_doublets=",
  total_doublets,
  "retained_singlets=",
  total_singlets,
  "primary_clusters=",
  length(
    unique(
      integrated_object$
        seurat_clusters
    )
  )
)

small_files <- unique(c(
  list.files(
    table_root,
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
      200 * 1024^2
]

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
      zip_error <<- conditionMessage(e)
      FALSE
    }
  )
}

cat("\n============================================================\n")
cat("STEP04E1J GSE183837 doublet/Harmony integration completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Pre-QC cells: ",
  format(50997L, big.mark = ","),
  "\n",
  sep = ""
)
cat(
  "Predicted doublets: ",
  format(total_doublets, big.mark = ","),
  " (",
  round(
    100 * total_doublets / 50997L,
    2
  ),
  "%)\n",
  sep = ""
)
cat(
  "Retained singlets: ",
  format(total_singlets, big.mark = ","),
  "\n",
  sep = ""
)
cat(
  "Primary resolution-0.4 clusters: ",
  length(
    unique(
      integrated_object$
        seurat_clusters
    )
  ),
  "\n",
  sep = ""
)
cat(
  "Count layers retained: ",
  length(count_layers),
  "\n",
  sep = ""
)
cat(
  "Joined data layers: ",
  length(data_layers_final),
  "\n",
  sep = ""
)
cat(
  "Final checkpoint: ",
  final_checkpoint_path,
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
    "\nThe integration outputs and QC report are complete.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 3: 65_step04E1K_GSE183837_detailed_integration_QC_v1_2_pathsafe
# ==============================================================================

# ======================================================================
# 65_step04E1K_GSE183837_detailed_integration_QC_v1_2_pathsafe.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Perform detailed, donor-aware QC of the completed GSE183837 Harmony
#   integration.
# - Quantify cluster donor/sample dominance, donor entropy, group support,
#   local graph-neighborhood mixing and resolution stability.
# - Review technical quality and residual-doublet indicators per cluster.
# - Specifically review clusters 7, 8, 13 and 15 without deleting cells.
# - Authorize memory-safe descriptive marker extraction only when the
#   detailed integration-QC Gate passes.
#
# UPSTREAM
# - PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS
# - integrated.rds with 47,161 singlets and 16 primary clusters.
#
# SCIENTIFIC GUARDRAILS
# - Sample mixing is assessed as a technical property, not as evidence
#   of disease biology.
# - Control-versus-RIF cell-level inference is prohibited.
# - Donor-dominant or rare clusters are retained for marker review.
# - No cluster is removed or biologically annotated by this script.
# - Marker authorization does not equal annotation lock.
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
dataset_id <- "GSE183837"

script_name <-
  "65_step04E1K_GSE183837_detailed_integration_QC_v1_2_pathsafe.R"

script_version <-
  "STEP04E1K_GSE183837_DETAILED_INTEGRATION_QC_v1.2_PATHSAFE"

integration_version <- "SCI183_v1_0"
QC_version <- "IQ183_v1_2"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1K_",
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

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  path_n <- nchar(
    normalize_slashes(path),
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
      normalize_slashes(path)
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
  }

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_csv_gz <- function(x, path) {
  assert_path_safe(path)

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  connection <- gzfile(
    path,
    open = "wt",
    encoding = "UTF-8"
  )

  on.exit(
    close(connection),
    add = TRUE
  )

  write.csv(
    x,
    connection,
    row.names = FALSE,
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
  }

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

safe_package_version <- function(package_name) {
  if (
    !requireNamespace(
      package_name,
      quietly = TRUE
    )
  ) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(package_name)
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
    as.numeric(as.character(x))
  )
}

choose2 <- function(x) {
  x * (x - 1) / 2
}

adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y)) {
    stop("ARI inputs have different lengths.")
  }

  contingency <- table(x, y)
  n <- sum(contingency)

  if (n < 2L) return(NA_real_)

  sum_comb_cells <- sum(
    choose2(contingency)
  )

  row_comb <- sum(
    choose2(
      rowSums(contingency)
    )
  )

  col_comb <- sum(
    choose2(
      colSums(contingency)
    )
  )

  total_comb <- choose2(n)

  expected <- (
    row_comb * col_comb /
      total_comb
  )

  maximum <- (
    row_comb + col_comb
  ) / 2

  denominator <- maximum - expected

  if (
    !is.finite(denominator) ||
      denominator == 0
  ) {
    return(NA_real_)
  }

  (sum_comb_cells - expected) /
    denominator
}

normalized_entropy <- function(
  counts,
  category_n
) {
  counts <- as.numeric(counts)
  counts <- counts[
    is.finite(counts) &
      counts > 0
  ]

  if (
    !length(counts) ||
      category_n <= 1
  ) {
    return(0)
  }

  probabilities <- counts /
    sum(counts)

  -sum(
    probabilities *
      log(probabilities)
  ) / log(category_n)
}

effective_category_n <- function(counts) {
  counts <- as.numeric(counts)
  counts <- counts[
    is.finite(counts) &
      counts > 0
  ]

  if (!length(counts)) return(0)

  probabilities <- counts /
    sum(counts)

  exp(
    -sum(
      probabilities *
        log(probabilities)
    )
  )
}

median_or_NA <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    median(x)
  }
}

quantile_or_NA <- function(
  x,
  probability
) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    as.numeric(
      stats::quantile(
        x,
        probs = probability,
        names = FALSE,
        type = 7
      )
    )
  }
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^2, 2)
  }
}

get_free_disk_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_LogicalDisk",
    "-Filter \"DeviceID='C:'\").FreeSpace"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(trimws(output[1]))
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(value / 1024^3, 2)
  }
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

make_review_flags <- function(
  cell_n,
  top_donor_fraction,
  donor_entropy,
  Control_donors_ge10,
  RIF_donors_ge10,
  median_percent_mt,
  median_nFeature,
  median_doublet_score
) {
  flags <- character()

  if (
    is.finite(cell_n) &&
      cell_n < 100
  ) {
    flags <- c(
      flags,
      "RARE_CLUSTER_BELOW_100_CELLS"
    )
  }

  if (
    is.finite(top_donor_fraction) &&
      top_donor_fraction >= 0.75
  ) {
    flags <- c(
      flags,
      "SEVERE_SINGLE_DONOR_DOMINANCE_GE_75_PERCENT"
    )
  } else if (
    is.finite(top_donor_fraction) &&
      top_donor_fraction >= 0.50
  ) {
    flags <- c(
      flags,
      "SINGLE_DONOR_ENRICHMENT_GE_50_PERCENT"
    )
  }

  if (
    is.finite(donor_entropy) &&
      donor_entropy < 0.50
  ) {
    flags <- c(
      flags,
      "LOW_DONOR_ENTROPY_BELOW_0_50"
    )
  } else if (
    is.finite(donor_entropy) &&
      donor_entropy < 0.70
  ) {
    flags <- c(
      flags,
      "MODERATE_DONOR_ENTROPY_BELOW_0_70"
    )
  }

  if (
    is.finite(Control_donors_ge10) &&
      Control_donors_ge10 < 2
  ) {
    flags <- c(
      flags,
      "LOW_CONTROL_DONOR_SUPPORT_GE10_CELLS"
    )
  }

  if (
    is.finite(RIF_donors_ge10) &&
      RIF_donors_ge10 < 3
  ) {
    flags <- c(
      flags,
      "LOW_RIF_DONOR_SUPPORT_GE10_CELLS"
    )
  }

  if (
    is.finite(median_percent_mt) &&
      is.finite(median_nFeature) &&
      median_percent_mt > 20 &&
      median_nFeature < 1000
  ) {
    flags <- c(
      flags,
      "POSSIBLE_LOW_QUALITY_CLUSTER"
    )
  }

  if (
    is.finite(median_doublet_score) &&
      median_doublet_score > 0.20
  ) {
    flags <- c(
      flags,
      "POSSIBLE_RESIDUAL_DOUBLET_CLUSTER"
    )
  }

  collapse_unique(flags)
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
# 2. Package preflight
# ----------------------------------------------------------------------

required_packages <- c(
  "SeuratObject",
  "Matrix",
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

parameter_root <- file.path(
  step04_root,
  "03_parameters"
)

environment_root <- file.path(
  step04_root,
  "04_environment"
)

integration_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

input_checkpoint <- file.path(
  integration_root,
  "checkpoints",
  dataset_id,
  "integrated.rds"
)

output_root <- file.path(
  integration_root,
  "iqc",
  "G183",
  QC_version
)

table_root <- file.path(
  output_root,
  "t"
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
  "m"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1K_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_dQC"
  )
)

figure_root <- file.path(
  run_qc_root,
  "figures"
)

for (path_value in c(
  execution_root,
  parameter_root,
  environment_root,
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
  artifact = c(
    "input_checkpoint",
    "output_root",
    "table_root",
    "source_data_root",
    "figure_output_root",
    "manifest_root",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    input_checkpoint,
    output_root,
    table_root,
    source_data_root,
    figure_output_root,
    manifest_root,
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

if (!all(planned_paths$path_safe)) {
  stop(
    "One or more planned paths are too long:\n",
    paste(
      planned_paths$path[
        !planned_paths$path_safe
      ],
      collapse = "\n"
    )
  )
}

write_csv_utf8(
  planned_paths,
  file.path(run_qc_root, "00_path_audit.csv")
)

planned_output_files <- data.frame(
  artifact = c(
    "cluster_merge_audit",
    "checkpoint_audit",
    "donor_identity",
    "cluster_sample",
    "cluster_QC",
    "review_clusters",
    "mixing_sample",
    "mixing_cluster",
    "graph_audit",
    "resolution_ARI",
    "resolution_02_to_04",
    "resolution_04_to_06",
    "cluster_stability",
    "target_review",
    "source_cluster_sample",
    "source_cell_mixing",
    "persistent_manifest",
    "detailed_gate",
    "detailed_card",
    "detailed_interpretation"
  ),
  path = normalize_slashes(c(
    file.path(table_root, "00_cluster_merge_count_audit.csv"),
    file.path(table_root, "01_checkpoint_audit.csv"),
    file.path(table_root, "02_donor_group_identity.csv"),
    file.path(table_root, "03_cluster_sample.csv"),
    file.path(table_root, "04_cluster_QC.csv"),
    file.path(table_root, "05_review_clusters.csv"),
    file.path(table_root, "06_mixing_sample.csv"),
    file.path(table_root, "07_mixing_cluster.csv"),
    file.path(table_root, "08_graph_audit.csv"),
    file.path(table_root, "09_resolution_ARI.csv"),
    file.path(table_root, "10_res02_to_04.csv"),
    file.path(table_root, "11_res04_to_06.csv"),
    file.path(table_root, "12_cluster_stability.csv"),
    file.path(table_root, "13_review_7_8_13_15.csv"),
    file.path(source_data_root, "SD1_cluster_sample.csv"),
    file.path(source_data_root, "SD2_cell_mixing.csv.gz"),
    file.path(manifest_root, "persistent_manifest.csv"),
    file.path(run_qc_root, "G183_DETAILED_QC_GATE.csv"),
    file.path(run_qc_root, "G183_DETAILED_QC_CARD.md"),
    file.path(run_qc_root, "G183_DETAILED_QC_INTERPRETATION.md")
  )),
  stringsAsFactors = FALSE
)

planned_output_files$character_n <- nchar(
  planned_output_files$path,
  type = "chars"
)

planned_output_files$path_safe <- (
  planned_output_files$character_n <= 235L
)

write_csv_utf8(
  planned_output_files,
  file.path(
    run_qc_root,
    "00_output_file_path_audit.csv"
  )
)

if (!all(planned_output_files$path_safe)) {
  stop(
    "One or more planned output files exceed the safe path limit:\n",
    paste(
      planned_output_files$path[
        !planned_output_files$path_safe
      ],
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream Gate
# ----------------------------------------------------------------------

integration_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

integration_qc_dirs <- integration_qc_dirs[
  grepl(
    "^STEP04E1J_[0-9]{8}_[0-9]{6}_integration$",
    basename(integration_qc_dirs),
    perl = TRUE
  )
]

if (!length(integration_qc_dirs)) {
  stop(
    "No STEP04E1J integration-QC directory was found."
  )
}

candidate_gate_rows <- list()

for (directory_path in integration_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "GSE183837_INTEGRATION_GATE.csv"
  )

  gate_value <- safe_read_csv(gate_path)

  if (
    !is.null(gate_value) &&
      nrow(gate_value) == 1L
  ) {
    gate_value$source_directory <-
      normalize_slashes(directory_path)

    gate_value$source_mtime <-
      as.character(
        file.info(directory_path)$mtime
      )

    candidate_gate_rows[[
      length(candidate_gate_rows) + 1L
    ]] <- gate_value
  }
}

if (!length(candidate_gate_rows)) {
  stop(
    "No readable GSE183837 integration Gate was found."
  )
}

candidate_gates <- as.data.frame(
  data.table::rbindlist(
    candidate_gate_rows,
    use.names = TRUE,
    fill = TRUE
  )
)

candidate_gates$source_mtime_POSIX <- as.POSIXct(
  candidate_gates$source_mtime
)

candidate_gates <- candidate_gates[
  order(
    candidate_gates$source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

authorized_gate_index <- which(
  candidate_gates$gate_status %in% c(
    "PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED",
    "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS"
  ) &
    vapply(
      candidate_gates$
        detailed_integration_QC_authorized,
      function(x) {
        isTRUE(x) ||
          toupper(as.character(x)) %in%
            c("TRUE", "T", "1", "YES")
      },
      logical(1)
    )
)

if (!length(authorized_gate_index)) {
  stop(
    "No STEP04E1J Gate authorizes detailed integration QC."
  )
}

upstream_gate <- as.data.frame(
  candidate_gates[
    authorized_gate_index[1],
    ,
    drop = FALSE
  ]
)

write_csv_utf8(
  upstream_gate,
  file.path(
    run_qc_root,
    "00_upstream_integration_gate.csv"
  )
)


# ----------------------------------------------------------------------
# 5. Load and validate checkpoint
# ----------------------------------------------------------------------

if (!file.exists(input_checkpoint)) {
  stop(
    "Integrated checkpoint is missing:\n",
    input_checkpoint
  )
}

free_memory_before_GB <- get_free_memory_gb()
free_disk_before_GB <- get_free_disk_gb()

if (
  is.finite(free_memory_before_GB) &&
    free_memory_before_GB < 10
) {
  stop(
    "Free physical memory is below 10 GB: ",
    free_memory_before_GB
  )
}

log_event(
  "Starting",
  run_id,
  "free_memory_GB=",
  free_memory_before_GB,
  "free_disk_GB=",
  free_disk_before_GB
)

integrated_object <- readRDS(
  input_checkpoint
)

if (!inherits(integrated_object, "Seurat")) {
  stop(
    "The integrated checkpoint is not a Seurat object."
  )
}

metadata <- integrated_object[[]]

required_metadata_columns <- c(
  "sample_accession",
  "donor_id",
  "analysis_group",
  "seurat_clusters",
  "harmony_snn_res.0.2",
  "harmony_snn_res.0.4",
  "harmony_snn_res.0.6",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "scDblFinder_score"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns)) {
  stop(
    "Integrated metadata is missing required columns: ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

if (
  ncol(integrated_object) != 47161L ||
    nrow(metadata) != 47161L ||
    length(
      unique(
        as.character(
          metadata$sample_accession
        )
      )
    ) != 9L ||
    length(
      unique(
        as.character(
          metadata$donor_id
        )
      )
    ) != 9L ||
    length(
      unique(
        as.character(
          metadata$seurat_clusters
        )
      )
    ) != 16L
) {
  stop(
    "Integrated checkpoint identity does not match the locked 47,161-cell, 9-donor, 16-cluster object."
  )
}

metadata$cell_barcode <- rownames(metadata)

metadata$sample_accession <- as.character(
  metadata$sample_accession
)

metadata$donor_id <- as.character(
  metadata$donor_id
)

metadata$analysis_group <- factor(
  as.character(
    metadata$analysis_group
  ),
  levels = c(
    "CONTROL",
    "RIF"
  )
)

metadata$seurat_clusters <- factor(
  as.character(
    metadata$seurat_clusters
  ),
  levels = as.character(0:15)
)

if (
  anyNA(metadata$analysis_group) ||
    anyNA(metadata$seurat_clusters)
) {
  stop(
    "Group or primary cluster metadata contains unexpected values."
  )
}

checkpoint_audit <- data.frame(
  metric = c(
    "feature_n",
    "cell_n",
    "sample_n",
    "donor_n",
    "group_n",
    "primary_cluster_n",
    "checkpoint_size_GB",
    "checkpoint_md5"
  ),
  value = c(
    nrow(integrated_object),
    ncol(integrated_object),
    length(unique(metadata$sample_accession)),
    length(unique(metadata$donor_id)),
    length(unique(metadata$analysis_group)),
    length(unique(metadata$seurat_clusters)),
    round(
      file.info(input_checkpoint)$size /
        1024^3,
      3
    ),
    unname(
      tools::md5sum(
        input_checkpoint
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  checkpoint_audit,
  file.path(
    table_root,
    "01_checkpoint_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 6. Donor/group identity and cluster composition
# ----------------------------------------------------------------------

metadata_DT <- data.table::as.data.table(
  metadata
)

donor_group_identity <- unique(
  metadata_DT[
    ,
    .(
      donor_id,
      sample_accession,
      analysis_group
    )
  ]
)

data.table::setorder(
  donor_group_identity,
  donor_id
)

if (
  nrow(donor_group_identity) != 9L ||
    anyDuplicated(
      donor_group_identity$donor_id
    ) > 0L
) {
  stop(
    "Donor/group identity is not one row per nine donors."
  )
}

write_csv_utf8(
  as.data.frame(donor_group_identity),
  file.path(
    table_root,
    "02_donor_group_identity.csv"
  )
)

cluster_sample_observed <- metadata_DT[
  ,
  .(
    cell_n = .N,
    median_nCount_RNA =
      median_or_NA(nCount_RNA),
    median_nFeature_RNA =
      median_or_NA(nFeature_RNA),
    median_percent_mt =
      median_or_NA(percent.mt),
    median_scDblFinder_score =
      median_or_NA(
        scDblFinder_score
      )
  ),
  by = .(
    seurat_clusters,
    sample_accession,
    donor_id,
    analysis_group
  )
]

cluster_levels <- as.character(0:15)
sample_levels <- donor_group_identity$
  sample_accession

cluster_sample_grid <- as.data.frame(
  data.table::CJ(
    seurat_clusters = cluster_levels,
    sample_accession = sample_levels,
    unique = TRUE
  )
)

cluster_sample_grid <- merge(
  cluster_sample_grid,
  as.data.frame(
    donor_group_identity
  ),
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

cluster_sample_full <- as.data.frame(
  merge(
    cluster_sample_grid,
    as.data.frame(
      cluster_sample_observed
    ),
    by = c(
      "seurat_clusters",
      "sample_accession",
      "donor_id",
      "analysis_group"
    ),
    all.x = TRUE,
    sort = FALSE
  )
)

cluster_sample_full$cell_n[
  is.na(cluster_sample_full$cell_n)
] <- 0L

cluster_sample_full$seurat_clusters <-
  factor(
    cluster_sample_full$
      seurat_clusters,
    levels = cluster_levels
  )

cluster_sample_full <- cluster_sample_full[
  order(
    as.integer(
      as.character(
        cluster_sample_full$
          seurat_clusters
      )
    ),
    match(
      cluster_sample_full$
        sample_accession,
      sample_levels
    )
  ),
  ,
  drop = FALSE
]

cluster_totals <- aggregate(
  cell_n ~ seurat_clusters,
  data = cluster_sample_full,
  FUN = sum
)

colnames(cluster_totals)[2] <-
  "cluster_cell_n"

cluster_sample_full <- merge(
  cluster_sample_full,
  cluster_totals,
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE
)

cluster_sample_full$
  within_cluster_fraction <- (
    cluster_sample_full$cell_n /
      cluster_sample_full$
        cluster_cell_n
  )

sample_totals <- aggregate(
  cell_n ~ sample_accession,
  data = cluster_sample_full,
  FUN = sum
)

colnames(sample_totals)[2] <-
  "sample_cell_n"

cluster_sample_full <- merge(
  cluster_sample_full,
  sample_totals,
  by = "sample_accession",
  all.x = TRUE,
  sort = FALSE
)

cluster_sample_full$
  within_sample_fraction <- (
    cluster_sample_full$cell_n /
      cluster_sample_full$
        sample_cell_n
  )

write_csv_utf8(
  cluster_sample_full,
  file.path(
    table_root,
    "03_cluster_sample.csv"
  )
)

write_csv_utf8(
  cluster_sample_full,
  file.path(
    source_data_root,
    "SD1_cluster_sample.csv"
  )
)

cluster_QC <- metadata_DT[
  ,
  .(
    technical_cell_n = .N,
    median_nCount_RNA =
      median_or_NA(nCount_RNA),
    Q1_nCount_RNA =
      quantile_or_NA(
        nCount_RNA,
        0.25
      ),
    Q3_nCount_RNA =
      quantile_or_NA(
        nCount_RNA,
        0.75
      ),
    median_nFeature_RNA =
      median_or_NA(nFeature_RNA),
    Q1_nFeature_RNA =
      quantile_or_NA(
        nFeature_RNA,
        0.25
      ),
    Q3_nFeature_RNA =
      quantile_or_NA(
        nFeature_RNA,
        0.75
      ),
    median_percent_mt =
      median_or_NA(percent.mt),
    Q1_percent_mt =
      quantile_or_NA(
        percent.mt,
        0.25
      ),
    Q3_percent_mt =
      quantile_or_NA(
        percent.mt,
        0.75
      ),
    median_scDblFinder_score =
      median_or_NA(
        scDblFinder_score
      ),
    Q90_scDblFinder_score =
      quantile_or_NA(
        scDblFinder_score,
        0.90
      ),
    Q95_scDblFinder_score =
      quantile_or_NA(
        scDblFinder_score,
        0.95
      )
  ),
  by = seurat_clusters
]

cluster_composition_rows <- list()

for (cluster_value in cluster_levels) {
  cluster_data <- cluster_sample_full[
    as.character(
      cluster_sample_full$
        seurat_clusters
    ) == cluster_value,
    ,
    drop = FALSE
  ]

  cluster_size <- sum(
    cluster_data$cell_n
  )

  top_index <- which.max(
    cluster_data$cell_n
  )

  group_cell_counts <- tapply(
    cluster_data$cell_n,
    cluster_data$analysis_group,
    sum
  )

  Control_cell_n <- if (
    "CONTROL" %in%
      names(group_cell_counts)
  ) {
    as.integer(
      group_cell_counts[
        "CONTROL"
      ]
    )
  } else {
    0L
  }

  RIF_cell_n <- if (
    "RIF" %in%
      names(group_cell_counts)
  ) {
    as.integer(
      group_cell_counts[
        "RIF"
      ]
    )
  } else {
    0L
  }

  cluster_composition_rows[[
    length(cluster_composition_rows) + 1L
  ]] <- data.frame(
    seurat_clusters =
      cluster_value,
    cell_n =
      cluster_size,
    cell_fraction =
      cluster_size /
        nrow(metadata),
    represented_donor_n =
      sum(
        cluster_data$cell_n > 0
      ),
    donors_ge5_cells =
      sum(
        cluster_data$cell_n >= 5
      ),
    donors_ge10_cells =
      sum(
        cluster_data$cell_n >= 10
      ),
    Control_donors_ge1_cell =
      sum(
        cluster_data$
          analysis_group ==
          "CONTROL" &
          cluster_data$cell_n > 0
      ),
    RIF_donors_ge1_cell =
      sum(
        cluster_data$
          analysis_group ==
          "RIF" &
          cluster_data$cell_n > 0
      ),
    Control_donors_ge10_cells =
      sum(
        cluster_data$
          analysis_group ==
          "CONTROL" &
          cluster_data$cell_n >= 10
      ),
    RIF_donors_ge10_cells =
      sum(
        cluster_data$
          analysis_group ==
          "RIF" &
          cluster_data$cell_n >= 10
      ),
    top_donor =
      cluster_data$
        sample_accession[
          top_index
        ],
    top_donor_group =
      as.character(
        cluster_data$
          analysis_group[
            top_index
          ]
      ),
    top_donor_cell_n =
      cluster_data$
        cell_n[
          top_index
        ],
    top_donor_fraction =
      cluster_data$
        cell_n[
          top_index
        ] /
        cluster_size,
    donor_entropy_normalized =
      normalized_entropy(
        cluster_data$cell_n,
        category_n = 9
      ),
    effective_donor_n =
      effective_category_n(
        cluster_data$cell_n
      ),
    Control_cell_n =
      Control_cell_n,
    RIF_cell_n =
      RIF_cell_n,
    Control_cell_fraction =
      Control_cell_n /
        cluster_size,
    RIF_cell_fraction =
      RIF_cell_n /
        cluster_size,
    stringsAsFactors = FALSE
  )
}

cluster_composition <- do.call(
  rbind,
  cluster_composition_rows
)

cluster_review <- merge(
  cluster_composition,
  as.data.frame(cluster_QC),
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE
)

required_cluster_review_columns <- c(
  "seurat_clusters",
  "cell_n",
  "technical_cell_n",
  "top_donor_fraction",
  "donor_entropy_normalized",
  "Control_donors_ge10_cells",
  "RIF_donors_ge10_cells",
  "median_percent_mt",
  "median_nFeature_RNA",
  "median_scDblFinder_score"
)

missing_cluster_review_columns <- setdiff(
  required_cluster_review_columns,
  colnames(cluster_review)
)

if (length(missing_cluster_review_columns)) {
  stop(
    "cluster_review merge is missing required columns: ",
    paste(
      missing_cluster_review_columns,
      collapse = ", "
    )
  )
}

if (
  nrow(cluster_review) != 16L ||
    anyNA(cluster_review$seurat_clusters) ||
    anyNA(cluster_review$cell_n) ||
    anyNA(cluster_review$technical_cell_n)
) {
  stop(
    "cluster_review merge did not produce 16 complete cluster rows."
  )
}

cluster_count_consistency <- data.frame(
  seurat_clusters =
    as.character(
      cluster_review$seurat_clusters
    ),
  composition_cell_n =
    as.integer(
      cluster_review$cell_n
    ),
  technical_cell_n =
    as.integer(
      cluster_review$technical_cell_n
    ),
  counts_identical =
    as.integer(
      cluster_review$cell_n
    ) ==
      as.integer(
        cluster_review$technical_cell_n
      ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  cluster_count_consistency,
  file.path(
    table_root,
    "00_cluster_merge_count_audit.csv"
  )
)

if (
  !all(
    cluster_count_consistency$
      counts_identical
  )
) {
  stop(
    "Cluster cell counts differ between composition and technical-QC tables."
  )
}

cluster_review$seurat_clusters <-
  as.character(
    cluster_review$seurat_clusters
  )

cluster_review <- cluster_review[
  match(
    cluster_levels,
    cluster_review$seurat_clusters
  ),
  ,
  drop = FALSE
]

review_flag_vector <- mapply(
  make_review_flags,
  cell_n =
    cluster_review$cell_n,
  top_donor_fraction =
    cluster_review$
      top_donor_fraction,
  donor_entropy =
    cluster_review$
      donor_entropy_normalized,
  Control_donors_ge10 =
    cluster_review$
      Control_donors_ge10_cells,
  RIF_donors_ge10 =
    cluster_review$
      RIF_donors_ge10_cells,
  median_percent_mt =
    cluster_review$
      median_percent_mt,
  median_nFeature =
    cluster_review$
      median_nFeature_RNA,
  median_doublet_score =
    cluster_review$
      median_scDblFinder_score,
  SIMPLIFY = TRUE,
  USE.NAMES = FALSE
)

if (
  length(review_flag_vector) !=
    nrow(cluster_review)
) {
  stop(
    "Review-flag vector length mismatch: expected ",
    nrow(cluster_review),
    ", obtained ",
    length(review_flag_vector),
    "."
  )
}

cluster_review$review_flags <-
  as.character(
    review_flag_vector
  )

cluster_review$review_required <-
  nzchar(
    cluster_review$
      review_flags
  )

cluster_review$critical_hold_flag <- (
  (
    cluster_review$
      top_donor_fraction >= 0.95 &
      cluster_review$cell_n >= 200
  ) |
    cluster_review$cell_n < 30 |
    (
      cluster_review$
        median_percent_mt > 22 &
        cluster_review$
          median_nFeature_RNA < 800
    ) |
    cluster_review$
      median_scDblFinder_score > 0.30
)

write_csv_utf8(
  cluster_review,
  file.path(
    table_root,
    "04_cluster_QC.csv"
  )
)

write_csv_utf8(
  cluster_review[
    cluster_review$review_required,
    ,
    drop = FALSE
  ],
  file.path(
    table_root,
    "05_review_clusters.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Harmony graph local-neighborhood mixing
# ----------------------------------------------------------------------

graph_names <- names(
  integrated_object@graphs
)

preferred_graph_names <- c(
  "harmony_nn",
  "harmony_snn"
)

graph_name <- preferred_graph_names[
  preferred_graph_names %in%
    graph_names
][1]

if (
  length(graph_name) != 1L ||
    is.na(graph_name)
) {
  stop(
    "Neither harmony_nn nor harmony_snn graph is present."
  )
}

mixing_graph <- integrated_object@graphs[[graph_name]]

if (
  nrow(mixing_graph) != nrow(metadata) ||
    ncol(mixing_graph) != nrow(metadata)
) {
  stop(
    "Harmony graph dimensions do not match metadata."
  )
}

graph_cells <- rownames(
  mixing_graph
)

if (is.null(graph_cells)) {
  graph_cells <- rownames(metadata)
}

if (
  !identical(
    graph_cells,
    rownames(metadata)
  )
) {
  if (
    setequal(
      graph_cells,
      rownames(metadata)
    )
  ) {
    metadata_graph <- metadata[
      graph_cells,
      ,
      drop = FALSE
    ]
  } else {
    stop(
      "Harmony graph cell identities do not match metadata."
    )
  }
} else {
  metadata_graph <- metadata
}

graph_edges <- data.table::as.data.table(
  Matrix::summary(
    mixing_graph
  )
)

required_edge_columns <- c(
  "i",
  "j",
  "x"
)

if (
  !all(
    required_edge_columns %in%
      colnames(graph_edges)
  )
) {
  stop(
    "Harmony graph summary lacks i, j and x columns."
  )
}

graph_edges <- graph_edges[
  i != j &
    is.finite(x) &
    x > 0
]

sample_code <- as.integer(
  factor(
    metadata_graph$
      sample_accession,
    levels =
      donor_group_identity$
        sample_accession
  )
)

group_code <- as.integer(
  factor(
    metadata_graph$
      analysis_group,
    levels = c(
      "CONTROL",
      "RIF"
    )
  )
)

cluster_code <- as.character(
  metadata_graph$
    seurat_clusters
)

graph_edges[
  ,
  `:=`(
    same_sample =
      sample_code[i] ==
        sample_code[j],
    same_group =
      group_code[i] ==
        group_code[j],
    neighbor_sample_code =
      sample_code[j]
  )
]

cell_mixing_basic <- graph_edges[
  ,
  .(
    neighbor_n = .N,
    weighted_degree =
      sum(x),
    same_sample_fraction =
      mean(same_sample),
    weighted_same_sample_fraction =
      sum(
        x * same_sample
      ) /
      sum(x),
    same_group_fraction =
      mean(same_group),
    weighted_same_group_fraction =
      sum(
        x * same_group
      ) /
      sum(x)
  ),
  by = i
]

neighbor_sample_counts <- graph_edges[
  ,
  .N,
  by = .(
    i,
    neighbor_sample_code
  )
]

neighbor_diversity <- neighbor_sample_counts[
  ,
  {
    probabilities <- N /
      sum(N)

    list(
      neighbor_donor_n =
        .N,
      neighbor_donor_entropy =
        -sum(
          probabilities *
            log(probabilities)
        ) / log(9),
      effective_neighbor_donor_n =
        exp(
          -sum(
            probabilities *
              log(probabilities)
          )
        )
    )
  },
  by = i
]

cell_mixing <- merge(
  cell_mixing_basic,
  neighbor_diversity,
  by = "i",
  all = TRUE,
  sort = FALSE
)

cell_identity <- data.frame(
  i = seq_len(
    nrow(metadata_graph)
  ),
  cell_barcode =
    rownames(metadata_graph),
  sample_accession =
    metadata_graph$
      sample_accession,
  donor_id =
    metadata_graph$
      donor_id,
  analysis_group =
    as.character(
      metadata_graph$
        analysis_group
    ),
  seurat_clusters =
    as.character(
      metadata_graph$
        seurat_clusters
    ),
  stringsAsFactors = FALSE
)

cell_mixing <- merge(
  cell_identity,
  cell_mixing,
  by = "i",
  all.x = TRUE,
  sort = FALSE
)

cell_mixing <- cell_mixing[
  order(cell_mixing$i),
  ,
  drop = FALSE
]

expected_sample_fraction <- cluster_sample_full[
  ,
  c(
    "seurat_clusters",
    "sample_accession",
    "cell_n",
    "cluster_cell_n"
  )
]

colnames(expected_sample_fraction)[
  colnames(expected_sample_fraction) ==
    "cell_n"
] <- "cluster_sample_cell_n"

cell_mixing <- merge(
  cell_mixing,
  expected_sample_fraction,
  by = c(
    "seurat_clusters",
    "sample_accession"
  ),
  all.x = TRUE,
  sort = FALSE
)

cell_mixing$
  expected_same_sample_fraction <- ifelse(
    cell_mixing$
      cluster_cell_n > 1,
    pmax(
      0,
      (
        cell_mixing$
          cluster_sample_cell_n - 1
      ) /
        (
          cell_mixing$
            cluster_cell_n - 1
        )
    ),
    0
  )

cell_mixing$
  same_sample_excess <- (
    cell_mixing$
      same_sample_fraction -
      cell_mixing$
        expected_same_sample_fraction
  )

cell_mixing$
  weighted_same_sample_excess <- (
    cell_mixing$
      weighted_same_sample_fraction -
      cell_mixing$
        expected_same_sample_fraction
  )

write_csv_gz(
  cell_mixing,
  file.path(
    source_data_root,
    "SD2_cell_mixing.csv.gz"
  )
)

mixing_sample_summary <- data.table::as.data.table(
  cell_mixing
)[
  ,
  .(
    cell_n = .N,
    median_neighbor_n =
      median_or_NA(neighbor_n),
    median_same_sample_fraction =
      median_or_NA(
        same_sample_fraction
      ),
    Q90_same_sample_fraction =
      quantile_or_NA(
        same_sample_fraction,
        0.90
      ),
    median_expected_same_sample_fraction =
      median_or_NA(
        expected_same_sample_fraction
      ),
    median_same_sample_excess =
      median_or_NA(
        same_sample_excess
      ),
    Q90_same_sample_excess =
      quantile_or_NA(
        same_sample_excess,
        0.90
      ),
    median_neighbor_donor_n =
      median_or_NA(
        neighbor_donor_n
      ),
    median_neighbor_donor_entropy =
      median_or_NA(
        neighbor_donor_entropy
      )
  ),
  by = .(
    sample_accession,
    donor_id,
    analysis_group
  )
]

mixing_cluster_summary <- data.table::as.data.table(
  cell_mixing
)[
  ,
  .(
    cell_n = .N,
    median_neighbor_n =
      median_or_NA(neighbor_n),
    median_same_sample_fraction =
      median_or_NA(
        same_sample_fraction
      ),
    Q90_same_sample_fraction =
      quantile_or_NA(
        same_sample_fraction,
        0.90
      ),
    median_expected_same_sample_fraction =
      median_or_NA(
        expected_same_sample_fraction
      ),
    median_same_sample_excess =
      median_or_NA(
        same_sample_excess
      ),
    Q90_same_sample_excess =
      quantile_or_NA(
        same_sample_excess,
        0.90
      ),
    median_neighbor_donor_n =
      median_or_NA(
        neighbor_donor_n
      ),
    median_neighbor_donor_entropy =
      median_or_NA(
        neighbor_donor_entropy
      )
  ),
  by = seurat_clusters
]

mixing_cluster_summary$
  local_mixing_review_flag <- (
    mixing_cluster_summary$
      median_same_sample_excess > 0.25 |
      mixing_cluster_summary$
        Q90_same_sample_excess > 0.50 |
      mixing_cluster_summary$
        median_neighbor_donor_n < 3
  )

write_csv_utf8(
  as.data.frame(
    mixing_sample_summary
  ),
  file.path(
    table_root,
    "06_mixing_sample.csv"
  )
)

write_csv_utf8(
  as.data.frame(
    mixing_cluster_summary
  ),
  file.path(
    table_root,
    "07_mixing_cluster.csv"
  )
)

graph_audit <- data.frame(
  metric = c(
    "graph_name",
    "graph_cell_n",
    "nonzero_edge_n_excluding_self",
    "median_neighbor_n",
    "minimum_neighbor_n",
    "maximum_neighbor_n",
    "overall_median_same_sample_fraction",
    "overall_median_expected_same_sample_fraction",
    "overall_median_same_sample_excess",
    "overall_median_neighbor_donor_n",
    "overall_median_neighbor_donor_entropy"
  ),
  value = c(
    graph_name,
    nrow(mixing_graph),
    nrow(graph_edges),
    median_or_NA(
      cell_mixing$neighbor_n
    ),
    min(
      cell_mixing$neighbor_n,
      na.rm = TRUE
    ),
    max(
      cell_mixing$neighbor_n,
      na.rm = TRUE
    ),
    median_or_NA(
      cell_mixing$
        same_sample_fraction
    ),
    median_or_NA(
      cell_mixing$
        expected_same_sample_fraction
    ),
    median_or_NA(
      cell_mixing$
        same_sample_excess
    ),
    median_or_NA(
      cell_mixing$
        neighbor_donor_n
    ),
    median_or_NA(
      cell_mixing$
        neighbor_donor_entropy
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  graph_audit,
  file.path(
    table_root,
    "08_graph_audit.csv"
  )
)

rm(
  graph_edges,
  neighbor_sample_counts,
  cell_mixing_basic,
  neighbor_diversity
)
gc()


# ----------------------------------------------------------------------
# 8. Resolution stability
# ----------------------------------------------------------------------

resolution_02 <- as.character(
  metadata$
    harmony_snn_res.0.2
)

resolution_04 <- as.character(
  metadata$
    harmony_snn_res.0.4
)

resolution_06 <- as.character(
  metadata$
    harmony_snn_res.0.6
)

resolution_ARI <- data.frame(
  comparison = c(
    "0.2_vs_0.4",
    "0.4_vs_0.6",
    "0.2_vs_0.6"
  ),
  adjusted_rand_index = c(
    adjusted_rand_index(
      resolution_02,
      resolution_04
    ),
    adjusted_rand_index(
      resolution_04,
      resolution_06
    ),
    adjusted_rand_index(
      resolution_02,
      resolution_06
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  resolution_ARI,
  file.path(
    table_root,
    "09_resolution_ARI.csv"
  )
)

make_transition_table <- function(
  source_values,
  target_values,
  source_name,
  target_name
) {
  contingency <- as.data.frame(
    table(
      source_values,
      target_values
    ),
    stringsAsFactors = FALSE
  )

  colnames(contingency) <- c(
    source_name,
    target_name,
    "cell_n"
  )

  source_totals <- aggregate(
    cell_n ~ .,
    data = contingency[
      ,
      c(
        source_name,
        "cell_n"
      )
    ],
    FUN = sum
  )

  colnames(source_totals)[2] <-
    "source_cluster_n"

  target_totals <- aggregate(
    cell_n ~ .,
    data = contingency[
      ,
      c(
        target_name,
        "cell_n"
      )
    ],
    FUN = sum
  )

  colnames(target_totals)[2] <-
    "target_cluster_n"

  output <- merge(
    contingency,
    source_totals,
    by = source_name,
    all.x = TRUE,
    sort = FALSE
  )

  output <- merge(
    output,
    target_totals,
    by = target_name,
    all.x = TRUE,
    sort = FALSE
  )

  output$within_source_fraction <- (
    output$cell_n /
      output$source_cluster_n
  )

  output$within_target_fraction <- (
    output$cell_n /
      output$target_cluster_n
  )

  output
}

transition_02_04 <- make_transition_table(
  resolution_02,
  resolution_04,
  "cluster_res_0_2",
  "cluster_res_0_4"
)

transition_04_06 <- make_transition_table(
  resolution_04,
  resolution_06,
  "cluster_res_0_4",
  "cluster_res_0_6"
)

write_csv_utf8(
  transition_02_04,
  file.path(
    table_root,
    "10_res02_to_04.csv"
  )
)

write_csv_utf8(
  transition_04_06,
  file.path(
    table_root,
    "11_res04_to_06.csv"
  )
)

primary_resolution_stability_rows <- list()

for (cluster_value in cluster_levels) {
  primary_index <- which(
    resolution_04 ==
      cluster_value
  )

  parent_table <- sort(
    table(
      resolution_02[
        primary_index
      ]
    ),
    decreasing = TRUE
  )

  child_table <- sort(
    table(
      resolution_06[
        primary_index
      ]
    ),
    decreasing = TRUE
  )

  primary_resolution_stability_rows[[
    length(primary_resolution_stability_rows) + 1L
  ]] <- data.frame(
    seurat_clusters =
      cluster_value,
    cell_n =
      length(primary_index),
    dominant_res_0_2_parent =
      names(parent_table)[1],
    dominant_parent_fraction =
      as.numeric(
        parent_table[1]
      ) /
        length(primary_index),
    res_0_2_parent_n =
      length(parent_table),
    dominant_res_0_6_child =
      names(child_table)[1],
    dominant_child_fraction =
      as.numeric(
        child_table[1]
      ) /
        length(primary_index),
    res_0_6_child_n =
      length(child_table),
    stringsAsFactors = FALSE
  )
}

primary_resolution_stability <- do.call(
  rbind,
  primary_resolution_stability_rows
)

write_csv_utf8(
  primary_resolution_stability,
  file.path(
    table_root,
    "12_cluster_stability.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Targeted review of clusters 7, 8, 13 and 15
# ----------------------------------------------------------------------

target_review_clusters <- c(
  "7",
  "8",
  "13",
  "15"
)

target_cluster_review <- cluster_review[
  cluster_review$
    seurat_clusters %in%
    target_review_clusters,
  ,
  drop = FALSE
]

target_cluster_review <- merge(
  target_cluster_review,
  as.data.frame(
    mixing_cluster_summary
  ),
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "",
    "_mixing"
  )
)

target_cluster_review <- merge(
  target_cluster_review,
  primary_resolution_stability,
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "",
    "_resolution"
  )
)

target_cluster_review <- target_cluster_review[
  match(
    target_review_clusters,
    target_cluster_review$
      seurat_clusters
  ),
  ,
  drop = FALSE
]

target_cluster_review$
  provisional_QC_interpretation <- c(
    "Donor-enriched cluster; retain for markers and determine whether it represents a donor-associated biological state.",
    "Donor-enriched cluster with low transcript complexity but low mitochondrial and doublet scores; do not label as low quality before marker review.",
    "Strongly dominated by one RIF donor; retain for marker extraction but prohibit disease interpretation unless the identity is reproduced across donors.",
    "Rare but represented across all nine donors and not overtly low quality or doublet-like; retain as a candidate rare population for marker review."
  )

write_csv_utf8(
  target_cluster_review,
  file.path(
    table_root,
    "13_review_7_8_13_15.csv"
  )
)


# ----------------------------------------------------------------------
# 10. QC figures
# ----------------------------------------------------------------------

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

heatmap_data <- cluster_sample_full

heatmap_data$seurat_clusters <-
  factor(
    as.character(
      heatmap_data$
        seurat_clusters
    ),
    levels = rev(cluster_levels)
  )

heatmap_data$sample_accession <-
  factor(
    heatmap_data$
      sample_accession,
    levels = sample_levels
  )

plot_donor_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = sample_accession,
    y = seurat_clusters,
    fill = within_cluster_fraction
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "GSE183837 donor contribution within each primary cluster",
    subtitle = "Each row sums to 100%; high single-donor fractions require marker-level review",
    x = "Donor/sample",
    y = "Primary cluster",
    fill = "Within-cluster fraction"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  )

save_plot_pair(
  plot_donor_heatmap,
  "QC1_cluster_donor_heatmap",
  width = 10,
  height = 7
)

dominance_plot_data <- cluster_review

dominance_plot_data$
  seurat_clusters <- factor(
    dominance_plot_data$
      seurat_clusters,
    levels = cluster_levels
  )

plot_dominance <- ggplot2::ggplot(
  dominance_plot_data,
  ggplot2::aes(
    x = top_donor_fraction,
    y = donor_entropy_normalized,
    size = cell_n,
    label = seurat_clusters,
    shape = review_required
  )
) +
  ggplot2::geom_vline(
    xintercept = c(
      0.50,
      0.75
    ),
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = c(
      0.50,
      0.70
    ),
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    alpha = 0.8
  ) +
  ggplot2::geom_text(
    size = 3,
    nudge_y = 0.025,
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Cluster donor dominance and donor entropy",
    subtitle = "Clusters 7, 8 and 13 show donor enrichment; cluster 15 is rare but donor-distributed",
    x = "Largest single-donor fraction",
    y = "Normalized donor entropy",
    size = "Cell count",
    shape = "Review flag"
  ) +
  theme_QC

save_plot_pair(
  plot_dominance,
  "QC2_cluster_donor_dominance_entropy",
  width = 8,
  height = 6
)

group_plot_data <- cluster_sample_full

group_plot_data <- aggregate(
  cell_n ~
    seurat_clusters +
    analysis_group,
  data = group_plot_data,
  FUN = sum
)

group_totals <- aggregate(
  cell_n ~ seurat_clusters,
  data = group_plot_data,
  FUN = sum
)

colnames(group_totals)[2] <-
  "cluster_cell_n"

group_plot_data <- merge(
  group_plot_data,
  group_totals,
  by = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE
)

group_plot_data$
  group_fraction <- (
    group_plot_data$cell_n /
      group_plot_data$
        cluster_cell_n
  )

group_plot_data$seurat_clusters <-
  factor(
    as.character(
      group_plot_data$
        seurat_clusters
    ),
    levels = cluster_levels
  )

plot_group <- ggplot2::ggplot(
  group_plot_data,
  ggplot2::aes(
    x = seurat_clusters,
    y = group_fraction,
    fill = analysis_group
  )
) +
  ggplot2::geom_col(
    position = "stack"
  ) +
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Control and RIF cell fractions by cluster",
    subtitle = "Descriptive composition only; no cell-level group inference",
    x = "Primary cluster",
    y = "Cell fraction",
    fill = "Group"
  ) +
  theme_QC

save_plot_pair(
  plot_group,
  "QC3_cluster_group_composition",
  width = 9,
  height = 5
)

technical_long <- data.frame(
  seurat_clusters = rep(
    cluster_review$
      seurat_clusters,
    times = 4
  ),
  metric = rep(
    c(
      "Median nFeature_RNA",
      "Median nCount_RNA",
      "Median percent.mt",
      "Median scDblFinder score"
    ),
    each = nrow(
      cluster_review
    )
  ),
  value = c(
    cluster_review$
      median_nFeature_RNA,
    cluster_review$
      median_nCount_RNA,
    cluster_review$
      median_percent_mt,
    cluster_review$
      median_scDblFinder_score
  ),
  stringsAsFactors = FALSE
)

technical_long$seurat_clusters <-
  factor(
    technical_long$
      seurat_clusters,
    levels = cluster_levels
  )

plot_technical <- ggplot2::ggplot(
  technical_long,
  ggplot2::aes(
    x = seurat_clusters,
    y = value,
    group = 1
  )
) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "Primary-cluster technical QC profiles",
    subtitle = "Low transcript complexity alone is not sufficient to remove an immune or rare cluster",
    x = "Primary cluster",
    y = "Median value"
  ) +
  theme_QC

save_plot_pair(
  plot_technical,
  "QC4_cluster_technical_profiles",
  width = 10,
  height = 7
)

mixing_plot_data <- cell_mixing

mixing_plot_data$seurat_clusters <-
  factor(
    mixing_plot_data$
      seurat_clusters,
    levels = cluster_levels
  )

plot_mixing <- ggplot2::ggplot(
  mixing_plot_data,
  ggplot2::aes(
    x = seurat_clusters,
    y = same_sample_excess
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_boxplot(
    outlier.shape = NA
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      quantile_or_NA(
        mixing_plot_data$
          same_sample_excess,
        0.01
      ),
      quantile_or_NA(
        mixing_plot_data$
          same_sample_excess,
        0.99
      )
    )
  ) +
  ggplot2::labs(
    title = "Local same-sample enrichment after Harmony",
    subtitle = "Observed same-sample neighbor fraction minus the expectation from cluster composition",
    x = "Primary cluster",
    y = "Same-sample neighborhood excess"
  ) +
  theme_QC

save_plot_pair(
  plot_mixing,
  "QC5_local_sample_mixing_by_cluster",
  width = 9,
  height = 5.5
)

transition_02_04_plot <- transition_02_04
transition_02_04_plot$
  cluster_res_0_2 <- factor(
    transition_02_04_plot$
      cluster_res_0_2,
    levels = rev(
      sort(
        unique(
          transition_02_04_plot$
            cluster_res_0_2
        )
      )
    )
  )

transition_02_04_plot$
  cluster_res_0_4 <- factor(
    transition_02_04_plot$
      cluster_res_0_4,
    levels = sort(
      unique(
        transition_02_04_plot$
          cluster_res_0_4
      )
    )
  )

plot_transition_02_04 <- ggplot2::ggplot(
  transition_02_04_plot,
  ggplot2::aes(
    x = cluster_res_0_4,
    y = cluster_res_0_2,
    fill = within_source_fraction
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Cluster transition from resolution 0.2 to 0.4",
    x = "Resolution 0.4 cluster",
    y = "Resolution 0.2 cluster",
    fill = "Fraction of 0.2 cluster"
  ) +
  theme_QC

save_plot_pair(
  plot_transition_02_04,
  "QC6_resolution_transition_0_2_to_0_4",
  width = 8,
  height = 6
)

transition_04_06_plot <- transition_04_06
transition_04_06_plot$
  cluster_res_0_4 <- factor(
    transition_04_06_plot$
      cluster_res_0_4,
    levels = rev(
      sort(
        unique(
          transition_04_06_plot$
            cluster_res_0_4
        )
      )
    )
  )

transition_04_06_plot$
  cluster_res_0_6 <- factor(
    transition_04_06_plot$
      cluster_res_0_6,
    levels = sort(
      unique(
        transition_04_06_plot$
          cluster_res_0_6
      )
    )
  )

plot_transition_04_06 <- ggplot2::ggplot(
  transition_04_06_plot,
  ggplot2::aes(
    x = cluster_res_0_6,
    y = cluster_res_0_4,
    fill = within_source_fraction
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Cluster transition from resolution 0.4 to 0.6",
    x = "Resolution 0.6 cluster",
    y = "Resolution 0.4 cluster",
    fill = "Fraction of 0.4 cluster"
  ) +
  theme_QC

save_plot_pair(
  plot_transition_04_06,
  "QC7_resolution_transition_0_4_to_0_6",
  width = 9,
  height = 6
)

cluster_size_plot_data <- cluster_review
cluster_size_plot_data$
  seurat_clusters <- factor(
    cluster_size_plot_data$
      seurat_clusters,
    levels = cluster_levels
  )

plot_cluster_size <- ggplot2::ggplot(
  cluster_size_plot_data,
  ggplot2::aes(
    x = seurat_clusters,
    y = cell_n,
    shape = review_required
  )
) +
  ggplot2::geom_hline(
    yintercept = 100,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_y_log10(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Primary-cluster sizes and review status",
    subtitle = "Cluster 15 is below 100 cells but is retained for marker review",
    x = "Primary cluster",
    y = "Cell count (log scale)",
    shape = "Review flag"
  ) +
  theme_QC

save_plot_pair(
  plot_cluster_size,
  "QC8_cluster_size_review",
  width = 8,
  height = 5
)

# Copy figure pairs into persistent result output.
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
# 11. Detailed integration-QC Gate
# ----------------------------------------------------------------------

critical_cluster_n <- sum(
  cluster_review$
    critical_hold_flag,
  na.rm = TRUE
)

review_cluster_ids <- cluster_review$
  seurat_clusters[
    cluster_review$
      review_required
  ]

severe_donor_dominant_clusters <-
  cluster_review$
    seurat_clusters[
      cluster_review$
        top_donor_fraction >=
        0.75
    ]

moderate_donor_enriched_clusters <-
  cluster_review$
    seurat_clusters[
      cluster_review$
        top_donor_fraction >=
        0.50 &
        cluster_review$
          top_donor_fraction <
          0.75
    ]

rare_clusters <- cluster_review$
  seurat_clusters[
    cluster_review$cell_n < 100
  ]

local_mixing_review_clusters <-
  as.character(
    mixing_cluster_summary$
      seurat_clusters[
        mixing_cluster_summary$
          local_mixing_review_flag
      ]
  )

ARI_02_04 <- resolution_ARI$
  adjusted_rand_index[
    resolution_ARI$comparison ==
      "0.2_vs_0.4"
  ]

ARI_04_06 <- resolution_ARI$
  adjusted_rand_index[
    resolution_ARI$comparison ==
      "0.4_vs_0.6"
  ]

resolution_review_flag <- (
  !is.finite(ARI_02_04) ||
    !is.finite(ARI_04_06) ||
    ARI_02_04 < 0.60 ||
    ARI_04_06 < 0.60
)

overall_review_flags <- c(
  if (length(severe_donor_dominant_clusters)) {
    paste0(
      "SEVERE_DONOR_DOMINANCE_CLUSTERS_",
      paste(
        severe_donor_dominant_clusters,
        collapse = "_"
      )
    )
  } else {
    character()
  },
  if (length(moderate_donor_enriched_clusters)) {
    paste0(
      "MODERATE_DONOR_ENRICHMENT_CLUSTERS_",
      paste(
        moderate_donor_enriched_clusters,
        collapse = "_"
      )
    )
  } else {
    character()
  },
  if (length(rare_clusters)) {
    paste0(
      "RARE_CLUSTERS_",
      paste(
        rare_clusters,
        collapse = "_"
      )
    )
  } else {
    character()
  },
  if (length(local_mixing_review_clusters)) {
    paste0(
      "LOCAL_MIXING_REVIEW_CLUSTERS_",
      paste(
        local_mixing_review_clusters,
        collapse = "_"
      )
    )
  } else {
    character()
  },
  if (resolution_review_flag) {
    "RESOLUTION_STABILITY_REVIEW"
  } else {
    character()
  }
)

overall_review_flags <- unique(
  overall_review_flags[
    nzchar(overall_review_flags)
  ]
)

gate_checks <- data.frame(
  check = c(
    "upstream_integration_gate_authorizes_detailed_QC",
    "checkpoint_cell_n_is_47161",
    "checkpoint_sample_n_is_9",
    "checkpoint_donor_n_is_9",
    "primary_cluster_n_is_16",
    "Harmony_graph_available",
    "Harmony_graph_dimensions_match_cells",
    "local_mixing_metrics_complete_for_all_cells",
    "all_clusters_represented_by_at_least_2_donors",
    "no_cluster_has_ge95_percent_single_donor_dominance",
    "no_cluster_has_fewer_than_30_cells",
    "no_cluster_meets_composite_low_quality_hold_rule",
    "no_cluster_has_median_doublet_score_above_0_30",
    "resolution_ARI_values_are_finite",
    "primary_resolution_stability_tables_complete",
    "target_clusters_7_8_13_15_reviewed",
    "marker_extraction_can_be_run_without_removing_review_clusters",
    "free_disk_after_run_at_least_8GB_or_unknown"
  ),
  pass = c(
    upstream_gate$gate_status[[1]] %in% c(
      "PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED",
      "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS"
    ),
    ncol(integrated_object) == 47161L,
    length(
      unique(
        metadata$sample_accession
      )
    ) == 9L,
    length(
      unique(
        metadata$donor_id
      )
    ) == 9L,
    length(
      unique(
        metadata$seurat_clusters
      )
    ) == 16L,
    graph_name %in%
      preferred_graph_names,
    nrow(mixing_graph) ==
      nrow(metadata) &&
      ncol(mixing_graph) ==
        nrow(metadata),
    nrow(cell_mixing) ==
      nrow(metadata) &&
      all(
        is.finite(
          cell_mixing$
            same_sample_fraction
        )
      ),
    all(
      cluster_review$
        represented_donor_n >= 2
    ),
    all(
      cluster_review$
        top_donor_fraction < 0.95
    ),
    all(
      cluster_review$cell_n >= 30
    ),
    !any(
      cluster_review$
        median_percent_mt > 22 &
        cluster_review$
          median_nFeature_RNA < 800
    ),
    all(
      cluster_review$
        median_scDblFinder_score <=
        0.30
    ),
    all(
      is.finite(
        resolution_ARI$
          adjusted_rand_index
      )
    ),
    nrow(
      primary_resolution_stability
    ) == 16L,
    all(
      target_review_clusters %in%
        target_cluster_review$
          seurat_clusters
    ),
    critical_cluster_n == 0L,
    is.na(
      get_free_disk_gb()
    ) ||
      get_free_disk_gb() >= 8
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

marker_analysis_authorized <- (
  all_gate_checks_pass &&
    critical_cluster_n == 0L
)

gate_status <- if (
  marker_analysis_authorized &&
    !length(overall_review_flags)
) {
  "PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED"
} else if (
  marker_analysis_authorized
) {
  "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED"
} else {
  "HOLD_DETAILED_INTEGRATION_QC_MARKERS_NOT_AUTHORIZED"
}

free_memory_after_GB <- get_free_memory_gb()
free_disk_after_GB <- get_free_disk_gb()

detailed_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  QC_version = QC_version,
  gate_status = gate_status,
  cell_n = nrow(metadata),
  donor_n = length(
    unique(
      metadata$donor_id
    )
  ),
  primary_cluster_n =
    length(
      unique(
        metadata$seurat_clusters
      )
    ),
  Harmony_graph = graph_name,
  severe_donor_dominant_clusters =
    collapse_unique(
      severe_donor_dominant_clusters
    ),
  moderate_donor_enriched_clusters =
    collapse_unique(
      moderate_donor_enriched_clusters
    ),
  rare_clusters =
    collapse_unique(
      rare_clusters
    ),
  local_mixing_review_clusters =
    collapse_unique(
      local_mixing_review_clusters
    ),
  ARI_res_0_2_vs_0_4 =
    ARI_02_04,
  ARI_res_0_4_vs_0_6 =
    ARI_04_06,
  critical_hold_cluster_n =
    critical_cluster_n,
  review_flags =
    collapse_unique(
      overall_review_flags
    ),
  descriptive_marker_analysis_authorized =
    marker_analysis_authorized,
  cluster_removal_authorized =
    FALSE,
  biological_annotation_locked =
    FALSE,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  detailed_gate,
  file.path(
    run_qc_root,
    "G183_DETAILED_QC_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1K GSE183837 Detailed Integration QC Card",
  "",
  paste0("- Dataset: ", dataset_id),
  paste0(
    "- Integration version: ",
    integration_version
  ),
  paste0("- QC version: ", QC_version),
  paste0("- Gate: ", gate_status),
  paste0(
    "- Cells: ",
    format(
      nrow(metadata),
      big.mark = ","
    )
  ),
  paste0(
    "- Donors: ",
    length(
      unique(metadata$donor_id)
    )
  ),
  paste0(
    "- Primary clusters: ",
    length(
      unique(
        metadata$seurat_clusters
      )
    )
  ),
  paste0(
    "- Severe donor-dominant clusters: ",
    if (
      length(
        severe_donor_dominant_clusters
      )
    ) {
      collapse_unique(
        severe_donor_dominant_clusters
      )
    } else {
      "NONE"
    }
  ),
  paste0(
    "- Moderate donor-enriched clusters: ",
    if (
      length(
        moderate_donor_enriched_clusters
      )
    ) {
      collapse_unique(
        moderate_donor_enriched_clusters
      )
    } else {
      "NONE"
    }
  ),
  paste0(
    "- Rare clusters: ",
    if (length(rare_clusters)) {
      collapse_unique(rare_clusters)
    } else {
      "NONE"
    }
  ),
  paste0(
    "- Resolution ARI 0.2 vs 0.4: ",
    round(ARI_02_04, 4)
  ),
  paste0(
    "- Resolution ARI 0.4 vs 0.6: ",
    round(ARI_04_06, 4)
  ),
  paste0(
    "- Marker analysis authorized: ",
    marker_analysis_authorized
  ),
  "- Cluster removal authorized: FALSE",
  "- Biological annotation locked: FALSE",
  "- Disease pseudobulk testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "G183_DETAILED_QC_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 detailed integration-QC interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Fixed interpretation rules",
  "",
  "- Cluster 13 is retained but must not support a disease claim unless its identity and signal are reproduced across multiple donors.",
  "- Clusters 7 and 8 are retained for marker extraction and reviewed for donor-associated states.",
  "- Cluster 15 is retained as a rare candidate population because rarity alone is not a deletion criterion.",
  "- Low RNA complexity alone is not treated as low quality when mitochondrial and doublet indicators are not concordantly abnormal.",
  "- Control-versus-RIF visualization is descriptive only; cell-level hypothesis testing remains prohibited.",
  "",
  "## Authorized next operation",
  "",
  "- Run memory-safe descriptive cluster marker extraction.",
  "- Review canonical lineage markers and the marker profiles of clusters 7, 8, 13 and 15.",
  "- Build a candidate broad-cell annotation map.",
  "",
  "## Not authorized yet",
  "",
  "- No cluster removal.",
  "- No biological annotation lock.",
  "- No disease pseudobulk testing.",
  "- No result lock."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "G183_DETAILED_QC_INTERPRETATION.md"
  )
)

QC_notes <- c(
  "# STEP04E1K QC notes",
  "",
  paste0("Gate: ", gate_status),
  "",
  "The integration checkpoint itself is unchanged.",
  "All review clusters remain in the object.",
  "Marker extraction is descriptive and cluster-wise.",
  "Disease inference remains donor-level and is deferred until annotation QC.",
  "",
  "Known pre-run composition from STEP04E1J:",
  "- Cluster 13: 402 cells; strongly dominated by GSM5572242.",
  "- Clusters 7 and 8: donor-enriched, particularly GSM5572241.",
  "- Cluster 15: 87 cells distributed across all nine donors.",
  "",
  paste0(
    "Detailed review flags: ",
    if (
      length(overall_review_flags)
    ) {
      collapse_unique(
        overall_review_flags
      )
    } else {
      "NONE"
    }
  )
)

write_text(
  QC_notes,
  file.path(
    run_qc_root,
    "QC_NOTES.md"
  )
)


# ----------------------------------------------------------------------
# 12. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_QC_versions = c(
    "INTQC183_v1_0_CELL_N_COLUMN_COLLISION",
    "INTQC183_v1_1_WINDOWS_PATH_LIMIT"
  ),
  merge_policy =
    "COMPOSITION_CELL_N_AND_TECHNICAL_CELL_N_EXPLICITLY_DISTINCT",
  path_policy =
    "SHORT_PERSISTENT_DIRECTORIES_AND_FILENAMES_MAX_235_CHARACTERS",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  QC_version = QC_version,
  upstream_gate =
    upstream_gate$gate_status[[1]],
  expected_cell_n = 47161L,
  expected_donor_n = 9L,
  expected_primary_cluster_n = 16L,
  graph_preference = c(
    "harmony_nn",
    "harmony_snn"
  ),
  donor_dominance_review_threshold = 0.50,
  donor_dominance_severe_threshold = 0.75,
  donor_dominance_hold_threshold = 0.95,
  rare_cluster_review_threshold = 100L,
  minimum_cluster_hold_threshold = 30L,
  donor_entropy_review_thresholds =
    c(0.50, 0.70),
  local_mixing_review_rules = list(
    median_same_sample_excess =
      0.25,
    Q90_same_sample_excess =
      0.50,
    median_neighbor_donor_n =
      3
  ),
  resolution_ARI_review_threshold =
    0.60,
  marker_analysis_policy =
    "DESCRIPTIVE_CLUSTER_MARKERS_ONLY",
  cluster_removal_authorized =
    FALSE,
  disease_cell_level_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  random_seed = 20260806L
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_QC_versions =
    "INTQC183_v1_0_CELL_N_COLUMN_COLLISION;INTQC183_v1_1_WINDOWS_PATH_LIMIT",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  QC_version = QC_version,
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
  cell_n = nrow(metadata),
  donor_n = length(
    unique(
      metadata$donor_id
    )
  ),
  primary_cluster_n =
    length(
      unique(
        metadata$seurat_clusters
      )
    ),
  severe_donor_dominant_clusters =
    collapse_unique(
      severe_donor_dominant_clusters
    ),
  moderate_donor_enriched_clusters =
    collapse_unique(
      moderate_donor_enriched_clusters
    ),
  rare_clusters =
    collapse_unique(
      rare_clusters
    ),
  marker_analysis_authorized =
    marker_analysis_authorized,
  cluster_removal_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  free_memory_before_GB =
    free_memory_before_GB,
  free_memory_after_GB =
    free_memory_after_GB,
  free_disk_before_GB =
    free_disk_before_GB,
  free_disk_after_GB =
    free_disk_after_GB,
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

persistent_files <- unique(c(
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
  )
))

persistent_files <- persistent_files[
  file.exists(persistent_files) &
    !dir.exists(persistent_files)
]

persistent_manifest <- data.frame(
  path = normalize_slashes(
    persistent_files
  ),
  size_bytes =
    file.info(persistent_files)$size,
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
    file.info(small_files)$size,
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
  "marker_analysis_authorized=",
  marker_analysis_authorized,
  "review_flags=",
  if (
    length(overall_review_flags)
  ) {
    collapse_unique(
      overall_review_flags
    )
  } else {
    "NONE"
  }
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
cat("STEP04E1K GSE183837 detailed integration QC completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Cells: ",
  format(
    nrow(metadata),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Donors: ",
  length(
    unique(
      metadata$donor_id
    )
  ),
  "\n",
  sep = ""
)
cat(
  "Primary clusters: ",
  length(
    unique(
      metadata$seurat_clusters
    )
  ),
  "\n",
  sep = ""
)
cat(
  "Severe donor-dominant clusters: ",
  if (
    length(
      severe_donor_dominant_clusters
    )
  ) {
    collapse_unique(
      severe_donor_dominant_clusters
    )
  } else {
    "NONE"
  },
  "\n",
  sep = ""
)
cat(
  "Moderate donor-enriched clusters: ",
  if (
    length(
      moderate_donor_enriched_clusters
    )
  ) {
    collapse_unique(
      moderate_donor_enriched_clusters
    )
  } else {
    "NONE"
  },
  "\n",
  sep = ""
)
cat(
  "Rare clusters: ",
  if (length(rare_clusters)) {
    collapse_unique(rare_clusters)
  } else {
    "NONE"
  },
  "\n",
  sep = ""
)
cat(
  "ARI 0.2 vs 0.4: ",
  round(ARI_02_04, 4),
  "\n",
  sep = ""
)
cat(
  "ARI 0.4 vs 0.6: ",
  round(ARI_04_06, 4),
  "\n",
  sep = ""
)
cat(
  "Marker analysis authorized: ",
  marker_analysis_authorized,
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
    "\nThe detailed QC outputs are complete in the project folders.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 4: 66_step04E1L_GSE183837_descriptive_markers_v1_1_panel_duplicate_fixed
# ==============================================================================

# ======================================================================
# 66_step04E1L_GSE183837_descriptive_markers_v1_1_panel_duplicate_fixed.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Extract descriptive positive markers for all 16 GSE183837 clusters.
# - Run each cluster sequentially with per-cluster checkpoints and resume.
# - Use a capped, reproducible Wilcoxon comparison for technical stability.
# - Quantify donor support for the leading markers of every cluster.
# - Summarize canonical lineage-marker expression across clusters.
# - Generate a candidate broad-lineage annotation table for manual review.
#
# UPSTREAM GATE
# - PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED
#
# SCIENTIFIC GUARDRAILS
# - Markers are descriptive cluster markers, not disease DE results.
# - Cell-level Control-versus-RIF testing is prohibited.
# - Cluster 13 markers cannot support a disease claim without donor
#   replication.
# - Clusters 7, 8, 13 and 15 remain mandatory review clusters.
# - Candidate annotations are not locked biological identities.
# - No cell or cluster is removed.
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
options(future.globals.maxSize = 8 * 1024^3)

set.seed(20260806)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE183837"

script_name <-
  "66_step04E1L_GSE183837_descriptive_markers_v1_1_panel_duplicate_fixed.R"

script_version <-
  "STEP04E1L_GSE183837_DESCRIPTIVE_MARKERS_v1.1_PANEL_DUPLICATE_FIXED"

integration_version <- "SCI183_v1_0"
marker_version <- "MK183_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1L_",
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

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  path_n <- nchar(
    normalize_slashes(path),
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
      normalize_slashes(path)
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
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

safe_package_version <- function(package_name) {
  if (
    !requireNamespace(
      package_name,
      quietly = TRUE
    )
  ) {
    return(NA_character_)
  }

  as.character(
    utils::packageVersion(package_name)
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

median_or_NA <- function(x) {
  x <- suppressWarnings(
    as.numeric(as.character(x))
  )

  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    median(x)
  }
}

safe_zscore <- function(x) {
  x <- as.numeric(x)

  if (
    length(x) < 2L ||
      !is.finite(stats::sd(x)) ||
      stats::sd(x) == 0
  ) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
}

detect_logFC_column <- function(column_names) {
  preferred <- c(
    "avg_log2FC",
    "avg_logFC",
    "avg_diff"
  )

  hit <- preferred[
    preferred %in% column_names
  ]

  if (length(hit)) {
    return(hit[1])
  }

  pattern_hit <- grep(
    "log.*fc|avg.*diff",
    column_names,
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(pattern_hit)) {
    return(pattern_hit[1])
  }

  ""
}

standardize_marker_table <- function(
  marker_table,
  cluster_id
) {
  if (
    is.null(marker_table) ||
      !nrow(marker_table)
  ) {
    return(
      data.frame(
        cluster = character(),
        gene = character(),
        p_val = numeric(),
        avg_log2FC = numeric(),
        pct.1 = numeric(),
        pct.2 = numeric(),
        p_val_adj = numeric(),
        rank = integer(),
        stringsAsFactors = FALSE
      )
    )
  }

  output <- as.data.frame(
    marker_table,
    stringsAsFactors = FALSE
  )

  output$gene <- rownames(output)
  rownames(output) <- NULL

  logFC_column <- detect_logFC_column(
    colnames(output)
  )

  if (!nzchar(logFC_column)) {
    stop(
      "Could not identify a log-fold-change column for cluster ",
      cluster_id,
      "."
    )
  }

  required_columns <- c(
    "p_val",
    "pct.1",
    "pct.2",
    "p_val_adj"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(output)
  )

  if (length(missing_columns)) {
    stop(
      "Marker table for cluster ",
      cluster_id,
      " is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  output$avg_log2FC <- as.numeric(
    output[[logFC_column]]
  )

  output$p_val <- as.numeric(
    output$p_val
  )

  output$pct.1 <- as.numeric(
    output$pct.1
  )

  output$pct.2 <- as.numeric(
    output$pct.2
  )

  output$p_val_adj <- as.numeric(
    output$p_val_adj
  )

  output$cluster <- as.character(
    cluster_id
  )

  output <- output[
    order(
      output$p_val_adj,
      -output$avg_log2FC,
      output$gene
    ),
    ,
    drop = FALSE
  ]

  output$rank <- seq_len(
    nrow(output)
  )

  preferred_order <- c(
    "cluster",
    "gene",
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "rank"
  )

  remaining_columns <- setdiff(
    colnames(output),
    preferred_order
  )

  output[
    ,
    c(
      preferred_order,
      remaining_columns
    ),
    drop = FALSE
  ]
}

resolve_features <- function(
  requested_features,
  available_features
) {
  # Preserve input length and order. The same canonical gene may
  # legitimately appear in more than one lineage panel.
  requested_features <- as.character(
    requested_features
  )

  available_features <- as.character(
    available_features
  )

  available_upper <- toupper(
    available_features
  )

  requested_upper <- toupper(
    requested_features
  )

  match_index <- match(
    requested_upper,
    available_upper
  )

  output <- data.frame(
    request_index =
      seq_along(requested_features),
    requested_gene =
      requested_features,
    resolved_feature = ifelse(
      is.na(match_index),
      NA_character_,
      available_features[match_index]
    ),
    available = !is.na(match_index),
    stringsAsFactors = FALSE
  )

  if (
    nrow(output) !=
      length(requested_features)
  ) {
    stop(
      "Feature-resolution output length does not match the requested panel length."
    )
  }

  output
}

get_data_matrix <- function(
  object,
  features
) {
  features <- unique(
    features[
      features %in%
        rownames(object)
    ]
  )

  if (!length(features)) {
    return(NULL)
  }

  output <- SeuratObject::LayerData(
    object = object,
    assay = "RNA",
    layer = "data",
    features = features
  )

  output[
    features,
    ,
    drop = FALSE
  ]
}

calculate_donor_marker_support <- function(
  object,
  metadata,
  cluster_id,
  marker_genes,
  minimum_cluster_cells_per_donor = 3L
) {
  marker_genes <- unique(
    marker_genes[
      marker_genes %in%
        rownames(object)
    ]
  )

  if (!length(marker_genes)) {
    return(
      data.frame(
        cluster = character(),
        gene = character(),
        donor_id = character(),
        analysis_group = character(),
        cluster_cell_n = integer(),
        outside_cell_n = integer(),
        pct_cluster = numeric(),
        pct_outside = numeric(),
        pct_delta = numeric(),
        avg_cluster = numeric(),
        avg_outside = numeric(),
        avg_delta = numeric(),
        eligible_donor = logical(),
        donor_supports_marker = logical(),
        stringsAsFactors = FALSE
      )
    )
  }

  expression_matrix <- get_data_matrix(
    object,
    marker_genes
  )

  donor_ids <- unique(
    as.character(
      metadata$donor_id
    )
  )

  output_rows <- list()
  row_index <- 0L

  for (donor_id_value in donor_ids) {
    donor_index <- which(
      metadata$donor_id ==
        donor_id_value
    )

    cluster_index <- donor_index[
      as.character(
        metadata$seurat_clusters[
          donor_index
        ]
      ) == as.character(cluster_id)
    ]

    outside_index <- donor_index[
      as.character(
        metadata$seurat_clusters[
          donor_index
        ]
      ) != as.character(cluster_id)
    ]

    cluster_cell_n <- length(
      cluster_index
    )

    outside_cell_n <- length(
      outside_index
    )

    analysis_group_value <- unique(
      as.character(
        metadata$analysis_group[
          donor_index
        ]
      )
    )

    if (
      length(analysis_group_value) != 1L
    ) {
      stop(
        "Donor/group identity is not unique for ",
        donor_id_value,
        "."
      )
    }

    if (cluster_cell_n > 0L) {
      pct_cluster <- Matrix::rowMeans(
        expression_matrix[
          ,
          cluster_index,
          drop = FALSE
        ] > 0
      )

      avg_cluster <- Matrix::rowMeans(
        expression_matrix[
          ,
          cluster_index,
          drop = FALSE
        ]
      )
    } else {
      pct_cluster <- rep(
        NA_real_,
        length(marker_genes)
      )

      avg_cluster <- rep(
        NA_real_,
        length(marker_genes)
      )
    }

    if (outside_cell_n > 0L) {
      pct_outside <- Matrix::rowMeans(
        expression_matrix[
          ,
          outside_index,
          drop = FALSE
        ] > 0
      )

      avg_outside <- Matrix::rowMeans(
        expression_matrix[
          ,
          outside_index,
          drop = FALSE
        ]
      )
    } else {
      pct_outside <- rep(
        NA_real_,
        length(marker_genes)
      )

      avg_outside <- rep(
        NA_real_,
        length(marker_genes)
      )
    }

    eligible_donor <- (
      cluster_cell_n >=
        minimum_cluster_cells_per_donor &&
        outside_cell_n >= 20L
    )

    donor_supports_marker <- (
      eligible_donor &
        is.finite(pct_cluster) &
        is.finite(pct_outside) &
        pct_cluster -
          pct_outside >= 0.10
    )

    row_index <- row_index + 1L

    output_rows[[row_index]] <- data.frame(
      cluster = as.character(
        cluster_id
      ),
      gene = marker_genes,
      donor_id = donor_id_value,
      analysis_group =
        analysis_group_value,
      cluster_cell_n =
        cluster_cell_n,
      outside_cell_n =
        outside_cell_n,
      pct_cluster =
        as.numeric(pct_cluster),
      pct_outside =
        as.numeric(pct_outside),
      pct_delta =
        as.numeric(
          pct_cluster -
            pct_outside
        ),
      avg_cluster =
        as.numeric(avg_cluster),
      avg_outside =
        as.numeric(avg_outside),
      avg_delta =
        as.numeric(
          avg_cluster -
            avg_outside
        ),
      eligible_donor =
        eligible_donor,
      donor_supports_marker =
        as.logical(
          donor_supports_marker
        ),
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    output_rows
  )
}

summarize_donor_marker_support <- function(
  support_long
) {
  if (
    is.null(support_long) ||
      !nrow(support_long)
  ) {
    return(
      data.frame(
        cluster = character(),
        gene = character(),
        eligible_donor_n = integer(),
        supporting_donor_n = integer(),
        supporting_donor_fraction = numeric(),
        median_pct_delta = numeric(),
        median_avg_delta = numeric(),
        Control_supporting_donor_n = integer(),
        RIF_supporting_donor_n = integer(),
        stringsAsFactors = FALSE
      )
    )
  }

  support_DT <- data.table::as.data.table(
    support_long
  )

  output <- support_DT[
    ,
    .(
      eligible_donor_n =
        sum(
          eligible_donor,
          na.rm = TRUE
        ),
      supporting_donor_n =
        sum(
          donor_supports_marker,
          na.rm = TRUE
        ),
      supporting_donor_fraction =
        if (
          sum(
            eligible_donor,
            na.rm = TRUE
          ) > 0L
        ) {
          sum(
            donor_supports_marker,
            na.rm = TRUE
          ) /
            sum(
              eligible_donor,
              na.rm = TRUE
            )
        } else {
          NA_real_
        },
      median_pct_delta =
        median_or_NA(
          pct_delta[
            eligible_donor
          ]
        ),
      median_avg_delta =
        median_or_NA(
          avg_delta[
            eligible_donor
          ]
        ),
      Control_supporting_donor_n =
        sum(
          donor_supports_marker &
            analysis_group ==
              "CONTROL",
          na.rm = TRUE
        ),
      RIF_supporting_donor_n =
        sum(
          donor_supports_marker &
            analysis_group ==
              "RIF",
          na.rm = TRUE
        )
    ),
    by = .(
      cluster,
      gene
    )
  ]

  as.data.frame(output)
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
  "Seurat",
  "SeuratObject",
  "Matrix",
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

parameter_root <- file.path(
  step04_root,
  "03_parameters"
)

environment_root <- file.path(
  step04_root,
  "04_environment"
)

integration_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

input_checkpoint <- file.path(
  integration_root,
  "checkpoints",
  dataset_id,
  "integrated.rds"
)

output_root <- file.path(
  integration_root,
  "mk",
  "G183",
  marker_version
)

model_root <- file.path(
  output_root,
  "m"
)

full_root <- file.path(
  output_root,
  "full"
)

top_root <- file.path(
  output_root,
  "top"
)

panel_root <- file.path(
  output_root,
  "panel"
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
    "STEP04E1L_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_markers"
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
  output_root,
  model_root,
  full_root,
  top_root,
  panel_root,
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

planned_files <- data.frame(
  artifact = c(
    "input_checkpoint",
    "output_root",
    "cluster_00_model",
    "cluster_15_model",
    "combined_markers",
    "combined_top100",
    "donor_support",
    "canonical_panel",
    "candidate_annotation",
    "gate",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    input_checkpoint,
    output_root,
    file.path(model_root, "C00.rds"),
    file.path(model_root, "C15.rds"),
    file.path(full_root, "all_markers.csv"),
    file.path(top_root, "all_top100.csv"),
    file.path(top_root, "donor_support_summary.csv"),
    file.path(panel_root, "canonical_panel.csv"),
    file.path(panel_root, "candidate_annotation.csv"),
    file.path(run_qc_root, "G183_MARKER_GATE.csv"),
    handoff_zip
  )),
  stringsAsFactors = FALSE
)

planned_files$character_n <- nchar(
  planned_files$path,
  type = "chars"
)

planned_files$path_safe <- (
  planned_files$character_n <= 235L
)

write_csv_utf8(
  planned_files,
  file.path(
    run_qc_root,
    "00_path_audit.csv"
  )
)

if (!all(planned_files$path_safe)) {
  stop(
    "One or more planned marker paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream detailed-QC Gate
# ----------------------------------------------------------------------

detailed_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

detailed_qc_dirs <- detailed_qc_dirs[
  grepl(
    "^STEP04E1K_[0-9]{8}_[0-9]{6}_dQC$",
    basename(detailed_qc_dirs),
    perl = TRUE
  )
]

if (!length(detailed_qc_dirs)) {
  stop(
    "No STEP04E1K detailed-QC directory was found."
  )
}

gate_candidates <- list()

for (directory_path in detailed_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "G183_DETAILED_QC_GATE.csv"
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

    gate_candidates[[
      length(gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(gate_candidates)) {
  stop(
    "No readable STEP04E1K detailed-QC Gate was found."
  )
}

gate_candidates <- as.data.frame(
  data.table::rbindlist(
    gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    gate_candidates$source_mtime
  )

gate_candidates <- gate_candidates[
  order(
    gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

authorized_index <- which(
  gate_candidates$gate_status %in% c(
    "PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED",
    "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED"
  ) &
    toupper(
      as.character(
        gate_candidates$
          descriptive_marker_analysis_authorized
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      )
)

if (!length(authorized_index)) {
  stop(
    "No detailed-QC Gate authorizes marker extraction."
  )
}

upstream_gate <- gate_candidates[
  authorized_index[1],
  ,
  drop = FALSE
]

write_csv_utf8(
  upstream_gate,
  file.path(
    run_qc_root,
    "00_upstream_gate.csv"
  )
)


# ----------------------------------------------------------------------
# 5. Load and validate integrated object
# ----------------------------------------------------------------------

if (!file.exists(input_checkpoint)) {
  stop(
    "Integrated checkpoint is missing:\n",
    input_checkpoint
  )
}

log_event(
  "Starting",
  run_id,
  "checkpoint=",
  normalize_slashes(
    input_checkpoint
  )
)

integrated_object <- readRDS(
  input_checkpoint
)

if (!inherits(integrated_object, "Seurat")) {
  stop(
    "Input checkpoint is not a Seurat object."
  )
}

SeuratObject::DefaultAssay(
  integrated_object
) <- "RNA"

metadata <- integrated_object[[]]

required_metadata <- c(
  "sample_accession",
  "donor_id",
  "analysis_group",
  "seurat_clusters"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(metadata)
)

if (length(missing_metadata)) {
  stop(
    "Integrated metadata is missing: ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}

if (
  ncol(integrated_object) != 47161L ||
    length(
      unique(
        as.character(
          metadata$donor_id
        )
      )
    ) != 9L ||
    length(
      unique(
        as.character(
          metadata$seurat_clusters
        )
      )
    ) != 16L
) {
  stop(
    "Integrated object identity does not match the locked 47,161-cell object."
  )
}

if (
  anyDuplicated(
    rownames(integrated_object)
  ) > 0L
) {
  stop(
    "Integrated feature names are not unique."
  )
}

data_layers <- SeuratObject::Layers(
  integrated_object,
  assay = "RNA",
  search = "^data"
)

count_layers <- SeuratObject::Layers(
  integrated_object,
  assay = "RNA",
  search = "^counts"
)

if (
  length(data_layers) != 1L ||
    data_layers[1] != "data" ||
    length(count_layers) != 9L
) {
  stop(
    "Expected one joined data layer and nine split count layers."
  )
}

metadata$cell_barcode <- rownames(
  metadata
)

metadata$donor_id <- as.character(
  metadata$donor_id
)

metadata$analysis_group <- as.character(
  metadata$analysis_group
)

metadata$seurat_clusters <- factor(
  as.character(
    metadata$seurat_clusters
  ),
  levels = as.character(0:15)
)

if (anyNA(metadata$seurat_clusters)) {
  stop(
    "Unexpected primary cluster identities."
  )
}

SeuratObject::Idents(
  integrated_object
) <- metadata$
  seurat_clusters

checkpoint_audit <- data.frame(
  metric = c(
    "feature_n",
    "cell_n",
    "donor_n",
    "cluster_n",
    "data_layer_n",
    "count_layer_n",
    "checkpoint_size_GB",
    "checkpoint_md5"
  ),
  value = c(
    nrow(integrated_object),
    ncol(integrated_object),
    length(
      unique(
        metadata$donor_id
      )
    ),
    length(
      unique(
        metadata$seurat_clusters
      )
    ),
    length(data_layers),
    length(count_layers),
    round(
      file.info(
        input_checkpoint
      )$size / 1024^3,
      3
    ),
    unname(
      tools::md5sum(
        input_checkpoint
      )
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  checkpoint_audit,
  file.path(
    run_qc_root,
    "01_checkpoint_audit.csv"
  )
)


# ----------------------------------------------------------------------
# 6. Sequential resumable marker extraction
# ----------------------------------------------------------------------

cluster_ids <- as.character(0:15)

marker_parameters <- list(
  test_use = "wilcox",
  only_positive = TRUE,
  minimum_pct = 0.10,
  logFC_threshold = 0.25,
  maximum_cells_per_identity = 5000L,
  minimum_cells_feature = 3L,
  minimum_cells_group = 3L,
  donor_support_top_n = 30L,
  top_output_n = 100L,
  random_seed = 20260806L
)

cluster_status_rows <- list()
cluster_marker_objects <- list()
donor_support_long_rows <- list()
donor_support_summary_rows <- list()

for (
  cluster_index in seq_along(
    cluster_ids
  )
) {
  cluster_id <- cluster_ids[
    cluster_index
  ]

  cluster_code <- sprintf(
    "C%02d",
    as.integer(cluster_id)
  )

  checkpoint_path <- file.path(
    model_root,
    paste0(
      cluster_code,
      ".rds"
    )
  )

  status_path <- file.path(
    model_root,
    paste0(
      cluster_code,
      "_status.csv"
    )
  )

  full_path <- file.path(
    full_root,
    paste0(
      cluster_code,
      ".csv"
    )
  )

  top_path <- file.path(
    top_root,
    paste0(
      cluster_code,
      "_top100.csv"
    )
  )

  donor_long_path <- file.path(
    top_root,
    paste0(
      cluster_code,
      "_donor_long.csv"
    )
  )

  donor_summary_path <- file.path(
    top_root,
    paste0(
      cluster_code,
      "_donor_summary.csv"
    )
  )

  existing_status <- safe_read_csv(
    status_path
  )

  resume_valid <- FALSE
  cluster_result <- NULL

  if (
    file.exists(checkpoint_path) &&
      !is.null(existing_status) &&
      nrow(existing_status) == 1L &&
      existing_status$status[[1]] ==
        "COMPLETE"
  ) {
    cluster_result <- tryCatch(
      readRDS(checkpoint_path),
      error = function(e) NULL
    )

    resume_valid <- (
      is.list(cluster_result) &&
        all(
          c(
            "markers",
            "donor_support_long",
            "donor_support_summary"
          ) %in%
            names(cluster_result)
        ) &&
        identical(
          as.character(
            cluster_result$cluster
          ),
          cluster_id
        )
    )
  }

  if (resume_valid) {
    log_event(
      "SKIP completed marker cluster",
      cluster_id,
      "marker_rows=",
      nrow(
        cluster_result$markers
      )
    )
  } else {
    cluster_start <- Sys.time()

    log_event(
      "START marker cluster",
      cluster_id,
      "(",
      cluster_index,
      "/",
      length(cluster_ids),
      ")"
    )

    set.seed(
      marker_parameters$
        random_seed +
        as.integer(cluster_id)
    )

    raw_marker_result <- tryCatch(
      Seurat::FindMarkers(
        object = integrated_object,
        ident.1 = cluster_id,
        ident.2 = NULL,
        assay = "RNA",
        slot = "data",
        test.use =
          marker_parameters$
            test_use,
        only.pos =
          marker_parameters$
            only_positive,
        min.pct =
          marker_parameters$
            minimum_pct,
        logfc.threshold =
          marker_parameters$
            logFC_threshold,
        min.cells.feature =
          marker_parameters$
            minimum_cells_feature,
        min.cells.group =
          marker_parameters$
            minimum_cells_group,
        max.cells.per.ident =
          marker_parameters$
            maximum_cells_per_identity,
        random.seed =
          marker_parameters$
            random_seed +
            as.integer(cluster_id),
        densify = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        error_status <- data.frame(
          cluster = cluster_id,
          status = "ERROR",
          error_message =
            conditionMessage(e),
          stringsAsFactors = FALSE
        )

        write_csv_utf8(
          error_status,
          status_path
        )

        stop(
          "FindMarkers failed for cluster ",
          cluster_id,
          ": ",
          conditionMessage(e)
        )
      }
    )

    standardized_markers <-
      standardize_marker_table(
        raw_marker_result,
        cluster_id
      )

    if (!nrow(standardized_markers)) {
      stop(
        "No positive marker rows were returned for cluster ",
        cluster_id,
        "."
      )
    }

    top100_markers <- head(
      standardized_markers,
      marker_parameters$
        top_output_n
    )

    donor_marker_genes <- head(
      standardized_markers$gene,
      marker_parameters$
        donor_support_top_n
    )

    donor_support_long <-
      calculate_donor_marker_support(
        object =
          integrated_object,
        metadata =
          metadata,
        cluster_id =
          cluster_id,
        marker_genes =
          donor_marker_genes,
        minimum_cluster_cells_per_donor =
          3L
      )

    donor_support_summary <-
      summarize_donor_marker_support(
        donor_support_long
      )

    standardized_markers <- merge(
      standardized_markers,
      donor_support_summary,
      by = c(
        "cluster",
        "gene"
      ),
      all.x = TRUE,
      sort = FALSE
    )

    standardized_markers <- standardized_markers[
      order(
        standardized_markers$rank
      ),
      ,
      drop = FALSE
    ]

    standardized_markers$
      robust_donor_supported <- (
        standardized_markers$
          p_val_adj < 0.05 &
          standardized_markers$
            avg_log2FC >= 0.50 &
          standardized_markers$
            pct.1 >= 0.20 &
          is.finite(
            standardized_markers$
              supporting_donor_fraction
          ) &
          standardized_markers$
            supporting_donor_fraction >=
            0.50 &
          standardized_markers$
            supporting_donor_n >= 2
      )

    top100_markers <- head(
      standardized_markers,
      marker_parameters$
        top_output_n
    )

    cluster_result <- list(
      cluster = cluster_id,
      markers =
        standardized_markers,
      top100 =
        top100_markers,
      donor_support_long =
        donor_support_long,
      donor_support_summary =
        donor_support_summary,
      parameters =
        marker_parameters
    )

    saveRDS(
      cluster_result,
      checkpoint_path,
      compress = "gzip"
    )

    write_csv_utf8(
      standardized_markers,
      full_path
    )

    write_csv_utf8(
      top100_markers,
      top_path
    )

    write_csv_utf8(
      donor_support_long,
      donor_long_path
    )

    write_csv_utf8(
      donor_support_summary,
      donor_summary_path
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

    cluster_status <- data.frame(
      cluster = cluster_id,
      status = "COMPLETE",
      marker_row_n =
        nrow(
          standardized_markers
        ),
      FDR05_marker_n =
        sum(
          standardized_markers$
            p_val_adj < 0.05,
          na.rm = TRUE
        ),
      robust_donor_supported_n =
        sum(
          standardized_markers$
            robust_donor_supported,
          na.rm = TRUE
        ),
      donor_support_gene_n =
        nrow(
          donor_support_summary
        ),
      elapsed_minutes =
        elapsed_minutes,
      checkpoint_md5 =
        unname(
          tools::md5sum(
            checkpoint_path
          )
        ),
      stringsAsFactors = FALSE
    )

    write_csv_utf8(
      cluster_status,
      status_path
    )

    existing_status <- cluster_status

    log_event(
      "DONE marker cluster",
      cluster_id,
      "rows=",
      nrow(
        standardized_markers
      ),
      "FDR05=",
      sum(
        standardized_markers$
          p_val_adj < 0.05,
        na.rm = TRUE
      ),
      "robust=",
      sum(
        standardized_markers$
          robust_donor_supported,
        na.rm = TRUE
      ),
      "elapsed_min=",
      elapsed_minutes
    )

    rm(
      raw_marker_result,
      standardized_markers,
      donor_support_long,
      donor_support_summary
    )

    gc()
  }

  if (is.null(existing_status)) {
    existing_status <- safe_read_csv(
      status_path
    )
  }

  cluster_status_rows[[cluster_index]] <-
    existing_status

  cluster_marker_objects[[cluster_index]] <-
    cluster_result$markers

  donor_support_long_rows[[cluster_index]] <-
    cluster_result$donor_support_long

  donor_support_summary_rows[[cluster_index]] <-
    cluster_result$donor_support_summary
}

cluster_status_all <- do.call(
  rbind,
  cluster_status_rows
)

cluster_status_all <- cluster_status_all[
  match(
    cluster_ids,
    cluster_status_all$cluster
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  cluster_status_all,
  file.path(
    output_root,
    "01_marker_status.csv"
  )
)

all_markers <- as.data.frame(
  data.table::rbindlist(
    cluster_marker_objects,
    use.names = TRUE,
    fill = TRUE
  )
)

all_markers$cluster <- as.character(
  all_markers$cluster
)

all_markers <- all_markers[
  order(
    as.integer(
      all_markers$cluster
    ),
    all_markers$rank
  ),
  ,
  drop = FALSE
]

all_top100 <- all_markers[
  all_markers$rank <=
    marker_parameters$
      top_output_n,
  ,
  drop = FALSE
]

donor_support_long_all <- as.data.frame(
  data.table::rbindlist(
    donor_support_long_rows,
    use.names = TRUE,
    fill = TRUE
  )
)

donor_support_summary_all <-
  as.data.frame(
    data.table::rbindlist(
      donor_support_summary_rows,
      use.names = TRUE,
      fill = TRUE
    )
  )

write_csv_utf8(
  all_markers,
  file.path(
    full_root,
    "all_markers.csv"
  )
)

write_csv_utf8(
  all_top100,
  file.path(
    top_root,
    "all_top100.csv"
  )
)

write_csv_utf8(
  donor_support_long_all,
  file.path(
    top_root,
    "donor_support_long.csv"
  )
)

write_csv_utf8(
  donor_support_summary_all,
  file.path(
    top_root,
    "donor_support_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Canonical lineage panel
# ----------------------------------------------------------------------

canonical_panels <- list(
  Unciliated_epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "PAX8",
    "MUC1"
  ),
  Ciliated_epithelial = c(
    "FOXJ1",
    "PIFO",
    "TPPP3",
    "CAPS",
    "C20orf85",
    "RSPH1"
  ),
  Stromal = c(
    "COL1A1",
    "COL1A2",
    "DCN",
    "COL3A1",
    "PDGFRA",
    "COL6A1"
  ),
  Blood_endothelial = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "RAMP2",
    "PLVAP"
  ),
  Lymphatic_endothelial = c(
    "CCL21",
    "PROX1",
    "LYVE1",
    "FLT4",
    "PDPN",
    "FABP4"
  ),
  NK = c(
    "NKG7",
    "GNLY",
    "KLRD1",
    "NCAM1",
    "PRF1",
    "FGFBP2"
  ),
  T = c(
    "CD3D",
    "CD3E",
    "TRBC1",
    "TRBC2",
    "IL7R",
    "LTB"
  ),
  B = c(
    "CD79A",
    "MS4A1",
    "CD37",
    "CD74",
    "CD22",
    "HLA-DRA"
  ),
  Plasma = c(
    "MZB1",
    "JCHAIN",
    "IGKC",
    "IGLC2",
    "CD79A",
    "SEC11C"
  ),
  Myeloid = c(
    "LYZ",
    "LST1",
    "FCER1G",
    "TYROBP",
    "CTSS",
    "HLA-DRA"
  ),
  pDC = c(
    "GZMB",
    "JCHAIN",
    "IL3RA",
    "CLEC4C",
    "TCF4",
    "SERPINF1"
  ),
  Neutrophil = c(
    "S100A8",
    "S100A9",
    "FCGR3B",
    "CSF3R",
    "FPR1",
    "MNDA"
  ),
  Mast = c(
    "TPSAB1",
    "TPSB2",
    "CPA3",
    "KIT",
    "HDC",
    "MS4A2"
  ),
  Cycling = c(
    "MKI67",
    "TOP2A",
    "TUBA1B",
    "STMN1",
    "CENPF",
    "UBE2C"
  )
)

panel_requested <- data.frame(
  lineage = rep(
    names(canonical_panels),
    lengths(canonical_panels)
  ),
  requested_gene = unlist(
    canonical_panels,
    use.names = FALSE
  ),
  stringsAsFactors = FALSE
)

panel_resolution <- resolve_features(
  panel_requested$requested_gene,
  rownames(integrated_object)
)

if (
  nrow(panel_resolution) !=
    nrow(panel_requested) ||
    !identical(
      panel_resolution$request_index,
      seq_len(nrow(panel_requested))
    ) ||
    !identical(
      panel_resolution$requested_gene,
      panel_requested$requested_gene
    )
) {
  stop(
    "Canonical-panel feature resolution is not one-to-one with the 84 requested panel rows."
  )
}

panel_requested$request_index <-
  panel_resolution$request_index

panel_requested$resolved_feature <-
  panel_resolution$resolved_feature

panel_requested$available <-
  panel_resolution$available

panel_duplicate_audit <- data.frame(
  metric = c(
    "requested_panel_row_n",
    "unique_requested_gene_n",
    "duplicated_cross_lineage_row_n",
    "resolved_panel_row_n",
    "available_panel_row_n",
    "unavailable_panel_row_n",
    "one_to_one_assignment_pass"
  ),
  value = c(
    nrow(panel_requested),
    length(
      unique(
        panel_requested$requested_gene
      )
    ),
    sum(
      duplicated(
        panel_requested$requested_gene
      )
    ),
    nrow(panel_resolution),
    sum(panel_requested$available),
    sum(!panel_requested$available),
    nrow(panel_resolution) ==
      nrow(panel_requested)
  ),
  stringsAsFactors = FALSE
)

panel_duplicate_rows <- panel_requested[
  duplicated(
    panel_requested$requested_gene
  ) |
    duplicated(
      panel_requested$requested_gene,
      fromLast = TRUE
    ),
  ,
  drop = FALSE
]

write_csv_utf8(
  panel_duplicate_audit,
  file.path(
    panel_root,
    "00_panel_row_audit.csv"
  )
)

write_csv_utf8(
  panel_duplicate_rows,
  file.path(
    panel_root,
    "00_cross_lineage_duplicate_genes.csv"
  )
)

write_csv_utf8(
  panel_requested,
  file.path(
    panel_root,
    "canonical_panel_manifest.csv"
  )
)

available_panel <- panel_requested[
  panel_requested$available,
  ,
  drop = FALSE
]

if (
  nrow(available_panel) <
    0.75 * nrow(panel_requested)
) {
  stop(
    "Canonical panel coverage is below 75%."
  )
}

panel_features <- unique(
  available_panel$
    resolved_feature
)

panel_expression <- get_data_matrix(
  integrated_object,
  panel_features
)

panel_rows <- list()
panel_row_index <- 0L

for (cluster_id in cluster_ids) {
  cluster_cells <- which(
    as.character(
      metadata$seurat_clusters
    ) == cluster_id
  )

  cluster_average <- Matrix::rowMeans(
    panel_expression[
      ,
      cluster_cells,
      drop = FALSE
    ]
  )

  cluster_pct <- Matrix::rowMeans(
    panel_expression[
      ,
      cluster_cells,
      drop = FALSE
    ] > 0
  )

  cluster_panel <- data.frame(
    cluster = cluster_id,
    resolved_feature =
      rownames(panel_expression),
    average_expression =
      as.numeric(
        cluster_average
      ),
    percent_expressing =
      as.numeric(
        cluster_pct
      ),
    stringsAsFactors = FALSE
  )

  cluster_panel <- merge(
    available_panel,
    cluster_panel,
    by = "resolved_feature",
    all.x = TRUE,
    sort = FALSE
  )

  panel_row_index <- panel_row_index + 1L

  panel_rows[[panel_row_index]] <-
    cluster_panel
}

canonical_panel_long <- do.call(
  rbind,
  panel_rows
)

canonical_panel_long$
  expression_for_scaling <- log1p(
    canonical_panel_long$
      average_expression
  )

canonical_panel_long$gene_z <- ave(
  canonical_panel_long$
    expression_for_scaling,
  canonical_panel_long$
    resolved_feature,
  FUN = safe_zscore
)

canonical_panel_long$gene_z <- pmax(
  -3,
  pmin(
    3,
    canonical_panel_long$gene_z
  )
)

lineage_scores <- data.table::as.data.table(
  canonical_panel_long
)[
  ,
  .(
    available_gene_n =
      data.table::uniqueN(
        resolved_feature
      ),
    lineage_score =
      mean(
        gene_z,
        na.rm = TRUE
      ),
    mean_percent_expressing =
      mean(
        percent_expressing,
        na.rm = TRUE
      ),
    mean_expression =
      mean(
        average_expression,
        na.rm = TRUE
      )
  ),
  by = .(
    cluster,
    lineage
  )
]

lineage_scores <- as.data.frame(
  lineage_scores
)

lineage_scores <- lineage_scores[
  order(
    as.integer(
      lineage_scores$cluster
    ),
    -lineage_scores$
      lineage_score
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  canonical_panel_long,
  file.path(
    panel_root,
    "canonical_panel.csv"
  )
)

write_csv_utf8(
  lineage_scores,
  file.path(
    panel_root,
    "lineage_scores.csv"
  )
)


# ----------------------------------------------------------------------
# 8. Candidate broad-lineage nomination
# ----------------------------------------------------------------------

annotation_lineages <- setdiff(
  names(canonical_panels),
  "Cycling"
)

top50_marker_panel_overlap <- merge(
  all_markers[
    all_markers$rank <= 50,
    c(
      "cluster",
      "gene"
    ),
    drop = FALSE
  ],
  available_panel[
    ,
    c(
      "lineage",
      "resolved_feature"
    ),
    drop = FALSE
  ],
  by.x = "gene",
  by.y = "resolved_feature",
  all = FALSE,
  sort = FALSE
)

top50_overlap_counts <- if (
  nrow(top50_marker_panel_overlap)
) {
  as.data.frame(
    data.table::as.data.table(
      top50_marker_panel_overlap
    )[
      ,
      .(
        top50_marker_overlap_n =
          data.table::uniqueN(gene)
      ),
      by = .(
        cluster,
        lineage
      )
    ]
  )
} else {
  data.frame(
    cluster = character(),
    lineage = character(),
    top50_marker_overlap_n =
      integer(),
    stringsAsFactors = FALSE
  )
}

candidate_rows <- list()

for (cluster_id in cluster_ids) {
  cluster_scores <- lineage_scores[
    lineage_scores$cluster ==
      cluster_id &
      lineage_scores$lineage %in%
        annotation_lineages,
    ,
    drop = FALSE
  ]

  cluster_scores <- cluster_scores[
    order(
      -cluster_scores$
        lineage_score,
      -cluster_scores$
        mean_percent_expressing
    ),
    ,
    drop = FALSE
  ]

  top_lineage <- cluster_scores$
    lineage[1]

  second_lineage <- cluster_scores$
    lineage[2]

  top_score <- cluster_scores$
    lineage_score[1]

  second_score <- cluster_scores$
    lineage_score[2]

  score_margin <- top_score -
    second_score

  overlap_value <- top50_overlap_counts[
    top50_overlap_counts$cluster ==
      cluster_id &
      top50_overlap_counts$lineage ==
        top_lineage,
    "top50_marker_overlap_n"
  ]

  if (!length(overlap_value)) {
    overlap_value <- 0L
  }

  cycling_score <- lineage_scores$
    lineage_score[
      lineage_scores$cluster ==
        cluster_id &
        lineage_scores$lineage ==
          "Cycling"
    ]

  if (!length(cycling_score)) {
    cycling_score <- NA_real_
  }

  review_cluster <- cluster_id %in%
    c(
      "7",
      "8",
      "13",
      "15"
    )

  support_status <- if (
    score_margin >= 0.50 &&
      overlap_value >= 2
  ) {
    "HIGH_SUPPORT_CANDIDATE"
  } else if (
    score_margin >= 0.25 ||
      overlap_value >= 2
  ) {
    "MODERATE_SUPPORT_CANDIDATE"
  } else {
    "AMBIGUOUS_REVIEW_REQUIRED"
  }

  if (review_cluster) {
    support_status <- paste0(
      support_status,
      "_MANDATORY_REVIEW"
    )
  }

  candidate_rows[[
    length(candidate_rows) + 1L
  ]] <- data.frame(
    cluster = cluster_id,
    candidate_broad_lineage =
      top_lineage,
    second_lineage =
      second_lineage,
    top_lineage_score =
      top_score,
    second_lineage_score =
      second_score,
    score_margin =
      score_margin,
    top50_marker_overlap_n =
      as.integer(overlap_value),
    cycling_score =
      cycling_score,
    candidate_support_status =
      support_status,
    mandatory_review_cluster =
      review_cluster,
    annotation_locked = FALSE,
    stringsAsFactors = FALSE
  )
}

candidate_annotation <- do.call(
  rbind,
  candidate_rows
)

candidate_annotation <- candidate_annotation[
  order(
    as.integer(
      candidate_annotation$cluster
    )
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  top50_overlap_counts,
  file.path(
    panel_root,
    "top50_panel_overlap.csv"
  )
)

write_csv_utf8(
  candidate_annotation,
  file.path(
    panel_root,
    "candidate_annotation.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Marker summaries and QC figures
# ----------------------------------------------------------------------

marker_summary <- data.table::as.data.table(
  all_markers
)[
  ,
  .(
    marker_row_n = .N,
    FDR05_marker_n =
      sum(
        p_val_adj < 0.05,
        na.rm = TRUE
      ),
    logFC05_FDR05_marker_n =
      sum(
        p_val_adj < 0.05 &
          avg_log2FC >= 0.50,
        na.rm = TRUE
      ),
    robust_donor_supported_n =
      sum(
        robust_donor_supported,
        na.rm = TRUE
      ),
    top_genes =
      paste(
        head(gene, 10),
        collapse = "; "
      ),
    top_robust_genes =
      paste(
        head(
          gene[
            robust_donor_supported
          ],
          10
        ),
        collapse = "; "
      )
  ),
  by = cluster
]

marker_summary <- as.data.frame(
  marker_summary
)

marker_summary <- merge(
  marker_summary,
  candidate_annotation,
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

marker_summary <- marker_summary[
  order(
    as.integer(
      marker_summary$cluster
    )
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  marker_summary,
  file.path(
    output_root,
    "02_marker_summary.csv"
  )
)

heatmap_gene_rows <- list()

for (cluster_id in cluster_ids) {
  cluster_markers <- all_markers[
    all_markers$cluster ==
      cluster_id,
    ,
    drop = FALSE
  ]

  robust_genes <- cluster_markers$gene[
    cluster_markers$
      robust_donor_supported
  ]

  selected_genes <- unique(
    c(
      head(
        robust_genes,
        5
      ),
      head(
        cluster_markers$gene,
        5
      )
    )
  )

  selected_genes <- head(
    selected_genes,
    5
  )

  heatmap_gene_rows[[
    length(heatmap_gene_rows) + 1L
  ]] <- data.frame(
    source_cluster =
      cluster_id,
    gene = selected_genes,
    stringsAsFactors = FALSE
  )
}

heatmap_gene_manifest <- unique(
  do.call(
    rbind,
    heatmap_gene_rows
  )
)

heatmap_genes <- unique(
  heatmap_gene_manifest$gene
)

heatmap_expression <- get_data_matrix(
  integrated_object,
  heatmap_genes
)

heatmap_rows <- list()
heatmap_row_index <- 0L

for (cluster_id in cluster_ids) {
  cluster_cells <- which(
    as.character(
      metadata$seurat_clusters
    ) == cluster_id
  )

  average_values <- Matrix::rowMeans(
    heatmap_expression[
      ,
      cluster_cells,
      drop = FALSE
    ]
  )

  heatmap_row_index <-
    heatmap_row_index + 1L

  heatmap_rows[[heatmap_row_index]] <-
    data.frame(
    cluster = cluster_id,
    gene =
      rownames(
        heatmap_expression
      ),
    average_expression =
      as.numeric(
        average_values
      ),
    stringsAsFactors = FALSE
  )
}

heatmap_data <- do.call(
  rbind,
  heatmap_rows
)

heatmap_data$gene_z <- ave(
  log1p(
    heatmap_data$
      average_expression
  ),
  heatmap_data$gene,
  FUN = safe_zscore
)

heatmap_data$gene_z <- pmax(
  -3,
  pmin(
    3,
    heatmap_data$gene_z
  )
)

write_csv_utf8(
  heatmap_gene_manifest,
  file.path(
    panel_root,
    "heatmap_gene_manifest.csv"
  )
)

write_csv_utf8(
  heatmap_data,
  file.path(
    panel_root,
    "heatmap_source_data.csv"
  )
)

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

marker_summary_plot <- marker_summary
marker_summary_plot$cluster <- factor(
  marker_summary_plot$cluster,
  levels = cluster_ids
)

plot_marker_counts <- ggplot2::ggplot(
  marker_summary_plot,
  ggplot2::aes(
    x = cluster,
    y = FDR05_marker_n,
    shape =
      mandatory_review_cluster
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "GSE183837 descriptive positive-marker counts",
    subtitle = "Wilcoxon cluster-versus-rest; descriptive only",
    x = "Cluster",
    y = "FDR < 0.05 marker rows",
    shape = "Mandatory review"
  ) +
  theme_QC

save_plot_pair(
  plot_marker_counts,
  "QC1_marker_counts",
  width = 8,
  height = 5
)

plot_robust_counts <- ggplot2::ggplot(
  marker_summary_plot,
  ggplot2::aes(
    x = cluster,
    y = robust_donor_supported_n,
    shape =
      mandatory_review_cluster
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Donor-supported descriptive markers",
    subtitle = "Requires FDR, effect size, detection and support in at least half of eligible donors",
    x = "Cluster",
    y = "Donor-supported marker count",
    shape = "Mandatory review"
  ) +
  theme_QC

save_plot_pair(
  plot_robust_counts,
  "QC2_donor_supported_marker_counts",
  width = 8,
  height = 5
)

heatmap_data$cluster <- factor(
  heatmap_data$cluster,
  levels = cluster_ids
)

heatmap_data$gene <- factor(
  heatmap_data$gene,
  levels = rev(
    unique(
      heatmap_gene_manifest$gene
    )
  )
)

plot_marker_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = cluster,
    y = gene,
    fill = gene_z
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    option = "C"
  ) +
  ggplot2::labs(
    title = "Top descriptive marker expression by cluster",
    subtitle = "Row-scaled log-normalized average expression",
    x = "Cluster",
    y = "Selected marker",
    fill = "Row z-score"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_text(
        size = 6
      )
  )

save_plot_pair(
  plot_marker_heatmap,
  "QC3_top_marker_heatmap",
  width = 10,
  height = 14
)

lineage_plot_data <- lineage_scores
lineage_plot_data$cluster <- factor(
  lineage_plot_data$cluster,
  levels = cluster_ids
)

lineage_plot_data$lineage <- factor(
  lineage_plot_data$lineage,
  levels = rev(
    names(canonical_panels)
  )
)

plot_lineage <- ggplot2::ggplot(
  lineage_plot_data,
  ggplot2::aes(
    x = cluster,
    y = lineage,
    fill = lineage_score
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    option = "C"
  ) +
  ggplot2::labs(
    title = "Canonical lineage-panel scores",
    subtitle = "Candidate broad identities require marker-level review",
    x = "Cluster",
    y = "Canonical panel",
    fill = "Mean gene z-score"
  ) +
  theme_QC

save_plot_pair(
  plot_lineage,
  "QC4_lineage_panel_heatmap",
  width = 10,
  height = 7
)

candidate_plot <- candidate_annotation
candidate_plot$cluster <- factor(
  candidate_plot$cluster,
  levels = cluster_ids
)

plot_margin <- ggplot2::ggplot(
  candidate_plot,
  ggplot2::aes(
    x = cluster,
    y = score_margin,
    shape =
      mandatory_review_cluster
  )
) +
  ggplot2::geom_hline(
    yintercept = c(
      0.25,
      0.50
    ),
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::labs(
    title = "Candidate broad-lineage score margins",
    subtitle = "Low margins indicate ambiguous lineage nomination",
    x = "Cluster",
    y = "Top minus second lineage score",
    shape = "Mandatory review"
  ) +
  theme_QC

save_plot_pair(
  plot_margin,
  "QC5_candidate_lineage_margin",
  width = 8,
  height = 5
)

review_support <- donor_support_summary_all[
  donor_support_summary_all$cluster %in%
    c(
      "7",
      "8",
      "13",
      "15"
    ),
  ,
  drop = FALSE
]

review_support <- merge(
  review_support,
  all_markers[
    ,
    c(
      "cluster",
      "gene",
      "rank",
      "avg_log2FC",
      "p_val_adj"
    ),
    drop = FALSE
  ],
  by = c(
    "cluster",
    "gene"
  ),
  all.x = TRUE,
  sort = FALSE
)

review_support <- review_support[
  review_support$rank <= 30,
  ,
  drop = FALSE
]

review_support$cluster <- factor(
  review_support$cluster,
  levels = c(
    "7",
    "8",
    "13",
    "15"
  )
)

write_csv_utf8(
  review_support,
  file.path(
    panel_root,
    "review_cluster_marker_support.csv"
  )
)

plot_review_support <- ggplot2::ggplot(
  review_support,
  ggplot2::aes(
    x = supporting_donor_fraction,
    y = avg_log2FC,
    shape = cluster
  )
) +
  ggplot2::geom_vline(
    xintercept = 0.50,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    alpha = 0.75
  ) +
  ggplot2::scale_x_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Marker donor support in mandatory-review clusters",
    subtitle = "Cluster 13 requires special caution because of severe donor dominance",
    x = "Supporting fraction of eligible donors",
    y = "Average log2 fold change",
    shape = "Cluster"
  ) +
  theme_QC

save_plot_pair(
  plot_review_support,
  "QC6_review_cluster_donor_support",
  width = 8,
  height = 6
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
# 10. Marker Gate
# ----------------------------------------------------------------------

ambiguous_clusters <- candidate_annotation$
  cluster[
    grepl(
      "^AMBIGUOUS",
      candidate_annotation$
        candidate_support_status
    )
  ]

review_clusters <- c(
  "7",
  "8",
  "13",
  "15"
)

gate_checks <- data.frame(
  check = c(
    "upstream_detailed_QC_authorizes_markers",
    "input_checkpoint_has_47161_cells",
    "all_16_cluster_marker_jobs_complete",
    "every_cluster_has_positive_marker_rows",
    "combined_marker_table_has_all_16_clusters",
    "all_top100_table_has_all_16_clusters",
    "donor_support_completed_for_all_16_clusters",
    "canonical_panel_coverage_at_least_75_percent",
    "lineage_scores_present_for_all_clusters",
    "candidate_annotation_has_16_clusters",
    "mandatory_review_clusters_7_8_13_15_present",
    "no_cluster_or_cell_was_removed",
    "candidate_annotation_review_can_start",
    "free_disk_at_least_8GB_or_unknown"
  ),
  pass = c(
    upstream_gate$gate_status[[1]] %in% c(
      "PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED",
      "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED"
    ),
    ncol(integrated_object) == 47161L,
    nrow(cluster_status_all) == 16L &&
      all(
        cluster_status_all$status ==
          "COMPLETE"
      ),
    all(
      cluster_status_all$
        marker_row_n > 0
    ),
    length(
      unique(
        all_markers$cluster
      )
    ) == 16L,
    length(
      unique(
        all_top100$cluster
      )
    ) == 16L,
    length(
      unique(
        donor_support_summary_all$
          cluster
      )
    ) == 16L,
    nrow(available_panel) >=
      0.75 * nrow(panel_requested),
    length(
      unique(
        lineage_scores$cluster
      )
    ) == 16L,
    nrow(candidate_annotation) ==
      16L,
    all(
      review_clusters %in%
        candidate_annotation$cluster
    ),
    ncol(integrated_object) ==
      47161L &&
      length(
        unique(
          metadata$seurat_clusters
        )
      ) == 16L,
    TRUE,
    is.na(
      suppressWarnings(
        as.numeric(
          system2(
            "powershell",
            args = c(
              "-NoProfile",
              "-Command",
              shQuote(
                "(Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='C:'\").FreeSpace / 1GB"
              )
            ),
            stdout = TRUE,
            stderr = FALSE
          )[1]
        )
      )
    ) ||
      suppressWarnings(
        as.numeric(
          system2(
            "powershell",
            args = c(
              "-NoProfile",
              "-Command",
              shQuote(
                "(Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='C:'\").FreeSpace / 1GB"
              )
            ),
            stdout = TRUE,
            stderr = FALSE
          )[1]
        )
      ) >= 8
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
  "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED"
} else {
  "HOLD_DESCRIPTIVE_MARKERS_INCOMPLETE"
}

review_flags <- c(
  "CLUSTER_13_SEVERE_SINGLE_DONOR_DOMINANCE",
  "CLUSTERS_7_8_MODERATE_DONOR_ENRICHMENT",
  "CLUSTER_15_RARE_POPULATION_REVIEW",
  if (length(ambiguous_clusters)) {
    paste0(
      "AMBIGUOUS_LINEAGE_CLUSTERS_",
      paste(
        ambiguous_clusters,
        collapse = "_"
      )
    )
  } else {
    character()
  }
)

marker_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version =
    marker_version,
  gate_status = gate_status,
  cell_n = ncol(integrated_object),
  cluster_n = 16L,
  total_marker_row_n =
    nrow(all_markers),
  FDR05_marker_row_n =
    sum(
      all_markers$p_val_adj < 0.05,
      na.rm = TRUE
    ),
  robust_donor_supported_marker_n =
    sum(
      all_markers$
        robust_donor_supported,
      na.rm = TRUE
    ),
  canonical_panel_available_n =
    nrow(available_panel),
  canonical_panel_requested_n =
    nrow(panel_requested),
  ambiguous_candidate_clusters =
    collapse_unique(
      ambiguous_clusters
    ),
  review_flags =
    collapse_unique(
      review_flags
    ),
  candidate_annotation_review_authorized =
    all_gate_checks_pass,
  broad_annotation_locked =
    FALSE,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  marker_gate,
  file.path(
    run_qc_root,
    "G183_MARKER_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1L GSE183837 Descriptive Marker Card",
  "",
  paste0("- Gate: ", gate_status),
  paste0(
    "- Cells retained: ",
    format(
      ncol(integrated_object),
      big.mark = ","
    )
  ),
  "- Primary clusters: 16",
  paste0(
    "- Total marker rows: ",
    format(
      nrow(all_markers),
      big.mark = ","
    )
  ),
  paste0(
    "- FDR < 0.05 marker rows: ",
    format(
      sum(
        all_markers$
          p_val_adj < 0.05,
        na.rm = TRUE
      ),
      big.mark = ","
    )
  ),
  paste0(
    "- Donor-supported marker rows: ",
    format(
      sum(
        all_markers$
          robust_donor_supported,
        na.rm = TRUE
      ),
      big.mark = ","
    )
  ),
  paste0(
    "- Ambiguous candidate clusters: ",
    if (
      length(ambiguous_clusters)
    ) {
      collapse_unique(
        ambiguous_clusters
      )
    } else {
      "NONE"
    }
  ),
  paste0(
    "- Candidate annotation review authorized: ",
    all_gate_checks_pass
  ),
  "- Broad annotation locked: FALSE",
  "- Disease pseudobulk testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "G183_MARKER_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 descriptive marker interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Interpretation constraints",
  "",
  "- Marker tests are cluster-versus-rest descriptions and are not Control-versus-RIF comparisons.",
  "- Candidate broad lineages are algorithmic nominations requiring manual review.",
  "- Cluster 13 cannot be described as an RIF-specific population because it is dominated by one RIF donor.",
  "- Clusters 7 and 8 require donor-robust marker confirmation.",
  "- Cluster 15 remains a retained rare candidate population.",
  "",
  "## Authorized next operation",
  "",
  "- Review top markers, canonical panels and donor support.",
  "- Construct and QC a broad-lineage annotation map.",
  "- Lock broad annotation only after marker and visual review.",
  "",
  "## Not authorized yet",
  "",
  "- No disease pseudobulk analysis.",
  "- No fine-subtype lock.",
  "- No result lock."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "G183_MARKER_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 11. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_execution =
    "V1_0_CANONICAL_PANEL_UNIQUE_COLLAPSED_84_ROWS_TO_81",
  panel_resolution_policy =
    "PRESERVE_CROSS_LINEAGE_DUPLICATE_GENES_AND_INPUT_ROW_ORDER",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version = marker_version,
  upstream_gate =
    upstream_gate$gate_status[[1]],
  marker_parameters =
    marker_parameters,
  marker_role =
    "DESCRIPTIVE_CLUSTER_VERSUS_REST",
  donor_support_rule =
    "PCT_DELTA_AT_LEAST_0_10_IN_AT_LEAST_HALF_OF_ELIGIBLE_DONORS",
  canonical_panel_role =
    "CANDIDATE_BROAD_LINEAGE_NOMINATION_ONLY",
  mandatory_review_clusters =
    c(
      "7",
      "8",
      "13",
      "15"
    ),
  broad_annotation_locked =
    FALSE,
  disease_cell_level_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_execution =
    "V1_0_CANONICAL_PANEL_UNIQUE_COLLAPSED_84_ROWS_TO_81",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version =
    marker_version,
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
  cell_n = ncol(integrated_object),
  cluster_n = 16L,
  total_marker_row_n =
    nrow(all_markers),
  FDR05_marker_row_n =
    sum(
      all_markers$p_val_adj < 0.05,
      na.rm = TRUE
    ),
  robust_donor_supported_marker_n =
    sum(
      all_markers$
        robust_donor_supported,
      na.rm = TRUE
    ),
  candidate_annotation_review_authorized =
    all_gate_checks_pass,
  broad_annotation_locked =
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

persistent_files <- unique(c(
  list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE
  )
))

persistent_files <- persistent_files[
  file.exists(persistent_files) &
    !dir.exists(persistent_files)
]

persistent_manifest <- data.frame(
  path = normalize_slashes(
    persistent_files
  ),
  size_bytes =
    file.info(persistent_files)$size,
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
  file.path(
    output_root,
    "01_marker_status.csv"
  ),
  file.path(
    output_root,
    "02_marker_summary.csv"
  ),
  list.files(
    top_root,
    pattern = "top100|donor_summary|donor_support_summary",
    full.names = TRUE
  ),
  list.files(
    panel_root,
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    figure_output_root,
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
    file.info(small_files)$size,
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
  "total_marker_rows=",
  nrow(all_markers),
  "robust_donor_supported=",
  sum(
    all_markers$
      robust_donor_supported,
    na.rm = TRUE
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
cat("STEP04E1L GSE183837 descriptive markers completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Cells: ",
  format(
    ncol(integrated_object),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat("Clusters completed: 16\n")
cat(
  "Total marker rows: ",
  format(
    nrow(all_markers),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "FDR < 0.05 marker rows: ",
  format(
    sum(
      all_markers$p_val_adj < 0.05,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Donor-supported marker rows: ",
  format(
    sum(
      all_markers$
        robust_donor_supported,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Ambiguous candidate clusters: ",
  if (
    length(ambiguous_clusters)
  ) {
    collapse_unique(
      ambiguous_clusters
    )
  } else {
    "NONE"
  },
  "\n",
  sep = ""
)
cat(
  "Candidate annotation review authorized: ",
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
    "\nMarker outputs are complete in the project folders.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 5: 67_step04E1M_GSE183837_broad_annotation_QC_v1
# ==============================================================================

# ======================================================================
# 67_step04E1M_GSE183837_broad_annotation_QC_v1.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Apply a manually reviewed, marker-supported broad annotation map to
#   the 16 GSE183837 primary clusters.
# - Preserve the ontology used for GSE250130 where possible:
#   Stromal cells, NK/T cells, Unciliated epithelial cells,
#   Ciliated epithelial cells, Myeloid cells, Endothelial cells.
# - Lock broad annotation while retaining cluster-level review flags.
# - Keep cluster 8 in the atlas as a low-confidence epithelial-like
#   population but exclude it from primary broad pseudobulk models.
# - Authorize VDR localization and pseudobulk preparation only after QC.
#
# UPSTREAM
# - PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED
# - 47,161 cells, 9 donors, 16 clusters.
#
# SCIENTIFIC GUARDRAILS
# - Broad annotation is not a disease result.
# - Cluster 13 is cycling stromal but strongly donor-dominant; it cannot
#   support an RIF-specific subtype claim.
# - Cluster 8 is low-confidence and sensitivity-only for downstream
#   broad-cell models.
# - No Control-versus-RIF test is performed here.
# - Original integrated.rds remains immutable.
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
dataset_id <- "GSE183837"

script_name <-
  "67_step04E1M_GSE183837_broad_annotation_QC_v1.R"

script_version <-
  "STEP04E1M_GSE183837_BROAD_ANNOTATION_QC_v1.0"

integration_version <- "SCI183_v1_0"
marker_version <- "MK183_v1_0"
annotation_version <- "ANN183_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1M_",
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

assert_path_safe <- function(
  path,
  max_chars = 235L
) {
  path_n <- nchar(
    normalize_slashes(path),
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
      normalize_slashes(path)
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

collapse_unique <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]

  if (!length(x)) {
    ""
  } else {
    paste(x, collapse = sep)
  }
}

split_semicolon <- function(x) {
  output <- unlist(
    strsplit(
      as.character(x),
      ";",
      fixed = TRUE
    ),
    use.names = FALSE
  )

  trimws(
    output[
      nzchar(trimws(output))
    ]
  )
}

safe_zscore <- function(x) {
  x <- as.numeric(x)

  if (
    length(x) < 2L ||
      !is.finite(stats::sd(x)) ||
      stats::sd(x) == 0
  ) {
    return(rep(0, length(x)))
  }

  as.numeric(scale(x))
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
  "SeuratObject",
  "Matrix",
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

parameter_root <- file.path(
  step04_root,
  "03_parameters"
)

environment_root <- file.path(
  step04_root,
  "04_environment"
)

integration_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

input_checkpoint <- file.path(
  integration_root,
  "checkpoints",
  dataset_id,
  "integrated.rds"
)

marker_root <- file.path(
  integration_root,
  "mk",
  "G183",
  marker_version
)

marker_summary_path <- file.path(
  marker_root,
  "02_marker_summary.csv"
)

all_top100_path <- file.path(
  marker_root,
  "top",
  "all_top100.csv"
)

donor_support_path <- file.path(
  marker_root,
  "top",
  "donor_support_summary.csv"
)

candidate_annotation_path <- file.path(
  marker_root,
  "panel",
  "candidate_annotation.csv"
)

detailed_qc_table_path <- file.path(
  integration_root,
  "iqc",
  "G183",
  "IQ183_v1_2",
  "t",
  "04_cluster_QC.csv"
)

output_root <- file.path(
  integration_root,
  "ann",
  "G183",
  annotation_version
)

table_root <- file.path(
  output_root,
  "t"
)

figure_output_root <- file.path(
  output_root,
  "fig"
)

manifest_root <- file.path(
  output_root,
  "man"
)

annotated_checkpoint <- file.path(
  output_root,
  "annotated.rds"
)

annotation_status_path <- file.path(
  output_root,
  "status.csv"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1M_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_annQC"
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
  output_root,
  table_root,
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
  artifact = c(
    "input_checkpoint",
    "marker_summary",
    "all_top100",
    "donor_support",
    "candidate_annotation",
    "detailed_cluster_QC",
    "output_root",
    "annotated_checkpoint",
    "annotation_status",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    input_checkpoint,
    marker_summary_path,
    all_top100_path,
    donor_support_path,
    candidate_annotation_path,
    detailed_qc_table_path,
    output_root,
    annotated_checkpoint,
    annotation_status_path,
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
    "One or more annotation paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream marker Gate
# ----------------------------------------------------------------------

marker_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

marker_qc_dirs <- marker_qc_dirs[
  grepl(
    "^STEP04E1L_[0-9]{8}_[0-9]{6}_markers$",
    basename(marker_qc_dirs),
    perl = TRUE
  )
]

if (!length(marker_qc_dirs)) {
  stop(
    "No STEP04E1L marker-QC directory was found."
  )
}

gate_candidates <- list()

for (directory_path in marker_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "G183_MARKER_GATE.csv"
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

    gate_candidates[[
      length(gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(gate_candidates)) {
  stop(
    "No readable STEP04E1L marker Gate was found."
  )
}

gate_candidates <- as.data.frame(
  data.table::rbindlist(
    gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    gate_candidates$source_mtime
  )

gate_candidates <- gate_candidates[
  order(
    gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

authorized_index <- which(
  gate_candidates$gate_status ==
    "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED" &
    toupper(
      as.character(
        gate_candidates$
          candidate_annotation_review_authorized
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      )
)

if (!length(authorized_index)) {
  stop(
    "No marker Gate authorizes annotation review."
  )
}

upstream_gate <- gate_candidates[
  authorized_index[1],
  ,
  drop = FALSE
]

write_csv_utf8(
  upstream_gate,
  file.path(
    run_qc_root,
    "00_upstream_marker_gate.csv"
  )
)


# ----------------------------------------------------------------------
# 5. Load required inputs
# ----------------------------------------------------------------------

required_input_paths <- c(
  input_checkpoint,
  marker_summary_path,
  all_top100_path,
  donor_support_path,
  candidate_annotation_path,
  detailed_qc_table_path
)

missing_input_paths <- required_input_paths[
  !file.exists(required_input_paths)
]

if (length(missing_input_paths)) {
  stop(
    "Missing required input file(s):\n",
    paste(
      normalize_slashes(
        missing_input_paths
      ),
      collapse = "\n"
    )
  )
}

marker_summary <- safe_read_csv(
  marker_summary_path
)

all_top100 <- safe_read_csv(
  all_top100_path
)

donor_support <- safe_read_csv(
  donor_support_path
)

algorithmic_candidates <- safe_read_csv(
  candidate_annotation_path
)

detailed_cluster_QC <- safe_read_csv(
  detailed_qc_table_path
)

integrated_object <- readRDS(
  input_checkpoint
)

if (!inherits(integrated_object, "Seurat")) {
  stop(
    "Input checkpoint is not a Seurat object."
  )
}

metadata <- integrated_object[[]]

if (
  ncol(integrated_object) != 47161L ||
    length(
      unique(
        as.character(
          metadata$donor_id
        )
      )
    ) != 9L ||
    length(
      unique(
        as.character(
          metadata$seurat_clusters
        )
      )
    ) != 16L
) {
  stop(
    "Integrated object does not match the locked 47,161-cell object."
  )
}

if (
  nrow(marker_summary) != 16L ||
    nrow(algorithmic_candidates) != 16L ||
    length(
      unique(
        all_top100$cluster
      )
    ) != 16L
) {
  stop(
    "Marker outputs do not contain all 16 clusters."
  )
}

log_event(
  "Starting",
  run_id,
  "cells=",
  ncol(integrated_object),
  "clusters=",
  length(
    unique(
      metadata$seurat_clusters
    )
  )
)


# ----------------------------------------------------------------------
# 6. Manually reviewed broad annotation map
# ----------------------------------------------------------------------

annotation_map <- data.frame(
  cluster = as.character(0:15),
  broad_cell_type = c(
    "Stromal cells",
    "Stromal cells",
    "Stromal cells",
    "Stromal cells",
    "Stromal cells",
    "Stromal cells",
    "NK/T cells",
    "Unciliated epithelial cells",
    "Unciliated epithelial cells",
    "NK/T cells",
    "Stromal cells",
    "Myeloid cells",
    "Endothelial cells",
    "Stromal cells",
    "NK/T cells",
    "Ciliated epithelial cells"
  ),
  subtype_label = c(
    "RARRES2_RAMP1_ribosomal_stromal",
    "GDF7_COL6A1_stromal",
    "CILP_APOD_PTGDS_stromal",
    "HOXA13_ITGA8_ECM_stromal",
    "MGP_PAMR1_MMP11_stromal",
    "KIF20A_TOP2A_cycling_stromal",
    "IL7R_CD3D_T",
    "PAX2_GABRP_SCGB_unciliated_epithelial",
    "SCGB1D2_low_complexity_epithelial_like",
    "GNLY_TRDC_KLRC1_NK_T",
    "RGS5_MYH11_NOTCH3_perivascular_stromal",
    "C1QC_LYZ_FOLR2_myeloid",
    "ACKR1_CLDN5_EMCN_blood_endothelial",
    "TOP2A_MKI67_MMP11_cycling_stromal",
    "FGFBP2_FCGR3A_cytotoxic_NK",
    "C20orf85_SNTN_ciliated_epithelial"
  ),
  annotation_confidence = c(
    "MODERATE",
    "MODERATE",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "LOW",
    "HIGH",
    "HIGH",
    "HIGH",
    "HIGH",
    "MODERATE",
    "HIGH",
    "HIGH"
  ),
  marker_evidence = c(
    "RARRES2;RAMP1;LGALS1",
    "GDF7;COL6A1;COL5A1",
    "CILP;APOD;PTGDS",
    "HOXA13;ITGA8;COL6A3",
    "MGP;PAMR1;MMP11",
    "KIF20A;TOP2A;MKI67",
    "IL7R;CD3D;TRAC",
    "PAX2;GABRP;KRT8;KRT18",
    "SCGB1D2;SCGB2A1",
    "GNLY;TRDC;KLRC1;NKG7",
    "RGS5;MYH11;NOTCH3;MCAM",
    "C1QC;LYZ;FOLR2;CD14",
    "ACKR1;CLDN5;EMCN;VWF",
    "TOP2A;MKI67;MMP11;COL1A1",
    "FGFBP2;FCGR3A;PRF1;NKG7",
    "C20orf85;SNTN;CFAP126;FOXJ1"
  ),
  review_flag = c(
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "MODERATE_DONOR_ENRICHMENT_IDENTITY_ROBUST",
    "LOW_COMPLEXITY_NO_ROBUST_DONOR_MARKERS",
    "",
    "",
    "",
    "",
    "SEVERE_SINGLE_DONOR_DOMINANCE_NO_DISEASE_SUBTYPE_CLAIM",
    "",
    "RARE_CLUSTER_DISTRIBUTED_ACROSS_DONORS"
  ),
  primary_VDR_localization_eligible = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  primary_broad_pseudobulk_eligible = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  ),
  subtype_disease_inference_eligible = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    FALSE
  ),
  annotation_locked = TRUE,
  stringsAsFactors = FALSE
)

if (
  nrow(annotation_map) != 16L ||
    anyDuplicated(
      annotation_map$cluster
    ) > 0L ||
    !identical(
      annotation_map$cluster,
      as.character(0:15)
    )
) {
  stop(
    "Annotation map must contain one ordered row for every cluster 0-15."
  )
}

allowed_broad_types <- c(
  "Stromal cells",
  "NK/T cells",
  "Unciliated epithelial cells",
  "Ciliated epithelial cells",
  "Myeloid cells",
  "Endothelial cells"
)

if (
  !all(
    annotation_map$broad_cell_type %in%
      allowed_broad_types
  )
) {
  stop(
    "Annotation map contains an unexpected broad-cell ontology."
  )
}

write_csv_utf8(
  annotation_map,
  file.path(
    table_root,
    "01_annotation_map.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Marker-evidence audit
# ----------------------------------------------------------------------

evidence_rows <- list()

for (row_index in seq_len(nrow(annotation_map))) {
  cluster_id <- annotation_map$
    cluster[row_index]

  evidence_genes <- split_semicolon(
    annotation_map$
      marker_evidence[row_index]
  )

  cluster_top <- all_top100[
    as.character(
      all_top100$cluster
    ) == cluster_id,
    ,
    drop = FALSE
  ]

  cluster_support <- donor_support[
    as.character(
      donor_support$cluster
    ) == cluster_id,
    ,
    drop = FALSE
  ]

  for (gene_value in evidence_genes) {
    top_match <- cluster_top[
      toupper(
        cluster_top$gene
      ) == toupper(gene_value),
      ,
      drop = FALSE
    ]

    support_match <- cluster_support[
      toupper(
        cluster_support$gene
      ) == toupper(gene_value),
      ,
      drop = FALSE
    ]

    evidence_rows[[
      length(evidence_rows) + 1L
    ]] <- data.frame(
      cluster = cluster_id,
      broad_cell_type =
        annotation_map$
          broad_cell_type[row_index],
      subtype_label =
        annotation_map$
          subtype_label[row_index],
      annotation_confidence =
        annotation_map$
          annotation_confidence[row_index],
      requested_marker =
        gene_value,
      present_in_top100 =
        nrow(top_match) > 0L,
      top100_rank =
        if (nrow(top_match)) {
          min(top_match$rank)
        } else {
          NA_integer_
        },
      avg_log2FC =
        if (nrow(top_match)) {
          top_match$avg_log2FC[1]
        } else {
          NA_real_
        },
      p_val_adj =
        if (nrow(top_match)) {
          top_match$p_val_adj[1]
        } else {
          NA_real_
        },
      robust_donor_supported =
        if (
          nrow(top_match) &&
            "robust_donor_supported" %in%
              colnames(top_match)
        ) {
          isTRUE(
            as.logical(
              top_match$
                robust_donor_supported[1]
            )
          )
        } else {
          FALSE
        },
      supporting_donor_n =
        if (nrow(support_match)) {
          support_match$
            supporting_donor_n[1]
        } else {
          NA_real_
        },
      supporting_donor_fraction =
        if (nrow(support_match)) {
          support_match$
            supporting_donor_fraction[1]
        } else {
          NA_real_
        },
      stringsAsFactors = FALSE
    )
  }
}

marker_evidence_audit <- do.call(
  rbind,
  evidence_rows
)

evidence_summary <- data.table::as.data.table(
  marker_evidence_audit
)[
  ,
  .(
    requested_marker_n = .N,
    top100_evidence_n =
      sum(
        present_in_top100,
        na.rm = TRUE
      ),
    robust_evidence_n =
      sum(
        robust_donor_supported,
        na.rm = TRUE
      ),
    median_supporting_donor_fraction =
      median(
        supporting_donor_fraction,
        na.rm = TRUE
      )
  ),
  by = .(
    cluster,
    broad_cell_type,
    subtype_label,
    annotation_confidence
  )
]

evidence_summary <- as.data.frame(
  evidence_summary
)

evidence_summary$
  median_supporting_donor_fraction[
    !is.finite(
      evidence_summary$
        median_supporting_donor_fraction
    )
  ] <- NA_real_

annotation_QC <- merge(
  annotation_map,
  evidence_summary,
  by = c(
    "cluster",
    "broad_cell_type",
    "subtype_label",
    "annotation_confidence"
  ),
  all.x = TRUE,
  sort = FALSE
)

annotation_QC <- merge(
  annotation_QC,
  marker_summary[
    ,
    c(
      "cluster",
      "marker_row_n",
      "FDR05_marker_n",
      "robust_donor_supported_n",
      "top_genes",
      "top_robust_genes"
    ),
    drop = FALSE
  ],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

annotation_QC <- merge(
  annotation_QC,
  algorithmic_candidates[
    ,
    c(
      "cluster",
      "candidate_broad_lineage",
      "second_lineage",
      "score_margin",
      "candidate_support_status"
    ),
    drop = FALSE
  ],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

annotation_QC <- merge(
  annotation_QC,
  detailed_cluster_QC[
    ,
    c(
      "seurat_clusters",
      "cell_n",
      "represented_donor_n",
      "top_donor",
      "top_donor_fraction",
      "donor_entropy_normalized",
      "median_nFeature_RNA",
      "median_percent_mt",
      "median_scDblFinder_score"
    ),
    drop = FALSE
  ],
  by.x = "cluster",
  by.y = "seurat_clusters",
  all.x = TRUE,
  sort = FALSE
)

annotation_QC <- annotation_QC[
  match(
    as.character(0:15),
    annotation_QC$cluster
  ),
  ,
  drop = FALSE
]

annotation_QC$
  evidence_pass <- ifelse(
    annotation_QC$
      annotation_confidence ==
      "LOW",
    annotation_QC$
      top100_evidence_n >= 1L,
    annotation_QC$
      top100_evidence_n >= 2L
  )

write_csv_utf8(
  marker_evidence_audit,
  file.path(
    table_root,
    "02_marker_evidence_audit.csv"
  )
)

write_csv_utf8(
  annotation_QC,
  file.path(
    table_root,
    "03_annotation_QC.csv"
  )
)


# ----------------------------------------------------------------------
# 8. Apply annotations without modifying the original checkpoint
# ----------------------------------------------------------------------

metadata$primary_cluster <- as.character(
  metadata$seurat_clusters
)

map_index <- match(
  metadata$primary_cluster,
  annotation_map$cluster
)

if (anyNA(map_index)) {
  stop(
    "One or more cells could not be matched to the annotation map."
  )
}

metadata$broad_cell_type <-
  annotation_map$
    broad_cell_type[map_index]

metadata$subtype_label <-
  annotation_map$
    subtype_label[map_index]

metadata$annotation_confidence <-
  annotation_map$
    annotation_confidence[map_index]

metadata$annotation_review_flag <-
  annotation_map$
    review_flag[map_index]

metadata$
  primary_VDR_localization_eligible <-
  annotation_map$
    primary_VDR_localization_eligible[
      map_index
    ]

metadata$
  primary_broad_pseudobulk_eligible <-
  annotation_map$
    primary_broad_pseudobulk_eligible[
      map_index
    ]

metadata$
  subtype_disease_inference_eligible <-
  annotation_map$
    subtype_disease_inference_eligible[
      map_index
    ]

metadata$annotation_version <-
  annotation_version

if (
  anyNA(metadata$broad_cell_type) ||
    anyNA(metadata$subtype_label) ||
    nrow(metadata) != 47161L
) {
  stop(
    "Annotation application produced missing values or cell loss."
  )
}

for (
  column_name in c(
    "primary_cluster",
    "broad_cell_type",
    "subtype_label",
    "annotation_confidence",
    "annotation_review_flag",
    "primary_VDR_localization_eligible",
    "primary_broad_pseudobulk_eligible",
    "subtype_disease_inference_eligible",
    "annotation_version"
  )
) {
  integrated_object[[
    column_name
  ]] <- metadata[
    ,
    column_name,
    drop = TRUE
  ]
}

integrated_object@misc$
  GSE183837_broad_annotation <- list(
    annotation_version =
      annotation_version,
    script_version =
      script_version,
    annotation_role =
      "MARKER_SUPPORTED_BROAD_LINEAGE",
    broad_annotation_locked =
      TRUE,
    cluster8_policy =
      "LOW_CONFIDENCE_EPITHELIAL_LIKE_ATLAS_RETAINED_PRIMARY_MODELS_EXCLUDED",
    cluster13_policy =
      "CYCLING_STROMAL_IDENTITY_RETAINED_NO_DISEASE_SUBTYPE_CLAIM",
    source_marker_version =
      marker_version,
    original_checkpoint_md5 =
      unname(
        tools::md5sum(
          input_checkpoint
        )
      )
  )

saveRDS(
  integrated_object,
  annotated_checkpoint,
  compress = "gzip"
)

annotated_checkpoint_md5 <- unname(
  tools::md5sum(
    annotated_checkpoint
  )
)

annotation_status <- data.frame(
  annotation_version =
    annotation_version,
  status = "COMPLETE",
  cell_n = ncol(
    integrated_object
  ),
  cluster_n = length(
    unique(
      integrated_object$
        primary_cluster
    )
  ),
  broad_cell_type_n = length(
    unique(
      integrated_object$
        broad_cell_type
    )
  ),
  checkpoint_path =
    normalize_slashes(
      annotated_checkpoint
    ),
  checkpoint_size_bytes =
    file.info(
      annotated_checkpoint
    )$size,
  checkpoint_md5 =
    annotated_checkpoint_md5,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  annotation_status,
  annotation_status_path
)


# ----------------------------------------------------------------------
# 9. Cell, donor and group composition tables
# ----------------------------------------------------------------------

annotated_metadata <- integrated_object[[]]

broad_cell_counts <- as.data.frame(
  table(
    broad_cell_type =
      annotated_metadata$
        broad_cell_type
  ),
  stringsAsFactors = FALSE
)

colnames(broad_cell_counts)[2] <-
  "cell_n"

broad_cell_counts$cell_fraction <- (
  broad_cell_counts$cell_n /
    sum(
      broad_cell_counts$cell_n
    )
)

cluster_cell_counts <- as.data.frame(
  table(
    cluster =
      annotated_metadata$
        primary_cluster,
    broad_cell_type =
      annotated_metadata$
        broad_cell_type,
    subtype_label =
      annotated_metadata$
        subtype_label
  ),
  stringsAsFactors = FALSE
)

cluster_cell_counts <- cluster_cell_counts[
  cluster_cell_counts$Freq > 0,
  ,
  drop = FALSE
]

colnames(cluster_cell_counts)[4] <-
  "cell_n"

donor_broad_counts <- as.data.frame(
  table(
    sample_accession =
      annotated_metadata$
        sample_accession,
    donor_id =
      annotated_metadata$
        donor_id,
    analysis_group =
      annotated_metadata$
        analysis_group,
    broad_cell_type =
      annotated_metadata$
        broad_cell_type
  ),
  stringsAsFactors = FALSE
)

donor_broad_counts <- donor_broad_counts[
  donor_broad_counts$Freq > 0,
  ,
  drop = FALSE
]

colnames(donor_broad_counts)[5] <-
  "cell_n"

donor_totals <- aggregate(
  cell_n ~
    sample_accession +
    donor_id +
    analysis_group,
  data = donor_broad_counts,
  FUN = sum
)

colnames(donor_totals)[4] <-
  "donor_total_cell_n"

donor_broad_counts <- merge(
  donor_broad_counts,
  donor_totals,
  by = c(
    "sample_accession",
    "donor_id",
    "analysis_group"
  ),
  all.x = TRUE,
  sort = FALSE
)

donor_broad_counts$
  within_donor_fraction <- (
    donor_broad_counts$cell_n /
      donor_broad_counts$
        donor_total_cell_n
  )

group_broad_summary <- data.table::as.data.table(
  donor_broad_counts
)[
  ,
  .(
    donor_n =
      data.table::uniqueN(
        donor_id
      ),
    median_within_donor_fraction =
      median(
        within_donor_fraction
      ),
    minimum_within_donor_fraction =
      min(
        within_donor_fraction
      ),
    maximum_within_donor_fraction =
      max(
        within_donor_fraction
      )
  ),
  by = .(
    analysis_group,
    broad_cell_type
  )
]

group_broad_summary <- as.data.frame(
  group_broad_summary
)

eligibility_counts <- as.data.frame(
  table(
    cluster =
      annotated_metadata$
        primary_cluster,
    broad_cell_type =
      annotated_metadata$
        broad_cell_type,
    VDR_localization_eligible =
      annotated_metadata$
        primary_VDR_localization_eligible,
    pseudobulk_eligible =
      annotated_metadata$
        primary_broad_pseudobulk_eligible
  ),
  stringsAsFactors = FALSE
)

eligibility_counts <- eligibility_counts[
  eligibility_counts$Freq > 0,
  ,
  drop = FALSE
]

colnames(eligibility_counts)[5] <-
  "cell_n"

write_csv_utf8(
  broad_cell_counts,
  file.path(
    table_root,
    "04_broad_cell_counts.csv"
  )
)

write_csv_utf8(
  cluster_cell_counts,
  file.path(
    table_root,
    "05_cluster_annotation_counts.csv"
  )
)

write_csv_utf8(
  donor_broad_counts,
  file.path(
    table_root,
    "06_donor_broad_composition.csv"
  )
)

write_csv_utf8(
  group_broad_summary,
  file.path(
    table_root,
    "07_group_broad_descriptive.csv"
  )
)

write_csv_utf8(
  eligibility_counts,
  file.path(
    table_root,
    "08_analysis_eligibility.csv"
  )
)


# ----------------------------------------------------------------------
# 10. Marker panel averages for annotation visualization
# ----------------------------------------------------------------------

visual_marker_list <- c(
  "COL1A1",
  "DCN",
  "APOD",
  "MKI67",
  "RGS5",
  "CD3D",
  "IL7R",
  "NKG7",
  "GNLY",
  "KRT8",
  "PAX8",
  "SCGB1D2",
  "FOXJ1",
  "C20orf85",
  "LYZ",
  "C1QC",
  "CLDN5",
  "VWF"
)

available_features <- rownames(
  integrated_object
)

feature_match <- match(
  toupper(
    visual_marker_list
  ),
  toupper(
    available_features
  )
)

visual_marker_manifest <- data.frame(
  requested_gene =
    visual_marker_list,
  resolved_feature = ifelse(
    is.na(feature_match),
    NA_character_,
    available_features[
      feature_match
    ]
  ),
  available =
    !is.na(feature_match),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  visual_marker_manifest,
  file.path(
    table_root,
    "09_visual_marker_manifest.csv"
  )
)

resolved_visual_features <- unique(
  visual_marker_manifest$
    resolved_feature[
      visual_marker_manifest$
        available
    ]
)

expression_matrix <- SeuratObject::LayerData(
  object = integrated_object,
  assay = "RNA",
  layer = "data",
  features =
    resolved_visual_features
)

marker_average_rows <- list()

for (cluster_id in as.character(0:15)) {
  cell_index <- which(
    annotated_metadata$
      primary_cluster ==
      cluster_id
  )

  average_expression <- Matrix::rowMeans(
    expression_matrix[
      ,
      cell_index,
      drop = FALSE
    ]
  )

  percent_expressing <- Matrix::rowMeans(
    expression_matrix[
      ,
      cell_index,
      drop = FALSE
    ] > 0
  )

  marker_average_rows[[
    length(marker_average_rows) + 1L
  ]] <- data.frame(
    cluster = cluster_id,
    gene =
      rownames(
        expression_matrix
      ),
    average_expression =
      as.numeric(
        average_expression
      ),
    percent_expressing =
      as.numeric(
        percent_expressing
      ),
    stringsAsFactors = FALSE
  )
}

marker_average <- do.call(
  rbind,
  marker_average_rows
)

marker_average$gene_z <- ave(
  log1p(
    marker_average$
      average_expression
  ),
  marker_average$gene,
  FUN = safe_zscore
)

marker_average$gene_z <- pmax(
  -3,
  pmin(
    3,
    marker_average$gene_z
  )
)

write_csv_utf8(
  marker_average,
  file.path(
    table_root,
    "10_annotation_marker_average.csv"
  )
)


# ----------------------------------------------------------------------
# 11. QC figures
# ----------------------------------------------------------------------

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

UMAP_embedding <- SeuratObject::Embeddings(
  integrated_object,
  reduction = "umap.harmony"
)

UMAP_metadata <- data.frame(
  cell_barcode =
    rownames(UMAP_embedding),
  UMAPH_1 =
    UMAP_embedding[, 1],
  UMAPH_2 =
    UMAP_embedding[, 2],
  broad_cell_type =
    integrated_object$
      broad_cell_type,
  subtype_label =
    integrated_object$
      subtype_label,
  primary_cluster =
    integrated_object$
      primary_cluster,
  annotation_confidence =
    integrated_object$
      annotation_confidence,
  stringsAsFactors = FALSE
)

plot_broad_UMAP <- ggplot2::ggplot(
  UMAP_metadata,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2,
    group = broad_cell_type
  )
) +
  ggplot2::geom_point(
    ggplot2::aes(
      shape = broad_cell_type
    ),
    size = 0.18,
    alpha = 0.55
  ) +
  ggplot2::labs(
    title = "GSE183837 broad cell annotation",
    subtitle = "Cluster 8 is retained as low-confidence epithelial-like and excluded from primary models",
    x = "Harmony UMAP 1",
    y = "Harmony UMAP 2",
    shape = "Broad cell type"
  ) +
  theme_QC

save_plot_pair(
  plot_broad_UMAP,
  "QC1_UMAP_broad_annotation",
  width = 10,
  height = 7
)

cluster_centroids <- aggregate(
  cbind(
    UMAPH_1,
    UMAPH_2
  ) ~
    primary_cluster +
    subtype_label,
  data = UMAP_metadata,
  FUN = median
)

plot_cluster_UMAP <- ggplot2::ggplot(
  UMAP_metadata,
  ggplot2::aes(
    x = UMAPH_1,
    y = UMAPH_2
  )
) +
  ggplot2::geom_point(
    size = 0.12,
    alpha = 0.30
  ) +
  ggplot2::geom_text(
    data = cluster_centroids,
    ggplot2::aes(
      label = primary_cluster
    ),
    size = 3.5
  ) +
  ggplot2::labs(
    title = "Primary clusters underlying the locked broad annotation",
    subtitle = "Numbers correspond to the reviewed cluster-to-lineage map",
    x = "Harmony UMAP 1",
    y = "Harmony UMAP 2"
  ) +
  theme_QC

save_plot_pair(
  plot_cluster_UMAP,
  "QC2_UMAP_cluster_labels",
  width = 9,
  height = 7
)

donor_plot_data <- donor_broad_counts

donor_plot_data$sample_accession <- factor(
  donor_plot_data$sample_accession,
  levels = unique(
    donor_plot_data$
      sample_accession
  )
)

plot_donor_composition <- ggplot2::ggplot(
  donor_plot_data,
  ggplot2::aes(
    x = sample_accession,
    y = within_donor_fraction,
    fill = broad_cell_type
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  ggplot2::labs(
    title = "Broad-cell composition by donor",
    subtitle = "Descriptive only; donor is the biological replicate",
    x = "Donor/sample",
    y = "Within-donor cell fraction",
    fill = "Broad cell type"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
  )

save_plot_pair(
  plot_donor_composition,
  "QC3_donor_broad_composition",
  width = 10,
  height = 6
)

marker_plot_data <- marker_average

marker_plot_data$cluster <- factor(
  marker_plot_data$cluster,
  levels = as.character(0:15)
)

marker_plot_data$gene <- factor(
  marker_plot_data$gene,
  levels = rev(
    resolved_visual_features
  )
)

plot_marker_heatmap <- ggplot2::ggplot(
  marker_plot_data,
  ggplot2::aes(
    x = cluster,
    y = gene,
    fill = gene_z
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    option = "C"
  ) +
  ggplot2::labs(
    title = "Canonical marker support for the broad annotation",
    subtitle = "Row-scaled log-normalized average expression",
    x = "Primary cluster",
    y = "Marker gene",
    fill = "Row z-score"
  ) +
  theme_QC

save_plot_pair(
  plot_marker_heatmap,
  "QC4_annotation_marker_heatmap",
  width = 10,
  height = 8
)

confidence_plot_data <- annotation_QC

confidence_plot_data$cluster <- factor(
  confidence_plot_data$cluster,
  levels = as.character(0:15)
)

plot_confidence <- ggplot2::ggplot(
  confidence_plot_data,
  ggplot2::aes(
    x = cluster,
    y = top100_evidence_n,
    shape = annotation_confidence,
    size = cell_n
  )
) +
  ggplot2::geom_point(
    alpha = 0.8
  ) +
  ggplot2::labs(
    title = "Annotation evidence and confidence by cluster",
    subtitle = "Cluster 8 remains low confidence and sensitivity-only",
    x = "Primary cluster",
    y = "Requested evidence markers present in top 100",
    shape = "Confidence",
    size = "Cell count"
  ) +
  theme_QC

save_plot_pair(
  plot_confidence,
  "QC5_annotation_confidence",
  width = 9,
  height = 5.5
)

eligibility_plot_data <- cluster_cell_counts

eligibility_plot_data <- merge(
  eligibility_plot_data,
  annotation_map[
    ,
    c(
      "cluster",
      "primary_broad_pseudobulk_eligible"
    ),
    drop = FALSE
  ],
  by = "cluster",
  all.x = TRUE,
  sort = FALSE
)

eligibility_plot_data$cluster <- factor(
  eligibility_plot_data$cluster,
  levels = as.character(0:15)
)

plot_eligibility <- ggplot2::ggplot(
  eligibility_plot_data,
  ggplot2::aes(
    x = cluster,
    y = cell_n,
    shape =
      primary_broad_pseudobulk_eligible
  )
) +
  ggplot2::geom_point(
    size = 3
  ) +
  ggplot2::scale_y_log10(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Cluster sizes and primary broad-pseudobulk eligibility",
    subtitle = "Cluster 8 is retained in the atlas but excluded from primary models",
    x = "Primary cluster",
    y = "Cell count (log scale)",
    shape = "Primary eligible"
  ) +
  theme_QC

save_plot_pair(
  plot_eligibility,
  "QC6_primary_model_eligibility",
  width = 8,
  height = 5
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
# 12. Annotation Gate
# ----------------------------------------------------------------------

cluster8_row <- annotation_map[
  annotation_map$cluster == "8",
  ,
  drop = FALSE
]

cluster13_row <- annotation_map[
  annotation_map$cluster == "13",
  ,
  drop = FALSE
]

gate_checks <- data.frame(
  check = c(
    "upstream_marker_gate_authorizes_annotation_review",
    "input_checkpoint_has_47161_cells",
    "annotation_map_has_16_unique_clusters",
    "all_cells_receive_broad_annotation",
    "all_cells_receive_subtype_label",
    "allowed_broad_ontology_only",
    "all_nonlow_confidence_clusters_have_at_least_2_top100_evidence_markers",
    "low_confidence_cluster8_has_at_least_1_top100_evidence_marker",
    "cluster8_is_retained_but_primary_models_excluded",
    "cluster13_is_cycling_stromal_with_donor_warning",
    "annotated_checkpoint_has_47161_cells",
    "annotated_checkpoint_exists",
    "annotated_checkpoint_md5_present",
    "donor_broad_composition_has_9_donors",
    "broad_annotation_can_be_locked",
    "VDR_localization_can_be_authorized",
    "disease_testing_remains_prohibited"
  ),
  pass = c(
    upstream_gate$gate_status[[1]] ==
      "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED",
    ncol(integrated_object) == 47161L,
    nrow(annotation_map) == 16L &&
      anyDuplicated(
        annotation_map$cluster
      ) == 0L,
    !anyNA(
      integrated_object$
        broad_cell_type
    ),
    !anyNA(
      integrated_object$
        subtype_label
    ),
    all(
      annotation_map$
        broad_cell_type %in%
        allowed_broad_types
    ),
    all(
      annotation_QC$
        evidence_pass[
          annotation_QC$
            annotation_confidence !=
            "LOW"
        ]
    ),
    isTRUE(
      annotation_QC$
        evidence_pass[
          annotation_QC$cluster ==
            "8"
        ]
    ),
    nrow(cluster8_row) == 1L &&
      !cluster8_row$
        primary_VDR_localization_eligible &&
      !cluster8_row$
        primary_broad_pseudobulk_eligible,
    nrow(cluster13_row) == 1L &&
      cluster13_row$
        broad_cell_type ==
        "Stromal cells" &&
      grepl(
        "DONOR",
        cluster13_row$review_flag
      ),
    ncol(integrated_object) ==
      47161L,
    file.exists(
      annotated_checkpoint
    ),
    nzchar(
      annotated_checkpoint_md5
    ),
    length(
      unique(
        donor_broad_counts$donor_id
      )
    ) == 9L,
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
  "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION"
} else {
  "HOLD_BROAD_ANNOTATION_NOT_LOCKED"
}

annotation_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version =
    marker_version,
  annotation_version =
    annotation_version,
  gate_status = gate_status,
  cell_n = ncol(
    integrated_object
  ),
  cluster_n = 16L,
  broad_cell_type_n =
    length(
      unique(
        integrated_object$
          broad_cell_type
      )
    ),
  broad_cell_types =
    collapse_unique(
      sort(
        unique(
          integrated_object$
            broad_cell_type
        )
      )
    ),
  low_confidence_clusters =
    collapse_unique(
      annotation_map$cluster[
        annotation_map$
          annotation_confidence ==
          "LOW"
      ]
    ),
  primary_excluded_clusters =
    collapse_unique(
      annotation_map$cluster[
        !annotation_map$
          primary_broad_pseudobulk_eligible
      ]
    ),
  subtype_disease_restricted_clusters =
    collapse_unique(
      annotation_map$cluster[
        !annotation_map$
          subtype_disease_inference_eligible
      ]
    ),
  review_flags =
    collapse_unique(
      annotation_map$review_flag
    ),
  broad_annotation_locked =
    all_gate_checks_pass,
  VDR_localization_authorized =
    all_gate_checks_pass,
  broad_pseudobulk_preparation_authorized =
    all_gate_checks_pass,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  annotated_checkpoint =
    normalize_slashes(
      annotated_checkpoint
    ),
  annotated_checkpoint_md5 =
    annotated_checkpoint_md5,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  annotation_gate,
  file.path(
    run_qc_root,
    "G183_ANNOTATION_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1M GSE183837 Broad Annotation Card",
  "",
  paste0("- Gate: ", gate_status),
  paste0(
    "- Cells retained: ",
    format(
      ncol(integrated_object),
      big.mark = ","
    )
  ),
  "- Primary clusters: 16",
  paste0(
    "- Broad cell types: ",
    length(
      unique(
        integrated_object$
          broad_cell_type
      )
    )
  ),
  paste0(
    "- Broad ontology: ",
    collapse_unique(
      sort(
        unique(
          integrated_object$
            broad_cell_type
        )
      )
    )
  ),
  "- Low-confidence cluster: 8",
  "- Primary model exclusion: cluster 8",
  "- Donor-dominant subtype restriction: cluster 13",
  paste0(
    "- Broad annotation locked: ",
    all_gate_checks_pass
  ),
  paste0(
    "- VDR localization authorized: ",
    all_gate_checks_pass
  ),
  "- Disease pseudobulk testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "G183_ANNOTATION_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 broad annotation interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Locked broad annotation",
  "",
  "- Stromal cells: clusters 0, 1, 2, 3, 4, 5, 10 and 13.",
  "- NK/T cells: clusters 6, 9 and 14.",
  "- Unciliated epithelial cells: clusters 7 and 8.",
  "- Myeloid cells: cluster 11.",
  "- Endothelial cells: cluster 12.",
  "- Ciliated epithelial cells: cluster 15.",
  "",
  "## Mandatory restrictions",
  "",
  "- Cluster 8 is retained as a low-confidence SCGB1D2-positive epithelial-like population but excluded from primary VDR localization and broad pseudobulk models.",
  "- Cluster 13 is retained as cycling stromal cells; its subtype cannot support an RIF-specific claim because of severe single-donor dominance.",
  "- Clusters 7 and 15 cannot support subtype-level disease inference at this stage.",
  "",
  "## Authorized next operation",
  "",
  "- Prepare eligible broad-cell pseudobulks by donor.",
  "- Localize VDR expression and the four locked VDR programs.",
  "- Run donor-level RIF comparisons only after pseudobulk QC and with explicit low-power wording.",
  "",
  "## Not authorized",
  "",
  "- No cell-level Control-versus-RIF inference.",
  "- No cluster-8 primary model.",
  "- No disease mechanism claim.",
  "- No final result lock."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "G183_ANNOTATION_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 13. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version = marker_version,
  annotation_version =
    annotation_version,
  upstream_marker_gate =
    upstream_gate$gate_status[[1]],
  broad_ontology =
    allowed_broad_types,
  cluster8_policy =
    "ATLAS_RETAINED_LOW_CONFIDENCE_PRIMARY_VDR_AND_PSEUDOBULK_EXCLUDED",
  cluster13_policy =
    "CYCLING_STROMAL_RETAINED_NO_DISEASE_SUBTYPE_CLAIM",
  broad_annotation_locked =
    all_gate_checks_pass,
  VDR_localization_authorized =
    all_gate_checks_pass,
  disease_pseudobulk_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(sessionInfo()),
  session_path
)

run_end <- Sys.time()

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  marker_version =
    marker_version,
  annotation_version =
    annotation_version,
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
  cell_n = ncol(
    integrated_object
  ),
  cluster_n = 16L,
  broad_cell_type_n =
    length(
      unique(
        integrated_object$
          broad_cell_type
      )
    ),
  low_confidence_clusters = "8",
  primary_excluded_clusters = "8",
  broad_annotation_locked =
    all_gate_checks_pass,
  VDR_localization_authorized =
    all_gate_checks_pass,
  disease_pseudobulk_testing_authorized =
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
  annotation_status_path,
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
  "cells=",
  ncol(integrated_object),
  "broad_types=",
  length(
    unique(
      integrated_object$
        broad_cell_type
    )
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
cat("STEP04E1M GSE183837 broad annotation QC completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Cells retained: ",
  format(
    ncol(integrated_object),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat("Primary clusters annotated: 16\n")
cat(
  "Broad cell types: ",
  length(
    unique(
      integrated_object$
        broad_cell_type
    )
  ),
  "\n",
  sep = ""
)
cat("Low-confidence cluster: 8\n")
cat("Primary model exclusion: 8\n")
cat(
  "Broad annotation locked: ",
  all_gate_checks_pass,
  "\n",
  sep = ""
)
cat(
  "VDR localization authorized: ",
  all_gate_checks_pass,
  "\n",
  sep = ""
)
cat(
  "Annotated checkpoint: ",
  annotated_checkpoint,
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
    "\nAnnotation outputs are complete in the project folders.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 6: 68_step04E1N_GSE183837_VDR_localization_pseudobulk_v1_1_unique_columns_fixed
# ==============================================================================

# ======================================================================
# 68_step04E1N_GSE183837_VDR_localization_pseudobulk_v1_1_unique_columns_fixed.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# 1. Aggregate raw counts by donor x locked broad cell type.
# 2. Preserve a primary matrix excluding low-confidence cluster 8.
# 3. Preserve an atlas-sensitivity matrix retaining cluster 8.
# 4. Localize VDR expression and four locked signed VDR programs.
# 5. Determine which broad cell types support donor-level
#    Control-versus-RIF models.
#
# LOCKED VDR PROGRAMS
# - RNA_ALL_436
# - OCCUPANCY_SUPPORTED_284
# - SHARED_OCCUPANCY_46
# - PROMOTER_OCCUPANCY_84
#
# SCORE DEFINITION
# - Raw pseudobulk counts are normalized to log2(CPM + 1).
# - Each candidate gene is standardized across primary donor x broad-type
#   profiles.
# - Standardized expression is multiplied by the sign of the locked
#   calcitriol-dependent VDR interaction logFC.
# - The program score is the mean signed standardized expression.
#
# PRIMARY / SENSITIVITY POLICY
# - Primary matrix: cluster 8 excluded.
# - Atlas sensitivity: cluster 8 retained.
# - Cluster 8 is never used in primary donor-level models.
# - Ciliated epithelial cells are descriptive when donor cell counts are
#   below the prespecified threshold.
#
# INTERPRETATION
# - Scores quantify relative alignment with the locked VDR interaction
#   direction.
# - Scores do not establish direct binding, causal activation or disease
#   differential expression.
# - No Control-versus-RIF statistical test is performed in this script.
#
# MEMORY / RESUME
# - Each of the nine split counts layers is aggregated independently.
# - Per-layer RDS files are written immediately.
# - Re-running skips validated completed layers.
# - Counts layers are never joined.
# - The large Seurat object is removed before program-score calculation.
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
options(future.globals.maxSize = 12 * 1024^3)

set.seed(20260806)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
dataset_id <- "GSE183837"

script_name <-
  "68_step04E1N_GSE183837_VDR_localization_pseudobulk_v1_1_unique_columns_fixed.R"

script_version <-
  "STEP04E1N_GSE183837_VDR_LOCALIZATION_PSEUDOBULK_v1.1_UNIQUE_COLUMNS_FIXED"

integration_version <- "SCI183_v1_0"
annotation_version <- "ANN183_v1_0"
pseudobulk_version <- "PB183_v1_0"

minimum_cells_primary <- 20L
expected_primary_cell_n <- 45630L
expected_atlas_cell_n <- 47161L
expected_cluster8_cell_n <- 1531L

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1N_",
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

  if (!dir.exists(dirname(path))) {
    stop(
      "Could not create output directory:\n",
      dirname(path)
    )
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

sanitize_token <- function(x) {
  output <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    as.character(x)
  )

  output <- gsub(
    "^_+|_+$",
    "",
    output
  )

  substr(output, 1L, 60L)
}

row_zscore_reference <- function(x) {
  x <- as.matrix(x)

  row_mean <- rowMeans(
    x,
    na.rm = TRUE
  )

  row_sd <- apply(
    x,
    1L,
    stats::sd,
    na.rm = TRUE
  )

  list(
    mean = row_mean,
    sd = row_sd,
    valid =
      is.finite(row_sd) &
      row_sd > 0
  )
}

apply_row_reference_z <- function(
  x,
  reference_mean,
  reference_sd
) {
  x <- as.matrix(x)

  if (
    !identical(
      rownames(x),
      names(reference_mean)
    ) ||
      !identical(
        rownames(x),
        names(reference_sd)
      )
  ) {
    stop(
      "Expression matrix and row-standardization reference are not aligned."
    )
  }

  output <- matrix(
    NA_real_,
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )

  valid <- is.finite(
    reference_sd
  ) &
    reference_sd > 0

  output[valid, ] <- sweep(
    sweep(
      x[
        valid,
        ,
        drop = FALSE
      ],
      MARGIN = 1,
      STATS =
        reference_mean[valid],
      FUN = "-"
    ),
    MARGIN = 1,
    STATS =
      reference_sd[valid],
    FUN = "/"
  )

  output
}

get_free_memory_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(
      trimws(output[1])
    )
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(
      value / 1024^2,
      2
    )
  }
}

get_free_disk_gb <- function() {
  if (.Platform$OS.type != "windows") {
    return(NA_real_)
  }

  command <- paste(
    "(Get-CimInstance Win32_LogicalDisk",
    "-Filter \"DeviceID='C:'\").FreeSpace"
  )

  output <- tryCatch(
    system2(
      "powershell",
      args = c(
        "-NoProfile",
        "-Command",
        shQuote(command)
      ),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character()
  )

  value <- suppressWarnings(
    as.numeric(
      trimws(output[1])
    )
  )

  if (
    length(value) != 1L ||
      !is.finite(value)
  ) {
    NA_real_
  } else {
    round(
      value / 1024^3,
      2
    )
  }
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
  "SeuratObject",
  "Matrix",
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

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
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

integration_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

annotated_checkpoint <- file.path(
  integration_root,
  "ann",
  "G183",
  annotation_version,
  "annotated.rds"
)

annotation_status_path <- file.path(
  integration_root,
  "ann",
  "G183",
  annotation_version,
  "status.csv"
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
  integration_root,
  "pb",
  "G183",
  pseudobulk_version
)

layer_root <- file.path(
  output_root,
  "lyr"
)

object_root <- file.path(
  output_root,
  "obj"
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

primary_counts_path <- file.path(
  object_root,
  "pb_primary_counts.rds"
)

atlas_counts_path <- file.path(
  object_root,
  "pb_atlas_counts.rds"
)

pb_status_path <- file.path(
  object_root,
  "pb_status.csv"
)

run_qc_root <- file.path(
  qc_base_root,
  paste0(
    "STEP04E1N_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_pbQC"
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
  layer_root,
  object_root,
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
  artifact = c(
    "annotated_checkpoint",
    "annotation_status",
    "candidate_table",
    "output_root",
    "layer_root",
    "primary_counts",
    "atlas_counts",
    "pb_status",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    annotated_checkpoint,
    annotation_status_path,
    candidate_path,
    output_root,
    layer_root,
    primary_counts_path,
    atlas_counts_path,
    pb_status_path,
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
    "One or more pseudobulk paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream annotation Gate and inputs
# ----------------------------------------------------------------------

annotation_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

annotation_qc_dirs <- annotation_qc_dirs[
  grepl(
    "^STEP04E1M_[0-9]{8}_[0-9]{6}_annQC$",
    basename(annotation_qc_dirs),
    perl = TRUE
  )
]

if (!length(annotation_qc_dirs)) {
  stop(
    "No STEP04E1M annotation-QC directory was found."
  )
}

gate_candidates <- list()

for (directory_path in annotation_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "G183_ANNOTATION_GATE.csv"
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

    gate_candidates[[
      length(gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(gate_candidates)) {
  stop(
    "No readable annotation Gate was found."
  )
}

gate_candidates <- as.data.frame(
  data.table::rbindlist(
    gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    gate_candidates$source_mtime
  )

gate_candidates <- gate_candidates[
  order(
    gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

authorized_index <- which(
  gate_candidates$gate_status ==
    "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION" &
    toupper(
      as.character(
        gate_candidates$
          broad_annotation_locked
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      ) &
    toupper(
      as.character(
        gate_candidates$
          VDR_localization_authorized
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      ) &
    toupper(
      as.character(
        gate_candidates$
          broad_pseudobulk_preparation_authorized
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      )
)

if (!length(authorized_index)) {
  stop(
    "No annotation Gate authorizes VDR localization and pseudobulk preparation."
  )
}

upstream_gate <- gate_candidates[
  authorized_index[1],
  ,
  drop = FALSE
]

write_csv_utf8(
  upstream_gate,
  file.path(
    run_qc_root,
    "00_upstream_annotation_gate.csv"
  )
)

required_inputs <- c(
  annotated_checkpoint,
  annotation_status_path,
  candidate_path
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(
      normalize_slashes(
        missing_inputs
      ),
      collapse = "\n"
    )
  )
}

annotation_status <- safe_read_csv(
  annotation_status_path
)

candidate_table <- safe_read_csv(
  candidate_path
)

if (
  is.null(annotation_status) ||
    nrow(annotation_status) != 1L ||
    annotation_status$status[[1]] !=
      "COMPLETE" ||
    annotation_status$cell_n[[1]] !=
      expected_atlas_cell_n
) {
  stop(
    "Annotated checkpoint status is absent or invalid."
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
    "Locked VDR candidate table is missing: ",
    paste(
      missing_candidate_columns,
      collapse = ", "
    )
  )
}

candidate_table$feature_id <- as.character(
  candidate_table$feature_id
)

candidate_table$interaction_logFC <-
  suppressWarnings(
    as.numeric(
      candidate_table$
        interaction_logFC
    )
  )

for (logical_column in c(
  "occupancy_supported",
  "shared_occupancy_supported",
  "any_promoter_overlap_support"
)) {
  candidate_table[[
    logical_column
  ]] <- as_logical_safe(
    candidate_table[[
      logical_column
    ]]
  )
}

if (
  nrow(candidate_table) != 436L ||
    anyDuplicated(
      candidate_table$feature_id
    ) > 0L ||
    sum(
      candidate_table$
        occupancy_supported
    ) != 284L ||
    sum(
      candidate_table$
        shared_occupancy_supported
    ) != 46L ||
    sum(
      candidate_table$
        any_promoter_overlap_support
    ) != 84L
) {
  stop(
    "Locked VDR program identity does not match 436/284/46/84."
  )
}


# ----------------------------------------------------------------------
# 5. Load and validate annotated object
# ----------------------------------------------------------------------

free_memory_before_GB <- get_free_memory_gb()
free_disk_before_GB <- get_free_disk_gb()

if (
  is.finite(
    free_memory_before_GB
  ) &&
    free_memory_before_GB < 10
) {
  stop(
    "Free physical memory is below 10 GB: ",
    free_memory_before_GB
  )
}

log_event(
  "Starting",
  run_id,
  "free_memory_GB=",
  free_memory_before_GB,
  "free_disk_GB=",
  free_disk_before_GB
)

annotated_object <- readRDS(
  annotated_checkpoint
)

if (!inherits(annotated_object, "Seurat")) {
  stop(
    "Annotated checkpoint is not a Seurat object."
  )
}

SeuratObject::DefaultAssay(
  annotated_object
) <- "RNA"

metadata <- annotated_object[[]]

required_metadata_columns <- c(
  "sample_accession",
  "donor_id",
  "analysis_group",
  "primary_cluster",
  "broad_cell_type",
  "annotation_confidence",
  "primary_VDR_localization_eligible",
  "primary_broad_pseudobulk_eligible"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns)) {
  stop(
    "Annotated metadata is missing: ",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

metadata$sample_accession <- as.character(
  metadata$sample_accession
)

metadata$donor_id <- as.character(
  metadata$donor_id
)

metadata$analysis_group <- as.character(
  metadata$analysis_group
)

metadata$primary_cluster <- as.character(
  metadata$primary_cluster
)

metadata$broad_cell_type <- as.character(
  metadata$broad_cell_type
)

metadata$
  primary_VDR_localization_eligible <-
  as_logical_safe(
    metadata$
      primary_VDR_localization_eligible
  )

metadata$
  primary_broad_pseudobulk_eligible <-
  as_logical_safe(
    metadata$
      primary_broad_pseudobulk_eligible
  )

if (
  ncol(annotated_object) !=
    expected_atlas_cell_n ||
    nrow(metadata) !=
      expected_atlas_cell_n ||
    length(
      unique(
        metadata$sample_accession
      )
    ) != 9L ||
    length(
      unique(
        metadata$donor_id
      )
    ) != 9L ||
    sum(
      metadata$
        primary_broad_pseudobulk_eligible
    ) != expected_primary_cell_n ||
    sum(
      metadata$primary_cluster ==
        "8"
    ) != expected_cluster8_cell_n
) {
  stop(
    "Annotated object identity or cluster-8 exclusion count is inconsistent."
  )
}

donor_group_identity <- unique(
  metadata[
    ,
    c(
      "sample_accession",
      "donor_id",
      "analysis_group"
    ),
    drop = FALSE
  ]
)

donor_group_identity <- donor_group_identity[
  order(
    donor_group_identity$
      sample_accession
  ),
  ,
  drop = FALSE
]

if (
  nrow(donor_group_identity) != 9L ||
    sum(
      donor_group_identity$
        analysis_group ==
        "CONTROL"
    ) != 3L ||
    sum(
      donor_group_identity$
        analysis_group ==
        "RIF"
    ) != 6L ||
    anyDuplicated(
      donor_group_identity$donor_id
    ) > 0L
) {
  stop(
    "Donor/group identity does not match 3 Control and 6 RIF donors."
  )
}

count_layers <- SeuratObject::Layers(
  annotated_object,
  assay = "RNA",
  search = "^counts"
)

data_layers <- SeuratObject::Layers(
  annotated_object,
  assay = "RNA",
  search = "^data"
)

if (
  length(count_layers) != 9L ||
    length(data_layers) != 1L ||
    data_layers[1] != "data"
) {
  stop(
    "Expected nine split count layers and one joined data layer."
  )
}

write_csv_utf8(
  donor_group_identity,
  file.path(
    table_root,
    "01_donor_group_identity.csv"
  )
)


# ----------------------------------------------------------------------
# 6. Per-layer resumable raw-count aggregation
# ----------------------------------------------------------------------

layer_status_rows <- list()
layer_result_paths <- character(
  length(count_layers)
)

for (
  layer_index in seq_along(
    count_layers
  )
) {
  layer_name <- count_layers[
    layer_index
  ]

  layer_code <- sprintf(
    "L%02d",
    layer_index
  )

  layer_result_path <- file.path(
    layer_root,
    paste0(
      layer_code,
      ".rds"
    )
  )

  layer_status_path <- file.path(
    layer_root,
    paste0(
      layer_code,
      "_status.csv"
    )
  )

  layer_result_paths[
    layer_index
  ] <- layer_result_path

  existing_status <- safe_read_csv(
    layer_status_path
  )

  layer_result <- NULL
  resume_valid <- FALSE

  if (
    file.exists(layer_result_path) &&
      !is.null(existing_status) &&
      nrow(existing_status) == 1L &&
      existing_status$status[[1]] ==
        "COMPLETE" &&
      existing_status$layer_name[[1]] ==
        layer_name
  ) {
    layer_result <- tryCatch(
      readRDS(layer_result_path),
      error = function(e) NULL
    )

    resume_valid <- (
      is.list(layer_result) &&
        all(
          c(
            "primary_counts",
            "atlas_counts",
            "primary_metadata",
            "atlas_metadata",
            "cell_audit",
            "raw_count_audit"
          ) %in%
            names(layer_result)
        )
    )
  }

  if (resume_valid) {
    log_event(
      "SKIP completed counts layer",
      layer_index,
      "/",
      length(count_layers),
      layer_name
    )
  } else {
    layer_start <- Sys.time()

    log_event(
      "START counts layer",
      layer_index,
      "/",
      length(count_layers),
      layer_name
    )

    layer_counts <- SeuratObject::LayerData(
      object = annotated_object,
      assay = "RNA",
      layer = layer_name
    )

    layer_cells <- colnames(
      layer_counts
    )

    if (
      is.null(layer_cells) ||
        !length(layer_cells) ||
        !all(
          layer_cells %in%
            rownames(metadata)
        )
    ) {
      stop(
        "Counts-layer cell identities are missing or invalid for ",
        layer_name,
        "."
      )
    }

    layer_metadata <- metadata[
      layer_cells,
      ,
      drop = FALSE
    ]

    layer_samples <- unique(
      layer_metadata$
        sample_accession
    )

    layer_donors <- unique(
      layer_metadata$donor_id
    )

    layer_groups <- unique(
      layer_metadata$
        analysis_group
    )

    if (
      length(layer_samples) != 1L ||
        length(layer_donors) != 1L ||
        length(layer_groups) != 1L
    ) {
      stop(
        "Counts layer does not map to exactly one sample/donor/group: ",
        layer_name
      )
    }

    aggregate_by_type <- function(
      counts_matrix,
      metadata_table,
      inclusion_vector,
      matrix_role
    ) {
      type_values <- sort(
        unique(
          metadata_table$
            broad_cell_type[
              inclusion_vector
            ]
        )
      )

      if (!length(type_values)) {
        stop(
          "No cells remain for ",
          matrix_role,
          " in ",
          layer_name,
          "."
        )
      }

      count_columns <- list()
      metadata_rows <- list()

      for (
        type_index in seq_along(
          type_values
        )
      ) {
        cell_type <- type_values[
          type_index
        ]

        selected_index <- which(
          inclusion_vector &
            metadata_table$
              broad_cell_type ==
              cell_type
        )

        if (!length(selected_index)) {
          next
        }

        pb_id <- paste(
          dataset_id,
          layer_samples,
          sanitize_token(cell_type),
          matrix_role,
          sep = "__"
        )

        aggregated_count <- Matrix::rowSums(
          counts_matrix[
            ,
            selected_index,
            drop = FALSE
          ]
        )

        count_columns[[
          length(count_columns) + 1L
        ]] <- Matrix::Matrix(
          aggregated_count,
          ncol = 1L,
          sparse = TRUE,
          dimnames = list(
            rownames(counts_matrix),
            pb_id
          )
        )

        metadata_rows[[
          length(metadata_rows) + 1L
        ]] <- data.frame(
          pb_id = pb_id,
          matrix_role =
            matrix_role,
          sample_accession =
            layer_samples,
          donor_id =
            layer_donors,
          analysis_group =
            layer_groups,
          broad_cell_type =
            cell_type,
          cell_n =
            length(selected_index),
          library_size =
            as.numeric(
              sum(aggregated_count)
            ),
          minimum_cell_threshold =
            minimum_cells_primary,
          minimum_cell_threshold_pass =
            length(selected_index) >=
              minimum_cells_primary,
          stringsAsFactors = FALSE
        )
      }

      list(
        counts = do.call(
          cbind,
          count_columns
        ),
        metadata = do.call(
          rbind,
          metadata_rows
        )
      )
    }

    primary_result <- aggregate_by_type(
      counts_matrix =
        layer_counts,
      metadata_table =
        layer_metadata,
      inclusion_vector =
        layer_metadata$
          primary_broad_pseudobulk_eligible,
      matrix_role =
        "PRIMARY"
    )

    atlas_result <- aggregate_by_type(
      counts_matrix =
        layer_counts,
      metadata_table =
        layer_metadata,
      inclusion_vector =
        rep(
          TRUE,
          nrow(layer_metadata)
        ),
      matrix_role =
        "ATLAS"
    )

    VDR_index <- match(
      "VDR",
      toupper(
        rownames(layer_counts)
      )
    )

    if (is.na(VDR_index)) {
      stop(
        "VDR feature is absent from counts layer ",
        layer_name,
        "."
      )
    }

    VDR_cell_rows <- list()

    for (cell_type in sort(
      unique(
        layer_metadata$
          broad_cell_type
      )
    )) {
      type_index <- which(
        layer_metadata$
          broad_cell_type ==
          cell_type
      )

      primary_index <- type_index[
        layer_metadata$
          primary_VDR_localization_eligible[
            type_index
          ]
      ]

      atlas_VDR_counts <- as.numeric(
        layer_counts[
          VDR_index,
          type_index,
          drop = TRUE
        ]
      )

      primary_VDR_counts <- if (
        length(primary_index)
      ) {
        as.numeric(
          layer_counts[
            VDR_index,
            primary_index,
            drop = TRUE
          ]
        )
      } else {
        numeric()
      }

      VDR_cell_rows[[
        length(VDR_cell_rows) + 1L
      ]] <- data.frame(
        sample_accession =
          layer_samples,
        donor_id =
          layer_donors,
        analysis_group =
          layer_groups,
        broad_cell_type =
          cell_type,
        atlas_cell_n =
          length(type_index),
        atlas_VDR_positive_n =
          sum(
            atlas_VDR_counts > 0
          ),
        atlas_VDR_positive_fraction =
          mean(
            atlas_VDR_counts > 0
          ),
        primary_cell_n =
          length(primary_index),
        primary_VDR_positive_n =
          sum(
            primary_VDR_counts > 0
          ),
        primary_VDR_positive_fraction =
          if (
            length(primary_index)
          ) {
            mean(
              primary_VDR_counts > 0
            )
          } else {
            NA_real_
          },
        stringsAsFactors = FALSE
      )
    }

    VDR_cell_audit <- do.call(
      rbind,
      VDR_cell_rows
    )

    raw_count_audit <- data.frame(
      layer_name = layer_name,
      sample_accession =
        layer_samples,
      donor_id =
        layer_donors,
      analysis_group =
        layer_groups,
      feature_n =
        nrow(layer_counts),
      atlas_cell_n =
        ncol(layer_counts),
      primary_cell_n =
        sum(
          layer_metadata$
            primary_broad_pseudobulk_eligible
        ),
      cluster8_excluded_cell_n =
        sum(
          layer_metadata$
            primary_cluster ==
            "8"
        ),
      raw_count_sum =
        as.numeric(
          sum(layer_counts)
        ),
      atlas_aggregated_count_sum =
        as.numeric(
          sum(
            atlas_result$counts
          )
        ),
      primary_input_count_sum =
        as.numeric(
          sum(
            layer_counts[
              ,
              layer_metadata$
                primary_broad_pseudobulk_eligible,
              drop = FALSE
            ]
          )
        ),
      primary_aggregated_count_sum =
        as.numeric(
          sum(
            primary_result$counts
          )
        ),
      stringsAsFactors = FALSE
    )

    if (
      raw_count_audit$
        raw_count_sum !=
        raw_count_audit$
          atlas_aggregated_count_sum ||
        raw_count_audit$
          primary_input_count_sum !=
          raw_count_audit$
            primary_aggregated_count_sum
    ) {
      stop(
        "Raw-count conservation failed for ",
        layer_name,
        "."
      )
    }

    layer_result <- list(
      layer_name = layer_name,
      primary_counts =
        primary_result$counts,
      atlas_counts =
        atlas_result$counts,
      primary_metadata =
        primary_result$metadata,
      atlas_metadata =
        atlas_result$metadata,
      VDR_cell_audit =
        VDR_cell_audit,
      cell_audit = data.frame(
        layer_name = layer_name,
        sample_accession =
          layer_samples,
        atlas_cell_n =
          nrow(layer_metadata),
        primary_cell_n =
          sum(
            layer_metadata$
              primary_broad_pseudobulk_eligible
          ),
        cluster8_excluded_cell_n =
          sum(
            layer_metadata$
              primary_cluster ==
              "8"
          ),
        stringsAsFactors = FALSE
      ),
      raw_count_audit =
        raw_count_audit
    )

    saveRDS(
      layer_result,
      layer_result_path,
      compress = "gzip"
    )

    layer_status <- data.frame(
      layer_name = layer_name,
      status = "COMPLETE",
      sample_accession =
        layer_samples,
      donor_id =
        layer_donors,
      atlas_cell_n =
        nrow(layer_metadata),
      primary_cell_n =
        sum(
          layer_metadata$
            primary_broad_pseudobulk_eligible
        ),
      cluster8_excluded_cell_n =
        sum(
          layer_metadata$
            primary_cluster ==
            "8"
        ),
      primary_profile_n =
        ncol(
          primary_result$counts
        ),
      atlas_profile_n =
        ncol(
          atlas_result$counts
        ),
      elapsed_minutes = round(
        as.numeric(
          difftime(
            Sys.time(),
            layer_start,
            units = "mins"
          )
        ),
        2
      ),
      checkpoint_md5 =
        unname(
          tools::md5sum(
            layer_result_path
          )
        ),
      stringsAsFactors = FALSE
    )

    write_csv_utf8(
      layer_status,
      layer_status_path
    )

    existing_status <- layer_status

    log_event(
      "DONE counts layer",
      layer_index,
      "/",
      length(count_layers),
      layer_name,
      "primary_cells=",
      layer_status$
        primary_cell_n,
      "atlas_cells=",
      layer_status$
        atlas_cell_n
    )

    rm(
      layer_counts,
      primary_result,
      atlas_result
    )

    gc()
  }

  layer_status_rows[[
    layer_index
  ]] <- existing_status
}

layer_status_all <- do.call(
  rbind,
  layer_status_rows
)

write_csv_utf8(
  layer_status_all,
  file.path(
    table_root,
    "02_layer_status.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Combine validated per-layer aggregates
# ----------------------------------------------------------------------

layer_results <- lapply(
  layer_result_paths,
  readRDS
)

feature_identities <- lapply(
  layer_results,
  function(x) {
    rownames(
      x$primary_counts
    )
  }
)

if (
  !all(
    vapply(
      feature_identities,
      identical,
      logical(1),
      feature_identities[[1]]
    )
  )
) {
  stop(
    "Per-layer pseudobulk feature identities are inconsistent."
  )
}

primary_counts <- do.call(
  cbind,
  lapply(
    layer_results,
    function(x) {
      x$primary_counts
    }
  )
)

atlas_counts <- do.call(
  cbind,
  lapply(
    layer_results,
    function(x) {
      x$atlas_counts
    }
  )
)

primary_metadata <- do.call(
  rbind,
  lapply(
    layer_results,
    function(x) {
      x$primary_metadata
    }
  )
)

atlas_metadata <- do.call(
  rbind,
  lapply(
    layer_results,
    function(x) {
      x$atlas_metadata
    }
  )
)

VDR_cell_audit <- do.call(
  rbind,
  lapply(
    layer_results,
    function(x) {
      x$VDR_cell_audit
    }
  )
)

cell_audit <- do.call(
  rbind,
  lapply(
    layer_results,
    function(x) {
      x$cell_audit
    }
  )
)

raw_count_audit <- do.call(
  rbind,
  lapply(
    layer_results,
    function(x) {
      x$raw_count_audit
    }
  )
)

if (
  anyDuplicated(
    colnames(primary_counts)
  ) > 0L ||
    anyDuplicated(
      colnames(atlas_counts)
    ) > 0L ||
    !identical(
      colnames(primary_counts),
      primary_metadata$pb_id
    ) ||
    !identical(
      colnames(atlas_counts),
      atlas_metadata$pb_id
    )
) {
  stop(
    "Combined pseudobulk columns and metadata are not aligned."
  )
}

if (
  sum(primary_metadata$cell_n) !=
    expected_primary_cell_n ||
    sum(atlas_metadata$cell_n) !=
      expected_atlas_cell_n ||
    sum(
      cell_audit$
        cluster8_excluded_cell_n
    ) != expected_cluster8_cell_n ||
    ncol(primary_counts) != 54L ||
    ncol(atlas_counts) != 54L
) {
  stop(
    "Combined pseudobulk identity does not match 45,630 primary cells, 47,161 atlas cells and 54 profiles."
  )
}

saveRDS(
  primary_counts,
  primary_counts_path,
  compress = "gzip"
)

saveRDS(
  atlas_counts,
  atlas_counts_path,
  compress = "gzip"
)

write_csv_utf8(
  primary_metadata,
  file.path(
    table_root,
    "03_primary_pb_metadata.csv"
  )
)

write_csv_utf8(
  atlas_metadata,
  file.path(
    table_root,
    "04_atlas_pb_metadata.csv"
  )
)

write_csv_utf8(
  VDR_cell_audit,
  file.path(
    table_root,
    "05_VDR_cell_detection.csv"
  )
)

write_csv_utf8(
  cell_audit,
  file.path(
    table_root,
    "06_cell_aggregation_audit.csv"
  )
)

write_csv_utf8(
  raw_count_audit,
  file.path(
    table_root,
    "07_raw_count_conservation.csv"
  )
)

pb_status <- data.frame(
  pseudobulk_version =
    pseudobulk_version,
  status = "COMPLETE",
  feature_n =
    nrow(primary_counts),
  primary_profile_n =
    ncol(primary_counts),
  atlas_profile_n =
    ncol(atlas_counts),
  primary_cell_n =
    sum(primary_metadata$cell_n),
  atlas_cell_n =
    sum(atlas_metadata$cell_n),
  cluster8_excluded_cell_n =
    sum(
      cell_audit$
        cluster8_excluded_cell_n
    ),
  count_layer_n =
    length(count_layers),
  minimum_cells_primary =
    minimum_cells_primary,
  primary_counts_path =
    normalize_slashes(
      primary_counts_path
    ),
  primary_counts_md5 =
    unname(
      tools::md5sum(
        primary_counts_path
      )
    ),
  atlas_counts_path =
    normalize_slashes(
      atlas_counts_path
    ),
  atlas_counts_md5 =
    unname(
      tools::md5sum(
        atlas_counts_path
      )
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  pb_status,
  pb_status_path
)

rm(
  annotated_object,
  metadata,
  layer_results
)

gc()


# ----------------------------------------------------------------------
# 8. VDR expression and locked-program scores
# ----------------------------------------------------------------------

feature_names <- rownames(
  primary_counts
)

candidate_match <- match(
  toupper(
    candidate_table$feature_id
  ),
  toupper(
    feature_names
  )
)

candidate_table$resolved_feature <- ifelse(
  is.na(candidate_match),
  NA_character_,
  feature_names[
    candidate_match
  ]
)

candidate_table$feature_available <-
  !is.na(candidate_match)

VDR_feature_index <- match(
  "VDR",
  toupper(feature_names)
)

if (is.na(VDR_feature_index)) {
  stop(
    "VDR is absent from pseudobulk counts."
  )
}

VDR_feature <- feature_names[
  VDR_feature_index
]

program_definitions <- list(
  RNA_ALL_436 =
    candidate_table$feature_id,
  OCCUPANCY_SUPPORTED_284 =
    candidate_table$feature_id[
      candidate_table$
        occupancy_supported
    ],
  SHARED_OCCUPANCY_46 =
    candidate_table$feature_id[
      candidate_table$
        shared_occupancy_supported
    ],
  PROMOTER_OCCUPANCY_84 =
    candidate_table$feature_id[
      candidate_table$
        any_promoter_overlap_support
    ]
)

program_sizes <- vapply(
  program_definitions,
  length,
  integer(1)
)

if (
  !identical(
    unname(program_sizes),
    c(
      436L,
      284L,
      46L,
      84L
    )
  )
) {
  stop(
    "Locked program sizes are not 436/284/46/84."
  )
}

primary_library_size <- Matrix::colSums(
  primary_counts
)

atlas_library_size <- Matrix::colSums(
  atlas_counts
)

if (
  any(primary_library_size <= 0) ||
    any(atlas_library_size <= 0)
) {
  stop(
    "One or more pseudobulk libraries are empty."
  )
}

primary_VDR_count <- as.numeric(
  primary_counts[
    VDR_feature,
    ,
    drop = TRUE
  ]
)

atlas_VDR_count <- as.numeric(
  atlas_counts[
    VDR_feature,
    ,
    drop = TRUE
  ]
)

primary_VDR_expression <- cbind(
  primary_metadata,
  data.frame(
    VDR_feature =
      VDR_feature,
    VDR_raw_count =
      primary_VDR_count,
    VDR_CPM =
      primary_VDR_count /
      primary_library_size *
      1e6,
    VDR_log2_CPM_plus_1 =
      log2(
        primary_VDR_count /
          primary_library_size *
          1e6 +
          1
      ),
    stringsAsFactors = FALSE
  )
)

atlas_VDR_expression <- cbind(
  atlas_metadata,
  data.frame(
    VDR_feature =
      VDR_feature,
    VDR_raw_count =
      atlas_VDR_count,
    VDR_CPM =
      atlas_VDR_count /
      atlas_library_size *
      1e6,
    VDR_log2_CPM_plus_1 =
      log2(
        atlas_VDR_count /
          atlas_library_size *
          1e6 +
          1
      ),
    stringsAsFactors = FALSE
  )
)

available_candidate_table <- candidate_table[
  candidate_table$feature_available &
    is.finite(
      candidate_table$
        interaction_logFC
    ) &
    candidate_table$
      interaction_logFC != 0,
  ,
  drop = FALSE
]

available_features <- unique(
  available_candidate_table$
    resolved_feature
)

primary_candidate_counts <- as.matrix(
  primary_counts[
    available_features,
    ,
    drop = FALSE
  ]
)

atlas_candidate_counts <- as.matrix(
  atlas_counts[
    available_features,
    ,
    drop = FALSE
  ]
)

primary_candidate_logCPM <- log2(
  sweep(
    primary_candidate_counts,
    MARGIN = 2,
    STATS =
      primary_library_size /
      1e6,
    FUN = "/"
  ) + 1
)

atlas_candidate_logCPM <- log2(
  sweep(
    atlas_candidate_counts,
    MARGIN = 2,
    STATS =
      atlas_library_size /
      1e6,
    FUN = "/"
  ) + 1
)

reference_statistics <- row_zscore_reference(
  primary_candidate_logCPM
)

names(
  reference_statistics$mean
) <- rownames(
  primary_candidate_logCPM
)

names(
  reference_statistics$sd
) <- rownames(
  primary_candidate_logCPM
)

primary_candidate_z <- apply_row_reference_z(
  primary_candidate_logCPM,
  reference_mean =
    reference_statistics$mean,
  reference_sd =
    reference_statistics$sd
)

atlas_candidate_z <- apply_row_reference_z(
  atlas_candidate_logCPM,
  reference_mean =
    reference_statistics$mean,
  reference_sd =
    reference_statistics$sd
)

score_programs <- function(
  standardized_matrix,
  pb_metadata,
  matrix_role
) {
  score_rows <- list()
  coverage_rows <- list()

  for (
    program_name in names(
      program_definitions
    )
  ) {
    frozen_genes <- unique(
      program_definitions[[
        program_name
      ]]
    )

    program_table <- candidate_table[
      candidate_table$feature_id %in%
        frozen_genes &
        candidate_table$
          feature_available,
      ,
      drop = FALSE
    ]

    measurable_features <- unique(
      program_table$
        resolved_feature[
          program_table$
            resolved_feature %in%
            rownames(
              standardized_matrix
            )
        ]
    )

    feature_rows <- match(
      measurable_features,
      program_table$
        resolved_feature
    )

    program_table <- program_table[
      feature_rows,
      ,
      drop = FALSE
    ]

    valid <- (
      is.finite(
        program_table$
          interaction_logFC
      ) &
        program_table$
          interaction_logFC != 0 &
        measurable_features %in%
          rownames(
            standardized_matrix
          )
    )

    measurable_features <-
      measurable_features[valid]

    program_table <-
      program_table[valid, , drop = FALSE]

    weights <- sign(
      program_table$
        interaction_logFC
    )

    variable_gene <- is.finite(
      reference_statistics$sd[
        measurable_features
      ]
    ) &
      reference_statistics$sd[
        measurable_features
      ] > 0

    measurable_features <-
      measurable_features[
        variable_gene
      ]

    weights <- weights[
      variable_gene
    ]

    if (length(measurable_features) < 10L) {
      stop(
        "Program ",
        program_name,
        " has fewer than 10 measurable variable genes."
      )
    }

    signed_matrix <- sweep(
      standardized_matrix[
        measurable_features,
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

    score_metadata <- as.data.frame(
      pb_metadata,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # pb_metadata already contains matrix_role from pseudobulk
    # aggregation. Overwrite it explicitly instead of appending a
    # second same-named column.
    score_metadata$matrix_role <-
      matrix_role

    score_statistics <- data.frame(
      program_name =
        program_name,
      score =
        as.numeric(
          score_value
        ),
      frozen_gene_n =
        length(frozen_genes),
      measurable_gene_n =
        length(
          measurable_features
        ),
      coverage_fraction =
        length(
          measurable_features
        ) /
        length(frozen_genes),
      score_definition =
        "MEAN_SIGNED_GENE_ZSCORE_USING_PRIMARY_REFERENCE",
      interpretation =
        "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    overlapping_columns <- intersect(
      colnames(score_metadata),
      colnames(score_statistics)
    )

    if (length(overlapping_columns)) {
      stop(
        "Program-score metadata/statistic columns overlap: ",
        paste(
          overlapping_columns,
          collapse = ", "
        )
      )
    }

    score_row <- cbind(
      score_metadata,
      score_statistics
    )

    if (
      anyDuplicated(
        colnames(score_row)
      ) > 0L
    ) {
      stop(
        "Program-score table contains duplicate column names for ",
        program_name,
        "."
      )
    }

    score_rows[[
      program_name
    ]] <- score_row

    coverage_rows[[
      program_name
    ]] <- data.frame(
      program_name =
        program_name,
      frozen_gene_n =
        length(frozen_genes),
      available_gene_n =
        sum(
          candidate_table$
            feature_id %in%
            frozen_genes &
            candidate_table$
              feature_available
        ),
      measurable_variable_gene_n =
        length(
          measurable_features
        ),
      coverage_fraction =
        length(
          measurable_features
        ) /
        length(frozen_genes),
      missing_genes =
        collapse_unique(
          setdiff(
            frozen_genes,
            candidate_table$
              feature_id[
                candidate_table$
                  feature_available
              ]
          )
        ),
      stringsAsFactors = FALSE
    )
  }

  list(
    scores = as.data.frame(
      data.table::rbindlist(
        score_rows,
        use.names = TRUE,
        fill = TRUE
      )
    ),
    coverage = do.call(
      rbind,
      coverage_rows
    )
  )
}

primary_score_result <- score_programs(
  standardized_matrix =
    primary_candidate_z,
  pb_metadata =
    primary_metadata,
  matrix_role =
    "PRIMARY"
)

atlas_score_result <- score_programs(
  standardized_matrix =
    atlas_candidate_z,
  pb_metadata =
    atlas_metadata,
  matrix_role =
    "ATLAS_SENSITIVITY"
)

program_coverage <- primary_score_result$
  coverage

if (
  any(
    program_coverage$
      coverage_fraction < 0.90
  )
) {
  stop(
    "One or more locked VDR programs have below 90% coverage."
  )
}

primary_program_scores <- primary_score_result$
  scores

atlas_program_scores <- atlas_score_result$
  scores

score_column_audit <- data.frame(
  matrix_role = c(
    "PRIMARY",
    "ATLAS_SENSITIVITY"
  ),
  row_n = c(
    nrow(primary_program_scores),
    nrow(atlas_program_scores)
  ),
  column_n = c(
    ncol(primary_program_scores),
    ncol(atlas_program_scores)
  ),
  duplicate_column_n = c(
    sum(
      duplicated(
        colnames(
          primary_program_scores
        )
      )
    ),
    sum(
      duplicated(
        colnames(
          atlas_program_scores
        )
      )
    )
  ),
  column_names = c(
    paste(
      colnames(
        primary_program_scores
      ),
      collapse = ";"
    ),
    paste(
      colnames(
        atlas_program_scores
      ),
      collapse = ";"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  score_column_audit,
  file.path(
    table_root,
    "00_program_score_column_audit.csv"
  )
)

if (
  any(
    score_column_audit$
      duplicate_column_n > 0L
  ) ||
    !all(
      c(
        "matrix_role",
        "program_name",
        "score",
        "broad_cell_type",
        "sample_accession",
        "minimum_cell_threshold_pass"
      ) %in%
        colnames(
          primary_program_scores
        )
    )
) {
  stop(
    "Program-score column-name audit failed."
  )
}

write_csv_utf8(
  candidate_table,
  file.path(
    table_root,
    "08_candidate_feature_audit.csv"
  )
)

write_csv_utf8(
  program_coverage,
  file.path(
    table_root,
    "09_program_coverage.csv"
  )
)

write_csv_utf8(
  primary_VDR_expression,
  file.path(
    table_root,
    "10_primary_VDR_expression.csv"
  )
)

write_csv_utf8(
  atlas_VDR_expression,
  file.path(
    table_root,
    "11_atlas_VDR_expression.csv"
  )
)

write_csv_utf8(
  primary_program_scores,
  file.path(
    table_root,
    "12_primary_program_scores.csv"
  )
)

write_csv_utf8(
  atlas_program_scores,
  file.path(
    table_root,
    "13_atlas_program_scores.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Model eligibility and cluster-8 sensitivity
# ----------------------------------------------------------------------

eligibility_rows <- list()

for (
  cell_type in sort(
    unique(
      primary_metadata$
        broad_cell_type
    )
  )
) {
  cell_type_metadata <- primary_metadata[
    primary_metadata$
      broad_cell_type ==
      cell_type,
    ,
    drop = FALSE
  ]

  Control_pass_n <- sum(
    cell_type_metadata$
      analysis_group ==
      "CONTROL" &
      cell_type_metadata$
        minimum_cell_threshold_pass
  )

  RIF_pass_n <- sum(
    cell_type_metadata$
      analysis_group ==
      "RIF" &
      cell_type_metadata$
        minimum_cell_threshold_pass
  )

  eligibility_rows[[
    length(eligibility_rows) + 1L
  ]] <- data.frame(
    broad_cell_type =
      cell_type,
    total_profile_n =
      nrow(cell_type_metadata),
    Control_total_n =
      sum(
        cell_type_metadata$
          analysis_group ==
          "CONTROL"
      ),
    RIF_total_n =
      sum(
        cell_type_metadata$
          analysis_group ==
          "RIF"
      ),
    Control_threshold_pass_n =
      Control_pass_n,
    RIF_threshold_pass_n =
      RIF_pass_n,
    minimum_cell_n =
      min(
        cell_type_metadata$
          cell_n
      ),
    median_cell_n =
      median(
        cell_type_metadata$
          cell_n
      ),
    maximum_cell_n =
      max(
        cell_type_metadata$
          cell_n
      ),
    donor_level_model_eligible =
      Control_pass_n >= 3L &&
      RIF_pass_n >= 3L,
    descriptive_localization_eligible =
      sum(
        cell_type_metadata$
          minimum_cell_threshold_pass
      ) >= 2L,
    stringsAsFactors = FALSE
  )
}

model_eligibility <- do.call(
  rbind,
  eligibility_rows
)

model_eligible_types <- model_eligibility$
  broad_cell_type[
    model_eligibility$
      donor_level_model_eligible
  ]

if (length(model_eligible_types) != 5L) {
  stop(
    "Expected five donor-level model-eligible broad cell types, obtained ",
    length(model_eligible_types),
    ": ",
    collapse_unique(
      model_eligible_types
    )
  )
}

primary_unciliated_scores <- primary_program_scores[
  primary_program_scores$
    broad_cell_type ==
    "Unciliated epithelial cells",
  ,
  drop = FALSE
]

atlas_unciliated_scores <- atlas_program_scores[
  atlas_program_scores$
    broad_cell_type ==
    "Unciliated epithelial cells",
  ,
  drop = FALSE
]

cluster8_program_sensitivity <- merge(
  primary_unciliated_scores[
    ,
    c(
      "sample_accession",
      "donor_id",
      "analysis_group",
      "program_name",
      "cell_n",
      "score"
    ),
    drop = FALSE
  ],
  atlas_unciliated_scores[
    ,
    c(
      "sample_accession",
      "donor_id",
      "analysis_group",
      "program_name",
      "cell_n",
      "score"
    ),
    drop = FALSE
  ],
  by = c(
    "sample_accession",
    "donor_id",
    "analysis_group",
    "program_name"
  ),
  all = TRUE,
  suffixes = c(
    "_primary",
    "_atlas"
  ),
  sort = FALSE
)

cluster8_program_sensitivity$
  cluster8_cell_n <- (
    cluster8_program_sensitivity$
      cell_n_atlas -
      cluster8_program_sensitivity$
        cell_n_primary
  )

cluster8_program_sensitivity$
  score_delta_atlas_minus_primary <- (
    cluster8_program_sensitivity$
      score_atlas -
      cluster8_program_sensitivity$
        score_primary
  )

cluster8_VDR_sensitivity <- merge(
  primary_VDR_expression[
    primary_VDR_expression$
      broad_cell_type ==
      "Unciliated epithelial cells",
    c(
      "sample_accession",
      "donor_id",
      "analysis_group",
      "cell_n",
      "VDR_log2_CPM_plus_1"
    ),
    drop = FALSE
  ],
  atlas_VDR_expression[
    atlas_VDR_expression$
      broad_cell_type ==
      "Unciliated epithelial cells",
    c(
      "sample_accession",
      "donor_id",
      "analysis_group",
      "cell_n",
      "VDR_log2_CPM_plus_1"
    ),
    drop = FALSE
  ],
  by = c(
    "sample_accession",
    "donor_id",
    "analysis_group"
  ),
  all = TRUE,
  suffixes = c(
    "_primary",
    "_atlas"
  ),
  sort = FALSE
)

cluster8_VDR_sensitivity$
  cluster8_cell_n <- (
    cluster8_VDR_sensitivity$
      cell_n_atlas -
      cluster8_VDR_sensitivity$
        cell_n_primary
  )

cluster8_VDR_sensitivity$
  VDR_delta_atlas_minus_primary <- (
    cluster8_VDR_sensitivity$
      VDR_log2_CPM_plus_1_atlas -
      cluster8_VDR_sensitivity$
        VDR_log2_CPM_plus_1_primary
  )

localization_summary <- data.table::as.data.table(
  primary_program_scores[
    primary_program_scores$
      minimum_cell_threshold_pass,
    ,
    drop = FALSE
  ]
)[
  ,
  .(
    donor_n =
      data.table::uniqueN(
        donor_id
      ),
    median_score =
      median(score),
    minimum_score =
      min(score),
    maximum_score =
      max(score)
  ),
  by = .(
    broad_cell_type,
    program_name
  )
]

localization_summary <- as.data.frame(
  localization_summary
)

VDR_localization_summary <- data.table::as.data.table(
  primary_VDR_expression[
    primary_VDR_expression$
      minimum_cell_threshold_pass,
    ,
    drop = FALSE
  ]
)[
  ,
  .(
    donor_n =
      data.table::uniqueN(
        donor_id
      ),
    median_VDR_log2_CPM_plus_1 =
      median(
        VDR_log2_CPM_plus_1
      ),
    minimum_VDR_log2_CPM_plus_1 =
      min(
        VDR_log2_CPM_plus_1
      ),
    maximum_VDR_log2_CPM_plus_1 =
      max(
        VDR_log2_CPM_plus_1
      ),
    median_VDR_positive_fraction =
      median(
        VDR_cell_audit$
          primary_VDR_positive_fraction[
            match(
              paste(
                donor_id,
                broad_cell_type
              ),
              paste(
                VDR_cell_audit$donor_id,
                VDR_cell_audit$
                  broad_cell_type
              )
            )
          ],
        na.rm = TRUE
      )
  ),
  by = broad_cell_type
]

VDR_localization_summary <- as.data.frame(
  VDR_localization_summary
)

write_csv_utf8(
  model_eligibility,
  file.path(
    table_root,
    "14_model_eligibility.csv"
  )
)

write_csv_utf8(
  cluster8_program_sensitivity,
  file.path(
    table_root,
    "15_cluster8_program_sensitivity.csv"
  )
)

write_csv_utf8(
  cluster8_VDR_sensitivity,
  file.path(
    table_root,
    "16_cluster8_VDR_sensitivity.csv"
  )
)

write_csv_utf8(
  localization_summary,
  file.path(
    table_root,
    "17_program_localization_summary.csv"
  )
)

write_csv_utf8(
  VDR_localization_summary,
  file.path(
    table_root,
    "18_VDR_localization_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 10. QC figures
# ----------------------------------------------------------------------

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

cell_count_plot_data <- primary_metadata

cell_count_plot_data$
  broad_cell_type <- factor(
    cell_count_plot_data$
      broad_cell_type,
    levels = sort(
      unique(
        cell_count_plot_data$
          broad_cell_type
      )
    )
  )

plot_cell_counts <- ggplot2::ggplot(
  cell_count_plot_data,
  ggplot2::aes(
    x = broad_cell_type,
    y = cell_n,
    shape = analysis_group
  )
) +
  ggplot2::geom_hline(
    yintercept =
      minimum_cells_primary,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.12,
        height = 0
      )
  ) +
  ggplot2::scale_y_log10(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "GSE183837 primary pseudobulk cell counts",
    subtitle = "Cluster 8 excluded; dashed line denotes the 20-cell model threshold",
    x = "Broad cell type",
    y = "Cells per donor pseudobulk (log scale)",
    shape = "Group"
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
  plot_cell_counts,
  "QC1_primary_pseudobulk_cell_counts",
  width = 10,
  height = 6
)

library_plot_data <- primary_metadata

plot_library <- ggplot2::ggplot(
  library_plot_data,
  ggplot2::aes(
    x = broad_cell_type,
    y = library_size,
    shape = analysis_group
  )
) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.12,
        height = 0
      )
  ) +
  ggplot2::scale_y_log10(
    labels = scales::comma
  ) +
  ggplot2::labs(
    title = "Primary pseudobulk library sizes",
    subtitle = "Raw counts aggregated by donor and locked broad cell type",
    x = "Broad cell type",
    y = "Library size (log scale)",
    shape = "Group"
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
  plot_library,
  "QC2_primary_library_sizes",
  width = 10,
  height = 6
)

VDR_plot_data <- primary_VDR_expression[
  primary_VDR_expression$
    minimum_cell_threshold_pass,
  ,
  drop = FALSE
]

plot_VDR <- ggplot2::ggplot(
  VDR_plot_data,
  ggplot2::aes(
    x = broad_cell_type,
    y = VDR_log2_CPM_plus_1,
    shape = analysis_group
  )
) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.12,
        height = 0
      )
  ) +
  ggplot2::labs(
    title = "Descriptive VDR expression localization",
    subtitle = "Threshold-passing primary pseudobulks; no group test in this step",
    x = "Broad cell type",
    y = "VDR log2(CPM + 1)",
    shape = "Group"
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
  plot_VDR,
  "QC3_VDR_localization",
  width = 10,
  height = 6
)

score_heatmap_data <- primary_program_scores[
  primary_program_scores$
    minimum_cell_threshold_pass,
  ,
  drop = FALSE
]

if (
  anyDuplicated(
    colnames(score_heatmap_data)
  ) > 0L
) {
  stop(
    "score_heatmap_data contains duplicate column names before plotting."
  )
}

score_heatmap_data$
  profile_label <- paste(
    score_heatmap_data$
      sample_accession,
    score_heatmap_data$
      broad_cell_type,
    sep = " | "
  )

score_heatmap_data$
  profile_label <- factor(
    score_heatmap_data$
      profile_label,
    levels = unique(
      score_heatmap_data$
        profile_label
    )
  )

score_heatmap_data$
  program_name <- factor(
    score_heatmap_data$
      program_name,
    levels = c(
      "RNA_ALL_436",
      "OCCUPANCY_SUPPORTED_284",
      "SHARED_OCCUPANCY_46",
      "PROMOTER_OCCUPANCY_84"
    )
  )

plot_score_heatmap <- ggplot2::ggplot(
  score_heatmap_data,
  ggplot2::aes(
    x = program_name,
    y = profile_label,
    fill = score
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2() +
  ggplot2::labs(
    title = "Locked VDR-program localization across donor pseudobulks",
    subtitle = "Scores use the primary matrix and signed interaction direction",
    x = "VDR program",
    y = "Donor | broad cell type",
    fill = "Signed score"
  ) +
  theme_QC +
  ggplot2::theme(
    axis.text.y =
      ggplot2::element_text(
        size = 5
      ),
    axis.text.x =
      ggplot2::element_text(
        angle = 35,
        hjust = 1
      )
  )

save_plot_pair(
  plot_score_heatmap,
  "QC4_program_score_heatmap",
  width = 10,
  height = 14
)

localization_plot_data <- localization_summary

localization_plot_data$
  broad_cell_type <- factor(
    localization_plot_data$
      broad_cell_type,
    levels = sort(
      unique(
        localization_plot_data$
          broad_cell_type
      )
    )
  )

plot_localization <- ggplot2::ggplot(
  localization_plot_data,
  ggplot2::aes(
    x = program_name,
    y = broad_cell_type,
    fill = median_score
  )
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2() +
  ggplot2::labs(
    title = "Median VDR-program alignment by broad cell type",
    subtitle = "Descriptive localization among threshold-passing donor profiles",
    x = "VDR program",
    y = "Broad cell type",
    fill = "Median score"
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
  plot_localization,
  "QC5_program_localization_summary",
  width = 9,
  height = 6
)

eligibility_plot_data <- model_eligibility

eligibility_long <- data.frame(
  broad_cell_type = rep(
    eligibility_plot_data$
      broad_cell_type,
    each = 2
  ),
  analysis_group = rep(
    c(
      "CONTROL",
      "RIF"
    ),
    times =
      nrow(
        eligibility_plot_data
      )
  ),
  threshold_pass_donor_n = as.vector(
    t(
      as.matrix(
        eligibility_plot_data[
          ,
          c(
            "Control_threshold_pass_n",
            "RIF_threshold_pass_n"
          )
        ]
      )
    )
  ),
  donor_level_model_eligible = rep(
    eligibility_plot_data$
      donor_level_model_eligible,
    each = 2
  ),
  stringsAsFactors = FALSE
)

plot_eligibility <- ggplot2::ggplot(
  eligibility_long,
  ggplot2::aes(
    x = analysis_group,
    y = broad_cell_type,
    fill = threshold_pass_donor_n
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(
      label =
        threshold_pass_donor_n
    )
  ) +
  ggplot2::labs(
    title = "Donor-level model eligibility",
    subtitle = "Primary model requires at least 3 threshold-passing donors in each group",
    x = "Group",
    y = "Broad cell type",
    fill = "Donors"
  ) +
  theme_QC

save_plot_pair(
  plot_eligibility,
  "QC6_model_eligibility",
  width = 7,
  height = 5.5
)

cluster8_plot_data <- cluster8_program_sensitivity

plot_cluster8 <- ggplot2::ggplot(
  cluster8_plot_data,
  ggplot2::aes(
    x = program_name,
    y =
      score_delta_atlas_minus_primary,
    shape = analysis_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.12,
        height = 0
      )
  ) +
  ggplot2::labs(
    title = "Cluster-8 sensitivity effect on unciliated epithelial scores",
    subtitle = "Atlas score minus primary score; sensitivity only",
    x = "VDR program",
    y = "Score difference",
    shape = "Group"
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
  plot_cluster8,
  "QC7_cluster8_program_sensitivity",
  width = 9,
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
# 11. Gate
# ----------------------------------------------------------------------

gate_checks <- data.frame(
  check = c(
    "upstream_annotation_gate_authorizes_pseudobulk",
    "annotated_checkpoint_has_47161_cells",
    "nine_count_layers_processed",
    "all_layer_jobs_complete",
    "primary_matrix_has_45630_cells",
    "atlas_matrix_has_47161_cells",
    "cluster8_exclusion_is_1531_cells",
    "primary_matrix_has_54_profiles",
    "atlas_matrix_has_54_profiles",
    "raw_count_conservation_passes_all_layers",
    "primary_count_conservation_passes_all_layers",
    "primary_metadata_aligned_to_counts",
    "atlas_metadata_aligned_to_counts",
    "donor_design_is_3_Control_6_RIF",
    "locked_program_sizes_are_436_284_46_84",
    "all_program_coverage_at_least_90_percent",
    "VDR_is_measurable",
    "five_broad_types_are_donor_model_eligible",
    "cluster8_is_sensitivity_only",
    "primary_and_atlas_count_checkpoints_exist",
    "donor_level_VDR_program_models_can_be_authorized",
    "disease_testing_not_performed_in_this_step",
    "free_disk_after_run_at_least_8GB_or_unknown"
  ),
  pass = c(
    upstream_gate$gate_status[[1]] ==
      "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION",
    expected_atlas_cell_n ==
      sum(atlas_metadata$cell_n),
    length(count_layers) == 9L,
    nrow(layer_status_all) == 9L &&
      all(
        layer_status_all$status ==
          "COMPLETE"
      ),
    sum(primary_metadata$cell_n) ==
      expected_primary_cell_n,
    sum(atlas_metadata$cell_n) ==
      expected_atlas_cell_n,
    sum(
      cell_audit$
        cluster8_excluded_cell_n
    ) == expected_cluster8_cell_n,
    ncol(primary_counts) == 54L,
    ncol(atlas_counts) == 54L,
    all(
      raw_count_audit$
        raw_count_sum ==
        raw_count_audit$
          atlas_aggregated_count_sum
    ),
    all(
      raw_count_audit$
        primary_input_count_sum ==
        raw_count_audit$
          primary_aggregated_count_sum
    ),
    identical(
      colnames(primary_counts),
      primary_metadata$pb_id
    ),
    identical(
      colnames(atlas_counts),
      atlas_metadata$pb_id
    ),
    nrow(donor_group_identity) == 9L &&
      sum(
        donor_group_identity$
          analysis_group ==
          "CONTROL"
      ) == 3L &&
      sum(
        donor_group_identity$
          analysis_group ==
          "RIF"
      ) == 6L,
    identical(
      unname(program_sizes),
      c(
        436L,
        284L,
        46L,
        84L
      )
    ),
    all(
      program_coverage$
        coverage_fraction >= 0.90
    ),
    VDR_feature %in%
      rownames(primary_counts),
    length(model_eligible_types) ==
      5L,
    all(
      primary_metadata$
        cell_n[
          primary_metadata$
            broad_cell_type ==
            "Unciliated epithelial cells"
        ] <=
        atlas_metadata$
          cell_n[
            match(
              paste(
                primary_metadata$
                  sample_accession[
                    primary_metadata$
                      broad_cell_type ==
                      "Unciliated epithelial cells"
                  ],
                "Unciliated epithelial cells"
              ),
              paste(
                atlas_metadata$
                  sample_accession,
                atlas_metadata$
                  broad_cell_type
              )
            )
          ]
    ),
    file.exists(
      primary_counts_path
    ) &&
      file.exists(
        atlas_counts_path
      ),
    TRUE,
    TRUE,
    is.na(
      get_free_disk_gb()
    ) ||
      get_free_disk_gb() >= 8
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
  "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED"
} else {
  "HOLD_GSE183837_PSEUDOBULK_OR_LOCALIZATION_FAILED"
}

review_flags <- c(
  "CLUSTER8_SENSITIVITY_ONLY",
  "CILIATED_EPITHELIAL_DESCRIPTIVE_LOW_CELL_COUNTS",
  "CONTROL_N3_RIF_N6_LOW_POWER_DONOR_MODELS",
  "NO_CELL_LEVEL_DISEASE_INFERENCE"
)

pb_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  annotation_version =
    annotation_version,
  pseudobulk_version =
    pseudobulk_version,
  gate_status = gate_status,
  feature_n =
    nrow(primary_counts),
  primary_profile_n =
    ncol(primary_counts),
  atlas_profile_n =
    ncol(atlas_counts),
  primary_cell_n =
    sum(primary_metadata$cell_n),
  atlas_cell_n =
    sum(atlas_metadata$cell_n),
  cluster8_excluded_cell_n =
    sum(
      cell_audit$
        cluster8_excluded_cell_n
    ),
  minimum_cells_primary =
    minimum_cells_primary,
  donor_model_eligible_type_n =
    length(model_eligible_types),
  donor_model_eligible_types =
    collapse_unique(
      model_eligible_types
    ),
  descriptive_only_types =
    collapse_unique(
      model_eligibility$
        broad_cell_type[
          !model_eligibility$
            donor_level_model_eligible
        ]
    ),
  review_flags =
    collapse_unique(
      review_flags
    ),
  donor_level_VDR_program_models_authorized =
    all_gate_checks_pass,
  gene_level_pseudobulk_DE_authorized =
    FALSE,
  cell_level_disease_testing_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  primary_counts_checkpoint =
    normalize_slashes(
      primary_counts_path
    ),
  primary_counts_md5 =
    unname(
      tools::md5sum(
        primary_counts_path
      )
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  pb_gate,
  file.path(
    run_qc_root,
    "G183_PBVDR_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1N GSE183837 Pseudobulk and VDR Localization Card",
  "",
  paste0("- Gate: ", gate_status),
  paste0(
    "- Primary cells: ",
    format(
      sum(primary_metadata$cell_n),
      big.mark = ","
    )
  ),
  paste0(
    "- Atlas cells: ",
    format(
      sum(atlas_metadata$cell_n),
      big.mark = ","
    )
  ),
  paste0(
    "- Cluster-8 cells excluded from primary: ",
    format(
      sum(
        cell_audit$
          cluster8_excluded_cell_n
      ),
      big.mark = ","
    )
  ),
  paste0(
    "- Primary pseudobulk profiles: ",
    ncol(primary_counts)
  ),
  paste0(
    "- Atlas sensitivity profiles: ",
    ncol(atlas_counts)
  ),
  paste0(
    "- Model-eligible broad cell types: ",
    collapse_unique(
      model_eligible_types
    )
  ),
  paste0(
    "- Descriptive-only broad cell types: ",
    collapse_unique(
      model_eligibility$
        broad_cell_type[
          !model_eligibility$
            donor_level_model_eligible
        ]
    )
  ),
  paste0(
    "- Donor-level VDR/program models authorized: ",
    all_gate_checks_pass
  ),
  "- Gene-level pseudobulk DE authorized: FALSE",
  "- Cell-level disease testing authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "G183_PBVDR_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 pseudobulk and VDR-localization interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Primary matrix",
  "",
  "- Uses donor x locked broad cell type raw-count aggregation.",
  "- Excludes all 1,531 low-confidence cluster-8 cells.",
  "- Retains 45,630 cells across 54 donor x broad-type profiles.",
  "",
  "## Atlas sensitivity matrix",
  "",
  "- Retains all 47,161 cells.",
  "- Used only to quantify how cluster 8 changes unciliated epithelial VDR and program scores.",
  "",
  "## Model authorization",
  "",
  paste0(
    "- Donor-level model-eligible broad cell types: ",
    collapse_unique(
      model_eligible_types
    ),
    "."
  ),
  "- Ciliated epithelial cells remain descriptive because too few donor profiles pass the 20-cell threshold.",
  "- Control n=3 and RIF n=6 imply low power; exact permutation and effect-size reporting are required in the next step.",
  "",
  "## Interpretation boundary",
  "",
  "- Program scores are relative alignment scores, not direct VDR target activity.",
  "- No group difference has been tested here.",
  "- No cell-level disease inference is allowed.",
  "- No result is locked."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "G183_PBVDR_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 12. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_execution =
    "V1_0_DUPLICATE_MATRIX_ROLE_COLUMN_BLOCKED_GGPLOT",
  score_table_column_policy =
    "OVERWRITE_EXISTING_MATRIX_ROLE_AND_REQUIRE_UNIQUE_COLUMN_NAMES",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  annotation_version =
    annotation_version,
  pseudobulk_version =
    pseudobulk_version,
  upstream_annotation_gate =
    upstream_gate$gate_status[[1]],
  minimum_cells_primary =
    minimum_cells_primary,
  expected_primary_cell_n =
    expected_primary_cell_n,
  expected_atlas_cell_n =
    expected_atlas_cell_n,
  cluster8_policy =
    "PRIMARY_EXCLUDED_ATLAS_SENSITIVITY_RETAINED",
  program_score_definition =
    "MEAN_SIGNED_GENE_ZSCORE_USING_PRIMARY_REFERENCE",
  program_interpretation =
    "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION",
  donor_model_eligibility =
    "AT_LEAST_3_CONTROL_AND_3_RIF_DONORS_WITH_AT_LEAST_20_CELLS",
  donor_level_VDR_program_models_authorized =
    all_gate_checks_pass,
  gene_level_pseudobulk_DE_authorized =
    FALSE,
  cell_level_disease_testing_authorized =
    FALSE,
  result_lock_authorized =
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
  supersedes_partial_execution =
    "V1_0_DUPLICATE_MATRIX_ROLE_COLUMN_BLOCKED_GGPLOT",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  annotation_version =
    annotation_version,
  pseudobulk_version =
    pseudobulk_version,
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
  primary_cell_n =
    sum(primary_metadata$cell_n),
  atlas_cell_n =
    sum(atlas_metadata$cell_n),
  primary_profile_n =
    ncol(primary_counts),
  atlas_profile_n =
    ncol(atlas_counts),
  donor_model_eligible_type_n =
    length(model_eligible_types),
  donor_model_eligible_types =
    collapse_unique(
      model_eligible_types
    ),
  donor_level_models_authorized =
    all_gate_checks_pass,
  result_lock_authorized =
    FALSE,
  free_memory_before_GB =
    free_memory_before_GB,
  free_disk_before_GB =
    free_disk_before_GB,
  free_disk_after_GB =
    get_free_disk_gb(),
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

methods_lines <- c(
  "# STEP04E1N Methods shell",
  "",
  "Raw counts were aggregated within each donor and locked broad cell type from the nine split Seurat v5 count layers. The primary matrix excluded the low-confidence cluster 8, whereas an atlas sensitivity matrix retained all cells. Profiles with at least 20 cells were considered threshold-passing. VDR expression was summarized as log2(CPM + 1). Four prospectively locked VDR interaction programs were scored using mean signed gene-wise z-scores, with signs determined by the calcitriol-dependent interaction log-fold changes. Program scores were interpreted as relative alignment with the locked interaction direction rather than evidence of direct VDR binding or causal activation."
)

results_lines <- c(
  "# STEP04E1N Results shell",
  "",
  paste0(
    "The primary donor-by-cell-type pseudobulk matrix retained ",
    format(
      sum(primary_metadata$cell_n),
      big.mark = ","
    ),
    " cells in ",
    ncol(primary_counts),
    " profiles after excluding ",
    format(
      sum(
        cell_audit$
          cluster8_excluded_cell_n
      ),
      big.mark = ","
    ),
    " low-confidence cluster-8 cells."
  ),
  paste0(
    "Five broad cell types met the prespecified donor-level model criterion: ",
    collapse_unique(
      model_eligible_types
    ),
    "."
  ),
  "Ciliated epithelial cells were retained for descriptive localization but did not satisfy the donor-level modelling threshold."
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

small_files <- unique(c(
  list.files(
    table_root,
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
  pb_status_path,
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
  "primary_profiles=",
  ncol(primary_counts),
  "eligible_types=",
  length(model_eligible_types)
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
cat("STEP04E1N GSE183837 pseudobulk/VDR localization completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Primary cells: ",
  format(
    sum(primary_metadata$cell_n),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Atlas cells: ",
  format(
    sum(atlas_metadata$cell_n),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Cluster-8 cells excluded: ",
  format(
    sum(
      cell_audit$
        cluster8_excluded_cell_n
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)
cat(
  "Primary pseudobulk profiles: ",
  ncol(primary_counts),
  "\n",
  sep = ""
)
cat(
  "Model-eligible broad cell types: ",
  length(model_eligible_types),
  "\n",
  sep = ""
)
cat(
  "Eligible types: ",
  collapse_unique(
    model_eligible_types
  ),
  "\n",
  sep = ""
)
cat(
  "Donor-level VDR/program models authorized: ",
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
    "\nPseudobulk outputs are complete in the project folders.\n",
    sep = ""
  )
}

# ==============================================================================
# Analysis stage 7: 69_step04E1O_GSE183837_donor_models_v1_1_type_safe_summary
# ==============================================================================

# ======================================================================
# 69_step04E1O_GSE183837_donor_models_v1_1_type_safe_summary.R
# Nutrition_Endometrial_VDR_v2
#
# PURPOSE
# - Run donor-level Control-versus-RIF models in the five broad cell
#   types authorized by STEP04E1N.
# - Model VDR expression and four prospectively locked VDR programs.
# - Use exact label-permutation tests as the primary inferential method.
# - Report Welch estimates, confidence intervals, Hedges g,
#   rank-biserial effects and leave-one-donor-out stability.
# - Compare primary and cluster-8-retaining atlas results for
#   unciliated epithelial cells.
#
# PRIMARY MODELS
# Five broad cell types x five endpoints = 25 planned models:
# - VDR_expression
# - RNA_ALL_436
# - OCCUPANCY_SUPPORTED_284
# - SHARED_OCCUPANCY_46
# - PROMOTER_OCCUPANCY_84
#
# DESIGN
# - Endothelial, Myeloid, NK/T and Stromal: 3 Control vs 6 RIF,
#   84 exact allocations.
# - Unciliated epithelial: 3 Control vs 4 RIF,
#   35 exact allocations.
# - Ciliated epithelial: descriptive only.
#
# MULTIPLICITY
# - Primary multiplicity: BH across all 25 planned models.
# - Secondary: BH within endpoint and within broad-cell-type families.
# - Nominal exact P < 0.05 without BH significance is exploratory only.
#
# INTERPRETATION
# - Positive effect = RIF greater than Control.
# - Negative VDR-program effect = reduced alignment in RIF with the
#   locked calcitriol-dependent interaction direction.
# - Program scores are not direct VDR activity, binding or causality.
# - Donor is the biological replicate.
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
dataset_id <- "GSE183837"

script_name <-
  "69_step04E1O_GSE183837_donor_models_v1_1_type_safe_summary.R"

script_version <-
  "STEP04E1O_GSE183837_DONOR_MODELS_v1.1_TYPE_SAFE_SUMMARY"

integration_version <- "SCI183_v1_0"
pseudobulk_version <- "PB183_v1_0"
model_version <- "MD183_v1_0"

run_start <- Sys.time()

run_id <- paste0(
  "STEP04E1O_",
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

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    mean(x)
  }
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    median(x)
  }
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2L) {
    NA_real_
  } else {
    stats::sd(x)
  }
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    min(x)
  }
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (!length(x)) {
    NA_real_
  } else {
    max(x)
  }
}

hedges_g <- function(
  Control_values,
  RIF_values
) {
  Control_values <- as.numeric(
    Control_values
  )

  RIF_values <- as.numeric(
    RIF_values
  )

  n_Control <- length(
    Control_values
  )

  n_RIF <- length(
    RIF_values
  )

  if (
    n_Control < 2L ||
      n_RIF < 2L
  ) {
    return(NA_real_)
  }

  sd_Control <- stats::sd(
    Control_values
  )

  sd_RIF <- stats::sd(
    RIF_values
  )

  df <- n_Control +
    n_RIF -
    2L

  pooled_variance <- (
    (
      n_Control - 1L
    ) *
      sd_Control^2 +
      (
        n_RIF - 1L
      ) *
      sd_RIF^2
  ) / df

  if (
    !is.finite(
      pooled_variance
    ) ||
      pooled_variance <= 0
  ) {
    return(NA_real_)
  }

  Cohen_d <- (
    mean(RIF_values) -
      mean(Control_values)
  ) / sqrt(pooled_variance)

  correction_J <- 1 -
    3 / (
      4 * df - 1
    )

  correction_J * Cohen_d
}

rank_biserial_RIF_minus_Control <- function(
  Control_values,
  RIF_values
) {
  Control_values <- as.numeric(
    Control_values
  )

  RIF_values <- as.numeric(
    RIF_values
  )

  n_Control <- length(
    Control_values
  )

  n_RIF <- length(
    RIF_values
  )

  if (
    n_Control < 1L ||
      n_RIF < 1L
  ) {
    return(NA_real_)
  }

  pooled_values <- c(
    RIF_values,
    Control_values
  )

  pooled_ranks <- rank(
    pooled_values,
    ties.method = "average"
  )

  RIF_rank_sum <- sum(
    pooled_ranks[
      seq_len(n_RIF)
    ]
  )

  U_RIF <- RIF_rank_sum -
    n_RIF *
    (
      n_RIF + 1
    ) / 2

  2 * U_RIF /
    (
      n_RIF *
        n_Control
    ) -
    1
}

exact_permutation_mean_difference <- function(
  values,
  groups,
  model_id,
  tolerance = 1e-12
) {
  values <- as.numeric(values)
  groups <- as.character(groups)

  if (
    length(values) !=
      length(groups) ||
      any(
        !is.finite(values)
      ) ||
      any(
        !groups %in%
          c(
            "CONTROL",
            "RIF"
          )
      )
  ) {
    stop(
      "Invalid values or groups for exact permutation model ",
      model_id,
      "."
    )
  }

  Control_index <- which(
    groups == "CONTROL"
  )

  n_total <- length(values)
  n_Control <- length(
    Control_index
  )

  n_RIF <- n_total -
    n_Control

  if (
    n_Control < 2L ||
      n_RIF < 2L
  ) {
    stop(
      "Exact permutation requires at least two donors per group for ",
      model_id,
      "."
    )
  }

  observed_difference <- mean(
    values[
      groups == "RIF"
    ]
  ) -
    mean(
      values[
        groups == "CONTROL"
      ]
    )

  Control_allocations <- utils::combn(
    n_total,
    n_Control
  )

  permutation_differences <- apply(
    Control_allocations,
    2,
    function(Control_allocation) {
      allocation_is_Control <- rep(
        FALSE,
        n_total
      )

      allocation_is_Control[
        Control_allocation
      ] <- TRUE

      mean(
        values[
          !allocation_is_Control
        ]
      ) -
        mean(
          values[
            allocation_is_Control
          ]
        )
    }
  )

  exact_P <- mean(
    abs(
      permutation_differences
    ) >=
      abs(
        observed_difference
      ) -
      tolerance
  )

  permutation_table <- data.frame(
    model_id = model_id,
    permutation_index =
      seq_along(
        permutation_differences
      ),
    permuted_mean_difference_RIF_minus_Control =
      as.numeric(
        permutation_differences
      ),
    absolute_permuted_difference =
      abs(
        permutation_differences
      ),
    at_least_as_extreme_as_observed =
      abs(
        permutation_differences
      ) >=
      abs(
        observed_difference
      ) -
      tolerance,
    stringsAsFactors = FALSE
  )

  list(
    observed_difference =
      observed_difference,
    exact_P =
      exact_P,
    permutation_n =
      length(
        permutation_differences
      ),
    minimum_resolvable_P =
      1 /
      length(
        permutation_differences
      ),
    permutation_table =
      permutation_table
  )
}

leave_one_donor_out <- function(
  model_data,
  full_difference,
  model_id
) {
  output_rows <- list()

  donors <- as.character(
    model_data$donor_id
  )

  for (
    donor_index in seq_along(
      donors
    )
  ) {
    excluded_donor <- donors[
      donor_index
    ]

    retained_data <- model_data[
      model_data$donor_id !=
        excluded_donor,
      ,
      drop = FALSE
    ]

    Control_n <- sum(
      retained_data$
        analysis_group ==
        "CONTROL"
    )

    RIF_n <- sum(
      retained_data$
        analysis_group ==
        "RIF"
    )

    if (
      Control_n < 2L ||
        RIF_n < 2L
    ) {
      LOO_difference <- NA_real_
    } else {
      LOO_difference <- mean(
        retained_data$value[
          retained_data$
            analysis_group ==
            "RIF"
        ]
      ) -
        mean(
          retained_data$value[
            retained_data$
              analysis_group ==
              "CONTROL"
          ]
        )
    }

    direction_concordant <- if (
      is.finite(
        LOO_difference
      ) &&
        is.finite(
          full_difference
        ) &&
        full_difference != 0
    ) {
      sign(
        LOO_difference
      ) ==
        sign(
          full_difference
        )
    } else {
      NA
    }

    output_rows[[
      donor_index
    ]] <- data.frame(
      model_id = model_id,
      excluded_donor =
        excluded_donor,
      excluded_group =
        model_data$
          analysis_group[
            donor_index
          ],
      retained_Control_n =
        Control_n,
      retained_RIF_n =
        RIF_n,
      LOO_mean_difference_RIF_minus_Control =
        LOO_difference,
      direction_concordant_with_full =
        direction_concordant,
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    output_rows
  )
}

fit_donor_model <- function(
  model_data,
  broad_cell_type,
  endpoint,
  endpoint_class,
  matrix_role
) {
  required_columns <- c(
    "donor_id",
    "sample_accession",
    "analysis_group",
    "broad_cell_type",
    "endpoint",
    "value",
    "cell_n"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(model_data)
  )

  if (length(missing_columns)) {
    stop(
      "Model data are missing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  model_data <- model_data[
    is.finite(
      model_data$value
    ),
    ,
    drop = FALSE
  ]

  model_data <- model_data[
    order(
      model_data$
        analysis_group,
      model_data$donor_id
    ),
    ,
    drop = FALSE
  ]

  if (
    anyDuplicated(
      model_data$donor_id
    ) > 0L
  ) {
    stop(
      "Duplicate donor rows in ",
      broad_cell_type,
      " / ",
      endpoint,
      " / ",
      matrix_role,
      "."
    )
  }

  Control_values <- model_data$value[
    model_data$analysis_group ==
      "CONTROL"
  ]

  RIF_values <- model_data$value[
    model_data$analysis_group ==
      "RIF"
  ]

  n_Control <- length(
    Control_values
  )

  n_RIF <- length(
    RIF_values
  )

  model_id <- paste(
    matrix_role,
    broad_cell_type,
    endpoint,
    sep = "__"
  )

  if (
    n_Control < 3L ||
      n_RIF < 3L
  ) {
    stop(
      "Model is not authorized by the 3-per-group rule: ",
      model_id,
      " (",
      n_Control,
      " Control, ",
      n_RIF,
      " RIF)."
    )
  }

  constant_endpoint <- (
    stats::sd(
      model_data$value
    ) == 0
  )

  permutation_result <-
    exact_permutation_mean_difference(
      values =
        model_data$value,
      groups =
        model_data$
          analysis_group,
      model_id =
        model_id
    )

  if (constant_endpoint) {
    Welch_P <- NA_real_
    Welch_statistic <- NA_real_
    Welch_df <- NA_real_
    CI_lower <- 0
    CI_upper <- 0
  } else {
    Welch_result <- tryCatch(
      stats::t.test(
        x = RIF_values,
        y = Control_values,
        alternative = "two.sided",
        var.equal = FALSE,
        conf.level = 0.95
      ),
      error = function(e) NULL
    )

    if (is.null(Welch_result)) {
      Welch_P <- NA_real_
      Welch_statistic <- NA_real_
      Welch_df <- NA_real_
      CI_lower <- NA_real_
      CI_upper <- NA_real_
    } else {
      Welch_P <- as.numeric(
        Welch_result$p.value
      )

      Welch_statistic <- as.numeric(
        Welch_result$statistic
      )

      Welch_df <- as.numeric(
        Welch_result$parameter
      )

      CI_lower <- as.numeric(
        Welch_result$conf.int[1]
      )

      CI_upper <- as.numeric(
        Welch_result$conf.int[2]
      )
    }
  }

  LOO_table <- leave_one_donor_out(
    model_data =
      model_data,
    full_difference =
      permutation_result$
        observed_difference,
    model_id =
      model_id
  )

  LOO_valid <- LOO_table[
    is.finite(
      LOO_table$
        LOO_mean_difference_RIF_minus_Control
    ),
    ,
    drop = FALSE
  ]

  LOO_direction_concordance <- if (
    nrow(LOO_valid) &&
      permutation_result$
        observed_difference != 0
  ) {
    mean(
      LOO_valid$
        direction_concordant_with_full,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  model_summary <- data.frame(
    model_id = model_id,
    matrix_role = matrix_role,
    broad_cell_type =
      broad_cell_type,
    endpoint = endpoint,
    endpoint_class =
      endpoint_class,
    n_Control = n_Control,
    n_RIF = n_RIF,
    total_donor_n =
      n_Control +
      n_RIF,
    Control_mean =
      safe_mean(
        Control_values
      ),
    RIF_mean =
      safe_mean(
        RIF_values
      ),
    Control_median =
      safe_median(
        Control_values
      ),
    RIF_median =
      safe_median(
        RIF_values
      ),
    Control_SD =
      safe_sd(
        Control_values
      ),
    RIF_SD =
      safe_sd(
        RIF_values
      ),
    mean_difference_RIF_minus_Control =
      permutation_result$
        observed_difference,
    median_difference_RIF_minus_Control =
      safe_median(
        RIF_values
      ) -
      safe_median(
        Control_values
      ),
    exact_permutation_P =
      permutation_result$
        exact_P,
    exact_permutation_n =
      permutation_result$
        permutation_n,
    minimum_resolvable_exact_P =
      permutation_result$
        minimum_resolvable_P,
    Welch_P = Welch_P,
    Welch_statistic =
      Welch_statistic,
    Welch_df = Welch_df,
    Welch_95CI_lower =
      CI_lower,
    Welch_95CI_upper =
      CI_upper,
    Hedges_g_RIF_minus_Control =
      hedges_g(
        Control_values,
        RIF_values
      ),
    rank_biserial_RIF_minus_Control =
      rank_biserial_RIF_minus_Control(
        Control_values,
        RIF_values
      ),
    constant_endpoint =
      constant_endpoint,
    LOO_valid_n =
      nrow(LOO_valid),
    LOO_direction_concordance =
      LOO_direction_concordance,
    LOO_minimum_difference =
      safe_min(
        LOO_valid$
          LOO_mean_difference_RIF_minus_Control
      ),
    LOO_maximum_difference =
      safe_max(
        LOO_valid$
          LOO_mean_difference_RIF_minus_Control
      ),
    effect_direction = if (
      permutation_result$
        observed_difference > 0
    ) {
      "RIF_HIGHER"
    } else if (
      permutation_result$
        observed_difference < 0
    ) {
      "RIF_LOWER"
    } else {
      "NO_DIFFERENCE"
    },
    biological_interpretation = if (
      endpoint_class ==
        "VDR_EXPRESSION"
    ) {
      "VDR_LOG2_CPM_PLUS_1"
    } else {
      "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION"
    },
    stringsAsFactors = FALSE
  )

  input_table <- model_data

  input_table$model_id <-
    model_id

  input_table$matrix_role <-
    matrix_role

  list(
    summary =
      model_summary,
    input =
      input_table,
    permutation =
      permutation_result$
        permutation_table,
    LOO =
      LOO_table
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

integration_root <- file.path(
  results_root,
  "S4E1_sc",
  integration_version
)

pb_root <- file.path(
  integration_root,
  "pb",
  "G183",
  pseudobulk_version
)

primary_VDR_path <- file.path(
  pb_root,
  "tab",
  "10_primary_VDR_expression.csv"
)

atlas_VDR_path <- file.path(
  pb_root,
  "tab",
  "11_atlas_VDR_expression.csv"
)

primary_program_path <- file.path(
  pb_root,
  "tab",
  "12_primary_program_scores.csv"
)

atlas_program_path <- file.path(
  pb_root,
  "tab",
  "13_atlas_program_scores.csv"
)

model_eligibility_path <- file.path(
  pb_root,
  "tab",
  "14_model_eligibility.csv"
)

cluster8_program_sensitivity_path <- file.path(
  pb_root,
  "tab",
  "15_cluster8_program_sensitivity.csv"
)

cluster8_VDR_sensitivity_path <- file.path(
  pb_root,
  "tab",
  "16_cluster8_VDR_sensitivity.csv"
)

output_root <- file.path(
  integration_root,
  "mdl",
  "G183",
  model_version
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
    "STEP04E1O_",
    format(run_start, "%Y%m%d_%H%M%S"),
    "_models"
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
  paste0(run_id, "_RUN_MANIFEST.csv")
)

parameter_path <- file.path(
  parameter_root,
  paste0(run_id, "_parameters.json")
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
    "primary_VDR",
    "atlas_VDR",
    "primary_program",
    "atlas_program",
    "model_eligibility",
    "output_root",
    "run_qc_root",
    "handoff_zip"
  ),
  path = normalize_slashes(c(
    primary_VDR_path,
    atlas_VDR_path,
    primary_program_path,
    atlas_program_path,
    model_eligibility_path,
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
    "One or more donor-model paths exceed the safe limit."
  )
}


# ----------------------------------------------------------------------
# 4. Validate upstream Gate
# ----------------------------------------------------------------------

pb_qc_dirs <- list.dirs(
  qc_base_root,
  recursive = FALSE,
  full.names = TRUE
)

pb_qc_dirs <- pb_qc_dirs[
  grepl(
    "^STEP04E1N_[0-9]{8}_[0-9]{6}_pbQC$",
    basename(pb_qc_dirs),
    perl = TRUE
  )
]

if (!length(pb_qc_dirs)) {
  stop(
    "No STEP04E1N pseudobulk-QC directory was found."
  )
}

gate_candidates <- list()

for (directory_path in pb_qc_dirs) {
  gate_path <- file.path(
    directory_path,
    "G183_PBVDR_GATE.csv"
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

    gate_candidates[[
      length(gate_candidates) + 1L
    ]] <- gate_value
  }
}

if (!length(gate_candidates)) {
  stop(
    "No readable STEP04E1N Gate was found."
  )
}

gate_candidates <- as.data.frame(
  data.table::rbindlist(
    gate_candidates,
    use.names = TRUE,
    fill = TRUE
  )
)

gate_candidates$source_mtime_POSIX <-
  as.POSIXct(
    gate_candidates$source_mtime
  )

gate_candidates <- gate_candidates[
  order(
    gate_candidates$
      source_mtime_POSIX,
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

authorized_index <- which(
  gate_candidates$gate_status ==
    "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED" &
    toupper(
      as.character(
        gate_candidates$
          donor_level_VDR_program_models_authorized
      )
    ) %in%
      c(
        "TRUE",
        "T",
        "1",
        "YES"
      )
)

if (!length(authorized_index)) {
  stop(
    "No pseudobulk Gate authorizes donor-level models."
  )
}

upstream_gate <- gate_candidates[
  authorized_index[1],
  ,
  drop = FALSE
]

write_csv_utf8(
  upstream_gate,
  file.path(
    run_qc_root,
    "00_upstream_PBVDR_gate.csv"
  )
)


# ----------------------------------------------------------------------
# 5. Load and validate model inputs
# ----------------------------------------------------------------------

required_inputs <- c(
  primary_VDR_path,
  atlas_VDR_path,
  primary_program_path,
  atlas_program_path,
  model_eligibility_path,
  cluster8_program_sensitivity_path,
  cluster8_VDR_sensitivity_path
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs)) {
  stop(
    "Missing required input file(s):\n",
    paste(
      normalize_slashes(
        missing_inputs
      ),
      collapse = "\n"
    )
  )
}

primary_VDR <- safe_read_csv(
  primary_VDR_path
)

atlas_VDR <- safe_read_csv(
  atlas_VDR_path
)

primary_program <- safe_read_csv(
  primary_program_path
)

atlas_program <- safe_read_csv(
  atlas_program_path
)

model_eligibility <- safe_read_csv(
  model_eligibility_path
)

cluster8_program_sensitivity <-
  safe_read_csv(
    cluster8_program_sensitivity_path
  )

cluster8_VDR_sensitivity <-
  safe_read_csv(
    cluster8_VDR_sensitivity_path
  )

for (
  data_name in c(
    "primary_VDR",
    "atlas_VDR",
    "primary_program",
    "atlas_program",
    "model_eligibility"
  )
) {
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

for (
  column_name in c(
    "minimum_cell_threshold_pass",
    "donor_level_model_eligible"
  )
) {
  if (
    column_name %in%
      colnames(primary_VDR)
  ) {
    primary_VDR[[
      column_name
    ]] <- as_logical_safe(
      primary_VDR[[
        column_name
      ]]
    )
  }

  if (
    column_name %in%
      colnames(atlas_VDR)
  ) {
    atlas_VDR[[
      column_name
    ]] <- as_logical_safe(
      atlas_VDR[[
        column_name
      ]]
    )
  }

  if (
    column_name %in%
      colnames(primary_program)
  ) {
    primary_program[[
      column_name
    ]] <- as_logical_safe(
      primary_program[[
        column_name
      ]]
    )
  }

  if (
    column_name %in%
      colnames(atlas_program)
  ) {
    atlas_program[[
      column_name
    ]] <- as_logical_safe(
      atlas_program[[
        column_name
      ]]
    )
  }

  if (
    column_name %in%
      colnames(model_eligibility)
  ) {
    model_eligibility[[
      column_name
    ]] <- as_logical_safe(
      model_eligibility[[
        column_name
      ]]
    )
  }
}

eligible_types <- model_eligibility$
  broad_cell_type[
    model_eligibility$
      donor_level_model_eligible
  ]

eligible_types <- sort(
  as.character(
    eligible_types
  )
)

expected_eligible_types <- sort(
  c(
    "Endothelial cells",
    "Myeloid cells",
    "NK/T cells",
    "Stromal cells",
    "Unciliated epithelial cells"
  )
)

if (
  !identical(
    eligible_types,
    expected_eligible_types
  )
) {
  stop(
    "Donor-model-eligible broad cell types do not match the locked five-type set."
  )
}

program_names <- c(
  "RNA_ALL_436",
  "OCCUPANCY_SUPPORTED_284",
  "SHARED_OCCUPANCY_46",
  "PROMOTER_OCCUPANCY_84"
)

if (
  !all(
    program_names %in%
      unique(
        primary_program$
          program_name
      )
  ) ||
    !all(
      program_names %in%
        unique(
          atlas_program$
            program_name
        )
    )
) {
  stop(
    "One or more locked VDR programs are absent."
  )
}

log_event(
  "Starting",
  run_id,
  "eligible_types=",
  collapse_unique(
    eligible_types
  )
)


# ----------------------------------------------------------------------
# 6. Construct primary donor-level model data
# ----------------------------------------------------------------------

primary_VDR_model_data <- primary_VDR[
  primary_VDR$
    broad_cell_type %in%
    eligible_types &
    primary_VDR$
      minimum_cell_threshold_pass,
  ,
  drop = FALSE
]

primary_VDR_model_data$endpoint <-
  "VDR_expression"

primary_VDR_model_data$
  endpoint_class <-
  "VDR_EXPRESSION"

primary_VDR_model_data$value <-
  as.numeric(
    primary_VDR_model_data$
      VDR_log2_CPM_plus_1
  )

primary_program_model_data <-
  primary_program[
    primary_program$
      broad_cell_type %in%
      eligible_types &
      primary_program$
        minimum_cell_threshold_pass &
      primary_program$
        program_name %in%
        program_names,
    ,
    drop = FALSE
  ]

primary_program_model_data$endpoint <-
  as.character(
    primary_program_model_data$
      program_name
  )

primary_program_model_data$
  endpoint_class <-
  "VDR_PROGRAM"

primary_program_model_data$value <-
  as.numeric(
    primary_program_model_data$
      score
  )

model_columns <- c(
  "donor_id",
  "sample_accession",
  "analysis_group",
  "broad_cell_type",
  "endpoint",
  "endpoint_class",
  "value",
  "cell_n"
)

primary_model_data <- rbind(
  primary_VDR_model_data[
    ,
    model_columns,
    drop = FALSE
  ],
  primary_program_model_data[
    ,
    model_columns,
    drop = FALSE
  ]
)

primary_model_data$donor_id <-
  as.character(
    primary_model_data$
      donor_id
  )

primary_model_data$sample_accession <-
  as.character(
    primary_model_data$
      sample_accession
  )

primary_model_data$analysis_group <-
  as.character(
    primary_model_data$
      analysis_group
  )

primary_model_data$broad_cell_type <-
  as.character(
    primary_model_data$
      broad_cell_type
  )

primary_model_data$endpoint <-
  as.character(
    primary_model_data$
      endpoint
  )

if (
  any(
    !primary_model_data$
      analysis_group %in%
      c(
        "CONTROL",
        "RIF"
      )
  ) ||
    any(
      !is.finite(
        primary_model_data$value
      )
    )
) {
  stop(
    "Primary model data contain invalid groups or values."
  )
}

write_csv_utf8(
  primary_model_data,
  file.path(
    source_data_root,
    "SD1_primary_model_data.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Run 25 primary exact-permutation models
# ----------------------------------------------------------------------

endpoint_order <- c(
  "VDR_expression",
  program_names
)

primary_summary_rows <- list()
primary_input_rows <- list()
primary_permutation_rows <- list()
primary_LOO_rows <- list()

model_index <- 0L

for (cell_type in eligible_types) {
  for (endpoint_value in endpoint_order) {
    model_index <- model_index + 1L

    model_data <- primary_model_data[
      primary_model_data$
        broad_cell_type ==
        cell_type &
        primary_model_data$
          endpoint ==
          endpoint_value,
      ,
      drop = FALSE
    ]

    endpoint_class_value <- unique(
      model_data$
        endpoint_class
    )

    if (
      length(endpoint_class_value) !=
        1L
    ) {
      stop(
        "Endpoint class is not unique for ",
        cell_type,
        " / ",
        endpoint_value,
        "."
      )
    }

    log_event(
      "MODEL",
      model_index,
      "/25",
      cell_type,
      endpoint_value,
      "donors=",
      nrow(model_data)
    )

    model_result <- fit_donor_model(
      model_data =
        model_data,
      broad_cell_type =
        cell_type,
      endpoint =
        endpoint_value,
      endpoint_class =
        endpoint_class_value,
      matrix_role =
        "PRIMARY"
    )

    primary_summary_rows[[
      model_index
    ]] <- model_result$summary

    primary_input_rows[[
      model_index
    ]] <- model_result$input

    primary_permutation_rows[[
      model_index
    ]] <- model_result$permutation

    primary_LOO_rows[[
      model_index
    ]] <- model_result$LOO
  }
}

primary_models <- do.call(
  rbind,
  primary_summary_rows
)

raw_primary_models_checkpoint <- file.path(
  source_data_root,
  "SD0_primary_models_pre_multiplicity.rds"
)

saveRDS(
  primary_models,
  raw_primary_models_checkpoint,
  compress = "gzip"
)

primary_model_inputs <- do.call(
  rbind,
  primary_input_rows
)

primary_permutations <- do.call(
  rbind,
  primary_permutation_rows
)

primary_LOO <- do.call(
  rbind,
  primary_LOO_rows
)

if (
  nrow(primary_models) != 25L ||
    length(
      unique(
        primary_models$
          broad_cell_type
      )
    ) != 5L ||
    length(
      unique(
        primary_models$
          endpoint
      )
    ) != 5L ||
    anyDuplicated(
      primary_models$model_id
    ) > 0L
) {
  stop(
    "Primary model grid is not the planned 5 x 5 design."
  )
}

primary_models$
  BH_global_25 <- stats::p.adjust(
    primary_models$
      exact_permutation_P,
    method = "BH"
  )

primary_models$
  BH_within_endpoint <- ave(
    primary_models$
      exact_permutation_P,
    primary_models$endpoint,
    FUN = function(x) {
      stats::p.adjust(
        x,
        method = "BH"
      )
    }
  )

primary_models$
  BH_within_cell_type <- ave(
    primary_models$
      exact_permutation_P,
    primary_models$
      broad_cell_type,
    FUN = function(x) {
      stats::p.adjust(
        x,
        method = "BH"
      )
    }
  )

primary_models$
  statistical_status <- ifelse(
    primary_models$
      constant_endpoint,
    "CONSTANT_UNTESTABLE_ENDPOINT",
    ifelse(
      primary_models$
        BH_global_25 < 0.05,
      "GLOBAL_BH_SIGNIFICANT",
      ifelse(
        primary_models$
          BH_within_endpoint < 0.05,
        "ENDPOINT_FAMILY_BH_SIGNIFICANT",
        ifelse(
          primary_models$
            exact_permutation_P < 0.05,
          "NOMINAL_EXACT_P_LT_0_05_EXPLORATORY",
          "NO_NOMINAL_EVIDENCE"
        )
      )
    )
  )

primary_models$
  LOO_review_flag <- ifelse(
    is.finite(
      primary_models$
        LOO_direction_concordance
    ) &
      primary_models$
        LOO_direction_concordance <
        0.75,
    "LOO_DIRECTION_UNSTABLE",
    ""
  )

primary_models$
  resolution_review_flag <- ifelse(
    primary_models$
      exact_permutation_n < 84,
    paste0(
      "LIMITED_EXACT_RESOLUTION_MIN_P_",
      format(
        primary_models$
          minimum_resolvable_exact_P,
        digits = 4
      )
    ),
    ""
  )

primary_models$
  interpretation_restriction <- ifelse(
    primary_models$
      endpoint_class ==
      "VDR_PROGRAM",
    "RELATIVE_ALIGNMENT_SCORE_NOT_DIRECT_VDR_ACTIVITY",
    "VDR_EXPRESSION_ONLY_NOT_CAUSAL_ACTIVITY"
  )

primary_models <- primary_models[
  order(
    match(
      primary_models$
        broad_cell_type,
      eligible_types
    ),
    match(
      primary_models$endpoint,
      endpoint_order
    )
  ),
  ,
  drop = FALSE
]

write_csv_utf8(
  primary_models,
  file.path(
    table_root,
    "01_primary_models.csv"
  )
)

write_csv_utf8(
  primary_model_inputs,
  file.path(
    source_data_root,
    "SD2_primary_model_inputs.csv"
  )
)

write_csv_utf8(
  primary_permutations,
  file.path(
    source_data_root,
    "SD3_exact_permutation_distributions.csv"
  )
)

write_csv_utf8(
  primary_LOO,
  file.path(
    source_data_root,
    "SD4_leave_one_donor_out.csv"
  )
)


# ----------------------------------------------------------------------
# 8. Primary result summaries
# ----------------------------------------------------------------------

group_summary <- data.table::as.data.table(
  primary_model_data
)[
  ,
  .(
    donor_n =
      as.integer(
        data.table::uniqueN(
          donor_id
        )
      ),
    mean_value =
      as.numeric(
        mean(value)
      ),
    median_value =
      as.numeric(
        median(value)
      ),
    SD_value =
      as.numeric(
        stats::sd(value)
      ),
    minimum_value =
      as.numeric(
        min(value)
      ),
    maximum_value =
      as.numeric(
        max(value)
      ),
    median_cell_n =
      as.numeric(
        median(cell_n)
      ),
    minimum_cell_n =
      as.numeric(
        min(cell_n)
      ),
    maximum_cell_n =
      as.numeric(
        max(cell_n)
      )
  ),
  by = .(
    broad_cell_type,
    endpoint,
    endpoint_class,
    analysis_group
  )
]

group_summary <- as.data.frame(
  group_summary
)

group_summary_type_audit <- data.frame(
  column_name =
    colnames(group_summary),
  column_class =
    vapply(
      group_summary,
      function(x) {
        paste(
          class(x),
          collapse = ";"
        )
      },
      character(1)
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  group_summary_type_audit,
  file.path(
    table_root,
    "00_group_summary_type_audit.csv"
  )
)

write_csv_utf8(
  group_summary,
  file.path(
    table_root,
    "02_group_summary.csv"
  )
)

multiplicity_summary <- data.frame(
  metric = c(
    "planned_primary_model_n",
    "constant_endpoint_n",
    "nominal_exact_P_lt_0_05_n",
    "global_BH_lt_0_05_n",
    "endpoint_family_BH_lt_0_05_n",
    "cell_type_family_BH_lt_0_05_n",
    "LOO_direction_unstable_n",
    "unciliated_models_with_35_permutations_n"
  ),
  value = c(
    nrow(primary_models),
    sum(
      primary_models$
        constant_endpoint
    ),
    sum(
      primary_models$
        exact_permutation_P <
        0.05
    ),
    sum(
      primary_models$
        BH_global_25 < 0.05
    ),
    sum(
      primary_models$
        BH_within_endpoint <
        0.05
    ),
    sum(
      primary_models$
        BH_within_cell_type <
        0.05
    ),
    sum(
      primary_models$
        LOO_review_flag ==
        "LOO_DIRECTION_UNSTABLE"
    ),
    sum(
      primary_models$
        exact_permutation_n ==
        35L
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  multiplicity_summary,
  file.path(
    table_root,
    "03_multiplicity_summary.csv"
  )
)

review_models <- primary_models[
  primary_models$
    statistical_status !=
    "NO_NOMINAL_EVIDENCE" |
    nzchar(
      primary_models$
        LOO_review_flag
    ) |
    nzchar(
      primary_models$
        resolution_review_flag
    ),
  ,
  drop = FALSE
]

write_csv_utf8(
  review_models,
  file.path(
    table_root,
    "04_models_requiring_review.csv"
  )
)


# ----------------------------------------------------------------------
# 9. Cluster-8 atlas sensitivity models
# ----------------------------------------------------------------------

atlas_VDR_model_data <- atlas_VDR[
  atlas_VDR$
    broad_cell_type ==
    "Unciliated epithelial cells" &
    atlas_VDR$
      minimum_cell_threshold_pass,
  ,
  drop = FALSE
]

atlas_VDR_model_data$endpoint <-
  "VDR_expression"

atlas_VDR_model_data$
  endpoint_class <-
  "VDR_EXPRESSION"

atlas_VDR_model_data$value <-
  as.numeric(
    atlas_VDR_model_data$
      VDR_log2_CPM_plus_1
  )

atlas_program_model_data <-
  atlas_program[
    atlas_program$
      broad_cell_type ==
      "Unciliated epithelial cells" &
      atlas_program$
        minimum_cell_threshold_pass &
      atlas_program$
        program_name %in%
        program_names,
    ,
    drop = FALSE
  ]

atlas_program_model_data$endpoint <-
  as.character(
    atlas_program_model_data$
      program_name
  )

atlas_program_model_data$
  endpoint_class <-
  "VDR_PROGRAM"

atlas_program_model_data$value <-
  as.numeric(
    atlas_program_model_data$
      score
  )

atlas_unciliated_data <- rbind(
  atlas_VDR_model_data[
    ,
    model_columns,
    drop = FALSE
  ],
  atlas_program_model_data[
    ,
    model_columns,
    drop = FALSE
  ]
)

sensitivity_summary_rows <- list()
sensitivity_input_rows <- list()
sensitivity_permutation_rows <- list()
sensitivity_LOO_rows <- list()

for (
  endpoint_index in seq_along(
    endpoint_order
  )
) {
  endpoint_value <- endpoint_order[
    endpoint_index
  ]

  model_data <- atlas_unciliated_data[
    atlas_unciliated_data$endpoint ==
      endpoint_value,
    ,
    drop = FALSE
  ]

  endpoint_class_value <- unique(
    model_data$
      endpoint_class
  )

  sensitivity_result <- fit_donor_model(
    model_data =
      model_data,
    broad_cell_type =
      "Unciliated epithelial cells",
    endpoint =
      endpoint_value,
    endpoint_class =
      endpoint_class_value,
    matrix_role =
      "ATLAS_SENSITIVITY"
  )

  sensitivity_summary_rows[[
    endpoint_index
  ]] <- sensitivity_result$summary

  sensitivity_input_rows[[
    endpoint_index
  ]] <- sensitivity_result$input

  sensitivity_permutation_rows[[
    endpoint_index
  ]] <- sensitivity_result$permutation

  sensitivity_LOO_rows[[
    endpoint_index
  ]] <- sensitivity_result$LOO
}

atlas_sensitivity_models <- do.call(
  rbind,
  sensitivity_summary_rows
)

atlas_sensitivity_inputs <- do.call(
  rbind,
  sensitivity_input_rows
)

atlas_sensitivity_permutations <- do.call(
  rbind,
  sensitivity_permutation_rows
)

atlas_sensitivity_LOO <- do.call(
  rbind,
  sensitivity_LOO_rows
)

atlas_sensitivity_models$
  BH_sensitivity_5 <- stats::p.adjust(
    atlas_sensitivity_models$
      exact_permutation_P,
    method = "BH"
  )

primary_unciliated_models <- primary_models[
  primary_models$broad_cell_type ==
    "Unciliated epithelial cells",
  ,
  drop = FALSE
]

cluster8_model_sensitivity <- merge(
  primary_unciliated_models[
    ,
    c(
      "endpoint",
      "n_Control",
      "n_RIF",
      "mean_difference_RIF_minus_Control",
      "exact_permutation_P",
      "minimum_resolvable_exact_P",
      "Hedges_g_RIF_minus_Control",
      "LOO_direction_concordance"
    ),
    drop = FALSE
  ],
  atlas_sensitivity_models[
    ,
    c(
      "endpoint",
      "n_Control",
      "n_RIF",
      "mean_difference_RIF_minus_Control",
      "exact_permutation_P",
      "minimum_resolvable_exact_P",
      "Hedges_g_RIF_minus_Control",
      "LOO_direction_concordance",
      "BH_sensitivity_5"
    ),
    drop = FALSE
  ],
  by = "endpoint",
  all = TRUE,
  suffixes = c(
    "_primary",
    "_atlas"
  ),
  sort = FALSE
)

cluster8_model_sensitivity$
  effect_delta_atlas_minus_primary <- (
    cluster8_model_sensitivity$
      mean_difference_RIF_minus_Control_atlas -
      cluster8_model_sensitivity$
        mean_difference_RIF_minus_Control_primary
  )

cluster8_model_sensitivity$
  direction_concordant <- ifelse(
    cluster8_model_sensitivity$
      mean_difference_RIF_minus_Control_primary ==
      0 |
      cluster8_model_sensitivity$
        mean_difference_RIF_minus_Control_atlas ==
        0,
    NA,
    sign(
      cluster8_model_sensitivity$
        mean_difference_RIF_minus_Control_primary
    ) ==
      sign(
        cluster8_model_sensitivity$
          mean_difference_RIF_minus_Control_atlas
      )
  )

cluster8_model_sensitivity$
  sensitivity_interpretation <- ifelse(
    cluster8_model_sensitivity$
      direction_concordant %in%
      TRUE,
    "DIRECTION_CONCORDANT",
    ifelse(
      cluster8_model_sensitivity$
        direction_concordant %in%
        FALSE,
      "DIRECTION_CHANGED_WITH_CLUSTER8",
      "ZERO_EFFECT_DIRECTION_NOT_CLASSIFIED"
    )
  )

write_csv_utf8(
  atlas_sensitivity_models,
  file.path(
    table_root,
    "05_atlas_unciliated_models.csv"
  )
)

write_csv_utf8(
  cluster8_model_sensitivity,
  file.path(
    table_root,
    "06_cluster8_model_sensitivity.csv"
  )
)

write_csv_utf8(
  atlas_sensitivity_inputs,
  file.path(
    source_data_root,
    "SD5_atlas_sensitivity_inputs.csv"
  )
)

write_csv_utf8(
  atlas_sensitivity_permutations,
  file.path(
    source_data_root,
    "SD6_atlas_sensitivity_permutations.csv"
  )
)

write_csv_utf8(
  atlas_sensitivity_LOO,
  file.path(
    source_data_root,
    "SD7_atlas_sensitivity_LOO.csv"
  )
)


# ----------------------------------------------------------------------
# 10. Ciliated epithelial descriptive-only summaries
# ----------------------------------------------------------------------

ciliated_VDR <- primary_VDR[
  primary_VDR$
    broad_cell_type ==
    "Ciliated epithelial cells",
  ,
  drop = FALSE
]

ciliated_VDR$endpoint <-
  "VDR_expression"

ciliated_VDR$value <-
  ciliated_VDR$
    VDR_log2_CPM_plus_1

ciliated_program <- primary_program[
  primary_program$
    broad_cell_type ==
    "Ciliated epithelial cells" &
    primary_program$
      program_name %in%
      program_names,
  ,
  drop = FALSE
]

ciliated_program$endpoint <-
  ciliated_program$
    program_name

ciliated_program$value <-
  ciliated_program$score

ciliated_columns <- c(
  "donor_id",
  "sample_accession",
  "analysis_group",
  "broad_cell_type",
  "endpoint",
  "value",
  "cell_n",
  "minimum_cell_threshold_pass"
)

ciliated_data <- rbind(
  ciliated_VDR[
    ,
    ciliated_columns,
    drop = FALSE
  ],
  ciliated_program[
    ,
    ciliated_columns,
    drop = FALSE
  ]
)

ciliated_descriptive <- data.table::as.data.table(
  ciliated_data
)[
  ,
  .(
    donor_profile_n =
      as.integer(.N),
    threshold_pass_profile_n =
      as.integer(
        sum(
          minimum_cell_threshold_pass
        )
      ),
    total_cell_n =
      as.numeric(
        sum(cell_n)
      ),
    median_cell_n =
      as.numeric(
        median(cell_n)
      ),
    mean_value =
      as.numeric(
        mean(value)
      ),
    median_value =
      as.numeric(
        median(value)
      ),
    minimum_value =
      as.numeric(
        min(value)
      ),
    maximum_value =
      as.numeric(
        max(value)
      )
  ),
  by = .(
    endpoint,
    analysis_group
  )
]

ciliated_descriptive <- as.data.frame(
  ciliated_descriptive
)

write_csv_utf8(
  ciliated_data,
  file.path(
    source_data_root,
    "SD8_ciliated_descriptive_data.csv"
  )
)

write_csv_utf8(
  ciliated_descriptive,
  file.path(
    table_root,
    "07_ciliated_descriptive_summary.csv"
  )
)


# ----------------------------------------------------------------------
# 11. Figures
# ----------------------------------------------------------------------

theme_QC <- ggplot2::theme_bw(
  base_size = 10
) +
  ggplot2::theme(
    panel.grid.minor =
      ggplot2::element_blank()
  )

forest_data <- primary_models

forest_data$broad_cell_type <- factor(
  forest_data$broad_cell_type,
  levels = rev(
    eligible_types
  )
)

forest_data$endpoint <- factor(
  forest_data$endpoint,
  levels = endpoint_order
)

plot_forest <- ggplot2::ggplot(
  forest_data,
  ggplot2::aes(
    x =
      mean_difference_RIF_minus_Control,
    y = broad_cell_type,
    shape =
      statistical_status
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = Welch_95CI_lower,
      xmax = Welch_95CI_upper
    ),
    height = 0.15,
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    size = 2.5
  ) +
  ggplot2::facet_wrap(
    ~ endpoint,
    scales = "free_x",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "GSE183837 donor-level Control-versus-RIF effects",
    subtitle = "Effect is RIF minus Control; intervals are secondary Welch 95% confidence intervals",
    x = "Mean difference",
    y = "Broad cell type",
    shape = "Statistical status"
  ) +
  theme_QC

save_plot_pair(
  plot_forest,
  "QC1_primary_effect_forest",
  width = 11,
  height = 9
)

pvalue_data <- primary_models

pvalue_data$
  minus_log10_exact_P <- -log10(
    pmax(
      pvalue_data$
        exact_permutation_P,
      1e-12
    )
  )

pvalue_data$
  minus_log10_global_BH <- -log10(
    pmax(
      pvalue_data$
        BH_global_25,
      1e-12
    )
  )

plot_pvalues <- ggplot2::ggplot(
  pvalue_data,
  ggplot2::aes(
    x =
      minus_log10_exact_P,
    y =
      minus_log10_global_BH,
    shape = endpoint
  )
) +
  ggplot2::geom_vline(
    xintercept =
      -log10(0.05),
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept =
      -log10(0.05),
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  ggplot2::labs(
    title = "Exact permutation evidence before and after global BH correction",
    subtitle = "Global multiplicity includes all 25 planned donor-level models",
    x = "-log10 exact permutation P",
    y = "-log10 global BH q",
    shape = "Endpoint"
  ) +
  theme_QC

save_plot_pair(
  plot_pvalues,
  "QC2_exact_P_and_global_BH",
  width = 8,
  height = 6
)

VDR_plot_data <- primary_model_data[
  primary_model_data$endpoint ==
    "VDR_expression",
  ,
  drop = FALSE
]

VDR_plot_data$broad_cell_type <- factor(
  VDR_plot_data$broad_cell_type,
  levels = eligible_types
)

plot_VDR_groups <- ggplot2::ggplot(
  VDR_plot_data,
  ggplot2::aes(
    x = analysis_group,
    y = value,
    shape = analysis_group
  )
) +
  ggplot2::geom_point(
    size = 2.5,
    alpha = 0.85,
    position =
      ggplot2::position_jitter(
        width = 0.08,
        height = 0
      )
  ) +
  ggplot2::facet_wrap(
    ~ broad_cell_type,
    scales = "free_y",
    ncol = 3
  ) +
  ggplot2::labs(
    title = "Donor-level VDR expression by group",
    subtitle = "Each point is one threshold-passing donor pseudobulk",
    x = "Group",
    y = "VDR log2(CPM + 1)",
    shape = "Group"
  ) +
  theme_QC

save_plot_pair(
  plot_VDR_groups,
  "QC3_VDR_group_profiles",
  width = 10,
  height = 7
)

program_plot_data <- primary_model_data[
  primary_model_data$
    endpoint_class ==
    "VDR_PROGRAM",
  ,
  drop = FALSE
]

program_plot_data$broad_cell_type <- factor(
  program_plot_data$
    broad_cell_type,
  levels = eligible_types
)

program_plot_data$endpoint <- factor(
  program_plot_data$endpoint,
  levels = program_names
)

plot_program_groups <- ggplot2::ggplot(
  program_plot_data,
  ggplot2::aes(
    x = analysis_group,
    y = value,
    shape = analysis_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    size = 1.9,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.07,
        height = 0
      )
  ) +
  ggplot2::facet_grid(
    broad_cell_type ~ endpoint,
    scales = "free_y"
  ) +
  ggplot2::labs(
    title = "Locked VDR-program scores by donor group",
    subtitle = "Positive scores indicate greater relative alignment with the locked interaction direction",
    x = "Group",
    y = "Signed program score",
    shape = "Group"
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
  plot_program_groups,
  "QC4_program_group_profiles",
  width = 13,
  height = 11
)

LOO_plot_data <- primary_models

LOO_plot_data$broad_cell_type <- factor(
  LOO_plot_data$broad_cell_type,
  levels = rev(
    eligible_types
  )
)

LOO_plot_data$endpoint <- factor(
  LOO_plot_data$endpoint,
  levels = endpoint_order
)

plot_LOO <- ggplot2::ggplot(
  LOO_plot_data,
  ggplot2::aes(
    x =
      mean_difference_RIF_minus_Control,
    y = broad_cell_type,
    shape =
      LOO_review_flag
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin =
        LOO_minimum_difference,
      xmax =
        LOO_maximum_difference
    ),
    height = 0.15
  ) +
  ggplot2::geom_point(
    size = 2.4
  ) +
  ggplot2::facet_wrap(
    ~ endpoint,
    scales = "free_x",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "Leave-one-donor-out effect stability",
    subtitle = "Intervals show the minimum and maximum RIF-minus-Control effect after excluding each donor",
    x = "Mean difference",
    y = "Broad cell type",
    shape = "LOO review"
  ) +
  theme_QC

save_plot_pair(
  plot_LOO,
  "QC5_leave_one_donor_out",
  width = 11,
  height = 9
)

sensitivity_plot_data <-
  cluster8_model_sensitivity

sensitivity_plot_long <- rbind(
  data.frame(
    endpoint =
      sensitivity_plot_data$
        endpoint,
    matrix_role = "PRIMARY",
    mean_difference =
      sensitivity_plot_data$
        mean_difference_RIF_minus_Control_primary,
    exact_P =
      sensitivity_plot_data$
        exact_permutation_P_primary,
    stringsAsFactors = FALSE
  ),
  data.frame(
    endpoint =
      sensitivity_plot_data$
        endpoint,
    matrix_role =
      "ATLAS_SENSITIVITY",
    mean_difference =
      sensitivity_plot_data$
        mean_difference_RIF_minus_Control_atlas,
    exact_P =
      sensitivity_plot_data$
        exact_permutation_P_atlas,
    stringsAsFactors = FALSE
  )
)

sensitivity_plot_long$endpoint <- factor(
  sensitivity_plot_long$endpoint,
  levels = endpoint_order
)

plot_sensitivity <- ggplot2::ggplot(
  sensitivity_plot_long,
  ggplot2::aes(
    x = matrix_role,
    y = mean_difference,
    group = endpoint,
    shape = endpoint
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(
    alpha = 0.6
  ) +
  ggplot2::geom_point(
    size = 2.7
  ) +
  ggplot2::labs(
    title = "Cluster-8 sensitivity of unciliated epithelial effects",
    subtitle = "Primary excludes cluster 8; atlas sensitivity retains it",
    x = "Analysis matrix",
    y = "RIF-minus-Control mean difference",
    shape = "Endpoint"
  ) +
  theme_QC

save_plot_pair(
  plot_sensitivity,
  "QC6_cluster8_model_sensitivity",
  width = 9,
  height = 6
)

ciliated_plot_data <- ciliated_data

ciliated_plot_data$endpoint <- factor(
  ciliated_plot_data$endpoint,
  levels = endpoint_order
)

plot_ciliated <- ggplot2::ggplot(
  ciliated_plot_data,
  ggplot2::aes(
    x = analysis_group,
    y = value,
    shape =
      minimum_cell_threshold_pass
  )
) +
  ggplot2::geom_point(
    size = 2.4,
    alpha = 0.8,
    position =
      ggplot2::position_jitter(
        width = 0.08,
        height = 0
      )
  ) +
  ggplot2::facet_wrap(
    ~ endpoint,
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::labs(
    title = "Ciliated epithelial descriptive profiles",
    subtitle = "Only two donor profiles pass the 20-cell threshold; no group test is performed",
    x = "Group",
    y = "Endpoint value",
    shape = "20-cell threshold pass"
  ) +
  theme_QC

save_plot_pair(
  plot_ciliated,
  "QC7_ciliated_descriptive_only",
  width = 10,
  height = 8
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
# 12. Gate
# ----------------------------------------------------------------------

nominal_models <- primary_models[
  primary_models$
    exact_permutation_P < 0.05,
  ,
  drop = FALSE
]

global_BH_models <- primary_models[
  primary_models$
    BH_global_25 < 0.05,
  ,
  drop = FALSE
]

constant_models <- primary_models[
  primary_models$
    constant_endpoint,
  ,
  drop = FALSE
]

direction_changed_sensitivity <- cluster8_model_sensitivity[
  cluster8_model_sensitivity$
    direction_concordant %in%
    FALSE,
  ,
  drop = FALSE
]

gate_checks <- data.frame(
  check = c(
    "upstream_gate_authorizes_donor_models",
    "five_locked_broad_cell_types_present",
    "five_locked_endpoints_present",
    "planned_primary_model_n_is_25",
    "all_primary_models_have_at_least_3_donors_per_group",
    "four_full_design_types_have_84_exact_allocations",
    "unciliated_models_have_35_exact_allocations",
    "endothelial_VDR_is_retained_as_constant_endpoint",
    "global_BH_values_complete",
    "within_endpoint_BH_values_complete",
    "within_cell_type_BH_values_complete",
    "LOO_results_exist_for_all_primary_models",
    "atlas_unciliated_sensitivity_has_5_models",
    "atlas_unciliated_uses_3_Control_6_RIF",
    "cluster8_primary_atlas_comparison_has_5_endpoints",
    "ciliated_epithelial_has_no_inferential_model",
    "cell_level_testing_not_performed",
    "gene_level_DE_not_performed",
    "independent_validation_review_can_be_authorized"
  ),
  pass = c(
    upstream_gate$gate_status[[1]] ==
      "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED",
    identical(
      sort(
        unique(
          primary_models$
            broad_cell_type
        )
      ),
      expected_eligible_types
    ),
    identical(
      unique(
        primary_models$endpoint[
          match(
            endpoint_order,
            primary_models$endpoint
          )
        ]
      ),
      endpoint_order
    ),
    nrow(primary_models) == 25L,
    all(
      primary_models$n_Control >= 3L &
        primary_models$n_RIF >= 3L
    ),
    all(
      primary_models$
        exact_permutation_n[
          primary_models$
            broad_cell_type !=
            "Unciliated epithelial cells"
        ] == 84L
    ),
    all(
      primary_models$
        exact_permutation_n[
          primary_models$
            broad_cell_type ==
            "Unciliated epithelial cells"
        ] == 35L
    ),
    nrow(constant_models) == 1L &&
      constant_models$
        broad_cell_type ==
        "Endothelial cells" &&
      constant_models$endpoint ==
        "VDR_expression",
    all(
      is.finite(
        primary_models$
          BH_global_25
      )
    ),
    all(
      is.finite(
        primary_models$
          BH_within_endpoint
      )
    ),
    all(
      is.finite(
        primary_models$
          BH_within_cell_type
      )
    ),
    length(
      unique(
        primary_LOO$model_id
      )
    ) == 25L,
    nrow(
      atlas_sensitivity_models
    ) == 5L,
    all(
      atlas_sensitivity_models$
        n_Control == 3L &
        atlas_sensitivity_models$
          n_RIF == 6L
    ),
    nrow(
      cluster8_model_sensitivity
    ) == 5L,
    !(
      "Ciliated epithelial cells" %in%
        primary_models$
          broad_cell_type
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
    "01_gate_checks.csv"
  )
)

all_gate_checks_pass <- all(
  gate_checks$pass
)

gate_status <- if (
  all_gate_checks_pass
) {
  "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED"
} else {
  "HOLD_GSE183837_DONOR_MODELS_INCOMPLETE"
}

review_flags <- c(
  "CONTROL_N3_RIF_N6_LOW_POWER",
  "UNCILIATED_PRIMARY_RIF_N4_EXACT_MIN_P_1_OVER_35",
  "ENDOTHELIAL_VDR_CONSTANT_ZERO",
  "CILIATED_EPITHELIAL_DESCRIPTIVE_ONLY",
  "CLUSTER8_ATLAS_SENSITIVITY_ONLY",
  if (nrow(nominal_models)) {
    paste0(
      "NOMINAL_EXACT_P_LT_0_05_MODELS_",
      nrow(nominal_models)
    )
  } else {
    "NO_NOMINAL_EXACT_P_LT_0_05"
  },
  if (nrow(global_BH_models)) {
    paste0(
      "GLOBAL_BH_SIGNIFICANT_MODELS_",
      nrow(global_BH_models)
    )
  } else {
    "NO_GLOBAL_BH_SIGNIFICANT_MODELS"
  },
  if (
    nrow(
      direction_changed_sensitivity
    )
  ) {
    paste0(
      "CLUSTER8_DIRECTION_CHANGED_ENDPOINTS_",
      collapse_unique(
        direction_changed_sensitivity$
          endpoint
      )
    )
  } else {
    "CLUSTER8_DIRECTIONS_CONCORDANT_OR_ZERO"
  }
)

model_gate <- data.frame(
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  pseudobulk_version =
    pseudobulk_version,
  model_version =
    model_version,
  gate_status = gate_status,
  primary_model_n =
    nrow(primary_models),
  testable_nonconstant_model_n =
    sum(
      !primary_models$
        constant_endpoint
    ),
  nominal_exact_P_lt_0_05_n =
    nrow(nominal_models),
  global_BH_lt_0_05_n =
    nrow(global_BH_models),
  endpoint_family_BH_lt_0_05_n =
    sum(
      primary_models$
        BH_within_endpoint <
        0.05
    ),
  LOO_direction_unstable_n =
    sum(
      primary_models$
        LOO_review_flag ==
        "LOO_DIRECTION_UNSTABLE"
    ),
  cluster8_direction_changed_n =
    nrow(
      direction_changed_sensitivity
    ),
  nominal_models =
    collapse_unique(
      paste(
        nominal_models$
          broad_cell_type,
        nominal_models$endpoint,
        sep = " | "
      )
    ),
  global_BH_models =
    collapse_unique(
      paste(
        global_BH_models$
          broad_cell_type,
        global_BH_models$endpoint,
        sep = " | "
      )
    ),
  review_flags =
    collapse_unique(
      review_flags
    ),
  independent_validation_review_authorized =
    all_gate_checks_pass,
  cross_dataset_concordance_review_authorized =
    all_gate_checks_pass,
  gene_level_pseudobulk_DE_authorized =
    FALSE,
  result_lock_authorized =
    FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  model_gate,
  file.path(
    run_qc_root,
    "G183_DONOR_MODEL_GATE.csv"
  )
)

card_lines <- c(
  "# STEP04E1O GSE183837 Donor Model Card",
  "",
  paste0("- Gate: ", gate_status),
  "- Primary models: 25",
  paste0(
    "- Nonconstant testable models: ",
    sum(
      !primary_models$
        constant_endpoint
    )
  ),
  paste0(
    "- Nominal exact P < 0.05: ",
    nrow(nominal_models)
  ),
  paste0(
    "- Global BH q < 0.05: ",
    nrow(global_BH_models)
  ),
  paste0(
    "- Endpoint-family BH q < 0.05: ",
    sum(
      primary_models$
        BH_within_endpoint <
        0.05
    )
  ),
  paste0(
    "- LOO direction-unstable models: ",
    sum(
      primary_models$
        LOO_review_flag ==
        "LOO_DIRECTION_UNSTABLE"
    )
  ),
  paste0(
    "- Cluster-8 sensitivity direction changes: ",
    nrow(
      direction_changed_sensitivity
    )
  ),
  paste0(
    "- Independent validation review authorized: ",
    all_gate_checks_pass
  ),
  "- Gene-level pseudobulk DE authorized: FALSE",
  "- Result lock authorized: FALSE"
)

write_text(
  card_lines,
  file.path(
    run_qc_root,
    "G183_DONOR_MODEL_CARD.md"
  )
)

interpretation_lines <- c(
  "# GSE183837 donor-model interpretation",
  "",
  paste0("- Run ID: `", run_id, "`"),
  paste0("- Gate: **", gate_status, "**"),
  "",
  "## Prespecified design",
  "",
  "- Twenty-five donor-level models were tested across five broad cell types and five endpoints.",
  "- Exact permutation tests were primary; Welch tests and confidence intervals were secondary.",
  "- Global BH correction covered all 25 planned models.",
  "",
  "## Interpretation rules",
  "",
  "- Nominal exact P < 0.05 without global BH significance is exploratory only.",
  "- A negative VDR-program effect means lower relative alignment in RIF with the locked calcitriol-dependent interaction direction.",
  "- Program effects do not establish direct binding, VDR activity or causality.",
  "- Ciliated epithelial cells remain descriptive only.",
  "- Cluster 8 is evaluated only through atlas sensitivity analysis.",
  "",
  "## Authorized next operation",
  "",
  "- Compare the independent GSE183837 effect directions with the locked GSE250130 donor-level results.",
  "- Determine convergent, discordant and non-informative evidence.",
  "- Lock the cross-dataset interpretation only after consistency QC.",
  "",
  "## Not authorized",
  "",
  "- No gene-level pseudobulk DE.",
  "- No cell-level disease testing.",
  "- No disease mechanism claim.",
  "- No final result lock."
)

write_text(
  interpretation_lines,
  file.path(
    run_qc_root,
    "G183_DONOR_MODEL_INTERPRETATION.md"
  )
)


# ----------------------------------------------------------------------
# 13. Manuscript shells
# ----------------------------------------------------------------------

nominal_description <- if (
  nrow(nominal_models)
) {
  paste(
    paste0(
      nominal_models$
        broad_cell_type,
      " / ",
      nominal_models$endpoint,
      " (exact P=",
      format(
        nominal_models$
          exact_permutation_P,
        digits = 4
      ),
      ", global BH q=",
      format(
        nominal_models$
          BH_global_25,
        digits = 4
      ),
      ")"
    ),
    collapse = "; "
  )
} else {
  "none"
}

methods_lines <- c(
  "# STEP04E1O Methods shell",
  "",
  "For each model-eligible broad cell type, donor-level VDR expression and four prospectively locked signed VDR-program scores were compared between Control and RIF donors. Two-sided exact label-permutation tests based on the RIF-minus-Control mean difference were used as the primary inferential procedure. All possible allocations were enumerated, yielding 84 allocations for the 3-versus-6 design and 35 allocations for the 3-versus-4 unciliated epithelial design. Multiplicity was controlled using the Benjamini-Hochberg procedure across all 25 planned primary models. Welch tests, 95% confidence intervals, Hedges g, rank-biserial effects and leave-one-donor-out analyses were reported as secondary estimates."
)

results_lines <- c(
  "# STEP04E1O Results shell",
  "",
  paste0(
    "Across 25 planned donor-level models, ",
    nrow(nominal_models),
    " met nominal exact P<0.05, whereas ",
    nrow(global_BH_models),
    " remained significant after global Benjamini-Hochberg correction."
  ),
  paste0(
    "The nominal models were: ",
    nominal_description,
    "."
  ),
  paste0(
    "Leave-one-donor-out direction instability was observed in ",
    sum(
      primary_models$
        LOO_review_flag ==
        "LOO_DIRECTION_UNSTABLE"
    ),
    " models."
  ),
  paste0(
    "Retaining cluster 8 changed the unciliated epithelial effect direction for ",
    nrow(
      direction_changed_sensitivity
    ),
    " of five endpoints."
  ),
  "Ciliated epithelial profiles were retained for descriptive localization only because only two donor profiles met the 20-cell threshold."
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


# ----------------------------------------------------------------------
# 14. Parameters, manifests and ZIP
# ----------------------------------------------------------------------

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  supersedes_partial_execution =
    "V1_0_DATA_TABLE_GROUPED_SUMMARY_INTEGER_DOUBLE_TYPE_MISMATCH",
  grouped_summary_policy =
    "EXPLICIT_INTEGER_COUNTS_AND_NUMERIC_CONTINUOUS_SUMMARIES",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  pseudobulk_version =
    pseudobulk_version,
  model_version =
    model_version,
  upstream_gate =
    upstream_gate$gate_status[[1]],
  eligible_broad_cell_types =
    eligible_types,
  endpoints =
    endpoint_order,
  primary_test =
    "TWO_SIDED_EXACT_LABEL_PERMUTATION_OF_RIF_MINUS_CONTROL_MEAN_DIFFERENCE",
  full_design_permutation_n =
    84L,
  unciliated_primary_permutation_n =
    35L,
  primary_multiplicity =
    "BH_ACROSS_ALL_25_PLANNED_MODELS",
  secondary_statistics = c(
    "WELCH_T_TEST",
    "WELCH_95_PERCENT_CI",
    "HEDGES_G",
    "RANK_BISERIAL",
    "LEAVE_ONE_DONOR_OUT"
  ),
  cluster8_policy =
    "ATLAS_SENSITIVITY_ONLY",
  ciliated_policy =
    "DESCRIPTIVE_ONLY_NO_GROUP_TEST",
  independent_validation_review_authorized =
    all_gate_checks_pass,
  gene_level_pseudobulk_DE_authorized =
    FALSE,
  result_lock_authorized =
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
  supersedes_partial_execution =
    "V1_0_DATA_TABLE_GROUPED_SUMMARY_INTEGER_DOUBLE_TYPE_MISMATCH",
  dataset_id = dataset_id,
  integration_version =
    integration_version,
  pseudobulk_version =
    pseudobulk_version,
  model_version =
    model_version,
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
  primary_model_n =
    nrow(primary_models),
  nominal_exact_P_lt_0_05_n =
    nrow(nominal_models),
  global_BH_lt_0_05_n =
    nrow(global_BH_models),
  LOO_direction_unstable_n =
    sum(
      primary_models$
        LOO_review_flag ==
        "LOO_DIRECTION_UNSTABLE"
    ),
  cluster8_direction_changed_n =
    nrow(
      direction_changed_sensitivity
    ),
  independent_validation_review_authorized =
    all_gate_checks_pass,
  result_lock_authorized =
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
  "primary_models=",
  nrow(primary_models),
  "nominal=",
  nrow(nominal_models),
  "global_BH=",
  nrow(global_BH_models)
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
cat("STEP04E1O GSE183837 donor-level models completed\n")
cat("============================================================\n")
cat("Gate status: ", gate_status, "\n", sep = "")
cat(
  "Primary models: ",
  nrow(primary_models),
  "\n",
  sep = ""
)
cat(
  "Nonconstant testable models: ",
  sum(
    !primary_models$
      constant_endpoint
  ),
  "\n",
  sep = ""
)
cat(
  "Nominal exact P < 0.05: ",
  nrow(nominal_models),
  "\n",
  sep = ""
)
cat(
  "Global BH q < 0.05: ",
  nrow(global_BH_models),
  "\n",
  sep = ""
)
cat(
  "LOO direction-unstable models: ",
  sum(
    primary_models$
      LOO_review_flag ==
      "LOO_DIRECTION_UNSTABLE"
  ),
  "\n",
  sep = ""
)
cat(
  "Cluster-8 sensitivity direction changes: ",
  nrow(
    direction_changed_sensitivity
  ),
  "\n",
  sep = ""
)
cat(
  "Independent validation review authorized: ",
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
    "\nDonor-model outputs are complete in the project folders.\n",
    sep = ""
  )
}
