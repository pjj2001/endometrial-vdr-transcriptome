#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D01: Prepare independent single-cell inputs

options(stringsAsFactors = FALSE, warn = 1)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "prepare_gse183837_inputs.R"

script_version <- "STEP04E1I_GSE183837_INTEGRATION_PREFLIGHT_v1.0"

planned_integration_version <- "SCI183_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1I_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_path_safe <- function(path, max_chars = 235L) {
    path_n <- nchar(normalize_slashes(path), type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Path exceeds the Windows/R safety limit (", path_n, " > ", max_chars, "):\n", normalize_slashes(path))
    }
    invisible(path)
}

write_csv_utf8 <- function(x, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    out <- tryCatch(.vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), 
        error = function(e) NULL)
    if (is.null(out)) {
        out <- .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    out
}

safe_package_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        return(NA_character_)
    }
    as.character(utils::packageVersion(pkg))
}

detect_first_column <- function(column_names, candidates) {
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
    }
    else {
        paste(x, collapse = sep)
    }
}

get_free_memory_gb <- function() {
    if (.Platform$OS.type != "windows") {
        return(NA_real_)
    }
    command <- paste("(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory")
    output <- tryCatch(system2("powershell", args = c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, 
        stderr = FALSE), error = function(e) character())
    value <- suppressWarnings(as.numeric(trimws(output[1])))
    if (length(value) != 1L || !is.finite(value)) {
        NA_real_
    }
    else {
        round(value/1024^2, 2)
    }
}

get_free_disk_gb <- function() {
    if (.Platform$OS.type != "windows") {
        return(NA_real_)
    }
    command <- paste("(Get-CimInstance Win32_LogicalDisk", "-Filter \"DeviceID='C:'\").FreeSpace")
    output <- tryCatch(system2("powershell", args = c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, 
        stderr = FALSE), error = function(e) character())
    value <- suppressWarnings(as.numeric(trimws(output[1])))
    if (length(value) != 1L || !is.finite(value)) {
        NA_real_
    }
    else {
        round(value/1024^3, 2)
    }
}

log_path <- NULL

log_event <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(line, "\n")
    if (!is.null(log_path)) {
        cat(line, "\n", file = log_path, append = TRUE)
    }
    flush.console()
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

design_path <- file.path(results_root, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0", "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv")

preQC_root <- file.path(results_root, "S4D_sc_preQC", "SCP_v1_0", dataset_id)

planned_output_root <- file.path(results_root, "S4E1_sc", planned_integration_version)

planned_checkpoint_root <- file.path(planned_output_root, "checkpoints", dataset_id)

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1I_", format(run_start, "%Y%m%d_%H%M%S"), "_preflight"))

for (path_value in c(execution_root, parameter_root, environment_root, run_qc_root)) {
    dir.create(path_value, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(execution_root, paste0(run_id, "_execution.log"))

run_manifest_path <- file.path(execution_root, paste0(run_id, "_RUN_MANIFEST.csv"))

parameter_path <- file.path(parameter_root, paste0(run_id, "_parameters.json"))

session_path <- file.path(environment_root, paste0("sessionInfo_", run_id, ".txt"))

desktop_candidates <- .vdr_report_root()

desktop_candidates <- .vdr_report_root()

zip_parent <- .vdr_report_root()

handoff_zip <- file.path(zip_parent, paste0(run_id, "_SMALL_REPORT.zip"))

DE_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

DE_qc_dirs <- DE_qc_dirs[grepl("^STEP04E1H_[0-9]{8}_[0-9]{6}_deQC$", basename(DE_qc_dirs), perl = TRUE)]

if (!length(DE_qc_dirs)) {
    stop("No completed STEP04E1H gene-level DE QC directory found.")
}

latest_DE_qc <- DE_qc_dirs[which.max(file.info(DE_qc_dirs)$mtime)]

DE_gate_path <- file.path(latest_DE_qc, "GENE_DE_GATE.csv")

DE_gate <- safe_read_csv(DE_gate_path)

if (is.null(DE_gate) || nrow(DE_gate) != 1L || DE_gate$gate_status[[1]] != "PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE" || 
    !isTRUE(DE_gate$GSE183837_independent_localization_authorized[[1]])) {
    stop("STEP04E1H Gate does not authorize GSE183837 independent localization.")
}

if (!file.exists(design_path)) {
    stop("Locked GSE183837 sample design is missing:\n", design_path)
}

design <- safe_read_csv(design_path)

if (is.null(design) || !"sample_accession" %in% colnames(design) || nrow(design) != 9L || anyDuplicated(design$sample_accession) > 
    0L) {
    stop("Locked GSE183837 sample design failed identity validation.")
}

group_column <- detect_first_column(colnames(design), c("clinical_group", "group", "condition", "phenotype", "case_control_status", 
    "case_control"))

donor_column <- detect_first_column(colnames(design), c("donor_id", "donor", "subject_id", "patient_id"))

if (!nzchar(group_column)) {
    stop("Could not identify the Control/RIF group column in the locked design.")
}

group_value <- toupper(trimws(as.character(design[[group_column]])))

group_value[group_value %in% c("CTRL", "CONTROL", "FERTILE", "NORMAL")] <- "CONTROL"

group_value[group_value %in% c("RIF", "RECURRENT_IMPLANTATION_FAILURE")] <- "RIF"

group_counts <- table(factor(group_value, levels = c("CONTROL", "RIF")))

design_identity_pass <- identical(as.integer(group_counts), c(3L, 6L))

design_audit <- data.frame(metric = c("locked_sample_n", "unique_sample_n", "group_column", "donor_column", "CONTROL_n", 
    "RIF_n", "design_identity_pass"), value = c(nrow(design), length(unique(design$sample_accession)), group_column, 
    donor_column, as.integer(group_counts["CONTROL"]), as.integer(group_counts["RIF"]), design_identity_pass), 
    stringsAsFactors = FALSE)

write_csv_utf8(design_audit, file.path(run_qc_root, "01_design_audit.csv"))

if (!design_identity_pass) {
    stop("Locked GSE183837 design is not exactly 3 Control and 6 RIF donors.")
}

if (!dir.exists(preQC_root)) {
    stop("GSE183837 pre-QC directory is missing:\n", preQC_root)
}

preQC_rows <- list()

for (sample_index in seq_len(nrow(design))) {
    sample_accession <- as.character(design$sample_accession[sample_index])
    object_path <- file.path(preQC_root, paste0(sample_accession, "_preQC.rds"))
    if (!file.exists(object_path)) {
        preQC_rows[[sample_index]] <- data.frame(sample_accession = sample_accession, object_exists = FALSE, object_size_MB = NA_real_, 
            feature_n = NA_integer_, retained_cell_n = NA_integer_, feature_table_row_n = NA_integer_, cell_metadata_row_n = NA_integer_, 
            required_components_pass = FALSE, matrix_dimnames_pass = FALSE, feature_alignment_pass = FALSE, cell_alignment_pass = FALSE, 
            sample_identity_pass = FALSE, md5 = NA_character_, object_path = normalize_slashes(object_path), stringsAsFactors = FALSE)
        next
    }
    object_value <- tryCatch(readRDS(object_path), error = function(e) NULL)
    required_components <- c("counts", "feature_table", "cell_metadata")
    components_pass <- (!is.null(object_value) && all(required_components %in% names(object_value)))
    if (!components_pass) {
        preQC_rows[[sample_index]] <- data.frame(sample_accession = sample_accession, object_exists = TRUE, object_size_MB = round(file.info(object_path)$size/1024^2, 
            2), feature_n = NA_integer_, retained_cell_n = NA_integer_, feature_table_row_n = NA_integer_, cell_metadata_row_n = NA_integer_, 
            required_components_pass = FALSE, matrix_dimnames_pass = FALSE, feature_alignment_pass = FALSE, cell_alignment_pass = FALSE, 
            sample_identity_pass = FALSE, md5 = unname(tools::md5sum(object_path)), object_path = normalize_slashes(object_path), 
            stringsAsFactors = FALSE)
        next
    }
    counts <- object_value$counts
    feature_table <- object_value$feature_table
    cell_metadata <- object_value$cell_metadata
    matrix_like <- (inherits(counts, "Matrix") || is.matrix(counts))
    matrix_dimnames_pass <- (matrix_like && !is.null(rownames(counts)) && !is.null(colnames(counts)) && !anyDuplicated(colnames(counts)))
    feature_alignment_pass <- (is.data.frame(feature_table) && nrow(feature_table) == nrow(counts))
    cell_alignment_pass <- (is.data.frame(cell_metadata) && nrow(cell_metadata) == ncol(counts))
    metadata_sample_pass <- TRUE
    if ("sample_accession" %in% colnames(cell_metadata)) {
        metadata_sample_pass <- (length(unique(as.character(cell_metadata$sample_accession))) == 1L && unique(as.character(cell_metadata$sample_accession)) == 
            sample_accession)
    }
    preQC_rows[[sample_index]] <- data.frame(sample_accession = sample_accession, object_exists = TRUE, object_size_MB = round(file.info(object_path)$size/1024^2, 
        2), feature_n = nrow(counts), retained_cell_n = ncol(counts), feature_table_row_n = nrow(feature_table), 
        cell_metadata_row_n = nrow(cell_metadata), required_components_pass = components_pass, matrix_dimnames_pass = matrix_dimnames_pass, 
        feature_alignment_pass = feature_alignment_pass, cell_alignment_pass = cell_alignment_pass, sample_identity_pass = metadata_sample_pass, 
        md5 = unname(tools::md5sum(object_path)), object_path = normalize_slashes(object_path), stringsAsFactors = FALSE)
    rm(object_value, counts, feature_table, cell_metadata)
    gc()
}

preQC_audit <- do.call(rbind, preQC_rows)

write_csv_utf8(preQC_audit, file.path(run_qc_root, "02_preQC_object_audit.csv"))

preQC_total_cells <- sum(preQC_audit$retained_cell_n, na.rm = TRUE)

preQC_all_pass <- (nrow(preQC_audit) == 9L && all(preQC_audit$object_exists) && all(preQC_audit$required_components_pass) && 
    all(preQC_audit$matrix_dimnames_pass) && all(preQC_audit$feature_alignment_pass) && all(preQC_audit$cell_alignment_pass) && 
    all(preQC_audit$sample_identity_pass) && preQC_total_cells == 50997L)

package_names <- c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "jsonlite", "future", "scDblFinder", 
    "SingleCellExperiment", "SummarizedExperiment", "zip")

package_audit <- data.frame(package = package_names, installed = vapply(package_names, requireNamespace, logical(1), 
    quietly = TRUE), version = vapply(package_names, safe_package_version, character(1)), required_for_integration = c(TRUE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE), stringsAsFactors = FALSE)

write_csv_utf8(package_audit, file.path(run_qc_root, "03_package_audit.csv"))

required_packages_pass <- all(package_audit$installed[package_audit$required_for_integration])

existing_GSE183837_files <- list.files(file.path(results_root, "S4E1_sc"), pattern = "GSE183837|183837", recursive = TRUE, 
    full.names = TRUE, ignore.case = TRUE)

checkpoint_inventory <- data.frame(path = normalize_slashes(existing_GSE183837_files), size_MB = if (length(existing_GSE183837_files)) {
    round(file.info(existing_GSE183837_files)$size/1024^2, 2)
} else {
    numeric()
}, modified_time = if (length(existing_GSE183837_files)) {
    as.character(file.info(existing_GSE183837_files)$mtime)
} else {
    character()
}, stringsAsFactors = FALSE)

write_csv_utf8(checkpoint_inventory, file.path(run_qc_root, "04_existing_output_inventory.csv"))

planned_paths <- data.frame(artifact = c("planned_output_root", "planned_checkpoint_root", "planned_integrated_checkpoint", 
    "planned_checkpoint_status", "planned_UMAP_metadata", "planned_doublet_QC"), path = normalize_slashes(c(planned_output_root, 
    planned_checkpoint_root, file.path(planned_checkpoint_root, "integrated.rds"), file.path(planned_checkpoint_root, 
        "status.csv"), file.path(planned_output_root, "tables", "UMAP_metadata.csv"), file.path(planned_output_root, 
        "tables", "doublet_QC.csv"))), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

write_csv_utf8(planned_paths, file.path(run_qc_root, "05_path_audit.csv"))

resource_audit <- data.frame(metric = c("free_physical_memory_GB", "free_C_drive_GB", "preQC_object_total_size_GB", 
    "preQC_total_cells", "preQC_sample_n"), value = c(get_free_memory_gb(), get_free_disk_gb(), round(sum(preQC_audit$object_size_MB)/1024, 
    3), preQC_total_cells, nrow(preQC_audit)), stringsAsFactors = FALSE)

write_csv_utf8(resource_audit, file.path(run_qc_root, "06_resource_audit.csv"))

integration_plan <- data.frame(parameter = c("dataset_processing", "biological_replicate", "doublet_method", "doublet_scope", 
    "normalization", "variable_features", "PCA_components", "Harmony_batch_variable", "Harmony_dimensions", "UMAP_dimensions", 
    "cluster_resolutions", "primary_resolution", "counts_layer_policy", "data_layer_policy", "scale_data_policy", 
    "annotation_policy", "disease_testing_policy"), value = c("GSE183837 processed independently from GSE250130", 
    "donor/sample", "scDblFinder", "run independently within each donor", "LogNormalize scale.factor=10000", "3000", 
    "50", "sample_accession", "1:30", "1:30", "0.2;0.4;0.6", "0.4", "retain 9 split raw-count layers", "join normalized data into one layer after integration", 
    "remove scale.data before final checkpoint", "cluster markers plus canonical panels and GSE250130 broad reference", 
    "no cell-level test; donor-level pseudobulk only after annotation QC"), stringsAsFactors = FALSE)

write_csv_utf8(integration_plan, file.path(run_qc_root, "07_integration_plan.csv"))

memory_value <- suppressWarnings(as.numeric(resource_audit$value[resource_audit$metric == "free_physical_memory_GB"]))

disk_value <- suppressWarnings(as.numeric(resource_audit$value[resource_audit$metric == "free_C_drive_GB"]))

memory_pass <- (is.na(memory_value) || memory_value >= 12)

disk_pass <- (is.na(disk_value) || disk_value >= 15)

gate_checks <- data.frame(check = c("STEP04E1H_authorizes_GSE183837", "locked_design_has_9_donors", "locked_design_is_3_Control_6_RIF", 
    "nine_preQC_objects_present", "all_preQC_object_structures_valid", "preQC_total_cells_is_50997", "required_integration_packages_installed", 
    "planned_paths_safe", "free_memory_at_least_12GB_or_unknown", "free_disk_at_least_15GB_or_unknown"), pass = c(DE_gate$GSE183837_independent_localization_authorized[[1]], 
    nrow(design) == 9L, design_identity_pass, nrow(preQC_audit) == 9L && all(preQC_audit$object_exists), all(preQC_audit$required_components_pass & 
        preQC_audit$matrix_dimnames_pass & preQC_audit$feature_alignment_pass & preQC_audit$cell_alignment_pass & 
        preQC_audit$sample_identity_pass), preQC_total_cells == 50997L, required_packages_pass, all(planned_paths$path_safe), 
    memory_pass, disk_pass), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "08_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_GSE183837_INTEGRATION_PREFLIGHT"
} else {
    "HOLD_GSE183837_INTEGRATION_PREFLIGHT_FAILED"
}

preflight_gate <- data.frame(dataset_id = dataset_id, script_version = script_version, planned_integration_version = planned_integration_version, 
    gate_status = gate_status, locked_sample_n = nrow(design), Control_n = as.integer(group_counts["CONTROL"]), 
    RIF_n = as.integer(group_counts["RIF"]), preQC_cell_n = preQC_total_cells, preQC_object_n = nrow(preQC_audit), 
    required_packages_pass = required_packages_pass, integration_authorized = all_gate_checks_pass, automatic_annotation_authorized = FALSE, 
    disease_pseudobulk_testing_authorized = FALSE, result_lock_authorized = FALSE, stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "GSE183837_PREFLIGHT_GATE.csv")

write_csv_utf8(preflight_gate, gate_path)

card_lines <- c("# STEP04E1I GSE183837 Integration Preflight Card", "", paste0("- Dataset: ", dataset_id), paste0("- Planned integration version: ", 
    planned_integration_version), paste0("- Gate: ", gate_status), paste0("- Locked donors: ", nrow(design), " (Control ", 
    as.integer(group_counts["CONTROL"]), "; RIF ", as.integer(group_counts["RIF"]), ")"), paste0("- Donor-level pre-QC objects: ", 
    nrow(preQC_audit)), paste0("- Pre-QC retained cells: ", format(preQC_total_cells, big.mark = ",")), paste0("- Required packages complete: ", 
    required_packages_pass), paste0("- Integration authorized: ", all_gate_checks_pass), "- Automatic annotation authorized: FALSE", 
    "- Disease pseudobulk testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "GSE183837_PREFLIGHT_CARD.md"))

interpretation_lines <- c("# GSE183837 integration preflight interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), "", "## Authorized next operation", "", "- Run sample-wise scDblFinder on the nine pre-QC donor objects.", 
    "- Integrate GSE183837 independently using sample accession as the Harmony batch variable.", "- Save staged checkpoints before marker extraction or biological annotation.", 
    "", "## Not authorized yet", "", "- No automatic broad-cell annotation is locked at this stage.", "- No cell-level Control-versus-RIF test is permitted.", 
    "- No donor-level pseudobulk disease test is permitted until broad-cell annotation QC passes.", "- No GSE250130 and GSE183837 cell pooling is permitted.")

write_text(interpretation_lines, file.path(run_qc_root, "GSE183837_PREFLIGHT_INTERPRETATION.md"))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to write the parameter manifest.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, dataset_id = dataset_id, 
    upstream_gate = DE_gate$gate_status[[1]], planned_integration_version = planned_integration_version, expected_sample_n = 9L, 
    expected_Control_n = 3L, expected_RIF_n = 6L, expected_preQC_cell_n = 50997L, minimum_free_memory_GB = 12, 
    minimum_free_disk_GB = 15, integration_plan = as.list(setNames(integration_plan$value, integration_plan$parameter)), 
    result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, dataset_id = dataset_id, 
    started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_seconds = round(as.numeric(difftime(run_end, 
        run_start, units = "secs")), 2), gate_status = gate_status, locked_sample_n = nrow(design), preQC_object_n = nrow(preQC_audit), 
    preQC_cell_n = preQC_total_cells, required_packages_pass = required_packages_pass, integration_authorized = all_gate_checks_pass, 
    planned_output_root = normalize_slashes(planned_output_root), qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

artifact_files <- unique(c(list.files(run_qc_root, recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, 
    parameter_path, session_path))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "preQC_objects=", nrow(preQC_audit), "preQC_cells=", 
    preQC_total_cells)

small_files <- unique(c(list.files(run_qc_root, recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, 
    parameter_path, session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files)]

if (!requireNamespace("zip", quietly = TRUE)) {
    try(install.packages("zip", repos = "https://cloud.r-project.org"), silent = TRUE)
}

zip_ok <- FALSE

zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
    zip_ok <- tryCatch({
        zip::zipr(zipfile = handoff_zip, files = small_files, recurse = FALSE, include_directories = FALSE)
        file.exists(handoff_zip)
    }, error = function(e) {
        zip_error <<- conditionMessage(e)
        FALSE
    })
}

cat("\n============================================================\n")

cat("STEP04E1I GSE183837 integration preflight completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Locked donors: ", nrow(design), " (Control ", as.integer(group_counts["CONTROL"]), "; RIF ", as.integer(group_counts["RIF"]), 
    ")\n", sep = "")

cat("Pre-QC objects: ", nrow(preQC_audit), "\n", sep = "")

cat("Pre-QC retained cells: ", format(preQC_total_cells, big.mark = ","), "\n", sep = "")

cat("Integration authorized: ", all_gate_checks_pass, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe preflight report is complete in the QC folder.\n", sep = "")
}

