#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D06: Aggregate primary and atlas-sensitivity pseudobulks

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

options(future.globals.maxSize = 12 * 1024^3)

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "aggregate_gse183837_pseudobulk.R"

script_version <- "STEP04E1N_GSE183837_VDR_LOCALIZATION_PSEUDOBULK_v1.1_UNIQUE_COLUMNS_FIXED"

integration_version <- "SCI183_v1_0"

annotation_version <- "ANN183_v1_0"

pseudobulk_version <- "PB183_v1_0"

minimum_cells_primary <- 20L

expected_primary_cell_n <- 45630L

expected_atlas_cell_n <- 47161L

expected_cluster8_cell_n <- 1531L

run_start <- Sys.time()

run_id <- paste0("STEP04E1N_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_path_safe <- function(path, max_chars = 235L) {
    normalized <- normalize_slashes(path)
    path_n <- nchar(normalized, type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Path exceeds the Windows/R safety limit (", path_n, " > ", max_chars, "):\n", normalized)
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
    writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    output <- tryCatch(.vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), 
        error = function(e) NULL)
    if (is.null(output)) {
        output <- .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    output
}

as_logical_safe <- function(x) {
    if (is.logical(x)) 
        return(x)
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
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

sanitize_token <- function(x) {
    output <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    output <- gsub("^_+|_+$", "", output)
    substr(output, 1L, 60L)
}

row_zscore_reference <- function(x) {
    x <- as.matrix(x)
    row_mean <- rowMeans(x, na.rm = TRUE)
    row_sd <- apply(x, 1L, stats::sd, na.rm = TRUE)
    list(mean = row_mean, sd = row_sd, valid = is.finite(row_sd) & row_sd > 0)
}

apply_row_reference_z <- function(x, reference_mean, reference_sd) {
    x <- as.matrix(x)
    if (!identical(rownames(x), names(reference_mean)) || !identical(rownames(x), names(reference_sd))) {
        stop("Expression matrix and row-standardization reference are not aligned.")
    }
    output <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    valid <- is.finite(reference_sd) & reference_sd > 0
    output[valid, ] <- sweep(sweep(x[valid, , drop = FALSE], MARGIN = 1, STATS = reference_mean[valid], FUN = "-"), 
        MARGIN = 1, STATS = reference_sd[valid], FUN = "/")
    output
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

save_plot_pair <- function(plot_object, stem, width, height) {
    png_path <- file.path(figure_root, paste0(stem, ".png"))
    pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))
    assert_path_safe(png_path)
    assert_path_safe(pdf_path)
    ggplot2::ggsave(filename = png_path, plot = plot_object, width = width, height = height, units = "in", dpi = 300, 
        bg = "white")
    ggplot2::ggsave(filename = pdf_path, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    c(png_path, pdf_path)
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

required_packages <- c("SeuratObject", "Matrix", "data.table", "ggplot2", "scales", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

integration_root <- file.path(results_root, "S4E1_sc", integration_version)

annotated_checkpoint <- file.path(integration_root, "ann", "G183", annotation_version, "annotated.rds")

annotation_status_path <- file.path(integration_root, "ann", "G183", annotation_version, "status.csv")

candidate_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

output_root <- file.path(integration_root, "pb", "G183", pseudobulk_version)

layer_root <- file.path(output_root, "lyr")

object_root <- file.path(output_root, "obj")

table_root <- file.path(output_root, "tab")

source_data_root <- file.path(output_root, "src")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "man")

primary_counts_path <- file.path(object_root, "pb_primary_counts.rds")

atlas_counts_path <- file.path(object_root, "pb_atlas_counts.rds")

pb_status_path <- file.path(object_root, "pb_status.csv")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1N_", format(run_start, "%Y%m%d_%H%M%S"), "_pbQC"))

figure_root <- file.path(run_qc_root, "fig")

for (path_value in c(execution_root, parameter_root, environment_root, manuscript_text_root, output_root, layer_root, 
    object_root, table_root, source_data_root, figure_output_root, manifest_root, run_qc_root, figure_root)) {
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

planned_paths <- data.frame(artifact = c("annotated_checkpoint", "annotation_status", "candidate_table", "output_root", 
    "layer_root", "primary_counts", "atlas_counts", "pb_status", "run_qc_root", "handoff_zip"), path = normalize_slashes(c(annotated_checkpoint, 
    annotation_status_path, candidate_path, output_root, layer_root, primary_counts_path, atlas_counts_path, pb_status_path, 
    run_qc_root, handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

if (!all(planned_paths$path_safe)) {
    stop("One or more pseudobulk paths exceed the safe limit.")
}

annotation_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

annotation_qc_dirs <- annotation_qc_dirs[grepl("^STEP04E1M_[0-9]{8}_[0-9]{6}_annQC$", basename(annotation_qc_dirs), 
    perl = TRUE)]

if (!length(annotation_qc_dirs)) {
    stop("No STEP04E1M annotation-QC directory was found.")
}

gate_candidates <- list()

for (directory_path in annotation_qc_dirs) {
    gate_path <- file.path(directory_path, "G183_ANNOTATION_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        gate_candidates[[length(gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(gate_candidates)) {
    stop("No readable annotation Gate was found.")
}

gate_candidates <- as.data.frame(data.table::rbindlist(gate_candidates, use.names = TRUE, fill = TRUE))

gate_candidates$source_mtime_POSIX <- as.POSIXct(gate_candidates$source_mtime)

gate_candidates <- gate_candidates[order(gate_candidates$source_mtime_POSIX, decreasing = TRUE), , drop = FALSE]

authorized_index <- which(gate_candidates$gate_status == "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION" & 
    toupper(as.character(gate_candidates$broad_annotation_locked)) %in% c("TRUE", "T", "1", "YES") & toupper(as.character(gate_candidates$VDR_localization_authorized)) %in% 
    c("TRUE", "T", "1", "YES") & toupper(as.character(gate_candidates$broad_pseudobulk_preparation_authorized)) %in% 
    c("TRUE", "T", "1", "YES"))

if (!length(authorized_index)) {
    stop("No annotation Gate authorizes VDR localization and pseudobulk preparation.")
}

upstream_gate <- gate_candidates[authorized_index[1], , drop = FALSE]

write_csv_utf8(upstream_gate, file.path(run_qc_root, "00_upstream_annotation_gate.csv"))

required_inputs <- c(annotated_checkpoint, annotation_status_path, candidate_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(normalize_slashes(missing_inputs), collapse = "\n"))
}

annotation_status <- safe_read_csv(annotation_status_path)

candidate_table <- safe_read_csv(candidate_path)

if (is.null(annotation_status) || nrow(annotation_status) != 1L || annotation_status$status[[1]] != "COMPLETE" || 
    annotation_status$cell_n[[1]] != expected_atlas_cell_n) {
    stop("Annotated checkpoint status is absent or invalid.")
}

required_candidate_columns <- c("feature_id", "interaction_logFC", "occupancy_supported", "shared_occupancy_supported", 
    "any_promoter_overlap_support")

missing_candidate_columns <- setdiff(required_candidate_columns, colnames(candidate_table))

if (length(missing_candidate_columns)) {
    stop("Locked VDR candidate table is missing: ", paste(missing_candidate_columns, collapse = ", "))
}

candidate_table$feature_id <- as.character(candidate_table$feature_id)

candidate_table$interaction_logFC <- suppressWarnings(as.numeric(candidate_table$interaction_logFC))

for (logical_column in c("occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support")) {
    candidate_table[[logical_column]] <- as_logical_safe(candidate_table[[logical_column]])
}

if (nrow(candidate_table) != 436L || anyDuplicated(candidate_table$feature_id) > 0L || sum(candidate_table$occupancy_supported) != 
    284L || sum(candidate_table$shared_occupancy_supported) != 46L || sum(candidate_table$any_promoter_overlap_support) != 
    84L) {
    stop("Locked VDR program identity does not match 436/284/46/84.")
}

free_memory_before_GB <- get_free_memory_gb()

free_disk_before_GB <- get_free_disk_gb()

if (is.finite(free_memory_before_GB) && free_memory_before_GB < 10) {
    stop("Free physical memory is below 10 GB: ", free_memory_before_GB)
}

log_event("Starting", run_id, "free_memory_GB=", free_memory_before_GB, "free_disk_GB=", free_disk_before_GB)

annotated_object <- readRDS(annotated_checkpoint)

if (!inherits(annotated_object, "Seurat")) {
    stop("Annotated checkpoint is not a Seurat object.")
}

SeuratObject::DefaultAssay(annotated_object) <- "RNA"

metadata <- annotated_object[[]]

required_metadata_columns <- c("sample_accession", "donor_id", "analysis_group", "primary_cluster", "broad_cell_type", 
    "annotation_confidence", "primary_VDR_localization_eligible", "primary_broad_pseudobulk_eligible")

missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata))

if (length(missing_metadata_columns)) {
    stop("Annotated metadata is missing: ", paste(missing_metadata_columns, collapse = ", "))
}

metadata$sample_accession <- as.character(metadata$sample_accession)

metadata$donor_id <- as.character(metadata$donor_id)

metadata$analysis_group <- as.character(metadata$analysis_group)

metadata$primary_cluster <- as.character(metadata$primary_cluster)

metadata$broad_cell_type <- as.character(metadata$broad_cell_type)

metadata$primary_VDR_localization_eligible <- as_logical_safe(metadata$primary_VDR_localization_eligible)

metadata$primary_broad_pseudobulk_eligible <- as_logical_safe(metadata$primary_broad_pseudobulk_eligible)

if (ncol(annotated_object) != expected_atlas_cell_n || nrow(metadata) != expected_atlas_cell_n || length(unique(metadata$sample_accession)) != 
    9L || length(unique(metadata$donor_id)) != 9L || sum(metadata$primary_broad_pseudobulk_eligible) != expected_primary_cell_n || 
    sum(metadata$primary_cluster == "8") != expected_cluster8_cell_n) {
    stop("Annotated object identity or cluster-8 exclusion count is inconsistent.")
}

donor_group_identity <- unique(metadata[, c("sample_accession", "donor_id", "analysis_group"), drop = FALSE])

donor_group_identity <- donor_group_identity[order(donor_group_identity$sample_accession), , drop = FALSE]

if (nrow(donor_group_identity) != 9L || sum(donor_group_identity$analysis_group == "CONTROL") != 3L || sum(donor_group_identity$analysis_group == 
    "RIF") != 6L || anyDuplicated(donor_group_identity$donor_id) > 0L) {
    stop("Donor/group identity does not match 3 Control and 6 RIF donors.")
}

count_layers <- SeuratObject::Layers(annotated_object, assay = "RNA", search = "^counts")

data_layers <- SeuratObject::Layers(annotated_object, assay = "RNA", search = "^data")

if (length(count_layers) != 9L || length(data_layers) != 1L || data_layers[1] != "data") {
    stop("Expected nine split count layers and one joined data layer.")
}

write_csv_utf8(donor_group_identity, file.path(table_root, "01_donor_group_identity.csv"))

layer_status_rows <- list()

layer_result_paths <- character(length(count_layers))

for (layer_index in seq_along(count_layers)) {
    layer_name <- count_layers[layer_index]
    layer_code <- sprintf("L%02d", layer_index)
    layer_result_path <- file.path(layer_root, paste0(layer_code, ".rds"))
    layer_status_path <- file.path(layer_root, paste0(layer_code, "_status.csv"))
    layer_result_paths[layer_index] <- layer_result_path
    existing_status <- safe_read_csv(layer_status_path)
    layer_result <- NULL
    resume_valid <- FALSE
    if (file.exists(layer_result_path) && !is.null(existing_status) && nrow(existing_status) == 1L && existing_status$status[[1]] == 
        "COMPLETE" && existing_status$layer_name[[1]] == layer_name) {
        layer_result <- tryCatch(readRDS(layer_result_path), error = function(e) NULL)
        resume_valid <- (is.list(layer_result) && all(c("primary_counts", "atlas_counts", "primary_metadata", "atlas_metadata", 
            "cell_audit", "raw_count_audit") %in% names(layer_result)))
    }
    if (resume_valid) {
        log_event("SKIP completed counts layer", layer_index, "/", length(count_layers), layer_name)
    }
    else {
        layer_start <- Sys.time()
        log_event("START counts layer", layer_index, "/", length(count_layers), layer_name)
        layer_counts <- SeuratObject::LayerData(object = annotated_object, assay = "RNA", layer = layer_name)
        layer_cells <- colnames(layer_counts)
        if (is.null(layer_cells) || !length(layer_cells) || !all(layer_cells %in% rownames(metadata))) {
            stop("Counts-layer cell identities are missing or invalid for ", layer_name, ".")
        }
        layer_metadata <- metadata[layer_cells, , drop = FALSE]
        layer_samples <- unique(layer_metadata$sample_accession)
        layer_donors <- unique(layer_metadata$donor_id)
        layer_groups <- unique(layer_metadata$analysis_group)
        if (length(layer_samples) != 1L || length(layer_donors) != 1L || length(layer_groups) != 1L) {
            stop("Counts layer does not map to exactly one sample/donor/group: ", layer_name)
        }
        aggregate_by_type <- function(counts_matrix, metadata_table, inclusion_vector, matrix_role) {
            type_values <- sort(unique(metadata_table$broad_cell_type[inclusion_vector]))
            if (!length(type_values)) {
                stop("No cells remain for ", matrix_role, " in ", layer_name, ".")
            }
            count_columns <- list()
            metadata_rows <- list()
            for (type_index in seq_along(type_values)) {
                cell_type <- type_values[type_index]
                selected_index <- which(inclusion_vector & metadata_table$broad_cell_type == cell_type)
                if (!length(selected_index)) {
                  next
                }
                pb_id <- paste(dataset_id, layer_samples, sanitize_token(cell_type), matrix_role, sep = "__")
                aggregated_count <- Matrix::rowSums(counts_matrix[, selected_index, drop = FALSE])
                count_columns[[length(count_columns) + 1L]] <- Matrix::Matrix(aggregated_count, ncol = 1L, sparse = TRUE, 
                  dimnames = list(rownames(counts_matrix), pb_id))
                metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(pb_id = pb_id, matrix_role = matrix_role, 
                  sample_accession = layer_samples, donor_id = layer_donors, analysis_group = layer_groups, broad_cell_type = cell_type, 
                  cell_n = length(selected_index), library_size = as.numeric(sum(aggregated_count)), minimum_cell_threshold = minimum_cells_primary, 
                  minimum_cell_threshold_pass = length(selected_index) >= minimum_cells_primary, stringsAsFactors = FALSE)
            }
            list(counts = do.call(cbind, count_columns), metadata = do.call(rbind, metadata_rows))
        }
        primary_result <- aggregate_by_type(counts_matrix = layer_counts, metadata_table = layer_metadata, inclusion_vector = layer_metadata$primary_broad_pseudobulk_eligible, 
            matrix_role = "PRIMARY")
        atlas_result <- aggregate_by_type(counts_matrix = layer_counts, metadata_table = layer_metadata, inclusion_vector = rep(TRUE, 
            nrow(layer_metadata)), matrix_role = "ATLAS")
        VDR_index <- match("VDR", toupper(rownames(layer_counts)))
        if (is.na(VDR_index)) {
            stop("VDR feature is absent from counts layer ", layer_name, ".")
        }
        VDR_cell_rows <- list()
        for (cell_type in sort(unique(layer_metadata$broad_cell_type))) {
            type_index <- which(layer_metadata$broad_cell_type == cell_type)
            primary_index <- type_index[layer_metadata$primary_VDR_localization_eligible[type_index]]
            atlas_VDR_counts <- as.numeric(layer_counts[VDR_index, type_index, drop = TRUE])
            primary_VDR_counts <- if (length(primary_index)) {
                as.numeric(layer_counts[VDR_index, primary_index, drop = TRUE])
            }
            else {
                numeric()
            }
            VDR_cell_rows[[length(VDR_cell_rows) + 1L]] <- data.frame(sample_accession = layer_samples, donor_id = layer_donors, 
                analysis_group = layer_groups, broad_cell_type = cell_type, atlas_cell_n = length(type_index), 
                atlas_VDR_positive_n = sum(atlas_VDR_counts > 0), atlas_VDR_positive_fraction = mean(atlas_VDR_counts > 
                  0), primary_cell_n = length(primary_index), primary_VDR_positive_n = sum(primary_VDR_counts > 
                  0), primary_VDR_positive_fraction = if (length(primary_index)) {
                  mean(primary_VDR_counts > 0)
                }
                else {
                  NA_real_
                }, stringsAsFactors = FALSE)
        }
        VDR_cell_audit <- do.call(rbind, VDR_cell_rows)
        raw_count_audit <- data.frame(layer_name = layer_name, sample_accession = layer_samples, donor_id = layer_donors, 
            analysis_group = layer_groups, feature_n = nrow(layer_counts), atlas_cell_n = ncol(layer_counts), primary_cell_n = sum(layer_metadata$primary_broad_pseudobulk_eligible), 
            cluster8_excluded_cell_n = sum(layer_metadata$primary_cluster == "8"), raw_count_sum = as.numeric(sum(layer_counts)), 
            atlas_aggregated_count_sum = as.numeric(sum(atlas_result$counts)), primary_input_count_sum = as.numeric(sum(layer_counts[, 
                layer_metadata$primary_broad_pseudobulk_eligible, drop = FALSE])), primary_aggregated_count_sum = as.numeric(sum(primary_result$counts)), 
            stringsAsFactors = FALSE)
        if (raw_count_audit$raw_count_sum != raw_count_audit$atlas_aggregated_count_sum || raw_count_audit$primary_input_count_sum != 
            raw_count_audit$primary_aggregated_count_sum) {
            stop("Raw-count conservation failed for ", layer_name, ".")
        }
        layer_result <- list(layer_name = layer_name, primary_counts = primary_result$counts, atlas_counts = atlas_result$counts, 
            primary_metadata = primary_result$metadata, atlas_metadata = atlas_result$metadata, VDR_cell_audit = VDR_cell_audit, 
            cell_audit = data.frame(layer_name = layer_name, sample_accession = layer_samples, atlas_cell_n = nrow(layer_metadata), 
                primary_cell_n = sum(layer_metadata$primary_broad_pseudobulk_eligible), cluster8_excluded_cell_n = sum(layer_metadata$primary_cluster == 
                  "8"), stringsAsFactors = FALSE), raw_count_audit = raw_count_audit)
        saveRDS(layer_result, layer_result_path, compress = "gzip")
        layer_status <- data.frame(layer_name = layer_name, status = "COMPLETE", sample_accession = layer_samples, 
            donor_id = layer_donors, atlas_cell_n = nrow(layer_metadata), primary_cell_n = sum(layer_metadata$primary_broad_pseudobulk_eligible), 
            cluster8_excluded_cell_n = sum(layer_metadata$primary_cluster == "8"), primary_profile_n = ncol(primary_result$counts), 
            atlas_profile_n = ncol(atlas_result$counts), elapsed_minutes = round(as.numeric(difftime(Sys.time(), 
                layer_start, units = "mins")), 2), checkpoint_md5 = unname(tools::md5sum(layer_result_path)), stringsAsFactors = FALSE)
        write_csv_utf8(layer_status, layer_status_path)
        existing_status <- layer_status
        log_event("DONE counts layer", layer_index, "/", length(count_layers), layer_name, "primary_cells=", layer_status$primary_cell_n, 
            "atlas_cells=", layer_status$atlas_cell_n)
        rm(layer_counts, primary_result, atlas_result)
        gc()
    }
    layer_status_rows[[layer_index]] <- existing_status
}

layer_status_all <- do.call(rbind, layer_status_rows)

write_csv_utf8(layer_status_all, file.path(table_root, "02_layer_status.csv"))

layer_results <- lapply(layer_result_paths, readRDS)

feature_identities <- lapply(layer_results, function(x) {
    rownames(x$primary_counts)
})

if (!all(vapply(feature_identities, identical, logical(1), feature_identities[[1]]))) {
    stop("Per-layer pseudobulk feature identities are inconsistent.")
}

primary_counts <- do.call(cbind, lapply(layer_results, function(x) {
    x$primary_counts
}))

atlas_counts <- do.call(cbind, lapply(layer_results, function(x) {
    x$atlas_counts
}))

primary_metadata <- do.call(rbind, lapply(layer_results, function(x) {
    x$primary_metadata
}))

atlas_metadata <- do.call(rbind, lapply(layer_results, function(x) {
    x$atlas_metadata
}))

VDR_cell_audit <- do.call(rbind, lapply(layer_results, function(x) {
    x$VDR_cell_audit
}))

cell_audit <- do.call(rbind, lapply(layer_results, function(x) {
    x$cell_audit
}))

raw_count_audit <- do.call(rbind, lapply(layer_results, function(x) {
    x$raw_count_audit
}))

if (anyDuplicated(colnames(primary_counts)) > 0L || anyDuplicated(colnames(atlas_counts)) > 0L || !identical(colnames(primary_counts), 
    primary_metadata$pb_id) || !identical(colnames(atlas_counts), atlas_metadata$pb_id)) {
    stop("Combined pseudobulk columns and metadata are not aligned.")
}

if (sum(primary_metadata$cell_n) != expected_primary_cell_n || sum(atlas_metadata$cell_n) != expected_atlas_cell_n || 
    sum(cell_audit$cluster8_excluded_cell_n) != expected_cluster8_cell_n || ncol(primary_counts) != 54L || ncol(atlas_counts) != 
    54L) {
    stop("Combined pseudobulk identity does not match 45,630 primary cells, 47,161 atlas cells and 54 profiles.")
}

saveRDS(primary_counts, primary_counts_path, compress = "gzip")

saveRDS(atlas_counts, atlas_counts_path, compress = "gzip")

write_csv_utf8(primary_metadata, file.path(table_root, "03_primary_pb_metadata.csv"))

write_csv_utf8(atlas_metadata, file.path(table_root, "04_atlas_pb_metadata.csv"))

write_csv_utf8(VDR_cell_audit, file.path(table_root, "05_VDR_cell_detection.csv"))

write_csv_utf8(cell_audit, file.path(table_root, "06_cell_aggregation_audit.csv"))

write_csv_utf8(raw_count_audit, file.path(table_root, "07_raw_count_conservation.csv"))

pb_status <- data.frame(pseudobulk_version = pseudobulk_version, status = "COMPLETE", feature_n = nrow(primary_counts), 
    primary_profile_n = ncol(primary_counts), atlas_profile_n = ncol(atlas_counts), primary_cell_n = sum(primary_metadata$cell_n), 
    atlas_cell_n = sum(atlas_metadata$cell_n), cluster8_excluded_cell_n = sum(cell_audit$cluster8_excluded_cell_n), 
    count_layer_n = length(count_layers), minimum_cells_primary = minimum_cells_primary, primary_counts_path = normalize_slashes(primary_counts_path), 
    primary_counts_md5 = unname(tools::md5sum(primary_counts_path)), atlas_counts_path = normalize_slashes(atlas_counts_path), 
    atlas_counts_md5 = unname(tools::md5sum(atlas_counts_path)), stringsAsFactors = FALSE)

write_csv_utf8(pb_status, pb_status_path)

rm(annotated_object, metadata, layer_results)

gc()

feature_names <- rownames(primary_counts)

candidate_match <- match(toupper(candidate_table$feature_id), toupper(feature_names))

candidate_table$resolved_feature <- ifelse(is.na(candidate_match), NA_character_, feature_names[candidate_match])

candidate_table$feature_available <- !is.na(candidate_match)

VDR_feature_index <- match("VDR", toupper(feature_names))

if (is.na(VDR_feature_index)) {
    stop("VDR is absent from pseudobulk counts.")
}

VDR_feature <- feature_names[VDR_feature_index]

program_definitions <- list(RNA_ALL_436 = candidate_table$feature_id, OCCUPANCY_SUPPORTED_284 = candidate_table$feature_id[candidate_table$occupancy_supported], 
    SHARED_OCCUPANCY_46 = candidate_table$feature_id[candidate_table$shared_occupancy_supported], PROMOTER_OCCUPANCY_84 = candidate_table$feature_id[candidate_table$any_promoter_overlap_support])

program_sizes <- vapply(program_definitions, length, integer(1))

if (!identical(unname(program_sizes), c(436L, 284L, 46L, 84L))) {
    stop("Locked program sizes are not 436/284/46/84.")
}

primary_library_size <- Matrix::colSums(primary_counts)

atlas_library_size <- Matrix::colSums(atlas_counts)

if (any(primary_library_size <= 0) || any(atlas_library_size <= 0)) {
    stop("One or more pseudobulk libraries are empty.")
}

primary_VDR_count <- as.numeric(primary_counts[VDR_feature, , drop = TRUE])

atlas_VDR_count <- as.numeric(atlas_counts[VDR_feature, , drop = TRUE])

primary_VDR_expression <- cbind(primary_metadata, data.frame(VDR_feature = VDR_feature, VDR_raw_count = primary_VDR_count, 
    VDR_CPM = primary_VDR_count/primary_library_size * 1000000, VDR_log2_CPM_plus_1 = log2(primary_VDR_count/primary_library_size * 
        1000000 + 1), stringsAsFactors = FALSE))

atlas_VDR_expression <- cbind(atlas_metadata, data.frame(VDR_feature = VDR_feature, VDR_raw_count = atlas_VDR_count, 
    VDR_CPM = atlas_VDR_count/atlas_library_size * 1000000, VDR_log2_CPM_plus_1 = log2(atlas_VDR_count/atlas_library_size * 
        1000000 + 1), stringsAsFactors = FALSE))

available_candidate_table <- candidate_table[candidate_table$feature_available & is.finite(candidate_table$interaction_logFC) & 
    candidate_table$interaction_logFC != 0, , drop = FALSE]

available_features <- unique(available_candidate_table$resolved_feature)

primary_candidate_counts <- as.matrix(primary_counts[available_features, , drop = FALSE])

atlas_candidate_counts <- as.matrix(atlas_counts[available_features, , drop = FALSE])

primary_candidate_logCPM <- log2(sweep(primary_candidate_counts, MARGIN = 2, STATS = primary_library_size/1000000, 
    FUN = "/") + 1)

atlas_candidate_logCPM <- log2(sweep(atlas_candidate_counts, MARGIN = 2, STATS = atlas_library_size/1000000, FUN = "/") + 
    1)

reference_statistics <- row_zscore_reference(primary_candidate_logCPM)

names(reference_statistics$mean) <- rownames(primary_candidate_logCPM)

names(reference_statistics$sd) <- rownames(primary_candidate_logCPM)

primary_candidate_z <- apply_row_reference_z(primary_candidate_logCPM, reference_mean = reference_statistics$mean, 
    reference_sd = reference_statistics$sd)

atlas_candidate_z <- apply_row_reference_z(atlas_candidate_logCPM, reference_mean = reference_statistics$mean, 
    reference_sd = reference_statistics$sd)

score_programs <- function(standardized_matrix, pb_metadata, matrix_role) {
    score_rows <- list()
    coverage_rows <- list()
    for (program_name in names(program_definitions)) {
        frozen_genes <- unique(program_definitions[[program_name]])
        program_table <- candidate_table[candidate_table$feature_id %in% frozen_genes & candidate_table$feature_available, 
            , drop = FALSE]
        measurable_features <- unique(program_table$resolved_feature[program_table$resolved_feature %in% rownames(standardized_matrix)])
        feature_rows <- match(measurable_features, program_table$resolved_feature)
        program_table <- program_table[feature_rows, , drop = FALSE]
        valid <- (is.finite(program_table$interaction_logFC) & program_table$interaction_logFC != 0 & measurable_features %in% 
            rownames(standardized_matrix))
        measurable_features <- measurable_features[valid]
        program_table <- program_table[valid, , drop = FALSE]
        weights <- sign(program_table$interaction_logFC)
        variable_gene <- is.finite(reference_statistics$sd[measurable_features]) & reference_statistics$sd[measurable_features] > 
            0
        measurable_features <- measurable_features[variable_gene]
        weights <- weights[variable_gene]
        if (length(measurable_features) < 10L) {
            stop("Program ", program_name, " has fewer than 10 measurable variable genes.")
        }
        signed_matrix <- sweep(standardized_matrix[measurable_features, , drop = FALSE], MARGIN = 1, STATS = weights, 
            FUN = "*")
        score_value <- colMeans(signed_matrix, na.rm = TRUE)
        score_metadata <- as.data.frame(pb_metadata, stringsAsFactors = FALSE, check.names = FALSE)
        score_metadata$matrix_role <- matrix_role
        score_statistics <- data.frame(program_name = program_name, score = as.numeric(score_value), frozen_gene_n = length(frozen_genes), 
            measurable_gene_n = length(measurable_features), coverage_fraction = length(measurable_features)/length(frozen_genes), 
            score_definition = "MEAN_SIGNED_GENE_ZSCORE_USING_PRIMARY_REFERENCE", interpretation = "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION", 
            stringsAsFactors = FALSE, check.names = FALSE)
        overlapping_columns <- intersect(colnames(score_metadata), colnames(score_statistics))
        if (length(overlapping_columns)) {
            stop("Program-score metadata/statistic columns overlap: ", paste(overlapping_columns, collapse = ", "))
        }
        score_row <- cbind(score_metadata, score_statistics)
        if (anyDuplicated(colnames(score_row)) > 0L) {
            stop("Program-score table contains duplicate column names for ", program_name, ".")
        }
        score_rows[[program_name]] <- score_row
        coverage_rows[[program_name]] <- data.frame(program_name = program_name, frozen_gene_n = length(frozen_genes), 
            available_gene_n = sum(candidate_table$feature_id %in% frozen_genes & candidate_table$feature_available), 
            measurable_variable_gene_n = length(measurable_features), coverage_fraction = length(measurable_features)/length(frozen_genes), 
            missing_genes = collapse_unique(setdiff(frozen_genes, candidate_table$feature_id[candidate_table$feature_available])), 
            stringsAsFactors = FALSE)
    }
    list(scores = as.data.frame(data.table::rbindlist(score_rows, use.names = TRUE, fill = TRUE)), coverage = do.call(rbind, 
        coverage_rows))
}

primary_score_result <- score_programs(standardized_matrix = primary_candidate_z, pb_metadata = primary_metadata, 
    matrix_role = "PRIMARY")

atlas_score_result <- score_programs(standardized_matrix = atlas_candidate_z, pb_metadata = atlas_metadata, matrix_role = "ATLAS_SENSITIVITY")

program_coverage <- primary_score_result$coverage

if (any(program_coverage$coverage_fraction < 0.90000000000000002)) {
    stop("One or more locked VDR programs have below 90% coverage.")
}

primary_program_scores <- primary_score_result$scores

atlas_program_scores <- atlas_score_result$scores

score_column_audit <- data.frame(matrix_role = c("PRIMARY", "ATLAS_SENSITIVITY"), row_n = c(nrow(primary_program_scores), 
    nrow(atlas_program_scores)), column_n = c(ncol(primary_program_scores), ncol(atlas_program_scores)), duplicate_column_n = c(sum(duplicated(colnames(primary_program_scores))), 
    sum(duplicated(colnames(atlas_program_scores)))), column_names = c(paste(colnames(primary_program_scores), 
    collapse = ";"), paste(colnames(atlas_program_scores), collapse = ";")), stringsAsFactors = FALSE)

write_csv_utf8(score_column_audit, file.path(table_root, "00_program_score_column_audit.csv"))

if (any(score_column_audit$duplicate_column_n > 0L) || !all(c("matrix_role", "program_name", "score", "broad_cell_type", 
    "sample_accession", "minimum_cell_threshold_pass") %in% colnames(primary_program_scores))) {
    stop("Program-score column-name audit failed.")
}

write_csv_utf8(candidate_table, file.path(table_root, "08_candidate_feature_audit.csv"))

write_csv_utf8(program_coverage, file.path(table_root, "09_program_coverage.csv"))

write_csv_utf8(primary_VDR_expression, file.path(table_root, "10_primary_VDR_expression.csv"))

write_csv_utf8(atlas_VDR_expression, file.path(table_root, "11_atlas_VDR_expression.csv"))

write_csv_utf8(primary_program_scores, file.path(table_root, "12_primary_program_scores.csv"))

write_csv_utf8(atlas_program_scores, file.path(table_root, "13_atlas_program_scores.csv"))

eligibility_rows <- list()

for (cell_type in sort(unique(primary_metadata$broad_cell_type))) {
    cell_type_metadata <- primary_metadata[primary_metadata$broad_cell_type == cell_type, , drop = FALSE]
    Control_pass_n <- sum(cell_type_metadata$analysis_group == "CONTROL" & cell_type_metadata$minimum_cell_threshold_pass)
    RIF_pass_n <- sum(cell_type_metadata$analysis_group == "RIF" & cell_type_metadata$minimum_cell_threshold_pass)
    eligibility_rows[[length(eligibility_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, total_profile_n = nrow(cell_type_metadata), 
        Control_total_n = sum(cell_type_metadata$analysis_group == "CONTROL"), RIF_total_n = sum(cell_type_metadata$analysis_group == 
            "RIF"), Control_threshold_pass_n = Control_pass_n, RIF_threshold_pass_n = RIF_pass_n, minimum_cell_n = min(cell_type_metadata$cell_n), 
        median_cell_n = median(cell_type_metadata$cell_n), maximum_cell_n = max(cell_type_metadata$cell_n), donor_level_model_eligible = Control_pass_n >= 
            3L && RIF_pass_n >= 3L, descriptive_localization_eligible = sum(cell_type_metadata$minimum_cell_threshold_pass) >= 
            2L, stringsAsFactors = FALSE)
}

model_eligibility <- do.call(rbind, eligibility_rows)

model_eligible_types <- model_eligibility$broad_cell_type[model_eligibility$donor_level_model_eligible]

if (length(model_eligible_types) != 5L) {
    stop("Expected five donor-level model-eligible broad cell types, obtained ", length(model_eligible_types), 
        ": ", collapse_unique(model_eligible_types))
}

primary_unciliated_scores <- primary_program_scores[primary_program_scores$broad_cell_type == "Unciliated epithelial cells", 
    , drop = FALSE]

atlas_unciliated_scores <- atlas_program_scores[atlas_program_scores$broad_cell_type == "Unciliated epithelial cells", 
    , drop = FALSE]

cluster8_program_sensitivity <- merge(primary_unciliated_scores[, c("sample_accession", "donor_id", "analysis_group", 
    "program_name", "cell_n", "score"), drop = FALSE], atlas_unciliated_scores[, c("sample_accession", "donor_id", 
    "analysis_group", "program_name", "cell_n", "score"), drop = FALSE], by = c("sample_accession", "donor_id", 
    "analysis_group", "program_name"), all = TRUE, suffixes = c("_primary", "_atlas"), sort = FALSE)

cluster8_program_sensitivity$cluster8_cell_n <- (cluster8_program_sensitivity$cell_n_atlas - cluster8_program_sensitivity$cell_n_primary)

cluster8_program_sensitivity$score_delta_atlas_minus_primary <- (cluster8_program_sensitivity$score_atlas - cluster8_program_sensitivity$score_primary)

cluster8_VDR_sensitivity <- merge(primary_VDR_expression[primary_VDR_expression$broad_cell_type == "Unciliated epithelial cells", 
    c("sample_accession", "donor_id", "analysis_group", "cell_n", "VDR_log2_CPM_plus_1"), drop = FALSE], atlas_VDR_expression[atlas_VDR_expression$broad_cell_type == 
    "Unciliated epithelial cells", c("sample_accession", "donor_id", "analysis_group", "cell_n", "VDR_log2_CPM_plus_1"), 
    drop = FALSE], by = c("sample_accession", "donor_id", "analysis_group"), all = TRUE, suffixes = c("_primary", 
    "_atlas"), sort = FALSE)

cluster8_VDR_sensitivity$cluster8_cell_n <- (cluster8_VDR_sensitivity$cell_n_atlas - cluster8_VDR_sensitivity$cell_n_primary)

cluster8_VDR_sensitivity$VDR_delta_atlas_minus_primary <- (cluster8_VDR_sensitivity$VDR_log2_CPM_plus_1_atlas - 
    cluster8_VDR_sensitivity$VDR_log2_CPM_plus_1_primary)

localization_summary <- data.table::as.data.table(primary_program_scores[primary_program_scores$minimum_cell_threshold_pass, 
    , drop = FALSE])[, .(donor_n = data.table::uniqueN(donor_id), median_score = median(score), minimum_score = min(score), 
    maximum_score = max(score)), by = .(broad_cell_type, program_name)]

localization_summary <- as.data.frame(localization_summary)

VDR_localization_summary <- data.table::as.data.table(primary_VDR_expression[primary_VDR_expression$minimum_cell_threshold_pass, 
    , drop = FALSE])[, .(donor_n = data.table::uniqueN(donor_id), median_VDR_log2_CPM_plus_1 = median(VDR_log2_CPM_plus_1), 
    minimum_VDR_log2_CPM_plus_1 = min(VDR_log2_CPM_plus_1), maximum_VDR_log2_CPM_plus_1 = max(VDR_log2_CPM_plus_1), 
    median_VDR_positive_fraction = median(VDR_cell_audit$primary_VDR_positive_fraction[match(paste(donor_id, broad_cell_type), 
        paste(VDR_cell_audit$donor_id, VDR_cell_audit$broad_cell_type))], na.rm = TRUE)), by = broad_cell_type]

VDR_localization_summary <- as.data.frame(VDR_localization_summary)

write_csv_utf8(model_eligibility, file.path(table_root, "14_model_eligibility.csv"))

write_csv_utf8(cluster8_program_sensitivity, file.path(table_root, "15_cluster8_program_sensitivity.csv"))

write_csv_utf8(cluster8_VDR_sensitivity, file.path(table_root, "16_cluster8_VDR_sensitivity.csv"))

write_csv_utf8(localization_summary, file.path(table_root, "17_program_localization_summary.csv"))

write_csv_utf8(VDR_localization_summary, file.path(table_root, "18_VDR_localization_summary.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

cell_count_plot_data <- primary_metadata

cell_count_plot_data$broad_cell_type <- factor(cell_count_plot_data$broad_cell_type, levels = sort(unique(cell_count_plot_data$broad_cell_type)))

plot_cell_counts <- ggplot2::ggplot(cell_count_plot_data, ggplot2::aes(x = broad_cell_type, y = cell_n, shape = analysis_group)) + 
    ggplot2::geom_hline(yintercept = minimum_cells_primary, linetype = "dashed") + ggplot2::geom_point(size = 2.5, 
    alpha = 0.80000000000000004, position = ggplot2::position_jitter(width = 0.12, height = 0)) + ggplot2::scale_y_log10(labels = scales::comma) + 
    ggplot2::labs(title = "GSE183837 primary pseudobulk cell counts", subtitle = "Cluster 8 excluded; dashed line denotes the 20-cell model threshold", 
        x = "Broad cell type", y = "Cells per donor pseudobulk (log scale)", shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_cell_counts, "QC1_primary_pseudobulk_cell_counts", width = 10, height = 6)

library_plot_data <- primary_metadata

plot_library <- ggplot2::ggplot(library_plot_data, ggplot2::aes(x = broad_cell_type, y = library_size, shape = analysis_group)) + 
    ggplot2::geom_point(size = 2.5, alpha = 0.80000000000000004, position = ggplot2::position_jitter(width = 0.12, 
        height = 0)) + ggplot2::scale_y_log10(labels = scales::comma) + ggplot2::labs(title = "Primary pseudobulk library sizes", 
    subtitle = "Raw counts aggregated by donor and locked broad cell type", x = "Broad cell type", y = "Library size (log scale)", 
    shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

save_plot_pair(plot_library, "QC2_primary_library_sizes", width = 10, height = 6)

VDR_plot_data <- primary_VDR_expression[primary_VDR_expression$minimum_cell_threshold_pass, , drop = FALSE]

plot_VDR <- ggplot2::ggplot(VDR_plot_data, ggplot2::aes(x = broad_cell_type, y = VDR_log2_CPM_plus_1, shape = analysis_group)) + 
    ggplot2::geom_point(size = 2.5, alpha = 0.80000000000000004, position = ggplot2::position_jitter(width = 0.12, 
        height = 0)) + ggplot2::labs(title = "Descriptive VDR expression localization", subtitle = "Threshold-passing primary pseudobulks; no group test in this step", 
    x = "Broad cell type", y = "VDR log2(CPM + 1)", shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_VDR, "QC3_VDR_localization", width = 10, height = 6)

score_heatmap_data <- primary_program_scores[primary_program_scores$minimum_cell_threshold_pass, , drop = FALSE]

if (anyDuplicated(colnames(score_heatmap_data)) > 0L) {
    stop("score_heatmap_data contains duplicate column names before plotting.")
}

score_heatmap_data$profile_label <- paste(score_heatmap_data$sample_accession, score_heatmap_data$broad_cell_type, 
    sep = " | ")

score_heatmap_data$profile_label <- factor(score_heatmap_data$profile_label, levels = unique(score_heatmap_data$profile_label))

score_heatmap_data$program_name <- factor(score_heatmap_data$program_name, levels = c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", 
    "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84"))

plot_score_heatmap <- ggplot2::ggplot(score_heatmap_data, ggplot2::aes(x = program_name, y = profile_label, fill = score)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2() + ggplot2::labs(title = "Locked VDR-program localization across donor pseudobulks", 
    subtitle = "Scores use the primary matrix and signed interaction direction", x = "VDR program", y = "Donor | broad cell type", 
    fill = "Signed score") + theme_QC + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5), axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_score_heatmap, "QC4_program_score_heatmap", width = 10, height = 14)

localization_plot_data <- localization_summary

localization_plot_data$broad_cell_type <- factor(localization_plot_data$broad_cell_type, levels = sort(unique(localization_plot_data$broad_cell_type)))

plot_localization <- ggplot2::ggplot(localization_plot_data, ggplot2::aes(x = program_name, y = broad_cell_type, 
    fill = median_score)) + ggplot2::geom_tile() + ggplot2::scale_fill_gradient2() + ggplot2::labs(title = "Median VDR-program alignment by broad cell type", 
    subtitle = "Descriptive localization among threshold-passing donor profiles", x = "VDR program", y = "Broad cell type", 
    fill = "Median score") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

save_plot_pair(plot_localization, "QC5_program_localization_summary", width = 9, height = 6)

eligibility_plot_data <- model_eligibility

eligibility_long <- data.frame(broad_cell_type = rep(eligibility_plot_data$broad_cell_type, each = 2), analysis_group = rep(c("CONTROL", 
    "RIF"), times = nrow(eligibility_plot_data)), threshold_pass_donor_n = as.vector(t(as.matrix(eligibility_plot_data[, 
    c("Control_threshold_pass_n", "RIF_threshold_pass_n")]))), donor_level_model_eligible = rep(eligibility_plot_data$donor_level_model_eligible, 
    each = 2), stringsAsFactors = FALSE)

plot_eligibility <- ggplot2::ggplot(eligibility_long, ggplot2::aes(x = analysis_group, y = broad_cell_type, fill = threshold_pass_donor_n)) + 
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = threshold_pass_donor_n)) + ggplot2::labs(title = "Donor-level model eligibility", 
    subtitle = "Primary model requires at least 3 threshold-passing donors in each group", x = "Group", y = "Broad cell type", 
    fill = "Donors") + theme_QC

save_plot_pair(plot_eligibility, "QC6_model_eligibility", width = 7, height = 5.5)

cluster8_plot_data <- cluster8_program_sensitivity

plot_cluster8 <- ggplot2::ggplot(cluster8_plot_data, ggplot2::aes(x = program_name, y = score_delta_atlas_minus_primary, 
    shape = analysis_group)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + ggplot2::geom_point(size = 2.5, 
    alpha = 0.80000000000000004, position = ggplot2::position_jitter(width = 0.12, height = 0)) + ggplot2::labs(title = "Cluster-8 sensitivity effect on unciliated epithelial scores", 
    subtitle = "Atlas score minus primary score; sensitivity only", x = "VDR program", y = "Score difference", 
    shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

save_plot_pair(plot_cluster8, "QC7_cluster8_program_sensitivity", width = 9, height = 5.5)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

gate_checks <- data.frame(check = c("upstream_annotation_gate_authorizes_pseudobulk", "annotated_checkpoint_has_47161_cells", 
    "nine_count_layers_processed", "all_layer_jobs_complete", "primary_matrix_has_45630_cells", "atlas_matrix_has_47161_cells", 
    "cluster8_exclusion_is_1531_cells", "primary_matrix_has_54_profiles", "atlas_matrix_has_54_profiles", "raw_count_conservation_passes_all_layers", 
    "primary_count_conservation_passes_all_layers", "primary_metadata_aligned_to_counts", "atlas_metadata_aligned_to_counts", 
    "donor_design_is_3_Control_6_RIF", "locked_program_sizes_are_436_284_46_84", "all_program_coverage_at_least_90_percent", 
    "VDR_is_measurable", "five_broad_types_are_donor_model_eligible", "cluster8_is_sensitivity_only", "primary_and_atlas_count_checkpoints_exist", 
    "donor_level_VDR_program_models_can_be_authorized", "disease_testing_not_performed_in_this_step", "free_disk_after_run_at_least_8GB_or_unknown"), 
    pass = c(upstream_gate$gate_status[[1]] == "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION", 
        expected_atlas_cell_n == sum(atlas_metadata$cell_n), length(count_layers) == 9L, nrow(layer_status_all) == 
            9L && all(layer_status_all$status == "COMPLETE"), sum(primary_metadata$cell_n) == expected_primary_cell_n, 
        sum(atlas_metadata$cell_n) == expected_atlas_cell_n, sum(cell_audit$cluster8_excluded_cell_n) == expected_cluster8_cell_n, 
        ncol(primary_counts) == 54L, ncol(atlas_counts) == 54L, all(raw_count_audit$raw_count_sum == raw_count_audit$atlas_aggregated_count_sum), 
        all(raw_count_audit$primary_input_count_sum == raw_count_audit$primary_aggregated_count_sum), identical(colnames(primary_counts), 
            primary_metadata$pb_id), identical(colnames(atlas_counts), atlas_metadata$pb_id), nrow(donor_group_identity) == 
            9L && sum(donor_group_identity$analysis_group == "CONTROL") == 3L && sum(donor_group_identity$analysis_group == 
            "RIF") == 6L, identical(unname(program_sizes), c(436L, 284L, 46L, 84L)), all(program_coverage$coverage_fraction >= 
            0.90000000000000002), VDR_feature %in% rownames(primary_counts), length(model_eligible_types) == 5L, 
        all(primary_metadata$cell_n[primary_metadata$broad_cell_type == "Unciliated epithelial cells"] <= atlas_metadata$cell_n[match(paste(primary_metadata$sample_accession[primary_metadata$broad_cell_type == 
            "Unciliated epithelial cells"], "Unciliated epithelial cells"), paste(atlas_metadata$sample_accession, 
            atlas_metadata$broad_cell_type))]), file.exists(primary_counts_path) && file.exists(atlas_counts_path), 
        TRUE, TRUE, is.na(get_free_disk_gb()) || get_free_disk_gb() >= 8), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED"
} else {
    "HOLD_GSE183837_PSEUDOBULK_OR_LOCALIZATION_FAILED"
}

review_flags <- c("CLUSTER8_SENSITIVITY_ONLY", "CILIATED_EPITHELIAL_DESCRIPTIVE_LOW_CELL_COUNTS", "CONTROL_N3_RIF_N6_LOW_POWER_DONOR_MODELS", 
    "NO_CELL_LEVEL_DISEASE_INFERENCE")

pb_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, annotation_version = annotation_version, 
    pseudobulk_version = pseudobulk_version, gate_status = gate_status, feature_n = nrow(primary_counts), primary_profile_n = ncol(primary_counts), 
    atlas_profile_n = ncol(atlas_counts), primary_cell_n = sum(primary_metadata$cell_n), atlas_cell_n = sum(atlas_metadata$cell_n), 
    cluster8_excluded_cell_n = sum(cell_audit$cluster8_excluded_cell_n), minimum_cells_primary = minimum_cells_primary, 
    donor_model_eligible_type_n = length(model_eligible_types), donor_model_eligible_types = collapse_unique(model_eligible_types), 
    descriptive_only_types = collapse_unique(model_eligibility$broad_cell_type[!model_eligibility$donor_level_model_eligible]), 
    review_flags = collapse_unique(review_flags), donor_level_VDR_program_models_authorized = all_gate_checks_pass, 
    gene_level_pseudobulk_DE_authorized = FALSE, cell_level_disease_testing_authorized = FALSE, result_lock_authorized = FALSE, 
    primary_counts_checkpoint = normalize_slashes(primary_counts_path), primary_counts_md5 = unname(tools::md5sum(primary_counts_path)), 
    stringsAsFactors = FALSE)

write_csv_utf8(pb_gate, file.path(run_qc_root, "G183_PBVDR_GATE.csv"))

card_lines <- c("# STEP04E1N GSE183837 Pseudobulk and VDR Localization Card", "", paste0("- Gate: ", gate_status), 
    paste0("- Primary cells: ", format(sum(primary_metadata$cell_n), big.mark = ",")), paste0("- Atlas cells: ", 
        format(sum(atlas_metadata$cell_n), big.mark = ",")), paste0("- Cluster-8 cells excluded from primary: ", 
        format(sum(cell_audit$cluster8_excluded_cell_n), big.mark = ",")), paste0("- Primary pseudobulk profiles: ", 
        ncol(primary_counts)), paste0("- Atlas sensitivity profiles: ", ncol(atlas_counts)), paste0("- Model-eligible broad cell types: ", 
        collapse_unique(model_eligible_types)), paste0("- Descriptive-only broad cell types: ", collapse_unique(model_eligibility$broad_cell_type[!model_eligibility$donor_level_model_eligible])), 
    paste0("- Donor-level VDR/program models authorized: ", all_gate_checks_pass), "- Gene-level pseudobulk DE authorized: FALSE", 
    "- Cell-level disease testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "G183_PBVDR_CARD.md"))

interpretation_lines <- c("# GSE183837 pseudobulk and VDR-localization interpretation", "", paste0("- Run ID: `", 
    run_id, "`"), paste0("- Gate: **", gate_status, "**"), "", "## Primary matrix", "", "- Uses donor x locked broad cell type raw-count aggregation.", 
    "- Excludes all 1,531 low-confidence cluster-8 cells.", "- Retains 45,630 cells across 54 donor x broad-type profiles.", 
    "", "## Atlas sensitivity matrix", "", "- Retains all 47,161 cells.", "- Used only to quantify how cluster 8 changes unciliated epithelial VDR and program scores.", 
    "", "## Model authorization", "", paste0("- Donor-level model-eligible broad cell types: ", collapse_unique(model_eligible_types), 
        "."), "- Ciliated epithelial cells remain descriptive because too few donor profiles pass the 20-cell threshold.", 
    "- Control n=3 and RIF n=6 imply low power; exact permutation and effect-size reporting are required in the next step.", 
    "", "## Interpretation boundary", "", "- Program scores are relative alignment scores, not direct VDR target activity.", 
    "- No group difference has been tested here.", "- No cell-level disease inference is allowed.", "- No result is locked.")

write_text(interpretation_lines, file.path(run_qc_root, "G183_PBVDR_INTERPRETATION.md"))

parameter_manifest <- list(script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_DUPLICATE_MATRIX_ROLE_COLUMN_BLOCKED_GGPLOT", 
    score_table_column_policy = "OVERWRITE_EXISTING_MATRIX_ROLE_AND_REQUIRE_UNIQUE_COLUMN_NAMES", dataset_id = dataset_id, 
    integration_version = integration_version, annotation_version = annotation_version, pseudobulk_version = pseudobulk_version, 
    upstream_annotation_gate = upstream_gate$gate_status[[1]], minimum_cells_primary = minimum_cells_primary, expected_primary_cell_n = expected_primary_cell_n, 
    expected_atlas_cell_n = expected_atlas_cell_n, cluster8_policy = "PRIMARY_EXCLUDED_ATLAS_SENSITIVITY_RETAINED", 
    program_score_definition = "MEAN_SIGNED_GENE_ZSCORE_USING_PRIMARY_REFERENCE", program_interpretation = "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION", 
    donor_model_eligibility = "AT_LEAST_3_CONTROL_AND_3_RIF_DONORS_WITH_AT_LEAST_20_CELLS", donor_level_VDR_program_models_authorized = all_gate_checks_pass, 
    gene_level_pseudobulk_DE_authorized = FALSE, cell_level_disease_testing_authorized = FALSE, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_DUPLICATE_MATRIX_ROLE_COLUMN_BLOCKED_GGPLOT", 
    dataset_id = dataset_id, integration_version = integration_version, annotation_version = annotation_version, 
    pseudobulk_version = pseudobulk_version, started_at = as.character(run_start), completed_at = as.character(run_end), 
    elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 2), gate_status = gate_status, 
    primary_cell_n = sum(primary_metadata$cell_n), atlas_cell_n = sum(atlas_metadata$cell_n), primary_profile_n = ncol(primary_counts), 
    atlas_profile_n = ncol(atlas_counts), donor_model_eligible_type_n = length(model_eligible_types), donor_model_eligible_types = collapse_unique(model_eligible_types), 
    donor_level_models_authorized = all_gate_checks_pass, result_lock_authorized = FALSE, free_memory_before_GB = free_memory_before_GB, 
    free_disk_before_GB = free_disk_before_GB, free_disk_after_GB = get_free_disk_gb(), output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

persistent_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)

persistent_files <- persistent_files[file.exists(persistent_files) & !dir.exists(persistent_files)]

persistent_manifest <- data.frame(path = normalize_slashes(persistent_files), size_bytes = file.info(persistent_files)$size, 
    md5 = unname(tools::md5sum(persistent_files)), stringsAsFactors = FALSE)

write_csv_utf8(persistent_manifest, file.path(manifest_root, "persistent_manifest.csv"))

methods_lines <- c("# STEP04E1N Methods shell", "", "Raw counts were aggregated within each donor and locked broad cell type from the nine split Seurat v5 count layers. The primary matrix excluded the low-confidence cluster 8, whereas an atlas sensitivity matrix retained all cells. Profiles with at least 20 cells were considered threshold-passing. VDR expression was summarized as log2(CPM + 1). Four prospectively locked VDR interaction programs were scored using mean signed gene-wise z-scores, with signs determined by the calcitriol-dependent interaction log-fold changes. Program scores were interpreted as relative alignment with the locked interaction direction rather than evidence of direct VDR binding or causal activation.")

results_lines <- c("# STEP04E1N Results shell", "", paste0("The primary donor-by-cell-type pseudobulk matrix retained ", 
    format(sum(primary_metadata$cell_n), big.mark = ","), " cells in ", ncol(primary_counts), " profiles after excluding ", 
    format(sum(cell_audit$cluster8_excluded_cell_n), big.mark = ","), " low-confidence cluster-8 cells."), paste0("Five broad cell types met the prespecified donor-level model criterion: ", 
    collapse_unique(model_eligible_types), "."), "Ciliated epithelial cells were retained for descriptive localization but did not satisfy the donor-level modelling threshold.")

write_text(methods_lines, file.path(manuscript_text_root, paste0(run_id, "_Methods.md")))

write_text(results_lines, file.path(manuscript_text_root, paste0(run_id, "_Results.md")))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(figure_output_root, 
    recursive = TRUE, full.names = TRUE), list.files(manifest_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), pb_status_path, log_path, run_manifest_path, parameter_path, session_path, 
    normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2]

artifact_manifest <- data.frame(path = normalize_slashes(small_files), size_bytes = file.info(small_files)$size, 
    md5 = unname(tools::md5sum(small_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "primary_profiles=", ncol(primary_counts), "eligible_types=", 
    length(model_eligible_types))

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

cat("STEP04E1N GSE183837 pseudobulk/VDR localization completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Primary cells: ", format(sum(primary_metadata$cell_n), big.mark = ","), "\n", sep = "")

cat("Atlas cells: ", format(sum(atlas_metadata$cell_n), big.mark = ","), "\n", sep = "")

cat("Cluster-8 cells excluded: ", format(sum(cell_audit$cluster8_excluded_cell_n), big.mark = ","), "\n", sep = "")

cat("Primary pseudobulk profiles: ", ncol(primary_counts), "\n", sep = "")

cat("Model-eligible broad cell types: ", length(model_eligible_types), "\n", sep = "")

cat("Eligible types: ", collapse_unique(model_eligible_types), "\n", sep = "")

cat("Donor-level VDR/program models authorized: ", all_gate_checks_pass, "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nPseudobulk outputs are complete in the project folders.\n", sep = "")
}

