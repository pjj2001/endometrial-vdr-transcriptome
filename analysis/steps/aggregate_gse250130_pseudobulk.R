#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C06: Aggregate raw counts and calculate program scores

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260805)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "aggregate_gse250130_pseudobulk.R"

script_version <- "STEP04E1F_GSE250130_VDR_LOCALIZATION_PSEUDOBULK_v1.0"

input_version <- "SCI_v1_1_memory_safe"

annotation_version <- "ANN_v1_2"

pseudobulk_version <- "PBVDR_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1F_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(path, max_chars = 245L) {
    normalized <- normalize_slashes(path)
    path_n <- nchar(normalized, type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Output path is too long for safe Windows/R writing (", path_n, " characters; safe limit ", max_chars, 
            "):\n", normalized)
    }
    invisible(path)
}

write_csv_utf8 <- function(x, path) {
    assert_windows_path_length(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
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
    if (!requireNamespace(pkg, quietly = TRUE)) 
        return(NA_character_)
    as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) 
        ""
    else paste(x, collapse = sep)
}

format_bytes <- function(bytes) {
    if (length(bytes) == 0L || is.na(bytes)) 
        return("")
    if (bytes < 1024) 
        return(sprintf("%.0f B", bytes))
    if (bytes < 1024^2) 
        return(sprintf("%.1f KB", bytes/1024))
    if (bytes < 1024^3) 
        return(sprintf("%.1f MB", bytes/1024^2))
    sprintf("%.2f GB", bytes/1024^3)
}

sanitize_file_token <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x <- gsub("_+", "_", x)
    sub("^_+|_+$", "", x)
}

get_free_memory_gb <- function() {
    if (.Platform$OS.type != "windows") 
        return(NA_real_)
    lines <- tryCatch(system(paste0("powershell -NoProfile -Command ", "\"$o=Get-CimInstance Win32_OperatingSystem; ", 
        "[math]::Round($o.FreePhysicalMemory/1MB,2)\""), intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
    if (!length(lines)) 
        return(NA_real_)
    suppressWarnings(as.numeric(trimws(lines[1])))
}

get_free_disk_gb <- function(drive_letter = "C") {
    if (.Platform$OS.type != "windows") 
        return(NA_real_)
    lines <- tryCatch(system(paste0("powershell -NoProfile -Command ", "\"[math]::Round((Get-PSDrive -Name ", drive_letter, 
        ").Free/1GB,2)\""), intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
    if (!length(lines)) 
        return(NA_real_)
    suppressWarnings(as.numeric(trimws(lines[1])))
}

row_zscore_dense <- function(x) {
    x <- as.matrix(x)
    row_mean <- rowMeans(x, na.rm = TRUE)
    row_sd <- apply(x, 1L, stats::sd, na.rm = TRUE)
    output <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    valid <- is.finite(row_sd) & row_sd > 0
    output[valid, ] <- (x[valid, , drop = FALSE] - row_mean[valid])/row_sd[valid]
    output
}

save_plot_pair <- function(plot_object, stem, width, height) {
    png_path <- file.path(figure_root, paste0(stem, ".png"))
    pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))
    assert_windows_path_length(png_path)
    assert_windows_path_length(pdf_path)
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

required_packages <- c("SeuratObject", "Matrix", "data.table", "ggplot2", "jsonlite", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
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

checkpoint_path <- file.path(results_root, "S4E1_sc", input_version, "checkpoints", dataset_id, paste0(dataset_id, 
    "_integrated_checkpoint.rds"))

checkpoint_status_path <- file.path(results_root, "S4E1_sc", input_version, "checkpoints", dataset_id, paste0(dataset_id, 
    "_checkpoint_status.csv"))

annotation_root <- file.path(results_root, "S4E1_sc", input_version, "annotations", dataset_id, annotation_version)

annotation_path <- file.path(annotation_root, "tables", "cell_annotation_metadata.csv")

candidate_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

output_root <- file.path(results_root, "S4E1_sc", input_version, "pb", dataset_id, pseudobulk_version)

per_layer_root <- file.path(output_root, "per_layer")

object_root <- file.path(output_root, "objects")

table_root <- file.path(output_root, "tables")

source_data_root <- file.path(output_root, "source_data")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1F_", format(run_start, "%Y%m%d_%H%M%S"), "_pbQC"))

figure_root <- file.path(run_qc_root, "figures")

for (p in c(execution_root, parameter_root, environment_root, manuscript_text_root, output_root, per_layer_root, 
    object_root, table_root, source_data_root, run_qc_root, figure_root)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(execution_root, paste0(run_id, "_execution.log"))

run_manifest_path <- file.path(execution_root, paste0(run_id, "_RUN_MANIFEST.csv"))

parameter_path <- file.path(parameter_root, paste0(run_id, "_parameters.json"))

session_path <- file.path(environment_root, paste0("sessionInfo_", run_id, ".txt"))

desktop_candidates <- .vdr_report_root()

desktop_candidates <- .vdr_report_root()

zip_parent <- .vdr_report_root()

handoff_zip <- file.path(zip_parent, paste0(run_id, "_SMALL_REPORT.zip"))

planned_paths <- data.frame(object = c("output_root", "per_layer_root", "object_root", "table_root", "source_data_root", 
    "run_qc_root", "figure_root", "log_path", "run_manifest_path", "parameter_path", "session_path", "handoff_zip"), 
    path = normalize_slashes(c(output_root, per_layer_root, object_root, table_root, source_data_root, run_qc_root, 
        figure_root, log_path, run_manifest_path, parameter_path, session_path, handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

if (.Platform$OS.type == "windows" && any(planned_paths$character_n > 245L)) {
    bad <- planned_paths[planned_paths$character_n > 245L, , drop = FALSE]
    stop("One or more planned paths are too long:\n", paste(paste0(bad$character_n, " chars: ", bad$path), collapse = "\n"))
}

required_inputs <- c(checkpoint_path, checkpoint_status_path, annotation_path, candidate_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

annotation_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

annotation_qc_dirs <- annotation_qc_dirs[grepl("^STEP04E1E_[0-9]{8}_[0-9]{6}_annQC$", basename(annotation_qc_dirs), 
    perl = TRUE)]

if (!length(annotation_qc_dirs)) {
    stop("No completed STEP04E1E annotation-QC folder found.")
}

latest_annotation_qc <- annotation_qc_dirs[which.max(file.info(annotation_qc_dirs)$mtime)]

annotation_gate_path <- file.path(latest_annotation_qc, "ANNOTATION_GATE.csv")

annotation_gate <- safe_read_csv(annotation_gate_path)

checkpoint_status <- safe_read_csv(checkpoint_status_path)

if (is.null(annotation_gate) || nrow(annotation_gate) != 1L || annotation_gate$gate_status[[1]] != "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW" || 
    !isTRUE(annotation_gate$broad_VDR_localization_authorized[[1]]) || !isTRUE(annotation_gate$broad_pseudobulk_aggregation_authorized[[1]])) {
    stop("STEP04E1E Gate does not authorize VDR localization and pseudobulk aggregation.")
}

if (is.null(checkpoint_status) || nrow(checkpoint_status) != 1L || checkpoint_status$checkpoint_status[[1]] != 
    "COMPLETE" || checkpoint_status$join_policy[[1]] != "DATA_ONLY_COUNTS_REMAIN_SPLIT") {
    stop("Integrated checkpoint status is absent or invalid.")
}

annotation <- .vdr_fread(annotation_path, data.table = FALSE, showProgress = FALSE)

candidate_table <- safe_read_csv(candidate_path)

required_annotation_columns <- c("cell_barcode", "sample_accession", "donor_id", "broad_cell_type", "clinical_group", 
    "LH_day", "temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "LH7_reference_group", "LH7_fertile_subgroup")

missing_annotation_columns <- setdiff(required_annotation_columns, colnames(annotation))

if (length(missing_annotation_columns)) {
    stop("Annotation table is missing columns: ", paste(missing_annotation_columns, collapse = ", "))
}

required_candidate_columns <- c("feature_id", "interaction_logFC", "occupancy_supported", "shared_occupancy_supported", 
    "any_promoter_overlap_support")

missing_candidate_columns <- setdiff(required_candidate_columns, colnames(candidate_table))

if (length(missing_candidate_columns)) {
    stop("Locked VDR candidate table is missing columns: ", paste(missing_candidate_columns, collapse = ", "))
}

if (nrow(annotation) != 236636L || length(unique(annotation$sample_accession)) != 28L || length(unique(annotation$broad_cell_type)) != 
    8L || anyDuplicated(annotation$cell_barcode) > 0L) {
    stop("Annotation identity checks failed.")
}

if (nrow(candidate_table) != 436L || length(unique(candidate_table$feature_id)) != 436L || sum(candidate_table$occupancy_supported) != 
    284L || sum(candidate_table$shared_occupancy_supported) != 46L || sum(candidate_table$any_promoter_overlap_support) != 
    84L) {
    stop("Locked 436/284/46/84 VDR candidate table failed validation.")
}

free_memory_gb <- get_free_memory_gb()

free_disk_gb <- get_free_disk_gb("C")

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
    stop("Free physical memory is below 18 GB. Restart RStudio and rerun in a clean session.")
}

if (is.finite(free_disk_gb) && free_disk_gb < 12) {
    stop("Free C: drive space is below 12 GB.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, input_version = input_version, 
    annotation_version = annotation_version, pseudobulk_version = pseudobulk_version, dataset_id = dataset_id, 
    aggregation_level = "sample_donor_x_broad_cell_type", raw_count_layer_policy = "PROCESS_SPLIT_COUNTS_LAYERS_SEPARATELY", 
    minimum_cells_for_score = 20L, minimum_cells_preferred = 50L, localization_normalization = "LOG2_LIBRARY_SIZE_CPM_PLUS_1", 
    program_score = "MEAN_GENE_ZSCORE_MULTIPLIED_BY_SIGN_INTERACTION_LOGFC", score_role = "DESCRIPTIVE_LOCALIZATION_NOT_CAUSAL_OR_INFERENTIAL", 
    temporal_model_eligibility = "AT_LEAST_2_FERTILE_DONORS_AT_EACH_LH3_LH5_LH7_LH9_LH11_AND_TOTAL_GE_12", LH7_model_eligibility = "AT_LEAST_3_FERTILE_LH7_AND_3_RIF_LH7_PSEUDOBULKS", 
    composition_testing_authorized = FALSE, subtype_analysis_authorized = FALSE, random_seed = 20260805L, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

resource_snapshot <- data.frame(metric = c("free_physical_memory_GB_before_load", "free_C_drive_GB_before_run", 
    "checkpoint_size_GB", "checkpoint_md5", "annotation_cell_n", "annotation_sample_n", "annotation_broad_type_n", 
    "candidate_gene_n", "occupancy_supported_gene_n", "shared_occupancy_gene_n", "promoter_occupancy_gene_n", "annotation_gate"), 
    value = c(free_memory_gb, free_disk_gb, round(file.info(checkpoint_path)$size/1024^3, 3), unname(tools::md5sum(checkpoint_path)), 
        nrow(annotation), length(unique(annotation$sample_accession)), length(unique(annotation$broad_cell_type)), 
        nrow(candidate_table), sum(candidate_table$occupancy_supported), sum(candidate_table$shared_occupancy_supported), 
        sum(candidate_table$any_promoter_overlap_support), annotation_gate$gate_status[[1]]), stringsAsFactors = FALSE)

write_csv_utf8(resource_snapshot, file.path(run_qc_root, "01_input_resource_snapshot.csv"))

log_event("Starting", run_id, "annotation_gate=", annotation_gate$gate_status[[1]], "free_memory_GB=", free_memory_gb, 
    "free_disk_GB=", free_disk_gb)

load_start <- Sys.time()

object <- readRDS(checkpoint_path)

log_event("Loaded checkpoint", "elapsed_min=", round(as.numeric(difftime(Sys.time(), load_start, units = "mins")), 
    2))

if (!inherits(object, "Seurat")) {
    stop("Checkpoint is not a Seurat object.")
}

if (ncol(object) != 236636L || nrow(object) != 36601L) {
    stop("Unexpected checkpoint dimensions: ", nrow(object), " genes x ", ncol(object), " cells")
}

object_cells <- colnames(object)

annotation_cells <- as.character(annotation$cell_barcode)

if (!setequal(object_cells, annotation_cells)) {
    stop("Checkpoint cells and annotation cells are not identical.")
}

annotation <- annotation[match(object_cells, annotation$cell_barcode), , drop = FALSE]

if (!identical(as.character(annotation$cell_barcode), object_cells)) {
    stop("Failed to align annotation rows to object cells.")
}

annotation_lookup <- annotation

rownames(annotation_lookup) <- annotation_lookup$cell_barcode

SeuratObject::DefaultAssay(object) <- "RNA"

count_layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts")

if (length(count_layers) != 28L) {
    stop("Expected 28 split counts layers; found ", length(count_layers), ".")
}

write_csv_utf8(data.frame(layer_index = seq_along(count_layers), layer = count_layers, stringsAsFactors = FALSE), 
    file.path(run_qc_root, "02_counts_layer_manifest.csv"))

layer_status_rows <- list()

processed_cell_lists <- list()

raw_count_total_from_layers <- 0

for (layer_index in seq_along(count_layers)) {
    layer_name <- count_layers[layer_index]
    layer_token <- sanitize_file_token(layer_name)
    layer_rds <- file.path(per_layer_root, paste0(sprintf("%02d", layer_index), "_", layer_token, ".rds"))
    layer_status_path <- file.path(per_layer_root, paste0(sprintf("%02d", layer_index), "_", layer_token, "_status.csv"))
    existing_status <- safe_read_csv(layer_status_path)
    if (file.exists(layer_rds) && !is.null(existing_status) && nrow(existing_status) == 1L && existing_status$status[[1]] == 
        "COMPLETE") {
        log_event("SKIP completed counts layer", layer_index, "/", length(count_layers), layer_name)
        saved <- readRDS(layer_rds)
        layer_status_rows[[length(layer_status_rows) + 1L]] <- existing_status
        processed_cell_lists[[length(processed_cell_lists) + 1L]] <- saved$cell_barcodes
        raw_count_total_from_layers <- raw_count_total_from_layers + as.numeric(saved$raw_count_total)
        rm(saved)
        gc()
        next
    }
    log_event("START counts layer", layer_index, "/", length(count_layers), layer_name)
    layer_start <- Sys.time()
    layer_error <- ""
    result <- tryCatch({
        count_matrix <- SeuratObject::LayerData(object = object, assay = "RNA", layer = layer_name)
        if (!inherits(count_matrix, "sparseMatrix")) {
            count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
        }
        layer_cells <- colnames(count_matrix)
        if (is.null(layer_cells) || !length(layer_cells)) {
            stop("Counts layer has no cell names.")
        }
        layer_annotation <- annotation_lookup[layer_cells, , drop = FALSE]
        if (nrow(layer_annotation) != length(layer_cells) || anyNA(layer_annotation$sample_accession) || anyNA(layer_annotation$broad_cell_type)) {
            stop("Counts-layer cells could not be fully annotated.")
        }
        group_id <- paste(layer_annotation$sample_accession, layer_annotation$broad_cell_type, sep = "||")
        group_factor <- factor(group_id, levels = unique(group_id))
        membership <- Matrix::sparseMatrix(i = seq_along(layer_cells), j = as.integer(group_factor), x = 1, dims = c(length(layer_cells), 
            nlevels(group_factor)), dimnames = list(layer_cells, levels(group_factor)))
        pb_counts <- count_matrix %*% membership
        colnames(pb_counts) <- levels(group_factor)
        cell_count_table <- as.data.frame(table(group_id), stringsAsFactors = FALSE)
        colnames(cell_count_table) <- c("pb_id", "cell_n")
        cell_count_table$cell_n <- as.integer(cell_count_table$cell_n)
        group_rows <- match(cell_count_table$pb_id, group_id)
        pb_metadata <- layer_annotation[group_rows, c("sample_accession", "donor_id", "broad_cell_type", "clinical_group", 
            "LH_day", "temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "LH7_reference_group", "LH7_fertile_subgroup", 
            "analysis_role"), drop = FALSE]
        pb_metadata$pb_id <- cell_count_table$pb_id
        pb_metadata$cell_n <- cell_count_table$cell_n
        pb_metadata$source_counts_layer <- layer_name
        pb_metadata <- pb_metadata[match(colnames(pb_counts), pb_metadata$pb_id), , drop = FALSE]
        if (anyNA(pb_metadata$pb_id) || !identical(as.character(pb_metadata$pb_id), colnames(pb_counts))) {
            stop("Pseudobulk metadata did not align with count columns.")
        }
        raw_count_total <- as.numeric(sum(count_matrix))
        aggregated_count_total <- as.numeric(sum(pb_counts))
        if (!isTRUE(all.equal(raw_count_total, aggregated_count_total, tolerance = 0))) {
            stop("Raw-count total changed during aggregation: ", raw_count_total, " vs ", aggregated_count_total)
        }
        saved_object <- list(layer_index = layer_index, layer_name = layer_name, counts = pb_counts, metadata = pb_metadata, 
            cell_barcodes = layer_cells, raw_count_total = raw_count_total, gene_n = nrow(pb_counts), pseudobulk_n = ncol(pb_counts))
        saveRDS(saved_object, layer_rds, compress = FALSE)
        list(saved_object = saved_object, raw_count_total = raw_count_total)
    }, error = function(e) {
        layer_error <<- conditionMessage(e)
        NULL
    })
    elapsed_minutes <- round(as.numeric(difftime(Sys.time(), layer_start, units = "mins")), 2)
    if (is.null(result)) {
        status <- data.frame(layer_index = layer_index, layer_name = layer_name, status = "FAILED", cell_n = NA_integer_, 
            sample_n = NA_integer_, broad_type_n = NA_integer_, pseudobulk_n = NA_integer_, raw_count_total = NA_real_, 
            elapsed_minutes = elapsed_minutes, error_message = layer_error, output_rds = normalize_slashes(layer_rds), 
            stringsAsFactors = FALSE)
        write_csv_utf8(status, layer_status_path)
        layer_status_rows[[length(layer_status_rows) + 1L]] <- status
        log_event("FAILED counts layer", layer_name, "error=", layer_error)
        stop("Counts-layer aggregation failed. Completed layers are safely saved. ", "Re-run the same script after reviewing the error.")
    }
    saved_object <- result$saved_object
    status <- data.frame(layer_index = layer_index, layer_name = layer_name, status = "COMPLETE", cell_n = length(saved_object$cell_barcodes), 
        sample_n = length(unique(saved_object$metadata$sample_accession)), broad_type_n = length(unique(saved_object$metadata$broad_cell_type)), 
        pseudobulk_n = ncol(saved_object$counts), raw_count_total = saved_object$raw_count_total, elapsed_minutes = elapsed_minutes, 
        error_message = "", output_rds = normalize_slashes(layer_rds), stringsAsFactors = FALSE)
    write_csv_utf8(status, layer_status_path)
    layer_status_rows[[length(layer_status_rows) + 1L]] <- status
    processed_cell_lists[[length(processed_cell_lists) + 1L]] <- saved_object$cell_barcodes
    raw_count_total_from_layers <- raw_count_total_from_layers + saved_object$raw_count_total
    log_event("DONE counts layer", layer_index, "/", length(count_layers), layer_name, "cells=", length(saved_object$cell_barcodes), 
        "pseudobulks=", ncol(saved_object$counts), "elapsed_min=", elapsed_minutes)
    rm(result, saved_object, count_matrix, membership, pb_counts)
    gc()
}

layer_status <- do.call(rbind, layer_status_rows)

layer_status <- layer_status[order(layer_status$layer_index), , drop = FALSE]

write_csv_utf8(layer_status, file.path(table_root, "01_layer_aggregation_status.csv"))

if (nrow(layer_status) != 28L || any(layer_status$status != "COMPLETE")) {
    stop("Not all 28 counts layers completed aggregation.")
}

processed_cells <- unlist(processed_cell_lists, use.names = FALSE)

cell_identity_audit <- data.frame(metric = c("processed_cell_n", "unique_processed_cell_n", "annotation_cell_n", 
    "duplicated_processed_cell_n", "missing_annotation_cell_n", "unexpected_processed_cell_n"), value = c(length(processed_cells), 
    length(unique(processed_cells)), nrow(annotation), sum(duplicated(processed_cells)), length(setdiff(annotation$cell_barcode, 
        processed_cells)), length(setdiff(processed_cells, annotation$cell_barcode))), stringsAsFactors = FALSE)

write_csv_utf8(cell_identity_audit, file.path(run_qc_root, "03_cell_identity_audit.csv"))

if (length(processed_cells) != 236636L || length(unique(processed_cells)) != 236636L || anyDuplicated(processed_cells) > 
    0L || !setequal(processed_cells, annotation$cell_barcode)) {
    stop("Counts-layer cell identity audit failed.")
}

rm(object, annotation_lookup, processed_cell_lists, processed_cells)

gc()

layer_rds_files <- layer_status$output_rds

saved_layers <- lapply(layer_rds_files, readRDS)

reference_genes <- rownames(saved_layers[[1]]$counts)

gene_identity_pass <- vapply(saved_layers, function(x) {
    identical(rownames(x$counts), reference_genes)
}, logical(1))

if (!all(gene_identity_pass)) {
    stop("Gene order differs across saved counts-layer pseudobulks.")
}

pb_matrices <- lapply(saved_layers, function(x) x$counts)

pb_metadata_list <- lapply(saved_layers, function(x) x$metadata)

pb_counts <- do.call(cbind, pb_matrices)

pb_metadata <- do.call(rbind, pb_metadata_list)

rownames(pb_metadata) <- NULL

if (anyDuplicated(colnames(pb_counts)) > 0L || anyDuplicated(pb_metadata$pb_id) > 0L) {
    stop("Duplicated pseudobulk IDs were created.")
}

pb_metadata <- pb_metadata[match(colnames(pb_counts), pb_metadata$pb_id), , drop = FALSE]

if (anyNA(pb_metadata$pb_id) || !identical(colnames(pb_counts), as.character(pb_metadata$pb_id))) {
    stop("Combined pseudobulk metadata alignment failed.")
}

library_size <- as.numeric(Matrix::colSums(pb_counts))

detected_gene_n <- as.numeric(Matrix::colSums(pb_counts > 0))

pb_metadata$library_size <- library_size

pb_metadata$detected_gene_n <- detected_gene_n

pb_metadata$minimum_cell_threshold_pass <- (pb_metadata$cell_n >= 20L)

pb_metadata$preferred_cell_threshold_pass <- (pb_metadata$cell_n >= 50L)

combined_raw_count_total <- as.numeric(sum(pb_counts))

if (!isTRUE(all.equal(raw_count_total_from_layers, combined_raw_count_total, tolerance = 0))) {
    stop("Combined raw-count total mismatch: ", raw_count_total_from_layers, " vs ", combined_raw_count_total)
}

pb_counts_path <- file.path(object_root, "pb_raw_counts.rds")

pb_metadata_path <- file.path(table_root, "02_pb_metadata.csv")

saveRDS(pb_counts, pb_counts_path, compress = FALSE)

write_csv_utf8(pb_metadata, pb_metadata_path)

combined_object_status <- data.frame(dataset_id = dataset_id, pseudobulk_version = pseudobulk_version, status = "COMPLETE", 
    gene_n = nrow(pb_counts), pseudobulk_n = ncol(pb_counts), sample_n = length(unique(pb_metadata$sample_accession)), 
    donor_n = length(unique(pb_metadata$donor_id)), broad_type_n = length(unique(pb_metadata$broad_cell_type)), 
    aggregated_cell_n = sum(pb_metadata$cell_n), raw_count_total = combined_raw_count_total, counts_file_size_bytes = file.info(pb_counts_path)$size, 
    counts_file_md5 = unname(tools::md5sum(pb_counts_path)), counts_path = normalize_slashes(pb_counts_path), metadata_path = normalize_slashes(pb_metadata_path), 
    stringsAsFactors = FALSE)

write_csv_utf8(combined_object_status, file.path(object_root, "pb_status.csv"))

rm(saved_layers, pb_matrices, pb_metadata_list)

gc()

broad_types <- sort(unique(pb_metadata$broad_cell_type))

pb_qc_by_type <- do.call(rbind, lapply(broad_types, function(cell_type) {
    one <- pb_metadata[pb_metadata$broad_cell_type == cell_type, , drop = FALSE]
    data.frame(broad_cell_type = cell_type, pseudobulk_n = nrow(one), sample_n = length(unique(one$sample_accession)), 
        donor_n = length(unique(one$donor_id)), median_cell_n = stats::median(one$cell_n), min_cell_n = min(one$cell_n), 
        median_library_size = stats::median(one$library_size), min_library_size = min(one$library_size), median_detected_gene_n = stats::median(one$detected_gene_n), 
        pseudobulk_n_cells_ge_20 = sum(one$cell_n >= 20L), pseudobulk_n_cells_ge_50 = sum(one$cell_n >= 50L), stringsAsFactors = FALSE)
}))

write_csv_utf8(pb_qc_by_type, file.path(table_root, "03_pb_QC_by_type.csv"))

eligibility_rows <- list()

for (cell_type in broad_types) {
    one <- pb_metadata[pb_metadata$broad_cell_type == cell_type & pb_metadata$cell_n >= 20L, , drop = FALSE]
    fertile_temporal <- one[one$clinical_group == "FERTILE" & one$temporal_fertile_eligible, , drop = FALSE]
    day_counts <- table(factor(fertile_temporal$LH_day, levels = c(3, 5, 7, 9, 11)))
    temporal_eligible <- (nrow(fertile_temporal) >= 12L && all(day_counts >= 2L))
    lh7 <- one[one$LH7_RIF_comparison_eligible, , drop = FALSE]
    fertile_lh7_n <- sum(lh7$LH7_reference_group == "FERTILE_LH7", na.rm = TRUE)
    rif_lh7_n <- sum(lh7$LH7_reference_group == "RIF", na.rm = TRUE)
    lh7_eligible <- (fertile_lh7_n >= 3L && rif_lh7_n >= 3L)
    localization_eligible <- (nrow(one) >= 5L && length(unique(one$donor_id)) >= 5L)
    eligibility_rows[[length(eligibility_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, adequate_pseudobulk_n = nrow(one), 
        adequate_donor_n = length(unique(one$donor_id)), fertile_LH3_n = as.integer(day_counts["3"]), fertile_LH5_n = as.integer(day_counts["5"]), 
        fertile_LH7_n = as.integer(day_counts["7"]), fertile_LH9_n = as.integer(day_counts["9"]), fertile_LH11_n = as.integer(day_counts["11"]), 
        RIF_LH7_n = rif_lh7_n, localization_eligible = localization_eligible, fertile_temporal_model_eligible = temporal_eligible, 
        LH7_RIF_vs_fertile_model_eligible = lh7_eligible, stringsAsFactors = FALSE)
}

model_eligibility <- do.call(rbind, eligibility_rows)

write_csv_utf8(model_eligibility, file.path(table_root, "04_model_eligibility.csv"))

program_definitions <- list(RNA_ALL_436 = candidate_table$feature_id, OCCUPANCY_SUPPORTED_284 = candidate_table$feature_id[candidate_table$occupancy_supported], 
    SHARED_OCCUPANCY_46 = candidate_table$feature_id[candidate_table$shared_occupancy_supported], PROMOTER_OCCUPANCY_84 = candidate_table$feature_id[candidate_table$any_promoter_overlap_support])

all_candidate_genes <- unique(candidate_table$feature_id)

present_candidate_genes <- intersect(all_candidate_genes, rownames(pb_counts))

candidate_counts <- as.matrix(pb_counts[present_candidate_genes, , drop = FALSE])

safe_library_size <- library_size

safe_library_size[safe_library_size <= 0] <- NA_real_

candidate_log_cpm <- log2(sweep(candidate_counts, MARGIN = 2, STATS = safe_library_size/1000000, FUN = "/") + 1)

candidate_z <- row_zscore_dense(candidate_log_cpm)

score_rows <- list()

coverage_rows <- list()

for (program_name in names(program_definitions)) {
    frozen_genes <- unique(program_definitions[[program_name]])
    measurable_genes <- intersect(frozen_genes, rownames(candidate_z))
    candidate_subset <- candidate_table[match(measurable_genes, candidate_table$feature_id), , drop = FALSE]
    weights <- sign(candidate_subset$interaction_logFC)
    valid <- is.finite(weights) & weights != 0
    measurable_genes <- measurable_genes[valid]
    weights <- weights[valid]
    if (length(measurable_genes) < 10L) {
        stop("Program ", program_name, " has fewer than 10 measurable genes.")
    }
    signed_matrix <- sweep(candidate_z[measurable_genes, , drop = FALSE], MARGIN = 1, STATS = weights, FUN = "*")
    score_value <- colMeans(signed_matrix, na.rm = TRUE)
    score_rows[[program_name]] <- data.frame(pb_id = colnames(pb_counts), program_name = program_name, score = as.numeric(score_value), 
        frozen_gene_n = length(frozen_genes), measurable_gene_n = length(measurable_genes), coverage_fraction = length(measurable_genes)/length(frozen_genes), 
        score_definition = "MEAN_SIGNED_GENE_ZSCORE_FROM_LOG2_CPM_PLUS_1", interpretation = "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION", 
        stringsAsFactors = FALSE)
    coverage_rows[[program_name]] <- data.frame(program_name = program_name, frozen_gene_n = length(frozen_genes), 
        measurable_gene_n = length(measurable_genes), coverage_fraction = length(measurable_genes)/length(frozen_genes), 
        missing_genes = collapse_unique(setdiff(frozen_genes, measurable_genes)), stringsAsFactors = FALSE)
}

program_scores <- do.call(rbind, score_rows)

program_coverage <- do.call(rbind, coverage_rows)

program_scores <- merge(program_scores, pb_metadata, by = "pb_id", all.x = TRUE, sort = FALSE)

if (anyNA(program_scores$broad_cell_type) || nrow(program_scores) != nrow(pb_metadata) * length(program_definitions)) {
    stop("Program-score metadata merge failed.")
}

program_scores$score_eligible <- (program_scores$cell_n >= 20L)

write_csv_utf8(program_scores, file.path(table_root, "05_VDR_program_scores.csv"))

write_csv_utf8(program_coverage, file.path(table_root, "06_VDR_program_coverage.csv"))

localization_summary <- do.call(rbind, lapply(split(program_scores[program_scores$score_eligible, , drop = FALSE], 
    list(program_scores$program_name[program_scores$score_eligible], program_scores$broad_cell_type[program_scores$score_eligible]), 
    drop = TRUE), function(one) {
    data.frame(program_name = one$program_name[1], broad_cell_type = one$broad_cell_type[1], pseudobulk_n = nrow(one), 
        donor_n = length(unique(one$donor_id)), mean_score = mean(one$score, na.rm = TRUE), median_score = stats::median(one$score, 
            na.rm = TRUE), score_sd = stats::sd(one$score, na.rm = TRUE), score_IQR = stats::IQR(one$score, na.rm = TRUE), 
        min_score = min(one$score, na.rm = TRUE), max_score = max(one$score, na.rm = TRUE), median_cell_n = stats::median(one$cell_n), 
        stringsAsFactors = FALSE)
}))

localization_summary <- localization_summary[order(localization_summary$program_name, -localization_summary$median_score), 
    , drop = FALSE]

write_csv_utf8(localization_summary, file.path(table_root, "07_VDR_localization_summary.csv"))

if ("VDR" %in% rownames(pb_counts)) {
    VDR_counts <- as.numeric(pb_counts["VDR", , drop = TRUE])
    VDR_log_cpm <- log2(VDR_counts/(safe_library_size/1000000) + 1)
    VDR_expression <- pb_metadata
    VDR_expression$VDR_raw_count <- VDR_counts
    VDR_expression$VDR_log2_CPM_plus_1 <- VDR_log_cpm
    write_csv_utf8(VDR_expression, file.path(table_root, "08_VDR_expression.csv"))
} else {
    VDR_expression <- data.frame()
    write_text("VDR was not present in the pseudobulk count matrix row names.", file.path(run_qc_root, "VDR_gene_absent_note.txt"))
}

theme_qc <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

plot_score_data <- program_scores[program_scores$score_eligible, , drop = FALSE]

plot_score_data$broad_cell_type <- factor(plot_score_data$broad_cell_type, levels = broad_types)

plot_scores <- ggplot2::ggplot(plot_score_data, ggplot2::aes(x = broad_cell_type, y = score)) + ggplot2::geom_boxplot(outlier.shape = NA, 
    width = 0.65000000000000002) + ggplot2::geom_jitter(ggplot2::aes(shape = clinical_group), width = 0.14999999999999999, 
    height = 0, size = 1.3999999999999999, alpha = 0.65000000000000002) + ggplot2::facet_wrap(~program_name, scales = "free_y", 
    ncol = 2) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + ggplot2::labs(title = "VDR-program localization across broad cell types", 
    subtitle = paste("Signed standardized pseudobulk scores;", "descriptive localization only"), x = "Broad cell type", 
    y = "Signed VDR-program score", shape = "Clinical group") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_scores, "QC1_VDR_program_scores", width = 12, height = 9)

heatmap_data <- localization_summary

heatmap_data$broad_cell_type <- factor(heatmap_data$broad_cell_type, levels = rev(broad_types))

plot_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = program_name, y = broad_cell_type, fill = median_score)) + 
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", median_score)), size = 3) + 
    ggplot2::scale_fill_gradient2(midpoint = 0) + ggplot2::labs(title = "Median VDR-program localization score", 
    subtitle = "Positive values align with the locked VDR interaction direction", x = "VDR program", y = "Broad cell type", 
    fill = "Median\nscore") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

save_plot_pair(plot_heatmap, "QC2_VDR_localization_heatmap", width = 9, height = 6)

plot_pb_cells <- ggplot2::ggplot(pb_metadata, ggplot2::aes(x = broad_cell_type, y = cell_n, shape = clinical_group)) + 
    ggplot2::geom_hline(yintercept = c(20, 50), linetype = c("solid", "dashed")) + ggplot2::geom_jitter(width = 0.14999999999999999, 
    height = 0, size = 1.8, alpha = 0.69999999999999996) + ggplot2::scale_y_log10(labels = scales::comma) + ggplot2::labs(title = "Pseudobulk cell-count QC", 
    subtitle = "Solid line: 20-cell minimum; dashed line: 50-cell preferred threshold", x = "Broad cell type", 
    y = "Aggregated cell count (log10 scale)", shape = "Clinical group") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_pb_cells, "QC3_pseudobulk_cell_counts", width = 10, height = 6)

plot_eligibility_data <- reshape(model_eligibility[, c("broad_cell_type", "localization_eligible", "fertile_temporal_model_eligible", 
    "LH7_RIF_vs_fertile_model_eligible")], varying = c("localization_eligible", "fertile_temporal_model_eligible", 
    "LH7_RIF_vs_fertile_model_eligible"), v.names = "eligible", timevar = "analysis", times = c("Localization", 
    "Fertile temporal", "LH7 RIF vs fertile"), direction = "long")

plot_eligibility <- ggplot2::ggplot(plot_eligibility_data, ggplot2::aes(x = analysis, y = broad_cell_type, fill = eligible)) + 
    ggplot2::geom_tile() + ggplot2::labs(title = "Broad-cell-type model eligibility", x = "Planned analysis", y = "Broad cell type", 
    fill = "Eligible") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

save_plot_pair(plot_eligibility, "QC4_model_eligibility", width = 8, height = 6)

if (nrow(VDR_expression)) {
    plot_VDR <- ggplot2::ggplot(VDR_expression[VDR_expression$cell_n >= 20L, , drop = FALSE], ggplot2::aes(x = broad_cell_type, 
        y = VDR_log2_CPM_plus_1, shape = clinical_group)) + ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65000000000000002) + 
        ggplot2::geom_jitter(width = 0.14999999999999999, height = 0, size = 1.5, alpha = 0.65000000000000002) + 
        ggplot2::labs(title = "VDR pseudobulk expression across broad cell types", subtitle = "Descriptive expression localization", 
            x = "Broad cell type", y = "VDR log2(CPM + 1)", shape = "Clinical group") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
        hjust = 1))
    save_plot_pair(plot_VDR, "QC5_VDR_expression", width = 10, height = 6)
}

technical_checks <- data.frame(check = c("annotation_gate_authorizes_broad_analysis", "checkpoint_complete", "counts_layer_n_is_28", 
    "all_counts_layers_aggregated", "all_236636_cells_aggregated_once", "raw_count_total_preserved", "sample_n_is_28", 
    "donor_n_is_28", "broad_type_n_is_8", "four_VDR_programs_scored", "all_program_coverage_ge_0_90", "at_least_one_type_eligible_for_temporal_model", 
    "at_least_one_type_eligible_for_LH7_model"), pass = c(annotation_gate$gate_status[[1]] == "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW", 
    checkpoint_status$checkpoint_status[[1]] == "COMPLETE", length(count_layers) == 28L, nrow(layer_status) == 
        28L && all(layer_status$status == "COMPLETE"), sum(pb_metadata$cell_n) == 236636L, isTRUE(all.equal(raw_count_total_from_layers, 
        combined_raw_count_total, tolerance = 0)), length(unique(pb_metadata$sample_accession)) == 28L, length(unique(pb_metadata$donor_id)) == 
        28L, length(unique(pb_metadata$broad_cell_type)) == 8L, length(unique(program_scores$program_name)) == 
        4L, all(program_coverage$coverage_fraction >= 0.90000000000000002), any(model_eligibility$fertile_temporal_model_eligible), 
    any(model_eligibility$LH7_RIF_vs_fertile_model_eligible)), stringsAsFactors = FALSE)

write_csv_utf8(technical_checks, file.path(run_qc_root, "04_gate_checks.csv"))

technical_pass <- all(technical_checks$pass)

gate_status <- if (technical_pass) {
    "PASS_PSEUDOBULK_AND_VDR_LOCALIZATION_READY_FOR_DONOR_MODELS"
} else {
    "HOLD_PSEUDOBULK_OR_LOCALIZATION_QC_FAILED"
}

gate_table <- data.frame(dataset_id = dataset_id, input_version = input_version, annotation_version = annotation_version, 
    pseudobulk_version = pseudobulk_version, gate_status = gate_status, gene_n = nrow(pb_counts), pseudobulk_n = ncol(pb_counts), 
    sample_n = length(unique(pb_metadata$sample_accession)), donor_n = length(unique(pb_metadata$donor_id)), broad_type_n = length(unique(pb_metadata$broad_cell_type)), 
    aggregated_cell_n = sum(pb_metadata$cell_n), low_cell_pseudobulk_n = sum(pb_metadata$cell_n < 20L), preferred_cell_threshold_fail_n = sum(pb_metadata$cell_n < 
        50L), program_n = length(unique(program_scores$program_name)), minimum_program_coverage = min(program_coverage$coverage_fraction), 
    temporal_model_eligible_type_n = sum(model_eligibility$fertile_temporal_model_eligible), LH7_model_eligible_type_n = sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible), 
    donor_level_modeling_authorized = technical_pass, cell_composition_inference_authorized = FALSE, subtype_modeling_authorized = FALSE, 
    result_lock_authorized = FALSE, stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "PBVDR_GATE.csv")

write_csv_utf8(gate_table, gate_path)

card_lines <- c("# STEP04E1F Pseudobulk and VDR Localization Result Card", "", paste0("- Dataset: ", dataset_id), 
    paste0("- Input version: ", input_version), paste0("- Annotation version: ", annotation_version), paste0("- Pseudobulk version: ", 
        pseudobulk_version), paste0("- Gate: ", gate_status), paste0("- Aggregated cells: ", format(sum(pb_metadata$cell_n), 
        big.mark = ",")), paste0("- Pseudobulks: ", ncol(pb_counts)), paste0("- Broad cell types: ", length(unique(pb_metadata$broad_cell_type))), 
    paste0("- Minimum VDR-program coverage: ", round(100 * min(program_coverage$coverage_fraction), 1), "%"), paste0("- Temporal-model eligible broad types: ", 
        sum(model_eligibility$fertile_temporal_model_eligible)), paste0("- LH7-model eligible broad types: ", sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible)), 
    paste0("- Donor-level modelling authorized: ", gate_table$donor_level_modeling_authorized), "- Cell-composition inference authorized: FALSE", 
    "- Subtype modelling authorized: FALSE", "- Result lock authorized: FALSE")

card_path <- file.path(run_qc_root, "PBVDR_RESULT_CARD.md")

write_text(card_lines, card_path)

interpretation_lines <- c("# STEP04E1F interpretation", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: **", 
    gate_status, "**"), "", "## Pseudobulk aggregation", "", paste0("- All ", format(sum(pb_metadata$cell_n), big.mark = ","), 
    " cells were aggregated exactly once across 28 samples and eight broad cell types."), "- Raw count totals were preserved during layer-wise aggregation.", 
    "- Raw pseudobulk counts were saved for subsequent donor-level modelling.", "", "## VDR localization", "", 
    "- Four locked VDR RNA programs were scored using signed standardized pseudobulk expression.", "- Scores provide relative cell-type localization of expression aligned with the locked VDR interaction direction.", 
    "- Scores do not establish direct VDR binding, causal activation or clinical differential expression.", "", 
    "## Inferential boundary", "", "- Only broad cell types passing donor/cell-count eligibility may enter temporal or LH+7 clinical models.", 
    "- Cell-composition hypothesis testing remains unauthorized.", "- Candidate subtype labels remain excluded from inferential analysis.", 
    "", "## Next step", "", "- Fit donor-level fertile temporal models for eligible broad cell types.", "- Fit LH+7 RIF-versus-fertile models for eligible broad cell types.", 
    "- Control multiple testing across broad cell types and programs.", "- Preserve VDR localization and donor-level clinical effects as distinct evidence layers.")

interpretation_path <- file.path(run_qc_root, "PBVDR_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

methods_lines <- c("# STEP04E1F Methods draft", "", paste("Raw RNA counts were aggregated separately within each split counts", 
    "layer using the accepted sample-level broad cell-type annotations."), paste("For each sample and broad cell type, a sparse cell-membership matrix", 
    "was multiplied by the raw gene-by-cell count matrix to obtain donor-level", "pseudobulk counts."), paste("Aggregation was performed layer by layer, and each intermediate", 
    "pseudobulk object was saved immediately to support reproducible restart."), paste("Pseudobulk profiles containing fewer than 20 cells were excluded from", 
    "program-score summaries and were retained with explicit QC flags."), paste("For descriptive VDR-program localization, pseudobulk counts were", 
    "converted to log2 counts per million plus one."), paste("Each program gene was standardized across all donor-by-cell-type", 
    "pseudobulk profiles, multiplied by the sign of its locked VDR", "interaction log fold-change, and averaged within each prespecified", 
    "program."), paste("The four programs comprised the complete 436-gene RNA interaction", "program, the 284-gene occupancy-supported subset, the 46-gene", 
    "shared-occupancy subset and the 84-gene promoter-occupancy subset."), paste("These scores were used as descriptive localization measures and not", 
    "as evidence of causal VDR activation or donor-level clinical effects."))

write_text(methods_lines, file.path(manuscript_text_root, "STEP04E1F_PBVDR_Methods_draft.md"))

results_lines <- c("# STEP04E1F Results shell", "", paste0("Raw counts from ", format(sum(pb_metadata$cell_n), 
    big.mark = ","), " cells were aggregated across ", length(unique(pb_metadata$sample_accession)), " samples and ", 
    length(unique(pb_metadata$broad_cell_type)), " broad cell populations."), paste0("The resulting pseudobulk matrix contained ", 
    format(nrow(pb_counts), big.mark = ","), " genes and ", ncol(pb_counts), " donor-by-cell-type profiles."), 
    paste0("All four locked VDR programs were measurable, with a minimum gene", " coverage of ", round(100 * min(program_coverage$coverage_fraction), 
        1), "%."), paste0(sum(model_eligibility$fertile_temporal_model_eligible), " broad cell types met the predefined fertile temporal modelling", 
        " criteria, and ", sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible), " met the LH+7 RIF-versus-fertile criteria."), 
    paste("Cell-type localization scores were retained as a descriptive evidence", "layer, whereas clinical and temporal comparisons were deferred to", 
        "donor-level models."))

write_text(results_lines, file.path(manuscript_text_root, "STEP04E1F_PBVDR_Results_shell.md"))

run_end <- Sys.time()

package_names <- c("SeuratObject", "Matrix", "data.table", "ggplot2", "jsonlite", "scales", "future")

package_table <- data.frame(package = package_names, installed = vapply(package_names, requireNamespace, logical(1), 
    quietly = TRUE), version = vapply(package_names, safe_package_version, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- unique(c(list.files(output_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1F_PBVDR_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1F_PBVDR_Results_shell.md")))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, input_version = input_version, 
    annotation_version = annotation_version, pseudobulk_version = pseudobulk_version, dataset_id = dataset_id, 
    started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, 
        run_start, units = "mins")), 2), gate_status = gate_status, gene_n = nrow(pb_counts), pseudobulk_n = ncol(pb_counts), 
    sample_n = length(unique(pb_metadata$sample_accession)), donor_n = length(unique(pb_metadata$donor_id)), broad_type_n = length(unique(pb_metadata$broad_cell_type)), 
    aggregated_cell_n = sum(pb_metadata$cell_n), counts_layer_n = length(count_layers), program_n = length(unique(program_scores$program_name)), 
    donor_level_modeling_authorized = gate_table$donor_level_modeling_authorized, composition_inference_authorized = FALSE, 
    subtype_modeling_authorized = FALSE, result_lock_authorized = FALSE, output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "pseudobulk_n=", ncol(pb_counts), "temporal_eligible_types=", 
    sum(model_eligibility$fertile_temporal_model_eligible), "LH7_eligible_types=", sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(source_data_root, 
    recursive = TRUE, full.names = TRUE), list.files(run_qc_root, recursive = TRUE, full.names = TRUE), file.path(object_root, 
    "pb_status.csv"), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1F_PBVDR_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1F_PBVDR_Results_shell.md"), normalizePath(.vdr_script_path(script_name), 
    winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    200 * 1024^2]

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

rm(candidate_counts, candidate_log_cpm, candidate_z, pb_counts)

gc()

cat("\n============================================================\n")

cat("STEP04E1F GSE250130 VDR localization/pseudobulk completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Aggregated cells: ", format(sum(pb_metadata$cell_n), big.mark = ","), "\n", sep = "")

cat("Pseudobulks: ", nrow(pb_metadata), "\n", sep = "")

cat("Broad cell types: ", length(unique(pb_metadata$broad_cell_type)), "\n", sep = "")

cat("Temporal-model eligible types: ", sum(model_eligibility$fertile_temporal_model_eligible), "\n", sep = "")

cat("LH7-model eligible types: ", sum(model_eligibility$LH7_RIF_vs_fertile_model_eligible), "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe output and QC folders are complete.\n", sep = "")
}

