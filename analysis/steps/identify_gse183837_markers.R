#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D04: Calculate descriptive cluster markers

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

options(future.globals.maxSize = 8 * 1024^3)

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "identify_gse183837_markers.R"

script_version <- "STEP04E1L_GSE183837_DESCRIPTIVE_MARKERS_v1.1_PANEL_DUPLICATE_FIXED"

integration_version <- "SCI183_v1_0"

marker_version <- "MK183_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1L_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

safe_package_version <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        return(NA_character_)
    }
    as.character(utils::packageVersion(package_name))
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

median_or_NA <- function(x) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        median(x)
    }
}

safe_zscore <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 2L || !is.finite(stats::sd(x)) || stats::sd(x) == 0) {
        return(rep(0, length(x)))
    }
    as.numeric(scale(x))
}

detect_logFC_column <- function(column_names) {
    preferred <- c("avg_log2FC", "avg_logFC", "avg_diff")
    hit <- preferred[preferred %in% column_names]
    if (length(hit)) {
        return(hit[1])
    }
    pattern_hit <- grep("log.*fc|avg.*diff", column_names, value = TRUE, ignore.case = TRUE)
    if (length(pattern_hit)) {
        return(pattern_hit[1])
    }
    ""
}

standardize_marker_table <- function(marker_table, cluster_id) {
    if (is.null(marker_table) || !nrow(marker_table)) {
        return(data.frame(cluster = character(), gene = character(), p_val = numeric(), avg_log2FC = numeric(), 
            pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), rank = integer(), stringsAsFactors = FALSE))
    }
    output <- as.data.frame(marker_table, stringsAsFactors = FALSE)
    output$gene <- rownames(output)
    rownames(output) <- NULL
    logFC_column <- detect_logFC_column(colnames(output))
    if (!nzchar(logFC_column)) {
        stop("Could not identify a log-fold-change column for cluster ", cluster_id, ".")
    }
    required_columns <- c("p_val", "pct.1", "pct.2", "p_val_adj")
    missing_columns <- setdiff(required_columns, colnames(output))
    if (length(missing_columns)) {
        stop("Marker table for cluster ", cluster_id, " is missing: ", paste(missing_columns, collapse = ", "))
    }
    output$avg_log2FC <- as.numeric(output[[logFC_column]])
    output$p_val <- as.numeric(output$p_val)
    output$pct.1 <- as.numeric(output$pct.1)
    output$pct.2 <- as.numeric(output$pct.2)
    output$p_val_adj <- as.numeric(output$p_val_adj)
    output$cluster <- as.character(cluster_id)
    output <- output[order(output$p_val_adj, -output$avg_log2FC, output$gene), , drop = FALSE]
    output$rank <- seq_len(nrow(output))
    preferred_order <- c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "rank")
    remaining_columns <- setdiff(colnames(output), preferred_order)
    output[, c(preferred_order, remaining_columns), drop = FALSE]
}

resolve_features <- function(requested_features, available_features) {
    requested_features <- as.character(requested_features)
    available_features <- as.character(available_features)
    available_upper <- toupper(available_features)
    requested_upper <- toupper(requested_features)
    match_index <- match(requested_upper, available_upper)
    output <- data.frame(request_index = seq_along(requested_features), requested_gene = requested_features, resolved_feature = ifelse(is.na(match_index), 
        NA_character_, available_features[match_index]), available = !is.na(match_index), stringsAsFactors = FALSE)
    if (nrow(output) != length(requested_features)) {
        stop("Feature-resolution output length does not match the requested panel length.")
    }
    output
}

get_data_matrix <- function(object, features) {
    features <- unique(features[features %in% rownames(object)])
    if (!length(features)) {
        return(NULL)
    }
    output <- SeuratObject::LayerData(object = object, assay = "RNA", layer = "data", features = features)
    output[features, , drop = FALSE]
}

calculate_donor_marker_support <- function(object, metadata, cluster_id, marker_genes, minimum_cluster_cells_per_donor = 3L) {
    marker_genes <- unique(marker_genes[marker_genes %in% rownames(object)])
    if (!length(marker_genes)) {
        return(data.frame(cluster = character(), gene = character(), donor_id = character(), analysis_group = character(), 
            cluster_cell_n = integer(), outside_cell_n = integer(), pct_cluster = numeric(), pct_outside = numeric(), 
            pct_delta = numeric(), avg_cluster = numeric(), avg_outside = numeric(), avg_delta = numeric(), eligible_donor = logical(), 
            donor_supports_marker = logical(), stringsAsFactors = FALSE))
    }
    expression_matrix <- get_data_matrix(object, marker_genes)
    donor_ids <- unique(as.character(metadata$donor_id))
    output_rows <- list()
    row_index <- 0L
    for (donor_id_value in donor_ids) {
        donor_index <- which(metadata$donor_id == donor_id_value)
        cluster_index <- donor_index[as.character(metadata$seurat_clusters[donor_index]) == as.character(cluster_id)]
        outside_index <- donor_index[as.character(metadata$seurat_clusters[donor_index]) != as.character(cluster_id)]
        cluster_cell_n <- length(cluster_index)
        outside_cell_n <- length(outside_index)
        analysis_group_value <- unique(as.character(metadata$analysis_group[donor_index]))
        if (length(analysis_group_value) != 1L) {
            stop("Donor/group identity is not unique for ", donor_id_value, ".")
        }
        if (cluster_cell_n > 0L) {
            pct_cluster <- Matrix::rowMeans(expression_matrix[, cluster_index, drop = FALSE] > 0)
            avg_cluster <- Matrix::rowMeans(expression_matrix[, cluster_index, drop = FALSE])
        }
        else {
            pct_cluster <- rep(NA_real_, length(marker_genes))
            avg_cluster <- rep(NA_real_, length(marker_genes))
        }
        if (outside_cell_n > 0L) {
            pct_outside <- Matrix::rowMeans(expression_matrix[, outside_index, drop = FALSE] > 0)
            avg_outside <- Matrix::rowMeans(expression_matrix[, outside_index, drop = FALSE])
        }
        else {
            pct_outside <- rep(NA_real_, length(marker_genes))
            avg_outside <- rep(NA_real_, length(marker_genes))
        }
        eligible_donor <- (cluster_cell_n >= minimum_cluster_cells_per_donor && outside_cell_n >= 20L)
        donor_supports_marker <- (eligible_donor & is.finite(pct_cluster) & is.finite(pct_outside) & pct_cluster - 
            pct_outside >= 0.10000000000000001)
        row_index <- row_index + 1L
        output_rows[[row_index]] <- data.frame(cluster = as.character(cluster_id), gene = marker_genes, donor_id = donor_id_value, 
            analysis_group = analysis_group_value, cluster_cell_n = cluster_cell_n, outside_cell_n = outside_cell_n, 
            pct_cluster = as.numeric(pct_cluster), pct_outside = as.numeric(pct_outside), pct_delta = as.numeric(pct_cluster - 
                pct_outside), avg_cluster = as.numeric(avg_cluster), avg_outside = as.numeric(avg_outside), avg_delta = as.numeric(avg_cluster - 
                avg_outside), eligible_donor = eligible_donor, donor_supports_marker = as.logical(donor_supports_marker), 
            stringsAsFactors = FALSE)
    }
    do.call(rbind, output_rows)
}

summarize_donor_marker_support <- function(support_long) {
    if (is.null(support_long) || !nrow(support_long)) {
        return(data.frame(cluster = character(), gene = character(), eligible_donor_n = integer(), supporting_donor_n = integer(), 
            supporting_donor_fraction = numeric(), median_pct_delta = numeric(), median_avg_delta = numeric(), 
            Control_supporting_donor_n = integer(), RIF_supporting_donor_n = integer(), stringsAsFactors = FALSE))
    }
    support_DT <- data.table::as.data.table(support_long)
    output <- support_DT[, .(eligible_donor_n = sum(eligible_donor, na.rm = TRUE), supporting_donor_n = sum(donor_supports_marker, 
        na.rm = TRUE), supporting_donor_fraction = if (sum(eligible_donor, na.rm = TRUE) > 0L) {
        sum(donor_supports_marker, na.rm = TRUE)/sum(eligible_donor, na.rm = TRUE)
    }
    else {
        NA_real_
    }, median_pct_delta = median_or_NA(pct_delta[eligible_donor]), median_avg_delta = median_or_NA(avg_delta[eligible_donor]), 
        Control_supporting_donor_n = sum(donor_supports_marker & analysis_group == "CONTROL", na.rm = TRUE), RIF_supporting_donor_n = sum(donor_supports_marker & 
            analysis_group == "RIF", na.rm = TRUE)), by = .(cluster, gene)]
    as.data.frame(output)
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

required_packages <- c("Seurat", "SeuratObject", "Matrix", "data.table", "ggplot2", "scales", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

integration_root <- file.path(results_root, "S4E1_sc", integration_version)

input_checkpoint <- file.path(integration_root, "checkpoints", dataset_id, "integrated.rds")

output_root <- file.path(integration_root, "mk", "G183", marker_version)

model_root <- file.path(output_root, "m")

full_root <- file.path(output_root, "full")

top_root <- file.path(output_root, "top")

panel_root <- file.path(output_root, "panel")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "man")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1L_", format(run_start, "%Y%m%d_%H%M%S"), "_markers"))

figure_root <- file.path(run_qc_root, "fig")

for (path_value in c(execution_root, parameter_root, environment_root, output_root, model_root, full_root, top_root, 
    panel_root, figure_output_root, manifest_root, run_qc_root, figure_root)) {
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

planned_files <- data.frame(artifact = c("input_checkpoint", "output_root", "cluster_00_model", "cluster_15_model", 
    "combined_markers", "combined_top100", "donor_support", "canonical_panel", "candidate_annotation", "gate", 
    "handoff_zip"), path = normalize_slashes(c(input_checkpoint, output_root, file.path(model_root, "C00.rds"), 
    file.path(model_root, "C15.rds"), file.path(full_root, "all_markers.csv"), file.path(top_root, "all_top100.csv"), 
    file.path(top_root, "donor_support_summary.csv"), file.path(panel_root, "canonical_panel.csv"), file.path(panel_root, 
        "candidate_annotation.csv"), file.path(run_qc_root, "G183_MARKER_GATE.csv"), handoff_zip)), stringsAsFactors = FALSE)

planned_files$character_n <- nchar(planned_files$path, type = "chars")

planned_files$path_safe <- (planned_files$character_n <= 235L)

write_csv_utf8(planned_files, file.path(run_qc_root, "00_path_audit.csv"))

if (!all(planned_files$path_safe)) {
    stop("One or more planned marker paths exceed the safe limit.")
}

detailed_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

detailed_qc_dirs <- detailed_qc_dirs[grepl("^STEP04E1K_[0-9]{8}_[0-9]{6}_dQC$", basename(detailed_qc_dirs), perl = TRUE)]

if (!length(detailed_qc_dirs)) {
    stop("No STEP04E1K detailed-QC directory was found.")
}

gate_candidates <- list()

for (directory_path in detailed_qc_dirs) {
    gate_path <- file.path(directory_path, "G183_DETAILED_QC_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        gate_candidates[[length(gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(gate_candidates)) {
    stop("No readable STEP04E1K detailed-QC Gate was found.")
}

gate_candidates <- as.data.frame(data.table::rbindlist(gate_candidates, use.names = TRUE, fill = TRUE))

gate_candidates$source_mtime_POSIX <- as.POSIXct(gate_candidates$source_mtime)

gate_candidates <- gate_candidates[order(gate_candidates$source_mtime_POSIX, decreasing = TRUE), , drop = FALSE]

authorized_index <- which(gate_candidates$gate_status %in% c("PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED", 
    "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED") & toupper(as.character(gate_candidates$descriptive_marker_analysis_authorized)) %in% 
    c("TRUE", "T", "1", "YES"))

if (!length(authorized_index)) {
    stop("No detailed-QC Gate authorizes marker extraction.")
}

upstream_gate <- gate_candidates[authorized_index[1], , drop = FALSE]

write_csv_utf8(upstream_gate, file.path(run_qc_root, "00_upstream_gate.csv"))

if (!file.exists(input_checkpoint)) {
    stop("Integrated checkpoint is missing:\n", input_checkpoint)
}

log_event("Starting", run_id, "checkpoint=", normalize_slashes(input_checkpoint))

integrated_object <- readRDS(input_checkpoint)

if (!inherits(integrated_object, "Seurat")) {
    stop("Input checkpoint is not a Seurat object.")
}

SeuratObject::DefaultAssay(integrated_object) <- "RNA"

metadata <- integrated_object[[]]

required_metadata <- c("sample_accession", "donor_id", "analysis_group", "seurat_clusters")

missing_metadata <- setdiff(required_metadata, colnames(metadata))

if (length(missing_metadata)) {
    stop("Integrated metadata is missing: ", paste(missing_metadata, collapse = ", "))
}

if (ncol(integrated_object) != 47161L || length(unique(as.character(metadata$donor_id))) != 9L || length(unique(as.character(metadata$seurat_clusters))) != 
    16L) {
    stop("Integrated object identity does not match the locked 47,161-cell object.")
}

if (anyDuplicated(rownames(integrated_object)) > 0L) {
    stop("Integrated feature names are not unique.")
}

data_layers <- SeuratObject::Layers(integrated_object, assay = "RNA", search = "^data")

count_layers <- SeuratObject::Layers(integrated_object, assay = "RNA", search = "^counts")

if (length(data_layers) != 1L || data_layers[1] != "data" || length(count_layers) != 9L) {
    stop("Expected one joined data layer and nine split count layers.")
}

metadata$cell_barcode <- rownames(metadata)

metadata$donor_id <- as.character(metadata$donor_id)

metadata$analysis_group <- as.character(metadata$analysis_group)

metadata$seurat_clusters <- factor(as.character(metadata$seurat_clusters), levels = as.character(0:15))

if (anyNA(metadata$seurat_clusters)) {
    stop("Unexpected primary cluster identities.")
}

SeuratObject::Idents(integrated_object) <- metadata$seurat_clusters

checkpoint_audit <- data.frame(metric = c("feature_n", "cell_n", "donor_n", "cluster_n", "data_layer_n", "count_layer_n", 
    "checkpoint_size_GB", "checkpoint_md5"), value = c(nrow(integrated_object), ncol(integrated_object), length(unique(metadata$donor_id)), 
    length(unique(metadata$seurat_clusters)), length(data_layers), length(count_layers), round(file.info(input_checkpoint)$size/1024^3, 
        3), unname(tools::md5sum(input_checkpoint))), stringsAsFactors = FALSE)

write_csv_utf8(checkpoint_audit, file.path(run_qc_root, "01_checkpoint_audit.csv"))

cluster_ids <- as.character(0:15)

marker_parameters <- list(test_use = "wilcox", only_positive = TRUE, minimum_pct = 0.10000000000000001, logFC_threshold = 0.25, 
    maximum_cells_per_identity = 5000L, minimum_cells_feature = 3L, minimum_cells_group = 3L, donor_support_top_n = 30L, 
    top_output_n = 100L, random_seed = 20260806L)

cluster_status_rows <- list()

cluster_marker_objects <- list()

donor_support_long_rows <- list()

donor_support_summary_rows <- list()

for (cluster_index in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[cluster_index]
    cluster_code <- sprintf("C%02d", as.integer(cluster_id))
    checkpoint_path <- file.path(model_root, paste0(cluster_code, ".rds"))
    status_path <- file.path(model_root, paste0(cluster_code, "_status.csv"))
    full_path <- file.path(full_root, paste0(cluster_code, ".csv"))
    top_path <- file.path(top_root, paste0(cluster_code, "_top100.csv"))
    donor_long_path <- file.path(top_root, paste0(cluster_code, "_donor_long.csv"))
    donor_summary_path <- file.path(top_root, paste0(cluster_code, "_donor_summary.csv"))
    existing_status <- safe_read_csv(status_path)
    resume_valid <- FALSE
    cluster_result <- NULL
    if (file.exists(checkpoint_path) && !is.null(existing_status) && nrow(existing_status) == 1L && existing_status$status[[1]] == 
        "COMPLETE") {
        cluster_result <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
        resume_valid <- (is.list(cluster_result) && all(c("markers", "donor_support_long", "donor_support_summary") %in% 
            names(cluster_result)) && identical(as.character(cluster_result$cluster), cluster_id))
    }
    if (resume_valid) {
        log_event("SKIP completed marker cluster", cluster_id, "marker_rows=", nrow(cluster_result$markers))
    }
    else {
        cluster_start <- Sys.time()
        log_event("START marker cluster", cluster_id, "(", cluster_index, "/", length(cluster_ids), ")")
        set.seed(marker_parameters$random_seed + as.integer(cluster_id))
        raw_marker_result <- tryCatch(Seurat::FindMarkers(object = integrated_object, ident.1 = cluster_id, ident.2 = NULL, 
            assay = "RNA", slot = "data", test.use = marker_parameters$test_use, only.pos = marker_parameters$only_positive, 
            min.pct = marker_parameters$minimum_pct, logfc.threshold = marker_parameters$logFC_threshold, min.cells.feature = marker_parameters$minimum_cells_feature, 
            min.cells.group = marker_parameters$minimum_cells_group, max.cells.per.ident = marker_parameters$maximum_cells_per_identity, 
            random.seed = marker_parameters$random_seed + as.integer(cluster_id), densify = FALSE, verbose = FALSE), 
            error = function(e) {
                error_status <- data.frame(cluster = cluster_id, status = "ERROR", error_message = conditionMessage(e), 
                  stringsAsFactors = FALSE)
                write_csv_utf8(error_status, status_path)
                stop("FindMarkers failed for cluster ", cluster_id, ": ", conditionMessage(e))
            })
        standardized_markers <- standardize_marker_table(raw_marker_result, cluster_id)
        if (!nrow(standardized_markers)) {
            stop("No positive marker rows were returned for cluster ", cluster_id, ".")
        }
        top100_markers <- head(standardized_markers, marker_parameters$top_output_n)
        donor_marker_genes <- head(standardized_markers$gene, marker_parameters$donor_support_top_n)
        donor_support_long <- calculate_donor_marker_support(object = integrated_object, metadata = metadata, cluster_id = cluster_id, 
            marker_genes = donor_marker_genes, minimum_cluster_cells_per_donor = 3L)
        donor_support_summary <- summarize_donor_marker_support(donor_support_long)
        standardized_markers <- merge(standardized_markers, donor_support_summary, by = c("cluster", "gene"), all.x = TRUE, 
            sort = FALSE)
        standardized_markers <- standardized_markers[order(standardized_markers$rank), , drop = FALSE]
        standardized_markers$robust_donor_supported <- (standardized_markers$p_val_adj < 0.050000000000000003 & 
            standardized_markers$avg_log2FC >= 0.5 & standardized_markers$pct.1 >= 0.20000000000000001 & is.finite(standardized_markers$supporting_donor_fraction) & 
            standardized_markers$supporting_donor_fraction >= 0.5 & standardized_markers$supporting_donor_n >= 
            2)
        top100_markers <- head(standardized_markers, marker_parameters$top_output_n)
        cluster_result <- list(cluster = cluster_id, markers = standardized_markers, top100 = top100_markers, donor_support_long = donor_support_long, 
            donor_support_summary = donor_support_summary, parameters = marker_parameters)
        saveRDS(cluster_result, checkpoint_path, compress = "gzip")
        write_csv_utf8(standardized_markers, full_path)
        write_csv_utf8(top100_markers, top_path)
        write_csv_utf8(donor_support_long, donor_long_path)
        write_csv_utf8(donor_support_summary, donor_summary_path)
        elapsed_minutes <- round(as.numeric(difftime(Sys.time(), cluster_start, units = "mins")), 2)
        cluster_status <- data.frame(cluster = cluster_id, status = "COMPLETE", marker_row_n = nrow(standardized_markers), 
            FDR05_marker_n = sum(standardized_markers$p_val_adj < 0.050000000000000003, na.rm = TRUE), robust_donor_supported_n = sum(standardized_markers$robust_donor_supported, 
                na.rm = TRUE), donor_support_gene_n = nrow(donor_support_summary), elapsed_minutes = elapsed_minutes, 
            checkpoint_md5 = unname(tools::md5sum(checkpoint_path)), stringsAsFactors = FALSE)
        write_csv_utf8(cluster_status, status_path)
        existing_status <- cluster_status
        log_event("DONE marker cluster", cluster_id, "rows=", nrow(standardized_markers), "FDR05=", sum(standardized_markers$p_val_adj < 
            0.050000000000000003, na.rm = TRUE), "robust=", sum(standardized_markers$robust_donor_supported, na.rm = TRUE), 
            "elapsed_min=", elapsed_minutes)
        rm(raw_marker_result, standardized_markers, donor_support_long, donor_support_summary)
        gc()
    }
    if (is.null(existing_status)) {
        existing_status <- safe_read_csv(status_path)
    }
    cluster_status_rows[[cluster_index]] <- existing_status
    cluster_marker_objects[[cluster_index]] <- cluster_result$markers
    donor_support_long_rows[[cluster_index]] <- cluster_result$donor_support_long
    donor_support_summary_rows[[cluster_index]] <- cluster_result$donor_support_summary
}

cluster_status_all <- do.call(rbind, cluster_status_rows)

cluster_status_all <- cluster_status_all[match(cluster_ids, cluster_status_all$cluster), , drop = FALSE]

write_csv_utf8(cluster_status_all, file.path(output_root, "01_marker_status.csv"))

all_markers <- as.data.frame(data.table::rbindlist(cluster_marker_objects, use.names = TRUE, fill = TRUE))

all_markers$cluster <- as.character(all_markers$cluster)

all_markers <- all_markers[order(as.integer(all_markers$cluster), all_markers$rank), , drop = FALSE]

all_top100 <- all_markers[all_markers$rank <= marker_parameters$top_output_n, , drop = FALSE]

donor_support_long_all <- as.data.frame(data.table::rbindlist(donor_support_long_rows, use.names = TRUE, fill = TRUE))

donor_support_summary_all <- as.data.frame(data.table::rbindlist(donor_support_summary_rows, use.names = TRUE, 
    fill = TRUE))

write_csv_utf8(all_markers, file.path(full_root, "all_markers.csv"))

write_csv_utf8(all_top100, file.path(top_root, "all_top100.csv"))

write_csv_utf8(donor_support_long_all, file.path(top_root, "donor_support_long.csv"))

write_csv_utf8(donor_support_summary_all, file.path(top_root, "donor_support_summary.csv"))

canonical_panels <- list(Unciliated_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "PAX8", "MUC1"), Ciliated_epithelial = c("FOXJ1", 
    "PIFO", "TPPP3", "CAPS", "C20orf85", "RSPH1"), Stromal = c("COL1A1", "COL1A2", "DCN", "COL3A1", "PDGFRA", "COL6A1"), 
    Blood_endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "RAMP2", "PLVAP"), Lymphatic_endothelial = c("CCL21", 
        "PROX1", "LYVE1", "FLT4", "PDPN", "FABP4"), NK = c("NKG7", "GNLY", "KLRD1", "NCAM1", "PRF1", "FGFBP2"), 
    T = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "LTB"), B = c("CD79A", "MS4A1", "CD37", "CD74", "CD22", "HLA-DRA"), 
    Plasma = c("MZB1", "JCHAIN", "IGKC", "IGLC2", "CD79A", "SEC11C"), Myeloid = c("LYZ", "LST1", "FCER1G", "TYROBP", 
        "CTSS", "HLA-DRA"), pDC = c("GZMB", "JCHAIN", "IL3RA", "CLEC4C", "TCF4", "SERPINF1"), Neutrophil = c("S100A8", 
        "S100A9", "FCGR3B", "CSF3R", "FPR1", "MNDA"), Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "HDC", "MS4A2"), 
    Cycling = c("MKI67", "TOP2A", "TUBA1B", "STMN1", "CENPF", "UBE2C"))

panel_requested <- data.frame(lineage = rep(names(canonical_panels), lengths(canonical_panels)), requested_gene = unlist(canonical_panels, 
    use.names = FALSE), stringsAsFactors = FALSE)

panel_resolution <- resolve_features(panel_requested$requested_gene, rownames(integrated_object))

if (nrow(panel_resolution) != nrow(panel_requested) || !identical(panel_resolution$request_index, seq_len(nrow(panel_requested))) || 
    !identical(panel_resolution$requested_gene, panel_requested$requested_gene)) {
    stop("Canonical-panel feature resolution is not one-to-one with the 84 requested panel rows.")
}

panel_requested$request_index <- panel_resolution$request_index

panel_requested$resolved_feature <- panel_resolution$resolved_feature

panel_requested$available <- panel_resolution$available

panel_duplicate_audit <- data.frame(metric = c("requested_panel_row_n", "unique_requested_gene_n", "duplicated_cross_lineage_row_n", 
    "resolved_panel_row_n", "available_panel_row_n", "unavailable_panel_row_n", "one_to_one_assignment_pass"), 
    value = c(nrow(panel_requested), length(unique(panel_requested$requested_gene)), sum(duplicated(panel_requested$requested_gene)), 
        nrow(panel_resolution), sum(panel_requested$available), sum(!panel_requested$available), nrow(panel_resolution) == 
            nrow(panel_requested)), stringsAsFactors = FALSE)

panel_duplicate_rows <- panel_requested[duplicated(panel_requested$requested_gene) | duplicated(panel_requested$requested_gene, 
    fromLast = TRUE), , drop = FALSE]

write_csv_utf8(panel_duplicate_audit, file.path(panel_root, "00_panel_row_audit.csv"))

write_csv_utf8(panel_duplicate_rows, file.path(panel_root, "00_cross_lineage_duplicate_genes.csv"))

write_csv_utf8(panel_requested, file.path(panel_root, "canonical_panel_manifest.csv"))

available_panel <- panel_requested[panel_requested$available, , drop = FALSE]

if (nrow(available_panel) < 0.75 * nrow(panel_requested)) {
    stop("Canonical panel coverage is below 75%.")
}

panel_features <- unique(available_panel$resolved_feature)

panel_expression <- get_data_matrix(integrated_object, panel_features)

panel_rows <- list()

panel_row_index <- 0L

for (cluster_id in cluster_ids) {
    cluster_cells <- which(as.character(metadata$seurat_clusters) == cluster_id)
    cluster_average <- Matrix::rowMeans(panel_expression[, cluster_cells, drop = FALSE])
    cluster_pct <- Matrix::rowMeans(panel_expression[, cluster_cells, drop = FALSE] > 0)
    cluster_panel <- data.frame(cluster = cluster_id, resolved_feature = rownames(panel_expression), average_expression = as.numeric(cluster_average), 
        percent_expressing = as.numeric(cluster_pct), stringsAsFactors = FALSE)
    cluster_panel <- merge(available_panel, cluster_panel, by = "resolved_feature", all.x = TRUE, sort = FALSE)
    panel_row_index <- panel_row_index + 1L
    panel_rows[[panel_row_index]] <- cluster_panel
}

canonical_panel_long <- do.call(rbind, panel_rows)

canonical_panel_long$expression_for_scaling <- log1p(canonical_panel_long$average_expression)

canonical_panel_long$gene_z <- ave(canonical_panel_long$expression_for_scaling, canonical_panel_long$resolved_feature, 
    FUN = safe_zscore)

canonical_panel_long$gene_z <- pmax(-3, pmin(3, canonical_panel_long$gene_z))

lineage_scores <- data.table::as.data.table(canonical_panel_long)[, .(available_gene_n = data.table::uniqueN(resolved_feature), 
    lineage_score = mean(gene_z, na.rm = TRUE), mean_percent_expressing = mean(percent_expressing, na.rm = TRUE), 
    mean_expression = mean(average_expression, na.rm = TRUE)), by = .(cluster, lineage)]

lineage_scores <- as.data.frame(lineage_scores)

lineage_scores <- lineage_scores[order(as.integer(lineage_scores$cluster), -lineage_scores$lineage_score), , drop = FALSE]

write_csv_utf8(canonical_panel_long, file.path(panel_root, "canonical_panel.csv"))

write_csv_utf8(lineage_scores, file.path(panel_root, "lineage_scores.csv"))

annotation_lineages <- setdiff(names(canonical_panels), "Cycling")

top50_marker_panel_overlap <- merge(all_markers[all_markers$rank <= 50, c("cluster", "gene"), drop = FALSE], available_panel[, 
    c("lineage", "resolved_feature"), drop = FALSE], by.x = "gene", by.y = "resolved_feature", all = FALSE, sort = FALSE)

top50_overlap_counts <- if (nrow(top50_marker_panel_overlap)) {
    as.data.frame(data.table::as.data.table(top50_marker_panel_overlap)[, .(top50_marker_overlap_n = data.table::uniqueN(gene)), 
        by = .(cluster, lineage)])
} else {
    data.frame(cluster = character(), lineage = character(), top50_marker_overlap_n = integer(), stringsAsFactors = FALSE)
}

candidate_rows <- list()

for (cluster_id in cluster_ids) {
    cluster_scores <- lineage_scores[lineage_scores$cluster == cluster_id & lineage_scores$lineage %in% annotation_lineages, 
        , drop = FALSE]
    cluster_scores <- cluster_scores[order(-cluster_scores$lineage_score, -cluster_scores$mean_percent_expressing), 
        , drop = FALSE]
    top_lineage <- cluster_scores$lineage[1]
    second_lineage <- cluster_scores$lineage[2]
    top_score <- cluster_scores$lineage_score[1]
    second_score <- cluster_scores$lineage_score[2]
    score_margin <- top_score - second_score
    overlap_value <- top50_overlap_counts[top50_overlap_counts$cluster == cluster_id & top50_overlap_counts$lineage == 
        top_lineage, "top50_marker_overlap_n"]
    if (!length(overlap_value)) {
        overlap_value <- 0L
    }
    cycling_score <- lineage_scores$lineage_score[lineage_scores$cluster == cluster_id & lineage_scores$lineage == 
        "Cycling"]
    if (!length(cycling_score)) {
        cycling_score <- NA_real_
    }
    review_cluster <- cluster_id %in% c("7", "8", "13", "15")
    support_status <- if (score_margin >= 0.5 && overlap_value >= 2) {
        "HIGH_SUPPORT_CANDIDATE"
    }
    else if (score_margin >= 0.25 || overlap_value >= 2) {
        "MODERATE_SUPPORT_CANDIDATE"
    }
    else {
        "AMBIGUOUS_REVIEW_REQUIRED"
    }
    if (review_cluster) {
        support_status <- paste0(support_status, "_MANDATORY_REVIEW")
    }
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(cluster = cluster_id, candidate_broad_lineage = top_lineage, 
        second_lineage = second_lineage, top_lineage_score = top_score, second_lineage_score = second_score, score_margin = score_margin, 
        top50_marker_overlap_n = as.integer(overlap_value), cycling_score = cycling_score, candidate_support_status = support_status, 
        mandatory_review_cluster = review_cluster, annotation_locked = FALSE, stringsAsFactors = FALSE)
}

candidate_annotation <- do.call(rbind, candidate_rows)

candidate_annotation <- candidate_annotation[order(as.integer(candidate_annotation$cluster)), , drop = FALSE]

write_csv_utf8(top50_overlap_counts, file.path(panel_root, "top50_panel_overlap.csv"))

write_csv_utf8(candidate_annotation, file.path(panel_root, "candidate_annotation.csv"))

marker_summary <- data.table::as.data.table(all_markers)[, .(marker_row_n = .N, FDR05_marker_n = sum(p_val_adj < 
    0.050000000000000003, na.rm = TRUE), logFC05_FDR05_marker_n = sum(p_val_adj < 0.050000000000000003 & avg_log2FC >= 
    0.5, na.rm = TRUE), robust_donor_supported_n = sum(robust_donor_supported, na.rm = TRUE), top_genes = paste(head(gene, 
    10), collapse = "; "), top_robust_genes = paste(head(gene[robust_donor_supported], 10), collapse = "; ")), 
    by = cluster]

marker_summary <- as.data.frame(marker_summary)

marker_summary <- merge(marker_summary, candidate_annotation, by = "cluster", all.x = TRUE, sort = FALSE)

marker_summary <- marker_summary[order(as.integer(marker_summary$cluster)), , drop = FALSE]

write_csv_utf8(marker_summary, file.path(output_root, "02_marker_summary.csv"))

heatmap_gene_rows <- list()

for (cluster_id in cluster_ids) {
    cluster_markers <- all_markers[all_markers$cluster == cluster_id, , drop = FALSE]
    robust_genes <- cluster_markers$gene[cluster_markers$robust_donor_supported]
    selected_genes <- unique(c(head(robust_genes, 5), head(cluster_markers$gene, 5)))
    selected_genes <- head(selected_genes, 5)
    heatmap_gene_rows[[length(heatmap_gene_rows) + 1L]] <- data.frame(source_cluster = cluster_id, gene = selected_genes, 
        stringsAsFactors = FALSE)
}

heatmap_gene_manifest <- unique(do.call(rbind, heatmap_gene_rows))

heatmap_genes <- unique(heatmap_gene_manifest$gene)

heatmap_expression <- get_data_matrix(integrated_object, heatmap_genes)

heatmap_rows <- list()

heatmap_row_index <- 0L

for (cluster_id in cluster_ids) {
    cluster_cells <- which(as.character(metadata$seurat_clusters) == cluster_id)
    average_values <- Matrix::rowMeans(heatmap_expression[, cluster_cells, drop = FALSE])
    heatmap_row_index <- heatmap_row_index + 1L
    heatmap_rows[[heatmap_row_index]] <- data.frame(cluster = cluster_id, gene = rownames(heatmap_expression), 
        average_expression = as.numeric(average_values), stringsAsFactors = FALSE)
}

heatmap_data <- do.call(rbind, heatmap_rows)

heatmap_data$gene_z <- ave(log1p(heatmap_data$average_expression), heatmap_data$gene, FUN = safe_zscore)

heatmap_data$gene_z <- pmax(-3, pmin(3, heatmap_data$gene_z))

write_csv_utf8(heatmap_gene_manifest, file.path(panel_root, "heatmap_gene_manifest.csv"))

write_csv_utf8(heatmap_data, file.path(panel_root, "heatmap_source_data.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

marker_summary_plot <- marker_summary

marker_summary_plot$cluster <- factor(marker_summary_plot$cluster, levels = cluster_ids)

plot_marker_counts <- ggplot2::ggplot(marker_summary_plot, ggplot2::aes(x = cluster, y = FDR05_marker_n, shape = mandatory_review_cluster)) + 
    ggplot2::geom_point(size = 3) + ggplot2::scale_y_continuous(labels = scales::comma) + ggplot2::labs(title = "GSE183837 descriptive positive-marker counts", 
    subtitle = "Wilcoxon cluster-versus-rest; descriptive only", x = "Cluster", y = "FDR < 0.05 marker rows", shape = "Mandatory review") + 
    theme_QC

save_plot_pair(plot_marker_counts, "QC1_marker_counts", width = 8, height = 5)

plot_robust_counts <- ggplot2::ggplot(marker_summary_plot, ggplot2::aes(x = cluster, y = robust_donor_supported_n, 
    shape = mandatory_review_cluster)) + ggplot2::geom_point(size = 3) + ggplot2::scale_y_continuous(labels = scales::comma) + 
    ggplot2::labs(title = "Donor-supported descriptive markers", subtitle = "Requires FDR, effect size, detection and support in at least half of eligible donors", 
        x = "Cluster", y = "Donor-supported marker count", shape = "Mandatory review") + theme_QC

save_plot_pair(plot_robust_counts, "QC2_donor_supported_marker_counts", width = 8, height = 5)

heatmap_data$cluster <- factor(heatmap_data$cluster, levels = cluster_ids)

heatmap_data$gene <- factor(heatmap_data$gene, levels = rev(unique(heatmap_gene_manifest$gene)))

plot_marker_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = cluster, y = gene, fill = gene_z)) + ggplot2::geom_tile() + 
    ggplot2::scale_fill_viridis_c(option = "C") + ggplot2::labs(title = "Top descriptive marker expression by cluster", 
    subtitle = "Row-scaled log-normalized average expression", x = "Cluster", y = "Selected marker", fill = "Row z-score") + 
    theme_QC + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))

save_plot_pair(plot_marker_heatmap, "QC3_top_marker_heatmap", width = 10, height = 14)

lineage_plot_data <- lineage_scores

lineage_plot_data$cluster <- factor(lineage_plot_data$cluster, levels = cluster_ids)

lineage_plot_data$lineage <- factor(lineage_plot_data$lineage, levels = rev(names(canonical_panels)))

plot_lineage <- ggplot2::ggplot(lineage_plot_data, ggplot2::aes(x = cluster, y = lineage, fill = lineage_score)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "C") + ggplot2::labs(title = "Canonical lineage-panel scores", 
    subtitle = "Candidate broad identities require marker-level review", x = "Cluster", y = "Canonical panel", 
    fill = "Mean gene z-score") + theme_QC

save_plot_pair(plot_lineage, "QC4_lineage_panel_heatmap", width = 10, height = 7)

candidate_plot <- candidate_annotation

candidate_plot$cluster <- factor(candidate_plot$cluster, levels = cluster_ids)

plot_margin <- ggplot2::ggplot(candidate_plot, ggplot2::aes(x = cluster, y = score_margin, shape = mandatory_review_cluster)) + 
    ggplot2::geom_hline(yintercept = c(0.25, 0.5), linetype = "dashed") + ggplot2::geom_point(size = 3) + ggplot2::labs(title = "Candidate broad-lineage score margins", 
    subtitle = "Low margins indicate ambiguous lineage nomination", x = "Cluster", y = "Top minus second lineage score", 
    shape = "Mandatory review") + theme_QC

save_plot_pair(plot_margin, "QC5_candidate_lineage_margin", width = 8, height = 5)

review_support <- donor_support_summary_all[donor_support_summary_all$cluster %in% c("7", "8", "13", "15"), , drop = FALSE]

review_support <- merge(review_support, all_markers[, c("cluster", "gene", "rank", "avg_log2FC", "p_val_adj"), 
    drop = FALSE], by = c("cluster", "gene"), all.x = TRUE, sort = FALSE)

review_support <- review_support[review_support$rank <= 30, , drop = FALSE]

review_support$cluster <- factor(review_support$cluster, levels = c("7", "8", "13", "15"))

write_csv_utf8(review_support, file.path(panel_root, "review_cluster_marker_support.csv"))

plot_review_support <- ggplot2::ggplot(review_support, ggplot2::aes(x = supporting_donor_fraction, y = avg_log2FC, 
    shape = cluster)) + ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed") + ggplot2::geom_point(alpha = 0.75) + 
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + ggplot2::labs(title = "Marker donor support in mandatory-review clusters", 
    subtitle = "Cluster 13 requires special caution because of severe donor dominance", x = "Supporting fraction of eligible donors", 
    y = "Average log2 fold change", shape = "Cluster") + theme_QC

save_plot_pair(plot_review_support, "QC6_review_cluster_donor_support", width = 8, height = 6)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

ambiguous_clusters <- candidate_annotation$cluster[grepl("^AMBIGUOUS", candidate_annotation$candidate_support_status)]

review_clusters <- c("7", "8", "13", "15")

gate_checks <- data.frame(check = c("upstream_detailed_QC_authorizes_markers", "input_checkpoint_has_47161_cells", 
    "all_16_cluster_marker_jobs_complete", "every_cluster_has_positive_marker_rows", "combined_marker_table_has_all_16_clusters", 
    "all_top100_table_has_all_16_clusters", "donor_support_completed_for_all_16_clusters", "canonical_panel_coverage_at_least_75_percent", 
    "lineage_scores_present_for_all_clusters", "candidate_annotation_has_16_clusters", "mandatory_review_clusters_7_8_13_15_present", 
    "no_cluster_or_cell_was_removed", "candidate_annotation_review_can_start", "free_disk_at_least_8GB_or_unknown"), 
    pass = c(upstream_gate$gate_status[[1]] %in% c("PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED", "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED"), 
        ncol(integrated_object) == 47161L, nrow(cluster_status_all) == 16L && all(cluster_status_all$status == 
            "COMPLETE"), all(cluster_status_all$marker_row_n > 0), length(unique(all_markers$cluster)) == 16L, 
        length(unique(all_top100$cluster)) == 16L, length(unique(donor_support_summary_all$cluster)) == 16L, nrow(available_panel) >= 
            0.75 * nrow(panel_requested), length(unique(lineage_scores$cluster)) == 16L, nrow(candidate_annotation) == 
            16L, all(review_clusters %in% candidate_annotation$cluster), ncol(integrated_object) == 47161L && length(unique(metadata$seurat_clusters)) == 
            16L, TRUE, is.na(suppressWarnings(as.numeric(system2("powershell", args = c("-NoProfile", "-Command", 
            shQuote("(Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='C:'\").FreeSpace / 1GB")), stdout = TRUE, 
            stderr = FALSE)[1]))) || suppressWarnings(as.numeric(system2("powershell", args = c("-NoProfile", "-Command", 
            shQuote("(Get-CimInstance Win32_LogicalDisk -Filter \"DeviceID='C:'\").FreeSpace / 1GB")), stdout = TRUE, 
            stderr = FALSE)[1])) >= 8), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED"
} else {
    "HOLD_DESCRIPTIVE_MARKERS_INCOMPLETE"
}

review_flags <- c("CLUSTER_13_SEVERE_SINGLE_DONOR_DOMINANCE", "CLUSTERS_7_8_MODERATE_DONOR_ENRICHMENT", "CLUSTER_15_RARE_POPULATION_REVIEW", 
    if (length(ambiguous_clusters)) {
        paste0("AMBIGUOUS_LINEAGE_CLUSTERS_", paste(ambiguous_clusters, collapse = "_"))
    } else {
        character()
    })

marker_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, marker_version = marker_version, 
    gate_status = gate_status, cell_n = ncol(integrated_object), cluster_n = 16L, total_marker_row_n = nrow(all_markers), 
    FDR05_marker_row_n = sum(all_markers$p_val_adj < 0.050000000000000003, na.rm = TRUE), robust_donor_supported_marker_n = sum(all_markers$robust_donor_supported, 
        na.rm = TRUE), canonical_panel_available_n = nrow(available_panel), canonical_panel_requested_n = nrow(panel_requested), 
    ambiguous_candidate_clusters = collapse_unique(ambiguous_clusters), review_flags = collapse_unique(review_flags), 
    candidate_annotation_review_authorized = all_gate_checks_pass, broad_annotation_locked = FALSE, disease_pseudobulk_testing_authorized = FALSE, 
    result_lock_authorized = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(marker_gate, file.path(run_qc_root, "G183_MARKER_GATE.csv"))

card_lines <- c("# STEP04E1L GSE183837 Descriptive Marker Card", "", paste0("- Gate: ", gate_status), paste0("- Cells retained: ", 
    format(ncol(integrated_object), big.mark = ",")), "- Primary clusters: 16", paste0("- Total marker rows: ", 
    format(nrow(all_markers), big.mark = ",")), paste0("- FDR < 0.05 marker rows: ", format(sum(all_markers$p_val_adj < 
    0.050000000000000003, na.rm = TRUE), big.mark = ",")), paste0("- Donor-supported marker rows: ", format(sum(all_markers$robust_donor_supported, 
    na.rm = TRUE), big.mark = ",")), paste0("- Ambiguous candidate clusters: ", if (length(ambiguous_clusters)) {
    collapse_unique(ambiguous_clusters)
} else {
    "NONE"
}), paste0("- Candidate annotation review authorized: ", all_gate_checks_pass), "- Broad annotation locked: FALSE", 
    "- Disease pseudobulk testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "G183_MARKER_CARD.md"))

interpretation_lines <- c("# GSE183837 descriptive marker interpretation", "", paste0("- Run ID: `", run_id, "`"), 
    paste0("- Gate: **", gate_status, "**"), "", "## Interpretation constraints", "", "- Marker tests are cluster-versus-rest descriptions and are not Control-versus-RIF comparisons.", 
    "- Candidate broad lineages are algorithmic nominations requiring manual review.", "- Cluster 13 cannot be described as an RIF-specific population because it is dominated by one RIF donor.", 
    "- Clusters 7 and 8 require donor-robust marker confirmation.", "- Cluster 15 remains a retained rare candidate population.", 
    "", "## Authorized next operation", "", "- Review top markers, canonical panels and donor support.", "- Construct and QC a broad-lineage annotation map.", 
    "- Lock broad annotation only after marker and visual review.", "", "## Not authorized yet", "", "- No disease pseudobulk analysis.", 
    "- No fine-subtype lock.", "- No result lock.")

write_text(interpretation_lines, file.path(run_qc_root, "G183_MARKER_INTERPRETATION.md"))

parameter_manifest <- list(script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_CANONICAL_PANEL_UNIQUE_COLLAPSED_84_ROWS_TO_81", 
    panel_resolution_policy = "PRESERVE_CROSS_LINEAGE_DUPLICATE_GENES_AND_INPUT_ROW_ORDER", dataset_id = dataset_id, 
    integration_version = integration_version, marker_version = marker_version, upstream_gate = upstream_gate$gate_status[[1]], 
    marker_parameters = marker_parameters, marker_role = "DESCRIPTIVE_CLUSTER_VERSUS_REST", donor_support_rule = "PCT_DELTA_AT_LEAST_0_10_IN_AT_LEAST_HALF_OF_ELIGIBLE_DONORS", 
    canonical_panel_role = "CANDIDATE_BROAD_LINEAGE_NOMINATION_ONLY", mandatory_review_clusters = c("7", "8", "13", 
        "15"), broad_annotation_locked = FALSE, disease_cell_level_testing_authorized = FALSE, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_CANONICAL_PANEL_UNIQUE_COLLAPSED_84_ROWS_TO_81", 
    dataset_id = dataset_id, integration_version = integration_version, marker_version = marker_version, started_at = as.character(run_start), 
    completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 
        2), gate_status = gate_status, cell_n = ncol(integrated_object), cluster_n = 16L, total_marker_row_n = nrow(all_markers), 
    FDR05_marker_row_n = sum(all_markers$p_val_adj < 0.050000000000000003, na.rm = TRUE), robust_donor_supported_marker_n = sum(all_markers$robust_donor_supported, 
        na.rm = TRUE), candidate_annotation_review_authorized = all_gate_checks_pass, broad_annotation_locked = FALSE, 
    result_lock_authorized = FALSE, output_root = normalize_slashes(output_root), qc_root = normalize_slashes(run_qc_root), 
    stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

persistent_files <- unique(c(list.files(output_root, recursive = TRUE, full.names = TRUE)))

persistent_files <- persistent_files[file.exists(persistent_files) & !dir.exists(persistent_files)]

persistent_manifest <- data.frame(path = normalize_slashes(persistent_files), size_bytes = file.info(persistent_files)$size, 
    md5 = unname(tools::md5sum(persistent_files)), stringsAsFactors = FALSE)

write_csv_utf8(persistent_manifest, file.path(manifest_root, "persistent_manifest.csv"))

small_files <- unique(c(file.path(output_root, "01_marker_status.csv"), file.path(output_root, "02_marker_summary.csv"), 
    list.files(top_root, pattern = "top100|donor_summary|donor_support_summary", full.names = TRUE), list.files(panel_root, 
        recursive = TRUE, full.names = TRUE), list.files(figure_output_root, recursive = TRUE, full.names = TRUE), 
    list.files(run_qc_root, recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, 
    session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2]

artifact_manifest <- data.frame(path = normalize_slashes(small_files), size_bytes = file.info(small_files)$size, 
    md5 = unname(tools::md5sum(small_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "total_marker_rows=", nrow(all_markers), "robust_donor_supported=", 
    sum(all_markers$robust_donor_supported, na.rm = TRUE))

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

cat("STEP04E1L GSE183837 descriptive markers completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Cells: ", format(ncol(integrated_object), big.mark = ","), "\n", sep = "")

cat("Clusters completed: 16\n")

cat("Total marker rows: ", format(nrow(all_markers), big.mark = ","), "\n", sep = "")

cat("FDR < 0.05 marker rows: ", format(sum(all_markers$p_val_adj < 0.050000000000000003, na.rm = TRUE), big.mark = ","), 
    "\n", sep = "")

cat("Donor-supported marker rows: ", format(sum(all_markers$robust_donor_supported, na.rm = TRUE), big.mark = ","), 
    "\n", sep = "")

cat("Ambiguous candidate clusters: ", if (length(ambiguous_clusters)) {
    collapse_unique(ambiguous_clusters)
} else {
    "NONE"
}, "\n", sep = "")

cat("Candidate annotation review authorized: ", all_gate_checks_pass, "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nMarker outputs are complete in the project folders.\n", sep = "")
}

