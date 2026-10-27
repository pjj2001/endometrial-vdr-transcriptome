#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C04: Calculate cluster marker expression

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260805)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "identify_gse250130_markers.R"

script_version <- "STEP04E1D_GSE250130_CLUSTER_MARKERS_v1.1_MEMORY_SAFE"

input_version <- "SCI_v1_1_memory_safe"

marker_version <- "MARKERS_v1_1_memory_safe"

cluster_column <- "RNA_snn_res.0.4"

run_start <- Sys.time()

run_id <- paste0("STEP04E1D_", dataset_id, "_v1_1_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

numeric_cluster_order <- function(x) {
    nx <- suppressWarnings(as.integer(as.character(x)))
    order(nx, as.character(x), na.last = TRUE)
}

save_plot_pair <- function(plot_object, stem, width, height) {
    png_path <- file.path(figure_root, paste0(stem, ".png"))
    pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))
    ggplot2::ggsave(filename = png_path, plot = plot_object, width = width, height = height, units = "in", dpi = 300, 
        bg = "white")
    ggplot2::ggsave(filename = pdf_path, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    c(png_path, pdf_path)
}

balanced_sample_cells <- function(cells, sample_by_cell, max_total) {
    cells <- unique(as.character(cells))
    cells <- cells[cells %in% names(sample_by_cell)]
    if (!length(cells)) 
        return(character())
    if (length(cells) <= max_total) 
        return(cells)
    sample_values <- as.character(sample_by_cell[cells])
    groups <- split(cells, sample_values, drop = TRUE)
    group_n <- length(groups)
    initial_quota <- max(1L, floor(max_total/group_n))
    chosen <- unlist(lapply(groups, function(one_group) {
        if (length(one_group) <= initial_quota) {
            one_group
        }
        else {
            sample(one_group, size = initial_quota, replace = FALSE)
        }
    }), use.names = FALSE)
    remaining_n <- max_total - length(chosen)
    if (remaining_n > 0L) {
        remaining_pool <- setdiff(cells, chosen)
        if (length(remaining_pool)) {
            extra <- sample(remaining_pool, size = min(remaining_n, length(remaining_pool)), replace = FALSE)
            chosen <- c(chosen, extra)
        }
    }
    if (length(chosen) > max_total) {
        chosen <- sample(chosen, size = max_total, replace = FALSE)
    }
    unique(chosen)
}

calculate_sample_breadth <- function(target_matrix, target_samples, detection_threshold = 0.10000000000000001) {
    target_samples <- factor(as.character(target_samples))
    sample_n <- nlevels(target_samples)
    if (sample_n == 0L) {
        return(list(breadth_fraction = rep(NA_real_, nrow(target_matrix)), sample_n = 0L))
    }
    membership <- Matrix::sparseMatrix(i = seq_along(target_samples), j = as.integer(target_samples), x = 1, dims = c(length(target_samples), 
        sample_n))
    detected <- target_matrix > 0
    detection_counts <- detected %*% membership
    cells_per_sample <- tabulate(as.integer(target_samples), nbins = sample_n)
    detection_fraction <- as.matrix(detection_counts)
    detection_fraction <- sweep(detection_fraction, MARGIN = 2, STATS = cells_per_sample, FUN = "/")
    breadth_fraction <- rowMeans(detection_fraction >= detection_threshold)
    list(breadth_fraction = as.numeric(breadth_fraction), sample_n = sample_n)
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

checkpoint_root <- file.path(results_root, "S4E1_sc", input_version, "checkpoints", dataset_id)

checkpoint_path <- file.path(checkpoint_root, paste0(dataset_id, "_integrated_checkpoint.rds"))

checkpoint_status_path <- file.path(checkpoint_root, paste0(dataset_id, "_checkpoint_status.csv"))

marker_result_root <- file.path(results_root, "S4E1_sc", input_version, "markers", dataset_id, marker_version)

per_cluster_root <- file.path(marker_result_root, "per_cluster")

table_root <- file.path(marker_result_root, "tables")

source_data_root <- file.path(marker_result_root, "source_data")

run_qc_root <- file.path(qc_base_root, paste0(run_id, "_cluster_markers_memory_safe"))

figure_root <- file.path(run_qc_root, "figures")

for (p in c(execution_root, parameter_root, environment_root, manuscript_text_root, marker_result_root, per_cluster_root, 
    table_root, source_data_root, run_qc_root, figure_root)) {
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

if (!file.exists(checkpoint_path)) {
    stop("Integrated checkpoint not found: ", checkpoint_path)
}

checkpoint_status <- safe_read_csv(checkpoint_status_path)

if (is.null(checkpoint_status) || nrow(checkpoint_status) != 1L || checkpoint_status$checkpoint_status[[1]] != 
    "COMPLETE" || checkpoint_status$join_policy[[1]] != "DATA_ONLY_COUNTS_REMAIN_SPLIT") {
    stop("Integrated checkpoint status is absent or invalid.")
}

qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

qc_dirs <- qc_dirs[grepl("^STEP04E1C_GSE250130_[0-9]{8}_[0-9]{6}_integration_QC$", basename(qc_dirs), perl = TRUE)]

if (!length(qc_dirs)) {
    stop("No STEP04E1C integration-QC folder found.")
}

latest_qc_dir <- qc_dirs[which.max(file.info(qc_dirs)$mtime)]

integration_gate_path <- file.path(latest_qc_dir, "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv")

integration_gate <- safe_read_csv(integration_gate_path)

if (is.null(integration_gate) || nrow(integration_gate) != 1L || !integration_gate$gate_status[[1]] %in% c("PASS", 
    "PASS_WITH_REVIEW_FLAGS") || !isTRUE(integration_gate$marker_analysis_authorized[[1]])) {
    stop("STEP04E1C Gate does not authorize marker analysis.")
}

review_clusters <- character()

if ("review_clusters" %in% colnames(integration_gate) && !is.na(integration_gate$review_clusters[[1]]) && nzchar(integration_gate$review_clusters[[1]])) {
    review_clusters <- trimws(unlist(strsplit(integration_gate$review_clusters[[1]], ";", fixed = TRUE)))
}

free_memory_gb <- get_free_memory_gb()

free_disk_gb <- get_free_disk_gb("C")

resource_snapshot <- data.frame(metric = c("free_physical_memory_GB_before_load", "free_C_drive_GB_before_run", 
    "checkpoint_size_GB", "checkpoint_md5", "integration_gate_status", "review_clusters", "failed_predecessor_marker_version"), 
    value = c(free_memory_gb, free_disk_gb, round(file.info(checkpoint_path)$size/1024^3, 3), unname(tools::md5sum(checkpoint_path)), 
        integration_gate$gate_status[[1]], collapse_unique(review_clusters), "MARKERS_v1_0_FINDMARKERS_MEMORY_FAILURE"), 
    stringsAsFactors = FALSE)

write_csv_utf8(resource_snapshot, file.path(run_qc_root, "01_RESOURCE_AND_INPUT_SNAPSHOT.csv"))

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
    stop("Free physical memory is below 18 GB. ", "Restart RStudio and rerun in a clean session.")
}

if (is.finite(free_disk_gb) && free_disk_gb < 10) {
    stop("Free C: drive space is below 10 GB.")
}

final_marker_path <- file.path(table_root, paste0(dataset_id, "_all_cluster_markers_descriptive.csv"))

final_gate_path <- file.path(run_qc_root, "STEP04E1D_GSE250130_MARKER_GATE.csv")

parameter_manifest <- list(script_name = script_name, script_version = script_version, input_version = input_version, 
    marker_version = marker_version, dataset_id = dataset_id, cluster_column = cluster_column, marker_method = "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING", 
    target_max_cells = 1500L, background_max_cells_per_other_cluster = 150L, target_min_detection_fraction = 0.10000000000000001, 
    minimum_detection_fraction_difference = 0.050000000000000003, minimum_mean_log_expression_difference = 0.20000000000000001, 
    minimum_target_sample_breadth = 0.25, fallback_target_min_detection_fraction = 0.050000000000000003, fallback_detection_fraction_difference = 0.02, 
    fallback_mean_log_expression_difference = 0.10000000000000001, primary_minimum_marker_n = 25L, sample_breadth_detection_threshold = 0.10000000000000001, 
    canonical_panel_max_cells_per_cluster = 500L, p_values_reported = FALSE, inferential_role = "NONE_ANNOTATION_SUPPORT_ONLY", 
    random_seed = 20260805L, integration_gate = integration_gate$gate_status[[1]], review_clusters = review_clusters, 
    output_status = "CANDIDATE_MARKERS_FOR_ANNOTATION")

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

log_event("Starting", run_id, "checkpoint=", checkpoint_path, "integration_gate=", integration_gate$gate_status[[1]], 
    "review_clusters=", collapse_unique(review_clusters))

load_start <- Sys.time()

object <- readRDS(checkpoint_path)

log_event("Loaded checkpoint", "elapsed_min=", round(as.numeric(difftime(Sys.time(), load_start, units = "mins")), 
    2))

if (!inherits(object, "Seurat")) {
    stop("Checkpoint is not a Seurat object.")
}

if (ncol(object) != 236636L || nrow(object) != 36601L) {
    stop("Unexpected object dimensions: ", nrow(object), " features x ", ncol(object), " cells")
}

if (!cluster_column %in% colnames(object@meta.data)) {
    stop("Cluster column missing: ", cluster_column)
}

if (!"sample_accession" %in% colnames(object@meta.data)) {
    stop("sample_accession is missing from object metadata.")
}

SeuratObject::DefaultAssay(object) <- "RNA"

data_layers <- SeuratObject::Layers(object[["RNA"]], search = "^data$")

if (!identical(data_layers, "data")) {
    stop("Expected exactly one normalized data layer named 'data'. Found: ", paste(data_layers, collapse = "; "))
}

cell_metadata <- data.frame(cell_barcode = rownames(object@meta.data), sample_accession = as.character(object$sample_accession), 
    cluster = as.character(object@meta.data[[cluster_column]]), stringsAsFactors = FALSE)

rownames(cell_metadata) <- cell_metadata$cell_barcode

sample_by_cell <- setNames(cell_metadata$sample_accession, cell_metadata$cell_barcode)

cells_by_cluster <- split(cell_metadata$cell_barcode, cell_metadata$cluster)

cluster_levels <- names(cells_by_cluster)

cluster_levels <- cluster_levels[numeric_cluster_order(cluster_levels)]

cluster_counts <- data.frame(cluster = cluster_levels, cell_n = vapply(cells_by_cluster[cluster_levels], length, 
    integer(1)), sample_n = vapply(cells_by_cluster[cluster_levels], function(cells) {
    length(unique(sample_by_cell[cells]))
}, integer(1)), review_cluster = (cluster_levels %in% review_clusters), stringsAsFactors = FALSE)

write_csv_utf8(cluster_counts, file.path(table_root, "01_CLUSTER_CELL_COUNTS.csv"))

if (length(cluster_levels) != 22L || sum(cluster_counts$cell_n) != 236636L) {
    stop("Cluster identity validation failed.")
}

layer_manifest <- data.frame(layer = SeuratObject::Layers(object[["RNA"]]), layer_type = ifelse(grepl("^counts", 
    SeuratObject::Layers(object[["RNA"]])), "counts", ifelse(grepl("^data", SeuratObject::Layers(object[["RNA"]])), 
    "data", "other")), stringsAsFactors = FALSE)

write_csv_utf8(layer_manifest, file.path(run_qc_root, "02_LAYER_MANIFEST.csv"))

marker_status_rows <- list()

for (cluster_id in cluster_levels) {
    cluster_file <- file.path(per_cluster_root, paste0(dataset_id, "_cluster_", cluster_id, "_descriptive_markers.csv"))
    status_file <- file.path(per_cluster_root, paste0(dataset_id, "_cluster_", cluster_id, "_marker_status.csv"))
    selected_cells_file <- file.path(per_cluster_root, paste0(dataset_id, "_cluster_", cluster_id, "_selected_cells.csv"))
    existing_status <- safe_read_csv(status_file)
    if (file.exists(cluster_file) && !is.null(existing_status) && nrow(existing_status) == 1L && existing_status$status[[1]] == 
        "COMPLETE") {
        log_event("SKIP completed cluster", cluster_id)
        marker_status_rows[[length(marker_status_rows) + 1L]] <- existing_status
        next
    }
    cluster_start <- Sys.time()
    cluster_error <- ""
    target_pool <- cells_by_cluster[[cluster_id]]
    target_cells <- balanced_sample_cells(cells = target_pool, sample_by_cell = sample_by_cell, max_total = 1500L)
    other_clusters <- setdiff(cluster_levels, cluster_id)
    background_cells <- unlist(lapply(other_clusters, function(other_cluster) {
        balanced_sample_cells(cells = cells_by_cluster[[other_cluster]], sample_by_cell = sample_by_cell, max_total = 150L)
    }), use.names = FALSE)
    background_cells <- unique(background_cells)
    background_cells <- setdiff(background_cells, target_cells)
    selected_cells <- c(target_cells, background_cells)
    selected_cell_manifest <- data.frame(cell_barcode = selected_cells, role = c(rep("TARGET", length(target_cells)), 
        rep("CLUSTER_BALANCED_BACKGROUND", length(background_cells))), sample_accession = sample_by_cell[selected_cells], 
        cluster = cell_metadata[selected_cells, "cluster"], target_cluster = cluster_id, stringsAsFactors = FALSE)
    write_csv_utf8(selected_cell_manifest, selected_cells_file)
    log_event("START descriptive marker cluster", cluster_id, "target_pool_n=", length(target_pool), "target_selected_n=", 
        length(target_cells), "background_selected_n=", length(background_cells))
    marker_result <- tryCatch({
        selected_matrix <- SeuratObject::LayerData(object = object, assay = "RNA", layer = "data", cells = selected_cells)
        if (!inherits(selected_matrix, "sparseMatrix")) {
            selected_matrix <- Matrix::Matrix(selected_matrix, sparse = TRUE)
        }
        target_index <- match(target_cells, colnames(selected_matrix))
        background_index <- match(background_cells, colnames(selected_matrix))
        if (anyNA(target_index) || anyNA(background_index)) {
            stop("Selected cells were not all returned by LayerData.")
        }
        target_matrix <- selected_matrix[, target_index, drop = FALSE]
        background_matrix <- selected_matrix[, background_index, drop = FALSE]
        mean_target <- Matrix::rowMeans(target_matrix)
        mean_background <- Matrix::rowMeans(background_matrix)
        pct_target <- Matrix::rowMeans(target_matrix > 0)
        pct_background <- Matrix::rowMeans(background_matrix > 0)
        mean_log_expression_difference <- (mean_target - mean_background)
        detection_fraction_difference <- (pct_target - pct_background)
        log2_detection_ratio <- log2((pct_target + 0.01)/(pct_background + 0.01))
        breadth <- calculate_sample_breadth(target_matrix = target_matrix, target_samples = sample_by_cell[target_cells], 
            detection_threshold = 0.10000000000000001)
        marker_score <- (pmax(mean_log_expression_difference, 0) * pmax(detection_fraction_difference, 0) * (0.5 + 
            breadth$breadth_fraction))
        all_gene_scores <- data.frame(cluster = cluster_id, gene = rownames(selected_matrix), mean_log_expression_target = as.numeric(mean_target), 
            mean_log_expression_background = as.numeric(mean_background), mean_log_expression_difference = as.numeric(mean_log_expression_difference), 
            pct_target = as.numeric(pct_target), pct_background = as.numeric(pct_background), detection_fraction_difference = as.numeric(detection_fraction_difference), 
            log2_detection_ratio = as.numeric(log2_detection_ratio), target_sample_breadth_fraction = as.numeric(breadth$breadth_fraction), 
            target_sample_n = as.integer(breadth$sample_n), target_pool_cell_n = as.integer(length(target_pool)), 
            target_selected_cell_n = as.integer(length(target_cells)), background_selected_cell_n = as.integer(length(background_cells)), 
            marker_score = as.numeric(marker_score), review_cluster = (cluster_id %in% review_clusters), stringsAsFactors = FALSE)
        primary_keep <- (all_gene_scores$pct_target >= 0.10000000000000001 & all_gene_scores$detection_fraction_difference >= 
            0.050000000000000003 & all_gene_scores$mean_log_expression_difference >= 0.20000000000000001 & all_gene_scores$target_sample_breadth_fraction >= 
            0.25)
        relaxed_keep <- (all_gene_scores$pct_target >= 0.050000000000000003 & all_gene_scores$detection_fraction_difference >= 
            0.02 & all_gene_scores$mean_log_expression_difference >= 0.10000000000000001 & all_gene_scores$target_sample_breadth_fraction >= 
            0.20000000000000001)
        all_gene_scores$selection_tier <- ifelse(primary_keep, "PRIMARY", ifelse(relaxed_keep, "RELAXED", "NOT_SELECTED"))
        primary_markers <- all_gene_scores[primary_keep, , drop = FALSE]
        if (nrow(primary_markers) >= 25L) {
            selected_markers <- primary_markers
        }
        else {
            selected_markers <- all_gene_scores[relaxed_keep, , drop = FALSE]
        }
        selected_markers <- selected_markers[order(-selected_markers$marker_score, -selected_markers$mean_log_expression_difference, 
            -selected_markers$detection_fraction_difference, -selected_markers$pct_target, selected_markers$gene), 
            , drop = FALSE]
        selected_markers$marker_role <- "CLUSTER_ANNOTATION_SUPPORT_ONLY"
        selected_markers$inferential_role <- "NONE_NO_P_VALUES_DONOR_INFERENCE_RESERVED_FOR_PSEUDOBULK"
        rm(selected_matrix, target_matrix, background_matrix, all_gene_scores, primary_markers)
        gc()
        selected_markers
    }, error = function(e) {
        cluster_error <<- conditionMessage(e)
        NULL
    })
    elapsed_minutes <- round(as.numeric(difftime(Sys.time(), cluster_start, units = "mins")), 2)
    if (is.null(marker_result) || !nrow(marker_result)) {
        marker_status <- data.frame(cluster = cluster_id, status = ifelse(nzchar(cluster_error), "FAILED", "NO_MARKERS_RETURNED"), 
            marker_n = 0L, primary_marker_n = 0L, relaxed_marker_n = 0L, target_pool_cell_n = length(target_pool), 
            target_selected_cell_n = length(target_cells), background_selected_cell_n = length(background_cells), 
            target_sample_n = length(unique(sample_by_cell[target_cells])), elapsed_minutes = elapsed_minutes, 
            error_message = cluster_error, output_file = normalize_slashes(cluster_file), stringsAsFactors = FALSE)
        write_csv_utf8(marker_status, status_file)
        marker_status_rows[[length(marker_status_rows) + 1L]] <- marker_status
        log_event(marker_status$status[[1]], "cluster=", cluster_id, "error=", cluster_error)
        rm(marker_result)
        gc()
        next
    }
    write_csv_utf8(marker_result, cluster_file)
    marker_status <- data.frame(cluster = cluster_id, status = "COMPLETE", marker_n = nrow(marker_result), primary_marker_n = sum(marker_result$selection_tier == 
        "PRIMARY"), relaxed_marker_n = sum(marker_result$selection_tier == "RELAXED"), target_pool_cell_n = length(target_pool), 
        target_selected_cell_n = length(target_cells), background_selected_cell_n = length(background_cells), target_sample_n = length(unique(sample_by_cell[target_cells])), 
        elapsed_minutes = elapsed_minutes, error_message = "", output_file = normalize_slashes(cluster_file), stringsAsFactors = FALSE)
    write_csv_utf8(marker_status, status_file)
    marker_status_rows[[length(marker_status_rows) + 1L]] <- marker_status
    log_event("DONE descriptive marker cluster", cluster_id, "marker_n=", nrow(marker_result), "primary_n=", sum(marker_result$selection_tier == 
        "PRIMARY"), "elapsed_min=", elapsed_minutes)
    rm(marker_result)
    gc()
}

marker_status_table <- do.call(rbind, marker_status_rows)

marker_status_table <- marker_status_table[numeric_cluster_order(marker_status_table$cluster), , drop = FALSE]

write_csv_utf8(marker_status_table, file.path(table_root, "02_CLUSTER_MARKER_RUN_STATUS.csv"))

failed_clusters <- marker_status_table$cluster[marker_status_table$status != "COMPLETE"]

completed_cluster_files <- marker_status_table$output_file[marker_status_table$status == "COMPLETE"]

completed_cluster_files <- completed_cluster_files[file.exists(completed_cluster_files)]

marker_tables <- lapply(completed_cluster_files, safe_read_csv)

marker_tables <- marker_tables[vapply(marker_tables, function(x) !is.null(x), logical(1))]

all_markers <- if (length(marker_tables)) {
    data.table::rbindlist(marker_tables, use.names = TRUE, fill = TRUE)
} else {
    data.table::data.table()
}

if (nrow(all_markers)) {
    all_markers$cluster <- as.character(all_markers$cluster)
    data.table::setorderv(all_markers, cols = c("cluster", "marker_score", "mean_log_expression_difference", "detection_fraction_difference"), 
        order = c(1L, -1L, -1L, -1L), na.last = TRUE)
    write_csv_utf8(as.data.frame(all_markers), final_marker_path)
    top25 <- all_markers[, head(.SD, 25L), by = cluster]
    write_csv_utf8(as.data.frame(top25), file.path(table_root, paste0(dataset_id, "_top25_descriptive_markers_by_cluster.csv")))
    top10_summary <- top25[, .(top_genes = paste(head(gene, 10L), collapse = "; "), median_marker_score = stats::median(marker_score, 
        na.rm = TRUE), minimum_target_sample_breadth = min(target_sample_breadth_fraction, na.rm = TRUE)), by = cluster]
    write_csv_utf8(as.data.frame(top10_summary), file.path(table_root, paste0(dataset_id, "_top10_marker_summary.csv")))
} else {
    top25 <- data.table::data.table()
    top10_summary <- data.table::data.table()
}

if (nrow(top25)) {
    quality_table <- top25[, {
        genes <- as.character(gene)
        top20_genes <- head(genes, 20L)
        mitochondrial <- grepl("^MT-", top20_genes, ignore.case = TRUE)
        ribosomal <- grepl("^RP[SL][0-9A-Z]", top20_genes, ignore.case = TRUE)
        stress <- top20_genes %in% c("FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "DDIT3", "HSPA1A", "HSPA1B", 
            "HSP90AA1", "NFKBIA", "IER2", "IER3")
        cycling <- top20_genes %in% c("MKI67", "TOP2A", "STMN1", "TUBB", "TYMS", "UBE2C", "CENPF", "CDK1", "PCNA")
        list(top20_gene_n = length(top20_genes), mitochondrial_fraction = mean(mitochondrial), ribosomal_fraction = mean(ribosomal), 
            stress_fraction = mean(stress), cycling_fraction = mean(cycling), minimum_sample_breadth = min(target_sample_breadth_fraction, 
                na.rm = TRUE), top20_genes = paste(top20_genes, collapse = "; "))
    }, by = cluster]
    quality_table$review_cluster <- (quality_table$cluster %in% review_clusters)
    quality_table$technical_marker_flag <- (quality_table$mitochondrial_fraction >= 0.29999999999999999 | quality_table$ribosomal_fraction >= 
        0.5 | quality_table$stress_fraction >= 0.29999999999999999)
    quality_table$cell_cycle_dominant_flag <- (quality_table$cycling_fraction >= 0.29999999999999999)
    write_csv_utf8(as.data.frame(quality_table), file.path(table_root, "03_CLUSTER_MARKER_QUALITY_AUDIT.csv"))
} else {
    quality_table <- data.frame()
}

marker_panels <- list(epithelial_general = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "MSLN", "WFDC2"), 
    ciliated_epithelial = c("FOXJ1", "PIFO", "TPPP3", "CAPS", "C20orf85", "CDC20B"), secretory_epithelial = c("PAEP", 
        "GPX3", "SCGB2A1", "DEFB1", "SLPI"), stromal_fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", 
        "COL6A1", "PDGFRA", "CFD"), decidual_stromal = c("PRL", "IGFBP1", "FOXO1", "HAND2", "IL15", "WNT4", "LEFTY2", 
        "PROK1"), endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "PLVAP"), lymphatic_endothelial = c("PROX1", 
        "PDPN", "LYVE1", "CCL21", "FLT4"), perivascular_smooth_muscle = c("RGS5", "MCAM", "CSPG4", "ACTA2", "TAGLN", 
        "MYH11", "DES", "NOTCH3"), immune_general = c("PTPRC", "TYROBP", "LST1", "FCER1G"), NK_cell = c("NKG7", 
        "GNLY", "KLRD1", "PRF1", "GZMB", "CCL5", "XCL1", "XCL2"), T_cell = c("CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", 
        "IL7R", "LTB"), B_cell = c("MS4A1", "CD79A", "CD74", "CD37", "CD22", "CD19"), plasma_cell = c("MZB1", "JCHAIN", 
        "SDC1", "IGHG1", "IGKC", "DERL3"), myeloid_macrophage = c("LST1", "FCER1G", "TYROBP", "CTSS", "LILRB1", 
        "CD68", "FCGR3A", "C1QA", "C1QB", "C1QC", "APOE", "CD163", "MRC1"), dendritic_cell = c("FCER1A", "CD1C", 
        "CLEC10A", "CLEC9A", "IRF8", "LILRA4"), mast_cell = c("TPSAB1", "TPSB2", "KIT", "HDC", "CPA3", "MS4A2"), 
    erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA"), cycling = c("MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C", 
        "CENPF"))

marker_panel_manifest <- do.call(rbind, lapply(names(marker_panels), function(panel_name) {
    data.frame(panel = panel_name, gene = marker_panels[[panel_name]], panel_role = ifelse(panel_name %in% c("immune_general", 
        "cycling"), "QC_OR_GENERAL_SUPPORT", "CANDIDATE_LINEAGE_SUPPORT"), evidence_basis = "CANONICAL_HUMAN_ENDOMETRIAL_BROAD_LINEAGE_MARKERS", 
        stringsAsFactors = FALSE)
}))

marker_panel_manifest$present_in_object <- (marker_panel_manifest$gene %in% rownames(object))

write_csv_utf8(marker_panel_manifest, file.path(table_root, "04_CANONICAL_MARKER_PANEL_MANIFEST.csv"))

panel_genes_present <- unique(marker_panel_manifest$gene[marker_panel_manifest$present_in_object])

panel_coverage <- aggregate(present_in_object ~ panel, data = marker_panel_manifest, FUN = function(x) {
    c(present_n = sum(x), total_n = length(x), fraction = mean(x))
})

panel_coverage_expanded <- data.frame(panel = panel_coverage$panel, present_n = panel_coverage$present_in_object[, 
    "present_n"], total_n = panel_coverage$present_in_object[, "total_n"], coverage_fraction = panel_coverage$present_in_object[, 
    "fraction"], stringsAsFactors = FALSE)

write_csv_utf8(panel_coverage_expanded, file.path(table_root, "05_CANONICAL_MARKER_PANEL_COVERAGE.csv"))

panel_cells <- unlist(lapply(cluster_levels, function(cluster_id) {
    balanced_sample_cells(cells = cells_by_cluster[[cluster_id]], sample_by_cell = sample_by_cell, max_total = 500L)
}), use.names = FALSE)

panel_matrix <- SeuratObject::LayerData(object = object, assay = "RNA", layer = "data", features = panel_genes_present, 
    cells = panel_cells)

if (!inherits(panel_matrix, "sparseMatrix")) {
    panel_matrix <- Matrix::Matrix(panel_matrix, sparse = TRUE)
}

panel_cell_cluster <- cell_metadata[colnames(panel_matrix), "cluster"]

panel_rows <- list()

for (cluster_id in cluster_levels) {
    idx <- which(panel_cell_cluster == cluster_id)
    one_matrix <- panel_matrix[, idx, drop = FALSE]
    panel_rows[[length(panel_rows) + 1L]] <- data.frame(cluster = cluster_id, gene = rownames(panel_matrix), average_log_normalized_expression = as.numeric(Matrix::rowMeans(one_matrix)), 
        detection_fraction = as.numeric(Matrix::rowMeans(one_matrix > 0)), sampled_cell_n = length(idx), stringsAsFactors = FALSE)
}

panel_long <- do.call(rbind, panel_rows)

panel_long <- merge(panel_long, marker_panel_manifest[marker_panel_manifest$present_in_object, c("panel", "gene", 
    "panel_role")], by = "gene", all.x = TRUE, sort = FALSE)

panel_long$cluster_numeric <- suppressWarnings(as.integer(panel_long$cluster))

panel_long <- panel_long[order(panel_long$cluster_numeric, panel_long$panel, panel_long$gene), , drop = FALSE]

scaled_rows <- list()

for (gene_id in unique(panel_long$gene)) {
    one <- panel_long[panel_long$gene == gene_id, , drop = FALSE]
    values <- one$average_log_normalized_expression
    value_sd <- stats::sd(values)
    one$scaled_average_expression <- if (!is.finite(value_sd) || value_sd == 0) {
        rep(0, length(values))
    }
    else {
        as.numeric(scale(values))
    }
    scaled_rows[[length(scaled_rows) + 1L]] <- one
}

panel_long <- do.call(rbind, scaled_rows)

write_csv_utf8(panel_long, file.path(source_data_root, "SOURCE_DATA_canonical_marker_panel_by_cluster.csv"))

candidate_panels <- setdiff(unique(marker_panel_manifest$panel), c("immune_general", "cycling"))

panel_support <- aggregate(cbind(scaled_average_expression, detection_fraction) ~ cluster + panel, data = panel_long[panel_long$panel %in% 
    candidate_panels, , drop = FALSE], FUN = mean)

colnames(panel_support)[colnames(panel_support) == "scaled_average_expression"] <- "mean_scaled_expression"

colnames(panel_support)[colnames(panel_support) == "detection_fraction"] <- "mean_detection_fraction"

panel_gene_counts <- aggregate(gene ~ cluster + panel, data = panel_long[panel_long$panel %in% candidate_panels, 
    , drop = FALSE], FUN = length)

colnames(panel_gene_counts)[3] <- "gene_n_present"

panel_support <- merge(panel_support, panel_gene_counts, by = c("cluster", "panel"), all.x = TRUE, sort = FALSE)

panel_support <- panel_support[order(suppressWarnings(as.integer(panel_support$cluster)), -panel_support$mean_scaled_expression), 
    , drop = FALSE]

write_csv_utf8(panel_support, file.path(table_root, "06_BROAD_LINEAGE_PANEL_SUPPORT.csv"))

top_panel_rows <- lapply(split(panel_support, panel_support$cluster), function(x) {
    x <- x[order(-x$mean_scaled_expression, -x$mean_detection_fraction), , drop = FALSE]
    first <- x[1, , drop = FALSE]
    second_score <- if (nrow(x) >= 2L) {
        x$mean_scaled_expression[2]
    }
    else {
        NA_real_
    }
    data.frame(cluster = first$cluster, top_support_panel = first$panel, top_support_score = first$mean_scaled_expression, 
        second_support_score = second_score, support_margin = (first$mean_scaled_expression - second_score), top_panel_gene_n_present = first$gene_n_present, 
        top_panel_detection_fraction = first$mean_detection_fraction, annotation_status = "ALGORITHMIC_SUPPORT_ONLY_REQUIRES_MANUAL_REVIEW", 
        stringsAsFactors = FALSE)
})

top_panel_summary <- do.call(rbind, top_panel_rows)

top_panel_summary$review_cluster <- (top_panel_summary$cluster %in% review_clusters)

top_panel_summary <- top_panel_summary[numeric_cluster_order(top_panel_summary$cluster), , drop = FALSE]

write_csv_utf8(top_panel_summary, file.path(table_root, "07_TOP_BROAD_LINEAGE_SUPPORT_BY_CLUSTER.csv"))

theme_qc <- ggplot2::theme_bw(base_size = 9) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95"))

panel_long$cluster <- factor(panel_long$cluster, levels = rev(cluster_levels))

panel_long$gene <- factor(panel_long$gene, levels = unique(marker_panel_manifest$gene[marker_panel_manifest$present_in_object]))

plot_dot <- ggplot2::ggplot(panel_long, ggplot2::aes(x = gene, y = cluster, size = detection_fraction, color = scaled_average_expression)) + 
    ggplot2::geom_point() + ggplot2::facet_grid(. ~ panel, scales = "free_x", space = "free_x") + ggplot2::scale_size_continuous(range = c(0.20000000000000001, 
    5), labels = scales::percent_format(accuracy = 1)) + ggplot2::scale_color_viridis_c(option = "C") + ggplot2::labs(title = "GSE250130 broad-lineage marker support", 
    subtitle = paste("Balanced sampled expression;", "not a locked cell-type annotation"), x = "Canonical marker", 
    y = "Cluster", size = "Detected", color = "Scaled\naverage") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
    vjust = 0.5, hjust = 1, size = 6), strip.text.x = ggplot2::element_text(angle = 90, size = 7))

save_plot_pair(plot_dot, "QC_Figure_1_broad_lineage_marker_dotplot", width = 20, height = 8)

if (nrow(top25)) {
    top5 <- top25[, head(.SD, 5L), by = cluster]
    top5$cluster <- factor(top5$cluster, levels = cluster_levels)
    top5$gene_cluster <- paste0(top5$gene, " [C", top5$cluster, "]")
    top5$gene_cluster <- factor(top5$gene_cluster, levels = rev(top5$gene_cluster))
    plot_top_markers <- ggplot2::ggplot(as.data.frame(top5), ggplot2::aes(x = cluster, y = gene_cluster, fill = marker_score)) + 
        ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "B") + ggplot2::labs(title = "Top five descriptive markers per cluster", 
        subtitle = paste("Balanced sparse marker scoring;", "annotation support only"), x = "Cluster", y = "Marker", 
        fill = "Marker\nscore") + theme_qc + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))
    save_plot_pair(plot_top_markers, "QC_Figure_2_top5_markers_per_cluster", width = 10, height = 18)
}

plot_support <- ggplot2::ggplot(panel_support, ggplot2::aes(x = panel, y = factor(cluster, levels = rev(cluster_levels)), 
    fill = mean_scaled_expression)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "D") + ggplot2::labs(title = "Broad-lineage support matrix", 
    subtitle = "Panel support requires marker review", x = "Marker panel", y = "Cluster", fill = "Mean scaled\nexpression") + 
    theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

save_plot_pair(plot_support, "QC_Figure_3_broad_lineage_support_matrix", width = 11, height = 7)

plot_marker_n <- ggplot2::ggplot(marker_status_table, ggplot2::aes(x = factor(cluster, levels = cluster_levels), 
    y = marker_n, fill = status)) + ggplot2::geom_col() + ggplot2::scale_y_continuous(labels = scales::comma) + 
    ggplot2::labs(title = "Descriptive marker yield by cluster", x = "Cluster", y = "Returned marker count", fill = "Status") + 
    theme_qc

save_plot_pair(plot_marker_n, "QC_Figure_4_marker_yield_by_cluster", width = 8, height = 5)

all_cluster_markers_complete <- (nrow(marker_status_table) == length(cluster_levels) && all(marker_status_table$status == 
    "COMPLETE"))

panel_coverage_pass <- (mean(marker_panel_manifest$present_in_object) >= 0.80000000000000004)

marker_yield_pass <- (all(marker_status_table$marker_n >= 10L))

gate_status <- if (all_cluster_markers_complete && panel_coverage_pass && marker_yield_pass) {
    "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED"
} else {
    "HOLD_DESCRIPTIVE_MARKER_OR_PANEL_QC_INCOMPLETE"
}

marker_gate <- data.frame(dataset_id = dataset_id, input_version = input_version, marker_version = marker_version, 
    supersedes_marker_version = "MARKERS_v1_0_FAILED_FINDMARKERS", gate_status = gate_status, marker_method = "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING", 
    p_values_reported = FALSE, cluster_n_expected = length(cluster_levels), cluster_marker_complete_n = sum(marker_status_table$status == 
        "COMPLETE"), cluster_marker_failed_n = sum(marker_status_table$status != "COMPLETE"), failed_clusters = collapse_unique(failed_clusters), 
    canonical_marker_present_fraction = mean(marker_panel_manifest$present_in_object), minimum_marker_n = if (nrow(marker_status_table)) {
        min(marker_status_table$marker_n)
    } else {
        NA_integer_
    }, review_clusters = collapse_unique(review_clusters), annotation_review_authorized = (gate_status == "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED"), 
    biological_annotation_locked = FALSE, VDR_localization_authorized = FALSE, pseudobulk_authorized = FALSE, result_lock_authorized = FALSE, 
    stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "STEP04E1D_GSE250130_MARKER_GATE.csv")

write_csv_utf8(marker_gate, gate_path)

card_lines <- c("# STEP04E1D Marker Result Card", "", paste0("- Dataset: ", dataset_id), paste0("- Input: ", input_version), 
    paste0("- Marker version: ", marker_version), paste0("- Gate: ", gate_status), paste0("- Completed clusters: ", 
        sum(marker_status_table$status == "COMPLETE"), "/", length(cluster_levels)), paste0("- Review clusters: ", 
        ifelse(length(review_clusters), collapse_unique(review_clusters), "none")), paste0("- Marker method: ", 
        marker_gate$marker_method), "- P-values reported: FALSE", "- Marker role: cluster annotation support only", 
    "- Donor-level inferential role: none", paste0("- Annotation review authorized: ", marker_gate$annotation_review_authorized), 
    "- Biological annotation locked: FALSE", "- VDR localization authorized: FALSE", "- Pseudobulk authorized: FALSE", 
    "- Result lock authorized: FALSE")

card_path <- file.path(run_qc_root, "STEP04E1D_GSE250130_MARKER_RESULT_CARD.md")

write_text(card_lines, card_path)

interpretation_lines <- c("# STEP04E1D GSE250130 descriptive cluster markers", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), paste0("- Completed clusters: ", sum(marker_status_table$status == 
    "COMPLETE"), "/", length(cluster_levels)), paste0("- Canonical marker-panel coverage: ", round(100 * mean(marker_panel_manifest$present_in_object), 
    1), "%"), "", "## Scientific interpretation rule", "", "- The failed MARKERS_v1_0 run is retained as a provenance record.", 
    "- MARKERS_v1_1 uses explicitly balanced, sampled sparse matrices and does not call FindMarkers.", "- Marker scores are descriptive and used only for cell-identity annotation.", 
    "- No cell-level p-value is reported or interpreted.", "- Clinical and temporal inference remains reserved for donor-aware pseudobulk.", 
    "- Clusters 12, 14, 20, and 21 require explicit review because of prior sample-mixing flags.", "", "## Next step", 
    "", "- Review top descriptive markers and canonical marker-panel support.", "- Construct an explicit cluster-to-cell-type mapping in STEP04E1E.", 
    "- Do not calculate final VDR localization or donor-aware pseudobulk until the annotation Gate passes.")

interpretation_path <- file.path(run_qc_root, "STEP04E1D_GSE250130_MARKER_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

methods_lines <- c("# STEP04E1D Methods draft", "", paste("Cluster markers were derived from the joined log-normalized RNA data", 
    "for the resolution-0.4 clustering solution using an explicitly", "memory-bounded descriptive scoring workflow."), 
    paste("For each cluster, up to 1,500 target cells were sampled while balancing", "representation across contributing samples."), 
    paste("The reference set was assembled by sampling up to 150 cells from each", "remaining cluster, also balancing cells across samples."), 
    paste("Only the selected cells were retrieved from the sparse normalized-data", "layer. For every gene, mean log-normalized expression, detection", 
        "fractions, detection-fraction difference, target-sample breadth, and a", "composite marker score were calculated."), 
    paste("Primary marker candidates required a target detection fraction of at", "least 10%, a detection-fraction difference of at least 5 percentage", 
        "points, a mean log-expression difference of at least 0.20, and detection", "in at least 25% of target-contributing samples."), 
    paste("A documented relaxed threshold was used only when fewer than 25 primary", "markers were available for a cluster."), 
    paste("These descriptive marker statistics were used solely for cell-identity", "annotation. No cell-level p-values were used for clinical or temporal", 
        "inference, which was reserved for subsequent donor-aware pseudobulk", "analysis."))

write_text(methods_lines, file.path(manuscript_text_root, "STEP04E1D_GSE250130_marker_Methods_draft_v1_1.md"))

run_end <- Sys.time()

package_names <- c("SeuratObject", "Matrix", "data.table", "ggplot2", "jsonlite", "scales", "future")

package_table <- data.frame(package = package_names, installed = vapply(package_names, requireNamespace, logical(1), 
    quietly = TRUE), version = vapply(package_names, safe_package_version, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- unique(c(list.files(marker_result_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1D_GSE250130_marker_Methods_draft_v1_1.md")))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, input_version = input_version, 
    marker_version = marker_version, supersedes_marker_version = "MARKERS_v1_0_FINDMARKERS_MEMORY_FAILURE", dataset_id = dataset_id, 
    started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, 
        run_start, units = "mins")), 2), gate_status = gate_status, cluster_n = length(cluster_levels), completed_cluster_n = sum(marker_status_table$status == 
        "COMPLETE"), failed_cluster_n = sum(marker_status_table$status != "COMPLETE"), marker_method = "BALANCED_SAMPLED_SPARSE_DESCRIPTIVE_SCORING", 
    p_values_reported = FALSE, inferential_role = "NONE_ANNOTATION_SUPPORT_ONLY", checkpoint_modified = FALSE, 
    biological_annotation_locked = FALSE, result_lock_authorized = FALSE, marker_result_root = normalize_slashes(marker_result_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "completed_clusters=", sum(marker_status_table$status == 
    "COMPLETE"), "failed_clusters=", collapse_unique(failed_clusters))

rm(panel_matrix, object)

gc()

small_files <- unique(c(final_marker_path, file.path(table_root, paste0(dataset_id, "_top25_descriptive_markers_by_cluster.csv")), 
    file.path(table_root, paste0(dataset_id, "_top10_marker_summary.csv")), list.files(table_root, recursive = TRUE, 
        full.names = TRUE), list.files(source_data_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
        recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, normalizePath(.vdr_script_path(script_name), 
        winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    150 * 1024^2]

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

cat("STEP04E1D GSE250130 memory-safe marker analysis completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Completed clusters: ", sum(marker_status_table$status == "COMPLETE"), "/", length(cluster_levels), "\n", sep = "")

cat("Failed clusters: ", ifelse(length(failed_clusters), collapse_unique(failed_clusters), "none"), "\n", sep = "")

cat("Marker result root: ", marker_result_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

cat("Run manifest: ", run_manifest_path, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe result and QC folders are complete.\n", sep = "")
}

