#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D03: Assess independent integration and donor representation

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "assess_gse183837_integration.R"

script_version <- "STEP04E1K_GSE183837_DETAILED_INTEGRATION_QC_v1.2_PATHSAFE"

integration_version <- "SCI183_v1_0"

QC_version <- "IQ183_v1_2"

run_start <- Sys.time()

run_id <- paste0("STEP04E1K_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

write_csv_gz <- function(x, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- gzfile(path, open = "wt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    write.csv(x, connection, row.names = FALSE, na = "")
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

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

choose2 <- function(x) {
    x * (x - 1)/2
}

adjusted_rand_index <- function(x, y) {
    if (length(x) != length(y)) {
        stop("ARI inputs have different lengths.")
    }
    contingency <- table(x, y)
    n <- sum(contingency)
    if (n < 2L) 
        return(NA_real_)
    sum_comb_cells <- sum(choose2(contingency))
    row_comb <- sum(choose2(rowSums(contingency)))
    col_comb <- sum(choose2(colSums(contingency)))
    total_comb <- choose2(n)
    expected <- (row_comb * col_comb/total_comb)
    maximum <- (row_comb + col_comb)/2
    denominator <- maximum - expected
    if (!is.finite(denominator) || denominator == 0) {
        return(NA_real_)
    }
    (sum_comb_cells - expected)/denominator
}

normalized_entropy <- function(counts, category_n) {
    counts <- as.numeric(counts)
    counts <- counts[is.finite(counts) & counts > 0]
    if (!length(counts) || category_n <= 1) {
        return(0)
    }
    probabilities <- counts/sum(counts)
    -sum(probabilities * log(probabilities))/log(category_n)
}

effective_category_n <- function(counts) {
    counts <- as.numeric(counts)
    counts <- counts[is.finite(counts) & counts > 0]
    if (!length(counts)) 
        return(0)
    probabilities <- counts/sum(counts)
    exp(-sum(probabilities * log(probabilities)))
}

median_or_NA <- function(x) {
    x <- safe_numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        median(x)
    }
}

quantile_or_NA <- function(x, probability) {
    x <- safe_numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        as.numeric(stats::quantile(x, probs = probability, names = FALSE, type = 7))
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

make_review_flags <- function(cell_n, top_donor_fraction, donor_entropy, Control_donors_ge10, RIF_donors_ge10, 
    median_percent_mt, median_nFeature, median_doublet_score) {
    flags <- character()
    if (is.finite(cell_n) && cell_n < 100) {
        flags <- c(flags, "RARE_CLUSTER_BELOW_100_CELLS")
    }
    if (is.finite(top_donor_fraction) && top_donor_fraction >= 0.75) {
        flags <- c(flags, "SEVERE_SINGLE_DONOR_DOMINANCE_GE_75_PERCENT")
    }
    else if (is.finite(top_donor_fraction) && top_donor_fraction >= 0.5) {
        flags <- c(flags, "SINGLE_DONOR_ENRICHMENT_GE_50_PERCENT")
    }
    if (is.finite(donor_entropy) && donor_entropy < 0.5) {
        flags <- c(flags, "LOW_DONOR_ENTROPY_BELOW_0_50")
    }
    else if (is.finite(donor_entropy) && donor_entropy < 0.69999999999999996) {
        flags <- c(flags, "MODERATE_DONOR_ENTROPY_BELOW_0_70")
    }
    if (is.finite(Control_donors_ge10) && Control_donors_ge10 < 2) {
        flags <- c(flags, "LOW_CONTROL_DONOR_SUPPORT_GE10_CELLS")
    }
    if (is.finite(RIF_donors_ge10) && RIF_donors_ge10 < 3) {
        flags <- c(flags, "LOW_RIF_DONOR_SUPPORT_GE10_CELLS")
    }
    if (is.finite(median_percent_mt) && is.finite(median_nFeature) && median_percent_mt > 20 && median_nFeature < 
        1000) {
        flags <- c(flags, "POSSIBLE_LOW_QUALITY_CLUSTER")
    }
    if (is.finite(median_doublet_score) && median_doublet_score > 0.20000000000000001) {
        flags <- c(flags, "POSSIBLE_RESIDUAL_DOUBLET_CLUSTER")
    }
    collapse_unique(flags)
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

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

integration_root <- file.path(results_root, "S4E1_sc", integration_version)

input_checkpoint <- file.path(integration_root, "checkpoints", dataset_id, "integrated.rds")

output_root <- file.path(integration_root, "iqc", "G183", QC_version)

table_root <- file.path(output_root, "t")

source_data_root <- file.path(output_root, "src")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "m")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1K_", format(run_start, "%Y%m%d_%H%M%S"), "_dQC"))

figure_root <- file.path(run_qc_root, "figures")

for (path_value in c(execution_root, parameter_root, environment_root, output_root, table_root, source_data_root, 
    figure_output_root, manifest_root, run_qc_root, figure_root)) {
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

planned_paths <- data.frame(artifact = c("input_checkpoint", "output_root", "table_root", "source_data_root", "figure_output_root", 
    "manifest_root", "run_qc_root", "handoff_zip"), path = normalize_slashes(c(input_checkpoint, output_root, table_root, 
    source_data_root, figure_output_root, manifest_root, run_qc_root, handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

if (!all(planned_paths$path_safe)) {
    stop("One or more planned paths are too long:\n", paste(planned_paths$path[!planned_paths$path_safe], collapse = "\n"))
}

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

planned_output_files <- data.frame(artifact = c("cluster_merge_audit", "checkpoint_audit", "donor_identity", "cluster_sample", 
    "cluster_QC", "review_clusters", "mixing_sample", "mixing_cluster", "graph_audit", "resolution_ARI", "resolution_02_to_04", 
    "resolution_04_to_06", "cluster_stability", "target_review", "source_cluster_sample", "source_cell_mixing", 
    "persistent_manifest", "detailed_gate", "detailed_card", "detailed_interpretation"), path = normalize_slashes(c(file.path(table_root, 
    "00_cluster_merge_count_audit.csv"), file.path(table_root, "01_checkpoint_audit.csv"), file.path(table_root, 
    "02_donor_group_identity.csv"), file.path(table_root, "03_cluster_sample.csv"), file.path(table_root, "04_cluster_QC.csv"), 
    file.path(table_root, "05_review_clusters.csv"), file.path(table_root, "06_mixing_sample.csv"), file.path(table_root, 
        "07_mixing_cluster.csv"), file.path(table_root, "08_graph_audit.csv"), file.path(table_root, "09_resolution_ARI.csv"), 
    file.path(table_root, "10_res02_to_04.csv"), file.path(table_root, "11_res04_to_06.csv"), file.path(table_root, 
        "12_cluster_stability.csv"), file.path(table_root, "13_review_7_8_13_15.csv"), file.path(source_data_root, 
        "SD1_cluster_sample.csv"), file.path(source_data_root, "SD2_cell_mixing.csv.gz"), file.path(manifest_root, 
        "persistent_manifest.csv"), file.path(run_qc_root, "G183_DETAILED_QC_GATE.csv"), file.path(run_qc_root, 
        "G183_DETAILED_QC_CARD.md"), file.path(run_qc_root, "G183_DETAILED_QC_INTERPRETATION.md"))), stringsAsFactors = FALSE)

planned_output_files$character_n <- nchar(planned_output_files$path, type = "chars")

planned_output_files$path_safe <- (planned_output_files$character_n <= 235L)

write_csv_utf8(planned_output_files, file.path(run_qc_root, "00_output_file_path_audit.csv"))

if (!all(planned_output_files$path_safe)) {
    stop("One or more planned output files exceed the safe path limit:\n", paste(planned_output_files$path[!planned_output_files$path_safe], 
        collapse = "\n"))
}

integration_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

integration_qc_dirs <- integration_qc_dirs[grepl("^STEP04E1J_[0-9]{8}_[0-9]{6}_integration$", basename(integration_qc_dirs), 
    perl = TRUE)]

if (!length(integration_qc_dirs)) {
    stop("No STEP04E1J integration-QC directory was found.")
}

candidate_gate_rows <- list()

for (directory_path in integration_qc_dirs) {
    gate_path <- file.path(directory_path, "GSE183837_INTEGRATION_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        candidate_gate_rows[[length(candidate_gate_rows) + 1L]] <- gate_value
    }
}

if (!length(candidate_gate_rows)) {
    stop("No readable GSE183837 integration Gate was found.")
}

candidate_gates <- as.data.frame(data.table::rbindlist(candidate_gate_rows, use.names = TRUE, fill = TRUE))

candidate_gates$source_mtime_POSIX <- as.POSIXct(candidate_gates$source_mtime)

candidate_gates <- candidate_gates[order(candidate_gates$source_mtime_POSIX, decreasing = TRUE), , drop = FALSE]

authorized_gate_index <- which(candidate_gates$gate_status %in% c("PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED", 
    "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS") & vapply(candidate_gates$detailed_integration_QC_authorized, 
    function(x) {
        isTRUE(x) || toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES")
    }, logical(1)))

if (!length(authorized_gate_index)) {
    stop("No STEP04E1J Gate authorizes detailed integration QC.")
}

upstream_gate <- as.data.frame(candidate_gates[authorized_gate_index[1], , drop = FALSE])

write_csv_utf8(upstream_gate, file.path(run_qc_root, "00_upstream_integration_gate.csv"))

if (!file.exists(input_checkpoint)) {
    stop("Integrated checkpoint is missing:\n", input_checkpoint)
}

free_memory_before_GB <- get_free_memory_gb()

free_disk_before_GB <- get_free_disk_gb()

if (is.finite(free_memory_before_GB) && free_memory_before_GB < 10) {
    stop("Free physical memory is below 10 GB: ", free_memory_before_GB)
}

log_event("Starting", run_id, "free_memory_GB=", free_memory_before_GB, "free_disk_GB=", free_disk_before_GB)

integrated_object <- readRDS(input_checkpoint)

if (!inherits(integrated_object, "Seurat")) {
    stop("The integrated checkpoint is not a Seurat object.")
}

metadata <- integrated_object[[]]

required_metadata_columns <- c("sample_accession", "donor_id", "analysis_group", "seurat_clusters", "harmony_snn_res.0.2", 
    "harmony_snn_res.0.4", "harmony_snn_res.0.6", "nCount_RNA", "nFeature_RNA", "percent.mt", "scDblFinder_score")

missing_metadata_columns <- setdiff(required_metadata_columns, colnames(metadata))

if (length(missing_metadata_columns)) {
    stop("Integrated metadata is missing required columns: ", paste(missing_metadata_columns, collapse = ", "))
}

if (ncol(integrated_object) != 47161L || nrow(metadata) != 47161L || length(unique(as.character(metadata$sample_accession))) != 
    9L || length(unique(as.character(metadata$donor_id))) != 9L || length(unique(as.character(metadata$seurat_clusters))) != 
    16L) {
    stop("Integrated checkpoint identity does not match the locked 47,161-cell, 9-donor, 16-cluster object.")
}

metadata$cell_barcode <- rownames(metadata)

metadata$sample_accession <- as.character(metadata$sample_accession)

metadata$donor_id <- as.character(metadata$donor_id)

metadata$analysis_group <- factor(as.character(metadata$analysis_group), levels = c("CONTROL", "RIF"))

metadata$seurat_clusters <- factor(as.character(metadata$seurat_clusters), levels = as.character(0:15))

if (anyNA(metadata$analysis_group) || anyNA(metadata$seurat_clusters)) {
    stop("Group or primary cluster metadata contains unexpected values.")
}

checkpoint_audit <- data.frame(metric = c("feature_n", "cell_n", "sample_n", "donor_n", "group_n", "primary_cluster_n", 
    "checkpoint_size_GB", "checkpoint_md5"), value = c(nrow(integrated_object), ncol(integrated_object), length(unique(metadata$sample_accession)), 
    length(unique(metadata$donor_id)), length(unique(metadata$analysis_group)), length(unique(metadata$seurat_clusters)), 
    round(file.info(input_checkpoint)$size/1024^3, 3), unname(tools::md5sum(input_checkpoint))), stringsAsFactors = FALSE)

write_csv_utf8(checkpoint_audit, file.path(table_root, "01_checkpoint_audit.csv"))

metadata_DT <- data.table::as.data.table(metadata)

donor_group_identity <- unique(metadata_DT[, .(donor_id, sample_accession, analysis_group)])

data.table::setorder(donor_group_identity, donor_id)

if (nrow(donor_group_identity) != 9L || anyDuplicated(donor_group_identity$donor_id) > 0L) {
    stop("Donor/group identity is not one row per nine donors.")
}

write_csv_utf8(as.data.frame(donor_group_identity), file.path(table_root, "02_donor_group_identity.csv"))

cluster_sample_observed <- metadata_DT[, .(cell_n = .N, median_nCount_RNA = median_or_NA(nCount_RNA), median_nFeature_RNA = median_or_NA(nFeature_RNA), 
    median_percent_mt = median_or_NA(percent.mt), median_scDblFinder_score = median_or_NA(scDblFinder_score)), 
    by = .(seurat_clusters, sample_accession, donor_id, analysis_group)]

cluster_levels <- as.character(0:15)

sample_levels <- donor_group_identity$sample_accession

cluster_sample_grid <- as.data.frame(data.table::CJ(seurat_clusters = cluster_levels, sample_accession = sample_levels, 
    unique = TRUE))

cluster_sample_grid <- merge(cluster_sample_grid, as.data.frame(donor_group_identity), by = "sample_accession", 
    all.x = TRUE, sort = FALSE)

cluster_sample_full <- as.data.frame(merge(cluster_sample_grid, as.data.frame(cluster_sample_observed), by = c("seurat_clusters", 
    "sample_accession", "donor_id", "analysis_group"), all.x = TRUE, sort = FALSE))

cluster_sample_full$cell_n[is.na(cluster_sample_full$cell_n)] <- 0L

cluster_sample_full$seurat_clusters <- factor(cluster_sample_full$seurat_clusters, levels = cluster_levels)

cluster_sample_full <- cluster_sample_full[order(as.integer(as.character(cluster_sample_full$seurat_clusters)), 
    match(cluster_sample_full$sample_accession, sample_levels)), , drop = FALSE]

cluster_totals <- aggregate(cell_n ~ seurat_clusters, data = cluster_sample_full, FUN = sum)

colnames(cluster_totals)[2] <- "cluster_cell_n"

cluster_sample_full <- merge(cluster_sample_full, cluster_totals, by = "seurat_clusters", all.x = TRUE, sort = FALSE)

cluster_sample_full$within_cluster_fraction <- (cluster_sample_full$cell_n/cluster_sample_full$cluster_cell_n)

sample_totals <- aggregate(cell_n ~ sample_accession, data = cluster_sample_full, FUN = sum)

colnames(sample_totals)[2] <- "sample_cell_n"

cluster_sample_full <- merge(cluster_sample_full, sample_totals, by = "sample_accession", all.x = TRUE, sort = FALSE)

cluster_sample_full$within_sample_fraction <- (cluster_sample_full$cell_n/cluster_sample_full$sample_cell_n)

write_csv_utf8(cluster_sample_full, file.path(table_root, "03_cluster_sample.csv"))

write_csv_utf8(cluster_sample_full, file.path(source_data_root, "SD1_cluster_sample.csv"))

cluster_QC <- metadata_DT[, .(technical_cell_n = .N, median_nCount_RNA = median_or_NA(nCount_RNA), Q1_nCount_RNA = quantile_or_NA(nCount_RNA, 
    0.25), Q3_nCount_RNA = quantile_or_NA(nCount_RNA, 0.75), median_nFeature_RNA = median_or_NA(nFeature_RNA), 
    Q1_nFeature_RNA = quantile_or_NA(nFeature_RNA, 0.25), Q3_nFeature_RNA = quantile_or_NA(nFeature_RNA, 0.75), 
    median_percent_mt = median_or_NA(percent.mt), Q1_percent_mt = quantile_or_NA(percent.mt, 0.25), Q3_percent_mt = quantile_or_NA(percent.mt, 
        0.75), median_scDblFinder_score = median_or_NA(scDblFinder_score), Q90_scDblFinder_score = quantile_or_NA(scDblFinder_score, 
        0.90000000000000002), Q95_scDblFinder_score = quantile_or_NA(scDblFinder_score, 0.94999999999999996)), 
    by = seurat_clusters]

cluster_composition_rows <- list()

for (cluster_value in cluster_levels) {
    cluster_data <- cluster_sample_full[as.character(cluster_sample_full$seurat_clusters) == cluster_value, , drop = FALSE]
    cluster_size <- sum(cluster_data$cell_n)
    top_index <- which.max(cluster_data$cell_n)
    group_cell_counts <- tapply(cluster_data$cell_n, cluster_data$analysis_group, sum)
    Control_cell_n <- if ("CONTROL" %in% names(group_cell_counts)) {
        as.integer(group_cell_counts["CONTROL"])
    }
    else {
        0L
    }
    RIF_cell_n <- if ("RIF" %in% names(group_cell_counts)) {
        as.integer(group_cell_counts["RIF"])
    }
    else {
        0L
    }
    cluster_composition_rows[[length(cluster_composition_rows) + 1L]] <- data.frame(seurat_clusters = cluster_value, 
        cell_n = cluster_size, cell_fraction = cluster_size/nrow(metadata), represented_donor_n = sum(cluster_data$cell_n > 
            0), donors_ge5_cells = sum(cluster_data$cell_n >= 5), donors_ge10_cells = sum(cluster_data$cell_n >= 
            10), Control_donors_ge1_cell = sum(cluster_data$analysis_group == "CONTROL" & cluster_data$cell_n > 
            0), RIF_donors_ge1_cell = sum(cluster_data$analysis_group == "RIF" & cluster_data$cell_n > 0), Control_donors_ge10_cells = sum(cluster_data$analysis_group == 
            "CONTROL" & cluster_data$cell_n >= 10), RIF_donors_ge10_cells = sum(cluster_data$analysis_group == 
            "RIF" & cluster_data$cell_n >= 10), top_donor = cluster_data$sample_accession[top_index], top_donor_group = as.character(cluster_data$analysis_group[top_index]), 
        top_donor_cell_n = cluster_data$cell_n[top_index], top_donor_fraction = cluster_data$cell_n[top_index]/cluster_size, 
        donor_entropy_normalized = normalized_entropy(cluster_data$cell_n, category_n = 9), effective_donor_n = effective_category_n(cluster_data$cell_n), 
        Control_cell_n = Control_cell_n, RIF_cell_n = RIF_cell_n, Control_cell_fraction = Control_cell_n/cluster_size, 
        RIF_cell_fraction = RIF_cell_n/cluster_size, stringsAsFactors = FALSE)
}

cluster_composition <- do.call(rbind, cluster_composition_rows)

cluster_review <- merge(cluster_composition, as.data.frame(cluster_QC), by = "seurat_clusters", all.x = TRUE, sort = FALSE)

required_cluster_review_columns <- c("seurat_clusters", "cell_n", "technical_cell_n", "top_donor_fraction", "donor_entropy_normalized", 
    "Control_donors_ge10_cells", "RIF_donors_ge10_cells", "median_percent_mt", "median_nFeature_RNA", "median_scDblFinder_score")

missing_cluster_review_columns <- setdiff(required_cluster_review_columns, colnames(cluster_review))

if (length(missing_cluster_review_columns)) {
    stop("cluster_review merge is missing required columns: ", paste(missing_cluster_review_columns, collapse = ", "))
}

if (nrow(cluster_review) != 16L || anyNA(cluster_review$seurat_clusters) || anyNA(cluster_review$cell_n) || anyNA(cluster_review$technical_cell_n)) {
    stop("cluster_review merge did not produce 16 complete cluster rows.")
}

cluster_count_consistency <- data.frame(seurat_clusters = as.character(cluster_review$seurat_clusters), composition_cell_n = as.integer(cluster_review$cell_n), 
    technical_cell_n = as.integer(cluster_review$technical_cell_n), counts_identical = as.integer(cluster_review$cell_n) == 
        as.integer(cluster_review$technical_cell_n), stringsAsFactors = FALSE)

write_csv_utf8(cluster_count_consistency, file.path(table_root, "00_cluster_merge_count_audit.csv"))

if (!all(cluster_count_consistency$counts_identical)) {
    stop("Cluster cell counts differ between composition and technical-QC tables.")
}

cluster_review$seurat_clusters <- as.character(cluster_review$seurat_clusters)

cluster_review <- cluster_review[match(cluster_levels, cluster_review$seurat_clusters), , drop = FALSE]

review_flag_vector <- mapply(make_review_flags, cell_n = cluster_review$cell_n, top_donor_fraction = cluster_review$top_donor_fraction, 
    donor_entropy = cluster_review$donor_entropy_normalized, Control_donors_ge10 = cluster_review$Control_donors_ge10_cells, 
    RIF_donors_ge10 = cluster_review$RIF_donors_ge10_cells, median_percent_mt = cluster_review$median_percent_mt, 
    median_nFeature = cluster_review$median_nFeature_RNA, median_doublet_score = cluster_review$median_scDblFinder_score, 
    SIMPLIFY = TRUE, USE.NAMES = FALSE)

if (length(review_flag_vector) != nrow(cluster_review)) {
    stop("Review-flag vector length mismatch: expected ", nrow(cluster_review), ", obtained ", length(review_flag_vector), 
        ".")
}

cluster_review$review_flags <- as.character(review_flag_vector)

cluster_review$review_required <- nzchar(cluster_review$review_flags)

cluster_review$critical_hold_flag <- ((cluster_review$top_donor_fraction >= 0.94999999999999996 & cluster_review$cell_n >= 
    200) | cluster_review$cell_n < 30 | (cluster_review$median_percent_mt > 22 & cluster_review$median_nFeature_RNA < 
    800) | cluster_review$median_scDblFinder_score > 0.29999999999999999)

write_csv_utf8(cluster_review, file.path(table_root, "04_cluster_QC.csv"))

write_csv_utf8(cluster_review[cluster_review$review_required, , drop = FALSE], file.path(table_root, "05_review_clusters.csv"))

graph_names <- names(integrated_object@graphs)

preferred_graph_names <- c("harmony_nn", "harmony_snn")

graph_name <- preferred_graph_names[preferred_graph_names %in% graph_names][1]

if (length(graph_name) != 1L || is.na(graph_name)) {
    stop("Neither harmony_nn nor harmony_snn graph is present.")
}

mixing_graph <- integrated_object@graphs[[graph_name]]

if (nrow(mixing_graph) != nrow(metadata) || ncol(mixing_graph) != nrow(metadata)) {
    stop("Harmony graph dimensions do not match metadata.")
}

graph_cells <- rownames(mixing_graph)

if (is.null(graph_cells)) {
    graph_cells <- rownames(metadata)
}

if (!identical(graph_cells, rownames(metadata))) {
    if (setequal(graph_cells, rownames(metadata))) {
        metadata_graph <- metadata[graph_cells, , drop = FALSE]
    }
    else {
        stop("Harmony graph cell identities do not match metadata.")
    }
} else {
    metadata_graph <- metadata
}

graph_edges <- data.table::as.data.table(Matrix::summary(mixing_graph))

required_edge_columns <- c("i", "j", "x")

if (!all(required_edge_columns %in% colnames(graph_edges))) {
    stop("Harmony graph summary lacks i, j and x columns.")
}

graph_edges <- graph_edges[i != j & is.finite(x) & x > 0]

sample_code <- as.integer(factor(metadata_graph$sample_accession, levels = donor_group_identity$sample_accession))

group_code <- as.integer(factor(metadata_graph$analysis_group, levels = c("CONTROL", "RIF")))

cluster_code <- as.character(metadata_graph$seurat_clusters)

graph_edges[, `:=`(same_sample = sample_code[i] == sample_code[j], same_group = group_code[i] == group_code[j], 
    neighbor_sample_code = sample_code[j])]

cell_mixing_basic <- graph_edges[, .(neighbor_n = .N, weighted_degree = sum(x), same_sample_fraction = mean(same_sample), 
    weighted_same_sample_fraction = sum(x * same_sample)/sum(x), same_group_fraction = mean(same_group), weighted_same_group_fraction = sum(x * 
        same_group)/sum(x)), by = i]

neighbor_sample_counts <- graph_edges[, .N, by = .(i, neighbor_sample_code)]

neighbor_diversity <- neighbor_sample_counts[, {
    probabilities <- N/sum(N)
    list(neighbor_donor_n = .N, neighbor_donor_entropy = -sum(probabilities * log(probabilities))/log(9), effective_neighbor_donor_n = exp(-sum(probabilities * 
        log(probabilities))))
}, by = i]

cell_mixing <- merge(cell_mixing_basic, neighbor_diversity, by = "i", all = TRUE, sort = FALSE)

cell_identity <- data.frame(i = seq_len(nrow(metadata_graph)), cell_barcode = rownames(metadata_graph), sample_accession = metadata_graph$sample_accession, 
    donor_id = metadata_graph$donor_id, analysis_group = as.character(metadata_graph$analysis_group), seurat_clusters = as.character(metadata_graph$seurat_clusters), 
    stringsAsFactors = FALSE)

cell_mixing <- merge(cell_identity, cell_mixing, by = "i", all.x = TRUE, sort = FALSE)

cell_mixing <- cell_mixing[order(cell_mixing$i), , drop = FALSE]

expected_sample_fraction <- cluster_sample_full[, c("seurat_clusters", "sample_accession", "cell_n", "cluster_cell_n")]

colnames(expected_sample_fraction)[colnames(expected_sample_fraction) == "cell_n"] <- "cluster_sample_cell_n"

cell_mixing <- merge(cell_mixing, expected_sample_fraction, by = c("seurat_clusters", "sample_accession"), all.x = TRUE, 
    sort = FALSE)

cell_mixing$expected_same_sample_fraction <- ifelse(cell_mixing$cluster_cell_n > 1, pmax(0, (cell_mixing$cluster_sample_cell_n - 
    1)/(cell_mixing$cluster_cell_n - 1)), 0)

cell_mixing$same_sample_excess <- (cell_mixing$same_sample_fraction - cell_mixing$expected_same_sample_fraction)

cell_mixing$weighted_same_sample_excess <- (cell_mixing$weighted_same_sample_fraction - cell_mixing$expected_same_sample_fraction)

write_csv_gz(cell_mixing, file.path(source_data_root, "SD2_cell_mixing.csv.gz"))

mixing_sample_summary <- data.table::as.data.table(cell_mixing)[, .(cell_n = .N, median_neighbor_n = median_or_NA(neighbor_n), 
    median_same_sample_fraction = median_or_NA(same_sample_fraction), Q90_same_sample_fraction = quantile_or_NA(same_sample_fraction, 
        0.90000000000000002), median_expected_same_sample_fraction = median_or_NA(expected_same_sample_fraction), 
    median_same_sample_excess = median_or_NA(same_sample_excess), Q90_same_sample_excess = quantile_or_NA(same_sample_excess, 
        0.90000000000000002), median_neighbor_donor_n = median_or_NA(neighbor_donor_n), median_neighbor_donor_entropy = median_or_NA(neighbor_donor_entropy)), 
    by = .(sample_accession, donor_id, analysis_group)]

mixing_cluster_summary <- data.table::as.data.table(cell_mixing)[, .(cell_n = .N, median_neighbor_n = median_or_NA(neighbor_n), 
    median_same_sample_fraction = median_or_NA(same_sample_fraction), Q90_same_sample_fraction = quantile_or_NA(same_sample_fraction, 
        0.90000000000000002), median_expected_same_sample_fraction = median_or_NA(expected_same_sample_fraction), 
    median_same_sample_excess = median_or_NA(same_sample_excess), Q90_same_sample_excess = quantile_or_NA(same_sample_excess, 
        0.90000000000000002), median_neighbor_donor_n = median_or_NA(neighbor_donor_n), median_neighbor_donor_entropy = median_or_NA(neighbor_donor_entropy)), 
    by = seurat_clusters]

mixing_cluster_summary$local_mixing_review_flag <- (mixing_cluster_summary$median_same_sample_excess > 0.25 | mixing_cluster_summary$Q90_same_sample_excess > 
    0.5 | mixing_cluster_summary$median_neighbor_donor_n < 3)

write_csv_utf8(as.data.frame(mixing_sample_summary), file.path(table_root, "06_mixing_sample.csv"))

write_csv_utf8(as.data.frame(mixing_cluster_summary), file.path(table_root, "07_mixing_cluster.csv"))

graph_audit <- data.frame(metric = c("graph_name", "graph_cell_n", "nonzero_edge_n_excluding_self", "median_neighbor_n", 
    "minimum_neighbor_n", "maximum_neighbor_n", "overall_median_same_sample_fraction", "overall_median_expected_same_sample_fraction", 
    "overall_median_same_sample_excess", "overall_median_neighbor_donor_n", "overall_median_neighbor_donor_entropy"), 
    value = c(graph_name, nrow(mixing_graph), nrow(graph_edges), median_or_NA(cell_mixing$neighbor_n), min(cell_mixing$neighbor_n, 
        na.rm = TRUE), max(cell_mixing$neighbor_n, na.rm = TRUE), median_or_NA(cell_mixing$same_sample_fraction), 
        median_or_NA(cell_mixing$expected_same_sample_fraction), median_or_NA(cell_mixing$same_sample_excess), 
        median_or_NA(cell_mixing$neighbor_donor_n), median_or_NA(cell_mixing$neighbor_donor_entropy)), stringsAsFactors = FALSE)

write_csv_utf8(graph_audit, file.path(table_root, "08_graph_audit.csv"))

rm(graph_edges, neighbor_sample_counts, cell_mixing_basic, neighbor_diversity)

gc()

resolution_02 <- as.character(metadata$harmony_snn_res.0.2)

resolution_04 <- as.character(metadata$harmony_snn_res.0.4)

resolution_06 <- as.character(metadata$harmony_snn_res.0.6)

resolution_ARI <- data.frame(comparison = c("0.2_vs_0.4", "0.4_vs_0.6", "0.2_vs_0.6"), adjusted_rand_index = c(adjusted_rand_index(resolution_02, 
    resolution_04), adjusted_rand_index(resolution_04, resolution_06), adjusted_rand_index(resolution_02, resolution_06)), 
    stringsAsFactors = FALSE)

write_csv_utf8(resolution_ARI, file.path(table_root, "09_resolution_ARI.csv"))

make_transition_table <- function(source_values, target_values, source_name, target_name) {
    contingency <- as.data.frame(table(source_values, target_values), stringsAsFactors = FALSE)
    colnames(contingency) <- c(source_name, target_name, "cell_n")
    source_totals <- aggregate(cell_n ~ ., data = contingency[, c(source_name, "cell_n")], FUN = sum)
    colnames(source_totals)[2] <- "source_cluster_n"
    target_totals <- aggregate(cell_n ~ ., data = contingency[, c(target_name, "cell_n")], FUN = sum)
    colnames(target_totals)[2] <- "target_cluster_n"
    output <- merge(contingency, source_totals, by = source_name, all.x = TRUE, sort = FALSE)
    output <- merge(output, target_totals, by = target_name, all.x = TRUE, sort = FALSE)
    output$within_source_fraction <- (output$cell_n/output$source_cluster_n)
    output$within_target_fraction <- (output$cell_n/output$target_cluster_n)
    output
}

transition_02_04 <- make_transition_table(resolution_02, resolution_04, "cluster_res_0_2", "cluster_res_0_4")

transition_04_06 <- make_transition_table(resolution_04, resolution_06, "cluster_res_0_4", "cluster_res_0_6")

write_csv_utf8(transition_02_04, file.path(table_root, "10_res02_to_04.csv"))

write_csv_utf8(transition_04_06, file.path(table_root, "11_res04_to_06.csv"))

primary_resolution_stability_rows <- list()

for (cluster_value in cluster_levels) {
    primary_index <- which(resolution_04 == cluster_value)
    parent_table <- sort(table(resolution_02[primary_index]), decreasing = TRUE)
    child_table <- sort(table(resolution_06[primary_index]), decreasing = TRUE)
    primary_resolution_stability_rows[[length(primary_resolution_stability_rows) + 1L]] <- data.frame(seurat_clusters = cluster_value, 
        cell_n = length(primary_index), dominant_res_0_2_parent = names(parent_table)[1], dominant_parent_fraction = as.numeric(parent_table[1])/length(primary_index), 
        res_0_2_parent_n = length(parent_table), dominant_res_0_6_child = names(child_table)[1], dominant_child_fraction = as.numeric(child_table[1])/length(primary_index), 
        res_0_6_child_n = length(child_table), stringsAsFactors = FALSE)
}

primary_resolution_stability <- do.call(rbind, primary_resolution_stability_rows)

write_csv_utf8(primary_resolution_stability, file.path(table_root, "12_cluster_stability.csv"))

target_review_clusters <- c("7", "8", "13", "15")

target_cluster_review <- cluster_review[cluster_review$seurat_clusters %in% target_review_clusters, , drop = FALSE]

target_cluster_review <- merge(target_cluster_review, as.data.frame(mixing_cluster_summary), by = "seurat_clusters", 
    all.x = TRUE, sort = FALSE, suffixes = c("", "_mixing"))

target_cluster_review <- merge(target_cluster_review, primary_resolution_stability, by = "seurat_clusters", all.x = TRUE, 
    sort = FALSE, suffixes = c("", "_resolution"))

target_cluster_review <- target_cluster_review[match(target_review_clusters, target_cluster_review$seurat_clusters), 
    , drop = FALSE]

target_cluster_review$provisional_QC_interpretation <- c("Donor-enriched cluster; retain for markers and determine whether it represents a donor-associated biological state.", 
    "Donor-enriched cluster with low transcript complexity but low mitochondrial and doublet scores; do not label as low quality before marker review.", 
    "Strongly dominated by one RIF donor; retain for marker extraction but prohibit disease interpretation unless the identity is reproduced across donors.", 
    "Rare but represented across all nine donors and not overtly low quality or doublet-like; retain as a candidate rare population for marker review.")

write_csv_utf8(target_cluster_review, file.path(table_root, "13_review_7_8_13_15.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

heatmap_data <- cluster_sample_full

heatmap_data$seurat_clusters <- factor(as.character(heatmap_data$seurat_clusters), levels = rev(cluster_levels))

heatmap_data$sample_accession <- factor(heatmap_data$sample_accession, levels = sample_levels)

plot_donor_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = sample_accession, y = seurat_clusters, fill = within_cluster_fraction)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) + ggplot2::labs(title = "GSE183837 donor contribution within each primary cluster", 
    subtitle = "Each row sums to 100%; high single-donor fractions require marker-level review", x = "Donor/sample", 
    y = "Primary cluster", fill = "Within-cluster fraction") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_donor_heatmap, "QC1_cluster_donor_heatmap", width = 10, height = 7)

dominance_plot_data <- cluster_review

dominance_plot_data$seurat_clusters <- factor(dominance_plot_data$seurat_clusters, levels = cluster_levels)

plot_dominance <- ggplot2::ggplot(dominance_plot_data, ggplot2::aes(x = top_donor_fraction, y = donor_entropy_normalized, 
    size = cell_n, label = seurat_clusters, shape = review_required)) + ggplot2::geom_vline(xintercept = c(0.5, 
    0.75), linetype = "dashed") + ggplot2::geom_hline(yintercept = c(0.5, 0.69999999999999996), linetype = "dotted") + 
    ggplot2::geom_point(alpha = 0.80000000000000004) + ggplot2::geom_text(size = 3, nudge_y = 0.025000000000000001, 
    show.legend = FALSE) + ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + ggplot2::labs(title = "Cluster donor dominance and donor entropy", 
    subtitle = "Clusters 7, 8 and 13 show donor enrichment; cluster 15 is rare but donor-distributed", x = "Largest single-donor fraction", 
    y = "Normalized donor entropy", size = "Cell count", shape = "Review flag") + theme_QC

save_plot_pair(plot_dominance, "QC2_cluster_donor_dominance_entropy", width = 8, height = 6)

group_plot_data <- cluster_sample_full

group_plot_data <- aggregate(cell_n ~ seurat_clusters + analysis_group, data = group_plot_data, FUN = sum)

group_totals <- aggregate(cell_n ~ seurat_clusters, data = group_plot_data, FUN = sum)

colnames(group_totals)[2] <- "cluster_cell_n"

group_plot_data <- merge(group_plot_data, group_totals, by = "seurat_clusters", all.x = TRUE, sort = FALSE)

group_plot_data$group_fraction <- (group_plot_data$cell_n/group_plot_data$cluster_cell_n)

group_plot_data$seurat_clusters <- factor(as.character(group_plot_data$seurat_clusters), levels = cluster_levels)

plot_group <- ggplot2::ggplot(group_plot_data, ggplot2::aes(x = seurat_clusters, y = group_fraction, fill = analysis_group)) + 
    ggplot2::geom_col(position = "stack") + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Control and RIF cell fractions by cluster", subtitle = "Descriptive composition only; no cell-level group inference", 
        x = "Primary cluster", y = "Cell fraction", fill = "Group") + theme_QC

save_plot_pair(plot_group, "QC3_cluster_group_composition", width = 9, height = 5)

technical_long <- data.frame(seurat_clusters = rep(cluster_review$seurat_clusters, times = 4), metric = rep(c("Median nFeature_RNA", 
    "Median nCount_RNA", "Median percent.mt", "Median scDblFinder score"), each = nrow(cluster_review)), value = c(cluster_review$median_nFeature_RNA, 
    cluster_review$median_nCount_RNA, cluster_review$median_percent_mt, cluster_review$median_scDblFinder_score), 
    stringsAsFactors = FALSE)

technical_long$seurat_clusters <- factor(technical_long$seurat_clusters, levels = cluster_levels)

plot_technical <- ggplot2::ggplot(technical_long, ggplot2::aes(x = seurat_clusters, y = value, group = 1)) + ggplot2::geom_point() + 
    ggplot2::geom_line() + ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) + ggplot2::labs(title = "Primary-cluster technical QC profiles", 
    subtitle = "Low transcript complexity alone is not sufficient to remove an immune or rare cluster", x = "Primary cluster", 
    y = "Median value") + theme_QC

save_plot_pair(plot_technical, "QC4_cluster_technical_profiles", width = 10, height = 7)

mixing_plot_data <- cell_mixing

mixing_plot_data$seurat_clusters <- factor(mixing_plot_data$seurat_clusters, levels = cluster_levels)

plot_mixing <- ggplot2::ggplot(mixing_plot_data, ggplot2::aes(x = seurat_clusters, y = same_sample_excess)) + ggplot2::geom_hline(yintercept = 0, 
    linetype = "dashed") + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::coord_cartesian(ylim = c(quantile_or_NA(mixing_plot_data$same_sample_excess, 
    0.01), quantile_or_NA(mixing_plot_data$same_sample_excess, 0.98999999999999999))) + ggplot2::labs(title = "Local same-sample enrichment after Harmony", 
    subtitle = "Observed same-sample neighbor fraction minus the expectation from cluster composition", x = "Primary cluster", 
    y = "Same-sample neighborhood excess") + theme_QC

save_plot_pair(plot_mixing, "QC5_local_sample_mixing_by_cluster", width = 9, height = 5.5)

transition_02_04_plot <- transition_02_04

transition_02_04_plot$cluster_res_0_2 <- factor(transition_02_04_plot$cluster_res_0_2, levels = rev(sort(unique(transition_02_04_plot$cluster_res_0_2))))

transition_02_04_plot$cluster_res_0_4 <- factor(transition_02_04_plot$cluster_res_0_4, levels = sort(unique(transition_02_04_plot$cluster_res_0_4)))

plot_transition_02_04 <- ggplot2::ggplot(transition_02_04_plot, ggplot2::aes(x = cluster_res_0_4, y = cluster_res_0_2, 
    fill = within_source_fraction)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Cluster transition from resolution 0.2 to 0.4", x = "Resolution 0.4 cluster", y = "Resolution 0.2 cluster", 
        fill = "Fraction of 0.2 cluster") + theme_QC

save_plot_pair(plot_transition_02_04, "QC6_resolution_transition_0_2_to_0_4", width = 8, height = 6)

transition_04_06_plot <- transition_04_06

transition_04_06_plot$cluster_res_0_4 <- factor(transition_04_06_plot$cluster_res_0_4, levels = rev(sort(unique(transition_04_06_plot$cluster_res_0_4))))

transition_04_06_plot$cluster_res_0_6 <- factor(transition_04_06_plot$cluster_res_0_6, levels = sort(unique(transition_04_06_plot$cluster_res_0_6)))

plot_transition_04_06 <- ggplot2::ggplot(transition_04_06_plot, ggplot2::aes(x = cluster_res_0_6, y = cluster_res_0_4, 
    fill = within_source_fraction)) + ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Cluster transition from resolution 0.4 to 0.6", x = "Resolution 0.6 cluster", y = "Resolution 0.4 cluster", 
        fill = "Fraction of 0.4 cluster") + theme_QC

save_plot_pair(plot_transition_04_06, "QC7_resolution_transition_0_4_to_0_6", width = 9, height = 6)

cluster_size_plot_data <- cluster_review

cluster_size_plot_data$seurat_clusters <- factor(cluster_size_plot_data$seurat_clusters, levels = cluster_levels)

plot_cluster_size <- ggplot2::ggplot(cluster_size_plot_data, ggplot2::aes(x = seurat_clusters, y = cell_n, shape = review_required)) + 
    ggplot2::geom_hline(yintercept = 100, linetype = "dashed") + ggplot2::geom_point(size = 3) + ggplot2::scale_y_log10(labels = scales::comma) + 
    ggplot2::labs(title = "Primary-cluster sizes and review status", subtitle = "Cluster 15 is below 100 cells but is retained for marker review", 
        x = "Primary cluster", y = "Cell count (log scale)", shape = "Review flag") + theme_QC

save_plot_pair(plot_cluster_size, "QC8_cluster_size_review", width = 8, height = 5)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

critical_cluster_n <- sum(cluster_review$critical_hold_flag, na.rm = TRUE)

review_cluster_ids <- cluster_review$seurat_clusters[cluster_review$review_required]

severe_donor_dominant_clusters <- cluster_review$seurat_clusters[cluster_review$top_donor_fraction >= 0.75]

moderate_donor_enriched_clusters <- cluster_review$seurat_clusters[cluster_review$top_donor_fraction >= 0.5 & cluster_review$top_donor_fraction < 
    0.75]

rare_clusters <- cluster_review$seurat_clusters[cluster_review$cell_n < 100]

local_mixing_review_clusters <- as.character(mixing_cluster_summary$seurat_clusters[mixing_cluster_summary$local_mixing_review_flag])

ARI_02_04 <- resolution_ARI$adjusted_rand_index[resolution_ARI$comparison == "0.2_vs_0.4"]

ARI_04_06 <- resolution_ARI$adjusted_rand_index[resolution_ARI$comparison == "0.4_vs_0.6"]

resolution_review_flag <- (!is.finite(ARI_02_04) || !is.finite(ARI_04_06) || ARI_02_04 < 0.59999999999999998 || 
    ARI_04_06 < 0.59999999999999998)

overall_review_flags <- c(if (length(severe_donor_dominant_clusters)) {
    paste0("SEVERE_DONOR_DOMINANCE_CLUSTERS_", paste(severe_donor_dominant_clusters, collapse = "_"))
} else {
    character()
}, if (length(moderate_donor_enriched_clusters)) {
    paste0("MODERATE_DONOR_ENRICHMENT_CLUSTERS_", paste(moderate_donor_enriched_clusters, collapse = "_"))
} else {
    character()
}, if (length(rare_clusters)) {
    paste0("RARE_CLUSTERS_", paste(rare_clusters, collapse = "_"))
} else {
    character()
}, if (length(local_mixing_review_clusters)) {
    paste0("LOCAL_MIXING_REVIEW_CLUSTERS_", paste(local_mixing_review_clusters, collapse = "_"))
} else {
    character()
}, if (resolution_review_flag) {
    "RESOLUTION_STABILITY_REVIEW"
} else {
    character()
})

overall_review_flags <- unique(overall_review_flags[nzchar(overall_review_flags)])

gate_checks <- data.frame(check = c("upstream_integration_gate_authorizes_detailed_QC", "checkpoint_cell_n_is_47161", 
    "checkpoint_sample_n_is_9", "checkpoint_donor_n_is_9", "primary_cluster_n_is_16", "Harmony_graph_available", 
    "Harmony_graph_dimensions_match_cells", "local_mixing_metrics_complete_for_all_cells", "all_clusters_represented_by_at_least_2_donors", 
    "no_cluster_has_ge95_percent_single_donor_dominance", "no_cluster_has_fewer_than_30_cells", "no_cluster_meets_composite_low_quality_hold_rule", 
    "no_cluster_has_median_doublet_score_above_0_30", "resolution_ARI_values_are_finite", "primary_resolution_stability_tables_complete", 
    "target_clusters_7_8_13_15_reviewed", "marker_extraction_can_be_run_without_removing_review_clusters", "free_disk_after_run_at_least_8GB_or_unknown"), 
    pass = c(upstream_gate$gate_status[[1]] %in% c("PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED", 
        "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS"), ncol(integrated_object) == 47161L, length(unique(metadata$sample_accession)) == 
        9L, length(unique(metadata$donor_id)) == 9L, length(unique(metadata$seurat_clusters)) == 16L, graph_name %in% 
        preferred_graph_names, nrow(mixing_graph) == nrow(metadata) && ncol(mixing_graph) == nrow(metadata), nrow(cell_mixing) == 
        nrow(metadata) && all(is.finite(cell_mixing$same_sample_fraction)), all(cluster_review$represented_donor_n >= 
        2), all(cluster_review$top_donor_fraction < 0.94999999999999996), all(cluster_review$cell_n >= 30), !any(cluster_review$median_percent_mt > 
        22 & cluster_review$median_nFeature_RNA < 800), all(cluster_review$median_scDblFinder_score <= 0.29999999999999999), 
        all(is.finite(resolution_ARI$adjusted_rand_index)), nrow(primary_resolution_stability) == 16L, all(target_review_clusters %in% 
            target_cluster_review$seurat_clusters), critical_cluster_n == 0L, is.na(get_free_disk_gb()) || get_free_disk_gb() >= 
            8), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

marker_analysis_authorized <- (all_gate_checks_pass && critical_cluster_n == 0L)

gate_status <- if (marker_analysis_authorized && !length(overall_review_flags)) {
    "PASS_DETAILED_INTEGRATION_QC_MARKERS_AUTHORIZED"
} else if (marker_analysis_authorized) {
    "PASS_DETAILED_INTEGRATION_QC_WITH_REVIEW_FLAGS_MARKERS_AUTHORIZED"
} else {
    "HOLD_DETAILED_INTEGRATION_QC_MARKERS_NOT_AUTHORIZED"
}

free_memory_after_GB <- get_free_memory_gb()

free_disk_after_GB <- get_free_disk_gb()

detailed_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, QC_version = QC_version, 
    gate_status = gate_status, cell_n = nrow(metadata), donor_n = length(unique(metadata$donor_id)), primary_cluster_n = length(unique(metadata$seurat_clusters)), 
    Harmony_graph = graph_name, severe_donor_dominant_clusters = collapse_unique(severe_donor_dominant_clusters), 
    moderate_donor_enriched_clusters = collapse_unique(moderate_donor_enriched_clusters), rare_clusters = collapse_unique(rare_clusters), 
    local_mixing_review_clusters = collapse_unique(local_mixing_review_clusters), ARI_res_0_2_vs_0_4 = ARI_02_04, 
    ARI_res_0_4_vs_0_6 = ARI_04_06, critical_hold_cluster_n = critical_cluster_n, review_flags = collapse_unique(overall_review_flags), 
    descriptive_marker_analysis_authorized = marker_analysis_authorized, cluster_removal_authorized = FALSE, biological_annotation_locked = FALSE, 
    disease_pseudobulk_testing_authorized = FALSE, result_lock_authorized = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(detailed_gate, file.path(run_qc_root, "G183_DETAILED_QC_GATE.csv"))

card_lines <- c("# STEP04E1K GSE183837 Detailed Integration QC Card", "", paste0("- Dataset: ", dataset_id), paste0("- Integration version: ", 
    integration_version), paste0("- QC version: ", QC_version), paste0("- Gate: ", gate_status), paste0("- Cells: ", 
    format(nrow(metadata), big.mark = ",")), paste0("- Donors: ", length(unique(metadata$donor_id))), paste0("- Primary clusters: ", 
    length(unique(metadata$seurat_clusters))), paste0("- Severe donor-dominant clusters: ", if (length(severe_donor_dominant_clusters)) {
    collapse_unique(severe_donor_dominant_clusters)
} else {
    "NONE"
}), paste0("- Moderate donor-enriched clusters: ", if (length(moderate_donor_enriched_clusters)) {
    collapse_unique(moderate_donor_enriched_clusters)
} else {
    "NONE"
}), paste0("- Rare clusters: ", if (length(rare_clusters)) {
    collapse_unique(rare_clusters)
} else {
    "NONE"
}), paste0("- Resolution ARI 0.2 vs 0.4: ", round(ARI_02_04, 4)), paste0("- Resolution ARI 0.4 vs 0.6: ", round(ARI_04_06, 
    4)), paste0("- Marker analysis authorized: ", marker_analysis_authorized), "- Cluster removal authorized: FALSE", 
    "- Biological annotation locked: FALSE", "- Disease pseudobulk testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "G183_DETAILED_QC_CARD.md"))

interpretation_lines <- c("# GSE183837 detailed integration-QC interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), "", "## Fixed interpretation rules", "", "- Cluster 13 is retained but must not support a disease claim unless its identity and signal are reproduced across multiple donors.", 
    "- Clusters 7 and 8 are retained for marker extraction and reviewed for donor-associated states.", "- Cluster 15 is retained as a rare candidate population because rarity alone is not a deletion criterion.", 
    "- Low RNA complexity alone is not treated as low quality when mitochondrial and doublet indicators are not concordantly abnormal.", 
    "- Control-versus-RIF visualization is descriptive only; cell-level hypothesis testing remains prohibited.", 
    "", "## Authorized next operation", "", "- Run memory-safe descriptive cluster marker extraction.", "- Review canonical lineage markers and the marker profiles of clusters 7, 8, 13 and 15.", 
    "- Build a candidate broad-cell annotation map.", "", "## Not authorized yet", "", "- No cluster removal.", 
    "- No biological annotation lock.", "- No disease pseudobulk testing.", "- No result lock.")

write_text(interpretation_lines, file.path(run_qc_root, "G183_DETAILED_QC_INTERPRETATION.md"))

QC_notes <- c("# STEP04E1K QC notes", "", paste0("Gate: ", gate_status), "", "The integration checkpoint itself is unchanged.", 
    "All review clusters remain in the object.", "Marker extraction is descriptive and cluster-wise.", "Disease inference remains donor-level and is deferred until annotation QC.", 
    "", "Known pre-run composition from STEP04E1J:", "- Cluster 13: 402 cells; strongly dominated by GSM5572242.", 
    "- Clusters 7 and 8: donor-enriched, particularly GSM5572241.", "- Cluster 15: 87 cells distributed across all nine donors.", 
    "", paste0("Detailed review flags: ", if (length(overall_review_flags)) {
        collapse_unique(overall_review_flags)
    } else {
        "NONE"
    }))

write_text(QC_notes, file.path(run_qc_root, "QC_NOTES.md"))

parameter_manifest <- list(script_name = script_name, script_version = script_version, supersedes_partial_QC_versions = c("INTQC183_v1_0_CELL_N_COLUMN_COLLISION", 
    "INTQC183_v1_1_WINDOWS_PATH_LIMIT"), merge_policy = "COMPOSITION_CELL_N_AND_TECHNICAL_CELL_N_EXPLICITLY_DISTINCT", 
    path_policy = "SHORT_PERSISTENT_DIRECTORIES_AND_FILENAMES_MAX_235_CHARACTERS", dataset_id = dataset_id, integration_version = integration_version, 
    QC_version = QC_version, upstream_gate = upstream_gate$gate_status[[1]], expected_cell_n = 47161L, expected_donor_n = 9L, 
    expected_primary_cluster_n = 16L, graph_preference = c("harmony_nn", "harmony_snn"), donor_dominance_review_threshold = 0.5, 
    donor_dominance_severe_threshold = 0.75, donor_dominance_hold_threshold = 0.94999999999999996, rare_cluster_review_threshold = 100L, 
    minimum_cluster_hold_threshold = 30L, donor_entropy_review_thresholds = c(0.5, 0.69999999999999996), local_mixing_review_rules = list(median_same_sample_excess = 0.25, 
        Q90_same_sample_excess = 0.5, median_neighbor_donor_n = 3), resolution_ARI_review_threshold = 0.59999999999999998, 
    marker_analysis_policy = "DESCRIPTIVE_CLUSTER_MARKERS_ONLY", cluster_removal_authorized = FALSE, disease_cell_level_testing_authorized = FALSE, 
    result_lock_authorized = FALSE, random_seed = 20260806L)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, supersedes_partial_QC_versions = "INTQC183_v1_0_CELL_N_COLUMN_COLLISION;INTQC183_v1_1_WINDOWS_PATH_LIMIT", 
    dataset_id = dataset_id, integration_version = integration_version, QC_version = QC_version, started_at = as.character(run_start), 
    completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 
        2), gate_status = gate_status, cell_n = nrow(metadata), donor_n = length(unique(metadata$donor_id)), primary_cluster_n = length(unique(metadata$seurat_clusters)), 
    severe_donor_dominant_clusters = collapse_unique(severe_donor_dominant_clusters), moderate_donor_enriched_clusters = collapse_unique(moderate_donor_enriched_clusters), 
    rare_clusters = collapse_unique(rare_clusters), marker_analysis_authorized = marker_analysis_authorized, cluster_removal_authorized = FALSE, 
    result_lock_authorized = FALSE, free_memory_before_GB = free_memory_before_GB, free_memory_after_GB = free_memory_after_GB, 
    free_disk_before_GB = free_disk_before_GB, free_disk_after_GB = free_disk_after_GB, output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

persistent_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(source_data_root, 
    recursive = TRUE, full.names = TRUE), list.files(figure_output_root, recursive = TRUE, full.names = TRUE)))

persistent_files <- persistent_files[file.exists(persistent_files) & !dir.exists(persistent_files)]

persistent_manifest <- data.frame(path = normalize_slashes(persistent_files), size_bytes = file.info(persistent_files)$size, 
    md5 = unname(tools::md5sum(persistent_files)), stringsAsFactors = FALSE)

write_csv_utf8(persistent_manifest, file.path(manifest_root, "persistent_manifest.csv"))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(figure_output_root, 
    recursive = TRUE, full.names = TRUE), list.files(manifest_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, normalizePath(.vdr_script_path(script_name), 
    winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2]

artifact_manifest <- data.frame(path = normalize_slashes(small_files), size_bytes = file.info(small_files)$size, 
    md5 = unname(tools::md5sum(small_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "marker_analysis_authorized=", marker_analysis_authorized, 
    "review_flags=", if (length(overall_review_flags)) {
        collapse_unique(overall_review_flags)
    } else {
        "NONE"
    })

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

cat("STEP04E1K GSE183837 detailed integration QC completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Cells: ", format(nrow(metadata), big.mark = ","), "\n", sep = "")

cat("Donors: ", length(unique(metadata$donor_id)), "\n", sep = "")

cat("Primary clusters: ", length(unique(metadata$seurat_clusters)), "\n", sep = "")

cat("Severe donor-dominant clusters: ", if (length(severe_donor_dominant_clusters)) {
    collapse_unique(severe_donor_dominant_clusters)
} else {
    "NONE"
}, "\n", sep = "")

cat("Moderate donor-enriched clusters: ", if (length(moderate_donor_enriched_clusters)) {
    collapse_unique(moderate_donor_enriched_clusters)
} else {
    "NONE"
}, "\n", sep = "")

cat("Rare clusters: ", if (length(rare_clusters)) {
    collapse_unique(rare_clusters)
} else {
    "NONE"
}, "\n", sep = "")

cat("ARI 0.2 vs 0.4: ", round(ARI_02_04, 4), "\n", sep = "")

cat("ARI 0.4 vs 0.6: ", round(ARI_04_06, 4), "\n", sep = "")

cat("Marker analysis authorized: ", marker_analysis_authorized, "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe detailed QC outputs are complete in the project folders.\n", 
        sep = "")
}

