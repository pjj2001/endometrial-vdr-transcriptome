#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C03: Assess integrated cell and donor representation

options(stringsAsFactors = FALSE, warn = 1)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "assess_gse250130_integration.R"

script_version <- "STEP04E1C_GSE250130_INTEGRATION_QC_v1.1_FIXED"

input_version <- "SCI_v1_1_memory_safe"

qc_version <- "GSE250130_INTEGRATION_QC_v1.1_FIXED"

set.seed(20260805)

run_start <- Sys.time()

run_id <- paste0("STEP04E1C_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

normalized_entropy <- function(counts) {
    counts <- counts[counts > 0]
    if (length(counts) <= 1L) 
        return(0)
    p <- counts/sum(counts)
    -sum(p * log(p))/log(length(p))
}

effective_sample_n <- function(counts) {
    counts <- counts[counts > 0]
    if (!length(counts)) 
        return(0)
    p <- counts/sum(counts)
    exp(-sum(p * log(p)))
}

detect_first_column <- function(column_names, candidates) {
    lower <- tolower(column_names)
    candidate_lower <- tolower(candidates)
    for (candidate in candidate_lower) {
        hit <- which(lower == candidate)
        if (length(hit)) 
            return(column_names[hit[1]])
    }
    ""
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

log_path <- NULL

log_event <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(line, "\n")
    if (!is.null(log_path)) {
        cat(line, "\n", file = log_path, append = TRUE)
    }
    flush.console()
}

required_packages <- c("data.table", "ggplot2", "jsonlite", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "), "\nInstall them first with:\ninstall.packages(c(", 
        paste0("\"", missing_packages, "\"", collapse = ", "), "))")
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

input_root <- file.path(results_root, "S4E1_sc", input_version)

input_table_root <- file.path(input_root, "tables", dataset_id)

input_checkpoint_root <- file.path(input_root, "checkpoints", dataset_id)

locked_design_path <- file.path(results_root, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0", "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv")

doublet_qc_path <- file.path(results_root, "S4E1_sc", "SCI_v1_0", "tables", dataset_id, paste0(dataset_id, "_sample_doublet_qc.csv"))

umap_path <- file.path(input_table_root, paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv"))

cluster_sample_path <- file.path(input_table_root, paste0(dataset_id, "_cluster_by_sample_prejoin.csv"))

checkpoint_status_path <- file.path(input_checkpoint_root, paste0(dataset_id, "_checkpoint_status.csv"))

run_root <- file.path(qc_base_root, paste0(run_id, "_integration_QC"))

table_root <- file.path(run_root, "tables")

figure_root <- file.path(run_root, "figures")

source_data_root <- file.path(run_root, "source_data")

for (p in c(execution_root, parameter_root, environment_root, run_root, table_root, figure_root, source_data_root)) {
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

required_inputs <- c(umap_path, cluster_sample_path, checkpoint_status_path, locked_design_path, doublet_qc_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

log_event("Starting", run_id)

umap <- .vdr_fread(umap_path, data.table = FALSE, showProgress = FALSE)

cluster_sample <- .vdr_fread(cluster_sample_path, data.table = FALSE, showProgress = FALSE)

checkpoint_status <- safe_read_csv(checkpoint_status_path)

locked_design <- safe_read_csv(locked_design_path)

doublet_qc <- safe_read_csv(doublet_qc_path)

required_umap_columns <- c("cell_barcode", "UMAPH_1", "UMAPH_2", "sample_accession", "seurat_clusters", "RNA_snn_res.0.2", 
    "RNA_snn_res.0.4", "RNA_snn_res.0.6")

missing_umap_columns <- setdiff(required_umap_columns, colnames(umap))

if (length(missing_umap_columns)) {
    stop("UMAP metadata is missing columns: ", paste(missing_umap_columns, collapse = ", "))
}

required_cluster_sample_columns <- c("sample_accession", "seurat_clusters", "cell_n")

missing_cluster_sample_columns <- setdiff(required_cluster_sample_columns, colnames(cluster_sample))

if (length(missing_cluster_sample_columns)) {
    stop("Cluster-by-sample table is missing columns: ", paste(missing_cluster_sample_columns, collapse = ", "))
}

if (is.null(checkpoint_status) || nrow(checkpoint_status) != 1L || checkpoint_status$checkpoint_status[[1]] != 
    "COMPLETE") {
    stop("Integrated checkpoint status is absent or not COMPLETE.")
}

if (is.null(locked_design) || !"sample_accession" %in% colnames(locked_design)) {
    stop("Locked GSE250130 design is absent or lacks sample_accession.")
}

if (is.null(doublet_qc) || !"sample_accession" %in% colnames(doublet_qc)) {
    stop("Doublet QC table is absent or lacks sample_accession.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, qc_version = qc_version, 
    input_version = input_version, dataset_id = dataset_id, expected_sample_n = 28L, expected_cell_n = 236636L, 
    expected_feature_n = 36601L, selected_resolution = "RNA_snn_res.0.4", expected_cluster_n_res_0_4 = 22L, rare_cluster_cell_threshold = 200L, 
    low_sample_coverage_threshold = 14L, moderate_sample_dominance_threshold = 0.29999999999999999, severe_sample_dominance_threshold = 0.5, 
    low_normalized_entropy_threshold = 0.65000000000000002, plot_max_cells_per_cluster = 5000L, random_seed = 20260805L, 
    output_status = "QC_ONLY_NOT_LOCKED")

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

umap$sample_accession <- as.character(umap$sample_accession)

umap$seurat_clusters <- as.character(umap$seurat_clusters)

cluster_sample$sample_accession <- as.character(cluster_sample$sample_accession)

cluster_sample$seurat_clusters <- as.character(cluster_sample$seurat_clusters)

locked_design$sample_accession <- as.character(locked_design$sample_accession)

design_duplicate_n <- sum(duplicated(locked_design$sample_accession))

if (design_duplicate_n > 0L) {
    stop("Locked design contains duplicated sample_accession rows: ", design_duplicate_n)
}

design_columns <- setdiff(colnames(locked_design), "sample_accession")

group_column <- detect_first_column(colnames(locked_design), c("group", "clinical_group", "phenotype", "disease_status", 
    "condition", "fertility_status", "patient_group", "case_control_status", "case_control"))

time_column <- detect_first_column(colnames(locked_design), c("timepoint", "time_point", "temporal_label", "cycle_day", 
    "lh_day", "day", "phase", "sampling_time"))

role_column <- detect_first_column(colnames(locked_design), c("analysis_role", "cohort_role", "sample_role", "role"))

donor_column <- detect_first_column(colnames(locked_design), c("donor_id", "subject_id", "patient_id", "individual_id"))

design_column_inventory <- data.frame(column = colnames(locked_design), detected_role = ifelse(colnames(locked_design) == 
    group_column, "GROUP", ifelse(colnames(locked_design) == time_column, "TIME", ifelse(colnames(locked_design) == 
    role_column, "ANALYSIS_ROLE", ifelse(colnames(locked_design) == donor_column, "DONOR", "OTHER")))), stringsAsFactors = FALSE)

write_csv_utf8(design_column_inventory, file.path(table_root, "01_LOCKED_DESIGN_COLUMN_INVENTORY.csv"))

umap_design <- merge(umap, locked_design, by = "sample_accession", all.x = TRUE, sort = FALSE)

unmatched_umap_sample_n <- length(setdiff(unique(umap$sample_accession), locked_design$sample_accession))

cluster_sample_design <- merge(cluster_sample, locked_design, by = "sample_accession", all.x = TRUE, sort = FALSE)

resolution_columns <- c("RNA_snn_res.0.2", "RNA_snn_res.0.4", "RNA_snn_res.0.6")

resolution_summary <- do.call(rbind, lapply(resolution_columns, function(column_name) {
    tab <- table(umap[[column_name]])
    data.frame(resolution_column = column_name, cluster_n = length(tab), min_cluster_cell_n = min(as.numeric(tab)), 
        median_cluster_cell_n = stats::median(as.numeric(tab)), max_cluster_cell_n = max(as.numeric(tab)), clusters_below_200_cells = sum(as.numeric(tab) < 
            200L), stringsAsFactors = FALSE)
}))

write_csv_utf8(resolution_summary, file.path(table_root, "02_RESOLUTION_SUMMARY.csv"))

finite_umap_fraction <- mean(is.finite(umap$UMAPH_1) & is.finite(umap$UMAPH_2))

global_qc <- data.frame(metric = c("checkpoint_status", "checkpoint_join_policy", "cell_n", "sample_n", "feature_n", 
    "cluster_n_resolution_0_4", "finite_UMAP_fraction", "locked_design_sample_n", "unmatched_UMAP_sample_n", "doublet_QC_sample_n", 
    "design_duplicate_sample_n", "group_column_detected", "time_column_detected", "role_column_detected", "donor_column_detected"), 
    observed = c(checkpoint_status$checkpoint_status[[1]], checkpoint_status$join_policy[[1]], nrow(umap), length(unique(umap$sample_accession)), 
        checkpoint_status$feature_n[[1]], length(unique(umap$seurat_clusters)), finite_umap_fraction, nrow(locked_design), 
        unmatched_umap_sample_n, nrow(doublet_qc), design_duplicate_n, group_column, time_column, role_column, 
        donor_column), expected_or_interpretation = c("COMPLETE", "DATA_ONLY_COUNTS_REMAIN_SPLIT", "236636", "28", 
        "36601", "22", "1.0", "28", "0", "28", "0", "Detected when available", "Detected when available", "Detected when available", 
        "Detected when available"), stringsAsFactors = FALSE)

write_csv_utf8(global_qc, file.path(table_root, "03_GLOBAL_INTEGRATION_QC.csv"))

cluster_sample_dt <- data.table::as.data.table(cluster_sample)

cluster_sample_dt[, `:=`(cell_n, as.numeric(cell_n))]

cluster_mixing <- cluster_sample_dt[, {
    counts <- as.numeric(cell_n)
    total <- as.numeric(sum(counts, na.rm = TRUE))
    proportions <- counts/total
    max_index <- which.max(proportions)
    list(cell_n = as.numeric(total), sample_n_present = as.integer(data.table::uniqueN(sample_accession)), min_sample_cell_n = as.numeric(min(counts, 
        na.rm = TRUE)), median_sample_cell_n = as.numeric(stats::median(counts, na.rm = TRUE)), max_sample_cell_n = as.numeric(max(counts, 
        na.rm = TRUE)), max_sample_accession = as.character(sample_accession[max_index]), max_sample_fraction = as.numeric(proportions[max_index]), 
        normalized_sample_entropy = as.numeric(normalized_entropy(counts)), effective_sample_n = as.numeric(effective_sample_n(counts)))
}, by = seurat_clusters]

cluster_mixing$cluster_numeric <- suppressWarnings(as.integer(cluster_mixing$seurat_clusters))

cluster_mixing <- cluster_mixing[order(cluster_mixing$cluster_numeric), ]

cluster_mixing$rare_cluster_flag <- (cluster_mixing$cell_n < 200L)

cluster_mixing$low_sample_coverage_flag <- (cluster_mixing$sample_n_present < 14L)

cluster_mixing$moderate_sample_dominance_flag <- (cluster_mixing$max_sample_fraction > 0.29999999999999999)

cluster_mixing$severe_sample_dominance_flag <- (cluster_mixing$max_sample_fraction > 0.5)

cluster_mixing$low_entropy_flag <- (cluster_mixing$normalized_sample_entropy < 0.65000000000000002)

cluster_mixing$review_flag <- (cluster_mixing$rare_cluster_flag | cluster_mixing$low_sample_coverage_flag | cluster_mixing$severe_sample_dominance_flag | 
    cluster_mixing$low_entropy_flag)

cluster_mixing$review_reason <- vapply(seq_len(nrow(cluster_mixing)), function(i) {
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
}, character(1))

write_csv_utf8(as.data.frame(cluster_mixing), file.path(table_root, "04_CLUSTER_SAMPLE_MIXING_METRICS.csv"))

cluster_total_lookup <- setNames(cluster_mixing$cell_n, cluster_mixing$seurat_clusters)

cluster_sample$cluster_fraction <- (cluster_sample$cell_n/cluster_total_lookup[cluster_sample$seurat_clusters])

sample_total_lookup <- tapply(cluster_sample$cell_n, cluster_sample$sample_accession, sum)

cluster_sample$within_sample_fraction <- (cluster_sample$cell_n/sample_total_lookup[cluster_sample$sample_accession])

write_csv_utf8(cluster_sample, file.path(source_data_root, "SOURCE_DATA_cluster_by_sample.csv"))

group_time_outputs <- list()

summarize_design_variable <- function(variable_name, output_stem) {
    if (!nzchar(variable_name)) 
        return(NULL)
    dt <- data.table::as.data.table(cluster_sample_design)
    dt[, `:=`(design_value, as.character(get(variable_name)))]
    out <- dt[, .(cell_n = as.numeric(sum(as.numeric(cell_n), na.rm = TRUE)), sample_n = as.integer(data.table::uniqueN(sample_accession))), 
        by = .(seurat_clusters, design_value)]
    out[, `:=`(cluster_fraction, cell_n/sum(cell_n)), by = seurat_clusters]
    write_csv_utf8(as.data.frame(out), file.path(table_root, paste0(output_stem, ".csv")))
    out
}

group_composition <- summarize_design_variable(group_column, "05_CLUSTER_BY_GROUP_COMPOSITION")

time_composition <- summarize_design_variable(time_column, "06_CLUSTER_BY_TIME_COMPOSITION")

role_composition <- summarize_design_variable(role_column, "07_CLUSTER_BY_ANALYSIS_ROLE_COMPOSITION")

theme_qc <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95"), 
    legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8))

umap_dt <- data.table::as.data.table(umap_design)

umap_plot_data <- umap_dt[, {
    n_keep <- min(.N, 5000L)
    .SD[sample.int(.N, n_keep)]
}, by = seurat_clusters]

umap_plot_data$seurat_clusters <- factor(umap_plot_data$seurat_clusters, levels = as.character(sort(unique(suppressWarnings(as.integer(umap_plot_data$seurat_clusters))))))

cluster_centroids <- umap_dt[, .(UMAPH_1 = stats::median(UMAPH_1), UMAPH_2 = stats::median(UMAPH_2)), by = seurat_clusters]

plot_umap_cluster <- ggplot2::ggplot(umap_plot_data, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = seurat_clusters)) + 
    ggplot2::geom_point(size = 0.080000000000000002, alpha = 0.55000000000000004, stroke = 0) + ggplot2::geom_text(data = cluster_centroids, 
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, label = seurat_clusters), inherit.aes = FALSE, size = 3, fontface = "bold") + 
    ggplot2::labs(title = "GSE250130 integration QC: Harmony UMAP", subtitle = "Resolution 0.4; deterministic plotting subset", 
        color = "Cluster", x = "UMAP 1", y = "UMAP 2") + theme_qc + ggplot2::theme(legend.position = "right")

save_plot_pair(plot_umap_cluster, "QC_Figure_1_UMAP_clusters_res0_4", width = 8.5, height = 7)

plot_cluster_size <- ggplot2::ggplot(as.data.frame(cluster_mixing), ggplot2::aes(x = reorder(factor(seurat_clusters), 
    cluster_numeric), y = cell_n, fill = review_flag)) + ggplot2::geom_col(width = 0.80000000000000004) + ggplot2::scale_y_continuous(labels = scales::comma) + 
    ggplot2::labs(title = "Cluster sizes at resolution 0.4", x = "Cluster", y = "Cell count", fill = "Review flag") + 
    theme_qc

save_plot_pair(plot_cluster_size, "QC_Figure_2_cluster_sizes", width = 8, height = 5)

plot_mixing <- ggplot2::ggplot(as.data.frame(cluster_mixing), ggplot2::aes(x = max_sample_fraction, y = normalized_sample_entropy, 
    size = cell_n, label = seurat_clusters)) + ggplot2::geom_vline(xintercept = c(0.29999999999999999, 0.5), linetype = c("dashed", 
    "solid")) + ggplot2::geom_hline(yintercept = 0.65000000000000002, linetype = "dashed") + ggplot2::geom_point(ggplot2::aes(shape = review_flag), 
    alpha = 0.80000000000000004) + ggplot2::geom_text(nudge_y = 0.025000000000000001, size = 3) + ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::scale_size_continuous(labels = scales::comma) + ggplot2::labs(title = "Cluster sample-mixing QC", 
    subtitle = "Lower dominance and higher entropy indicate broader sample mixing", x = "Largest single-sample fraction", 
    y = "Normalized sample entropy", size = "Cell count", shape = "Review flag") + theme_qc

save_plot_pair(plot_mixing, "QC_Figure_3_cluster_sample_mixing", width = 8, height = 6)

heatmap_data <- cluster_sample

heatmap_data$seurat_clusters <- factor(heatmap_data$seurat_clusters, levels = as.character(sort(unique(suppressWarnings(as.integer(heatmap_data$seurat_clusters))))))

sample_order <- names(sort(tapply(heatmap_data$cell_n, heatmap_data$sample_accession, sum), decreasing = TRUE))

heatmap_data$sample_accession <- factor(heatmap_data$sample_accession, levels = rev(sample_order))

plot_heatmap <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = seurat_clusters, y = sample_accession, fill = within_sample_fraction)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(labels = scales::percent_format(accuracy = 1), option = "C") + 
    ggplot2::labs(title = "Within-sample cluster composition", x = "Cluster", y = "Sample", fill = "Fraction") + 
    theme_qc + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))

save_plot_pair(plot_heatmap, "QC_Figure_4_sample_cluster_composition_heatmap", width = 9, height = 8)

resolution_long <- do.call(rbind, lapply(resolution_columns, function(column_name) {
    tab <- as.data.frame(table(umap[[column_name]]))
    colnames(tab) <- c("cluster", "cell_n")
    tab$resolution <- sub("RNA_snn_res.", "", column_name, fixed = TRUE)
    tab
}))

resolution_long$cluster <- as.character(resolution_long$cluster)

write_csv_utf8(resolution_long, file.path(source_data_root, "SOURCE_DATA_resolution_cluster_sizes.csv"))

plot_resolution <- ggplot2::ggplot(resolution_long, ggplot2::aes(x = reorder(cluster, cell_n), y = cell_n)) + ggplot2::geom_col() + 
    ggplot2::facet_wrap(~resolution, scales = "free_x", nrow = 1) + ggplot2::scale_y_continuous(labels = scales::comma) + 
    ggplot2::labs(title = "Cluster-size distribution across tested resolutions", x = "Cluster", y = "Cell count") + 
    theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))

save_plot_pair(plot_resolution, "QC_Figure_5_resolution_cluster_sizes", width = 11, height = 5)

technical_checks <- data.frame(check = c("checkpoint_complete", "join_policy_correct", "cell_count_locked", "sample_count_locked", 
    "feature_count_locked", "cluster_count_res0_4_expected", "all_UMAP_coordinates_finite", "all_UMAP_samples_in_locked_design", 
    "locked_design_has_28_samples", "doublet_QC_has_28_samples"), pass = c(checkpoint_status$checkpoint_status[[1]] == 
    "COMPLETE", checkpoint_status$join_policy[[1]] == "DATA_ONLY_COUNTS_REMAIN_SPLIT", nrow(umap) == 236636L, length(unique(umap$sample_accession)) == 
    28L, checkpoint_status$feature_n[[1]] == 36601L, length(unique(umap$seurat_clusters)) == 22L, identical(finite_umap_fraction, 
    1), unmatched_umap_sample_n == 0L, nrow(locked_design) == 28L, nrow(doublet_qc) == 28L), stringsAsFactors = FALSE)

write_csv_utf8(technical_checks, file.path(table_root, "08_TECHNICAL_GATE_CHECKS.csv"))

technical_pass <- all(technical_checks$pass)

review_clusters <- cluster_mixing$seurat_clusters[cluster_mixing$review_flag]

gate_status <- if (!technical_pass) {
    "HOLD_TECHNICAL_INTEGRATION_QC_FAILED"
} else if (length(review_clusters)) {
    "PASS_WITH_REVIEW_FLAGS"
} else {
    "PASS"
}

gate_table <- data.frame(dataset_id = dataset_id, qc_version = qc_version, input_version = input_version, gate_status = gate_status, 
    technical_check_n = nrow(technical_checks), technical_pass_n = sum(technical_checks$pass), technical_fail_n = sum(!technical_checks$pass), 
    cluster_n = nrow(cluster_mixing), review_cluster_n = length(review_clusters), review_clusters = collapse_unique(review_clusters), 
    group_column = group_column, time_column = time_column, role_column = role_column, donor_column = donor_column, 
    marker_analysis_authorized = technical_pass, biological_annotation_locked = FALSE, result_lock_authorized = FALSE, 
    stringsAsFactors = FALSE)

gate_path <- file.path(run_root, "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv")

write_csv_utf8(gate_table, gate_path)

review_detail <- cluster_mixing[cluster_mixing$review_flag, c("seurat_clusters", "cell_n", "sample_n_present", 
    "max_sample_accession", "max_sample_fraction", "normalized_sample_entropy", "review_reason")]

interpretation_lines <- c("# STEP04E1C GSE250130 integration QC", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: **", 
    gate_status, "**"), paste0("- Cells: ", format(nrow(umap), big.mark = ",")), paste0("- Samples: ", length(unique(umap$sample_accession))), 
    paste0("- Resolution 0.4 clusters: ", length(unique(umap$seurat_clusters))), paste0("- Finite UMAP fraction: ", 
        finite_umap_fraction), paste0("- Review clusters: ", ifelse(length(review_clusters), collapse_unique(review_clusters), 
        "none")), "", "## Interpretation", "", if (technical_pass) {
        "- The GSE250130 integrated object passes the predefined technical identity and completeness checks."
    } else {
        "- One or more technical identity/completeness checks failed; marker analysis must not proceed."
    }, if (length(review_clusters)) {
        "- Review flags identify clusters that are rare, present in few samples, dominated by one sample, or show low normalized sample entropy."
    } else {
        "- No cluster crossed the predefined review thresholds."
    }, "- A review flag does not prove a batch artifact. These clusters must be evaluated using marker genes, sample metadata, and biological plausibility.", 
    "- The current Gate authorizes marker analysis when technical checks pass, but it does not authorize final annotation or result locking.", 
    "", "## Required next analysis", "", "1. Compute cluster markers using the joined normalized data layer.", 
    "2. Apply marker-based broad-cell annotation with explicit marker coverage.", "3. Review flagged clusters before accepting labels.", 
    "4. Only after annotation QC, calculate VDR program localization and donor-aware pseudobulk.")

interpretation_path <- file.path(run_root, "STEP04E1C_GSE250130_QC_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

card_lines <- c("# STEP04E1C Result Card", "", paste0("- Dataset: ", dataset_id), paste0("- Input: ", input_version), 
    paste0("- QC version: ", qc_version), paste0("- Status: ", gate_status), paste0("- Marker analysis authorized: ", 
        technical_pass), "- Biological annotation locked: FALSE", "- Result lock authorized: FALSE", paste0("- Review clusters: ", 
        ifelse(length(review_clusters), collapse_unique(review_clusters), "none")), "", "This card documents integration QC only. It is not a scientific result lock.")

card_path <- file.path(run_root, "STEP04E1C_GSE250130_RESULT_CARD.md")

write_text(card_lines, card_path)

run_end <- Sys.time()

package_table <- data.frame(package = c("data.table", "ggplot2", "jsonlite", "scales"), version = vapply(c("data.table", 
    "ggplot2", "jsonlite", "scales"), safe_package_version, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- list.files(run_root, recursive = TRUE, full.names = TRUE)

artifact_manifest <- data.frame(relative_path = normalize_slashes(substring(artifact_files, nchar(run_root) + 2L)), 
    size_bytes = file.info(artifact_files)$size, md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, qc_version = qc_version, 
    input_version = input_version, dataset_id = dataset_id, started_at = as.character(run_start), completed_at = as.character(run_end), 
    elapsed_seconds = round(as.numeric(difftime(run_end, run_start, units = "secs")), 2), gate_status = gate_status, 
    large_Seurat_checkpoint_loaded = FALSE, project_analysis_modified = FALSE, marker_analysis_authorized = technical_pass, 
    result_lock_authorized = FALSE, output_root = normalize_slashes(run_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

small_handoff_root <- file.path(run_root, "small_handoff")

dir.create(small_handoff_root, recursive = TRUE, showWarnings = FALSE)

small_files <- unique(c(list.files(run_root, recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, 
    parameter_path, session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2 & !startsWith(normalize_slashes(small_files), normalize_slashes(small_handoff_root))]

copy_manifest <- list()

for (i in seq_along(small_files)) {
    source_path <- small_files[i]
    destination_path <- file.path(small_handoff_root, paste0(sprintf("%03d", i), "_", basename(source_path)))
    copied <- file.copy(source_path, destination_path, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    copy_manifest[[i]] <- data.frame(source = normalize_slashes(source_path), destination = normalize_slashes(destination_path), 
        copied = isTRUE(copied), stringsAsFactors = FALSE)
}

if (length(copy_manifest)) {
    write_csv_utf8(do.call(rbind, copy_manifest), file.path(run_root, "SMALL_HANDOFF_COPY_MANIFEST.csv"))
}

if (!requireNamespace("zip", quietly = TRUE)) {
    try(install.packages("zip", repos = "https://cloud.r-project.org"), silent = TRUE)
}

zip_ok <- FALSE

zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
    zip_ok <- tryCatch({
        zip::zipr(zipfile = handoff_zip, files = c(gate_path, card_path, interpretation_path, table_root, figure_root, 
            source_data_root, log_path, run_manifest_path, parameter_path, session_path, file.path(run_root, "PACKAGE_VERSIONS.csv"), 
            file.path(run_root, "ARTIFACT_MANIFEST.csv"), file.path(run_root, "SMALL_HANDOFF_COPY_MANIFEST.csv")), 
            recurse = TRUE, include_directories = TRUE)
        file.exists(handoff_zip)
    }, error = function(e) {
        zip_error <<- conditionMessage(e)
        FALSE
    })
}

log_event("Completed", run_id, "gate_status=", gate_status, "review_clusters=", collapse_unique(review_clusters))

cat("\n============================================================\n")

cat("STEP04E1C GSE250130 integration QC completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("QC folder: ", run_root, "\n", sep = "")

cat("Run manifest: ", run_manifest_path, "\n", sep = "")

cat("Review clusters: ", ifelse(length(review_clusters), collapse_unique(review_clusters), "none"), "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe QC folder is complete; compress it manually if needed.\n", sep = "")
}

