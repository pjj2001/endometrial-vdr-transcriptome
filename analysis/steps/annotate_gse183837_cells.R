#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D05: Apply accepted broad-cell annotations

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "annotate_gse183837_cells.R"

script_version <- "STEP04E1M_GSE183837_BROAD_ANNOTATION_QC_v1.0"

integration_version <- "SCI183_v1_0"

marker_version <- "MK183_v1_0"

annotation_version <- "ANN183_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1M_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

split_semicolon <- function(x) {
    output <- unlist(strsplit(as.character(x), ";", fixed = TRUE), use.names = FALSE)
    trimws(output[nzchar(trimws(output))])
}

safe_zscore <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 2L || !is.finite(stats::sd(x)) || stats::sd(x) == 0) {
        return(rep(0, length(x)))
    }
    as.numeric(scale(x))
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

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

integration_root <- file.path(results_root, "S4E1_sc", integration_version)

input_checkpoint <- file.path(integration_root, "checkpoints", dataset_id, "integrated.rds")

marker_root <- file.path(integration_root, "mk", "G183", marker_version)

marker_summary_path <- file.path(marker_root, "02_marker_summary.csv")

all_top100_path <- file.path(marker_root, "top", "all_top100.csv")

donor_support_path <- file.path(marker_root, "top", "donor_support_summary.csv")

candidate_annotation_path <- file.path(marker_root, "panel", "candidate_annotation.csv")

detailed_qc_table_path <- file.path(integration_root, "iqc", "G183", "IQ183_v1_2", "t", "04_cluster_QC.csv")

output_root <- file.path(integration_root, "ann", "G183", annotation_version)

table_root <- file.path(output_root, "t")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "man")

annotated_checkpoint <- file.path(output_root, "annotated.rds")

annotation_status_path <- file.path(output_root, "status.csv")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1M_", format(run_start, "%Y%m%d_%H%M%S"), "_annQC"))

figure_root <- file.path(run_qc_root, "fig")

for (path_value in c(execution_root, parameter_root, environment_root, output_root, table_root, figure_output_root, 
    manifest_root, run_qc_root, figure_root)) {
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

planned_paths <- data.frame(artifact = c("input_checkpoint", "marker_summary", "all_top100", "donor_support", "candidate_annotation", 
    "detailed_cluster_QC", "output_root", "annotated_checkpoint", "annotation_status", "run_qc_root", "handoff_zip"), 
    path = normalize_slashes(c(input_checkpoint, marker_summary_path, all_top100_path, donor_support_path, candidate_annotation_path, 
        detailed_qc_table_path, output_root, annotated_checkpoint, annotation_status_path, run_qc_root, handoff_zip)), 
    stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

if (!all(planned_paths$path_safe)) {
    stop("One or more annotation paths exceed the safe limit.")
}

marker_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

marker_qc_dirs <- marker_qc_dirs[grepl("^STEP04E1L_[0-9]{8}_[0-9]{6}_markers$", basename(marker_qc_dirs), perl = TRUE)]

if (!length(marker_qc_dirs)) {
    stop("No STEP04E1L marker-QC directory was found.")
}

gate_candidates <- list()

for (directory_path in marker_qc_dirs) {
    gate_path <- file.path(directory_path, "G183_MARKER_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        gate_candidates[[length(gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(gate_candidates)) {
    stop("No readable STEP04E1L marker Gate was found.")
}

gate_candidates <- as.data.frame(data.table::rbindlist(gate_candidates, use.names = TRUE, fill = TRUE))

gate_candidates$source_mtime_POSIX <- as.POSIXct(gate_candidates$source_mtime)

gate_candidates <- gate_candidates[order(gate_candidates$source_mtime_POSIX, decreasing = TRUE), , drop = FALSE]

authorized_index <- which(gate_candidates$gate_status == "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED" & 
    toupper(as.character(gate_candidates$candidate_annotation_review_authorized)) %in% c("TRUE", "T", "1", "YES"))

if (!length(authorized_index)) {
    stop("No marker Gate authorizes annotation review.")
}

upstream_gate <- gate_candidates[authorized_index[1], , drop = FALSE]

write_csv_utf8(upstream_gate, file.path(run_qc_root, "00_upstream_marker_gate.csv"))

required_input_paths <- c(input_checkpoint, marker_summary_path, all_top100_path, donor_support_path, candidate_annotation_path, 
    detailed_qc_table_path)

missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]

if (length(missing_input_paths)) {
    stop("Missing required input file(s):\n", paste(normalize_slashes(missing_input_paths), collapse = "\n"))
}

marker_summary <- safe_read_csv(marker_summary_path)

all_top100 <- safe_read_csv(all_top100_path)

donor_support <- safe_read_csv(donor_support_path)

algorithmic_candidates <- safe_read_csv(candidate_annotation_path)

detailed_cluster_QC <- safe_read_csv(detailed_qc_table_path)

integrated_object <- readRDS(input_checkpoint)

if (!inherits(integrated_object, "Seurat")) {
    stop("Input checkpoint is not a Seurat object.")
}

metadata <- integrated_object[[]]

if (ncol(integrated_object) != 47161L || length(unique(as.character(metadata$donor_id))) != 9L || length(unique(as.character(metadata$seurat_clusters))) != 
    16L) {
    stop("Integrated object does not match the locked 47,161-cell object.")
}

if (nrow(marker_summary) != 16L || nrow(algorithmic_candidates) != 16L || length(unique(all_top100$cluster)) != 
    16L) {
    stop("Marker outputs do not contain all 16 clusters.")
}

log_event("Starting", run_id, "cells=", ncol(integrated_object), "clusters=", length(unique(metadata$seurat_clusters)))

annotation_map <- data.frame(cluster = as.character(0:15), broad_cell_type = c("Stromal cells", "Stromal cells", 
    "Stromal cells", "Stromal cells", "Stromal cells", "Stromal cells", "NK/T cells", "Unciliated epithelial cells", 
    "Unciliated epithelial cells", "NK/T cells", "Stromal cells", "Myeloid cells", "Endothelial cells", "Stromal cells", 
    "NK/T cells", "Ciliated epithelial cells"), subtype_label = c("RARRES2_RAMP1_ribosomal_stromal", "GDF7_COL6A1_stromal", 
    "CILP_APOD_PTGDS_stromal", "HOXA13_ITGA8_ECM_stromal", "MGP_PAMR1_MMP11_stromal", "KIF20A_TOP2A_cycling_stromal", 
    "IL7R_CD3D_T", "PAX2_GABRP_SCGB_unciliated_epithelial", "SCGB1D2_low_complexity_epithelial_like", "GNLY_TRDC_KLRC1_NK_T", 
    "RGS5_MYH11_NOTCH3_perivascular_stromal", "C1QC_LYZ_FOLR2_myeloid", "ACKR1_CLDN5_EMCN_blood_endothelial", "TOP2A_MKI67_MMP11_cycling_stromal", 
    "FGFBP2_FCGR3A_cytotoxic_NK", "C20orf85_SNTN_ciliated_epithelial"), annotation_confidence = c("MODERATE", "MODERATE", 
    "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "LOW", "HIGH", "HIGH", "HIGH", "HIGH", "MODERATE", "HIGH", 
    "HIGH"), marker_evidence = c("RARRES2;RAMP1;LGALS1", "GDF7;COL6A1;COL5A1", "CILP;APOD;PTGDS", "HOXA13;ITGA8;COL6A3", 
    "MGP;PAMR1;MMP11", "KIF20A;TOP2A;MKI67", "IL7R;CD3D;TRAC", "PAX2;GABRP;KRT8;KRT18", "SCGB1D2;SCGB2A1", "GNLY;TRDC;KLRC1;NKG7", 
    "RGS5;MYH11;NOTCH3;MCAM", "C1QC;LYZ;FOLR2;CD14", "ACKR1;CLDN5;EMCN;VWF", "TOP2A;MKI67;MMP11;COL1A1", "FGFBP2;FCGR3A;PRF1;NKG7", 
    "C20orf85;SNTN;CFAP126;FOXJ1"), review_flag = c("", "", "", "", "", "", "", "MODERATE_DONOR_ENRICHMENT_IDENTITY_ROBUST", 
    "LOW_COMPLEXITY_NO_ROBUST_DONOR_MARKERS", "", "", "", "", "SEVERE_SINGLE_DONOR_DOMINANCE_NO_DISEASE_SUBTYPE_CLAIM", 
    "", "RARE_CLUSTER_DISTRIBUTED_ACROSS_DONORS"), primary_VDR_localization_eligible = c(TRUE, TRUE, TRUE, TRUE, 
    TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), primary_broad_pseudobulk_eligible = c(TRUE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), subtype_disease_inference_eligible = c(TRUE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE), annotation_locked = TRUE, 
    stringsAsFactors = FALSE)

if (nrow(annotation_map) != 16L || anyDuplicated(annotation_map$cluster) > 0L || !identical(annotation_map$cluster, 
    as.character(0:15))) {
    stop("Annotation map must contain one ordered row for every cluster 0-15.")
}

allowed_broad_types <- c("Stromal cells", "NK/T cells", "Unciliated epithelial cells", "Ciliated epithelial cells", 
    "Myeloid cells", "Endothelial cells")

if (!all(annotation_map$broad_cell_type %in% allowed_broad_types)) {
    stop("Annotation map contains an unexpected broad-cell ontology.")
}

write_csv_utf8(annotation_map, file.path(table_root, "01_annotation_map.csv"))

evidence_rows <- list()

for (row_index in seq_len(nrow(annotation_map))) {
    cluster_id <- annotation_map$cluster[row_index]
    evidence_genes <- split_semicolon(annotation_map$marker_evidence[row_index])
    cluster_top <- all_top100[as.character(all_top100$cluster) == cluster_id, , drop = FALSE]
    cluster_support <- donor_support[as.character(donor_support$cluster) == cluster_id, , drop = FALSE]
    for (gene_value in evidence_genes) {
        top_match <- cluster_top[toupper(cluster_top$gene) == toupper(gene_value), , drop = FALSE]
        support_match <- cluster_support[toupper(cluster_support$gene) == toupper(gene_value), , drop = FALSE]
        evidence_rows[[length(evidence_rows) + 1L]] <- data.frame(cluster = cluster_id, broad_cell_type = annotation_map$broad_cell_type[row_index], 
            subtype_label = annotation_map$subtype_label[row_index], annotation_confidence = annotation_map$annotation_confidence[row_index], 
            requested_marker = gene_value, present_in_top100 = nrow(top_match) > 0L, top100_rank = if (nrow(top_match)) {
                min(top_match$rank)
            }
            else {
                NA_integer_
            }, avg_log2FC = if (nrow(top_match)) {
                top_match$avg_log2FC[1]
            }
            else {
                NA_real_
            }, p_val_adj = if (nrow(top_match)) {
                top_match$p_val_adj[1]
            }
            else {
                NA_real_
            }, robust_donor_supported = if (nrow(top_match) && "robust_donor_supported" %in% colnames(top_match)) {
                isTRUE(as.logical(top_match$robust_donor_supported[1]))
            }
            else {
                FALSE
            }, supporting_donor_n = if (nrow(support_match)) {
                support_match$supporting_donor_n[1]
            }
            else {
                NA_real_
            }, supporting_donor_fraction = if (nrow(support_match)) {
                support_match$supporting_donor_fraction[1]
            }
            else {
                NA_real_
            }, stringsAsFactors = FALSE)
    }
}

marker_evidence_audit <- do.call(rbind, evidence_rows)

evidence_summary <- data.table::as.data.table(marker_evidence_audit)[, .(requested_marker_n = .N, top100_evidence_n = sum(present_in_top100, 
    na.rm = TRUE), robust_evidence_n = sum(robust_donor_supported, na.rm = TRUE), median_supporting_donor_fraction = median(supporting_donor_fraction, 
    na.rm = TRUE)), by = .(cluster, broad_cell_type, subtype_label, annotation_confidence)]

evidence_summary <- as.data.frame(evidence_summary)

evidence_summary$median_supporting_donor_fraction[!is.finite(evidence_summary$median_supporting_donor_fraction)] <- NA_real_

annotation_QC <- merge(annotation_map, evidence_summary, by = c("cluster", "broad_cell_type", "subtype_label", 
    "annotation_confidence"), all.x = TRUE, sort = FALSE)

annotation_QC <- merge(annotation_QC, marker_summary[, c("cluster", "marker_row_n", "FDR05_marker_n", "robust_donor_supported_n", 
    "top_genes", "top_robust_genes"), drop = FALSE], by = "cluster", all.x = TRUE, sort = FALSE)

annotation_QC <- merge(annotation_QC, algorithmic_candidates[, c("cluster", "candidate_broad_lineage", "second_lineage", 
    "score_margin", "candidate_support_status"), drop = FALSE], by = "cluster", all.x = TRUE, sort = FALSE)

annotation_QC <- merge(annotation_QC, detailed_cluster_QC[, c("seurat_clusters", "cell_n", "represented_donor_n", 
    "top_donor", "top_donor_fraction", "donor_entropy_normalized", "median_nFeature_RNA", "median_percent_mt", 
    "median_scDblFinder_score"), drop = FALSE], by.x = "cluster", by.y = "seurat_clusters", all.x = TRUE, sort = FALSE)

annotation_QC <- annotation_QC[match(as.character(0:15), annotation_QC$cluster), , drop = FALSE]

annotation_QC$evidence_pass <- ifelse(annotation_QC$annotation_confidence == "LOW", annotation_QC$top100_evidence_n >= 
    1L, annotation_QC$top100_evidence_n >= 2L)

write_csv_utf8(marker_evidence_audit, file.path(table_root, "02_marker_evidence_audit.csv"))

write_csv_utf8(annotation_QC, file.path(table_root, "03_annotation_QC.csv"))

metadata$primary_cluster <- as.character(metadata$seurat_clusters)

map_index <- match(metadata$primary_cluster, annotation_map$cluster)

if (anyNA(map_index)) {
    stop("One or more cells could not be matched to the annotation map.")
}

metadata$broad_cell_type <- annotation_map$broad_cell_type[map_index]

metadata$subtype_label <- annotation_map$subtype_label[map_index]

metadata$annotation_confidence <- annotation_map$annotation_confidence[map_index]

metadata$annotation_review_flag <- annotation_map$review_flag[map_index]

metadata$primary_VDR_localization_eligible <- annotation_map$primary_VDR_localization_eligible[map_index]

metadata$primary_broad_pseudobulk_eligible <- annotation_map$primary_broad_pseudobulk_eligible[map_index]

metadata$subtype_disease_inference_eligible <- annotation_map$subtype_disease_inference_eligible[map_index]

metadata$annotation_version <- annotation_version

if (anyNA(metadata$broad_cell_type) || anyNA(metadata$subtype_label) || nrow(metadata) != 47161L) {
    stop("Annotation application produced missing values or cell loss.")
}

for (column_name in c("primary_cluster", "broad_cell_type", "subtype_label", "annotation_confidence", "annotation_review_flag", 
    "primary_VDR_localization_eligible", "primary_broad_pseudobulk_eligible", "subtype_disease_inference_eligible", 
    "annotation_version")) {
    integrated_object[[column_name]] <- metadata[, column_name, drop = TRUE]
}

integrated_object@misc$GSE183837_broad_annotation <- list(annotation_version = annotation_version, script_version = script_version, 
    annotation_role = "MARKER_SUPPORTED_BROAD_LINEAGE", broad_annotation_locked = TRUE, cluster8_policy = "LOW_CONFIDENCE_EPITHELIAL_LIKE_ATLAS_RETAINED_PRIMARY_MODELS_EXCLUDED", 
    cluster13_policy = "CYCLING_STROMAL_IDENTITY_RETAINED_NO_DISEASE_SUBTYPE_CLAIM", source_marker_version = marker_version, 
    original_checkpoint_md5 = unname(tools::md5sum(input_checkpoint)))

saveRDS(integrated_object, annotated_checkpoint, compress = "gzip")

annotated_checkpoint_md5 <- unname(tools::md5sum(annotated_checkpoint))

annotation_status <- data.frame(annotation_version = annotation_version, status = "COMPLETE", cell_n = ncol(integrated_object), 
    cluster_n = length(unique(integrated_object$primary_cluster)), broad_cell_type_n = length(unique(integrated_object$broad_cell_type)), 
    checkpoint_path = normalize_slashes(annotated_checkpoint), checkpoint_size_bytes = file.info(annotated_checkpoint)$size, 
    checkpoint_md5 = annotated_checkpoint_md5, stringsAsFactors = FALSE)

write_csv_utf8(annotation_status, annotation_status_path)

annotated_metadata <- integrated_object[[]]

broad_cell_counts <- as.data.frame(table(broad_cell_type = annotated_metadata$broad_cell_type), stringsAsFactors = FALSE)

colnames(broad_cell_counts)[2] <- "cell_n"

broad_cell_counts$cell_fraction <- (broad_cell_counts$cell_n/sum(broad_cell_counts$cell_n))

cluster_cell_counts <- as.data.frame(table(cluster = annotated_metadata$primary_cluster, broad_cell_type = annotated_metadata$broad_cell_type, 
    subtype_label = annotated_metadata$subtype_label), stringsAsFactors = FALSE)

cluster_cell_counts <- cluster_cell_counts[cluster_cell_counts$Freq > 0, , drop = FALSE]

colnames(cluster_cell_counts)[4] <- "cell_n"

donor_broad_counts <- as.data.frame(table(sample_accession = annotated_metadata$sample_accession, donor_id = annotated_metadata$donor_id, 
    analysis_group = annotated_metadata$analysis_group, broad_cell_type = annotated_metadata$broad_cell_type), 
    stringsAsFactors = FALSE)

donor_broad_counts <- donor_broad_counts[donor_broad_counts$Freq > 0, , drop = FALSE]

colnames(donor_broad_counts)[5] <- "cell_n"

donor_totals <- aggregate(cell_n ~ sample_accession + donor_id + analysis_group, data = donor_broad_counts, FUN = sum)

colnames(donor_totals)[4] <- "donor_total_cell_n"

donor_broad_counts <- merge(donor_broad_counts, donor_totals, by = c("sample_accession", "donor_id", "analysis_group"), 
    all.x = TRUE, sort = FALSE)

donor_broad_counts$within_donor_fraction <- (donor_broad_counts$cell_n/donor_broad_counts$donor_total_cell_n)

group_broad_summary <- data.table::as.data.table(donor_broad_counts)[, .(donor_n = data.table::uniqueN(donor_id), 
    median_within_donor_fraction = median(within_donor_fraction), minimum_within_donor_fraction = min(within_donor_fraction), 
    maximum_within_donor_fraction = max(within_donor_fraction)), by = .(analysis_group, broad_cell_type)]

group_broad_summary <- as.data.frame(group_broad_summary)

eligibility_counts <- as.data.frame(table(cluster = annotated_metadata$primary_cluster, broad_cell_type = annotated_metadata$broad_cell_type, 
    VDR_localization_eligible = annotated_metadata$primary_VDR_localization_eligible, pseudobulk_eligible = annotated_metadata$primary_broad_pseudobulk_eligible), 
    stringsAsFactors = FALSE)

eligibility_counts <- eligibility_counts[eligibility_counts$Freq > 0, , drop = FALSE]

colnames(eligibility_counts)[5] <- "cell_n"

write_csv_utf8(broad_cell_counts, file.path(table_root, "04_broad_cell_counts.csv"))

write_csv_utf8(cluster_cell_counts, file.path(table_root, "05_cluster_annotation_counts.csv"))

write_csv_utf8(donor_broad_counts, file.path(table_root, "06_donor_broad_composition.csv"))

write_csv_utf8(group_broad_summary, file.path(table_root, "07_group_broad_descriptive.csv"))

write_csv_utf8(eligibility_counts, file.path(table_root, "08_analysis_eligibility.csv"))

visual_marker_list <- c("COL1A1", "DCN", "APOD", "MKI67", "RGS5", "CD3D", "IL7R", "NKG7", "GNLY", "KRT8", "PAX8", 
    "SCGB1D2", "FOXJ1", "C20orf85", "LYZ", "C1QC", "CLDN5", "VWF")

available_features <- rownames(integrated_object)

feature_match <- match(toupper(visual_marker_list), toupper(available_features))

visual_marker_manifest <- data.frame(requested_gene = visual_marker_list, resolved_feature = ifelse(is.na(feature_match), 
    NA_character_, available_features[feature_match]), available = !is.na(feature_match), stringsAsFactors = FALSE)

write_csv_utf8(visual_marker_manifest, file.path(table_root, "09_visual_marker_manifest.csv"))

resolved_visual_features <- unique(visual_marker_manifest$resolved_feature[visual_marker_manifest$available])

expression_matrix <- SeuratObject::LayerData(object = integrated_object, assay = "RNA", layer = "data", features = resolved_visual_features)

marker_average_rows <- list()

for (cluster_id in as.character(0:15)) {
    cell_index <- which(annotated_metadata$primary_cluster == cluster_id)
    average_expression <- Matrix::rowMeans(expression_matrix[, cell_index, drop = FALSE])
    percent_expressing <- Matrix::rowMeans(expression_matrix[, cell_index, drop = FALSE] > 0)
    marker_average_rows[[length(marker_average_rows) + 1L]] <- data.frame(cluster = cluster_id, gene = rownames(expression_matrix), 
        average_expression = as.numeric(average_expression), percent_expressing = as.numeric(percent_expressing), 
        stringsAsFactors = FALSE)
}

marker_average <- do.call(rbind, marker_average_rows)

marker_average$gene_z <- ave(log1p(marker_average$average_expression), marker_average$gene, FUN = safe_zscore)

marker_average$gene_z <- pmax(-3, pmin(3, marker_average$gene_z))

write_csv_utf8(marker_average, file.path(table_root, "10_annotation_marker_average.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

UMAP_embedding <- SeuratObject::Embeddings(integrated_object, reduction = "umap.harmony")

UMAP_metadata <- data.frame(cell_barcode = rownames(UMAP_embedding), UMAPH_1 = UMAP_embedding[, 1], UMAPH_2 = UMAP_embedding[, 
    2], broad_cell_type = integrated_object$broad_cell_type, subtype_label = integrated_object$subtype_label, primary_cluster = integrated_object$primary_cluster, 
    annotation_confidence = integrated_object$annotation_confidence, stringsAsFactors = FALSE)

plot_broad_UMAP <- ggplot2::ggplot(UMAP_metadata, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, group = broad_cell_type)) + 
    ggplot2::geom_point(ggplot2::aes(shape = broad_cell_type), size = 0.17999999999999999, alpha = 0.55000000000000004) + 
    ggplot2::labs(title = "GSE183837 broad cell annotation", subtitle = "Cluster 8 is retained as low-confidence epithelial-like and excluded from primary models", 
        x = "Harmony UMAP 1", y = "Harmony UMAP 2", shape = "Broad cell type") + theme_QC

save_plot_pair(plot_broad_UMAP, "QC1_UMAP_broad_annotation", width = 10, height = 7)

cluster_centroids <- aggregate(cbind(UMAPH_1, UMAPH_2) ~ primary_cluster + subtype_label, data = UMAP_metadata, 
    FUN = median)

plot_cluster_UMAP <- ggplot2::ggplot(UMAP_metadata, ggplot2::aes(x = UMAPH_1, y = UMAPH_2)) + ggplot2::geom_point(size = 0.12, 
    alpha = 0.29999999999999999) + ggplot2::geom_text(data = cluster_centroids, ggplot2::aes(label = primary_cluster), 
    size = 3.5) + ggplot2::labs(title = "Primary clusters underlying the locked broad annotation", subtitle = "Numbers correspond to the reviewed cluster-to-lineage map", 
    x = "Harmony UMAP 1", y = "Harmony UMAP 2") + theme_QC

save_plot_pair(plot_cluster_UMAP, "QC2_UMAP_cluster_labels", width = 9, height = 7)

donor_plot_data <- donor_broad_counts

donor_plot_data$sample_accession <- factor(donor_plot_data$sample_accession, levels = unique(donor_plot_data$sample_accession))

plot_donor_composition <- ggplot2::ggplot(donor_plot_data, ggplot2::aes(x = sample_accession, y = within_donor_fraction, 
    fill = broad_cell_type)) + ggplot2::geom_col() + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Broad-cell composition by donor", subtitle = "Descriptive only; donor is the biological replicate", 
        x = "Donor/sample", y = "Within-donor cell fraction", fill = "Broad cell type") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_donor_composition, "QC3_donor_broad_composition", width = 10, height = 6)

marker_plot_data <- marker_average

marker_plot_data$cluster <- factor(marker_plot_data$cluster, levels = as.character(0:15))

marker_plot_data$gene <- factor(marker_plot_data$gene, levels = rev(resolved_visual_features))

plot_marker_heatmap <- ggplot2::ggplot(marker_plot_data, ggplot2::aes(x = cluster, y = gene, fill = gene_z)) + 
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = "C") + ggplot2::labs(title = "Canonical marker support for the broad annotation", 
    subtitle = "Row-scaled log-normalized average expression", x = "Primary cluster", y = "Marker gene", fill = "Row z-score") + 
    theme_QC

save_plot_pair(plot_marker_heatmap, "QC4_annotation_marker_heatmap", width = 10, height = 8)

confidence_plot_data <- annotation_QC

confidence_plot_data$cluster <- factor(confidence_plot_data$cluster, levels = as.character(0:15))

plot_confidence <- ggplot2::ggplot(confidence_plot_data, ggplot2::aes(x = cluster, y = top100_evidence_n, shape = annotation_confidence, 
    size = cell_n)) + ggplot2::geom_point(alpha = 0.80000000000000004) + ggplot2::labs(title = "Annotation evidence and confidence by cluster", 
    subtitle = "Cluster 8 remains low confidence and sensitivity-only", x = "Primary cluster", y = "Requested evidence markers present in top 100", 
    shape = "Confidence", size = "Cell count") + theme_QC

save_plot_pair(plot_confidence, "QC5_annotation_confidence", width = 9, height = 5.5)

eligibility_plot_data <- cluster_cell_counts

eligibility_plot_data <- merge(eligibility_plot_data, annotation_map[, c("cluster", "primary_broad_pseudobulk_eligible"), 
    drop = FALSE], by = "cluster", all.x = TRUE, sort = FALSE)

eligibility_plot_data$cluster <- factor(eligibility_plot_data$cluster, levels = as.character(0:15))

plot_eligibility <- ggplot2::ggplot(eligibility_plot_data, ggplot2::aes(x = cluster, y = cell_n, shape = primary_broad_pseudobulk_eligible)) + 
    ggplot2::geom_point(size = 3) + ggplot2::scale_y_log10(labels = scales::comma) + ggplot2::labs(title = "Cluster sizes and primary broad-pseudobulk eligibility", 
    subtitle = "Cluster 8 is retained in the atlas but excluded from primary models", x = "Primary cluster", y = "Cell count (log scale)", 
    shape = "Primary eligible") + theme_QC

save_plot_pair(plot_eligibility, "QC6_primary_model_eligibility", width = 8, height = 5)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

cluster8_row <- annotation_map[annotation_map$cluster == "8", , drop = FALSE]

cluster13_row <- annotation_map[annotation_map$cluster == "13", , drop = FALSE]

gate_checks <- data.frame(check = c("upstream_marker_gate_authorizes_annotation_review", "input_checkpoint_has_47161_cells", 
    "annotation_map_has_16_unique_clusters", "all_cells_receive_broad_annotation", "all_cells_receive_subtype_label", 
    "allowed_broad_ontology_only", "all_nonlow_confidence_clusters_have_at_least_2_top100_evidence_markers", "low_confidence_cluster8_has_at_least_1_top100_evidence_marker", 
    "cluster8_is_retained_but_primary_models_excluded", "cluster13_is_cycling_stromal_with_donor_warning", "annotated_checkpoint_has_47161_cells", 
    "annotated_checkpoint_exists", "annotated_checkpoint_md5_present", "donor_broad_composition_has_9_donors", 
    "broad_annotation_can_be_locked", "VDR_localization_can_be_authorized", "disease_testing_remains_prohibited"), 
    pass = c(upstream_gate$gate_status[[1]] == "PASS_DESCRIPTIVE_MARKERS_COMPLETE_CANDIDATE_ANNOTATION_REVIEW_AUTHORIZED", 
        ncol(integrated_object) == 47161L, nrow(annotation_map) == 16L && anyDuplicated(annotation_map$cluster) == 
            0L, !anyNA(integrated_object$broad_cell_type), !anyNA(integrated_object$subtype_label), all(annotation_map$broad_cell_type %in% 
            allowed_broad_types), all(annotation_QC$evidence_pass[annotation_QC$annotation_confidence != "LOW"]), 
        isTRUE(annotation_QC$evidence_pass[annotation_QC$cluster == "8"]), nrow(cluster8_row) == 1L && !cluster8_row$primary_VDR_localization_eligible && 
            !cluster8_row$primary_broad_pseudobulk_eligible, nrow(cluster13_row) == 1L && cluster13_row$broad_cell_type == 
            "Stromal cells" && grepl("DONOR", cluster13_row$review_flag), ncol(integrated_object) == 47161L, file.exists(annotated_checkpoint), 
        nzchar(annotated_checkpoint_md5), length(unique(donor_broad_counts$donor_id)) == 9L, TRUE, TRUE, TRUE), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_BROAD_ANNOTATION_LOCKED_WITH_CLUSTER8_PRIMARY_EXCLUSION"
} else {
    "HOLD_BROAD_ANNOTATION_NOT_LOCKED"
}

annotation_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, marker_version = marker_version, 
    annotation_version = annotation_version, gate_status = gate_status, cell_n = ncol(integrated_object), cluster_n = 16L, 
    broad_cell_type_n = length(unique(integrated_object$broad_cell_type)), broad_cell_types = collapse_unique(sort(unique(integrated_object$broad_cell_type))), 
    low_confidence_clusters = collapse_unique(annotation_map$cluster[annotation_map$annotation_confidence == "LOW"]), 
    primary_excluded_clusters = collapse_unique(annotation_map$cluster[!annotation_map$primary_broad_pseudobulk_eligible]), 
    subtype_disease_restricted_clusters = collapse_unique(annotation_map$cluster[!annotation_map$subtype_disease_inference_eligible]), 
    review_flags = collapse_unique(annotation_map$review_flag), broad_annotation_locked = all_gate_checks_pass, 
    VDR_localization_authorized = all_gate_checks_pass, broad_pseudobulk_preparation_authorized = all_gate_checks_pass, 
    disease_pseudobulk_testing_authorized = FALSE, result_lock_authorized = FALSE, annotated_checkpoint = normalize_slashes(annotated_checkpoint), 
    annotated_checkpoint_md5 = annotated_checkpoint_md5, stringsAsFactors = FALSE)

write_csv_utf8(annotation_gate, file.path(run_qc_root, "G183_ANNOTATION_GATE.csv"))

card_lines <- c("# STEP04E1M GSE183837 Broad Annotation Card", "", paste0("- Gate: ", gate_status), paste0("- Cells retained: ", 
    format(ncol(integrated_object), big.mark = ",")), "- Primary clusters: 16", paste0("- Broad cell types: ", 
    length(unique(integrated_object$broad_cell_type))), paste0("- Broad ontology: ", collapse_unique(sort(unique(integrated_object$broad_cell_type)))), 
    "- Low-confidence cluster: 8", "- Primary model exclusion: cluster 8", "- Donor-dominant subtype restriction: cluster 13", 
    paste0("- Broad annotation locked: ", all_gate_checks_pass), paste0("- VDR localization authorized: ", all_gate_checks_pass), 
    "- Disease pseudobulk testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "G183_ANNOTATION_CARD.md"))

interpretation_lines <- c("# GSE183837 broad annotation interpretation", "", paste0("- Run ID: `", run_id, "`"), 
    paste0("- Gate: **", gate_status, "**"), "", "## Locked broad annotation", "", "- Stromal cells: clusters 0, 1, 2, 3, 4, 5, 10 and 13.", 
    "- NK/T cells: clusters 6, 9 and 14.", "- Unciliated epithelial cells: clusters 7 and 8.", "- Myeloid cells: cluster 11.", 
    "- Endothelial cells: cluster 12.", "- Ciliated epithelial cells: cluster 15.", "", "## Mandatory restrictions", 
    "", "- Cluster 8 is retained as a low-confidence SCGB1D2-positive epithelial-like population but excluded from primary VDR localization and broad pseudobulk models.", 
    "- Cluster 13 is retained as cycling stromal cells; its subtype cannot support an RIF-specific claim because of severe single-donor dominance.", 
    "- Clusters 7 and 15 cannot support subtype-level disease inference at this stage.", "", "## Authorized next operation", 
    "", "- Prepare eligible broad-cell pseudobulks by donor.", "- Localize VDR expression and the four locked VDR programs.", 
    "- Run donor-level RIF comparisons only after pseudobulk QC and with explicit low-power wording.", "", "## Not authorized", 
    "", "- No cell-level Control-versus-RIF inference.", "- No cluster-8 primary model.", "- No disease mechanism claim.", 
    "- No final result lock.")

write_text(interpretation_lines, file.path(run_qc_root, "G183_ANNOTATION_INTERPRETATION.md"))

parameter_manifest <- list(script_name = script_name, script_version = script_version, dataset_id = dataset_id, 
    integration_version = integration_version, marker_version = marker_version, annotation_version = annotation_version, 
    upstream_marker_gate = upstream_gate$gate_status[[1]], broad_ontology = allowed_broad_types, cluster8_policy = "ATLAS_RETAINED_LOW_CONFIDENCE_PRIMARY_VDR_AND_PSEUDOBULK_EXCLUDED", 
    cluster13_policy = "CYCLING_STROMAL_RETAINED_NO_DISEASE_SUBTYPE_CLAIM", broad_annotation_locked = all_gate_checks_pass, 
    VDR_localization_authorized = all_gate_checks_pass, disease_pseudobulk_testing_authorized = FALSE, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, dataset_id = dataset_id, 
    integration_version = integration_version, marker_version = marker_version, annotation_version = annotation_version, 
    started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, 
        run_start, units = "mins")), 2), gate_status = gate_status, cell_n = ncol(integrated_object), cluster_n = 16L, 
    broad_cell_type_n = length(unique(integrated_object$broad_cell_type)), low_confidence_clusters = "8", primary_excluded_clusters = "8", 
    broad_annotation_locked = all_gate_checks_pass, VDR_localization_authorized = all_gate_checks_pass, disease_pseudobulk_testing_authorized = FALSE, 
    output_root = normalize_slashes(output_root), qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

persistent_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)

persistent_files <- persistent_files[file.exists(persistent_files) & !dir.exists(persistent_files)]

persistent_manifest <- data.frame(path = normalize_slashes(persistent_files), size_bytes = file.info(persistent_files)$size, 
    md5 = unname(tools::md5sum(persistent_files)), stringsAsFactors = FALSE)

write_csv_utf8(persistent_manifest, file.path(manifest_root, "persistent_manifest.csv"))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(figure_output_root, 
    recursive = TRUE, full.names = TRUE), list.files(manifest_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), annotation_status_path, log_path, run_manifest_path, parameter_path, 
    session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2]

artifact_manifest <- data.frame(path = normalize_slashes(small_files), size_bytes = file.info(small_files)$size, 
    md5 = unname(tools::md5sum(small_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "cells=", ncol(integrated_object), "broad_types=", 
    length(unique(integrated_object$broad_cell_type)))

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

cat("STEP04E1M GSE183837 broad annotation QC completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Cells retained: ", format(ncol(integrated_object), big.mark = ","), "\n", sep = "")

cat("Primary clusters annotated: 16\n")

cat("Broad cell types: ", length(unique(integrated_object$broad_cell_type)), "\n", sep = "")

cat("Low-confidence cluster: 8\n")

cat("Primary model exclusion: 8\n")

cat("Broad annotation locked: ", all_gate_checks_pass, "\n", sep = "")

cat("VDR localization authorized: ", all_gate_checks_pass, "\n", sep = "")

cat("Annotated checkpoint: ", annotated_checkpoint, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nAnnotation outputs are complete in the project folders.\n", sep = "")
}

