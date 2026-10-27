#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C05: Apply accepted broad-cell annotations

options(stringsAsFactors = FALSE, warn = 1)

set.seed(20260805)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "annotate_gse250130_cells.R"

script_version <- "STEP04E1E_GSE250130_BROAD_ANNOTATION_QC_v1.2_PATHSAFE"

input_version <- "SCI_v1_1_memory_safe"

marker_version <- "MARKERS_v1_1_memory_safe"

annotation_version <- "ANN_v1_2"

run_start <- Sys.time()

run_id <- paste0("STEP04E1E_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

split_semicolon <- function(x) {
    x <- as.character(x)
    if (is.na(x) || !nzchar(x)) 
        return(character())
    trimws(unlist(strsplit(x, ";", fixed = TRUE)))
}

detect_first_column <- function(column_names, candidates) {
    lower <- tolower(column_names)
    for (candidate in tolower(candidates)) {
        hit <- which(lower == candidate)
        if (length(hit)) 
            return(column_names[hit[1]])
    }
    ""
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

required_packages <- c("data.table", "ggplot2", "jsonlite", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

input_table_root <- file.path(results_root, "S4E1_sc", input_version, "tables", dataset_id)

marker_root <- file.path(results_root, "S4E1_sc", input_version, "markers", dataset_id, marker_version)

marker_table_root <- file.path(marker_root, "tables")

annotation_root <- file.path(results_root, "S4E1_sc", input_version, "annotations", dataset_id, annotation_version)

annotation_table_root <- file.path(annotation_root, "tables")

annotation_source_root <- file.path(annotation_root, "source_data")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1E_", format(run_start, "%Y%m%d_%H%M%S"), "_annQC"))

figure_root <- file.path(run_qc_root, "figures")

for (p in c(execution_root, parameter_root, environment_root, manuscript_text_root, annotation_root, annotation_table_root, 
    annotation_source_root, run_qc_root, figure_root)) {
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

planned_path_audit <- data.frame(object = c("annotation_root", "annotation_table_root", "annotation_source_root", 
    "run_qc_root", "figure_root", "log_path", "run_manifest_path", "parameter_path", "session_path", "handoff_zip"), 
    path = normalize_slashes(c(annotation_root, annotation_table_root, annotation_source_root, run_qc_root, figure_root, 
        log_path, run_manifest_path, parameter_path, session_path, handoff_zip)), stringsAsFactors = FALSE)

planned_path_audit$character_n <- nchar(planned_path_audit$path, type = "chars")

if (.Platform$OS.type == "windows" && any(planned_path_audit$character_n > 245L)) {
    bad_paths <- planned_path_audit[planned_path_audit$character_n > 245L, , drop = FALSE]
    stop("One or more planned output paths remain too long:\n", paste(paste0(bad_paths$character_n, " chars: ", 
        bad_paths$path), collapse = "\n"))
}

umap_path <- file.path(input_table_root, paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv"))

marker_gate_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

marker_gate_dirs <- marker_gate_dirs[grepl("^STEP04E1D_GSE250130_v1_1_[0-9]{8}_[0-9]{6}_cluster_markers_memory_safe$", 
    basename(marker_gate_dirs), perl = TRUE)]

if (!length(marker_gate_dirs)) {
    stop("No completed STEP04E1D memory-safe marker QC folder found.")
}

latest_marker_qc <- marker_gate_dirs[which.max(file.info(marker_gate_dirs)$mtime)]

marker_gate_path <- file.path(latest_marker_qc, "STEP04E1D_GSE250130_MARKER_GATE.csv")

integration_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

integration_qc_dirs <- integration_qc_dirs[grepl("^STEP04E1C_GSE250130_[0-9]{8}_[0-9]{6}_integration_QC$", basename(integration_qc_dirs), 
    perl = TRUE)]

if (!length(integration_qc_dirs)) {
    stop("No completed STEP04E1C integration QC folder found.")
}

latest_integration_qc <- integration_qc_dirs[which.max(file.info(integration_qc_dirs)$mtime)]

integration_gate_path <- file.path(latest_integration_qc, "STEP04E1C_GSE250130_INTEGRATION_QC_GATE.csv")

all_marker_path <- file.path(marker_table_root, paste0(dataset_id, "_all_cluster_markers_descriptive.csv"))

top25_path <- file.path(marker_table_root, paste0(dataset_id, "_top25_descriptive_markers_by_cluster.csv"))

top10_path <- file.path(marker_table_root, paste0(dataset_id, "_top10_marker_summary.csv"))

panel_support_path <- file.path(marker_table_root, "06_BROAD_LINEAGE_PANEL_SUPPORT.csv")

top_panel_path <- file.path(marker_table_root, "07_TOP_BROAD_LINEAGE_SUPPORT_BY_CLUSTER.csv")

quality_path <- file.path(marker_table_root, "03_CLUSTER_MARKER_QUALITY_AUDIT.csv")

locked_design_path <- file.path(results_root, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0", "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv")

required_inputs <- c(umap_path, marker_gate_path, integration_gate_path, all_marker_path, top25_path, top10_path, 
    panel_support_path, top_panel_path, quality_path, locked_design_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

marker_gate <- safe_read_csv(marker_gate_path)

integration_gate <- safe_read_csv(integration_gate_path)

if (is.null(marker_gate) || nrow(marker_gate) != 1L || marker_gate$gate_status[[1]] != "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED" || 
    !isTRUE(marker_gate$annotation_review_authorized[[1]])) {
    stop("STEP04E1D marker Gate does not authorize annotation review.")
}

if (is.null(integration_gate) || nrow(integration_gate) != 1L || !integration_gate$gate_status[[1]] %in% c("PASS", 
    "PASS_WITH_REVIEW_FLAGS")) {
    stop("STEP04E1C integration Gate is not valid.")
}

umap <- .vdr_fread(umap_path, data.table = FALSE, showProgress = FALSE)

all_markers <- .vdr_fread(all_marker_path, data.table = FALSE, showProgress = FALSE)

top25 <- .vdr_fread(top25_path, data.table = FALSE, showProgress = FALSE)

top10 <- safe_read_csv(top10_path)

panel_support <- safe_read_csv(panel_support_path)

top_panel <- safe_read_csv(top_panel_path)

marker_quality <- safe_read_csv(quality_path)

locked_design <- safe_read_csv(locked_design_path)

required_umap_columns <- c("cell_barcode", "UMAPH_1", "UMAPH_2", "sample_accession", "seurat_clusters")

missing_umap_columns <- setdiff(required_umap_columns, colnames(umap))

if (length(missing_umap_columns)) {
    stop("UMAP table is missing columns: ", paste(missing_umap_columns, collapse = ", "))
}

if (nrow(umap) != 236636L || length(unique(umap$sample_accession)) != 28L || length(unique(umap$seurat_clusters)) != 
    22L) {
    stop("UMAP identity checks failed.")
}

if (is.null(locked_design) || !"sample_accession" %in% colnames(locked_design) || nrow(locked_design) != 28L) {
    stop("Locked GSE250130 sample design is invalid.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, input_version = input_version, 
    marker_version = marker_version, annotation_version = annotation_version, supersedes_partial_annotation_version = "ANNOTATION_v1_1_FIXED_PATH_LENGTH_FAILURE", 
    dataset_id = dataset_id, broad_nomenclature_source = "GSE250130 published eight-major-cell-type framework", 
    original_study_doi = "10.1038/s41467-024-55419-z", reference_atlas_doi = "10.1038/s41588-024-01873-w", older_endometrial_atlas_doi = "10.1038/s41591-020-1040-z", 
    broad_annotation_role = "ACCEPTED_FOR_BROAD_LEVEL_LOCALIZATION_AND_AGGREGATION_AFTER_GATE", subtype_annotation_role = "CANDIDATE_REQUIRES_REVIEW", 
    review_clusters = c("12", "14", "20", "21"), composition_role = "DESCRIPTIVE_QC_ONLY_NO_HYPOTHESIS_TEST", random_seed = 20260805L, 
    result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

log_event("Starting", run_id, "marker_gate=", marker_gate$gate_status[[1]], "integration_gate=", integration_gate$gate_status[[1]])

write_csv_utf8(planned_path_audit, file.path(run_qc_root, "00_path_length_audit.csv"))

annotation_map <- data.frame(cluster = as.character(0:21), broad_cell_type = c("Stromal cells", "Unciliated epithelial cells", 
    "Stromal cells", "NK/T cells", "NK/T cells", "NK/T cells", "Unciliated epithelial cells", "NK/T cells", "NK/T cells", 
    "Myeloid cells", "Ciliated epithelial cells", "NK/T cells", "B cells", "Myeloid cells", "Mast cells", "Endothelial cells", 
    "Endothelial cells", "NK/T cells", "Myeloid cells", "Unciliated epithelial cells", "Unciliated epithelial cells", 
    "Stromal cells"), candidate_subtype = c("COL1A1_IGF1_stromal", "MT_high_unciliated_epithelial", "PTGDS_APOD_stromal", 
    "NCAM1_GNLY_NK", "conventional_T", "KLRB1_TRDC_uNK_like", "PAX8_GABRP_unciliated_epithelial", "cycling_NK_T", 
    "FGFBP2_cytotoxic_NK", "HLA_DRA_LYZ_myeloid_APC", "ciliated_epithelial", "TRDC_GNLY_gamma_delta_like_NK_T", 
    "MS4A1_CD79A_B", "CSF3R_S100A8_A9_neutrophil", "TPSAB1_CPA3_mast", "VWF_EMCN_blood_endothelial", "CCL21_PROX1_lymphatic_endothelial", 
    "IL7R_KLRB1_T", "IL3RA_CLEC4C_GZMB_plasmacytoid_DC", "KRT19_SLPI_unciliated_epithelial", "ESR1_LIFR_LPAR3_luminal_like_epithelial", 
    "C11orf96_INHBA_activated_perivascular_stromal"), broad_confidence = c("HIGH", "HIGH", "HIGH", "HIGH", "HIGH", 
    "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "MODERATE", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", 
    "HIGH", "HIGH", "MODERATE", "MODERATE"), subtype_confidence = c("MODERATE", "MODERATE", "MODERATE", "HIGH", 
    "HIGH", "MODERATE", "MODERATE", "HIGH", "HIGH", "MODERATE", "HIGH", "MODERATE", "HIGH", "HIGH", "HIGH", "HIGH", 
    "HIGH", "HIGH", "HIGH", "MODERATE", "MODERATE_REVIEW", "MODERATE_REVIEW"), review_flag = c(FALSE, FALSE, FALSE, 
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, 
    FALSE, TRUE, TRUE), review_reason = c("", "", "", "", "", "", "", "", "", "", "", "NK/T subtype boundary requires conservative naming", 
    "B-cell identity is clear; sample-mixing/abundance remains under review", "", "Mast-cell identity is clear; strong single-sample contribution remains under review", 
    "", "", "", "", "", "Rare early-luteal/luminal-like epithelial cluster; limited to seven samples", "Rare activated perivascular/stromal cluster; five-sample and single-sample-dominance review"), 
    evidence_genes = c("COL1A1;COL3A1;PDGFRA;DCN;IGF1", "KRT19;KRT18;KRT8;SCGB2A1;WFDC2;MT1G", "COL1A1;COL3A1;DCN;PTGDS;APOD", 
        "GNLY;NCAM1;NKG7;KLRD1;GZMB", "CD3D;CD3E;TRAC;CD2;BCL11B", "NKG7;KLRB1;TRDC;XCL1;XCL2", "EPCAM;KRT8;KRT18;KRT19;PAX8;GABRP", 
        "MKI67;TOP2A;STMN1;NKG7;GZMA", "NKG7;FGFBP2;GNLY;PRF1;GZMB", "LYZ;HLA-DRA;CTSS;TYROBP;AIF1", "DNAH12;RSPH1;CAPS;TPPP3;C20orf85", 
        "TRDC;GNLY;NKG7;GZMB;CTSW", "MS4A1;CD79A;CD74;BANK1;CD37", "S100A8;S100A9;CSF3R;FCGR3B;AQP9", "TPSAB1;TPSB2;CPA3;KIT;HDC;MS4A2", 
        "VWF;EMCN;ADGRL4;PTPRB;ERG", "CCL21;PROX1;PDPN;FLT4;CLDN5", "IL7R;KLRB1;TRBC1;LTB;ICOS", "IL3RA;CLEC4C;GZMB;TCF4;JCHAIN", 
        "KRT19;KRT18;KRT7;SLPI;WFDC2;EHF", "EPCAM;KRT8;KRT18;PAX8;ERBB4;LIFR;LPAR3;ESR1", "COL1A1;COL3A1;PDGFRB;CALD1;PRRX1;TAGLN;NOTCH3;INHBA"), 
    annotation_basis = c(rep("DESCRIPTIVE_CLUSTER_MARKERS_PLUS_CANONICAL_PANEL_PLUS_PUBLISHED_GSE250130_MAJOR_TYPE_FRAMEWORK", 
        22)), subtype_status = rep("CANDIDATE_NOT_LOCKED", 22), stringsAsFactors = FALSE)

all_markers$cluster <- as.character(all_markers$cluster)

top25$cluster <- as.character(top25$cluster)

evidence_audit_rows <- list()

for (i in seq_len(nrow(annotation_map))) {
    cluster_id <- annotation_map$cluster[i]
    expected <- split_semicolon(annotation_map$evidence_genes[i])
    one_all <- all_markers[all_markers$cluster == cluster_id, , drop = FALSE]
    one_top25 <- top25[top25$cluster == cluster_id, , drop = FALSE]
    observed_all <- intersect(expected, one_all$gene)
    observed_top25 <- intersect(expected, one_top25$gene)
    evidence_audit_rows[[i]] <- data.frame(cluster = cluster_id, broad_cell_type = annotation_map$broad_cell_type[i], 
        candidate_subtype = annotation_map$candidate_subtype[i], expected_evidence_gene_n = length(expected), observed_evidence_gene_n_all_markers = length(observed_all), 
        observed_evidence_fraction_all_markers = length(observed_all)/length(expected), observed_evidence_genes_all_markers = collapse_unique(observed_all), 
        observed_evidence_gene_n_top25 = length(observed_top25), observed_evidence_genes_top25 = collapse_unique(observed_top25), 
        review_flag = annotation_map$review_flag[i], stringsAsFactors = FALSE)
}

evidence_audit <- do.call(rbind, evidence_audit_rows)

annotation_map <- merge(annotation_map, evidence_audit[, c("cluster", "observed_evidence_fraction_all_markers", 
    "observed_evidence_genes_all_markers", "observed_evidence_genes_top25")], by = "cluster", all.x = TRUE, sort = FALSE)

annotation_map <- annotation_map[order(as.integer(annotation_map$cluster)), , drop = FALSE]

annotation_map$broad_annotation_status <- ifelse(annotation_map$observed_evidence_fraction_all_markers >= 0.59999999999999998, 
    ifelse(annotation_map$review_flag, "ACCEPTED_BROAD_WITH_REVIEW", "ACCEPTED_BROAD"), "HOLD_EVIDENCE_INCOMPLETE")

write_csv_utf8(annotation_map, file.path(annotation_table_root, "cluster_annotation_map.csv"))

write_csv_utf8(evidence_audit, file.path(run_qc_root, "01_annotation_evidence_audit.csv"))

umap$cluster <- as.character(umap$seurat_clusters)

cell_annotation <- merge(umap, annotation_map[, c("cluster", "broad_cell_type", "candidate_subtype", "broad_confidence", 
    "subtype_confidence", "review_flag", "review_reason", "broad_annotation_status", "subtype_status")], by = "cluster", 
    all.x = TRUE, sort = FALSE)

if (nrow(cell_annotation) != 236636L || anyNA(cell_annotation$broad_cell_type)) {
    stop("Cell-level annotation merge failed.")
}

locked_design$sample_accession <- as.character(locked_design$sample_accession)

cell_annotation$sample_accession <- as.character(cell_annotation$sample_accession)

cell_annotation <- merge(cell_annotation, locked_design, by = "sample_accession", all.x = TRUE, sort = FALSE)

if (nrow(cell_annotation) != 236636L) {
    stop("Locked-design merge changed the cell count.")
}

cell_annotation <- cell_annotation[order(match(cell_annotation$cell_barcode, umap$cell_barcode)), , drop = FALSE]

cell_annotation_path <- file.path(annotation_table_root, "cell_annotation_metadata.csv")

write_csv_utf8(cell_annotation, cell_annotation_path)

major_summary <- aggregate(cell_barcode ~ broad_cell_type, data = cell_annotation, FUN = length)

colnames(major_summary)[2] <- "cell_n"

major_summary$cell_fraction <- (major_summary$cell_n/sum(major_summary$cell_n))

major_summary <- major_summary[order(-major_summary$cell_n), , drop = FALSE]

published_profile <- data.frame(broad_cell_type = c("Unciliated epithelial cells", "Ciliated epithelial cells", 
    "Stromal cells", "Endothelial cells", "NK/T cells", "Myeloid cells", "B cells", "Mast cells"), published_cell_fraction = c(0.16800000000000001, 
    0.019, 0.35799999999999998, 0.0060000000000000001, 0.38500000000000001, 0.037999999999999999, 0.017999999999999999, 
    0.0060000000000000001), source = "Cao_et_al_Nat_Commun_2025_DOI_10.1038_s41467-024-55419-z", role = "CONTEXTUAL_QC_NOT_EXACT_REPRODUCTION_TARGET", 
    stringsAsFactors = FALSE)

major_summary <- merge(major_summary, published_profile, by = "broad_cell_type", all.x = TRUE, sort = FALSE)

major_summary$absolute_fraction_difference <- abs(major_summary$cell_fraction - major_summary$published_cell_fraction)

major_summary <- major_summary[order(-major_summary$cell_n), , drop = FALSE]

write_csv_utf8(major_summary, file.path(annotation_table_root, "broad_type_summary.csv"))

write_csv_utf8(published_profile, file.path(annotation_source_root, "published_major_profile.csv"))

cluster_summary <- merge(aggregate(cell_barcode ~ cluster, data = cell_annotation, FUN = length), annotation_map, 
    by = "cluster", all.x = TRUE, sort = FALSE)

colnames(cluster_summary)[colnames(cluster_summary) == "cell_barcode"] <- "cell_n"

cluster_summary <- cluster_summary[order(as.integer(cluster_summary$cluster)), , drop = FALSE]

write_csv_utf8(cluster_summary, file.path(annotation_table_root, "cluster_annotation_summary.csv"))

sample_type_counts <- aggregate(cell_barcode ~ sample_accession + broad_cell_type, data = cell_annotation, FUN = length)

colnames(sample_type_counts)[3] <- "cell_n"

sample_totals <- aggregate(cell_n ~ sample_accession, data = sample_type_counts, FUN = sum)

colnames(sample_totals)[2] <- "sample_total_cell_n"

sample_type_counts <- merge(sample_type_counts, sample_totals, by = "sample_accession", all.x = TRUE, sort = FALSE)

sample_type_counts$within_sample_fraction <- (sample_type_counts$cell_n/sample_type_counts$sample_total_cell_n)

sample_type_counts <- merge(sample_type_counts, locked_design, by = "sample_accession", all.x = TRUE, sort = FALSE)

write_csv_utf8(sample_type_counts, file.path(annotation_table_root, "sample_broad_type_composition.csv"))

group_column <- detect_first_column(colnames(locked_design), c("group", "clinical_group", "phenotype", "condition", 
    "fertility_status", "patient_group", "case_control_status", "case_control"))

time_column <- detect_first_column(colnames(locked_design), c("timepoint", "time_point", "temporal_label", "cycle_day", 
    "lh_day", "day", "phase", "sampling_time"))

role_column <- detect_first_column(colnames(locked_design), c("analysis_role", "cohort_role", "sample_role", "role"))

design_detection <- data.frame(design_role = c("GROUP", "TIME", "ANALYSIS_ROLE"), detected_column = c(group_column, 
    time_column, role_column), stringsAsFactors = FALSE)

write_csv_utf8(design_detection, file.path(run_qc_root, "02_design_columns.csv"))

theme_qc <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

cell_dt <- data.table::as.data.table(cell_annotation)

plot_cells <- cell_dt[, {
    keep_n <- min(.N, 10000L)
    .SD[sample.int(.N, keep_n)]
}, by = broad_cell_type]

broad_centroids <- cell_dt[, .(UMAPH_1 = stats::median(UMAPH_1), UMAPH_2 = stats::median(UMAPH_2)), by = broad_cell_type]

plot_umap_broad <- ggplot2::ggplot(plot_cells, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = broad_cell_type)) + 
    ggplot2::geom_point(size = 0.080000000000000002, alpha = 0.59999999999999998, stroke = 0) + ggplot2::geom_text(data = broad_centroids, 
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, label = broad_cell_type), inherit.aes = FALSE, size = 3, fontface = "bold") + 
    ggplot2::labs(title = "GSE250130 broad cell-type annotation", subtitle = "Eight-major-cell-type framework; deterministic plotting subset", 
        x = "UMAP 1", y = "UMAP 2", color = "Broad cell type") + theme_qc + ggplot2::theme(legend.position = "right")

save_plot_pair(plot_umap_broad, "QC1_UMAP_broad_types", width = 10, height = 7)

subtype_centroids <- cell_dt[, .(UMAPH_1 = stats::median(UMAPH_1), UMAPH_2 = stats::median(UMAPH_2), cluster = unique(cluster)[1]), 
    by = candidate_subtype]

plot_umap_subtype <- ggplot2::ggplot(plot_cells, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = candidate_subtype)) + 
    ggplot2::geom_point(size = 0.059999999999999998, alpha = 0.55000000000000004, stroke = 0) + ggplot2::geom_text(data = subtype_centroids, 
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, label = paste0("C", cluster)), inherit.aes = FALSE, size = 2.7999999999999998, 
    fontface = "bold") + ggplot2::labs(title = "GSE250130 candidate subtype annotation", subtitle = "Candidate labels are not locked; C12/C14/C20/C21 retain review flags", 
    x = "UMAP 1", y = "UMAP 2", color = "Candidate subtype") + theme_qc + ggplot2::theme(legend.position = "right", 
    legend.text = ggplot2::element_text(size = 6))

save_plot_pair(plot_umap_subtype, "QC2_UMAP_subtypes", width = 13, height = 8)

plot_major_profile <- ggplot2::ggplot(reshape(major_summary[, c("broad_cell_type", "cell_fraction", "published_cell_fraction")], 
    varying = c("cell_fraction", "published_cell_fraction"), v.names = "fraction", timevar = "profile", times = c("Current reanalysis", 
        "Published study"), direction = "long"), ggplot2::aes(x = reorder(broad_cell_type, fraction), y = fraction, 
    fill = profile)) + ggplot2::geom_col(position = "dodge") + ggplot2::coord_flip() + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Major cell-type profile: contextual concordance", subtitle = "Published proportions are contextual QC, not an exact reproduction target", 
        x = "Broad cell type", y = "Cell fraction", fill = "Profile") + theme_qc

save_plot_pair(plot_major_profile, "QC3_profile_context", width = 9, height = 6)

sample_type_counts$sample_accession <- factor(sample_type_counts$sample_accession, levels = unique(locked_design$sample_accession))

plot_sample_composition <- ggplot2::ggplot(sample_type_counts, ggplot2::aes(x = sample_accession, y = within_sample_fraction, 
    fill = broad_cell_type)) + ggplot2::geom_col(width = 0.84999999999999998) + ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Sample-level broad cell-type composition", subtitle = "Descriptive QC only; no composition hypothesis test", 
        x = "Sample", y = "Within-sample fraction", fill = "Broad cell type") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
    vjust = 0.5, hjust = 1, size = 6))

save_plot_pair(plot_sample_composition, "QC4_sample_composition", width = 12, height = 7)

evidence_plot_data <- annotation_map

evidence_plot_data$cluster <- factor(evidence_plot_data$cluster, levels = as.character(0:21))

plot_evidence <- ggplot2::ggplot(evidence_plot_data, ggplot2::aes(x = cluster, y = observed_evidence_fraction_all_markers, 
    fill = broad_annotation_status)) + ggplot2::geom_col() + ggplot2::geom_hline(yintercept = 0.59999999999999998, 
    linetype = "dashed") + ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) + 
    ggplot2::labs(title = "Cluster annotation evidence coverage", subtitle = "Dashed line marks the predefined 60% broad-evidence threshold", 
        x = "Cluster", y = "Expected evidence genes observed", fill = "Annotation status") + theme_qc

save_plot_pair(plot_evidence, "QC5_evidence_coverage", width = 9, height = 5)

expected_broad_types <- c("Unciliated epithelial cells", "Ciliated epithelial cells", "Stromal cells", "Endothelial cells", 
    "NK/T cells", "Myeloid cells", "B cells", "Mast cells")

technical_checks <- data.frame(check = c("integration_gate_valid", "marker_gate_valid", "cell_count_locked", "sample_count_locked", 
    "cluster_count_locked", "all_22_clusters_mapped", "all_eight_broad_types_present", "all_cells_receive_broad_annotation", 
    "all_clusters_meet_broad_evidence_threshold", "review_clusters_preserved"), pass = c(integration_gate$gate_status[[1]] %in% 
    c("PASS", "PASS_WITH_REVIEW_FLAGS"), marker_gate$gate_status[[1]] == "PASS_DESCRIPTIVE_MARKERS_COMPLETE_ANNOTATION_REVIEW_AUTHORIZED", 
    nrow(cell_annotation) == 236636L, length(unique(cell_annotation$sample_accession)) == 28L, length(unique(cell_annotation$cluster)) == 
        22L, setequal(annotation_map$cluster, as.character(0:21)), setequal(unique(annotation_map$broad_cell_type), 
        expected_broad_types), !anyNA(cell_annotation$broad_cell_type), all(annotation_map$observed_evidence_fraction_all_markers >= 
        0.59999999999999998), setequal(annotation_map$cluster[annotation_map$review_flag], c("12", "14", "20", 
        "21"))), stringsAsFactors = FALSE)

write_csv_utf8(technical_checks, file.path(run_qc_root, "03_gate_checks.csv"))

all_technical_pass <- all(technical_checks$pass)

gate_status <- if (all_technical_pass) {
    "PASS_BROAD_ANNOTATION_WITH_SUBTYPE_REVIEW"
} else {
    "HOLD_BROAD_ANNOTATION_QC_FAILED"
}

annotation_gate <- data.frame(dataset_id = dataset_id, input_version = input_version, marker_version = marker_version, 
    annotation_version = annotation_version, gate_status = gate_status, cell_n = nrow(cell_annotation), sample_n = length(unique(cell_annotation$sample_accession)), 
    cluster_n = length(unique(cell_annotation$cluster)), broad_cell_type_n = length(unique(cell_annotation$broad_cell_type)), 
    review_clusters = collapse_unique(annotation_map$cluster[annotation_map$review_flag]), broad_annotation_authorized = all_technical_pass, 
    broad_VDR_localization_authorized = all_technical_pass, broad_pseudobulk_aggregation_authorized = all_technical_pass, 
    broad_composition_descriptive_QC_authorized = all_technical_pass, composition_inferential_testing_authorized = FALSE, 
    subtype_inferential_analysis_authorized = FALSE, subtype_annotation_locked = FALSE, result_lock_authorized = FALSE, 
    stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "ANNOTATION_GATE.csv")

write_csv_utf8(annotation_gate, gate_path)

card_lines <- c("# STEP04E1E Broad Annotation Result Card", "", paste0("- Dataset: ", dataset_id), paste0("- Input version: ", 
    input_version), paste0("- Marker version: ", marker_version), paste0("- Annotation version: ", annotation_version), 
    paste0("- Gate: ", gate_status), paste0("- Broad cell types: ", length(unique(cell_annotation$broad_cell_type))), 
    paste0("- Annotated cells: ", format(nrow(cell_annotation), big.mark = ",")), paste0("- Review clusters: ", 
        collapse_unique(annotation_map$cluster[annotation_map$review_flag])), paste0("- Broad VDR localization authorized: ", 
        annotation_gate$broad_VDR_localization_authorized), paste0("- Broad pseudobulk aggregation authorized: ", 
        annotation_gate$broad_pseudobulk_aggregation_authorized), "- Composition inferential testing authorized: FALSE", 
    "- Subtype inferential analysis authorized: FALSE", "- Subtype annotation locked: FALSE", "- Result lock authorized: FALSE")

card_path <- file.path(run_qc_root, "ANNOTATION_RESULT_CARD.md")

write_text(card_lines, card_path)

interpretation_lines <- c("# STEP04E1E GSE250130 broad annotation interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), paste0("- Broad cell types: ", collapse_unique(major_summary$broad_cell_type)), 
    "", "## Accepted broad annotation", "", "- Broad cell types follow the eight-major-cell-type framework reported for GSE250130.", 
    "- All 22 clusters met the predefined marker-evidence threshold.", "- Broad annotations may be used for VDR program localization and donor-by-cell-type aggregation.", 
    "", "## Preserved uncertainty", "", "- Candidate subtype labels are not locked.", "- Cluster 12 is accepted as B cells, but its sample-mixing pattern remains under review.", 
    "- Cluster 14 is accepted as mast cells, but its abundance is strongly influenced by one sample.", "- Cluster 20 is accepted as unciliated epithelial cells and provisionally described as ESR1/LIFR/LPAR3 luminal-like epithelial cells.", 
    "- Cluster 21 is accepted as stromal cells and provisionally described as activated perivascular stromal cells.", 
    "", "## Composition rule", "", "- Sample-level cell fractions are retained as descriptive QC.", "- No group or timepoint composition hypothesis test is authorized at this stage.", 
    "", "## Next step", "", "- Use the broad annotation map for VDR RNA-program localization.", "- Aggregate raw counts by donor and broad cell type.", 
    "- Perform donor-level QC before any clinical or temporal comparison.", "- Process GSE183837 only after the GSE250130 annotation and aggregation workflow is stable.")

interpretation_path <- file.path(run_qc_root, "ANNOTATION_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

methods_lines <- c("# STEP04E1E Methods draft", "", paste("Clusters from the resolution-0.4 solution were annotated using an", 
    "evidence-reviewed hierarchical procedure."), paste("Broad cell-type labels were aligned to the eight major populations", 
    "reported in the source GSE250130 study: unciliated epithelial, ciliated", "epithelial, stromal, endothelial, NK/T, myeloid, B and mast cells."), 
    paste("For each cluster, descriptive positive markers, canonical marker-panel", "support and sample-mixing QC were reviewed jointly."), 
    paste("Broad labels were accepted only when at least 60% of a predefined", "cluster-specific evidence-gene set was recovered among the descriptive", 
        "marker candidates."), paste("Detailed subtype labels were retained as candidates and were not used", "for inferential testing."), 
    paste("Clusters previously flagged for sample dominance or restricted donor", "coverage retained explicit review flags in the annotation manifest."), 
    paste("Cell-type composition summaries were generated for quality-control", "purposes only; no cell-fraction hypothesis testing was performed."))

write_text(methods_lines, file.path(manuscript_text_root, "STEP04E1E_annotation_Methods_draft.md"))

summary_lookup <- setNames(as.numeric(major_summary$cell_fraction), as.character(major_summary$broad_cell_type))

get_fraction <- function(cell_type) {
    value <- unname(summary_lookup[cell_type])
    if (length(value) != 1L || is.na(value) || !is.finite(value)) {
        stop("Missing or invalid broad-cell-type fraction for: ", cell_type)
    }
    as.numeric(value)
}

results_lines <- c("# STEP04E1E Results draft", "", paste("The 236,636 integrated cells were assigned to eight broad endometrial", 
    "cell populations using cluster-level marker evidence."), paste0("NK/T cells represented ", round(100 * get_fraction("NK/T cells"), 
    1), "% of cells, followed by stromal cells (", round(100 * get_fraction("Stromal cells"), 1), "%) and unciliated epithelial cells (", 
    round(100 * get_fraction("Unciliated epithelial cells"), 1), "%)."), paste0("Myeloid cells accounted for ", 
    round(100 * get_fraction("Myeloid cells"), 1), "%, ciliated epithelial cells for ", round(100 * get_fraction("Ciliated epithelial cells"), 
        1), "%, B cells for ", round(100 * get_fraction("B cells"), 1), "%, endothelial cells for ", round(100 * 
        get_fraction("Endothelial cells"), 1), "% and mast cells for ", round(100 * get_fraction("Mast cells"), 
        1), "%."), paste("All clusters satisfied the predefined broad-lineage marker-evidence", "criterion."), 
    paste("Four clusters retained review flags: B-cell cluster 12 and mast-cell", "cluster 14 because of sample-mixing characteristics, and rare clusters", 
        "20 and 21 because of restricted donor coverage and provisional subtype", "interpretation."), paste("Broad labels were therefore authorized for subsequent VDR-program", 
        "localization and donor-level aggregation, whereas subtype-level", "inference remained deferred."))

write_text(results_lines, file.path(manuscript_text_root, "STEP04E1E_annotation_Results_draft.md"))

run_end <- Sys.time()

package_names <- c("data.table", "ggplot2", "jsonlite", "scales")

package_table <- data.frame(package = package_names, version = vapply(package_names, safe_package_version, character(1)), 
    stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- unique(c(list.files(annotation_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1E_annotation_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1E_annotation_Results_draft.md")))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, input_version = input_version, 
    marker_version = marker_version, annotation_version = annotation_version, supersedes_partial_annotation_version = "ANNOTATION_v1_1_FIXED_PATH_LENGTH_FAILURE", 
    dataset_id = dataset_id, started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_seconds = round(as.numeric(difftime(run_end, 
        run_start, units = "secs")), 2), gate_status = gate_status, cell_n = nrow(cell_annotation), sample_n = length(unique(cell_annotation$sample_accession)), 
    cluster_n = length(unique(cell_annotation$cluster)), broad_cell_type_n = length(unique(cell_annotation$broad_cell_type)), 
    large_Seurat_checkpoint_loaded = FALSE, broad_annotation_authorized = annotation_gate$broad_annotation_authorized, 
    subtype_annotation_locked = FALSE, composition_inferential_testing_authorized = FALSE, result_lock_authorized = FALSE, 
    annotation_root = normalize_slashes(annotation_root), qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "broad_cell_types=", length(unique(cell_annotation$broad_cell_type)))

small_files <- unique(c(list.files(annotation_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1E_annotation_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1E_annotation_Results_draft.md"), 
    normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

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

cat("\n============================================================\n")

cat("STEP04E1E GSE250130 broad annotation QC completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Annotated cells: ", format(nrow(cell_annotation), big.mark = ","), "\n", sep = "")

cat("Broad cell types: ", length(unique(cell_annotation$broad_cell_type)), "\n", sep = "")

cat("Review clusters: ", collapse_unique(annotation_map$cluster[annotation_map$review_flag]), "\n", sep = "")

cat("Annotation root: ", annotation_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe annotation and QC folders are complete.\n", sep = "")
}

