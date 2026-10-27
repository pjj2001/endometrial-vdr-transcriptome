#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A05: Fit the original RPKM factorial limma model

options(stringsAsFactors = FALSE)

options(timeout = 7200)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03B1_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

rpkm_filter_threshold <- 1

minimum_samples_above_threshold <- 3L

log_pseudocount <- 1

program_fdr <- 0.050000000000000003

program_abs_logfc <- 0.5

stringent_abs_logfc <- 1

top_variable_features <- 1000L

top_interaction_heatmap <- 30L

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

invisible(NULL)

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive = file.path(module_root, "10_archive_work", run_id))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

qc_figure_dir <- file.path(dirs$figures, "A_QC_only", "STEP03B1_RNA_QC")

candidate_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v1")

for (d in c(qc_figure_dir, file.path(candidate_root, "figures"), file.path(candidate_root, "source_data"), file.path(candidate_root, 
    "code"), file.path(candidate_root, "legends"), file.path(candidate_root, "qc"), file.path(candidate_root, "manifests"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

cran_required <- c("data.table", "ggplot2", "digest", "jsonlite")

cran_missing <- cran_required[!vapply(cran_required, requireNamespace, logical(1), quietly = TRUE)]

if (length(cran_missing) > 0L) {
    install.packages(cran_missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma", ask = FALSE, update = FALSE)
}

required_packages <- c(cran_required, "limma")

if (!all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required packages could not be loaded.")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

md5_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    unname(tools::md5sum(path))
}

sha256_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

normalize_label <- function(x) {
    gsub("[^a-z0-9]+", "", tolower(x))
}

extract_matrix_replicate <- function(x) {
    hit <- regexec("(0[1-3]|[1-3])$", x, perl = TRUE)
    parts <- regmatches(x, hit)[[1]]
    if (length(parts) >= 2L) {
        as.integer(parts[2])
    }
    else {
        NA_integer_
    }
}

infer_matrix_genetic <- function(x) {
    y <- tolower(x)
    if (grepl("sivdr|si.?vdr|vdr.?kd|vdrknock|human.?vdr", y, perl = TRUE)) {
        return("VDR_knockdown")
    }
    if (grepl("sint|si.?nt|nontarget|non.?target|scrambl|control", y, perl = TRUE)) {
        return("siRNA_control")
    }
    NA_character_
}

infer_matrix_ligand <- function(x) {
    y <- tolower(x)
    if (grepl("calcitriol|1.?25|vitamin.?d|ligand|2.?nm", y, perl = TRUE)) {
        return("calcitriol_2nM")
    }
    if (grepl("vehicle|veh|etoh|ethanol", y, perl = TRUE)) {
        return("vehicle")
    }
    NA_character_
}

safe_median <- function(x) {
    median(x, na.rm = TRUE)
}

save_candidate_plot <- function(plot_object, stem, width, height) {
    ggplot2::ggsave(filename = paste0(stem, ".pdf"), plot = plot_object, width = width, height = height, units = "in")
    ggplot2::ggsave(filename = paste0(stem, ".png"), plot = plot_object, width = width, height = height, units = "in", 
        dpi = 300)
}

save_qc_plot <- function(plot_object, path, width, height) {
    ggplot2::ggsave(filename = path, plot = plot_object, width = width, height = height, units = "in", dpi = 200)
}

log_event("Starting", run_id)

gate_path <- file.path(dirs$results, "STEP03A_v5_GATE_TABLE.csv")

design_path <- file.path(dirs$results, "STEP03A_v4_sample_design_LOCKED.csv")

contrast_plan_path <- file.path(dirs$results, "STEP03A_v4_contrast_plan_LOCKED.csv")

download_manifest_path <- file.path(dirs$results, "STEP03A_v5_targeted_download_manifest.csv")

explicit_mapping_path <- file.path(dirs$inputs, "STEP03B1_explicit_matrix_column_to_GSM_lock.csv")

for (p in c(gate_path, design_path, contrast_plan_path, download_manifest_path, explicit_mapping_path)) {
    if (!file.exists(p)) {
        stop("Required input is missing:\n", p)
    }
}

gate <- .vdr_read_csv(gate_path, stringsAsFactors = FALSE, check.names = FALSE)

gate_value <- gate$status[gate$gate_component == "STEP03A v5 gate"]

if (length(gate_value) != 1L || gate_value != "READY_FOR_STEP03B_RNA_ANALYSIS") {
    stop("STEP03A v5 gate does not authorize RNA analysis: ", paste(gate_value, collapse = ", "))
}

locked_design <- .vdr_read_csv(design_path, stringsAsFactors = FALSE, check.names = FALSE)

rna_design <- locked_design[locked_design$accession == "GSE254251", , drop = FALSE]

if (nrow(rna_design) != 12L) {
    stop("Expected 12 locked RNA samples but found ", nrow(rna_design), ".")
}

download_manifest <- .vdr_read_csv(download_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

rna_file_row <- download_manifest[download_manifest$file_role == "processed_expression_matrix", , drop = FALSE]

if (nrow(rna_file_row) != 1L) {
    stop("Expected exactly one processed RNA expression matrix.")
}

expression_path <- rna_file_row$destination[1]

if (!file.exists(expression_path)) {
    stop("Downloaded RPKM matrix is missing:\n", expression_path)
}

if (md5_file(expression_path) != rna_file_row$md5[1]) {
    stop("RPKM matrix MD5 does not match the locked download manifest.")
}

if (sha256_file(expression_path) != rna_file_row$sha256[1]) {
    stop("RPKM matrix SHA-256 does not match the locked download manifest.")
}

raw <- .vdr_fread(expression_path, data.table = FALSE, check.names = FALSE, showProgress = TRUE)

if (ncol(raw) != 13L) {
    write_csv_utf8(data.frame(column_index = seq_along(names(raw)), column_name = names(raw), stringsAsFactors = FALSE), 
        file.path(dirs$qc, "STEP03B1_matrix_header_UNEXPECTED.csv"))
    stop("Expected 13 matrix columns (one feature column plus 12 samples), found ", ncol(raw), ".")
}

feature_column <- names(raw)[1]

matrix_sample_columns <- names(raw)[-1]

explicit_mapping <- .vdr_read_csv(explicit_mapping_path, stringsAsFactors = FALSE, check.names = FALSE)

required_mapping_columns <- c("mapping_order", "matrix_column", "gsm", "genetic_context", "ligand_status_locked", 
    "replicate_locked", "verification_basis", "verification_status")

missing_mapping_columns <- setdiff(required_mapping_columns, names(explicit_mapping))

if (length(missing_mapping_columns) > 0L) {
    stop("Explicit mapping lock is missing columns: ", paste(missing_mapping_columns, collapse = ", "))
}

if (nrow(explicit_mapping) != 12L || anyDuplicated(explicit_mapping$matrix_column) || anyDuplicated(explicit_mapping$gsm)) {
    stop("Explicit mapping lock must contain exactly 12 unique matrix columns ", "and 12 unique GSM accessions.")
}

if (!setequal(explicit_mapping$matrix_column, matrix_sample_columns)) {
    write_csv_utf8(data.frame(matrix_column_observed = sort(matrix_sample_columns), matrix_column_locked = sort(explicit_mapping$matrix_column), 
        stringsAsFactors = FALSE), file.path(dirs$qc, "STEP03B1_explicit_mapping_header_mismatch.csv"))
    stop("The explicit mapping lock does not match the 12 matrix headers.")
}

if (!setequal(explicit_mapping$gsm, rna_design$gsm)) {
    stop("The explicit mapping lock does not contain the same 12 GSM accessions ", "as the STEP03A_v4 locked design.")
}

mapping <- data.frame(matrix_column = matrix_sample_columns, matrix_column_normalized = normalize_label(matrix_sample_columns), 
    inferred_genetic = vapply(matrix_sample_columns, infer_matrix_genetic, character(1)), inferred_ligand = vapply(matrix_sample_columns, 
        infer_matrix_ligand, character(1)), inferred_replicate = vapply(matrix_sample_columns, extract_matrix_replicate, 
        integer(1)), stringsAsFactors = FALSE)

mapping <- merge(mapping, explicit_mapping[c("matrix_column", "gsm", "genetic_context", "ligand_status_locked", 
    "replicate_locked", "verification_basis", "verification_status")], by = "matrix_column", all.x = TRUE, sort = FALSE)

mapping <- mapping[match(matrix_sample_columns, mapping$matrix_column), , drop = FALSE]

mapping <- merge(mapping, rna_design[c("gsm", "genetic_context", "ligand_status_locked", "replicate_locked")], 
    by = "gsm", all.x = TRUE, suffixes = c("_mapping_lock", "_GEO_lock"), sort = FALSE)

mapping <- mapping[match(matrix_sample_columns, mapping$matrix_column), , drop = FALSE]

mapping$header_vs_mapping_genetic <- mapping$inferred_genetic == mapping$genetic_context_mapping_lock

mapping$header_vs_mapping_ligand <- mapping$inferred_ligand == mapping$ligand_status_locked_mapping_lock

mapping$header_vs_mapping_replicate <- mapping$inferred_replicate == mapping$replicate_locked_mapping_lock

mapping$mapping_vs_GEO_genetic <- mapping$genetic_context_mapping_lock == mapping$genetic_context_GEO_lock

mapping$mapping_vs_GEO_ligand <- mapping$ligand_status_locked_mapping_lock == mapping$ligand_status_locked_GEO_lock

mapping$mapping_vs_GEO_replicate <- mapping$replicate_locked_mapping_lock == mapping$replicate_locked_GEO_lock

mapping$mapped_gsm <- mapping$gsm

mapping$mapping_method <- "explicit_LOCKED_column_to_GSM_table"

mapping$mapping_complete <- with(mapping, !is.na(mapped_gsm) & mapped_gsm %in% rna_design$gsm & header_vs_mapping_genetic & 
    header_vs_mapping_ligand & header_vs_mapping_replicate & mapping_vs_GEO_genetic & mapping_vs_GEO_ligand & mapping_vs_GEO_replicate & 
    verification_status == "LOCKED")

write_csv_utf8(mapping, file.path(dirs$qc, "STEP03B1_matrix_column_mapping_audit.csv"))

write_csv_utf8(explicit_mapping, file.path(dirs$results, "STEP03B1_explicit_matrix_column_to_GSM_lock_USED.csv"))

if (!all(mapping$mapping_complete) || anyDuplicated(mapping$mapped_gsm) || !setequal(mapping$mapped_gsm, rna_design$gsm)) {
    write_text(c("# STEP03B1 explicit mapping validation hold", "", "- The explicit matrix-column-to-GSM lock failed one or more checks.", 
        "- Review STEP03B1_matrix_column_mapping_audit.csv.", "- No expression analysis was performed."), file.path(dirs$qc, 
        "STEP03B1_EXPLICIT_MAPPING_VALIDATION_HOLD.md"))
    stop("Explicit sample-column mapping failed validation. ", "No expression analysis was performed.")
}

locked_order <- rna_design$gsm

matrix_col_for_gsm <- mapping$matrix_column[match(locked_order, mapping$mapped_gsm)]

sample_manifest <- merge(rna_design, mapping[c("matrix_column", "mapped_gsm", "mapping_method")], by.x = "gsm", 
    by.y = "mapped_gsm", all.x = TRUE, sort = FALSE)

sample_manifest <- sample_manifest[match(locked_order, sample_manifest$gsm), , drop = FALSE]

sample_manifest$analysis_group <- factor(paste(sample_manifest$genetic_context, sample_manifest$ligand_status_locked, 
    sep = "__"), levels = c("siRNA_control__vehicle", "siRNA_control__calcitriol_2nM", "VDR_knockdown__vehicle", 
    "VDR_knockdown__calcitriol_2nM"))

sample_manifest$short_label <- paste0(ifelse(sample_manifest$genetic_context == "siRNA_control", "NT", "KD"), "_", 
    ifelse(sample_manifest$ligand_status_locked == "vehicle", "VEH", "CAL"), "_", sample_manifest$replicate_locked)

write_csv_utf8(sample_manifest, file.path(dirs$results, "STEP03B1_sample_manifest_FINAL.csv"))

feature_id <- trimws(as.character(raw[[feature_column]]))

expr_raw_df <- raw[matrix_col_for_gsm]

expr_numeric <- suppressWarnings(as.matrix(data.frame(lapply(expr_raw_df, as.numeric), check.names = FALSE)))

colnames(expr_numeric) <- locked_order

rownames(expr_numeric) <- feature_id

non_numeric_count <- sum(is.na(expr_numeric) & !is.na(as.matrix(expr_raw_df)) & nzchar(trimws(as.character(as.matrix(expr_raw_df)))))

if (non_numeric_count > 0L) {
    stop("The expression matrix contains ", non_numeric_count, " nonnumeric nonempty values.")
}

valid_feature <- !is.na(feature_id) & nzchar(feature_id) & apply(expr_numeric, 1, function(x) all(is.finite(x)))

expr_numeric <- expr_numeric[valid_feature, , drop = FALSE]

feature_id <- rownames(expr_numeric)

if (any(expr_numeric < 0)) {
    stop("Negative values were found in the RPKM matrix.")
}

duplicate_feature_n <- sum(duplicated(feature_id))

if (duplicate_feature_n > 0L) {
    feature_counts <- table(feature_id)
    collapsed_sum <- rowsum(expr_numeric, group = feature_id, reorder = FALSE)
    expr_numeric <- collapsed_sum/as.numeric(feature_counts[rownames(collapsed_sum)])
    feature_id <- rownames(expr_numeric)
}

feature_manifest <- data.frame(feature_id = feature_id, identifier_type = ifelse(grepl("^ENSG[0-9]+", feature_id), 
    "Ensembl_gene_id", ifelse(grepl("^[A-Za-z0-9.-]+$", feature_id), "gene_symbol_or_compact_identifier", "other")), 
    mean_RPKM = rowMeans(expr_numeric), median_RPKM = apply(expr_numeric, 1, median), max_RPKM = apply(expr_numeric, 
        1, max), n_samples_RPKM_gt_1 = rowSums(expr_numeric > rpkm_filter_threshold), stringsAsFactors = FALSE)

keep <- rowSums(expr_numeric > rpkm_filter_threshold) >= minimum_samples_above_threshold

expr_filtered <- expr_numeric[keep, , drop = FALSE]

log_expr <- log2(expr_filtered + log_pseudocount)

write_csv_utf8(feature_manifest, file.path(dirs$results, "STEP03B1_feature_manifest.csv"))

saveRDS(expr_numeric, file.path(dirs$results, "STEP03B1_RPKM_matrix_all_features.rds"), compress = "xz")

saveRDS(log_expr, file.path(dirs$results, "STEP03B1_log2_RPKM_filtered_matrix.rds"), compress = "xz")

sample_qc <- data.frame(gsm = colnames(expr_numeric), short_label = sample_manifest$short_label, group = as.character(sample_manifest$analysis_group), 
    total_RPKM = colSums(expr_numeric), median_RPKM = apply(expr_numeric, 2, median), zero_fraction = colMeans(expr_numeric == 
        0), detected_RPKM_gt_1 = colSums(expr_numeric > rpkm_filter_threshold), stringsAsFactors = FALSE)

sample_cor <- cor(log_expr, method = "pearson")

sample_qc$median_sample_correlation <- vapply(seq_len(ncol(sample_cor)), function(i) {
    median(sample_cor[i, -i], na.rm = TRUE)
}, numeric(1))

cor_median <- median(sample_qc$median_sample_correlation)

cor_mad <- mad(sample_qc$median_sample_correlation, constant = 1)

sample_qc$correlation_outlier_flag <- if (is.finite(cor_mad) && cor_mad > 0) {
    sample_qc$median_sample_correlation < cor_median - 3 * cor_mad
} else {
    FALSE
}

write_csv_utf8(sample_qc, file.path(dirs$qc, "STEP03B1_sample_QC_metrics.csv"))

cor_long <- as.data.frame(as.table(sample_cor), stringsAsFactors = FALSE)

names(cor_long) <- c("sample_x", "sample_y", "correlation")

cor_long$sample_x <- factor(cor_long$sample_x, levels = locked_order, labels = sample_manifest$short_label)

cor_long$sample_y <- factor(cor_long$sample_y, levels = rev(locked_order), labels = rev(sample_manifest$short_label))

write_csv_utf8(cor_long, file.path(dirs$qc, "SourceData_STEP03B1_sample_correlation.csv"))

p_cor <- ggplot2::ggplot(cor_long, ggplot2::aes(x = sample_x, y = sample_y, fill = correlation)) + ggplot2::geom_tile() + 
    ggplot2::coord_equal() + ggplot2::labs(title = "GSE254251 sample correlation", x = NULL, y = NULL, fill = "Pearson r") + 
    ggplot2::theme_bw(base_size = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

save_qc_plot(p_cor, file.path(qc_figure_dir, "QC_STEP03B1_sample_correlation.png"), width = 7.5, height = 6.5)

distribution_df <- data.frame(value = as.vector(log_expr), gsm = rep(colnames(log_expr), each = nrow(log_expr)), 
    stringsAsFactors = FALSE)

distribution_df$short_label <- sample_manifest$short_label[match(distribution_df$gsm, sample_manifest$gsm)]

distribution_df$short_label <- factor(distribution_df$short_label, levels = sample_manifest$short_label)

p_distribution <- ggplot2::ggplot(distribution_df, ggplot2::aes(x = short_label, y = value)) + ggplot2::geom_boxplot(outlier.shape = NA) + 
    ggplot2::labs(title = "GSE254251 log2(RPKM + 1) distributions", x = NULL, y = "log2(RPKM + 1)") + ggplot2::theme_bw(base_size = 10) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

save_qc_plot(p_distribution, file.path(qc_figure_dir, "QC_STEP03B1_expression_distributions.png"), width = 8.1999999999999993, 
    height = 5.5)

feature_variance <- apply(log_expr, 1, var)

pca_n <- min(top_variable_features, sum(is.finite(feature_variance)))

pca_features <- names(sort(feature_variance, decreasing = TRUE))[seq_len(pca_n)]

pca <- prcomp(t(log_expr[pca_features, , drop = FALSE]), center = TRUE, scale. = FALSE)

pca_variance <- (pca$sdev^2/sum(pca$sdev^2)) * 100

pca_df <- data.frame(gsm = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)

pca_df <- merge(pca_df, sample_manifest[c("gsm", "short_label", "genetic_context", "ligand_status_locked", "replicate_locked", 
    "analysis_group")], by = "gsm", all.x = TRUE, sort = FALSE)

pca_df <- pca_df[match(locked_order, pca_df$gsm), , drop = FALSE]

write_csv_utf8(pca_df, file.path(candidate_root, "source_data", "SourceData_STEP03B1_PCA.csv"))

p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = analysis_group, shape = genetic_context)) + 
    ggplot2::geom_point(size = 3.2000000000000002, alpha = 0.84999999999999998) + ggplot2::geom_text(ggplot2::aes(label = short_label), 
    vjust = -0.80000000000000004, size = 3, show.legend = FALSE) + ggplot2::labs(title = "GSE254251 RNA-seq sample structure", 
    subtitle = paste0("Top ", pca_n, " variable expressed features"), x = paste0("PC1 (", round(pca_variance[1], 
        1), "%)"), y = paste0("PC2 (", round(pca_variance[2], 1), "%)"), color = "Condition", shape = "siRNA") + 
    ggplot2::theme_bw(base_size = 11) + ggplot2::theme(plot.title.position = "plot")

save_candidate_plot(p_pca, file.path(candidate_root, "figures", "Figure_STEP03B1_PCA_candidate_v1"), width = 8, 
    height = 6)

group <- sample_manifest$analysis_group

design <- model.matrix(~0 + group)

colnames(design) <- c("NT_VEH", "NT_CAL", "KD_VEH", "KD_CAL")

rownames(design) <- sample_manifest$gsm

if (qr(design)$rank != ncol(design)) {
    stop("RNA design matrix is rank deficient.")
}

contrast_matrix <- limma::makeContrasts(RNA_C1_ligand_in_siNT = NT_CAL - NT_VEH, RNA_C2_ligand_in_siVDR = KD_CAL - 
    KD_VEH, RNA_C3_VDRKD_in_vehicle = KD_VEH - NT_VEH, RNA_C4_VDRKD_in_ligand = KD_CAL - NT_CAL, RNA_C5_siRNA_by_ligand_interaction = (KD_CAL - 
    KD_VEH) - (NT_CAL - NT_VEH), levels = design)

fit <- limma::lmFit(log_expr, design)

fit2 <- limma::contrasts.fit(fit, contrast_matrix)

fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)

saveRDS(fit2, file.path(dirs$results, "STEP03B1_limma_fit_factorial.rds"), compress = "xz")

contrast_names <- colnames(contrast_matrix)

all_tables <- list()

for (contrast_name in contrast_names) {
    tab <- limma::topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P")
    tab$feature_id <- rownames(tab)
    rownames(tab) <- NULL
    tab$contrast_id <- contrast_name
    tab$program_candidate <- (tab$adj.P.Val < program_fdr & abs(tab$logFC) >= program_abs_logfc)
    tab$stringent_DE <- (tab$adj.P.Val < program_fdr & abs(tab$logFC) >= stringent_abs_logfc)
    tab$direction <- ifelse(tab$logFC > 0, "up", ifelse(tab$logFC < 0, "down", "zero"))
    tab <- tab[c("feature_id", "contrast_id", setdiff(names(tab), c("feature_id", "contrast_id")))]
    all_tables[[contrast_name]] <- tab
    write_csv_utf8(tab, file.path(dirs$results, paste0("STEP03B1_limma_", contrast_name, "_all_features.csv")))
}

all_results <- do.call(rbind, all_tables)

rownames(all_results) <- NULL

write_csv_utf8(all_results, file.path(dirs$results, "STEP03B1_limma_all_contrasts_long.csv"))

deg_summary <- do.call(rbind, lapply(all_tables, function(tab) {
    data.frame(contrast_id = unique(tab$contrast_id), n_tested = nrow(tab), fdr_only = sum(tab$adj.P.Val < program_fdr, 
        na.rm = TRUE), program_up = sum(tab$program_candidate & tab$logFC > 0, na.rm = TRUE), program_down = sum(tab$program_candidate & 
        tab$logFC < 0, na.rm = TRUE), stringent_up = sum(tab$stringent_DE & tab$logFC > 0, na.rm = TRUE), stringent_down = sum(tab$stringent_DE & 
        tab$logFC < 0, na.rm = TRUE), stringsAsFactors = FALSE)
}))

write_csv_utf8(deg_summary, file.path(dirs$results, "STEP03B1_contrast_DE_summary.csv"))

program_candidates <- all_results[all_results$program_candidate, c("feature_id", "contrast_id", "logFC", "AveExpr", 
    "P.Value", "adj.P.Val", "direction"), drop = FALSE]

write_csv_utf8(program_candidates, file.path(dirs$results, "STEP03B1_program_gene_candidates_NOT_FROZEN.csv"))

vdr_indices <- which(toupper(gsub("\\.[0-9]+$", "", rownames(log_expr))) %in% c("VDR", "ENSG00000111424"))

vdr_found <- length(vdr_indices) >= 1L

if (vdr_found) {
    vdr_index <- vdr_indices[which.max(rowMeans(log_expr[vdr_indices, , drop = FALSE]))]
    vdr_df <- data.frame(gsm = colnames(log_expr), log2_RPKM_plus_1 = as.numeric(log_expr[vdr_index, ]), stringsAsFactors = FALSE)
    vdr_df <- merge(vdr_df, sample_manifest[c("gsm", "short_label", "genetic_context", "ligand_status_locked", 
        "analysis_group")], by = "gsm", all.x = TRUE, sort = FALSE)
    vdr_df <- vdr_df[match(locked_order, vdr_df$gsm), , drop = FALSE]
    vdr_group_means <- aggregate(log2_RPKM_plus_1 ~ genetic_context + ligand_status_locked, data = vdr_df, FUN = mean)
    nt_vehicle <- vdr_group_means$log2_RPKM_plus_1[vdr_group_means$genetic_context == "siRNA_control" & vdr_group_means$ligand_status_locked == 
        "vehicle"]
    kd_vehicle <- vdr_group_means$log2_RPKM_plus_1[vdr_group_means$genetic_context == "VDR_knockdown" & vdr_group_means$ligand_status_locked == 
        "vehicle"]
    nt_ligand <- vdr_group_means$log2_RPKM_plus_1[vdr_group_means$genetic_context == "siRNA_control" & vdr_group_means$ligand_status_locked == 
        "calcitriol_2nM"]
    kd_ligand <- vdr_group_means$log2_RPKM_plus_1[vdr_group_means$genetic_context == "VDR_knockdown" & vdr_group_means$ligand_status_locked == 
        "calcitriol_2nM"]
    vdr_direction_vehicle <- length(nt_vehicle) == 1L && length(kd_vehicle) == 1L && kd_vehicle < nt_vehicle
    vdr_direction_ligand <- length(nt_ligand) == 1L && length(kd_ligand) == 1L && kd_ligand < nt_ligand
    vdr_direction_pass <- vdr_direction_vehicle && vdr_direction_ligand
    write_csv_utf8(vdr_df, file.path(candidate_root, "source_data", "SourceData_STEP03B1_VDR_expression.csv"))
    write_csv_utf8(vdr_group_means, file.path(dirs$qc, "STEP03B1_VDR_group_means.csv"))
    p_vdr <- ggplot2::ggplot(vdr_df, ggplot2::aes(x = analysis_group, y = log2_RPKM_plus_1, shape = genetic_context)) + 
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55000000000000004) + ggplot2::geom_jitter(width = 0.080000000000000002, 
        height = 0, size = 2.7999999999999998) + ggplot2::labs(title = "VDR expression in the GSE254251 perturbation design", 
        x = NULL, y = "VDR log2(RPKM + 1)", shape = "siRNA") + ggplot2::theme_bw(base_size = 11) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, 
        hjust = 1), plot.title.position = "plot")
    save_candidate_plot(p_vdr, file.path(candidate_root, "figures", "Figure_STEP03B1_VDR_expression_candidate_v1"), 
        width = 7.5, height = 5.5)
} else {
    vdr_direction_vehicle <- NA
    vdr_direction_ligand <- NA
    vdr_direction_pass <- NA
    write_text(c("# VDR identifier review required", "", "- No exact VDR or ENSG00000111424 feature was found after filtering.", 
        "- Inspect the matrix feature identifier type.", "- Do not interpret knockdown efficacy until VDR is resolved."), 
        file.path(dirs$qc, "STEP03B1_VDR_IDENTIFIER_REVIEW.md"))
}

deg_plot_df <- deg_summary[rep(seq_len(nrow(deg_summary)), each = 2L), c("contrast_id", "program_up", "program_down")]

deg_plot_df$direction <- rep(c("up", "down"), times = nrow(deg_summary))

deg_plot_df$n_genes <- c(rbind(deg_summary$program_up, deg_summary$program_down))

deg_plot_df$signed_n <- ifelse(deg_plot_df$direction == "up", deg_plot_df$n_genes, -deg_plot_df$n_genes)

write_csv_utf8(deg_plot_df, file.path(candidate_root, "source_data", "SourceData_STEP03B1_contrast_summary.csv"))

p_deg_summary <- ggplot2::ggplot(deg_plot_df, ggplot2::aes(x = contrast_id, y = signed_n, fill = direction)) + 
    ggplot2::geom_col() + ggplot2::geom_hline(yintercept = 0, linewidth = 0.40000000000000002) + ggplot2::coord_flip() + 
    ggplot2::labs(title = "Candidate VDR and ligand-response gene sets", subtitle = paste0("FDR < ", program_fdr, 
        " and absolute log2FC >= ", program_abs_logfc), x = NULL, y = "Candidate genes (down shown below zero)", 
        fill = "Direction") + ggplot2::theme_bw(base_size = 11) + ggplot2::theme(plot.title.position = "plot")

save_candidate_plot(p_deg_summary, file.path(candidate_root, "figures", "Figure_STEP03B1_contrast_summary_candidate_v1"), 
    width = 8, height = 5.7999999999999998)

c1 <- all_tables[["RNA_C1_ligand_in_siNT"]][, c("feature_id", "logFC", "adj.P.Val"), drop = FALSE]

names(c1) <- c("feature_id", "ligand_logFC_siNT", "ligand_FDR_siNT")

c2 <- all_tables[["RNA_C2_ligand_in_siVDR"]][, c("feature_id", "logFC", "adj.P.Val"), drop = FALSE]

names(c2) <- c("feature_id", "ligand_logFC_siVDR", "ligand_FDR_siVDR")

c5 <- all_tables[["RNA_C5_siRNA_by_ligand_interaction"]][, c("feature_id", "logFC", "adj.P.Val"), drop = FALSE]

names(c5) <- c("feature_id", "interaction_logFC", "interaction_FDR")

interaction_df <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = FALSE), list(c1, c2, c5))

interaction_df$interaction_category <- ifelse(interaction_df$interaction_FDR < program_fdr & interaction_df$interaction_logFC > 
    0, "positive interaction", ifelse(interaction_df$interaction_FDR < program_fdr & interaction_df$interaction_logFC < 
    0, "negative interaction", "not significant"))

write_csv_utf8(interaction_df, file.path(candidate_root, "source_data", "SourceData_STEP03B1_ligand_response_interaction.csv"))

p_interaction <- ggplot2::ggplot(interaction_df, ggplot2::aes(x = ligand_logFC_siNT, y = ligand_logFC_siVDR, color = interaction_category)) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.5) + ggplot2::geom_hline(yintercept = 0, 
    linewidth = 0.29999999999999999) + ggplot2::geom_vline(xintercept = 0, linewidth = 0.29999999999999999) + ggplot2::geom_point(alpha = 0.55000000000000004, 
    size = 1.3999999999999999) + ggplot2::labs(title = "Ligand-response effects with and without VDR knockdown", 
    subtitle = "Interaction categories are based on the formal siRNA-by-ligand contrast", x = "Ligand log2FC in siRNA control", 
    y = "Ligand log2FC after VDR knockdown", color = "Interaction") + ggplot2::theme_bw(base_size = 11) + ggplot2::theme(plot.title.position = "plot")

save_candidate_plot(p_interaction, file.path(candidate_root, "figures", "Figure_STEP03B1_ligand_interaction_candidate_v1"), 
    width = 7.2999999999999998, height = 6.2000000000000002)

interaction_rank <- all_tables[["RNA_C5_siRNA_by_ligand_interaction"]]

interaction_rank <- interaction_rank[order(interaction_rank$adj.P.Val, -abs(interaction_rank$logFC)), , drop = FALSE]

heatmap_n <- min(top_interaction_heatmap, nrow(interaction_rank))

heatmap_features <- interaction_rank$feature_id[seq_len(heatmap_n)]

heatmap_matrix <- log_expr[heatmap_features, , drop = FALSE]

row_center <- rowMeans(heatmap_matrix)

row_sd <- apply(heatmap_matrix, 1, sd)

row_sd[!is.finite(row_sd) | row_sd == 0] <- 1

heatmap_z <- (heatmap_matrix - row_center)/row_sd

heatmap_long <- data.frame(feature_id = rep(rownames(heatmap_z), times = ncol(heatmap_z)), gsm = rep(colnames(heatmap_z), 
    each = nrow(heatmap_z)), z_score = as.vector(heatmap_z), stringsAsFactors = FALSE)

heatmap_long$short_label <- sample_manifest$short_label[match(heatmap_long$gsm, sample_manifest$gsm)]

heatmap_long$feature_id <- factor(heatmap_long$feature_id, levels = rev(heatmap_features))

heatmap_long$short_label <- factor(heatmap_long$short_label, levels = sample_manifest$short_label)

write_csv_utf8(heatmap_long, file.path(candidate_root, "source_data", "SourceData_STEP03B1_top_interaction_heatmap.csv"))

p_heatmap <- ggplot2::ggplot(heatmap_long, ggplot2::aes(x = short_label, y = feature_id, fill = z_score)) + ggplot2::geom_tile() + 
    ggplot2::labs(title = "Top siRNA-by-ligand interaction features", subtitle = paste0("Top ", heatmap_n, " features ranked by interaction FDR"), 
        x = NULL, y = NULL, fill = "Row z-score") + ggplot2::theme_bw(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1), plot.title.position = "plot")

save_candidate_plot(p_heatmap, file.path(candidate_root, "figures", "Figure_STEP03B1_top_interaction_heatmap_candidate_v1"), 
    width = 8.5, height = 8)

qc_summary <- data.frame(metric = c("raw_feature_rows", "valid_feature_rows", "duplicate_feature_ids_collapsed", 
    "filtered_feature_rows", "samples", "matrix_mapping_complete", "design_rank", "design_columns", "sample_correlation_outliers", 
    "VDR_identifier_found", "VDR_lower_after_knockdown_vehicle", "VDR_lower_after_knockdown_ligand"), value = c(nrow(raw), 
    nrow(expr_numeric), duplicate_feature_n, nrow(log_expr), ncol(log_expr), all(mapping$mapping_complete), qr(design)$rank, 
    ncol(design), sum(sample_qc$correlation_outlier_flag), vdr_found, vdr_direction_vehicle, vdr_direction_ligand), 
    stringsAsFactors = FALSE)

write_csv_utf8(qc_summary, file.path(dirs$qc, "STEP03B1_QC_summary.csv"))

structural_pass <- (ncol(log_expr) == 12L && all(mapping$mapping_complete) && qr(design)$rank == ncol(design) && 
    nrow(log_expr) > 1000L)

biological_qc_status <- if (!vdr_found) {
    "REVIEW_VDR_IDENTIFIER"
} else if (isTRUE(vdr_direction_pass)) {
    "PASS"
} else {
    "HOLD_VDR_KNOCKDOWN_DIRECTION"
}

automated_gate <- if (structural_pass && biological_qc_status == "PASS") {
    "PROCEED_TO_MANUAL_RESULT_AND_FIGURE_QC"
} else if (structural_pass && biological_qc_status == "REVIEW_VDR_IDENTIFIER") {
    "HOLD_VDR_IDENTIFIER_REVIEW"
} else if (structural_pass) {
    "HOLD_BIOLOGICAL_PERTURBATION_QC"
} else {
    "HOLD_MATRIX_OR_DESIGN_QC"
}

gate_table <- data.frame(gate_component = c("STEP03A v5 acquisition gate", "Matrix column mapping", "Expression matrix structure", 
    "Filtering", "Design rank", "Sample-correlation outlier flag", "VDR identifier", "VDR knockdown direction", 
    "Program candidates", "Figure class", "STEP03B1 gate"), status = c("PASS", ifelse(all(mapping$mapping_complete), 
    "PASS", "FAIL"), paste0(nrow(expr_numeric), " features x ", ncol(expr_numeric), " samples"), paste0(nrow(log_expr), 
    " retained"), paste0(qr(design)$rank, "/", ncol(design)), paste0("N=", sum(sample_qc$correlation_outlier_flag)), 
    ifelse(vdr_found, "PASS", "REVIEW"), biological_qc_status, "NOT_FROZEN", "A_QC_ONLY_AND_B_SUBMISSION_CANDIDATE", 
    automated_gate), required_action = c("Retain the acquisition manifest.", "Do not use column order as a silent fallback.", 
    "Review feature identifier and matrix dimensions.", "Retain the prespecified RPKM filter.", "Do not add unsupported covariates.", 
    "Review any flagged sample before exclusion.", "Resolve feature identifier if absent.", "Confirm effective perturbation before program freezing.", 
    "Freeze programs only after content/statistical/visual QC.", "Do not promote figures to class C.", "Upload results and actual rendered figures for review."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B1_GATE_TABLE.csv"))

methods_lines <- c("# STEP03B1 Methods draft", "", "The GSE254251 processed RPKM matrix was verified against the locked download manifest using MD5 and SHA-256 hashes. The 12 expression columns were mapped to GEO samples using GSM accessions embedded in the matrix header or unique treatment-by-replicate parsing; matrix-column order was not used as an implicit fallback.", 
    "", paste0("Features were retained when RPKM exceeded ", rpkm_filter_threshold, " in at least ", minimum_samples_above_threshold, 
        " samples. Values were transformed as log2(RPKM + ", log_pseudocount, "). A no-intercept four-group limma model represented siRNA-control/vehicle, siRNA-control/calcitriol, VDR-knockdown/vehicle and VDR-knockdown/calcitriol conditions. Empirical Bayes moderation used trend=TRUE and robust=TRUE."), 
    "", "Five contrasts were prespecified: ligand response under siRNA control, ligand response after VDR knockdown, VDR-knockdown effect under vehicle, VDR-knockdown effect under ligand and the formal siRNA-by-ligand interaction. P values were adjusted by the Benjamini-Hochberg method.")

results_lines <- c("# STEP03B1 Results draft", "", paste0("The matrix contained ", nrow(expr_numeric), " valid features across 12 samples; ", 
    nrow(log_expr), " features passed the prespecified expression filter."), "", paste0("No. of correlation-outlier flags: ", 
    sum(sample_qc$correlation_outlier_flag), "."), "", paste0("Biological VDR-knockdown QC status: ", biological_qc_status, 
    "."), "", "Candidate program counts by contrast are reported in STEP03B1_contrast_DE_summary.csv. These gene sets remain provisional until scientific, statistical and visual QC are completed.")

discussion_lines <- c("# STEP03B1 Discussion draft", "", "The factorial design separates ligand response, VDR-knockdown effects and their interaction. The interaction contrast is essential because a gene can respond to calcitriol in both siRNA backgrounds while still showing a quantitatively altered response after VDR knockdown.", 
    "", "The processed matrix contains RPKM values rather than raw counts. Accordingly, the analysis uses log-transformed continuous-expression modeling and does not claim count-model differential expression. Candidate programs are experimental T-HESC perturbation programs and should be interpreted as mechanistic anchors for later patient-cohort projection, not as decidualization programs.")

limitation_lines <- c("# STEP03B1 Limitation note", "", "- Three biological wells per condition provide limited precision.", 
    "- RPKM is a processed continuous-expression measure, not raw count data.", "- No sequencing sample was experimentally decidualized.", 
    "- Results are derived from the T-HESC cell model.", "- Candidate programs are not frozen until cross-omics and visual QC.")

transition_lines <- c("# STEP03B1 Transition note", "", if (automated_gate == "PROCEED_TO_MANUAL_RESULT_AND_FIGURE_QC") {
    "Proceed to manual review of all five contrasts and the actual candidate figure renders. After QC, freeze the RNA program definitions and begin STEP03B2 ATAC-seq analysis."
} else {
    paste0("Do not proceed to program freezing. Resolve gate: ", automated_gate, ".")
})

write_text(methods_lines, file.path(dirs$text, "STEP03B1_Methods_draft.md"))

write_text(results_lines, file.path(dirs$text, "STEP03B1_Results_draft.md"))

write_text(discussion_lines, file.path(dirs$text, "STEP03B1_Discussion_draft.md"))

write_text(limitation_lines, file.path(dirs$text, "STEP03B1_Limitation_note.md"))

write_text(transition_lines, file.path(dirs$text, "STEP03B1_Transition_to_next_substep.md"))

legend_lines <- c("# STEP03B1 candidate figure legends", "", "**PCA candidate.** Principal-component analysis of the 12 GSE254251 samples using the most variable expressed features after the prespecified RPKM filter. Point color indicates the four siRNA-by-ligand conditions and point shape indicates siRNA status.", 
    "", "**VDR-expression candidate.** VDR log2(RPKM + 1) values across the four experimental conditions. Individual points represent biological wells.", 
    "", "**Contrast-summary candidate.** Numbers of provisional up- and downregulated genes for the five prespecified contrasts using FDR < 0.05 and absolute log2 fold change >= 0.5. These counts are descriptive and the programs are not yet frozen.", 
    "", "**Ligand-interaction candidate.** Comparison of ligand-response log2 fold changes under siRNA control and VDR knockdown. Interaction categories are defined by the formal siRNA-by-ligand interaction FDR and sign.", 
    "", "**Interaction-heatmap candidate.** Row-standardized expression of the top features ranked by the siRNA-by-ligand interaction FDR.")

write_text(legend_lines, file.path(candidate_root, "legends", "STEP03B1_candidate_legends.md"))

blueprint_lines <- c("# STEP03B1 Figure Blueprint Card", "", "- Figure class: B_submission_candidate", "- Final manuscript role: NOT_LOCKED", 
    "- Scientific question: How do VDR knockdown and calcitriol jointly alter the T-HESC transcriptome?", "- Authoritative source: STEP03B1 limma tables and source-data CSV files.", 
    "- Statistical unit: biological well.", "- Main candidate panels: VDR expression; ligand-response interaction; top interaction heatmap.", 
    "- Supporting candidate panel: PCA.", "- Descriptive candidate panel: contrast gene-count summary.", "- Allowed claim: experimental VDR and ligand-response programs were identified in a 2 x 2 T-HESC perturbation design.", 
    "- Prohibited claim: these data directly define decidualization or prove patient-disease causality.", "- Final composite layout: DEFERRED.", 
    "- Cross-format QC: DEFERRED UNTIL CLASS C.", "- Manuscript-context QC: DEFERRED UNTIL CLASS C.")

write_text(blueprint_lines, file.path(candidate_root, "FIGURE_BLUEPRINT_CARD.md"))

candidate_qc_lines <- c("# STEP03B1 candidate QC status", "", paste0("- Analysis gate: ", automated_gate), "- Scientific QC: PENDING MANUAL REVIEW", 
    "- Statistical QC: PENDING MANUAL REVIEW", "- Actual-render QC: PENDING UPLOAD AND REVIEW", "- Source data: PRESENT", 
    "- Generating code: PRESENT", "- Draft legends: PRESENT", "- Submission-ready: NO", "- Program freeze: NOT PERMITTED YET")

write_text(candidate_qc_lines, file.path(candidate_root, "qc", "CANDIDATE_QC_STATUS.md"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "fit_experimental_rpkm_model.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
    file.copy(source_file, file.path(candidate_root, "code", "fit_experimental_rpkm_model.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
}

candidate_files <- list.files(candidate_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

candidate_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", candidate_files), nchar(gsub("\\\\", 
    "/", candidate_root)) + 2L), size_bytes = as.numeric(file.info(candidate_files)$size), md5 = vapply(candidate_files, 
    md5_file, character(1)), sha256 = vapply(candidate_files, sha256_file, character(1)), figure_class = "B_submission_candidate", 
    release_status = "NOT_SUBMISSION_READY", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest, file.path(candidate_root, "manifests", "STEP03B1_candidate_package_manifest.csv"))

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03B1", 
    run_id = run_id, input_accession = "GSE254251", input_type = "RPKM", rpkm_filter_threshold = rpkm_filter_threshold, 
    minimum_samples_above_threshold = minimum_samples_above_threshold, transformation = paste0("log2(RPKM+", log_pseudocount, 
        ")"), limma = list(trend = TRUE, robust = TRUE), program_threshold = list(FDR = program_fdr, absolute_logFC = program_abs_logfc), 
    stringent_threshold = list(FDR = program_fdr, absolute_logFC = stringent_abs_logfc), candidate_programs_frozen = FALSE, 
    figure_classes = c("A_QC_only", "B_submission_candidate"), class_C_release = FALSE)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B1_parameters_v0.1.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B1.txt"))

for (d in c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", "06_results", 
    "07_figures_tables", "08_qc", "09_manuscript_text")) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_groups <- list(`01_inputs` = c(gate_path, design_path, contrast_plan_path, download_manifest_path, explicit_mapping_path, 
    file.path(dirs$results, "STEP03B1_explicit_matrix_column_to_GSM_lock_USED.csv"), file.path(dirs$results, "STEP03B1_sample_manifest_FINAL.csv")), 
    `02_code` = list.files(dirs$code, pattern = "09_step03B1_v2_", full.names = TRUE), `03_parameters` = list.files(dirs$parameters, 
        pattern = "^STEP03B1_", full.names = TRUE), `04_environment` = list.files(dirs$environment, pattern = "STEP03B1", 
        full.names = TRUE), `05_execution` = log_path, `06_results` = list.files(dirs$results, pattern = "^STEP03B1_", 
        full.names = TRUE), `08_qc` = list.files(dirs$qc, pattern = "^STEP03B1_|^SourceData_STEP03B1_", full.names = TRUE), 
    `09_manuscript_text` = list.files(dirs$text, pattern = "^STEP03B1_", full.names = TRUE))

for (subdir in names(copy_groups)) {
    paths <- copy_groups[[subdir]]
    paths <- paths[file.exists(paths)]
    for (p in paths) {
        file.copy(p, file.path(dirs$archive, subdir, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

archive_qc_fig <- file.path(dirs$archive, "07_figures_tables", "A_QC_only")

archive_candidate <- file.path(dirs$archive, "07_figures_tables", "B_submission_candidate")

dir.create(archive_qc_fig, recursive = TRUE, showWarnings = FALSE)

dir.create(archive_candidate, recursive = TRUE, showWarnings = FALSE)

qc_files <- list.files(qc_figure_dir, full.names = TRUE)

for (p in qc_files) {
    file.copy(p, file.path(archive_qc_fig, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

copy_candidate_files <- list.files(candidate_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

for (p in copy_candidate_files) {
    rel <- substring(gsub("\\\\", "/", p), nchar(gsub("\\\\", "/", candidate_root)) + 2L)
    destination <- file.path(archive_candidate, rel)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    file.copy(p, destination, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03B1_archive_manifest.csv"))

write_text(c("# STEP03B1 archive status", "", paste0("- Run ID: ", run_id), paste0("- Gate: ", automated_gate), 
    "- Matrix and sample mapping: CAPTURED", "- Limma result tables: CAPTURED", "- Candidate programs: NOT_FROZEN", 
    "- A-class QC figures: CAPTURED", "- B-class candidate figures: CAPTURED_NOT_RELEASED", "- Scientific/statistical/visual QC: PENDING", 
    "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B1", paste0("- Run ID: ", run_id), paste0("- Gate: ", 
    automated_gate), "- RNA-seq factorial analysis: COMPLETED_PENDING_QC", "- RNA program candidates: NOT_FROZEN", 
    "- ATAC analysis: NOT_STARTED", "- CUT&RUN analysis: NOT_STARTED", "- Figure candidates: B_CLASS_PENDING_RENDER_QC", 
    "- Full module status: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03B1 v3 RNA factorial analysis completed.\n")

cat("Gate: ", automated_gate, "\n")

cat("Samples: ", ncol(log_expr), "\n")

cat("Filtered features: ", nrow(log_expr), "\n")

cat("VDR biological QC: ", biological_qc_status, "\n")

cat("Candidate programs remain NOT FROZEN.\n")

cat("B-class figures require actual-render QC.\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B1_GATE_TABLE.csv\n")

cat("  06_results/STEP03B1_contrast_DE_summary.csv\n")

cat("  06_results/STEP03B1_program_gene_candidates_NOT_FROZEN.csv\n")

cat("  08_qc/STEP03B1_QC_summary.csv\n")

cat("  08_qc/STEP03B1_sample_QC_metrics.csv\n")

cat("  08_qc/STEP03B1_matrix_column_mapping_audit.csv\n")

cat("  07_figures_tables/A_QC_only/STEP03B1_RNA_QC/\n")

cat("  07_figures_tables/B_submission_candidate/STEP03B1_RNA_candidates_v1/\n")

cat("============================================================\n")

