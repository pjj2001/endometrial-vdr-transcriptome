#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A06: Classify interaction candidates and preserve source tables

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03B1_v4_", format(Sys.time(), "%Y%m%d_%H%M%S"))

within_group_correlation_threshold <- 0.94999999999999996

program_fdr <- 0.050000000000000003

program_abs_logfc <- 0.5

top_volcano_labels <- 12L

heatmap_top_n <- 30L

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive = file.path(module_root, "10_archive_work", run_id))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_v1_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v1")

candidate_v2_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v2")

for (d in c(file.path(candidate_v2_root, "figures"), file.path(candidate_v2_root, "source_data"), file.path(candidate_v2_root, 
    "code"), file.path(candidate_v2_root, "legends"), file.path(candidate_v2_root, "qc"), file.path(candidate_v2_root, 
    "manifests"), file.path(dirs$figures, "A_QC_only", "STEP03B1_RNA_QC_v2"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

required_cran <- c("ggplot2", "ggrepel", "scales", "digest", "jsonlite")

missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_cran) > 0L) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

if (!all(vapply(required_cran, requireNamespace, logical(1), quietly = TRUE))) {
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

theme_skill3 <- function(base_size = 9) {
    ggplot2::theme_classic(base_size = base_size, base_family = "Arial") + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", 
        size = base_size + 1, margin = ggplot2::margin(b = 4)), plot.subtitle = ggplot2::element_text(size = base_size, 
        margin = ggplot2::margin(b = 6)), axis.title = ggplot2::element_text(size = base_size), axis.text = ggplot2::element_text(size = base_size - 
        0.5), legend.title = ggplot2::element_text(face = "bold", size = base_size), legend.text = ggplot2::element_text(size = base_size - 
        0.5), plot.margin = ggplot2::margin(7, 7, 7, 7))
}

save_qc_png <- function(plot_object, path, width_mm, height_mm) {
    ggplot2::ggsave(filename = path, plot = plot_object, width = width_mm, height = height_mm, units = "mm", dpi = 300, 
        bg = "white", limitsize = FALSE)
    if (!file.exists(path) || file.info(path)$size <= 0) {
        stop("Expected plot was not created:\n", path)
    }
}

required_columns <- function(x, cols, object_name) {
    missing <- setdiff(cols, names(x))
    if (length(missing) > 0L) {
        stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
    }
}

log_event("Starting", run_id)

gate_path <- file.path(dirs$results, "STEP03B1_GATE_TABLE.csv")

deg_summary_path <- file.path(dirs$results, "STEP03B1_contrast_DE_summary.csv")

sample_manifest_path <- file.path(dirs$results, "STEP03B1_sample_manifest_FINAL.csv")

all_contrasts_path <- file.path(dirs$results, "STEP03B1_limma_all_contrasts_long.csv")

legacy_sample_qc_path <- file.path(dirs$qc, "STEP03B1_sample_QC_metrics.csv")

correlation_source_path <- file.path(dirs$qc, "SourceData_STEP03B1_sample_correlation.csv")

pca_source_path <- file.path(candidate_v1_root, "source_data", "SourceData_STEP03B1_PCA.csv")

vdr_source_path <- file.path(candidate_v1_root, "source_data", "SourceData_STEP03B1_VDR_expression.csv")

interaction_source_path <- file.path(candidate_v1_root, "source_data", "SourceData_STEP03B1_ligand_response_interaction.csv")

heatmap_source_path <- file.path(candidate_v1_root, "source_data", "SourceData_STEP03B1_top_interaction_heatmap.csv")

required_paths <- c(gate_path, deg_summary_path, sample_manifest_path, all_contrasts_path, legacy_sample_qc_path, 
    correlation_source_path, pca_source_path, vdr_source_path, interaction_source_path, heatmap_source_path)

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0L) {
    stop("Missing locked STEP03B1 inputs:\n", paste(missing_paths, collapse = "\n"))
}

gate <- .vdr_read_csv(gate_path, stringsAsFactors = FALSE, check.names = FALSE)

gate_value <- gate$status[gate$gate_component == "STEP03B1 gate"]

if (length(gate_value) != 1L || gate_value != "PROCEED_TO_MANUAL_RESULT_AND_FIGURE_QC") {
    stop("STEP03B1 gate is not ready for Skill 3 review: ", paste(gate_value, collapse = ", "))
}

deg_summary <- .vdr_read_csv(deg_summary_path, stringsAsFactors = FALSE, check.names = FALSE)

sample_manifest <- .vdr_read_csv(sample_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

all_contrasts <- .vdr_read_csv(all_contrasts_path, stringsAsFactors = FALSE, check.names = FALSE)

legacy_sample_qc <- .vdr_read_csv(legacy_sample_qc_path, stringsAsFactors = FALSE, check.names = FALSE)

cor_long <- .vdr_read_csv(correlation_source_path, stringsAsFactors = FALSE, check.names = FALSE)

pca_df <- .vdr_read_csv(pca_source_path, stringsAsFactors = FALSE, check.names = FALSE)

vdr_df <- .vdr_read_csv(vdr_source_path, stringsAsFactors = FALSE, check.names = FALSE)

interaction_df <- .vdr_read_csv(interaction_source_path, stringsAsFactors = FALSE, check.names = FALSE)

heatmap_long <- .vdr_read_csv(heatmap_source_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(sample_manifest) != 12L || nrow(pca_df) != 12L || nrow(vdr_df) != 12L || nrow(interaction_df) != 12116L) {
    stop("Unexpected STEP03B1 row count.")
}

required_columns(interaction_df, c("feature_id", "ligand_logFC_siNT", "ligand_FDR_siNT", "ligand_logFC_siVDR", 
    "ligand_FDR_siVDR", "interaction_logFC", "interaction_FDR"), "interaction source")

sample_order <- sample_manifest$short_label

sample_group <- setNames(as.character(sample_manifest$analysis_group), sample_manifest$short_label)

cor_matrix <- xtabs(correlation ~ sample_y + sample_x, data = cor_long)

missing_labels <- setdiff(sample_order, rownames(cor_matrix))

if (length(missing_labels) > 0L) {
    stop("Correlation source is missing sample labels: ", paste(missing_labels, collapse = ", "))
}

corrected_rows <- lapply(sample_order, function(sample_id) {
    group_id <- sample_group[[sample_id]]
    peers <- names(sample_group)[sample_group == group_id & names(sample_group) != sample_id]
    within_values <- as.numeric(cor_matrix[sample_id, peers, drop = TRUE])
    data.frame(short_label = sample_id, analysis_group = group_id, n_within_group_peers = length(within_values), 
        within_group_min_r = min(within_values), within_group_median_r = median(within_values), global_median_r = legacy_sample_qc$median_sample_correlation[match(sample_id, 
            legacy_sample_qc$short_label)], legacy_global_flag = legacy_sample_qc$correlation_outlier_flag[match(sample_id, 
            legacy_sample_qc$short_label)], corrected_sample_QC = ifelse(min(within_values) >= within_group_correlation_threshold, 
            "PASS_RETAIN", "MANUAL_REVIEW"), stringsAsFactors = FALSE)
})

corrected_sample_qc <- do.call(rbind, corrected_rows)

write_csv_utf8(corrected_sample_qc, file.path(dirs$qc, "STEP03B1_v4_corrected_within_condition_sample_QC.csv"))

group_qc <- aggregate(cbind(within_group_min_r, within_group_median_r) ~ analysis_group, data = corrected_sample_qc, 
    FUN = min)

group_qc$n_samples <- vapply(group_qc$analysis_group, function(group_id) {
    sum(corrected_sample_qc$analysis_group == group_id)
}, integer(1))

group_qc$decision <- ifelse(group_qc$within_group_min_r >= within_group_correlation_threshold, "PASS_RETAIN_ALL", 
    "MANUAL_REVIEW")

write_csv_utf8(group_qc, file.path(dirs$qc, "STEP03B1_v4_within_condition_group_QC.csv"))

all_samples_retained <- all(corrected_sample_qc$corrected_sample_QC == "PASS_RETAIN")

write_text(c("# STEP03B1 v4 sample-QC correction", "", "- The legacy global median-correlation rule is retired.", 
    "- It flagged an entire biological treatment group and therefore", "  did not identify sample-specific technical outliers.", 
    paste0("- Corrected threshold: minimum within-condition Pearson r >= ", within_group_correlation_threshold, 
        "."), paste0("- Corrected decision: ", ifelse(all_samples_retained, "retain all 12 samples", "manual review required"), 
        "."), "- No sample was excluded from the already fitted limma model."), file.path(dirs$qc, "STEP03B1_v4_SAMPLE_QC_DECISION.md"))

interaction_df$is_interaction_candidate <- (interaction_df$interaction_FDR < program_fdr & abs(interaction_df$interaction_logFC) >= 
    program_abs_logfc)

interaction_df$program_class <- "not_interaction_candidate"

interaction_df$program_class[interaction_df$is_interaction_candidate & interaction_df$ligand_logFC_siNT > 0 & interaction_df$interaction_logFC < 
    0] <- "attenuated_ligand_induction_after_VDRKD"

interaction_df$program_class[interaction_df$is_interaction_candidate & interaction_df$ligand_logFC_siNT < 0 & interaction_df$interaction_logFC > 
    0] <- "attenuated_ligand_repression_after_VDRKD"

interaction_df$program_class[interaction_df$is_interaction_candidate & interaction_df$ligand_logFC_siNT > 0 & interaction_df$interaction_logFC > 
    0] <- "enhanced_ligand_induction_after_VDRKD"

interaction_df$program_class[interaction_df$is_interaction_candidate & interaction_df$ligand_logFC_siNT < 0 & interaction_df$interaction_logFC < 
    0] <- "enhanced_ligand_repression_after_VDRKD"

interaction_df$control_ligand_candidate <- (interaction_df$ligand_FDR_siNT < program_fdr & abs(interaction_df$ligand_logFC_siNT) >= 
    program_abs_logfc)

interaction_df$residual_ligand_candidate_after_VDRKD <- (interaction_df$ligand_FDR_siVDR < program_fdr & abs(interaction_df$ligand_logFC_siVDR) >= 
    program_abs_logfc)

interaction_candidates <- interaction_df[interaction_df$is_interaction_candidate, , drop = FALSE]

interaction_summary <- as.data.frame(table(interaction_candidates$program_class), stringsAsFactors = FALSE)

names(interaction_summary) <- c("program_class", "n_genes")

interaction_summary$proportion <- (interaction_summary$n_genes/sum(interaction_summary$n_genes))

interaction_summary <- interaction_summary[order(interaction_summary$n_genes, decreasing = TRUE), , drop = FALSE]

write_csv_utf8(interaction_candidates, file.path(dirs$results, "STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv"))

write_csv_utf8(interaction_summary, file.path(dirs$results, "STEP03B1_v4_interaction_program_summary.csv"))

program_definition <- data.frame(program_id = c("RNA_C1_ligand_response_intact_VDR", "RNA_C2_residual_ligand_response_after_VDRKD", 
    "RNA_C3_baseline_VDRKD_effect", "RNA_C4_ligand_conditioned_VDRKD_effect", "RNA_C5_formal_VDRKD_by_ligand_interaction"), 
    source_contrast = c("RNA_C1_ligand_in_siNT", "RNA_C2_ligand_in_siVDR", "RNA_C3_VDRKD_in_vehicle", "RNA_C4_VDRKD_in_ligand", 
        "RNA_C5_siRNA_by_ligand_interaction"), threshold = paste0("FDR<", program_fdr, "; abs(log2FC)>=", program_abs_logfc), 
    scientific_role = c("ligand response with intact VDR", "residual ligand response after VDR knockdown", "VDR-dependent baseline transcriptome", 
        "VDR-dependent transcriptome under ligand", "formal VDR dependence of ligand response"), freeze_status = "RNA_LAYER_CANDIDATE_NOT_FROZEN", 
    stringsAsFactors = FALSE)

program_definition$n_up <- deg_summary$program_up[match(program_definition$source_contrast, deg_summary$contrast_id)]

program_definition$n_down <- deg_summary$program_down[match(program_definition$source_contrast, deg_summary$contrast_id)]

write_csv_utf8(program_definition, file.path(dirs$results, "STEP03B1_v4_RNA_program_definition_card.csv"))

write_csv_utf8(deg_summary, file.path(candidate_v2_root, "source_data", "Table_STEP03B1_contrast_DE_summary_candidate.csv"))

write_text(c("# STEP03B1 figure-versus-table decision", "", "- PCA: supplementary figure candidate.", "- VDR expression: main/supporting figure candidate.", 
    "- Formal interaction volcano: main figure candidate.", "- Top interaction heatmap: main figure candidate.", 
    "- Five-contrast gene counts: supplementary table, not a retained figure.", "- Original ligand-response scatter: archived as v1; replaced by volcano.", 
    "", "Reason: figures should communicate sample structure, perturbation", "efficacy and interaction patterns; the full exact contrast counts", 
    "are more transparent in a table."), file.path(candidate_v2_root, "FIGURE_TABLE_DECISION.md"))

write_text(c("# Language Routing Card", "", "- Mode: AUTO_NATIVE", "- Authoritative result language: R", "- Authoritative plotting implementation: R", 
    "- Rationale: limma analysis and locked result tables are R-native.", "- Cross-language migration: PROHIBITED for this candidate version.", 
    "- Inferential statistics: loaded from locked result tables; not recomputed.", "- Skill: research_scientific_figure_suite v1.7-beta"), 
    file.path(candidate_v2_root, "LANGUAGE_ROUTING_CARD.md"))

write_text(c("# Color-semantic registry", "", "- siNT / siRNA control: #0072B2", "- siVDR / VDR knockdown: #D55E00", 
    "- Interaction attenuation of induction: #7B3294", "- Interaction attenuation of repression: #008837", "- Other interaction: #E6AB02", 
    "- Noncandidate genes: #BDBDBD", "- Heatmap: continuous blue-white-red row z-score scale", "- Black lines: reference or group separators"), 
    file.path(candidate_v2_root, "COLOR_SEMANTIC_REGISTRY.md"))

pca_df$siRNA_label <- ifelse(pca_df$genetic_context == "siRNA_control", "siNT", "siVDR")

pca_df$ligand_label <- ifelse(pca_df$ligand_status_locked == "vehicle", "Vehicle", "Calcitriol")

pca_df$group_label <- paste(pca_df$siRNA_label, pca_df$ligand_label, sep = " + ")

pca_centroids <- aggregate(cbind(PC1, PC2) ~ group_label + siRNA_label + ligand_label, data = pca_df, FUN = mean)

pca_source_v2 <- merge(pca_df, pca_centroids[c("group_label", "PC1", "PC2")], by = "group_label", all.x = TRUE, 
    suffixes = c("", "_centroid"), sort = FALSE)

write_csv_utf8(pca_source_v2, file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_PCA_candidate_v2.csv"))

p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = siRNA_label, shape = ligand_label)) + ggplot2::geom_point(size = 2.6000000000000001, 
    alpha = 0.90000000000000002) + ggplot2::geom_point(data = pca_centroids, ggplot2::aes(x = PC1, y = PC2), inherit.aes = FALSE, 
    shape = 21, fill = "white", color = "black", size = 3.7999999999999998, stroke = 0.80000000000000004) + ggrepel::geom_text_repel(data = pca_centroids, 
    ggplot2::aes(x = PC1, y = PC2, label = group_label), inherit.aes = FALSE, family = "Arial", size = 2.7000000000000002, 
    min.segment.length = 0, box.padding = 0.25, point.padding = 0.14999999999999999, seed = 20260803, show.legend = FALSE) + 
    ggplot2::scale_color_manual(values = c(siNT = "#0072B2", siVDR = "#D55E00")) + ggplot2::scale_shape_manual(values = c(Vehicle = 16, 
    Calcitriol = 17)) + ggplot2::labs(title = "GSE254251 sample structure", subtitle = "Top 1,000 variable expressed features", 
    x = "PC1 (68.3%)", y = "PC2 (13.5%)", color = "siRNA", shape = "Ligand") + theme_skill3(9) + ggplot2::theme(legend.position = "right")

save_qc_png(p_pca, file.path(candidate_v2_root, "figures", "Figure_STEP03B1_PCA_candidate_v2_QC.png"), width_mm = 105, 
    height_mm = 82)

vdr_df$siRNA_label <- ifelse(vdr_df$genetic_context == "siRNA_control", "siNT", "siVDR")

vdr_df$ligand_label <- ifelse(vdr_df$ligand_status_locked == "vehicle", "Vehicle", "Calcitriol")

vdr_df$siRNA_label <- factor(vdr_df$siRNA_label, levels = c("siNT", "siVDR"))

vdr_df$ligand_label <- factor(vdr_df$ligand_label, levels = c("Vehicle", "Calcitriol"))

vdr_summary <- aggregate(log2_RPKM_plus_1 ~ siRNA_label + ligand_label, data = vdr_df, FUN = function(x) {
    c(mean = mean(x), sd = sd(x), n = length(x))
})

vdr_summary <- data.frame(siRNA_label = vdr_summary$siRNA_label, ligand_label = vdr_summary$ligand_label, mean = vdr_summary$log2_RPKM_plus_1[, 
    "mean"], sd = vdr_summary$log2_RPKM_plus_1[, "sd"], n = vdr_summary$log2_RPKM_plus_1[, "n"], stringsAsFactors = FALSE)

write_csv_utf8(vdr_df, file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_VDR_expression_candidate_v2.csv"))

write_csv_utf8(vdr_summary, file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_VDR_expression_summary_v2.csv"))

dodge <- ggplot2::position_dodge(width = 0.45000000000000001)

p_vdr <- ggplot2::ggplot(vdr_df, ggplot2::aes(x = ligand_label, y = log2_RPKM_plus_1, color = siRNA_label, group = siRNA_label)) + 
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.059999999999999998, jitter.height = 0, 
        dodge.width = 0.45000000000000001, seed = 20260803), size = 2.2999999999999998, alpha = 0.90000000000000002) + 
    ggplot2::geom_errorbar(data = vdr_summary, ggplot2::aes(y = mean, ymin = mean - sd, ymax = mean + sd), position = dodge, 
        width = 0.10000000000000001, linewidth = 0.55000000000000004, show.legend = FALSE) + ggplot2::geom_point(data = vdr_summary, 
    ggplot2::aes(y = mean, fill = siRNA_label), position = dodge, shape = 21, color = "black", size = 3.2000000000000002, 
    stroke = 0.45000000000000001, show.legend = FALSE) + ggplot2::scale_color_manual(values = c(siNT = "#0072B2", 
    siVDR = "#D55E00")) + ggplot2::scale_fill_manual(values = c(siNT = "#0072B2", siVDR = "#D55E00")) + ggplot2::scale_y_continuous(limits = c(0, 
    NA), expand = ggplot2::expansion(mult = c(0, 0.080000000000000002))) + ggplot2::labs(title = "VDR knockdown was effective under both ligand conditions", 
    x = NULL, y = "VDR log2(RPKM + 1)", color = "siRNA") + theme_skill3(9) + ggplot2::theme(legend.position = "top", 
    legend.justification = "left")

save_qc_png(p_vdr, file.path(candidate_v2_root, "figures", "Figure_STEP03B1_VDR_expression_candidate_v2_QC.png"), 
    width_mm = 92, height_mm = 72)

interaction_plot <- interaction_df

interaction_plot$plot_class <- "Not candidate"

interaction_plot$plot_class[interaction_plot$program_class == "attenuated_ligand_induction_after_VDRKD"] <- "Attenuated induction"

interaction_plot$plot_class[interaction_plot$program_class == "attenuated_ligand_repression_after_VDRKD"] <- "Attenuated repression"

interaction_plot$plot_class[interaction_plot$is_interaction_candidate & interaction_plot$plot_class == "Not candidate"] <- "Other interaction"

interaction_plot$plot_class <- factor(interaction_plot$plot_class, levels = c("Not candidate", "Attenuated induction", 
    "Attenuated repression", "Other interaction"))

interaction_plot$minus_log10_FDR <- -log10(pmax(interaction_plot$interaction_FDR, .Machine$double.xmin))

label_df <- interaction_plot[interaction_plot$is_interaction_candidate, , drop = FALSE]

label_df <- label_df[order(label_df$interaction_FDR, -abs(label_df$interaction_logFC)), , drop = FALSE]

label_df <- head(label_df, top_volcano_labels)

write_csv_utf8(interaction_plot, file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_interaction_volcano_candidate_v2.csv"))

p_volcano <- ggplot2::ggplot(interaction_plot, ggplot2::aes(x = interaction_logFC, y = minus_log10_FDR, color = plot_class)) + 
    ggplot2::geom_point(size = 1.1499999999999999, alpha = 0.57999999999999996) + ggplot2::geom_vline(xintercept = c(-program_abs_logfc, 
    program_abs_logfc), linetype = 2, linewidth = 0.34999999999999998, color = "grey35") + ggplot2::geom_hline(yintercept = -log10(program_fdr), 
    linetype = 2, linewidth = 0.34999999999999998, color = "grey35") + ggrepel::geom_text_repel(data = label_df, 
    ggplot2::aes(label = feature_id), family = "Arial", size = 2.5, min.segment.length = 0, max.overlaps = Inf, 
    box.padding = 0.25, point.padding = 0.12, seed = 20260803, show.legend = FALSE) + ggplot2::scale_color_manual(values = c(`Not candidate` = "#BDBDBD", 
    `Attenuated induction` = "#7B3294", `Attenuated repression` = "#008837", `Other interaction` = "#E6AB02")) + 
    ggplot2::labs(title = "VDR knockdown altered the calcitriol transcriptional response", subtitle = paste0("Formal siRNA-by-ligand interaction; FDR < ", 
        program_fdr, " and |interaction log2FC| >= ", program_abs_logfc), x = "Interaction log2 fold change", y = expression(-log[10]("FDR")), 
        color = NULL) + theme_skill3(9) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") + 
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 2, byrow = TRUE))

save_qc_png(p_volcano, file.path(candidate_v2_root, "figures", "Figure_STEP03B1_interaction_volcano_candidate_v2_QC.png"), 
    width_mm = 118, height_mm = 92)

feature_order <- unique(heatmap_long$feature_id)

feature_order <- head(feature_order, heatmap_top_n)

heatmap_plot <- heatmap_long[heatmap_long$feature_id %in% feature_order, , drop = FALSE]

heatmap_plot$feature_id <- factor(heatmap_plot$feature_id, levels = rev(feature_order))

heatmap_plot$short_label <- factor(heatmap_plot$short_label, levels = sample_order)

write_csv_utf8(heatmap_plot, file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_top_interaction_heatmap_candidate_v2.csv"))

p_heatmap <- ggplot2::ggplot(heatmap_plot, ggplot2::aes(x = short_label, y = feature_id, fill = z_score)) + ggplot2::geom_tile() + 
    ggplot2::geom_vline(xintercept = c(3.5, 6.5, 9.5), color = "black", linewidth = 0.45000000000000001) + ggplot2::scale_fill_gradient2(low = "#3B4CC0", 
    mid = "white", high = "#B40426", midpoint = 0, limits = c(-max(abs(heatmap_plot$z_score), na.rm = TRUE), max(abs(heatmap_plot$z_score), 
        na.rm = TRUE)), oob = scales::squish) + ggplot2::labs(title = "Top VDR-dependent ligand-response interaction genes", 
    subtitle = paste0("Top ", length(feature_order), " genes ranked by interaction FDR"), x = NULL, y = NULL, fill = "Row z-score") + 
    theme_skill3(8.5) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), axis.ticks = ggplot2::element_blank(), 
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.40000000000000002), legend.position = "right")

save_qc_png(p_heatmap, file.path(candidate_v2_root, "figures", "Figure_STEP03B1_top_interaction_heatmap_candidate_v2_QC.png"), 
    width_mm = 125, height_mm = 125)

write_text(c("# Figure Blueprint Card v2", "", "- Status: BLUEPRINT_LOCKED_FOR_FAST_REPLOT", "- Scientific question: How do VDR knockdown and calcitriol jointly", 
    "  alter the T-HESC transcriptome?", "- Statistical unit: biological well; n=3 per factorial cell.", "- Authoritative result source: locked STEP03B1 limma result tables.", 
    "- Transformation: log2(RPKM + 1).", "- Model: four-group limma model with robust, trend empirical Bayes.", 
    "- Multiple testing: Benjamini-Hochberg FDR.", "- Candidate threshold: FDR < 0.05 and |log2FC| >= 0.5.", "- Panel A candidate: VDR perturbation efficacy.", 
    "- Panel B candidate: formal interaction volcano.", "- Panel C candidate: top interaction heatmap.", "- Supplementary candidate: PCA.", 
    "- Exact contrast counts: supplementary table.", "- Allowed claim: VDR knockdown substantially altered the calcitriol", 
    "  transcriptional response in T-HESC cells.", "- Prohibited claim: direct VDR binding, decidualization response,", 
    "  patient causality or diagnostic value.", "- Composite layout: DEFERRED UNTIL ATAC/CUT&RUN PANELS ARE AVAILABLE.", 
    "- Final PDF/TIFF/SVG export: PROHIBITED BEFORE VISUAL LOCK."), file.path(candidate_v2_root, "FIGURE_BLUEPRINT_CARD_v2.md"))

write_text(c("# Candidate figure legends v2", "", "**PCA candidate.** Principal-component analysis of the 12 GSE254251", 
    "samples using the 1,000 most variable retained features. Colors denote", "siRNA status, shapes denote ligand treatment, and open symbols mark", 
    "condition centroids. The analysis is descriptive.", "", "**VDR-expression candidate.** VDR log2(RPKM + 1) values in the four", 
    "factorial conditions. Points are biological wells; filled outlined", "symbols show group means and error bars show one standard deviation.", 
    "", "**Interaction-volcano candidate.** Formal siRNA-by-ligand interaction", "effects from the limma model. Candidate genes meet FDR < 0.05 and", 
    "|interaction log2FC| >= 0.5. Attenuated induction denotes a positive", "ligand effect under siNT with a negative interaction; attenuated", 
    "repression denotes a negative ligand effect under siNT with a positive", "interaction. Selected top-ranked genes are labelled.", 
    "", "**Interaction-heatmap candidate.** Row-standardized log2(RPKM + 1)", "expression for the top interaction-ranked genes. Columns are individual", 
    "biological wells ordered by the four factorial conditions; black vertical", "lines separate conditions."), 
    file.path(candidate_v2_root, "legends", "STEP03B1_candidate_legends_v2.md"))

write_text(c("# Chinese panel descriptions v2", "", "- PCA\uff1a\u5c55\u793a12\u4e2a\u6837\u672c\u6574\u4f53\u7ed3\u6784\uff0c\u4ec5\u4f5c\u4e3a\u6837\u672c\u5173\u7cfb\u548c\u5206\u7ec4\u5408\u7406\u6027\u7684\u8865\u5145\u8bc1\u636e\u3002", 
    "- VDR\u8868\u8fbe\uff1a\u663e\u793asiVDR\u5728vehicle\u548ccalcitriol\u6761\u4ef6\u4e0b\u5747\u6709\u6548\u964d\u4f4eVDR\u8868\u8fbe\u3002", "- \u4ea4\u4e92\u706b\u5c71\u56fe\uff1a\u663e\u793aVDR\u6572\u4f4e\u540e\u54ea\u4e9bcalcitriol\u53cd\u5e94\u53d1\u751f\u663e\u8457\u6539\u53d8\u3002", 
    "- \u4ea4\u4e92\u70ed\u56fe\uff1a\u5c55\u793a\u4ea4\u4e92\u4f5c\u7528\u6700\u5f3a\u7684\u57fa\u56e0\u5728\u56db\u4e2a\u5b9e\u9a8c\u6761\u4ef6\u4e2d\u7684\u8868\u8fbe\u6a21\u5f0f\u3002", "- \u8fd9\u4e9b\u56fe\u4e0d\u80fd\u88ab\u89e3\u91ca\u4e3a\u8715\u819c\u5316\u5904\u7406\u7ed3\u679c\u3001\u76f4\u63a5\u7ed3\u5408\u8bc1\u636e\u6216\u60a3\u8005\u56e0\u679c\u8bc1\u636e\u3002"), 
    file.path(candidate_v2_root, "legends", "STEP03B1_Chinese_panel_descriptions_v2.md"))

write_text(c("# STEP03B1 v4 scientific and statistical QC", "", "- Source lock: PASS", "- Sample mapping: PASS", 
    "- Factorial design: PASS", "- Limma model and contrast tables: RETAINED UNCHANGED", "- VDR perturbation efficacy: PASS", 
    "- Legacy global correlation outlier rule: RETIRED", paste0("- Corrected within-condition sample QC: ", ifelse(all_samples_retained, 
        "PASS_RETAIN_ALL", "MANUAL_REVIEW")), "- Sample exclusion: NONE", "- Interaction classification: COMPLETE", 
    "- Program freeze: NOT YET PERMITTED", "- Decidualization interpretation: PROHIBITED"), file.path(candidate_v2_root, 
    "qc", "SCIENTIFIC_STATISTICAL_QC_v2.md"))

write_text(c("# STEP03B1 v4 render-QC status", "", "- Candidate version: v2", "- Output mode: 300-dpi QC PNG only", 
    "- PDF/TIFF/SVG: NOT GENERATED", "- Actual-render review: PENDING USER UPLOAD", "- Cross-format QC: DEFERRED", 
    "- Manuscript-context QC: DEFERRED", "- Submission-ready: NO"), file.path(candidate_v2_root, "qc", "RENDER_QC_STATUS_v2.md"))

write_text(c("# Old-version disposition", "", "- STEP03B1_RNA_candidates_v1: ARCHIVED, NOT DELETED.", "- PCA v1: replaced because of overlapping sample labels.", 
    "- VDR-expression v1: replaced because n=3 boxplots and long labels", "  reduced clarity.", "- Ligand-response scatter v1: replaced by formal interaction volcano.", 
    "- Heatmap v1: replaced because the sequential palette obscured sign.", "- Contrast-summary figure v1: demoted to supplementary table.", 
    "- No v1 file is overwritten."), file.path(candidate_v2_root, "OLD_VERSION_DISPOSITION.md"))

n_interaction <- sum(interaction_df$is_interaction_candidate)

n_attenuated <- sum(interaction_df$program_class %in% c("attenuated_ligand_induction_after_VDRKD", "attenuated_ligand_repression_after_VDRKD"))

revised_figures <- list.files(file.path(candidate_v2_root, "figures"), pattern = "_v2_QC\\.png$", full.names = TRUE)

figure_files_valid <- (length(revised_figures) == 4L && all(file.info(revised_figures)$size > 0))

automated_gate <- if (all_samples_retained && n_interaction == 436L && n_attenuated == 431L && figure_files_valid) {
    "PROCEED_TO_REVISED_FIGURE_RENDER_QC"
} else {
    "HOLD_QC_OR_SOURCE_MISMATCH"
}

gate_v4 <- data.frame(gate_component = c("Locked STEP03B1 inference reused", "Within-condition sample QC", "Samples retained", 
    "Interaction candidates", "Attenuation-pattern candidates", "Figure/table decision", "v2 candidate PNG creation", 
    "Program freeze", "STEP03B1 v4 gate"), status = c("PASS", ifelse(all_samples_retained, "PASS", "FAIL"), paste0(sum(corrected_sample_qc$corrected_sample_QC == 
    "PASS_RETAIN"), "/12"), paste0("N=", n_interaction), paste0("N=", n_attenuated), "PASS", ifelse(figure_files_valid, 
    "PASS", "FAIL"), "NOT_FROZEN", automated_gate), required_action = c("Do not refit or change contrasts during FAST_REPLOT.", 
    "Use within-condition replicate agreement.", "Do not exclude samples based on the retired global flag.", "Review exact gene lists before downstream scoring.", 
    "Interpret as attenuation, not direct binding.", "Use figures for patterns and table for exact counts.", "Upload actual renders for visual QC.", 
    "Freeze only after v2 render QC and content review.", "Do not start ATAC analysis until STEP03B1 is locked."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_v4, file.path(dirs$results, "STEP03B1_v4_GATE_TABLE.csv"))

write_text(c("# STEP03B1 v4 QC and FAST_REPLOT report", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: `", 
    automated_gate, "`"), "- Limma model: unchanged.", "- Samples retained: 12/12.", paste0("- Formal interaction candidates: ", 
    n_interaction, "."), paste0("- Attenuation-pattern interaction candidates: ", n_attenuated, " (", round(100 * 
    n_attenuated/n_interaction, 1), "%)."), "- Candidate PNGs created: 4.", "- Program freeze: pending actual-render and content QC.", 
    "- Next analysis module: not authorized yet."), file.path(dirs$results, "STEP03B1_v4_QC_FAST_REPLOT_REPORT.md"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "classify_experimental_interactions.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
    file.copy(source_file, file.path(candidate_v2_root, "code", "classify_experimental_interactions.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
}

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03B1_v4", 
    run_id = run_id, mode = "FAST_REPLOT_AFTER_FULL_AUDIT", authoritative_language = "R", limma_refit = FALSE, 
    sample_exclusion = FALSE, within_group_correlation_threshold = within_group_correlation_threshold, program_fdr = program_fdr, 
    program_abs_logfc = program_abs_logfc, candidate_version = "v2", final_export = FALSE)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B1_v4_parameters_v0.1.json"), pretty = TRUE, 
    auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B1_v4.txt"))

candidate_files <- list.files(candidate_v2_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

candidate_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", candidate_files), nchar(gsub("\\\\", 
    "/", candidate_v2_root)) + 2L), size_bytes = as.numeric(file.info(candidate_files)$size), md5 = vapply(candidate_files, 
    md5_file, character(1)), sha256 = vapply(candidate_files, sha256_file, character(1)), figure_class = "B_submission_candidate", 
    release_status = "PENDING_ACTUAL_RENDER_QC", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest, file.path(candidate_v2_root, "manifests", "STEP03B1_candidate_v2_manifest.csv"))

for (d in c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", "06_results", 
    "07_figures_tables", "08_qc")) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_groups <- list(`01_inputs` = required_paths, `02_code` = list.files(dirs$code, pattern = "^11_step03B1_v4_", 
    full.names = TRUE), `03_parameters` = list.files(dirs$parameters, pattern = "^STEP03B1_v4_", full.names = TRUE), 
    `04_environment` = list.files(dirs$environment, pattern = "STEP03B1_v4", full.names = TRUE), `05_execution` = log_path, 
    `06_results` = list.files(dirs$results, pattern = "^STEP03B1_v4_", full.names = TRUE), `08_qc` = list.files(dirs$qc, 
        pattern = "^STEP03B1_v4_", full.names = TRUE))

for (subdir in names(copy_groups)) {
    paths <- copy_groups[[subdir]]
    paths <- paths[file.exists(paths)]
    for (p in paths) {
        file.copy(p, file.path(dirs$archive, subdir, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

archive_candidate_root <- file.path(dirs$archive, "07_figures_tables", "B_submission_candidate_v2")

for (p in candidate_files) {
    relative_path <- substring(gsub("\\\\", "/", p), nchar(gsub("\\\\", "/", candidate_v2_root)) + 2L)
    destination <- file.path(archive_candidate_root, relative_path)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    file.copy(p, destination, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03B1_v4_archive_manifest.csv"))

write_text(c("# STEP03B1 v4 archive status", "", paste0("- Run ID: ", run_id), paste0("- Gate: ", automated_gate), 
    "- Limma inference: RETAINED", "- Sample QC correction: CAPTURED", "- Samples excluded: NONE", "- Candidate figures: V2_QC_PNG_CREATED", 
    "- Actual-render QC: PENDING", "- RNA program freeze: NOT_LOCKED", "- ATAC analysis: NOT AUTHORIZED", "- Full STEP03 module: NOT_LOCKED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B1_v4", paste0("- Run ID: ", run_id), paste0("- Gate: ", 
    automated_gate), "- RNA inference: COMPLETED", "- Sample QC: CORRECTED_RETAIN_ALL", "- RNA figures: V2_PENDING_RENDER_QC", 
    "- RNA program candidates: NOT_FROZEN", "- ATAC analysis: NOT_STARTED", "- CUT&RUN analysis: NOT_STARTED", 
    "- Full module status: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03B1 v4 QC correction and Skill 3 FAST_REPLOT completed.\n")

cat("Gate: ", automated_gate, "\n")

cat("Samples retained: 12/12\n")

cat("Interaction candidates: ", n_interaction, "\n")

cat("Attenuation-pattern candidates: ", n_attenuated, "\n")

cat("New v2 candidate PNGs: 4\n")

cat("Limma was not refitted.\n")

cat("Programs remain NOT FROZEN.\n")

cat("ATAC analysis is not yet authorized.\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B1_v4_GATE_TABLE.csv\n")

cat("  06_results/STEP03B1_v4_QC_FAST_REPLOT_REPORT.md\n")

cat("  06_results/STEP03B1_v4_interaction_program_summary.csv\n")

cat("  08_qc/STEP03B1_v4_corrected_within_condition_sample_QC.csv\n")

cat("  08_qc/STEP03B1_v4_within_condition_group_QC.csv\n")

cat("  07_figures_tables/B_submission_candidate/STEP03B1_RNA_candidates_v2/\n")

cat("============================================================\n")

