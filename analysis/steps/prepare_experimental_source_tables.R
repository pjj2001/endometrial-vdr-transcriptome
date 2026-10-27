#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A06b: Prepare schema-consistent experimental source tables

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B1_v7_", format(Sys.time(), "%Y%m%d_%H%M%S"))

program_fdr <- 0.050000000000000003

program_abs_logfc <- 0.5

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

dirs <- list(control = file.path(module_root, "00_module_control"), code = file.path(module_root, "02_code"), parameters = file.path(module_root, 
    "03_parameters"), environment = file.path(module_root, "04_environment"), execution = file.path(module_root, 
    "05_execution"), results = file.path(module_root, "06_results"), figures = file.path(module_root, "07_figures_tables"), 
    qc = file.path(module_root, "08_qc"), archive = file.path(module_root, "10_archive_work", run_id))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_v2_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v2")

candidate_v3_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v3")

for (d in c(file.path(candidate_v3_root, "figures"), file.path(candidate_v3_root, "source_data"), file.path(candidate_v3_root, 
    "source_data", "upstream_refs"), file.path(candidate_v3_root, "code"), file.path(candidate_v3_root, "legends"), 
    file.path(candidate_v3_root, "qc"), file.path(candidate_v3_root, "manifests"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(candidate_v2_root)) {
    stop("STEP03B1 v2 candidate package not found:\n", candidate_v2_root)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

required_cran <- c("ggplot2", "ggrepel", "digest", "jsonlite")

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
        size = base_size + 1, margin = ggplot2::margin(b = 3)), plot.subtitle = ggplot2::element_text(size = base_size, 
        lineheight = 1.05, margin = ggplot2::margin(b = 6)), axis.title = ggplot2::element_text(size = base_size), 
        axis.text = ggplot2::element_text(size = base_size - 0.5), legend.title = ggplot2::element_text(face = "bold", 
            size = base_size), legend.text = ggplot2::element_text(size = base_size - 0.5), plot.margin = ggplot2::margin(7, 
            9, 7, 7))
}

save_qc_png <- function(plot_object, path, width_mm, height_mm) {
    ggplot2::ggsave(filename = path, plot = plot_object, width = width_mm, height = height_mm, units = "mm", dpi = 300, 
        bg = "white")
}

copy_validate <- function(source, destination) {
    if (!file.exists(source)) {
        stop("Source file missing:\n", source)
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    copied <- file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    if (!isTRUE(copied) || !file.exists(destination) || md5_file(source) != md5_file(destination) || sha256_file(source) != 
        sha256_file(destination)) {
        stop("Copy/hash validation failed:\n", source, "\n->\n", destination)
    }
    TRUE
}

log_event("Starting", run_id)

v4_gate_path <- file.path(dirs$results, "STEP03B1_v4_GATE_TABLE.csv")

upstream_result_paths <- c(interaction_program_summary = file.path(dirs$results, "STEP03B1_v4_interaction_program_summary.csv"), 
    program_definition_card = file.path(dirs$results, "STEP03B1_v4_RNA_program_definition_card.csv"), interaction_candidates = file.path(dirs$results, 
        "STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv"), sample_manifest = file.path(dirs$results, 
        "STEP03B1_sample_manifest_FINAL.csv"), limma_all_contrasts = file.path(dirs$results, "STEP03B1_limma_all_contrasts_long.csv"))

required_source_paths <- c(PCA = file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_PCA_candidate_v2.csv"), 
    VDR = file.path(candidate_v2_root, "source_data", "SourceData_STEP03B1_VDR_expression_candidate_v2.csv"), VDR_summary = file.path(candidate_v2_root, 
        "source_data", "SourceData_STEP03B1_VDR_expression_summary_v2.csv"), volcano = file.path(candidate_v2_root, 
        "source_data", "SourceData_STEP03B1_interaction_volcano_candidate_v2.csv"), heatmap = file.path(candidate_v2_root, 
        "source_data", "SourceData_STEP03B1_top_interaction_heatmap_candidate_v2.csv"), heatmap_png = file.path(candidate_v2_root, 
        "figures", "Figure_STEP03B1_top_interaction_heatmap_candidate_v2_QC.png"))

if (!file.exists(v4_gate_path)) {
    stop("STEP03B1 v4 gate is missing.")
}

missing_inputs <- c(required_source_paths[!file.exists(required_source_paths)], upstream_result_paths[!file.exists(upstream_result_paths)])

if (length(missing_inputs) > 0L) {
    stop("Missing locked source or upstream result files:\n", paste(missing_inputs, collapse = "\n"))
}

v4_gate <- .vdr_read_csv(v4_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

v4_gate_value <- v4_gate$status[v4_gate$gate_component == "STEP03B1 v4 gate"]

if (length(v4_gate_value) != 1L || v4_gate_value != "PROCEED_TO_REVISED_FIGURE_RENDER_QC") {
    stop("STEP03B1 v4 gate does not authorize this targeted replot.")
}

pca_df <- .vdr_read_csv(required_source_paths["PCA"], stringsAsFactors = FALSE, check.names = FALSE)

vdr_df <- .vdr_read_csv(required_source_paths["VDR"], stringsAsFactors = FALSE, check.names = FALSE)

vdr_summary <- .vdr_read_csv(required_source_paths["VDR_summary"], stringsAsFactors = FALSE, check.names = FALSE)

interaction_plot <- .vdr_read_csv(required_source_paths["volcano"], stringsAsFactors = FALSE, check.names = FALSE)

heatmap_source <- .vdr_read_csv(required_source_paths["heatmap"], stringsAsFactors = FALSE, check.names = FALSE)

required_pca_columns <- c("gsm", "PC1", "PC2", "short_label", "genetic_context", "ligand_status_locked", "replicate_locked", 
    "analysis_group")

missing_pca_columns <- setdiff(required_pca_columns, names(pca_df))

if (length(missing_pca_columns) > 0L) {
    stop("PCA source table is missing required columns: ", paste(missing_pca_columns, collapse = ", "))
}

pca_core <- pca_df[, required_pca_columns, drop = FALSE]

if (nrow(pca_core) != 12L || anyDuplicated(pca_core$gsm)) {
    stop("PCA source must contain exactly 12 unique GSM sample rows.")
}

pca_core$siRNA_label <- ifelse(pca_core$genetic_context == "siRNA_control", "siNT", "siVDR")

pca_core$ligand_label <- ifelse(pca_core$ligand_status_locked == "vehicle", "Vehicle", "Calcitriol")

pca_core$group_label <- paste(pca_core$siRNA_label, pca_core$ligand_label, sep = " + ")

pca_centroids <- aggregate(cbind(PC1, PC2) ~ group_label + siRNA_label + ligand_label, data = pca_core, FUN = mean)

if (nrow(pca_centroids) != 4L || any(table(pca_core$group_label) != 3L)) {
    stop("PCA source does not contain four groups with three samples each.")
}

pca_centroids$short_group_label <- c(`siNT + Calcitriol` = "NT + Cal", `siNT + Vehicle` = "NT + Veh", `siVDR + Calcitriol` = "KD + Cal", 
    `siVDR + Vehicle` = "KD + Veh")[pca_centroids$group_label]

pca_centroids$label_x <- pca_centroids$PC1

pca_centroids$label_y <- pca_centroids$PC2

for (i in seq_len(nrow(pca_centroids))) {
    label <- pca_centroids$short_group_label[i]
    if (label == "NT + Veh") {
        pca_centroids$label_x[i] <- pca_centroids$PC1[i] - 2.2000000000000002
        pca_centroids$label_y[i] <- pca_centroids$PC2[i] + 0.69999999999999996
    }
    if (label == "NT + Cal") {
        pca_centroids$label_x[i] <- pca_centroids$PC1[i] + 1
        pca_centroids$label_y[i] <- pca_centroids$PC2[i]
    }
    if (label == "KD + Veh") {
        pca_centroids$label_x[i] <- pca_centroids$PC1[i] - 1.7
        pca_centroids$label_y[i] <- pca_centroids$PC2[i] + 1
    }
    if (label == "KD + Cal") {
        pca_centroids$label_x[i] <- pca_centroids$PC1[i] - 1.7
        pca_centroids$label_y[i] <- pca_centroids$PC2[i] - 1
    }
}

pca_centroids_for_merge <- pca_centroids

names(pca_centroids_for_merge)[names(pca_centroids_for_merge) == "PC1"] <- "PC1_centroid"

names(pca_centroids_for_merge)[names(pca_centroids_for_merge) == "PC2"] <- "PC2_centroid"

pca_source_v3 <- merge(pca_core, pca_centroids_for_merge, by = c("group_label", "siRNA_label", "ligand_label"), 
    all.x = TRUE, sort = FALSE)

pca_source_v3 <- pca_source_v3[match(pca_core$gsm, pca_source_v3$gsm), , drop = FALSE]

write_csv_utf8(pca_source_v3, file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_PCA_candidate_v3.csv"))

p_pca <- ggplot2::ggplot(pca_core, ggplot2::aes(x = PC1, y = PC2, color = siRNA_label, shape = ligand_label)) + 
    ggplot2::geom_point(size = 2.5, alpha = 0.90000000000000002) + ggplot2::geom_segment(data = pca_centroids, 
    ggplot2::aes(x = PC1, y = PC2, xend = label_x, yend = label_y), inherit.aes = FALSE, linewidth = 0.34999999999999998, 
    color = "grey35") + ggplot2::geom_point(data = pca_centroids, ggplot2::aes(x = PC1, y = PC2), inherit.aes = FALSE, 
    shape = 21, fill = "white", color = "black", size = 3.6000000000000001, stroke = 0.80000000000000004) + ggplot2::geom_text(data = pca_centroids, 
    ggplot2::aes(x = label_x, y = label_y, label = short_group_label), inherit.aes = FALSE, family = "Arial", size = 2.7000000000000002, 
    fontface = "plain") + ggplot2::scale_color_manual(values = c(siNT = "#0072B2", siVDR = "#D55E00")) + ggplot2::scale_shape_manual(values = c(Vehicle = 16, 
    Calcitriol = 17)) + ggplot2::labs(title = "GSE254251 sample structure", subtitle = "Top 1,000 variable expressed features", 
    x = "PC1 (68.3%)", y = "PC2 (13.5%)", color = "siRNA", shape = "Ligand") + theme_skill3(9) + ggplot2::theme(legend.position = "right")

save_qc_png(p_pca, file.path(candidate_v3_root, "figures", "Figure_STEP03B1_PCA_candidate_v3_QC.png"), width_mm = 110, 
    height_mm = 82)

required_vdr_columns <- c("gsm", "log2_RPKM_plus_1", "short_label", "genetic_context", "ligand_status_locked", 
    "analysis_group")

missing_vdr_columns <- setdiff(required_vdr_columns, names(vdr_df))

if (length(missing_vdr_columns) > 0L) {
    stop("VDR source table is missing required columns: ", paste(missing_vdr_columns, collapse = ", "))
}

vdr_core <- vdr_df[, required_vdr_columns, drop = FALSE]

if (nrow(vdr_core) != 12L || anyDuplicated(vdr_core$gsm)) {
    stop("VDR source must contain exactly 12 unique GSM sample rows.")
}

vdr_core$replicate_locked <- suppressWarnings(as.integer(sub("^.*_([1-3])$", "\\1", vdr_core$short_label, perl = TRUE)))

if (any(!is.finite(vdr_core$replicate_locked)) || any(!vdr_core$replicate_locked %in% 1:3)) {
    stop("Could not recover replicate 1/2/3 from every VDR short_label.")
}

vdr_core$siRNA_label <- factor(ifelse(vdr_core$genetic_context == "siRNA_control", "siNT", "siVDR"), levels = c("siNT", 
    "siVDR"))

vdr_core$ligand_label <- factor(ifelse(vdr_core$ligand_status_locked == "vehicle", "Vehicle", "Calcitriol"), levels = c("Vehicle", 
    "Calcitriol"))

vdr_core$plot_group <- interaction(vdr_core$siRNA_label, vdr_core$ligand_label, drop = TRUE, lex.order = TRUE)

vdr_group_counts <- table(vdr_core$siRNA_label, vdr_core$ligand_label)

if (!all(dim(vdr_group_counts) == c(2L, 2L)) || any(vdr_group_counts != 3L)) {
    stop("VDR source does not contain four siRNA-by-ligand groups with n=3.")
}

replicate_check <- vapply(split(vdr_core$replicate_locked, vdr_core$plot_group), function(x) {
    identical(sort(as.integer(x)), 1:3)
}, logical(1))

if (!all(replicate_check)) {
    stop("At least one VDR group does not contain replicates 1, 2 and 3.")
}

vdr_summary <- aggregate(log2_RPKM_plus_1 ~ siRNA_label + ligand_label, data = vdr_core, FUN = function(x) {
    c(mean = mean(x), sd = stats::sd(x), n = length(x))
})

vdr_summary <- data.frame(siRNA_label = vdr_summary$siRNA_label, ligand_label = vdr_summary$ligand_label, mean = vdr_summary$log2_RPKM_plus_1[, 
    "mean"], sd = vdr_summary$log2_RPKM_plus_1[, "sd"], n = vdr_summary$log2_RPKM_plus_1[, "n"], stringsAsFactors = FALSE)

vdr_summary$siRNA_label <- factor(vdr_summary$siRNA_label, levels = c("siNT", "siVDR"))

vdr_summary$ligand_label <- factor(vdr_summary$ligand_label, levels = c("Vehicle", "Calcitriol"))

vdr_core$x_base <- ifelse(vdr_core$ligand_label == "Vehicle", 1, 2)

vdr_core$x_group_offset <- ifelse(vdr_core$siRNA_label == "siNT", -0.14000000000000001, 0.14000000000000001)

replicate_offsets <- c(`1` = -0.044999999999999998, `2` = 0, `3` = 0.044999999999999998)

vdr_core$x_plot <- vdr_core$x_base + vdr_core$x_group_offset + unname(replicate_offsets[as.character(vdr_core$replicate_locked)])

if (length(vdr_core$x_plot) != nrow(vdr_core) || any(!is.finite(vdr_core$x_plot))) {
    stop("VDR x-position construction failed.")
}

vdr_summary$x_base <- ifelse(vdr_summary$ligand_label == "Vehicle", 1, 2)

vdr_summary$x_group_offset <- ifelse(vdr_summary$siRNA_label == "siNT", -0.14000000000000001, 0.14000000000000001)

vdr_summary$x_plot <- vdr_summary$x_base + vdr_summary$x_group_offset

write_csv_utf8(vdr_core, file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_VDR_expression_candidate_v3.csv"))

write_csv_utf8(vdr_summary, file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_VDR_expression_summary_v3.csv"))

p_vdr <- ggplot2::ggplot() + ggplot2::geom_errorbar(data = vdr_summary, ggplot2::aes(x = x_plot, ymin = mean - 
    sd, ymax = mean + sd, color = siRNA_label), width = 0.055, linewidth = 0.55000000000000004, show.legend = FALSE) + 
    ggplot2::geom_point(data = vdr_core, ggplot2::aes(x = x_plot, y = log2_RPKM_plus_1, color = siRNA_label), size = 2.1499999999999999, 
        alpha = 0.90000000000000002) + ggplot2::geom_crossbar(data = vdr_summary, ggplot2::aes(x = x_plot, y = mean, 
    ymin = mean, ymax = mean, color = siRNA_label), width = 0.13, linewidth = 0.65000000000000002, show.legend = FALSE) + 
    ggplot2::scale_color_manual(values = c(siNT = "#0072B2", siVDR = "#D55E00")) + ggplot2::scale_x_continuous(breaks = c(1, 
    2), labels = c("Vehicle", "Calcitriol"), limits = c(0.65000000000000002, 2.3500000000000001)) + ggplot2::scale_y_continuous(limits = c(0, 
    NA), expand = ggplot2::expansion(mult = c(0, 0.080000000000000002))) + ggplot2::labs(title = "VDR knockdown reduced VDR expression", 
    subtitle = "Vehicle and calcitriol conditions", x = NULL, y = "VDR log2(RPKM + 1)", color = "siRNA") + theme_skill3(9) + 
    ggplot2::theme(legend.position = "top", legend.justification = "left")

save_qc_png(p_vdr, file.path(candidate_v3_root, "figures", "Figure_STEP03B1_VDR_expression_candidate_v3_QC.png"), 
    width_mm = 105, height_mm = 75)

interaction_plot$plot_class <- factor(interaction_plot$plot_class, levels = c("Not candidate", "Attenuated induction", 
    "Attenuated repression", "Other interaction"))

select_top_labels <- function(x, class_name, n) {
    y <- x[x$plot_class == class_name & x$is_interaction_candidate, , drop = FALSE]
    y <- y[order(y$interaction_FDR, -abs(y$interaction_logFC)), , drop = FALSE]
    head(y, n)
}

label_df <- rbind(select_top_labels(interaction_plot, "Attenuated induction", 5L), select_top_labels(interaction_plot, 
    "Attenuated repression", 5L), select_top_labels(interaction_plot, "Other interaction", 2L))

write_csv_utf8(interaction_plot, file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_interaction_volcano_candidate_v3.csv"))

write_csv_utf8(label_df, file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_interaction_volcano_labels_v3.csv"))

p_volcano <- ggplot2::ggplot(interaction_plot, ggplot2::aes(x = interaction_logFC, y = minus_log10_FDR, color = plot_class)) + 
    ggplot2::geom_point(size = 1.1200000000000001, alpha = 0.57999999999999996) + ggplot2::geom_vline(xintercept = c(-program_abs_logfc, 
    program_abs_logfc), linetype = 2, linewidth = 0.34999999999999998, color = "grey35") + ggplot2::geom_hline(yintercept = -log10(program_fdr), 
    linetype = 2, linewidth = 0.34999999999999998, color = "grey35") + ggrepel::geom_text_repel(data = label_df, 
    ggplot2::aes(label = feature_id), family = "Arial", size = 2.3500000000000001, min.segment.length = 0, max.overlaps = Inf, 
    box.padding = 0.29999999999999999, point.padding = 0.12, seed = 20260803, force = 1.2, show.legend = FALSE) + 
    ggplot2::scale_color_manual(values = c(`Not candidate` = "#BDBDBD", `Attenuated induction` = "#7B3294", `Attenuated repression` = "#008837", 
        `Other interaction` = "#E6AB02")) + ggplot2::labs(title = "VDR knockdown attenuated calcitriol-responsive transcription", 
    subtitle = paste0("Formal siRNA x ligand interaction\n", "FDR < ", program_fdr, "; |interaction log2FC| >= ", 
        program_abs_logfc), x = "Interaction log2 fold change", y = expression(-log[10]("FDR")), color = NULL) + 
    theme_skill3(9) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical") + ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2.3999999999999999, 
    alpha = 1), nrow = 2, byrow = TRUE))

save_qc_png(p_volcano, file.path(candidate_v3_root, "figures", "Figure_STEP03B1_interaction_volcano_candidate_v3_QC.png"), 
    width_mm = 132, height_mm = 94)

copy_validate(required_source_paths["heatmap"], file.path(candidate_v3_root, "source_data", "SourceData_STEP03B1_top_interaction_heatmap_ACCEPTED_v2.csv"))

copy_validate(required_source_paths["heatmap_png"], file.path(candidate_v3_root, "figures", "Figure_STEP03B1_top_interaction_heatmap_ACCEPTED_v2_QC.png"))

for (source_path in upstream_result_paths) {
    copy_validate(source_path, file.path(candidate_v3_root, "source_data", "upstream_refs", basename(source_path)))
}

figure_source_map <- data.frame(figure_file = c("figures/Figure_STEP03B1_PCA_candidate_v3_QC.png", "figures/Figure_STEP03B1_VDR_expression_candidate_v3_QC.png", 
    "figures/Figure_STEP03B1_interaction_volcano_candidate_v3_QC.png", "figures/Figure_STEP03B1_top_interaction_heatmap_ACCEPTED_v2_QC.png"), 
    figure_role = c("Supplementary candidate", "Strong mechanistic-QC candidate", "Strong main-figure candidate", 
        "Strong main-figure candidate"), direct_source_data = c("source_data/SourceData_STEP03B1_PCA_candidate_v3.csv", 
        paste("source_data/SourceData_STEP03B1_VDR_expression_candidate_v3.csv;", "source_data/SourceData_STEP03B1_VDR_expression_summary_v3.csv"), 
        paste("source_data/SourceData_STEP03B1_interaction_volcano_candidate_v3.csv;", "source_data/SourceData_STEP03B1_interaction_volcano_labels_v3.csv"), 
        "source_data/SourceData_STEP03B1_top_interaction_heatmap_ACCEPTED_v2.csv"), upstream_result_reference = c("source_data/upstream_refs/STEP03B1_sample_manifest_FINAL.csv", 
        paste("source_data/upstream_refs/STEP03B1_sample_manifest_FINAL.csv;", "source_data/upstream_refs/STEP03B1_limma_all_contrasts_long.csv"), 
        paste("source_data/upstream_refs/STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv;", "source_data/upstream_refs/STEP03B1_v4_interaction_program_summary.csv;", 
            "source_data/upstream_refs/STEP03B1_limma_all_contrasts_long.csv"), paste("source_data/upstream_refs/STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv;", 
            "source_data/upstream_refs/STEP03B1_v4_RNA_program_definition_card.csv")), authoritative_generator_code = "code/15_step03B1_v7_Skill3_FAST_REPLOT_schema_safe.R", 
    authoritative_plotting_language = "R", statistical_unit = "biological well", analysis_recomputed_during_plotting = "NO", 
    current_figure_class = "B_submission_candidate", current_release_status = "PENDING_FINAL_RENDER_QC", stringsAsFactors = FALSE)

write_csv_utf8(figure_source_map, file.path(candidate_v3_root, "manifests", "STEP03B1_v3_FIGURE_SOURCE_MAP.csv"))

source_data_readme <- c("# STEP03B1 v3 figure source-data guide", "", "This folder is the self-contained source package for later figure QC.", 
    "No differential-expression analysis needs to be rerun to revise these plots.", "", "## Direct figure sources", 
    "", "- PCA: `SourceData_STEP03B1_PCA_candidate_v3.csv`", "- VDR raw points: `SourceData_STEP03B1_VDR_expression_candidate_v3.csv`", 
    "- VDR mean/SD: `SourceData_STEP03B1_VDR_expression_summary_v3.csv`", "- Interaction volcano: `SourceData_STEP03B1_interaction_volcano_candidate_v3.csv`", 
    "- Volcano labels: `SourceData_STEP03B1_interaction_volcano_labels_v3.csv`", "- Accepted interaction heatmap: `SourceData_STEP03B1_top_interaction_heatmap_ACCEPTED_v2.csv`", 
    "", "## Upstream references", "", "The `upstream_refs/` subfolder contains the exact locked result tables", 
    "needed to verify figure meaning, sample identity, contrast statistics and", "program classification without searching elsewhere in the project.", 
    "", "## Routing", "", "See `../manifests/STEP03B1_v3_FIGURE_SOURCE_MAP.csv` for the exact", "figure-to-source-to-code relationship.", 
    "", "## Status", "", "- Figure class: B submission candidate", "- Final PDF/TIFF/SVG release: not yet generated", 
    "- Program freeze: not yet authorized")

write_text(source_data_readme, file.path(candidate_v3_root, "source_data", "README_SOURCE_DATA.md"))

language_routing_card <- c("# STEP03B1 v3 language routing card", "", "- Governing skill: research_scientific_figure_suite v1.7-beta", 
    "- Mode: FAST_REPLOT", "- Authoritative analysis language: R", "- Authoritative plotting language: R", "- Cross-language migration: NO", 
    "- Inferential statistics recomputed during plotting: NO", "- Direct source tables: LOCKED", "- Actual-render QC: PENDING", 
    "- Final release: NOT AUTHORIZED")

write_text(language_routing_card, file.path(candidate_v3_root, "LANGUAGE_ROUTING_CARD.md"))

legend_lines <- c("# STEP03B1 v3 candidate figure legends", "", "**PCA candidate.** Principal-component analysis of the 12 GSE254251 samples using the top 1,000 variable expressed features. Individual points are biological wells; open circles mark condition centroids.", 
    "", "**VDR-expression candidate.** VDR log2(RPKM + 1) values under vehicle and calcitriol conditions. Small points represent biological wells, horizontal bars show means and error bars show standard deviations.", 
    "", "**Interaction-volcano candidate.** Formal siRNA-by-ligand interaction effects. Candidate interactions satisfy FDR < 0.05 and absolute interaction log2 fold change >= 0.5. Labels were selected in a balanced manner across interaction classes.", 
    "", "**Interaction heatmap.** Row-standardized expression of the top 30 genes ranked by interaction FDR. The accepted v2 heatmap was copied without regeneration.")

write_text(legend_lines, file.path(candidate_v3_root, "legends", "STEP03B1_v3_candidate_legends.md"))

qc_lines <- c("# STEP03B1 Skill 3 targeted FAST_REPLOT round 2", "", "- PCA: REVISED_PENDING_ACTUAL_RENDER_QC", 
    "- VDR expression: REVISED_PENDING_ACTUAL_RENDER_QC", "- Interaction volcano: REVISED_PENDING_ACTUAL_RENDER_QC", 
    "- Interaction heatmap: PASS_AS_B_CLASS_CANDIDATE", "- Limma model: UNCHANGED", "- Samples: UNCHANGED", "- Thresholds: UNCHANGED", 
    "- RNA programs: NOT_FROZEN", "- ATAC analysis: NOT_AUTHORIZED", "- Submission-ready: NO")

write_text(qc_lines, file.path(candidate_v3_root, "qc", "STEP03B1_v3_RENDER_QC_STATUS.md"))

gate_table <- data.frame(gate_component = c("Locked v2 source data", "PCA targeted replot", "VDR-expression targeted replot", 
    "Interaction-volcano targeted replot", "Accepted heatmap copy and hash validation", "Limma/model changes", 
    "Program freeze", "STEP03B1 v5 gate"), status = c("PASS", "PASS_CREATED", "PASS_CREATED", "PASS_CREATED", "PASS", 
    "NONE", "NOT_FROZEN", "PROCEED_TO_FINAL_RENDER_QC"), required_action = c("Retain source hashes.", "Inspect actual v3 PNG.", 
    "Inspect raw-point visibility and clipping.", "Inspect title, subtitle, labels and legend.", "Retain accepted v2 render unchanged.", 
    "Do not refit or alter statistics.", "Freeze only after final render and content QC.", "Upload the complete v3 candidate package."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B1_v7_GATE_TABLE.csv"))

write_text(c("# STEP03B1 v7 schema-safe FAST_REPLOT and source-routing report", "", paste0("- Run ID: `", run_id, 
    "`"), "- Gate: `PROCEED_TO_FINAL_RENDER_QC`", "- Figure-source map: CREATED.", "- PCA: revised.", "- VDR expression: revised.", 
    "- Interaction volcano: revised.", "- Interaction heatmap: accepted v2 copied unchanged.", "- PCA source schema normalized before centroid recomputation.", 
    "- VDR replicate numbers recovered from locked short labels.", "- VDR descriptive mean/SD rebuilt from locked raw points.", 
    "- Limma refit: NO.", "- Threshold change: NO.", "- Program freeze: NO.", "- ATAC authorization: NO."), file.path(dirs$results, 
    "STEP03B1_v7_TARGETED_FAST_REPLOT_REPORT.md"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "prepare_experimental_source_tables.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
    file.copy(source_file, file.path(candidate_v3_root, "code", "prepare_experimental_source_tables.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
}

parameters <- list(schema_version = "1.0", run_id = run_id, source_candidate_version = "v2", source_routing_manifest = "STEP03B1_v3_FIGURE_SOURCE_MAP.csv", 
    output_candidate_version = "v3", limma_refit = FALSE, sample_change = FALSE, threshold_change = FALSE, program_frozen = FALSE, 
    heatmap_regenerated = FALSE, heatmap_status = "accepted_v2", gate = "PROCEED_TO_FINAL_RENDER_QC")

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B1_v7_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B1_v7.txt"))

candidate_files <- list.files(candidate_v3_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

candidate_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", candidate_files), nchar(gsub("\\\\", 
    "/", candidate_v3_root)) + 2L), size_bytes = as.numeric(file.info(candidate_files)$size), md5 = vapply(candidate_files, 
    md5_file, character(1)), sha256 = vapply(candidate_files, sha256_file, character(1)), figure_class = "B_submission_candidate", 
    release_status = "PENDING_FINAL_RENDER_QC", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest, file.path(candidate_v3_root, "manifests", "STEP03B1_v3_candidate_package_manifest.csv"))

archive_bv3 <- file.path(dirs$archive, "07_figures_tables", "Bv3")

dir.create(archive_bv3, recursive = TRUE, showWarnings = FALSE)

candidate_files <- list.files(candidate_v3_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_rows <- vector("list", length(candidate_files))

for (i in seq_along(candidate_files)) {
    source_path <- candidate_files[i]
    relative_path <- substring(gsub("\\\\", "/", source_path), nchar(gsub("\\\\", "/", candidate_v3_root)) + 2L)
    destination_path <- file.path(archive_bv3, relative_path)
    if (nchar(destination_path) >= 250L) {
        stop("Archive destination path remains too long:\n", destination_path)
    }
    copy_validate(source_path, destination_path)
    archive_rows[[i]] <- data.frame(relative_path = relative_path, source_md5 = md5_file(source_path), archived_md5 = md5_file(destination_path), 
        source_sha256 = sha256_file(source_path), archived_sha256 = sha256_file(destination_path), validation_pass = md5_file(source_path) == 
            md5_file(destination_path) && sha256_file(source_path) == sha256_file(destination_path), stringsAsFactors = FALSE)
}

archive_copy_manifest <- do.call(rbind, archive_rows)

write_csv_utf8(archive_copy_manifest, file.path(dirs$archive, "00_archive_control", "STEP03B1_v7_Bv3_archive_copy_manifest.csv"))

if (!all(archive_copy_manifest$validation_pass)) {
    stop("Bv3 archive copy/hash validation failed.")
}

write_text(c("# STEP03B1 v7 archive status", "", paste0("- Run ID: ", run_id), "- Bv3 candidate archive: PASS", 
    paste0("- Candidate files validated: ", nrow(archive_copy_manifest), "/", nrow(archive_copy_manifest)), "- Authoritative archive path: `07_figures_tables/Bv3/`", 
    "- Limma/model: UNCHANGED", "- RNA program: NOT_FROZEN", "- Gate: PROCEED_TO_FINAL_RENDER_QC", "- ATAC analysis: NOT_AUTHORIZED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B1_v7", paste0("- Run ID: ", run_id), "- RNA analysis: COMPLETED_PENDING_FINAL_RENDER_QC", 
    "- RNA program: NOT_FROZEN", "- Candidate figure package: v3", "- ATAC analysis: NOT_STARTED", "- Full module: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate: PROCEED_TO_FINAL_RENDER_QC")

cat("\n============================================================\n")

cat("STEP03B1 v7 Skill 3 schema-safe FAST_REPLOT completed.\n")

cat("Gate: PROCEED_TO_FINAL_RENDER_QC\n")

cat("PCA v3: created\n")

cat("VDR expression v3: created\n")

cat("Interaction volcano v3: created\n")

cat("Interaction heatmap: accepted v2 copied unchanged\n")

cat("Limma was not refitted.\n")

cat("RNA programs remain NOT FROZEN.\n")

cat("ATAC analysis remains unauthorized.\n")

cat("Figure-source routing manifest: CREATED\n")

cat("Self-contained source-data package: CREATED\n")

cat("Bv3 archive copy/hash validation: PASS\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B1_v7_GATE_TABLE.csv\n")

cat("  06_results/STEP03B1_v7_TARGETED_FAST_REPLOT_REPORT.md\n")

cat("  07_figures_tables/B_submission_candidate/STEP03B1_RNA_candidates_v3/\n")

cat("============================================================\n")

