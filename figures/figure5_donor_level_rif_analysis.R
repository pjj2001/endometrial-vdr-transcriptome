#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="figure5_donor_level_rif_analysis.R" & file.exists(.paths)]
if (length(.paths)) {.script <- tail(.paths,1)} else {
  .arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if (!length(.arg)) stop('Run this file with Rscript or source its complete path.')
  .script <- sub('^--file=','',.arg[[1]])
}
repo_root <- normalizePath(file.path(dirname(.script),'..'),winslash='/',mustWork=TRUE)
source(file.path(repo_root,'figures','figure_common.R'),local=TRUE)
figure_data_root <- normalizePath(Sys.getenv('VDR_FIGURE_DATA',unset=file.path(repo_root,'data','figures')),winslash='/',mustWork=TRUE)
figure_index <- read.csv(file.path(figure_data_root,'inputs.csv'),stringsAsFactors=FALSE,check.names=FALSE)
figure_font <- 'Arial'
figure_preflight('Figure5',c("data.table","ggplot2","patchwork","gtable","ragg","systemfonts","scales","R.utils"))
if ('--check' %in% commandArgs(TRUE)) {message('Figure5: inputs, packages and font available');quit(status=0)}
source_paths <- setNames(vapply(figure_index$key[figure_index$figure=="Figure5"],function(k)figure_input("Figure5",k),character(1)),figure_index$key[figure_index$figure=="Figure5"])

cell_levels <- c("Ciliated epithelial cells", "Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells")

cell_labels <- c("Ciliated epi.", "Endothelial", "Myeloid", "NK/T", "Stromal", "Unciliated epi.")

names(cell_labels) <- cell_levels

cell_colors <- c(`Ciliated epithelial cells` = "#E69F00", `Endothelial cells` = "#009E73", `Myeloid cells` = "#00A6C7", `NK/T cells` = "#007DB3", 
    `Stromal cells` = "#7A5ABB", `Unciliated epithelial cells` = "#CC79A7")

group_colors <- c(CONTROL = "#5B7896", RIF = "#C65D68")

endpoint_levels <- c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84", "VDR_expression")

endpoint_labels <- c("RNA program\n(n=436)", "Occupancy-supported\n(n=284)", "Shared occupancy\n(n=46)", "Promoter-associated\n(n=84)", 
    "VDR expression")

names(endpoint_labels) <- endpoint_levels

endpoint_short_labels <- c("RNA\n436", "Occupancy\n284", "Shared\n46", "Promoter\n84", "VDR")

names(endpoint_short_labels) <- endpoint_levels

plot_data <- capture_stage("read_and_prepare_frozen_sources", {
    g250_program <- data.table::fread(source_paths[["g250_program"]])
    g250_vdr <- data.table::fread(source_paths[["g250_vdr"]])
    g250_primary_de <- data.table::fread(source_paths[["g250_primary_de"]])
    g250_sensitivity_de <- data.table::fread(source_paths[["g250_sensitivity_de"]])
    g183_annotation <- data.table::fread(source_paths[["g183_annotation"]])
    g183_umap <- data.table::fread(source_paths[["g183_umap"]], select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters", 
        "analysis_group", "sample_accession"))
    g183_models <- data.table::fread(source_paths[["g183_models"]])
    g183_group_summary <- data.table::fread(source_paths[["g183_group_summary"]])
    g183_model_data <- data.table::fread(source_paths[["g183_model_data"]])
    g183_ciliated_summary <- data.table::fread(source_paths[["g183_ciliated_summary"]])
    g183_ciliated_data <- data.table::fread(source_paths[["g183_ciliated_data"]])
    a_program <- g250_program[, .(broad_cell_type, endpoint = program_name, estimate = RIF_minus_fertile, ci_low = Welch_ci_low, 
        ci_high = Welch_ci_high, FDR = FDR_BH, exact_P = exact_permutation_p, fertile_n, RIF_n)]
    a_vdr <- g250_vdr[, .(broad_cell_type, endpoint = "VDR_expression", estimate = RIF_minus_fertile, ci_low = Welch_ci_low, 
        ci_high = Welch_ci_high, FDR = FDR_BH, exact_P = exact_permutation_p, fertile_n, RIF_n)]
    panel_A <- data.table::rbindlist(list(a_program, a_vdr), use.names = TRUE)
    panel_A[, `:=`(endpoint_display, factor(endpoint_labels[endpoint], levels = endpoint_labels))]
    panel_A[, `:=`(endpoint_display_short, factor(endpoint_short_labels[endpoint], levels = endpoint_short_labels))]
    panel_A[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(cell_levels)))]
    primary_counts <- g250_primary_de[, .(broad_cell_type, count = whole_genome_up_n + whole_genome_down_n)]
    primary_counts[, `:=`(scope, "Primary LH+7\nRIF vs fertile")]
    sensitivity_counts <- g250_sensitivity_de[, .N, by = broad_cell_type]
    data.table::setnames(sensitivity_counts, "N", "count")
    sensitivity_counts[, `:=`(scope, "Fertile-reference\nsensitivity")]
    all_cells <- data.table::data.table(broad_cell_type = cell_levels)
    primary_counts <- merge(all_cells, primary_counts, by = "broad_cell_type", all.x = TRUE)
    primary_counts[is.na(count), `:=`(count, 0)]
    primary_counts[is.na(scope), `:=`(scope, "Primary LH+7\nRIF vs fertile")]
    sensitivity_counts <- merge(all_cells, sensitivity_counts, by = "broad_cell_type", all.x = TRUE)
    sensitivity_counts[is.na(count), `:=`(count, 0)]
    sensitivity_counts[is.na(scope), `:=`(scope, "Fertile-reference\nsensitivity")]
    panel_B <- data.table::rbindlist(list(primary_counts, sensitivity_counts), use.names = TRUE)
    panel_B[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(cell_levels)))]
    panel_B[, `:=`(scope, factor(scope, levels = c("Primary LH+7\nRIF vs fertile", "Fertile-reference\nsensitivity")))]
    g183_umap[, `:=`(cluster, as.character(seurat_clusters))]
    g183_annotation[, `:=`(cluster, as.character(cluster))]
    panel_C_all <- merge(g183_umap, g183_annotation[, .(cluster, broad_cell_type)], by = "cluster", all.x = TRUE, sort = FALSE)
    if (anyNA(panel_C_all$broad_cell_type)) 
        stop("GSE183837 UMAP-to-annotation join is incomplete.")
    data.table::setorder(panel_C_all, broad_cell_type, cell_barcode)
    panel_C <- panel_C_all[, .SD[seq_len(min(.N, 6000L))], by = broad_cell_type]
    centroids <- panel_C[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)), by = broad_cell_type]
    centroids[, `:=`(label_x = UMAPH_1 + data.table::fcase(broad_cell_type == "Unciliated epithelial cells", -2.2000000000000002, 
        broad_cell_type == "Endothelial cells", 2.2000000000000002, default = 0), label_y = UMAPH_2)]
    panel_D <- g183_models[, .(broad_cell_type, endpoint, estimate = mean_difference_RIF_minus_Control, ci_low = Welch_95CI_lower, 
        ci_high = Welch_95CI_upper, exact_P = exact_permutation_P, global_BH = BH_global_25, statistical_status, constant_endpoint)]
    panel_D[!is.finite(ci_low), `:=`(ci_low, estimate)]
    panel_D[!is.finite(ci_high), `:=`(ci_high, estimate)]
    panel_D[, `:=`(endpoint_display, factor(endpoint_labels[endpoint], levels = endpoint_labels))]
    panel_D[, `:=`(endpoint_display_short, factor(endpoint_short_labels[endpoint], levels = endpoint_short_labels))]
    panel_D[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(cell_levels)))]
    panel_D[, `:=`(nominal_marker, ifelse(exact_P < 0.050000000000000003, "*", ""))]
    panel_E <- data.table::copy(g183_ciliated_data)
    panel_E[, `:=`(endpoint_display, factor(endpoint_labels[endpoint], levels = endpoint_labels))]
    panel_E[, `:=`(analysis_group, factor(analysis_group, levels = c("CONTROL", "RIF")))]
    panel_E[, `:=`(threshold_status, factor(ifelse(minimum_cell_threshold_pass, "Threshold pass", "Below threshold"), levels = c("Threshold pass", 
        "Below threshold")))]
    list(g250_program = g250_program, g250_vdr = g250_vdr, g250_primary_de = g250_primary_de, g250_sensitivity_de = g250_sensitivity_de, 
        g183_models = g183_models, g183_group_summary = g183_group_summary, g183_model_data = g183_model_data, g183_ciliated_summary = g183_ciliated_summary, 
        A = panel_A, B = panel_B, C_all = panel_C_all, C = panel_C, C_centroids = centroids, D = panel_D, E = panel_E)
})

scientific_checks <- data.frame(check_id = sprintf("S%02d", 1:15), check = c("GSE250130 LH+7 program models", "GSE250130 LH+7 program FDR-significant models", 
    "GSE250130 LH+7 VDR models", "GSE250130 LH+7 VDR FDR-significant models", "GSE250130 primary whole-genome DE rows", "GSE250130 sensitivity whole-genome DE rows", 
    "GSE183837 frozen UMAP cells", "GSE183837 locked broad cell types", "GSE183837 donors", "GSE183837 primary donor models", 
    "GSE183837 nominal exact-P signals", "GSE183837 global-BH signals", "GSE183837 ciliated endpoints", "GSE183837 ciliated donors per endpoint", 
    "Ciliated primary-threshold-pass RIF donors"), observed = c(nrow(plot_data$g250_program), sum(plot_data$g250_program$FDR_BH < 
    0.050000000000000003), nrow(plot_data$g250_vdr), sum(plot_data$g250_vdr$FDR_BH < 0.050000000000000003), sum(plot_data$g250_primary_de$whole_genome_up_n + 
    plot_data$g250_primary_de$whole_genome_down_n), nrow(plot_data$g250_sensitivity_de), nrow(plot_data$C_all), data.table::uniqueN(plot_data$C_all$broad_cell_type), 
    data.table::uniqueN(plot_data$C_all$sample_accession), nrow(plot_data$g183_models), sum(plot_data$g183_models$exact_permutation_P < 
        0.050000000000000003), sum(plot_data$g183_models$BH_global_25 < 0.050000000000000003), data.table::uniqueN(plot_data$E$endpoint), 
    min(plot_data$E[, data.table::uniqueN(donor_id), by = endpoint]$V1), data.table::uniqueN(plot_data$E[analysis_group == 
        "RIF" & minimum_cell_threshold_pass, donor_id])), expected = c(24, 0, 6, 0, 0, 28, 47161, 6, 9, 25, 1, 0, 5, 9, 2), 
    stringsAsFactors = FALSE)

scientific_checks$pass <- scientific_checks$observed == scientific_checks$expected

if (!all(scientific_checks$pass)) stop("Figure 5 scientific QC failed.")

color_text <- "#161616"

color_muted <- "#555555"

color_grid <- "#D9D9D9"

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.0999999999999996) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major = ggplot2::element_line(linewidth = 0.23000000000000001, colour = color_grid), axis.title = ggplot2::element_text(size = 6.7000000000000002, 
        colour = color_text), axis.text = ggplot2::element_text(size = 5.5999999999999996, colour = color_text), strip.background = ggplot2::element_rect(fill = "#F2F2F2", 
        colour = NA), strip.text = ggplot2::element_text(size = 5.5, face = "bold", colour = color_text, lineheight = 0.90000000000000002), 
    legend.title = ggplot2::element_text(size = 5.7000000000000002, face = "bold"), legend.text = ggplot2::element_text(size = 5.4000000000000004), 
    plot.margin = ggplot2::margin(1.2, 1.5, 1.2, 1.5, unit = "mm"))

panel_A <- capture_stage("panel_A_construction", {
    ggplot2::ggplot(plot_data$A, ggplot2::aes(x = estimate, y = broad_cell_type, colour = broad_cell_type)) + ggplot2::geom_vline(xintercept = 0, 
        linetype = 2, linewidth = 0.28000000000000003, colour = "#777777") + ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, 
        xmax = ci_high), orientation = "y", width = 0.17999999999999999, linewidth = 0.41999999999999998) + ggplot2::geom_point(size = 1.45) + 
        ggplot2::facet_grid(~endpoint_display_short, scales = "free_x") + ggplot2::scale_colour_manual(values = cell_colors, 
        drop = FALSE) + ggplot2::scale_y_discrete(labels = cell_labels) + ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3), 
        labels = scales::label_number(accuracy = 0.01), guide = ggplot2::guide_axis(check.overlap = TRUE)) + ggplot2::labs(x = "Mean difference (RIF - fertile; independent facet scales)", 
        y = NULL) + base_theme + ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_text(size = 4.7999999999999998), axis.text.y = ggplot2::element_text(size = 5.7000000000000002), 
        axis.title.x = ggplot2::element_text(size = 6.2999999999999998), panel.spacing.x = grid::unit(1.3, "mm"))
})

panel_B <- capture_stage("panel_B_construction", {
    ggplot2::ggplot(plot_data$B, ggplot2::aes(x = count, y = broad_cell_type, fill = broad_cell_type)) + ggplot2::geom_blank(data = data.frame(count = c(1, 
        16), broad_cell_type = factor(rep("Ciliated epithelial cells", 2), levels = rev(cell_levels)), scope = factor(c("Primary LH+7\nRIF vs fertile", 
        "Fertile-reference\nsensitivity"), levels = c("Primary LH+7\nRIF vs fertile", "Fertile-reference\nsensitivity")))) + 
        ggplot2::geom_col(width = 0.62) + ggplot2::geom_point(size = 0.90000000000000002, colour = color_text) + ggplot2::geom_text(ggplot2::aes(label = count), 
        hjust = -0.25, size = 2, family = figure_font, colour = color_text) + ggplot2::facet_grid(~scope, scales = "free_x") + 
        ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE) + ggplot2::scale_y_discrete(labels = cell_labels) + 
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.16)), breaks = scales::breaks_pretty(n = 4)) + 
        ggplot2::labs(x = "FDR-significant result rows (independent facet axes)", y = NULL) + base_theme + ggplot2::theme(legend.position = "none", 
        panel.grid.major.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 5.4000000000000004), axis.text.x = ggplot2::element_text(size = 5.0999999999999996), 
        panel.spacing.x = grid::unit(3, "mm"))
})

panel_C <- capture_stage("panel_C_construction", {
    ggplot2::ggplot(plot_data$C, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, colour = broad_cell_type)) + ggplot2::geom_point(size = 0.080000000000000002, 
        alpha = 0.47999999999999998, stroke = 0) + ggplot2::geom_label(data = plot_data$C_centroids, ggplot2::aes(x = label_x, 
        y = label_y, label = broad_cell_type), inherit.aes = FALSE, family = figure_font, fontface = "bold", size = 1.75, 
        linewidth = 0, label.padding = grid::unit(0.55000000000000004, "mm"), fill = scales::alpha("white", 0.81999999999999995), 
        colour = color_text) + ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) + ggplot2::coord_fixed(clip = "off") + 
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2") + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
        axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 5.9000000000000004), legend.position = "none")
})

panel_D <- capture_stage("panel_D_construction", {
    ggplot2::ggplot(plot_data$D, ggplot2::aes(x = estimate, y = broad_cell_type, colour = broad_cell_type)) + ggplot2::geom_vline(xintercept = 0, 
        linetype = 2, linewidth = 0.28000000000000003, colour = "#777777") + ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, 
        xmax = ci_high), orientation = "y", width = 0.17999999999999999, linewidth = 0.41999999999999998) + ggplot2::geom_point(size = 1.3500000000000001) + 
        ggplot2::geom_text(ggplot2::aes(label = nominal_marker), nudge_y = 0.23999999999999999, family = figure_font, fontface = "bold", 
            size = 2.7000000000000002, colour = color_text) + ggplot2::facet_grid(~endpoint_display_short, scales = "free_x") + 
        ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) + ggplot2::scale_y_discrete(labels = cell_labels) + 
        ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(accuracy = 0.01), 
            guide = ggplot2::guide_axis(check.overlap = TRUE)) + ggplot2::labs(x = "Mean difference (RIF - control; independent facet scales)", 
        y = NULL) + base_theme + ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_text(size = 4.2999999999999998), axis.text.y = ggplot2::element_text(size = 5.2000000000000002), 
        axis.title.x = ggplot2::element_text(size = 5.9000000000000004), panel.spacing.x = grid::unit(0.80000000000000004, 
            "mm"))
})

panel_E <- capture_stage("panel_E_construction", {
    ggplot2::ggplot(plot_data$E, ggplot2::aes(x = analysis_group, y = value, colour = analysis_group)) + ggplot2::geom_point(ggplot2::aes(shape = threshold_status), 
        position = ggplot2::position_jitter(width = 0.074999999999999997, height = 0, seed = 2505), size = 1.55, stroke = 0.55000000000000004) + 
        ggplot2::stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.41999999999999998, 
            linewidth = 0.47999999999999998) + ggplot2::facet_grid(~endpoint_display, scales = "free_y") + ggplot2::scale_colour_manual(values = group_colors) + 
        ggplot2::scale_shape_manual(values = c(`Threshold pass` = 16, `Below threshold` = 1)) + ggplot2::labs(x = NULL, y = "Donor-level value", 
        colour = "Group", shape = "Ciliated profile") + base_theme + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_text(size = 5.5), axis.text.y = ggplot2::element_text(size = 5), legend.position = "bottom", 
        legend.box = "horizontal", legend.margin = ggplot2::margin(t = 0.5, b = 0), panel.spacing.x = grid::unit(1.3999999999999999, 
            "mm"))
})

build_body_gtable <- function(plot_object) {
    ggplot2::ggplotGrob(plot_object + ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) + ggplot2::theme(plot.tag = ggplot2::element_blank(), 
        plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank()))
}

header_title_row_height_mm <- 5.0999999999999996

header_subtitle_row_height_mm <- 3.2000000000000002

header_tag_to_title_gap_mm <- 1.8

header_tag_size_pt <- 8.1999999999999993

header_title_size_pt <- 8

header_subtitle_size_pt <- 6.4000000000000004

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle, header_shift_mm = 0) {
    panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
    if (!nrow(panel_cells)) 
        stop("Could not locate panel cells for ", panel_tag)
    panel_left <- min(panel_cells$l)
    panel_right <- max(panel_cells$r)
    header_gt <- gtable::gtable(widths = body_gtable$widths, heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), 
        "mm"))
    header_shift <- grid::unit(header_shift_mm, "mm")
    title_grob <- grid::textGrob(panel_title, x = grid::unit(0, "npc") + header_shift, y = grid::unit(0.52000000000000002, 
        "npc"), hjust = 0, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_title_size_pt, 
        col = color_text))
    tag_grob <- grid::textGrob(panel_tag, x = grid::unit(0, "npc") + header_shift - grid::unit(header_tag_to_title_gap_mm, 
        "mm"), y = grid::unit(0.52000000000000002, "npc"), hjust = 1, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", 
        fontsize = header_tag_size_pt, col = color_text))
    subtitle_grob <- grid::textGrob(panel_subtitle, x = grid::unit(0, "npc") + header_shift, y = grid::unit(0.56000000000000005, 
        "npc"), hjust = 0, gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, col = color_muted))
    header_gt <- gtable::gtable_add_grob(header_gt, grid::grobTree(tag_grob, title_grob), t = 1, b = 1, l = panel_left, r = panel_right, 
        clip = "off", name = paste0("F5-header-title-", panel_tag))
    header_gt <- gtable::gtable_add_grob(header_gt, subtitle_grob, t = 2, b = 2, l = panel_left, r = panel_right, clip = "off", 
        name = paste0("F5-header-subtitle-", panel_tag))
    list(gtable = header_gt, geometry = data.frame(panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE", panel_left_column = panel_left, 
        panel_right_column = panel_right, tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt, subtitle_size_pt = header_subtitle_size_pt, 
        tag_title_gap_mm = header_tag_to_title_gap_mm, header_shift_mm = header_shift_mm, stringsAsFactors = FALSE))
}

figure_width_mm <- 180

figure_height_mm <- 250

metrics_device_path <- tempfile(pattern = "FIG5_metrics_", fileext = ".png")

ragg::agg_png(metrics_device_path, width = figure_width_mm, height = figure_height_mm, units = "mm", res = 300, background = "white")

composite_parts <- capture_stage("gtable_and_composite_assembly", {
    body_A <- build_body_gtable(panel_A)
    body_B <- build_body_gtable(panel_B)
    body_C <- build_body_gtable(panel_C)
    body_D <- build_body_gtable(panel_D)
    body_E <- build_body_gtable(panel_E)
    header_A <- make_header_gtable(body_A, "A", "GSE250130 LH+7 donor effects", "0/24 program and 0/6 VDR models at FDR < 0.05")
    header_B <- make_header_gtable(body_B, "B", "GSE250130 LH+7 differential expression", "Primary: 0 rows; fertile-reference sensitivity: 28 rows")
    header_C <- make_header_gtable(body_C, "C", "GSE183837 broad-cell landscape", "47,161 cells \u00b7 9 donors \u00b7 deterministic display subset", 
        header_shift_mm = 6)
    header_D <- make_header_gtable(body_D, "D", "GSE183837 donor-level effects", "25 prespecified models \u00b7 1 nominal exact P < 0.05 \u00b7 0 global-BH")
    header_E <- make_header_gtable(body_E, "E", "Ciliated epithelial profiles", "Descriptive only \u00b7 0/3 control and 2/6 RIF donors pass the primary cell threshold")
    design <- "\nKK\nAA\nBB\n##\nCC\nDD\n##\nEF\nGH\n##\nII\nJJ\nLL\n"
    figure_object <- (patchwork::wrap_elements(full = header_A$gtable) + patchwork::wrap_elements(full = body_A) + patchwork::wrap_elements(full = header_B$gtable) + 
        patchwork::wrap_elements(full = body_B) + patchwork::wrap_elements(full = header_C$gtable) + patchwork::wrap_elements(full = header_D$gtable) + 
        patchwork::wrap_elements(full = body_C) + patchwork::wrap_elements(full = body_D) + patchwork::wrap_elements(full = header_E$gtable) + 
        patchwork::wrap_elements(full = body_E) + patchwork::plot_spacer() + patchwork::plot_spacer()) + patchwork::plot_layout(design = design, 
        widths = c(0.81999999999999995, 1.1799999999999999), heights = grid::unit(c(2, 8.3000000000000007, 44, 2.5, 8.3000000000000007, 
            21, 2.5, 8.3000000000000007, 50, 2.5, 8.3000000000000007, 43, 2), "mm")) + patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 
        3, 2, 3, unit = "mm")))
    list(figure = figure_object, header_geometry = do.call(rbind, lapply(list(header_A, header_B, header_C, header_D, header_E), 
        `[[`, "geometry")))
})

grDevices::dev.off()
figure_export(composite_parts$figure,'Figure5',180,250)
