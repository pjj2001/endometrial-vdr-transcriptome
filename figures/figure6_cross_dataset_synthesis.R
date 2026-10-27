#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="figure6_cross_dataset_synthesis.R" & file.exists(.paths)]
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
figure_preflight('Figure6',c("data.table","ggplot2","patchwork","gtable","ragg","systemfonts","scales","R.utils"))
if ('--check' %in% commandArgs(TRUE)) {message('Figure6: inputs, packages and font available');quit(status=0)}
source_paths <- setNames(vapply(figure_index$key[figure_index$figure=="Figure6"],function(k)figure_input("Figure6",k),character(1)),figure_index$key[figure_index$figure=="Figure6"])

read_frozen <- function(path) data.table::fread(path, check.names = FALSE, encoding = "UTF-8")

frozen <- lapply(source_paths, read_frozen)

cell_display_order <- c("Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells")

endpoint_order <- c("VDR_expression", "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")

endpoint_labels <- c(VDR_expression = "VDR\nexpression", RNA_ALL_436 = "RNA all\n436", OCCUPANCY_SUPPORTED_284 = "Occupancy\nsupported 284", 
    SHARED_OCCUPANCY_46 = "Shared\noccupancy 46", PROMOTER_OCCUPANCY_84 = "Promoter\noccupancy 84")

status_levels <- c("CONCORDANT", "DISCORDANT", "NONINFORMATIVE")

status_labels <- c(CONCORDANT = "Concordant direction", DISCORDANT = "Discordant direction", NONINFORMATIVE = "Non-informative")

status_colors <- c(CONCORDANT = "#3E7C73", DISCORDANT = "#D58655", NONINFORMATIVE = "#D7D7D7")

dir_symbol <- function(x) {
    out <- rep("0", length(x))
    out[x == "RIF_HIGHER"] <- "+"
    out[x == "RIF_LOWER"] <- "-"
    out
}

plot_data <- capture_stage("frozen_source_transform", {
    concordance <- data.table::copy(frozen$concordance)
    numeric_cols <- c("Hedges_g_G250", "Hedges_g_G183", "exact_P_G250", "primary_q_G250", "exact_P_G183", "primary_q_G183")
    concordance[, `:=`((numeric_cols), lapply(.SD, as.numeric)), .SDcols = numeric_cols]
    concordance[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(cell_display_order)))]
    concordance[, `:=`(endpoint, factor(endpoint, levels = endpoint_order))]
    concordance[, `:=`(direction_status, factor(direction_status, levels = status_levels))]
    concordance[, `:=`(direction_pair, paste0("250 ", dir_symbol(direction_G250), "\n183 ", dir_symbol(direction_G183)))]
    concordance[, `:=`(label_colour, ifelse(direction_status == "NONINFORMATIVE", "dark", "light"))]
    concordance[, `:=`(nominal_any, nominal_significant_G250 | nominal_significant_G183)]
    concordance[, `:=`(endpoint_kind, ifelse(endpoint == "VDR_expression", "VDR expression", "VDR program"))]
    scatter <- concordance[!is.na(Hedges_g_G250) & !is.na(Hedges_g_G183)]
    label_points <- scatter[model_key %in% c("Myeloid cells | VDR_expression", "Stromal cells | VDR_expression")]
    label_points[, `:=`(label_x = ifelse(grepl("Myeloid", model_key), 0.56000000000000005, -1.4199999999999999), label_y = ifelse(grepl("Myeloid", 
        model_key), -1.4199999999999999, 0.40000000000000002), plot_label = ifelse(grepl("Myeloid", model_key), "Myeloid VDR\nnominal in G183", 
        "Stromal VDR\nnominal in G250"))]
    myeloid <- data.table::copy(frozen$myeloid_vdr)
    myeloid_plot <- data.table::data.table(dataset = factor(c("GSE183837", "GSE250130"), levels = c("GSE183837", "GSE250130")), 
        y = c(1, 2), hedges_g = c(as.numeric(myeloid$Hedges_g_G183[[1]]), as.numeric(myeloid$Hedges_g_G250[[1]])), exact_p = c(as.numeric(myeloid$exact_P_G183[[1]]), 
            as.numeric(myeloid$exact_P_G250[[1]])), primary_q = c(as.numeric(myeloid$primary_q_G183[[1]]), as.numeric(myeloid$primary_q_G250[[1]])), 
        dataset_colour = c("G183", "G250"))
    myeloid_plot[, `:=`(annotation, c("P=0.0357; q=0.893", "P=0.570; q=0.570"))]
    facts <- setNames(as.character(frozen$locked_facts$value), frozen$locked_facts$fact_id)
    evidence_boxes <- data.table::data.table(box = c("anchor", "temporal", "bulk", "disease"), xmin = c(1, 26, 51, 76), xmax = c(24, 
        49, 74, 99), ymin = 58, ymax = 99, fill = c("#EEF1F4", "#DCEFE8", "#F8E8CF", "#F3DDDA"), title = c("Mechanistic anchor", 
        "Normal remodeling", "Bulk projection", "Disease replication"), detail = c(paste0(facts[["F001"]], " RNA candidates\n", 
        facts[["F002"]], " occupancy-supported\n", facts[["F003"]], " shared occupancy"), paste0(facts[["F009"]], "/", facts[["F008"]], 
        " program models FDR+\n", facts[["F010"]], " temporal DE rows\nStrong positive evidence"), paste0(facts[["F007"]], 
        " FDR+ program rows\nacross ", facts[["F006"]], " contrasts\nContext-specific"), paste0("G250: 0 corrected\nG183: 0 global BH\nCross-cohort: ", 
        facts[["F022"]], " replicated")))
    final_box <- data.table::data.table(xmin = 8, xmax = 92, ymin = 3, ymax = 32, title = "Integrated evidence boundary", 
        detail = "No broad, reproducible RIF-associated disruption of the VDR programs\nNo pooling; no direct-target, causal, or biomarker claim")
    list(A = concordance, B = scatter, B_labels = label_points, C = myeloid_plot, D = evidence_boxes, D_final = final_box, 
        facts = facts)
})

metric <- setNames(as.character(frozen$overall$value), frozen$overall$metric)

interpretation <- setNames(as.character(frozen$interpretation_lock$value), frozen$interpretation_lock$item)

scientific_checks <- data.frame(check_id = sprintf("S%02d", 1:17), check = c("Planned common models", "Informative models", 
    "Direction-concordant models", "Direction-discordant models", "Non-informative models", "G250 nominal models", "G250 corrected models", 
    "G183 nominal models", "G183 corrected models", "Corrected replicated models", "Exploratory directional-support models", 
    "Myeloid VDR G250 direction", "Myeloid VDR G183 direction", "Myeloid VDR G183 exact P", "Myeloid VDR G183 global q", 
    "Broad reproducible dysregulation supported", "Pooled meta-analysis performed"), observed = c(nrow(plot_data$A), nrow(plot_data$B), 
    sum(plot_data$A$direction_status == "CONCORDANT"), sum(plot_data$A$direction_status == "DISCORDANT"), sum(plot_data$A$direction_status == 
        "NONINFORMATIVE"), sum(plot_data$A$nominal_significant_G250), sum(plot_data$A$corrected_significant_G250), sum(plot_data$A$nominal_significant_G183), 
    sum(plot_data$A$corrected_significant_G183), metric[["corrected_replicated_concordant_n"]], metric[["exploratory_directional_support_n"]], 
    as.character(frozen$myeloid_vdr$direction_G250[[1]]), as.character(frozen$myeloid_vdr$direction_G183[[1]]), sprintf("%.7f", 
        as.numeric(frozen$myeloid_vdr$exact_P_G183[[1]])), sprintf("%.6f", as.numeric(frozen$myeloid_vdr$primary_q_G183[[1]])), 
    interpretation[["broad_RIF_VDR_program_dysregulation_supported"]], interpretation[["pooled_meta_analysis_performed"]]), 
    expected = c(25, 24, 12, 12, 1, 1, 0, 1, 0, "0", "0", "RIF_HIGHER", "RIF_LOWER", "0.0357143", "0.892857", "FALSE", "FALSE"), 
    stringsAsFactors = FALSE)

scientific_checks$pass <- as.character(scientific_checks$observed) == as.character(scientific_checks$expected)

if (!all(scientific_checks$pass)) stop("Figure 6 scientific QC failed.")

color_text <- "#161616"

color_muted <- "#555555"

color_grid <- "#D9D9D9"

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2000000000000002) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major = ggplot2::element_line(linewidth = 0.23000000000000001, colour = color_grid), axis.title = ggplot2::element_text(size = 6.7000000000000002, 
        colour = color_text), axis.text = ggplot2::element_text(size = 5.7000000000000002, colour = color_text), legend.title = ggplot2::element_text(size = 5.5999999999999996, 
        face = "bold"), legend.text = ggplot2::element_text(size = 5.2000000000000002), plot.margin = ggplot2::margin(1.2, 
        1.5, 1.2, 1.5, unit = "mm"))

panel_A <- capture_stage("panel_A_construction", {
    ggplot2::ggplot(plot_data$A, ggplot2::aes(x = endpoint, y = broad_cell_type, fill = direction_status)) + ggplot2::geom_tile(colour = "white", 
        linewidth = 0.80000000000000004, width = 0.97999999999999998, height = 0.95999999999999996) + ggplot2::geom_text(ggplot2::aes(label = direction_pair, 
        colour = label_colour), family = figure_font, fontface = "bold", size = 2.0499999999999998, lineheight = 0.90000000000000002) + 
        ggplot2::scale_fill_manual(values = status_colors, labels = status_labels, drop = FALSE) + ggplot2::scale_colour_manual(values = c(dark = color_text, 
        light = "white"), guide = "none") + ggplot2::scale_x_discrete(labels = endpoint_labels, expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 
        0)) + ggplot2::labs(x = NULL, y = NULL, fill = "Cross-dataset status") + ggplot2::theme_void(base_family = figure_font, 
        base_size = 7.2000000000000002) + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5.7000000000000002, colour = color_text, 
        lineheight = 0.90000000000000002, margin = ggplot2::margin(t = 2)), axis.text.y = ggplot2::element_text(size = 6, 
        colour = color_text, margin = ggplot2::margin(r = 2)), legend.position = "bottom", legend.direction = "horizontal", 
        legend.title.position = "left", legend.key.width = grid::unit(4.2000000000000002, "mm"), legend.key.height = grid::unit(2.7999999999999998, 
            "mm"), legend.margin = ggplot2::margin(t = 0.5, b = 0), plot.margin = ggplot2::margin(1, 1.5, 0.5, 1.5, unit = "mm"))
})

panel_B <- capture_stage("panel_B_construction", {
    ggplot2::ggplot(plot_data$B, ggplot2::aes(x = Hedges_g_G250, y = Hedges_g_G183)) + ggplot2::geom_abline(slope = 1, intercept = 0, 
        colour = "#A6A6A6", linewidth = 0.34999999999999998, linetype = 3) + ggplot2::geom_vline(xintercept = 0, colour = "#777777", 
        linewidth = 0.29999999999999999, linetype = 2) + ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.29999999999999999, 
        linetype = 2) + ggplot2::geom_point(ggplot2::aes(fill = direction_status, shape = nominal_any), size = 2, stroke = 0.45000000000000001, 
        colour = color_text) + ggplot2::geom_segment(data = plot_data$B_labels, ggplot2::aes(x = Hedges_g_G250, y = Hedges_g_G183, 
        xend = label_x, yend = label_y), inherit.aes = FALSE, linewidth = 0.28000000000000003, colour = "#6D6D6D") + ggplot2::geom_label(data = plot_data$B_labels, 
        ggplot2::aes(x = label_x, y = label_y, label = plot_label), inherit.aes = FALSE, family = figure_font, size = 1.75, 
        lineheight = 0.90000000000000002, linewidth = 0.17999999999999999, label.padding = grid::unit(0.55000000000000004, 
            "mm"), fill = scales::alpha("white", 0.92000000000000004), colour = color_text) + ggplot2::scale_fill_manual(values = status_colors, 
        guide = "none") + ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") + ggplot2::scale_x_continuous(limits = c(-1.8, 
        0.80000000000000004), breaks = c(-1.5, -1, -0.5, 0, 0.5)) + ggplot2::scale_y_continuous(limits = c(-1.8, 0.80000000000000004), 
        breaks = c(-1.5, -1, -0.5, 0, 0.5)) + ggplot2::labs(x = "GSE250130 Hedges g", y = "GSE183837 Hedges g") + ggplot2::coord_cartesian(clip = "off") + 
        base_theme + ggplot2::theme(legend.position = "none", panel.grid.major = ggplot2::element_line(linewidth = 0.20000000000000001, 
        colour = color_grid), axis.text = ggplot2::element_text(size = 5.5), axis.title = ggplot2::element_text(size = 6.2000000000000002), 
        plot.margin = ggplot2::margin(1.5, 2.5, 1.2, 1.5, unit = "mm"))
})

panel_C <- capture_stage("panel_C_construction", {
    ggplot2::ggplot(plot_data$C, ggplot2::aes(x = hedges_g, y = y, colour = dataset_colour)) + ggplot2::geom_vline(xintercept = 0, 
        colour = "#777777", linewidth = 0.32000000000000001, linetype = 2) + ggplot2::geom_segment(ggplot2::aes(x = 0, xend = hedges_g, 
        yend = y), linewidth = 0.75, alpha = 0.80000000000000004) + ggplot2::geom_point(size = 2.7999999999999998) + ggplot2::geom_text(ggplot2::aes(label = annotation), 
        nudge_y = 0.20000000000000001, family = figure_font, size = 1.8999999999999999, colour = color_text) + ggplot2::annotate("label", 
        x = -0.68999999999999995, y = 1.5, label = "Opposite directions\nNo cross-dataset support", family = figure_font, 
        fontface = "bold", size = 1.8999999999999999, lineheight = 0.92000000000000004, linewidth = 0.20000000000000001, 
        label.padding = grid::unit(0.69999999999999996, "mm"), fill = "#F7F7F7", colour = color_text) + ggplot2::scale_colour_manual(values = c(G183 = "#D58655", 
        G250 = "#3E7C73"), guide = "none") + ggplot2::scale_y_continuous(breaks = c(1, 2), labels = c("GSE183837", "GSE250130"), 
        limits = c(0.71999999999999997, 2.2799999999999998)) + ggplot2::scale_x_continuous(limits = c(-2.0499999999999998, 
        0.65000000000000002), breaks = c(-2, -1.5, -1, -0.5, 0, 0.5)) + ggplot2::labs(x = "Hedges g (RIF - reference)", y = NULL) + 
        base_theme + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "none", axis.text.x = ggplot2::element_text(size = 5.2000000000000002), 
        axis.text.y = ggplot2::element_text(size = 5.7999999999999998), axis.title.x = ggplot2::element_text(size = 6.0999999999999996), 
        plot.margin = ggplot2::margin(1.5, 2, 1.2, 1.5, unit = "mm"))
})

panel_D <- capture_stage("panel_D_construction", {
    arrow_data <- data.frame(x = c(12.5, 37.5, 62.5, 87.5), y = 57, xend = c(34, 45, 55, 66), yend = 33)
    ggplot2::ggplot() + ggplot2::geom_segment(data = arrow_data, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.38, 
        colour = "#8A8A8A", arrow = grid::arrow(length = grid::unit(1.6000000000000001, "mm"), type = "closed")) + ggplot2::geom_rect(data = plot_data$D, 
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), colour = "#A6A6A6", linewidth = 0.38) + 
        ggplot2::scale_fill_identity() + ggplot2::geom_text(data = plot_data$D, ggplot2::aes(x = (xmin + xmax)/2, y = 85, 
        label = title), family = figure_font, fontface = "bold", size = 2.3500000000000001, colour = color_text) + ggplot2::geom_text(data = plot_data$D, 
        ggplot2::aes(x = (xmin + xmax)/2, y = 70, label = detail), family = figure_font, size = 1.95, lineheight = 1.02, 
        colour = color_text) + ggplot2::geom_rect(data = plot_data$D_final, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, 
        ymax = ymax), fill = "#E8EEF3", colour = "#526778", linewidth = 0.55000000000000004) + ggplot2::geom_text(data = plot_data$D_final, 
        ggplot2::aes(x = (xmin + xmax)/2, y = 24, label = title), family = figure_font, fontface = "bold", size = 2.5499999999999998, 
        colour = "#2A4152") + ggplot2::geom_text(data = plot_data$D_final, ggplot2::aes(x = (xmin + xmax)/2, y = 13.5, label = detail), 
        family = figure_font, size = 2.1000000000000001, lineheight = 1.04, colour = color_text) + ggplot2::coord_cartesian(xlim = c(0, 
        100), ylim = c(0, 100), clip = "off") + ggplot2::theme_void(base_family = figure_font) + ggplot2::theme(plot.margin = ggplot2::margin(-1.3999999999999999, 
        1.5, 1.2, 1.5, unit = "mm"))
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

header_subtitle_size_pt <- 6.2999999999999998

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
        clip = "off", name = paste0("F6-header-title-", panel_tag))
    header_gt <- gtable::gtable_add_grob(header_gt, subtitle_grob, t = 2, b = 2, l = panel_left, r = panel_right, clip = "off", 
        name = paste0("F6-header-subtitle-", panel_tag))
    list(gtable = header_gt, geometry = data.frame(panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE", panel_left_column = panel_left, 
        panel_right_column = panel_right, tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt, subtitle_size_pt = header_subtitle_size_pt, 
        tag_title_gap_mm = header_tag_to_title_gap_mm, header_shift_mm = header_shift_mm, stringsAsFactors = FALSE))
}

figure_width_mm <- 180

figure_height_mm <- 212

metrics_device_path <- tempfile(pattern = "FIG6_metrics_", fileext = ".png")

ragg::agg_png(metrics_device_path, width = figure_width_mm, height = figure_height_mm, units = "mm", res = 300, background = "white")

composite_parts <- capture_stage("gtable_and_composite_assembly", {
    body_A <- build_body_gtable(panel_A)
    body_B <- build_body_gtable(panel_B)
    body_C <- build_body_gtable(panel_C)
    body_D <- build_body_gtable(panel_D)
    header_A <- make_header_gtable(body_A, "A", "Direction concordance across 25 shared models", "250/183: effect direction in GSE250130/GSE183837; 12 concordant, 12 discordant, 1 non-informative")
    header_B <- make_header_gtable(body_B, "B", "Standardized effects across cohorts", "24 informative models; color follows Panel A; triangles mark nominal P < 0.05 in one cohort")
    header_C <- make_header_gtable(body_C, "C", "Myeloid VDR is directionally discordant", "Nominal only in GSE183837; no corrected support")
    header_D <- make_header_gtable(body_D, "D", "Integrated evidence states and inference boundary", "Positive normal remodeling and context-specific bulk signals remain distinct from negative disease replication", 
        header_shift_mm = 8.4000000000000004)
    design <- "\nAA\nBB\n##\nCD\nEF\n##\nGG\nHH\n"
    figure_object <- (patchwork::wrap_elements(full = header_A$gtable) + patchwork::wrap_elements(full = body_A) + patchwork::wrap_elements(full = header_B$gtable) + 
        patchwork::wrap_elements(full = header_C$gtable) + patchwork::wrap_elements(full = body_B) + patchwork::wrap_elements(full = body_C) + 
        patchwork::wrap_elements(full = header_D$gtable) + patchwork::wrap_elements(full = body_D)) + patchwork::plot_layout(design = design, 
        widths = c(1.1799999999999999, 0.81999999999999995), heights = grid::unit(c(8.3000000000000007, 52, 2.5, 8.3000000000000007, 
            50, 2.5, 8.3000000000000007, 56), "mm")) + patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 
        3, 2, 3, unit = "mm")))
    list(figure = figure_object, header_geometry = do.call(rbind, lapply(list(header_A, header_B, header_C, header_D), `[[`, 
        "geometry")))
})

grDevices::dev.off()
figure_export(composite_parts$figure,'Figure6',180,212)
