#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="figure3_bulk_projection.R" & file.exists(.paths)]
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
figure_preflight('Figure3',c("data.table","ggplot2","patchwork","gtable","ragg","systemfonts","scales","R.utils"))
if ('--check' %in% commandArgs(TRUE)) {message('Figure3: inputs, packages and font available');quit(status=0)}
source_root <- file.path(figure_data_root,'Figure3')
.metric_file <- tempfile(fileext='.png')
ragg::agg_png(.metric_file,width=180,height=184,units='mm',res=300,background='white')

score_models <- read.csv(file.path(source_root, "Figure3_panelA_program_effect_heatmap.csv"), stringsAsFactors = FALSE, check.names = FALSE, 
    fileEncoding = "UTF-8-BOM")

coverage <- read.csv(file.path(source_root, "Figure3_panelC_candidate_context_coverage.csv"), stringsAsFactors = FALSE, check.names = FALSE, 
    fileEncoding = "UTF-8-BOM")

panel_D_data <- read.csv(file.path(source_root, "Figure3_panelD_context_effects_long.csv"), stringsAsFactors = FALSE, check.names = FALSE, 
    fileEncoding = "UTF-8-BOM")

selected_examples <- read.csv(file.path(source_root, "Figure3_panelD_selected_examples.csv"), stringsAsFactors = FALSE, check.names = FALSE, 
    fileEncoding = "UTF-8-BOM")

if (nrow(score_models) != 16L) stop("Expected 16 program-by-contrast rows.")

if (sum(as.logical(score_models$FDR_significant)) != 4L) stop("Expected 4 FDR-significant program rows.")

if (nrow(coverage) != 5L || sum(coverage$candidate_n) != 436L) stop("Coverage table mismatch.")

if (nrow(selected_examples) != 10L || nrow(panel_D_data) != 40L) stop("Panel D source mismatch.")

program_levels <- c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")

program_labels <- c(RNA_ALL_436 = "RNA program\n(n=436)", OCCUPANCY_SUPPORTED_284 = "Occupancy-supported\n(n=284)", SHARED_OCCUPANCY_46 = "Shared occupancy\n(n=46)", 
    PROMOTER_OCCUPANCY_84 = "Promoter-associated\n(n=84)")

context_ids <- c("GSE165004__RPL_vs_CONTROL", "GSE165004__UIF_vs_CONTROL", "GSE92324__IF_vs_CONTROL", "GSE111974__RIF_vs_CONTROL")

context_labels <- c(GSE165004__RPL_vs_CONTROL = "GSE165004\nRPL", GSE165004__UIF_vs_CONTROL = "GSE165004\nUIF", GSE92324__IF_vs_CONTROL = "GSE92324\nIF", 
    GSE111974__RIF_vs_CONTROL = "GSE111974\nRIF")

score_models$program_label <- factor(unname(program_labels[score_models$program_name]), levels = unname(program_labels[rev(program_levels)]))

score_models$context_label <- factor(unname(context_labels[score_models$context_id]), levels = unname(context_labels[context_ids]))

score_models$logFC <- as.numeric(score_models$logFC)

score_models$adj.P.Val <- as.numeric(score_models$adj.P.Val)

score_models$FDR_significant <- as.logical(score_models$FDR_significant)

score_models$neg_log10_FDR <- -log10(pmax(score_models$adj.P.Val, 1e-300))

score_models$display_label <- paste0(formatC(score_models$logFC, format = "f", digits = 2), ifelse(score_models$FDR_significant, 
    "*", ""))

coverage$coverage_label <- factor(coverage$measured_context_n, levels = 0:4)

panel_D_data$context_label <- factor(unname(context_labels[panel_D_data$context_id]), levels = unname(context_labels[context_ids]))

ordered_genes <- selected_examples$feature_id[order(selected_examples$selection_rank)]

panel_D_data$feature_id <- factor(panel_D_data$feature_id, levels = rev(ordered_genes))

panel_D_data$logFC <- as.numeric(panel_D_data$logFC)

panel_D_data$FDR_significant <- as.logical(panel_D_data$FDR_significant)

if (!"display_label" %in% names(panel_D_data)) {
    panel_D_data$display_label <- paste0(formatC(panel_D_data$logFC, format = "f", digits = 2), ifelse(panel_D_data$FDR_significant, 
        "*", ""))
}

header_tag_size_pt <- 11

header_title_size_pt <- 9.1999999999999993

header_subtitle_size_pt <- 7.2000000000000002

header_tag_to_title_gap_mm <- 1.6000000000000001

header_title_row_height_mm <- 5.4000000000000004

header_subtitle_row_height_mm <- 4

figure_width_mm <- 180

figure_height_mm <- 184

color_negative <- "#3F6FA8"

color_neutral <- "#F5F5F2"

color_positive <- "#B65B5B"

color_bar <- "#6C8198"

color_text <- "#161616"

color_muted <- "#505050"

color_grid <- "#D7D7D7"

program_effect_limit <- max(0.5, ceiling(max(abs(score_models$logFC), na.rm = TRUE) * 10)/10)

gene_effect_limit <- max(1, ceiling(max(abs(panel_D_data$logFC), na.rm = TRUE) * 10)/10)

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2000000000000002) + ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.25, 
    colour = color_grid), panel.grid.minor = ggplot2::element_blank(), axis.title = ggplot2::element_text(family = figure_font, 
    size = 8, colour = color_text), axis.text = ggplot2::element_text(family = figure_font, size = 6.7999999999999998, colour = color_text), 
    legend.title = ggplot2::element_text(family = figure_font, size = 7, face = "bold", colour = color_text), legend.text = ggplot2::element_text(family = figure_font, 
        size = 6.7999999999999998, colour = color_text), plot.margin = ggplot2::margin(3, 4, 3, 4))

panel_A <- ggplot2::ggplot(score_models, ggplot2::aes(x = context_label, y = program_label, fill = logFC)) + ggplot2::geom_tile(linewidth = 0.45000000000000001, 
    colour = "white") + ggplot2::geom_text(ggplot2::aes(label = display_label), family = figure_font, size = 2.4500000000000002, 
    colour = color_text) + ggplot2::scale_fill_gradient2(low = color_negative, mid = color_neutral, high = color_positive, 
    midpoint = 0, limits = c(-program_effect_limit, program_effect_limit), oob = scales::squish, name = "Program logFC") + 
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::coord_cartesian(clip = "off") + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), 
    axis.text.x = ggplot2::element_text(family = figure_font, size = 6.5999999999999996, lineheight = 0.92000000000000004, 
        margin = ggplot2::margin(t = 3)), axis.text.y = ggplot2::element_text(family = figure_font, size = 6.5, lineheight = 0.93999999999999995, 
        margin = ggplot2::margin(r = 3)), legend.position = "none")

panel_B <- ggplot2::ggplot(score_models, ggplot2::aes(x = context_label, y = program_label)) + ggplot2::geom_point(ggplot2::aes(size = neg_log10_FDR, 
    fill = logFC), shape = 21, colour = "#9B9B9B", stroke = 0.34999999999999998) + ggplot2::geom_point(data = score_models[score_models$FDR_significant, 
    , drop = FALSE], ggplot2::aes(size = neg_log10_FDR, fill = logFC), shape = 21, colour = "#171717", stroke = 0.90000000000000002) + 
    ggplot2::scale_fill_gradient2(low = color_negative, mid = color_neutral, high = color_positive, midpoint = 0, limits = c(-program_effect_limit, 
        program_effect_limit), oob = scales::squish, name = "Program logFC") + ggplot2::scale_size_continuous(range = c(2.1000000000000001, 
    6.7000000000000002), breaks = c(1, 2, 4), name = expression(-log[10](FDR))) + ggplot2::labs(x = NULL, y = NULL) + base_theme + 
    ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.34999999999999998, colour = "#E1E1E1"), axis.text.x = ggplot2::element_text(family = figure_font, 
        size = 6.5999999999999996, lineheight = 0.92000000000000004, margin = ggplot2::margin(t = 3)), axis.text.y = ggplot2::element_text(family = figure_font, 
        size = 6.5, lineheight = 0.93999999999999995, margin = ggplot2::margin(r = 3)), legend.position = "bottom", legend.box = "horizontal", 
        legend.justification = "center", legend.spacing.x = grid::unit(1.2, "mm"), legend.key.width = grid::unit(7, "mm"), 
        legend.margin = ggplot2::margin(t = 1, b = 0)) + ggplot2::guides(fill = "none", size = ggplot2::guide_legend(title.position = "left", 
    nrow = 1, byrow = TRUE, override.aes = list(fill = "#D6D6D6", colour = "#171717")))

panel_C <- ggplot2::ggplot(coverage, ggplot2::aes(x = coverage_label, y = candidate_n)) + ggplot2::geom_col(width = 0.68000000000000005, 
    fill = color_bar) + ggplot2::geom_text(ggplot2::aes(label = format(candidate_n, big.mark = ",", scientific = FALSE, trim = TRUE)), 
    family = figure_font, size = 2.5499999999999998, vjust = -0.45000000000000001, colour = color_text) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 
    0.14000000000000001)), labels = scales::label_number(big.mark = ",")) + ggplot2::labs(x = "Number of patient bulk contexts measured", 
    y = "Candidate genes") + base_theme + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(family = figure_font, 
    size = 7))

panel_D <- ggplot2::ggplot(panel_D_data, ggplot2::aes(x = context_label, y = feature_id, fill = logFC)) + ggplot2::geom_tile(linewidth = 0.45000000000000001, 
    colour = "white") + ggplot2::geom_text(ggplot2::aes(label = display_label), family = figure_font, size = 2.25, colour = color_text) + 
    ggplot2::scale_fill_gradient2(low = color_negative, mid = color_neutral, high = color_positive, midpoint = 0, limits = c(-gene_effect_limit, 
        gene_effect_limit), oob = scales::squish, name = "Gene logFC") + ggplot2::labs(x = NULL, y = NULL) + base_theme + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(family = figure_font, size = 6.5, 
        lineheight = 0.92000000000000004, margin = ggplot2::margin(t = 3)), axis.text.y = ggplot2::element_text(family = figure_font, 
        size = 6.7000000000000002, face = "italic", margin = ggplot2::margin(r = 3)), legend.position = "none")

panel_A_guide_reference <- panel_A + ggplot2::theme(legend.position = "bottom", legend.key.width = grid::unit(12, "mm"), 
    legend.margin = ggplot2::margin(t = 1, b = 0)) + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "left", 
    title.hjust = 0.5, barheight = grid::unit(2.6000000000000001, "mm")))

panel_D_guide_reference <- panel_D + ggplot2::theme(legend.position = "bottom", legend.key.width = grid::unit(12, "mm")) + 
    ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "left", barheight = grid::unit(2.6000000000000001, "mm")))

find_colorbar_guide_gtable <- function(g) {
    if (inherits(g, "gtable") && any(g$layout$name == "bar") && any(g$layout$name == "title")) {
        return(g)
    }
    if (!is.null(g$grobs)) {
        for (child in g$grobs) {
            ans <- find_colorbar_guide_gtable(child)
            if (!is.null(ans)) 
                return(ans)
        }
    }
    NULL
}

measure_data_body_right_fraction <- function(plot_object, panel_tag) {
    built <- ggplot2::ggplot_build(plot_object)
    if (!length(built$layout$panel_params)) {
        stop("No panel parameters found for Panel ", panel_tag)
    }
    x_range <- built$layout$panel_params[[1]]$x.range
    if (length(x_range) != 2L || !all(is.finite(x_range))) {
        stop("Could not resolve x-range for Panel ", panel_tag)
    }
    tile_layer_idx <- which(vapply(built$data, function(d) all(c("xmin", "xmax") %in% names(d)), logical(1)))
    if (!length(tile_layer_idx)) {
        stop("Could not identify geom_tile layer for Panel ", panel_tag)
    }
    tile_dat <- built$data[[tile_layer_idx[[1]]]]
    tile_right <- max(tile_dat$xmax, na.rm = TRUE)
    tile_left <- min(tile_dat$xmin, na.rm = TRUE)
    right_frac <- (tile_right - x_range[[1]])/(x_range[[2]] - x_range[[1]])
    left_frac <- (tile_left - x_range[[1]])/(x_range[[2]] - x_range[[1]])
    if (!is.finite(right_frac) || right_frac <= 0 || right_frac >= 1 || !is.finite(left_frac) || left_frac <= 0 || left_frac >= 
        right_frac) {
        stop("Invalid tile-body fractions for Panel ", panel_tag)
    }
    data.frame(panel = panel_tag, x_range_min = x_range[[1]], x_range_max = x_range[[2]], tile_left = tile_left, tile_right = tile_right, 
        tile_left_fraction = left_frac, tile_right_fraction = right_frac, stringsAsFactors = FALSE)
}

measure_color_scale <- function(plot_object, panel_tag, approved_bar_width_mm) {
    gt <- ggplot2::ggplotGrob(plot_object)
    guide_idx <- which(gt$layout$name == "guide-box-bottom")
    if (!length(guide_idx)) {
        guide_idx <- which(grepl("^guide-box", gt$layout$name))
    }
    if (!length(guide_idx)) {
        stop("Could not locate guide box for Panel ", panel_tag)
    }
    guide_box <- gt$grobs[[guide_idx[[1]]]]
    guide_gt <- find_colorbar_guide_gtable(guide_box)
    if (is.null(guide_gt)) {
        stop("Could not locate colorbar guide gtable for Panel ", panel_tag)
    }
    bar_cell <- guide_gt$layout[guide_gt$layout$name == "bar", , drop = FALSE]
    if (nrow(bar_cell) != 1L) {
        stop("Expected exactly one colorbar body for Panel ", panel_tag)
    }
    bar_height_mm <- grid::convertHeight(sum(guide_gt$heights[bar_cell$t:bar_cell$b]), "mm", valueOnly = TRUE)
    if (!is.finite(bar_height_mm) || bar_height_mm <= 0) {
        stop("Invalid colorbar height for Panel ", panel_tag)
    }
    built <- ggplot2::ggplot_build(plot_object)
    fill_scale <- built$plot$scales$get_scales("fill")
    if (is.null(fill_scale)) {
        stop("Could not locate fill scale for Panel ", panel_tag)
    }
    limits <- fill_scale$get_limits()
    breaks <- fill_scale$get_breaks()
    breaks <- breaks[is.finite(breaks) & breaks >= min(limits) & breaks <= max(limits)]
    labels <- fill_scale$get_labels(breaks)
    raster_values <- seq(limits[[1]], limits[[2]], length.out = 512L)
    raster_colors <- fill_scale$map(raster_values)
    list(panel = panel_tag, approved_bar_width_mm = approved_bar_width_mm, bar_height_mm = bar_height_mm, limits = limits, 
        breaks = breaks, labels = labels, raster_colors = raster_colors)
}

make_data_body_aligned_colorbar_grob <- function(title, metrics, tile_right_fraction, figure_font, title_size_pt = 7, label_size_pt = 6.7999999999999998, 
    color_text = "#161616", title_to_bar_gap_mm = 2) {
    bar_w <- grid::unit(metrics$approved_bar_width_mm, "mm")
    bar_h <- grid::unit(metrics$bar_height_mm, "mm")
    bar_bottom <- grid::unit(4.0999999999999996, "mm")
    bar_top <- bar_bottom + bar_h
    bar_right <- grid::unit(tile_right_fraction, "npc")
    bar_left <- bar_right - bar_w
    gradient_matrix <- matrix(metrics$raster_colors, nrow = 1L)
    bar <- grid::rasterGrob(image = gradient_matrix, x = bar_left + 0.5 * bar_w, y = bar_bottom + 0.5 * bar_h, width = bar_w, 
        height = bar_h, interpolate = TRUE)
    frac <- (metrics$breaks - metrics$limits[[1]])/(metrics$limits[[2]] - metrics$limits[[1]])
    tick_x <- bar_left + frac * bar_w
    tick_grobs <- lapply(seq_along(tick_x), function(i) {
        grid::segmentsGrob(x0 = tick_x[[i]], x1 = tick_x[[i]], y0 = bar_bottom, y1 = bar_top, gp = grid::gpar(col = "white", 
            lwd = 0.45000000000000001))
    })
    label_grobs <- lapply(seq_along(tick_x), function(i) {
        if (i == 1L) {
            just_x <- "left"
        }
        else if (i == length(tick_x)) {
            just_x <- "right"
        }
        else {
            just_x <- "center"
        }
        grid::textGrob(as.character(metrics$labels[[i]]), x = tick_x[[i]], y = grid::unit(1.6499999999999999, "mm"), just = c(just_x, 
            "center"), gp = grid::gpar(fontfamily = figure_font, fontsize = label_size_pt, col = color_text))
    })
    title_grob <- grid::textGrob(title, x = bar_left - grid::unit(title_to_bar_gap_mm, "mm"), y = bar_bottom, just = c("right", 
        "bottom"), gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = title_size_pt, col = color_text))
    grid::grobTree(title_grob, bar, do.call(grid::grobTree, tick_grobs), do.call(grid::grobTree, label_grobs))
}

add_data_body_colorbar_to_body <- function(body_gtable, legend_grob, legend_height_mm, panel_tag) {
    panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
    if (!nrow(panel_cells)) {
        stop("Could not locate scientific panel cell for Panel ", panel_tag)
    }
    panel_left <- min(panel_cells$l)
    panel_right <- max(panel_cells$r)
    body_gtable <- gtable::gtable_add_rows(body_gtable, heights = grid::unit(legend_height_mm, "mm"), pos = -1)
    new_row <- length(body_gtable$heights)
    body_gtable <- gtable::gtable_add_grob(body_gtable, legend_grob, t = new_row, b = new_row, l = panel_left, r = panel_right, 
        clip = "off", name = paste0("F3v148-data-body-colorbar-", panel_tag))
    body_gtable
}

A_approved_bar_width_mm <- 50.880000000000003

D_approved_bar_width_mm <- 59.939999999999998

A_tile_geometry <- measure_data_body_right_fraction(panel_A, "A")

D_tile_geometry <- measure_data_body_right_fraction(panel_D, "D")

A_colorbar_metrics <- measure_color_scale(panel_A_guide_reference, "A", approved_bar_width_mm = A_approved_bar_width_mm)

D_colorbar_metrics <- measure_color_scale(panel_D_guide_reference, "D", approved_bar_width_mm = D_approved_bar_width_mm)

A_constraint_legend <- make_data_body_aligned_colorbar_grob(title = "Program logFC", metrics = A_colorbar_metrics, tile_right_fraction = A_tile_geometry$tile_right_fraction[[1]], 
    figure_font = figure_font, color_text = color_text)

D_constraint_legend <- make_data_body_aligned_colorbar_grob(title = "Gene logFC", metrics = D_colorbar_metrics, tile_right_fraction = D_tile_geometry$tile_right_fraction[[1]], 
    figure_font = figure_font, color_text = color_text)

build_body_gtable <- function(plot_object) {
    ggplot2::ggplotGrob(plot_object + ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) + ggplot2::theme(plot.tag = ggplot2::element_blank(), 
        plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank()))
}

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle) {
    panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
    if (!nrow(panel_cells)) 
        stop("Could not locate scientific panel cell for ", panel_tag)
    panel_left <- min(panel_cells$l)
    panel_right <- max(panel_cells$r)
    header_gt <- gtable::gtable(widths = body_gtable$widths, heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), 
        "mm"))
    title_grob <- grid::textGrob(panel_title, x = grid::unit(0, "npc"), y = grid::unit(0.52000000000000002, "npc"), hjust = 0, 
        vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_title_size_pt, lineheight = 0.93999999999999995, 
            col = color_text))
    tag_grob <- grid::textGrob(panel_tag, x = grid::unit(0, "npc") - grid::unit(header_tag_to_title_gap_mm, "mm"), y = grid::unit(0.52000000000000002, 
        "npc"), hjust = 1, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_tag_size_pt, 
        col = color_text))
    subtitle_grob <- grid::textGrob(panel_subtitle, x = grid::unit(0, "npc"), y = grid::unit(0.56000000000000005, "npc"), 
        hjust = 0, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, col = color_muted))
    header_gt <- gtable::gtable_add_grob(header_gt, grid::grobTree(tag_grob, title_grob), t = 1, b = 1, l = panel_left, r = panel_right, 
        clip = "off", name = paste0("F3v147-title-", panel_tag))
    header_gt <- gtable::gtable_add_grob(header_gt, subtitle_grob, t = 2, b = 2, l = panel_left, r = panel_right, clip = "off", 
        name = paste0("F3v147-subtitle-", panel_tag))
    list(gtable = header_gt, geometry = data.frame(panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE", panel_left_column = panel_left, 
        panel_right_column = panel_right, tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt, subtitle_size_pt = header_subtitle_size_pt, 
        tag_title_gap_mm = header_tag_to_title_gap_mm, header_composition = "DEDICATED_SHARED_ROW", stringsAsFactors = FALSE))
}

body_A <- build_body_gtable(panel_A)

body_B <- build_body_gtable(panel_B)

body_C <- build_body_gtable(panel_C)

body_D <- build_body_gtable(panel_D)

body_A <- add_data_body_colorbar_to_body(body_A, A_constraint_legend, legend_height_mm = 8, panel_tag = "A")

body_D <- add_data_body_colorbar_to_body(body_D, D_constraint_legend, legend_height_mm = 8, panel_tag = "D")

header_A <- make_header_gtable(body_A, "A", "Program effects across patient contrasts", "All contrasts vs Control \u00b7 cell = logFC; * FDR < 0.05")

header_B <- make_header_gtable(body_B, "B", "Effect magnitude and multiplicity", "vs Control \u00b7 fill = logFC; size = -log10(FDR)")

header_C <- make_header_gtable(body_C, "C", "Candidate measurement coverage", "Patient contexts with an estimable gene-level effect")

header_D <- make_header_gtable(body_D, "D", "Context-dependent gene examples", "Top 10 mixed-direction genes \u00b7 * FDR < 0.05")

header_geometry <- do.call(rbind, lapply(list(header_A, header_B, header_C, header_D), `[[`, "geometry"))

design <- "\nAB\nCD\n##\nEF\nGH\n##\n"

figure_object <- (patchwork::wrap_elements(full = header_A$gtable) + patchwork::wrap_elements(full = header_B$gtable) + patchwork::wrap_elements(full = body_A) + 
    patchwork::wrap_elements(full = body_B) + patchwork::wrap_elements(full = header_C$gtable) + patchwork::wrap_elements(full = header_D$gtable) + 
    patchwork::wrap_elements(full = body_C) + patchwork::wrap_elements(full = body_D)) + patchwork::plot_layout(design = design, 
    widths = c(1, 1), heights = grid::unit(c(10, 72, 4, 10, 70, 8), "mm")) + patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 
    3, 2, 3, unit = "mm")))

grDevices::dev.off()
figure_export(figure_object,'Figure3',180,184)
