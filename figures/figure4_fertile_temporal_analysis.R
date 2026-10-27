#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="figure4_fertile_temporal_analysis.R" & file.exists(.paths)]
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
figure_preflight('Figure4',c("data.table","ggplot2","patchwork","gtable","ragg","systemfonts","scales","R.utils"))
if ('--check' %in% commandArgs(TRUE)) {message('Figure4: inputs, packages and font available');quit(status=0)}
source_paths <- setNames(vapply(figure_index$key[figure_index$figure=="Figure4"],function(k)figure_input("Figure4",k),character(1)),figure_index$key[figure_index$figure=="Figure4"])

plot_data <- capture_stage("locked_table_read_and_plotting_transform", {
    broad_summary <- data.table::fread(source_paths[["broad_type_summary"]])
    cluster_map <- data.table::fread(source_paths[["cluster_annotation_map"]])
    composition <- data.table::fread(source_paths[["sample_composition"]])
    program_summary <- data.table::fread(source_paths[["temporal_program_summary"]])
    program_models <- data.table::fread(source_paths[["temporal_program_models"]])
    VDR_summary <- data.table::fread(source_paths[["temporal_VDR_summary"]])
    VDR_models <- data.table::fread(source_paths[["temporal_VDR_models"]])
    omnibus_DE <- data.table::fread(source_paths[["temporal_omnibus_DE_summary"]])
    pairwise_DE <- data.table::fread(source_paths[["temporal_pairwise_DE_summary"]])
    cell_meta <- data.table::fread(source_paths[["cell_annotation_metadata"]], select = c("cell_barcode", "UMAPH_1", "UMAPH_2", 
        "broad_cell_type"))
    if (nrow(cell_meta) != 236636L || data.table::uniqueN(cell_meta$cell_barcode) != 236636L) {
        stop("Locked cell annotation metadata row/ID count mismatch.")
    }
    if (any(!is.finite(cell_meta$UMAPH_1)) || any(!is.finite(cell_meta$UMAPH_2))) {
        stop("Non-finite locked UMAP coordinates detected.")
    }
    full_type_order <- c("B cells", "Ciliated epithelial cells", "Endothelial cells", "Mast cells", "Myeloid cells", "NK/T cells", 
        "Stromal cells", "Unciliated epithelial cells")
    temporal_type_order <- c("Ciliated epithelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells")
    program_order <- c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")
    data.table::setorder(cell_meta, broad_cell_type, cell_barcode)
    umap_subset <- cell_meta[, {
        keep_n <- min(.N, 10000L)
        idx <- unique(as.integer(round(seq.int(1L, .N, length.out = keep_n))))
        .SD[idx]
    }, by = broad_cell_type]
    type_sizes <- setNames(broad_summary$cell_n, broad_summary$broad_cell_type)
    umap_subset[, `:=`(full_cell_n, as.integer(type_sizes[broad_cell_type]))]
    type_display_order <- broad_summary[order(-cell_n), broad_cell_type]
    umap_subset[, `:=`(display_order, match(broad_cell_type, type_display_order))]
    data.table::setorder(umap_subset, display_order, cell_barcode)
    umap_centroids <- cell_meta[, .(UMAPH_1 = stats::median(UMAPH_1), UMAPH_2 = stats::median(UMAPH_2), full_cell_n = .N), 
        by = broad_cell_type]
    umap_centroids[, `:=`(broad_cell_type, factor(broad_cell_type, levels = full_type_order))]
    composition[, `:=`(temporal_fertile_eligible, as.logical(temporal_fertile_eligible))]
    comp_fertile <- composition[temporal_fertile_eligible == TRUE]
    sample_meta <- unique(comp_fertile[, .(sample_accession, LH_day)])
    data.table::setorder(sample_meta, LH_day, sample_accession)
    sample_meta[, `:=`(sample_order, seq_len(.N))]
    sample_meta[, `:=`(donor_label, substr(sample_accession, nchar(sample_accession) - 3L, nchar(sample_accession)))]
    sample_meta[, `:=`(stage_label, factor(paste0("LH+", LH_day), levels = paste0("LH+", c(3, 5, 7, 9, 11))))]
    comp_fertile <- merge(comp_fertile, sample_meta, by = c("sample_accession", "LH_day"), all.x = TRUE)
    comp_fertile[, `:=`(broad_cell_type, factor(broad_cell_type, levels = full_type_order))]
    comp_fertile[, `:=`(donor_label, factor(donor_label, levels = sample_meta$donor_label))]
    program_plot <- merge(program_summary, program_models[, .(broad_cell_type, program_name, donor_n, FDR_BH)], by = c("broad_cell_type", 
        "program_name"), all.x = TRUE, sort = FALSE)
    program_labels <- c(RNA_ALL_436 = "RNA program\n(n=436)", OCCUPANCY_SUPPORTED_284 = "Occupancy-supported\n(n=284)", SHARED_OCCUPANCY_46 = "Shared occupancy\n(n=46)", 
        PROMOTER_OCCUPANCY_84 = "Promoter-associated\n(n=84)")
    cell_short <- c(`Ciliated epithelial cells` = "Ciliated epi.", `Myeloid cells` = "Myeloid", `NK/T cells` = "NK/T", `Stromal cells` = "Stromal", 
        `Unciliated epithelial cells` = "Unciliated epi.")
    program_plot[, `:=`(program_label, factor(unname(program_labels[program_name]), levels = unname(program_labels[program_order])))]
    program_plot[, `:=`(cell_label, factor(unname(cell_short[broad_cell_type]), levels = unname(cell_short[temporal_type_order])))]
    program_plot[, `:=`(stage_index, match(LH_day, c(3, 5, 7, 9, 11)))]
    program_sig <- unique(program_plot[, .(broad_cell_type, program_name, program_label, cell_label, FDR_BH, significance_label = ifelse(FDR_BH < 
        0.050000000000000003, "*", ""))])
    VDR_plot <- merge(VDR_summary, VDR_models[, .(broad_cell_type, donor_n, FDR_BH)], by = "broad_cell_type", all.x = TRUE, 
        sort = FALSE)
    VDR_plot[, `:=`(stage_label, factor(paste0("LH+", LH_day), levels = paste0("LH+", c(3, 5, 7, 9, 11))))]
    VDR_plot[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(temporal_type_order)))]
    VDR_plot[, `:=`(display_label, formatC(mean, format = "f", digits = 1))]
    VDR_plot[, `:=`(text_color, ifelse(mean >= 3.2999999999999998, "white", "#161616"))]
    pairwise_count <- pairwise_DE[, .(significant_row_n = sum(whole_genome_up_n + whole_genome_down_n)), by = broad_cell_type]
    omnibus_count <- omnibus_DE[, .(broad_cell_type, significant_row_n = whole_genome_FDR_lt_0_05_n)]
    omnibus_count[, `:=`(result_scope, "Omnibus\ngene \u00d7 cell type")]
    pairwise_count[, `:=`(result_scope, "Planned pairwise\ngene \u00d7 cell type \u00d7 contrast")]
    DE_count_plot <- data.table::rbindlist(list(omnibus_count, pairwise_count), use.names = TRUE)
    DE_count_plot[, `:=`(result_scope, factor(result_scope, levels = c("Omnibus\ngene \u00d7 cell type", "Planned pairwise\ngene \u00d7 cell type \u00d7 contrast")))]
    DE_count_plot[, `:=`(broad_cell_type, factor(broad_cell_type, levels = rev(temporal_type_order)))]
    DE_count_plot[, `:=`(count_label, format(significant_row_n, big.mark = ",", scientific = FALSE, trim = TRUE))]
    DE_count_plot[, `:=`(facet_max, max(significant_row_n)), by = result_scope]
    DE_count_plot[, `:=`(label_inside, significant_row_n >= 0.17999999999999999 * facet_max)]
    DE_count_plot[, `:=`(label_x, data.table::fifelse(label_inside, significant_row_n - 0.025000000000000001 * facet_max, 
        significant_row_n + 0.025000000000000001 * facet_max))]
    DE_count_plot[, `:=`(label_hjust, data.table::fifelse(label_inside, 1, 0))]
    DE_count_plot[, `:=`(label_color, data.table::fifelse(label_inside, "white", "#161616"))]
    list(broad_summary = broad_summary, full_cell_n = nrow(cell_meta), cluster_map = cluster_map, composition = composition, 
        sample_meta = sample_meta, umap_subset = umap_subset, umap_centroids = umap_centroids, comp_fertile = comp_fertile, 
        program_summary = program_summary, program_models = program_models, program_plot = program_plot, program_sig = program_sig, 
        VDR_summary = VDR_summary, VDR_models = VDR_models, VDR_plot = VDR_plot, omnibus_DE = omnibus_DE, pairwise_DE = pairwise_DE, 
        DE_count_plot = DE_count_plot, full_type_order = full_type_order, temporal_type_order = temporal_type_order, program_order = program_order)
})

composition_sum <- plot_data$comp_fertile[, .(fraction_sum = sum(within_sample_fraction)), by = sample_accession]

program_wide <- data.table::dcast(plot_data$program_summary, broad_cell_type + program_name ~ LH_day, value.var = "mean")

program_compare <- merge(plot_data$program_models, program_wide, by = c("broad_cell_type", "program_name"), all = FALSE)

program_mean_diff <- max(abs(c(program_compare$LH3_mean - program_compare[["3"]], program_compare$LH5_mean - program_compare[["5"]], 
    program_compare$LH7_mean - program_compare[["7"]], program_compare$LH9_mean - program_compare[["9"]], program_compare$LH11_mean - 
        program_compare[["11"]])), na.rm = TRUE)

VDR_wide <- data.table::dcast(plot_data$VDR_summary, broad_cell_type ~ LH_day, value.var = "mean")

VDR_compare <- merge(plot_data$VDR_models, VDR_wide, by = "broad_cell_type", all = FALSE)

VDR_mean_diff <- max(abs(c(VDR_compare$LH3_mean - VDR_compare[["3"]], VDR_compare$LH5_mean - VDR_compare[["5"]], VDR_compare$LH7_mean - 
    VDR_compare[["7"]], VDR_compare$LH9_mean - VDR_compare[["9"]], VDR_compare$LH11_mean - VDR_compare[["11"]])), na.rm = TRUE)

omnibus_total <- sum(plot_data$omnibus_DE$whole_genome_FDR_lt_0_05_n)

pairwise_total <- sum(plot_data$pairwise_DE$whole_genome_up_n + plot_data$pairwise_DE$whole_genome_down_n)

statistical_qc <- data.frame(check_id = sprintf("S%02d", 1:12), check = c("Locked embedding cell rows", "Locked broad cell types", 
    "All GSE250130 samples", "Fertile temporal donors used in Panel B", "Panel B within-donor fractions sum to one", "Temporal program models", 
    "Temporal program models with FDR < 0.05", "Program summary means match frozen model means", "Temporal VDR models with FDR < 0.05", 
    "VDR summary means match frozen model means", "Temporal omnibus significant gene-by-cell-type rows", "Planned-pairwise significant gene-by-cell-type-contrast rows"), 
    observed = c(plot_data$full_cell_n, nrow(plot_data$broad_summary), data.table::uniqueN(plot_data$composition$sample_accession), 
        nrow(plot_data$sample_meta), max(abs(composition_sum$fraction_sum - 1)), nrow(plot_data$program_models), sum(plot_data$program_models$FDR_BH < 
            0.050000000000000003), program_mean_diff, sum(plot_data$VDR_models$FDR_BH < 0.050000000000000003), VDR_mean_diff, 
        omnibus_total, pairwise_total), expected = c(236636, 8, 28, 18, 0, 20, 13, 0, 0, 0, 24168, 55539), tolerance = c(0, 
        0, 0, 0, 1e-10, 0, 0, 9.9999999999999998e-13, 0, 9.9999999999999998e-13, 0, 0), stringsAsFactors = FALSE)

statistical_qc$pass <- abs(as.numeric(statistical_qc$observed) - statistical_qc$expected) <= statistical_qc$tolerance

if (!all(statistical_qc$pass)) stop("Figure 4 statistical/source consistency QC failed.")

cell_colors <- c(`B cells` = "#D55E00", `Ciliated epithelial cells` = "#E69F00", `Endothelial cells` = "#009E73", `Mast cells` = "#A6761D", 
    `Myeloid cells` = "#00A6D6", `NK/T cells` = "#0072B2", `Stromal cells` = "#7A5AB5", `Unciliated epithelial cells` = "#CC79A7")

header_tag_size_pt <- 11

header_title_size_pt <- 9.1999999999999993

header_subtitle_size_pt <- 7.2000000000000002

header_tag_to_title_gap_mm <- 1.6000000000000001

header_title_row_height_mm <- 5.4000000000000004

header_subtitle_row_height_mm <- 4

figure_width_mm <- 180

figure_height_mm <- 250

color_text <- "#161616"

color_muted <- "#535353"

color_grid <- "#D9D9D9"

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2000000000000002) + ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.23999999999999999, 
    colour = color_grid), panel.grid.minor = ggplot2::element_blank(), axis.title = ggplot2::element_text(family = figure_font, 
    size = 7.5, colour = color_text), axis.text = ggplot2::element_text(family = figure_font, size = 6.2000000000000002, 
    colour = color_text), legend.title = ggplot2::element_text(family = figure_font, size = 6.5, face = "bold", colour = color_text), 
    legend.text = ggplot2::element_text(family = figure_font, size = 6.2000000000000002, colour = color_text), strip.text = ggplot2::element_text(family = figure_font, 
        size = 6.2000000000000002, face = "bold", colour = color_text), plot.margin = ggplot2::margin(2, 2.5, 2, 2.5, unit = "mm"))

panel_A <- capture_stage("panel_A_construction", {
    ggplot2::ggplot(plot_data$umap_subset, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, colour = broad_cell_type)) + ggplot2::geom_point(size = 0.074999999999999997, 
        alpha = 0.57999999999999996, stroke = 0) + ggplot2::geom_label(data = plot_data$umap_centroids, ggplot2::aes(x = UMAPH_1, 
        y = UMAPH_2, label = broad_cell_type), inherit.aes = FALSE, family = figure_font, fontface = "bold", size = 2.1499999999999999, 
        linewidth = 0, label.padding = grid::unit(0.69999999999999996, "mm"), fill = scales::alpha("white", 0.78000000000000003), 
        colour = color_text) + ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) + ggplot2::coord_fixed(clip = "off") + 
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2") + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), 
        axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 6.2999999999999998), legend.position = "none", 
        plot.margin = ggplot2::margin(1.5, 1.5, 1.5, 1.5, unit = "mm"))
})

panel_B <- capture_stage("panel_B_construction", {
    ggplot2::ggplot(plot_data$comp_fertile, ggplot2::aes(x = donor_label, y = within_sample_fraction, fill = broad_cell_type)) + 
        ggplot2::geom_col(width = 0.92000000000000004) + ggplot2::facet_grid(~stage_label, scales = "free_x", space = "free_x") + 
        ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE, name = "Broad cell type") + ggplot2::scale_y_continuous(expand = c(0, 
        0), labels = scales::label_percent(accuracy = 25)) + ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) + ggplot2::labs(x = "Donor (GSM suffix)", 
        y = "Within-donor fraction") + base_theme + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5.4000000000000004, 
        angle = 90, hjust = 1, vjust = 0.5), axis.text.y = ggplot2::element_text(size = 5.7999999999999998), strip.background = ggplot2::element_rect(fill = "#F2F2F2", 
        colour = NA), strip.text.x = ggplot2::element_text(size = 6.2000000000000002, face = "bold"), panel.spacing.x = grid::unit(0.59999999999999998, 
        "mm"), legend.position = "bottom", legend.box = "vertical", legend.margin = ggplot2::margin(t = 0.5, b = 0), legend.spacing.y = grid::unit(0.20000000000000001, 
        "mm"), legend.key.height = grid::unit(2.5, "mm"), legend.key.width = grid::unit(3.2000000000000002, "mm"), plot.margin = ggplot2::margin(1.5, 
        2, 1, 1.5, unit = "mm")) + ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE, title.position = "top"))
})

panel_C <- capture_stage("panel_C_construction", {
    ggplot2::ggplot(plot_data$program_plot, ggplot2::aes(x = stage_index, y = mean, group = 1, colour = broad_cell_type, 
        fill = broad_cell_type)) + ggplot2::geom_hline(yintercept = 0, linewidth = 0.22, linetype = "dashed", colour = "#8A8A8A") + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high), alpha = 0.14999999999999999, colour = NA) + ggplot2::geom_line(linewidth = 0.47999999999999998) + 
        ggplot2::geom_point(size = 0.90000000000000002) + ggplot2::geom_text(data = plot_data$program_sig, ggplot2::aes(x = 5, 
        y = Inf, label = significance_label), inherit.aes = FALSE, family = figure_font, fontface = "bold", size = 2.7000000000000002, 
        hjust = 1.1000000000000001, vjust = 1, colour = color_text) + ggplot2::facet_grid(program_label ~ cell_label, scales = "free_y", 
        switch = "y") + ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) + ggplot2::scale_fill_manual(values = cell_colors, 
        drop = FALSE) + ggplot2::scale_x_continuous(breaks = 1:5, labels = paste0("LH+", c(3, 5, 7, 9, 11)), expand = ggplot2::expansion(mult = c(0.080000000000000002, 
        0.080000000000000002))) + ggplot2::labs(x = "Biopsy timing", y = "Signed VDR-program score") + base_theme + ggplot2::theme(legend.position = "none", 
        axis.text.x = ggplot2::element_text(size = 5.2000000000000002, angle = 35, hjust = 1), axis.text.y = ggplot2::element_text(size = 5.2000000000000002), 
        axis.title = ggplot2::element_text(size = 6.5999999999999996), strip.background = ggplot2::element_rect(fill = "#F2F2F2", 
            colour = NA), strip.placement = "outside", strip.text.x = ggplot2::element_text(size = 5.7999999999999998, lineheight = 0.92000000000000004), 
        strip.text.y.left = ggplot2::element_text(size = 5.5, angle = 0, hjust = 1, lineheight = 0.90000000000000002), panel.spacing = grid::unit(0.65000000000000002, 
            "mm"), plot.margin = ggplot2::margin(1, 1.5, 1.2, 1.5, unit = "mm"))
})

panel_D <- capture_stage("panel_D_construction", {
    ggplot2::ggplot(plot_data$VDR_plot, ggplot2::aes(x = stage_label, y = broad_cell_type, fill = mean)) + ggplot2::geom_tile(linewidth = 0.45000000000000001, 
        colour = "white") + ggplot2::geom_text(ggplot2::aes(label = display_label, colour = text_color), family = figure_font, 
        size = 2.0499999999999998) + ggplot2::scale_colour_identity() + ggplot2::scale_fill_gradient(low = "#F4F3F0", high = "#405466", 
        name = "Mean VDR\nlog2(CPM + 1)", limits = c(0, 6), breaks = c(0, 2, 4, 6), oob = scales::squish) + ggplot2::labs(x = "Biopsy timing", 
        y = NULL) + base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5.9000000000000004, 
        angle = 35, hjust = 1), axis.text.y = ggplot2::element_text(size = 6, lineheight = 0.92000000000000004), legend.position = "bottom", 
        legend.direction = "horizontal", legend.title = ggplot2::element_text(size = 6.2000000000000002, face = "bold", margin = ggplot2::margin(r = 2, 
            unit = "mm")), legend.text = ggplot2::element_text(size = 5.7999999999999998), legend.key.height = grid::unit(2.6000000000000001, 
            "mm"), legend.margin = ggplot2::margin(t = 0.5, b = 0), plot.margin = ggplot2::margin(1.5, 2, 1, 1.5, unit = "mm")) + 
        ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "left", title.hjust = 0.5, barwidth = grid::unit(34, 
            "mm"), barheight = grid::unit(2.6000000000000001, "mm")))
})

panel_E <- capture_stage("panel_E_construction", {
    ggplot2::ggplot(plot_data$DE_count_plot, ggplot2::aes(x = significant_row_n, y = broad_cell_type, fill = broad_cell_type)) + 
        ggplot2::geom_col(width = 0.68000000000000005) + ggplot2::geom_text(ggplot2::aes(x = label_x, label = count_label, 
        hjust = label_hjust, colour = label_color), family = figure_font, size = 2, show.legend = FALSE) + ggplot2::facet_grid(~result_scope, 
        scales = "free_x", space = "free_x") + ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE) + ggplot2::scale_colour_identity() + 
        ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 3), labels = scales::label_number(scale_cut = scales::cut_short_scale()), 
            expand = ggplot2::expansion(mult = c(0, 0.080000000000000002))) + ggplot2::labs(x = "FDR-significant rows", y = NULL) + 
        base_theme + ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5.5), 
        axis.text.y = ggplot2::element_text(size = 5.9000000000000004, lineheight = 0.92000000000000004), axis.title.x = ggplot2::element_text(size = 6.5), 
        strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA), strip.text.x = ggplot2::element_text(size = 5.7999999999999998, 
            lineheight = 0.92000000000000004), panel.spacing.x = grid::unit(2, "mm"), plot.margin = ggplot2::margin(1.5, 
            2, 1, 1.5, unit = "mm"))
})

build_body_gtable <- function(plot_object) {
    ggplot2::ggplotGrob(plot_object + ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) + ggplot2::theme(plot.tag = ggplot2::element_blank(), 
        plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank()))
}

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle, horizontal_nudge_mm = 0) {
    panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
    if (!nrow(panel_cells)) 
        stop("Could not locate scientific panel cell for Panel ", panel_tag)
    panel_left <- min(panel_cells$l)
    panel_right <- max(panel_cells$r)
    header_gt <- gtable::gtable(widths = body_gtable$widths, heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), 
        "mm"))
    title_grob <- grid::textGrob(panel_title, x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm"), y = grid::unit(0.52000000000000002, 
        "npc"), hjust = 0, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_title_size_pt, 
        lineheight = 0.93999999999999995, col = color_text))
    tag_grob <- grid::textGrob(panel_tag, x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm") - grid::unit(header_tag_to_title_gap_mm, 
        "mm"), y = grid::unit(0.52000000000000002, "npc"), hjust = 1, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, 
        fontface = "bold", fontsize = header_tag_size_pt, col = color_text))
    subtitle_grob <- grid::textGrob(panel_subtitle, x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm"), y = grid::unit(0.56000000000000005, 
        "npc"), hjust = 0, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, lineheight = 0.93999999999999995, 
        col = color_muted))
    header_gt <- gtable::gtable_add_grob(header_gt, grid::grobTree(tag_grob, title_grob), t = 1, b = 1, l = panel_left, r = panel_right, 
        clip = "off", name = paste0("F4-header-title-", panel_tag))
    header_gt <- gtable::gtable_add_grob(header_gt, subtitle_grob, t = 2, b = 2, l = panel_left, r = panel_right, clip = "off", 
        name = paste0("F4-header-subtitle-", panel_tag))
    list(gtable = header_gt, geometry = data.frame(panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE", panel_left_column = panel_left, 
        panel_right_column = panel_right, tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt, subtitle_size_pt = header_subtitle_size_pt, 
        tag_title_gap_mm = header_tag_to_title_gap_mm, horizontal_nudge_mm = horizontal_nudge_mm, header_composition = "DEDICATED_SHARED_ROW", 
        stringsAsFactors = FALSE))
}

metrics_device_path <- tempfile(pattern = "FIG4_ragg_metrics_", fileext = ".png")

ragg::agg_png(filename = metrics_device_path, width = figure_width_mm, height = figure_height_mm, units = "mm", res = 300, 
    background = "white")

composite_parts <- capture_stage("gtable_and_composite_assembly", {
    body_A <- build_body_gtable(panel_A)
    body_B <- build_body_gtable(panel_B)
    body_C <- build_body_gtable(panel_C)
    body_D <- build_body_gtable(panel_D)
    body_E <- build_body_gtable(panel_E)
    header_A <- make_header_gtable(body_A, "A", "Broad cell-type landscape", "8 broad types \u00b7 236,636 cells (display subset)", 
        horizontal_nudge_mm = 14)
    header_B <- make_header_gtable(body_B, "B", "Fertile donor composition", "18 fertile donors \u00b7 descriptive only")
    header_C <- make_header_gtable(body_C, "C", "Stage-associated VDR-program remodeling", "Cross-sectional donor means \u00b1 95% CI \u00b7 * omnibus FDR < 0.05")
    header_D <- make_header_gtable(body_D, "D", "VDR expression across stages", "Mean log2(CPM + 1) \u00b7 0/5 omnibus FDR < 0.05")
    header_E <- make_header_gtable(body_E, "E", "Transcriptome-wide temporal DE", "Genome-wide FDR < 0.05 \u00b7 independent axes")
    design <- "\n##\nAB\nCD\n##\nEE\nFF\n##\nGH\nIJ\n##\n"
    figure_object <- (patchwork::wrap_elements(full = header_A$gtable) + patchwork::wrap_elements(full = header_B$gtable) + 
        patchwork::wrap_elements(full = body_A) + patchwork::wrap_elements(full = body_B) + patchwork::wrap_elements(full = header_C$gtable) + 
        patchwork::wrap_elements(full = body_C) + patchwork::wrap_elements(full = header_D$gtable) + patchwork::wrap_elements(full = header_E$gtable) + 
        patchwork::wrap_elements(full = body_D) + patchwork::wrap_elements(full = body_E)) + patchwork::plot_layout(design = design, 
        widths = c(0.92000000000000004, 1.0800000000000001), heights = grid::unit(c(1, 10, 52, 3, 10, 57, 3, 10, 74, 18), 
            "mm")) + patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm")))
    list(figure = figure_object, header_geometry = do.call(rbind, lapply(list(header_A, header_B, header_C, header_D, header_E), 
        `[[`, "geometry")))
})

grDevices::dev.off()
figure_export(composite_parts$figure,'Figure4',180,250)
