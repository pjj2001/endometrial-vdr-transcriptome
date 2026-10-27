#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="figure2_perturbation_occupancy.R" & file.exists(.paths)]
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
figure_preflight('Figure2',c("data.table","ggplot2","patchwork","gtable","ragg","systemfonts","scales","R.utils","ggrepel"))
if ('--check' %in% commandArgs(TRUE)) {message('Figure2: inputs, packages and font available');quit(status=0)}
source_paths <- setNames(vapply(figure_index$key[figure_index$figure=="Figure2"],function(k)figure_input("Figure2",k),character(1)),figure_index$key[figure_index$figure=="Figure2"])
.metric_file <- tempfile(fileext='.png')
ragg::agg_png(.metric_file,width=180,height=170,units='mm',res=96,background='white')

target_width_mm <- 180

target_height_mm <- 170

stress_test_width_mm <- 170

top_block_tag_x <- 0

top_block_tag_y <- 1

top_block_tag_to_title_gap_pt <- 0

top_block_anchor_reference <- "NATIVE_TOP_BLOCK_STRIPPED"

header_anchor_reference <- "GTABLE_PANEL_CELL_LEFT_EDGE"

header_tag_size_pt <- 9

header_title_size_pt <- 8.5

header_subtitle_size_pt <- 7.2000000000000002

header_tag_to_title_gap_mm <- 2

header_title_row_height_mm <- 5

header_subtitle_row_height_mm <- 4

header_title_lineheight <- 0.93999999999999995

panel_B_axis_text_y_size <- 5.4000000000000004

panel_C_x_limit <- 310

panel_D_x_limit <- 205

panel_A_y_tick_to_axis_gap_pt <- 0.75

panel_A_y_title_to_tick_gap_pt <- 0.69999999999999996

panel_A_horizontal_geometry_state <- "LABEL_WIDTH_DECOUPLED_LOCK"

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    output <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), error = function(e) NULL)
    if (is.null(output)) {
        output <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    output
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

as_logical_safe <- function(x) {
    if (is.logical(x)) 
        return(x)
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

first_existing <- function(x, candidates, label) {
    hit <- candidates[candidates %in% names(x)]
    if (!length(hit)) {
        stop("Could not identify ", label, ". Expected one of: ", paste(candidates, collapse = ", "))
    }
    hit[[1]]
}

require_columns <- function(x, columns, label) {
    missing_columns <- setdiff(columns, names(x))
    if (length(missing_columns)) {
        stop(label, " is missing column(s): ", paste(missing_columns, collapse = ", "))
    }
    invisible(TRUE)
}

volcano <- safe_read_csv(source_paths[["volcano"]])

volcano_labels <- safe_read_csv(source_paths[["volcano_labels"]])

heatmap <- safe_read_csv(source_paths[["heatmap"]])

sample_manifest <- safe_read_csv(source_paths[["sample_manifest"]])

RNA_frozen <- safe_read_csv(source_paths[["RNA_frozen"]])

integrated <- safe_read_csv(source_paths[["integrated_436"]])

routing <- safe_read_csv(source_paths[["routing_436"]])

require_columns(volcano, c("feature_id", "interaction_logFC", "interaction_FDR"), "Volcano source")

require_columns(heatmap, c("feature_id", "short_label", "z_score"), "Heatmap source")

require_columns(integrated, c("feature_id", "program_class", "interaction_logFC", "interaction_FDR", "occupancy_supported", 
    "shared_occupancy_supported", "any_promoter_overlap_support", "vehicle_peak_n", "calcitriol_peak_n", "shared_peak_n", 
    "occupancy_evidence_group"), "RNA/occupancy evidence table")

require_columns(routing, c("feature_id", "downstream_route_class", "direct_target_status"), "Downstream routing table")

if (nrow(volcano) != 12116L || nrow(integrated) != 436L || nrow(routing) != 436L || length(unique(heatmap$feature_id)) != 
    30L) {
    stop("Figure 2 source row counts do not match the locked design.")
}

integrated$occupancy_supported <- as_logical_safe(integrated$occupancy_supported)

integrated$shared_occupancy_supported <- as_logical_safe(integrated$shared_occupancy_supported)

integrated$any_promoter_overlap_support <- as_logical_safe(integrated$any_promoter_overlap_support)

if (sum(integrated$occupancy_supported) != 284L || sum(!integrated$occupancy_supported) != 152L || sum(integrated$shared_occupancy_supported) != 
    46L || sum(integrated$any_promoter_overlap_support) != 84L || sum(safe_numeric(integrated$vehicle_peak_n) > 0) != 77L || 
    sum(safe_numeric(integrated$calcitriol_peak_n) > 0) != 275L) {
    stop("RNA/occupancy evidence counts do not match the locked values.")
}

if (any(routing$direct_target_status != "NOT_ASSIGNED")) {
    stop("One or more routing rows incorrectly assign direct-target status.")
}

volcano$interaction_logFC <- safe_numeric(volcano$interaction_logFC)

volcano$interaction_FDR <- safe_numeric(volcano$interaction_FDR)

if (!"minus_log10_FDR" %in% names(volcano)) {
    volcano$minus_log10_FDR <- -log10(pmax(volcano$interaction_FDR, .Machine$double.xmin))
}

if (!"is_interaction_candidate" %in% names(volcano)) {
    volcano$is_interaction_candidate <- (volcano$interaction_FDR < 0.050000000000000003 & abs(volcano$interaction_logFC) >= 
        0.5)
} else {
    volcano$is_interaction_candidate <- as_logical_safe(volcano$is_interaction_candidate)
}

if (!"plot_class" %in% names(volcano)) {
    program_match <- integrated$program_class[match(volcano$feature_id, integrated$feature_id)]
    volcano$plot_class <- "Not candidate"
    volcano$plot_class[program_match == "attenuated_ligand_induction_after_VDRKD"] <- "Attenuated induction"
    volcano$plot_class[program_match == "attenuated_ligand_repression_after_VDRKD"] <- "Attenuated repression"
    volcano$plot_class[volcano$is_interaction_candidate & volcano$plot_class == "Not candidate"] <- "Other interaction"
}

volcano$plot_class <- factor(as.character(volcano$plot_class), levels = c("Not candidate", "Attenuated induction", "Attenuated repression", 
    "Other interaction"))

label_gene_column <- first_existing(volcano_labels, c("feature_id", "gene_symbol", "label"), "volcano label identifier")

label_ids <- unique(as.character(volcano_labels[[label_gene_column]]))

label_data <- volcano[volcano$feature_id %in% label_ids, , drop = FALSE]

label_data$display_label <- as.character(label_data$feature_id)

manual_label_layout <- data.frame(feature_id = c("CYP24A1", "COL14A1", "RASSF2", "EFL1", "VEGFA", "IL11", "PPP1R3C", "SLC16A6", 
    "NFKBIZ", "F2RL2", "RGMB", "FOS"), label_x = c(-5.1200000000000001, -2.8199999999999998, -2.52, -1.55, -0.92000000000000004, 
    -1.52, 1.48, 1.5600000000000001, 1.55, 1.55, 1.24, 1.1599999999999999), label_y = c(13.539999999999999, 8.4199999999999999, 
    7.1500000000000004, 8.8200000000000003, 8.1199999999999992, 4.4500000000000002, 7.5800000000000001, 6.9000000000000004, 
    6.25, 5.6200000000000001, 4.9500000000000002, 3.6000000000000001), hjust = c(0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0), use_connector = c(FALSE, 
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), manual_priority = TRUE, stringsAsFactors = FALSE)

label_data <- merge(label_data, manual_label_layout, by = "feature_id", all.x = TRUE, sort = FALSE)

if (any(is.na(label_data$label_x)) || any(is.na(label_data$label_y)) || any(is.na(label_data$hjust))) {
    stop("The manual volcano-label layout is incomplete.")
}

label_data$connector_start_x <- ifelse(label_data$hjust == 0, label_data$label_x - 0.055, label_data$label_x + 0.055)

label_data$connector_start_y <- label_data$label_y

label_data$connector_end_x <- label_data$interaction_logFC

label_data$connector_end_y <- label_data$minus_log10_FDR

label_data$connector_crosses_zero_side <- (sign(label_data$label_x) != sign(label_data$interaction_logFC))

if (any(label_data$connector_crosses_zero_side)) {
    stop("A volcano label was placed on the opposite interaction side ", "from its source point.")
}

connector_data <- label_data[label_data$use_connector, , drop = FALSE]

direct_label_data <- label_data[!label_data$use_connector, , drop = FALSE]

sample_order <- if ("short_label" %in% names(sample_manifest)) {
    as.character(sample_manifest$short_label)
} else {
    unique(as.character(heatmap$short_label))
}

sample_order <- sample_order[sample_order %in% unique(as.character(heatmap$short_label))]

if (length(sample_order) != length(unique(heatmap$short_label))) {
    sample_order <- unique(as.character(heatmap$short_label))
}

feature_order <- unique(as.character(heatmap$feature_id))

heatmap$feature_id <- factor(heatmap$feature_id, levels = rev(feature_order))

heatmap$short_label <- factor(heatmap$short_label, levels = sample_order)

heatmap$z_score <- safe_numeric(heatmap$z_score)

heatmap_limit <- max(abs(heatmap$z_score), na.rm = TRUE)

support_counts <- data.frame(category = c("Occupancy-supported", "RNA only", "Calcitriol occupancy", "Vehicle occupancy", 
    "Shared occupancy", "Promoter-associated"), count = c(284L, 152L, 275L, 77L, 46L, 84L), group = c("Partition", "Partition", 
    "Overlapping attributes", "Overlapping attributes", "Overlapping attributes", "Overlapping attributes"), color_key = c("Occupancy-supported", 
    "RNA only", "Calcitriol occupancy", "Vehicle occupancy", "Shared occupancy", "Promoter-associated"), stringsAsFactors = FALSE)

support_counts$category <- factor(support_counts$category, levels = rev(support_counts$category))

support_counts$value_label <- paste0(support_counts$count, " (", round(100 * support_counts$count/436, 1), "%)")

support_counts$label_inside <- (support_counts$count >= 250L)

route_counts <- as.data.frame(table(routing$downstream_route_class), stringsAsFactors = FALSE)

names(route_counts) <- c("route_class", "candidate_n")

route_label_map <- c(ROUTE_1_SHARED_AND_PROMOTER_SUPPORT = "Shared + promoter", ROUTE_2_SHARED_NONPROMOTER_SUPPORT = "Shared, non-promoter", 
    ROUTE_3_NONSHARED_PROMOTER_SUPPORT = "Non-shared promoter", ROUTE_4_OTHER_OCCUPANCY_SUPPORT = "Other occupancy support", 
    ROUTE_5_RNA_ONLY = "RNA only")

route_order <- names(route_label_map)

if (!setequal(route_counts$route_class, route_order)) {
    stop("Downstream route classes do not match the locked five-class design.")
}

route_counts <- route_counts[match(route_order, route_counts$route_class), , drop = FALSE]

route_counts$route_label <- unname(route_label_map[route_counts$route_class])

route_counts$route_label <- factor(route_counts$route_label, levels = rev(unname(route_label_map)))

if (!identical(as.integer(route_counts$candidate_n), c(21L, 25L, 63L, 175L, 152L))) {
    stop("Downstream route counts do not match 21/25/63/175/152.")
}

strip_native_top_block <- function(plot_object) {
    plot_object + ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) + ggplot2::theme(plot.tag = ggplot2::element_blank(), 
        plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank())
}

build_body_gtable <- function(plot_object, free_left_labels = FALSE) {
    stripped_plot <- strip_native_top_block(plot_object)
    body_gtable <- ggplot2::ggplotGrob(stripped_plot)
    attr(body_gtable, "independent_left_label_width") <- TRUE
    attr(body_gtable, "legacy_free_left_requested") <- isTRUE(free_left_labels)
    body_gtable
}

add_explicit_panel_header <- function(plot_object, panel_tag, panel_title, panel_subtitle, free_left_labels = FALSE) {
    body_gtable <- build_body_gtable(plot_object, free_left_labels)
    panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
    if (!nrow(panel_cells)) 
        stop("Could not locate panel cell for Panel ", panel_tag)
    panel_left_column <- min(panel_cells$l)
    panel_right_column <- max(panel_cells$r)
    body_gtable <- gtable::gtable_add_rows(body_gtable, heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), 
        "mm"), pos = 0)
    title_grob <- grid::textGrob(panel_title, x = grid::unit(0, "npc"), y = grid::unit(0.52000000000000002, "npc"), hjust = 0, 
        vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_title_size_pt, lineheight = header_title_lineheight, 
            col = "#111111"))
    tag_grob <- grid::textGrob(panel_tag, x = grid::unit(0, "npc") - grid::unit(header_tag_to_title_gap_mm, "mm"), y = grid::unit(0.52000000000000002, 
        "npc"), hjust = 1, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_tag_size_pt, 
        col = "#111111"))
    subtitle_grob <- grid::textGrob(panel_subtitle, x = grid::unit(0, "npc"), y = grid::unit(0.56000000000000005, "npc"), 
        hjust = 0, vjust = 0.5, gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, col = "#444444"))
    body_gtable <- gtable::gtable_add_grob(body_gtable, grid::grobTree(tag_grob, title_grob), t = 1, b = 1, l = panel_left_column, 
        r = panel_right_column, clip = "off", name = paste0("explicit-title-row-", panel_tag))
    body_gtable <- gtable::gtable_add_grob(body_gtable, subtitle_grob, t = 2, b = 2, l = panel_left_column, r = panel_right_column, 
        clip = "off", name = paste0("explicit-subtitle-row-", panel_tag))
    list(gtable = body_gtable, geometry = data.frame(panel = panel_tag, title_left_reference = header_anchor_reference, detected_panel_left_column = panel_left_column, 
        detected_panel_right_column = panel_right_column, title_x_npc = 0, tag_right_offset_mm = -header_tag_to_title_gap_mm, 
        tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt, subtitle_size_pt = header_subtitle_size_pt, 
        native_title_tag_stripped = TRUE, left_label_width_strategy = "INDEPENDENT_GGPLOT_GTABLE", patchwork_free_used = FALSE, 
        patchworkGrob_used = FALSE, legacy_free_left_requested = isTRUE(free_left_labels), stringsAsFactors = FALSE))
}

palette <- c(`Not candidate` = "#BDBDBD", `Attenuated induction` = "#7B3294", `Attenuated repression` = "#008837", `Other interaction` = "#E6AB02", 
    `Occupancy-supported` = "#0072B2", `RNA only` = "#A6A6A6", `Calcitriol occupancy` = "#009E73", `Vehicle occupancy` = "#56B4E9", 
    `Shared occupancy` = "#CC79A7", `Promoter-associated` = "#E69F00")

route_palette <- c(`Shared + promoter` = "#004C6D", `Shared, non-promoter` = "#1677A6", `Non-shared promoter` = "#E69F00", 
    `Other occupancy support` = "#56B4E9", `RNA only` = "#A6A6A6")

theme_pub <- ggplot2::theme_classic(base_size = 8.4000000000000004, base_family = figure_font) + ggplot2::theme(plot.title.position = "panel", 
    plot.title = ggplot2::element_text(size = header_title_size_pt, face = "bold", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, 
        b = 3)), plot.subtitle = ggplot2::element_text(size = 7.2000000000000002, color = "#444444", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, 
        b = 4)), plot.tag = ggplot2::element_text(size = header_tag_size_pt, face = "bold", family = figure_font, hjust = 0, 
        vjust = 1), plot.tag.position = c(top_block_tag_x, top_block_tag_y), plot.tag.location = "panel", axis.title = ggplot2::element_text(size = 8), 
    axis.text = ggplot2::element_text(size = 7, color = "#222222"), legend.title = ggplot2::element_text(size = 7.2000000000000002, 
        face = "bold"), legend.text = ggplot2::element_text(size = 6.7999999999999998), plot.margin = ggplot2::margin(t = 4, 
        r = 4, b = 4, l = 2))

panel_A <- ggplot2::ggplot(volcano, ggplot2::aes(x = interaction_logFC, y = minus_log10_FDR, color = plot_class)) + ggplot2::geom_point(size = 0.75, 
    alpha = 0.54000000000000004, stroke = 0) + ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.28000000000000003, 
    color = "#555555") + ggplot2::geom_hline(yintercept = -log10(0.050000000000000003), linetype = 2, linewidth = 0.28000000000000003, 
    color = "#555555") + ggplot2::geom_segment(data = connector_data, ggplot2::aes(x = connector_start_x, y = connector_start_y, 
    xend = connector_end_x, yend = connector_end_y), inherit.aes = FALSE, linewidth = 0.20000000000000001, color = "#888888", 
    lineend = "round", show.legend = FALSE) + ggplot2::geom_point(data = label_data, ggplot2::aes(x = interaction_logFC, 
    y = minus_log10_FDR), inherit.aes = FALSE, shape = 21, size = 1.45, stroke = 0.29999999999999999, color = "#333333", 
    fill = NA, show.legend = FALSE) + ggplot2::geom_label(data = label_data, ggplot2::aes(x = label_x, y = label_y, label = display_label, 
    hjust = hjust), inherit.aes = FALSE, size = 2.2200000000000002, family = figure_font, color = "#222222", fill = scales::alpha("white", 
    0.81999999999999995), linewidth = 0, label.padding = grid::unit(0.029999999999999999, "lines"), label.r = grid::unit(0.014999999999999999, 
    "lines"), show.legend = FALSE) + ggplot2::scale_color_manual(values = palette[c("Not candidate", "Attenuated induction", 
    "Attenuated repression", "Other interaction")], labels = c("Not candidate", "Induction attenuated", "Repression attenuated", 
    "Other"), drop = FALSE) + ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top", 
    title.hjust = 0)) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.10000000000000001))) + ggplot2::coord_cartesian(xlim = c(-5.7800000000000002, 
    3.2000000000000002), clip = "off") + ggplot2::labs(tag = "A", title = "VDR RNA interaction program", subtitle = "436 candidates \u00b7 FDR < 0.05 \u00b7 |log2FC| \u2265 0.5", 
    x = "siRNA \u00d7 ligand interaction log2FC", y = expression(-log[10]("FDR")), color = "RNA program class") + theme_pub + 
    ggplot2::theme(plot.title = ggplot2::element_text(size = header_title_size_pt, face = "bold", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, 
        b = 1.2)), plot.subtitle = ggplot2::element_text(size = 7.2000000000000002, color = "#444444", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, 
        b = 1)), legend.position = "bottom", legend.justification = "left", legend.direction = "horizontal", legend.box = "vertical", 
        legend.key.width = grid::unit(0.75, "lines"), legend.key.height = grid::unit(0.5, "lines"), legend.spacing.x = grid::unit(0.20000000000000001, 
            "lines"), legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.050000000000000003, 
            "lines"), plot.tag.position = c(top_block_tag_x, top_block_tag_y), plot.tag.location = "panel", axis.title.y = ggplot2::element_text(size = 8, 
            angle = 90, vjust = 0.5, margin = ggplot2::margin(r = panel_A_y_title_to_tick_gap_pt)), axis.text.y = ggplot2::element_text(size = 7, 
            margin = ggplot2::margin(r = panel_A_y_tick_to_axis_gap_pt)), axis.title.x = ggplot2::element_text(size = 8, 
            margin = ggplot2::margin(t = 0.80000000000000004)), panel.border = ggplot2::element_rect(fill = NA, color = "#333333", 
            linewidth = 0.34999999999999998), plot.margin = ggplot2::margin(t = 4, r = 4, b = 0, l = 1))

separator_positions <- if (length(sample_order) == 12L) {
    c(3.5, 6.5, 9.5)
} else {
    numeric(0)
}

panel_B <- ggplot2::ggplot(heatmap, ggplot2::aes(x = short_label, y = feature_id, fill = z_score)) + ggplot2::geom_tile() + 
    ggplot2::geom_vline(xintercept = separator_positions, linewidth = 0.32000000000000001, color = "#222222") + ggplot2::scale_fill_gradient2(low = "#3B4CC0", 
    mid = "white", high = "#B40426", midpoint = 0, limits = c(-heatmap_limit, heatmap_limit), oob = scales::squish) + ggplot2::labs(tag = "B", 
    title = "Top interaction candidates\nacross 12 samples", subtitle = "Top 30 candidates from the perturbation analysis", 
    x = NULL, y = NULL, fill = "Row z-score") + theme_pub + ggplot2::theme(plot.title = ggplot2::element_text(size = header_title_size_pt, 
    face = "bold", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 2)), plot.subtitle = ggplot2::element_text(size = 7.2000000000000002, 
    color = "#444444", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 0.80000000000000004)), 
    plot.tag.position = c(top_block_tag_x, top_block_tag_y), plot.tag.location = "panel", axis.text.x = ggplot2::element_text(angle = 45, 
        hjust = 1, vjust = 1, size = 6.5999999999999996), axis.text.y = ggplot2::element_text(size = panel_B_axis_text_y_size, 
        margin = ggplot2::margin(r = 1)), axis.ticks = ggplot2::element_blank(), legend.position = "right", panel.border = ggplot2::element_rect(fill = NA, 
        color = "#333333", linewidth = 0.34999999999999998), plot.margin = ggplot2::margin(t = 3, r = 4, b = 4, l = 2))

panel_C <- ggplot2::ggplot(support_counts, ggplot2::aes(x = count, y = category, fill = color_key)) + ggplot2::geom_col(width = 0.68000000000000005, 
    color = "#333333", linewidth = 0.25) + ggplot2::geom_text(data = support_counts[support_counts$label_inside, , drop = FALSE], 
    ggplot2::aes(label = value_label), hjust = 1.0800000000000001, size = 2.3799999999999999, family = figure_font, fontface = "bold", 
    color = "white", show.legend = FALSE) + ggplot2::geom_text(data = support_counts[!support_counts$label_inside, , drop = FALSE], 
    ggplot2::aes(label = value_label), hjust = -0.10000000000000001, size = 2.3799999999999999, family = figure_font, color = "#222222", 
    show.legend = FALSE) + ggplot2::geom_hline(yintercept = 4.5, linetype = 3, linewidth = 0.32000000000000001, color = "#777777") + 
    ggplot2::scale_fill_manual(values = palette[c("Occupancy-supported", "RNA only", "Calcitriol occupancy", "Vehicle occupancy", 
        "Shared occupancy", "Promoter-associated")], guide = "none") + ggplot2::scale_y_discrete(labels = c(`Occupancy-supported` = "Occupancy-\nsupported", 
    `RNA only` = "RNA only", `Calcitriol occupancy` = "Calcitriol\noccupancy", `Vehicle occupancy` = "Vehicle\noccupancy", 
    `Shared occupancy` = "Shared\noccupancy", `Promoter-associated` = "Promoter-\nassociated")) + ggplot2::scale_x_continuous(limits = c(0, 
    panel_C_x_limit), expand = ggplot2::expansion(mult = c(0, 0))) + ggplot2::labs(tag = "C", title = "Occupancy support across\n436 RNA candidates", 
    subtitle = "Top two partition; lower attributes overlap", x = "RNA candidates", y = NULL) + theme_pub + ggplot2::theme(plot.title = ggplot2::element_text(size = header_title_size_pt, 
    face = "bold", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 1.2)), plot.subtitle = ggplot2::element_text(size = 7.2000000000000002, 
    color = "#444444", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 1)), plot.tag.position = c(top_block_tag_x, 
    top_block_tag_y), plot.tag.location = "panel", panel.grid.major.x = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.25), 
    panel.grid.major.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
    axis.text.y = ggplot2::element_text(size = 6.4000000000000004, lineheight = 0.90000000000000002, margin = ggplot2::margin(r = 1.2)), 
    plot.margin = ggplot2::margin(t = 3, r = 3, b = 4, l = 2))

panel_D <- ggplot2::ggplot(route_counts, ggplot2::aes(x = candidate_n, y = route_label, fill = route_label)) + ggplot2::geom_col(width = 0.68000000000000005, 
    color = "#333333", linewidth = 0.25) + ggplot2::geom_text(ggplot2::aes(label = candidate_n), hjust = -0.16, size = 2.3799999999999999, 
    family = figure_font, color = "#222222") + ggplot2::scale_fill_manual(values = route_palette, guide = "none") + ggplot2::scale_y_discrete(labels = c(`Shared + promoter` = "Shared +\npromoter", 
    `Shared, non-promoter` = "Shared,\nnon-promoter", `Non-shared promoter` = "Non-shared\npromoter", `Other occupancy support` = "Other occupancy\nsupport", 
    `RNA only` = "RNA only")) + ggplot2::scale_x_continuous(limits = c(0, panel_D_x_limit), expand = ggplot2::expansion(mult = c(0, 
    0))) + ggplot2::labs(tag = "D", title = "Mutually exclusive\nevidence routes", subtitle = "Routing categories; no direct-target assignment", 
    x = "RNA candidates", y = NULL) + theme_pub + ggplot2::theme(plot.title = ggplot2::element_text(size = header_title_size_pt, 
    face = "bold", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 1.2)), plot.subtitle = ggplot2::element_text(size = 7.2000000000000002, 
    color = "#444444", hjust = 0, margin = ggplot2::margin(l = top_block_tag_to_title_gap_pt, b = 1)), plot.tag.position = c(top_block_tag_x, 
    top_block_tag_y), plot.tag.location = "panel", panel.grid.major.x = ggplot2::element_line(color = "#E5E5E5", linewidth = 0.25), 
    panel.grid.major.y = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), 
    axis.text.y = ggplot2::element_text(size = 6.4000000000000004, lineheight = 0.90000000000000002, margin = ggplot2::margin(r = 1.2)), 
    plot.margin = ggplot2::margin(t = 3, r = 2, b = 4, l = 1))

composite_design <- "\nAB\n##\nCD\n"

panel_header_spec <- data.frame(panel = c("A", "B", "C", "D"), title = c("VDR RNA interaction program", "Top interaction candidates (12 samples)", 
    "Occupancy support across 436 RNA candidates", "Mutually exclusive evidence routes"), subtitle = c("436 candidates \u00b7 FDR < 0.05 \u00b7 |log2FC| \u2265 0.5", 
    "Top 30 candidates from the perturbation analysis", "Top two partition; lower attributes overlap", "Routing categories; no direct-target assignment"), 
    body_left_edge = c("volcano_y_axis_line", "heatmap_left_border", "bar_fill_left_edge", "bar_fill_left_edge"), stringsAsFactors = FALSE)

header_A <- add_explicit_panel_header(panel_A, "A", panel_header_spec$title[panel_header_spec$panel == "A"], panel_header_spec$subtitle[panel_header_spec$panel == 
    "A"], TRUE)

header_B <- add_explicit_panel_header(panel_B, "B", panel_header_spec$title[panel_header_spec$panel == "B"], panel_header_spec$subtitle[panel_header_spec$panel == 
    "B"])

header_C <- add_explicit_panel_header(panel_C, "C", panel_header_spec$title[panel_header_spec$panel == "C"], panel_header_spec$subtitle[panel_header_spec$panel == 
    "C"])

header_D <- add_explicit_panel_header(panel_D, "D", panel_header_spec$title[panel_header_spec$panel == "D"], panel_header_spec$subtitle[panel_header_spec$panel == 
    "D"])

explicit_header_geometry <- do.call(rbind, list(header_A$geometry, header_B$geometry, header_C$geometry, header_D$geometry))

explicit_header_geometry <- merge(explicit_header_geometry, panel_header_spec[, c("panel", "body_left_edge")], by = "panel", 
    all.x = TRUE, sort = FALSE)

panel_A_layout <- patchwork::wrap_elements(full = header_A$gtable)

panel_B_layout <- patchwork::wrap_elements(full = header_B$gtable)

panel_C_layout <- patchwork::wrap_elements(full = header_C$gtable)

panel_D_layout <- patchwork::wrap_elements(full = header_D$gtable)

figure_object <- (panel_A_layout + panel_B_layout + panel_C_layout + panel_D_layout) + patchwork::plot_layout(design = composite_design, 
    widths = c(1, 1), heights = c(1.1000000000000001, 0.044999999999999998, 0.85499999999999998), guides = "keep") + patchwork::plot_annotation(theme = ggplot2::theme(plot.margin = ggplot2::margin(3, 
    3, 3, 3)))

grDevices::dev.off()
figure_export(figure_object,'Figure2',180,170)
