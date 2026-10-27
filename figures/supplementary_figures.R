#!/usr/bin/env Rscript
# Plot existing results; no alignment, model fitting or hypothesis testing.
options(stringsAsFactors=FALSE,warn=1)
.frames <- sys.frames()
.paths <- vapply(.frames,function(f)if(is.null(f$ofile)) '' else as.character(f$ofile),character(1))
.paths <- .paths[basename(.paths)=="supplementary_figures.R" & file.exists(.paths)]
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
target <- setdiff(commandArgs(TRUE),'--check')
target <- if(length(target)) toupper(target[[1]]) else 'ALL'
if(!target %in% c('ALL','S1','S2','S3','S4')) stop('Choose ALL, S1, S2, S3 or S4.')
ids <- if(target=='ALL') paste0('S',1:4) else target
for(id in ids) figure_preflight(id,c('data.table','ggplot2','patchwork','gtable','ragg','systemfonts','scales','ggrepel','ggrastr','R.utils'))
if('--check' %in% commandArgs(TRUE)){message('Supplementary inputs, packages and font available');quit(status=0)}
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel);library(scales)})
assert_file <- function(path) {
    if (!file.exists(path)) 
        stop("Missing authority file: ", path)
    invisible(path)
}

assert_columns <- function(x, cols, label) {
    missing <- setdiff(cols, names(x))
    if (length(missing) > 0L) {
        stop(label, " is missing columns: ", paste(missing, collapse = ", "))
    }
    invisible(TRUE)
}

theme_pub <- function(base_size = 7.5) {
    theme_classic(base_family = "Arial", base_size = base_size) + theme(panel.border = element_rect(colour = "#333333", fill = NA, 
        linewidth = 0.34999999999999998), axis.line = element_blank(), axis.ticks = element_line(linewidth = 0.29999999999999999, 
        colour = "#333333"), axis.title = element_text(size = 7.5, colour = "#222222"), axis.text = element_text(size = 6.5, 
        colour = "#333333"), legend.title = element_text(size = 7, face = "bold"), legend.text = element_text(size = 6.5), 
        legend.key = element_blank(), plot.title = element_text(size = 8.5, face = "bold", hjust = 0, margin = margin(b = 1, 
            unit = "mm")), plot.subtitle = element_text(size = 7.2000000000000002, hjust = 0, colour = "#444444", margin = margin(b = 1, 
            unit = "mm")), plot.title.position = "panel", plot.margin = margin(0.59999999999999998, 3, 3, 3, unit = "mm"))
}

panel_title <- function(tag, title) paste0(tag, "   ", title)

ARTICLE_FONT <- "Arial"

PANEL_TAG_PT <- 9

PANEL_TITLE_PT <- 8.5

PANEL_SUBTITLE_PT <- 7.2000000000000002

TAG_TITLE_GAP_MM <- 2

panel_module <- function(plot, tag) {
    full_title <- as.character(plot$labels$title)
    if (length(full_title) != 1L || !nzchar(full_title)) {
        stop("Panel ", tag, " is missing its locked title")
    }
    title_text <- sub(paste0("^", tag, "[[:space:]]+"), "", full_title)
    subtitle_text <- as.character(plot$labels$subtitle)
    if (length(subtitle_text) == 0L || is.na(subtitle_text)) 
        subtitle_text <- ""
    body <- plot + labs(title = NULL, subtitle = NULL) + theme(plot.title = element_blank(), plot.subtitle = element_blank())
    gt <- ggplotGrob(body)
    panel_rows <- gt$layout[grepl("^panel", gt$layout$name), , drop = FALSE]
    if (nrow(panel_rows) == 0L) 
        stop("Could not locate scientific body for panel ", tag)
    body_left_col <- min(panel_rows$l)
    title_lines <- length(strsplit(title_text, "\n", fixed = TRUE)[[1]])
    subtitle_lines <- if (nzchar(subtitle_text)) {
        length(strsplit(subtitle_text, "\n", fixed = TRUE)[[1]])
    }
    else 0L
    title_row_mm <- 3.5499999999999998 * title_lines + 0.55000000000000004
    subtitle_row_mm <- if (subtitle_lines > 0L) 
        3.0499999999999998 * subtitle_lines + 0.75
    else 0.45000000000000001
    gt <- gtable::gtable_add_rows(gt, grid::unit(c(title_row_mm, subtitle_row_mm), "mm"), pos = 0)
    baseline_y <- grid::unit(0.71999999999999997, "mm")
    tag_grob <- grid::textGrob(tag, x = grid::unit(-TAG_TITLE_GAP_MM, "mm"), y = baseline_y, hjust = 1, vjust = 0, gp = grid::gpar(fontfamily = ARTICLE_FONT, 
        fontsize = PANEL_TAG_PT, fontface = "bold", col = "#111111"))
    title_grob <- grid::textGrob(title_text, x = grid::unit(0, "mm"), y = baseline_y, hjust = 0, vjust = 0, gp = grid::gpar(fontfamily = ARTICLE_FONT, 
        fontsize = PANEL_TITLE_PT, fontface = "bold", col = "#111111", lineheight = 0.97999999999999998))
    subtitle_grob <- grid::textGrob(subtitle_text, x = grid::unit(0, "mm"), y = grid::unit(0.62, "npc"), hjust = 0, vjust = 0.5, 
        gp = grid::gpar(fontfamily = ARTICLE_FONT, fontsize = PANEL_SUBTITLE_PT, col = "#444444", lineheight = 0.95999999999999996))
    gt <- gtable::gtable_add_grob(gt, grobs = list(tag_grob, title_grob), t = c(1, 1), l = c(body_left_col, body_left_col), 
        z = c(Inf, Inf), clip = "off", name = c(paste0("panel-tag-", tag), paste0("panel-title-", tag)))
    if (nzchar(subtitle_text)) {
        gt <- gtable::gtable_add_grob(gt, grobs = subtitle_grob, t = 2, l = body_left_col, z = Inf, clip = "off", name = paste0("panel-subtitle-", 
            tag))
    }
    patchwork::wrap_elements(full = gt, clip = FALSE)
}

group_palette <- c(Control = "#0072B2", RPL = "#D55E00", UI = "#009E73")

broad_palette <- c(`B cells` = "#CC79A7", `Ciliated epithelial cells` = "#E69F00", `Endothelial cells` = "#009E73", `Mast cells` = "#8B6C00", 
    `Myeloid cells` = "#56B4E9", `NK/T cells` = "#0072B2", `Stromal cells` = "#6A3D9A", `Unciliated epithelial cells` = "#D55E00", 
    `Low-confidence epithelial-like` = "#7F7F7F")

short_endpoint <- c(VDR_expression = "VDR expression", RNA_ALL_436 = "RNA 436", OCCUPANCY_SUPPORTED_284 = "Occupancy 284", 
    SHARED_OCCUPANCY_46 = "Shared 46", PROMOTER_OCCUPANCY_84 = "Promoter 84")

endpoint_order <- names(short_endpoint)

cell_type_order <- c("Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells")

short_cell_type <- c(`Endothelial cells` = "Endothelial", `Myeloid cells` = "Myeloid", `NK/T cells` = "NK/T", `Stromal cells` = "Stromal", 
    `Unciliated epithelial cells` = "Unciliated epi.")

deterministic_subset <- function(dt, group_col, id_col, max_per_group) {
    x <- copy(dt)
    setorderv(x, c(group_col, id_col))
    x[, {
        n_take <- min(.N, max_per_group)
        idx <- unique(as.integer(round(seq(1, .N, length.out = n_take))))
        .SD[idx]
    }, by = group_col]
}

build_s1 <- function() {
pca_path <- figure_input("S1", "pca_path")

cent_path <- figure_input("S1", "cent_path")

alloc_path <- figure_input("S1", "alloc_path")

loo_path <- figure_input("S1", "loo_path")

lapply(c(pca_path, cent_path, alloc_path, loo_path), assert_file)

pca <- fread(pca_path)

cent <- fread(cent_path)

alloc <- fread(alloc_path)

loo <- fread(loo_path)

assert_columns(pca, c("geo_accession", "PC1", "PC2", "group", "slide_short"), "STEP00 PCA")

assert_columns(cent, c("PC1", "PC2", "group", "slide_label"), "STEP00 centroids")

assert_columns(alloc, c("slide_barcode", "group", "n_samples"), "STEP00 allocation")

assert_columns(loo, c("comparison", "spearman_logFC", "omitted_group"), "STEP00 LOO")

stopifnot(nrow(pca) == 72L, nrow(cent) == 9L, nrow(alloc) == 9L, nrow(loo) == 27L)

pca[, `:=`(group, factor(group, levels = names(group_palette)))]

cent[, `:=`(group, factor(group, levels = names(group_palette)))]

p_a <- ggplot(pca, aes(PC1, PC2, colour = group, shape = group)) + geom_point(size = 1.25, alpha = 0.71999999999999997) + 
  ggrepel::geom_text_repel(data = cent, aes(label = slide_label, colour = group), seed = 180801, size = 2.0499999999999998, 
      fontface = "bold", show.legend = FALSE, min.segment.length = 0, box.padding = 0.25, point.padding = 0.20000000000000001, 
      segment.size = 0.25, max.overlaps = Inf) + scale_colour_manual(values = group_palette, drop = FALSE) + scale_shape_manual(values = c(Control = 16, 
  RPL = 17, UI = 15), drop = FALSE) + labs(title = panel_title("A", "PCA and physical-slide centroids"), subtitle = "72 samples; labels identify nine slide centroids", 
  x = "PC1 (33.3%)", y = "PC2 (23.3%)", colour = "Group", shape = "Group") + theme_pub() + theme(axis.title.x = element_text(margin = margin(t = 0.34999999999999998, 
  unit = "mm")), axis.title.y = element_text(margin = margin(r = 0.34999999999999998, unit = "mm")), legend.position = "bottom", 
  legend.box = "horizontal", legend.box.spacing = unit(0.14999999999999999, "mm"), legend.margin = margin(0, 0, 0, 0, unit = "mm"), 
  plot.margin = margin(0.59999999999999998, 3, 0.80000000000000004, 3, unit = "mm")) + guides(colour = guide_legend(nrow = 1), 
  shape = guide_legend(nrow = 1))

alloc[, `:=`(group, factor(group, levels = names(group_palette)))]

alloc[, `:=`(slide_barcode_char, as.character(slide_barcode))]

alloc[, `:=`(slide_short, substr(slide_barcode_char, nchar(slide_barcode_char) - 3L, nchar(slide_barcode_char)))]

setorder(alloc, group, slide_short)

alloc[, `:=`(slide_short, factor(slide_short, levels = rev(slide_short)))]

p_b <- ggplot(alloc, aes(n_samples, slide_short, fill = group)) + geom_col(width = 0.71999999999999997) + geom_text(aes(label = paste0("n=", 
  n_samples)), hjust = 1.1499999999999999, colour = "white", size = 2.2000000000000002, family = "Arial") + scale_fill_manual(values = group_palette, 
  drop = FALSE) + scale_x_continuous(limits = c(0, 8.4000000000000004), breaks = c(0, 4, 8), expand = c(0, 0)) + labs(title = panel_title("B", 
  "Phenotype allocation across slides"), subtitle = "Each physical slide contains one phenotype group", x = "Samples", 
  y = "Slide suffix", fill = "Group") + theme_pub() + theme(legend.position = "none")

loo[, `:=`(comparison, factor(comparison, levels = c("UI_vs_Control", "RPL_vs_Control", "UI_vs_RPL"), labels = c("UI vs Control", 
  "RPL vs Control", "UI vs RPL")))]

loo[, `:=`(omitted_group, factor(omitted_group, levels = names(group_palette)))]

p_c <- ggplot(loo, aes(comparison, spearman_logFC, shape = omitted_group)) + geom_hline(yintercept = 0.84999999999999998, 
  linetype = 2, linewidth = 0.34999999999999998) + geom_point(position = position_jitter(width = 0.11, height = 0, seed = 180802), 
  size = 1.5, alpha = 0.81999999999999995, colour = "#2B2B2B") + scale_shape_manual(values = c(Control = 16, RPL = 17, 
  UI = 15), drop = FALSE) + scale_y_continuous(limits = c(0.84499999999999997, 1.0049999999999999), breaks = c(0.84999999999999998, 
  0.90000000000000002, 0.94999999999999996, 1)) + labs(title = panel_title("C", "Leave-one-slide-out robustness"), subtitle = "Dashed line: prespecified Spearman threshold of 0.85", 
  x = NULL, y = "Spearman correlation", shape = "Omitted group") + theme_pub() + theme(axis.text.x = element_text(angle = 18, 
  hjust = 1), legend.position = "bottom", legend.box.spacing = unit(0.14999999999999999, "mm"), legend.margin = margin(0, 
  0, 0, 0, unit = "mm"), plot.margin = margin(0.59999999999999998, 3, 0.80000000000000004, 3, unit = "mm")) + guides(shape = guide_legend(nrow = 1))

p_a <- panel_module(p_a, "A")

p_b <- panel_module(p_b, "B")

p_c <- panel_module(p_c, "C")

design <- c(patchwork::area(t = 1, l = 1, b = 2, r = 1), patchwork::area(t = 1, l = 2, b = 1, r = 2), patchwork::area(t = 2, 
  l = 2, b = 2, r = 2))

composite <- wrap_plots(p_a, p_b, p_c, design = design, widths = c(1.28, 1), heights = c(1, 1))
figure_export(composite,'FigS1',183,140)
}

build_s2 <- function() {
cell_path <- figure_input("S2", "cell_path")

doublet_path <- figure_input("S2", "doublet_path")

mixing_path <- figure_input("S2", "mixing_path")

marker_path <- figure_input("S2", "marker_path")

annotation_path <- figure_input("S2", "annotation_path")

lapply(c(cell_path, doublet_path, mixing_path, marker_path, annotation_path), assert_file)

cells <- fread(cell_path, select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters", "broad_cell_type", "sample_accession", 
  "clinical_group", "LH_day"))

assert_columns(cells, c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters", "broad_cell_type", "sample_accession"), 
  "GSE250130 cell metadata")

stopifnot(nrow(cells) == 236636L, uniqueN(cells$cell_barcode) == nrow(cells))

cells[, `:=`(cluster, as.character(seurat_clusters))]

display_cells <- deterministic_subset(cells, "cluster", "cell_barcode", 2500L)

broad_counts <- cells[, .N, by = broad_cell_type][order(-N)]

display_cells[, `:=`(draw_order, match(broad_cell_type, broad_counts$broad_cell_type))]

setorder(display_cells, draw_order, cell_barcode)

display_cells[, `:=`(broad_cell_type, factor(broad_cell_type, levels = names(broad_palette)))]

centres <- display_cells[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)), by = cluster]

cluster_palette <- setNames(grDevices::hcl.colors(uniqueN(display_cells$cluster), "Dynamic"), sort(unique(display_cells$cluster)))

p_a <- ggplot(display_cells, aes(UMAPH_1, UMAPH_2, colour = cluster)) + geom_point(size = 0.10000000000000001, alpha = 0.41999999999999998, 
  stroke = 0) + ggrepel::geom_text_repel(data = centres, aes(x = UMAPH_1, y = UMAPH_2, label = cluster), inherit.aes = FALSE, 
  seed = 180803, size = 2.2000000000000002, fontface = "bold", colour = "#1A1A1A", min.segment.length = 0, segment.size = 0.20000000000000001, 
  box.padding = 0.17999999999999999, point.padding = 0.12, max.overlaps = Inf) + scale_colour_manual(values = cluster_palette, 
  guide = "none") + coord_equal() + labs(title = panel_title("A", "Integrated cluster map"), subtitle = paste0(format(nrow(display_cells), 
  big.mark = ","), "-cell subset; direct cluster labels"), x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + theme_pub() + 
  theme(legend.position = "none")

sample_design <- unique(cells[, .(sample_accession, clinical_group, LH_day)])

doublet <- fread(doublet_path)

doublet <- merge(doublet, sample_design, by = "sample_accession", all.x = TRUE, sort = FALSE)

stopifnot(nrow(doublet) == 28L, all(!is.na(doublet$clinical_group)))

doublet[, `:=`(sample_class, ifelse(toupper(clinical_group) == "RIF", "RIF LH+7", paste0("Fertile LH+", LH_day)))]

class_order <- c("Fertile LH+3", "Fertile LH+5", "Fertile LH+7", "Fertile LH+9", "Fertile LH+11", "RIF LH+7")

doublet[, `:=`(sample_class, factor(sample_class, levels = class_order))]

doublet[, `:=`(analysis_group, ifelse(grepl("^RIF", sample_class), "RIF", "Fertile"))]

p_c <- ggplot(doublet, aes(sample_class, predicted_doublet_fraction, shape = analysis_group, colour = analysis_group)) + 
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 180804), size = 1.5, alpha = 0.84999999999999998) + 
  scale_colour_manual(values = c(Fertile = "#0072B2", RIF = "#D55E00")) + scale_shape_manual(values = c(Fertile = 16, RIF = 17)) + 
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.029999999999999999, 0.080000000000000002))) + 
  labs(title = panel_title("C", "Sample-level doublet removal"), subtitle = "Blue circles: fertile; orange triangles: RIF", 
      x = NULL, y = "Predicted doublet fraction", colour = NULL, shape = NULL) + theme_pub() + theme(axis.text.x = element_text(angle = 32, 
  hjust = 1), legend.position = "none", plot.margin = margin(0.59999999999999998, 3, 0.59999999999999998, 3, unit = "mm"))

mixing <- fread(mixing_path)

stopifnot(nrow(mixing) == 22L)

mixing[, `:=`(review_flag, as.logical(review_flag))]

label_mix <- mixing[review_flag %in% TRUE]

p_b <- ggplot(mixing, aes(max_sample_fraction, normalized_sample_entropy)) + geom_vline(xintercept = c(0.29999999999999999, 
  0.5), linetype = c(2, 1), linewidth = 0.34999999999999998, colour = "#555555") + geom_hline(yintercept = 0.65000000000000002, 
  linetype = 2, linewidth = 0.34999999999999998, colour = "#555555") + geom_point(aes(size = cell_n, shape = review_flag), 
  alpha = 0.81999999999999995, colour = "#2B2B2B") + ggrepel::geom_text_repel(data = label_mix, aes(label = seurat_clusters), 
  seed = 180805, size = 2.2000000000000002, min.segment.length = 0, segment.size = 0.25, box.padding = 0.25, max.overlaps = Inf) + 
  scale_x_continuous(labels = percent_format(accuracy = 1)) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) + 
  scale_size_continuous(range = c(1.2, 4), labels = comma) + labs(title = panel_title("B", "Sample-mixing QC"), subtitle = "Prespecified review clusters are labeled", 
  x = "Largest single-sample fraction", y = "Normalized sample entropy", shape = "Review", size = "Cells") + theme_pub() + 
  theme(legend.position = "bottom", legend.box = "horizontal", legend.box.spacing = unit(0.14999999999999999, "mm"), legend.margin = margin(0, 
      0, 0, 0, unit = "mm"), plot.margin = margin(0.59999999999999998, 3, 1, 3, unit = "mm")) + guides(size = "none", shape = guide_legend(nrow = 1))

marker <- fread(marker_path)

stopifnot(nrow(marker) == 352L, uniqueN(marker$cluster) == 22L)

lineage_labels <- c(B_cell = "B", ciliated_epithelial = "Ciliated", cycling = "Cycling", decidual_stromal = "Decidual", dendritic_cell = "Dendritic", 
  endothelial = "Endothelial", epithelial_general = "Epithelial", erythroid = "Erythroid", immune_general = "Immune", lymphatic_endothelial = "Lymphatic endo.", 
  mast_cell = "Mast", myeloid_macrophage = "Myeloid", NK_cell = "NK", perivascular_smooth_muscle = "Perivascular", plasma_cell = "Plasma", 
  secretory_epithelial = "Secretory", stromal_fibroblast = "Stromal", T_cell = "T")

marker[, `:=`(panel_label, lineage_labels[panel])]

marker[, `:=`(cluster, factor(cluster, levels = rev(as.character(0:21))))]

marker[, `:=`(panel_label, factor(panel_label, levels = unname(lineage_labels)))]

p_d <- ggplot(marker, aes(panel_label, cluster, fill = mean_scaled_expression)) + geom_tile(colour = "white", linewidth = 0.12) + 
  scale_fill_gradientn(colours = c("#3B0F70", "#2C7BB6", "#22A884", "#FDE725"), name = "Mean scaled\nexpression") + labs(title = panel_title("D", 
  "Broad-lineage marker support"), subtitle = "Cluster-by-panel mean scaled expression; annotation support, not cell-level inference", 
  x = "Canonical marker panel", y = "Cluster") + theme_pub() + theme(axis.text.x = element_text(angle = 42, hjust = 1, 
  vjust = 1, size = 5.7999999999999998), axis.text.y = element_text(size = 5.7999999999999998), legend.position = "right", 
  plot.margin = margin(0.20000000000000001, 3, 3, 3, unit = "mm"))

p_a_pdf <- ggrastr::rasterise(p_a, layers = "Point", dpi = 600, dev = "ragg")

p_a <- panel_module(p_a, "A")

p_a_pdf <- panel_module(p_a_pdf, "A")

p_b <- panel_module(p_b, "B")

p_c <- panel_module(p_c, "C")

p_d <- panel_module(p_d, "D")

design <- c(patchwork::area(t = 1, l = 1, b = 1, r = 1), patchwork::area(t = 1, l = 2, b = 2, r = 2), patchwork::area(t = 2, 
  l = 1, b = 2, r = 1), patchwork::area(t = 3, l = 1, b = 3, r = 2))

composite <- wrap_plots(p_a, p_b, p_c, p_d, design = design, heights = c(1, 0.94999999999999996, 1.23), widths = c(1.1200000000000001, 
  0.88))

composite_pdf <- wrap_plots(p_a_pdf, p_b, p_c, p_d, design = design, heights = c(1, 0.94999999999999996, 1.23), widths = c(1.1200000000000001, 
  0.88))
figure_export(composite_pdf,'FigS2',183,188)
}

build_s3 <- function() {
umap_path <- figure_input("S3", "umap_path")

annotation_map_path <- figure_input("S3", "annotation_map_path")

mixing_path <- figure_input("S3", "mixing_path")

annotation_qc_path <- figure_input("S3", "annotation_qc_path")

sensitivity_path <- figure_input("S3", "sensitivity_path")

lapply(c(umap_path, annotation_map_path, mixing_path, annotation_qc_path, sensitivity_path), assert_file)

umap <- fread(umap_path, select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "sample_accession", "analysis_group", "seurat_clusters"))

ann_map <- fread(annotation_map_path)

umap <- merge(umap, ann_map[, .(cluster, broad_cell_type, annotation_confidence, review_flag)], by.x = "seurat_clusters", 
  by.y = "cluster", all.x = TRUE, sort = FALSE)

stopifnot(nrow(umap) == 47161L, all(!is.na(umap$broad_cell_type)))

umap[, `:=`(cluster, as.character(seurat_clusters))]

display_umap <- deterministic_subset(umap, "broad_cell_type", "cell_barcode", 6000L)

broad_counts <- umap[, .N, by = broad_cell_type][order(-N)]

display_umap[, `:=`(draw_order, match(broad_cell_type, broad_counts$broad_cell_type))]

setorder(display_umap, draw_order, cell_barcode)

display_umap[, `:=`(broad_cell_type, factor(broad_cell_type, levels = names(broad_palette)))]

centres <- display_umap[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)), by = cluster]

centres[, `:=`(label_colour, ifelse(cluster == "8", "#D55E00", "#1A1A1A"))]

p_a <- ggplot(display_umap, aes(UMAPH_1, UMAPH_2, colour = broad_cell_type)) + geom_point(size = 0.14000000000000001, alpha = 0.46000000000000002, 
  stroke = 0) + ggrepel::geom_text_repel(data = centres, aes(x = UMAPH_1, y = UMAPH_2, label = cluster, colour = I(label_colour)), 
  inherit.aes = FALSE, seed = 180806, size = 2.25, fontface = "bold", min.segment.length = 0, segment.size = 0.20000000000000001, 
  box.padding = 0.20000000000000001, point.padding = 0.12, max.overlaps = Inf, show.legend = FALSE) + scale_colour_manual(values = broad_palette, 
  drop = FALSE) + coord_equal() + labs(title = panel_title("A", "Integrated UMAP and annotations"), subtitle = "Cluster 8 (orange label) is low confidence\nand excluded from primary models", 
  x = "Harmony UMAP 1", y = "Harmony UMAP 2", colour = "Broad cell type") + theme_pub() + theme(legend.position = "none")

mixing <- fread(mixing_path)

stopifnot(nrow(mixing) == 16L)

mixing[, `:=`(review_required, as.logical(review_required))]

p_b <- ggplot(mixing, aes(top_donor_fraction, donor_entropy_normalized)) + geom_vline(xintercept = c(0.5, 0.75), linetype = 2, 
  linewidth = 0.34999999999999998, colour = "#555555") + geom_hline(yintercept = c(0.5, 0.69999999999999996), linetype = 3, 
  linewidth = 0.34999999999999998, colour = "#555555") + geom_point(aes(size = cell_n, shape = review_required), alpha = 0.81999999999999995, 
  colour = "#2B2B2B") + ggrepel::geom_text_repel(data = mixing[review_required %in% TRUE], aes(label = seurat_clusters), 
  seed = 180807, size = 2.2000000000000002, min.segment.length = 0, segment.size = 0.25, box.padding = 0.25, max.overlaps = Inf) + 
  scale_x_continuous(labels = percent_format(accuracy = 1)) + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) + 
  scale_size_continuous(range = c(1.2, 4), labels = comma) + labs(title = panel_title("B", "Donor-mixing diagnostics"), 
  subtitle = "Labels identify clusters requiring\ndonor/technical review", x = "Largest single-donor fraction", y = "Normalized donor entropy", 
  shape = "Review", size = "Cells") + guides(size = "none", shape = guide_legend(nrow = 1)) + theme_pub() + theme(legend.position = "bottom", 
  legend.box = "vertical", legend.box.spacing = unit(1.25, "mm"), legend.margin = margin(0, 0, 0, 0, unit = "mm"), legend.box.margin = margin(t = 0, 
      unit = "mm"), plot.margin = margin(0.59999999999999998, 3, 0.80000000000000004, 3, unit = "mm"))

ann_qc <- fread(annotation_qc_path)

stopifnot(nrow(ann_qc) == 16L)

ann_qc[, `:=`(cluster_factor, factor(cluster, levels = as.character(0:15)))]

p_c <- ggplot(ann_qc, aes(cluster_factor, top100_evidence_n, shape = annotation_confidence, size = cell_n)) + geom_point(colour = "#2B2B2B", 
  alpha = 0.81999999999999995) + geom_point(data = ann_qc[cluster == 8], colour = "#D55E00", stroke = 0.90000000000000002) + 
  ggrepel::geom_text_repel(data = ann_qc[cluster == 8], aes(label = "cluster 8"), seed = 180808, size = 2.2000000000000002, 
      colour = "#D55E00", min.segment.length = 0, segment.size = 0.25) + scale_shape_manual(values = c(HIGH = 16, MODERATE = 15, 
  LOW = 17)) + scale_size_continuous(range = c(1.2, 4.2000000000000002), labels = comma) + scale_y_continuous(breaks = 0:4, 
  limits = c(0, 4.25)) + labs(title = panel_title("C", "Annotation evidence and confidence"), subtitle = "Requested evidence markers among the top-100 marker summaries", 
  x = "Primary cluster", y = "Evidence markers present", shape = "Confidence", size = "Cells") + guides(size = "none", 
  shape = guide_legend(nrow = 1)) + theme_pub() + theme(axis.title.x = element_text(margin = margin(t = 1.25, unit = "mm")), 
  legend.position = "bottom", legend.box = "vertical", legend.box.spacing = unit(0.050000000000000003, "mm"), legend.margin = margin(0, 
      0, 0, 0, unit = "mm"), legend.box.margin = margin(t = -0.25, unit = "mm"), plot.margin = margin(0.59999999999999998, 
      3, 0.69999999999999996, 3, unit = "mm"))

sensitivity <- fread(sensitivity_path)

stopifnot(nrow(sensitivity) == 5L, all(as.logical(sensitivity$direction_concordant)))

sens_long <- rbindlist(list(sensitivity[, .(endpoint, analysis_matrix = "Primary (cluster 8 excluded)", effect = mean_difference_RIF_minus_Control_primary, 
  exact_P = exact_permutation_P_primary, loo = LOO_direction_concordance_primary)], sensitivity[, .(endpoint, analysis_matrix = "Atlas sensitivity (cluster 8 retained)", 
  effect = mean_difference_RIF_minus_Control_atlas, exact_P = exact_permutation_P_atlas, loo = LOO_direction_concordance_atlas)]))

sens_long[, `:=`(endpoint_label, factor(short_endpoint[endpoint], levels = rev(unname(short_endpoint))))]

sens_long[, `:=`(analysis_matrix, factor(analysis_matrix, levels = c("Primary (cluster 8 excluded)", "Atlas sensitivity (cluster 8 retained)"), 
  labels = c("Primary: cluster 8 excluded", "Sensitivity: cluster 8 retained")))]

p_d <- ggplot(sens_long, aes(effect, endpoint_label, shape = analysis_matrix, colour = analysis_matrix)) + geom_vline(xintercept = 0, 
  linetype = 2, linewidth = 0.34999999999999998) + geom_point(size = 2.1000000000000001, stroke = 0.80000000000000004, 
  position = position_dodge(width = 0.45000000000000001)) + scale_colour_manual(values = c(`Primary: cluster 8 excluded` = "#0072B2", 
  `Sensitivity: cluster 8 retained` = "#D55E00")) + scale_shape_manual(values = c(`Primary: cluster 8 excluded` = 16, `Sensitivity: cluster 8 retained` = 1)) + 
  labs(title = panel_title("D", "Cluster-8 sensitivity"), subtitle = "All five directions are concordant;\nalternative definitions are not connected", 
      x = "RIF minus Control mean difference", y = NULL, colour = "Analysis", shape = "Analysis") + guides(colour = guide_legend(nrow = 2, 
  byrow = TRUE, keyheight = unit(2.7000000000000002, "mm")), shape = guide_legend(nrow = 2, byrow = TRUE, keyheight = unit(2.7000000000000002, 
  "mm"))) + theme_pub() + theme(legend.position = "bottom", legend.box = "vertical", legend.box.spacing = unit(1.1499999999999999, 
  "mm"), legend.key.height = unit(2.7000000000000002, "mm"), legend.spacing.y = unit(0.050000000000000003, "mm"), legend.margin = margin(0, 
  0, 0, 0, unit = "mm"), legend.box.margin = margin(t = 0, unit = "mm"), plot.margin = margin(0.59999999999999998, 3, 0.69999999999999996, 
  3, unit = "mm"))

p_a_pdf <- ggrastr::rasterise(p_a, layers = "Point", dpi = 600, dev = "ragg")

p_a <- panel_module(p_a, "A")

p_a_pdf <- panel_module(p_a_pdf, "A")

p_b <- panel_module(p_b, "B")

p_c <- panel_module(p_c, "C")

p_d <- panel_module(p_d, "D")

composite <- (p_a + p_b)/(p_c + p_d) + plot_layout(widths = c(0.95999999999999996, 1.04), heights = c(1.03, 0.96999999999999997))

composite_pdf <- (p_a_pdf + p_b)/(p_c + p_d) + plot_layout(widths = c(0.95999999999999996, 1.04), heights = c(1.03, 0.96999999999999997))
figure_export(composite_pdf,'FigS3',183,180)
}

build_s4 <- function() {
cross_path <- figure_input("S4", "cross_path")

overall_path <- figure_input("S4", "overall_path")

evidence_path <- figure_input("S4", "evidence_path")

lapply(c(cross_path, overall_path, evidence_path), assert_file)

cross <- fread(cross_path)

overall <- fread(overall_path)

evidence <- fread(evidence_path)

stopifnot(nrow(cross) == 25L, nrow(overall) == 13L, nrow(evidence) == 4L)

long <- rbindlist(list(cross[, .(broad_cell_type, endpoint, dataset = "GSE250130", exact_P = exact_P_G250, q_value = primary_q_G250, 
  q_scope = primary_q_scope_G250, loo = LOO_direction_concordance_G250, constant_endpoint = constant_endpoint_G250)], cross[, 
  .(broad_cell_type, endpoint, dataset = "GSE183837", exact_P = exact_P_G183, q_value = primary_q_G183, q_scope = primary_q_scope_G183, 
      loo = LOO_direction_concordance_G183, constant_endpoint = constant_endpoint_G183)]))

long[, `:=`(cell_label, factor(short_cell_type[broad_cell_type], levels = unname(short_cell_type[cell_type_order])))]

long[, `:=`(endpoint_label, factor(short_endpoint[endpoint], levels = rev(unname(short_endpoint[endpoint_order]))))]

long[, `:=`(dataset, factor(dataset, levels = c("GSE250130", "GSE183837")))]

long[, `:=`(p_text, ifelse(is.na(exact_P), "NA", formatC(exact_P, format = "f", digits = 3)))]

long[, `:=`(q_text, ifelse(is.na(q_value), "NA", formatC(q_value, format = "f", digits = 3)))]

long[, `:=`(loo_text, ifelse(is.na(loo), "NA", formatC(loo, format = "f", digits = 2)))]

long[, `:=`(p_fill, -log10(pmax(exact_P, 9.9999999999999995e-07)))]

long[, `:=`(q_fill, -log10(pmax(q_value, 9.9999999999999995e-07)))]

heat_theme <- theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6), axis.text.y = element_text(size = 6), 
  strip.text = element_text(size = 7.2000000000000002, face = "bold"), strip.background = element_rect(fill = "#F2F2F2", 
      colour = "#666666", linewidth = 0.29999999999999999), legend.position = "right", legend.justification = c(0.5, 1), 
  legend.box.just = "top", legend.margin = margin(0, 0, 0, 1, unit = "mm"))

p_a <- ggplot(long, aes(cell_label, endpoint_label, fill = p_fill)) + geom_tile(colour = "white", linewidth = 0.34999999999999998) + 
  geom_text(aes(label = p_text), family = "Arial", size = 2) + facet_grid(. ~ dataset) + scale_fill_gradientn(colours = c("#F3F3F3", 
  "#92C5DE", "#0571B0"), limits = c(0, max(long$p_fill, na.rm = TRUE)), name = expression(-log[10](italic(P))), guide = guide_colourbar(barheight = unit(27.300000000000001, 
  "mm"), barwidth = unit(4, "mm"), title.position = "top")) + labs(title = panel_title("A", "Exact-permutation P values for all 25 models"), 
  subtitle = "Tiles give exact P values; no cross-cohort combination is performed", x = NULL, y = NULL) + heat_theme

p_b <- ggplot(long, aes(cell_label, endpoint_label, fill = q_fill)) + geom_tile(colour = "white", linewidth = 0.34999999999999998) + 
  geom_text(aes(label = q_text), family = "Arial", size = 2) + facet_grid(. ~ dataset) + scale_fill_gradientn(colours = c("#F3F3F3", 
  "#FDB863", "#B2182B"), limits = c(0, max(long$q_fill, na.rm = TRUE)), name = expression(-log[10](italic(q))), guide = guide_colourbar(barheight = unit(27.300000000000001, 
  "mm"), barwidth = unit(4, "mm"), title.position = "top")) + labs(title = panel_title("B", "Prespecified within-cohort adjusted q values"), 
  subtitle = "GSE250130 uses endpoint-family BH; GSE183837 uses global BH across 25 models", x = NULL, y = NULL) + heat_theme

p_c <- ggplot(long, aes(cell_label, endpoint_label, fill = loo)) + geom_tile(colour = "white", linewidth = 0.34999999999999998) + 
  geom_text(aes(label = loo_text), family = "Arial", size = 2) + facet_grid(. ~ dataset) + scale_fill_gradientn(colours = c("#B2182B", 
  "#F7F7F7", "#2166AC"), limits = c(0, 1), na.value = "#E5E5E5", name = "LOO direction\nconcordance", guide = guide_colourbar(barheight = unit(25.800000000000001, 
  "mm"), barwidth = unit(4, "mm"), title.position = "top")) + labs(title = panel_title("C", "Leave-one-donor-out direction concordance"), 
  subtitle = "NA denotes a constant or non-informative endpoint; values are model diagnostics", x = NULL, y = NULL) + heat_theme

evidence[, `:=`(evidence_label, factor(evidence_class, levels = c("CONCORDANT_NO_STATISTICAL_SUPPORT", "DISCORDANT_NO_STATISTICAL_SUPPORT", 
  "DISCORDANT_WITH_NOMINAL_OR_CORRECTED_SIGNAL", "NONINFORMATIVE_CONSTANT_OR_ZERO"), labels = c("Concordant; no statistical support", 
  "Discordant; no statistical support", "Discordant; nominal signal in one cohort", "Non-informative constant/zero")))]

evidence[, `:=`(evidence_label, factor(evidence_label, levels = rev(levels(evidence_label))))]

evidence_cols <- c(`Concordant; no statistical support` = "#4DAF4A", `Discordant; no statistical support` = "#984EA3", `Discordant; nominal signal in one cohort` = "#E41A1C", 
  `Non-informative constant/zero` = "#999999")

p_d <- ggplot(evidence, aes(model_n, evidence_label, fill = evidence_label)) + geom_col(width = 0.68000000000000005) + geom_text(aes(label = model_n), 
  hjust = -0.34999999999999998, family = "Arial", size = 2.3999999999999999) + scale_fill_manual(values = evidence_cols, 
  guide = "none") + scale_x_continuous(limits = c(0, 13), breaks = c(0, 4, 8, 12), expand = c(0, 0)) + labs(title = panel_title("D", 
  "Cross-dataset evidence classes"), subtitle = "24 informative models: 12 concordant and 12 discordant; corrected replication = 0", 
  x = "Models", y = NULL) + theme_pub() + theme(axis.text.y = element_text(size = 6.2000000000000002))

p_a <- panel_module(p_a, "A")

p_b <- panel_module(p_b, "B")

p_c <- panel_module(p_c, "C")

p_d <- panel_module(p_d, "D")

composite <- p_a/p_b/p_c/p_d + plot_layout(heights = c(1, 1, 1, 0.71999999999999997))
figure_export(composite,'FigS4',183,220)
}

.metric_file <- tempfile(fileext='.png')
ragg::agg_png(.metric_file,width=800,height=800,units='px',res=150,background='white')
tryCatch({for(id in ids)get(paste0('build_',tolower(id)))()},finally=grDevices::dev.off())
