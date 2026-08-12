
# Publication figure stage 1
# ======================================================================
# 108_step05C_Figure3_data_body_aligned_colorbars_v1_4_8_QC.R
# Nutrition_Endometrial_VDR_v2
#
# Skill 3 v2.2-beta.6 GEOMETRY-CONSTRAINED legend alignment candidate.
# Reads ONLY frozen Figure 3 source CSVs; no statistics are recomputed.
# Produces ONE 300-dpi QC PNG and a SMALL_REPORT ZIP on Desktop.
#
# v1.4.8 RULE-DRIVEN CHANGE ONLY:
# 1) Panels A/D use custom code-layer horizontal colorbars, not raster patches.
# 2) IMPORTANT: alignment target is the ACTUAL DATA/TILE BODY right boundary,
#    not the outer panel/gtable/axis/whitespace boundary.
# 3) Colorbar BODY right edge is constrained to the right edge of the last tile.
# 4) Legend-title TEXT bottom edge is constrained to the colorbar BODY bottom edge.
# 5) Bar-body widths are preserved from the previously approved v1.3 actual render:
#      Panel A = 50.88 mm; Panel D = 59.94 mm.
#    This prevents Panel A from being stretched to the generic 60-mm width.
# 6) Endpoint tick labels are inward-justified so the maximum label cannot clip.
# 7) No guessed visual pixel shifts are used for placement.

#
# Frozen repair from v1.2 actual-render QC:
# 1) BLOCKING: the lower-row canvas ended too close to the device boundary.
#    Panel C x-axis title and Panel D colorbar tick labels were clipped at the
#    bottom edge in the actual 300-dpi PNG.
#    Repair: reserve a dedicated 8-mm empty bottom row inside the 184-mm canvas
#    and slightly rebalance body-row heights. No scientific/statistical content
#    or typography sizes are changed.
#
# SAVE AS:
# <repository module script>
#
# RUN:
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "108_step05C_Figure3_data_body_aligned_colorbars_v1_4_8_QC.R"
script_path <- file.path(project_root, "04_code", "R", script_name)
figure_version <- "F3_v1_4_8_QC"
source_version <- "F3_v1_2_QC"
source_lock <- "STEP04_LOCK_v1.2"
bulk_projection_version <- "BULK_VDR_PROGRAM_PROJECTION_v1.3"
governing_skill <- "research_scientific_figure_suite v2.2-beta.6"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

required_packages <- c(
  "ggplot2", "patchwork", "gtable", "grid", "scales",
  "systemfonts", "ragg", "digest", "zip"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

normalize_slashes <- function(x) gsub("\\\\", "/", x)
write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_text <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}
sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
copy_exact <- function(src, dst) {
  if (!file.exists(src)) stop("Missing source file: ", src)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dst) || sha256_file(src) != sha256_file(dst)) {
    stop("Validated copy failed: ", src, " -> ", dst)
  }
  invisible(TRUE)
}
find_desktop <- function() {
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
  )
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) stop("Could not locate Desktop.")
  candidates[[1]]
}

# ----------------------------------------------------------------------
# 1. Governance / Skill preflight
# ----------------------------------------------------------------------

skill_path <- file.path(
  "<FIGURE_STYLE_REFERENCE>",
  "research_scientific_figure_suite",
  "SKILL.md"
)
if (!file.exists(skill_path)) stop("Active Skill 3 file is missing: ", skill_path)
skill_text <- paste(readLines(skill_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.6", skill_text, fixed = TRUE)) {
  stop("Active Skill 3 is not v2.2-beta.6.")
}
if (!file.exists(script_path)) {
  stop(
    "Please save this code first as:\n", normalize_slashes(script_path),
    "\nThen run it with source(...)."
  )
}

figure_font <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}
font_match <- systemfonts::match_fonts(figure_font)
font_info <- systemfonts::font_info(font_match$path)
font_ok <- nrow(font_match) >= 1L && file.exists(font_match$path[[1]]) &&
  grepl("Arial", font_info$family[[1]], ignore.case = TRUE)
if (!font_ok) stop("Arial exact-family font preflight failed.")

# ----------------------------------------------------------------------
# 2. Read v1.1 figure-source tables only (NO statistical recomputation)
# ----------------------------------------------------------------------

step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")
fig3_parent <- file.path(step05_root, "06_results", "FIG3")
source_root <- file.path(fig3_parent, source_version, "source_data")
text_root_v11 <- file.path(fig3_parent, source_version, "manuscript_text")

out_root <- file.path(fig3_parent, figure_version)
fig_root <- file.path(out_root, "figures")
qc_root <- file.path(out_root, "qc")
code_root <- file.path(out_root, "code")
source_copy_root <- file.path(out_root, "source_data")
text_copy_root <- file.path(out_root, "manuscript_text")

if (dir.exists(out_root)) {
  stop(
    "Refusing to overwrite reviewed candidate folder:\n",
    normalize_slashes(out_root),
    "\nIf this is an incomplete failed run, inspect it manually before removing it."
  )
}
for (p in c(fig_root, qc_root, code_root, source_copy_root, text_copy_root)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

source_files_needed <- c(
  "Figure3_panelA_program_effect_heatmap.csv",
  "Figure3_panelB_FDR_effect_matrix.csv",
  "Figure3_panelC_candidate_context_coverage.csv",
  "Figure3_panelD_context_effects_long.csv",
  "Figure3_panelD_selected_examples.csv",
  "Figure3_locked_blueprint_row.csv",
  "Figure3_score_significance_summary.csv",
  "Figure3_candidate_direction_summary.csv"
)
source_paths <- file.path(source_root, source_files_needed)
if (!all(file.exists(source_paths))) {
  stop("One or more Figure 3 v1.2 QC source tables are missing.")
}

for (p in source_paths) {
  copy_exact(p, file.path(source_copy_root, basename(p)))
}

score_models <- read.csv(
  file.path(source_root, "Figure3_panelA_program_effect_heatmap.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
coverage <- read.csv(
  file.path(source_root, "Figure3_panelC_candidate_context_coverage.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
panel_D_data <- read.csv(
  file.path(source_root, "Figure3_panelD_context_effects_long.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
selected_examples <- read.csv(
  file.path(source_root, "Figure3_panelD_selected_examples.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)

# Locked structure checks only.
if (nrow(score_models) != 16L) stop("Expected 16 program-by-contrast rows.")
if (sum(as.logical(score_models$FDR_significant)) != 4L) stop("Expected 4 FDR-significant program rows.")
if (nrow(coverage) != 5L || sum(coverage$candidate_n) != 436L) stop("Coverage table mismatch.")
if (nrow(selected_examples) != 10L || nrow(panel_D_data) != 40L) stop("Panel D source mismatch.")

program_levels <- c(
  "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284",
  "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84"
)
program_labels <- c(
  "RNA_ALL_436" = "RNA program\n(n=436)",
  "OCCUPANCY_SUPPORTED_284" = "Occupancy-supported\n(n=284)",
  "SHARED_OCCUPANCY_46" = "Shared occupancy\n(n=46)",
  "PROMOTER_OCCUPANCY_84" = "Promoter-associated\n(n=84)"
)
context_ids <- c(
  "GSE165004__RPL_vs_CONTROL",
  "GSE165004__UIF_vs_CONTROL",
  "GSE92324__IF_vs_CONTROL",
  "GSE111974__RIF_vs_CONTROL"
)
context_labels <- c(
  "GSE165004__RPL_vs_CONTROL" = "GSE165004\nRPL",
  "GSE165004__UIF_vs_CONTROL" = "GSE165004\nUIF",
  "GSE92324__IF_vs_CONTROL" = "GSE92324\nIF",
  "GSE111974__RIF_vs_CONTROL" = "GSE111974\nRIF"
)

score_models$program_label <- factor(
  unname(program_labels[score_models$program_name]),
  levels = unname(program_labels[rev(program_levels)])
)
score_models$context_label <- factor(
  unname(context_labels[score_models$context_id]),
  levels = unname(context_labels[context_ids])
)
score_models$logFC <- as.numeric(score_models$logFC)
score_models$adj.P.Val <- as.numeric(score_models$adj.P.Val)
score_models$FDR_significant <- as.logical(score_models$FDR_significant)
score_models$neg_log10_FDR <- -log10(pmax(score_models$adj.P.Val, 1e-300))
score_models$display_label <- paste0(
  formatC(score_models$logFC, format = "f", digits = 2),
  ifelse(score_models$FDR_significant, "*", "")
)

coverage$coverage_label <- factor(coverage$measured_context_n, levels = 0:4)

panel_D_data$context_label <- factor(
  unname(context_labels[panel_D_data$context_id]),
  levels = unname(context_labels[context_ids])
)
ordered_genes <- selected_examples$feature_id[order(selected_examples$selection_rank)]
panel_D_data$feature_id <- factor(panel_D_data$feature_id, levels = rev(ordered_genes))
panel_D_data$logFC <- as.numeric(panel_D_data$logFC)
panel_D_data$FDR_significant <- as.logical(panel_D_data$FDR_significant)
if (!"display_label" %in% names(panel_D_data)) {
  panel_D_data$display_label <- paste0(
    formatC(panel_D_data$logFC, format = "f", digits = 2),
    ifelse(panel_D_data$FDR_significant, "*", "")
  )
}

# ----------------------------------------------------------------------
# 3. Locked article typography + documented Figure 3 height exception
# ----------------------------------------------------------------------

header_tag_size_pt <- 11.0
header_title_size_pt <- 9.2
header_subtitle_size_pt <- 7.2
header_tag_to_title_gap_mm <- 1.6
header_title_row_height_mm <- 5.4
header_subtitle_row_height_mm <- 4.0
figure_width_mm <- 180
figure_height_mm <- 184

color_negative <- "#3F6FA8"
color_neutral <- "#F5F5F2"
color_positive <- "#B65B5B"
color_bar <- "#6C8198"
color_text <- "#161616"
color_muted <- "#505050"
color_grid <- "#D7D7D7"

program_effect_limit <- max(
  0.5,
  ceiling(max(abs(score_models$logFC), na.rm = TRUE) * 10) / 10
)
gene_effect_limit <- max(
  1.0,
  ceiling(max(abs(panel_D_data$logFC), na.rm = TRUE) * 10) / 10
)

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = color_grid),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(family = figure_font, size = 8.0, colour = color_text),
    axis.text = ggplot2::element_text(family = figure_font, size = 6.8, colour = color_text),
    legend.title = ggplot2::element_text(family = figure_font, size = 7.0, face = "bold", colour = color_text),
    legend.text = ggplot2::element_text(family = figure_font, size = 6.8, colour = color_text),
    plot.margin = ggplot2::margin(3, 4, 3, 4)
  )

# ----------------------------------------------------------------------
# 4. Panels: identical science/statistics to v1.2
# ----------------------------------------------------------------------

panel_A <- ggplot2::ggplot(
  score_models,
  ggplot2::aes(x = context_label, y = program_label, fill = logFC)
) +
  ggplot2::geom_tile(linewidth = 0.45, colour = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = display_label),
    family = figure_font, size = 2.45, colour = color_text
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-program_effect_limit, program_effect_limit),
    oob = scales::squish,
    name = "Program logFC"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::coord_cartesian(clip = "off") +
  base_theme +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.6, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.94,
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "none"
  )

panel_B <- ggplot2::ggplot(
  score_models,
  ggplot2::aes(x = context_label, y = program_label)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = neg_log10_FDR, fill = logFC),
    shape = 21, colour = "#9B9B9B", stroke = 0.35
  ) +
  ggplot2::geom_point(
    data = score_models[score_models$FDR_significant, , drop = FALSE],
    ggplot2::aes(size = neg_log10_FDR, fill = logFC),
    shape = 21, colour = "#171717", stroke = 0.9
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-program_effect_limit, program_effect_limit),
    oob = scales::squish,
    name = "Program logFC"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.1, 6.7), breaks = c(1, 2, 4),
    name = expression(-log[10](FDR))
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.35, colour = "#E1E1E1"),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.6, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.94,
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.key.width = grid::unit(7, "mm"),
    legend.margin = ggplot2::margin(t = 1, b = 0)
  ) +
  ggplot2::guides(
    fill = "none",
    size = ggplot2::guide_legend(
      title.position = "left", nrow = 1, byrow = TRUE,
      override.aes = list(fill = "#D6D6D6", colour = "#171717")
    )
  )

panel_C <- ggplot2::ggplot(
  coverage,
  ggplot2::aes(x = coverage_label, y = candidate_n)
) +
  ggplot2::geom_col(width = 0.68, fill = color_bar) +
  ggplot2::geom_text(
    ggplot2::aes(label = format(candidate_n, big.mark = ",", scientific = FALSE, trim = TRUE)),
    family = figure_font, size = 2.55, vjust = -0.45, colour = color_text
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.14)),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::labs(
    x = "Number of patient bulk contexts measured",
    y = "Candidate genes"
  ) +
  base_theme +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(family = figure_font, size = 7.0)
  )

panel_D <- ggplot2::ggplot(
  panel_D_data,
  ggplot2::aes(x = context_label, y = feature_id, fill = logFC)
) +
  ggplot2::geom_tile(linewidth = 0.45, colour = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = display_label),
    family = figure_font, size = 2.25, colour = color_text
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-gene_effect_limit, gene_effect_limit),
    oob = scales::squish,
    name = "Gene logFC"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.7, face = "italic",
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "none"
  )


# ----------------------------------------------------------------------
# 4B. Data-body-constrained horizontal colorbar helpers
# ----------------------------------------------------------------------

# Recreate guide-bearing reference plots ONLY for color scale / height extraction.
# They are never directly rendered into Figure 3.
panel_A_guide_reference <- panel_A +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(12, "mm"),
    legend.margin = ggplot2::margin(t = 1, b = 0)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "left", title.hjust = 0.5,
      barheight = grid::unit(2.6, "mm")
    )
  )

panel_D_guide_reference <- panel_D +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(12, "mm")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "left",
      barheight = grid::unit(2.6, "mm")
    )
  )

find_colorbar_guide_gtable <- function(g) {
  if (inherits(g, "gtable") &&
      any(g$layout$name == "bar") &&
      any(g$layout$name == "title")) {
    return(g)
  }
  if (!is.null(g$grobs)) {
    for (child in g$grobs) {
      ans <- find_colorbar_guide_gtable(child)
      if (!is.null(ans)) return(ans)
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

  tile_layer_idx <- which(vapply(
    built$data,
    function(d) all(c("xmin", "xmax") %in% names(d)),
    logical(1)
  ))

  if (!length(tile_layer_idx)) {
    stop("Could not identify geom_tile layer for Panel ", panel_tag)
  }

  tile_dat <- built$data[[tile_layer_idx[[1]]]]
  tile_right <- max(tile_dat$xmax, na.rm = TRUE)
  tile_left <- min(tile_dat$xmin, na.rm = TRUE)

  right_frac <- (tile_right - x_range[[1]]) /
    (x_range[[2]] - x_range[[1]])
  left_frac <- (tile_left - x_range[[1]]) /
    (x_range[[2]] - x_range[[1]])

  if (!is.finite(right_frac) || right_frac <= 0 || right_frac >= 1 ||
      !is.finite(left_frac) || left_frac <= 0 || left_frac >= right_frac) {
    stop("Invalid tile-body fractions for Panel ", panel_tag)
  }

  data.frame(
    panel = panel_tag,
    x_range_min = x_range[[1]],
    x_range_max = x_range[[2]],
    tile_left = tile_left,
    tile_right = tile_right,
    tile_left_fraction = left_frac,
    tile_right_fraction = right_frac,
    stringsAsFactors = FALSE
  )
}

measure_color_scale <- function(
  plot_object,
  panel_tag,
  approved_bar_width_mm
) {
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

  bar_height_mm <- grid::convertHeight(
    sum(guide_gt$heights[bar_cell$t:bar_cell$b]),
    "mm", valueOnly = TRUE
  )
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
  breaks <- breaks[
    is.finite(breaks) &
      breaks >= min(limits) &
      breaks <= max(limits)
  ]
  labels <- fill_scale$get_labels(breaks)

  raster_values <- seq(limits[[1]], limits[[2]], length.out = 512L)
  raster_colors <- fill_scale$map(raster_values)

  list(
    panel = panel_tag,
    approved_bar_width_mm = approved_bar_width_mm,
    bar_height_mm = bar_height_mm,
    limits = limits,
    breaks = breaks,
    labels = labels,
    raster_colors = raster_colors
  )
}

make_data_body_aligned_colorbar_grob <- function(
  title,
  metrics,
  tile_right_fraction,
  figure_font,
  title_size_pt = 7.0,
  label_size_pt = 6.8,
  color_text = "#161616",
  title_to_bar_gap_mm = 2.0
) {
  bar_w <- grid::unit(metrics$approved_bar_width_mm, "mm")
  bar_h <- grid::unit(metrics$bar_height_mm, "mm")

  # Vertical anchor:
  # text bottom == colorbar body bottom.
  bar_bottom <- grid::unit(4.1, "mm")
  bar_top <- bar_bottom + bar_h

  # Horizontal anchor:
  # bar right == LAST TILE / DATA BODY right boundary.
  bar_right <- grid::unit(tile_right_fraction, "npc")
  bar_left <- bar_right - bar_w

  gradient_matrix <- matrix(metrics$raster_colors, nrow = 1L)

  bar <- grid::rasterGrob(
    image = gradient_matrix,
    x = bar_left + 0.5 * bar_w,
    y = bar_bottom + 0.5 * bar_h,
    width = bar_w,
    height = bar_h,
    interpolate = TRUE
  )

  frac <- (metrics$breaks - metrics$limits[[1]]) /
    (metrics$limits[[2]] - metrics$limits[[1]])
  tick_x <- bar_left + frac * bar_w

  tick_grobs <- lapply(seq_along(tick_x), function(i) {
    grid::segmentsGrob(
      x0 = tick_x[[i]], x1 = tick_x[[i]],
      y0 = bar_bottom,
      y1 = bar_top,
      gp = grid::gpar(col = "white", lwd = 0.45)
    )
  })

  # Endpoint labels are inward-justified:
  # first break extends rightward from left edge;
  # last break extends leftward from right edge.
  label_grobs <- lapply(seq_along(tick_x), function(i) {
    if (i == 1L) {
      just_x <- "left"
    } else if (i == length(tick_x)) {
      just_x <- "right"
    } else {
      just_x <- "center"
    }

    grid::textGrob(
      as.character(metrics$labels[[i]]),
      x = tick_x[[i]],
      y = grid::unit(1.65, "mm"),
      just = c(just_x, "center"),
      gp = grid::gpar(
        fontfamily = figure_font,
        fontsize = label_size_pt,
        col = color_text
      )
    )
  })

  title_grob <- grid::textGrob(
    title,
    x = bar_left - grid::unit(title_to_bar_gap_mm, "mm"),
    y = bar_bottom,
    just = c("right", "bottom"),
    gp = grid::gpar(
      fontfamily = figure_font,
      fontface = "bold",
      fontsize = title_size_pt,
      col = color_text
    )
  )

  grid::grobTree(
    title_grob,
    bar,
    do.call(grid::grobTree, tick_grobs),
    do.call(grid::grobTree, label_grobs)
  )
}

add_data_body_colorbar_to_body <- function(
  body_gtable,
  legend_grob,
  legend_height_mm,
  panel_tag
) {
  panel_cells <- body_gtable$layout[
    grepl("^panel", body_gtable$layout$name), , drop = FALSE
  ]
  if (!nrow(panel_cells)) {
    stop("Could not locate scientific panel cell for Panel ", panel_tag)
  }

  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)

  body_gtable <- gtable::gtable_add_rows(
    body_gtable,
    heights = grid::unit(legend_height_mm, "mm"),
    pos = -1
  )
  new_row <- length(body_gtable$heights)

  # The custom legend's 0..1 npc coordinates now map ONLY to the
  # scientific panel column(s), not to axis/title/whitespace columns.
  body_gtable <- gtable::gtable_add_grob(
    body_gtable,
    legend_grob,
    t = new_row, b = new_row,
    l = panel_left, r = panel_right,
    clip = "off",
    name = paste0("F3v148-data-body-colorbar-", panel_tag)
  )

  body_gtable
}

# Frozen bar-body widths measured from the previously approved v1.3
# actual 300-dpi render, not guessed:
#   Panel A: 601 px at 300 dpi = 50.88 mm
#   Panel D: 708 px at 300 dpi = 59.94 mm
A_approved_bar_width_mm <- 50.88
D_approved_bar_width_mm <- 59.94

A_tile_geometry <- measure_data_body_right_fraction(panel_A, "A")
D_tile_geometry <- measure_data_body_right_fraction(panel_D, "D")

A_colorbar_metrics <- measure_color_scale(
  panel_A_guide_reference,
  "A",
  approved_bar_width_mm = A_approved_bar_width_mm
)
D_colorbar_metrics <- measure_color_scale(
  panel_D_guide_reference,
  "D",
  approved_bar_width_mm = D_approved_bar_width_mm
)

A_constraint_legend <- make_data_body_aligned_colorbar_grob(
  title = "Program logFC",
  metrics = A_colorbar_metrics,
  tile_right_fraction = A_tile_geometry$tile_right_fraction[[1]],
  figure_font = figure_font,
  color_text = color_text
)
D_constraint_legend <- make_data_body_aligned_colorbar_grob(
  title = "Gene logFC",
  metrics = D_colorbar_metrics,
  tile_right_fraction = D_tile_geometry$tile_right_fraction[[1]],
  figure_font = figure_font,
  color_text = color_text
)

# 5. v1.3 bottom-canvas repair: shared header rows, body rows separate
# ----------------------------------------------------------------------

build_body_gtable <- function(plot_object) {
  ggplot2::ggplotGrob(
    plot_object +
      ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) +
      ggplot2::theme(
        plot.tag = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank()
      )
  )
}

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle) {
  panel_cells <- body_gtable$layout[
    grepl("^panel", body_gtable$layout$name), , drop = FALSE
  ]
  if (!nrow(panel_cells)) stop("Could not locate scientific panel cell for ", panel_tag)
  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)

  header_gt <- gtable::gtable(
    widths = body_gtable$widths,
    heights = grid::unit(
      c(header_title_row_height_mm, header_subtitle_row_height_mm),
      "mm"
    )
  )

  title_grob <- grid::textGrob(
    panel_title,
    x = grid::unit(0, "npc"), y = grid::unit(0.52, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_title_size_pt, lineheight = 0.94,
      col = color_text
    )
  )
  tag_grob <- grid::textGrob(
    panel_tag,
    x = grid::unit(0, "npc") - grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"),
    hjust = 1, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_tag_size_pt, col = color_text
    )
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle,
    x = grid::unit(0, "npc"), y = grid::unit(0.56, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font,
      fontsize = header_subtitle_size_pt,
      col = color_muted
    )
  )

  header_gt <- gtable::gtable_add_grob(
    header_gt, grid::grobTree(tag_grob, title_grob),
    t = 1, b = 1, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F3v147-title-", panel_tag)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, subtitle_grob,
    t = 2, b = 2, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F3v147-subtitle-", panel_tag)
  )

  list(
    gtable = header_gt,
    geometry = data.frame(
      panel = panel_tag,
      anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE",
      panel_left_column = panel_left,
      panel_right_column = panel_right,
      tag_size_pt = header_tag_size_pt,
      title_size_pt = header_title_size_pt,
      subtitle_size_pt = header_subtitle_size_pt,
      tag_title_gap_mm = header_tag_to_title_gap_mm,
      header_composition = "DEDICATED_SHARED_ROW",
      stringsAsFactors = FALSE
    )
  )
}

body_A <- build_body_gtable(panel_A)
body_B <- build_body_gtable(panel_B)
body_C <- build_body_gtable(panel_C)
body_D <- build_body_gtable(panel_D)

body_A <- add_data_body_colorbar_to_body(
  body_A,
  A_constraint_legend,
  legend_height_mm = 8.0,
  panel_tag = "A"
)
body_D <- add_data_body_colorbar_to_body(
  body_D,
  D_constraint_legend,
  legend_height_mm = 8.0,
  panel_tag = "D"
)

header_A <- make_header_gtable(
  body_A, "A", "Program effects across patient contrasts",
  "All contrasts vs Control · cell = logFC; * FDR < 0.05"
)
header_B <- make_header_gtable(
  body_B, "B", "Effect magnitude and multiplicity",
  "vs Control · fill = logFC; size = -log10(FDR)"
)
header_C <- make_header_gtable(
  body_C, "C", "Candidate measurement coverage",
  "Patient contexts with an estimable gene-level effect"
)
header_D <- make_header_gtable(
  body_D, "D", "Context-dependent gene examples",
  "Top 10 mixed-direction genes · * FDR < 0.05"
)

header_geometry <- do.call(
  rbind,
  lapply(list(header_A, header_B, header_C, header_D), `[[`, "geometry")
)

design <- "
AB
CD
##
EF
GH
##
"

figure_object <- (
  patchwork::wrap_elements(full = header_A$gtable) +
  patchwork::wrap_elements(full = header_B$gtable) +
  patchwork::wrap_elements(full = body_A) +
  patchwork::wrap_elements(full = body_B) +
  patchwork::wrap_elements(full = header_C$gtable) +
  patchwork::wrap_elements(full = header_D$gtable) +
  patchwork::wrap_elements(full = body_C) +
  patchwork::wrap_elements(full = body_D)
) +
  patchwork::plot_layout(
    design = design,
    widths = c(1, 1),
    heights = grid::unit(c(10, 72, 4, 10, 70, 8), "mm")
  ) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm")
    )
  )


# ----------------------------------------------------------------------
# 5B. Data-body alignment geometry QC
# ----------------------------------------------------------------------

constraint_geometry_qc <- rbind(
  data.frame(
    panel = "A",
    approved_bar_width_mm = A_approved_bar_width_mm,
    tile_left_fraction = A_tile_geometry$tile_left_fraction[[1]],
    tile_right_fraction = A_tile_geometry$tile_right_fraction[[1]],
    rule_right_edge = "COLORBAR_BODY_RIGHT == LAST_TILE_RIGHT_EDGE",
    rule_bottom_edge = "LEGEND_TITLE_TEXT_BOTTOM == COLORBAR_BODY_BOTTOM",
    endpoint_label_rule = "FIRST_LEFT_JUSTIFIED_LAST_RIGHT_JUSTIFIED",
    stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "D",
    approved_bar_width_mm = D_approved_bar_width_mm,
    tile_left_fraction = D_tile_geometry$tile_left_fraction[[1]],
    tile_right_fraction = D_tile_geometry$tile_right_fraction[[1]],
    rule_right_edge = "COLORBAR_BODY_RIGHT == LAST_TILE_RIGHT_EDGE",
    rule_bottom_edge = "LEGEND_TITLE_TEXT_BOTTOM == COLORBAR_BODY_BOTTOM",
    endpoint_label_rule = "FIRST_LEFT_JUSTIFIED_LAST_RIGHT_JUSTIFIED",
    stringsAsFactors = FALSE
  )
)

constraint_geometry_qc$alignment_target <-
  "ACTUAL_DATA_TILE_BODY_NOT_OUTER_PANEL"
constraint_geometry_qc$alignment_method <-
  "CODE_LAYER_DATA_COORDINATE_CONSTRAINT"
constraint_geometry_qc$pixel_shift_guessing <- "FORBIDDEN"

write_csv_utf8(
  constraint_geometry_qc,
  file.path(qc_root, "FIG3_DATA_BODY_COLORBAR_CONSTRAINTS_v1_4_8.csv")
)

skill3_rule_update <- data.frame(
  rule_id = c(
    "SKILL3-DATA-BODY-RIGHT-EDGE-CONSTRAINT-001",
    "SKILL3-COLORBAR-WIDTH-PRESERVATION-001",
    "SKILL3-LEGEND-TITLE-BOTTOM-EDGE-CONSTRAINT-001",
    "SKILL3-COLORBAR-ENDPOINT-LABEL-SAFETY-001"
  ),
  rule = c(
    paste(
      "For heatmap/tile panels, the right-edge alignment target is the actual",
      "data-body boundary formed by the last tile edge, not the outer panel,",
      "gtable, axis, tick, or whitespace boundary."
    ),
    paste(
      "When polishing an already approved figure, preserve the approved actual-render",
      "colorbar-body width unless a deliberate redesign is authorized; do not allow",
      "a generic legend container width to stretch the bar."
    ),
    paste(
      "For a left-side horizontal-colorbar title, constrain the title text bottom",
      "edge to the colorbar-body bottom edge using code/grob geometry."
    ),
    paste(
      "Colorbar endpoint labels must remain fully visible. The first label may be",
      "left-justified and the last right-justified at their endpoint ticks to avoid",
      "clipping without moving the bar away from the data-body boundary."
    )
  ),
  status = "CARRY_FORWARD_TO_NEXT_SKILL3_REVISION",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  skill3_rule_update,
  file.path(qc_root, "SKILL3_RULE_UPDATE_DATA_BODY_COLORBAR_v1_4_8.csv")
)

# 6. BEFORE LOCK: one QC PNG only
# ----------------------------------------------------------------------

figure_png <- file.path(fig_root, "Figure3_v1_4_8_QC_180mm_300dpi.png")
ragg::agg_png(
  filename = figure_png,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = 300,
  background = "white"
)
print(figure_object)
grDevices::dev.off()

if (!file.exists(figure_png) || file.info(figure_png)$size < 50000) {
  stop("QC PNG export failed or is unexpectedly small.")
}

# ----------------------------------------------------------------------
# 7. QC evidence / frozen revision scope
# ----------------------------------------------------------------------

revision_scope <- data.frame(
  issue_id = "F3V12-01",
  severity = "BLOCKING",
  v1_2_finding = paste(
    "Actual 300-dpi PNG showed bottom-edge clipping:",
    "Panel C x-axis title and Panel D colorbar tick labels extended beyond the visible canvas."
  ),
  v1_3_action = paste(
    "Reserve an explicit 8-mm empty bottom row within the existing 184-mm canvas;",
    "rebalance body rows from 75/75 mm to 72/70 mm;",
    "retain all data, scales, labels, colors and typography sizes."
  ),
  allowed_change_scope = "LAYOUT_ONLY_BOTTOM_CANVAS_RESERVE",
  scientific_change = "NONE",
  statistical_change = "NONE",
  stringsAsFactors = FALSE
)
write_csv_utf8(revision_scope, file.path(qc_root, "FIG3_V1_2_TO_V1_3_REVISION_SCOPE_LOCK.csv"))

layout_exception <- data.frame(
  component = "figure_height",
  article_profile_value = "178 mm",
  figure3_candidate_value = "184 mm",
  rationale = paste(
    "Figure 3 is a dense 2x2 composite. The 184-mm candidate preserves the locked",
    "11.0/9.2/7.2-pt header system and final-size readability without shrinking text."
  ),
  status = "PENDING_USER_APPROVAL",
  scientific_effect = "NONE",
  statistical_effect = "NONE",
  stringsAsFactors = FALSE
)
write_csv_utf8(layout_exception, file.path(qc_root, "FIG3_LAYOUT_EXCEPTION_v1_3.csv"))

font_preflight <- data.frame(
  requested_font_family = figure_font,
  resolved_font_family = font_info$family[[1]],
  resolved_font_path = normalize_slashes(font_match$path[[1]]),
  exact_family_match = font_ok,
  skill_version = "2.2-beta.6",
  stringsAsFactors = FALSE
)
write_csv_utf8(font_preflight, file.path(qc_root, "FIG3_FONT_PREFLIGHT_v1_3.csv"))
write_csv_utf8(header_geometry, file.path(qc_root, "FIG3_PANEL_HEADER_GEOMETRY_v1_3.csv"))

source_manifest <- data.frame(
  file = basename(source_paths),
  source_version = source_version,
  sha256 = vapply(source_paths, sha256_file, character(1)),
  bytes = as.numeric(file.info(source_paths)$size),
  stringsAsFactors = FALSE
)
write_csv_utf8(source_manifest, file.path(qc_root, "FIG3_SOURCE_HASH_MANIFEST_v1_3.csv"))

candidate_status <- data.frame(
  figure_id = "Figure 3",
  figure_version = figure_version,
  status = "HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW",
  candidate_class = "QC_PNG_ONLY_BEFORE_LOCK",
  governing_skill = governing_skill,
  source_lock = source_lock,
  bulk_projection_version = bulk_projection_version,
  scientific_QC = "PASS_INHERITED_LOCKED_SOURCE",
  statistical_QC = "PASS_NO_RECOMPUTATION",
  actual_render_QC = "PENDING_USER_AND_CHAT_REVIEW",
  cross_format_QC = "NOT_RUN_BEFORE_NATIVE_RENDER_PASS",
  manuscript_ready = FALSE,
  immutable_release = FALSE,
  stringsAsFactors = FALSE
)
write_csv_utf8(candidate_status, file.path(qc_root, "FIGURE_STATUS_v1_3_QC.csv"))

manual_qc <- data.frame(
  check_id = sprintf("FIG3V13-%02d", 1:12),
  check = c(
    "A and B panel tags/titles share one top baseline",
    "C and D panel tags/titles share one top baseline",
    "Panel B subtitle is fully visible and not clipped",
    "All titles start at each panel scientific-body left edge",
    "Tag-to-title gap is visually identical across A/B/C/D",
    "Panel A heatmap values and asterisks are legible",
    "Panel B size differences and significant outlines are legible",
    "Panel C labels 69, 1, 0, 62 and 304 are all visible",
    "Panel D gene names, four contexts, values and asterisks are legible",
    "Blue-neutral-red direction semantics remain consistent in A/B/D",
    "No panel implies pooled replication, causality or definitive direct targets",
    "Panel C x-axis title and Panel D colorbar ticks are fully visible above the bottom canvas edge"
  ),
  status = "PENDING_ACTUAL_RENDER_REVIEW",
  stringsAsFactors = FALSE
)
write_csv_utf8(manual_qc, file.path(qc_root, "FIG3_STRICT_RENDER_QC_v1_3_PENDING.csv"))

write_text(
  c(
    "# Figure 3 v1.4.8 QC handoff",
    "",
    paste0("- Governing skill: ", governing_skill),
    paste0("- Source lock: ", source_lock),
    paste0("- Bulk projection: ", bulk_projection_version),
    "- Statistics recomputed: NO",
    "- Export before visual lock: ONE 300-dpi QC PNG only",
    "- Current status: HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW",
    "- Cross-format export/verification: deferred until native render passes and the user approves the candidate",
    "- Skill 8 immutable lock: NOT AUTHORIZED YET",
    "",
    "Frozen v1.3 repair scope:",
    "1. Correct only the bottom-canvas clipping seen in the v1.2 actual render.",
    "2. Add an explicit 8-mm empty bottom row and rebalance body-row heights within the same 184-mm canvas.",
    "3. Preserve all scientific/statistical content, scales, values, labels, colors and typography sizes.",
    "4. 184-mm Figure-3-specific height exception remains PENDING_USER_APPROVAL."
  ),
  file.path(qc_root, "FIG3_V1_3_QC_HANDOFF.md")
)

# Copy unchanged manuscript-facing text for review context.
if (dir.exists(text_root_v11)) {
  txt_files <- list.files(text_root_v11, full.names = TRUE)
  for (p in txt_files[file.exists(txt_files) & !dir.exists(txt_files)]) {
    copy_exact(p, file.path(text_copy_root, basename(p)))
  }
}
copy_exact(script_path, file.path(code_root, script_name))
writeLines(capture.output(sessionInfo()), file.path(qc_root, "sessionInfo_FIG3_v1_4_8.txt"))

write_text(
  c(
    "# Figure 3 v1.4.8 data-body-aligned legend handoff",
    "",
    "- Statistical recomputation: NO.",
    "- Raster patching: NO.",
    "- Pixel-shift guessing: NO.",
    "- Panel A/D default guides removed from rendered figure.",
    "- Fill scales/heights are measured from reference ggplots; bar widths are preserved from approved v1.3 actual render.",
    "- Constraint 1: colorbar body right edge == actual last-tile/data-body right edge.",
    "- Constraint 2: legend-title text bottom edge == colorbar body bottom edge.
    - Endpoint labels use inward justification to prevent clipping.",
    "- Next gate: strict actual-render QC, followed by final user visual confirmation."
  ),
  file.path(qc_root, "FIG3_V1_4_8_DATA_BODY_ALIGNED_HANDOFF.md")
)

# ----------------------------------------------------------------------
# 8. Desktop SMALL_REPORT ZIP for upload
# ----------------------------------------------------------------------

desktop_root <- find_desktop()
report_root <- file.path(desktop_root, paste0("FIG3_V1_4_8_QC_REPORT_", run_stamp))
zip_path <- file.path(desktop_root, paste0("STEP05C_FIG3_V1_4_8_QC_", run_stamp, "_SMALL_REPORT.zip"))
if (dir.exists(report_root) || file.exists(zip_path)) stop("Desktop report target already exists.")

report_fig3 <- file.path(report_root, "FIGURE3")
for (subdir in c("figures", "qc", "source_data", "code", "manuscript_text")) {
  dir.create(file.path(report_fig3, subdir), recursive = TRUE, showWarnings = FALSE)
}

copy_exact(figure_png, file.path(report_fig3, "figures", basename(figure_png)))
for (p in list.files(qc_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "qc", basename(p)))
}
for (p in list.files(source_copy_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "source_data", basename(p)))
}
for (p in list.files(code_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "code", basename(p)))
}
for (p in list.files(text_copy_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "manuscript_text", basename(p)))
}

report_files <- list.files(report_root, recursive = TRUE, full.names = TRUE)
report_files <- report_files[file.exists(report_files) & !dir.exists(report_files)]
report_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(report_files),
    nchar(normalize_slashes(report_root)) + 2L
  ),
  bytes = as.numeric(file.info(report_files)$size),
  sha256 = vapply(report_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))

zip::zipr(
  zipfile = zip_path,
  files = basename(report_root),
  root = dirname(report_root),
  include_directories = TRUE
)
if (!file.exists(zip_path) || file.info(zip_path)$size < 50000) stop("Desktop ZIP creation failed.")

cat("\n============================================================\n")
cat("Figure 3 v1.4.8 QC candidate generated.\n")
cat("Gate: HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW\n")
cat("Scientific/statistical recomputation: NONE\n")
cat("QC PNG:\n", normalize_slashes(figure_png), "\n", sep = "")
cat("Upload this ZIP:\n", normalize_slashes(zip_path), "\n", sep = "")
cat("Do NOT run Skill 8 lock yet.\n")
cat("============================================================\n")

# Publication figure stage 2
# ======================================================================
# 113_step05C_Figure3_submission_formats_600ppi_pdf_v1_4_8_QC3.R
# Nutrition_Endometrial_VDR_v2
#
# Skill 3 v2.2-beta.6 SUBMISSION-FORMAT cross-format candidate.
# Reads ONLY frozen Figure 3 source CSVs; no statistics are recomputed.
# Produces true 600-ppi PNG + true 600-ppi TIFF + vector PDF at 180 x 184 mm,
# plus the approved 300-ppi PNG as a QC reference and a SMALL_REPORT ZIP.
#
# v1.4.8 RULE-DRIVEN CHANGE ONLY:
# 1) Panels A/D use custom code-layer horizontal colorbars, not raster patches.
# 2) IMPORTANT: alignment target is the ACTUAL DATA/TILE BODY right boundary,
#    not the outer panel/gtable/axis/whitespace boundary.
# 3) Colorbar BODY right edge is constrained to the right edge of the last tile.
# 4) Legend-title TEXT bottom edge is constrained to the colorbar BODY bottom edge.
# 5) Bar-body widths are preserved from the previously approved v1.3 actual render:
#      Panel A = 50.88 mm; Panel D = 59.94 mm.
#    This prevents Panel A from being stretched to the generic 60-mm width.
# 6) Endpoint tick labels are inward-justified so the maximum label cannot clip.
# 7) No guessed visual pixel shifts are used for placement.

#
# Frozen repair from v1.2 actual-render QC:
# 1) BLOCKING: the lower-row canvas ended too close to the device boundary.
#    Panel C x-axis title and Panel D colorbar tick labels were clipped at the
#    bottom edge in the actual 300-dpi PNG.
#    Repair: reserve a dedicated 8-mm empty bottom row inside the 184-mm canvas
#    and slightly rebalance body-row heights. No scientific/statistical content
#    or typography sizes are changed.
#
# SAVE AS:
# <repository module script>
#
# RUN:
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "113_step05C_Figure3_submission_formats_600ppi_pdf_v1_4_8_QC3.R"
script_path <- file.path(project_root, "04_code", "R", script_name)
figure_version <- "F3_v1_4_8_SUBMISSION_FORMATS_QC3"
source_version <- "F3_v1_2_QC"
source_lock <- "STEP04_LOCK_v1.2"
bulk_projection_version <- "BULK_VDR_PROGRAM_PROJECTION_v1.3"
governing_skill <- "research_scientific_figure_suite v2.2-beta.6"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

required_packages <- c(
  "ggplot2", "patchwork", "gtable", "grid", "scales",
  "systemfonts", "ragg", "digest", "zip"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

normalize_slashes <- function(x) gsub("\\\\", "/", x)
write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_text <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(x), path, useBytes = TRUE)
}
sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
copy_exact <- function(src, dst) {
  if (!file.exists(src)) stop("Missing source file: ", src)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dst) || sha256_file(src) != sha256_file(dst)) {
    stop("Validated copy failed: ", src, " -> ", dst)
  }
  invisible(TRUE)
}
find_desktop <- function() {
  candidates <- c(
    file.path(Sys.getenv("USERPROFILE"), "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "桌面")
  )
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) stop("Could not locate Desktop.")
  candidates[[1]]
}

# ----------------------------------------------------------------------
# 1. Governance / Skill preflight
# ----------------------------------------------------------------------

skill_path <- file.path(
  "<FIGURE_STYLE_REFERENCE>",
  "research_scientific_figure_suite",
  "SKILL.md"
)
if (!file.exists(skill_path)) stop("Active Skill 3 file is missing: ", skill_path)
skill_text <- paste(readLines(skill_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.6", skill_text, fixed = TRUE)) {
  stop("Active Skill 3 is not v2.2-beta.6.")
}
if (!file.exists(script_path)) {
  stop(
    "Please save this code first as:\n", normalize_slashes(script_path),
    "\nThen run it with source(...)."
  )
}

figure_font <- "Arial"
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}
font_match <- systemfonts::match_fonts(figure_font)
font_info <- systemfonts::font_info(font_match$path)
font_ok <- nrow(font_match) >= 1L && file.exists(font_match$path[[1]]) &&
  grepl("Arial", font_info$family[[1]], ignore.case = TRUE)
if (!font_ok) stop("Arial exact-family font preflight failed.")

# ----------------------------------------------------------------------
# 2. Read v1.1 figure-source tables only (NO statistical recomputation)
# ----------------------------------------------------------------------

step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")
fig3_parent <- file.path(step05_root, "06_results", "FIG3")
source_root <- file.path(fig3_parent, source_version, "source_data")
text_root_v11 <- file.path(fig3_parent, source_version, "manuscript_text")

out_root <- file.path(fig3_parent, figure_version)
fig_root <- file.path(out_root, "figures")
qc_root <- file.path(out_root, "qc")
code_root <- file.path(out_root, "code")
source_copy_root <- file.path(out_root, "source_data")
text_copy_root <- file.path(out_root, "manuscript_text")

if (dir.exists(out_root)) {
  stop(
    "Refusing to overwrite reviewed candidate folder:\n",
    normalize_slashes(out_root),
    "\nIf this is an incomplete failed run, inspect it manually before removing it."
  )
}
for (p in c(fig_root, qc_root, code_root, source_copy_root, text_copy_root)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

source_files_needed <- c(
  "Figure3_panelA_program_effect_heatmap.csv",
  "Figure3_panelB_FDR_effect_matrix.csv",
  "Figure3_panelC_candidate_context_coverage.csv",
  "Figure3_panelD_context_effects_long.csv",
  "Figure3_panelD_selected_examples.csv",
  "Figure3_locked_blueprint_row.csv",
  "Figure3_score_significance_summary.csv",
  "Figure3_candidate_direction_summary.csv"
)
source_paths <- file.path(source_root, source_files_needed)
if (!all(file.exists(source_paths))) {
  stop("One or more Figure 3 v1.2 QC source tables are missing.")
}

for (p in source_paths) {
  copy_exact(p, file.path(source_copy_root, basename(p)))
}

score_models <- read.csv(
  file.path(source_root, "Figure3_panelA_program_effect_heatmap.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
coverage <- read.csv(
  file.path(source_root, "Figure3_panelC_candidate_context_coverage.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
panel_D_data <- read.csv(
  file.path(source_root, "Figure3_panelD_context_effects_long.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)
selected_examples <- read.csv(
  file.path(source_root, "Figure3_panelD_selected_examples.csv"),
  stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"
)

# Locked structure checks only.
if (nrow(score_models) != 16L) stop("Expected 16 program-by-contrast rows.")
if (sum(as.logical(score_models$FDR_significant)) != 4L) stop("Expected 4 FDR-significant program rows.")
if (nrow(coverage) != 5L || sum(coverage$candidate_n) != 436L) stop("Coverage table mismatch.")
if (nrow(selected_examples) != 10L || nrow(panel_D_data) != 40L) stop("Panel D source mismatch.")

program_levels <- c(
  "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284",
  "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84"
)
program_labels <- c(
  "RNA_ALL_436" = "RNA program\n(n=436)",
  "OCCUPANCY_SUPPORTED_284" = "Occupancy-supported\n(n=284)",
  "SHARED_OCCUPANCY_46" = "Shared occupancy\n(n=46)",
  "PROMOTER_OCCUPANCY_84" = "Promoter-associated\n(n=84)"
)
context_ids <- c(
  "GSE165004__RPL_vs_CONTROL",
  "GSE165004__UIF_vs_CONTROL",
  "GSE92324__IF_vs_CONTROL",
  "GSE111974__RIF_vs_CONTROL"
)
context_labels <- c(
  "GSE165004__RPL_vs_CONTROL" = "GSE165004\nRPL",
  "GSE165004__UIF_vs_CONTROL" = "GSE165004\nUIF",
  "GSE92324__IF_vs_CONTROL" = "GSE92324\nIF",
  "GSE111974__RIF_vs_CONTROL" = "GSE111974\nRIF"
)

score_models$program_label <- factor(
  unname(program_labels[score_models$program_name]),
  levels = unname(program_labels[rev(program_levels)])
)
score_models$context_label <- factor(
  unname(context_labels[score_models$context_id]),
  levels = unname(context_labels[context_ids])
)
score_models$logFC <- as.numeric(score_models$logFC)
score_models$adj.P.Val <- as.numeric(score_models$adj.P.Val)
score_models$FDR_significant <- as.logical(score_models$FDR_significant)
score_models$neg_log10_FDR <- -log10(pmax(score_models$adj.P.Val, 1e-300))
score_models$display_label <- paste0(
  formatC(score_models$logFC, format = "f", digits = 2),
  ifelse(score_models$FDR_significant, "*", "")
)

coverage$coverage_label <- factor(coverage$measured_context_n, levels = 0:4)

panel_D_data$context_label <- factor(
  unname(context_labels[panel_D_data$context_id]),
  levels = unname(context_labels[context_ids])
)
ordered_genes <- selected_examples$feature_id[order(selected_examples$selection_rank)]
panel_D_data$feature_id <- factor(panel_D_data$feature_id, levels = rev(ordered_genes))
panel_D_data$logFC <- as.numeric(panel_D_data$logFC)
panel_D_data$FDR_significant <- as.logical(panel_D_data$FDR_significant)
if (!"display_label" %in% names(panel_D_data)) {
  panel_D_data$display_label <- paste0(
    formatC(panel_D_data$logFC, format = "f", digits = 2),
    ifelse(panel_D_data$FDR_significant, "*", "")
  )
}

# ----------------------------------------------------------------------
# 3. Locked article typography + documented Figure 3 height exception
# ----------------------------------------------------------------------

header_tag_size_pt <- 11.0
header_title_size_pt <- 9.2
header_subtitle_size_pt <- 7.2
header_tag_to_title_gap_mm <- 1.6
header_title_row_height_mm <- 5.4
header_subtitle_row_height_mm <- 4.0
figure_width_mm <- 180
figure_height_mm <- 184

color_negative <- "#3F6FA8"
color_neutral <- "#F5F5F2"
color_positive <- "#B65B5B"
color_bar <- "#6C8198"
color_text <- "#161616"
color_muted <- "#505050"
color_grid <- "#D7D7D7"

program_effect_limit <- max(
  0.5,
  ceiling(max(abs(score_models$logFC), na.rm = TRUE) * 10) / 10
)
gene_effect_limit <- max(
  1.0,
  ceiling(max(abs(panel_D_data$logFC), na.rm = TRUE) * 10) / 10
)

base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = color_grid),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(family = figure_font, size = 8.0, colour = color_text),
    axis.text = ggplot2::element_text(family = figure_font, size = 6.8, colour = color_text),
    legend.title = ggplot2::element_text(family = figure_font, size = 7.0, face = "bold", colour = color_text),
    legend.text = ggplot2::element_text(family = figure_font, size = 6.8, colour = color_text),
    plot.margin = ggplot2::margin(3, 4, 3, 4)
  )

# ----------------------------------------------------------------------
# 4. Panels: identical science/statistics to v1.2
# ----------------------------------------------------------------------

panel_A <- ggplot2::ggplot(
  score_models,
  ggplot2::aes(x = context_label, y = program_label, fill = logFC)
) +
  ggplot2::geom_tile(linewidth = 0.45, colour = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = display_label),
    family = figure_font, size = 2.45, colour = color_text
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-program_effect_limit, program_effect_limit),
    oob = scales::squish,
    name = "Program logFC"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::coord_cartesian(clip = "off") +
  base_theme +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.6, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.94,
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "none"
  )

panel_B <- ggplot2::ggplot(
  score_models,
  ggplot2::aes(x = context_label, y = program_label)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = neg_log10_FDR, fill = logFC),
    shape = 21, colour = "#9B9B9B", stroke = 0.35
  ) +
  ggplot2::geom_point(
    data = score_models[score_models$FDR_significant, , drop = FALSE],
    ggplot2::aes(size = neg_log10_FDR, fill = logFC),
    shape = 21, colour = "#171717", stroke = 0.9
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-program_effect_limit, program_effect_limit),
    oob = scales::squish,
    name = "Program logFC"
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.1, 6.7), breaks = c(1, 2, 4),
    name = expression(-log[10](FDR))
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.35, colour = "#E1E1E1"),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.6, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.94,
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.spacing.x = grid::unit(1.2, "mm"),
    legend.key.width = grid::unit(7, "mm"),
    legend.margin = ggplot2::margin(t = 1, b = 0)
  ) +
  ggplot2::guides(
    fill = "none",
    size = ggplot2::guide_legend(
      title.position = "left", nrow = 1, byrow = TRUE,
      override.aes = list(fill = "#D6D6D6", colour = "#171717")
    )
  )

panel_C <- ggplot2::ggplot(
  coverage,
  ggplot2::aes(x = coverage_label, y = candidate_n)
) +
  ggplot2::geom_col(width = 0.68, fill = color_bar) +
  ggplot2::geom_text(
    ggplot2::aes(label = format(candidate_n, big.mark = ",", scientific = FALSE, trim = TRUE)),
    family = figure_font, size = 2.55, vjust = -0.45, colour = color_text
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.14)),
    labels = scales::label_number(big.mark = ",")
  ) +
  ggplot2::labs(
    x = "Number of patient bulk contexts measured",
    y = "Candidate genes"
  ) +
  base_theme +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(family = figure_font, size = 7.0)
  )

panel_D <- ggplot2::ggplot(
  panel_D_data,
  ggplot2::aes(x = context_label, y = feature_id, fill = logFC)
) +
  ggplot2::geom_tile(linewidth = 0.45, colour = "white") +
  ggplot2::geom_text(
    ggplot2::aes(label = display_label),
    family = figure_font, size = 2.25, colour = color_text
  ) +
  ggplot2::scale_fill_gradient2(
    low = color_negative, mid = color_neutral, high = color_positive,
    midpoint = 0,
    limits = c(-gene_effect_limit, gene_effect_limit),
    oob = scales::squish,
    name = "Gene logFC"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  base_theme +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      family = figure_font, size = 6.5, lineheight = 0.92,
      margin = ggplot2::margin(t = 3)
    ),
    axis.text.y = ggplot2::element_text(
      family = figure_font, size = 6.7, face = "italic",
      margin = ggplot2::margin(r = 3)
    ),
    legend.position = "none"
  )


# ----------------------------------------------------------------------
# 4B. Data-body-constrained horizontal colorbar helpers
# ----------------------------------------------------------------------

# Recreate guide-bearing reference plots ONLY for color scale / height extraction.
# They are never directly rendered into Figure 3.
panel_A_guide_reference <- panel_A +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(12, "mm"),
    legend.margin = ggplot2::margin(t = 1, b = 0)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "left", title.hjust = 0.5,
      barheight = grid::unit(2.6, "mm")
    )
  )

panel_D_guide_reference <- panel_D +
  ggplot2::theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(12, "mm")
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "left",
      barheight = grid::unit(2.6, "mm")
    )
  )

find_colorbar_guide_gtable <- function(g) {
  if (inherits(g, "gtable") &&
      any(g$layout$name == "bar") &&
      any(g$layout$name == "title")) {
    return(g)
  }
  if (!is.null(g$grobs)) {
    for (child in g$grobs) {
      ans <- find_colorbar_guide_gtable(child)
      if (!is.null(ans)) return(ans)
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

  tile_layer_idx <- which(vapply(
    built$data,
    function(d) all(c("xmin", "xmax") %in% names(d)),
    logical(1)
  ))

  if (!length(tile_layer_idx)) {
    stop("Could not identify geom_tile layer for Panel ", panel_tag)
  }

  tile_dat <- built$data[[tile_layer_idx[[1]]]]
  tile_right <- max(tile_dat$xmax, na.rm = TRUE)
  tile_left <- min(tile_dat$xmin, na.rm = TRUE)

  right_frac <- (tile_right - x_range[[1]]) /
    (x_range[[2]] - x_range[[1]])
  left_frac <- (tile_left - x_range[[1]]) /
    (x_range[[2]] - x_range[[1]])

  if (!is.finite(right_frac) || right_frac <= 0 || right_frac >= 1 ||
      !is.finite(left_frac) || left_frac <= 0 || left_frac >= right_frac) {
    stop("Invalid tile-body fractions for Panel ", panel_tag)
  }

  data.frame(
    panel = panel_tag,
    x_range_min = x_range[[1]],
    x_range_max = x_range[[2]],
    tile_left = tile_left,
    tile_right = tile_right,
    tile_left_fraction = left_frac,
    tile_right_fraction = right_frac,
    stringsAsFactors = FALSE
  )
}

measure_color_scale <- function(
  plot_object,
  panel_tag,
  approved_bar_width_mm
) {
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

  bar_height_mm <- grid::convertHeight(
    sum(guide_gt$heights[bar_cell$t:bar_cell$b]),
    "mm", valueOnly = TRUE
  )
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
  breaks <- breaks[
    is.finite(breaks) &
      breaks >= min(limits) &
      breaks <= max(limits)
  ]
  labels <- fill_scale$get_labels(breaks)

  raster_values <- seq(limits[[1]], limits[[2]], length.out = 512L)
  raster_colors <- fill_scale$map(raster_values)

  list(
    panel = panel_tag,
    approved_bar_width_mm = approved_bar_width_mm,
    bar_height_mm = bar_height_mm,
    limits = limits,
    breaks = breaks,
    labels = labels,
    raster_colors = raster_colors
  )
}

make_data_body_aligned_colorbar_grob <- function(
  title,
  metrics,
  tile_right_fraction,
  figure_font,
  title_size_pt = 7.0,
  label_size_pt = 6.8,
  color_text = "#161616",
  title_to_bar_gap_mm = 2.0
) {
  bar_w <- grid::unit(metrics$approved_bar_width_mm, "mm")
  bar_h <- grid::unit(metrics$bar_height_mm, "mm")

  # Vertical anchor:
  # text bottom == colorbar body bottom.
  bar_bottom <- grid::unit(4.1, "mm")
  bar_top <- bar_bottom + bar_h

  # Horizontal anchor:
  # bar right == LAST TILE / DATA BODY right boundary.
  bar_right <- grid::unit(tile_right_fraction, "npc")
  bar_left <- bar_right - bar_w

  gradient_matrix <- matrix(metrics$raster_colors, nrow = 1L)

  bar <- grid::rasterGrob(
    image = gradient_matrix,
    x = bar_left + 0.5 * bar_w,
    y = bar_bottom + 0.5 * bar_h,
    width = bar_w,
    height = bar_h,
    interpolate = TRUE
  )

  frac <- (metrics$breaks - metrics$limits[[1]]) /
    (metrics$limits[[2]] - metrics$limits[[1]])
  tick_x <- bar_left + frac * bar_w

  tick_grobs <- lapply(seq_along(tick_x), function(i) {
    grid::segmentsGrob(
      x0 = tick_x[[i]], x1 = tick_x[[i]],
      y0 = bar_bottom,
      y1 = bar_top,
      gp = grid::gpar(col = "white", lwd = 0.45)
    )
  })

  # Endpoint labels are inward-justified:
  # first break extends rightward from left edge;
  # last break extends leftward from right edge.
  label_grobs <- lapply(seq_along(tick_x), function(i) {
    if (i == 1L) {
      just_x <- "left"
    } else if (i == length(tick_x)) {
      just_x <- "right"
    } else {
      just_x <- "center"
    }

    grid::textGrob(
      as.character(metrics$labels[[i]]),
      x = tick_x[[i]],
      y = grid::unit(1.65, "mm"),
      just = c(just_x, "center"),
      gp = grid::gpar(
        fontfamily = figure_font,
        fontsize = label_size_pt,
        col = color_text
      )
    )
  })

  title_grob <- grid::textGrob(
    title,
    x = bar_left - grid::unit(title_to_bar_gap_mm, "mm"),
    y = bar_bottom,
    just = c("right", "bottom"),
    gp = grid::gpar(
      fontfamily = figure_font,
      fontface = "bold",
      fontsize = title_size_pt,
      col = color_text
    )
  )

  grid::grobTree(
    title_grob,
    bar,
    do.call(grid::grobTree, tick_grobs),
    do.call(grid::grobTree, label_grobs)
  )
}

add_data_body_colorbar_to_body <- function(
  body_gtable,
  legend_grob,
  legend_height_mm,
  panel_tag
) {
  panel_cells <- body_gtable$layout[
    grepl("^panel", body_gtable$layout$name), , drop = FALSE
  ]
  if (!nrow(panel_cells)) {
    stop("Could not locate scientific panel cell for Panel ", panel_tag)
  }

  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)

  body_gtable <- gtable::gtable_add_rows(
    body_gtable,
    heights = grid::unit(legend_height_mm, "mm"),
    pos = -1
  )
  new_row <- length(body_gtable$heights)

  # The custom legend's 0..1 npc coordinates now map ONLY to the
  # scientific panel column(s), not to axis/title/whitespace columns.
  body_gtable <- gtable::gtable_add_grob(
    body_gtable,
    legend_grob,
    t = new_row, b = new_row,
    l = panel_left, r = panel_right,
    clip = "off",
    name = paste0("F3v148-data-body-colorbar-", panel_tag)
  )

  body_gtable
}

# Frozen bar-body widths measured from the previously approved v1.3
# actual 300-dpi render, not guessed:
#   Panel A: 601 px at 300 dpi = 50.88 mm
#   Panel D: 708 px at 300 dpi = 59.94 mm
A_approved_bar_width_mm <- 50.88
D_approved_bar_width_mm <- 59.94

A_tile_geometry <- measure_data_body_right_fraction(panel_A, "A")
D_tile_geometry <- measure_data_body_right_fraction(panel_D, "D")

A_colorbar_metrics <- measure_color_scale(
  panel_A_guide_reference,
  "A",
  approved_bar_width_mm = A_approved_bar_width_mm
)
D_colorbar_metrics <- measure_color_scale(
  panel_D_guide_reference,
  "D",
  approved_bar_width_mm = D_approved_bar_width_mm
)

A_constraint_legend <- make_data_body_aligned_colorbar_grob(
  title = "Program logFC",
  metrics = A_colorbar_metrics,
  tile_right_fraction = A_tile_geometry$tile_right_fraction[[1]],
  figure_font = figure_font,
  color_text = color_text
)
D_constraint_legend <- make_data_body_aligned_colorbar_grob(
  title = "Gene logFC",
  metrics = D_colorbar_metrics,
  tile_right_fraction = D_tile_geometry$tile_right_fraction[[1]],
  figure_font = figure_font,
  color_text = color_text
)

# 5. v1.3 bottom-canvas repair: shared header rows, body rows separate
# ----------------------------------------------------------------------

build_body_gtable <- function(plot_object) {
  ggplot2::ggplotGrob(
    plot_object +
      ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) +
      ggplot2::theme(
        plot.tag = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank()
      )
  )
}

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle) {
  panel_cells <- body_gtable$layout[
    grepl("^panel", body_gtable$layout$name), , drop = FALSE
  ]
  if (!nrow(panel_cells)) stop("Could not locate scientific panel cell for ", panel_tag)
  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)

  header_gt <- gtable::gtable(
    widths = body_gtable$widths,
    heights = grid::unit(
      c(header_title_row_height_mm, header_subtitle_row_height_mm),
      "mm"
    )
  )

  title_grob <- grid::textGrob(
    panel_title,
    x = grid::unit(0, "npc"), y = grid::unit(0.52, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_title_size_pt, lineheight = 0.94,
      col = color_text
    )
  )
  tag_grob <- grid::textGrob(
    panel_tag,
    x = grid::unit(0, "npc") - grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"),
    hjust = 1, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_tag_size_pt, col = color_text
    )
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle,
    x = grid::unit(0, "npc"), y = grid::unit(0.56, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font,
      fontsize = header_subtitle_size_pt,
      col = color_muted
    )
  )

  header_gt <- gtable::gtable_add_grob(
    header_gt, grid::grobTree(tag_grob, title_grob),
    t = 1, b = 1, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F3v147-title-", panel_tag)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, subtitle_grob,
    t = 2, b = 2, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F3v147-subtitle-", panel_tag)
  )

  list(
    gtable = header_gt,
    geometry = data.frame(
      panel = panel_tag,
      anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE",
      panel_left_column = panel_left,
      panel_right_column = panel_right,
      tag_size_pt = header_tag_size_pt,
      title_size_pt = header_title_size_pt,
      subtitle_size_pt = header_subtitle_size_pt,
      tag_title_gap_mm = header_tag_to_title_gap_mm,
      header_composition = "DEDICATED_SHARED_ROW",
      stringsAsFactors = FALSE
    )
  )
}

body_A <- build_body_gtable(panel_A)
body_B <- build_body_gtable(panel_B)
body_C <- build_body_gtable(panel_C)
body_D <- build_body_gtable(panel_D)

body_A <- add_data_body_colorbar_to_body(
  body_A,
  A_constraint_legend,
  legend_height_mm = 8.0,
  panel_tag = "A"
)
body_D <- add_data_body_colorbar_to_body(
  body_D,
  D_constraint_legend,
  legend_height_mm = 8.0,
  panel_tag = "D"
)

header_A <- make_header_gtable(
  body_A, "A", "Program effects across patient contrasts",
  "All contrasts vs Control · cell = logFC; * FDR < 0.05"
)
header_B <- make_header_gtable(
  body_B, "B", "Effect magnitude and multiplicity",
  "vs Control · fill = logFC; size = -log10(FDR)"
)
header_C <- make_header_gtable(
  body_C, "C", "Candidate measurement coverage",
  "Patient contexts with an estimable gene-level effect"
)
header_D <- make_header_gtable(
  body_D, "D", "Context-dependent gene examples",
  "Top 10 mixed-direction genes · * FDR < 0.05"
)

header_geometry <- do.call(
  rbind,
  lapply(list(header_A, header_B, header_C, header_D), `[[`, "geometry")
)

design <- "
AB
CD
##
EF
GH
##
"

figure_object <- (
  patchwork::wrap_elements(full = header_A$gtable) +
  patchwork::wrap_elements(full = header_B$gtable) +
  patchwork::wrap_elements(full = body_A) +
  patchwork::wrap_elements(full = body_B) +
  patchwork::wrap_elements(full = header_C$gtable) +
  patchwork::wrap_elements(full = header_D$gtable) +
  patchwork::wrap_elements(full = body_C) +
  patchwork::wrap_elements(full = body_D)
) +
  patchwork::plot_layout(
    design = design,
    widths = c(1, 1),
    heights = grid::unit(c(10, 72, 4, 10, 70, 8), "mm")
  ) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm")
    )
  )


# ----------------------------------------------------------------------
# 5B. Data-body alignment geometry QC
# ----------------------------------------------------------------------

constraint_geometry_qc <- rbind(
  data.frame(
    panel = "A",
    approved_bar_width_mm = A_approved_bar_width_mm,
    tile_left_fraction = A_tile_geometry$tile_left_fraction[[1]],
    tile_right_fraction = A_tile_geometry$tile_right_fraction[[1]],
    rule_right_edge = "COLORBAR_BODY_RIGHT == LAST_TILE_RIGHT_EDGE",
    rule_bottom_edge = "LEGEND_TITLE_TEXT_BOTTOM == COLORBAR_BODY_BOTTOM",
    endpoint_label_rule = "FIRST_LEFT_JUSTIFIED_LAST_RIGHT_JUSTIFIED",
    stringsAsFactors = FALSE
  ),
  data.frame(
    panel = "D",
    approved_bar_width_mm = D_approved_bar_width_mm,
    tile_left_fraction = D_tile_geometry$tile_left_fraction[[1]],
    tile_right_fraction = D_tile_geometry$tile_right_fraction[[1]],
    rule_right_edge = "COLORBAR_BODY_RIGHT == LAST_TILE_RIGHT_EDGE",
    rule_bottom_edge = "LEGEND_TITLE_TEXT_BOTTOM == COLORBAR_BODY_BOTTOM",
    endpoint_label_rule = "FIRST_LEFT_JUSTIFIED_LAST_RIGHT_JUSTIFIED",
    stringsAsFactors = FALSE
  )
)

constraint_geometry_qc$alignment_target <-
  "ACTUAL_DATA_TILE_BODY_NOT_OUTER_PANEL"
constraint_geometry_qc$alignment_method <-
  "CODE_LAYER_DATA_COORDINATE_CONSTRAINT"
constraint_geometry_qc$pixel_shift_guessing <- "FORBIDDEN"

write_csv_utf8(
  constraint_geometry_qc,
  file.path(qc_root, "FIG3_DATA_BODY_COLORBAR_CONSTRAINTS_v1_4_8.csv")
)

skill3_rule_update <- data.frame(
  rule_id = c(
    "SKILL3-DATA-BODY-RIGHT-EDGE-CONSTRAINT-001",
    "SKILL3-COLORBAR-WIDTH-PRESERVATION-001",
    "SKILL3-LEGEND-TITLE-BOTTOM-EDGE-CONSTRAINT-001",
    "SKILL3-COLORBAR-ENDPOINT-LABEL-SAFETY-001"
  ),
  rule = c(
    paste(
      "For heatmap/tile panels, the right-edge alignment target is the actual",
      "data-body boundary formed by the last tile edge, not the outer panel,",
      "gtable, axis, tick, or whitespace boundary."
    ),
    paste(
      "When polishing an already approved figure, preserve the approved actual-render",
      "colorbar-body width unless a deliberate redesign is authorized; do not allow",
      "a generic legend container width to stretch the bar."
    ),
    paste(
      "For a left-side horizontal-colorbar title, constrain the title text bottom",
      "edge to the colorbar-body bottom edge using code/grob geometry."
    ),
    paste(
      "Colorbar endpoint labels must remain fully visible. The first label may be",
      "left-justified and the last right-justified at their endpoint ticks to avoid",
      "clipping without moving the bar away from the data-body boundary."
    )
  ),
  status = "CARRY_FORWARD_TO_NEXT_SKILL3_REVISION",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  skill3_rule_update,
  file.path(qc_root, "SKILL3_RULE_UPDATE_DATA_BODY_COLORBAR_v1_4_8.csv")
)

# 6. BEFORE LOCK: submission-format export candidate
# ----------------------------------------------------------------------

# Submission geometry
submission_width_mm <- figure_width_mm
submission_height_mm <- figure_height_mm
submission_width_in <- submission_width_mm / 25.4
submission_height_in <- submission_height_mm / 25.4

# Expected pixel geometry at 600 ppi
expected_width_px_600 <- round(submission_width_in * 600)
expected_height_px_600 <- round(submission_height_in * 600)

approved_reference_png <- file.path(
  fig3_parent,
  "F3_v1_4_8_QC",
  "figures",
  "Figure3_v1_4_8_QC_180mm_300dpi.png"
)
if (!file.exists(approved_reference_png)) {
  stop("Approved v1.4.8 300-dpi reference PNG is missing.")
}

reference_root <- file.path(out_root, "reference")
dir.create(reference_root, recursive = TRUE, showWarnings = FALSE)
reference_copy <- file.path(
  reference_root,
  "Figure3_v1_4_8_APPROVED_REFERENCE_180mm_300dpi.png"
)
copy_exact(approved_reference_png, reference_copy)

figure_png_600 <- file.path(
  fig_root,
  "Figure3_v1_4_8_180mm_600ppi.png"
)
figure_tiff_600 <- file.path(
  fig_root,
  "Figure3_v1_4_8_180mm_600ppi.tiff"
)
figure_pdf <- file.path(
  fig_root,
  "Figure3_v1_4_8_180mm.pdf"
)

# TRUE 600-ppi PNG: rerender from the frozen code-layer figure object.
ragg::agg_png(
  filename = figure_png_600,
  width = expected_width_px_600,
  height = expected_height_px_600,
  units = "px",
  res = 600,
  background = "white"
)
print(figure_object)
grDevices::dev.off()

# TRUE 600-ppi TIFF: rerender from the same frozen figure object.
ragg::agg_tiff(
  filename = figure_tiff_600,
  width = expected_width_px_600,
  height = expected_height_px_600,
  units = "px",
  res = 600,
  compression = "lzw",
  background = "white"
)
print(figure_object)
grDevices::dev.off()

# Vector PDF at the exact same physical dimensions.
if (!isTRUE(base::capabilities("cairo"))) {
  stop("Cairo PDF support is required for the submission PDF export.")
}

grDevices::cairo_pdf(
  filename = figure_pdf,
  width = submission_width_in,
  height = submission_height_in,
  family = figure_font,
  bg = "white",
  onefile = TRUE
)
print(figure_object)
grDevices::dev.off()

for (p in c(figure_png_600, figure_tiff_600, figure_pdf)) {
  if (!file.exists(p) || file.info(p)$size < 50000) {
    stop("Submission-format export failed or is unexpectedly small: ", p)
  }
}

# Raster geometry verification.
png600_info <- magick::image_info(magick::image_read(figure_png_600))
tif600_info <- magick::image_info(magick::image_read(figure_tiff_600))

if (png600_info$width[[1]] != expected_width_px_600 ||
    png600_info$height[[1]] != expected_height_px_600) {
  stop(
    "600-ppi PNG pixel geometry mismatch. Expected ",
    expected_width_px_600, " x ", expected_height_px_600,
    "; observed ", png600_info$width[[1]], " x ", png600_info$height[[1]]
  )
}

if (tif600_info$width[[1]] != expected_width_px_600 ||
    tif600_info$height[[1]] != expected_height_px_600) {
  stop(
    "600-ppi TIFF pixel geometry mismatch. Expected ",
    expected_width_px_600, " x ", expected_height_px_600,
    "; observed ", tif600_info$width[[1]], " x ", tif600_info$height[[1]]
  )
}

# PDF structural preflight.
pdf_con <- file(figure_pdf, "rb")
pdf_raw <- readBin(pdf_con, what = "raw", n = 5L)
close(pdf_con)
pdf_header <- rawToChar(pdf_raw)
if (!identical(pdf_header, "%PDF-")) {
  stop("PDF export does not have a valid %PDF- header.")
}

submission_format_qc <- data.frame(
  format = c("PNG", "TIFF", "PDF"),
  file = c(
    basename(figure_png_600),
    basename(figure_tiff_600),
    basename(figure_pdf)
  ),
  physical_width_mm = c(
    submission_width_mm,
    submission_width_mm,
    submission_width_mm
  ),
  physical_height_mm = c(
    submission_height_mm,
    submission_height_mm,
    submission_height_mm
  ),
  ppi = c(600, 600, NA),
  expected_width_px = c(expected_width_px_600, expected_width_px_600, NA),
  expected_height_px = c(expected_height_px_600, expected_height_px_600, NA),
  observed_width_px = c(
    png600_info$width[[1]],
    tif600_info$width[[1]],
    NA
  ),
  observed_height_px = c(
    png600_info$height[[1]],
    tif600_info$height[[1]],
    NA
  ),
  structural_check = c(
    "PASS_TRUE_600PPI_GEOMETRY",
    "PASS_TRUE_600PPI_GEOMETRY",
    "PASS_VALID_PDF_HEADER"
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  submission_format_qc,
  file.path(qc_root, "FIG3_SUBMISSION_FORMAT_PREFLIGHT_v1_4_8.csv")
)

submission_hash_manifest <- data.frame(
  role = c(
    "approved_reference_300ppi_png",
    "submission_600ppi_png",
    "submission_600ppi_tiff",
    "submission_vector_pdf"
  ),
  file = c(
    basename(reference_copy),
    basename(figure_png_600),
    basename(figure_tiff_600),
    basename(figure_pdf)
  ),
  sha256 = c(
    sha256_file(reference_copy),
    sha256_file(figure_png_600),
    sha256_file(figure_tiff_600),
    sha256_file(figure_pdf)
  ),
  bytes = c(
    as.numeric(file.info(reference_copy)$size),
    as.numeric(file.info(figure_png_600)$size),
    as.numeric(file.info(figure_tiff_600)$size),
    as.numeric(file.info(figure_pdf)$size)
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  submission_hash_manifest,
  file.path(qc_root, "FIG3_SUBMISSION_FORMAT_HASH_MANIFEST_v1_4_8.csv")
)

write_text(
  c(
    "# Figure 3 v1.4.8 submission-format QC handoff",
    "",
    "- Scientific/statistical content: FROZEN.",
    "- Analysis/statistics recomputed: NO.",
    "- Figure rebuilt only from the already frozen Figure-3 source tables/code layer.",
    "- Approved 300-ppi PNG retained as QC reference only.",
    "- Submission PNG: TRUE 600 ppi at explicit 4252 x 4346 px (~180 x 184 mm).",
    "- Submission TIFF: TRUE 600 ppi at explicit 4252 x 4346 px (~180 x 184 mm).",
    "- Submission PDF: vector/Cairo PDF at 180 x 184 mm.",
    "- Final lock is NOT yet authorized until cross-format actual-render QC passes."
  ),
  file.path(qc_root, "FIG3_SUBMISSION_FORMATS_HANDOFF_v1_4_8.md")
)

# 7. QC evidence / frozen revision scope
# ----------------------------------------------------------------------

revision_scope <- data.frame(
  issue_id = "F3V12-01",
  severity = "BLOCKING",
  v1_2_finding = paste(
    "Actual 300-dpi PNG showed bottom-edge clipping:",
    "Panel C x-axis title and Panel D colorbar tick labels extended beyond the visible canvas."
  ),
  v1_3_action = paste(
    "Reserve an explicit 8-mm empty bottom row within the existing 184-mm canvas;",
    "rebalance body rows from 75/75 mm to 72/70 mm;",
    "retain all data, scales, labels, colors and typography sizes."
  ),
  allowed_change_scope = "LAYOUT_ONLY_BOTTOM_CANVAS_RESERVE",
  scientific_change = "NONE",
  statistical_change = "NONE",
  stringsAsFactors = FALSE
)
write_csv_utf8(revision_scope, file.path(qc_root, "FIG3_V1_2_TO_V1_3_REVISION_SCOPE_LOCK.csv"))

layout_exception <- data.frame(
  component = "figure_height",
  article_profile_value = "178 mm",
  figure3_candidate_value = "184 mm",
  rationale = paste(
    "Figure 3 is a dense 2x2 composite. The 184-mm candidate preserves the locked",
    "11.0/9.2/7.2-pt header system and final-size readability without shrinking text."
  ),
  status = "PENDING_USER_APPROVAL",
  scientific_effect = "NONE",
  statistical_effect = "NONE",
  stringsAsFactors = FALSE
)
write_csv_utf8(layout_exception, file.path(qc_root, "FIG3_LAYOUT_EXCEPTION_v1_3.csv"))

font_preflight <- data.frame(
  requested_font_family = figure_font,
  resolved_font_family = font_info$family[[1]],
  resolved_font_path = normalize_slashes(font_match$path[[1]]),
  exact_family_match = font_ok,
  skill_version = "2.2-beta.6",
  stringsAsFactors = FALSE
)
write_csv_utf8(font_preflight, file.path(qc_root, "FIG3_FONT_PREFLIGHT_v1_3.csv"))
write_csv_utf8(header_geometry, file.path(qc_root, "FIG3_PANEL_HEADER_GEOMETRY_v1_3.csv"))

source_manifest <- data.frame(
  file = basename(source_paths),
  source_version = source_version,
  sha256 = vapply(source_paths, sha256_file, character(1)),
  bytes = as.numeric(file.info(source_paths)$size),
  stringsAsFactors = FALSE
)
write_csv_utf8(source_manifest, file.path(qc_root, "FIG3_SOURCE_HASH_MANIFEST_v1_3.csv"))

candidate_status <- data.frame(
  figure_id = "Figure 3",
  figure_version = figure_version,
  status = "HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW",
  candidate_class = "QC_PNG_ONLY_BEFORE_LOCK",
  governing_skill = governing_skill,
  source_lock = source_lock,
  bulk_projection_version = bulk_projection_version,
  scientific_QC = "PASS_INHERITED_LOCKED_SOURCE",
  statistical_QC = "PASS_NO_RECOMPUTATION",
  actual_render_QC = "PENDING_USER_AND_CHAT_REVIEW",
  cross_format_QC = "NOT_RUN_BEFORE_NATIVE_RENDER_PASS",
  manuscript_ready = FALSE,
  immutable_release = FALSE,
  stringsAsFactors = FALSE
)
write_csv_utf8(candidate_status, file.path(qc_root, "FIGURE_STATUS_v1_3_QC.csv"))

manual_qc <- data.frame(
  check_id = sprintf("FIG3V13-%02d", 1:12),
  check = c(
    "A and B panel tags/titles share one top baseline",
    "C and D panel tags/titles share one top baseline",
    "Panel B subtitle is fully visible and not clipped",
    "All titles start at each panel scientific-body left edge",
    "Tag-to-title gap is visually identical across A/B/C/D",
    "Panel A heatmap values and asterisks are legible",
    "Panel B size differences and significant outlines are legible",
    "Panel C labels 69, 1, 0, 62 and 304 are all visible",
    "Panel D gene names, four contexts, values and asterisks are legible",
    "Blue-neutral-red direction semantics remain consistent in A/B/D",
    "No panel implies pooled replication, causality or definitive direct targets",
    "Panel C x-axis title and Panel D colorbar ticks are fully visible above the bottom canvas edge"
  ),
  status = "PENDING_ACTUAL_RENDER_REVIEW",
  stringsAsFactors = FALSE
)
write_csv_utf8(manual_qc, file.path(qc_root, "FIG3_STRICT_RENDER_QC_v1_3_PENDING.csv"))

write_text(
  c(
    "# Figure 3 v1.4.8 QC handoff",
    "",
    paste0("- Governing skill: ", governing_skill),
    paste0("- Source lock: ", source_lock),
    paste0("- Bulk projection: ", bulk_projection_version),
    "- Statistics recomputed: NO",
    "- Export before visual lock: ONE 300-dpi QC PNG only",
    "- Current status: HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW",
    "- Cross-format export/verification: deferred until native render passes and the user approves the candidate",
    "- Skill 8 immutable lock: NOT AUTHORIZED YET",
    "",
    "Frozen v1.3 repair scope:",
    "1. Correct only the bottom-canvas clipping seen in the v1.2 actual render.",
    "2. Add an explicit 8-mm empty bottom row and rebalance body-row heights within the same 184-mm canvas.",
    "3. Preserve all scientific/statistical content, scales, values, labels, colors and typography sizes.",
    "4. 184-mm Figure-3-specific height exception remains PENDING_USER_APPROVAL."
  ),
  file.path(qc_root, "FIG3_V1_3_QC_HANDOFF.md")
)

# Copy unchanged manuscript-facing text for review context.
if (dir.exists(text_root_v11)) {
  txt_files <- list.files(text_root_v11, full.names = TRUE)
  for (p in txt_files[file.exists(txt_files) & !dir.exists(txt_files)]) {
    copy_exact(p, file.path(text_copy_root, basename(p)))
  }
}
copy_exact(script_path, file.path(code_root, script_name))
writeLines(capture.output(sessionInfo()), file.path(qc_root, "sessionInfo_FIG3_v1_4_8.txt"))

write_text(
  c(
    "# Figure 3 v1.4.8 data-body-aligned legend handoff",
    "",
    "- Statistical recomputation: NO.",
    "- Raster patching: NO.",
    "- Pixel-shift guessing: NO.",
    "- Panel A/D default guides removed from rendered figure.",
    "- Fill scales/heights are measured from reference ggplots; bar widths are preserved from approved v1.3 actual render.",
    "- Constraint 1: colorbar body right edge == actual last-tile/data-body right edge.",
    "- Constraint 2: legend-title text bottom edge == colorbar body bottom edge.
    - Endpoint labels use inward justification to prevent clipping.",
    "- Next gate: strict actual-render QC, followed by final user visual confirmation."
  ),
  file.path(qc_root, "FIG3_V1_4_8_DATA_BODY_ALIGNED_HANDOFF.md")
)

# ----------------------------------------------------------------------
# 8. Desktop SMALL_REPORT ZIP for upload
# ----------------------------------------------------------------------

desktop_root <- find_desktop()
report_root <- file.path(desktop_root, paste0("FIG3_V1_4_8_SUBMISSION_FORMATS_QC3_REPORT_", run_stamp))
zip_path <- file.path(desktop_root, paste0("STEP05C_FIG3_V1_4_8_SUBMISSION_FORMATS_QC3_", run_stamp, "_SMALL_REPORT.zip"))
if (dir.exists(report_root) || file.exists(zip_path)) stop("Desktop report target already exists.")

report_fig3 <- file.path(report_root, "FIGURE3")
for (subdir in c("figures", "qc", "source_data", "code", "manuscript_text")) {
  dir.create(file.path(report_fig3, subdir), recursive = TRUE, showWarnings = FALSE)
}

copy_exact(figure_png, file.path(report_fig3, "figures", basename(figure_png)))
for (p in list.files(qc_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "qc", basename(p)))
}
for (p in list.files(source_copy_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "source_data", basename(p)))
}
for (p in list.files(code_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "code", basename(p)))
}
for (p in list.files(text_copy_root, full.names = TRUE)) {
  if (file.exists(p) && !dir.exists(p)) copy_exact(p, file.path(report_fig3, "manuscript_text", basename(p)))
}

report_files <- list.files(report_root, recursive = TRUE, full.names = TRUE)
report_files <- report_files[file.exists(report_files) & !dir.exists(report_files)]
report_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(report_files),
    nchar(normalize_slashes(report_root)) + 2L
  ),
  bytes = as.numeric(file.info(report_files)$size),
  sha256 = vapply(report_files, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))

zip::zipr(
  zipfile = zip_path,
  files = basename(report_root),
  root = dirname(report_root),
  include_directories = TRUE
)
if (!file.exists(zip_path) || file.info(zip_path)$size < 50000) stop("Desktop ZIP creation failed.")

cat("\n============================================================\n")
cat("Figure 3 v1.4.8 QC candidate generated.\n")
cat("Gate: HOLD_FOR_STRICT_ACTUAL_RENDER_REVIEW\n")
cat("Scientific/statistical recomputation: NONE\n")
cat("QC PNG:\n", normalize_slashes(figure_png), "\n", sep = "")
cat("Upload this ZIP:\n", normalize_slashes(zip_path), "\n", sep = "")
cat("Do NOT run Skill 8 lock yet.\n")
cat("============================================================\n")
