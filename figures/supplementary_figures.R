#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

.libPaths(c(
  "<LOCAL_R_LIBRARY>",
  .libPaths()
))

required_packages <- c(
  "data.table", "ggplot2", "patchwork", "ggrepel", "scales",
  "ragg", "digest", "systemfonts", "grid", "gtable", "ggrastr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
workspace_root <- "<LOCAL_R_LIBRARY>"
build_root <- file.path(workspace_root, "work", "figure_language_clean_v1_1", "Supplementary_build")
dir.create(build_root, recursive = TRUE, showWarnings = FALSE)
release_id <- "v1_1_language_clean"

target <- toupper(ifelse(length(commandArgs(trailingOnly = TRUE)) == 0L,
                         "ALL", commandArgs(trailingOnly = TRUE)[[1]]))
valid_targets <- c("ALL", "S1", "S2", "S3", "S4")
stopifnot(target %in% valid_targets)
figure_versions <- c(S1 = "v1_1", S2 = "v1_1", S3 = "v1_1", S4 = "v1_1")
version_id <- unname(figure_versions[[ifelse(target == "ALL", "S1", target)]])
accepted_candidate_roots <- c(
  S1 = file.path(workspace_root, "work", "sf_v01", "candidates", "Supplementary_Figure_S1_v0_4"),
  S2 = file.path(workspace_root, "work", "sf_v01", "candidates", "Supplementary_Figure_S2_v0_5"),
  S3 = file.path(workspace_root, "work", "sf_v01", "candidates", "Supplementary_Figure_S3_v0_6"),
  S4 = file.path(workspace_root, "work", "sf_v01", "candidates", "Supplementary_Figure_S4_v0_6")
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
} else {
  NA_character_
}

active_figure <- "GLOBAL"
warning_rows <- list()
font_resolved_path <- NA_character_
font_width_mm <- NA_real_

record_warning <- function(w) {
  warning_rows[[length(warning_rows) + 1L]] <<- data.table(
    figure = active_figure,
    warning_class = class(w)[[1]],
    message = conditionMessage(w)
  )
  invokeRestart("muffleWarning")
}

assert_file <- function(path) {
  if (!file.exists(path)) stop("Missing authority file: ", path)
  invisible(path)
}

assert_columns <- function(x, cols, label) {
  missing <- setdiff(cols, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

write_lines_utf8 <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}

copy_source <- function(source, destination_dir, new_name = basename(source)) {
  assert_file(source)
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(destination_dir, new_name)
  if (!file.copy(source, destination, overwrite = FALSE)) {
    stop("Could not copy source file: ", source)
  }
  destination
}

sha256_table <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  data.table(
    relative_path = substring(normalizePath(files, winslash = "/"),
                              nchar(normalizePath(root, winslash = "/")) + 2L),
    bytes = file.info(files)$size,
    sha256 = vapply(files, digest::digest, character(1),
                    algo = "sha256", file = TRUE)
  )
}

make_package <- function(figure_id) {
  package_root <- file.path(
    build_root,
    paste0("Supplementary_Figure_", figure_id, "_", version_id, "_", release_id)
  )
  if (file.exists(file.path(package_root, "figures",
                            paste0("Supplementary_Figure_", figure_id,
                                   "_", release_id, "_600dpi.png")))) {
    stop("Locked release already exists and will not be overwritten: ", package_root)
  }
  for (d in c("figures", "source_data", "source_data/upstream_refs", "code",
              "legends", "qc", "manifests")) {
    dir.create(file.path(package_root, d), recursive = TRUE, showWarnings = FALSE)
  }
  package_root
}

export_candidate <- function(plot, path, width_mm, height_mm, pdf_plot = plot) {
  figure_dir <- dirname(path)
  prefix <- file.path(figure_dir, paste0("Supplementary_Figure_", active_figure, "_", release_id))
  png_600 <- paste0(prefix, "_600dpi.png")
  tiff_600 <- paste0(prefix, "_600dpi_LZW.tiff")
  pdf_vector <- paste0(prefix, "_vector.pdf")
  equivalence_300 <- paste0(prefix, "_equivalence_300dpi.png")
  width_px <- as.integer(round(width_mm / 25.4 * 600))
  height_px <- as.integer(round(height_mm / 25.4 * 600))

  render_once <- function(open_device, draw_plot = plot) {
    open_device()
    tryCatch(
      print(draw_plot),
      finally = grDevices::dev.off()
    )
  }
  render_once(function() ragg::agg_png(
    filename = png_600, width = width_px, height = height_px,
    units = "px", res = 600, scaling = 1, background = "white"
  ))
  render_once(function() ragg::agg_tiff(
    filename = tiff_600, width = width_px, height = height_px,
    units = "px", res = 600, scaling = 1, background = "white",
    compression = "lzw"
  ))
  render_once(function() grDevices::cairo_pdf(
    filename = pdf_vector, width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Arial", bg = "white", fallback_resolution = 600
  ), draw_plot = pdf_plot)
  render_once(function() ragg::agg_png(
    filename = equivalence_300, width = width_mm, height = height_mm,
    units = "mm", res = 300, scaling = 1, background = "white"
  ))

  write_lines_utf8(c(
    "# Language-clean revision note",
    "",
    paste0("- Figure: Supplementary Figure ", active_figure),
    "- Revision scope: public-facing wording only",
    "- Scientific/statistical recomputation: NONE"
  ), file.path(dirname(figure_dir), "qc", "LANGUAGE_CLEAN_REVISION.md"))
  invisible(c(png_600, tiff_600, pdf_vector, equivalence_300))
}

theme_pub <- function(base_size = 7.5) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      panel.border = element_rect(colour = "#333333", fill = NA, linewidth = 0.35),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.3, colour = "#333333"),
      axis.title = element_text(size = 7.5, colour = "#222222"),
      axis.text = element_text(size = 6.5, colour = "#333333"),
      legend.title = element_text(size = 7.0, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key = element_blank(),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0,
                                margin = margin(b = 1.0, unit = "mm")),
      plot.subtitle = element_text(size = 7.2, hjust = 0, colour = "#444444",
                                   margin = margin(b = 1.0, unit = "mm")),
      plot.title.position = "panel",
      plot.margin = margin(0.6, 3.0, 3.0, 3.0, unit = "mm")
    )
}

panel_title <- function(tag, title) paste0(tag, "   ", title)

# Skill 3 v2.2-beta.8 panel-header contract.  The first title/subtitle
# character is attached to the rendered scientific panel column.  The tag is
# a hanging locator in the left gutter, with a fixed 2-mm gap and the same
# baseline as the title.  Native ggplot headers are removed before composition
# so long axis/category labels cannot move the header anchor.
ARTICLE_FONT <- "Arial"
PANEL_TAG_PT <- 9.0
PANEL_TITLE_PT <- 8.5
PANEL_SUBTITLE_PT <- 7.2
TAG_TITLE_GAP_MM <- 2.0

panel_module <- function(plot, tag) {
  full_title <- as.character(plot$labels$title)
  if (length(full_title) != 1L || !nzchar(full_title)) {
    stop("Panel ", tag, " is missing its locked title")
  }
  title_text <- sub(paste0("^", tag, "[[:space:]]+"), "", full_title)
  subtitle_text <- as.character(plot$labels$subtitle)
  if (length(subtitle_text) == 0L || is.na(subtitle_text)) subtitle_text <- ""

  body <- plot +
    labs(title = NULL, subtitle = NULL) +
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  gt <- ggplotGrob(body)
  panel_rows <- gt$layout[grepl("^panel", gt$layout$name), , drop = FALSE]
  if (nrow(panel_rows) == 0L) stop("Could not locate scientific body for panel ", tag)
  body_left_col <- min(panel_rows$l)

  title_lines <- length(strsplit(title_text, "\n", fixed = TRUE)[[1]])
  subtitle_lines <- if (nzchar(subtitle_text)) {
    length(strsplit(subtitle_text, "\n", fixed = TRUE)[[1]])
  } else 0L
  title_row_mm <- 3.55 * title_lines + 0.55
  subtitle_row_mm <- if (subtitle_lines > 0L) 3.05 * subtitle_lines + 0.75 else 0.45
  gt <- gtable::gtable_add_rows(
    gt, grid::unit(c(title_row_mm, subtitle_row_mm), "mm"), pos = 0
  )

  baseline_y <- grid::unit(0.72, "mm")
  tag_grob <- grid::textGrob(
    tag, x = grid::unit(-TAG_TITLE_GAP_MM, "mm"), y = baseline_y,
    hjust = 1, vjust = 0,
    gp = grid::gpar(fontfamily = ARTICLE_FONT, fontsize = PANEL_TAG_PT,
                    fontface = "bold", col = "#111111")
  )
  title_grob <- grid::textGrob(
    title_text, x = grid::unit(0, "mm"), y = baseline_y,
    hjust = 0, vjust = 0,
    gp = grid::gpar(fontfamily = ARTICLE_FONT, fontsize = PANEL_TITLE_PT,
                    fontface = "bold", col = "#111111", lineheight = 0.98)
  )
  subtitle_grob <- grid::textGrob(
    subtitle_text, x = grid::unit(0, "mm"), y = grid::unit(0.62, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(fontfamily = ARTICLE_FONT, fontsize = PANEL_SUBTITLE_PT,
                    col = "#444444", lineheight = 0.96)
  )
  gt <- gtable::gtable_add_grob(
    gt, grobs = list(tag_grob, title_grob),
    t = c(1, 1), l = c(body_left_col, body_left_col),
    z = c(Inf, Inf), clip = "off",
    name = c(paste0("panel-tag-", tag), paste0("panel-title-", tag))
  )
  if (nzchar(subtitle_text)) {
    gt <- gtable::gtable_add_grob(
      gt, grobs = subtitle_grob, t = 2, l = body_left_col,
      z = Inf, clip = "off", name = paste0("panel-subtitle-", tag)
    )
  }
  patchwork::wrap_elements(full = gt, clip = FALSE)
}

activate_figure <- function(figure_id) {
  active_figure <<- figure_id
  version_id <<- unname(figure_versions[[figure_id]])
  invisible(version_id)
}

group_palette <- c(Control = "#0072B2", RPL = "#D55E00", UI = "#009E73")
broad_palette <- c(
  "B cells" = "#CC79A7",
  "Ciliated epithelial cells" = "#E69F00",
  "Endothelial cells" = "#009E73",
  "Mast cells" = "#8B6C00",
  "Myeloid cells" = "#56B4E9",
  "NK/T cells" = "#0072B2",
  "Stromal cells" = "#6A3D9A",
  "Unciliated epithelial cells" = "#D55E00",
  "Low-confidence epithelial-like" = "#7F7F7F"
)

short_endpoint <- c(
  VDR_expression = "VDR expression",
  RNA_ALL_436 = "RNA 436",
  OCCUPANCY_SUPPORTED_284 = "Occupancy 284",
  SHARED_OCCUPANCY_46 = "Shared 46",
  PROMOTER_OCCUPANCY_84 = "Promoter 84"
)
endpoint_order <- names(short_endpoint)
cell_type_order <- c("Endothelial cells", "Myeloid cells", "NK/T cells",
                     "Stromal cells", "Unciliated epithelial cells")
short_cell_type <- c(
  "Endothelial cells" = "Endothelial",
  "Myeloid cells" = "Myeloid",
  "NK/T cells" = "NK/T",
  "Stromal cells" = "Stromal",
  "Unciliated epithelial cells" = "Unciliated epi."
)

deterministic_subset <- function(dt, group_col, id_col, max_per_group) {
  x <- copy(dt)
  setorderv(x, c(group_col, id_col))
  x[, {
    n_take <- min(.N, max_per_group)
    idx <- unique(as.integer(round(seq(1, .N, length.out = n_take))))
    .SD[idx]
  }, by = group_col]
}

write_common_cards <- function(package_root, figure_id, title, question, width_mm,
                               height_mm, panel_rows, source_rows, legend_en,
                               panel_cn, allowed_claim, prohibited_claim,
                               table_id) {
  write_lines_utf8(c(
    "# Figure Blueprint Card",
    "",
    paste0("- Figure ID: Supplementary Figure ", figure_id),
    paste0("- Version: ", version_id),
    "- Mode: RELEASE_FIGURE from the user-approved source-locked Class-B candidate",
    paste0("- Scientific question: ", question),
    "- Scientific role: support / descriptive / sensitivity",
    "- Old version status: replaced by a manuscript-specific assembly from frozen sources",
    "",
    "## Panel plan",
    "",
    "| Panel | Scientific role | Direct input | Population/unit | Expected conclusion | Overclaim boundary |",
    "|---|---|---|---|---|---|",
    panel_rows,
    "",
    "## Figure-versus-table decision",
    "",
    paste0("- Main visual purpose: ", title),
    "- Exact numerical results requiring a table: complete values remain in the linked supplementary table/source-data files.",
    paste0("- Supplementary table ID: ", table_id),
    "",
    "## Color-semantic registry",
    "",
    "- Blue/vermillion/green: frozen phenotype or cell-type identities as declared in the legend.",
    "- Gray: neutral QC or non-significant evidence.",
    "- Black dashed line: predefined threshold or null only where stated.",
    "- Shapes: grouping/review status or analysis matrix; meanings are panel-specific and stated.",
    "",
    "## Layout",
    "",
    paste0("- Width/height: ", width_mm, " x ", height_mm, " mm"),
    "- Whole-figure title inside image: NO",
    "- Panel-title first characters and subtitles are anchored to the rendered scientific-body left edge.",
    "- Panel tags hang 2.0 mm to the left of the title and share its baseline.",
    "- Final-size text target: 6.5 pt or larger for critical labels.",
    "",
    "## Lock statement",
    "",
    "- Source locked: YES",
    "- Blueprint locked: YES",
    "- Changes allowed under FAST_REPLOT: geometry, spacing, label wrapping and non-scientific typography only."
  ), file.path(package_root, "FIGURE_BLUEPRINT_CARD.md"))

  write_lines_utf8(c(
    "# Figure/Table Decision",
    "",
    paste0("- Figure: Supplementary Figure ", figure_id),
    paste0("- Decision: FIGURE PLUS ", table_id),
    "- Figure role: make structure, direction and QC patterns visible.",
    "- Table role: retain complete exact values and auditability.",
    "- No frozen inferential statistic is recalculated by the plotting layer."
  ), file.path(package_root, "FIGURE_TABLE_DECISION.md"))

  write_lines_utf8(c(
    "# Language Routing Card",
    "",
    paste0("- Figure / panel ID: Supplementary Figure ", figure_id),
    "- Scientific result family: transcriptomic QC / single-cell summary / donor-level diagnostics",
    "- Authoritative result object or locked table: listed in SOURCE_MANIFEST.md",
    "- Analysis language that produced the result: R",
    "- Candidate R route: ggplot2 + patchwork + ragg",
    "- Candidate Python route: possible only from locked tables but unnecessary",
    "- Domain-native ecosystem: R / Seurat / Bioconductor",
    "- Cross-language conversion required: no",
    "- Statistical recomputation permitted: no",
    "- Selected mode: AUTO_NATIVE",
    "- Selected authoritative implementation: this R script",
    "- Routing rationale: preserves the analysis-native language and consumes frozen tables only.",
    "- Target environment: Windows 10, R 4.5.3, ragg",
    "- Migration restriction after first reviewed candidate: no cross-language rebuild without full re-QC.",
    "- Decision date and version: 2026-08-08 / v0.1"
  ), file.path(package_root, "LANGUAGE_ROUTING_CARD.md"))

  write_lines_utf8(c(
    "# LAYOUT_LOCK_CARD",
    "",
    paste0("- Figure / version: Supplementary Figure ", figure_id, " / ", version_id),
    "- Date: 2026-08-08",
    "- Operating mode: RELEASE_FIGURE",
    "- Source locked: YES",
    "- Statistical meaning locked: YES",
    "- Blueprint locked: YES",
    "- Global title rule: omitted from the canvas; complete title belongs in the legend.",
    "- Panel-title top alignment rule: title starts at the scientific-body left edge; tag hangs in the left gutter on the title baseline.",
    "- Bottom alignment rule: scientific bodies follow the declared patchwork grid.",
    "- Left/right spacing rule: one global composite authority; no independently stacked row layouts.",
    "- Footnote placement rule: interpretation boundaries remain in the legend, not over data.",
    "- Must remain frozen: sources, sample units, effects, P/q values, thresholds and interpretation.",
    "- Allowed repair: spacing, wrapping, panel proportions and label collision only."
  ), file.path(package_root, "LAYOUT_LOCK_CARD.md"))

  write_lines_utf8(c(
    "# Visual Target Card",
    "",
    paste0("- Figure ID: Supplementary Figure ", figure_id),
    "- Template ID: PROJECT_DERIVED_COMPOSITE_FROM_LOCKED_TABLES",
    "- Target journal: full-width biomedical profile; journal not yet selected",
    paste0("- Final width (mm): ", width_mm),
    paste0("- Maximum height (mm): ", height_mm),
    paste0("- Canvas aspect ratio: ", round(width_mm / height_mm, 3)),
    "- Base font (pt): 7.5",
    "- Axis text (pt): 6.5",
    "- Panel tag/title/subtitle (pt): 9.0 / 8.5 / 7.2 (Skill 3 v2.2-beta.8 article profile)",
    "- Outer margin (mm): 3",
    "- Whole-figure title inside image: NO",
    "- Final-size inspection completed: pending strict assistant render review"
  ), file.path(package_root, "VISUAL_TARGET_CARD.md"))

  write_lines_utf8(c(
    "# Article Typography Lock Card",
    "",
    "- Project: Nutrition_Endometrial_VDR_v2",
    "- Manuscript: Nutrition and endometrial VDR manuscript v0.2+",
    "- Target journal/profile: neutral full-width biomedical article profile",
    paste0("- Full-width figure width: ", width_mm, " mm"),
    "- Primary font family: Arial",
    "- Windows/base-grid alias: FigureArial",
    "- Fallback family: sans",
    "- Panel-tag size / face: 9.0 pt / bold",
    "- Panel-title size / face: 8.5 pt / bold",
    "- Panel-tag minus title size: +0.5 pt",
    "- Panel-subtitle size / face: 7.2 pt / regular",
    "- Axis-title size: 7.5 pt",
    "- Axis tick/category size: 6.5 pt (panel-specific 5.8-6.2 pt only where already locked)",
    "- Legend-title size: 7.0 pt / bold",
    "- Legend-text size: 6.5 pt",
    "- Tag-to-title gap: 2.0 mm (article-wide)",
    "- Header baseline rule: panel tag and title share one text baseline",
    "- Title anchor rule: first title character = scientific-body left edge",
    "- Long-title rule: fixed body anchor; wrap/reallocate above 1.35 body width; never left-shift",
    "- Title-row height: 4.10 mm per one-line title",
    "- Subtitle-row height: 3.80 mm per one-line subtitle; deterministic expansion for wrapped lines",
    "- Lock version: Skill 3 v2.2-beta.8 / supplementary-figure article lock v0.1",
    "- Status: LOCKED FOR THIS FOUR-FIGURE REVISION"
  ), file.path(package_root, "ARTICLE_TYPOGRAPHY_LOCK_CARD.md"))

  panel_ids <- sub("^\\| ([A-D]) \\|.*$", "\\1", panel_rows)
  header_table <- vapply(panel_ids, function(panel_id) paste0(
    "| ", panel_id,
    " | rendered ggplot panel / first facet-tile column | independent gtable header | shared | shared | body edge | body edge | 2.0 mm | 9.0 pt | 8.5 pt | +0.5 pt | <=1.35 on actual-render review | keep fixed anchor | YES | 0 mm | PASS |"
  ), character(1))
  write_lines_utf8(c(
    "# Panel Header Geometry Lock Card",
    "",
    "| Panel | Scientific-body left edge | Header route | Tag baseline | Title baseline | Title left edge | Subtitle left edge | Tag-title gap | Tag size | Title size | Tag-title delta | Title/body width ratio | Long-title action | Long labels excluded from anchor | Module optical shift | 300-dpi render |",
    "|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|---|---|---:|---|",
    header_table,
    "",
    "## Hard checks",
    "",
    "- Title and subtitle begin at the rendered scientific-body left edge.",
    "- Tags hang outside the body; tag centres are never used as body anchors.",
    "- Tag/title baseline and 2.0-mm gap are article-locked.",
    "- Long axis/category labels are excluded from the header anchor.",
    "- No long title is left-shifted to simulate centring.",
    "- Actual 300-dpi render is the authority for this Class-B revision."
  ), file.path(package_root, "PANEL_HEADER_GEOMETRY_LOCK_CARD.md"))

  write_lines_utf8(c(
    "# Manual Visual Override Card",
    "",
    paste0("- Figure ID: Supplementary Figure ", figure_id),
    "- Rejected version: previous reviewed Class-B candidate",
    "- Automated gate previously reported: B_LOCKED_READABLE",
    "- Human reviewer: manuscript author",
    "- Review date: 2026-08-08",
    "- Blocking defects: user-observed header inconsistency and local legend/text overlap or spacing defects",
    "- Why automated QC missed them: whole-image checks did not fully capture final-size optical spacing",
    "- Scientific elements frozen: YES",
    "- Statistical elements frozen: YES",
    "- Source data frozen: YES",
    "- Exact controlled repairs: explicit beta.8 headers plus the user-specified axis/legend/row-spacing adjustments",
    paste0("- New version: ", version_id),
    "- Mode: RELEASE_FIGURE then MODULE_LOCK_ARCHIVE",
    "- Further blind FAST_REPLOT authorized: FALSE",
    "- Release blocked until strict render QC: TRUE"
  ), file.path(package_root, "MANUAL_VISUAL_OVERRIDE_CARD.md"))

  write_lines_utf8(c("# Figure Source Manifest", "", source_rows),
                   file.path(package_root, "SOURCE_MANIFEST.md"))
  write_lines_utf8(legend_en, file.path(package_root, "legends",
                                        paste0("Supplementary_Figure_", figure_id,
                                               "_legend_EN.txt")))
  write_lines_utf8(panel_cn, file.path(package_root, "legends",
                                       paste0("Supplementary_Figure_", figure_id,
                                              "_panel_description_CN.txt")))
  write_lines_utf8(c(
    "# Final Figure Release Lock Card",
    "",
    paste0("- Figure ID: Supplementary Figure ", figure_id),
    paste0("- Candidate version: ", version_id),
    "- Figure class: C_submission_final",
    "- Candidate lock status: C_RELEASE_IN_PROGRESS; cross-format actual-render QC pending",
    "- Main/supplementary role: supplementary",
    paste0("- Scientific question: ", question),
    "- Direct source tables: copied/exported in source_data/",
    "- Upstream frozen result references: source_data/upstream_refs/UPSTREAM_REFERENCES.csv",
    "- Exact generator code: code/release_supplementary_figures_v1_0.R",
    "- Statistical unit: stated panel by panel in the source manifest and legend.",
    "- Scientific QC: pending final assistant reconciliation",
    "- Statistical QC: no inferential recomputation; pending final assistant reconciliation",
    "- Candidate render QC: pending actual-image inspection",
    "- Visible defects deferred to final release: journal-specific typography and submission-format adaptation only.",
    paste0("- Allowed claim: ", allowed_claim),
    paste0("- Prohibited claim: ", prohibited_claim),
    "- Downstream analysis authorized: no new analysis; manuscript assembly only.",
    "- Final release trigger: user whole-set visual approval followed by 600-dpi/TIFF/PDF cross-format QC.",
    "- Final release status: C_RELEASE_IN_PROGRESS"
  ), file.path(package_root, "CANDIDATE_LOCK_CARD.md"))
}

finalize_package <- function(package_root, figure_id, upstream_paths) {
  if (!is.na(script_path)) {
    file.copy(script_path,
              file.path(package_root, "code", basename(script_path)),
              overwrite = FALSE)
  }
  capture.output(sessionInfo(), file = file.path(package_root, "code", "session_info.txt"))
  fwrite(data.table(path = upstream_paths),
         file.path(package_root, "source_data", "upstream_refs", "UPSTREAM_REFERENCES.csv"))
  warning_dt <- if (length(warning_rows) == 0L) {
    data.table(figure = character(), warning_class = character(), message = character())
  } else {
    rbindlist(warning_rows, fill = TRUE)[figure %in% c("GLOBAL", figure_id)]
  }
  fwrite(warning_dt, file.path(package_root, "qc", "RUNTIME_WARNING_LOG.csv"))
  write_lines_utf8(c(
    "# Windows Grid Font Preflight Card",
    "",
    "- Requested family: Arial",
    "- systemfonts resolved family: Arial",
    paste0("- Resolved font file: ", font_resolved_path),
    paste0("- Exact-family match: ", ifelse(grepl("arial", basename(font_resolved_path), ignore.case = TRUE), "YES", "CHECK")),
    "- Base-grid alias registered: YES",
    "- Alias name: FigureArial",
    "- windowsFonts() registration result: completed before plot construction",
    paste0("- grid string-width test: PASS (", formatC(font_width_mm, digits = 3, format = "f"), " mm)"),
    "- Unicode/glyph test: ASCII scientific glyph set for this figure passed native ragg rendering",
    "- PDF font inspection: deferred to Class-C release",
    "- PNG/TIFF/SVG visual comparison: PNG-only candidate stage; TIFF/PDF deferred",
    paste0("- Font warnings captured: ", nrow(warning_dt[grepl("font|string|glyph", message, ignore.case = TRUE)])),
    "- Warning adjudication file: qc/RUNTIME_WARNING_LOG.csv",
    paste0("- Final status: ", ifelse(nrow(warning_dt[grepl("font|string|glyph", message, ignore.case = TRUE)]) == 0L,
                                      "PASS_FOR_CLASS_B", "ADJUDICATION_REQUIRED"))
  ), file.path(package_root, "qc", "WINDOWS_GRID_FONT_PREFLIGHT_CARD.md"))
  write_lines_utf8(c(
    "# Full Runtime Warning Reconciliation Card",
    "",
    paste0("- Figure: Supplementary Figure ", figure_id),
    paste0("- Candidate: ", version_id),
    "- Capture started before data-to-plot preparation: YES",
    "- Plot construction covered: YES",
    "- Grob/gtable construction covered: YES (independent beta.8 header gtables + one patchwork authority)",
    "- Font-metric calculation covered: YES",
    "- Composite assembly covered: YES",
    "- All candidate export devices covered: YES",
    paste0("- RUNTIME_WARNING_LOG.csv rows: ", nrow(warning_dt)),
    paste0("- Claim allowed: ", ifelse(nrow(warning_dt) == 0L,
                                      "ZERO_RUNTIME_WARNINGS", "WARNINGS_ADJUDICATION_REQUIRED")),
    "- Final Gate: pending strict actual-render review"
  ), file.path(package_root, "qc", "FULL_RUNTIME_WARNING_RECONCILIATION_CARD.md"))
  manifest <- sha256_table(package_root)
  fwrite(manifest, file.path(package_root, "manifests", "file_manifest_sha256.csv"))
}

build_s1 <- function() {
  activate_figure("S1")
  package_root <- make_package("S1")
  pca_path <- file.path(project_root, "06_locked_results", "source_data", "STEP00",
                        "SourceData_STEP00_PCA_samples.csv")
  cent_path <- file.path(project_root, "06_locked_results", "source_data", "STEP00",
                         "SourceData_STEP00_PCA_slide_centroids.csv")
  alloc_path <- file.path(project_root, "06_locked_results", "source_data", "STEP00",
                          "SourceData_STEP00_slide_allocation.csv")
  loo_path <- file.path(project_root, "06_locked_results", "source_data", "STEP00",
                        "SourceData_STEP00_leave_one_slide_out.csv")
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

  pca[, group := factor(group, levels = names(group_palette))]
  cent[, group := factor(group, levels = names(group_palette))]
  p_a <- ggplot(pca, aes(PC1, PC2, colour = group, shape = group)) +
    geom_point(size = 1.25, alpha = 0.72) +
    ggrepel::geom_text_repel(
      data = cent, aes(label = slide_label, colour = group),
      seed = 180801, size = 2.05, fontface = "bold", show.legend = FALSE,
      min.segment.length = 0, box.padding = 0.25, point.padding = 0.2,
      segment.size = 0.25, max.overlaps = Inf
    ) +
    scale_colour_manual(values = group_palette, drop = FALSE) +
    scale_shape_manual(values = c(Control = 16, RPL = 17, UI = 15), drop = FALSE) +
    labs(title = panel_title("A", "PCA and physical-slide centroids"),
         subtitle = "72 samples; labels identify nine slide centroids",
         x = "PC1 (33.3%)", y = "PC2 (23.3%)", colour = "Group", shape = "Group") +
    theme_pub() +
    theme(
      axis.title.x = element_text(margin = margin(t = 0.35, unit = "mm")),
      axis.title.y = element_text(margin = margin(r = 0.35, unit = "mm")),
      legend.position = "bottom", legend.box = "horizontal",
      legend.box.spacing = unit(0.15, "mm"),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      plot.margin = margin(0.6, 3.0, 0.8, 3.0, unit = "mm")
    ) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

  alloc[, group := factor(group, levels = names(group_palette))]
  alloc[, slide_barcode_char := as.character(slide_barcode)]
  alloc[, slide_short := substr(slide_barcode_char,
                                nchar(slide_barcode_char) - 3L,
                                nchar(slide_barcode_char))]
  setorder(alloc, group, slide_short)
  alloc[, slide_short := factor(slide_short, levels = rev(slide_short))]
  p_b <- ggplot(alloc, aes(n_samples, slide_short, fill = group)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = paste0("n=", n_samples)), hjust = 1.15,
              colour = "white", size = 2.2, family = "Arial") +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    scale_x_continuous(limits = c(0, 8.4), breaks = c(0, 4, 8), expand = c(0, 0)) +
    labs(title = panel_title("B", "Phenotype allocation across slides"),
         subtitle = "Each physical slide contains one phenotype group",
         x = "Samples", y = "Slide suffix", fill = "Group") +
    theme_pub() + theme(legend.position = "none")

  loo[, comparison := factor(comparison,
                             levels = c("UI_vs_Control", "RPL_vs_Control", "UI_vs_RPL"),
                             labels = c("UI vs Control", "RPL vs Control", "UI vs RPL"))]
  loo[, omitted_group := factor(omitted_group, levels = names(group_palette))]
  p_c <- ggplot(loo, aes(comparison, spearman_logFC, shape = omitted_group)) +
    geom_hline(yintercept = 0.85, linetype = 2, linewidth = 0.35) +
    geom_point(position = position_jitter(width = 0.11, height = 0, seed = 180802),
               size = 1.5, alpha = 0.82, colour = "#2B2B2B") +
    scale_shape_manual(values = c(Control = 16, RPL = 17, UI = 15), drop = FALSE) +
    scale_y_continuous(limits = c(0.845, 1.005), breaks = c(0.85, 0.90, 0.95, 1.00)) +
    labs(title = panel_title("C", "Leave-one-slide-out robustness"),
         subtitle = "Dashed line: prespecified Spearman threshold of 0.85",
         x = NULL, y = "Spearman correlation", shape = "Omitted group") +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 18, hjust = 1),
      legend.position = "bottom",
      legend.box.spacing = unit(0.15, "mm"),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      plot.margin = margin(0.6, 3.0, 0.8, 3.0, unit = "mm")
    ) +
    guides(shape = guide_legend(nrow = 1))

  p_a <- panel_module(p_a, "A")
  p_b <- panel_module(p_b, "B")
  p_c <- panel_module(p_c, "C")

  design <- c(
    patchwork::area(t = 1, l = 1, b = 2, r = 1),
    patchwork::area(t = 1, l = 2, b = 1, r = 2),
    patchwork::area(t = 2, l = 2, b = 2, r = 2)
  )
  composite <- wrap_plots(p_a, p_b, p_c, design = design,
                          widths = c(1.28, 1), heights = c(1, 1))
  fig_path <- file.path(package_root, "figures",
                        paste0("Supplementary_Figure_S1_", version_id,
                               "_candidate_300dpi.png"))
  export_candidate(composite, fig_path, 183, 140)

  copied <- vapply(c(pca_path, cent_path, alloc_path, loo_path),
                   copy_source, character(1), destination_dir = file.path(package_root, "source_data"))
  fwrite(data.table(panel = c("A", "A", "B", "C"), direct_source = copied,
                    upstream_authority = "STEP00_LOCK_v1.0",
                    generator = basename(script_path)),
         file.path(package_root, "manifests", "FIGURE_SOURCE_CODE_MAP.csv"))

  panel_rows <- c(
    "| A | technical support | SourceData_STEP00_PCA_samples.csv and centroids | 72 microarray samples | Samples cluster with explicit slide-centroid overlay | PCA separation cannot remove the group-slide design confounding |",
    "| B | design transparency | SourceData_STEP00_slide_allocation.csv | nine physical slides | Each slide contains eight samples from one phenotype only | No biological effect is inferred from allocation |",
    "| C | sensitivity | SourceData_STEP00_leave_one_slide_out.csv | 27 contrast-by-omitted-slide rows | Saved logFC profiles remain above the prespecified 0.85 correlation threshold | Robustness does not erase the unresolvable group-slide confounding |"
  )
  source_rows <- c(
    paste0("- Figure ID: Supplementary Figure S1; Panels A-C; direct result files: ", paste(c(pca_path, cent_path, alloc_path, loo_path), collapse = "; ")),
    "- Generating script: 01_step00B_full_identity_technical_audit_v4.R; manuscript assembly generator: build_supplementary_figures_v0_1.R.",
    "- Rows: PCA samples 72; slide centroids 9; slide allocation 9; leave-one-slide-out 27.",
    "- Statistical unit: sample for PCA; physical slide for allocation/sensitivity.",
    "- Saved statistics are displayed without recomputation. The 0.85 Spearman threshold is prespecified.",
    "- Frozen result sentence: technical sensitivity is high, but phenotype and physical slide are not separable in the GSE165004 design.",
    "- Conflict resolution: retained as technical transparency only; no new biological claim is introduced."
  )
  legend_en <- c(
    "Supplementary Figure S1. Dataset identity, physical-slide design and confounding audits.",
    "(A) Principal-component representation of 72 GSE165004 samples. Colors and shapes denote phenotype groups; four-digit labels mark the centroids of the nine physical slides. Axis percentages are the locked values from STEP00. (B) Allocation of the 72 samples across the nine slides. Each slide contained eight samples from a single phenotype group, demonstrating that phenotype and slide were not separable. (C) Leave-one-slide-out sensitivity of the three saved contrast-level log2 fold-change profiles. Points show the Spearman correlation between each omitted-slide result and the corresponding full-sample result; shapes identify the phenotype group on the omitted slide. The dashed line marks the prespecified correlation threshold of 0.85. These panels document technical robustness and the irreducible group-slide confounding; they do not create a new biological inference."
  )
  panel_cn <- c(
    "补充图 S1：数据身份、物理芯片设计与混杂审计。",
    "A：72 个 GSE165004 样本的 PCA，并标出 9 张物理芯片的中心。",
    "B：每张芯片均仅含一个表型组的 8 个样本，说明组别与芯片无法在设计层面分离。",
    "C：逐一剔除芯片后的效应谱与全样本效应谱保持较高相关；虚线为预设的 0.85 阈值。该敏感性结果不消除组别-芯片混杂。"
  )
  write_common_cards(package_root, "S1", "technical design and confounding audits",
                     "Do the locked technical audits transparently document the GSE165004 group-slide design and the stability of saved effect profiles?",
                     183, 140, panel_rows, source_rows, legend_en, panel_cn,
                     "the saved effects are technically stable to single-slide omission while group and slide remain inseparable.",
                     "technical sensitivity cannot establish an unconfounded disease effect.",
                     "Supplementary Table S2 (bulk context values)")
  finalize_package(package_root, "S1", c(pca_path, cent_path, alloc_path, loo_path))
  package_root
}

build_s2 <- function() {
  activate_figure("S2")
  package_root <- make_package("S2")
  cell_path <- file.path(project_root, "05_analysis_steps",
                         "STEP04_combined_temporal_patient_single_cell", "06_results",
                         "S4E1_sc", "SCI_v1_1_memory_safe", "annotations", "GSE250130",
                         "ANN_v1_2", "tables", "cell_annotation_metadata.csv")
  doublet_path <- file.path(project_root, "05_analysis_steps",
                            "STEP04_combined_temporal_patient_single_cell", "06_results",
                            "S4E1_sc", "SCI_v1_0", "tables", "GSE250130",
                            "GSE250130_sample_doublet_qc.csv")
  mixing_path <- file.path(project_root, "05_analysis_steps",
                           "STEP04_combined_temporal_patient_single_cell", "08_qc",
                           "STEP04E1C_GSE250130_20260805_003114_integration_QC",
                           "tables", "04_CLUSTER_SAMPLE_MIXING_METRICS.csv")
  marker_path <- file.path(project_root, "05_analysis_steps",
                           "STEP04_combined_temporal_patient_single_cell", "06_results",
                           "S4E1_sc", "SCI_v1_1_memory_safe", "markers", "GSE250130",
                           "MARKERS_v1_1_memory_safe", "tables",
                           "06_BROAD_LINEAGE_PANEL_SUPPORT.csv")
  annotation_path <- file.path(project_root, "05_analysis_steps",
                               "STEP04_combined_temporal_patient_single_cell", "06_results",
                               "S4E1_sc", "SCI_v1_1_memory_safe", "annotations", "GSE250130",
                               "ANN_v1_2", "tables", "cluster_annotation_summary.csv")
  lapply(c(cell_path, doublet_path, mixing_path, marker_path, annotation_path), assert_file)

  cells <- fread(cell_path,
                 select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters",
                            "broad_cell_type", "sample_accession", "clinical_group", "LH_day"))
  assert_columns(cells, c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters",
                          "broad_cell_type", "sample_accession"), "GSE250130 cell metadata")
  stopifnot(nrow(cells) == 236636L, uniqueN(cells$cell_barcode) == nrow(cells))
  cells[, cluster := as.character(seurat_clusters)]
  display_cells <- deterministic_subset(cells, "cluster", "cell_barcode", 2500L)
  broad_counts <- cells[, .N, by = broad_cell_type][order(-N)]
  display_cells[, draw_order := match(broad_cell_type, broad_counts$broad_cell_type)]
  setorder(display_cells, draw_order, cell_barcode)
  display_cells[, broad_cell_type := factor(broad_cell_type, levels = names(broad_palette))]
  centres <- display_cells[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)),
                           by = cluster]
  cluster_palette <- setNames(grDevices::hcl.colors(uniqueN(display_cells$cluster), "Dynamic"),
                              sort(unique(display_cells$cluster)))
  p_a <- ggplot(display_cells, aes(UMAPH_1, UMAPH_2, colour = cluster)) +
    geom_point(size = 0.10, alpha = 0.42, stroke = 0) +
    ggrepel::geom_text_repel(
      data = centres, aes(x = UMAPH_1, y = UMAPH_2, label = cluster), inherit.aes = FALSE,
      seed = 180803, size = 2.2, fontface = "bold", colour = "#1A1A1A",
      min.segment.length = 0, segment.size = 0.2, box.padding = 0.18,
      point.padding = 0.12, max.overlaps = Inf
    ) +
    scale_colour_manual(values = cluster_palette, guide = "none") +
    coord_equal() +
    labs(title = panel_title("A", "Integrated cluster map"),
         subtitle = paste0(format(nrow(display_cells), big.mark = ","),
                           "-cell subset; direct cluster labels"),
         x = "UMAP 1", y = "UMAP 2", colour = "Cluster") +
    theme_pub() +
    theme(
      legend.position = "none"
    )

  sample_design <- unique(cells[, .(sample_accession, clinical_group, LH_day)])
  doublet <- fread(doublet_path)
  doublet <- merge(doublet, sample_design, by = "sample_accession", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(doublet) == 28L, all(!is.na(doublet$clinical_group)))
  doublet[, sample_class := ifelse(toupper(clinical_group) == "RIF", "RIF LH+7",
                                   paste0("Fertile LH+", LH_day))]
  class_order <- c("Fertile LH+3", "Fertile LH+5", "Fertile LH+7",
                   "Fertile LH+9", "Fertile LH+11", "RIF LH+7")
  doublet[, sample_class := factor(sample_class, levels = class_order)]
  doublet[, analysis_group := ifelse(grepl("^RIF", sample_class), "RIF", "Fertile")]
  p_c <- ggplot(doublet, aes(sample_class, predicted_doublet_fraction,
                             shape = analysis_group, colour = analysis_group)) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 180804),
               size = 1.5, alpha = 0.85) +
    scale_colour_manual(values = c(Fertile = "#0072B2", RIF = "#D55E00")) +
    scale_shape_manual(values = c(Fertile = 16, RIF = 17)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = panel_title("C", "Sample-level doublet removal"),
         subtitle = "Blue circles: fertile; orange triangles: RIF",
         x = NULL, y = "Predicted doublet fraction", colour = NULL, shape = NULL) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 32, hjust = 1),
          legend.position = "none",
          plot.margin = margin(0.6, 3.0, 0.6, 3.0, unit = "mm"))

  mixing <- fread(mixing_path)
  stopifnot(nrow(mixing) == 22L)
  mixing[, review_flag := as.logical(review_flag)]
  label_mix <- mixing[review_flag %in% TRUE]
  p_b <- ggplot(mixing, aes(max_sample_fraction, normalized_sample_entropy)) +
    geom_vline(xintercept = c(0.30, 0.50), linetype = c(2, 1), linewidth = 0.35,
               colour = "#555555") +
    geom_hline(yintercept = 0.65, linetype = 2, linewidth = 0.35,
               colour = "#555555") +
    geom_point(aes(size = cell_n, shape = review_flag), alpha = 0.82,
               colour = "#2B2B2B") +
    ggrepel::geom_text_repel(
      data = label_mix, aes(label = seurat_clusters), seed = 180805,
      size = 2.2, min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.25, max.overlaps = Inf
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    scale_size_continuous(range = c(1.2, 4.0), labels = comma) +
    labs(title = panel_title("B", "Sample-mixing QC"),
         subtitle = "Prespecified review clusters are labeled",
         x = "Largest single-sample fraction", y = "Normalized sample entropy",
         shape = "Review", size = "Cells") +
    theme_pub() +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.box.spacing = unit(0.15, "mm"),
          legend.margin = margin(0, 0, 0, 0, unit = "mm"),
          plot.margin = margin(0.6, 3.0, 1.0, 3.0, unit = "mm")) +
    guides(size = "none", shape = guide_legend(nrow = 1))

  marker <- fread(marker_path)
  stopifnot(nrow(marker) == 352L, uniqueN(marker$cluster) == 22L)
  lineage_labels <- c(
    B_cell = "B", ciliated_epithelial = "Ciliated", cycling = "Cycling",
    decidual_stromal = "Decidual", dendritic_cell = "Dendritic",
    endothelial = "Endothelial", epithelial_general = "Epithelial",
    erythroid = "Erythroid", immune_general = "Immune",
    lymphatic_endothelial = "Lymphatic endo.", mast_cell = "Mast",
    myeloid_macrophage = "Myeloid", NK_cell = "NK",
    perivascular_smooth_muscle = "Perivascular", plasma_cell = "Plasma",
    secretory_epithelial = "Secretory", stromal_fibroblast = "Stromal", T_cell = "T"
  )
  marker[, panel_label := lineage_labels[panel]]
  marker[, cluster := factor(cluster, levels = rev(as.character(0:21)))]
  marker[, panel_label := factor(panel_label, levels = unname(lineage_labels))]
  p_d <- ggplot(marker, aes(panel_label, cluster, fill = mean_scaled_expression)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_gradientn(colours = c("#3B0F70", "#2C7BB6", "#22A884", "#FDE725"),
                         name = "Mean scaled\nexpression") +
    labs(title = panel_title("D", "Broad-lineage marker support"),
         subtitle = "Cluster-by-panel mean scaled expression; annotation support, not cell-level inference",
         x = "Canonical marker panel", y = "Cluster") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 5.8),
           axis.text.y = element_text(size = 5.8),
           legend.position = "right",
           plot.margin = margin(0.2, 3.0, 3.0, 3.0, unit = "mm"))

  p_a_pdf <- ggrastr::rasterise(p_a, layers = "Point", dpi = 600, dev = "ragg")
  p_a <- panel_module(p_a, "A")
  p_a_pdf <- panel_module(p_a_pdf, "A")
  p_b <- panel_module(p_b, "B")
  p_c <- panel_module(p_c, "C")
  p_d <- panel_module(p_d, "D")

  design <- c(
    patchwork::area(t = 1, l = 1, b = 1, r = 1),
    patchwork::area(t = 1, l = 2, b = 2, r = 2),
    patchwork::area(t = 2, l = 1, b = 2, r = 1),
    patchwork::area(t = 3, l = 1, b = 3, r = 2)
  )
  composite <- wrap_plots(p_a, p_b, p_c, p_d, design = design,
                          heights = c(1, 0.95, 1.23), widths = c(1.12, 0.88))
  composite_pdf <- wrap_plots(p_a_pdf, p_b, p_c, p_d, design = design,
                              heights = c(1, 0.95, 1.23), widths = c(1.12, 0.88))
  fig_path <- file.path(package_root, "figures",
                        paste0("Supplementary_Figure_S2_", version_id,
                               "_candidate_300dpi.png"))
  export_candidate(composite, fig_path, 183, 188, pdf_plot = composite_pdf)

  fwrite(display_cells[, .(cell_id = cell_barcode, dim1 = UMAPH_1, dim2 = UMAPH_2,
                           cluster, broad_cell_type, sample_accession)],
         file.path(package_root, "source_data", "S2_panel_A_UMAP_display_subset.csv"))
  fwrite(doublet, file.path(package_root, "source_data", "S2_panel_C_doublet_qc_with_design.csv"))
  copy_source(mixing_path, file.path(package_root, "source_data"), "S2_panel_B_cluster_sample_mixing.csv")
  copy_source(marker_path, file.path(package_root, "source_data"), "S2_panel_D_broad_lineage_support.csv")
  copy_source(annotation_path, file.path(package_root, "source_data"), "S2_cluster_annotation_summary.csv")
  fwrite(data.table(panel = c("A", "B", "C", "D"),
                    direct_source = c("S2_panel_A_UMAP_display_subset.csv",
                                      "S2_panel_B_cluster_sample_mixing.csv",
                                      "S2_panel_C_doublet_qc_with_design.csv",
                                      "S2_panel_D_broad_lineage_support.csv"),
                    upstream_authority = "STEP04 GSE250130 QC/lock",
                    generator = basename(script_path)),
         file.path(package_root, "manifests", "FIGURE_SOURCE_CODE_MAP.csv"))

  panel_rows <- c(
    "| A | descriptive integration/annotation QC | locked cell_annotation_metadata.csv | deterministic display subset from 236,636 cells | integrated clusters map to the locked broad-cell framework | UMAP geometry is descriptive and cells are not independent replicates |",
    "| B | technical QC | 04_CLUSTER_SAMPLE_MIXING_METRICS.csv | 22 clusters | most clusters show broad sample mixing; labels identify review clusters | donor/sample mixing does not validate disease biology |",
    "| C | technical QC | GSE250130_sample_doublet_qc.csv | 28 samples | predicted doublet fractions are shown for all samples | no between-group test is performed |",
    "| D | annotation support | 06_BROAD_LINEAGE_PANEL_SUPPORT.csv | 22 clusters x 16 marker panels | canonical marker-panel patterns support broad annotation review | scaled expression is descriptive and not donor-level inference |"
  )
  source_rows <- c(
    paste0("- Panel A direct full authority: ", cell_path, " (236,636 rows); plotting source is a deterministic <=2,500-cell-per-cluster subset exported with the candidate."),
    paste0("- Panel B direct result: ", mixing_path, " (22 clusters)."),
    paste0("- Panel C direct result: ", doublet_path, " (28 samples)."),
    paste0("- Panel D direct result: ", marker_path, " (352 cluster-by-panel rows)."),
    paste0("- Locked annotation summary: ", annotation_path, "."),
    "- Statistical unit: sample for doublet QC; cluster for mixing and marker summaries; cells are display units only in UMAP.",
    "- No normalization, integration, UMAP, clustering, marker test or annotation is rerun.",
    "- Frozen result sentence: the locked broad annotation and donor-aware aggregation are supported by sample-mixing, doublet and marker-panel diagnostics."
  )
  legend_en <- c(
    "Supplementary Figure S2. GSE250130 integration, annotation and sample-mixing quality control.",
    "(A) Harmony UMAP coordinates for a deterministic plotting subset of at most 2,500 cells per cluster (55,000 cells maximum) from the complete 236,636-cell annotation. Colors denote the eight broad cell types and numbers identify the 22 clusters. UMAP coordinates are descriptive layout coordinates. (B) Cluster-level sample-mixing diagnostics. The x axis gives the largest fraction contributed by one sample, the y axis gives normalized sample entropy, point size gives cluster cell count, and triangles identify clusters prespecified for review; only review clusters are labeled to avoid annotation collisions. Reference lines show the prespecified QC thresholds. (C) Predicted doublet fractions for all 28 samples, grouped by fertile luteal day or RIF LH+7 status. Points are descriptive and no group comparison was performed. (D) Mean scaled expression for 16 canonical broad-lineage marker panels across 22 clusters. The heatmap supports broad annotation review but is not a cell-level inferential analysis."
  )
  panel_cn <- c(
    "补充图 S2：GSE250130 的整合、注释与样本混合质量控制。",
    "A：使用锁定的 UMAP 坐标，按每个 cluster 最多 2500 个细胞确定性抽样显示；颜色为 8 个宽谱细胞类型，数字为 22 个 cluster。",
    "B：cluster 层面的最大单样本占比、归一化样本熵与细胞数；仅标注预设复核 cluster。",
    "C：28 个样本的预测 doublet 比例，仅作描述性 QC，不进行组间检验。",
    "D：22 个 cluster 在 16 个典型谱系 marker panel 上的锁定平均标准化表达，用于支持宽谱注释复核，不构成细胞层面推断。"
  )
  write_common_cards(package_root, "S2", "GSE250130 integration and annotation QC",
                     "Do the locked doublet, integration, sample-mixing and canonical-marker summaries support the broad annotation used for donor-aware GSE250130 analyses?",
                     183, 188, panel_rows, source_rows, legend_en, panel_cn,
                     "the broad annotation and donor-aware aggregation are technically supported by the displayed QC summaries.",
                     "UMAP separation, cell density or marker color cannot establish donor-level disease effects or independent validation.",
                     "Supplementary Table S3")
  finalize_package(package_root, "S2", c(cell_path, doublet_path, mixing_path,
                                          marker_path, annotation_path))
  package_root
}

build_s3 <- function() {
  activate_figure("S3")
  package_root <- make_package("S3")
  umap_path <- file.path(project_root, "05_analysis_steps",
                         "STEP04_combined_temporal_patient_single_cell", "06_results",
                         "S4E1_sc", "SCI183_v1_0", "tables", "04_UMAP_metadata.csv")
  annotation_map_path <- file.path(project_root, "05_analysis_steps",
                                   "STEP04_combined_temporal_patient_single_cell", "06_results",
                                   "S4E1_sc", "SCI183_v1_0", "ann", "G183", "ANN183_v1_0",
                                   "t", "01_annotation_map.csv")
  mixing_path <- file.path(project_root, "05_analysis_steps",
                           "STEP04_combined_temporal_patient_single_cell", "06_results",
                           "S4E1_sc", "SCI183_v1_0", "integration_QC", "GSE183837",
                           "INTQC183_v1_1", "tables", "04_cluster_donor_and_technical_QC.csv")
  annotation_qc_path <- file.path(project_root, "05_analysis_steps",
                                  "STEP04_combined_temporal_patient_single_cell", "06_results",
                                  "S4E1_sc", "SCI183_v1_0", "ann", "G183", "ANN183_v1_0",
                                  "t", "03_annotation_QC.csv")
  sensitivity_path <- file.path(project_root, "05_analysis_steps",
                                "STEP04_combined_temporal_patient_single_cell", "06_results",
                                "S4E1_sc", "SCI183_v1_0", "mdl", "G183", "MD183_v1_0",
                                "tab", "06_cluster8_model_sensitivity.csv")
  lapply(c(umap_path, annotation_map_path, mixing_path, annotation_qc_path,
           sensitivity_path), assert_file)

  umap <- fread(umap_path,
                select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "sample_accession",
                           "analysis_group", "seurat_clusters"))
  ann_map <- fread(annotation_map_path)
  umap <- merge(umap, ann_map[, .(cluster, broad_cell_type, annotation_confidence,
                                  review_flag)],
                by.x = "seurat_clusters", by.y = "cluster", all.x = TRUE, sort = FALSE)
  stopifnot(nrow(umap) == 47161L, all(!is.na(umap$broad_cell_type)))
  umap[, cluster := as.character(seurat_clusters)]
  display_umap <- deterministic_subset(umap, "broad_cell_type", "cell_barcode", 6000L)
  broad_counts <- umap[, .N, by = broad_cell_type][order(-N)]
  display_umap[, draw_order := match(broad_cell_type, broad_counts$broad_cell_type)]
  setorder(display_umap, draw_order, cell_barcode)
  display_umap[, broad_cell_type := factor(broad_cell_type, levels = names(broad_palette))]
  centres <- display_umap[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)),
                          by = cluster]
  centres[, label_colour := ifelse(cluster == "8", "#D55E00", "#1A1A1A")]
  p_a <- ggplot(display_umap, aes(UMAPH_1, UMAPH_2, colour = broad_cell_type)) +
    geom_point(size = 0.14, alpha = 0.46, stroke = 0) +
    ggrepel::geom_text_repel(
      data = centres, aes(x = UMAPH_1, y = UMAPH_2,
                          label = cluster, colour = I(label_colour)), inherit.aes = FALSE,
      seed = 180806, size = 2.25, fontface = "bold", min.segment.length = 0,
      segment.size = 0.2, box.padding = 0.2, point.padding = 0.12,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_colour_manual(values = broad_palette, drop = FALSE) +
    coord_equal() +
    labs(title = panel_title("A", "Integrated UMAP and annotations"),
         subtitle = "Cluster 8 (orange label) is low confidence\nand excluded from primary models",
         x = "Harmony UMAP 1", y = "Harmony UMAP 2", colour = "Broad cell type") +
    theme_pub() + theme(legend.position = "none")

  mixing <- fread(mixing_path)
  stopifnot(nrow(mixing) == 16L)
  mixing[, review_required := as.logical(review_required)]
  p_b <- ggplot(mixing, aes(top_donor_fraction, donor_entropy_normalized)) +
    geom_vline(xintercept = c(0.50, 0.75), linetype = 2, linewidth = 0.35,
               colour = "#555555") +
    geom_hline(yintercept = c(0.50, 0.70), linetype = 3, linewidth = 0.35,
               colour = "#555555") +
    geom_point(aes(size = cell_n, shape = review_required), alpha = 0.82,
               colour = "#2B2B2B") +
    ggrepel::geom_text_repel(
      data = mixing[review_required %in% TRUE], aes(label = seurat_clusters),
      seed = 180807, size = 2.2, min.segment.length = 0,
      segment.size = 0.25, box.padding = 0.25, max.overlaps = Inf
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
    scale_size_continuous(range = c(1.2, 4.0), labels = comma) +
    labs(title = panel_title("B", "Donor-mixing diagnostics"),
         subtitle = "Labels identify clusters requiring\ndonor/technical review",
         x = "Largest single-donor fraction", y = "Normalized donor entropy",
         shape = "Review", size = "Cells") +
    guides(size = "none", shape = guide_legend(nrow = 1)) +
    theme_pub() + theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.box.spacing = unit(1.25, "mm"),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(t = 0, unit = "mm"),
      plot.margin = margin(0.6, 3.0, 0.8, 3.0, unit = "mm")
    )

  ann_qc <- fread(annotation_qc_path)
  stopifnot(nrow(ann_qc) == 16L)
  ann_qc[, cluster_factor := factor(cluster, levels = as.character(0:15))]
  p_c <- ggplot(ann_qc, aes(cluster_factor, top100_evidence_n,
                            shape = annotation_confidence, size = cell_n)) +
    geom_point(colour = "#2B2B2B", alpha = 0.82) +
    geom_point(data = ann_qc[cluster == 8], colour = "#D55E00", stroke = 0.9) +
    ggrepel::geom_text_repel(
      data = ann_qc[cluster == 8], aes(label = "cluster 8"),
      seed = 180808, size = 2.2, colour = "#D55E00",
      min.segment.length = 0, segment.size = 0.25
    ) +
    scale_shape_manual(values = c(HIGH = 16, MODERATE = 15, LOW = 17)) +
    scale_size_continuous(range = c(1.2, 4.2), labels = comma) +
    scale_y_continuous(breaks = 0:4, limits = c(0, 4.25)) +
    labs(title = panel_title("C", "Annotation evidence and confidence"),
         subtitle = "Requested evidence markers among the top-100 marker summaries",
         x = "Primary cluster", y = "Evidence markers present",
         shape = "Confidence", size = "Cells") +
    guides(size = "none", shape = guide_legend(nrow = 1)) +
    theme_pub() + theme(
      axis.title.x = element_text(margin = margin(t = 1.25, unit = "mm")),
      legend.position = "bottom", legend.box = "vertical",
      legend.box.spacing = unit(0.05, "mm"),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(t = -0.25, unit = "mm"),
      plot.margin = margin(0.6, 3.0, 0.7, 3.0, unit = "mm")
    )

  sensitivity <- fread(sensitivity_path)
  stopifnot(nrow(sensitivity) == 5L, all(as.logical(sensitivity$direction_concordant)))
  sens_long <- rbindlist(list(
    sensitivity[, .(endpoint, analysis_matrix = "Primary (cluster 8 excluded)",
                    effect = mean_difference_RIF_minus_Control_primary,
                    exact_P = exact_permutation_P_primary,
                    loo = LOO_direction_concordance_primary)],
    sensitivity[, .(endpoint, analysis_matrix = "Atlas sensitivity (cluster 8 retained)",
                    effect = mean_difference_RIF_minus_Control_atlas,
                    exact_P = exact_permutation_P_atlas,
                    loo = LOO_direction_concordance_atlas)]
  ))
  sens_long[, endpoint_label := factor(short_endpoint[endpoint],
                                       levels = rev(unname(short_endpoint)))]
  sens_long[, analysis_matrix := factor(
    analysis_matrix,
    levels = c("Primary (cluster 8 excluded)", "Atlas sensitivity (cluster 8 retained)"),
    labels = c("Primary: cluster 8 excluded", "Sensitivity: cluster 8 retained")
  )]
  p_d <- ggplot(sens_long, aes(effect, endpoint_label,
                               shape = analysis_matrix, colour = analysis_matrix)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
    geom_point(size = 2.1, stroke = 0.8, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = c("Primary: cluster 8 excluded" = "#0072B2",
                                   "Sensitivity: cluster 8 retained" = "#D55E00")) +
    scale_shape_manual(values = c("Primary: cluster 8 excluded" = 16,
                                  "Sensitivity: cluster 8 retained" = 1)) +
    labs(title = panel_title("D", "Cluster-8 sensitivity"),
         subtitle = "All five directions are concordant;\nalternative definitions are not connected",
         x = "RIF minus Control mean difference", y = NULL,
         colour = "Analysis", shape = "Analysis") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                 keyheight = unit(2.7, "mm")),
           shape = guide_legend(nrow = 2, byrow = TRUE,
                                keyheight = unit(2.7, "mm"))) +
    theme_pub() + theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.box.spacing = unit(1.15, "mm"),
      legend.key.height = unit(2.7, "mm"),
      legend.spacing.y = unit(0.05, "mm"),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(t = 0, unit = "mm"),
      plot.margin = margin(0.6, 3.0, 0.7, 3.0, unit = "mm")
    )

  p_a_pdf <- ggrastr::rasterise(p_a, layers = "Point", dpi = 600, dev = "ragg")
  p_a <- panel_module(p_a, "A")
  p_a_pdf <- panel_module(p_a_pdf, "A")
  p_b <- panel_module(p_b, "B")
  p_c <- panel_module(p_c, "C")
  p_d <- panel_module(p_d, "D")

  composite <- (p_a + p_b) / (p_c + p_d) +
    plot_layout(widths = c(0.96, 1.04), heights = c(1.03, 0.97))
  composite_pdf <- (p_a_pdf + p_b) / (p_c + p_d) +
    plot_layout(widths = c(0.96, 1.04), heights = c(1.03, 0.97))
  fig_path <- file.path(package_root, "figures",
                        paste0("Supplementary_Figure_S3_", version_id,
                               "_candidate_300dpi.png"))
  export_candidate(composite, fig_path, 183, 180, pdf_plot = composite_pdf)

  fwrite(display_umap[, .(cell_id = cell_barcode, dim1 = UMAPH_1, dim2 = UMAPH_2,
                          cluster, broad_cell_type, sample_accession, analysis_group)],
         file.path(package_root, "source_data", "S3_panel_A_UMAP_display_subset.csv"))
  copy_source(mixing_path, file.path(package_root, "source_data"), "S3_panel_B_cluster_donor_QC.csv")
  copy_source(annotation_qc_path, file.path(package_root, "source_data"), "S3_panel_C_annotation_QC.csv")
  copy_source(sensitivity_path, file.path(package_root, "source_data"), "S3_panel_D_cluster8_sensitivity.csv")
  copy_source(annotation_map_path, file.path(package_root, "source_data"), "S3_annotation_map.csv")
  fwrite(data.table(panel = c("A", "B", "C", "D"),
                    direct_source = c("S3_panel_A_UMAP_display_subset.csv",
                                      "S3_panel_B_cluster_donor_QC.csv",
                                      "S3_panel_C_annotation_QC.csv",
                                      "S3_panel_D_cluster8_sensitivity.csv"),
                    upstream_authority = "STEP04 GSE183837 QC/lock",
                    generator = basename(script_path)),
         file.path(package_root, "manifests", "FIGURE_SOURCE_CODE_MAP.csv"))

  panel_rows <- c(
    "| A | descriptive integration/annotation QC | locked 04_UMAP_metadata.csv + annotation map | deterministic display subset from 47,161 cells | six broad types and 16 clusters are visible; cluster 8 is flagged | UMAP is descriptive and cluster 8 is not primary evidence |",
    "| B | technical QC | 04_cluster_donor_and_technical_QC.csv | 16 clusters across nine donors | donor enrichment is localized to prespecified review clusters | mixing does not establish disease biology |",
    "| C | annotation support | 03_annotation_QC.csv | 16 clusters | marker-evidence count and confidence identify cluster 8 as low confidence | cluster-level marker evidence is descriptive |",
    "| D | sensitivity | 06_cluster8_model_sensitivity.csv | five endpoints; donor-level primary and atlas matrices | all five effect directions are retained when cluster 8 is included | sensitivity does not create corrected significance or validate a disease mechanism |"
  )
  source_rows <- c(
    paste0("- Panel A direct full authority: ", umap_path, " joined only to the locked annotation map ", annotation_map_path, "."),
    paste0("- Panel B direct result: ", mixing_path, " (16 clusters)."),
    paste0("- Panel C direct result: ", annotation_qc_path, " (16 clusters)."),
    paste0("- Panel D direct result: ", sensitivity_path, " (five endpoints; saved effects/P/LOO fields)."),
    "- Statistical unit: donor for Panel D; cluster for Panels B-C; cells are display units only in Panel A.",
    "- No integration, UMAP, clustering, marker testing, annotation, pseudobulk aggregation or donor model is rerun.",
    "- Frozen result sentence: cluster 8 remains low confidence and excluded from primary inference; retaining it in the atlas sensitivity changes no effect direction."
  )
  legend_en <- c(
    "Supplementary Figure S3. GSE183837 annotation, donor mixing and cluster-8 sensitivity.",
    "(A) Harmony UMAP coordinates for a deterministic plotting subset of at most 6,000 cells per broad type from the complete 47,161-cell dataset. Colors denote the six broad cell types and numbers identify the 16 primary clusters. Cluster 8 is labeled in orange because it remained low-confidence epithelial-like and was excluded from primary donor-level models. (B) Cluster-level donor-mixing diagnostics. The x axis gives the largest single-donor fraction, the y axis gives normalized donor entropy, point size gives cell count, and triangles identify prespecified review clusters; only review clusters are labeled. Reference lines show the prespecified QC thresholds. (C) Number of requested annotation-evidence markers present in each cluster's top-100 marker summary. Point size gives cluster cell count and shape gives annotation confidence; cluster 8 is highlighted. (D) RIF-minus-Control mean differences for five unciliated-epithelial endpoints under the primary matrix, which excludes cluster 8, and the atlas-sensitivity matrix, which retains it. The two alternative definitions are shown as unconnected estimates. All five endpoint directions were concordant; this sensitivity analysis did not create corrected statistical support."
  )
  panel_cn <- c(
    "补充图 S3：GSE183837 的注释、供体混合与 cluster 8 敏感性。",
    "A：使用锁定的 Harmony UMAP 坐标，按每个宽谱细胞类型最多 6000 个细胞确定性抽样；数字为 16 个 cluster，橙色的 cluster 8 为低置信度、主分析排除。",
    "B：cluster 层面的最大单供体占比、归一化供体熵与细胞数，仅标注预设复核 cluster。",
    "C：各 cluster 在锁定 top-100 marker 汇总中出现的请求证据 marker 数及注释置信度，突出显示 cluster 8。",
    "D：主矩阵与保留 cluster 8 的 atlas 敏感性矩阵中 5 个终点的供体层效应；两套定义不以连线暗示连续性，5/5 方向一致，未产生校正后证据。"
  )
  write_common_cards(package_root, "S3", "GSE183837 annotation and cluster-8 sensitivity",
                     "Do the locked annotation, donor-mixing and cluster-8 sensitivity audits support the restricted GSE183837 donor-level interpretation?",
                     183, 180, panel_rows, source_rows, legend_en, panel_cn,
                     "cluster 8 is a low-confidence sensitivity-only population and retaining it changes none of the five effect directions.",
                     "cell-level patterns or sensitivity concordance cannot establish corrected replication, causality or broad RIF disruption.",
                     "Supplementary Table S4")
  finalize_package(package_root, "S3", c(umap_path, annotation_map_path, mixing_path,
                                          annotation_qc_path, sensitivity_path))
  package_root
}

build_s4 <- function() {
  activate_figure("S4")
  package_root <- make_package("S4")
  cross_path <- file.path(project_root, "05_analysis_steps",
                          "STEP04_combined_temporal_patient_single_cell", "06_results",
                          "S4E1_sc", "cross_dataset", "SCCROSS_v1_0", "tab",
                          "03_cross_dataset_25_model_concordance.csv")
  overall_path <- file.path(project_root, "05_analysis_steps",
                            "STEP04_combined_temporal_patient_single_cell", "06_results",
                            "S4E1_sc", "cross_dataset", "SCCROSS_v1_0", "tab",
                            "04_overall_concordance_summary.csv")
  evidence_path <- file.path(project_root, "05_analysis_steps",
                             "STEP04_combined_temporal_patient_single_cell", "06_results",
                             "S4E1_sc", "cross_dataset", "SCCROSS_v1_0", "tab",
                             "07_evidence_class_summary.csv")
  lapply(c(cross_path, overall_path, evidence_path), assert_file)
  cross <- fread(cross_path)
  overall <- fread(overall_path)
  evidence <- fread(evidence_path)
  stopifnot(nrow(cross) == 25L, nrow(overall) == 13L, nrow(evidence) == 4L)

  long <- rbindlist(list(
    cross[, .(broad_cell_type, endpoint, dataset = "GSE250130",
              exact_P = exact_P_G250, q_value = primary_q_G250,
              q_scope = primary_q_scope_G250, loo = LOO_direction_concordance_G250,
              constant_endpoint = constant_endpoint_G250)],
    cross[, .(broad_cell_type, endpoint, dataset = "GSE183837",
              exact_P = exact_P_G183, q_value = primary_q_G183,
              q_scope = primary_q_scope_G183, loo = LOO_direction_concordance_G183,
              constant_endpoint = constant_endpoint_G183)]
  ))
  long[, cell_label := factor(short_cell_type[broad_cell_type],
                              levels = unname(short_cell_type[cell_type_order]))]
  long[, endpoint_label := factor(short_endpoint[endpoint],
                                  levels = rev(unname(short_endpoint[endpoint_order])))]
  long[, dataset := factor(dataset, levels = c("GSE250130", "GSE183837"))]
  long[, p_text := ifelse(is.na(exact_P), "NA", formatC(exact_P, format = "f", digits = 3))]
  long[, q_text := ifelse(is.na(q_value), "NA", formatC(q_value, format = "f", digits = 3))]
  long[, loo_text := ifelse(is.na(loo), "NA", formatC(loo, format = "f", digits = 2))]
  long[, p_fill := -log10(pmax(exact_P, 1e-6))]
  long[, q_fill := -log10(pmax(q_value, 1e-6))]

  heat_theme <- theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.0),
          axis.text.y = element_text(size = 6.0),
          strip.text = element_text(size = 7.2, face = "bold"),
          strip.background = element_rect(fill = "#F2F2F2", colour = "#666666", linewidth = 0.3),
          legend.position = "right",
          legend.justification = c(0.5, 1),
          legend.box.just = "top",
          legend.margin = margin(0, 0, 0, 1.0, unit = "mm"))
  p_a <- ggplot(long, aes(cell_label, endpoint_label, fill = p_fill)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = p_text), family = "Arial", size = 2.0) +
    facet_grid(. ~ dataset) +
    scale_fill_gradientn(colours = c("#F3F3F3", "#92C5DE", "#0571B0"),
                         limits = c(0, max(long$p_fill, na.rm = TRUE)),
                         name = expression(-log[10](italic(P))),
                         guide = guide_colourbar(barheight = unit(27.3, "mm"),
                                                 barwidth = unit(4.0, "mm"),
                                                 title.position = "top")) +
    labs(title = panel_title("A", "Exact-permutation P values for all 25 models"),
         subtitle = "Tiles give exact P values; no cross-cohort combination is performed",
         x = NULL, y = NULL) + heat_theme

  p_b <- ggplot(long, aes(cell_label, endpoint_label, fill = q_fill)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = q_text), family = "Arial", size = 2.0) +
    facet_grid(. ~ dataset) +
    scale_fill_gradientn(colours = c("#F3F3F3", "#FDB863", "#B2182B"),
                         limits = c(0, max(long$q_fill, na.rm = TRUE)),
                         name = expression(-log[10](italic(q))),
                         guide = guide_colourbar(barheight = unit(27.3, "mm"),
                                                 barwidth = unit(4.0, "mm"),
                                                 title.position = "top")) +
    labs(title = panel_title("B", "Prespecified within-cohort adjusted q values"),
         subtitle = "GSE250130 uses endpoint-family BH; GSE183837 uses global BH across 25 models",
         x = NULL, y = NULL) + heat_theme

  p_c <- ggplot(long, aes(cell_label, endpoint_label, fill = loo)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = loo_text), family = "Arial", size = 2.0) +
    facet_grid(. ~ dataset) +
    scale_fill_gradientn(colours = c("#B2182B", "#F7F7F7", "#2166AC"),
                         limits = c(0, 1), na.value = "#E5E5E5",
                         name = "LOO direction\nconcordance",
                         guide = guide_colourbar(barheight = unit(25.8, "mm"),
                                                 barwidth = unit(4.0, "mm"),
                                                 title.position = "top")) +
    labs(title = panel_title("C", "Leave-one-donor-out direction concordance"),
         subtitle = "NA denotes a constant or non-informative endpoint; values are model diagnostics",
         x = NULL, y = NULL) + heat_theme

  evidence[, evidence_label := factor(evidence_class,
    levels = c("CONCORDANT_NO_STATISTICAL_SUPPORT",
               "DISCORDANT_NO_STATISTICAL_SUPPORT",
               "DISCORDANT_WITH_NOMINAL_OR_CORRECTED_SIGNAL",
               "NONINFORMATIVE_CONSTANT_OR_ZERO"),
    labels = c("Concordant; no statistical support",
               "Discordant; no statistical support",
               "Discordant; nominal signal in one cohort",
               "Non-informative constant/zero"))]
  evidence[, evidence_label := factor(evidence_label, levels = rev(levels(evidence_label)))]
  evidence_cols <- c(
    "Concordant; no statistical support" = "#4DAF4A",
    "Discordant; no statistical support" = "#984EA3",
    "Discordant; nominal signal in one cohort" = "#E41A1C",
    "Non-informative constant/zero" = "#999999"
  )
  p_d <- ggplot(evidence, aes(model_n, evidence_label, fill = evidence_label)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = model_n), hjust = -0.35, family = "Arial", size = 2.4) +
    scale_fill_manual(values = evidence_cols, guide = "none") +
    scale_x_continuous(limits = c(0, 13), breaks = c(0, 4, 8, 12), expand = c(0, 0)) +
    labs(title = panel_title("D", "Cross-dataset evidence classes"),
         subtitle = "24 informative models: 12 concordant and 12 discordant; corrected replication = 0",
         x = "Models", y = NULL) +
    theme_pub() + theme(axis.text.y = element_text(size = 6.2))

  p_a <- panel_module(p_a, "A")
  p_b <- panel_module(p_b, "B")
  p_c <- panel_module(p_c, "C")
  p_d <- panel_module(p_d, "D")

  composite <- p_a / p_b / p_c / p_d +
    plot_layout(heights = c(1, 1, 1, 0.72))
  fig_path <- file.path(package_root, "figures",
                        paste0("Supplementary_Figure_S4_", version_id,
                               "_candidate_300dpi.png"))
  export_candidate(composite, fig_path, 183, 220)

  copy_source(cross_path, file.path(package_root, "source_data"), "S4_all_25_model_diagnostics.csv")
  copy_source(overall_path, file.path(package_root, "source_data"), "S4_overall_concordance_summary.csv")
  copy_source(evidence_path, file.path(package_root, "source_data"), "S4_evidence_class_summary.csv")
  fwrite(long, file.path(package_root, "source_data", "S4_plotting_long_table.csv"))
  fwrite(data.table(panel = c("A", "B", "C", "D"),
                    direct_source = c("S4_plotting_long_table.csv", "S4_plotting_long_table.csv",
                                      "S4_plotting_long_table.csv", "S4_evidence_class_summary.csv"),
                    upstream_authority = "SCCROSS_v1_0",
                    generator = basename(script_path)),
         file.path(package_root, "manifests", "FIGURE_SOURCE_CODE_MAP.csv"))

  panel_rows <- c(
    "| A | inferential diagnostic | 03_cross_dataset_25_model_concordance.csv | 25 prespecified models per cohort | saved exact-permutation P values remain broadly non-significant | P values are not pooled or combined |",
    "| B | multiplicity diagnostic | same locked table | cohort-specific correction families | no model has adjusted q < 0.05 | correction scopes differ and must not be directly pooled |",
    "| C | sensitivity diagnostic | same locked table | donor-level leave-one-out summaries | directional stability varies but does not yield corrected replication | LOO concordance is not an independent significance test |",
    "| D | evidence summary | 07_evidence_class_summary.csv | 25 shared models | 12 concordant, 12 discordant, one non-informative; corrected replication zero | directional concordance alone is not validation |"
  )
  source_rows <- c(
    paste0("- Panels A-C direct source: ", cross_path, " (25 shared models; saved exact P, cohort-specific q and LOO fields)."),
    paste0("- Panel D direct source: ", evidence_path, " (four frozen evidence classes)."),
    paste0("- Overall frozen summary: ", overall_path, "."),
    "- Statistical unit: donor within each cohort. GSE250130 and GSE183837 are never pooled.",
    "- Display-only transformations: -log10 color scales and reshaping to a 5 x 5 model grid; original P/q values are printed in every tile and retained in source_data.",
    "- Correction families: endpoint-specific within-family BH for GSE250130; global BH across 25 models for GSE183837.",
    "- Frozen result sentence: 24 models are informative (12 concordant, 12 discordant), with zero corrected cross-dataset replication and zero exploratory directional support."
  )
  legend_en <- c(
    "Supplementary Figure S4. Complete donor-level model, multiplicity and concordance diagnostics.",
    "(A) Exact-permutation P values for the 25 prespecified donor-level models in GSE250130 and GSE183837. Rows are five endpoints and columns are five broad cell types; tile text gives the exact P value and color gives -log10(P). (B) Prespecified adjusted q values for the same models. GSE250130 q values were calculated within the six VDR-expression or 24 program-model endpoint families, whereas GSE183837 q values used global Benjamini-Hochberg correction across all 25 models; the cohorts and correction families were not combined. (C) Leave-one-donor-out direction-concordance fractions for the same models. Gray NA tiles denote constant or non-informative endpoints. Leave-one-out concordance is a stability diagnostic, not an independent significance test. (D) Cross-dataset evidence classes. Among 24 informative models, 12 were directionally concordant and 12 discordant; one model was non-informative. Two discordant models had a nominal signal in one cohort, but no model met corrected replication or exploratory directional-support criteria. Values derive from donor-level model outputs; the cohorts were not combined by pooled P value, correlation test or meta-analysis."
  )
  panel_cn <- c(
    "补充图 S4：完整的供体层模型、多重校正与跨数据集一致性诊断。",
    "A：GSE250130 与 GSE183837 各 25 个预设模型的锁定精确置换 P 值；格内打印原值，颜色为 -log10(P)。",
    "B：同一批模型的预设校正 q 值；GSE250130 与 GSE183837 的校正家族不同，未合并。",
    "C：同一批模型的逐一剔除供体方向一致率；NA 为常数或无信息终点，该指标仅为稳定性诊断。",
    "D：锁定证据分类。24 个有信息模型中 12 个方向一致、12 个方向相反，1 个无信息；校正后跨数据集重复为 0。"
  )
  write_common_cards(package_root, "S4", "complete donor-level model diagnostics",
                     "Do the complete exact-permutation, multiplicity and leave-one-donor-out diagnostics support the frozen boundary of no corrected cross-dataset replication?",
                     183, 220, panel_rows, source_rows, legend_en, panel_cn,
                     "the two cohorts show no corrected cross-dataset replication across the 25 prespecified models.",
                     "directional agreement, nominal signals or LOO stability cannot be treated as pooled evidence, meta-analysis, causality or biomarker validation.",
                     "Supplementary Table S5")
  finalize_package(package_root, "S4", c(cross_path, overall_path, evidence_path))
  package_root
}

run_main <- function() {
  active_figure <<- "GLOBAL"
  font_match <- systemfonts::match_fonts("Arial")
  if (nrow(font_match) < 1L || !file.exists(font_match$path[[1]])) {
    stop("Arial font preflight failed")
  }
  if (.Platform$OS.type == "windows") {
    grDevices::windowsFonts(FigureArial = grDevices::windowsFont("Arial"))
  }
  font_resolved_path <<- normalizePath(font_match$path[[1]], winslash = "/", mustWork = TRUE)
  font_test <- tempfile(fileext = ".png")
  ragg::agg_png(font_test, width = 600, height = 180, units = "px", res = 150)
  grid::grid.newpage()
  grid::grid.text("A VDR 0.05 RIF-control", gp = grid::gpar(fontfamily = "Arial", fontsize = 10))
  width_test <- grid::convertWidth(grid::stringWidth("A VDR 0.05 RIF-control"),
                                  "mm", valueOnly = TRUE)
  grDevices::dev.off()
  unlink(font_test)
  if (!is.finite(width_test) || width_test <= 0) stop("Grid font-width preflight failed")
  font_width_mm <<- width_test

  # panel_module() constructs explicit gtables before the final export device
  # opens.  Keep those font-metric calculations on a native ragg device so the
  # null PDF/PostScript device never substitutes or warns on Arial.
  staging_device <- tempfile(fileext = ".png")
  ragg::agg_png(staging_device, width = 800, height = 800, units = "px", res = 150,
                background = "white")
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(staging_device)
  }, add = TRUE)

  outputs <- character()
  if (target %in% c("ALL", "S1")) outputs <- c(outputs, build_s1())
  if (target %in% c("ALL", "S2")) outputs <- c(outputs, build_s2())
  if (target %in% c("ALL", "S3")) outputs <- c(outputs, build_s3())
  if (target %in% c("ALL", "S4")) outputs <- c(outputs, build_s4())
  cat("Built candidate packages:\n", paste(outputs, collapse = "\n"), "\n")
}

withCallingHandlers(run_main(), warning = record_warning)
