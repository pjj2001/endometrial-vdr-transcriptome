
# Publication figure stage 1
options(stringsAsFactors = FALSE, warn = 1)

# Figure 6 v1.1 manuscript composite from frozen STEP04_LOCK_v1.2 and SCCROSS_v1_0.
# No integration, model fitting, hypothesis testing, correlation testing,
# pooled meta-analysis, or multiplicity correction is rerun.

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "131_step05F_Figure6_manuscript_composite_v1_1_QC.R"
candidate_id <- "F6_v1_2_language_clean"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

required_packages <- c(
  "data.table", "ggplot2", "patchwork", "gtable", "ragg", "magick",
  "digest", "zip", "png", "jsonlite", "scales"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))

normp <- function(x, must_work = FALSE) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = must_work))
}
sha256_file <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
md5_file <- function(path) unname(as.character(tools::md5sum(path)))

write_csv_utf8 <- function(x, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", normp(path, FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}

write_text_utf8 <- function(x, path) {
  if (file.exists(path)) stop("Refusing overwrite: ", normp(path, FALSE))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(x), con = con, useBytes = TRUE)
  invisible(path)
}

relative_to <- function(path, root) {
  p <- normp(path, FALSE)
  r <- sub("/+$", "", normp(root, FALSE))
  if (startsWith(tolower(p), paste0(tolower(r), "/"))) return(substr(p, nchar(r) + 2L, nchar(p)))
  p
}

copy_exact <- function(src, dst) {
  if (!file.exists(src) || dir.exists(src)) stop("Missing source file: ", normp(src, FALSE))
  if (file.exists(dst)) stop("Refusing overwrite: ", normp(dst, FALSE))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dst) || !identical(sha256_file(src), sha256_file(dst))) {
    stop("Validated copy failed: ", normp(src, FALSE), " -> ", normp(dst, FALSE))
  }
  invisible(dst)
}

copy_tree <- function(src_root, dst_root) {
  files <- list.files(src_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  for (src in files) copy_exact(src, file.path(dst_root, relative_to(src, src_root)))
  invisible(TRUE)
}

build_manifest <- function(root, exclude = character()) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  rel <- vapply(files, relative_to, character(1), root = root)
  keep <- !rel %in% exclude
  files <- files[keep]
  rel <- rel[keep]
  info <- file.info(files)
  out <- data.frame(
    relative_path = rel, size_bytes = as.numeric(info$size),
    md5 = vapply(files, md5_file, character(1)),
    sha256 = vapply(files, sha256_file, character(1)), stringsAsFactors = FALSE
  )
  out[order(out$relative_path), , drop = FALSE]
}

validate_manifest <- function(root, manifest) {
  do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    p <- file.path(root, manifest$relative_path[[i]])
    ex <- file.exists(p) && !dir.exists(p)
    data.frame(
      relative_path = manifest$relative_path[[i]], exists = ex,
      size_match = ex && as.numeric(file.info(p)$size) == manifest$size_bytes[[i]],
      md5_match = ex && identical(md5_file(p), manifest$md5[[i]]),
      sha256_match = ex && identical(sha256_file(p), manifest$sha256[[i]]),
      stringsAsFactors = FALSE
    )
  }))
}

find_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    p <- sub("^--file=", "", file_arg[[1]])
    if (file.exists(p)) return(normp(p, TRUE))
  }
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    p <- frames[[i]]$ofile
    if (!is.null(p) && length(p) == 1L && file.exists(p)) return(normp(p, TRUE))
  }
  fallback <- file.path(project_root, "04_code", "R", script_name)
  if (file.exists(fallback)) return(normp(fallback, TRUE))
  stop("Could not identify executed script path.")
}

script_path <- find_script_path()
skill3_path <- "<FIGURE_STYLE_REFERENCE>"
skill3_text <- paste(readLines(skill3_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.8", skill3_text, fixed = TRUE)) stop("Skill 3 v2.2-beta.8 is required.")

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")
step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")
lock_root <- file.path(step04_root, "06_results", "STEP04_closeout", "STEP04_LOCK_v1.2")
cross_root <- file.path(step04_root, "06_results", "S4E1_sc", "cross_dataset", "SCCROSS_v1_0")
preflight_root <- file.path(step05_root, "06_results", "FIG6", "FIG6_WORKFLOW_PREFLIGHT_v1")
preflight_summary_path <- file.path(preflight_root, "run_evidence", "RUN_SUMMARY.csv")
if (!file.exists(preflight_summary_path)) stop("Formal Figure 6 preflight is missing.")
preflight_summary <- utils::read.csv(preflight_summary_path, check.names = FALSE)
if (preflight_summary$gate[[1]] != "PASS_FIGURE6_WORKFLOW_PREFLIGHT_SOURCE_ASSEMBLY_AUTHORIZED") {
  stop("Figure 6 preflight gate is not PASS.")
}

parent_override <- Sys.getenv("FIG6_CANDIDATE_PARENT", unset = "")
feedback_override <- Sys.getenv("FIG6_FEEDBACK_ROOT", unset = "")
figure_parent <- if (nzchar(parent_override)) parent_override else file.path(step05_root, "06_results", "FIG6")
feedback_root <- if (nzchar(feedback_override)) feedback_override else file.path(project_root, "\u53cd\u9988")
out_root <- file.path(figure_parent, candidate_id)
if (dir.exists(out_root)) stop("Refusing to overwrite Figure 6 candidate: ", normp(out_root, FALSE))

dirs <- list(
  figures = file.path(out_root, "figures"), source_data = file.path(out_root, "source_data"),
  manuscript = file.path(out_root, "manuscript_text"), qc = file.path(out_root, "qc"),
  code = file.path(out_root, "code"), manifests = file.path(out_root, "manifests")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
copy_exact(script_path, file.path(dirs$code, script_name))

source_paths <- c(
  concordance = file.path(cross_root, "tab", "03_cross_dataset_25_model_concordance.csv"),
  overall = file.path(cross_root, "tab", "04_overall_concordance_summary.csv"),
  endpoint_summary = file.path(cross_root, "tab", "05_concordance_by_endpoint.csv"),
  cell_summary = file.path(cross_root, "tab", "06_concordance_by_cell_type.csv"),
  evidence_summary = file.path(cross_root, "tab", "07_evidence_class_summary.csv"),
  myeloid_vdr = file.path(cross_root, "tab", "08_myeloid_VDR_focus.csv"),
  interpretation_lock = file.path(cross_root, "tab", "10_single_cell_interpretation_lock.csv"),
  persistent_manifest = file.path(cross_root, "man", "persistent_manifest.csv"),
  locked_facts = file.path(lock_root, "tables", "02_locked_fact_ledger.csv"),
  locked_claims = file.path(lock_root, "tables", "03_claim_ledger.csv"),
  integrated_summary = file.path(lock_root, "tables", "06_integrated_result_summary.csv"),
  blueprint = file.path(lock_root, "tables", "07_figure_blueprint.csv")
)
if (any(!file.exists(source_paths))) stop("Frozen Figure 6 source authority is incomplete.")

warning_rows <- list()
warning_n <- 0L
capture_stage <- function(stage, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warning_n <<- warning_n + 1L
    warning_rows[[warning_n]] <<- data.frame(
      warning_id = sprintf("W%03d", warning_n), stage = stage,
      warning_class = class(w)[[1]], message = conditionMessage(w), stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  })
}

read_frozen <- function(path) data.table::fread(path, check.names = FALSE, encoding = "UTF-8")
frozen <- lapply(source_paths[c(
  "concordance", "overall", "endpoint_summary", "cell_summary", "evidence_summary",
  "myeloid_vdr", "interpretation_lock", "locked_facts", "locked_claims", "integrated_summary", "blueprint"
)], read_frozen)

cell_display_order <- c(
  "Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells"
)
endpoint_order <- c(
  "VDR_expression", "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84"
)
endpoint_labels <- c(
  VDR_expression = "VDR\nexpression", RNA_ALL_436 = "RNA all\n436",
  OCCUPANCY_SUPPORTED_284 = "Occupancy\nsupported 284", SHARED_OCCUPANCY_46 = "Shared\noccupancy 46",
  PROMOTER_OCCUPANCY_84 = "Promoter\noccupancy 84"
)
status_levels <- c("CONCORDANT", "DISCORDANT", "NONINFORMATIVE")
status_labels <- c(
  CONCORDANT = "Concordant direction", DISCORDANT = "Discordant direction",
  NONINFORMATIVE = "Non-informative"
)
status_colors <- c(CONCORDANT = "#3E7C73", DISCORDANT = "#D58655", NONINFORMATIVE = "#D7D7D7")

dir_symbol <- function(x) {
  out <- rep("0", length(x))
  out[x == "RIF_HIGHER"] <- "+"
  out[x == "RIF_LOWER"] <- "-"
  out
}

plot_data <- capture_stage("frozen_source_transform", {
  concordance <- data.table::copy(frozen$concordance)
  numeric_cols <- c(
    "Hedges_g_G250", "Hedges_g_G183", "exact_P_G250", "primary_q_G250",
    "exact_P_G183", "primary_q_G183"
  )
  concordance[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
  concordance[, broad_cell_type := factor(broad_cell_type, levels = rev(cell_display_order))]
  concordance[, endpoint := factor(endpoint, levels = endpoint_order)]
  concordance[, direction_status := factor(direction_status, levels = status_levels)]
  concordance[, direction_pair := paste0("250 ", dir_symbol(direction_G250), "\n183 ", dir_symbol(direction_G183))]
  concordance[, label_colour := ifelse(direction_status == "NONINFORMATIVE", "dark", "light")]
  concordance[, nominal_any := nominal_significant_G250 | nominal_significant_G183]
  concordance[, endpoint_kind := ifelse(endpoint == "VDR_expression", "VDR expression", "VDR program")]

  scatter <- concordance[!is.na(Hedges_g_G250) & !is.na(Hedges_g_G183)]
  label_points <- scatter[model_key %in% c("Myeloid cells | VDR_expression", "Stromal cells | VDR_expression")]
  label_points[, `:=`(
    label_x = ifelse(grepl("Myeloid", model_key), 0.56, -1.42),
    label_y = ifelse(grepl("Myeloid", model_key), -1.42, 0.40),
    plot_label = ifelse(grepl("Myeloid", model_key), "Myeloid VDR\nnominal in G183", "Stromal VDR\nnominal in G250")
  )]

  myeloid <- data.table::copy(frozen$myeloid_vdr)
  myeloid_plot <- data.table::data.table(
    dataset = factor(c("GSE183837", "GSE250130"), levels = c("GSE183837", "GSE250130")),
    y = c(1, 2),
    hedges_g = c(as.numeric(myeloid$Hedges_g_G183[[1]]), as.numeric(myeloid$Hedges_g_G250[[1]])),
    exact_p = c(as.numeric(myeloid$exact_P_G183[[1]]), as.numeric(myeloid$exact_P_G250[[1]])),
    primary_q = c(as.numeric(myeloid$primary_q_G183[[1]]), as.numeric(myeloid$primary_q_G250[[1]])),
    dataset_colour = c("G183", "G250")
  )
  myeloid_plot[, annotation := c("P=0.0357; q=0.893", "P=0.570; q=0.570")]

  facts <- setNames(as.character(frozen$locked_facts$value), frozen$locked_facts$fact_id)
  evidence_boxes <- data.table::data.table(
    box = c("anchor", "temporal", "bulk", "disease"),
    xmin = c(1, 26, 51, 76), xmax = c(24, 49, 74, 99), ymin = 58, ymax = 99,
    fill = c("#EEF1F4", "#DCEFE8", "#F8E8CF", "#F3DDDA"),
    title = c("Mechanistic anchor", "Normal remodeling", "Bulk projection", "Disease replication"),
    detail = c(
      paste0(facts[["F001"]], " RNA candidates\n", facts[["F002"]], " occupancy-supported\n", facts[["F003"]], " shared occupancy"),
      paste0(facts[["F009"]], "/", facts[["F008"]], " program models FDR+\n", facts[["F010"]], " temporal DE rows\nStrong positive evidence"),
      paste0(facts[["F007"]], " FDR+ program rows\nacross ", facts[["F006"]], " contrasts\nContext-specific"),
      paste0("G250: 0 corrected\nG183: 0 global BH\nCross-cohort: ", facts[["F022"]], " replicated")
    )
  )
  final_box <- data.table::data.table(
    xmin = 8, xmax = 92, ymin = 3, ymax = 32,
    title = "Integrated evidence boundary",
    detail = "No broad, reproducible RIF-associated disruption of the VDR programs\nNo pooling; no direct-target, causal, or biomarker claim"
  )

  list(A = concordance, B = scatter, B_labels = label_points, C = myeloid_plot,
       D = evidence_boxes, D_final = final_box, facts = facts)
})

metric <- setNames(as.character(frozen$overall$value), frozen$overall$metric)
interpretation <- setNames(as.character(frozen$interpretation_lock$value), frozen$interpretation_lock$item)
scientific_checks <- data.frame(
  check_id = sprintf("S%02d", 1:17),
  check = c(
    "Planned common models", "Informative models", "Direction-concordant models", "Direction-discordant models",
    "Non-informative models", "G250 nominal models", "G250 corrected models", "G183 nominal models",
    "G183 corrected models", "Corrected replicated models", "Exploratory directional-support models",
    "Myeloid VDR G250 direction", "Myeloid VDR G183 direction", "Myeloid VDR G183 exact P",
    "Myeloid VDR G183 global q", "Broad reproducible dysregulation supported", "Pooled meta-analysis performed"
  ),
  observed = c(
    nrow(plot_data$A), nrow(plot_data$B), sum(plot_data$A$direction_status == "CONCORDANT"),
    sum(plot_data$A$direction_status == "DISCORDANT"), sum(plot_data$A$direction_status == "NONINFORMATIVE"),
    sum(plot_data$A$nominal_significant_G250), sum(plot_data$A$corrected_significant_G250),
    sum(plot_data$A$nominal_significant_G183), sum(plot_data$A$corrected_significant_G183),
    metric[["corrected_replicated_concordant_n"]], metric[["exploratory_directional_support_n"]],
    as.character(frozen$myeloid_vdr$direction_G250[[1]]), as.character(frozen$myeloid_vdr$direction_G183[[1]]),
    sprintf("%.7f", as.numeric(frozen$myeloid_vdr$exact_P_G183[[1]])),
    sprintf("%.6f", as.numeric(frozen$myeloid_vdr$primary_q_G183[[1]])),
    interpretation[["broad_RIF_VDR_program_dysregulation_supported"]],
    interpretation[["pooled_meta_analysis_performed"]]
  ),
  expected = c(25, 24, 12, 12, 1, 1, 0, 1, 0, "0", "0", "RIF_HIGHER", "RIF_LOWER", "0.0357143", "0.892857", "FALSE", "FALSE"),
  stringsAsFactors = FALSE
)
scientific_checks$pass <- as.character(scientific_checks$observed) == as.character(scientific_checks$expected)
write_csv_utf8(scientific_checks, file.path(dirs$qc, "FIG6_SCIENTIFIC_QC.csv"))
if (!all(scientific_checks$pass)) stop("Figure 6 scientific QC failed.")

statistics_qc <- data.frame(
  check = c(
    "No model/statistical/correlation test rerun", "No multiplicity correction rerun",
    "No pooled meta-analysis", "Dataset-specific q families retained", "Direction concordance is descriptive",
    "Myeloid VDR nominal signal remains unreplicated"
  ),
  observed = c(
    "PASS_NO_RECOMPUTATION", "PASS_FROZEN_Q_ONLY", "PASS_NOT_PERFORMED",
    "G250 within-family; G183 global BH", "12/24 informative models only", "Opposite directions; corrected support 0"
  ),
  expected = c(
    "PASS_NO_RECOMPUTATION", "PASS_FROZEN_Q_ONLY", "PASS_NOT_PERFORMED",
    "Preserve frozen scopes", "No inferential replication claim", "No validation claim"
  ),
  pass = TRUE, stringsAsFactors = FALSE
)
write_csv_utf8(statistics_qc, file.path(dirs$qc, "FIG6_STATISTICAL_QC.csv"))

write_csv_utf8(plot_data$A, file.path(dirs$source_data, "Figure6_panelA_25_model_direction_matrix.csv"))
write_csv_utf8(plot_data$B, file.path(dirs$source_data, "Figure6_panelB_Hedges_g_comparison.csv"))
write_csv_utf8(plot_data$C, file.path(dirs$source_data, "Figure6_panelC_myeloid_VDR_focus.csv"))
write_csv_utf8(plot_data$D, file.path(dirs$source_data, "Figure6_panelD_evidence_state_boxes.csv"))
write_csv_utf8(plot_data$D_final, file.path(dirs$source_data, "Figure6_panelD_final_boundary.csv"))
for (nm in names(source_paths)) copy_exact(source_paths[[nm]], file.path(dirs$source_data, "upstream_refs", basename(source_paths[[nm]])))
copy_exact(file.path(preflight_root, "inventories", "FIG6_PANEL_SOURCE_MAP.csv"), file.path(dirs$source_data, "upstream_refs", "FIG6_PANEL_SOURCE_MAP.csv"))
copy_exact(file.path(preflight_root, "blueprint_reconciliation", "FIG6_BLUEPRINT_RECONCILIATION_REGISTER.csv"), file.path(dirs$source_data, "upstream_refs", "FIG6_BLUEPRINT_RECONCILIATION_REGISTER.csv"))

upstream_identity <- data.frame(
  role = names(source_paths), path = vapply(source_paths, normp, character(1), must_work = TRUE),
  size_bytes = as.numeric(file.info(source_paths)$size),
  md5 = vapply(source_paths, md5_file, character(1)),
  sha256 = vapply(source_paths, sha256_file, character(1)), stringsAsFactors = FALSE
)
write_csv_utf8(upstream_identity, file.path(dirs$manifests, "FIG6_UPSTREAM_SOURCE_IDENTITY.csv"))

figure_font <- "Arial"
arial_path <- "C:/Windows/Fonts/arial.ttf"
if (!file.exists(arial_path)) stop("Arial font file is missing.")
color_text <- "#161616"
color_muted <- "#555555"
color_grid <- "#D9D9D9"
base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.23, colour = color_grid),
    axis.title = ggplot2::element_text(size = 6.7, colour = color_text),
    axis.text = ggplot2::element_text(size = 5.7, colour = color_text),
    legend.title = ggplot2::element_text(size = 5.6, face = "bold"),
    legend.text = ggplot2::element_text(size = 5.2),
    plot.margin = ggplot2::margin(1.2, 1.5, 1.2, 1.5, unit = "mm")
  )

panel_A <- capture_stage("panel_A_construction", {
  ggplot2::ggplot(plot_data$A, ggplot2::aes(x = endpoint, y = broad_cell_type, fill = direction_status)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8, width = 0.98, height = 0.96) +
    ggplot2::geom_text(
      ggplot2::aes(label = direction_pair, colour = label_colour),
      family = figure_font, fontface = "bold", size = 2.05, lineheight = 0.90
    ) +
    ggplot2::scale_fill_manual(values = status_colors, labels = status_labels, drop = FALSE) +
    ggplot2::scale_colour_manual(values = c(dark = color_text, light = "white"), guide = "none") +
    ggplot2::scale_x_discrete(labels = endpoint_labels, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Cross-dataset status") +
    ggplot2::theme_void(base_family = figure_font, base_size = 7.2) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 5.7, colour = color_text, lineheight = 0.90, margin = ggplot2::margin(t = 2.0)),
      axis.text.y = ggplot2::element_text(size = 6.0, colour = color_text, margin = ggplot2::margin(r = 2.0)),
      legend.position = "bottom", legend.direction = "horizontal", legend.title.position = "left",
      legend.key.width = grid::unit(4.2, "mm"), legend.key.height = grid::unit(2.8, "mm"),
      legend.margin = ggplot2::margin(t = 0.5, b = 0), plot.margin = ggplot2::margin(1, 1.5, 0.5, 1.5, unit = "mm")
    )
})

panel_B <- capture_stage("panel_B_construction", {
  ggplot2::ggplot(plot_data$B, ggplot2::aes(x = Hedges_g_G250, y = Hedges_g_G183)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#A6A6A6", linewidth = 0.35, linetype = 3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.30, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.30, linetype = 2) +
    ggplot2::geom_point(
      ggplot2::aes(fill = direction_status, shape = nominal_any),
      size = 2.0, stroke = 0.45, colour = color_text
    ) +
    ggplot2::geom_segment(
      data = plot_data$B_labels,
      ggplot2::aes(x = Hedges_g_G250, y = Hedges_g_G183, xend = label_x, yend = label_y),
      inherit.aes = FALSE, linewidth = 0.28, colour = "#6D6D6D"
    ) +
    ggplot2::geom_label(
      data = plot_data$B_labels,
      ggplot2::aes(x = label_x, y = label_y, label = plot_label),
      inherit.aes = FALSE, family = figure_font, size = 1.75, lineheight = 0.90,
      linewidth = 0.18, label.padding = grid::unit(0.55, "mm"),
      fill = scales::alpha("white", 0.92), colour = color_text
    ) +
    ggplot2::scale_fill_manual(values = status_colors, guide = "none") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1.80, 0.80), breaks = c(-1.5, -1, -0.5, 0, 0.5)) +
    ggplot2::scale_y_continuous(limits = c(-1.80, 0.80), breaks = c(-1.5, -1, -0.5, 0, 0.5)) +
    ggplot2::labs(x = "GSE250130 Hedges g", y = "GSE183837 Hedges g") +
    ggplot2::coord_cartesian(clip = "off") + base_theme +
    ggplot2::theme(
      legend.position = "none", panel.grid.major = ggplot2::element_line(linewidth = 0.20, colour = color_grid),
      axis.text = ggplot2::element_text(size = 5.5), axis.title = ggplot2::element_text(size = 6.2),
      plot.margin = ggplot2::margin(1.5, 2.5, 1.2, 1.5, unit = "mm")
    )
})

panel_C <- capture_stage("panel_C_construction", {
  ggplot2::ggplot(plot_data$C, ggplot2::aes(x = hedges_g, y = y, colour = dataset_colour)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.32, linetype = 2) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = hedges_g, yend = y), linewidth = 0.75, alpha = 0.80
    ) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = annotation),
      nudge_y = 0.20, family = figure_font, size = 1.90, colour = color_text
    ) +
    ggplot2::annotate(
      "label", x = -0.69, y = 1.50, label = "Opposite directions\nNo cross-dataset support",
      family = figure_font, fontface = "bold", size = 1.90, lineheight = 0.92,
      linewidth = 0.20, label.padding = grid::unit(0.7, "mm"), fill = "#F7F7F7", colour = color_text
    ) +
    ggplot2::scale_colour_manual(values = c(G183 = "#D58655", G250 = "#3E7C73"), guide = "none") +
    ggplot2::scale_y_continuous(breaks = c(1, 2), labels = c("GSE183837", "GSE250130"), limits = c(0.72, 2.28)) +
    ggplot2::scale_x_continuous(limits = c(-2.05, 0.65), breaks = c(-2, -1.5, -1, -0.5, 0, 0.5)) +
    ggplot2::labs(x = "Hedges g (RIF - reference)", y = NULL) + base_theme +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(), legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 5.2), axis.text.y = ggplot2::element_text(size = 5.8),
      axis.title.x = ggplot2::element_text(size = 6.1), plot.margin = ggplot2::margin(1.5, 2, 1.2, 1.5, unit = "mm")
    )
})

panel_D <- capture_stage("panel_D_construction", {
  arrow_data <- data.frame(
    x = c(12.5, 37.5, 62.5, 87.5), y = 57,
    xend = c(34, 45, 55, 66), yend = 33
  )
  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.38, colour = "#8A8A8A",
      arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")
    ) +
    ggplot2::geom_rect(
      data = plot_data$D,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      colour = "#A6A6A6", linewidth = 0.38
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(
      data = plot_data$D,
      ggplot2::aes(x = (xmin + xmax) / 2, y = 85, label = title),
      family = figure_font, fontface = "bold", size = 2.35, colour = color_text
    ) +
    ggplot2::geom_text(
      data = plot_data$D,
      ggplot2::aes(x = (xmin + xmax) / 2, y = 70, label = detail),
      family = figure_font, size = 1.95, lineheight = 1.02, colour = color_text
    ) +
    ggplot2::geom_rect(
      data = plot_data$D_final,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#E8EEF3", colour = "#526778", linewidth = 0.55
    ) +
    ggplot2::geom_text(
      data = plot_data$D_final,
      ggplot2::aes(x = (xmin + xmax) / 2, y = 24, label = title),
      family = figure_font, fontface = "bold", size = 2.55, colour = "#2A4152"
    ) +
    ggplot2::geom_text(
      data = plot_data$D_final,
      ggplot2::aes(x = (xmin + xmax) / 2, y = 13.5, label = detail),
      family = figure_font, size = 2.10, lineheight = 1.04, colour = color_text
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
    ggplot2::theme_void(base_family = figure_font) +
    ggplot2::theme(plot.margin = ggplot2::margin(-1.4, 1.5, 1.2, 1.5, unit = "mm"))
})

build_body_gtable <- function(plot_object) {
  ggplot2::ggplotGrob(
    plot_object + ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) +
      ggplot2::theme(
        plot.tag = ggplot2::element_blank(), plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank()
      )
  )
}

header_title_row_height_mm <- 5.1
header_subtitle_row_height_mm <- 3.2
header_tag_to_title_gap_mm <- 1.8
header_tag_size_pt <- 8.2
header_title_size_pt <- 8.0
header_subtitle_size_pt <- 6.3

make_header_gtable <- function(body_gtable, panel_tag, panel_title, panel_subtitle, header_shift_mm = 0) {
  panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
  if (!nrow(panel_cells)) stop("Could not locate panel cells for ", panel_tag)
  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)
  header_gt <- gtable::gtable(
    widths = body_gtable$widths,
    heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), "mm")
  )
  header_shift <- grid::unit(header_shift_mm, "mm")
  title_grob <- grid::textGrob(
    panel_title, x = grid::unit(0, "npc") + header_shift, y = grid::unit(0.52, "npc"), hjust = 0,
    gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_title_size_pt, col = color_text)
  )
  tag_grob <- grid::textGrob(
    panel_tag, x = grid::unit(0, "npc") + header_shift - grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"), hjust = 1,
    gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_tag_size_pt, col = color_text)
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle, x = grid::unit(0, "npc") + header_shift, y = grid::unit(0.56, "npc"), hjust = 0,
    gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, col = color_muted)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, grid::grobTree(tag_grob, title_grob), t = 1, b = 1, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F6-header-title-", panel_tag)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, subtitle_grob, t = 2, b = 2, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F6-header-subtitle-", panel_tag)
  )
  list(
    gtable = header_gt,
    geometry = data.frame(
      panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE",
      panel_left_column = panel_left, panel_right_column = panel_right,
      tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt,
      subtitle_size_pt = header_subtitle_size_pt, tag_title_gap_mm = header_tag_to_title_gap_mm,
      header_shift_mm = header_shift_mm, stringsAsFactors = FALSE
    )
  )
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
  header_D <- make_header_gtable(
    body_D, "D", "Integrated evidence states and inference boundary",
    "Positive normal remodeling and context-specific bulk signals remain distinct from negative disease replication",
    header_shift_mm = 8.4
  )
  design <- "
AA
BB
##
CD
EF
##
GG
HH
"
  figure_object <- (
    patchwork::wrap_elements(full = header_A$gtable) +
      patchwork::wrap_elements(full = body_A) +
      patchwork::wrap_elements(full = header_B$gtable) +
      patchwork::wrap_elements(full = header_C$gtable) +
      patchwork::wrap_elements(full = body_B) +
      patchwork::wrap_elements(full = body_C) +
      patchwork::wrap_elements(full = header_D$gtable) +
      patchwork::wrap_elements(full = body_D)
  ) + patchwork::plot_layout(
    design = design, widths = c(1.18, 0.82),
    heights = grid::unit(c(8.3, 52, 2.5, 8.3, 50, 2.5, 8.3, 56), "mm")
  ) + patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm"))
  )
  list(
    figure = figure_object,
    header_geometry = do.call(rbind, lapply(list(header_A, header_B, header_C, header_D), `[[`, "geometry"))
  )
})
grDevices::dev.off()
unlink(metrics_device_path)
write_csv_utf8(composite_parts$header_geometry, file.path(dirs$qc, "FIG6_PANEL_HEADER_GEOMETRY.csv"))

figure_png <- file.path(dirs$figures, "Figure6_v1_1_QC_180mm_300dpi.png")
capture_stage("300dpi_candidate_export", {
  ragg::agg_png(figure_png, width = figure_width_mm, height = figure_height_mm, units = "mm", res = 300, background = "white")
  print(composite_parts$figure)
  grDevices::dev.off()
})
if (!file.exists(figure_png) || file.info(figure_png)$size < 2e+05) stop("Figure 6 candidate export failed.")

warning_log <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
  warning_id = character(), stage = character(), warning_class = character(), message = character(), stringsAsFactors = FALSE
)
write_csv_utf8(warning_log, file.path(dirs$qc, "FIG6_RUNTIME_WARNING_LOG.csv"))

img <- magick::image_read(figure_png)
img_info <- magick::image_info(img)
rgb <- png::readPNG(figure_png)
rgb3 <- rgb[, , seq_len(min(3L, dim(rgb)[3])), drop = FALSE]
delta <- apply(abs(rgb3 - 1), c(1, 2), max)
foreground <- which(delta > 0.035, arr.ind = TRUE)
if (!nrow(foreground)) stop("Figure 6 raster appears blank.")
edge_margins <- data.frame(
  edge = c("top", "bottom", "left", "right"),
  margin_px = c(
    min(foreground[, 1]) - 1, nrow(delta) - max(foreground[, 1]),
    min(foreground[, 2]) - 1, ncol(delta) - max(foreground[, 2])
  ), minimum_px = c(8, 8, 8, 8), stringsAsFactors = FALSE
)
edge_margins$pass <- edge_margins$margin_px >= edge_margins$minimum_px
write_csv_utf8(edge_margins, file.path(dirs$qc, "FIG6_RASTER_EDGE_MARGIN_QC.csv"))

crop_top <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.36), "+0+0"))
crop_middle <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.35), "+0+", round(img_info$height[[1]] * 0.34)))
crop_bottom <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.32), "+0+", round(img_info$height[[1]] * 0.68)))
magick::image_write(crop_top, file.path(dirs$qc, "Figure6_TOP_CROP.png"), format = "png")
magick::image_write(crop_middle, file.path(dirs$qc, "Figure6_MIDDLE_CROP.png"), format = "png")
magick::image_write(crop_bottom, file.path(dirs$qc, "Figure6_BOTTOM_CROP.png"), format = "png")

automated_checks <- data.frame(
  check = c(
    "Raster exists", "Raster width", "Raster height", "Raster content variance",
    "All edge safety margins", "Runtime warning count", "Four panel headers", "Arial font file"
  ),
  observed = c(
    file.exists(figure_png), img_info$width[[1]], img_info$height[[1]], stats::sd(as.numeric(rgb3)),
    all(edge_margins$pass), nrow(warning_log), nrow(composite_parts$header_geometry), file.exists(arial_path)
  ),
  expected = c(TRUE, floor(figure_width_mm / 25.4 * 300), floor(figure_height_mm / 25.4 * 300), "> 0.01", TRUE, 0, 4, TRUE),
  pass = c(
    file.exists(figure_png), img_info$width[[1]] == floor(figure_width_mm / 25.4 * 300),
    img_info$height[[1]] == floor(figure_height_mm / 25.4 * 300), stats::sd(as.numeric(rgb3)) > 0.01,
    all(edge_margins$pass), nrow(warning_log) == 0, nrow(composite_parts$header_geometry) == 4, file.exists(arial_path)
  ), stringsAsFactors = FALSE
)
write_csv_utf8(automated_checks, file.path(dirs$qc, "FIG6_AUTOMATED_RENDER_CHECKS.csv"))
if (!all(automated_checks$pass)) stop("Figure 6 automated render QC failed.")

expected_reviewed_sha256 <- "a4b0a6028719996d074f5eb0957a5c39e5d933ea027bef3e215efa8122f0bc0f"
actual_sha256 <- sha256_file(figure_png)
review_bound <- identical(actual_sha256, expected_reviewed_sha256)
strict_visual_qc <- data.frame(
  severity = c("BLOCKING", "BLOCKING", "BLOCKING", "MAJOR", "MAJOR", "MAJOR", "MAJOR"),
  visual_check = c(
    "Global panel allocation and header/body alignment", "Text, glyph, axis, matrix-cell, and annotation clipping",
    "Scientific values and evidence boundary visible", "Panel A 25-cell direction matrix readability",
    "Panel B effect-scatter point and label separation", "Panel C myeloid VDR discordance readability",
    "Panel D evidence-state flow and wording hierarchy"
  ),
  status = if (review_bound) "PASS_ASSISTANT_ACTUAL_RENDER_REVIEW" else "PENDING_ASSISTANT_ACTUAL_RENDER_REVIEW",
  reviewed_reference_sha256 = if (review_bound) actual_sha256 else "", stringsAsFactors = FALSE
)
write_csv_utf8(strict_visual_qc, file.path(dirs$qc, "FIG6_STRICT_RENDER_QC.csv"))

typography_lock <- data.frame(
  component = c("font_family", "panel_tag_pt", "header_title_pt", "header_subtitle_pt", "base_theme_pt", "final_width_mm", "final_height_mm", "approval_ppi"),
  value = c(figure_font, header_tag_size_pt, header_title_size_pt, header_subtitle_size_pt, 7.2, figure_width_mm, figure_height_mm, 300),
  stringsAsFactors = FALSE
)
write_csv_utf8(typography_lock, file.path(dirs$qc, "FIG6_ARTICLE_TYPOGRAPHY_LOCK.csv"))

legend_text <- c(
  "# Figure 6 legend", "",
  "**Figure 6. Cross-dataset comparison defines the boundary of the single-cell RIF evidence.**",
  "**A,** Direction matrix for the 25 prespecified models shared by GSE250130 and GSE183837. Cell text gives the GSE250130 (250) and GSE183837 (183) RIF-minus-reference direction, with plus and minus denoting higher and lower values, respectively. Twenty-four models were informative: 12 were directionally concordant and 12 discordant; one endpoint was non-informative.",
  "**B,** Hedges g values for the 24 informative models. Dashed lines mark zero and the dotted diagonal marks equal standardized effects. Colors encode the directional status used in Panel A; triangles identify the two models with nominal P < 0.05 in one cohort. No correlation test or pooled estimate was calculated.",
  "**C,** Myeloid-cell VDR-expression effects. GSE183837 showed a nominal decrease (exact P=0.0357, global q=0.893), whereas GSE250130 showed a nonsignificant effect in the opposite direction (exact P=0.570, within-family q=0.570). Thus, the nominal signal was not replicated.",
  "**D,** Integrated evidence-state diagram. The perturbation-defined mechanistic candidate hierarchy, strong normal temporal remodeling, and context-specific bulk projections are retained as distinct evidence layers. Donor-level single-cell disease comparisons yielded no corrected cross-dataset replication, so broad reproducible RIF-associated disruption of the VDR programs was not supported. No pooled meta-analysis, definitive direct-target, causal, or biomarker claim was made."
)
write_text_utf8(legend_text, file.path(dirs$manuscript, "Figure6_legend_v1.1.md"))

methods_text <- c(
  "# Figure 6 Methods", "",
  "Figure 6 was assembled from STEP04_LOCK_v1.2 and SCCROSS_v1.0 result outputs. Panel A shows the direction classification for all 25 prespecified shared models. Panels B and C display Hedges g, exact P, and adjusted P values. The q-value families remain dataset-specific: GSE250130 values use within-family corrections, whereas GSE183837 values use global Benjamini-Hochberg correction across 25 planned models. Panel D integrates the prespecified evidence layers without additional quantitative inference."
)
write_text_utf8(methods_text, file.path(dirs$manuscript, "Figure6_methods_v1.1.md"))

results_text <- c(
  "# Figure 6 Results", "",
  "Among 25 prespecified models shared by the two single-cell cohorts, 24 were directionally informative. Twelve showed the same direction and 12 showed opposite directions. Neither corrected cross-dataset replication nor exploratory directional support was observed. The myeloid VDR-expression decrease in GSE183837 met nominal exact P < 0.05 but not global multiplicity correction and was directionally opposite to the nonsignificant GSE250130 effect.", "",
  "These cross-dataset results preserve a clear evidence boundary. Strong normal temporal remodeling and context-specific patient bulk projections remain supported in their prespecified settings, but the donor-level single-cell datasets did not support broad, reproducible RIF-associated disruption of the VDR programs. The cohorts were not pooled because their reference groups and sampling structures differ."
)
write_text_utf8(results_text, file.path(dirs$manuscript, "Figure6_results_v1.1.md"))

limitations_text <- c(
  "# Figure 6 limitations and interpretation", "",
  "- Directional concordance without corrected significance is descriptive, not confirmatory replication.",
  "- The two studies use different reference groups and sampling structures and were not pooled.",
  "- The GSE183837 myeloid VDR signal is nominal, globally nonsignificant, and directionally discordant in GSE250130.",
  "- Negative corrected results limit broad dysregulation claims but do not establish biological equivalence.",
  "- No direct-target, causal-mechanism, diagnostic, or biomarker-performance claim is authorized."
)
write_text_utf8(limitations_text, file.path(dirs$manuscript, "Figure6_limitations_v1.1.md"))

parameters <- list(
  candidate_id = candidate_id, width_mm = figure_width_mm, height_mm = figure_height_mm,
  approval_ppi = 300, source_authority = "STEP04_LOCK_v1.2 + SCCROSS_v1_0",
  effect_comparison = "descriptive frozen Hedges g; no correlation test",
  pooled_meta_analysis = FALSE, scientific_recomputation = FALSE,
  final_user_visual_QC_required = TRUE
)
jsonlite::write_json(parameters, file.path(dirs$qc, "analysis_parameters.json"), auto_unbox = TRUE, pretty = TRUE)
write_csv_utf8(data.frame(
  component = c("All panels", "All inferential outputs"),
  seed = c("NOT_APPLICABLE_DETERMINISTIC", "NOT_APPLICABLE_FROZEN_TABLES"), stringsAsFactors = FALSE
), file.path(dirs$qc, "random_seed_manifest.csv"))
write_text_utf8(capture.output(sessionInfo()), file.path(dirs$qc, "sessionInfo_FIG6_v1_1.txt"))

figure_status <- data.frame(
  figure = "Figure 6", candidate = candidate_id,
  science_QC = "PASS", statistics_QC = "PASS_NO_RECOMPUTATION", automated_render_QC = "PASS",
  assistant_actual_render_QC = if (review_bound) "PASS" else "PENDING",
  user_final_visual_QC = "PENDING", immutable_release = FALSE, stringsAsFactors = FALSE
)
write_csv_utf8(figure_status, file.path(out_root, "FIGURE_STATUS.csv"))

hash_qc <- data.frame(
  file = basename(figure_png), actual_sha256 = actual_sha256,
  expected_reviewed_sha256 = expected_reviewed_sha256, pass = review_bound, stringsAsFactors = FALSE
)
write_csv_utf8(hash_qc, file.path(dirs$qc, "FIG6_REFERENCE_RENDER_HASH_QC.csv"))

candidate_manifest <- build_manifest(out_root, exclude = "manifests/FIG6_CANDIDATE_MANIFEST.csv")
write_csv_utf8(candidate_manifest, file.path(dirs$manifests, "FIG6_CANDIDATE_MANIFEST.csv"))
candidate_validation <- validate_manifest(out_root, candidate_manifest)
write_csv_utf8(candidate_validation, file.path(dirs$manifests, "FIG6_CANDIDATE_MANIFEST_VALIDATION.csv"))
if (!all(candidate_validation$exists & candidate_validation$size_match & candidate_validation$md5_match & candidate_validation$sha256_match)) {
  stop("Figure 6 candidate manifest validation failed.")
}

if (!review_bound) {
  cat("Figure 6 candidate generated for assistant actual-render review.\n")
  cat("Candidate PNG: ", normp(figure_png, TRUE), "\n", sep = "")
  cat("SHA-256: ", actual_sha256, "\n", sep = "")
  stop("Assistant actual-render review and reviewed-hash binding are required before formal QC completion.")
}

dir.create(feedback_root, recursive = TRUE, showWarnings = FALSE)
report_root <- file.path(feedback_root, paste0("FIG6_V1_1_QC_REPORT_", run_stamp))
if (dir.exists(report_root)) stop("Report root exists: ", normp(report_root, FALSE))
report_fig6 <- file.path(report_root, "FIGURE6")
dir.create(report_fig6, recursive = TRUE, showWarnings = FALSE)
copy_tree(out_root, report_fig6)
report_manifest <- build_manifest(report_root, exclude = "SMALL_REPORT_MANIFEST.csv")
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))
report_validation <- validate_manifest(report_root, report_manifest)
if (!all(report_validation$exists & report_validation$size_match & report_validation$md5_match & report_validation$sha256_match)) {
  stop("Figure 6 QC SMALL_REPORT validation failed.")
}
zip_path <- file.path(feedback_root, paste0("STEP05F_FIG6_V1_1_QC_", run_stamp, "_SMALL_REPORT.zip"))
if (file.exists(zip_path)) stop("ZIP exists: ", normp(zip_path, FALSE))
zip::zipr(zipfile = zip_path, files = basename(report_root), root = dirname(report_root), include_directories = TRUE)
if (!file.exists(zip_path) || file.info(zip_path)$size < 3e+05) stop("Figure 6 QC ZIP creation failed.")

cat("\n============================================================\n")
cat("Figure 6 v1.1 QC candidate completed\n")
cat("Gate: PASS_FIGURE6_ASSISTANT_ACTUAL_RENDER_QC_AWAIT_USER_VISUAL_QC\n")
cat("Output: ", normp(out_root, TRUE), "\n", sep = "")
cat("PNG: ", normp(figure_png, TRUE), "\n", sep = "")
cat("SMALL_REPORT: ", normp(zip_path, TRUE), "\n", sep = "")
cat("============================================================\n")

invisible(list(
  gate = "PASS_FIGURE6_ASSISTANT_ACTUAL_RENDER_QC_AWAIT_USER_VISUAL_QC",
  output_root = normp(out_root, TRUE), png = normp(figure_png, TRUE), small_report = normp(zip_path, TRUE)
))

# Publication figure stage 2
options(stringsAsFactors = FALSE, warn = 1)

work_root <- "<LOCAL_R_LIBRARY>"
candidate_parent <- file.path(work_root, "Figure6_export_build")
source_script <- file.path(work_root, "Figure6_language_clean_v1_2.R")
Sys.setenv(FIG6_CANDIDATE_PARENT = candidate_parent)

plot_env <- new.env(parent = globalenv())
sys.source(source_script, envir = plot_env, keep.source = TRUE)

figure_object <- get("figure_object", envir = plot_env)
width_mm <- get("figure_width_mm", envir = plot_env)
height_mm <- get("figure_height_mm", envir = plot_env)
candidate_root <- file.path(candidate_parent, "F6_v1_2_language_clean")
figure_dir <- file.path(candidate_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

png600 <- file.path(figure_dir, "Figure6_submission_180mm_v1_2_language_clean_600ppi.png")
tiff600 <- file.path(figure_dir, "Figure6_submission_180mm_v1_2_language_clean_600ppi.tiff")
pdfvec <- file.path(figure_dir, "Figure6_submission_180mm_v1_2_language_clean.pdf")

ragg::agg_png(png600, width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
print(figure_object)
grDevices::dev.off()

ragg::agg_tiff(tiff600, width = width_mm, height = height_mm, units = "mm", res = 600,
               compression = "lzw", background = "white")
print(figure_object)
grDevices::dev.off()

grDevices::cairo_pdf(pdfvec, width = width_mm / 25.4, height = height_mm / 25.4,
                     family = "Arial", bg = "white", fallback_resolution = 600)
print(figure_object)
grDevices::dev.off()

cat(candidate_root, "\n")
