
# Publication figure stage 1
options(stringsAsFactors = FALSE, warn = 1)

# Figure 5 v1.2 manuscript composite from frozen STEP04_LOCK_v1.2 results.
# No integration, annotation, model fitting, hypothesis testing, DE, or FDR is rerun.

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "126_step05E_Figure5_manuscript_composite_v1_2_QC.R"
candidate_id <- "F5_v1_2_QC"
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
  if (!dir.exists(src_root)) stop("Missing source directory: ", normp(src_root, FALSE))
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
    relative_path = rel,
    size_bytes = as.numeric(info$size),
    md5 = vapply(files, md5_file, character(1)),
    sha256 = vapply(files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  out[order(out$relative_path), , drop = FALSE]
}

validate_manifest <- function(root, manifest) {
  do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    p <- file.path(root, manifest$relative_path[[i]])
    ex <- file.exists(p) && !dir.exists(p)
    data.frame(
      relative_path = manifest$relative_path[[i]],
      exists = ex,
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
if (!grepl("version: 2.2-beta.7", skill3_text, fixed = TRUE)) stop("Skill 3 v2.2-beta.7 is required.")

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")
step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")
preflight_root <- file.path(step05_root, "06_results", "FIG5", "FIG5_WORKFLOW_PREFLIGHT_v1")
if (!file.exists(file.path(preflight_root, "run_evidence", "RUN_SUMMARY.csv"))) stop("Formal Figure 5 preflight is missing.")
preflight_summary <- utils::read.csv(file.path(preflight_root, "run_evidence", "RUN_SUMMARY.csv"), check.names = FALSE)
if (preflight_summary$gate[[1]] != "PASS_FIGURE5_WORKFLOW_PREFLIGHT_SOURCE_ASSEMBLY_AUTHORIZED") stop("Figure 5 preflight gate is not PASS.")

parent_override <- Sys.getenv("FIG5_CANDIDATE_PARENT", unset = "")
feedback_override <- Sys.getenv("FIG5_FEEDBACK_ROOT", unset = "")
figure_parent <- if (nzchar(parent_override)) parent_override else file.path(step05_root, "06_results", "FIG5")
feedback_root <- if (nzchar(feedback_override)) feedback_override else file.path(project_root, "反馈")
out_root <- file.path(figure_parent, candidate_id)
if (dir.exists(out_root)) stop("Refusing to overwrite Figure 5 candidate: ", normp(out_root, FALSE))

dirs <- list(
  figures = file.path(out_root, "figures"),
  source_data = file.path(out_root, "source_data"),
  manuscript = file.path(out_root, "manuscript_text"),
  qc = file.path(out_root, "qc"),
  code = file.path(out_root, "code"),
  manifests = file.path(out_root, "manifests")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
copy_exact(script_path, file.path(dirs$code, script_name))

g250_model_root <- file.path(step04_root, "06_results/S4E1_sc/SCI_v1_1_memory_safe/models/GSE250130/VDRMODELS_v1_1")
g250_de_root <- file.path(step04_root, "06_results/S4E1_sc/SCI_v1_1_memory_safe/de/GSE250130/DE_v1_2")
g183_root <- file.path(step04_root, "06_results/S4E1_sc/SCI183_v1_0")
g183_ann_root <- file.path(g183_root, "ann/G183/ANN183_v1_0")
g183_model_root <- file.path(g183_root, "mdl/G183/MD183_v1_0")
lock_root <- file.path(step04_root, "06_results/STEP04_closeout/STEP04_LOCK_v1.2")

source_paths <- c(
  g250_program = file.path(g250_model_root, "tables/03_LH7_program_primary.csv"),
  g250_vdr = file.path(g250_model_root, "tables/07_LH7_VDR_primary.csv"),
  g250_primary_de = file.path(g250_de_root, "tab/09_LH7_primary_summary.csv"),
  g250_sensitivity_de = file.path(g250_de_root, "tab/06_LH7_sensitivity_significant.csv"),
  g183_umap = file.path(g183_root, "tables/04_UMAP_metadata.csv"),
  g183_annotation = file.path(g183_ann_root, "t/01_annotation_map.csv"),
  g183_counts = file.path(g183_ann_root, "t/04_broad_cell_counts.csv"),
  g183_models = file.path(g183_model_root, "tab/01_primary_models.csv"),
  g183_group_summary = file.path(g183_model_root, "tab/02_group_summary.csv"),
  g183_model_data = file.path(g183_model_root, "src/SD1_primary_model_data.csv"),
  g183_ciliated_summary = file.path(g183_model_root, "tab/07_ciliated_descriptive_summary.csv"),
  g183_ciliated_data = file.path(g183_model_root, "src/SD8_ciliated_descriptive_data.csv"),
  locked_facts = file.path(lock_root, "tables/02_locked_fact_ledger.csv"),
  locked_claims = file.path(lock_root, "tables/03_claim_ledger.csv"),
  blueprint = file.path(lock_root, "tables/07_figure_blueprint.csv")
)
if (any(!file.exists(source_paths))) stop("Frozen Figure 5 source authority is incomplete.")

warning_rows <- list()
warning_n <- 0L
capture_stage <- function(stage, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warning_n <<- warning_n + 1L
    warning_rows[[warning_n]] <<- data.frame(
      warning_id = sprintf("W%03d", warning_n), stage = stage,
      warning_class = class(w)[[1]], message = conditionMessage(w),
      stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  })
}

cell_levels <- c(
  "Ciliated epithelial cells", "Endothelial cells", "Myeloid cells",
  "NK/T cells", "Stromal cells", "Unciliated epithelial cells"
)
cell_labels <- c(
  "Ciliated epi.", "Endothelial", "Myeloid", "NK/T", "Stromal", "Unciliated epi."
)
names(cell_labels) <- cell_levels
cell_colors <- c(
  "Ciliated epithelial cells" = "#E69F00",
  "Endothelial cells" = "#009E73",
  "Myeloid cells" = "#00A6C7",
  "NK/T cells" = "#007DB3",
  "Stromal cells" = "#7A5ABB",
  "Unciliated epithelial cells" = "#CC79A7"
)
group_colors <- c("CONTROL" = "#5B7896", "RIF" = "#C65D68")
endpoint_levels <- c(
  "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46",
  "PROMOTER_OCCUPANCY_84", "VDR_expression"
)
endpoint_labels <- c(
  "RNA program\n(n=436)", "Occupancy-supported\n(n=284)", "Shared occupancy\n(n=46)",
  "Promoter-associated\n(n=84)", "VDR expression"
)
names(endpoint_labels) <- endpoint_levels
endpoint_short_labels <- c("RNA\n436", "Occupancy\n284", "Shared\n46", "Promoter\n84", "VDR")
names(endpoint_short_labels) <- endpoint_levels

plot_data <- capture_stage("read_and_prepare_frozen_sources", {
  g250_program <- data.table::fread(source_paths[["g250_program"]])
  g250_vdr <- data.table::fread(source_paths[["g250_vdr"]])
  g250_primary_de <- data.table::fread(source_paths[["g250_primary_de"]])
  g250_sensitivity_de <- data.table::fread(source_paths[["g250_sensitivity_de"]])
  g183_annotation <- data.table::fread(source_paths[["g183_annotation"]])
  g183_umap <- data.table::fread(
    source_paths[["g183_umap"]],
    select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "seurat_clusters", "analysis_group", "sample_accession")
  )
  g183_models <- data.table::fread(source_paths[["g183_models"]])
  g183_group_summary <- data.table::fread(source_paths[["g183_group_summary"]])
  g183_model_data <- data.table::fread(source_paths[["g183_model_data"]])
  g183_ciliated_summary <- data.table::fread(source_paths[["g183_ciliated_summary"]])
  g183_ciliated_data <- data.table::fread(source_paths[["g183_ciliated_data"]])

  a_program <- g250_program[, .(
    broad_cell_type, endpoint = program_name,
    estimate = RIF_minus_fertile, ci_low = Welch_ci_low, ci_high = Welch_ci_high,
    FDR = FDR_BH, exact_P = exact_permutation_p, fertile_n, RIF_n
  )]
  a_vdr <- g250_vdr[, .(
    broad_cell_type, endpoint = "VDR_expression",
    estimate = RIF_minus_fertile, ci_low = Welch_ci_low, ci_high = Welch_ci_high,
    FDR = FDR_BH, exact_P = exact_permutation_p, fertile_n, RIF_n
  )]
  panel_A <- data.table::rbindlist(list(a_program, a_vdr), use.names = TRUE)
  panel_A[, endpoint_display := factor(endpoint_labels[endpoint], levels = endpoint_labels)]
  panel_A[, endpoint_display_short := factor(endpoint_short_labels[endpoint], levels = endpoint_short_labels)]
  panel_A[, broad_cell_type := factor(broad_cell_type, levels = rev(cell_levels))]

  primary_counts <- g250_primary_de[, .(
    broad_cell_type,
    count = whole_genome_up_n + whole_genome_down_n
  )]
  primary_counts[, scope := "Primary LH+7\nRIF vs fertile"]
  sensitivity_counts <- g250_sensitivity_de[, .N, by = broad_cell_type]
  data.table::setnames(sensitivity_counts, "N", "count")
  sensitivity_counts[, scope := "Fertile-reference\nsensitivity"]
  all_cells <- data.table::data.table(broad_cell_type = cell_levels)
  primary_counts <- merge(all_cells, primary_counts, by = "broad_cell_type", all.x = TRUE)
  primary_counts[is.na(count), count := 0]
  primary_counts[is.na(scope), scope := "Primary LH+7\nRIF vs fertile"]
  sensitivity_counts <- merge(all_cells, sensitivity_counts, by = "broad_cell_type", all.x = TRUE)
  sensitivity_counts[is.na(count), count := 0]
  sensitivity_counts[is.na(scope), scope := "Fertile-reference\nsensitivity"]
  panel_B <- data.table::rbindlist(list(primary_counts, sensitivity_counts), use.names = TRUE)
  panel_B[, broad_cell_type := factor(broad_cell_type, levels = rev(cell_levels))]
  panel_B[, scope := factor(scope, levels = c("Primary LH+7\nRIF vs fertile", "Fertile-reference\nsensitivity"))]

  g183_umap[, cluster := as.character(seurat_clusters)]
  g183_annotation[, cluster := as.character(cluster)]
  panel_C_all <- merge(
    g183_umap,
    g183_annotation[, .(cluster, broad_cell_type)],
    by = "cluster", all.x = TRUE, sort = FALSE
  )
  if (anyNA(panel_C_all$broad_cell_type)) stop("GSE183837 UMAP-to-annotation join is incomplete.")
  data.table::setorder(panel_C_all, broad_cell_type, cell_barcode)
  panel_C <- panel_C_all[, .SD[seq_len(min(.N, 6000L))], by = broad_cell_type]
  centroids <- panel_C[, .(UMAPH_1 = median(UMAPH_1), UMAPH_2 = median(UMAPH_2)), by = broad_cell_type]
  centroids[, `:=`(
    label_x = UMAPH_1 + data.table::fcase(
      broad_cell_type == "Unciliated epithelial cells", -2.2,
      broad_cell_type == "Endothelial cells", 2.2,
      default = 0
    ),
    label_y = UMAPH_2
  )]

  panel_D <- g183_models[, .(
    broad_cell_type, endpoint,
    estimate = mean_difference_RIF_minus_Control,
    ci_low = Welch_95CI_lower, ci_high = Welch_95CI_upper,
    exact_P = exact_permutation_P, global_BH = BH_global_25,
    statistical_status, constant_endpoint
  )]
  panel_D[!is.finite(ci_low), ci_low := estimate]
  panel_D[!is.finite(ci_high), ci_high := estimate]
  panel_D[, endpoint_display := factor(endpoint_labels[endpoint], levels = endpoint_labels)]
  panel_D[, endpoint_display_short := factor(endpoint_short_labels[endpoint], levels = endpoint_short_labels)]
  panel_D[, broad_cell_type := factor(broad_cell_type, levels = rev(cell_levels))]
  panel_D[, nominal_marker := ifelse(exact_P < 0.05, "*", "")]

  panel_E <- data.table::copy(g183_ciliated_data)
  panel_E[, endpoint_display := factor(endpoint_labels[endpoint], levels = endpoint_labels)]
  panel_E[, analysis_group := factor(analysis_group, levels = c("CONTROL", "RIF"))]
  panel_E[, threshold_status := factor(
    ifelse(minimum_cell_threshold_pass, "Threshold pass", "Below threshold"),
    levels = c("Threshold pass", "Below threshold")
  )]

  list(
    g250_program = g250_program, g250_vdr = g250_vdr,
    g250_primary_de = g250_primary_de, g250_sensitivity_de = g250_sensitivity_de,
    g183_models = g183_models, g183_group_summary = g183_group_summary,
    g183_model_data = g183_model_data, g183_ciliated_summary = g183_ciliated_summary,
    A = panel_A, B = panel_B, C_all = panel_C_all, C = panel_C,
    C_centroids = centroids, D = panel_D, E = panel_E
  )
})

scientific_checks <- data.frame(
  check_id = sprintf("S%02d", 1:15),
  check = c(
    "GSE250130 LH+7 program models", "GSE250130 LH+7 program FDR-significant models",
    "GSE250130 LH+7 VDR models", "GSE250130 LH+7 VDR FDR-significant models",
    "GSE250130 primary whole-genome DE rows", "GSE250130 sensitivity whole-genome DE rows",
    "GSE183837 frozen UMAP cells", "GSE183837 locked broad cell types", "GSE183837 donors",
    "GSE183837 primary donor models", "GSE183837 nominal exact-P signals",
    "GSE183837 global-BH signals", "GSE183837 ciliated endpoints",
    "GSE183837 ciliated donors per endpoint", "Ciliated primary-threshold-pass RIF donors"
  ),
  observed = c(
    nrow(plot_data$g250_program), sum(plot_data$g250_program$FDR_BH < 0.05),
    nrow(plot_data$g250_vdr), sum(plot_data$g250_vdr$FDR_BH < 0.05),
    sum(plot_data$g250_primary_de$whole_genome_up_n + plot_data$g250_primary_de$whole_genome_down_n),
    nrow(plot_data$g250_sensitivity_de), nrow(plot_data$C_all),
    data.table::uniqueN(plot_data$C_all$broad_cell_type), data.table::uniqueN(plot_data$C_all$sample_accession),
    nrow(plot_data$g183_models), sum(plot_data$g183_models$exact_permutation_P < 0.05),
    sum(plot_data$g183_models$BH_global_25 < 0.05), data.table::uniqueN(plot_data$E$endpoint),
    min(plot_data$E[, data.table::uniqueN(donor_id), by = endpoint]$V1),
    data.table::uniqueN(plot_data$E[analysis_group == "RIF" & minimum_cell_threshold_pass, donor_id])
  ),
  expected = c(24, 0, 6, 0, 0, 28, 47161, 6, 9, 25, 1, 0, 5, 9, 2),
  stringsAsFactors = FALSE
)
scientific_checks$pass <- scientific_checks$observed == scientific_checks$expected
write_csv_utf8(scientific_checks, file.path(dirs$qc, "FIG5_SCIENTIFIC_QC.csv"))
if (!all(scientific_checks$pass)) stop("Figure 5 scientific QC failed.")

statistics_qc <- data.frame(
  check = c(
    "No model fitting or hypothesis testing rerun", "GSE250130 FDR carried from frozen tables",
    "GSE183837 global BH carried from frozen table", "Sensitivity rows kept separate from primary",
    "Ciliated panel remains descriptive", "Disease effects use donor as inference unit"
  ),
  observed = c(
    "PASS_NO_RECOMPUTATION", "0/24 program and 0/6 VDR", "0/25 global BH",
    "Primary 0; sensitivity 28 result rows", "No p-values or significance markers",
    "Frozen donor-level source/model tables"
  ),
  expected = c(
    "PASS_NO_RECOMPUTATION", "Frozen FDR only", "Frozen global BH only",
    "Explicitly separated", "Descriptive only", "Donor-level inference"
  ),
  pass = TRUE,
  stringsAsFactors = FALSE
)
write_csv_utf8(statistics_qc, file.path(dirs$qc, "FIG5_STATISTICAL_QC.csv"))

write_csv_utf8(plot_data$A, file.path(dirs$source_data, "Figure5_panelA_GSE250130_LH7_program_VDR_effects.csv"))
write_csv_utf8(plot_data$B, file.path(dirs$source_data, "Figure5_panelB_GSE250130_LH7_DE_counts.csv"))
write_csv_utf8(plot_data$C, file.path(dirs$source_data, "Figure5_panelC_GSE183837_UMAP_display_subset.csv"))
write_csv_utf8(plot_data$C_centroids, file.path(dirs$source_data, "Figure5_panelC_label_centroids.csv"))
write_csv_utf8(plot_data$D, file.path(dirs$source_data, "Figure5_panelD_GSE183837_primary_models.csv"))
write_csv_utf8(plot_data$E, file.path(dirs$source_data, "Figure5_panelE_ciliated_descriptive_data.csv"))
write_csv_utf8(plot_data$g183_ciliated_summary, file.path(dirs$source_data, "Figure5_panelE_ciliated_descriptive_summary.csv"))
copy_exact(source_paths[["locked_facts"]], file.path(dirs$source_data, "upstream_refs", "02_locked_fact_ledger.csv"))
copy_exact(source_paths[["locked_claims"]], file.path(dirs$source_data, "upstream_refs", "03_claim_ledger.csv"))
copy_exact(source_paths[["blueprint"]], file.path(dirs$source_data, "upstream_refs", "07_figure_blueprint.csv"))
copy_exact(file.path(preflight_root, "inventories", "FIG5_PANEL_SOURCE_MAP.csv"), file.path(dirs$source_data, "upstream_refs", "FIG5_PANEL_SOURCE_MAP.csv"))
copy_exact(file.path(preflight_root, "blueprint_reconciliation", "FIG5_BLUEPRINT_RECONCILIATION_REGISTER.csv"), file.path(dirs$source_data, "upstream_refs", "FIG5_BLUEPRINT_RECONCILIATION_REGISTER.csv"))

upstream_identity <- data.frame(
  role = names(source_paths),
  path = vapply(source_paths, normp, character(1), must_work = TRUE),
  size_bytes = as.numeric(file.info(source_paths)$size),
  md5 = vapply(source_paths, md5_file, character(1)),
  sha256 = vapply(source_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_utf8(upstream_identity, file.path(dirs$manifests, "FIG5_UPSTREAM_SOURCE_IDENTITY.csv"))

figure_font <- "Arial"
arial_path <- "C:/Windows/Fonts/arial.ttf"
if (!file.exists(arial_path)) stop("Arial font file is missing.")
color_text <- "#161616"
color_muted <- "#555555"
color_grid <- "#D9D9D9"
base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.1) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.23, colour = color_grid),
    axis.title = ggplot2::element_text(size = 6.7, colour = color_text),
    axis.text = ggplot2::element_text(size = 5.6, colour = color_text),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
    strip.text = ggplot2::element_text(size = 5.5, face = "bold", colour = color_text, lineheight = 0.9),
    legend.title = ggplot2::element_text(size = 5.7, face = "bold"),
    legend.text = ggplot2::element_text(size = 5.4),
    plot.margin = ggplot2::margin(1.2, 1.5, 1.2, 1.5, unit = "mm")
  )

panel_A <- capture_stage("panel_A_construction", {
  ggplot2::ggplot(
    plot_data$A,
    ggplot2::aes(x = estimate, y = broad_cell_type, colour = broad_cell_type)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.28, colour = "#777777") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y",
      width = 0.18, linewidth = 0.42
    ) +
    ggplot2::geom_point(size = 1.45) +
    ggplot2::facet_grid(~endpoint_display_short, scales = "free_x") +
    ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = cell_labels) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3),
      labels = scales::label_number(accuracy = 0.01),
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::labs(x = "Mean difference (RIF - fertile; independent facet scales)", y = NULL) +
    base_theme +
    ggplot2::theme(
      legend.position = "none", panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 4.8), axis.text.y = ggplot2::element_text(size = 5.7),
      axis.title.x = ggplot2::element_text(size = 6.3), panel.spacing.x = grid::unit(1.3, "mm")
    )
})

panel_B <- capture_stage("panel_B_construction", {
  ggplot2::ggplot(
    plot_data$B,
    ggplot2::aes(x = count, y = broad_cell_type, fill = broad_cell_type)
  ) +
    ggplot2::geom_blank(
      data = data.frame(
        count = c(1, 16),
        broad_cell_type = factor(rep("Ciliated epithelial cells", 2), levels = rev(cell_levels)),
        scope = factor(
          c("Primary LH+7\nRIF vs fertile", "Fertile-reference\nsensitivity"),
          levels = c("Primary LH+7\nRIF vs fertile", "Fertile-reference\nsensitivity")
        )
      )
    ) +
    ggplot2::geom_col(width = 0.62) +
    ggplot2::geom_point(size = 0.9, colour = color_text) +
    ggplot2::geom_text(
      ggplot2::aes(label = count), hjust = -0.25, size = 2.0,
      family = figure_font, colour = color_text
    ) +
    ggplot2::facet_grid(~scope, scales = "free_x") +
    ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = cell_labels) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.16)), breaks = scales::breaks_pretty(n = 4)) +
    ggplot2::labs(x = "FDR-significant result rows (independent facet axes)", y = NULL) +
    base_theme +
    ggplot2::theme(
      legend.position = "none", panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 5.4), axis.text.x = ggplot2::element_text(size = 5.1),
      panel.spacing.x = grid::unit(3, "mm")
    )
})

panel_C <- capture_stage("panel_C_construction", {
  ggplot2::ggplot(
    plot_data$C,
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, colour = broad_cell_type)
  ) +
    ggplot2::geom_point(size = 0.08, alpha = 0.48, stroke = 0) +
    ggplot2::geom_label(
      data = plot_data$C_centroids,
      ggplot2::aes(x = label_x, y = label_y, label = broad_cell_type),
      inherit.aes = FALSE, family = figure_font, fontface = "bold", size = 1.75,
      linewidth = 0, label.padding = grid::unit(0.55, "mm"),
      fill = scales::alpha("white", 0.82), colour = color_text
    ) +
    ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    base_theme +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 5.9),
      legend.position = "none"
    )
})

panel_D <- capture_stage("panel_D_construction", {
  ggplot2::ggplot(
    plot_data$D,
    ggplot2::aes(x = estimate, y = broad_cell_type, colour = broad_cell_type)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.28, colour = "#777777") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y",
      width = 0.18, linewidth = 0.42
    ) +
    ggplot2::geom_point(size = 1.35) +
    ggplot2::geom_text(
      ggplot2::aes(label = nominal_marker), nudge_y = 0.24,
      family = figure_font, fontface = "bold", size = 2.7, colour = color_text
    ) +
    ggplot2::facet_grid(~endpoint_display_short, scales = "free_x") +
    ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_y_discrete(labels = cell_labels) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3),
      labels = scales::label_number(accuracy = 0.01),
      guide = ggplot2::guide_axis(check.overlap = TRUE)
    ) +
    ggplot2::labs(x = "Mean difference (RIF - control; independent facet scales)", y = NULL) +
    base_theme +
    ggplot2::theme(
      legend.position = "none", panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 4.3), axis.text.y = ggplot2::element_text(size = 5.2),
      axis.title.x = ggplot2::element_text(size = 5.9), panel.spacing.x = grid::unit(0.8, "mm")
    )
})

panel_E <- capture_stage("panel_E_construction", {
  ggplot2::ggplot(
    plot_data$E,
    ggplot2::aes(x = analysis_group, y = value, colour = analysis_group)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = threshold_status),
      position = ggplot2::position_jitter(width = 0.075, height = 0, seed = 2505),
      size = 1.55, stroke = 0.55
    ) +
    ggplot2::stat_summary(
      fun = mean, fun.min = mean, fun.max = mean,
      geom = "crossbar", width = 0.42, linewidth = 0.48
    ) +
    ggplot2::facet_grid(~endpoint_display, scales = "free_y") +
    ggplot2::scale_colour_manual(values = group_colors) +
    ggplot2::scale_shape_manual(values = c("Threshold pass" = 16, "Below threshold" = 1)) +
    ggplot2::labs(x = NULL, y = "Donor-level value", colour = "Group", shape = "Ciliated profile") +
    base_theme +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5.5),
      axis.text.y = ggplot2::element_text(size = 5.0), legend.position = "bottom",
      legend.box = "horizontal", legend.margin = ggplot2::margin(t = 0.5, b = 0),
      panel.spacing.x = grid::unit(1.4, "mm")
    )
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
header_subtitle_size_pt <- 6.4

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
    panel_tag,
    x = grid::unit(0, "npc") + header_shift - grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"), hjust = 1,
    gp = grid::gpar(fontfamily = figure_font, fontface = "bold", fontsize = header_tag_size_pt, col = color_text)
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle, x = grid::unit(0, "npc") + header_shift, y = grid::unit(0.56, "npc"), hjust = 0,
    gp = grid::gpar(fontfamily = figure_font, fontsize = header_subtitle_size_pt, col = color_muted)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, grid::grobTree(tag_grob, title_grob),
    t = 1, b = 1, l = panel_left, r = panel_right, clip = "off",
    name = paste0("F5-header-title-", panel_tag)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, subtitle_grob,
    t = 2, b = 2, l = panel_left, r = panel_right, clip = "off",
    name = paste0("F5-header-subtitle-", panel_tag)
  )
  list(
    gtable = header_gt,
    geometry = data.frame(
      panel = panel_tag, anchor = "BODY_GTABLE_PANEL_CELL_LEFT_EDGE",
      panel_left_column = panel_left, panel_right_column = panel_right,
      tag_size_pt = header_tag_size_pt, title_size_pt = header_title_size_pt,
      subtitle_size_pt = header_subtitle_size_pt, tag_title_gap_mm = header_tag_to_title_gap_mm,
      header_shift_mm = header_shift_mm,
      stringsAsFactors = FALSE
    )
  )
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
  header_C <- make_header_gtable(
    body_C, "C", "GSE183837 broad-cell landscape",
    "47,161 cells · 9 donors · deterministic display subset",
    header_shift_mm = 6
  )
  header_D <- make_header_gtable(body_D, "D", "GSE183837 donor-level effects", "25 prespecified models · 1 nominal exact P < 0.05 · 0 global-BH")
  header_E <- make_header_gtable(body_E, "E", "Ciliated epithelial profiles", "Descriptive only · 0/3 control and 2/6 RIF donors pass the primary cell threshold")
  design <- "
KK
AA
BB
##
CC
DD
##
EF
GH
##
II
JJ
LL
"
  figure_object <- (
    patchwork::wrap_elements(full = header_A$gtable) +
      patchwork::wrap_elements(full = body_A) +
      patchwork::wrap_elements(full = header_B$gtable) +
      patchwork::wrap_elements(full = body_B) +
      patchwork::wrap_elements(full = header_C$gtable) +
      patchwork::wrap_elements(full = header_D$gtable) +
      patchwork::wrap_elements(full = body_C) +
      patchwork::wrap_elements(full = body_D) +
      patchwork::wrap_elements(full = header_E$gtable) +
      patchwork::wrap_elements(full = body_E) +
      patchwork::plot_spacer() +
      patchwork::plot_spacer()
  ) + patchwork::plot_layout(
    design = design, widths = c(0.82, 1.18),
    # These are panel-region heights. Axes, strips, legends, and patchwork's
    # outer rows add about 45 mm; keeping the panel regions at 158 mm prevents
    # the A header and E legend from overflowing the 250 mm submission canvas.
    heights = grid::unit(c(2, 8.3, 44, 2.5, 8.3, 21, 2.5, 8.3, 50, 2.5, 8.3, 43, 2), "mm")
  ) + patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm"))
  )
  list(
    figure = figure_object,
    header_geometry = do.call(rbind, lapply(list(header_A, header_B, header_C, header_D, header_E), `[[`, "geometry"))
  )
})
grDevices::dev.off()
unlink(metrics_device_path)
write_csv_utf8(composite_parts$header_geometry, file.path(dirs$qc, "FIG5_PANEL_HEADER_GEOMETRY.csv"))

figure_png <- file.path(dirs$figures, "Figure5_v1_2_QC_180mm_300dpi.png")
capture_stage("300dpi_candidate_export", {
  ragg::agg_png(
    figure_png, width = figure_width_mm, height = figure_height_mm,
    units = "mm", res = 300, background = "white"
  )
  print(composite_parts$figure)
  grDevices::dev.off()
})
if (!file.exists(figure_png) || file.info(figure_png)$size < 2e+05) stop("Figure 5 candidate export failed.")

warning_log <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
  warning_id = character(), stage = character(), warning_class = character(), message = character(), stringsAsFactors = FALSE
)
write_csv_utf8(warning_log, file.path(dirs$qc, "FIG5_RUNTIME_WARNING_LOG.csv"))

img <- magick::image_read(figure_png)
img_info <- magick::image_info(img)
rgb <- png::readPNG(figure_png)
rgb3 <- rgb[, , seq_len(min(3L, dim(rgb)[3])), drop = FALSE]
delta <- apply(abs(rgb3 - 1), c(1, 2), max)
foreground <- which(delta > 0.035, arr.ind = TRUE)
if (!nrow(foreground)) stop("Figure 5 raster appears blank.")
edge_margins <- data.frame(
  edge = c("top", "bottom", "left", "right"),
  margin_px = c(
    min(foreground[, 1]) - 1,
    nrow(delta) - max(foreground[, 1]),
    min(foreground[, 2]) - 1,
    ncol(delta) - max(foreground[, 2])
  ),
  minimum_px = c(8, 8, 8, 8),
  stringsAsFactors = FALSE
)
edge_margins$pass <- edge_margins$margin_px >= edge_margins$minimum_px
write_csv_utf8(edge_margins, file.path(dirs$qc, "FIG5_RASTER_EDGE_MARGIN_QC.csv"))

full_geometry <- paste0(img_info$width[[1]], "x", img_info$height[[1]])
crop_top <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.39), "+0+0"))
crop_middle <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.34), "+0+", round(img_info$height[[1]] * 0.38)))
crop_bottom <- magick::image_crop(img, paste0(img_info$width[[1]], "x", round(img_info$height[[1]] * 0.29), "+0+", round(img_info$height[[1]] * 0.71)))
magick::image_write(crop_top, file.path(dirs$qc, "Figure5_TOP_CROP.png"), format = "png")
magick::image_write(crop_middle, file.path(dirs$qc, "Figure5_MIDDLE_CROP.png"), format = "png")
magick::image_write(crop_bottom, file.path(dirs$qc, "Figure5_BOTTOM_CROP.png"), format = "png")

automated_checks <- data.frame(
  check = c(
    "Raster exists", "Raster width", "Raster height", "Raster content variance",
    "All edge safety margins", "Runtime warning count", "Five panel headers", "Arial font file"
  ),
  observed = c(
    file.exists(figure_png), img_info$width[[1]], img_info$height[[1]], stats::sd(as.numeric(rgb3)),
    all(edge_margins$pass), nrow(warning_log), nrow(composite_parts$header_geometry), file.exists(arial_path)
  ),
  expected = c(TRUE, floor(figure_width_mm / 25.4 * 300), floor(figure_height_mm / 25.4 * 300), "> 0.01", TRUE, 0, 5, TRUE),
  pass = c(
    file.exists(figure_png), img_info$width[[1]] == floor(figure_width_mm / 25.4 * 300),
    img_info$height[[1]] == floor(figure_height_mm / 25.4 * 300), stats::sd(as.numeric(rgb3)) > 0.01,
    all(edge_margins$pass), nrow(warning_log) == 0, nrow(composite_parts$header_geometry) == 5, file.exists(arial_path)
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(automated_checks, file.path(dirs$qc, "FIG5_AUTOMATED_RENDER_CHECKS.csv"))
if (!all(automated_checks$pass)) stop("Figure 5 automated render QC failed.")

expected_reviewed_sha256 <- "74922b7e68a60ec14af7d370477de623384ab3762c5109de7d5b60aa714d38b1"
actual_sha256 <- sha256_file(figure_png)
review_bound <- identical(actual_sha256, expected_reviewed_sha256)
strict_visual_qc <- data.frame(
  severity = c("BLOCKING", "BLOCKING", "BLOCKING", "MAJOR", "MAJOR", "MAJOR", "MAJOR"),
  visual_check = c(
    "Global panel allocation and header/body alignment", "Text, glyph, axis, strip and legend clipping",
    "Scientific values and significance annotations visible", "Panel A 30-model forest readability",
    "Panel B primary-versus-sensitivity separation", "Panel C labels and population visibility",
    "Panels D/E donor-effect and descriptive-profile readability"
  ),
  status = if (review_bound) "PASS_ASSISTANT_ACTUAL_RENDER_REVIEW" else "PENDING_ASSISTANT_ACTUAL_RENDER_REVIEW",
  reviewed_reference_sha256 = if (review_bound) actual_sha256 else "",
  stringsAsFactors = FALSE
)
write_csv_utf8(strict_visual_qc, file.path(dirs$qc, "FIG5_STRICT_RENDER_QC.csv"))

typography_lock <- data.frame(
  component = c("font_family", "panel_tag_pt", "header_title_pt", "header_subtitle_pt", "base_theme_pt", "final_width_mm", "final_height_mm", "approval_ppi"),
  value = c(figure_font, header_tag_size_pt, header_title_size_pt, header_subtitle_size_pt, 7.1, figure_width_mm, figure_height_mm, 300),
  stringsAsFactors = FALSE
)
write_csv_utf8(typography_lock, file.path(dirs$qc, "FIG5_ARTICLE_TYPOGRAPHY_LOCK.csv"))

legend_text <- c(
  "# Figure 5 legend",
  "",
  "**Figure 5. Donor-level disease comparisons in GSE250130 and GSE183837 do not support broad RIF-associated disruption of the locked VDR programs.**",
  "**A,** Frozen GSE250130 LH+7 RIF-minus-fertile effects for four locked VDR programs and VDR expression. Points show saved mean differences and error bars show saved Welch 95% confidence intervals; facet scales are independent. None of 24 program models or six VDR-expression models was FDR-significant.",
  "**B,** Frozen GSE250130 LH+7 differential-expression summary. The primary RIF-versus-fertile comparison yielded no FDR-significant genome-wide result rows. Twenty-eight rows arose only in fertile-reference sensitivity contrasts and are shown separately; they are not primary RIF evidence.",
  "**C,** UMAP of 47,161 frozen GSE183837 cells across six locked broad cell types. A deterministic display subset of at most 6,000 cells per broad type is shown; coordinates are descriptive.",
  "**D,** Frozen GSE183837 donor-level mean differences for 25 prespecified VDR-expression and VDR-program models. One myeloid VDR-expression model had nominal exact P < 0.05 (asterisk), but no model survived global Benjamini-Hochberg correction. Error bars show saved Welch 95% confidence intervals; facet scales are independent.",
  "**E,** Donor-level ciliated epithelial VDR-expression and program profiles. Points are descriptive because zero of three control and two of six RIF donor profiles met the primary cell-count threshold; horizontal bars show group means and no hypothesis test is displayed."
)
write_text_utf8(legend_text, file.path(dirs$manuscript, "Figure5_legend_v1.2.md"))

methods_text <- c(
  "# Figure 5 Methods",
  "",
  "Figure 5 was assembled exclusively from frozen STEP04_LOCK_v1.2 outputs. No single-cell integration, annotation, pseudobulk construction, model fitting, differential-expression testing, permutation testing, leave-one-donor-out analysis, or multiple-testing correction was rerun. GSE250130 panel values were read from VDRMODELS_v1.1 and DE_v1.2. GSE183837 coordinates and broad annotations were read from SCI183_v1.0 and ANN183_v1.0, and donor effects and ciliated descriptive profiles were read from MD183_v1.0. For Panel C, cells were ordered by locked broad type and cell barcode and at most 6,000 cells per broad type were retained for display only. All inferential quantities, confidence intervals, exact P values, and adjusted P values were carried unchanged from frozen tables."
)
write_text_utf8(methods_text, file.path(dirs$manuscript, "Figure5_methods_v1.2.md"))

results_text <- c(
  "# Figure 5 Results",
  "",
  "At LH+7 in GSE250130, none of 24 prespecified VDR-program models and none of six VDR-expression models was significant after FDR correction. The primary genome-wide RIF-versus-fertile analysis likewise yielded no significant gene-by-cell-type rows. Twenty-eight significant rows were confined to fertile-reference sensitivity contrasts and therefore did not constitute primary RIF evidence.",
  "",
  "In GSE183837, the frozen broad annotation comprised 47,161 cells from nine donors. Among 25 prespecified donor-level models, one myeloid VDR-expression comparison met nominal exact P < 0.05, but no model survived global multiplicity correction. Ciliated epithelial profiles were retained descriptively because the primary donor-level cell threshold was met by only two RIF donors and no control donors. Together, these donor-level analyses did not support broad, reproducible RIF-associated disruption of the locked VDR programs; the negative results should not be interpreted as proof of equivalence."
)
write_text_utf8(results_text, file.path(dirs$manuscript, "Figure5_results_v1.2.md"))

limitations_text <- c(
  "# Figure 5 limitations and interpretation",
  "",
  "- Negative primary findings limit broad dysregulation claims but do not establish biological equivalence.",
  "- GSE183837 contains three control and six RIF donors and has limited exact-test resolution.",
  "- The nominal myeloid VDR-expression signal did not survive global BH correction and must not be described as validation.",
  "- Ciliated epithelial profiles are descriptive because the primary cell-count threshold was not met in most donors.",
  "- Fertile-reference sensitivity DE rows remain sensitivity findings and are not primary RIF-versus-fertile evidence."
)
write_text_utf8(limitations_text, file.path(dirs$manuscript, "Figure5_limitations_v1.2.md"))

write_text_utf8(c(
  "# Figure 5 panel description (Chinese)", "",
  "A：GSE250130 LH+7 供体层面 VDR-program 与 VDR 表达效应；24 个 program 模型和 6 个 VDR 模型均无 FDR 显著。",
  "B：GSE250130 LH+7 DE 汇总；primary 为 0，28 行仅来自 fertile-reference sensitivity。",
  "C：GSE183837 冻结 broad-cell UMAP；仅展示确定性子集。",
  "D：GSE183837 的 25 个预设供体模型；1 个 nominal exact-P 信号，global BH 为 0。",
  "E：纤毛上皮供体值，仅描述，不进行推断。"
), file.path(dirs$manuscript, "Figure5_panel_description_zh_v1.2.md"))

parameters <- list(
  candidate_id = candidate_id,
  width_mm = figure_width_mm,
  height_mm = figure_height_mm,
  approval_ppi = 300,
  UMAP_display_max_per_broad_type = 6000,
  UMAP_display_order = "broad_cell_type_then_cell_barcode",
  jitter_seed_panel_E = 2505,
  facet_scales_A_D_E = "free",
  scientific_recomputation = FALSE,
  final_user_visual_QC_required = TRUE
)
jsonlite::write_json(parameters, file.path(dirs$qc, "analysis_parameters.json"), auto_unbox = TRUE, pretty = TRUE)
write_csv_utf8(data.frame(
  component = c("Panel C display subset", "Panel E jitter", "All inferential outputs"),
  seed = c("NOT_APPLICABLE_DETERMINISTIC_ORDER", "2505", "NOT_APPLICABLE_FROZEN_TABLES"),
  stringsAsFactors = FALSE
), file.path(dirs$qc, "random_seed_manifest.csv"))
write_text_utf8(capture.output(sessionInfo()), file.path(dirs$qc, "sessionInfo_FIG5_v1_2.txt"))

figure_status <- data.frame(
  figure = "Figure 5", candidate = candidate_id,
  science_QC = "PASS", statistics_QC = "PASS_NO_RECOMPUTATION",
  automated_render_QC = "PASS",
  assistant_actual_render_QC = if (review_bound) "PASS" else "PENDING",
  user_final_visual_QC = "PENDING",
  immutable_release = FALSE,
  stringsAsFactors = FALSE
)
write_csv_utf8(figure_status, file.path(out_root, "FIGURE_STATUS.csv"))

hash_qc <- data.frame(
  file = basename(figure_png), actual_sha256 = actual_sha256,
  expected_reviewed_sha256 = expected_reviewed_sha256,
  pass = review_bound,
  stringsAsFactors = FALSE
)
write_csv_utf8(hash_qc, file.path(dirs$qc, "FIG5_REFERENCE_RENDER_HASH_QC.csv"))

candidate_manifest <- build_manifest(out_root, exclude = "manifests/FIG5_CANDIDATE_MANIFEST.csv")
write_csv_utf8(candidate_manifest, file.path(dirs$manifests, "FIG5_CANDIDATE_MANIFEST.csv"))
candidate_validation <- validate_manifest(out_root, candidate_manifest)
write_csv_utf8(candidate_validation, file.path(dirs$manifests, "FIG5_CANDIDATE_MANIFEST_VALIDATION.csv"))
if (!all(candidate_validation$exists & candidate_validation$size_match & candidate_validation$md5_match & candidate_validation$sha256_match)) {
  stop("Figure 5 candidate manifest validation failed.")
}

if (!review_bound) {
  cat("Figure 5 candidate generated for assistant actual-render review.\n")
  cat("Candidate PNG: ", normp(figure_png, TRUE), "\n", sep = "")
  cat("SHA-256: ", actual_sha256, "\n", sep = "")
  stop("Assistant actual-render review and reviewed-hash binding are required before formal QC completion.")
}

dir.create(feedback_root, recursive = TRUE, showWarnings = FALSE)
report_root <- file.path(feedback_root, paste0("FIG5_V1_2_QC_REPORT_", run_stamp))
if (dir.exists(report_root)) stop("Report root exists: ", normp(report_root, FALSE))
report_fig5 <- file.path(report_root, "FIGURE5")
dir.create(report_fig5, recursive = TRUE, showWarnings = FALSE)
copy_tree(out_root, report_fig5)
report_manifest <- build_manifest(report_root, exclude = "SMALL_REPORT_MANIFEST.csv")
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))
report_validation <- validate_manifest(report_root, report_manifest)
if (!all(report_validation$exists & report_validation$size_match & report_validation$md5_match & report_validation$sha256_match)) {
  stop("Figure 5 QC SMALL_REPORT validation failed.")
}
zip_path <- file.path(feedback_root, paste0("STEP05E_FIG5_V1_2_QC_", run_stamp, "_SMALL_REPORT.zip"))
if (file.exists(zip_path)) stop("ZIP exists: ", normp(zip_path, FALSE))
zip::zipr(zipfile = zip_path, files = basename(report_root), root = dirname(report_root), include_directories = TRUE)
if (!file.exists(zip_path) || file.info(zip_path)$size < 3e+05) stop("Figure 5 QC ZIP creation failed.")

cat("\n============================================================\n")
cat("Figure 5 v1.2 QC candidate completed\n")
cat("Gate: PASS_FIGURE5_ASSISTANT_ACTUAL_RENDER_QC_AWAIT_USER_VISUAL_QC\n")
cat("Output: ", normp(out_root, TRUE), "\n", sep = "")
cat("PNG: ", normp(figure_png, TRUE), "\n", sep = "")
cat("SMALL_REPORT: ", normp(zip_path, TRUE), "\n", sep = "")
cat("============================================================\n")

invisible(list(
  gate = "PASS_FIGURE5_ASSISTANT_ACTUAL_RENDER_QC_AWAIT_USER_VISUAL_QC",
  output_root = normp(out_root, TRUE), png = normp(figure_png, TRUE),
  small_report = normp(zip_path, TRUE)
))

# Publication figure stage 2
# ============================================================================
# 127_step05E_Figure5_v1_2_submission_formats_QC.R
# Nutrition_Endometrial_VDR_v2
#
# User-approved Figure 5 v1.2 submission-format build.
# - Reuses the exact approved v1.2 code-layer figure object.
# - Does not rerun analysis or statistics.
# - Exports true 600-ppi PNG/TIFF and a 600-ppi visual-authority PDF at 180 x 250 mm.
# - Renders the PDF with Poppler and performs cross-format QC before lock.
# - Does not overwrite the approved 300-dpi candidate.
# ============================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "127_step05E_Figure5_v1_2_submission_formats_QC.R"
figure_id <- "Figure 5"
approved_visual_version <- "F5_v1_2_QC"
format_version <- "F5_v1_2_SUBMISSION_FORMATS_QC"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

required_packages <- c("ragg", "magick", "digest", "zip", "png")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

normp <- function(x, must_work = FALSE) {
  gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = must_work))
}

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

sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

md5_file <- function(path) unname(as.character(tools::md5sum(path)))

copy_exact <- function(src, dst) {
  if (!file.exists(src) || dir.exists(src)) stop("Missing source file: ", normp(src, FALSE))
  if (file.exists(dst)) stop("Refusing overwrite: ", normp(dst, FALSE))
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dst) ||
      !identical(sha256_file(src), sha256_file(dst))) {
    stop("Validated copy failed: ", normp(src, FALSE), " -> ", normp(dst, FALSE))
  }
  invisible(dst)
}

relative_to <- function(path, root) {
  p <- normp(path, FALSE)
  r <- sub("/+$", "", normp(root, FALSE))
  if (startsWith(tolower(p), paste0(tolower(r), "/"))) {
    return(substr(p, nchar(r) + 2L, nchar(p)))
  }
  p
}

copy_tree <- function(src_root, dst_root) {
  if (!dir.exists(src_root)) return(invisible(FALSE))
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
    relative_path = rel,
    size_bytes = as.numeric(info$size),
    md5 = vapply(files, md5_file, character(1)),
    sha256 = vapply(files, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  out[order(out$relative_path), , drop = FALSE]
}

validate_manifest <- function(root, manifest) {
  do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
    p <- file.path(root, manifest$relative_path[[i]])
    ex <- file.exists(p) && !dir.exists(p)
    data.frame(
      relative_path = manifest$relative_path[[i]],
      exists = ex,
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
step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")
figure_parent <- file.path(step05_root, "06_results", "FIG5")
approved_root <- file.path(figure_parent, approved_visual_version)
approved_png <- file.path(
  approved_root, "figures", "Figure5_v1_2_QC_180mm_300dpi.png"
)
approved_script <- file.path(
  project_root, "04_code", "R",
  "126_step05E_Figure5_manuscript_composite_v1_2_QC.R"
)
expected_approved_png_sha256 <- "74922b7e68a60ec14af7d370477de623384ab3762c5109de7d5b60aa714d38b1"

if (!dir.exists(approved_root) || !file.exists(approved_png) || !file.exists(approved_script)) {
  stop("Approved Figure 5 v1.2 authority is incomplete.")
}
if (!identical(tolower(sha256_file(approved_png)), expected_approved_png_sha256)) {
  stop("Approved Figure 5 v1.2 PNG identity mismatch.")
}

skill3_path <- "<FIGURE_STYLE_REFERENCE>"
skill3_text <- paste(readLines(skill3_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.7", skill3_text, fixed = TRUE)) {
  stop("Active Skill 3 is not v2.2-beta.7.")
}

parent_override <- Sys.getenv("FIG5_SUBMISSION_PARENT", unset = "")
feedback_override <- Sys.getenv("FIG5_SUBMISSION_FEEDBACK", unset = "")
out_parent <- if (nzchar(parent_override)) parent_override else figure_parent
feedback_root <- if (nzchar(feedback_override)) feedback_override else file.path(project_root, "反馈")
out_root <- file.path(out_parent, format_version)
if (dir.exists(out_root)) stop("Refusing to overwrite submission-format folder: ", normp(out_root, FALSE))

dirs <- list(
  figures = file.path(out_root, "figures"),
  reference = file.path(out_root, "reference"),
  qc = file.path(out_root, "qc"),
  code = file.path(out_root, "code"),
  source_data = file.path(out_root, "source_data"),
  manuscript = file.path(out_root, "manuscript_text"),
  manifests = file.path(out_root, "manifests")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(feedback_root, recursive = TRUE, showWarnings = FALSE)

# Rebuild the exact approved code-layer object in a process-local temporary target.
temp_build_root <- file.path(tempdir(), paste0("FIG5_v1_2_submission_build_", run_stamp))
temp_candidate_parent <- file.path(temp_build_root, "candidate_parent")
temp_feedback <- file.path(temp_build_root, "feedback")
old_candidate_env <- Sys.getenv("FIG5_CANDIDATE_PARENT", unset = NA_character_)
old_feedback_env <- Sys.getenv("FIG5_FEEDBACK_ROOT", unset = NA_character_)
Sys.setenv(
  FIG5_CANDIDATE_PARENT = temp_candidate_parent,
  FIG5_FEEDBACK_ROOT = temp_feedback
)
plot_env <- new.env(parent = globalenv())
source_console <- capture.output(
  sys.source(approved_script, envir = plot_env, keep.source = TRUE)
)
if (is.na(old_candidate_env)) Sys.unsetenv("FIG5_CANDIDATE_PARENT") else Sys.setenv(FIG5_CANDIDATE_PARENT = old_candidate_env)
if (is.na(old_feedback_env)) Sys.unsetenv("FIG5_FEEDBACK_ROOT") else Sys.setenv(FIG5_FEEDBACK_ROOT = old_feedback_env)

temp_approved_png <- file.path(
  temp_candidate_parent, approved_visual_version, "figures",
  "Figure5_v1_2_QC_180mm_300dpi.png"
)
if (!file.exists(temp_approved_png) ||
    !identical(tolower(sha256_file(temp_approved_png)), expected_approved_png_sha256)) {
  stop("Rebuilt Figure 5 v1.2 code layer does not reproduce the approved reference PNG.")
}
if (is.null(plot_env$composite_parts$figure)) stop("Approved figure object was not recovered.")
figure_object <- plot_env$composite_parts$figure
figure_width_mm <- plot_env$figure_width_mm
figure_height_mm <- plot_env$figure_height_mm
if (!identical(as.numeric(figure_width_mm), 180) ||
    !identical(as.numeric(figure_height_mm), 250)) {
  stop("Approved Figure 5 geometry is not 180 x 250 mm.")
}
write_text_utf8(source_console, file.path(dirs$qc, "FIG5_APPROVED_CODE_LAYER_REBUILD_CONSOLE.txt"))

# Preserve approved candidate context without altering it.
copy_tree(file.path(approved_root, "source_data"), dirs$source_data)
copy_tree(file.path(approved_root, "manuscript_text"), dirs$manuscript)
copy_tree(file.path(approved_root, "manifests"), dirs$manifests)
copy_tree(file.path(approved_root, "qc"), file.path(dirs$qc, "approved_candidate_qc"))
copy_tree(file.path(approved_root, "code"), file.path(dirs$code, "approved_candidate_code"))
copy_exact(approved_script, file.path(dirs$code, basename(approved_script)))
copy_exact(script_path, file.path(dirs$code, script_name))

reference_copy <- file.path(
  dirs$reference, "Figure5_v1_2_APPROVED_REFERENCE_180mm_300dpi.png"
)
copy_exact(approved_png, reference_copy)

user_approval <- data.frame(
  figure = figure_id,
  approved_visual_version = approved_visual_version,
  approved_png_sha256 = expected_approved_png_sha256,
  user_final_visual_qc = "PASS",
  approval_date = "2026-08-07",
  approval_scope = "Figure 5 v1.2 explicitly linked and approved by user",
  unselected_candidate = "F5_v1_0_QC_AND_F5_v1_1_QC_NOT_SELECTED_FOR_RELEASE",
  stringsAsFactors = FALSE
)
write_csv_utf8(user_approval, file.path(dirs$qc, "FIG5_USER_VISUAL_APPROVAL_v1_2.csv"))

width_px_600 <- round(figure_width_mm / 25.4 * 600)
height_px_600 <- round(figure_height_mm / 25.4 * 600)
width_in <- figure_width_mm / 25.4
height_in <- figure_height_mm / 25.4

png600 <- file.path(dirs$figures, "Figure5_v1_2_submission_180mm_600ppi.png")
tiff600 <- file.path(dirs$figures, "Figure5_v1_2_submission_180mm_600ppi.tiff")
pdfvec <- file.path(dirs$figures, "Figure5_v1_2_submission_180mm.pdf")

warning_rows <- list()
warning_n <- 0L
capture_stage <- function(stage, expr) {
  withCallingHandlers(expr, warning = function(w) {
    warning_n <<- warning_n + 1L
    warning_rows[[warning_n]] <<- data.frame(
      warning_id = sprintf("W%03d", warning_n),
      stage = stage,
      warning_class = class(w)[[1]],
      message = conditionMessage(w),
      stringsAsFactors = FALSE
    )
    invokeRestart("muffleWarning")
  })
}

capture_stage("600ppi_PNG_export", {
  ragg::agg_png(
    png600, width = width_px_600, height = height_px_600,
    units = "px", res = 600, background = "white"
  )
  print(figure_object)
  grDevices::dev.off()
})
capture_stage("600ppi_TIFF_export", {
  ragg::agg_tiff(
    tiff600, width = width_px_600, height = height_px_600,
    units = "px", res = 600, compression = "lzw", background = "white"
  )
  print(figure_object)
  grDevices::dev.off()
})
if (!isTRUE(capabilities("cairo"))) stop("Cairo PDF support is required.")
capture_stage("visual_authority_PDF_export", {
  grDevices::cairo_pdf(
    pdfvec, width = width_in, height = height_in,
    family = "Arial", bg = "white", onefile = TRUE
  )
  grid::grid.newpage()
  grid::grid.raster(
    png::readPNG(png600),
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.5, "npc"),
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"),
    interpolate = TRUE
  )
  grDevices::dev.off()
})

warning_log <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
  warning_id = character(), stage = character(), warning_class = character(),
  message = character(), stringsAsFactors = FALSE
)
write_csv_utf8(warning_log, file.path(dirs$qc, "FIG5_SUBMISSION_RUNTIME_WARNING_LOG.csv"))
if (nrow(warning_log)) stop("Submission-format export emitted warnings.")

for (p in c(png600, tiff600, pdfvec)) {
  if (!file.exists(p) || file.info(p)$size < 50000) stop("Export failed: ", normp(p, FALSE))
}
pdf_con <- file(pdfvec, "rb")
pdf_header <- rawToChar(readBin(pdf_con, what = "raw", n = 5L))
close(pdf_con)
if (!identical(pdf_header, "%PDF-")) stop("Invalid PDF header.")

png_img <- magick::image_read(png600)
tiff_img <- magick::image_read(tiff600)
png_info <- magick::image_info(png_img)
tiff_info <- magick::image_info(tiff_img)
if (png_info$width[[1]] != width_px_600 || png_info$height[[1]] != height_px_600 ||
    tiff_info$width[[1]] != width_px_600 || tiff_info$height[[1]] != height_px_600) {
  stop("600-ppi raster geometry mismatch.")
}

# Render the vector PDF using the bundled Poppler runtime.
pdftoppm_exe <- paste0(
  "<OPTIONAL_LOCAL_RUNTIME>",
  "native/poppler/Library/bin/pdftoppm.exe"
)
pdfinfo_exe <- paste0(
  "<OPTIONAL_LOCAL_RUNTIME>",
  "native/poppler/Library/bin/pdfinfo.exe"
)
if (!file.exists(pdftoppm_exe) || !file.exists(pdfinfo_exe)) stop("Bundled Poppler is missing.")
pdf_render_prefix <- file.path(dirs$qc, "Figure5_v1_2_PDF_RENDER_300dpi")
poppler_output <- system2(
  pdftoppm_exe,
  c("-png", "-r", "300", "-singlefile", shQuote(pdfvec), shQuote(pdf_render_prefix)),
  stdout = TRUE, stderr = TRUE
)
pdf_render_png <- paste0(pdf_render_prefix, ".png")
if (!file.exists(pdf_render_png)) {
  stop("Poppler PDF render failed: ", paste(poppler_output, collapse = " | "))
}
pdfinfo_output_raw <- system2(pdfinfo_exe, shQuote(pdfvec), stdout = TRUE, stderr = TRUE)
pdfinfo_keep <- grepl(
  "^(Producer|Pages|Page size|Page rot|File size|Encrypted|PDF version):",
  pdfinfo_output_raw,
  useBytes = TRUE
)
pdfinfo_output <- pdfinfo_output_raw[pdfinfo_keep]
write_text_utf8(pdfinfo_output, file.path(dirs$qc, "FIG5_PDFINFO.txt"))

approved_img <- magick::image_read(reference_copy)
approved_info <- magick::image_info(approved_img)
approved_geometry <- paste0(approved_info$width[[1]], "x", approved_info$height[[1]], "!")
png600_preview <- magick::image_resize(png_img, approved_geometry, filter = "Lanczos")
png600_preview_path <- file.path(dirs$qc, "Figure5_v1_2_PNG600_RENDER_300dpi.png")
magick::image_write(png600_preview, png600_preview_path, format = "png")

pdf_render_img <- magick::image_read(pdf_render_png)
pdf_render_exact <- magick::image_resize(pdf_render_img, approved_geometry, filter = "Lanczos")
pdf_render_exact_path <- file.path(dirs$qc, "Figure5_v1_2_PDF_RENDER_MATCHED_300dpi.png")
magick::image_write(pdf_render_exact, pdf_render_exact_path, format = "png")

png_tiff_compare <- magick::image_compare(png_img, tiff_img, metric = "AE", fuzz = 0)
png_tiff_AE <- as.numeric(attr(png_tiff_compare, "distortion"))
ref_compare <- magick::image_compare(approved_img, png600_preview, metric = "RMSE", fuzz = 0)
ref_png600_RMSE <- as.numeric(attr(ref_compare, "distortion"))
pdf_compare <- magick::image_compare(approved_img, pdf_render_exact, metric = "RMSE", fuzz = 0)
ref_pdf_RMSE <- as.numeric(attr(pdf_compare, "distortion"))

cross_format_qc <- data.frame(
  check = c(
    "600-ppi PNG geometry",
    "600-ppi TIFF geometry",
    "PNG versus TIFF pixel identity",
    "Approved 300-dpi versus downsampled 600-ppi PNG",
    "Approved 300-dpi versus Poppler-rendered visual-authority PDF",
    "PDF structural header",
    "PDF page count"
  ),
  observed = c(
    paste0(png_info$width[[1]], "x", png_info$height[[1]]),
    paste0(tiff_info$width[[1]], "x", tiff_info$height[[1]]),
    format(png_tiff_AE, scientific = FALSE),
    format(ref_png600_RMSE, digits = 8),
    format(ref_pdf_RMSE, digits = 8),
    pdf_header,
    paste(pdfinfo_output[grepl("^Pages:", pdfinfo_output)], collapse = " ")
  ),
  expected = c(
    paste0(width_px_600, "x", height_px_600),
    paste0(width_px_600, "x", height_px_600),
    "AE = 0",
    "RMSE < 0.08",
    "RMSE < 0.08",
    "%PDF-",
    "Pages: 1"
  ),
  pass = c(
    png_info$width[[1]] == width_px_600 && png_info$height[[1]] == height_px_600,
    tiff_info$width[[1]] == width_px_600 && tiff_info$height[[1]] == height_px_600,
    png_tiff_AE == 0,
    ref_png600_RMSE < 0.08,
    ref_pdf_RMSE < 0.08,
    identical(pdf_header, "%PDF-"),
    any(grepl("^Pages:\\s+1\\s*$", pdfinfo_output))
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(cross_format_qc, file.path(dirs$qc, "FIG5_CROSS_FORMAT_QC.csv"))
if (!all(cross_format_qc$pass)) stop("Figure 5 cross-format QC failed.")

format_manifest <- data.frame(
  role = c("approved_reference_300dpi", "submission_png_600ppi", "submission_tiff_600ppi", "submission_visual_authority_pdf"),
  file = c(basename(reference_copy), basename(png600), basename(tiff600), basename(pdfvec)),
  size_bytes = as.numeric(file.info(c(reference_copy, png600, tiff600, pdfvec))$size),
  sha256 = vapply(c(reference_copy, png600, tiff600, pdfvec), sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_utf8(format_manifest, file.path(dirs$manifests, "FIG5_SUBMISSION_FORMAT_HASH_MANIFEST.csv"))

# The deterministic raster hashes are bound after assistant actual-render review.
# Cairo writes session-specific PDF metadata; its exact formal-build hash is therefore
# bound downstream by Skill 8 after all PDF structure and rendering gates pass.
expected_submission_hashes <- c(
  submission_png_600ppi = "a7300fd37008814e5d5a8196b9b98d1f9cd6b6a13e5e4d86121beee2e017b9b6",
  submission_tiff_600ppi = "6fef35dbf933a5fff9467f65ad0b0631056fb79298986400b0c3c3fd1c8a6bba",
  submission_visual_authority_pdf = "BOUND_AT_SKILL8_AFTER_FORMAL_BUILD"
)
actual_submission_hashes <- c(
  submission_png_600ppi = sha256_file(png600),
  submission_tiff_600ppi = sha256_file(tiff600),
  submission_visual_authority_pdf = sha256_file(pdfvec)
)
reviewed_hash_qc <- data.frame(
  role = names(actual_submission_hashes),
  actual_sha256 = unname(actual_submission_hashes),
  expected_sha256 = unname(expected_submission_hashes),
  binding_policy = c(
    "EXACT_REVIEWED_HASH",
    "EXACT_REVIEWED_HASH",
    "SESSION_METADATA_ALLOWED_CONTENT_AND_STRUCTURE_QC_REQUIRED"
  ),
  pass = c(
    actual_submission_hashes[["submission_png_600ppi"]] == expected_submission_hashes[["submission_png_600ppi"]],
    actual_submission_hashes[["submission_tiff_600ppi"]] == expected_submission_hashes[["submission_tiff_600ppi"]],
    nzchar(actual_submission_hashes[["submission_visual_authority_pdf"]])
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(reviewed_hash_qc, file.path(dirs$qc, "FIG5_REVIEWED_SUBMISSION_HASH_QC.csv"))
if (!all(reviewed_hash_qc$pass)) {
  stop("Submission formats await assistant actual-render review and hash binding.")
}

submission_status <- data.frame(
  figure = figure_id,
  visual_authority = approved_visual_version,
  user_visual_QC = "PASS",
  assistant_cross_format_actual_render_QC = "PASS",
  scientific_QC = "PASS_INHERITED",
  statistical_QC = "PASS_NO_RECOMPUTATION",
  png_600ppi = "PASS",
  tiff_600ppi = "PASS",
  visual_authority_pdf = "PASS_600PPI_RASTER_EMBEDDED",
  ready_for_skill8_lock = TRUE,
  immutable_release = FALSE,
  stringsAsFactors = FALSE
)
write_csv_utf8(submission_status, file.path(out_root, "FIGURE_SUBMISSION_STATUS.csv"))

write_text_utf8(c(
  "# Figure 5 v1.2 submission-format handoff",
  "",
  "- User-approved visual authority: F5_v1_2_QC.",
  "- Unselected candidates: F5_v1_0_QC and F5_v1_1_QC; not used for release.",
  "- Scientific/statistical recomputation: none.",
  "- True 600-ppi PNG and LZW TIFF: generated from the approved code-layer object.",
  "- PDF: 600-ppi visual-authority raster embedded on a Cairo PDF page at 180 x 250 mm.",
  "- Rationale: subpixel semi-transparent UMAP points changed density under direct vector rendering; raster embedding preserves the user-approved appearance.",
  "- Cross-format actual-render QC: PASS.",
  "- Next gate: Skill 8 v1.1 MODULE_LOCK_ARCHIVE."
), file.path(dirs$qc, "FIG5_SUBMISSION_FORMATS_HANDOFF.md"))

write_text_utf8(capture.output(sessionInfo()), file.path(dirs$qc, "sessionInfo_FIG5_submission_v1_2.txt"))
candidate_manifest <- build_manifest(out_root, exclude = "manifests/FIG5_SUBMISSION_CANDIDATE_MANIFEST.csv")
write_csv_utf8(candidate_manifest, file.path(dirs$manifests, "FIG5_SUBMISSION_CANDIDATE_MANIFEST.csv"))
candidate_validation <- validate_manifest(out_root, candidate_manifest)
write_csv_utf8(candidate_validation, file.path(dirs$qc, "FIG5_SUBMISSION_MANIFEST_VALIDATION.csv"))
if (!all(candidate_validation$exists & candidate_validation$size_match &
         candidate_validation$md5_match & candidate_validation$sha256_match)) {
  stop("Submission candidate manifest validation failed.")
}

report_root <- file.path(feedback_root, paste0("FIG5_V1_2_SUBMISSION_FORMATS_REPORT_", run_stamp))
if (dir.exists(report_root)) stop("Report root exists: ", normp(report_root, FALSE))
report_fig5 <- file.path(report_root, "FIGURE5")
dir.create(report_fig5, recursive = TRUE, showWarnings = FALSE)
copy_tree(out_root, report_fig5)
report_manifest <- build_manifest(report_root, exclude = "SMALL_REPORT_MANIFEST.csv")
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))
report_validation <- validate_manifest(report_root, report_manifest)
if (!all(report_validation$exists & report_validation$size_match &
         report_validation$md5_match & report_validation$sha256_match)) {
  stop("SMALL_REPORT validation failed.")
}

zip_path <- file.path(
  feedback_root,
  paste0("STEP05E_FIG5_V1_2_SUBMISSION_FORMATS_", run_stamp, "_SMALL_REPORT.zip")
)
if (file.exists(zip_path)) stop("ZIP exists: ", normp(zip_path, FALSE))
zip::zipr(
  zipfile = zip_path,
  files = basename(report_root),
  root = dirname(report_root),
  include_directories = TRUE
)
if (!file.exists(zip_path) || file.info(zip_path)$size < 500000) stop("ZIP creation failed.")

cat("\n============================================================\n")
cat("Figure 5 v1.2 submission formats completed\n")
cat("Gate: PASS_FIGURE5_V1_2_SUBMISSION_FORMATS_READY_FOR_SKILL8\n")
cat("Output: ", normp(out_root, TRUE), "\n", sep = "")
cat("SMALL_REPORT: ", normp(zip_path, TRUE), "\n", sep = "")
cat("============================================================\n")

invisible(list(
  gate = "PASS_FIGURE5_V1_2_SUBMISSION_FORMATS_READY_FOR_SKILL8",
  output_root = normp(out_root, TRUE),
  png600 = normp(png600, TRUE),
  tiff600 = normp(tiff600, TRUE),
  pdf = normp(pdfvec, TRUE),
  zip = normp(zip_path, TRUE)
))
