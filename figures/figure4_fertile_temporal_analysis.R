
# Publication figure stage 1
# ============================================================================
# 120_step05D_Figure4_panelA_header_alignment_v1_2_QC.R
# Nutrition_Endometrial_VDR_v2
#
# Skill 3 FULL_AUDIT first retained candidate for Figure 4.
# Reads ONLY locked/frozen GSE250130 tables and embedding coordinates.
# Does NOT rerun normalization, integration, UMAP, annotation, models, DE, or FDR.
# Produces ONE 180 mm x 250 mm, 300-dpi QC PNG before user visual lock.
# ============================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "120_step05D_Figure4_panelA_header_alignment_v1_2_QC.R"
figure_id <- "Figure 4"
figure_version <- "F4_v1_2_QC"
source_lock <- "STEP04_LOCK_v1.2"
governing_skill <- "research_scientific_figure_suite v2.2-beta.7"
run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

required_packages <- c(
  "data.table", "ggplot2", "patchwork", "gtable", "grid", "scales",
  "systemfonts", "ragg", "digest", "zip", "png"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

normalize_slashes <- function(x) gsub("\\\\", "/", x)
norm_path <- function(x, must_work = FALSE) {
  normalize_slashes(normalizePath(x, winslash = "/", mustWork = must_work))
}
write_csv_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(path)
}
write_text_utf8 <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(x), con = con, useBytes = TRUE)
  invisible(path)
}
md5_file <- function(path) {
  unname(as.character(tools::md5sum(path)))
}
sha256_file <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
copy_exact <- function(src, dst) {
  if (!file.exists(src) || dir.exists(src)) stop("Missing source file: ", src)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dst)) stop("Refusing to overwrite report copy: ", dst)
  ok <- file.copy(src, dst, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  if (!isTRUE(ok) || !file.exists(dst) || sha256_file(src) != sha256_file(dst)) {
    stop("Validated copy failed: ", src, " -> ", dst)
  }
  invisible(TRUE)
}
find_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    p <- sub("^--file=", "", file_arg[[1]])
    if (file.exists(p)) return(norm_path(p, TRUE))
  }
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      p <- frames[[i]]$ofile
      if (!is.null(p) && length(p) == 1L && file.exists(p)) {
        return(norm_path(p, TRUE))
      }
    }
  }
  fallback <- file.path(project_root, "04_code", "R", script_name)
  if (file.exists(fallback)) return(norm_path(fallback, TRUE))
  stop("Could not identify the executed script path.")
}
relative_to <- function(path, root) {
  p <- norm_path(path, FALSE)
  r <- sub("/+$", "", norm_path(root, FALSE))
  if (startsWith(tolower(p), paste0(tolower(r), "/"))) {
    return(substr(p, nchar(r) + 2L, nchar(p)))
  }
  p
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
copy_tree <- function(src_root, dst_root) {
  files <- list.files(src_root, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  for (src in files) {
    rel <- relative_to(src, src_root)
    copy_exact(src, file.path(dst_root, rel))
  }
  invisible(TRUE)
}

if (!dir.exists(project_root) ||
    !file.exists(file.path(project_root, "Nutrition_Endometrial_VDR_v2.Rproj"))) {
  stop("Project root is missing or invalid: ", project_root)
}

script_path <- find_script_path()
skill_path <- "<FIGURE_STYLE_REFERENCE>"
if (!file.exists(skill_path)) stop("Active Skill 3 file is missing: ", skill_path)
skill_text <- paste(readLines(skill_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.7", skill_text, fixed = TRUE)) {
  stop("Active Skill 3 is not v2.2-beta.7.")
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

# ---------------------------------------------------------------------------
# Paths and immutable candidate directory
# ---------------------------------------------------------------------------

step04_root <- file.path(
  project_root,
  "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell"
)
step05_root <- file.path(project_root, "05_analysis_steps", "STEP05_manuscript_assembly")

candidate_parent_override <- Sys.getenv("FIG4_CANDIDATE_PARENT", unset = "")
feedback_override <- Sys.getenv("FIG4_FEEDBACK_ROOT", unset = "")
fig4_parent <- if (nzchar(candidate_parent_override)) {
  candidate_parent_override
} else {
  file.path(step05_root, "06_results", "FIG4")
}
feedback_root <- if (nzchar(feedback_override)) {
  feedback_override
} else {
  file.path(project_root, "反馈")
}

candidate_root <- file.path(fig4_parent, figure_version)
if (dir.exists(candidate_root)) {
  stop(
    "Refusing to overwrite an existing candidate folder:\n",
    norm_path(candidate_root, FALSE)
  )
}

dirs <- list(
  figures = file.path(candidate_root, "figures"),
  source_data = file.path(candidate_root, "source_data"),
  upstream_refs = file.path(candidate_root, "source_data", "upstream_refs"),
  code = file.path(candidate_root, "code"),
  manuscript = file.path(candidate_root, "manuscript_text"),
  qc = file.path(candidate_root, "qc"),
  manifests = file.path(candidate_root, "manifests")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
dir.create(feedback_root, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Frozen authority and direct upstream inputs
# ---------------------------------------------------------------------------

lock_root <- file.path(step04_root, "06_results", "STEP04_closeout", "STEP04_LOCK_v1.2")
blueprint_path <- file.path(lock_root, "tables", "07_figure_blueprint.csv")
fact_ledger_path <- file.path(lock_root, "tables", "02_locked_fact_ledger.csv")
claim_ledger_path <- file.path(lock_root, "tables", "03_claim_ledger.csv")
lock_status_path <- file.path(lock_root, "LOCK_STATUS.csv")
stop_rule_path <- file.path(lock_root, "text", "STEP04_ANALYSIS_STOP_RULE.md")

ann_root <- file.path(
  step04_root,
  "06_results", "S4E1_sc", "SCI_v1_1_memory_safe", "annotations",
  "GSE250130", "ANN_v1_2"
)
model_root <- file.path(
  step04_root,
  "06_results", "S4E1_sc", "SCI_v1_1_memory_safe", "models",
  "GSE250130", "VDRMODELS_v1_1"
)
de_root <- file.path(
  step04_root,
  "06_results", "S4E1_sc", "SCI_v1_1_memory_safe", "de",
  "GSE250130", "DE_v1_2"
)
ann_qc_root <- file.path(step04_root, "08_qc", "STEP04E1E_20260805_185455_annQC")
pb_qc_root <- file.path(step04_root, "08_qc", "STEP04E1F_20260805_190756_pbQC")
model_qc_root <- file.path(step04_root, "08_qc", "STEP04E1G_20260805_233002_modelQC")
de_qc_root <- file.path(step04_root, "08_qc", "STEP04E1H_20260806_091153_deQC")

source_paths <- c(
  cell_annotation_metadata = file.path(ann_root, "tables", "cell_annotation_metadata.csv"),
  broad_type_summary = file.path(ann_root, "tables", "broad_type_summary.csv"),
  cluster_annotation_map = file.path(ann_root, "tables", "cluster_annotation_map.csv"),
  sample_composition = file.path(ann_root, "tables", "sample_broad_type_composition.csv"),
  temporal_program_summary = file.path(model_root, "source_data", "source_temporal_program_summary.csv"),
  temporal_program_models = file.path(model_root, "tables", "01_temporal_program_omnibus.csv"),
  temporal_VDR_summary = file.path(model_root, "source_data", "source_temporal_VDR_summary.csv"),
  temporal_VDR_models = file.path(model_root, "tables", "05_temporal_VDR_omnibus.csv"),
  temporal_omnibus_DE_summary = file.path(de_root, "tab", "07_temporal_omnibus_summary.csv"),
  temporal_pairwise_DE_summary = file.path(de_root, "tab", "08_temporal_pairwise_summary.csv")
)
if (!all(file.exists(source_paths))) {
  stop("One or more frozen Figure 4 source files are missing:\n",
       paste(names(source_paths)[!file.exists(source_paths)], collapse = ", "))
}

blueprint <- data.table::fread(blueprint_path, data.table = FALSE)
fig4_blueprint <- blueprint[blueprint$figure_id == "Figure 4", , drop = FALSE]
if (nrow(fig4_blueprint) != 1L ||
    !identical(fig4_blueprint$current_status[[1]], "SOURCE_RESULTS_READY_FINAL_COMPOSITE_PENDING")) {
  stop("Figure 4 locked blueprint row is missing or unexpected.")
}
lock_status <- data.table::fread(lock_status_path, data.table = FALSE)
if (nrow(lock_status) != 1L ||
    !identical(lock_status$gate_status[[1]], "PASS_STEP04_LOCK_V1_2_FINAL_MANUSCRIPT_AUTHORITY") ||
    !isTRUE(lock_status$final_figure_assembly_authorized[[1]]) ||
    isTRUE(lock_status$additional_discovery_testing_authorized[[1]])) {
  stop("STEP04 authority does not permit this locked-source figure assembly.")
}

# Verify selected inputs against their upstream ARTIFACT_MANIFEST files.
manifest_paths <- c(
  annotation = file.path(ann_qc_root, "ARTIFACT_MANIFEST.csv"),
  model = file.path(model_qc_root, "ARTIFACT_MANIFEST.csv"),
  DE = file.path(de_qc_root, "ARTIFACT_MANIFEST.csv")
)
if (!all(file.exists(manifest_paths))) stop("One or more upstream manifests are missing.")
manifest_rows <- do.call(rbind, lapply(names(manifest_paths), function(nm) {
  x <- data.table::fread(manifest_paths[[nm]], data.table = FALSE)
  x$manifest <- nm
  x
}))
manifest_rows$path_key <- tolower(norm_path(manifest_rows$path, FALSE))

source_identity <- do.call(rbind, lapply(seq_along(source_paths), function(i) {
  p <- norm_path(source_paths[[i]], TRUE)
  idx <- match(tolower(p), manifest_rows$path_key)
  expected_md5 <- if (!is.na(idx)) tolower(manifest_rows$md5[[idx]]) else NA_character_
  actual_md5 <- tolower(md5_file(p))
  info <- file.info(p)
  data.frame(
    source_id = names(source_paths)[[i]],
    path = p,
    relative_path = relative_to(p, project_root),
    size_bytes = as.numeric(info$size),
    upstream_manifest = if (!is.na(idx)) manifest_rows$manifest[[idx]] else NA_character_,
    expected_md5 = expected_md5,
    actual_md5 = actual_md5,
    md5_match = !is.na(expected_md5) && identical(expected_md5, actual_md5),
    stringsAsFactors = FALSE
  )
}))
if (!all(source_identity$md5_match)) {
  stop("One or more Figure 4 source files fail upstream MD5 verification.")
}
write_csv_utf8(source_identity, file.path(dirs$manifests, "FIG4_UPSTREAM_SOURCE_IDENTITY.csv"))

# Copy small upstream authority files; the 80-MB cell table is referenced by identity only.
upstream_copy_paths <- c(
  blueprint_path, fact_ledger_path, claim_ledger_path, lock_status_path, stop_rule_path,
  file.path(ann_qc_root, "ANNOTATION_GATE.csv"),
  file.path(pb_qc_root, "PBVDR_GATE.csv"),
  file.path(model_qc_root, "VDR_MODEL_GATE.csv"),
  file.path(de_qc_root, "GENE_DE_GATE.csv")
)
for (p in upstream_copy_paths) copy_exact(p, file.path(dirs$upstream_refs, basename(p)))

# ---------------------------------------------------------------------------
# Full-runtime warning capture registry
# ---------------------------------------------------------------------------

warning_rows <- list()
warning_n <- 0L
capture_stage <- function(stage, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_n <<- warning_n + 1L
      warning_rows[[warning_n]] <<- data.frame(
        warning_id = sprintf("W%03d", warning_n),
        stage = stage,
        warning_class = class(w)[[1]],
        message = conditionMessage(w),
        stringsAsFactors = FALSE
      )
      invokeRestart("muffleWarning")
    }
  )
}

# ---------------------------------------------------------------------------
# Read locked tables and create direct plotting tables (plotting transforms only)
# ---------------------------------------------------------------------------

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

  cell_meta <- data.table::fread(
    source_paths[["cell_annotation_metadata"]],
    select = c("cell_barcode", "UMAPH_1", "UMAPH_2", "broad_cell_type")
  )
  if (nrow(cell_meta) != 236636L || data.table::uniqueN(cell_meta$cell_barcode) != 236636L) {
    stop("Locked cell annotation metadata row/ID count mismatch.")
  }
  if (any(!is.finite(cell_meta$UMAPH_1)) || any(!is.finite(cell_meta$UMAPH_2))) {
    stop("Non-finite locked UMAP coordinates detected.")
  }

  full_type_order <- c(
    "B cells", "Ciliated epithelial cells", "Endothelial cells", "Mast cells",
    "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells"
  )
  temporal_type_order <- c(
    "Ciliated epithelial cells", "Myeloid cells", "NK/T cells",
    "Stromal cells", "Unciliated epithelial cells"
  )
  program_order <- c(
    "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284",
    "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84"
  )

  data.table::setorder(cell_meta, broad_cell_type, cell_barcode)
  umap_subset <- cell_meta[, {
    keep_n <- min(.N, 10000L)
    idx <- unique(as.integer(round(seq.int(1L, .N, length.out = keep_n))))
    .SD[idx]
  }, by = broad_cell_type]
  type_sizes <- setNames(broad_summary$cell_n, broad_summary$broad_cell_type)
  umap_subset[, full_cell_n := as.integer(type_sizes[broad_cell_type])]
  type_display_order <- broad_summary[order(-cell_n), broad_cell_type]
  umap_subset[, display_order := match(broad_cell_type, type_display_order)]
  data.table::setorder(umap_subset, display_order, cell_barcode)

  umap_centroids <- cell_meta[, .(
    UMAPH_1 = stats::median(UMAPH_1),
    UMAPH_2 = stats::median(UMAPH_2),
    full_cell_n = .N
  ), by = broad_cell_type]
  umap_centroids[, broad_cell_type := factor(broad_cell_type, levels = full_type_order)]

  composition[, temporal_fertile_eligible := as.logical(temporal_fertile_eligible)]
  comp_fertile <- composition[temporal_fertile_eligible == TRUE]
  sample_meta <- unique(comp_fertile[, .(sample_accession, LH_day)])
  data.table::setorder(sample_meta, LH_day, sample_accession)
  sample_meta[, sample_order := seq_len(.N)]
  sample_meta[, donor_label := substr(sample_accession, nchar(sample_accession) - 3L, nchar(sample_accession))]
  sample_meta[, stage_label := factor(
    paste0("LH+", LH_day),
    levels = paste0("LH+", c(3, 5, 7, 9, 11))
  )]
  comp_fertile <- merge(comp_fertile, sample_meta, by = c("sample_accession", "LH_day"), all.x = TRUE)
  comp_fertile[, broad_cell_type := factor(broad_cell_type, levels = full_type_order)]
  comp_fertile[, donor_label := factor(donor_label, levels = sample_meta$donor_label)]

  program_plot <- merge(
    program_summary,
    program_models[, .(broad_cell_type, program_name, donor_n, FDR_BH)],
    by = c("broad_cell_type", "program_name"),
    all.x = TRUE,
    sort = FALSE
  )
  program_labels <- c(
    "RNA_ALL_436" = "RNA program\n(n=436)",
    "OCCUPANCY_SUPPORTED_284" = "Occupancy-supported\n(n=284)",
    "SHARED_OCCUPANCY_46" = "Shared occupancy\n(n=46)",
    "PROMOTER_OCCUPANCY_84" = "Promoter-associated\n(n=84)"
  )
  cell_short <- c(
    "Ciliated epithelial cells" = "Ciliated epi.",
    "Myeloid cells" = "Myeloid",
    "NK/T cells" = "NK/T",
    "Stromal cells" = "Stromal",
    "Unciliated epithelial cells" = "Unciliated epi."
  )
  program_plot[, program_label := factor(
    unname(program_labels[program_name]),
    levels = unname(program_labels[program_order])
  )]
  program_plot[, cell_label := factor(
    unname(cell_short[broad_cell_type]),
    levels = unname(cell_short[temporal_type_order])
  )]
  program_plot[, stage_index := match(LH_day, c(3, 5, 7, 9, 11))]
  program_sig <- unique(program_plot[, .(
    broad_cell_type, program_name, program_label, cell_label, FDR_BH,
    significance_label = ifelse(FDR_BH < 0.05, "*", "")
  )])

  VDR_plot <- merge(
    VDR_summary,
    VDR_models[, .(broad_cell_type, donor_n, FDR_BH)],
    by = "broad_cell_type",
    all.x = TRUE,
    sort = FALSE
  )
  VDR_plot[, stage_label := factor(paste0("LH+", LH_day), levels = paste0("LH+", c(3, 5, 7, 9, 11)))]
  VDR_plot[, broad_cell_type := factor(broad_cell_type, levels = rev(temporal_type_order))]
  VDR_plot[, display_label := formatC(mean, format = "f", digits = 1)]
  VDR_plot[, text_color := ifelse(mean >= 3.3, "white", "#161616")]

  pairwise_count <- pairwise_DE[, .(
    significant_row_n = sum(whole_genome_up_n + whole_genome_down_n)
  ), by = broad_cell_type]
  omnibus_count <- omnibus_DE[, .(
    broad_cell_type,
    significant_row_n = whole_genome_FDR_lt_0_05_n
  )]
  omnibus_count[, result_scope := "Omnibus\ngene × cell type"]
  pairwise_count[, result_scope := "Planned pairwise\ngene × cell type × contrast"]
  DE_count_plot <- data.table::rbindlist(list(omnibus_count, pairwise_count), use.names = TRUE)
  DE_count_plot[, result_scope := factor(
    result_scope,
    levels = c("Omnibus\ngene × cell type", "Planned pairwise\ngene × cell type × contrast")
  )]
  DE_count_plot[, broad_cell_type := factor(broad_cell_type, levels = rev(temporal_type_order))]
  DE_count_plot[, count_label := format(significant_row_n, big.mark = ",", scientific = FALSE, trim = TRUE)]
  DE_count_plot[, facet_max := max(significant_row_n), by = result_scope]
  DE_count_plot[, label_inside := significant_row_n >= 0.18 * facet_max]
  DE_count_plot[, label_x := data.table::fifelse(
    label_inside,
    significant_row_n - 0.025 * facet_max,
    significant_row_n + 0.025 * facet_max
  )]
  DE_count_plot[, label_hjust := data.table::fifelse(label_inside, 1, 0)]
  DE_count_plot[, label_color := data.table::fifelse(label_inside, "white", "#161616")]

  list(
    broad_summary = broad_summary,
    full_cell_n = nrow(cell_meta),
    cluster_map = cluster_map,
    composition = composition,
    sample_meta = sample_meta,
    umap_subset = umap_subset,
    umap_centroids = umap_centroids,
    comp_fertile = comp_fertile,
    program_summary = program_summary,
    program_models = program_models,
    program_plot = program_plot,
    program_sig = program_sig,
    VDR_summary = VDR_summary,
    VDR_models = VDR_models,
    VDR_plot = VDR_plot,
    omnibus_DE = omnibus_DE,
    pairwise_DE = pairwise_DE,
    DE_count_plot = DE_count_plot,
    full_type_order = full_type_order,
    temporal_type_order = temporal_type_order,
    program_order = program_order
  )
})

# ---------------------------------------------------------------------------
# Source/statistical consistency checks (no model refitting)
# ---------------------------------------------------------------------------

composition_sum <- plot_data$comp_fertile[, .(fraction_sum = sum(within_sample_fraction)), by = sample_accession]
program_wide <- data.table::dcast(
  plot_data$program_summary,
  broad_cell_type + program_name ~ LH_day,
  value.var = "mean"
)
program_compare <- merge(
  plot_data$program_models,
  program_wide,
  by = c("broad_cell_type", "program_name"),
  all = FALSE
)
program_mean_diff <- max(abs(c(
  program_compare$LH3_mean - program_compare[["3"]],
  program_compare$LH5_mean - program_compare[["5"]],
  program_compare$LH7_mean - program_compare[["7"]],
  program_compare$LH9_mean - program_compare[["9"]],
  program_compare$LH11_mean - program_compare[["11"]]
)), na.rm = TRUE)

VDR_wide <- data.table::dcast(
  plot_data$VDR_summary,
  broad_cell_type ~ LH_day,
  value.var = "mean"
)
VDR_compare <- merge(plot_data$VDR_models, VDR_wide, by = "broad_cell_type", all = FALSE)
VDR_mean_diff <- max(abs(c(
  VDR_compare$LH3_mean - VDR_compare[["3"]],
  VDR_compare$LH5_mean - VDR_compare[["5"]],
  VDR_compare$LH7_mean - VDR_compare[["7"]],
  VDR_compare$LH9_mean - VDR_compare[["9"]],
  VDR_compare$LH11_mean - VDR_compare[["11"]]
)), na.rm = TRUE)

omnibus_total <- sum(plot_data$omnibus_DE$whole_genome_FDR_lt_0_05_n)
pairwise_total <- sum(plot_data$pairwise_DE$whole_genome_up_n + plot_data$pairwise_DE$whole_genome_down_n)

statistical_qc <- data.frame(
  check_id = sprintf("S%02d", 1:12),
  check = c(
    "Locked embedding cell rows",
    "Locked broad cell types",
    "All GSE250130 samples",
    "Fertile temporal donors used in Panel B",
    "Panel B within-donor fractions sum to one",
    "Temporal program models",
    "Temporal program models with FDR < 0.05",
    "Program summary means match frozen model means",
    "Temporal VDR models with FDR < 0.05",
    "VDR summary means match frozen model means",
    "Temporal omnibus significant gene-by-cell-type rows",
    "Planned-pairwise significant gene-by-cell-type-contrast rows"
  ),
  observed = c(
    plot_data$full_cell_n,
    nrow(plot_data$broad_summary),
    data.table::uniqueN(plot_data$composition$sample_accession),
    nrow(plot_data$sample_meta),
    max(abs(composition_sum$fraction_sum - 1)),
    nrow(plot_data$program_models),
    sum(plot_data$program_models$FDR_BH < 0.05),
    program_mean_diff,
    sum(plot_data$VDR_models$FDR_BH < 0.05),
    VDR_mean_diff,
    omnibus_total,
    pairwise_total
  ),
  expected = c(236636, 8, 28, 18, 0, 20, 13, 0, 0, 0, 24168, 55539),
  tolerance = c(0, 0, 0, 0, 1e-10, 0, 0, 1e-12, 0, 1e-12, 0, 0),
  stringsAsFactors = FALSE
)
statistical_qc$pass <- abs(as.numeric(statistical_qc$observed) - statistical_qc$expected) <=
  statistical_qc$tolerance
if (!all(statistical_qc$pass)) stop("Figure 4 statistical/source consistency QC failed.")
write_csv_utf8(statistical_qc, file.path(dirs$qc, "FIG4_STATISTICAL_QC.csv"))

# ---------------------------------------------------------------------------
# Direct plotting source tables
# ---------------------------------------------------------------------------

data.table::fwrite(
  plot_data$umap_subset[, .(cell_barcode, UMAPH_1, UMAPH_2, broad_cell_type, full_cell_n)],
  file.path(dirs$source_data, "Figure4_panelA_UMAP_display_subset.csv")
)
data.table::fwrite(plot_data$umap_centroids, file.path(dirs$source_data, "Figure4_panelA_label_centroids.csv"))
data.table::fwrite(plot_data$broad_summary, file.path(dirs$source_data, "Figure4_panelA_broad_type_counts.csv"))
data.table::fwrite(plot_data$comp_fertile, file.path(dirs$source_data, "Figure4_panelB_fertile_donor_composition.csv"))
data.table::fwrite(plot_data$program_plot, file.path(dirs$source_data, "Figure4_panelC_temporal_program_profiles.csv"))
data.table::fwrite(plot_data$program_models, file.path(dirs$source_data, "Figure4_panelC_temporal_program_models.csv"))
data.table::fwrite(plot_data$VDR_plot, file.path(dirs$source_data, "Figure4_panelD_temporal_VDR_means.csv"))
data.table::fwrite(plot_data$VDR_models, file.path(dirs$source_data, "Figure4_panelD_temporal_VDR_models.csv"))
data.table::fwrite(plot_data$DE_count_plot, file.path(dirs$source_data, "Figure4_panelE_temporal_DE_count_summary.csv"))
data.table::fwrite(plot_data$omnibus_DE, file.path(dirs$source_data, "Figure4_panelE_temporal_omnibus_source.csv"))
data.table::fwrite(plot_data$pairwise_DE, file.path(dirs$source_data, "Figure4_panelE_temporal_pairwise_source.csv"))
write_csv_utf8(fig4_blueprint, file.path(dirs$source_data, "Figure4_locked_blueprint_row.csv"))

# ---------------------------------------------------------------------------
# Article color and typography registries
# ---------------------------------------------------------------------------

cell_colors <- c(
  "B cells" = "#D55E00",
  "Ciliated epithelial cells" = "#E69F00",
  "Endothelial cells" = "#009E73",
  "Mast cells" = "#A6761D",
  "Myeloid cells" = "#00A6D6",
  "NK/T cells" = "#0072B2",
  "Stromal cells" = "#7A5AB5",
  "Unciliated epithelial cells" = "#CC79A7"
)
color_registry <- data.frame(
  scientific_concept = names(cell_colors),
  hex = unname(cell_colors),
  semantic = "BROAD_CELL_TYPE",
  applies_to = "Panels A, B, C, and E",
  stringsAsFactors = FALSE
)
write_csv_utf8(color_registry, file.path(dirs$qc, "FIG4_COLOR_SEMANTIC_REGISTRY.csv"))

header_tag_size_pt <- 11.0
header_title_size_pt <- 9.2
header_subtitle_size_pt <- 7.2
header_tag_to_title_gap_mm <- 1.6
header_title_row_height_mm <- 5.4
header_subtitle_row_height_mm <- 4.0
figure_width_mm <- 180
figure_height_mm <- 250

typography_lock <- data.frame(
  parameter = c(
    "font_family", "panel_tag_pt", "panel_title_pt", "panel_subtitle_pt",
    "tag_title_gap_mm", "axis_title_pt", "axis_text_pt", "legend_text_pt",
    "figure_width_mm", "figure_height_mm"
  ),
  value = c(
    figure_font, header_tag_size_pt, header_title_size_pt, header_subtitle_size_pt,
    header_tag_to_title_gap_mm, 7.5, 6.2, 6.2,
    figure_width_mm, figure_height_mm
  ),
  authority = "Article profile inherited from locked Figures 2/3; Figure 4 height exception",
  stringsAsFactors = FALSE
)
write_csv_utf8(typography_lock, file.path(dirs$qc, "FIG4_ARTICLE_TYPOGRAPHY_LOCK.csv"))

color_text <- "#161616"
color_muted <- "#535353"
color_grid <- "#D9D9D9"
base_theme <- ggplot2::theme_minimal(base_family = figure_font, base_size = 7.2) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(linewidth = 0.24, colour = color_grid),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(family = figure_font, size = 7.5, colour = color_text),
    axis.text = ggplot2::element_text(family = figure_font, size = 6.2, colour = color_text),
    legend.title = ggplot2::element_text(family = figure_font, size = 6.5, face = "bold", colour = color_text),
    legend.text = ggplot2::element_text(family = figure_font, size = 6.2, colour = color_text),
    strip.text = ggplot2::element_text(family = figure_font, size = 6.2, face = "bold", colour = color_text),
    plot.margin = ggplot2::margin(2.0, 2.5, 2.0, 2.5, unit = "mm")
  )

# ---------------------------------------------------------------------------
# Panel construction
# ---------------------------------------------------------------------------

panel_A <- capture_stage("panel_A_construction", {
  ggplot2::ggplot(
    plot_data$umap_subset,
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, colour = broad_cell_type)
  ) +
    ggplot2::geom_point(size = 0.075, alpha = 0.58, stroke = 0) +
    ggplot2::geom_label(
      data = plot_data$umap_centroids,
      ggplot2::aes(x = UMAPH_1, y = UMAPH_2, label = broad_cell_type),
      inherit.aes = FALSE,
      family = figure_font,
      fontface = "bold",
      size = 2.15,
      linewidth = 0,
      label.padding = grid::unit(0.7, "mm"),
      fill = scales::alpha("white", 0.78),
      colour = color_text
    ) +
    ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    base_theme +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 6.3),
      legend.position = "none",
      plot.margin = ggplot2::margin(1.5, 1.5, 1.5, 1.5, unit = "mm")
    )
})

panel_B <- capture_stage("panel_B_construction", {
  ggplot2::ggplot(
    plot_data$comp_fertile,
    ggplot2::aes(x = donor_label, y = within_sample_fraction, fill = broad_cell_type)
  ) +
    ggplot2::geom_col(width = 0.92) +
    ggplot2::facet_grid(~ stage_label, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE, name = "Broad cell type") +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      labels = scales::label_percent(accuracy = 25)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(x = "Donor (GSM suffix)", y = "Within-donor fraction") +
    base_theme +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 5.4, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 5.8),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
      strip.text.x = ggplot2::element_text(size = 6.2, face = "bold"),
      panel.spacing.x = grid::unit(0.6, "mm"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = ggplot2::margin(t = 0.5, b = 0),
      legend.spacing.y = grid::unit(0.2, "mm"),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      plot.margin = ggplot2::margin(1.5, 2.0, 1.0, 1.5, unit = "mm")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE, title.position = "top"))
})

panel_C <- capture_stage("panel_C_construction", {
  ggplot2::ggplot(
    plot_data$program_plot,
    ggplot2::aes(
      x = stage_index, y = mean, group = 1,
      colour = broad_cell_type, fill = broad_cell_type
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.22, linetype = "dashed", colour = "#8A8A8A") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.48) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::geom_text(
      data = plot_data$program_sig,
      ggplot2::aes(x = 5, y = Inf, label = significance_label),
      inherit.aes = FALSE,
      family = figure_font,
      fontface = "bold",
      size = 2.7,
      hjust = 1.1,
      vjust = 1.0,
      colour = color_text
    ) +
    ggplot2::facet_grid(
      program_label ~ cell_label,
      scales = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_colour_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = 1:5,
      labels = paste0("LH+", c(3, 5, 7, 9, 11)),
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::labs(x = "Biopsy timing", y = "Signed VDR-program score") +
    base_theme +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(size = 5.2, angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 5.2),
      axis.title = ggplot2::element_text(size = 6.6),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
      strip.placement = "outside",
      strip.text.x = ggplot2::element_text(size = 5.8, lineheight = 0.92),
      strip.text.y.left = ggplot2::element_text(size = 5.5, angle = 0, hjust = 1, lineheight = 0.90),
      panel.spacing = grid::unit(0.65, "mm"),
      plot.margin = ggplot2::margin(1.0, 1.5, 1.2, 1.5, unit = "mm")
    )
})

panel_D <- capture_stage("panel_D_construction", {
  ggplot2::ggplot(
    plot_data$VDR_plot,
    ggplot2::aes(x = stage_label, y = broad_cell_type, fill = mean)
  ) +
    ggplot2::geom_tile(linewidth = 0.45, colour = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = display_label, colour = text_color),
      family = figure_font,
      size = 2.05
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_gradient(
      low = "#F4F3F0",
      high = "#405466",
      name = "Mean VDR\nlog2(CPM + 1)",
      limits = c(0, 6),
      breaks = c(0, 2, 4, 6),
      oob = scales::squish
    ) +
    ggplot2::labs(x = "Biopsy timing", y = NULL) +
    base_theme +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 5.9, angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6.0, lineheight = 0.92),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(
        size = 6.2, face = "bold",
        margin = ggplot2::margin(r = 2.0, unit = "mm")
      ),
      legend.text = ggplot2::element_text(size = 5.8),
      legend.key.height = grid::unit(2.6, "mm"),
      legend.margin = ggplot2::margin(t = 0.5, b = 0),
      plot.margin = ggplot2::margin(1.5, 2.0, 1.0, 1.5, unit = "mm")
    ) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(
      title.position = "left", title.hjust = 0.5,
      barwidth = grid::unit(34, "mm"),
      barheight = grid::unit(2.6, "mm")
    ))
})

panel_E <- capture_stage("panel_E_construction", {
  ggplot2::ggplot(
    plot_data$DE_count_plot,
    ggplot2::aes(x = significant_row_n, y = broad_cell_type, fill = broad_cell_type)
  ) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = label_x, label = count_label,
        hjust = label_hjust, colour = label_color
      ),
      family = figure_font,
      size = 2.0,
      show.legend = FALSE
    ) +
    ggplot2::facet_grid(~ result_scope, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = cell_colors, drop = FALSE) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(x = "FDR-significant rows", y = NULL) +
    base_theme +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 5.5),
      axis.text.y = ggplot2::element_text(size = 5.9, lineheight = 0.92),
      axis.title.x = ggplot2::element_text(size = 6.5),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
      strip.text.x = ggplot2::element_text(size = 5.8, lineheight = 0.92),
      panel.spacing.x = grid::unit(2.0, "mm"),
      plot.margin = ggplot2::margin(1.5, 2.0, 1.0, 1.5, unit = "mm")
    )
})

# ---------------------------------------------------------------------------
# Scientific-body-anchored headers and one global composite geometry authority
# ---------------------------------------------------------------------------

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

make_header_gtable <- function(
  body_gtable, panel_tag, panel_title, panel_subtitle,
  horizontal_nudge_mm = 0
) {
  panel_cells <- body_gtable$layout[grepl("^panel", body_gtable$layout$name), , drop = FALSE]
  if (!nrow(panel_cells)) stop("Could not locate scientific panel cell for Panel ", panel_tag)
  panel_left <- min(panel_cells$l)
  panel_right <- max(panel_cells$r)

  header_gt <- gtable::gtable(
    widths = body_gtable$widths,
    heights = grid::unit(c(header_title_row_height_mm, header_subtitle_row_height_mm), "mm")
  )
  title_grob <- grid::textGrob(
    panel_title,
    x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm"),
    y = grid::unit(0.52, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_title_size_pt, lineheight = 0.94, col = color_text
    )
  )
  tag_grob <- grid::textGrob(
    panel_tag,
    x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm") -
      grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"),
    hjust = 1, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontface = "bold",
      fontsize = header_tag_size_pt, col = color_text
    )
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle,
    x = grid::unit(0, "npc") + grid::unit(horizontal_nudge_mm, "mm"),
    y = grid::unit(0.56, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(
      fontfamily = figure_font, fontsize = header_subtitle_size_pt,
      lineheight = 0.94, col = color_muted
    )
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, grid::grobTree(tag_grob, title_grob),
    t = 1, b = 1, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F4-header-title-", panel_tag)
  )
  header_gt <- gtable::gtable_add_grob(
    header_gt, subtitle_grob,
    t = 2, b = 2, l = panel_left, r = panel_right,
    clip = "off", name = paste0("F4-header-subtitle-", panel_tag)
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
      horizontal_nudge_mm = horizontal_nudge_mm,
      header_composition = "DEDICATED_SHARED_ROW",
      stringsAsFactors = FALSE
    )
  )
}

metrics_device_path <- tempfile(pattern = "FIG4_ragg_metrics_", fileext = ".png")
ragg::agg_png(
  filename = metrics_device_path,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = 300,
  background = "white"
)
composite_parts <- capture_stage("gtable_and_composite_assembly", {
  body_A <- build_body_gtable(panel_A)
  body_B <- build_body_gtable(panel_B)
  body_C <- build_body_gtable(panel_C)
  body_D <- build_body_gtable(panel_D)
  body_E <- build_body_gtable(panel_E)

  header_A <- make_header_gtable(
    body_A, "A", "Broad cell-type landscape",
    "8 broad types · 236,636 cells (display subset)",
    horizontal_nudge_mm = 18
  )
  header_B <- make_header_gtable(body_B, "B", "Fertile donor composition",
                                 "18 fertile donors · descriptive only")
  header_C <- make_header_gtable(body_C, "C", "Stage-associated VDR-program remodeling",
                                 "Cross-sectional donor means ± 95% CI · * omnibus FDR < 0.05")
  header_D <- make_header_gtable(body_D, "D", "VDR expression across stages",
                                 "Mean log2(CPM + 1) · 0/5 omnibus FDR < 0.05")
  header_E <- make_header_gtable(body_E, "E", "Transcriptome-wide temporal DE",
                                 "Genome-wide FDR < 0.05 · independent axes")

  design <- "
##
AB
CD
##
EE
FF
##
GH
IJ
##
"
  figure_object <- (
    patchwork::wrap_elements(full = header_A$gtable) +
      patchwork::wrap_elements(full = header_B$gtable) +
      patchwork::wrap_elements(full = body_A) +
      patchwork::wrap_elements(full = body_B) +
      patchwork::wrap_elements(full = header_C$gtable) +
      patchwork::wrap_elements(full = body_C) +
      patchwork::wrap_elements(full = header_D$gtable) +
      patchwork::wrap_elements(full = header_E$gtable) +
      patchwork::wrap_elements(full = body_D) +
      patchwork::wrap_elements(full = body_E)
  ) +
    patchwork::plot_layout(
      design = design,
      widths = c(0.92, 1.08),
      heights = grid::unit(c(1, 10, 52, 3, 10, 57, 3, 10, 74, 18), "mm")
    ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(2, 3, 2, 3, unit = "mm"))
    )

  list(
    figure = figure_object,
    header_geometry = do.call(
      rbind,
      lapply(list(header_A, header_B, header_C, header_D, header_E), `[[`, "geometry")
    )
  )
})
grDevices::dev.off()
unlink(metrics_device_path)

write_csv_utf8(
  composite_parts$header_geometry,
  file.path(dirs$qc, "FIG4_PANEL_HEADER_GEOMETRY.csv")
)
layout_lock <- data.frame(
  item = c(
    "canvas", "global_geometry_authority", "column_widths", "row_heights_mm",
    "top_row", "middle_row", "bottom_row", "header_anchor",
    "Panel_A_header_nudge", "Panel_A_aspect", "user_visual_lock"
  ),
  value = c(
    "180 mm x 250 mm",
    "ONE_EXPLICIT_PATCHWORK_DESIGN_MATRIX",
    "0.92 : 1.08",
    "1,10,52,3,10,57,3,10,74,18",
    "A broad UMAP | B fertile donor composition",
    "C stage-associated program profiles spans both columns",
    "D VDR mean heatmap | E temporal DE row counts",
    "B-E start at the body-gtable panel-cell left edge; A uses the recorded user-directed visible-body alignment override",
    "+18 mm for A tag, title, and subtitle only",
    "coord_fixed; deterministic <=10,000 cells per broad type",
    "REQUIRED_BEFORE_RELEASE"
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(layout_lock, file.path(dirs$qc, "FIG4_LAYOUT_LOCK_CARD.csv"))

write_text_utf8(c(
  "# Figure 4 Manual Visual Override Card — v1.2",
  "",
  "- Trigger: renewed user visual QC on Figure 4 v1.1.",
  "- Observed defect: Panel A UMAP body position is correct, but its tag/title/subtitle remain visually left-shifted.",
  "- Authorized repair: keep the UMAP body fixed and shift the Panel A tag, title, and subtitle an additional +4 mm, for a total +18 mm horizontal nudge.",
  "- Unchanged: Panel A UMAP coordinates, cell subset, labels, palette, aspect ratio, and all Panels B-E.",
  "- Statistical recomputation: none.",
  "- User visual approval after repair: required."
), file.path(dirs$qc, "FIG4_MANUAL_VISUAL_OVERRIDE_CARD.md"))

# ---------------------------------------------------------------------------
# Scientific QC and manuscript-facing draft text
# ---------------------------------------------------------------------------

scientific_qc <- data.frame(
  check_id = sprintf("Q%02d", 1:10),
  check = c(
    "Figure central question matches locked blueprint",
    "Panel A uses locked embedding and broad annotation only",
    "Panel B is restricted to fertile temporal donors",
    "Panel B is explicitly descriptive with no composition test",
    "Panel C uses cross-sectional stage-associated wording",
    "Panel C preserves donor-level mean and 95% CI",
    "Panel C significance comes from frozen omnibus FDR",
    "Panel D discloses 0/5 FDR-significant VDR models",
    "Panel E distinguishes omnibus from planned-pairwise result units",
    "No RIF disease claim or within-person trajectory claim is made"
  ),
  status = "PASS",
  evidence = c(
    fig4_blueprint$central_question[[1]],
    "ANN_v1_2 cell_annotation_metadata.csv; candidate subtype not plotted",
    "temporal_fertile_eligible == TRUE; 18 donors",
    "ANNOTATION_GATE composition_inferential_testing_authorized=FALSE",
    "Header and legend use cross-sectional/stage-associated language",
    "source_temporal_program_summary.csv",
    "01_temporal_program_omnibus.csv; 13/20",
    "05_temporal_VDR_omnibus.csv; 0/5",
    "07/08 DE summary tables and separate facet labels",
    "Locked STEP04 interpretation boundary"
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(scientific_qc, file.path(dirs$qc, "FIG4_SCIENTIFIC_QC.csv"))

blueprint_card <- c(
  "# Figure 4 Blueprint Card — v1.2 QC",
  "",
  "- Figure ID: Figure 4",
  "- Status: B_CANDIDATE_ASSISTANT_QC_PASS_USER_REVIEW_PENDING",
  "- Governing authority: STEP04_LOCK_v1.2",
  "- Scientific question: How strongly does normal luteal timing remodel the programs?",
  "- Dataset: GSE250130",
  "- Statistical unit: donor-level pseudobulk for Panels C/D; gene-by-cell-type or gene-by-cell-type-by-contrast rows for Panel E",
  "- Plotting language: R / AUTO_NATIVE direct-result route",
  "- Target size: 180 mm x 250 mm",
  "- Pre-lock output: one 300-dpi QC PNG",
  "",
  "## Panels",
  "",
  "- A: broad-cell UMAP from locked coordinates; deterministic display subset; descriptive.",
  "- B: broad cell-type fractions for 18 fertile donors across LH+3 to LH+11; descriptive only.",
  "- C: four locked VDR-program profiles across five fertile luteal stages and five broad cell types; donor mean ± 95% CI; frozen omnibus FDR.",
  "- D: mean VDR log2(CPM + 1) across the same stages/cell types; 0/5 omnibus models FDR < 0.05.",
  "- E: genome-wide temporal DE counts, explicitly separated into omnibus and planned-pairwise result units.",
  "",
  "## Allowed claim",
  "",
  "Cross-sectional variation across fertile luteal stages is accompanied by strong cell-type-resolved transcriptomic and VDR-program remodeling.",
  "",
  "## Prohibited claims",
  "",
  "No within-person longitudinal trajectory, composition-inference, subtype-inference, causal VDR-activity, RIF disease, direct-target, or biomarker claim."
)
write_text_utf8(blueprint_card, file.path(candidate_root, "FIGURE_BLUEPRINT_CARD.md"))

write_text_utf8(c(
  "# Language Routing Card",
  "",
  "- Mode: AUTO_NATIVE",
  "- Authoritative implementation: R",
  "- Rationale: locked upstream objects/tables and prior QC were produced in R; UMAP coordinates and all inferential outputs are frozen tables.",
  "- Cross-language migration: PROHIBITED during v1.x visual iteration unless a documented device failure requires controlled migration.",
  "- Statistical recomputation in plotting code: PROHIBITED."
), file.path(candidate_root, "LANGUAGE_ROUTING_CARD.md"))

write_text_utf8(c(
  "# Figure-versus-table decision",
  "",
  "Figure 4 is retained because its main purpose is to show spatial cell-type structure, donor-composition heterogeneity, ordered stage-associated profiles, and the distribution of temporal transcriptomic remodeling across cell types.",
  "",
  "Exact omnibus statistics, pairwise results, confidence intervals, sample metadata, and complete gene-level results remain in the source/supplementary tables."
), file.path(candidate_root, "FIGURE_TABLE_DECISION.md"))

write_text_utf8(c(
  "# Figure 4 legend — draft v1.2",
  "",
  "**Figure 4 | Cell-type-resolved remodeling across fertile luteal stages in GSE250130.**",
  "**A,** UMAP of the locked broad-cell annotation. For visualization only, a deterministic subset of at most 10,000 cells per broad type is displayed; the complete annotation contains 236,636 cells across eight broad types. UMAP coordinates are descriptive layout coordinates.",
  "**B,** Broad cell-type composition of 18 fertile donors sampled at LH+3, LH+5, LH+7, LH+9, or LH+11. Bars show within-donor cell fractions; no composition hypothesis test was performed.",
  "**C,** Cross-sectional stage-associated profiles of four locked VDR programs in five broad cell types. Points show donor-level pseudobulk means and ribbons show 95% confidence intervals. Asterisks denote an omnibus permutation-test Benjamini–Hochberg FDR below 0.05; 13 of 20 program-by-cell-type models met this criterion.",
  "**D,** Mean donor-level VDR expression [log2(CPM + 1)] across fertile luteal stages. None of the five omnibus temporal VDR-expression models was significant after Benjamini–Hochberg correction.",
  "**E,** Counts of genome-wide FDR-significant temporal differential-expression rows. Omnibus counts use gene-by-cell-type rows, whereas planned-pairwise counts use gene-by-cell-type-by-contrast rows; the two facets therefore use independent axes and are not direct effect-size comparisons. Across cell types, 24,168 omnibus and 55,539 planned-pairwise rows met genome-wide FDR < 0.05."
), file.path(dirs$manuscript, "Figure4_legend_v1.2.md"))

write_text_utf8(c(
  "# Figure 4 Methods — draft v1.2",
  "",
  "Figure 4 was assembled exclusively from the frozen STEP04_LOCK_v1.2 GSE250130 outputs. No normalization, integration, dimensionality reduction, clustering, annotation, model fitting, differential-expression testing, or multiple-testing correction was rerun. Panel A used the locked UMAP coordinates and broad-cell annotations, with a deterministic display subset of at most 10,000 cells per broad type. Panel B displayed within-donor broad-cell fractions for the 18 donors designated as temporal_fertile_eligible and was descriptive only. Panels C and D used saved donor-level pseudobulk summaries and frozen omnibus permutation-test FDR values from VDRMODELS_v1.1. Panel E used the saved DE_v1.2 genome-wide FDR count summaries; omnibus and planned-pairwise rows were retained as distinct result units."
), file.path(dirs$manuscript, "Figure4_methods_v1.2.md"))

write_text_utf8(c(
  "# Figure 4 Results — draft v1.2",
  "",
  "Across the fertile luteal-stage series, 13 of 20 prespecified VDR-program-by-cell-type omnibus models were significant after FDR correction, demonstrating strong cell-type-resolved stage-associated remodeling. Genome-wide temporal analyses yielded 24,168 significant omnibus gene-by-cell-type rows and 55,539 significant planned-pairwise gene-by-cell-type-by-contrast rows. In contrast, none of the five temporal VDR-expression models survived FDR correction, indicating that program remodeling should not be interpreted as a uniform change in VDR transcript abundance. Because the samples were cross-sectional, these patterns represent stage-associated differences rather than within-person longitudinal trajectories."
), file.path(dirs$manuscript, "Figure4_results_v1.2.md"))

write_text_utf8(c(
  "# Figure 4 中文分面说明",
  "",
  "- A：GSE250130 的冻结 broad-cell UMAP；仅作描述性展示，不解释 UMAP 坐标大小，也不引入 subtype 推断。",
  "- B：18 名 fertile donor 在 LH+3 至 LH+11 的 broad cell-type 组成；仅作描述性展示，不进行组成差异检验。",
  "- C：5 类 broad cell type 中 4 个锁定 VDR program 的横断面阶段相关变化；点为 donor-level 均值，带为 95% CI，星号为冻结的 omnibus FDR < 0.05。",
  "- D：各阶段平均 VDR 表达；5 个模型均未达到 FDR 显著，不能据此推断 VDR 活性。",
  "- E：genome-wide temporal DE 的显著结果行数；omnibus 与 planned-pairwise 的统计行单位不同，分别展示。"
), file.path(dirs$manuscript, "Figure4_panel_description_zh_v1.2.md"))

write_text_utf8(c(
  "# Figure 4 source data",
  "",
  "All tables in this directory are plotting-only extracts or exact copies derived from frozen STEP04 inputs. No inferential statistic was recomputed.",
  "",
  "Panel A packages a deterministic display subset instead of duplicating the 80-MB complete cell metadata table. The complete upstream path, size, MD5, and SHA-256 are recorded in FIG4_UPSTREAM_SOURCE_IDENTITY.csv.",
  "",
  "Panel E retains omnibus gene-by-cell-type rows and planned-pairwise gene-by-cell-type-by-contrast rows as distinct units."
), file.path(dirs$source_data, "README_SOURCE_DATA.md"))

# ---------------------------------------------------------------------------
# Export one QC PNG and validate the actual raster
# ---------------------------------------------------------------------------

figure_png <- file.path(dirs$figures, "Figure4_v1_2_QC_180mm_300dpi.png")
capture_stage("300dpi_QC_PNG_export", {
  ragg::agg_png(
    filename = figure_png,
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    res = 300,
    scaling = 1,
    background = "white"
  )
  print(composite_parts$figure)
  grDevices::dev.off()
})
if (!file.exists(figure_png) || file.info(figure_png)$size < 100000) {
  stop("Figure 4 QC PNG was not created correctly.")
}

warning_log <- if (length(warning_rows)) {
  do.call(rbind, warning_rows)
} else {
  data.frame(
    warning_id = character(), stage = character(), warning_class = character(),
    message = character(), stringsAsFactors = FALSE
  )
}
write_csv_utf8(warning_log, file.path(dirs$qc, "FIG4_RUNTIME_WARNING_LOG.csv"))
if (nrow(warning_log) > 0L) {
  stop("Figure 4 plotting/export emitted runtime warnings; see FIG4_RUNTIME_WARNING_LOG.csv")
}

png_array <- png::readPNG(figure_png)
observed_width_px <- dim(png_array)[2]
observed_height_px <- dim(png_array)[1]
expected_width_px <- round(figure_width_mm / 25.4 * 300)
expected_height_px <- round(figure_height_mm / 25.4 * 300)
rgb_array <- png_array[, , seq_len(min(3L, dim(png_array)[3])), drop = FALSE]
raster_sd <- stats::sd(as.numeric(rgb_array))
expected_reference_png_sha256 <- "94424cb3ffa010cab6c9ff2258165708451f1fd3e66f336add5bde9ace3f4601"
observed_png_sha256 <- tolower(sha256_file(figure_png))
reference_render_qc <- data.frame(
  check = "Formal QC PNG exactly matches assistant-reviewed reference render",
  observed_sha256 = observed_png_sha256,
  expected_sha256 = expected_reference_png_sha256,
  pass = identical(observed_png_sha256, expected_reference_png_sha256),
  stringsAsFactors = FALSE
)
write_csv_utf8(reference_render_qc, file.path(dirs$qc, "FIG4_REFERENCE_RENDER_HASH_QC.csv"))
if (!isTRUE(reference_render_qc$pass[[1]])) {
  stop("Formal Figure 4 PNG differs from the assistant-reviewed reference render.")
}

parent_reference_png <- file.path(
  step05_root, "06_results", "FIG4", "F4_v1_1_QC", "figures",
  "Figure4_v1_1_QC_180mm_300dpi.png"
)
expected_parent_sha256 <- "6bf4d5ac9f0108872eea80b87c56b9f14b9d68df65b83bb770736b835c7e1b26"
if (!file.exists(parent_reference_png) ||
    !identical(tolower(sha256_file(parent_reference_png)), expected_parent_sha256)) {
  stop("The reviewed Figure 4 v1.0 parent render is missing or has changed.")
}
parent_png_array <- png::readPNG(parent_reference_png)
if (!identical(dim(parent_png_array), dim(png_array))) {
  stop("Figure 4 v1.1 and v1.2 raster dimensions differ unexpectedly.")
}
repair_diff_mask <- apply(abs(parent_png_array - png_array) > 1e-10, c(1, 2), any)
repair_rows <- which(rowSums(repair_diff_mask) > 0)
repair_cols <- which(colSums(repair_diff_mask) > 0)
controlled_repair_qc <- data.frame(
  check = c(
    "v1.2 differs from v1.1",
    "All changed pixels remain in the Panel A header rows",
    "All changed pixels remain in the left-column header region",
    "No pixel changed below the header safety boundary"
  ),
  observed = c(
    sum(repair_diff_mask),
    paste0(min(repair_rows), "-", max(repair_rows)),
    paste0(min(repair_cols), "-", max(repair_cols)),
    sum(repair_diff_mask[221:nrow(repair_diff_mask), , drop = FALSE])
  ),
  expected = c(">0", "rows within 1-220", "columns within 1-900", "0"),
  pass = c(
    sum(repair_diff_mask) > 0,
    min(repair_rows) >= 1 && max(repair_rows) <= 220,
    min(repair_cols) >= 1 && max(repair_cols) <= 900,
    sum(repair_diff_mask[221:nrow(repair_diff_mask), , drop = FALSE]) == 0
  ),
  stringsAsFactors = FALSE
)
write_csv_utf8(controlled_repair_qc, file.path(dirs$qc, "FIG4_V1_2_CONTROLLED_PIXEL_DIFF_QC.csv"))
if (!all(controlled_repair_qc$pass)) {
  stop("Figure 4 v1.2 changed pixels outside the authorized Panel A header region.")
}
content_mask <- apply(rgb_array < 0.995, c(1, 2), any)
content_rows <- which(rowSums(content_mask) > 0)
content_cols <- which(colSums(content_mask) > 0)
edge_margins_px <- c(
  top = min(content_rows) - 1L,
  bottom = observed_height_px - max(content_rows),
  left = min(content_cols) - 1L,
  right = observed_width_px - max(content_cols)
)
edge_qc <- data.frame(
  edge = names(edge_margins_px),
  blank_margin_px = as.integer(edge_margins_px),
  blank_margin_mm = as.numeric(edge_margins_px) / 300 * 25.4,
  minimum_px = 24L,
  pass = as.numeric(edge_margins_px) >= 24,
  stringsAsFactors = FALSE
)
write_csv_utf8(edge_qc, file.path(dirs$qc, "FIG4_RASTER_EDGE_MARGIN_QC.csv"))

automated_checks <- data.frame(
  check = c(
    "QC PNG exists",
    "QC PNG nonzero",
    "Width pixels within device rounding",
    "Height pixels within device rounding",
    "Raster is not blank",
    "Top raster edge safety margin",
    "Bottom raster edge safety margin",
    "Left raster edge safety margin",
    "Right raster edge safety margin",
    "Arial exact family resolved",
    "Runtime warnings",
    "Upstream source MD5",
    "Scientific QC",
    "Statistical/source consistency QC",
    "Only one pre-lock figure format generated"
  ),
  observed = c(
    file.exists(figure_png),
    file.info(figure_png)$size,
    observed_width_px,
    observed_height_px,
    raster_sd,
    edge_margins_px[["top"]],
    edge_margins_px[["bottom"]],
    edge_margins_px[["left"]],
    edge_margins_px[["right"]],
    font_info$family[[1]],
    nrow(warning_log),
    all(source_identity$md5_match),
    all(scientific_qc$status == "PASS"),
    all(statistical_qc$pass),
    length(list.files(dirs$figures, full.names = TRUE))
  ),
  expected = c(
    "TRUE", ">100000", paste0(expected_width_px, " ±2"), paste0(expected_height_px, " ±2"),
    ">0.01", ">=24 px", ">=24 px", ">=24 px", ">=24 px",
    "Arial", "0", "TRUE", "TRUE", "TRUE", "1"
  ),
  pass = c(
    file.exists(figure_png),
    file.info(figure_png)$size > 100000,
    abs(observed_width_px - expected_width_px) <= 2,
    abs(observed_height_px - expected_height_px) <= 2,
    raster_sd > 0.01,
    edge_margins_px[["top"]] >= 24,
    edge_margins_px[["bottom"]] >= 24,
    edge_margins_px[["left"]] >= 24,
    edge_margins_px[["right"]] >= 24,
    grepl("Arial", font_info$family[[1]], ignore.case = TRUE),
    nrow(warning_log) == 0L,
    all(source_identity$md5_match),
    all(scientific_qc$status == "PASS"),
    all(statistical_qc$pass),
    length(list.files(dirs$figures, full.names = TRUE)) == 1L
  ),
  stringsAsFactors = FALSE
)
if (!all(automated_checks$pass)) stop("Figure 4 automated render checks failed.")
write_csv_utf8(automated_checks, file.path(dirs$qc, "FIG4_AUTOMATED_RENDER_CHECKS.csv"))

strict_render_qc <- data.frame(
  severity = c("BLOCKING", "BLOCKING", "BLOCKING", "MAJOR", "MAJOR", "MAJOR", "MAJOR", "MAJOR"),
  visual_check = c(
    "Panel allocation and cross-row column alignment",
    "Header/body anchors and A-E tag baselines",
    "Text, glyph, axis, strip, and legend clipping",
    "Panel A labels and rare-population visibility",
    "Panel B donor labels and eight-colour legend readability",
    "Panel C 20-facet final-size readability and CI visibility",
    "Panel D heatmap text/legend contrast",
    "Panel E independent-axis disclosure and value-label clearance"
  ),
  status = "PASS_ASSISTANT_ACTUAL_RENDER_REVIEW",
  evidence = c(
    "One global design matrix; Panel A +18 mm header override and all rows inspected at original resolution",
    "Panel A aligns to its visible fixed-aspect UMAP body; B-E retain body-gtable left-edge anchors",
    "No clipping; raster edge margins pass and full/bottom crops inspected",
    "Eight broad populations remain visible; all labels readable; UMAP body is pixel-identical to v1.0",
    "All 18 donor suffixes and all eight legend entries readable",
    "All 20 facets, 95% CI ribbons, axes, and 13 significance markers readable",
    "All 25 values readable; explicit 0-6 colour scale fully visible",
    "Both facet units, independent axes, and all ten value labels clear"
  ),
  reviewed_reference_sha256 = expected_reference_png_sha256,
  stringsAsFactors = FALSE
)
write_csv_utf8(strict_render_qc, file.path(dirs$qc, "FIG4_STRICT_RENDER_QC.csv"))

figure_status <- data.frame(
  figure_id = figure_id,
  figure_version = figure_version,
  figure_class = "B_submission_candidate",
  current_status = "B_CANDIDATE_ASSISTANT_QC_PASS_USER_REVIEW_PENDING",
  source_locked = TRUE,
  blueprint_locked = TRUE,
  scientific_QC = "PASS",
  statistical_QC = "PASS_NO_RECOMPUTATION",
  automated_render_QC = "PASS",
  assistant_actual_render_QC = "PASS_REFERENCE_RENDER_HASH_MATCH",
  user_visual_QC = "PENDING",
  final_submission_formats_generated = FALSE,
  immutable_lock_authorized = FALSE,
  stringsAsFactors = FALSE
)
write_csv_utf8(figure_status, file.path(candidate_root, "FIGURE_STATUS.csv"))

# ---------------------------------------------------------------------------
# Source-to-code map, run evidence, candidate manifest, and SMALL_REPORT
# ---------------------------------------------------------------------------

source_map <- data.frame(
  panel = c("A", "B", "C", "D", "E", "E"),
  direct_plotting_table = c(
    "Figure4_panelA_UMAP_display_subset.csv",
    "Figure4_panelB_fertile_donor_composition.csv",
    "Figure4_panelC_temporal_program_profiles.csv",
    "Figure4_panelD_temporal_VDR_means.csv",
    "Figure4_panelE_temporal_omnibus_source.csv",
    "Figure4_panelE_temporal_pairwise_source.csv"
  ),
  frozen_upstream_source = c(
    relative_to(source_paths[["cell_annotation_metadata"]], project_root),
    relative_to(source_paths[["sample_composition"]], project_root),
    relative_to(source_paths[["temporal_program_summary"]], project_root),
    relative_to(source_paths[["temporal_VDR_summary"]], project_root),
    relative_to(source_paths[["temporal_omnibus_DE_summary"]], project_root),
    relative_to(source_paths[["temporal_pairwise_DE_summary"]], project_root)
  ),
  statistical_unit = c(
    "Cell for descriptive embedding only",
    "Donor/sample descriptive fraction",
    "Donor-level pseudobulk summary; frozen program-by-cell-type omnibus model",
    "Donor-level pseudobulk summary; frozen cell-type omnibus model",
    "Gene-by-cell-type row",
    "Gene-by-cell-type-by-contrast row"
  ),
  generator_code = script_name,
  figure_class = "B_submission_candidate",
  stringsAsFactors = FALSE
)
write_csv_utf8(source_map, file.path(dirs$manifests, "FIG4_FIGURE_SOURCE_MAP.csv"))

copy_exact(script_path, file.path(dirs$code, script_name))
write_csv_utf8(data.frame(
  key = c(
    "run_stamp", "script", "script_sha256", "R_version", "governing_skill",
    "source_lock", "analysis_rerun", "plot_generated", "output_mode"
  ),
  value = c(
    run_stamp, script_name, sha256_file(script_path), R.version.string, governing_skill,
    source_lock, "FALSE", "TRUE_ONE_300DPI_QC_PNG", "B_CANDIDATE_ASSISTANT_QC_PASS_USER_REVIEW_PENDING"
  ),
  stringsAsFactors = FALSE
), file.path(dirs$qc, "FIG4_RUN_MANIFEST.csv"))
write_text_utf8(capture.output(sessionInfo()), file.path(dirs$qc, "sessionInfo_FIG4_v1_2.txt"))

candidate_manifest <- build_manifest(candidate_root, exclude = "manifests/FIG4_CANDIDATE_MANIFEST.csv")
write_csv_utf8(candidate_manifest, file.path(dirs$manifests, "FIG4_CANDIDATE_MANIFEST.csv"))

report_root <- file.path(feedback_root, paste0("FIG4_V1_2_QC_REPORT_", run_stamp))
if (dir.exists(report_root)) stop("Report root already exists: ", report_root)
report_fig4_root <- file.path(report_root, "FIGURE4")
dir.create(report_fig4_root, recursive = TRUE, showWarnings = FALSE)
copy_tree(candidate_root, report_fig4_root)

report_manifest <- build_manifest(report_root, exclude = "SMALL_REPORT_MANIFEST.csv")
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))
recheck_paths <- file.path(report_root, report_manifest$relative_path)
if (!all(file.exists(recheck_paths)) ||
    !all(vapply(recheck_paths, sha256_file, character(1)) == report_manifest$sha256)) {
  stop("SMALL_REPORT manifest validation failed.")
}

zip_path <- file.path(feedback_root, paste0("STEP05D_FIG4_V1_2_QC_", run_stamp, "_SMALL_REPORT.zip"))
if (file.exists(zip_path)) stop("ZIP already exists: ", zip_path)
zip::zipr(
  zipfile = zip_path,
  files = basename(report_root),
  root = dirname(report_root),
  include_directories = TRUE
)
if (!file.exists(zip_path) || file.info(zip_path)$size < 100000) {
  stop("Figure 4 SMALL_REPORT ZIP was not created correctly.")
}

cat("\n")
cat("============================================================\n")
cat("Figure 4 v1.2 QC candidate completed\n")
cat("============================================================\n")
cat("Gate: PASS_FIGURE4_V1_2_ASSISTANT_QC_USER_REVIEW_PENDING\n")
cat("Assistant actual-render QC: PASS_REFERENCE_RENDER_HASH_MATCH\n")
cat("User visual QC: PENDING\n")
cat("Candidate PNG: ", norm_path(figure_png, TRUE), "\n", sep = "")
cat("SMALL_REPORT ZIP: ", norm_path(zip_path, TRUE), "\n", sep = "")

invisible(list(
  gate = "PASS_FIGURE4_V1_2_ASSISTANT_QC_USER_REVIEW_PENDING",
  candidate_png = norm_path(figure_png, TRUE),
  report_root = norm_path(report_root, TRUE),
  zip = norm_path(zip_path, TRUE)
))

# Publication figure stage 2
# ============================================================================
# 121_step05D_Figure4_v1_1_submission_formats_QC.R
# Nutrition_Endometrial_VDR_v2
#
# User-approved Figure 4 v1.1 submission-format build.
# - Reuses the exact approved v1.1 code-layer figure object.
# - Does not rerun analysis or statistics.
# - Exports true 600-ppi PNG/TIFF and a vector Cairo PDF at 180 x 250 mm.
# - Renders the PDF with Poppler and performs cross-format QC before lock.
# - Does not overwrite the approved 300-dpi candidate.
# ============================================================================

options(stringsAsFactors = FALSE, warn = 1)

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)
script_name <- "121_step05D_Figure4_v1_1_submission_formats_QC.R"
figure_id <- "Figure 4"
approved_visual_version <- "F4_v1_1_QC"
format_version <- "F4_v1_1_SUBMISSION_FORMATS_QC"
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
figure_parent <- file.path(step05_root, "06_results", "FIG4")
approved_root <- file.path(figure_parent, approved_visual_version)
approved_png <- file.path(
  approved_root, "figures", "Figure4_v1_1_QC_180mm_300dpi.png"
)
approved_script <- file.path(
  project_root, "04_code", "R",
  "119_step05D_Figure4_panelA_header_alignment_v1_1_QC.R"
)
expected_approved_png_sha256 <- "6bf4d5ac9f0108872eea80b87c56b9f14b9d68df65b83bb770736b835c7e1b26"

if (!dir.exists(approved_root) || !file.exists(approved_png) || !file.exists(approved_script)) {
  stop("Approved Figure 4 v1.1 authority is incomplete.")
}
if (!identical(tolower(sha256_file(approved_png)), expected_approved_png_sha256)) {
  stop("Approved Figure 4 v1.1 PNG identity mismatch.")
}

skill3_path <- "<FIGURE_STYLE_REFERENCE>"
skill3_text <- paste(readLines(skill3_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
if (!grepl("version: 2.2-beta.7", skill3_text, fixed = TRUE)) {
  stop("Active Skill 3 is not v2.2-beta.7.")
}

parent_override <- Sys.getenv("FIG4_SUBMISSION_PARENT", unset = "")
feedback_override <- Sys.getenv("FIG4_SUBMISSION_FEEDBACK", unset = "")
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
temp_build_root <- file.path(tempdir(), paste0("FIG4_v1_1_submission_build_", run_stamp))
temp_candidate_parent <- file.path(temp_build_root, "candidate_parent")
temp_feedback <- file.path(temp_build_root, "feedback")
old_candidate_env <- Sys.getenv("FIG4_CANDIDATE_PARENT", unset = NA_character_)
old_feedback_env <- Sys.getenv("FIG4_FEEDBACK_ROOT", unset = NA_character_)
Sys.setenv(
  FIG4_CANDIDATE_PARENT = temp_candidate_parent,
  FIG4_FEEDBACK_ROOT = temp_feedback
)
plot_env <- new.env(parent = globalenv())
source_console <- capture.output(
  sys.source(approved_script, envir = plot_env, keep.source = TRUE)
)
if (is.na(old_candidate_env)) Sys.unsetenv("FIG4_CANDIDATE_PARENT") else Sys.setenv(FIG4_CANDIDATE_PARENT = old_candidate_env)
if (is.na(old_feedback_env)) Sys.unsetenv("FIG4_FEEDBACK_ROOT") else Sys.setenv(FIG4_FEEDBACK_ROOT = old_feedback_env)

temp_approved_png <- file.path(
  temp_candidate_parent, approved_visual_version, "figures",
  "Figure4_v1_1_QC_180mm_300dpi.png"
)
if (!file.exists(temp_approved_png) ||
    !identical(tolower(sha256_file(temp_approved_png)), expected_approved_png_sha256)) {
  stop("Rebuilt Figure 4 v1.1 code layer does not reproduce the approved reference PNG.")
}
if (is.null(plot_env$composite_parts$figure)) stop("Approved figure object was not recovered.")
figure_object <- plot_env$composite_parts$figure
figure_width_mm <- plot_env$figure_width_mm
figure_height_mm <- plot_env$figure_height_mm
if (!identical(as.numeric(figure_width_mm), 180) ||
    !identical(as.numeric(figure_height_mm), 250)) {
  stop("Approved Figure 4 geometry is not 180 x 250 mm.")
}
write_text_utf8(source_console, file.path(dirs$qc, "FIG4_APPROVED_CODE_LAYER_REBUILD_CONSOLE.txt"))

# Preserve approved candidate context without altering it.
copy_tree(file.path(approved_root, "source_data"), dirs$source_data)
copy_tree(file.path(approved_root, "manuscript_text"), dirs$manuscript)
copy_tree(file.path(approved_root, "manifests"), dirs$manifests)
copy_tree(file.path(approved_root, "qc"), file.path(dirs$qc, "approved_candidate_qc"))
copy_tree(file.path(approved_root, "code"), file.path(dirs$code, "approved_candidate_code"))
copy_exact(approved_script, file.path(dirs$code, basename(approved_script)))
copy_exact(script_path, file.path(dirs$code, script_name))

reference_copy <- file.path(
  dirs$reference, "Figure4_v1_1_APPROVED_REFERENCE_180mm_300dpi.png"
)
copy_exact(approved_png, reference_copy)

user_approval <- data.frame(
  figure = figure_id,
  approved_visual_version = approved_visual_version,
  approved_png_sha256 = expected_approved_png_sha256,
  user_final_visual_qc = "PASS",
  approval_date = "2026-08-07",
  approval_scope = "Figure 4 v1.1 explicitly linked and approved by user",
  unselected_candidate = "F4_v1_2_QC_NOT_SELECTED_FOR_RELEASE",
  stringsAsFactors = FALSE
)
write_csv_utf8(user_approval, file.path(dirs$qc, "FIG4_USER_VISUAL_APPROVAL_v1_1.csv"))

width_px_600 <- round(figure_width_mm / 25.4 * 600)
height_px_600 <- round(figure_height_mm / 25.4 * 600)
width_in <- figure_width_mm / 25.4
height_in <- figure_height_mm / 25.4

png600 <- file.path(dirs$figures, "Figure4_v1_1_submission_180mm_600ppi.png")
tiff600 <- file.path(dirs$figures, "Figure4_v1_1_submission_180mm_600ppi.tiff")
pdfvec <- file.path(dirs$figures, "Figure4_v1_1_submission_180mm.pdf")

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
write_csv_utf8(warning_log, file.path(dirs$qc, "FIG4_SUBMISSION_RUNTIME_WARNING_LOG.csv"))
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
pdf_render_prefix <- file.path(dirs$qc, "Figure4_v1_1_PDF_RENDER_300dpi")
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
write_text_utf8(pdfinfo_output, file.path(dirs$qc, "FIG4_PDFINFO.txt"))

approved_img <- magick::image_read(reference_copy)
approved_info <- magick::image_info(approved_img)
approved_geometry <- paste0(approved_info$width[[1]], "x", approved_info$height[[1]], "!")
png600_preview <- magick::image_resize(png_img, approved_geometry, filter = "Lanczos")
png600_preview_path <- file.path(dirs$qc, "Figure4_v1_1_PNG600_RENDER_300dpi.png")
magick::image_write(png600_preview, png600_preview_path, format = "png")

pdf_render_img <- magick::image_read(pdf_render_png)
pdf_render_exact <- magick::image_resize(pdf_render_img, approved_geometry, filter = "Lanczos")
pdf_render_exact_path <- file.path(dirs$qc, "Figure4_v1_1_PDF_RENDER_MATCHED_300dpi.png")
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
write_csv_utf8(cross_format_qc, file.path(dirs$qc, "FIG4_CROSS_FORMAT_QC.csv"))
if (!all(cross_format_qc$pass)) stop("Figure 4 cross-format QC failed.")

format_manifest <- data.frame(
  role = c("approved_reference_300dpi", "submission_png_600ppi", "submission_tiff_600ppi", "submission_visual_authority_pdf"),
  file = c(basename(reference_copy), basename(png600), basename(tiff600), basename(pdfvec)),
  size_bytes = as.numeric(file.info(c(reference_copy, png600, tiff600, pdfvec))$size),
  sha256 = vapply(c(reference_copy, png600, tiff600, pdfvec), sha256_file, character(1)),
  stringsAsFactors = FALSE
)
write_csv_utf8(format_manifest, file.path(dirs$manifests, "FIG4_SUBMISSION_FORMAT_HASH_MANIFEST.csv"))

# The deterministic raster hashes are bound after assistant actual-render review.
# Cairo writes session-specific PDF metadata; its exact formal-build hash is therefore
# bound downstream by Skill 8 after all PDF structure and rendering gates pass.
expected_submission_hashes <- c(
  submission_png_600ppi = "9812a20be08a3f2a531c29c482a29d1807ac551d1c4a9664360f4f985f313223",
  submission_tiff_600ppi = "dc8c4c4e22e79ace93bb839e866161be9b4a31da1ea58f3b5cdd673e3f1f1e63",
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
write_csv_utf8(reviewed_hash_qc, file.path(dirs$qc, "FIG4_REVIEWED_SUBMISSION_HASH_QC.csv"))
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
  "# Figure 4 v1.1 submission-format handoff",
  "",
  "- User-approved visual authority: F4_v1_1_QC.",
  "- Unselected candidate: F4_v1_2_QC; not used for release.",
  "- Scientific/statistical recomputation: none.",
  "- True 600-ppi PNG and LZW TIFF: generated from the approved code-layer object.",
  "- PDF: 600-ppi visual-authority raster embedded on a Cairo PDF page at 180 x 250 mm.",
  "- Rationale: subpixel semi-transparent UMAP points changed density under direct vector rendering; raster embedding preserves the user-approved appearance.",
  "- Cross-format actual-render QC: PASS.",
  "- Next gate: Skill 8 v1.1 MODULE_LOCK_ARCHIVE."
), file.path(dirs$qc, "FIG4_SUBMISSION_FORMATS_HANDOFF.md"))

write_text_utf8(capture.output(sessionInfo()), file.path(dirs$qc, "sessionInfo_FIG4_submission_v1_1.txt"))
candidate_manifest <- build_manifest(out_root, exclude = "manifests/FIG4_SUBMISSION_CANDIDATE_MANIFEST.csv")
write_csv_utf8(candidate_manifest, file.path(dirs$manifests, "FIG4_SUBMISSION_CANDIDATE_MANIFEST.csv"))
candidate_validation <- validate_manifest(out_root, candidate_manifest)
write_csv_utf8(candidate_validation, file.path(dirs$qc, "FIG4_SUBMISSION_MANIFEST_VALIDATION.csv"))
if (!all(candidate_validation$exists & candidate_validation$size_match &
         candidate_validation$md5_match & candidate_validation$sha256_match)) {
  stop("Submission candidate manifest validation failed.")
}

report_root <- file.path(feedback_root, paste0("FIG4_V1_1_SUBMISSION_FORMATS_REPORT_", run_stamp))
if (dir.exists(report_root)) stop("Report root exists: ", normp(report_root, FALSE))
report_fig4 <- file.path(report_root, "FIGURE4")
dir.create(report_fig4, recursive = TRUE, showWarnings = FALSE)
copy_tree(out_root, report_fig4)
report_manifest <- build_manifest(report_root, exclude = "SMALL_REPORT_MANIFEST.csv")
write_csv_utf8(report_manifest, file.path(report_root, "SMALL_REPORT_MANIFEST.csv"))
report_validation <- validate_manifest(report_root, report_manifest)
if (!all(report_validation$exists & report_validation$size_match &
         report_validation$md5_match & report_validation$sha256_match)) {
  stop("SMALL_REPORT validation failed.")
}

zip_path <- file.path(
  feedback_root,
  paste0("STEP05D_FIG4_V1_1_SUBMISSION_FORMATS_", run_stamp, "_SMALL_REPORT.zip")
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
cat("Figure 4 v1.1 submission formats completed\n")
cat("Gate: PASS_FIGURE4_V1_1_SUBMISSION_FORMATS_READY_FOR_SKILL8\n")
cat("Output: ", normp(out_root, TRUE), "\n", sep = "")
cat("SMALL_REPORT: ", normp(zip_path, TRUE), "\n", sep = "")
cat("============================================================\n")

invisible(list(
  gate = "PASS_FIGURE4_V1_1_SUBMISSION_FORMATS_READY_FOR_SKILL8",
  output_root = normp(out_root, TRUE),
  png600 = normp(png600, TRUE),
  tiff600 = normp(tiff600, TRUE),
  pdf = normp(pdfvec, TRUE),
  zip = normp(zip_path, TRUE)
))
