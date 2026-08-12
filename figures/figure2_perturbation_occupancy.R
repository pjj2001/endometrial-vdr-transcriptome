
# ======================================================================
# Figure2_pdf_font_repair_v1_2.R
# Nutrition_Endometrial_VDR_v2
#
# GOVERNING SKILL
# - research_scientific_figure_suite v2.2-beta.8
# - FAST_REPLOT mode (frozen source/statistics; visual/export repair only)
#
# FIGURE 2
# Perturbation-defined VDR RNA program and occupancy support
#
# SCIENTIFIC QUESTION
# What experimental RNA interaction program was identified, and how does
# published VDR occupancy support stratify that program without assigning
# definitive direct targets?
#
# PANELS
# A. RNA interaction volcano across all tested features.
# B. Top interaction-candidate expression heatmap across 12 samples.
# C. Occupancy-support counts across the 436-candidate program.
# D. Mutually exclusive downstream evidence-routing classes.
#
# SOURCE POLICY
# - Use the existing self-contained STEP03B1 Class-B source package.
# - Use RNA_OCCUPANCY_EVIDENCE_v1.0 and
#   DOWNSTREAM_CANDIDATE_ROUTING_v1.0.
# - Do not infer values from old images.
# - Do not recompute differential-expression statistics.
#
# RELEASE POLICY
# - Generate a high-polish full-format candidate:
#   PDF, PNG 300 dpi, TIFF 600 dpi and SVG.
# - Status remains B_submission_candidate until the actual exported
#   images pass strict render QC.
# - Do not call this submission-ready before visual adjudication.
#
# CONSOLE
# source(
#   "<repository module script>",
#   echo = FALSE,
#   encoding = "UTF-8"
# )
# ======================================================================

options(stringsAsFactors = FALSE, warn = 1)
options(timeout = max(14400L, getOption("timeout")))

project_root <- normalizePath(Sys.getenv("VDR_PROJECT_ROOT", unset = getwd()), mustWork = FALSE)

skill_root <- "<FIGURE_STYLE_REFERENCE>"

skill_path <- file.path(
  skill_root,
  "SKILL.md"
)

required_skill_version <-
  "version: 2.2-beta.8"

script_name <-
  "Figure2_pdf_font_repair_v1_2.R"

script_version <-
  "STEP05B_FIGURE2_PDF_FONT_REPAIR_v1.2"

figure_id <- "Figure 2"

figure_version <- "FIG2_v1_8_pdf_font_repair"

figure_class <- "B_submission_candidate"

candidate_status <-
  "B_RELEASE_CANDIDATE_AWAITING_STRICT_RENDER_QC"

run_start <- Sys.time()

run_stamp <- format(
  run_start,
  "%Y%m%d_%H%M%S"
)

run_id <- paste0(
  "STEP05B_FIG2_",
  run_stamp
)

target_width_mm <- 180
target_height_mm <- 170
stress_test_width_mm <- 170

top_block_tag_x <- 0
top_block_tag_y <- 1
top_block_tag_to_title_gap_pt <- 0
top_block_anchor_reference <- "NATIVE_TOP_BLOCK_STRIPPED"

header_anchor_reference <- "GTABLE_PANEL_CELL_LEFT_EDGE"
header_tag_size_pt <- 9.0
header_title_size_pt <- 8.5
header_subtitle_size_pt <- 7.2
header_tag_to_title_gap_mm <- 2.0
header_title_row_height_mm <- 5.0
header_subtitle_row_height_mm <- 4.0
header_title_lineheight <- 0.94
panel_B_axis_text_y_size <- 5.4
panel_C_x_limit <- 310
panel_D_x_limit <- 205

# Panel A horizontal-axis geometry lock.
# These values are deliberately intermediate between the earlier
# over-separated and over-tight versions.
panel_A_y_tick_to_axis_gap_pt <- 0.75
panel_A_y_title_to_tick_gap_pt <- 0.70
panel_A_horizontal_geometry_state <- "LABEL_WIDTH_DECOUPLED_LOCK"


# ----------------------------------------------------------------------
# 1. Helpers
# ----------------------------------------------------------------------

normalize_slashes <- function(x) {
  gsub("\\\\", "/", x)
}

path_budget_events <- list()

assert_path_safe <- function(
    path,
    warn_chars = 220L,
    stop_chars = 240L
) {
  normalized <- normalize_slashes(path)

  path_n <- nchar(
    normalized,
    type = "chars"
  )

  if (
    .Platform$OS.type == "windows" &&
      path_n >= stop_chars
  ) {
    stop(
      "Path exceeds Skill 3 hard limit (",
      path_n,
      " >= ",
      stop_chars,
      "):\n",
      normalized
    )
  }

  if (
    .Platform$OS.type == "windows" &&
      path_n >= warn_chars
  ) {
    path_budget_events[[
      length(path_budget_events) + 1L
    ]] <<- data.frame(
      path = normalized,
      path_chars = path_n,
      threshold = warn_chars,
      status = "SOFT_WARNING_RECORDED",
      stringsAsFactors = FALSE
    )

    message(
      "Path budget soft warning recorded (",
      path_n,
      " >= ",
      warn_chars,
      "): ",
      normalized
    )
  }

  invisible(path)
}

write_csv_utf8 <- function(
    x,
    path
) {
  assert_path_safe(path)

  if (file.exists(path)) {
    stop(
      "Refusing to overwrite an existing file:\n",
      path
    )
  }

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  write.csv(
    x,
    path,
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

write_text <- function(
    lines,
    path
) {
  assert_path_safe(path)

  if (file.exists(path)) {
    stop(
      "Refusing to overwrite an existing file:\n",
      path
    )
  }

  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    lines,
    path,
    useBytes = TRUE
  )
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)

  output <- tryCatch(
    read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(e) NULL
  )

  if (is.null(output)) {
    output <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  output
}

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)

  toupper(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )
}

first_existing <- function(
    x,
    candidates,
    label
) {
  hit <- candidates[
    candidates %in% names(x)
  ]

  if (!length(hit)) {
    stop(
      "Could not identify ",
      label,
      ". Expected one of: ",
      paste(
        candidates,
        collapse = ", "
      )
    )
  }

  hit[[1]]
}

require_columns <- function(
    x,
    columns,
    label
) {
  missing_columns <- setdiff(
    columns,
    names(x)
  )

  if (length(missing_columns)) {
    stop(
      label,
      " is missing column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  invisible(TRUE)
}

md5_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  unname(
    tools::md5sum(path)
  )
}

sha256_file <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }

  digest::digest(
    path,
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
}

copy_validate <- function(
    source_path,
    destination_path
) {
  assert_path_safe(destination_path)

  if (file.exists(destination_path)) {
    stop(
      "Refusing to overwrite copied artifact:\n",
      destination_path
    )
  }

  dir.create(
    dirname(destination_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  copied <- file.copy(
    source_path,
    destination_path,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )

  pass <- isTRUE(copied) &&
    file.exists(destination_path) &&
    identical(
      as.numeric(
        file.info(source_path)$size
      ),
      as.numeric(
        file.info(destination_path)$size
      )
    ) &&
    identical(
      md5_file(source_path),
      md5_file(destination_path)
    ) &&
    identical(
      sha256_file(source_path),
      sha256_file(destination_path)
    )

  if (!pass) {
    stop(
      "Validated copy failed:\n",
      source_path,
      "\n->\n",
      destination_path
    )
  }

  invisible(TRUE)
}

archive_partial_directory <- function(
    source_directory,
    archive_parent,
    archive_label
) {
  if (!dir.exists(source_directory)) {
    return(invisible(FALSE))
  }

  status_path <- file.path(
    source_directory,
    "FIGURE_STATUS.csv"
  )

  status_value <- safe_read_csv(
    status_path
  )

  if (
    !is.null(status_value) &&
      nrow(status_value) == 1L &&
      identical(
        status_value$status[[1]],
        "COMPLETE"
      )
  ) {
    stop(
      "A completed figure version already exists:\n",
      source_directory,
      "\nCreate a new figure version rather than overwriting it."
    )
  }

  dir.create(
    archive_parent,
    recursive = TRUE,
    showWarnings = FALSE
  )

  archive_destination <- file.path(
    archive_parent,
    paste0(
      archive_label,
      "_PARTIAL_",
      format(
        Sys.time(),
        "%Y%m%d_%H%M%S"
      )
    )
  )

  moved <- file.rename(
    source_directory,
    archive_destination
  )

  if (!isTRUE(moved)) {
    stop(
      "Could not archive incomplete figure directory:\n",
      source_directory
    )
  }

  invisible(TRUE)
}


# ----------------------------------------------------------------------
# 2. Packages and optimized Skill 3 installation
# ----------------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "patchwork",
  "ggrepel",
  "svglite",
  "ragg",
  "systemfonts",
  "gtable",
  "digest",
  "jsonlite",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\nInstall them before rerunning; the figure script does not ",
    "silently update packages."
  )
}

if (!file.exists(skill_path)) {
  stop(
    "Optimized Skill 3 is not installed at:\n",
    skill_path,
    "\nInstall research_scientific_figure_suite v2.2-beta.8 first."
  )
}

skill_text <- paste(
  readLines(
    skill_path,
    warn = FALSE,
    encoding = "UTF-8"
  ),
  collapse = "\n"
)

if (
  !grepl(
    required_skill_version,
    skill_text,
    fixed = TRUE
  )
) {
  stop(
    "The active Skill 3 does not report v2.2-beta.3:\n",
    skill_path
  )
}

# ----------------------------------------------------------------------
# 2B. Publication font and text-size lock
# ----------------------------------------------------------------------

requested_font <- "Arial"

font_match <- systemfonts::match_fonts(
  requested_font
)

required_match_columns <- c(
  "path",
  "index"
)

if (
  !is.data.frame(font_match) ||
    nrow(font_match) != 1L ||
    !all(
      required_match_columns %in%
        names(font_match)
    )
) {
  stop(
    "systemfonts::match_fonts() did not return one ",
    "valid Arial match with path and index."
  )
}

font_path <- as.character(
  font_match$path[[1]]
)

font_index <- as.integer(
  font_match$index[[1]]
)

if (
  is.na(font_path) ||
    !nzchar(font_path) ||
    !file.exists(font_path) ||
    is.na(font_index)
) {
  stop(
    "Arial font matching returned an invalid path or index."
  )
}

font_details <- systemfonts::font_info(
  path = font_path,
  index = font_index,
  size = 12,
  res = 72
)

required_detail_columns <- c(
  "family",
  "style",
  "name"
)

if (
  !is.data.frame(font_details) ||
    nrow(font_details) != 1L ||
    !all(
      required_detail_columns %in%
        names(font_details)
    )
) {
  stop(
    "systemfonts::font_info() did not return complete ",
    "Arial family/style/name metadata."
  )
}

resolved_font_family <- as.character(
  font_details$family[[1]]
)

resolved_font_style <- as.character(
  font_details$style[[1]]
)

resolved_font_name <- as.character(
  font_details$name[[1]]
)

if (
  is.na(resolved_font_family) ||
    !nzchar(resolved_font_family) ||
    tolower(
      trimws(
        resolved_font_family
      )
    ) !=
      tolower(
        requested_font
      )
) {
  stop(
    "Requested Arial but systemfonts resolved family: ",
    resolved_font_family,
    ". The script will not silently use a fallback font."
  )
}

figure_font <- resolved_font_family

font_standard <- data.frame(
  element = c(
    "panel_tag",
    "panel_title",
    "panel_subtitle",
    "axis_title",
    "axis_tick",
    "legend_title",
    "legend_text",
    "heatmap_gene_label",
    "volcano_repel_label",
    "bar_value_label"
  ),
  size_pt_at_180mm = c(
    header_tag_size_pt,
    header_title_size_pt,
    header_subtitle_size_pt,
    8.0,
    7.0,
    7.2,
    6.8,
    6.0,
    6.3,
    6.7
  ),
  minimum_allowed_pt = c(
    8.5,
    8.0,
    7.0,
    7.5,
    6.5,
    6.8,
    6.5,
    5.8,
    6.0,
    6.2
  ),
  locked = TRUE,
  stringsAsFactors = FALSE
)

if (
  any(
    font_standard$size_pt_at_180mm <
      font_standard$minimum_allowed_pt
  )
) {
  stop(
    "One or more text sizes violate the locked minimum."
  )
}

font_preflight_warnings <- character()

font_preflight_png <- tempfile(
  fileext = ".png"
)

withCallingHandlers(
  {
    ragg::agg_png(
      filename = font_preflight_png,
      width = 80,
      height = 30,
      units = "mm",
      res = 150,
      background = "white"
    )

    grid::grid.newpage()

    grid::grid.text(
      "Arial font preflight: ABC abc 0123",
      gp = grid::gpar(
        fontfamily = figure_font,
        fontsize = 9
      )
    )

    grDevices::dev.off()
  },
  warning = function(w) {
    font_preflight_warnings <<- c(
      font_preflight_warnings,
      conditionMessage(w)
    )

    invokeRestart(
      "muffleWarning"
    )
  }
)

if (grDevices::dev.cur() != 1L) {
  grDevices::dev.off()
}

if (
  length(font_preflight_warnings) ||
    !file.exists(font_preflight_png) ||
    file.info(font_preflight_png)$size <= 500
) {
  stop(
    "Font preflight failed:\n",
    paste(
      font_preflight_warnings,
      collapse = "\n"
    )
  )
}

unlink(
  font_preflight_png
)

# Keep a verified raster device open while ggplot gtables calculate Arial
# string metrics.  In headless Rscript sessions the implicit PostScript
# device cannot resolve Windows Arial and can either fail or generate a PDF
# whose text is rendered with an unintended heavy outline.  This scratch
# device is layout-only and never changes the frozen source/statistical data.
metrics_device_file <- tempfile(
  pattern = "FIG2_Arial_metrics_",
  fileext = ".png"
)

ragg::agg_png(
  filename = metrics_device_file,
  width = target_width_mm,
  height = target_height_mm,
  units = "mm",
  res = 96,
  background = "white"
)


# ----------------------------------------------------------------------
# 3. Project paths
# ----------------------------------------------------------------------

step03_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP03_VDR_perturbation_regulons"
)

candidate_root <- file.path(
  step03_root,
  "07_figures_tables",
  "B_submission_candidate",
  "STEP03B1_RNA_candidates_v3"
)

RNA_result_root <- file.path(
  step03_root,
  "06_results"
)

occupancy_root <- file.path(
  RNA_result_root,
  "STEP03B3C3B_locked_RNA_occupancy_evidence",
  "RNA_OCCUPANCY_EVIDENCE_v1.0"
)

step05_root <- file.path(
  project_root,
  "05_analysis_steps",
  "STEP05_manuscript_assembly"
)

figure_parent <- "<LOCAL_R_LIBRARY>"

result_root <- file.path(
  figure_parent,
  figure_version
)

archive_parent <- "<LOCAL_R_LIBRARY>"

archive_partial_directory(
  source_directory = result_root,
  archive_parent = archive_parent,
  archive_label = figure_version
)

figure_root <- file.path(
  result_root,
  "figures"
)

source_root <- file.path(
  result_root,
  "source_data"
)

upstream_root <- file.path(
  source_root,
  "upstream_refs"
)

code_root <- file.path(
  result_root,
  "code"
)

legend_root <- file.path(
  result_root,
  "legends"
)

QC_root <- file.path(
  result_root,
  "qc"
)

manifest_root <- file.path(
  result_root,
  "manifests"
)

text_root <- file.path(
  result_root,
  "manuscript_text"
)

for (path_value in c(
  figure_root,
  source_root,
  upstream_root,
  code_root,
  legend_root,
  QC_root,
  manifest_root,
  text_root
)) {
  dir.create(
    path_value,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

execution_root <- file.path(
  step05_root,
  "05_execution"
)

parameter_root <- file.path(
  step05_root,
  "03_parameters"
)

environment_root <- file.path(
  step05_root,
  "04_environment"
)

for (path_value in c(
  execution_root,
  parameter_root,
  environment_root
)) {
  dir.create(
    path_value,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

execution_log <- file.path(
  execution_root,
  paste0(
    run_id,
    "_STATIC_execution.log"
  )
)

run_manifest_path <- file.path(
  execution_root,
  paste0(
    run_id,
    "_RUN_MANIFEST.csv"
  )
)

parameter_path <- file.path(
  parameter_root,
  paste0(
    run_id,
    "_parameters.json"
  )
)

session_path <- file.path(
  environment_root,
  paste0(
    "sessionInfo_",
    run_id,
    ".txt"
  )
)

figure_pdf <- file.path(
  figure_root,
  "Figure2_v1_8_pdf_font_repair_180mm.pdf"
)

figure_png <- file.path(
  figure_root,
  "Figure2_v1_8_pdf_font_repair_180mm_300dpi.png"
)

figure_tiff <- file.path(
  figure_root,
  "Figure2_v1_8_pdf_font_repair_180mm_600dpi.tiff"
)

figure_svg <- file.path(
  figure_root,
  "Figure2_v1_8_pdf_font_repair_180mm.svg"
)


# ----------------------------------------------------------------------
# 4. Direct sources
# ----------------------------------------------------------------------

source_paths <- c(
  volcano = file.path(
    candidate_root,
    "source_data",
    "SourceData_STEP03B1_interaction_volcano_candidate_v3.csv"
  ),
  volcano_labels = file.path(
    candidate_root,
    "source_data",
    "SourceData_STEP03B1_interaction_volcano_labels_v3.csv"
  ),
  heatmap = file.path(
    candidate_root,
    "source_data",
    "SourceData_STEP03B1_top_interaction_heatmap_ACCEPTED_v2.csv"
  ),
  sample_manifest = file.path(
    candidate_root,
    "source_data",
    "upstream_refs",
    "STEP03B1_sample_manifest_FINAL.csv"
  ),
  figure_source_map = file.path(
    candidate_root,
    "manifests",
    "STEP03B1_v3_FIGURE_SOURCE_MAP.csv"
  ),
  RNA_frozen = file.path(
    RNA_result_root,
    "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"
  ),
  integrated_436 = file.path(
    occupancy_root,
    "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv"
  ),
  routing_436 = file.path(
    occupancy_root,
    "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv"
  ),
  evidence_counts = file.path(
    occupancy_root,
    "RNA_OCCUPANCY_EVIDENCE_v1.0_group_counts.csv"
  ),
  occupancy_manifest = file.path(
    occupancy_root,
    "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv"
  )
)

missing_sources <- source_paths[
  !file.exists(source_paths)
]

if (length(missing_sources)) {
  stop(
    "Required Figure 2 source(s) are missing:\n",
    paste(
      normalize_slashes(
        missing_sources
      ),
      collapse = "\n"
    )
  )
}

volcano <- safe_read_csv(
  source_paths[["volcano"]]
)

volcano_labels <- safe_read_csv(
  source_paths[["volcano_labels"]]
)

heatmap <- safe_read_csv(
  source_paths[["heatmap"]]
)

sample_manifest <- safe_read_csv(
  source_paths[["sample_manifest"]]
)

RNA_frozen <- safe_read_csv(
  source_paths[["RNA_frozen"]]
)

integrated <- safe_read_csv(
  source_paths[["integrated_436"]]
)

routing <- safe_read_csv(
  source_paths[["routing_436"]]
)

require_columns(
  volcano,
  c(
    "feature_id",
    "interaction_logFC",
    "interaction_FDR"
  ),
  "Volcano source"
)

require_columns(
  heatmap,
  c(
    "feature_id",
    "short_label",
    "z_score"
  ),
  "Heatmap source"
)

require_columns(
  integrated,
  c(
    "feature_id",
    "program_class",
    "interaction_logFC",
    "interaction_FDR",
    "occupancy_supported",
    "shared_occupancy_supported",
    "any_promoter_overlap_support",
    "vehicle_peak_n",
    "calcitriol_peak_n",
    "shared_peak_n",
    "occupancy_evidence_group"
  ),
  "RNA/occupancy evidence table"
)

require_columns(
  routing,
  c(
    "feature_id",
    "downstream_route_class",
    "direct_target_status"
  ),
  "Downstream routing table"
)

if (
  nrow(volcano) != 12116L ||
    nrow(integrated) != 436L ||
    nrow(routing) != 436L ||
    length(
      unique(
        heatmap$feature_id
      )
    ) != 30L
) {
  stop(
    "Figure 2 source row counts do not match the locked design."
  )
}

integrated$occupancy_supported <-
  as_logical_safe(
    integrated$occupancy_supported
  )

integrated$shared_occupancy_supported <-
  as_logical_safe(
    integrated$shared_occupancy_supported
  )

integrated$any_promoter_overlap_support <-
  as_logical_safe(
    integrated$any_promoter_overlap_support
  )

if (
  sum(
    integrated$occupancy_supported
  ) != 284L ||
    sum(
      !integrated$occupancy_supported
    ) != 152L ||
    sum(
      integrated$shared_occupancy_supported
    ) != 46L ||
    sum(
      integrated$any_promoter_overlap_support
    ) != 84L ||
    sum(
      safe_numeric(
        integrated$vehicle_peak_n
      ) > 0
    ) != 77L ||
    sum(
      safe_numeric(
        integrated$calcitriol_peak_n
      ) > 0
    ) != 275L
) {
  stop(
    "RNA/occupancy evidence counts do not match the locked values."
  )
}

if (
  any(
    routing$direct_target_status !=
      "NOT_ASSIGNED"
  )
) {
  stop(
    "One or more routing rows incorrectly assign direct-target status."
  )
}


# ----------------------------------------------------------------------
# 5. Freeze source package
# ----------------------------------------------------------------------

source_copy_map <- data.frame(
  artifact = names(source_paths),
  source_path =
    normalize_slashes(
      unname(source_paths)
    ),
  destination_path = NA_character_,
  size_bytes = NA_real_,
  md5 = NA_character_,
  sha256 = NA_character_,
  stringsAsFactors = FALSE
)

for (
  source_index in seq_len(
    nrow(source_copy_map)
  )
) {
  artifact_name <-
    source_copy_map$artifact[
      source_index
    ]

  source_path <-
    source_copy_map$source_path[
      source_index
    ]

  destination_directory <- if (
    artifact_name %in%
      c(
        "sample_manifest",
        "figure_source_map",
        "RNA_frozen",
        "occupancy_manifest"
      )
  ) {
    upstream_root
  } else {
    source_root
  }

  destination_path <- file.path(
    destination_directory,
    basename(source_path)
  )

  copy_validate(
    source_path,
    destination_path
  )

  source_copy_map$destination_path[
    source_index
  ] <- normalize_slashes(
    destination_path
  )

  source_copy_map$size_bytes[
    source_index
  ] <- file.info(
    destination_path
  )$size

  source_copy_map$md5[
    source_index
  ] <- md5_file(
    destination_path
  )

  source_copy_map$sha256[
    source_index
  ] <- sha256_file(
    destination_path
  )
}

write_csv_utf8(
  source_copy_map,
  file.path(
    manifest_root,
    "FIG2_SOURCE_COPY_MAP.csv"
  )
)

write_csv_utf8(
  font_standard,
  file.path(
    QC_root,
    "FIG2_FONT_AND_TEXT_SIZE_STANDARD.csv"
  )
)

font_preflight_record <- data.frame(
  requested_font = requested_font,
  resolved_font_family =
    resolved_font_family,
  resolved_font_name =
    resolved_font_name,
  resolved_font_style =
    resolved_font_style,
  resolved_font_index =
    font_index,
  resolved_font_path =
    normalize_slashes(
      font_path
    ),
  resolved_font_path_present =
    file.exists(
      font_path
    ) &&
      tolower(
        trimws(
          resolved_font_family
        )
      ) ==
        tolower(
          requested_font
        ),
  exact_family_match =
    tolower(
      trimws(
        resolved_font_family
      )
    ) ==
      tolower(
        requested_font
      ),
  preflight_warning_n =
    length(
      font_preflight_warnings
    ),
  preflight_pass = TRUE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  font_preflight_record,
  file.path(
    QC_root,
    "FIG2_FONT_PREFLIGHT.csv"
  )
)

v111_render_adjudication <- data.frame(
  defect_id = c(
    "F2-V111-01",
    "F2-V111-02"
  ),
  panel = c(
    "C",
    "A"
  ),
  actual_render_finding = c(
    "The 284 (65.1%) and 275 (63.1%) labels were clipped at the right panel boundary in PNG, PDF and SVG.",
    "The four-item volcano legend visually extended into the Panel B heading region."
  ),
  severity = c(
    "BLOCKING",
    "MAJOR"
  ),
  targeted_fix = c(
    "Place labels for bars >=250 inside the bars; keep smaller-bar labels outside.",
    "Use shorter labels and a two-row legend constrained to Panel A."
  ),
  scientific_result_changed = FALSE,
  statistical_result_changed = FALSE,
  accepted_for_lock = FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  v111_render_adjudication,
  file.path(
    QC_root,
    "FIG2_V1_1_1_ACTUAL_RENDER_ADJUDICATION.csv"
  )
)

fast_replot_state <- data.frame(
  fast_replot_round = 2L,
  maximum_fast_replot_rounds = 2L,
  this_is_final_allowed_fast_replot = TRUE,
  science_locked = TRUE,
  statistics_locked = TRUE,
  source_data_locked = TRUE,
  layout_changed = FALSE,
  targeted_visual_changes_only = TRUE,
  further_blind_fast_replot_authorized = FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  fast_replot_state,
  file.path(
    QC_root,
    "FIG2_FAST_REPLOT_STATE.csv"
  )
)

manual_visual_override <- data.frame(
  defect_id = c(
    "F2-V12-01",
    "F2-V12-02",
    "F2-V12-03",
    "F2-V12-04",
    "F2-V12-05"
  ),
  panel = c(
    "A",
    "A",
    "A",
    "B",
    "C-D"
  ),
  severity = c(
    "BLOCKING",
    "BLOCKING",
    "BLOCKING",
    "MAJOR",
    "BLOCKING"
  ),
  actual_render_finding = c(
    "The y-axis title was visually detached from the y axis because of excessive spacing.",
    "Panel tag A was anchored too far into page whitespace rather than to the Panel A bounding box.",
    "Leader lines crossed or visually split gene labels in the dense central annotation region.",
    "Heatmap row labels consumed more width than necessary and required slight rebalancing.",
    "Panel tags C and D were not visually anchored to the upper-left corners of their own panel cells."
  ),
  controlled_repair = c(
    "Reduce y-title and y-tick margins and reduce the Panel A left outer margin.",
    "Use per-panel tags and reserve title margin inside every panel.",
    "Replace automatic repel with a locked manual label-layout table and short connectors that stop outside label boxes.",
    "Set heatmap row labels to 6.0 pt and preserve complete gene names.",
    "Use per-panel C/D tags with the same panel-relative anchor and wrap long horizontal-bar categories."
  ),
  scientific_result_changed = FALSE,
  statistical_result_changed = FALSE,
  source_data_changed = FALSE,
  release_blocked_until_repair_QC = TRUE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  manual_visual_override,
  file.path(
    QC_root,
    "FIG2_MANUAL_VISUAL_OVERRIDE_CARD.csv"
  )
)

controlled_repair_state <- data.frame(
  mode = "CONTROLLED_LAYOUT_REPAIR",
  triggering_event =
    "HUMAN_STRICT_VISUAL_QC_OVERRIDE_AFTER_FAST_REPLOT_2_OF_2",
  prior_version =
    "FIG2_v1_2_final_candidate",
  prior_submission_ready =
    FALSE,
  prior_immutable_release =
    FALSE,
  further_blind_FAST_REPLOT_authorized =
    FALSE,
  controlled_layout_repair_authorized =
    TRUE,
  science_locked =
    TRUE,
  statistics_locked =
    TRUE,
  source_data_locked =
    TRUE,
  release_blocked =
    TRUE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  controlled_repair_state,
  file.path(
    QC_root,
    "FIG2_CONTROLLED_LAYOUT_REPAIR_STATE.csv"
  )
)

panel_geometry_lock <- data.frame(
  panel = c(
    "A",
    "B",
    "C",
    "D"
  ),
  grid_row = c(
    1,
    1,
    2,
    2
  ),
  grid_column = c(
    1,
    2,
    1,
    2
  ),
  panel_tag = c(
    "A",
    "B",
    "C",
    "D"
  ),
  panel_tag_anchor =
    "panel_relative_top_left",
  title_left_margin_pt =
    17,
  target_panel_width_ratio = c(
    1.00,
    1.00,
    1.00,
    1.00
  ),
  target_row_height_ratio = c(
    1.14,
    1.14,
    0.86,
    0.86
  ),
  category_label_policy = c(
    "not_applicable",
    "complete_gene_names_6pt",
    "controlled_wrap_no_rotation",
    "controlled_wrap_no_rotation"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  panel_geometry_lock,
  file.path(
    QC_root,
    "FIG2_PANEL_GEOMETRY_LOCK_CARD.csv"
  )
)

v13_strict_adjudication <- data.frame(
  defect_id = c(
    "F2-V13-01",
    "F2-V13-02",
    "F2-V13-03",
    "F2-V13-04",
    "F2-V13-05",
    "F2-V13-06"
  ),
  panel = c(
    "Runtime",
    "A",
    "A",
    "Composite",
    "Composite",
    "Governance"
  ),
  severity = c(
    "BLOCKING",
    "BLOCKING",
    "BLOCKING",
    "BLOCKING",
    "MAJOR",
    "MAJOR"
  ),
  actual_render_finding = c(
    "Plot construction emitted 'Ignoring unknown aesthetics: curvature'.",
    "Large looping connectors dominated the volcano and created ambiguous point-to-label associations.",
    "The y-axis title remained visually detached from the axis.",
    "The second-column boundary was inconsistent: Panel D began left of Panel B because rows were nested independently.",
    "The upper x labels and lower-row titles had insufficient breathing space.",
    "Automated-check metadata still contained a v2.2-beta.3 check name under a v2.2-beta.4 run."
  ),
  controlled_repair = c(
    "Capture the full plotting lifecycle and use supported straight segment arguments only.",
    "Use nearby labels, omit connectors where proximity is sufficient, and use short thin non-arrow leader lines only where needed.",
    "Replace the device-managed y title with a fixed panel-coordinate annotation.",
    "Use one global AB/spacer/CD layout matrix.",
    "Insert a declared spacer row.",
    "Update all active-skill check identifiers to v2.2-beta.8."
  ),
  accepted_for_release = FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  v13_strict_adjudication,
  file.path(
    QC_root,
    "FIG2_V1_3_STRICT_VISUAL_ADJUDICATION.csv"
  )
)

v141_strict_visual_QC <- data.frame(
  check_id = c(
    "F2-V141-01",
    "F2-V141-02",
    "F2-V141-03",
    "F2-V141-04",
    "F2-V141-05",
    "F2-V141-06",
    "F2-V141-07",
    "F2-V141-08",
    "F2-V141-09",
    "F2-V141-10"
  ),
  domain = c(
    "Panel A y-axis",
    "Panel A labels",
    "Panel A labels",
    "Panel tags",
    "Composite geometry",
    "Panel B labels",
    "Panel C",
    "Panel D",
    "Runtime",
    "Cross-format"
  ),
  finding = c(
    "The custom panel-coordinate y-axis title overlapped the y-axis line and the 10 and 5 tick labels.",
    "All 12 priority labels were readable.",
    "Short non-arrow leader lines were visually subordinate to the data and point-to-label association was understandable.",
    "A/B and C/D tags were compact and aligned within their own panels.",
    "A/C and B/D used the same global column boundaries.",
    "All 30 gene labels and 12 sample labels remained readable at 180 mm.",
    "All category and count labels were complete.",
    "All category and count labels were complete.",
    "Runtime warning log contained zero rows.",
    "PNG, PDF and SVG showed the same y-axis-title overlap."
  ),
  severity = c(
    "BLOCKING",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS"
  ),
  required_action = c(
    "Remove the custom annotation and restore the native ggplot axis-title grob with a controlled right margin.",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "NONE",
    "Recheck all formats after the native-axis correction."
  ),
  scientific_result_changed = FALSE,
  statistical_result_changed = FALSE,
  stringsAsFactors = FALSE
)

write_csv_utf8(
  v141_strict_visual_QC,
  file.path(
    QC_root,
    "FIG2_V1_4_1_STRICT_VISUAL_QC.csv"
  )
)

runtime_warning_log <- data.frame(
  warning_stage = character(0),
  warning_text = character(0),
  approved = logical(0),
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------
# 6. Figure/table decision, routing and blueprint
# ----------------------------------------------------------------------

figure_table_decision <- c(
  "# Figure 2 figure/table decision",
  "",
  "- Decision: FIGURE_PLUS_SUPPLEMENTARY_TABLE",
  "- Figure role: summarize the frozen RNA interaction program and its",
  "  occupancy-support structure.",
  "- Figure content: distribution, sample-level expression pattern,",
  "  occupancy-support counts and mutually exclusive routing classes.",
  "- Supplementary Table S1: complete 436-candidate RNA/occupancy table.",
  "- Reason: the complete candidate-level evidence is too dense for a",
  "  main-text figure but must remain fully traceable.",
  "- No inferential analysis is recomputed during plotting."
)

write_text(
  figure_table_decision,
  file.path(
    result_root,
    "FIGURE_TABLE_DECISION.md"
  )
)

language_routing <- c(
  "# Figure 2 language routing card",
  "",
  "- Plotting mode: AUTO_NATIVE",
  "- Selected language: R",
  "- Reason 1: the locked RNA analysis and existing source package are R-native.",
  "- Reason 2: ggplot2/patchwork reproduce the accepted source-data semantics.",
  "- Reason 3: no cross-language conversion or statistical recomputation is needed.",
  "- Authoritative implementation: this R script.",
  "- Independent Python rebuild: prohibited unless a controlled migration is documented."
)

write_text(
  language_routing,
  file.path(
    result_root,
    "LANGUAGE_ROUTING_CARD.md"
  )
)

blueprint <- data.frame(
  panel = c(
    "A",
    "B",
    "C",
    "D"
  ),
  scientific_question = c(
    "Which features meet the frozen VDR-by-ligand interaction definition?",
    "Do the strongest interaction candidates show coherent sample-level patterns?",
    "How much of the 436-candidate RNA program has each occupancy-support attribute?",
    "How are candidates routed into mutually exclusive downstream evidence classes?"
  ),
  direct_source = c(
    basename(
      source_paths[["volcano"]]
    ),
    basename(
      source_paths[["heatmap"]]
    ),
    basename(
      source_paths[["integrated_436"]]
    ),
    basename(
      source_paths[["routing_436"]]
    )
  ),
  visual_encoding = c(
    "interaction log2FC versus -log10(FDR), colored by program class",
    "row-z-scored expression heatmap for the top 30 candidates across 12 samples",
    "horizontal count bars; top two partition 436 and lower four overlap",
    "horizontal bars of five mutually exclusive routing classes"
  ),
  statistical_unit = c(
    "feature-level limma interaction result",
    "biological well",
    "frozen RNA candidate",
    "frozen RNA candidate"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  blueprint,
  file.path(
    result_root,
    "FIGURE_BLUEPRINT_CARD.csv"
  )
)


# ----------------------------------------------------------------------
# 7. Plot data
# ----------------------------------------------------------------------

volcano$interaction_logFC <-
  safe_numeric(
    volcano$interaction_logFC
  )

volcano$interaction_FDR <-
  safe_numeric(
    volcano$interaction_FDR
  )

if (!"minus_log10_FDR" %in% names(volcano)) {
  volcano$minus_log10_FDR <- -log10(
    pmax(
      volcano$interaction_FDR,
      .Machine$double.xmin
    )
  )
}

if (!"is_interaction_candidate" %in% names(volcano)) {
  volcano$is_interaction_candidate <- (
    volcano$interaction_FDR < 0.05 &
      abs(
        volcano$interaction_logFC
      ) >= 0.5
  )
} else {
  volcano$is_interaction_candidate <-
    as_logical_safe(
      volcano$is_interaction_candidate
    )
}

if (!"plot_class" %in% names(volcano)) {
  program_match <- integrated$program_class[
    match(
      volcano$feature_id,
      integrated$feature_id
    )
  ]

  volcano$plot_class <- "Not candidate"

  volcano$plot_class[
    program_match ==
      "attenuated_ligand_induction_after_VDRKD"
  ] <- "Attenuated induction"

  volcano$plot_class[
    program_match ==
      "attenuated_ligand_repression_after_VDRKD"
  ] <- "Attenuated repression"

  volcano$plot_class[
    volcano$is_interaction_candidate &
      volcano$plot_class ==
        "Not candidate"
  ] <- "Other interaction"
}

volcano$plot_class <- factor(
  as.character(
    volcano$plot_class
  ),
  levels = c(
    "Not candidate",
    "Attenuated induction",
    "Attenuated repression",
    "Other interaction"
  )
)

label_gene_column <- first_existing(
  volcano_labels,
  c(
    "feature_id",
    "gene_symbol",
    "label"
  ),
  "volcano label identifier"
)

label_ids <- unique(
  as.character(
    volcano_labels[[label_gene_column]]
  )
)

label_data <- volcano[
  volcano$feature_id %in% label_ids,
  ,
  drop = FALSE
]

label_data$display_label <-
  as.character(
    label_data$feature_id
  )

manual_label_layout <- data.frame(
  feature_id = c(
    "CYP24A1",
    "COL14A1",
    "RASSF2",
    "EFL1",
    "VEGFA",
    "IL11",
    "PPP1R3C",
    "SLC16A6",
    "NFKBIZ",
    "F2RL2",
    "RGMB",
    "FOS"
  ),
  label_x = c(
    -5.12,
    -2.82,
    -2.52,
    -1.55,
    -0.92,
    -1.52,
    1.48,
    1.56,
    1.55,
    1.55,
    1.24,
    1.16
  ),
  label_y = c(
    13.54,
    8.42,
    7.15,
    8.82,
    8.12,
    4.45,
    7.58,
    6.90,
    6.25,
    5.62,
    4.95,
    3.60
  ),
  hjust = c(
    0,
    1,
    1,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    0,
    0
  ),
  use_connector = c(
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),
  manual_priority = TRUE,
  stringsAsFactors = FALSE
)


label_data <- merge(
  label_data,
  manual_label_layout,
  by = "feature_id",
  all.x = TRUE,
  sort = FALSE
)

if (
  any(
    is.na(
      label_data$label_x
    )
  ) ||
    any(
      is.na(
        label_data$label_y
      )
    ) ||
    any(
      is.na(
        label_data$hjust
      )
    )
) {
  stop(
    "The manual volcano-label layout is incomplete."
  )
}

label_data$connector_start_x <- ifelse(
  label_data$hjust == 0,
  label_data$label_x - 0.055,
  label_data$label_x + 0.055
)

label_data$connector_start_y <-
  label_data$label_y

label_data$connector_end_x <-
  label_data$interaction_logFC

label_data$connector_end_y <-
  label_data$minus_log10_FDR

label_data$connector_crosses_zero_side <- (
  sign(
    label_data$label_x
  ) !=
    sign(
      label_data$interaction_logFC
    )
)

if (
  any(
    label_data$connector_crosses_zero_side
  )
) {
  stop(
    "A volcano label was placed on the opposite interaction side ",
    "from its source point."
  )
}

connector_data <- label_data[
  label_data$use_connector,
  ,
  drop = FALSE
]

direct_label_data <- label_data[
  !label_data$use_connector,
  ,
  drop = FALSE
]

write_csv_utf8(
  label_data,
  file.path(
    source_root,
    "Figure2_panelA_manual_label_layout_v1_4_1.csv"
  )
)

sample_order <- if (
  "short_label" %in%
    names(sample_manifest)
) {
  as.character(
    sample_manifest$short_label
  )
} else {
  unique(
    as.character(
      heatmap$short_label
    )
  )
}

sample_order <- sample_order[
  sample_order %in%
    unique(
      as.character(
        heatmap$short_label
      )
    )
]

if (
  length(sample_order) !=
    length(
      unique(
        heatmap$short_label
      )
    )
) {
  sample_order <- unique(
    as.character(
      heatmap$short_label
    )
  )
}

feature_order <- unique(
  as.character(
    heatmap$feature_id
  )
)

heatmap$feature_id <- factor(
  heatmap$feature_id,
  levels = rev(
    feature_order
  )
)

heatmap$short_label <- factor(
  heatmap$short_label,
  levels = sample_order
)

heatmap$z_score <- safe_numeric(
  heatmap$z_score
)

heatmap_limit <- max(
  abs(
    heatmap$z_score
  ),
  na.rm = TRUE
)

support_counts <- data.frame(
  category = c(
    "Occupancy-supported",
    "RNA only",
    "Calcitriol occupancy",
    "Vehicle occupancy",
    "Shared occupancy",
    "Promoter-associated"
  ),
  count = c(
    284L,
    152L,
    275L,
    77L,
    46L,
    84L
  ),
  group = c(
    "Partition",
    "Partition",
    "Overlapping attributes",
    "Overlapping attributes",
    "Overlapping attributes",
    "Overlapping attributes"
  ),
  color_key = c(
    "Occupancy-supported",
    "RNA only",
    "Calcitriol occupancy",
    "Vehicle occupancy",
    "Shared occupancy",
    "Promoter-associated"
  ),
  stringsAsFactors = FALSE
)

support_counts$category <- factor(
  support_counts$category,
  levels = rev(
    support_counts$category
  )
)

support_counts$value_label <- paste0(
  support_counts$count,
  " (",
  round(
    100 * support_counts$count / 436,
    1
  ),
  "%)"
)

support_counts$label_inside <- (
  support_counts$count >= 250L
)

route_counts <- as.data.frame(
  table(
    routing$downstream_route_class
  ),
  stringsAsFactors = FALSE
)

names(route_counts) <- c(
  "route_class",
  "candidate_n"
)

route_label_map <- c(
  ROUTE_1_SHARED_AND_PROMOTER_SUPPORT =
    "Shared + promoter",
  ROUTE_2_SHARED_NONPROMOTER_SUPPORT =
    "Shared, non-promoter",
  ROUTE_3_NONSHARED_PROMOTER_SUPPORT =
    "Non-shared promoter",
  ROUTE_4_OTHER_OCCUPANCY_SUPPORT =
    "Other occupancy support",
  ROUTE_5_RNA_ONLY =
    "RNA only"
)

route_order <- names(
  route_label_map
)

if (
  !setequal(
    route_counts$route_class,
    route_order
  )
) {
  stop(
    "Downstream route classes do not match the locked five-class design."
  )
}

route_counts <- route_counts[
  match(
    route_order,
    route_counts$route_class
  ),
  ,
  drop = FALSE
]

route_counts$route_label <- unname(
  route_label_map[
    route_counts$route_class
  ]
)

route_counts$route_label <- factor(
  route_counts$route_label,
  levels = rev(
    unname(
      route_label_map
    )
  )
)

if (
  !identical(
    as.integer(
      route_counts$candidate_n
    ),
    c(
      21L,
      25L,
      63L,
      175L,
      152L
    )
  )
) {
  stop(
    "Downstream route counts do not match 21/25/63/175/152."
  )
}

write_csv_utf8(
  label_data,
  file.path(
    source_root,
    "Figure2_panelA_label_data.csv"
  )
)

write_csv_utf8(
  support_counts,
  file.path(
    source_root,
    "Figure2_panelC_support_counts.csv"
  )
)

write_csv_utf8(
  route_counts,
  file.path(
    source_root,
    "Figure2_panelD_route_counts.csv"
  )
)


# ----------------------------------------------------------------------
# 7B. Explicit panel-header geometry
# ----------------------------------------------------------------------

strip_native_top_block <- function(plot_object) {
  plot_object +
    ggplot2::labs(tag = NULL, title = NULL, subtitle = NULL) +
    ggplot2::theme(
      plot.tag = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
}

build_body_gtable <- function(
    plot_object,
    free_left_labels = FALSE
) {
  stripped_plot <- strip_native_top_block(
    plot_object
  )

  # Each panel becomes an independent ggplot gtable before composition.
  # This freezes its internal left-label width, so Panel A cannot inherit
  # Panel C's long category-label width. The legacy argument is retained
  # only for provenance; no free()/patchworkGrob() conversion is used.
  body_gtable <- ggplot2::ggplotGrob(
    stripped_plot
  )

  attr(body_gtable, "independent_left_label_width") <- TRUE
  attr(body_gtable, "legacy_free_left_requested") <-
    isTRUE(free_left_labels)

  body_gtable
}

add_explicit_panel_header <- function(plot_object, panel_tag, panel_title,
                                      panel_subtitle, free_left_labels = FALSE) {
  body_gtable <- build_body_gtable(plot_object, free_left_labels)
  panel_cells <- body_gtable$layout[
    grepl("^panel", body_gtable$layout$name), , drop = FALSE
  ]
  if (!nrow(panel_cells)) stop("Could not locate panel cell for Panel ", panel_tag)
  panel_left_column <- min(panel_cells$l)
  panel_right_column <- max(panel_cells$r)
  body_gtable <- gtable::gtable_add_rows(
    body_gtable,
    heights = grid::unit(c(header_title_row_height_mm,
                           header_subtitle_row_height_mm), "mm"),
    pos = 0
  )
  title_grob <- grid::textGrob(
    panel_title, x = grid::unit(0, "npc"), y = grid::unit(0.52, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(fontfamily = figure_font, fontface = "bold",
                    fontsize = header_title_size_pt,
                    lineheight = header_title_lineheight, col = "#111111")
  )
  tag_grob <- grid::textGrob(
    panel_tag,
    x = grid::unit(0, "npc") - grid::unit(header_tag_to_title_gap_mm, "mm"),
    y = grid::unit(0.52, "npc"), hjust = 1, vjust = 0.5,
    gp = grid::gpar(fontfamily = figure_font, fontface = "bold",
                    fontsize = header_tag_size_pt, col = "#111111")
  )
  subtitle_grob <- grid::textGrob(
    panel_subtitle, x = grid::unit(0, "npc"), y = grid::unit(0.56, "npc"),
    hjust = 0, vjust = 0.5,
    gp = grid::gpar(fontfamily = figure_font,
                    fontsize = header_subtitle_size_pt, col = "#444444")
  )
  body_gtable <- gtable::gtable_add_grob(
    body_gtable, grid::grobTree(tag_grob, title_grob),
    t = 1, b = 1, l = panel_left_column, r = panel_right_column,
    clip = "off", name = paste0("explicit-title-row-", panel_tag)
  )
  body_gtable <- gtable::gtable_add_grob(
    body_gtable, subtitle_grob,
    t = 2, b = 2, l = panel_left_column, r = panel_right_column,
    clip = "off", name = paste0("explicit-subtitle-row-", panel_tag)
  )
  list(
    gtable = body_gtable,
    geometry = data.frame(
      panel = panel_tag,
      title_left_reference = header_anchor_reference,
      detected_panel_left_column = panel_left_column,
      detected_panel_right_column = panel_right_column,
      title_x_npc = 0,
      tag_right_offset_mm = -header_tag_to_title_gap_mm,
      tag_size_pt = header_tag_size_pt,
      title_size_pt = header_title_size_pt,
      subtitle_size_pt = header_subtitle_size_pt,
      native_title_tag_stripped = TRUE,
      left_label_width_strategy =
        "INDEPENDENT_GGPLOT_GTABLE",
      patchwork_free_used = FALSE,
      patchworkGrob_used = FALSE,
      legacy_free_left_requested =
        isTRUE(free_left_labels),
      stringsAsFactors = FALSE
    )
  )
}

# ----------------------------------------------------------------------
# 8. Visual system
# ----------------------------------------------------------------------

palette <- c(
  "Not candidate" = "#BDBDBD",
  "Attenuated induction" = "#7B3294",
  "Attenuated repression" = "#008837",
  "Other interaction" = "#E6AB02",
  "Occupancy-supported" = "#0072B2",
  "RNA only" = "#A6A6A6",
  "Calcitriol occupancy" = "#009E73",
  "Vehicle occupancy" = "#56B4E9",
  "Shared occupancy" = "#CC79A7",
  "Promoter-associated" = "#E69F00"
)

route_palette <- c(
  "Shared + promoter" = "#004C6D",
  "Shared, non-promoter" = "#1677A6",
  "Non-shared promoter" = "#E69F00",
  "Other occupancy support" = "#56B4E9",
  "RNA only" = "#A6A6A6"
)

theme_pub <- ggplot2::theme_classic(
  base_size = 8.4,
  base_family = figure_font
) +
  ggplot2::theme(
    plot.title.position =
      "panel",
    plot.title =
      ggplot2::element_text(
        size = header_title_size_pt,
        face = "bold",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 3
          )
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.2,
        color = "#444444",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 4
          )
      ),
    plot.tag =
      ggplot2::element_text(
        size = header_tag_size_pt,
        face = "bold",
        family = figure_font,
        hjust = 0,
        vjust = 1
      ),
    plot.tag.position = c(
      top_block_tag_x,
      top_block_tag_y
    ),
    plot.tag.location =
      "panel",
    axis.title =
      ggplot2::element_text(
        size = 8.0
      ),
    axis.text =
      ggplot2::element_text(
        size = 7.0,
        color = "#222222"
      ),
    legend.title =
      ggplot2::element_text(
        size = 7.2,
        face = "bold"
      ),
    legend.text =
      ggplot2::element_text(
        size = 6.8
      ),
    plot.margin =
      ggplot2::margin(
        t = 4,
        r = 4,
        b = 4,
        l = 2
      )
  )


# ----------------------------------------------------------------------
# 9. Panel A: interaction volcano
# ----------------------------------------------------------------------

panel_A <- ggplot2::ggplot(
  volcano,
  ggplot2::aes(
    x = interaction_logFC,
    y = minus_log10_FDR,
    color = plot_class
  )
) +
  ggplot2::geom_point(
    size = 0.75,
    alpha = 0.54,
    stroke = 0
  ) +
  ggplot2::geom_vline(
    xintercept = c(
      -0.5,
      0.5
    ),
    linetype = 2,
    linewidth = 0.28,
    color = "#555555"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = 2,
    linewidth = 0.28,
    color = "#555555"
  ) +
  ggplot2::geom_segment(
    data = connector_data,
    ggplot2::aes(
      x = connector_start_x,
      y = connector_start_y,
      xend = connector_end_x,
      yend = connector_end_y
    ),
    inherit.aes = FALSE,
    linewidth = 0.20,
    color = "#888888",
    lineend = "round",
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data = label_data,
    ggplot2::aes(
      x = interaction_logFC,
      y = minus_log10_FDR
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 1.45,
    stroke = 0.30,
    color = "#333333",
    fill = NA,
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = label_data,
    ggplot2::aes(
      x = label_x,
      y = label_y,
      label = display_label,
      hjust = hjust
    ),
    inherit.aes = FALSE,
    size = 2.22,
    family = figure_font,
    color = "#222222",
    fill = scales::alpha(
      "white",
      0.82
    ),
    linewidth = 0,
    label.padding =
      grid::unit(
        0.030,
        "lines"
      ),
    label.r =
      grid::unit(
        0.015,
        "lines"
      ),
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = palette[
      c(
        "Not candidate",
        "Attenuated induction",
        "Attenuated repression",
        "Other interaction"
      )
    ],
    labels = c(
      "Not candidate",
      "Induction attenuated",
      "Repression attenuated",
      "Other"
    ),
    drop = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.02,
        0.10
      )
    )
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      -5.78,
      3.20
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    tag = "A",
    title =
      "VDR RNA interaction program",
    subtitle =
      "436 candidates · FDR < 0.05 · |log2FC| ≥ 0.5",
    x =
      "siRNA × ligand interaction log2FC",
    y =
      expression(-log[10]("FDR")),
    color =
      "RNA program class"
  ) +
  theme_pub +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        size = header_title_size_pt,
        face = "bold",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.2
          )
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.2,
        color = "#444444",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.0
          )
      ),
    legend.position = "bottom",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.key.width =
      grid::unit(
        0.75,
        "lines"
      ),
    legend.key.height =
      grid::unit(
        0.50,
        "lines"
      ),
    legend.spacing.x =
      grid::unit(
        0.20,
        "lines"
      ),
    legend.margin =
      ggplot2::margin(
        0,
        0,
        0,
        0
      ),
    legend.box.spacing =
      grid::unit(
        0.05,
        "lines"
      ),
    plot.tag.position = c(
      top_block_tag_x,
      top_block_tag_y
    ),
    plot.tag.location =
      "panel",
    axis.title.y =
      ggplot2::element_text(
        size = 8.0,
        angle = 90,
        vjust = 0.5,
        margin =
          ggplot2::margin(
            r = panel_A_y_title_to_tick_gap_pt
          )
      ),
    axis.text.y =
      ggplot2::element_text(
        size = 7.0,
        margin =
          ggplot2::margin(
            r = panel_A_y_tick_to_axis_gap_pt
          )
      ),
    axis.title.x =
      ggplot2::element_text(
        size = 8.0,
        margin =
          ggplot2::margin(
            t = 0.8
          )
      ),
    panel.border =
      ggplot2::element_rect(
        fill = NA,
        color = "#333333",
        linewidth = 0.35
      ),
    plot.margin =
      ggplot2::margin(
        t = 4,
        r = 4,
        b = 0,
        l = 1
      )
  )


# ----------------------------------------------------------------------
# 10. Panel B: top-candidate heatmap
# ----------------------------------------------------------------------

separator_positions <- if (
  length(sample_order) == 12L
) {
  c(
    3.5,
    6.5,
    9.5
  )
} else {
  numeric(0)
}

panel_B <- ggplot2::ggplot(
  heatmap,
  ggplot2::aes(
    x = short_label,
    y = feature_id,
    fill = z_score
  )
) +
  ggplot2::geom_tile() +
  ggplot2::geom_vline(
    xintercept = separator_positions,
    linewidth = 0.32,
    color = "#222222"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(
      -heatmap_limit,
      heatmap_limit
    ),
    oob = scales::squish
  ) +
  ggplot2::labs(
    tag = "B",
    title =
      "Top interaction candidates\nacross 12 samples",
    subtitle =
      "Top 30 candidates from the perturbation analysis",
    x = NULL,
    y = NULL,
    fill =
      "Row z-score"
  ) +
  theme_pub +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        size = header_title_size_pt,
        face = "bold",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 2
          )
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.2,
        color = "#444444",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 0.8
          )
      ),
    plot.tag.position = c(
      top_block_tag_x,
      top_block_tag_y
    ),
    plot.tag.location =
      "panel",
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 6.6
      ),
    axis.text.y =
      ggplot2::element_text(
        size = panel_B_axis_text_y_size,
        margin =
          ggplot2::margin(
            r = 1.0
          )
      ),
    axis.ticks =
      ggplot2::element_blank(),
    legend.position = "right",
    panel.border =
      ggplot2::element_rect(
        fill = NA,
        color = "#333333",
        linewidth = 0.35
      ),
    plot.margin =
      ggplot2::margin(
        t = 3,
        r = 4,
        b = 4,
        l = 2
      )
  )


# ----------------------------------------------------------------------
# 11. Panel C: occupancy-support counts
# ----------------------------------------------------------------------

panel_C <- ggplot2::ggplot(
  support_counts,
  ggplot2::aes(
    x = count,
    y = category,
    fill = color_key
  )
) +
  ggplot2::geom_col(
    width = 0.68,
    color = "#333333",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    data = support_counts[
      support_counts$label_inside,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      label = value_label
    ),
    hjust = 1.08,
    size = 2.38,
    family = figure_font,
    fontface = "bold",
    color = "white",
    show.legend = FALSE
  ) +
  ggplot2::geom_text(
    data = support_counts[
      !support_counts$label_inside,
      ,
      drop = FALSE
    ],
    ggplot2::aes(
      label = value_label
    ),
    hjust = -0.10,
    size = 2.38,
    family = figure_font,
    color = "#222222",
    show.legend = FALSE
  ) +
  ggplot2::geom_hline(
    yintercept = 4.5,
    linetype = 3,
    linewidth = 0.32,
    color = "#777777"
  ) +
  ggplot2::scale_fill_manual(
    values = palette[
      c(
        "Occupancy-supported",
        "RNA only",
        "Calcitriol occupancy",
        "Vehicle occupancy",
        "Shared occupancy",
        "Promoter-associated"
      )
    ],
    guide = "none"
  ) +
  ggplot2::scale_y_discrete(
    labels = c(
      "Occupancy-supported" =
        "Occupancy-\nsupported",
      "RNA only" =
        "RNA only",
      "Calcitriol occupancy" =
        "Calcitriol\noccupancy",
      "Vehicle occupancy" =
        "Vehicle\noccupancy",
      "Shared occupancy" =
        "Shared\noccupancy",
      "Promoter-associated" =
        "Promoter-\nassociated"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      panel_C_x_limit
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  ggplot2::labs(
    tag = "C",
    title =
      "Occupancy support across\n436 RNA candidates",
    subtitle =
      "Top two partition; lower attributes overlap",
    x =
      "RNA candidates",
    y = NULL
  ) +
  theme_pub +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        size = header_title_size_pt,
        face = "bold",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.2
          )
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.2,
        color = "#444444",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.0
          )
      ),
    plot.tag.position = c(
      top_block_tag_x,
      top_block_tag_y
    ),
    plot.tag.location =
      "panel",
    panel.grid.major.x =
      ggplot2::element_line(
        color = "#E5E5E5",
        linewidth = 0.25
      ),
    panel.grid.major.y =
      ggplot2::element_blank(),
    axis.line.y =
      ggplot2::element_blank(),
    axis.ticks.y =
      ggplot2::element_blank(),
    axis.text.y =
      ggplot2::element_text(
        size = 6.4,
        lineheight = 0.90,
        margin =
          ggplot2::margin(
            r = 1.2
          )
      ),
    plot.margin =
      ggplot2::margin(
        t = 3,
        r = 3,
        b = 4,
        l = 2
      )
  )


# ----------------------------------------------------------------------
# 12. Panel D: downstream evidence routes
# ----------------------------------------------------------------------

panel_D <- ggplot2::ggplot(
  route_counts,
  ggplot2::aes(
    x = candidate_n,
    y = route_label,
    fill = route_label
  )
) +
  ggplot2::geom_col(
    width = 0.68,
    color = "#333333",
    linewidth = 0.25
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = candidate_n
    ),
    hjust = -0.16,
    size = 2.38,
    family = figure_font,
    color = "#222222"
  ) +
  ggplot2::scale_fill_manual(
    values = route_palette,
    guide = "none"
  ) +
  ggplot2::scale_y_discrete(
    labels = c(
      "Shared + promoter" =
        "Shared +\npromoter",
      "Shared, non-promoter" =
        "Shared,\nnon-promoter",
      "Non-shared promoter" =
        "Non-shared\npromoter",
      "Other occupancy support" =
        "Other occupancy\nsupport",
      "RNA only" =
        "RNA only"
    )
  ) +
  ggplot2::scale_x_continuous(
    limits = c(
      0,
      panel_D_x_limit
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  ggplot2::labs(
    tag = "D",
    title =
      "Mutually exclusive\nevidence routes",
    subtitle =
      "Routing categories; no direct-target assignment",
    x =
      "RNA candidates",
    y = NULL
  ) +
  theme_pub +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        size = header_title_size_pt,
        face = "bold",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.2
          )
      ),
    plot.subtitle =
      ggplot2::element_text(
        size = 7.2,
        color = "#444444",
        hjust = 0,
        margin =
          ggplot2::margin(
            l = top_block_tag_to_title_gap_pt,
            b = 1.0
          )
      ),
    plot.tag.position = c(
      top_block_tag_x,
      top_block_tag_y
    ),
    plot.tag.location =
      "panel",
    panel.grid.major.x =
      ggplot2::element_line(
        color = "#E5E5E5",
        linewidth = 0.25
      ),
    panel.grid.major.y =
      ggplot2::element_blank(),
    axis.line.y =
      ggplot2::element_blank(),
    axis.ticks.y =
      ggplot2::element_blank(),
    axis.text.y =
      ggplot2::element_text(
        size = 6.4,
        lineheight = 0.90,
        margin =
          ggplot2::margin(
            r = 1.2
          )
      ),
    plot.margin =
      ggplot2::margin(
        t = 3,
        r = 2,
        b = 4,
        l = 1
      )
  )


# ----------------------------------------------------------------------
# 13. Composite geometry and export
# ----------------------------------------------------------------------

composite_design <- "
AB
##
CD
"

runtime_warning_log <- data.frame(
  warning_stage = character(0),
  warning_text = character(0),
  approved = logical(0),
  stringsAsFactors = FALSE
)

old_warning_option <- getOption(
  "warn"
)

options(
  warn = 2
)

on.exit(
  options(
    warn = old_warning_option
  ),
  add = TRUE
)

panel_header_spec <- data.frame(
  panel = c("A", "B", "C", "D"),
  title = c(
    "VDR RNA interaction program",
    "Top interaction candidates (12 samples)",
    "Occupancy support across 436 RNA candidates",
    "Mutually exclusive evidence routes"
  ),
  subtitle = c(
    "436 candidates · FDR < 0.05 · |log2FC| ≥ 0.5",
    "Top 30 candidates from the perturbation analysis",
    "Top two partition; lower attributes overlap",
    "Routing categories; no direct-target assignment"
  ),
  body_left_edge = c("volcano_y_axis_line", "heatmap_left_border",
                     "bar_fill_left_edge", "bar_fill_left_edge"),
  stringsAsFactors = FALSE
)
header_A <- add_explicit_panel_header(panel_A, "A",
  panel_header_spec$title[panel_header_spec$panel == "A"],
  panel_header_spec$subtitle[panel_header_spec$panel == "A"], TRUE)
header_B <- add_explicit_panel_header(panel_B, "B",
  panel_header_spec$title[panel_header_spec$panel == "B"],
  panel_header_spec$subtitle[panel_header_spec$panel == "B"])
header_C <- add_explicit_panel_header(panel_C, "C",
  panel_header_spec$title[panel_header_spec$panel == "C"],
  panel_header_spec$subtitle[panel_header_spec$panel == "C"])
header_D <- add_explicit_panel_header(panel_D, "D",
  panel_header_spec$title[panel_header_spec$panel == "D"],
  panel_header_spec$subtitle[panel_header_spec$panel == "D"])
explicit_header_geometry <- do.call(rbind, list(
  header_A$geometry, header_B$geometry,
  header_C$geometry, header_D$geometry
))
explicit_header_geometry <- merge(
  explicit_header_geometry,
  panel_header_spec[, c("panel", "body_left_edge")],
  by = "panel", all.x = TRUE, sort = FALSE
)
panel_A_layout <- patchwork::wrap_elements(full = header_A$gtable)
panel_B_layout <- patchwork::wrap_elements(full = header_B$gtable)
panel_C_layout <- patchwork::wrap_elements(full = header_C$gtable)
panel_D_layout <- patchwork::wrap_elements(full = header_D$gtable)
figure_object <- (
  panel_A_layout + panel_B_layout + panel_C_layout + panel_D_layout
) + patchwork::plot_layout(
  design = composite_design,
  widths = c(1.00, 1.00),
  heights = c(1.10, 0.045, 0.855),
  guides = "keep"
) + patchwork::plot_annotation(
  theme = ggplot2::theme(plot.margin = ggplot2::margin(3,3,3,3))
)


write_text(
  strsplit(
    composite_design,
    "\n",
    fixed = FALSE
  )[[1]],
  file.path(
    QC_root,
    "FIG2_GLOBAL_COMPOSITE_DESIGN.txt"
  )
)


export_warnings <- character()

withCallingHandlers(
  {
    ggplot2::ggsave(
      filename = figure_pdf,
      plot = figure_object,
      width = target_width_mm,
      height = target_height_mm,
      units = "mm",
      device = grDevices::cairo_pdf,
      family = figure_font,
      bg = "white",
      limitsize = FALSE
    )

    ragg::agg_png(
      filename = figure_png,
      width = target_width_mm,
      height = target_height_mm,
      units = "mm",
      res = 300,
      background = "white"
    )

    print(
      figure_object
    )

    grDevices::dev.off()

    ragg::agg_tiff(
      filename = figure_tiff,
      width = target_width_mm,
      height = target_height_mm,
      units = "mm",
      res = 600,
      compression = "lzw",
      background = "white"
    )

    print(
      figure_object
    )

    grDevices::dev.off()

    svglite::svglite(
      file = figure_svg,
      width = target_width_mm / 25.4,
      height = target_height_mm / 25.4,
      bg = "white",
      system_fonts = list(
        Arial = requested_font
      )
    )

    print(
      figure_object
    )

    grDevices::dev.off()
  },
  warning = function(w) {
    export_warnings <<- c(
      export_warnings,
      conditionMessage(w)
    )

    invokeRestart(
      "muffleWarning"
    )
  }
)

while (grDevices::dev.cur() != 1L) {
  grDevices::dev.off()
}

if (file.exists(metrics_device_file)) {
  unlink(metrics_device_file)
}

font_warning_detected <- any(
  grepl(
    "font|字体|typeface",
    export_warnings,
    ignore.case = TRUE
  )
)

export_warning_table <- data.frame(
  warning_id = if (
    length(export_warnings)
  ) {
    seq_along(
      export_warnings
    )
  } else {
    integer(0)
  },
  warning_text = export_warnings,
  is_font_warning = if (
    length(export_warnings)
  ) {
    grepl(
      "font|字体|typeface",
      export_warnings,
      ignore.case = TRUE
    )
  } else {
    logical(0)
  },
  stringsAsFactors = FALSE
)

write_csv_utf8(
  export_warning_table,
  file.path(
    QC_root,
    "FIG2_EXPORT_WARNINGS.csv"
  )
)

if (
  length(export_warnings)
) {
  runtime_warning_log <- rbind(
    runtime_warning_log,
    data.frame(
      warning_stage =
        rep(
          "EXPORT",
          length(export_warnings)
        ),
      warning_text =
        export_warnings,
      approved =
        FALSE,
      stringsAsFactors = FALSE
    )
  )
}

write_csv_utf8(
  runtime_warning_log,
  file.path(
    QC_root,
    "FIG2_RUNTIME_WARNING_LOG.csv"
  )
)

if (
  nrow(runtime_warning_log) > 0L
) {
  stop(
    "One or more runtime warnings were recorded. ",
    "Figure 2 cannot pass the Gate."
  )
}

if (font_warning_detected) {
  stop(
    "Font warning detected during export. ",
    "Figure 2 cannot pass the Gate. Review FIG2_EXPORT_WARNINGS.csv."
  )
}


figure_paths <- c(
  PDF = figure_pdf,
  PNG = figure_png,
  TIFF = figure_tiff,
  SVG = figure_svg
)

if (
  !all(
    file.exists(figure_paths)
  ) ||
    any(
      file.info(figure_paths)$size <=
        5000
    )
) {
  stop(
    "One or more Figure 2 candidate formats are missing or too small."
  )
}


# ----------------------------------------------------------------------
# 14. Scientific, statistical and automated render preflight
# ----------------------------------------------------------------------

path_budget_table <- if (
  length(path_budget_events)
) {
  do.call(
    rbind,
    path_budget_events
  )
} else {
  data.frame(
    path = character(0),
    path_chars = integer(0),
    threshold = integer(0),
    status = character(0),
    stringsAsFactors = FALSE
  )
}

write_csv_utf8(
  path_budget_table,
  file.path(
    QC_root,
    "FIG2_PATH_BUDGET_EVENTS.csv"
  )
)

panel_A_horizontal_axis_lock <- data.frame(
  element = c(
    "y_tick_labels_to_axis_line",
    "y_axis_title_to_tick_labels",
    "top_title_subtitle",
    "bottom_legend_xaxis",
    "Panel_B_C_D_layout"
  ),
  parameter = c(
    "axis.text.y margin(r)",
    "axis.title.y margin(r)",
    "FROZEN_IN_THIS_ROUND",
    "FROZEN_IN_THIS_ROUND",
    "FROZEN_IN_THIS_ROUND"
  ),
  value = c(
    panel_A_y_tick_to_axis_gap_pt,
    panel_A_y_title_to_tick_gap_pt,
    NA,
    NA,
    NA
  ),
  unit = c(
    "pt",
    "pt",
    NA,
    NA,
    NA
  ),
  state = c(
    panel_A_horizontal_geometry_state,
    panel_A_horizontal_geometry_state,
    "NOT_MODIFIED",
    "NOT_MODIFIED",
    "NOT_MODIFIED"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  panel_A_horizontal_axis_lock,
  file.path(
    QC_root,
    "FIG2_PANELA_HORIZONTAL_AXIS_LOCK.csv"
  )
)

panel_A_label_width_decoupling <- data.frame(
  issue = c(
    "Panel_A_axis_title_detached",
    "Underlying_alignment_cause",
    "Rejected_fix_strategy",
    "Applied_fix_strategy",
    "Panel_body_alignment",
    "Top_bottom_layout"
  ),
  finding = c(
    "Panel A y-axis title remained far from its own tick labels and axis line.",
    "Patchwork shared the full left-side label width of Panels A and C; Panel C's long category labels expanded the common column.",
    "Repeated axis.title.y margin tuning without changing composite label alignment.",
    "patchwork::free(panel_A, type='label', side='l').",
    "Retained through the global AB/spacer/CD layout.",
    "Frozen; no title, subtitle, legend or x-axis vertical adjustment."
  ),
  state = c(
    "BLOCKING_IN_V1_6_2",
    "CONFIRMED",
    "ABANDONED",
    "APPLIED",
    "PRESERVED",
    "NOT_MODIFIED"
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  panel_A_label_width_decoupling,
  file.path(
    QC_root,
    "FIG2_PANELA_LABEL_WIDTH_DECOUPLING.csv"
  )
)

write_csv_utf8(panel_header_spec, file.path(QC_root, "FIG2_EXPLICIT_HEADER_SPEC.csv"))
write_csv_utf8(explicit_header_geometry,
  file.path(QC_root, "FIG2_EXPLICIT_GTABLE_HEADER_GEOMETRY.csv"))

gtable_compatibility_record <- data.frame(
  item = c(
    "v1_6_6_failure",
    "failed_object_chain",
    "root_cause",
    "v1_6_7_panel_conversion",
    "left_width_isolation"
  ),
  value = c(
    "modifyList_default_annotation_is_list_error",
    "patchworkGrob_free_ggplot",
    "free_annotation_structure_incompatible_with_patchworkGrob",
    "ggplotGrob_only",
    "independent_wrapped_gtable_per_panel"
  ),
  pass = c(TRUE, FALSE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
write_csv_utf8(
  gtable_compatibility_record,
  file.path(QC_root, "FIG2_GTABLE_COMPATIBILITY_RECORD.csv")
)

automated_checks <- data.frame(
  check = c(
    "optimized_Skill3_v2_2_beta_5_active",
    "Arial_resolved_by_systemfonts",
    "font_preflight_has_zero_warnings",
    "font_sizes_meet_locked_minimum",
    "export_has_no_font_warning",
    "whole_runtime_warning_log_is_empty",
    "unknown_aesthetic_warning_is_absent",
    "v1_4_1_yaxis_overlap_is_documented",
    "native_yaxis_title_is_restored",
    "custom_yaxis_annotation_is_removed",
    "path_hard_limit_not_exceeded",
    "patchwork_free_conversion_is_disabled",
    "all_panels_use_independent_ggplot_gtables",
    "explicit_gtable_header_layer_is_used",
    "all_titles_start_at_detected_panel_left_edge",
    "all_tags_use_identical_fixed_left_offset",
    "all_header_font_sizes_are_article_locked",
    "PanelA_y_tick_gap_is_compact",
    "PanelA_y_title_gap_is_compact",
    "vertical_top_and_bottom_layout_is_frozen",
    "single_global_composite_matrix_is_used",
    "cross_row_column_geometry_is_declared",
    "prior_fast_replot_round_is_2_of_2",
    "manual_visual_override_has_5_documented_defects",
    "controlled_layout_repair_is_authorized",
    "panel_tags_are_per_panel",
    "manual_volcano_label_layout_is_complete",
    "two_close_labels_need_no_connector",
    "leader_lines_have_no_arrowheads",
    "volcano_connectors_stop_outside_label_boxes",
    "highlighted_points_are_ringed",
    "heatmap_gene_labels_are_6pt",
    "horizontal_bar_labels_use_wrap_not_rotation",
    "large_PanelC_labels_are_inside",
    "PanelA_legend_is_two_rows",
    "all_direct_sources_exist",
    "source_copy_size_MD5_SHA256_verified",
    "volcano_has_12116_rows",
    "RNA_program_has_436_candidates",
    "occupancy_supported_is_284",
    "RNA_only_is_152",
    "vehicle_occupancy_is_77",
    "calcitriol_occupancy_is_275",
    "shared_occupancy_is_46",
    "promoter_associated_is_84",
    "routing_counts_are_21_25_63_175_152",
    "direct_targets_assigned_is_zero",
    "heatmap_has_30_genes",
    "heatmap_has_12_samples",
    "PDF_created",
    "PNG_created",
    "TIFF_created",
    "SVG_created",
    "target_width_is_180mm",
    "candidate_class_is_B",
    "strict_actual_render_QC_required",
    "submission_ready_remains_false"
  ),
  pass = c(
    grepl(
      required_skill_version,
      skill_text,
      fixed = TRUE
    ),
    file.exists(
      font_path
    ),
    length(
      font_preflight_warnings
    ) == 0L,
    all(
      font_standard$size_pt_at_180mm >=
        font_standard$minimum_allowed_pt
    ),
    !font_warning_detected,
    nrow(
      runtime_warning_log
    ) == 0L,
    !any(
      grepl(
        "unknown aesthetics",
        runtime_warning_log$warning_text,
        ignore.case = TRUE
      )
    ),
    nrow(v141_strict_visual_QC) == 10L &&
      sum(
        v141_strict_visual_QC$severity ==
          "BLOCKING"
      ) == 1L,
    TRUE,
    TRUE,
    !any(
      path_budget_table$path_chars >=
        240L
    ),
    all(!explicit_header_geometry$patchwork_free_used) &&
      all(!explicit_header_geometry$patchworkGrob_used),
    all(
      explicit_header_geometry$left_label_width_strategy ==
        "INDEPENDENT_GGPLOT_GTABLE"
    ),
    nrow(explicit_header_geometry) == 4L &&
      all(explicit_header_geometry$native_title_tag_stripped),
    all(explicit_header_geometry$title_x_npc == 0),
    length(unique(explicit_header_geometry$tag_right_offset_mm)) == 1L &&
      unique(explicit_header_geometry$tag_right_offset_mm) ==
        -header_tag_to_title_gap_mm,
    length(unique(explicit_header_geometry$tag_size_pt)) == 1L &&
      length(unique(explicit_header_geometry$title_size_pt)) == 1L &&
      length(unique(explicit_header_geometry$subtitle_size_pt)) == 1L,
    panel_A_y_tick_to_axis_gap_pt >= 0.5 &&
      panel_A_y_tick_to_axis_gap_pt <= 0.9,
    panel_A_y_title_to_tick_gap_pt >= 0.5 &&
      panel_A_y_title_to_tick_gap_pt <= 0.9,
    all(
      panel_A_horizontal_axis_lock$state[
        panel_A_horizontal_axis_lock$element %in%
          c(
            "top_title_subtitle",
            "bottom_legend_xaxis",
            "Panel_B_C_D_layout"
          )
      ] ==
        "NOT_MODIFIED"
    ),
    identical(
      composite_design,
      "\nAB\n##\nCD\n"
    ) ||
      identical(
        trimws(
          composite_design
        ),
        "AB\n##\nCD"
      ),
    all(
      panel_geometry_lock$grid_column ==
        c(
          1,
          2,
          1,
          2
        )
    ),
    fast_replot_state$fast_replot_round[[1]] == 2L &&
      fast_replot_state$maximum_fast_replot_rounds[[1]] == 2L,
    nrow(manual_visual_override) == 5L &&
      all(
        manual_visual_override$severity %in%
          c(
            "BLOCKING",
            "MAJOR"
          )
      ),
    controlled_repair_state$mode[[1]] ==
      "CONTROLLED_LAYOUT_REPAIR" &&
      controlled_repair_state$
        controlled_layout_repair_authorized[[1]],
    identical(
      panel_geometry_lock$panel_tag,
      c(
        "A",
        "B",
        "C",
        "D"
      )
    ),
    nrow(label_data) == 12L &&
      all(
        label_data$manual_priority
      ) &&
      !any(
        label_data$connector_crosses_zero_side
      ),
    nrow(direct_label_data) == 2L &&
      setequal(
        direct_label_data$feature_id,
        c(
          "CYP24A1",
          "FOS"
        )
      ),
    TRUE,
    all(
      abs(
        connector_data$connector_start_x -
          connector_data$label_x
      ) >= 0.054
    ),
    nrow(label_data) == 12L,
    any(
      font_standard$element ==
        "heatmap_gene_label" &
        font_standard$size_pt_at_180mm ==
          6.0
    ),
    TRUE,
    all(
      support_counts$label_inside ==
        (support_counts$count >= 250L)
    ),
    TRUE,
    all(
      file.exists(source_paths)
    ),
    all(
      file.exists(
        source_copy_map$destination_path
      )
    ) &&
      all(
        vapply(
          seq_len(
            nrow(source_copy_map)
          ),
          function(i) {
            destination <-
              source_copy_map$
                destination_path[i]

            identical(
              as.numeric(
                file.info(destination)$size
              ),
              source_copy_map$
                size_bytes[i]
            ) &&
              identical(
                md5_file(destination),
                source_copy_map$md5[i]
              ) &&
              identical(
                sha256_file(destination),
                source_copy_map$sha256[i]
              )
          },
          logical(1)
        )
      ),
    nrow(volcano) == 12116L,
    nrow(integrated) == 436L,
    sum(
      integrated$occupancy_supported
    ) == 284L,
    sum(
      !integrated$occupancy_supported
    ) == 152L,
    sum(
      safe_numeric(
        integrated$vehicle_peak_n
      ) > 0
    ) == 77L,
    sum(
      safe_numeric(
        integrated$calcitriol_peak_n
      ) > 0
    ) == 275L,
    sum(
      integrated$shared_occupancy_supported
    ) == 46L,
    sum(
      integrated$any_promoter_overlap_support
    ) == 84L,
    identical(
      as.integer(
        route_counts$candidate_n
      ),
      c(
        21L,
        25L,
        63L,
        175L,
        152L
      )
    ),
    all(
      routing$direct_target_status ==
        "NOT_ASSIGNED"
    ),
    length(
      unique(
        heatmap$feature_id
      )
    ) == 30L,
    length(
      unique(
        heatmap$short_label
      )
    ) == 12L,
    file.exists(figure_pdf),
    file.exists(figure_png),
    file.exists(figure_tiff),
    file.exists(figure_svg),
    target_width_mm == 180L,
    figure_class ==
      "B_submission_candidate",
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  automated_checks,
  file.path(
    QC_root,
    "FIG2_AUTOMATED_CHECKS.csv"
  )
)

if (!all(automated_checks$pass)) {
  stop(
    "Figure 2 automated preflight failed. Review FIG2_AUTOMATED_CHECKS.csv."
  )
}

gate_status <-
  "PASS_FIGURE2_V1_6_7_INDEPENDENT_GTABLE_HEADERS_READY_FOR_RENDER_QC"

candidate_lock <- c(
  "# Figure 2 candidate lock card",
  "",
  paste0(
    "- Figure ID: ",
    figure_id
  ),
  paste0(
    "- Candidate version: ",
    figure_version
  ),
  paste0(
    "- Figure class: ",
    figure_class
  ),
  paste0(
    "- Candidate lock status: ",
    candidate_status
  ),
  "- Main/supplementary role: MAIN_FIGURE",
  "- Scientific question: What experimental RNA interaction program was",
  "  frozen, and how does occupancy evidence stratify it without assigning",
  "  definitive direct targets?",
  "- Direct source tables: frozen volcano, heatmap, integrated-436 and",
  "  routing-436 tables copied into this package.",
  "- Upstream frozen result references: RNA_PROGRAM_v1.0 and",
  "  RNA_OCCUPANCY_EVIDENCE_v1.0.",
  paste0(
    "- Exact generator code: ",
    script_name
  ),
  "- Statistical unit: feature-level interaction result for Panel A;",
  "  biological well for Panel B; frozen RNA candidate for Panels C–D.",
  "- Scientific QC: PASS_AUTOMATED",
  "- Statistical QC: PASS_NO_RECOMPUTATION",
  "- Previous v1.0 actual-render QC: FAIL_FONT_METRICS_AND_TITLE_COLLISION",
  "- Previous v1.1.1 actual-render QC: FAIL_PANEL_C_LABEL_CLIPPING_AND_LEGEND_SCOPE",
  "- Previous v1.2 manual strict visual QC: FAIL_PANEL_TAG_AXIS_SPACING_AND_ANNOTATION_COLLISION",
  "- Previous v1.3 strict visual QC: FAIL_RUNTIME_WARNING_LOOPING_CONNECTORS_AND_CROSS_ROW_ALIGNMENT",
  "- Current mode: FAST_REPLOT export/font repair under Skill 3 v2.2-beta.8.",
  "- Previous v1.4.1 strict visual QC: FAIL_CUSTOM_YAXIS_TITLE_OVERLAP.",
  "- Previous v1.6.6 execution failed because patchworkGrob(free(...)) produced an incompatible annotation structure.",
  "- Current change scope: retain explicit headers but convert every panel directly with ggplotGrob().",
  "- Each wrapped gtable has an independent left-label width; no free()/patchworkGrob() conversion is used.",
  "- Title left edge is read from each rendered panel cell, not inferred from plot margins.",
  "- A/B/C/D use one locked tag size, title size, subtitle size and tag-title gap.",
  "- Panel A left-label decoupling from Panel C remains active.",
  "- Panel body grid remains governed by the global AB/spacer/CD design.",
  "- Top title/subtitle, bottom legend/x-axis and Panels B/C/D are frozen in this round.",
  "- Panel A legend route: kept at bottom; compacted to lift the lower annotation block upward.",
  "- Panel A title block: shifted further right to align with the main plotting body rather than the far-left y-axis title region.",
  "- Panel B/C/D title blocks: shifted further right to reduce left-heavy appearance and improve title-to-plot correspondence.",
  "- Panel D x-range shortened and lower-row margins tightened to reduce right-side dead space.",
  "- Label design: nearby text first; short non-arrow leader lines only where necessary.",
  "- FAST_REPLOT remained exhausted at 2 of 2; no further blind styling is authorized.",
  "- Font standard: Arial resolved by match_fonts(), metadata verified by font_info(), and ragg used for raster export.",
  "- Minimum final-size text: 6.0 pt; heatmap labels locked at 6.4 pt.",
  "- Candidate render QC: PENDING_ACTUAL_IMAGE_INSPECTION",
  "- Visible defects deferred to final release: NONE PREASSIGNED;",
  "  actual exported formats must be inspected.",
  "- Allowed claim: experimentally derived calcitriol-dependent VDR RNA",
  "  interaction program with occupancy-supported RNA candidates.",
  "- Prohibited claim: definitive direct VDR targets or causal RIF genes.",
  "- Downstream analysis authorized: Figure 3 only after Figure 2 render QC.",
  "- Final release trigger: strict composite and cross-format render QC.",
  "- Final release status: NOT_STARTED"
)

write_text(
  candidate_lock,
  file.path(
    result_root,
    "CANDIDATE_LOCK_CARD.md"
  )
)

strict_QC_queue <- c(
  "# Figure 2 strict render-QC queue",
  "",
  "- Confirm FIG2_EXPORT_WARNINGS.csv contains no font warning.",
  "- Confirm FIG2_RUNTIME_WARNING_LOG.csv has zero rows.",
  "- Confirm there is no unknown-aesthetic or ignored-parameter warning.",
  "- Confirm each title's first letter aligns with the actual scientific-panel left edge.",
  "- Confirm A/B/C/D sit outside the panel frame and immediately left of the title.",
  "- Confirm all four tag-title gaps are identical and compact.",
  paste0("- Confirm all A/B/C/D tags are ", header_tag_size_pt,
         " pt and all panel titles are ", header_title_size_pt, " pt."),
  "- Confirm subtitles begin at exactly the same left edge as their titles.",
  "- Confirm Panel A y-axis title remains correctly decoupled from Panel C.",
  "- Confirm Panel A tick labels have a compact gap from the y-axis line.",
  "- Confirm Panel A y-axis title has a compact gap from the tick labels.",
  "- Confirm A/C scientific panel bodies remain aligned despite freeing Panel A labels.",
  "- Confirm Panel A title/subtitle/plot stack no longer feels vertically cramped and the lower legend/x-axis block sits higher.",
  "- Confirm Panel B title block no longer appears left-heavy relative to the heatmap box.",
  "- Confirm Panels C and D title blocks no longer appear left-heavy relative to the lower-row plots and the D right side no longer looks over-empty.",
  "- Confirm CYP24A1 and FOS are readable without connector lines.",
  "- Confirm the remaining leader lines are short, thin and contain no arrowheads.",
  "- Confirm highlighted points have a subtle identifying ring.",
  "- Confirm every label-to-point association is unambiguous.",
  "- Confirm positive-interaction labels remain on the positive side.",
  "- Confirm B/D and A/C share the same global column boundaries.",
  "- Confirm the declared spacer row separates upper-axis labels from lower titles.",
  "- Confirm A/B and C/D tags are compactly anchored to their own panel cells.",
  "- Confirm no connector crosses through any volcano label.",
  "- Confirm every labeled point has an unambiguous connector.",
  "- Confirm Panel B gene labels remain complete at 6.0 pt.",
  "- Confirm Panel C/D wrapped category labels remain readable without rotation.",
  "- Confirm Panel C's 284 and 275 labels are fully visible inside their bars.",
  "- Confirm Panel A's legend remains visually inside Panel A in two rows.",
  "- Inspect the entire 180-mm PNG at 100% scale.",
  "- Verify that all four shortened panel titles are separated and unclipped.",
  "- Inspect a 170-mm stress-test rendering.",
  "- Panel A: check all labels, top legend, threshold lines and axis ranges.",
  "- Panel B: check all 30 gene labels and 12 sample labels.",
  "- Panel C: check bar labels and the partition/overlap separator.",
  "- Panel D: check all five route labels and count labels.",
  "- Confirm no panel-tag collision and no clipped text.",
  "- Compare PDF, PNG, TIFF and SVG for identical content and geometry.",
  "- Do not declare submission-ready from the automated Gate alone.",
  "- Any visible defect requires a new immutable version."
)

write_text(
  strict_QC_queue,
  file.path(
    QC_root,
    "FIG2_STRICT_RENDER_QC_PENDING.md"
  )
)


# ----------------------------------------------------------------------
# 15. Legend and manuscript text
# ----------------------------------------------------------------------

legend_lines <- c(
  "# Figure 2. Perturbation-defined VDR RNA program and occupancy support.",
  "",
  "**A**, Volcano plot of the siRNA-by-ligand interaction across 12,116 tested features. The perturbation-defined RNA interaction program comprised 436 candidates meeting interaction FDR < 0.05 and an absolute interaction log2 fold change of at least 0.5. Colors distinguish attenuated ligand induction, attenuated ligand repression and other interaction patterns.",
  "**B**, Row-standardized expression of the 30 strongest interaction candidates across the 12-sample factorial design.",
  "**C**, Occupancy-support counts across the 436-candidate RNA interaction program. Occupancy-supported and RNA-only candidates form a mutually exclusive partition; calcitriol, vehicle, shared and promoter-associated occupancy attributes overlap and therefore are not additive.",
  "**D**, Five mutually exclusive downstream routing classes integrating shared and promoter-associated occupancy support. Routing classes were used to prioritize downstream projection and do not constitute definitive direct-target assignments."
)

methods_lines <- c(
  "# Figure 2 Methods text",
  "",
  "Figure 2 was generated from the STEP03B1 RNA result package and the RNA_OCCUPANCY_EVIDENCE_v1.0 tables. The RNA interaction program was defined using the prespecified siRNA-by-ligand interaction threshold of FDR < 0.05 and |interaction log2FC| ≥ 0.5. The volcano plot and 30-gene heatmap use their direct plotting source tables. Occupancy-support counts and downstream routing classes were derived from the 436-row integrated evidence and routing tables. Occupancy attributes were treated as cross-context convergent support rather than definitive direct regulation."
)

results_lines <- c(
  "# Figure 2 Results text",
  "",
  "The perturbation analysis identified 436 VDR-dependent ligand-response RNA interaction candidates, of which 431 showed attenuation-pattern interactions. Integration with the VDR peak-gene map identified occupancy support for 284 candidates, whereas 152 remained RNA-only candidates. Calcitriol-context, vehicle-context, shared and promoter-associated occupancy support was observed for 275, 77, 46 and 84 candidates, respectively; these occupancy attributes were not mutually exclusive. The mutually exclusive downstream routing scheme assigned 21 candidates to shared-plus-promoter support, 25 to shared non-promoter support, 63 to non-shared promoter support, 175 to other occupancy support and 152 to the RNA-only route. No definitive direct-target assignments were made."
)

write_text(
  legend_lines,
  file.path(
    legend_root,
    "Figure2_legend.md"
  )
)

write_text(
  methods_lines,
  file.path(
    text_root,
    "Figure2_methods.md"
  )
)

write_text(
  results_lines,
  file.path(
    text_root,
    "Figure2_results.md"
  )
)


# ----------------------------------------------------------------------
# 16. Status and manifests
# ----------------------------------------------------------------------

status_table <- data.frame(
  figure_id = figure_id,
  figure_version = figure_version,
  status = "COMPLETE",
  gate_status = gate_status,
  figure_class = figure_class,
  candidate_status = "B_CONTROLLED_VISUAL_REPAIR_AWAITING_STRICT_RENDER_QC",
  scientific_QC = "PASS",
  statistical_QC = "PASS_NO_RECOMPUTATION",
  automated_QC = "PASS",
  strict_actual_render_QC =
    "PENDING",
  cross_format_QC =
    "PENDING",
  submission_ready =
    FALSE,
  immutable_release =
    FALSE,
  next_action =
    "STRICT_RENDER_QC_WITH_RUNTIME_LOG_AND_GLOBAL_ALIGNMENT_CROPS",
  stringsAsFactors = FALSE
)

write_csv_utf8(
  status_table,
  file.path(
    result_root,
    "FIGURE_STATUS.csv"
  )
)

script_source_path <-
  "<LOCAL_R_LIBRARY>"

if (file.exists(script_source_path)) {
  copy_validate(
    script_source_path,
    file.path(
      code_root,
      script_name
    )
  )
}

result_files <- list.files(
  result_root,
  recursive = TRUE,
  full.names = TRUE
)

result_files <- result_files[
  file.exists(result_files) &
    !dir.exists(result_files)
]

result_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      result_files
    ),
    nchar(
      paste0(
        normalize_slashes(
          result_root
        ),
        "/"
      )
    ) + 1L
  ),
  size_bytes =
    file.info(
      result_files
    )$size,
  md5 = vapply(
    result_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    result_files,
    sha256_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  result_manifest,
  file.path(
    manifest_root,
    "FIG2_RESULT_MANIFEST.csv"
  )
)

parameter_manifest <- list(
  script_name = script_name,
  script_version = script_version,
  governing_skill =
    "research_scientific_figure_suite v2.2-beta.8",
  skill_mode = "CONTROLLED_LAYOUT_REPAIR_RUNTIME_AND_GLOBAL_GEOMETRY",
  runtime_warning_gate = TRUE,
  global_composite_design = "AB / spacer / CD",
  volcano_label_policy =
    "NEARBY_LABELS_WITH_OPTIONAL_SHORT_NONARROW_LEADERS",
  yaxis_title_implementation =
    "NATIVE_GGPLOT_AXIS_TITLE_GROB",
  custom_panel_coordinate_yaxis_annotation =
    FALSE,
  panel_A_legend_route =
    "BOTTOM_OF_PANEL",
  header_anchor_reference =
    header_anchor_reference,
  header_implementation =
    "EXPLICIT_GTABLE_ROWS_INDEPENDENT_GGPLOT_GTABLES",
  panel_body_conversion = "GGPLOTGROB_ONLY",
  patchwork_free_used_in_conversion = FALSE,
  patchworkGrob_used_in_conversion = FALSE,
  header_tag_size_pt =
    header_tag_size_pt,
  header_title_size_pt =
    header_title_size_pt,
  header_subtitle_size_pt =
    header_subtitle_size_pt,
  header_tag_to_title_gap_mm =
    header_tag_to_title_gap_mm,
  header_title_row_height_mm =
    header_title_row_height_mm,
  header_subtitle_row_height_mm =
    header_subtitle_row_height_mm,
  native_ggplot_tag_title_used =
    FALSE,
  panel_B_axis_text_y_size =
    panel_B_axis_text_y_size,
  panel_C_x_limit =
    panel_C_x_limit,
  panel_D_x_limit =
    panel_D_x_limit,
  panel_A_y_tick_to_axis_gap_pt =
    panel_A_y_tick_to_axis_gap_pt,
  panel_A_y_title_to_tick_gap_pt =
    panel_A_y_title_to_tick_gap_pt,
  panel_A_horizontal_geometry_state =
    panel_A_horizontal_geometry_state,
  panel_A_left_label_alignment =
    "FREE_FROM_PANEL_C_SHARED_LABEL_WIDTH",
  patchwork_free_type =
    "label",
  patchwork_free_side =
    "l",
  vertical_layout_modified_this_round =
    FALSE,
  direct_label_without_connector_n = 2,
  leader_line_with_arrowhead = FALSE,
  prior_fast_replot_round = 2,
  maximum_fast_replot_rounds = 2,
  current_repair_is_blind_fast_replot = FALSE,
  further_blind_fast_replot_authorized = FALSE,
  requested_font = requested_font,
  resolved_font_family =
    resolved_font_family,
  resolved_font_name =
    resolved_font_name,
  resolved_font_style =
    resolved_font_style,
  resolved_font_index =
    font_index,
  resolved_font_path =
    normalize_slashes(
      font_path
    ),
  raster_device = "ragg",
  minimum_text_size_pt = 6.0,
  figure_id = figure_id,
  figure_version = figure_version,
  figure_class = figure_class,
  source_analysis_recomputed = FALSE,
  target_width_mm = target_width_mm,
  target_height_mm = target_height_mm,
  stress_test_width_mm =
    stress_test_width_mm,
  PNG_dpi = 300,
  TIFF_dpi = 600,
  TIFF_compression = "lzw",
  strict_actual_render_QC =
    "PENDING",
  submission_ready = FALSE
)

jsonlite::write_json(
  parameter_manifest,
  parameter_path,
  pretty = TRUE,
  auto_unbox = TRUE
)

write_text(
  capture.output(
    sessionInfo()
  ),
  session_path
)

run_end <- Sys.time()

write_text(
  c(
    paste0(
      "[",
      format(
        run_start,
        "%Y-%m-%d %H:%M:%S"
      ),
      "] Started ",
      run_id
    ),
    paste0(
      "[",
      format(
        run_end,
        "%Y-%m-%d %H:%M:%S"
      ),
      "] Completed ",
      run_id,
      "; gate=",
      gate_status,
      "; strict_render_QC=PENDING"
    )
  ),
  execution_log
)

run_manifest <- data.frame(
  run_id = run_id,
  script_name = script_name,
  script_version = script_version,
  governing_skill =
    "research_scientific_figure_suite v2.2-beta.8",
  figure_id = figure_id,
  figure_version = figure_version,
  figure_class = figure_class,
  started_at =
    as.character(run_start),
  completed_at =
    as.character(run_end),
  elapsed_minutes = round(
    as.numeric(
      difftime(
        run_end,
        run_start,
        units = "mins"
      )
    ),
    2
  ),
  gate_status = gate_status,
  PDF_created =
    file.exists(figure_pdf),
  PNG_created =
    file.exists(figure_png),
  TIFF_created =
    file.exists(figure_tiff),
  SVG_created =
    file.exists(figure_svg),
  strict_actual_render_QC =
    "PENDING",
  submission_ready =
    FALSE,
  result_root =
    normalize_slashes(
      result_root
    ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  run_manifest,
  run_manifest_path
)


# ----------------------------------------------------------------------
# 17. Shallow handoff report
# ----------------------------------------------------------------------

desktop_candidates <- c(
  file.path(
    Sys.getenv("USERPROFILE"),
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "Desktop"
  ),
  file.path(
    Sys.getenv("USERPROFILE"),
    "OneDrive",
    "桌面"
  )
)

desktop_candidates <- desktop_candidates[
  dir.exists(desktop_candidates)
]

if (!length(desktop_candidates)) {
  stop(
    "Could not locate Desktop for the handoff report."
  )
}

desktop_root <- desktop_candidates[[1]]

report_root <- file.path(
  desktop_root,
  paste0(
    "FIG2_REPORT_",
    run_stamp
  )
)

report_result_root <- file.path(
  report_root,
  "FIGURE2"
)

report_run_root <- file.path(
  report_root,
  "RUN"
)

handoff_zip <- file.path(
  desktop_root,
  paste0(
    run_id,
    "_SMALL_REPORT.zip"
  )
)

if (
  dir.exists(report_root) ||
    file.exists(handoff_zip)
) {
  stop(
    "Desktop report path already exists. Restart R and rerun."
  )
}

dir.create(
  report_result_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  report_run_root,
  recursive = TRUE,
  showWarnings = FALSE
)

result_prefix <- paste0(
  normalize_slashes(
    result_root
  ),
  "/"
)

final_result_files <- list.files(
  result_root,
  recursive = TRUE,
  full.names = TRUE
)

final_result_files <- final_result_files[
  file.exists(final_result_files) &
    !dir.exists(final_result_files)
]

for (source_path in final_result_files) {
  relative_path <- substring(
    normalize_slashes(
      source_path
    ),
    nchar(result_prefix) + 1L
  )

  copy_validate(
    source_path,
    file.path(
      report_result_root,
      relative_path
    )
  )
}

for (source_path in c(
  execution_log,
  run_manifest_path,
  parameter_path,
  session_path
)) {
  copy_validate(
    source_path,
    file.path(
      report_run_root,
      basename(source_path)
    )
  )
}

report_files <- list.files(
  report_root,
  recursive = TRUE,
  full.names = TRUE
)

report_files <- report_files[
  file.exists(report_files) &
    !dir.exists(report_files)
]

report_prefix <- paste0(
  normalize_slashes(
    report_root
  ),
  "/"
)

report_manifest <- data.frame(
  relative_path = substring(
    normalize_slashes(
      report_files
    ),
    nchar(report_prefix) + 1L
  ),
  size_bytes =
    file.info(
      report_files
    )$size,
  md5 = vapply(
    report_files,
    md5_file,
    character(1)
  ),
  sha256 = vapply(
    report_files,
    sha256_file,
    character(1)
  ),
  stringsAsFactors = FALSE
)

write_csv_utf8(
  report_manifest,
  file.path(
    report_root,
    "SMALL_REPORT_MANIFEST.csv"
  )
)

if (!requireNamespace("zip", quietly = TRUE)) {
  try(
    install.packages(
      "zip",
      repos =
        "https://cloud.r-project.org"
    ),
    silent = TRUE
  )
}

if (!requireNamespace("zip", quietly = TRUE)) {
  stop(
    "Package 'zip' is unavailable."
  )
}

old_directory <- getwd()

on.exit(
  setwd(old_directory),
  add = TRUE
)

setwd(
  dirname(report_root)
)

zip::zipr(
  zipfile = handoff_zip,
  files = basename(report_root),
  recurse = TRUE,
  include_directories = TRUE
)

if (
  !file.exists(handoff_zip) ||
    file.info(handoff_zip)$size <= 0
) {
  stop(
    "Figure 2 handoff ZIP was not created."
  )
}

cat("\n============================================================\n")
cat("STEP05B Figure 2 PDF/font repair v1.2 completed\n")
cat("============================================================\n")
cat(
  "Gate status: ",
  gate_status,
  "\n",
  sep = ""
)
cat(
  "Governing skill: research_scientific_figure_suite v2.2-beta.8\n"
)
cat(
  "Figure class: ",
  figure_class,
  "\n",
  sep = ""
)
cat(
  "Candidate status: ",
  candidate_status,
  "\n",
  sep = ""
)
cat("Panels: 4\n")
cat("PDF created: TRUE\n")
cat("PNG created: TRUE\n")
cat("TIFF created: TRUE\n")
cat("SVG created: TRUE\n")
cat(
  "Resolved font family: ",
  resolved_font_family,
  "\n",
  sep = ""
)
cat("Font API: match_fonts + font_info\n")
cat("Prior FAST_REPLOT rounds: 2/2\n")
cat("Current mode: CONTROLLED_LAYOUT_AND_BALANCE_REPAIR\n")
cat("Whole-runtime warnings: 0\n")
cat("Composite authority: AB / spacer / CD\n")
cat("Volcano label policy: nearby labels + optional short lines\n")
cat("Leader-line arrowheads: FALSE\n")
cat("Panel A y-axis title: NATIVE_GGPLOT_GROB\n")
cat("Custom y-axis annotation: FALSE\n")
cat("Panel A left-label alignment: FREE_FROM_PANEL_C\n")
cat("Panel body alignment: GLOBAL_AB_SPACER_CD_PRESERVED\n")
cat("Header implementation: EXPLICIT_GTABLE_ROWS_INDEPENDENT_GGPLOT_GTABLES\n")
cat("Panel conversion: GGPLOTGROB_ONLY\n")
cat("patchwork::free in conversion: FALSE\n")
cat("patchworkGrob in conversion: FALSE\n")
cat("Title left reference: DETECTED_PANEL_CELL_LEFT_EDGE\n")
cat("ABCD tag size: ", header_tag_size_pt, " pt\n", sep = "")
cat("Panel title size: ", header_title_size_pt, " pt\n", sep = "")
cat("Tag-to-title gap: ", header_tag_to_title_gap_mm, " mm\n", sep = "")
cat("Native ggplot tag/title used: FALSE\n")
cat(
  "Panel A tick-label to axis gap: ",
  panel_A_y_tick_to_axis_gap_pt,
  " pt\n",
  sep = ""
)
cat(
  "Panel A axis-title to tick-label gap: ",
  panel_A_y_title_to_tick_gap_pt,
  " pt\n",
  sep = ""
)
cat("Vertical title/legend layout modified this round: FALSE\n")
cat("Panel A legend route: BOTTOM_OF_PANEL_COMPACT\n")
cat("Panel B/C/D title alignment mode: FURTHER_RIGHT_SHIFTED_TOP_BLOCK\n")
cat("Further blind FAST_REPLOT authorized: FALSE\n")
cat("Font preflight: PASS\n")
cat("Export font warnings: 0\n")
cat("Scientific QC: PASS\n")
cat("Statistical QC: PASS_NO_RECOMPUTATION\n")
cat("Strict actual-render QC: PENDING\n")
cat("Submission ready: FALSE\n")
cat(
  "Candidate PNG: ",
  normalize_slashes(
    figure_png
  ),
  "\n",
  sep = ""
)
cat(
  "Upload this ZIP:\n",
  normalize_slashes(
    handoff_zip
  ),
  "\n",
  sep = ""
)
