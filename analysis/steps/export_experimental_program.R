#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A07: Export the accepted RNA-program definition

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B1_v10_", format(Sys.time(), "%Y%m%d_%H%M%S"))

freeze_version <- "RNA_PROGRAM_v1.0"

candidate_status <- "B_LOCKED_DEFERRED_POLISH"

dirs <- list(control = file.path(module_root, "00_module_control"), code = file.path(module_root, "02_code"), parameters = file.path(module_root, 
    "03_parameters"), environment = file.path(module_root, "04_environment"), execution = file.path(module_root, 
    "05_execution"), results = file.path(module_root, "06_results"), figures = file.path(module_root, "07_figures_tables"), 
    qc = file.path(module_root, "08_qc"), text = file.path(module_root, "09_manuscript_text"), archive_parent = file.path(module_root, 
        "10_archive_work"))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_root <- file.path(dirs$figures, "B_submission_candidate", "STEP03B1_RNA_candidates_v3")

if (!dir.exists(candidate_root)) {
    stop("Candidate package not found:\n", candidate_root)
}

if (!requireNamespace("digest", quietly = TRUE)) {
    install.packages("digest", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos = "https://cloud.r-project.org")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

md5_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    unname(tools::md5sum(path))
}

sha256_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

copy_validate <- function(source_path, destination_path) {
    dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
    if (nchar(destination_path) >= 240L) {
        stop("Destination path exceeds the 240-character hard stop:\n", destination_path)
    }
    copied <- file.copy(source_path, destination_path, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    pass <- isTRUE(copied) && file.exists(destination_path) && identical(as.numeric(file.info(source_path)$size), 
        as.numeric(file.info(destination_path)$size)) && identical(md5_file(source_path), md5_file(destination_path)) && 
        identical(sha256_file(source_path), sha256_file(destination_path))
    if (!pass) {
        stop("Copy/hash validation failed:\n", source_path)
    }
    TRUE
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

log_event("Starting", run_id)

interaction_path <- file.path(dirs$results, "STEP03B1_v4_interaction_program_candidates_NOT_FROZEN.csv")

interaction_summary_path <- file.path(dirs$results, "STEP03B1_v4_interaction_program_summary.csv")

all_programs_path <- file.path(dirs$results, "STEP03B1_program_gene_candidates_NOT_FROZEN.csv")

definition_path <- file.path(dirs$results, "STEP03B1_v4_RNA_program_definition_card.csv")

required_analysis_files <- c(interaction_path, interaction_summary_path, all_programs_path, definition_path)

missing_analysis <- required_analysis_files[!file.exists(required_analysis_files)]

if (length(missing_analysis) > 0L) {
    stop("Missing analysis-lock inputs:\n", paste(missing_analysis, collapse = "\n"))
}

interaction <- .vdr_read_csv(interaction_path, stringsAsFactors = FALSE, check.names = FALSE)

interaction_summary <- .vdr_read_csv(interaction_summary_path, stringsAsFactors = FALSE, check.names = FALSE)

all_programs <- .vdr_read_csv(all_programs_path, stringsAsFactors = FALSE, check.names = FALSE)

definition <- .vdr_read_csv(definition_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(interaction) != 436L) {
    stop("Expected 436 interaction candidates, found ", nrow(interaction))
}

required_summary_columns <- c("program_class", "n_genes")

missing_summary_columns <- setdiff(required_summary_columns, names(interaction_summary))

if (length(missing_summary_columns) > 0L) {
    stop("Interaction summary is missing columns: ", paste(missing_summary_columns, collapse = ", "))
}

interaction_summary$n_genes <- suppressWarnings(as.integer(as.character(interaction_summary$n_genes)))

if (any(!is.finite(interaction_summary$n_genes)) || any(interaction_summary$n_genes < 0)) {
    stop("Interaction summary contains invalid n_genes values.")
}

summary_total <- sum(interaction_summary$n_genes)

summary_attenuated <- sum(interaction_summary$n_genes[grepl("^attenuated_ligand_", interaction_summary$program_class)])

candidate_total <- nrow(interaction)

if (!"program_class" %in% names(interaction)) {
    stop("Interaction candidate table is missing program_class.")
}

candidate_attenuated <- sum(grepl("^attenuated_ligand_", interaction$program_class))

if (candidate_total != 436L || candidate_attenuated != 431L) {
    stop("Candidate-table totals differ from the reviewed lock: total=", candidate_total, "; attenuated=", candidate_attenuated, 
        ".")
}

if (summary_total != candidate_total || summary_attenuated != candidate_attenuated) {
    stop("Interaction summary does not reconcile with the candidate table: ", "summary total=", summary_total, 
        ", candidate total=", candidate_total, "; summary attenuated=", summary_attenuated, ", candidate attenuated=", 
        candidate_attenuated, ".")
}

interaction_count_reconciliation <- data.frame(metric = c("candidate_total", "summary_total", "candidate_attenuated", 
    "summary_attenuated"), value = c(candidate_total, summary_total, candidate_attenuated, summary_attenuated), 
    expected = c(436L, 436L, 431L, 431L), pass = c(candidate_total == 436L, summary_total == 436L, candidate_attenuated == 
        431L, summary_attenuated == 431L), stringsAsFactors = FALSE)

write_csv_utf8(interaction_count_reconciliation, file.path(dirs$results, "STEP03B1_v10_INTERACTION_COUNT_RECONCILIATION.csv"))

required_candidate_files <- c(file.path(candidate_root, "figures", "Figure_STEP03B1_interaction_volcano_candidate_v3_QC.png"), 
    file.path(candidate_root, "figures", "Figure_STEP03B1_PCA_candidate_v3_QC.png"), file.path(candidate_root, 
        "figures", "Figure_STEP03B1_top_interaction_heatmap_ACCEPTED_v2_QC.png"), file.path(candidate_root, "figures", 
        "Figure_STEP03B1_VDR_expression_candidate_v3_QC.png"), file.path(candidate_root, "manifests", "STEP03B1_v3_FIGURE_SOURCE_MAP.csv"), 
    file.path(candidate_root, "source_data", "README_SOURCE_DATA.md"), file.path(candidate_root, "LANGUAGE_ROUTING_CARD.md"), 
    file.path(candidate_root, "legends", "STEP03B1_v3_candidate_legends.md"))

missing_candidate <- required_candidate_files[!file.exists(required_candidate_files)]

if (length(missing_candidate) > 0L) {
    stop("Class B candidate package is incomplete:\n", paste(missing_candidate, collapse = "\n"))
}

direct_sources <- list.files(file.path(candidate_root, "source_data"), pattern = "\\.csv$", full.names = TRUE, 
    recursive = FALSE)

upstream_refs <- list.files(file.path(candidate_root, "source_data", "upstream_refs"), pattern = "\\.csv$", full.names = TRUE)

generator_code <- list.files(file.path(candidate_root, "code"), pattern = "\\.R$", full.names = TRUE)

if (length(direct_sources) < 6L || length(upstream_refs) < 4L || length(generator_code) < 1L) {
    stop("The self-contained source package does not meet the Class B minimum.")
}

freeze_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

add_freeze <- function(x, scope) {
    cbind(data.frame(program_freeze_version = rep(freeze_version, nrow(x)), program_freeze_run = rep(run_id, nrow(x)), 
        program_freeze_date = rep(freeze_date, nrow(x)), membership_status = rep("FROZEN_AT_RNA_STATISTICAL_LAYER", 
            nrow(x)), program_scope = rep(scope, nrow(x)), direct_target_tier = rep("PENDING_ATAC_AND_CUTRUN", 
            nrow(x)), stringsAsFactors = FALSE), x)
}

interaction_frozen <- add_freeze(interaction, paste("Formal siRNA-by-ligand interaction;", "FDR < 0.05 and absolute interaction log2FC >= 0.5"))

all_programs_frozen <- add_freeze(all_programs, "Five prespecified GSE254251 RNA contrasts")

definition$program_freeze_version <- freeze_version

definition$program_freeze_run <- run_id

definition$program_freeze_date <- freeze_date

definition$membership_status <- "FROZEN_AT_RNA_STATISTICAL_LAYER"

definition$direct_target_tier <- "PENDING_ATAC_AND_CUTRUN"

write_csv_utf8(interaction_frozen, file.path(dirs$results, "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"))

write_csv_utf8(all_programs_frozen, file.path(dirs$results, "STEP03B1_v10_RNA_all_contrast_programs_FROZEN.csv"))

write_csv_utf8(definition, file.path(dirs$results, "STEP03B1_v10_RNA_program_definition_FROZEN.csv"))

figure_paths <- required_candidate_files[1:4]

figure_lock <- data.frame(figure_id = c("STEP03B1-RNA-VOLCANO", "STEP03B1-RNA-PCA", "STEP03B1-RNA-HEATMAP", "STEP03B1-RNA-VDR"), 
    figure_file = basename(figure_paths), current_class = "B_submission_candidate", candidate_lock_status = candidate_status, 
    scientific_qc = "PASS", statistical_qc = "PASS", candidate_readability_qc = "PASS_FOR_RESULT_VERIFICATION", 
    final_visual_qc = "DEFERRED", final_release_status = "NOT_STARTED", proposed_role = c("MAIN_CANDIDATE", "SUPPLEMENTARY_CANDIDATE", 
        "MAIN_CANDIDATE", "SUPPORTING_MECHANISTIC_CANDIDATE"), deferred_defect = c(paste("Rebalance final labels and typography", 
        "after composite width is known."), paste("Reposition centroid labels after", "supplementary-panel geometry is known."), 
        paste("Recheck font size and color bar", "inside the final composite."), paste("Refine point/summary geometry", 
            "at final physical panel width.")), md5 = vapply(figure_paths, md5_file, character(1)), sha256 = vapply(figure_paths, 
        sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(figure_lock, file.path(dirs$results, "STEP03B1_v10_B_CANDIDATE_LOCK_DECISION.csv"))

write_text(c("# STEP03B1 deferred final figure release", "", "- Current class: B_submission_candidate", paste0("- Candidate lock status: ", 
    candidate_status), "- Scientific QC: PASS", "- Statistical QC: PASS", "- Candidate readability: sufficient to verify the result", 
    "- Final publication visual QC: DEFERRED", "- PDF/TIFF/SVG release: NOT_STARTED", "- Main/supplementary final numbering: DEFERRED", 
    "- Composite geometry: DEFERRED", "", "The current QC PNGs are not final submission figures.", "All source tables, upstream references, exact code, legends, QC cards", 
    "and manifests are retained for manuscript-wide release without", "rerunning the RNA analysis."), file.path(candidate_root, 
    "qc", "STEP03B1_v9_DEFERRED_FINAL_RELEASE_CARD.md"))

write_text(c("# STEP03B1 RNA program lock card", "", paste0("- Version: ", freeze_version), paste0("- Run: ", run_id), 
    "- Samples: 12/12 retained", paste0("- Interaction candidates: ", candidate_total), paste0("- Attenuation-pattern candidates: ", 
        candidate_attenuated), "- RNA statistical membership: FROZEN", "- Direct-target tier: PENDING ATAC/CUT&RUN", 
    "- Figure source package: LOCKED", paste0("- Candidate visual status: ", candidate_status), "- Final submission artwork: DEFERRED", 
    "", "ATAC and CUT&RUN may add support tiers but may not silently change", "RNA membership."), file.path(dirs$results, 
    "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

gate <- data.frame(gate_component = c("RNA scientific/statistical QC", "RNA interaction membership", "Class B source package", 
    "Candidate readability", "Final publication visual QC", "C-class release", "Direct-target tier", "STEP03B1 v10 gate"), 
    status = c("PASS", "FROZEN_RNA_PROGRAM_v1.0", "PASS_COMPLETE", "PASS_FOR_RESULT_VERIFICATION", "DEFERRED_TO_MANUSCRIPT_RELEASE", 
        "NOT_STARTED", "PENDING_ATAC_AND_CUTRUN", "READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED"), 
    required_action = c("Retain accepted analysis outputs.", "Do not silently change membership.", "Use source map for later figure release.", 
        "No more standalone candidate polishing now.", "Perform once after manuscript figure plan is stable.", 
        "Generate PDF/TIFF/PNG/SVG only in RELEASE_FIGURE.", "Assign only after cross-omics integration.", "Begin STEP03B2 ATAC."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate, file.path(dirs$results, "STEP03B1_v10_GATE_TABLE.csv"))

write_text(c("# STEP03B2 ATAC authorization", "", paste0("- Authorized by: ", run_id), "- Authorization: YES", 
    "- RNA anchor: RNA_PROGRAM_v1.0", "- RNA membership: IMMUTABLE", "- Figure C-class release: NOT REQUIRED FOR ATAC START", 
    "- ATAC role: chromatin-accessibility support and mechanistic tiering", "- Prohibited: count-based differential-accessibility claims from", 
    "  normalized BigWig and condition-level merged BED alone", "- Required: self-contained source package for every retained ATAC figure"), 
    file.path(dirs$control, "STEP03B2_AUTHORIZATION.md"))

write_text(c("# STEP03B1 transition", "", "RNA scientific/statistical results and program membership are locked.", 
    "The four current plots are reproducible B-class candidates with final", "visual polish deferred to manuscript-wide figure release.", 
    "Proceed to STEP03B2 ATAC without rerunning RNA analysis or continuing", "standalone cosmetic iterations."), 
    file.path(dirs$text, "STEP03B1_v10_TRANSITION_TO_STEP03B2.md"))

archive_root <- file.path(dirs$archive_parent, paste0("S3B1v10_", format(Sys.time(), "%Y%m%d_%H%M%S")))

for (d in c(file.path(archive_root, "B3"), file.path(archive_root, "res"), file.path(archive_root, "ctl"), file.path(archive_root, 
    "txt"), file.path(archive_root, "code"), file.path(archive_root, "env"), file.path(archive_root, "log"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

candidate_files <- list.files(candidate_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

copy_rows <- vector("list", length(candidate_files))

for (i in seq_along(candidate_files)) {
    source_path <- candidate_files[i]
    relative_path <- substring(normalize_slashes(source_path), nchar(normalize_slashes(candidate_root)) + 2L)
    destination_path <- file.path(archive_root, "B3", relative_path)
    copy_validate(source_path, destination_path)
    copy_rows[[i]] <- data.frame(relative_path = relative_path, source_md5 = md5_file(source_path), archive_md5 = md5_file(destination_path), 
        source_sha256 = sha256_file(source_path), archive_sha256 = sha256_file(destination_path), validation_pass = TRUE, 
        stringsAsFactors = FALSE)
}

copy_manifest <- do.call(rbind, copy_rows)

write_csv_utf8(copy_manifest, file.path(archive_root, "ctl", "STEP03B1_v10_B3_COPY_MANIFEST.csv"))

same_normalized_path <- function(path_a, path_b) {
    if (is.null(path_a) || is.null(path_b) || length(path_a) != 1L || length(path_b) != 1L || is.na(path_a) || 
        is.na(path_b)) {
        return(FALSE)
    }
    normalized_a <- tryCatch(normalizePath(path_a, winslash = "/", mustWork = FALSE), error = function(error_object) {
        gsub("\\\\", "/", path_a)
    })
    normalized_b <- tryCatch(normalizePath(path_b, winslash = "/", mustWork = FALSE), error = function(error_object) {
        gsub("\\\\", "/", path_b)
    })
    identical(tolower(normalized_a), tolower(normalized_b))
}

safe_file_copy <- function(from, to, overwrite = FALSE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE) {
    if (same_normalized_path(from, to)) {
        return(TRUE)
    }
    base::file.copy(from = from, to = to, overwrite = overwrite, recursive = recursive, copy.mode = copy.mode, 
        copy.date = copy.date)
}

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    project_code <- .vdr_script_path("20_step03B1_v10_lock_RNA_defer_final_figures_authorize_ATAC.R")
    module_code <- file.path(dirs$code, "export_experimental_program.R")
    safe_file_copy(source_file, project_code, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, module_code, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

parameters <- list(schema_version = "1.0", run_id = run_id, policy = "candidate-first source-locked release-deferred", 
    RNA_program_version = freeze_version, RNA_membership_frozen = TRUE, figure_class = "B_submission_candidate", 
    candidate_lock_status = candidate_status, final_figure_release = "DEFERRED", ATAC_authorized = TRUE)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B1_v10_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B1_v10.txt"))

archive_artifacts <- c(file.path(dirs$results, "STEP03B1_v10_INTERACTION_COUNT_RECONCILIATION.csv"), file.path(dirs$results, 
    "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"), file.path(dirs$results, "STEP03B1_v10_RNA_all_contrast_programs_FROZEN.csv"), 
    file.path(dirs$results, "STEP03B1_v10_RNA_program_definition_FROZEN.csv"), file.path(dirs$results, "STEP03B1_v10_B_CANDIDATE_LOCK_DECISION.csv"), 
    file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"), file.path(dirs$results, "STEP03B1_v10_GATE_TABLE.csv"), 
    file.path(dirs$text, "STEP03B1_v10_TRANSITION_TO_STEP03B2.md"), file.path(dirs$control, "STEP03B2_AUTHORIZATION.md"), 
    file.path(dirs$parameters, "STEP03B1_v10_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B1_v10.txt"), 
    log_path)

for (p in archive_artifacts[file.exists(archive_artifacts)]) {
    section <- if (grepl("\\.csv$|LOCK_CARD", basename(p))) {
        "res"
    }
    else if (grepl("TRANSITION|AUTHORIZATION", basename(p))) {
        "txt"
    }
    else if (grepl("sessionInfo|parameters", basename(p))) {
        "env"
    }
    else {
        "log"
    }
    copy_validate(p, file.path(archive_root, section, basename(p)))
}

archive_files <- list.files(archive_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(normalize_slashes(archive_files), nchar(normalize_slashes(archive_root)) + 
    2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, character(1)), 
    sha256 = vapply(archive_files, sha256_file, character(1)), path_length = nchar(archive_files), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B1_v10_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B1 v10 archive status", "", paste0("- Run ID: ", run_id), "- RNA result: LOCKED", paste0("- Candidate figure status: ", 
    candidate_status), "- Final publication visual QC: DEFERRED", "- C-class release: NOT_STARTED", paste0("- Candidate package files copied: ", 
    nrow(copy_manifest)), paste0("- Maximum archive path length: ", max(archive_manifest$path_length)), "- STEP03B2 ATAC: AUTHORIZED"), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Completed substep: STEP03B1", paste0("- Accepted run: ", run_id), 
    "- RNA program: FROZEN RNA_PROGRAM_v1.0", paste0("- RNA figures: ", candidate_status), "- Final figure release: DEFERRED", 
    "- Active next substep: STEP03B2 ATAC", "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate: READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED")

cat("\n============================================================\n")

cat("STEP03B1 v10 candidate-first closeout completed.\n")

cat("RNA program frozen: RNA_PROGRAM_v1.0\n")

cat("Interaction candidates: ", candidate_total, "\n")

cat("Attenuation-pattern candidates: ", candidate_attenuated, "\n")

cat("Current figure class: B_submission_candidate\n")

cat("Candidate lock: B_LOCKED_DEFERRED_POLISH\n")

cat("Final submission figure release: DEFERRED\n")

cat("Gate: READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED\n")

cat("Authoritative archive:\n  ", archive_root, "\n")

cat("============================================================\n")

