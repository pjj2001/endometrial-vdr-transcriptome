#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B05: Export bulk results and single-cell input checks

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04D_QC1_", format(Sys.time(), "%Y%m%d_%H%M%S"))

repair_version <- "STEP04D_QC_REPAIR_v1.1"

projection_version <- "BULK_VDR_PROGRAM_PROJECTION_v1.3"

required_packages <- c("digest", "jsonlite", "limma")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(step04_root, "00_module_control"), code = file.path(step04_root, "02_code"), parameters = file.path(step04_root, 
    "03_parameters"), environment = file.path(step04_root, "04_environment"), execution = file.path(step04_root, 
    "05_execution"), results = file.path(step04_root, "06_results"), qc = file.path(step04_root, "08_qc"), text = file.path(step04_root, 
    "09_manuscript_text"), archive_parent = file.path(step04_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory_path)) {
        stop("Required directory could not be created: ", directory_path)
    }
}

output_root <- file.path(dirs$results, "S4D_qc_lock", "BVP_v1_3")

archive_root <- file.path(dirs$archive_parent, paste0("S4DQC_", format(Sys.time(), "%Y%m%d_%H%M%S")))

for (directory_path in c(output_root, file.path(archive_root, "res"), file.path(archive_root, "qc"), file.path(archive_root, 
    "txt"), file.path(archive_root, "code"), file.path(archive_root, "env"), file.path(archive_root, "log"), file.path(archive_root, 
    "ctl"))) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory_path)) {
        stop("Required output/archive directory could not be created: ", directory_path)
    }
    if (nchar(directory_path) >= 220L) {
        stop("Directory exceeds the 220-character soft limit: ", directory_path)
    }
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    message_text <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(message_text, "\n")
    cat(message_text, "\n", file = log_path, append = TRUE)
}

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
    if (nchar(path) >= 220L) {
        stop("Output path exceeds the 220-character soft limit: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Output directory could not be created: ", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

write_text <- function(lines, path) {
    if (nchar(path) >= 220L) {
        stop("Output path exceeds the 220-character soft limit: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

md5_file <- function(path) {
    if (!file.exists(path)) {
        return(NA_character_)
    }
    unname(tools::md5sum(path))
}

sha256_file <- function(path) {
    if (!file.exists(path)) {
        return(NA_character_)
    }
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

same_normalized_path <- function(path_a, path_b) {
    if (is.null(path_a) || is.null(path_b) || length(path_a) != 1L || length(path_b) != 1L || is.na(path_a) || 
        is.na(path_b)) {
        return(FALSE)
    }
    normalized_a <- tryCatch(normalizePath(path_a, winslash = "/", mustWork = FALSE), error = function(error_object) {
        normalize_slashes(path_a)
    })
    normalized_b <- tryCatch(normalizePath(path_b, winslash = "/", mustWork = FALSE), error = function(error_object) {
        normalize_slashes(path_b)
    })
    identical(tolower(normalized_a), tolower(normalized_b))
}

copy_if_different <- function(source_path, destination_path) {
    if (same_normalized_path(source_path, destination_path)) {
        return(TRUE)
    }
    dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(source_path, destination_path, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

copy_validate <- function(source_path, destination_path) {
    copied <- copy_if_different(source_path, destination_path)
    pass <- isTRUE(copied) && file.exists(destination_path) && identical(as.numeric(file.info(source_path)$size), 
        as.numeric(file.info(destination_path)$size)) && identical(md5_file(source_path), md5_file(destination_path)) && 
        identical(sha256_file(source_path), sha256_file(destination_path))
    if (!pass) {
        stop("Copy/hash validation failed: ", source_path)
    }
    TRUE
}

log_event("Starting", run_id)

v6_projection_root <- file.path(dirs$results, "S4D_bulk", "BVP_v1_2")

v6_gate_path <- file.path(dirs$results, "STEP04D_GATE_TABLE.csv")

v6_manifest_path <- file.path(dirs$results, "S4D_artifact_manifest.csv")

sample_scores_path <- file.path(v6_projection_root, "sample_scores.csv")

candidate_gene_results_path <- file.path(v6_projection_root, "candidate_gene_results.csv")

sc_sample_qc_path <- file.path(dirs$results, "S4D_sc_sample_qc.csv")

sc_dataset_qc_path <- file.path(dirs$results, "S4D_sc_dataset_qc.csv")

locked_design_root <- file.path(dirs$results, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0")

candidate_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

required_paths <- c(v6_gate_path, v6_manifest_path, sample_scores_path, candidate_gene_results_path, sc_sample_qc_path, 
    sc_dataset_qc_path, candidate_path, file.path(locked_design_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE165004.csv"), 
    file.path(locked_design_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE92324.csv"), file.path(locked_design_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE111974.csv"))

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0L) {
    stop("Required STEP04D-QC1 input is missing:\n", paste(missing_paths, collapse = "\n"))
}

v6_gate <- .vdr_read_csv(v6_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

v6_gate_status <- v6_gate$status[v6_gate$gate_component == "STEP04D gate"]

if (length(v6_gate_status) != 1L || v6_gate_status != "READY_FOR_SINGLE_CELL_INTEGRATION_AND_FINAL_MULTILAYER_CANDIDATE_PRIORITIZATION") {
    stop("STEP04D v6 did not finish with the expected READY Gate.")
}

v6_manifest <- .vdr_read_csv(v6_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

repair_input_paths <- c(sample_scores_path, candidate_gene_results_path, sc_sample_qc_path, sc_dataset_qc_path)

repair_input_names <- c("sample_scores", "candidate_gene_results", "sc_sample_qc", "sc_dataset_qc")

repair_input_normalized <- normalize_slashes(repair_input_paths)

manifest_path_normalized <- normalize_slashes(v6_manifest$path)

manifest_match_index <- match(repair_input_normalized, manifest_path_normalized)

input_identity_audit <- data.frame(input_name = repair_input_names, path = repair_input_normalized, file_exists = file.exists(repair_input_paths), 
    size_bytes = as.numeric(file.info(repair_input_paths)$size), listed_in_v6_manifest = !is.na(manifest_match_index), 
    manifest_md5 = NA_character_, observed_md5 = vapply(repair_input_paths, md5_file, character(1)), manifest_sha256 = NA_character_, 
    observed_sha256 = vapply(repair_input_paths, sha256_file, character(1)), identity_status = NA_character_, stringsAsFactors = FALSE)

listed_rows <- which(input_identity_audit$listed_in_v6_manifest)

if (length(listed_rows) > 0L) {
    matched_manifest_rows <- manifest_match_index[listed_rows]
    input_identity_audit$manifest_md5[listed_rows] <- as.character(v6_manifest$md5[matched_manifest_rows])
    input_identity_audit$manifest_sha256[listed_rows] <- as.character(v6_manifest$sha256[matched_manifest_rows])
}

input_identity_audit$identity_status <- ifelse(!input_identity_audit$file_exists, "FAIL_FILE_MISSING", ifelse(input_identity_audit$size_bytes <= 
    0 | is.na(input_identity_audit$size_bytes), "FAIL_EMPTY_FILE", ifelse(input_identity_audit$listed_in_v6_manifest & 
    input_identity_audit$manifest_md5 == input_identity_audit$observed_md5 & input_identity_audit$manifest_sha256 == 
    input_identity_audit$observed_sha256, "PASS_V6_MANIFEST_HASH_MATCH", ifelse(input_identity_audit$listed_in_v6_manifest, 
    "FAIL_V6_MANIFEST_HASH_MISMATCH", "PASS_DIRECT_HASH_CAPTURE_V6_MANIFEST_OMISSION"))))

write_csv_utf8(input_identity_audit, file.path(output_root, "STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv"))

if (any(grepl("^FAIL_", input_identity_audit$identity_status))) {
    stop("At least one STEP04D-QC1 repair input failed identity validation. Review STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv.")
}

manifest_omission_n <- sum(!input_identity_audit$listed_in_v6_manifest)

log_event("Repair input identity audit passed.", "Inputs listed in v6 manifest:", sum(input_identity_audit$listed_in_v6_manifest), 
    "/4;", "direct hash captures for manifest omissions:", manifest_omission_n)

read_locked_design <- function(dataset_id) {
    .vdr_read_csv(file.path(locked_design_root, paste0("COMBINED_SAMPLE_DESIGN_v1.0_", dataset_id, ".csv")), stringsAsFactors = FALSE, 
        check.names = FALSE)
}

design_165004 <- read_locked_design("GSE165004")

design_92324 <- read_locked_design("GSE92324")

design_111974 <- read_locked_design("GSE111974")

sample_scores <- .vdr_read_csv(sample_scores_path, stringsAsFactors = FALSE, check.names = FALSE)

candidate_gene_results <- .vdr_read_csv(candidate_gene_results_path, stringsAsFactors = FALSE, check.names = FALSE)

candidate_table <- .vdr_read_csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)

sc_sample_qc <- .vdr_read_csv(sc_sample_qc_path, stringsAsFactors = FALSE, check.names = FALSE)

sc_dataset_qc <- .vdr_read_csv(sc_dataset_qc_path, stringsAsFactors = FALSE, check.names = FALSE)

expected_programs <- c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")

if (nrow(candidate_table) != 436L || length(unique(candidate_table$feature_id)) != 436L) {
    stop("The locked 436-candidate table failed validation.")
}

if (!setequal(unique(sample_scores$program_name), expected_programs)) {
    stop("The sample-score file does not contain the four expected program names.")
}

if (nrow(sc_sample_qc) != 37L || sum(sc_sample_qc$dataset_id == "GSE250130") != 28L || sum(sc_sample_qc$dataset_id == 
    "GSE183837") != 9L) {
    stop("The 37-sample single-cell QC table failed validation.")
}

build_score_matrix <- function(dataset_id, sample_design) {
    score_subset <- sample_scores[sample_scores$dataset_id == dataset_id, , drop = FALSE]
    score_matrix <- matrix(NA_real_, nrow = length(expected_programs), ncol = nrow(sample_design), dimnames = list(expected_programs, 
        sample_design$sample_accession))
    for (program_name in expected_programs) {
        program_rows <- score_subset[score_subset$program_name == program_name, , drop = FALSE]
        program_rows <- program_rows[match(sample_design$sample_accession, program_rows$sample_accession), , drop = FALSE]
        if (nrow(program_rows) != nrow(sample_design) || any(is.na(program_rows$score))) {
            stop(dataset_id, " sample-score alignment failed for ", program_name)
        }
        score_matrix[program_name, ] <- program_rows$score
    }
    score_matrix
}

fit_score_matrix <- function(score_matrix, design_matrix, contrast_matrix, dataset_id) {
    if (ncol(score_matrix) != nrow(design_matrix)) {
        stop(dataset_id, " score/design dimensions do not match.")
    }
    fit <- limma::lmFit(score_matrix, design_matrix)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit, trend = FALSE, robust = FALSE)
    result_rows <- list()
    for (contrast_name in colnames(contrast_matrix)) {
        table_value <- limma::topTable(fit, coef = contrast_name, number = Inf, sort.by = "none", adjust.method = "BH")
        table_value$program_name <- rownames(table_value)
        table_value$dataset_id <- dataset_id
        table_value$contrast <- contrast_name
        table_value$score_model <- "FOUR_PRESPECIFIED_PROGRAMS_MODELED_TOGETHER"
        result_rows[[contrast_name]] <- table_value
    }
    do.call(rbind, result_rows)
}

score_matrix_165004 <- build_score_matrix("GSE165004", design_165004)

group_165004 <- factor(design_165004$group, levels = c("CONTROL", "RPL", "UIF"))

design_matrix_165004 <- model.matrix(~0 + group_165004)

colnames(design_matrix_165004) <- levels(group_165004)

contrast_matrix_165004 <- limma::makeContrasts(RPL_vs_CONTROL = RPL - CONTROL, UIF_vs_CONTROL = UIF - CONTROL, 
    levels = design_matrix_165004)

score_result_165004 <- fit_score_matrix(score_matrix_165004, design_matrix_165004, contrast_matrix_165004, "GSE165004")

score_matrix_92324 <- build_score_matrix("GSE92324", design_92324)

group_92324 <- factor(design_92324$group, levels = c("CONTROL", "IMPLANTATION_FAILURE"))

day_92324 <- factor(design_92324$collection_day)

design_matrix_92324 <- model.matrix(~0 + group_92324 + day_92324)

colnames(design_matrix_92324) <- make.names(colnames(design_matrix_92324))

contrast_matrix_92324 <- limma::makeContrasts(IF_vs_CONTROL = group_92324IMPLANTATION_FAILURE - group_92324CONTROL, 
    levels = design_matrix_92324)

score_result_92324 <- fit_score_matrix(score_matrix_92324, design_matrix_92324, contrast_matrix_92324, "GSE92324")

score_matrix_111974 <- build_score_matrix("GSE111974", design_111974)

group_111974 <- factor(design_111974$group, levels = c("CONTROL", "RIF"))

design_matrix_111974 <- model.matrix(~0 + group_111974)

colnames(design_matrix_111974) <- levels(group_111974)

contrast_matrix_111974 <- limma::makeContrasts(RIF_vs_CONTROL = RIF - CONTROL, levels = design_matrix_111974)

score_result_111974 <- fit_score_matrix(score_matrix_111974, design_matrix_111974, contrast_matrix_111974, "GSE111974")

score_models_corrected <- rbind(score_result_165004, score_result_92324, score_result_111974)

score_models_corrected <- score_models_corrected[, c("dataset_id", "contrast", "program_name", "logFC", "AveExpr", 
    "t", "P.Value", "adj.P.Val", "B", "score_model"), drop = FALSE]

if (nrow(score_models_corrected) != 16L || length(unique(score_models_corrected$program_name)) != 4L || any(score_models_corrected$program_name == 
    "1")) {
    stop("Corrected program-score models failed the 16-row/4-program validation.")
}

write_csv_utf8(score_models_corrected, file.path(output_root, "score_models_v1_3.csv"))

candidate_base_columns <- c("feature_id", "program_class", "interaction_logFC", "interaction_FDR", "occupancy_supported", 
    "shared_occupancy_supported", "any_promoter_overlap_support")

missing_candidate_columns <- setdiff(candidate_base_columns, names(candidate_table))

if (length(missing_candidate_columns) > 0L) {
    stop("Locked candidate table lacks required column(s): ", paste(missing_candidate_columns, collapse = ", "))
}

candidate_convergence <- candidate_table[, candidate_base_columns, drop = FALSE]

expected_contexts <- data.frame(dataset_id = c("GSE165004", "GSE165004", "GSE92324", "GSE111974"), contrast = c("UIF_vs_CONTROL", 
    "RPL_vs_CONTROL", "IF_vs_CONTROL", "RIF_vs_CONTROL"), context_id = c("GSE165004__UIF_vs_CONTROL", "GSE165004__RPL_vs_CONTROL", 
    "GSE92324__IF_vs_CONTROL", "GSE111974__RIF_vs_CONTROL"), stringsAsFactors = FALSE)

for (context_index in seq_len(nrow(expected_contexts))) {
    dataset_id <- expected_contexts$dataset_id[context_index]
    contrast_name <- expected_contexts$contrast[context_index]
    context_id <- expected_contexts$context_id[context_index]
    context_rows <- candidate_gene_results[!is.na(candidate_gene_results$dataset_id) & !is.na(candidate_gene_results$contrast) & 
        candidate_gene_results$dataset_id == dataset_id & candidate_gene_results$contrast == contrast_name, c("feature_id", 
        "logFC", "P.Value", "adj.P.Val"), drop = FALSE]
    if (anyDuplicated(context_rows$feature_id)) {
        stop("Duplicate candidate/context rows detected for ", context_id)
    }
    names(context_rows)[-1] <- paste0(context_id, "__", c("logFC", "P.Value", "adj.P.Val"))
    candidate_convergence <- merge(candidate_convergence, context_rows, by = "feature_id", all.x = TRUE, sort = FALSE)
}

logFC_columns <- paste0(expected_contexts$context_id, "__logFC")

adj_columns <- paste0(expected_contexts$context_id, "__adj.P.Val")

logFC_matrix <- as.matrix(candidate_convergence[, logFC_columns, drop = FALSE])

adj_matrix <- as.matrix(candidate_convergence[, adj_columns, drop = FALSE])

candidate_convergence$measured_context_n <- rowSums(is.finite(logFC_matrix))

candidate_convergence$positive_context_n <- rowSums(logFC_matrix > 0, na.rm = TRUE)

candidate_convergence$negative_context_n <- rowSums(logFC_matrix < 0, na.rm = TRUE)

candidate_convergence$FDR05_context_n <- rowSums(adj_matrix < 0.050000000000000003, na.rm = TRUE)

candidate_convergence$cross_context_direction_label <- ifelse(candidate_convergence$measured_context_n == 0L, "NOT_MEASURED_IN_ANY_BULK_CONTEXT", 
    ifelse(candidate_convergence$positive_context_n == candidate_convergence$measured_context_n, "POSITIVE_ALL_MEASURED_CONTEXTS", 
        ifelse(candidate_convergence$negative_context_n == candidate_convergence$measured_context_n, "NEGATIVE_ALL_MEASURED_CONTEXTS", 
            "CONTEXT_DEPENDENT_OR_MIXED_DIRECTION")))

if (nrow(candidate_convergence) != 436L || any(grepl("^NA__NA", names(candidate_convergence))) || any(candidate_convergence$measured_context_n == 
    0L & candidate_convergence$cross_context_direction_label != "NOT_MEASURED_IN_ANY_BULK_CONTEXT")) {
    stop("Corrected cross-context candidate table failed semantic validation.")
}

write_csv_utf8(candidate_convergence, file.path(output_root, "cross_context_candidates_v1_3.csv"))

context_coverage <- data.frame(measured_context_n = 0:4, candidate_n = as.integer(table(factor(candidate_convergence$measured_context_n, 
    levels = 0:4))), stringsAsFactors = FALSE)

write_csv_utf8(context_coverage, file.path(output_root, "candidate_context_coverage.csv"))

sc_sample_qc$QC_review_tier <- ifelse(sc_sample_qc$retained_fraction < 0.65000000000000002 | sc_sample_qc$median_percent_mito > 
    20, "REVIEW_REQUIRED", ifelse(sc_sample_qc$retained_fraction < 0.75 | sc_sample_qc$median_percent_mito > 15, 
    "CAUTION", "PASS"))

sc_sample_qc$integration_inclusion <- ifelse(sc_sample_qc$QC_review_tier == "PASS", "INCLUDE", "INCLUDE_WITH_DONOR_LEVEL_QC_REVIEW")

sc_watchlist <- sc_sample_qc[sc_sample_qc$QC_review_tier != "PASS", , drop = FALSE]

write_csv_utf8(sc_sample_qc, file.path(output_root, "sc_sample_qc_reviewed.csv"))

write_csv_utf8(sc_watchlist, file.path(output_root, "sc_sample_watchlist.csv"))

watchlist_summary <- aggregate(sample_accession ~ dataset_id + QC_review_tier, data = sc_sample_qc, FUN = length)

names(watchlist_summary)[names(watchlist_summary) == "sample_accession"] <- "sample_n"

write_csv_utf8(watchlist_summary, file.path(output_root, "sc_watchlist_summary.csv"))

score_significant_summary <- aggregate(adj.P.Val ~ dataset_id + contrast, data = transform(score_models_corrected, 
    adj.P.Val = adj.P.Val < 0.050000000000000003), FUN = sum)

names(score_significant_summary)[names(score_significant_summary) == "adj.P.Val"] <- "FDR05_program_n"

direction_counts <- as.data.frame(table(candidate_convergence$cross_context_direction_label), stringsAsFactors = FALSE)

names(direction_counts) <- c("direction_label", "candidate_n")

write_csv_utf8(score_significant_summary, file.path(output_root, "score_significance_summary.csv"))

write_csv_utf8(direction_counts, file.path(output_root, "candidate_direction_summary.csv"))

review_required_n <- sum(sc_sample_qc$QC_review_tier == "REVIEW_REQUIRED")

caution_n <- sum(sc_sample_qc$QC_review_tier == "CAUTION")

gate_status <- "READY_FOR_SINGLE_CELL_INTEGRATION_WITH_GSE183837_SAMPLE_WATCHLIST"

gate_table <- data.frame(gate_component = c("STEP04D v6 computational completion", "STEP04D v6 repair-input identity audit", 
    "Program-score name integrity", "Program-score multiplicity correction", "Cross-context table NA__NA removal", 
    "Unmeasured-candidate direction semantics", "Candidate row integrity", "Single-cell sample count", "Single-cell review-required samples", 
    "Single-cell caution samples", "Automatic sample exclusion", "Definitive direct-target assignment", "Scientific figures", 
    "Corrected bulk projection version", "STEP04D-QC1 lock", "STEP04D-QC1 gate"), status = c("PASS", paste0("PASS_", 
    4L - manifest_omission_n, "_MANIFEST_HASH_MATCH_", manifest_omission_n, "_DIRECT_HASH_CAPTURE"), "PASS_4_NAMED_PROGRAMS", 
    "BH_WITHIN_EACH_CONTRAST_ACROSS_4_PROGRAMS", "PASS_NO_NA_NA_COLUMNS", "PASS_EXPLICIT_NOT_MEASURED_LABEL", "PASS_436_UNIQUE_CANDIDATES", 
    "PASS_37_OF_37", paste0("WATCHLIST_", review_required_n), paste0("CAUTION_", caution_n), "NONE", "NOT_PERFORMED", 
    "NONE_CREATED", paste0(projection_version, "_LOCKED"), paste0(repair_version, "_LOCKED"), gate_status), required_action = c("Retain the v6 bulk matrices and single-cell pre-QC objects.", 
    "Retain STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv; direct hashes repair the v6 manifest omission without altering source files.", 
    "Use only the corrected v1.3 score table.", "Interpret FDR within each prespecified four-program contrast family.", 
    "Use only the corrected v1.3 candidate table.", "Do not classify unmeasured candidates as positive or negative.", 
    "Preserve the locked 436-candidate universe.", "Proceed donor-aware.", "Inspect these donors after normalization, doublet scoring and cell-type annotation.", 
    "Retain and inspect these donors during integration.", "Do not exclude donors solely from this pre-QC watchlist.", 
    "Continue occupancy-supported candidate terminology.", "Generate figures only after integrated biological QC.", 
    "Use v1.3 for downstream prioritization.", "Retain repair manifest and archive.", "Proceed to single-cell integration with the watchlist active."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(output_root, "STEP04D_QC1_GATE.csv"))

result_card <- c("# STEP04D-QC1 corrected handoff card", "", paste0("- Run ID: ", run_id), paste0("- Repair version: ", 
    repair_version, " LOCKED"), paste0("- Corrected bulk projection: ", projection_version, " LOCKED"), "- Bulk matrices retained: BULK_GENE_MATRIX_v1.2", 
    "- Single-cell pre-QC objects retained: SINGLE_CELL_PREQC_OBJECTS_v1.0", paste0("- Repair inputs listed in v6 manifest/direct-hash captured: ", 
        4L - manifest_omission_n, "/", manifest_omission_n), "- Program-score rows: 16", "- Named programs: 4", 
    "- Multiple-testing rule: BH across four prespecified programs within each contrast", "- Cross-context candidates: 436", 
    paste0("- Candidates measured in 0/1/2/3/4 contexts: ", paste(context_coverage$candidate_n, collapse = "/")), 
    paste0("- Single-cell review-required samples: ", review_required_n), paste0("- Single-cell caution samples: ", 
        caution_n), "- Automatic donor exclusions: 0", "- GSE111974 independent external validation: NO", "- Definitive direct-target assignments: 0", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status))

write_text(result_card, file.path(output_root, "STEP04D_QC1_CARD.md"))

write_text(c("# STEP04D-QC1 interpretation note", "", "The STEP04D v6 computational run completed successfully, but manual output review identified two table-construction defects: program names were lost in the score-model table, and the cross-context candidate table contained an artificial NA__NA context and misclassified candidates with zero measurable bulk contexts.", 
    "", "The corrected v1.3 score models preserve all four prespecified program names and control the false-discovery rate across the four programs within each dataset contrast. The corrected candidate table contains only the four authorized biological contexts and explicitly labels candidates not measured in any bulk cohort.", 
    "", "Single-cell preprocessing is technically adequate for integration. Samples on the review watchlist remain included; they require donor-level review after normalization, doublet assessment and cell-type annotation. No donor is excluded at this stage."), 
    file.path(dirs$text, "STEP04D_QC1_interpretation_note.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, repair_version = repair_version, projection_version = projection_version, 
    candidate_n = nrow(candidate_convergence), score_model_row_n = nrow(score_models_corrected), score_program_n = length(unique(score_models_corrected$program_name)), 
    sc_sample_n = nrow(sc_sample_qc), sc_review_required_n = review_required_n, sc_caution_n = caution_n, automatic_donor_exclusion_n = 0, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04D_QC1_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04D_QC1.txt"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) {
    NULL
})

if (!is.null(source_file) && file.exists(source_file)) {
    copy_if_different(source_file, file.path(dirs$code, "export_patient_bulk_results.R"))
}

artifact_paths <- c(file.path(output_root, "STEP04D_QC1_INPUT_IDENTITY_AUDIT.csv"), file.path(output_root, "score_models_v1_3.csv"), 
    file.path(output_root, "cross_context_candidates_v1_3.csv"), file.path(output_root, "candidate_context_coverage.csv"), 
    file.path(output_root, "sc_sample_qc_reviewed.csv"), file.path(output_root, "sc_sample_watchlist.csv"), file.path(output_root, 
        "sc_watchlist_summary.csv"), file.path(output_root, "score_significance_summary.csv"), file.path(output_root, 
        "candidate_direction_summary.csv"), file.path(output_root, "STEP04D_QC1_GATE.csv"), file.path(output_root, 
        "STEP04D_QC1_CARD.md"), file.path(dirs$text, "STEP04D_QC1_interpretation_note.md"), file.path(dirs$parameters, 
        "STEP04D_QC1_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04D_QC1.txt"), file.path(dirs$code, 
        "export_patient_bulk_results.R"), log_path)

missing_artifacts <- artifact_paths[!file.exists(artifact_paths)]

if (length(missing_artifacts) > 0L) {
    stop("Expected repair artifact is missing:\n", paste(missing_artifacts, collapse = "\n"))
}

artifact_manifest <- data.frame(path = normalize_slashes(artifact_paths), size_bytes = as.numeric(file.info(artifact_paths)$size), 
    md5 = vapply(artifact_paths, md5_file, character(1)), sha256 = vapply(artifact_paths, sha256_file, character(1)), 
    artifact_status = "LOCKED", stringsAsFactors = FALSE)

manifest_path <- file.path(output_root, "STEP04D_QC1_LOCKED_MANIFEST.csv")

write_csv_utf8(artifact_manifest, manifest_path)

archive_artifacts <- c(artifact_paths, manifest_path)

for (path_value in archive_artifacts) {
    section <- if (grepl("\\.R$", basename(path_value))) {
        "code"
    }
    else if (grepl("sessionInfo|parameters", basename(path_value))) {
        "env"
    }
    else if (grepl("\\.log$", basename(path_value))) {
        "log"
    }
    else if (grepl("\\.md$", basename(path_value))) {
        "txt"
    }
    else if (grepl("GATE|CARD|MANIFEST|SUMMARY|WATCHLIST|COVERAGE", basename(path_value), ignore.case = TRUE)) {
        "qc"
    }
    else {
        "res"
    }
    copy_validate(path_value, file.path(archive_root, section, basename(path_value)))
}

archive_files <- list.files(archive_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(normalize_slashes(archive_files), nchar(normalize_slashes(archive_root)) + 
    2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, character(1)), 
    sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04D_QC1_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04D-QC1 archive status", "", paste0("- Run ID: ", run_id), paste0("- Repair version: ", repair_version, 
    " LOCKED"), paste0("- Projection version: ", projection_version, " LOCKED"), "- Bulk matrices rerun: NO", "- Single-cell pre-QC rerun: NO", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Accepted correction run: ", run_id), "- Bulk matrices: BULK_GENE_MATRIX_v1.2 RETAINED", 
    paste0("- Bulk projection: ", projection_version, " LOCKED"), "- Single-cell pre-QC: SINGLE_CELL_PREQC_OBJECTS_v1.0 RETAINED_WITH_WATCHLIST", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04D-QC1 v2 repair and lock completed.\n")

cat("Repair version: ", repair_version, " LOCKED\n")

cat("Corrected projection: ", projection_version, " LOCKED\n")

cat("Program-score rows/names: ", nrow(score_models_corrected), " / ", length(unique(score_models_corrected$program_name)), 
    "\n", sep = "")

cat("Candidate rows: ", nrow(candidate_convergence), "\n")

cat("Candidate context coverage 0/1/2/3/4: ", paste(context_coverage$candidate_n, collapse = "/"), "\n")

cat("Single-cell REVIEW_REQUIRED / CAUTION: ", review_required_n, " / ", caution_n, "\n")

cat("Automatic donor exclusions: 0\n")

cat("Definitive direct-target assignments: 0\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/S4D_qc_lock/BVP_v1_3/STEP04D_QC1_GATE.csv\n")

cat("  06_results/S4D_qc_lock/BVP_v1_3/score_models_v1_3.csv\n")

cat("  06_results/S4D_qc_lock/BVP_v1_3/cross_context_candidates_v1_3.csv\n")

cat("  06_results/S4D_qc_lock/BVP_v1_3/sc_sample_watchlist.csv\n")

cat("  06_results/S4D_qc_lock/BVP_v1_3/STEP04D_QC1_CARD.md\n")

cat("============================================================\n")

