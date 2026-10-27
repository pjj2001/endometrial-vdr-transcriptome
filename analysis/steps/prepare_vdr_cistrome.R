#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A10: Prepare accepted published VDR cistrome sources

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3B_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"

RNA_anchor <- "RNA_PROGRAM_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

required_packages <- c("digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(module_root, "00_module_control"), code = file.path(module_root, "02_code"), parameters = file.path(module_root, 
    "03_parameters"), environment = file.path(module_root, "04_environment"), execution = file.path(module_root, 
    "05_execution"), results = file.path(module_root, "06_results"), qc = file.path(module_root, "08_qc"), text = file.path(module_root, 
    "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

archive_root <- file.path(dirs$archive_parent, paste0("S3B3Bv3_", format(Sys.time(), "%Y%m%d_%H%M%S")))

for (directory_path in c(file.path(archive_root, "res"), file.path(archive_root, "qc"), file.path(archive_root, 
    "txt"), file.path(archive_root, "code"), file.path(archive_root, "env"), file.path(archive_root, "log"), file.path(archive_root, 
    "ctl"))) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    message_text <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(message_text, "\n")
    cat(message_text, "\n", file = log_path, append = TRUE)
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
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

copy_validate <- function(source_path, destination_path) {
    dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
    if (nchar(destination_path) >= 240L) {
        stop("Archive destination exceeds 240 characters:\n", destination_path)
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

is_nonempty_value <- function(x) {
    !is.na(x) & nzchar(trimws(as.character(x)))
}

log_event("Starting", run_id)

required_paths <- c(gate_v2 = file.path(dirs$results, "STEP03B3B_v2_GATE_TABLE.csv"), source_candidate_v2 = file.path(dirs$results, 
    "STEP03B3B_v2_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"), coordinate_summary_v2 = file.path(dirs$qc, "STEP03B3B_v2_COORDINATE_SCHEMA_SUMMARY.csv"), 
    chromosome_style_v2 = file.path(dirs$qc, "STEP03B3B_v2_CHROMOSOME_STYLE_SUMMARY.csv"), invalid_rows_v2 = file.path(dirs$qc, 
        "STEP03B3B_v2_INVALID_COORDINATE_ROWS.csv"), reconciliation_v2 = file.path(dirs$qc, "STEP03B3B_v2_PUBLISHED_TABLE_RECONCILIATION.csv"), 
    RNA_lock_card = file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"), ATAC_role_policy = file.path(dirs$results, 
        "STEP03B2A2_ATAC_EVIDENCE_ROLE_POLICY.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required source-lock input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v2 <- .vdr_read_csv(required_paths[["gate_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

source_candidate_v2 <- .vdr_read_csv(required_paths[["source_candidate_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

coordinate_summary_v2 <- .vdr_read_csv(required_paths[["coordinate_summary_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

chromosome_style_v2 <- .vdr_read_csv(required_paths[["chromosome_style_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

invalid_rows_v2 <- .vdr_read_csv(required_paths[["invalid_rows_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

reconciliation_v2 <- .vdr_read_csv(required_paths[["reconciliation_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

ATAC_role_policy <- .vdr_read_csv(required_paths[["ATAC_role_policy"]], stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- gate_v2$status[gate_v2$gate_component == "STEP03B3B v2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_PUBLISHED_PEAK_SOURCE_LOCK") {
    stop("STEP03B3B v2 does not authorize source lock. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

if (!all(ATAC_role_policy$status == ATAC_role)) {
    stop("ATAC role policy is not uniformly restricted descriptive support.")
}

expected_audits <- c("SHARED_TO_VEHICLE", "SHARED_TO_CALCITRIOL", "VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME", "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME")

if (nrow(reconciliation_v2) != 4L || !setequal(reconciliation_v2$audit_id, expected_audits) || !all(reconciliation_v2$any_overlap_fraction == 
    1) || !all(reconciliation_v2$exact_match_fraction == 1)) {
    stop("Published-table reconciliation is not exact for all four locked audits.")
}

if (nrow(source_candidate_v2) != 5L || !all(source_candidate_v2$source_lock_status == "CANDIDATE_FOR_MANUAL_SOURCE_LOCK") || 
    !all(source_candidate_v2$chromosome_style == "canonical_UCSC_chr")) {
    stop("The five source-lock candidates are incomplete.")
}

source_candidate_v2$file_exists <- file.exists(source_candidate_v2$clean_table_path)

source_candidate_v2$observed_md5 <- vapply(source_candidate_v2$clean_table_path, md5_file, character(1))

source_candidate_v2$observed_sha256 <- vapply(source_candidate_v2$clean_table_path, sha256_file, character(1))

source_candidate_v2$identity_pass <- with(source_candidate_v2, file_exists & clean_table_md5 == observed_md5 & 
    clean_table_sha256 == observed_sha256)

if (!all(source_candidate_v2$identity_pass)) {
    stop("One or more canonical clean source tables failed hash revalidation.")
}

if (nrow(invalid_rows_v2) != 16L || !all(invalid_rows_v2$source_id_from_list == "VDR_CALCITRIOL_CISTROME") || !all(!invalid_rows_v2$valid_coordinate)) {
    stop("The invalid-row audit is not the expected 16-row calcitriol appendix.")
}

coordinate_fields_empty <- (!is_nonempty_value(invalid_rows_v2$chromosome_source)) & is.na(invalid_rows_v2$start0) & 
    is.na(invalid_rows_v2$end)

if (!all(coordinate_fields_empty)) {
    stop("At least one invalid row contains a partial genomic coordinate and requires manual review.")
}

has_summary_label <- is_nonempty_value(invalid_rows_v2$raw_geneOverlap)

has_summary_value <- if ("raw_promoterOverlap" %in% names(invalid_rows_v2)) {
    is_nonempty_value(invalid_rows_v2$raw_promoterOverlap)
} else {
    rep(FALSE, nrow(invalid_rows_v2))
}

invalid_rows_v2$manual_row_class <- ifelse(has_summary_label | has_summary_value, "NON_COORDINATE_APPENDIX_SUMMARY_ROW", 
    "BLANK_SEPARATOR_ROW")

invalid_rows_v2$peak_membership_status <- "NOT_A_PEAK_ROW"

invalid_rows_v2$coordinate_exclusion_reason <- ifelse(invalid_rows_v2$manual_row_class == "NON_COORDINATE_APPENDIX_SUMMARY_ROW", 
    "AUTHOR_APPENDED_ANNOTATION_SUMMARY", "AUTHOR_APPENDED_BLANK_SEPARATOR")

invalid_rows_v2$manual_decision <- "RETAIN_IN_AUDIT_EXCLUDE_FROM_COORDINATE_SOURCE"

write_csv_utf8(invalid_rows_v2, file.path(dirs$qc, "STEP03B3B_v3_NON_COORDINATE_APPENDIX_ROW_DECISIONS.csv"))

appendix_summary <- as.data.frame(table(invalid_rows_v2$manual_row_class), stringsAsFactors = FALSE)

names(appendix_summary) <- c("manual_row_class", "n_rows")

appendix_summary$source_id <- "VDR_CALCITRIOL_CISTROME"

appendix_summary$peak_rows_removed <- 0L

appendix_summary$decision <- "NON_PEAK_APPENDIX_ROWS_EXCLUDED_FROM_COORDINATE_SOURCE"

write_csv_utf8(appendix_summary, file.path(dirs$qc, "STEP03B3B_v3_NON_COORDINATE_APPENDIX_SUMMARY.csv"))

if (sum(invalid_rows_v2$manual_row_class == "NON_COORDINATE_APPENDIX_SUMMARY_ROW") != 13L || sum(invalid_rows_v2$manual_row_class == 
    "BLANK_SEPARATOR_ROW") != 3L) {
    stop("The 16-row appendix did not resolve to 13 summary rows and 3 blank separators.")
}

expected_counts <- data.frame(source_id = c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME", "VDR_SHARED_CISTROME", 
    "VEHICLE_CISTROME_RNA_OVERLAP", "CALCITRIOL_CISTROME_RNA_OVERLAP"), total_source_rows = c(885L, 6108L, 519L, 
    134L, 2022L), valid_coordinate_rows = c(885L, 6092L, 519L, 134L, 2022L), primary_nuclear_rows = c(884L, 6091L, 
    518L, 134L, 2022L), nonprimary_standard_rows = c(1L, 1L, 1L, 0L, 0L), stringsAsFactors = FALSE)

summary_source_map <- c(MEDIA6_pLentiVeh = "VDR_VEHICLE_CISTROME", MEDIA6_pLenti125 = "VDR_CALCITRIOL_CISTROME", 
    MEDIA6_pLentiVeh_125overlap = "VDR_SHARED_CISTROME", MEDIA2_OverlapPeaksVeh_ReceptorDEG = "VEHICLE_CISTROME_RNA_OVERLAP", 
    MEDIA2_OverlapPeaks125_125DEG = "CALCITRIOL_CISTROME_RNA_OVERLAP")

observed_counts <- coordinate_summary_v2

observed_counts$source_id_locked <- unname(summary_source_map[observed_counts$source_id])

observed_counts$nonprimary_standard_rows <- observed_counts$valid_coordinate_rows - observed_counts$primary_coordinate_rows

count_check <- merge(expected_counts, observed_counts[c("source_id_locked", "total_rows", "valid_coordinate_rows", 
    "primary_coordinate_rows", "nonprimary_standard_rows")], by.x = "source_id", by.y = "source_id_locked", all = TRUE, 
    sort = FALSE, suffixes = c("_expected", "_observed"))

count_check$count_pass <- with(count_check, total_source_rows == total_rows & valid_coordinate_rows_expected == 
    valid_coordinate_rows_observed & primary_nuclear_rows == primary_coordinate_rows & nonprimary_standard_rows_expected == 
    nonprimary_standard_rows_observed)

write_csv_utf8(count_check, file.path(dirs$qc, "STEP03B3B_v3_SOURCE_COUNT_LOCK_CHECK.csv"))

if (!all(count_check$count_pass)) {
    stop("One or more published source counts differ from the reviewed lock.")
}

source_lock <- source_candidate_v2

source_lock$source_lock_version <- source_lock_version

source_lock$source_class <- ifelse(source_lock$source_id %in% c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME", 
    "VDR_SHARED_CISTROME"), "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_REFERENCE")

source_lock$allowed_project_role <- ifelse(source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", 
    "PROJECT_SPECIFIC_PEAK_TO_GENE_MAPPING_AND_RNA_OCCUPANCY_SUPPORT", "RECONCILIATION_AND_AUTHOR_RESULT_CROSS_CHECK_ONLY")

source_lock$prohibited_project_role <- ifelse(source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", 
    "CANNOT_ALONE_ESTABLISH_CAUSAL_OR_DEFINITIVE_DIRECT_TARGET_STATUS", "CANNOT_REPLACE_RNA_PROGRAM_v1.0_OR_SERVE_AS_INDEPENDENT_VALIDATION")

source_lock$nuclear_integration_policy <- "USE_chr1_TO_chr22_chrX_chrY;RETAIN_chrM_SEPARATELY"

source_lock$local_peak_calling_status <- "DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS"

source_lock$source_lock_status <- "LOCKED"

write_csv_utf8(source_lock, file.path(dirs$results, "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"))

source_role_summary <- data.frame(source_class = c("PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_REFERENCE"), 
    n_sources = c(3L, 2L), project_authority = c("PRIMARY_CHROMATIN_OCCUPANCY_SOURCE", "DERIVED_CROSS_CHECK_ONLY"), 
    independent_validation = c(FALSE, FALSE), changes_RNA_membership = c(FALSE, FALSE), stringsAsFactors = FALSE)

write_csv_utf8(source_role_summary, file.path(dirs$results, "STEP03B3B_v3_SOURCE_ROLE_SUMMARY.csv"))

gate_table <- data.frame(gate_component = c("STEP03B3B v2 gate", "Canonical source identities", "Published-table reconciliation", 
    "Chromosome style", "Non-coordinate calcitriol rows", "Primary VDR cistrome sources", "Derived author overlap sources", 
    "chrM handling", "Local peak calling", "RNA membership", "ATAC evidence role", "Definitive direct-target assignment", 
    "Scientific figures", "STEP03B3C authorization", "STEP03B3B v3 gate"), status = c("PASS", "PASS_5_OF_5", "PASS_4_OF_4_EXACT", 
    "LOCKED_CANONICAL_UCSC_chr", "13_SUMMARY_PLUS_3_BLANK_NOT_PEAK_ROWS", "LOCKED_3_OF_3", "LOCKED_2_OF_2_DERIVED_ONLY", 
    "RETAINED_SEPARATELY_EXCLUDED_FROM_PRIMARY_NUCLEAR_JOIN", "DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS", paste0("IMMUTABLE_", 
        RNA_anchor), ATAC_role, "NOT_PERFORMED", "NONE_CREATED", "AUTHORIZED", "READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP"), 
    required_action = c("Retain the v2 repair record.", "Use only hash-locked canonical tables.", "Preserve the exact-match reconciliation record.", 
        "Use chr-prefixed names in all later genomic joins.", "Retain these rows in audit only; do not count them as failed peaks.", 
        "Use media-6 primary sources for project-specific occupancy support.", "Use media-2 tables only as author-derived cross-checks.", 
        "Do not silently mix mitochondrial and nuclear integration results.", "Do not install or run MACS3 for the primary analysis.", 
        "CUT&RUN may add evidence tiers but cannot change RNA membership.", "Use ATAC only as optional accessible context.", 
        "Use the term direct-target candidate only after project-specific integration.", "No Scientific Figure Skill 3 workflow is triggered.", 
        "Audit peak annotation fields before creating gene-level evidence tiers.", "Proceed to STEP03B3C."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3B_v3_GATE_TABLE.csv"))

write_text(c("# Published VDR cistrome source lock", "", paste0("- Source-lock version: ", source_lock_version), 
    paste0("- Run ID: ", run_id), "- Genome build: hg38 as reported by the source article", "- Chromosome style: canonical UCSC chr", 
    "- Primary source 1: media-6 pLentiVeh, 885 valid coordinate rows", "- Primary source 2: media-6 pLenti125, 6,092 valid coordinate rows", 
    "- Primary source 3: media-6 pLentiVeh_125overlap, 519 valid coordinate rows", "- Nuclear primary rows: 884 vehicle, 6,091 calcitriol, 518 shared", 
    "- chrM records: retained separately, one in each primary source", "- Derived media-2 overlap sources: 134 vehicle and 2,022 calcitriol rows", 
    "- Reconciliation: four of four audits are 100% exact matches", "- Non-coordinate appendix: 13 summary rows plus 3 blank separator rows", 
    "- Peak rows removed by appendix classification: 0", "- Local peak calling: deferred; not required for primary analysis", 
    paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE"), paste0("- ATAC role: ", ATAC_role), "- Definitive direct-target assignment: not performed", 
    "- Scientific figures: none", "- Next: STEP03B3C peak-to-gene schema audit and integration preparation"), file.path(dirs$results, 
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK_CARD.md"))

write_text(c("# STEP03B3C authorization", "", paste0("- Authorized by: ", run_id), paste0("- Published source lock: ", 
    source_lock_version), "- Authorization: YES", "- Authorized task 1: audit author-provided peak annotation fields", 
    "- Authorized task 2: build a transparent peak-to-gene mapping table", "- Authorized task 3: prepare project-specific RNA/CUT&RUN evidence tiers", 
    "- RNA membership changes: PROHIBITED", "- Definitive causal or direct-target claims: PROHIBITED", "- media-2 author overlap tables: cross-check only", 
    "- ATAC: optional accessible context only", "- Local MACS peak calling: not required for primary analysis", 
    "- Scientific figures: create only after gene-level evidence tables pass QC"), file.path(dirs$control, "STEP03B3C_AUTHORIZATION.md"))

write_text(c("# STEP03B3B v3 Methods addendum", "", "Following chromosome-style harmonization, all four published-table reconciliation tests showed complete exact coordinate matching. The author-processed vehicle, calcitriol and shared VDR cistrome coordinate tables were therefore locked as the primary chromatin-occupancy sources for subsequent project-specific peak-to-gene mapping.", 
    "", "Sixteen trailing rows in the calcitriol worksheet contained no genomic coordinates. Manual inspection classified 13 as appended annotation-summary rows and three as blank separator rows; none represented failed peak coordinates. These rows were retained in the audit record but excluded from the coordinate source.", 
    "", "One mitochondrial interval was retained in each primary cistrome table. Mitochondrial intervals remain in the source archive but will be reported separately from the primary nuclear-genome integration."), 
    file.path(dirs$text, "STEP03B3B_v3_Methods_addendum.md"))

write_text(c("# STEP03B3B v3 Results draft", "", "Canonical chromosome harmonization resolved the apparent failure of the shared-peak reconciliation. All 518 primary-chromosome shared regions were exact matches in both the vehicle and calcitriol cistrome tables. The two author-derived cistrome-transcriptome overlap tables likewise showed complete exact reconciliation with their corresponding primary cistrome sources (134 of 134 vehicle regions and 2,022 of 2,022 calcitriol regions).", 
    "", "The primary author-processed sources contained 885 valid vehicle regions, 6,092 valid calcitriol regions and 519 valid shared regions. Sixteen additional rows in the calcitriol worksheet were non-coordinate appendix content rather than peak records."), 
    file.path(dirs$text, "STEP03B3B_v3_Results_draft.md"))

parameters <- list(schema_version = "3.0", run_id = run_id, source_lock_version = source_lock_version, RNA_anchor = RNA_anchor, 
    ATAC_role = ATAC_role, primary_sources = 3, derived_sources = 2, exact_reconciliation_audits = 4, non_coordinate_summary_rows = 13, 
    blank_separator_rows = 3, peak_rows_removed = 0, local_peak_calling_required_for_primary_analysis = FALSE, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, STEP03B3C_authorized = TRUE, 
    gate = "READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP")

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3B_v3_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B_v3.txt"))

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

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    safe_file_copy(source_file, .vdr_script_path("35_step03B3B_v3_lock_published_VDR_cistrome_source.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "prepare_vdr_cistrome.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3B_v3_SOURCE_ROLE_SUMMARY.csv"), file.path(dirs$results, "STEP03B3B_v3_GATE_TABLE.csv"), file.path(dirs$results, 
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK_CARD.md"), file.path(dirs$qc, "STEP03B3B_v3_NON_COORDINATE_APPENDIX_ROW_DECISIONS.csv"), 
    file.path(dirs$qc, "STEP03B3B_v3_NON_COORDINATE_APPENDIX_SUMMARY.csv"), file.path(dirs$qc, "STEP03B3B_v3_SOURCE_COUNT_LOCK_CHECK.csv"), 
    file.path(dirs$control, "STEP03B3C_AUTHORIZATION.md"), file.path(dirs$text, "STEP03B3B_v3_Methods_addendum.md"), 
    file.path(dirs$text, "STEP03B3B_v3_Results_draft.md"), file.path(dirs$parameters, "STEP03B3B_v3_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B_v3.txt"), file.path(dirs$code, "prepare_vdr_cistrome.R"), 
    log_path)

for (path_value in archive_artifacts[file.exists(archive_artifacts)]) {
    section <- if (grepl("\\.R$", basename(path_value))) {
        "code"
    }
    else if (grepl("sessionInfo|parameters", basename(path_value))) {
        "env"
    }
    else if (grepl("\\.log$", basename(path_value))) {
        "log"
    }
    else if (grepl("Methods|Results|AUTHORIZATION", basename(path_value))) {
        "txt"
    }
    else if (grepl("SUMMARY|DECISION|CHECK|LOCK_CARD", basename(path_value), ignore.case = TRUE)) {
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
    sha256 = vapply(archive_files, sha256_file, character(1)), path_length = nchar(archive_files), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3B_v3_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3B v3 archive status", "", paste0("- Run ID: ", run_id), paste0("- Source lock: ", source_lock_version), 
    "- Primary VDR cistrome sources: 3", "- Derived author overlap references: 2", "- Exact reconciliation audits: 4/4", 
    "- Non-coordinate appendix rows: 16", "- Peak rows removed: 0", "- Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS", 
    "- STEP03B3C authorization: YES", paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Completed substep: STEP03B3B", paste0("- Accepted run: ", run_id), 
    paste0("- Published VDR cistrome source: ", source_lock_version, " LOCKED"), paste0("- RNA anchor: ", RNA_anchor, 
        " IMMUTABLE"), paste0("- ATAC role: ", ATAC_role), "- Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS", 
    "- Definitive direct-target assignment: NOT_PERFORMED", "- Active next substep: STEP03B3C", "- Scientific figures: NONE", 
    "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate: READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP")

cat("\n============================================================\n")

cat("STEP03B3B v3 published VDR cistrome source lock completed.\n")

cat("Source lock: ", source_lock_version, "\n")

cat("Primary VDR cistrome sources: 3/3 LOCKED\n")

cat("Derived author overlap references: 2/2 LOCKED AS DERIVED\n")

cat("Published reconciliation: 4/4 EXACT\n")

cat("Non-coordinate calcitriol rows: 13 summary + 3 blank\n")

cat("Peak rows removed: 0\n")

cat("Local peak calling: DEFERRED_NOT_REQUIRED_FOR_PRIMARY_ANALYSIS\n")

cat("RNA anchor: ", RNA_anchor, " IMMUTABLE\n")

cat("ATAC role: ", ATAC_role, "\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("STEP03B3C authorization: YES\n")

cat("Gate: READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP\n")

cat("Authoritative archive:\n  ", archive_root, "\n")

cat("============================================================\n")

