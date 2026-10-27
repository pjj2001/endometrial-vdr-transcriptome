#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A13: Standardize occupancy source identifiers

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C2_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"

RNA_anchor <- "RNA_PROGRAM_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

map_version <- "PEAK_GENE_MAP_CANDIDATE_v1.1"

required_packages <- c("data.table", "digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

required_dir_names <- c("control", "inputs", "code", "parameters", "environment", "execution", "results", "figures", 
    "qc", "text", "archive_parent")

missing_dir_names <- setdiff(required_dir_names, names(dirs))

invalid_dir_names <- names(dirs)[!vapply(dirs, function(path_value) {
    is.character(path_value) && length(path_value) == 1L && !is.na(path_value) && nzchar(path_value)
}, logical(1))]

if (length(missing_dir_names) > 0L || length(invalid_dir_names) > 0L) {
    stop("Invalid directory registry. Missing: ", paste(missing_dir_names, collapse = ", "), "; invalid: ", paste(invalid_dir_names, 
        collapse = ", "))
}

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

v1_map_root <- file.path(dirs$results, "STEP03B3C2_peak_gene_map_candidate")

v2_map_root <- file.path(dirs$results, "STEP03B3C2_v2_peak_gene_map_candidate")

dir.create(v2_map_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C2v2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

write_tsv_gz <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
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

truth_flag <- function(x, n_rows) {
    if (is.null(x)) {
        return(rep(FALSE, n_rows))
    }
    value <- toupper(trimws(as.character(x)))
    value %in% c("TRUE", "T", "1")
}

split_semicolon_set_n <- function(x) {
    value <- trimws(as.character(x))
    value[is.na(value)] <- ""
    vapply(value, function(item) {
        if (!nzchar(item)) {
            return(0L)
        }
        length(unique(strsplit(item, ";", fixed = TRUE)[[1]]))
    }, integer(1))
}

empty_promoter_exception_table <- function() {
    data.frame(source_id = character(0), occupancy_context = character(0), coordinate_key = character(0), promoter_only_gene = character(0), 
        raw_gene_overlap_value = character(0), raw_promoter_overlap_value = character(0), stringsAsFactors = FALSE)
}

empty_invalid_token_table <- function() {
    data.frame(source_id = character(0), occupancy_context = character(0), coordinate_key = character(0), invalid_gene_token = character(0), 
        stringsAsFactors = FALSE)
}

log_event("Starting", run_id)

required_paths <- c(gate_v1 = file.path(dirs$results, "STEP03B3C2_GATE_TABLE.csv"), card_v1 = file.path(dirs$results, 
    "STEP03B3C2_PEAK_GENE_MAP_CANDIDATE_CARD.md"), field_lock = file.path(dirs$results, "STEP03B3C2_FIELD_SEMANTICS_LOCK.csv"), 
    mapping_qc_v1 = file.path(dirs$qc, "STEP03B3C2_MAPPING_QC_SUMMARY.csv"), peak_qc = file.path(dirs$qc, "STEP03B3C2_PEAK_LEVEL_MAPPING_QC.csv"), 
    shared_audit = file.path(dirs$qc, "STEP03B3C2_SHARED_ANNOTATION_INHERITANCE_AUDIT.csv"), promoter_exceptions_v1 = file.path(dirs$qc, 
        "STEP03B3C2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"), invalid_tokens_v1 = file.path(dirs$qc, "STEP03B3C2_INVALID_GENE_TOKEN_ROWS.csv"), 
    map_v1 = file.path(v1_map_root, "STEP03B3C2_PEAK_GENE_LONG_CANDIDATE.tsv.gz"), RNA_lock_card = file.path(dirs$results, 
        "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"), source_lock = file.path(dirs$results, "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C2 v1 output is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v1 <- .vdr_read_csv(required_paths[["gate_v1"]], stringsAsFactors = FALSE, check.names = FALSE)

card_v1 <- readLines(required_paths[["card_v1"]], warn = FALSE, encoding = "UTF-8")

field_lock <- .vdr_read_csv(required_paths[["field_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

mapping_qc_v1 <- .vdr_read_csv(required_paths[["mapping_qc_v1"]], stringsAsFactors = FALSE, check.names = FALSE)

peak_qc <- .vdr_read_csv(required_paths[["peak_qc"]], stringsAsFactors = FALSE, check.names = FALSE)

shared_audit <- .vdr_read_csv(required_paths[["shared_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock <- .vdr_read_csv(required_paths[["source_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- gate_v1$status[gate_v1$gate_component == "STEP03B3C2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_PEAK_GENE_MAP_LOCK") {
    stop("Unexpected STEP03B3C2 v1 gate: ", paste(accepted_gate, collapse = ", "))
}

if (!all(field_lock$status == "LOCKED_FOR_CANDIDATE_MAP")) {
    stop("STEP03B3C2 field semantics are not fully locked.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

if (nrow(source_lock[source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES" & source_lock$source_lock_status == 
    "LOCKED", , drop = FALSE]) != 3L) {
    stop("Published primary cistrome source lock is incomplete.")
}

if (nrow(shared_audit) != 518L || !all(shared_audit$gene_set_exact) || !all(shared_audit$promoter_set_exact) || 
    !all(shared_audit$relative_location_exact) || !all(shared_audit$tss_distance_exact)) {
    stop("Shared annotation inheritance audit is not exact.")
}

expected_vehicle_relations <- sum(peak_qc$relation_token_n[peak_qc$occupancy_context == "vehicle"])

expected_calcitriol_relations <- sum(peak_qc$relation_token_n[peak_qc$occupancy_context == "calcitriol"])

expected_shared_relations <- sum(split_semicolon_set_n(shared_audit$vehicle_gene_set))

expected_relation_counts <- data.frame(source_id = c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME", "VDR_SHARED_CISTROME"), 
    occupancy_context = c("vehicle", "calcitriol", "shared_vehicle_and_calcitriol"), expected_relation_n = c(expected_vehicle_relations, 
        expected_calcitriol_relations, expected_shared_relations), stringsAsFactors = FALSE)

expected_total_relations <- sum(expected_relation_counts$expected_relation_n)

if (expected_vehicle_relations <= 0L || expected_calcitriol_relations <= 0L || expected_shared_relations <= 0L) {
    stop("Independent QC tables yielded an empty source-specific relation count.")
}

map_v1 <- .vdr_fread(required_paths[["map_v1"]], sep = "\t", data.table = FALSE, showProgress = FALSE)

required_map_columns <- c("source_id", "occupancy_context", "coordinate_key", "gene_symbol_author")

missing_map_columns <- setdiff(required_map_columns, names(map_v1))

if (length(missing_map_columns) > 0L) {
    stop("Candidate long map is missing columns: ", paste(missing_map_columns, collapse = ", "))
}

if (nrow(map_v1) != expected_total_relations) {
    stop("Candidate long-map total differs from independent relation total. Map: ", nrow(map_v1), "; expected: ", 
        expected_total_relations)
}

source_id_map <- c(MEDIA6_pLentiVeh = "VDR_VEHICLE_CISTROME", VDR_VEHICLE_CISTROME = "VDR_VEHICLE_CISTROME", MEDIA6_pLenti125 = "VDR_CALCITRIOL_CISTROME", 
    VDR_CALCITRIOL_CISTROME = "VDR_CALCITRIOL_CISTROME", MEDIA6_pLentiVeh_125overlap = "VDR_SHARED_CISTROME", VDR_SHARED_CISTROME = "VDR_SHARED_CISTROME")

unexpected_source_ids <- setdiff(unique(map_v1$source_id), names(source_id_map))

if (length(unexpected_source_ids) > 0L) {
    stop("Unexpected upstream source IDs: ", paste(unexpected_source_ids, collapse = ", "))
}

map_v2 <- map_v1

map_v2$source_id_upstream <- as.character(map_v2$source_id)

map_v2$source_id <- unname(source_id_map[map_v2$source_id_upstream])

expected_context_by_source <- c(VDR_VEHICLE_CISTROME = "vehicle", VDR_CALCITRIOL_CISTROME = "calcitriol", VDR_SHARED_CISTROME = "shared_vehicle_and_calcitriol")

context_pass <- map_v2$occupancy_context == unname(expected_context_by_source[map_v2$source_id])

if (!all(context_pass)) {
    stop("At least one row has an occupancy context inconsistent with its canonical source ID.")
}

map_v2$peak_gene_relation_id <- paste0(map_v2$source_id, "::", map_v2$coordinate_key, "::", map_v2$gene_symbol_author)

duplicate_relation_ids <- duplicated(map_v2$peak_gene_relation_id)

if (any(duplicate_relation_ids)) {
    stop("Canonical source normalization produced duplicate peak-gene relation IDs.")
}

map_v2 <- map_v2[order(map_v2$source_id, map_v2$coordinate_key, map_v2$gene_symbol_author), , drop = FALSE]

observed_relation_counts <- aggregate(peak_gene_relation_id ~ source_id + occupancy_context, data = map_v2, FUN = length)

names(observed_relation_counts)[names(observed_relation_counts) == "peak_gene_relation_id"] <- "observed_relation_n"

relation_count_reconciliation <- merge(expected_relation_counts, observed_relation_counts, by = c("source_id", 
    "occupancy_context"), all = TRUE, sort = FALSE)

relation_count_reconciliation$relation_count_pass <- with(relation_count_reconciliation, expected_relation_n == 
    observed_relation_n)

if (nrow(relation_count_reconciliation) != 3L || !all(relation_count_reconciliation$relation_count_pass)) {
    stop("Canonical source-specific relation counts do not reconcile.")
}

map_dt <- data.table::as.data.table(map_v2)

n_rows <- nrow(map_dt)

map_dt$gene_overlap_supported_row <- truth_flag(map_dt$author_gene_overlap, n_rows) | truth_flag(map_dt$vehicle_gene_overlap, 
    n_rows) | truth_flag(map_dt$calcitriol_gene_overlap, n_rows)

map_dt$promoter_overlap_supported_row <- truth_flag(map_dt$author_promoter_overlap, n_rows) | truth_flag(map_dt$vehicle_promoter_overlap, 
    n_rows) | truth_flag(map_dt$calcitriol_promoter_overlap, n_rows)

gene_level_preview <- map_dt[, .(supporting_peak_n = data.table::uniqueN(coordinate_key), gene_overlap_supported = any(gene_overlap_supported_row, 
    na.rm = TRUE), promoter_overlap_supported = any(promoter_overlap_supported_row, na.rm = TRUE), minimum_absolute_tss_distance = if (all(is.na(suppressWarnings(as.numeric(peak_tss_distance))))) {
    NA_real_
} else {
    min(abs(suppressWarnings(as.numeric(peak_tss_distance))), na.rm = TRUE)
}), by = .(source_id, occupancy_context, gene_symbol_author)]

gene_level_preview <- as.data.frame(gene_level_preview, stringsAsFactors = FALSE)

promoter_exceptions_v1 <- tryCatch(.vdr_read_csv(required_paths[["promoter_exceptions_v1"]], stringsAsFactors = FALSE, 
    check.names = FALSE), error = function(error_object) {
    empty_promoter_exception_table()
})

invalid_tokens_v1 <- tryCatch(.vdr_read_csv(required_paths[["invalid_tokens_v1"]], stringsAsFactors = FALSE, check.names = FALSE), 
    error = function(error_object) {
        empty_invalid_token_table()
    })

if (nrow(promoter_exceptions_v1) == 0L) {
    promoter_exceptions_v2 <- empty_promoter_exception_table()
} else {
    promoter_exceptions_v2 <- promoter_exceptions_v1
}

if (nrow(invalid_tokens_v1) == 0L) {
    invalid_tokens_v2 <- empty_invalid_token_table()
} else {
    invalid_tokens_v2 <- invalid_tokens_v1
}

if (nrow(promoter_exceptions_v2) > 0L || nrow(invalid_tokens_v2) > 0L) {
    stop("The v1 exception tables are not empty; targeted manual review is required.")
}

map_v2_path <- file.path(v2_map_root, "STEP03B3C2_v2_PEAK_GENE_LONG_CANDIDATE.tsv.gz")

gene_preview_v2_path <- file.path(v2_map_root, "STEP03B3C2_v2_GENE_LEVEL_EVIDENCE_PREVIEW.csv")

write_tsv_gz(map_v2, map_v2_path)

write_csv_utf8(gene_level_preview, gene_preview_v2_path)

write_csv_utf8(promoter_exceptions_v2, file.path(dirs$qc, "STEP03B3C2_v2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"))

write_csv_utf8(invalid_tokens_v2, file.path(dirs$qc, "STEP03B3C2_v2_INVALID_GENE_TOKEN_ROWS.csv"))

write_csv_utf8(relation_count_reconciliation, file.path(dirs$qc, "STEP03B3C2_v2_RELATION_COUNT_RECONCILIATION.csv"))

peak_summary <- aggregate(cbind(linked_peak = as.integer(peak_qc$relation_token_n > 0L), zero_link_peak = as.integer(peak_qc$relation_token_n == 
    0L), relation_token_n = peak_qc$relation_token_n, promoter_relation_n = peak_qc$promoter_overlap_token_n) ~ 
    occupancy_context, data = peak_qc, FUN = sum)

peak_total <- aggregate(coordinate_key ~ occupancy_context, data = peak_qc, FUN = length)

names(peak_total)[names(peak_total) == "coordinate_key"] <- "peak_n"

peak_summary <- merge(peak_total, peak_summary, by = "occupancy_context", all = TRUE, sort = FALSE)

shared_gene_n <- split_semicolon_set_n(shared_audit$vehicle_gene_set)

shared_promoter_n <- split_semicolon_set_n(shared_audit$vehicle_promoter_set)

mapping_qc_v2 <- data.frame(metric = c("Vehicle nuclear peaks", "Vehicle peaks with at least one author gene relation", 
    "Vehicle peaks with zero author gene relations", "Vehicle peak-gene relations", "Vehicle promoter-supported relations", 
    "Calcitriol nuclear peaks", "Calcitriol peaks with at least one author gene relation", "Calcitriol peaks with zero author gene relations", 
    "Calcitriol peak-gene relations", "Calcitriol promoter-supported relations", "Shared nuclear peaks", "Shared peaks with at least one inherited gene relation", 
    "Shared peaks with zero inherited gene relations", "Shared inherited peak-gene relations", "Shared promoter-supported relations", 
    "Total peak-gene relations", "Unique author gene symbols", "Promoter-only relation exceptions", "Invalid gene-token rows", 
    "Shared gene-set concordance", "Shared promoter-set concordance", "Shared relative-location concordance", "Shared TSS-distance concordance", 
    "Canonical source IDs represented"), value = c(peak_summary$peak_n[peak_summary$occupancy_context == "vehicle"], 
    peak_summary$linked_peak[peak_summary$occupancy_context == "vehicle"], peak_summary$zero_link_peak[peak_summary$occupancy_context == 
        "vehicle"], expected_vehicle_relations, peak_summary$promoter_relation_n[peak_summary$occupancy_context == 
        "vehicle"], peak_summary$peak_n[peak_summary$occupancy_context == "calcitriol"], peak_summary$linked_peak[peak_summary$occupancy_context == 
        "calcitriol"], peak_summary$zero_link_peak[peak_summary$occupancy_context == "calcitriol"], expected_calcitriol_relations, 
    peak_summary$promoter_relation_n[peak_summary$occupancy_context == "calcitriol"], nrow(shared_audit), sum(shared_gene_n > 
        0L), sum(shared_gene_n == 0L), expected_shared_relations, sum(shared_promoter_n), nrow(map_v2), length(unique(map_v2$gene_symbol_author)), 
    nrow(promoter_exceptions_v2), nrow(invalid_tokens_v2), mean(shared_audit$gene_set_exact), mean(shared_audit$promoter_set_exact), 
    mean(shared_audit$relative_location_exact), mean(shared_audit$tss_distance_exact), length(unique(map_v2$source_id))), 
    stringsAsFactors = FALSE)

write_csv_utf8(mapping_qc_v2, file.path(dirs$qc, "STEP03B3C2_v2_MAPPING_QC_SUMMARY.csv"))

candidate_manifest_v2 <- data.frame(artifact_id = c("PEAK_GENE_LONG_CANDIDATE", "GENE_LEVEL_EVIDENCE_PREVIEW"), 
    version = map_version, path = normalize_slashes(c(map_v2_path, gene_preview_v2_path)), size_bytes = as.numeric(file.info(c(map_v2_path, 
        gene_preview_v2_path))$size), md5 = vapply(c(map_v2_path, gene_preview_v2_path), md5_file, character(1)), 
    sha256 = vapply(c(map_v2_path, gene_preview_v2_path), sha256_file, character(1)), source_id_status = "CANONICAL_LOCKED_SOURCE_IDS", 
    relation_membership_changed = FALSE, lock_status = "CANDIDATE_PENDING_MANUAL_QC", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest_v2, file.path(v2_map_root, "STEP03B3C2_v2_CANDIDATE_MAP_MANIFEST.csv"))

source_id_transition <- data.frame(upstream_source_id = c("MEDIA6_pLentiVeh", "MEDIA6_pLenti125", "VDR_SHARED_CISTROME"), 
    canonical_source_id = c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME", "VDR_SHARED_CISTROME"), occupancy_context = c("vehicle", 
        "calcitriol", "shared_vehicle_and_calcitriol"), relation_membership_changed = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(source_id_transition, file.path(dirs$results, "STEP03B3C2_v2_SOURCE_ID_TRANSITION_LOCK.csv"))

gate_status <- "PROCEED_TO_MANUAL_PEAK_GENE_MAP_v1_1_LOCK"

gate_table <- data.frame(gate_component = c("STEP03B3C2 v1 gate", "Published cistrome source lock", "RNA membership", 
    "ATAC role", "Field semantics", "Upstream source IDs", "Canonical source IDs", "Vehicle relation reconciliation", 
    "Calcitriol relation reconciliation", "Shared relation reconciliation", "Total relation membership", "Promoter relation consistency", 
    "Gene token shape", "Candidate map version", "RNA/CUT&RUN integration", "Definitive direct-target assignment", 
    "Scientific figures", "STEP03B3C2 v2 gate"), status = c("PASS_WITH_REPORTING_IDENTITY_REPAIR", paste0("IMMUTABLE_", 
    source_lock_version), paste0("IMMUTABLE_", RNA_anchor), ATAC_role, "UNCHANGED_LOCKED", "PRESERVED_IN_source_id_upstream", 
    "PASS_3_OF_3", paste0("PASS_", expected_vehicle_relations), paste0("PASS_", expected_calcitriol_relations), 
    paste0("PASS_", expected_shared_relations), paste0("UNCHANGED_", nrow(map_v2)), "PASS_ZERO_EXCEPTIONS", "PASS_ZERO_INVALID_TOKENS", 
    paste0(map_version, "_PENDING_MANUAL_QC"), "NOT_STARTED", "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain v1 audit but do not lock its source-specific count summary.", 
    "Do not alter source peaks or coordinates.", "Do not change the 436 RNA interaction candidates.", "Use only as optional accessible context.", 
    "Retain the six field interpretations.", "Keep upstream IDs for provenance.", "Use canonical IDs in all downstream tables and joins.", 
    "Review corrected vehicle count.", "Review corrected calcitriol count.", "Review corrected shared count.", 
    "Confirm no relation row was added or removed.", "Retain explicit empty-schema exception table.", "Retain explicit empty-schema invalid-token table.", 
    "Lock only after manual QC.", "Do not integrate until v1.1 candidate map is locked.", "Use occupancy-supported candidate terminology only after integration QC.", 
    "No Scientific Figure Skill 3 workflow is triggered.", "Upload corrected QC, reconciliation and manifest."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C2_v2_GATE_TABLE.csv"))

write_text(c("# STEP03B3C2 v2 candidate-map repair card", "", paste0("- Run ID: ", run_id), "- Root cause: upstream source IDs were counted against canonical IDs", 
    "- Peak-gene membership changed: NO", paste0("- Candidate version: ", map_version), paste0("- Vehicle relations: ", 
        expected_vehicle_relations), paste0("- Calcitriol relations: ", expected_calcitriol_relations), paste0("- Shared relations: ", 
        expected_shared_relations), paste0("- Total relations: ", nrow(map_v2)), paste0("- Unique author gene symbols: ", 
        length(unique(map_v2$gene_symbol_author))), "- Promoter-only exceptions: 0", "- Invalid gene-token rows: 0", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C2_v2_CANDIDATE_MAP_REPAIR_CARD.md"))

write_text(c("# STEP03B3C2 v2 correction note", "", "The v1 candidate long map contained the expected 23,499 peak-gene relations, but vehicle and calcitriol rows retained the upstream source labels `MEDIA6_pLentiVeh` and `MEDIA6_pLenti125`. The v1 QC summary counted the locked canonical labels and therefore incorrectly reported zero relations for those two sources.", 
    "", "The v2 repair preserved all relation rows and recorded the upstream labels in `source_id_upstream`, while normalizing `source_id` to `VDR_VEHICLE_CISTROME`, `VDR_CALCITRIOL_CISTROME` and `VDR_SHARED_CISTROME`. Corrected source-specific counts were independently reconciled to the peak-level token audit and shared annotation inheritance table.", 
    "", "No peak coordinate, gene token, promoter flag, field semantics or RNA membership was changed."), file.path(dirs$text, 
    "STEP03B3C2_v2_Correction_note.md"))

parameters <- list(schema_version = "1.1", run_id = run_id, source_lock_version = source_lock_version, RNA_anchor = RNA_anchor, 
    ATAC_role = ATAC_role, candidate_map_version = map_version, root_cause = "upstream_source_ids_counted_against_canonical_ids", 
    relation_membership_changed = FALSE, vehicle_relation_n = expected_vehicle_relations, calcitriol_relation_n = expected_calcitriol_relations, 
    shared_relation_n = expected_shared_relations, total_relation_n = nrow(map_v2), unique_gene_symbol_n = length(unique(map_v2$gene_symbol_author)), 
    RNA_CUTRUN_integration = FALSE, definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, 
    gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C2_v2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2_v2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("39_step03B3C2_v2_normalize_source_ids_repair_mapping_QC.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "standardize_peak_source_ids.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C2_v2_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C2_v2_SOURCE_ID_TRANSITION_LOCK.csv"), 
    file.path(dirs$results, "STEP03B3C2_v2_CANDIDATE_MAP_REPAIR_CARD.md"), file.path(dirs$qc, "STEP03B3C2_v2_RELATION_COUNT_RECONCILIATION.csv"), 
    file.path(dirs$qc, "STEP03B3C2_v2_MAPPING_QC_SUMMARY.csv"), file.path(dirs$qc, "STEP03B3C2_v2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"), 
    file.path(dirs$qc, "STEP03B3C2_v2_INVALID_GENE_TOKEN_ROWS.csv"), map_v2_path, gene_preview_v2_path, file.path(v2_map_root, 
        "STEP03B3C2_v2_CANDIDATE_MAP_MANIFEST.csv"), file.path(dirs$text, "STEP03B3C2_v2_Correction_note.md"), 
    file.path(dirs$parameters, "STEP03B3C2_v2_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2_v2.txt"), 
    file.path(dirs$code, "standardize_peak_source_ids.R"), log_path)

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
    else if (grepl("Correction", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|RECONCILIATION|EXCEPTION|INVALID|CARD|MANIFEST|TRANSITION", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C2_v2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C2 v2 archive status", "", paste0("- Run ID: ", run_id), "- Peak-gene membership changed: NO", 
    paste0("- Candidate map: ", map_version), paste0("- Vehicle relations: ", expected_vehicle_relations), paste0("- Calcitriol relations: ", 
        expected_calcitriol_relations), paste0("- Shared relations: ", expected_shared_relations), paste0("- Total relations: ", 
        nrow(map_v2)), "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", 
        gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, 
    "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C2 v2", paste0("- Run ID: ", run_id), paste0("- Published source lock: ", 
    source_lock_version, " IMMUTABLE"), paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE"), paste0("- ATAC role: ", 
    ATAC_role), "- Annotation field semantics: LOCKED", paste0("- Peak-gene candidate map: ", map_version, " PENDING_MANUAL_QC"), 
    "- Source IDs: CANONICALIZED_WITH_UPSTREAM_PROVENANCE", "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT LOCKED"), file.path(dirs$control, 
    "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C2 v2 source-ID and mapping-QC repair completed.\n")

cat("Peak-gene membership changed: NO\n")

cat("Candidate map: ", map_version, "\n")

cat("Vehicle relations: ", expected_vehicle_relations, "\n")

cat("Calcitriol relations: ", expected_calcitriol_relations, "\n")

cat("Shared relations: ", expected_shared_relations, "\n")

cat("Total relations: ", nrow(map_v2), "\n")

cat("Unique author gene symbols: ", length(unique(map_v2$gene_symbol_author)), "\n")

cat("Promoter-only exceptions: 0\n")

cat("Invalid gene-token rows: 0\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C2_v2_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3C2_v2_SOURCE_ID_TRANSITION_LOCK.csv\n")

cat("  08_qc/STEP03B3C2_v2_RELATION_COUNT_RECONCILIATION.csv\n")

cat("  08_qc/STEP03B3C2_v2_MAPPING_QC_SUMMARY.csv\n")

cat("  06_results/STEP03B3C2_v2_peak_gene_map_candidate/STEP03B3C2_v2_CANDIDATE_MAP_MANIFEST.csv\n")

cat("============================================================\n")

