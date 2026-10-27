#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A14: Export the accepted peak-gene map

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C2_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"

RNA_anchor <- "RNA_PROGRAM_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

expected_relation_counts <- c(VDR_VEHICLE_CISTROME = 2848L, VDR_CALCITRIOL_CISTROME = 19037L, VDR_SHARED_CISTROME = 1614L)

expected_total_relations <- 23499L

expected_unique_gene_symbols <- 10042L

expected_primary_source_n <- 3L

expected_shared_peak_n <- 518L

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

candidate_root <- file.path(dirs$results, "STEP03B3C2_v2_peak_gene_map_candidate")

locked_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

dir.create(locked_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C2v3_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

parse_bool <- function(x) {
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

read_empty_schema_csv <- function(path, expected_columns) {
    table_value <- .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(table_value) != 0L || !identical(names(table_value), expected_columns)) {
        stop("Zero-row exception table failed schema validation: ", path)
    }
    table_value
}

log_event("Starting", run_id)

required_paths <- c(gate_v2 = file.path(dirs$results, "STEP03B3C2_v2_GATE_TABLE.csv"), transition_lock = file.path(dirs$results, 
    "STEP03B3C2_v2_SOURCE_ID_TRANSITION_LOCK.csv"), field_semantics = file.path(dirs$results, "STEP03B3C2_FIELD_SEMANTICS_LOCK.csv"), 
    relation_reconciliation = file.path(dirs$qc, "STEP03B3C2_v2_RELATION_COUNT_RECONCILIATION.csv"), mapping_qc = file.path(dirs$qc, 
        "STEP03B3C2_v2_MAPPING_QC_SUMMARY.csv"), promoter_exceptions = file.path(dirs$qc, "STEP03B3C2_v2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"), 
    invalid_tokens = file.path(dirs$qc, "STEP03B3C2_v2_INVALID_GENE_TOKEN_ROWS.csv"), candidate_manifest = file.path(candidate_root, 
        "STEP03B3C2_v2_CANDIDATE_MAP_MANIFEST.csv"), candidate_long_map = file.path(candidate_root, "STEP03B3C2_v2_PEAK_GENE_LONG_CANDIDATE.tsv.gz"), 
    candidate_gene_preview = file.path(candidate_root, "STEP03B3C2_v2_GENE_LEVEL_EVIDENCE_PREVIEW.csv"), source_lock = file.path(dirs$results, 
        "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"), RNA_lock_card = file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C2 v3 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v2 <- .vdr_read_csv(required_paths[["gate_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

transition_lock <- .vdr_read_csv(required_paths[["transition_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

field_semantics <- .vdr_read_csv(required_paths[["field_semantics"]], stringsAsFactors = FALSE, check.names = FALSE)

relation_reconciliation <- .vdr_read_csv(required_paths[["relation_reconciliation"]], stringsAsFactors = FALSE, 
    check.names = FALSE)

mapping_qc <- .vdr_read_csv(required_paths[["mapping_qc"]], stringsAsFactors = FALSE, check.names = FALSE)

candidate_manifest <- .vdr_read_csv(required_paths[["candidate_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock <- .vdr_read_csv(required_paths[["source_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- gate_v2$status[gate_v2$gate_component == "STEP03B3C2 v2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_PEAK_GENE_MAP_v1_1_LOCK") {
    stop("STEP03B3C2 v2 does not authorize map lock. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (nrow(field_semantics) != 6L || !all(field_semantics$status == "LOCKED_FOR_CANDIDATE_MAP")) {
    stop("The six annotation-field semantics are not fully locked.")
}

if (nrow(transition_lock) != 3L || !setequal(transition_lock$canonical_source_id, names(expected_relation_counts)) || 
    any(parse_bool(transition_lock$relation_membership_changed))) {
    stop("Source-ID transition lock failed.")
}

primary_source_lock <- source_lock[source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", 
    , drop = FALSE]

if (nrow(primary_source_lock) != expected_primary_source_n || !all(primary_source_lock$source_lock_version == source_lock_version) || 
    !all(primary_source_lock$source_lock_status == "LOCKED")) {
    stop("Published VDR cistrome source lock is incomplete.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

promoter_exception_columns <- c("source_id", "occupancy_context", "coordinate_key", "promoter_only_gene", "raw_gene_overlap_value", 
    "raw_promoter_overlap_value")

invalid_token_columns <- c("source_id", "occupancy_context", "coordinate_key", "invalid_gene_token")

promoter_exceptions <- read_empty_schema_csv(required_paths[["promoter_exceptions"]], promoter_exception_columns)

invalid_tokens <- read_empty_schema_csv(required_paths[["invalid_tokens"]], invalid_token_columns)

if (nrow(relation_reconciliation) != 3L || !setequal(relation_reconciliation$source_id, names(expected_relation_counts)) || 
    !all(relation_reconciliation$relation_count_pass)) {
    stop("Source-specific relation reconciliation is incomplete.")
}

relation_reconciliation$expected_locked <- unname(expected_relation_counts[relation_reconciliation$source_id])

if (!all(relation_reconciliation$expected_relation_n == relation_reconciliation$expected_locked) || !all(relation_reconciliation$observed_relation_n == 
    relation_reconciliation$expected_locked)) {
    stop("Source-specific relation counts differ from the reviewed lock.")
}

qc_value <- setNames(mapping_qc$value, mapping_qc$metric)

required_qc_values <- c(`Vehicle peak-gene relations` = 2848, `Calcitriol peak-gene relations` = 19037, `Shared inherited peak-gene relations` = 1614, 
    `Total peak-gene relations` = 23499, `Unique author gene symbols` = 10042, `Promoter-only relation exceptions` = 0, 
    `Invalid gene-token rows` = 0, `Shared gene-set concordance` = 1, `Shared promoter-set concordance` = 1, `Shared relative-location concordance` = 1, 
    `Shared TSS-distance concordance` = 1, `Canonical source IDs represented` = 3)

missing_qc_metrics <- setdiff(names(required_qc_values), names(qc_value))

if (length(missing_qc_metrics) > 0L) {
    stop("Required QC metrics are missing: ", paste(missing_qc_metrics, collapse = ", "))
}

observed_qc_values <- as.numeric(qc_value[names(required_qc_values)])

if (!all(observed_qc_values == as.numeric(required_qc_values))) {
    stop("Mapping QC values differ from the reviewed lock.")
}

if (as.numeric(qc_value[["Shared nuclear peaks"]]) != expected_shared_peak_n) {
    stop("Shared nuclear peak count differs from 518.")
}

if (nrow(candidate_manifest) != 2L || !setequal(candidate_manifest$artifact_id, c("PEAK_GENE_LONG_CANDIDATE", "GENE_LEVEL_EVIDENCE_PREVIEW")) || 
    !all(candidate_manifest$version == "PEAK_GENE_MAP_CANDIDATE_v1.1") || !all(candidate_manifest$source_id_status == 
    "CANONICAL_LOCKED_SOURCE_IDS") || any(parse_bool(candidate_manifest$relation_membership_changed))) {
    stop("Candidate-map manifest metadata failed validation.")
}

candidate_manifest$file_exists <- file.exists(candidate_manifest$path)

candidate_manifest$observed_size_bytes <- as.numeric(file.info(candidate_manifest$path)$size)

candidate_manifest$observed_md5 <- vapply(candidate_manifest$path, md5_file, character(1))

candidate_manifest$observed_sha256 <- vapply(candidate_manifest$path, sha256_file, character(1))

candidate_manifest$identity_pass <- with(candidate_manifest, file_exists & size_bytes == observed_size_bytes & 
    md5 == observed_md5 & sha256 == observed_sha256)

if (!all(candidate_manifest$identity_pass)) {
    stop("One or more candidate-map artifacts failed hash revalidation.")
}

write_csv_utf8(candidate_manifest, file.path(dirs$qc, "STEP03B3C2_v3_CANDIDATE_ARTIFACT_IDENTITY_RECHECK.csv"))

peak_gene_map <- .vdr_fread(required_paths[["candidate_long_map"]], sep = "\t", data.table = FALSE, showProgress = FALSE)

required_map_columns <- c("source_id", "source_id_upstream", "occupancy_context", "coordinate_key", "gene_symbol_author", 
    "peak_gene_relation_id")

missing_map_columns <- setdiff(required_map_columns, names(peak_gene_map))

if (length(missing_map_columns) > 0L) {
    stop("Peak-gene map is missing columns: ", paste(missing_map_columns, collapse = ", "))
}

if (nrow(peak_gene_map) != expected_total_relations) {
    stop("Peak-gene map row count differs from 23,499.")
}

if (any(is.na(peak_gene_map$source_id) | !nzchar(trimws(peak_gene_map$source_id))) || any(is.na(peak_gene_map$coordinate_key) | 
    !nzchar(trimws(peak_gene_map$coordinate_key))) || any(is.na(peak_gene_map$gene_symbol_author) | !nzchar(trimws(peak_gene_map$gene_symbol_author)))) {
    stop("Peak-gene map contains missing required identifiers.")
}

if (length(unique(peak_gene_map$gene_symbol_author)) != expected_unique_gene_symbols) {
    stop("Unique author gene-symbol count differs from 10,042.")
}

if (anyDuplicated(peak_gene_map$peak_gene_relation_id) > 0L) {
    stop("Peak-gene relation IDs are not unique.")
}

expected_relation_id <- paste0(peak_gene_map$source_id, "::", peak_gene_map$coordinate_key, "::", peak_gene_map$gene_symbol_author)

if (!identical(as.character(peak_gene_map$peak_gene_relation_id), as.character(expected_relation_id))) {
    stop("At least one peak-gene relation ID is inconsistent with its fields.")
}

observed_source_counts <- table(peak_gene_map$source_id)

if (!setequal(names(observed_source_counts), names(expected_relation_counts)) || !all(as.integer(observed_source_counts[names(expected_relation_counts)]) == 
    as.integer(expected_relation_counts))) {
    stop("Deep source-specific relation counts failed.")
}

context_by_source <- c(VDR_VEHICLE_CISTROME = "vehicle", VDR_CALCITRIOL_CISTROME = "calcitriol", VDR_SHARED_CISTROME = "shared_vehicle_and_calcitriol")

if (!all(peak_gene_map$occupancy_context == unname(context_by_source[peak_gene_map$source_id]))) {
    stop("Occupancy contexts are inconsistent with canonical source IDs.")
}

upstream_by_source <- c(VDR_VEHICLE_CISTROME = "MEDIA6_pLentiVeh", VDR_CALCITRIOL_CISTROME = "MEDIA6_pLenti125", 
    VDR_SHARED_CISTROME = "VDR_SHARED_CISTROME")

if (!all(peak_gene_map$source_id_upstream == unname(upstream_by_source[peak_gene_map$source_id]))) {
    stop("Upstream source provenance is inconsistent.")
}

gene_preview <- .vdr_read_csv(required_paths[["candidate_gene_preview"]], stringsAsFactors = FALSE, check.names = FALSE)

required_preview_columns <- c("source_id", "occupancy_context", "gene_symbol_author", "supporting_peak_n", "gene_overlap_supported", 
    "promoter_overlap_supported")

missing_preview_columns <- setdiff(required_preview_columns, names(gene_preview))

if (length(missing_preview_columns) > 0L) {
    stop("Gene-level preview is missing columns: ", paste(missing_preview_columns, collapse = ", "))
}

preview_key <- paste(gene_preview$source_id, gene_preview$occupancy_context, gene_preview$gene_symbol_author, sep = "::")

if (anyDuplicated(preview_key) > 0L || any(gene_preview$supporting_peak_n < 1L)) {
    stop("Gene-level evidence preview failed uniqueness or supporting-peak checks.")
}

locked_long_map_path <- file.path(locked_root, "PEAK_GENE_MAP_v1.1_peak_gene_long.tsv.gz")

locked_gene_preview_path <- file.path(locked_root, "PEAK_GENE_MAP_v1.1_gene_level_evidence_preview.csv")

copy_validate(required_paths[["candidate_long_map"]], locked_long_map_path)

copy_validate(required_paths[["candidate_gene_preview"]], locked_gene_preview_path)

locked_manifest <- data.frame(artifact_id = c("PEAK_GENE_LONG", "GENE_LEVEL_EVIDENCE_PREVIEW"), lock_version = peak_gene_map_lock_version, 
    path = normalize_slashes(c(locked_long_map_path, locked_gene_preview_path)), size_bytes = as.numeric(file.info(c(locked_long_map_path, 
        locked_gene_preview_path))$size), md5 = vapply(c(locked_long_map_path, locked_gene_preview_path), md5_file, 
        character(1)), sha256 = vapply(c(locked_long_map_path, locked_gene_preview_path), sha256_file, character(1)), 
    source_lock_version = source_lock_version, RNA_anchor = RNA_anchor, relation_membership_status = "FROZEN_23499_RELATIONS", 
    lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(locked_manifest, file.path(locked_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"))

source_count_lock <- data.frame(source_id = names(expected_relation_counts), occupancy_context = unname(context_by_source[names(expected_relation_counts)]), 
    relation_n = as.integer(expected_relation_counts), map_total_relation_n = expected_total_relations, unique_author_gene_symbol_n = expected_unique_gene_symbols, 
    relation_membership_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(source_count_lock, file.path(dirs$results, "STEP03B3C2_v3_SOURCE_RELATION_COUNT_LOCK.csv"))

map_lock_record <- data.frame(lock_component = c("Lock version", "Long-map rows", "Unique author gene symbols", 
    "Canonical source IDs", "Vehicle relations", "Calcitriol relations", "Shared relations", "Promoter-only exceptions", 
    "Invalid gene-token rows", "Field semantics", "RNA membership", "ATAC role", "Direct-target status", "Scientific figures"), 
    locked_value = c(peak_gene_map_lock_version, as.character(expected_total_relations), as.character(expected_unique_gene_symbols), 
        as.character(expected_primary_source_n), as.character(expected_relation_counts[["VDR_VEHICLE_CISTROME"]]), 
        as.character(expected_relation_counts[["VDR_CALCITRIOL_CISTROME"]]), as.character(expected_relation_counts[["VDR_SHARED_CISTROME"]]), 
        "0", "0", "LOCKED_SIX_FIELD_INTERPRETATIONS", paste0(RNA_anchor, "_IMMUTABLE_436_CANDIDATES"), ATAC_role, 
        "NOT_PERFORMED", "NONE_CREATED"), lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(map_lock_record, file.path(dirs$results, "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK.csv"))

gate_status <- "READY_FOR_STEP03B3C3_RNA_GENE_ID_HARMONIZATION_AND_OCCUPANCY_INTEGRATION"

gate_table <- data.frame(gate_component = c("STEP03B3C2 v2 gate", "Published cistrome source lock", "Field semantics", 
    "Canonical source IDs", "Relation-count reconciliation", "Peak-gene relation IDs", "Peak-gene membership", 
    "Unique author gene symbols", "Promoter relation consistency", "Gene token shape", "Gene-level preview", "Peak-gene map lock", 
    "RNA membership", "ATAC role", "RNA/CUT&RUN integration", "Definitive direct-target assignment", "Scientific figures", 
    "STEP03B3C3 authorization", "STEP03B3C2 v3 gate"), status = c("PASS", paste0("IMMUTABLE_", source_lock_version), 
    "LOCKED_SIX_OF_SIX", "PASS_3_OF_3", "PASS_3_OF_3", "UNIQUE_23499_OF_23499", "FROZEN_23499", "LOCKED_10042", 
    "PASS_ZERO_EXCEPTIONS", "PASS_ZERO_INVALID_TOKENS", "PASS_UNIQUE_SOURCE_CONTEXT_GENE_KEYS", paste0("LOCKED_", 
        peak_gene_map_lock_version), paste0("IMMUTABLE_", RNA_anchor, "_436_CANDIDATES"), ATAC_role, "NOT_STARTED", 
    "NOT_PERFORMED", "NONE_CREATED", "AUTHORIZED", gate_status), required_action = c("Retain the repaired v2 audit record.", 
    "Do not alter source peaks or coordinates.", "Do not reinterpret the six fields without a new version.", "Use canonical IDs in every downstream join.", 
    "Retain source-specific counts.", "Use locked relation IDs as row authority.", "Do not add or remove peak-gene relations.", 
    "Audit RNA-to-cistrome symbol matching next.", "Retain explicit zero-row exception table.", "Retain explicit zero-row invalid-token table.", 
    "Use as a summary only; long map remains authoritative.", "Use only the locked versioned files and hashes.", 
    "Do not change the RNA interaction set.", "Use only as optional accessible context.", "Begin only in STEP03B3C3 after gene-ID audit.", 
    "Use occupancy-supported candidate terminology after integration QC.", "No Scientific Figure Skill 3 workflow is triggered.", 
    "Proceed to RNA gene-ID harmonization and evidence integration.", "STEP03B3C2 is substep-locked."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C2_v3_GATE_TABLE.csv"))

write_text(c("# PEAK_GENE_MAP_v1.1 lock card", "", paste0("- Run ID: ", run_id), paste0("- Lock version: ", peak_gene_map_lock_version), 
    paste0("- Upstream cistrome source: ", source_lock_version, " IMMUTABLE"), "- Peak-gene relations: 23,499", 
    "- Unique author gene symbols: 10,042", "- Vehicle relations: 2,848", "- Calcitriol relations: 19,037", "- Shared inherited relations: 1,614", 
    "- Canonical source IDs: 3/3", "- Promoter-only exceptions: 0", "- Invalid gene-token rows: 0", "- Shared annotation concordance: 100%", 
    "- Relation IDs: unique 23,499/23,499", paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE; 436 interaction candidates"), 
    paste0("- ATAC role: ", ATAC_role), "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", "- Next: STEP03B3C3 RNA gene-ID harmonization and occupancy integration"), file.path(dirs$results, 
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK_CARD.md"))

write_text(c("# STEP03B3C3 authorization", "", paste0("- Authorized by: ", run_id), paste0("- Peak-gene map: ", 
    peak_gene_map_lock_version, " LOCKED"), paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE"), "- Authorization: YES", 
    "- Authorized task 1: identify the authoritative 436-gene RNA interaction table and revalidate its hash.", 
    "- Authorized task 2: audit exact RNA-symbol matches to the locked 10,042-symbol cistrome universe.", "- Authorized task 3: preserve unmatched and ambiguous symbols explicitly; do not silently alias-map.", 
    "- Authorized task 4: create gene-level occupancy evidence by vehicle, calcitriol and shared context.", "- Authorized task 5: retain promoter support, supporting-peak counts and relation provenance separately.", 
    "- Authorized task 6: construct RNA-anchored occupancy-supported candidate evidence tiers.", "- RNA membership changes: PROHIBITED", 
    "- Peak-gene relation membership changes: PROHIBITED", "- media-2 derived overlap tables: CROSS-CHECK ONLY", 
    "- ATAC: OPTIONAL ACCESSIBLE_CONTEXT ONLY", "- Definitive causal/direct-target claims: PROHIBITED", "- Scientific figures: create only after integrated evidence tables pass content and visual QC"), 
    file.path(dirs$control, "STEP03B3C3_AUTHORIZATION.md"))

write_text(c("# STEP03B3C2 v3 Methods addendum", "", "After correction of the source-identifier reporting mismatch, the peak-by-gene candidate map was revalidated against source-specific relation counts, unique relation identifiers and explicit exception schemas. The accepted map comprised 23,499 unique peak-gene relations across three canonical occupancy contexts and 10,042 unique author-reported gene symbols.", 
    "", "The validated long-form map and gene-level evidence preview were copied into a versioned locked directory and recorded with MD5 and SHA-256 hashes. Peak coordinates, gene tokens, promoter flags, field semantics and RNA interaction membership were unchanged."), 
    file.path(dirs$text, "STEP03B3C2_v3_Methods_addendum.md"))

write_text(c("# STEP03B3C2 v3 Results draft", "", "The locked peak-gene map contained 2,848 vehicle, 19,037 calcitriol and 1,614 shared occupancy-supported relations, for a total of 23,499 unique peak-gene relations involving 10,042 author-reported gene symbols. No promoter-only relation exceptions or malformed gene tokens were detected.", 
    "", "All relation identifiers were unique and source-specific counts reconciled exactly with the independent peak-level token audit. The shared annotations retained complete concordance between vehicle and calcitriol source tables."), 
    file.path(dirs$text, "STEP03B3C2_v3_Results_draft.md"))

parameters <- list(schema_version = "1.1-lock", run_id = run_id, source_lock_version = source_lock_version, peak_gene_map_lock_version = peak_gene_map_lock_version, 
    RNA_anchor = RNA_anchor, RNA_interaction_candidate_n = 436, ATAC_role = ATAC_role, vehicle_relation_n = 2848, 
    calcitriol_relation_n = 19037, shared_relation_n = 1614, total_relation_n = 23499, unique_author_gene_symbol_n = 10042, 
    relation_ids_unique = TRUE, promoter_only_exceptions = 0, invalid_gene_tokens = 0, RNA_CUTRUN_integration = FALSE, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, STEP03B3C3_authorized = TRUE, 
    gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C2_v3_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2_v3.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("40_step03B3C2_v3_lock_PEAK_GENE_MAP_v1_1.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "export_peak_gene_map.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3C2_v3_SOURCE_RELATION_COUNT_LOCK.csv"), file.path(dirs$results, "STEP03B3C2_v3_GATE_TABLE.csv"), file.path(dirs$results, 
    "STEP03B3C2_v3_PEAK_GENE_MAP_LOCK_CARD.md"), file.path(dirs$qc, "STEP03B3C2_v3_CANDIDATE_ARTIFACT_IDENTITY_RECHECK.csv"), 
    file.path(locked_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"), locked_long_map_path, locked_gene_preview_path, 
    file.path(dirs$control, "STEP03B3C3_AUTHORIZATION.md"), file.path(dirs$text, "STEP03B3C2_v3_Methods_addendum.md"), 
    file.path(dirs$text, "STEP03B3C2_v3_Results_draft.md"), file.path(dirs$parameters, "STEP03B3C2_v3_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2_v3.txt"), file.path(dirs$code, "export_peak_gene_map.R"), 
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
    else if (grepl("Methods|Results|AUTHORIZATION", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|RECHECK|MANIFEST|LOCK_CARD|COUNT_LOCK", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C2_v3_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C2 v3 archive status", "", paste0("- Run ID: ", run_id), paste0("- Peak-gene map lock: ", 
    peak_gene_map_lock_version), "- Peak-gene relations: 23,499", "- Unique author gene symbols: 10,042", "- Canonical source IDs: 3", 
    "- Relation IDs unique: YES", "- Relation membership changed: NO", "- RNA/CUT&RUN integration: NOT STARTED", 
    "- Scientific figures: NONE", "- STEP03B3C3 authorization: YES", paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Completed substep: STEP03B3C2", paste0("- Accepted run: ", run_id), 
    paste0("- Published cistrome source: ", source_lock_version, " IMMUTABLE"), paste0("- Peak-gene map: ", peak_gene_map_lock_version, 
        " LOCKED"), paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE; 436 candidates"), paste0("- ATAC role: ", 
        ATAC_role), "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", "- Active next substep: STEP03B3C3", "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, 
    "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C2 v3 PEAK_GENE_MAP_v1.1 lock completed.\n")

cat("Peak-gene map: ", peak_gene_map_lock_version, " LOCKED\n")

cat("Peak-gene relations: 23,499\n")

cat("Unique author gene symbols: 10,042\n")

cat("Vehicle relations: 2,848\n")

cat("Calcitriol relations: 19,037\n")

cat("Shared relations: 1,614\n")

cat("Canonical source IDs: 3/3\n")

cat("Relation IDs unique: 23,499/23,499\n")

cat("Promoter-only exceptions: 0\n")

cat("Invalid gene-token rows: 0\n")

cat("RNA anchor: ", RNA_anchor, " IMMUTABLE; 436 candidates\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("STEP03B3C3 authorization: YES\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n")

cat("============================================================\n")

