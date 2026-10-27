#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A18: Calculate the 436-gene occupancy integration

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3B_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

RNA_source_lock_version <- "RNA_INTERACTION_SOURCE_v1.0"

gene_ID_policy_lock_version <- "RNA_CISTROME_GENE_ID_POLICY_v1.0"

gene_ID_map_lock_version <- "RNA_CISTROME_GENE_ID_MAP_v1.0"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

integration_candidate_version <- "RNA_OCCUPANCY_EVIDENCE_CANDIDATE_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_RNA_n <- 436L

expected_exact_join_n <- 279L

expected_manual_join_n <- 5L

expected_join_authorized_n <- 284L

expected_RNA_only_n <- 152L

expected_peak_gene_relation_n <- 23499L

expected_cistrome_gene_n <- 10042L

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

integration_root <- file.path(dirs$results, "STEP03B3C3B_RNA_occupancy_integration_candidate")

dir.create(integration_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3B_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

parse_bool <- function(x) {
    value <- toupper(trimws(as.character(x)))
    value %in% c("TRUE", "T", "1")
}

safe_bool_column <- function(data_frame, column_name) {
    if (!column_name %in% names(data_frame)) {
        return(rep(FALSE, nrow(data_frame)))
    }
    parse_bool(data_frame[[column_name]])
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

log_event("Starting", run_id)

RNA_gene_ID_root <- file.path(dirs$results, "STEP03B3C3A_locked_RNA_gene_ID", gene_ID_policy_lock_version)

peak_gene_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

required_paths <- c(gate_c3a = file.path(dirs$results, "STEP03B3C3A_v4_GATE_TABLE.csv"), authorization = file.path(dirs$control, 
    "STEP03B3C3B_AUTHORIZATION.md"), RNA_gene_ID_manifest = file.path(RNA_gene_ID_root, "RNA_CISTROME_GENE_ID_POLICY_v1.0_LOCKED_MANIFEST.csv"), 
    RNA_source_lock = file.path(RNA_gene_ID_root, "RNA_INTERACTION_SOURCE_v1.0_LOCK.csv"), RNA_gene_ID_map = file.path(RNA_gene_ID_root, 
        "RNA_CISTROME_GENE_ID_MAP_v1.0.csv"), peak_gene_manifest = file.path(peak_gene_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3B input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_c3a <- .vdr_read_csv(required_paths[["gate_c3a"]], stringsAsFactors = FALSE, check.names = FALSE)

authorization_text <- readLines(required_paths[["authorization"]], warn = FALSE, encoding = "UTF-8")

RNA_gene_ID_manifest <- .vdr_read_csv(required_paths[["RNA_gene_ID_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_source_lock <- .vdr_read_csv(required_paths[["RNA_source_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_gene_ID_map <- .vdr_read_csv(required_paths[["RNA_gene_ID_map"]], stringsAsFactors = FALSE, check.names = FALSE)

peak_gene_manifest <- .vdr_read_csv(required_paths[["peak_gene_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- gate_c3a$status[gate_c3a$gate_component == "STEP03B3C3A v4 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "READY_FOR_STEP03B3C3B_RNA_OCCUPANCY_EVIDENCE_INTEGRATION") {
    stop("STEP03B3C3A v4 does not authorize integration. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (!any(grepl("Authorization: YES", authorization_text, fixed = TRUE))) {
    stop("STEP03B3C3B authorization is incomplete.")
}

if (nrow(RNA_gene_ID_manifest) != 4L || !all(RNA_gene_ID_manifest$lock_status == "LOCKED")) {
    stop("RNA/cistrome gene-ID locked manifest is incomplete.")
}

RNA_gene_ID_manifest$file_exists <- file.exists(RNA_gene_ID_manifest$path)

RNA_gene_ID_manifest$observed_md5 <- vapply(RNA_gene_ID_manifest$path, md5_file, character(1))

RNA_gene_ID_manifest$observed_sha256 <- vapply(RNA_gene_ID_manifest$path, sha256_file, character(1))

RNA_gene_ID_manifest$identity_pass <- with(RNA_gene_ID_manifest, file_exists & md5 == observed_md5 & sha256 == 
    observed_sha256)

if (!all(RNA_gene_ID_manifest$identity_pass)) {
    stop("A locked RNA/gene-ID artifact failed hash validation.")
}

if (nrow(peak_gene_manifest) != 2L || !all(peak_gene_manifest$lock_version == peak_gene_map_lock_version) || !all(peak_gene_manifest$lock_status == 
    "LOCKED")) {
    stop("Peak-gene locked manifest is incomplete.")
}

peak_gene_manifest$file_exists <- file.exists(peak_gene_manifest$path)

peak_gene_manifest$observed_md5 <- vapply(peak_gene_manifest$path, md5_file, character(1))

peak_gene_manifest$observed_sha256 <- vapply(peak_gene_manifest$path, sha256_file, character(1))

peak_gene_manifest$identity_pass <- with(peak_gene_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(peak_gene_manifest$identity_pass)) {
    stop("A locked peak-gene artifact failed hash validation.")
}

if (nrow(RNA_source_lock) != 1L || RNA_source_lock$lock_version != RNA_source_lock_version || RNA_source_lock$gene_symbol_column != 
    "feature_id" || RNA_source_lock$n_rows != expected_RNA_n || RNA_source_lock$unique_gene_n != expected_RNA_n || 
    RNA_source_lock$lock_status != "LOCKED") {
    stop("RNA source lock metadata failed validation.")
}

RNA_source_path <- RNA_source_lock$source_file

if (!file.exists(RNA_source_path) || md5_file(RNA_source_path) != RNA_source_lock$md5 || sha256_file(RNA_source_path) != 
    RNA_source_lock$sha256) {
    stop("Frozen RNA source failed file or hash validation.")
}

RNA_table <- .vdr_read_csv(RNA_source_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(RNA_table) != expected_RNA_n || !"feature_id" %in% names(RNA_table) || length(unique(RNA_table$feature_id)) != 
    expected_RNA_n) {
    stop("Frozen RNA table failed row or feature_id validation.")
}

required_gene_ID_columns <- c("RNA_row", "RNA_gene_symbol_original", "RNA_gene_symbol_cistrome_join", "mapping_status", 
    "mapping_rule", "exact_original_match", "manual_mapping_accepted", "cistrome_join_authorized", "RNA_membership_status")

missing_gene_ID_columns <- setdiff(required_gene_ID_columns, names(RNA_gene_ID_map))

if (length(missing_gene_ID_columns) > 0L) {
    stop("Gene-ID map is missing columns: ", paste(missing_gene_ID_columns, collapse = ", "))
}

if (nrow(RNA_gene_ID_map) != expected_RNA_n || length(unique(RNA_gene_ID_map$RNA_row)) != expected_RNA_n || length(unique(RNA_gene_ID_map$RNA_gene_symbol_original)) != 
    expected_RNA_n || !identical(as.character(RNA_gene_ID_map$RNA_gene_symbol_original), as.character(RNA_table$feature_id))) {
    stop("Gene-ID map does not exactly track the frozen RNA source.")
}

RNA_gene_ID_map$exact_original_match <- parse_bool(RNA_gene_ID_map$exact_original_match)

RNA_gene_ID_map$manual_mapping_accepted <- parse_bool(RNA_gene_ID_map$manual_mapping_accepted)

RNA_gene_ID_map$cistrome_join_authorized <- parse_bool(RNA_gene_ID_map$cistrome_join_authorized)

if (sum(RNA_gene_ID_map$exact_original_match) != expected_exact_join_n || sum(RNA_gene_ID_map$manual_mapping_accepted) != 
    expected_manual_join_n || sum(RNA_gene_ID_map$cistrome_join_authorized) != expected_join_authorized_n || sum(!RNA_gene_ID_map$cistrome_join_authorized) != 
    expected_RNA_only_n) {
    stop("Gene-ID lock counts failed reconciliation.")
}

peak_gene_long_path <- peak_gene_manifest$path[peak_gene_manifest$artifact_id == "PEAK_GENE_LONG"]

if (length(peak_gene_long_path) != 1L) {
    stop("Locked peak-gene long map could not be uniquely identified.")
}

peak_gene_map <- .vdr_fread(peak_gene_long_path, sep = "\t", data.table = FALSE, showProgress = FALSE)

required_peak_gene_columns <- c("source_id", "occupancy_context", "coordinate_key", "gene_symbol_author", "peak_gene_relation_id")

missing_peak_gene_columns <- setdiff(required_peak_gene_columns, names(peak_gene_map))

if (length(missing_peak_gene_columns) > 0L) {
    stop("Peak-gene map is missing columns: ", paste(missing_peak_gene_columns, collapse = ", "))
}

if (nrow(peak_gene_map) != expected_peak_gene_relation_n || anyDuplicated(peak_gene_map$peak_gene_relation_id) > 
    0L) {
    stop("Peak-gene map row count or relation-ID uniqueness failed.")
}

cistrome_symbols <- unique(peak_gene_map$gene_symbol_author)

if (length(cistrome_symbols) != expected_cistrome_gene_n) {
    stop("Cistrome gene-symbol universe differs from 10,042.")
}

authorized_ID_map <- RNA_gene_ID_map[RNA_gene_ID_map$cistrome_join_authorized, , drop = FALSE]

relation_join <- merge(authorized_ID_map, peak_gene_map, by.x = "RNA_gene_symbol_cistrome_join", by.y = "gene_symbol_author", 
    all.x = TRUE, sort = FALSE)

if (nrow(relation_join) == 0L || any(is.na(relation_join$peak_gene_relation_id))) {
    stop("At least one cistrome-authorized RNA symbol lacks a peak-gene relation.")
}

joined_authorized_symbols <- unique(relation_join$RNA_gene_symbol_original)

if (length(joined_authorized_symbols) != expected_join_authorized_n || !setequal(joined_authorized_symbols, authorized_ID_map$RNA_gene_symbol_original)) {
    stop("The relation-level join does not cover all 284 authorized RNA symbols.")
}

relation_join$author_gene_overlap_supported <- safe_bool_column(relation_join, "author_gene_overlap") | safe_bool_column(relation_join, 
    "vehicle_gene_overlap") | safe_bool_column(relation_join, "calcitriol_gene_overlap")

relation_join$promoter_overlap_supported <- safe_bool_column(relation_join, "author_promoter_overlap") | safe_bool_column(relation_join, 
    "vehicle_promoter_overlap") | safe_bool_column(relation_join, "calcitriol_promoter_overlap")

shared_coordinate_keys <- unique(peak_gene_map$coordinate_key[peak_gene_map$source_id == "VDR_SHARED_CISTROME"])

relation_join$is_shared_coordinate <- relation_join$coordinate_key %in% shared_coordinate_keys

relation_join$is_vehicle_nonshared_relation <- relation_join$source_id == "VDR_VEHICLE_CISTROME" & !relation_join$is_shared_coordinate

relation_join$is_calcitriol_nonshared_relation <- relation_join$source_id == "VDR_CALCITRIOL_CISTROME" & !relation_join$is_shared_coordinate

relation_join$integration_candidate_version <- integration_candidate_version

relation_join$scientific_role <- "RNA_INTERACTION_CANDIDATE_WITH_CROSS_CONTEXT_VDR_OCCUPANCY_SUPPORT"

relation_join$definitive_direct_target <- FALSE

relation_dt <- data.table::as.data.table(relation_join)

gene_occupancy_summary <- relation_dt[, .(total_unique_peak_n = data.table::uniqueN(coordinate_key), total_peak_gene_relation_n = data.table::uniqueN(peak_gene_relation_id), 
    vehicle_peak_n = data.table::uniqueN(coordinate_key[source_id == "VDR_VEHICLE_CISTROME"]), calcitriol_peak_n = data.table::uniqueN(coordinate_key[source_id == 
        "VDR_CALCITRIOL_CISTROME"]), shared_peak_n = data.table::uniqueN(coordinate_key[source_id == "VDR_SHARED_CISTROME"]), 
    vehicle_nonshared_peak_n = data.table::uniqueN(coordinate_key[is_vehicle_nonshared_relation]), calcitriol_nonshared_peak_n = data.table::uniqueN(coordinate_key[is_calcitriol_nonshared_relation]), 
    any_gene_overlap_support = any(author_gene_overlap_supported, na.rm = TRUE), any_promoter_overlap_support = any(promoter_overlap_supported, 
        na.rm = TRUE), vehicle_promoter_peak_n = data.table::uniqueN(coordinate_key[source_id == "VDR_VEHICLE_CISTROME" & 
        promoter_overlap_supported]), calcitriol_promoter_peak_n = data.table::uniqueN(coordinate_key[source_id == 
        "VDR_CALCITRIOL_CISTROME" & promoter_overlap_supported]), shared_promoter_peak_n = data.table::uniqueN(coordinate_key[source_id == 
        "VDR_SHARED_CISTROME" & promoter_overlap_supported])), by = .(RNA_row, RNA_gene_symbol_original)]

gene_occupancy_summary <- as.data.frame(gene_occupancy_summary, stringsAsFactors = FALSE)

integer_summary_columns <- c("total_unique_peak_n", "total_peak_gene_relation_n", "vehicle_peak_n", "calcitriol_peak_n", 
    "shared_peak_n", "vehicle_nonshared_peak_n", "calcitriol_nonshared_peak_n", "vehicle_promoter_peak_n", "calcitriol_promoter_peak_n", 
    "shared_promoter_peak_n")

for (column_name in integer_summary_columns) {
    gene_occupancy_summary[[column_name]][is.na(gene_occupancy_summary[[column_name]])] <- 0L
}

RNA_with_ID <- cbind(RNA_table, RNA_gene_ID_map[setdiff(names(RNA_gene_ID_map), c("RNA_row", "RNA_gene_symbol_original"))])

if (nrow(RNA_with_ID) != expected_RNA_n || !identical(as.character(RNA_with_ID$feature_id), as.character(RNA_gene_ID_map$RNA_gene_symbol_original))) {
    stop("RNA source and gene-ID map alignment failed.")
}

integration_candidate <- merge(RNA_with_ID, gene_occupancy_summary, by.x = "feature_id", by.y = "RNA_gene_symbol_original", 
    all.x = TRUE, sort = FALSE)

integration_candidate <- integration_candidate[match(RNA_table$feature_id, integration_candidate$feature_id), , 
    drop = FALSE]

unexpected_merge_suffix_columns <- grep("\\.(x|y)$", names(integration_candidate), value = TRUE)

if (length(unexpected_merge_suffix_columns) > 0L) {
    stop("Unexpected duplicate columns after RNA/occupancy merge: ", paste(unexpected_merge_suffix_columns, collapse = ", "))
}

if (nrow(integration_candidate) != expected_RNA_n || anyDuplicated(integration_candidate$feature_id) > 0L || !identical(as.character(integration_candidate$feature_id), 
    as.character(RNA_table$feature_id))) {
    stop("The 436-row integration candidate failed row identity or ordering.")
}

for (column_name in integer_summary_columns) {
    integration_candidate[[column_name]][is.na(integration_candidate[[column_name]])] <- 0L
}

logical_summary_columns <- c("any_gene_overlap_support", "any_promoter_overlap_support")

for (column_name in logical_summary_columns) {
    integration_candidate[[column_name]][is.na(integration_candidate[[column_name]])] <- FALSE
}

integration_candidate$occupancy_supported <- integration_candidate$total_peak_gene_relation_n > 0L

integration_candidate$shared_occupancy_supported <- integration_candidate$shared_peak_n > 0L

integration_candidate$vehicle_occupancy_supported <- integration_candidate$vehicle_peak_n > 0L

integration_candidate$calcitriol_occupancy_supported <- integration_candidate$calcitriol_peak_n > 0L

integration_candidate$occupancy_evidence_group <- ifelse(integration_candidate$shared_peak_n > 0L, "RNA_INTERACTION_PLUS_SHARED_VDR_OCCUPANCY", 
    ifelse(integration_candidate$vehicle_peak_n > 0L & integration_candidate$calcitriol_peak_n > 0L, "RNA_INTERACTION_PLUS_BOTH_CONDITION_OCCUPANCY_NO_SHARED_PEAK", 
        ifelse(integration_candidate$calcitriol_peak_n > 0L, "RNA_INTERACTION_PLUS_CALCITRIOL_OCCUPANCY", ifelse(integration_candidate$vehicle_peak_n > 
            0L, "RNA_INTERACTION_PLUS_VEHICLE_OCCUPANCY", "RNA_INTERACTION_ONLY_NO_CURRENT_CISTROME_SUPPORT"))))

integration_candidate$occupancy_positional_support <- ifelse(integration_candidate$any_promoter_overlap_support, 
    "PROMOTER_OCCUPANCY_PRESENT", ifelse(integration_candidate$occupancy_supported, "NONPROMOTER_OR_UNSPECIFIED_GENE_OVERLAP_ONLY", 
        "NO_CURRENT_OCCUPANCY_SUPPORT"))

integration_candidate$integration_candidate_version <- integration_candidate_version

integration_candidate$integrated_scientific_label <- ifelse(integration_candidate$occupancy_supported, "VDR_OCCUPANCY_SUPPORTED_RNA_INTERACTION_CANDIDATE", 
    "RNA_INTERACTION_CANDIDATE_WITHOUT_CURRENT_CISTROME_SUPPORT")

integration_candidate$direct_target_assignment <- "NOT_PERFORMED"

integration_candidate$causal_claim_authorized <- FALSE

integration_candidate$ATAC_evidence_status <- "NOT_ADDED_RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

if (sum(integration_candidate$occupancy_supported) != expected_join_authorized_n || sum(!integration_candidate$occupancy_supported) != 
    expected_RNA_only_n) {
    stop("Integrated occupancy-supported/RNA-only counts differ from 284/152.")
}

if (any(integration_candidate$cistrome_join_authorized & !integration_candidate$occupancy_supported) || any(!integration_candidate$cistrome_join_authorized & 
    integration_candidate$occupancy_supported)) {
    stop("Gene-ID join authorization and observed occupancy support are inconsistent.")
}

manual_mapping_join_audit <- integration_candidate[integration_candidate$manual_mapping_accepted, c("feature_id", 
    "RNA_gene_symbol_cistrome_join", "mapping_status", "total_unique_peak_n", "vehicle_peak_n", "calcitriol_peak_n", 
    "shared_peak_n", "any_promoter_overlap_support", "occupancy_evidence_group"), drop = FALSE]

if (nrow(manual_mapping_join_audit) != expected_manual_join_n || any(manual_mapping_join_audit$total_unique_peak_n < 
    1L)) {
    stop("At least one manually accepted symbol failed occupancy integration.")
}

RNA_only_audit <- integration_candidate[!integration_candidate$occupancy_supported, c("feature_id", "RNA_gene_symbol_cistrome_join", 
    "mapping_status", "program_class", "interaction_logFC", "interaction_FDR", "occupancy_evidence_group"), drop = FALSE]

if (nrow(RNA_only_audit) != expected_RNA_only_n || any(!is.na(RNA_only_audit$RNA_gene_symbol_cistrome_join) & nzchar(trimws(RNA_only_audit$RNA_gene_symbol_cistrome_join)))) {
    stop("The RNA-only 152-row audit failed.")
}

evidence_group_counts <- as.data.frame(table(integration_candidate$occupancy_evidence_group), stringsAsFactors = FALSE)

names(evidence_group_counts) <- c("occupancy_evidence_group", "RNA_candidate_n")

evidence_group_counts$RNA_candidate_fraction <- evidence_group_counts$RNA_candidate_n/expected_RNA_n

program_class_by_evidence <- as.data.frame(table(integration_candidate$program_class, integration_candidate$occupancy_evidence_group), 
    stringsAsFactors = FALSE)

names(program_class_by_evidence) <- c("program_class", "occupancy_evidence_group", "RNA_candidate_n")

occupancy_source_summary <- data.frame(occupancy_metric = c("RNA candidates with any occupancy support", "RNA candidates with vehicle occupancy", 
    "RNA candidates with calcitriol occupancy", "RNA candidates with shared occupancy", "RNA candidates with promoter occupancy", 
    "RNA-only candidates without current cistrome support", "Manually harmonized candidates with occupancy support"), 
    RNA_candidate_n = c(sum(integration_candidate$occupancy_supported), sum(integration_candidate$vehicle_occupancy_supported), 
        sum(integration_candidate$calcitriol_occupancy_supported), sum(integration_candidate$shared_occupancy_supported), 
        sum(integration_candidate$any_promoter_overlap_support), sum(!integration_candidate$occupancy_supported), 
        nrow(manual_mapping_join_audit)), stringsAsFactors = FALSE)

relation_source_summary <- as.data.frame(table(relation_join$source_id), stringsAsFactors = FALSE)

names(relation_source_summary) <- c("source_id", "joined_peak_gene_relation_n")

relation_source_gene_counts <- aggregate(RNA_gene_symbol_original ~ source_id, data = relation_join, FUN = function(x) {
    length(unique(x))
})

names(relation_source_gene_counts)[names(relation_source_gene_counts) == "RNA_gene_symbol_original"] <- "joined_RNA_candidate_n"

relation_source_summary <- merge(relation_source_summary, relation_source_gene_counts, by = "source_id", all = TRUE, 
    sort = FALSE)

mapping_status_summary <- as.data.frame(table(integration_candidate$mapping_status, integration_candidate$occupancy_supported), 
    stringsAsFactors = FALSE)

names(mapping_status_summary) <- c("mapping_status", "occupancy_supported", "RNA_candidate_n")

integration_QC_summary <- data.frame(metric = c("Frozen RNA rows", "Unique frozen RNA feature_id", "Exact gene-ID joins", 
    "Manual punctuation joins", "Cistrome-authorized RNA candidates", "RNA-only candidates", "Integrated candidate rows", 
    "Unique integrated feature_id", "Occupancy-supported RNA candidates", "RNA-only integrated rows", "Relation-level joined rows", 
    "Joined unique RNA candidates", "Manual mappings successfully joined", "Locked peak-gene relation rows", "Peak-gene membership changed", 
    "RNA membership changed", "Definitive direct-target assignments", "Scientific figures"), value = c(nrow(RNA_table), 
    length(unique(RNA_table$feature_id)), sum(RNA_gene_ID_map$exact_original_match), sum(RNA_gene_ID_map$manual_mapping_accepted), 
    sum(RNA_gene_ID_map$cistrome_join_authorized), sum(!RNA_gene_ID_map$cistrome_join_authorized), nrow(integration_candidate), 
    length(unique(integration_candidate$feature_id)), sum(integration_candidate$occupancy_supported), sum(!integration_candidate$occupancy_supported), 
    nrow(relation_join), length(unique(relation_join$RNA_gene_symbol_original)), nrow(manual_mapping_join_audit), 
    nrow(peak_gene_map), 0L, 0L, 0L, 0L), stringsAsFactors = FALSE)

integration_candidate_path <- file.path(integration_root, "STEP03B3C3B_RNA_OCCUPANCY_INTEGRATED_CANDIDATE.csv")

relation_join_path <- file.path(integration_root, "STEP03B3C3B_RNA_OCCUPANCY_RELATION_LONG.tsv.gz")

write_csv_utf8(integration_candidate, integration_candidate_path)

write_tsv_gz(relation_join, relation_join_path)

write_csv_utf8(evidence_group_counts, file.path(integration_root, "STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv"))

write_csv_utf8(program_class_by_evidence, file.path(integration_root, "STEP03B3C3B_PROGRAM_CLASS_BY_EVIDENCE_GROUP.csv"))

write_csv_utf8(occupancy_source_summary, file.path(integration_root, "STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv"))

write_csv_utf8(relation_source_summary, file.path(integration_root, "STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv"))

write_csv_utf8(mapping_status_summary, file.path(integration_root, "STEP03B3C3B_MAPPING_STATUS_BY_OCCUPANCY.csv"))

write_csv_utf8(integration_QC_summary, file.path(dirs$qc, "STEP03B3C3B_INTEGRATION_QC_SUMMARY.csv"))

write_csv_utf8(manual_mapping_join_audit, file.path(dirs$qc, "STEP03B3C3B_MANUAL_MAPPING_JOIN_AUDIT.csv"))

write_csv_utf8(RNA_only_audit, file.path(dirs$qc, "STEP03B3C3B_RNA_ONLY_152_AUDIT.csv"))

candidate_artifact_paths <- c(integration_candidate_path, relation_join_path, file.path(integration_root, "STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv"), 
    file.path(integration_root, "STEP03B3C3B_PROGRAM_CLASS_BY_EVIDENCE_GROUP.csv"), file.path(integration_root, 
        "STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv"), file.path(integration_root, "STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv"))

candidate_manifest <- data.frame(artifact_id = c("RNA_OCCUPANCY_INTEGRATED_CANDIDATE", "RNA_OCCUPANCY_RELATION_LONG", 
    "EVIDENCE_GROUP_COUNTS", "PROGRAM_CLASS_BY_EVIDENCE_GROUP", "OCCUPANCY_SOURCE_SUMMARY", "RELATION_SOURCE_SUMMARY"), 
    candidate_version = integration_candidate_version, path = normalize_slashes(candidate_artifact_paths), size_bytes = as.numeric(file.info(candidate_artifact_paths)$size), 
    md5 = vapply(candidate_artifact_paths, md5_file, character(1)), sha256 = vapply(candidate_artifact_paths, sha256_file, 
        character(1)), RNA_membership_status = "IMMUTABLE_436", peak_gene_membership_status = "IMMUTABLE_23499", 
    lock_status = "CANDIDATE_PENDING_MANUAL_QC", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest, file.path(integration_root, "STEP03B3C3B_INTEGRATION_CANDIDATE_MANIFEST.csv"))

gate_status <- "PROCEED_TO_MANUAL_RNA_OCCUPANCY_INTEGRATION_LOCK"

gate_table <- data.frame(gate_component = c("STEP03B3C3A authorization", "RNA source lock", "Gene-ID policy lock", 
    "Peak-gene map lock", "Frozen RNA rows", "Cistrome-authorized RNA candidates", "RNA-only candidates retained", 
    "Relation-level join completeness", "Manual mapping integration", "Integrated candidate uniqueness", "Occupancy evidence groups", 
    "Promoter evidence", "RNA membership", "Peak-gene membership", "ATAC role", "Definitive direct-target assignment", 
    "Scientific figures", "STEP03B3C3B gate"), status = c("PASS", paste0("LOCKED_", RNA_source_lock_version), paste0("LOCKED_", 
    gene_ID_policy_lock_version), paste0("LOCKED_", peak_gene_map_lock_version), "PASS_436", "PASS_284", "PASS_152", 
    "PASS_284_OF_284", "PASS_5_OF_5", "PASS_436_UNIQUE", "CREATED_NONCAUSAL", "RETAINED_SEPARATELY", "UNCHANGED_436", 
    "UNCHANGED_23499", ATAC_role, "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain authorization.", 
    "Use only the locked RNA source.", "Use only the locked 436-row gene-ID map.", "Use only the locked 23,499-row relation map.", 
    "Do not add or remove RNA rows.", "Review all occupancy-supported rows.", "Retain as RNA-only candidates.", 
    "Confirm all authorized joins have relations.", "Review five manually harmonized joins.", "Use feature_id as integrated row authority.", 
    "Use evidence groups, not causal tiers.", "Do not collapse promoter and generic overlap evidence.", "Membership changes are prohibited.", 
    "Membership changes are prohibited.", "Do not add ATAC until selected-locus support is prespecified.", "Use occupancy-supported candidate terminology only.", 
    "No Scientific Figure Skill 3 workflow is triggered.", "Upload candidate tables and QC summaries."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3B_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3B RNA/occupancy integration candidate card", "", paste0("- Run ID: ", run_id), paste0("- Candidate version: ", 
    integration_candidate_version), "- Frozen RNA candidates: 436", "- Cistrome-authorized RNA candidates: 284", 
    "- RNA-only candidates retained: 152", paste0("- Relation-level joined rows: ", nrow(relation_join)), paste0("- Occupancy-supported RNA candidates: ", 
        sum(integration_candidate$occupancy_supported)), paste0("- Shared-occupancy-supported RNA candidates: ", 
        sum(integration_candidate$shared_occupancy_supported)), paste0("- Promoter-occupancy-supported RNA candidates: ", 
        sum(integration_candidate$any_promoter_overlap_support)), "- Manual punctuation mappings integrated: 5/5", 
    "- RNA membership changed: NO", "- Peak-gene membership changed: NO", "- ATAC evidence: NOT ADDED", "- Definitive direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C3B_RNA_OCCUPANCY_INTEGRATION_CANDIDATE_CARD.md"))

write_text(c("# STEP03B3C3B Methods draft", "", "The frozen 436-gene RNA interaction program was integrated with the locked VDR peak-gene map through the locked RNA/cistrome gene-ID policy. The original RNA `feature_id` remained the row authority, while 279 exact and five manually reviewed punctuation-equivalent cistrome join symbols were used for relation-level matching.", 
    "", "For each RNA candidate, unique supporting peaks and peak-gene relations were summarized separately for vehicle, calcitriol and shared VDR occupancy. Condition-level totals were distinguished from nonshared condition-specific peaks, and promoter overlap was retained separately from generic author-reported gene overlap.", 
    "", "All 152 RNA candidates without a current cistrome match were retained as RNA-only rows. Occupancy evidence was interpreted as cross-context convergent support rather than definitive direct regulation. No ATAC evidence or scientific figure was added in this substep."), 
    file.path(dirs$text, "STEP03B3C3B_Methods_draft.md"))

write_text(c("# STEP03B3C3B interpretation constraints", "", "- The RNA program is derived from siVDR +/- calcitriol interaction analysis.", 
    "- The VDR cistrome is derived from VDR-overexpression occupancy contexts.", "- Their integration is convergent cross-context evidence.", 
    "- Shared occupancy support does not establish causality.", "- Promoter overlap is stronger positional support than generic overlap but is not proof of direct regulation.", 
    "- A missing cistrome match does not invalidate an RNA interaction candidate.", "- The 436-member RNA set and 23,499 peak-gene relations remain immutable.", 
    "- ATAC may later be used only as selected-locus accessible context.", "- Terms such as definitive direct target, causal target or diagnostic marker remain prohibited."), 
    file.path(dirs$text, "STEP03B3C3B_Interpretation_constraints.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, integration_candidate_version = integration_candidate_version, 
    RNA_anchor = RNA_anchor, RNA_source_lock_version = RNA_source_lock_version, gene_ID_policy_lock_version = gene_ID_policy_lock_version, 
    peak_gene_map_lock_version = peak_gene_map_lock_version, RNA_n = expected_RNA_n, cistrome_authorized_RNA_n = expected_join_authorized_n, 
    RNA_only_n = expected_RNA_only_n, relation_level_join_n = nrow(relation_join), occupancy_supported_RNA_n = sum(integration_candidate$occupancy_supported), 
    shared_occupancy_supported_RNA_n = sum(integration_candidate$shared_occupancy_supported), promoter_occupancy_supported_RNA_n = sum(integration_candidate$any_promoter_overlap_support), 
    RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, ATAC_evidence_added = FALSE, definitive_direct_target_assignment = FALSE, 
    scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3B_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3B.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("45_step03B3C3B_integrate_RNA_occupancy_evidence.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "integrate_rna_occupancy.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3B_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3B_RNA_OCCUPANCY_INTEGRATION_CANDIDATE_CARD.md"), 
    integration_candidate_path, relation_join_path, file.path(integration_root, "STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv"), 
    file.path(integration_root, "STEP03B3C3B_PROGRAM_CLASS_BY_EVIDENCE_GROUP.csv"), file.path(integration_root, 
        "STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv"), file.path(integration_root, "STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv"), 
    file.path(integration_root, "STEP03B3C3B_MAPPING_STATUS_BY_OCCUPANCY.csv"), file.path(integration_root, "STEP03B3C3B_INTEGRATION_CANDIDATE_MANIFEST.csv"), 
    file.path(dirs$qc, "STEP03B3C3B_INTEGRATION_QC_SUMMARY.csv"), file.path(dirs$qc, "STEP03B3C3B_MANUAL_MAPPING_JOIN_AUDIT.csv"), 
    file.path(dirs$qc, "STEP03B3C3B_RNA_ONLY_152_AUDIT.csv"), file.path(dirs$text, "STEP03B3C3B_Methods_draft.md"), 
    file.path(dirs$text, "STEP03B3C3B_Interpretation_constraints.md"), file.path(dirs$parameters, "STEP03B3C3B_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3B.txt"), file.path(dirs$code, "integrate_rna_occupancy.R"), 
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
    else if (grepl("Methods|Interpretation", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|AUDIT|SUMMARY|CARD|MANIFEST|COUNTS|MAPPING_STATUS", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3B_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3B archive status", "", paste0("- Run ID: ", run_id), paste0("- Integration candidate: ", 
    integration_candidate_version), "- Frozen RNA candidates: 436", "- Occupancy-supported RNA candidates: 284", 
    "- RNA-only candidates: 152", paste0("- Relation-level joined rows: ", nrow(relation_join)), "- RNA membership changed: NO", 
    "- Peak-gene membership changed: NO", "- ATAC evidence added: NO", "- Direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C3B", paste0("- Run ID: ", run_id), paste0("- RNA source: ", 
    RNA_source_lock_version, " IMMUTABLE"), paste0("- Gene-ID policy: ", gene_ID_policy_lock_version, " IMMUTABLE"), 
    paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE"), paste0("- Integration candidate: ", 
        integration_candidate_version, " PENDING_MANUAL_QC"), "- RNA membership: 436 UNCHANGED", "- Occupancy-supported RNA candidates: 284", 
    "- RNA-only candidates: 152", "- ATAC evidence: NOT ADDED", "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3B RNA/occupancy evidence integration completed.\n")

cat("Integration candidate: ", integration_candidate_version, "\n")

cat("Frozen RNA candidates: 436\n")

cat("Cistrome-authorized RNA candidates: 284\n")

cat("RNA-only candidates retained: 152\n")

cat("Relation-level joined rows: ", nrow(relation_join), "\n")

cat("Shared-occupancy-supported RNA candidates: ", sum(integration_candidate$shared_occupancy_supported), "\n")

cat("Promoter-occupancy-supported RNA candidates: ", sum(integration_candidate$any_promoter_overlap_support), "\n")

cat("Manual punctuation mappings integrated: 5/5\n")

cat("RNA membership changed: NO\n")

cat("Peak-gene membership changed: NO\n")

cat("ATAC evidence added: NO\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C3B_GATE_TABLE.csv\n")

cat("  08_qc/STEP03B3C3B_INTEGRATION_QC_SUMMARY.csv\n")

cat("  08_qc/STEP03B3C3B_MANUAL_MAPPING_JOIN_AUDIT.csv\n")

cat("  06_results/STEP03B3C3B_RNA_occupancy_integration_candidate/STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv\n")

cat("  06_results/STEP03B3C3B_RNA_occupancy_integration_candidate/STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv\n")

cat("  06_results/STEP03B3C3B_RNA_occupancy_integration_candidate/STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv\n")

cat("  06_results/STEP03B3C3B_RNA_occupancy_integration_candidate/STEP03B3C3B_INTEGRATION_CANDIDATE_MANIFEST.csv\n")

cat("============================================================\n")

