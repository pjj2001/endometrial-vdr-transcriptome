#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A12: Build peak-gene relations from published annotations

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"

RNA_anchor <- "RNA_PROGRAM_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

peak_gene_map_candidate_version <- "PEAK_GENE_MAP_CANDIDATE_v1.0"

minimum_shared_annotation_concordance <- 1

maximum_unparseable_gene_token_fraction <- 0.001

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

map_root <- file.path(dirs$results, "STEP03B3C2_peak_gene_map_candidate")

dir.create(map_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

parse_logical_text <- function(x) {
    value <- toupper(trimws(as.character(x)))
    value == "TRUE"
}

clean_author_field <- function(x) {
    value <- trimws(as.character(x))
    value[is.na(value) | !nzchar(value) | value %in% c(".", "NA", "N/A", "NULL", "None", "none")] <- NA_character_
    value
}

split_gene_tokens <- function(x) {
    value <- clean_author_field(x)
    if (length(value) == 0L || is.na(value)) {
        return(character(0))
    }
    tokens <- unlist(strsplit(value, "\\s+", perl = TRUE), use.names = FALSE)
    tokens <- trimws(tokens)
    tokens <- tokens[nzchar(tokens) & tokens != "."]
    unique(tokens)
}

valid_gene_symbol_shape <- function(x) {
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", x, perl = TRUE)
}

canonical_gene_set_string <- function(x) {
    tokens <- split_gene_tokens(x)
    if (length(tokens) == 0L) {
        return("")
    }
    paste(sort(tokens), collapse = ";")
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))
}

required_annotation_columns <- c("source_id", "source_row", "chromosome", "start0", "end", "width", "valid_coordinate", 
    "primary_chromosome", "coordinate_key", "raw_distance", "raw_tssDistance", "raw_relativeLocation", "raw_geneOverlap", 
    "raw_promoterOverlap")

build_condition_peak_gene_map <- function(source_table, occupancy_context) {
    missing_columns <- setdiff(required_annotation_columns, names(source_table))
    if (length(missing_columns) > 0L) {
        stop(occupancy_context, " source table is missing required columns: ", paste(missing_columns, collapse = ", "))
    }
    valid_nuclear <- parse_logical_text(source_table$valid_coordinate) & parse_logical_text(source_table$primary_chromosome)
    peak_table <- source_table[valid_nuclear, , drop = FALSE]
    long_rows <- vector("list", nrow(peak_table))
    promoter_exceptions <- vector("list", nrow(peak_table))
    invalid_token_rows <- vector("list", nrow(peak_table))
    peak_qc_rows <- vector("list", nrow(peak_table))
    for (row_index in seq_len(nrow(peak_table))) {
        peak_row <- peak_table[row_index, , drop = FALSE]
        gene_tokens <- split_gene_tokens(peak_row$raw_geneOverlap)
        promoter_tokens <- split_gene_tokens(peak_row$raw_promoterOverlap)
        relation_tokens <- unique(c(gene_tokens, promoter_tokens))
        promoter_only <- setdiff(promoter_tokens, gene_tokens)
        invalid_tokens <- relation_tokens[!valid_gene_symbol_shape(relation_tokens)]
        if (length(relation_tokens) > 0L) {
            long_rows[[row_index]] <- data.frame(source_id = as.character(peak_row$source_id), occupancy_context = occupancy_context, 
                source_row = as.character(peak_row$source_row), coordinate_key = as.character(peak_row$coordinate_key), 
                chromosome = as.character(peak_row$chromosome), start0 = safe_numeric(peak_row$start0), end = safe_numeric(peak_row$end), 
                width = safe_numeric(peak_row$width), gene_symbol_author = relation_tokens, author_gene_overlap = relation_tokens %in% 
                  gene_tokens, author_promoter_overlap = relation_tokens %in% promoter_tokens, peak_relative_location = as.character(peak_row$raw_relativeLocation), 
                peak_tss_distance = safe_numeric(peak_row$raw_tssDistance), peak_distance_author = safe_numeric(peak_row$raw_distance), 
                raw_gene_overlap_value = as.character(peak_row$raw_geneOverlap), raw_promoter_overlap_value = as.character(peak_row$raw_promoterOverlap), 
                gene_linkage_provenance = ifelse(relation_tokens %in% promoter_tokens & relation_tokens %in% gene_tokens, 
                  "AUTHOR_GENE_OVERLAP_AND_PROMOTER_OVERLAP", ifelse(relation_tokens %in% promoter_tokens, "AUTHOR_PROMOTER_OVERLAP_ONLY", 
                    "AUTHOR_GENE_OVERLAP_ONLY")), stringsAsFactors = FALSE)
        }
        if (length(promoter_only) > 0L) {
            promoter_exceptions[[row_index]] <- data.frame(source_id = as.character(peak_row$source_id), occupancy_context = occupancy_context, 
                coordinate_key = as.character(peak_row$coordinate_key), promoter_only_gene = promoter_only, raw_gene_overlap_value = as.character(peak_row$raw_geneOverlap), 
                raw_promoter_overlap_value = as.character(peak_row$raw_promoterOverlap), stringsAsFactors = FALSE)
        }
        if (length(invalid_tokens) > 0L) {
            invalid_token_rows[[row_index]] <- data.frame(source_id = as.character(peak_row$source_id), occupancy_context = occupancy_context, 
                coordinate_key = as.character(peak_row$coordinate_key), invalid_gene_token = invalid_tokens, stringsAsFactors = FALSE)
        }
        peak_qc_rows[[row_index]] <- data.frame(source_id = as.character(peak_row$source_id), occupancy_context = occupancy_context, 
            coordinate_key = as.character(peak_row$coordinate_key), gene_overlap_token_n = length(gene_tokens), 
            promoter_overlap_token_n = length(promoter_tokens), relation_token_n = length(relation_tokens), promoter_only_token_n = length(promoter_only), 
            invalid_token_n = length(invalid_tokens), stringsAsFactors = FALSE)
    }
    long_table <- data.table::rbindlist(long_rows, fill = TRUE)
    promoter_exception_table <- data.table::rbindlist(promoter_exceptions, fill = TRUE)
    invalid_token_table <- data.table::rbindlist(invalid_token_rows, fill = TRUE)
    peak_qc_table <- data.table::rbindlist(peak_qc_rows, fill = TRUE)
    list(peak_source_table = peak_table, long_table = as.data.frame(long_table), promoter_exceptions = as.data.frame(promoter_exception_table), 
        invalid_tokens = as.data.frame(invalid_token_table), peak_qc = as.data.frame(peak_qc_table))
}

inherit_shared_annotations <- function(shared_table, vehicle_table, calcitriol_table) {
    shared_valid <- parse_logical_text(shared_table$valid_coordinate) & parse_logical_text(shared_table$primary_chromosome)
    vehicle_valid <- parse_logical_text(vehicle_table$valid_coordinate) & parse_logical_text(vehicle_table$primary_chromosome)
    calcitriol_valid <- parse_logical_text(calcitriol_table$valid_coordinate) & parse_logical_text(calcitriol_table$primary_chromosome)
    shared_nuclear <- shared_table[shared_valid, , drop = FALSE]
    vehicle_nuclear <- vehicle_table[vehicle_valid, , drop = FALSE]
    calcitriol_nuclear <- calcitriol_table[calcitriol_valid, , drop = FALSE]
    vehicle_match <- match(shared_nuclear$coordinate_key, vehicle_nuclear$coordinate_key)
    calcitriol_match <- match(shared_nuclear$coordinate_key, calcitriol_nuclear$coordinate_key)
    if (any(is.na(vehicle_match)) || any(is.na(calcitriol_match))) {
        stop("At least one shared nuclear peak lacks an exact coordinate match in vehicle or calcitriol.")
    }
    vehicle_annotation <- vehicle_nuclear[vehicle_match, , drop = FALSE]
    calcitriol_annotation <- calcitriol_nuclear[calcitriol_match, , drop = FALSE]
    gene_set_vehicle <- vapply(vehicle_annotation$raw_geneOverlap, canonical_gene_set_string, character(1))
    gene_set_calcitriol <- vapply(calcitriol_annotation$raw_geneOverlap, canonical_gene_set_string, character(1))
    promoter_set_vehicle <- vapply(vehicle_annotation$raw_promoterOverlap, canonical_gene_set_string, character(1))
    promoter_set_calcitriol <- vapply(calcitriol_annotation$raw_promoterOverlap, canonical_gene_set_string, character(1))
    relative_location_vehicle <- trimws(as.character(vehicle_annotation$raw_relativeLocation))
    relative_location_calcitriol <- trimws(as.character(calcitriol_annotation$raw_relativeLocation))
    tss_vehicle <- safe_numeric(vehicle_annotation$raw_tssDistance)
    tss_calcitriol <- safe_numeric(calcitriol_annotation$raw_tssDistance)
    audit <- data.frame(coordinate_key = shared_nuclear$coordinate_key, vehicle_gene_set = gene_set_vehicle, calcitriol_gene_set = gene_set_calcitriol, 
        gene_set_exact = gene_set_vehicle == gene_set_calcitriol, vehicle_promoter_set = promoter_set_vehicle, 
        calcitriol_promoter_set = promoter_set_calcitriol, promoter_set_exact = promoter_set_vehicle == promoter_set_calcitriol, 
        vehicle_relative_location = relative_location_vehicle, calcitriol_relative_location = relative_location_calcitriol, 
        relative_location_exact = relative_location_vehicle == relative_location_calcitriol, vehicle_tss_distance = tss_vehicle, 
        calcitriol_tss_distance = tss_calcitriol, tss_distance_exact = (is.na(tss_vehicle) & is.na(tss_calcitriol)) | 
            (is.finite(tss_vehicle) & is.finite(tss_calcitriol) & tss_vehicle == tss_calcitriol), stringsAsFactors = FALSE)
    long_rows <- vector("list", nrow(shared_nuclear))
    for (row_index in seq_len(nrow(shared_nuclear))) {
        vehicle_gene_tokens <- split_gene_tokens(vehicle_annotation$raw_geneOverlap[row_index])
        calcitriol_gene_tokens <- split_gene_tokens(calcitriol_annotation$raw_geneOverlap[row_index])
        vehicle_promoter_tokens <- split_gene_tokens(vehicle_annotation$raw_promoterOverlap[row_index])
        calcitriol_promoter_tokens <- split_gene_tokens(calcitriol_annotation$raw_promoterOverlap[row_index])
        relation_tokens <- unique(c(vehicle_gene_tokens, calcitriol_gene_tokens, vehicle_promoter_tokens, calcitriol_promoter_tokens))
        if (length(relation_tokens) == 0L) {
            next
        }
        long_rows[[row_index]] <- data.frame(source_id = "VDR_SHARED_CISTROME", occupancy_context = "shared_vehicle_and_calcitriol", 
            source_row = as.character(shared_nuclear$source_row[row_index]), coordinate_key = as.character(shared_nuclear$coordinate_key[row_index]), 
            chromosome = as.character(shared_nuclear$chromosome[row_index]), start0 = safe_numeric(shared_nuclear$start0[row_index]), 
            end = safe_numeric(shared_nuclear$end[row_index]), width = safe_numeric(shared_nuclear$width[row_index]), 
            gene_symbol_author = relation_tokens, vehicle_gene_overlap = relation_tokens %in% vehicle_gene_tokens, 
            calcitriol_gene_overlap = relation_tokens %in% calcitriol_gene_tokens, vehicle_promoter_overlap = relation_tokens %in% 
                vehicle_promoter_tokens, calcitriol_promoter_overlap = relation_tokens %in% calcitriol_promoter_tokens, 
            inherited_gene_set_exact = audit$gene_set_exact[row_index], inherited_promoter_set_exact = audit$promoter_set_exact[row_index], 
            peak_relative_location = vehicle_annotation$raw_relativeLocation[row_index], peak_tss_distance = safe_numeric(vehicle_annotation$raw_tssDistance[row_index]), 
            annotation_inheritance_provenance = "EXACT_COORDINATE_MATCH_TO_VEHICLE_AND_CALCITRIOL", stringsAsFactors = FALSE)
    }
    shared_long <- data.table::rbindlist(long_rows, fill = TRUE)
    list(audit = audit, long_table = as.data.frame(shared_long))
}

log_event("Starting", run_id)

required_paths <- c(gate_c1 = file.path(dirs$results, "STEP03B3C1_GATE_TABLE.csv"), gene_decision_c1 = file.path(dirs$results, 
    "STEP03B3C1_GENE_FIELD_DECISION_TABLE.csv"), coverage_c1 = file.path(dirs$results, "STEP03B3C1_FIELD_ROLE_COVERAGE_MATRIX.csv"), 
    source_summary_c1 = file.path(dirs$qc, "STEP03B3C1_SOURCE_FIELD_SUMMARY.csv"), column_audit_c1 = file.path(dirs$qc, 
        "STEP03B3C1_ANNOTATION_COLUMN_AUDIT.csv"), source_lock = file.path(dirs$results, "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"), 
    RNA_lock_card = file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C2 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_c1 <- .vdr_read_csv(required_paths[["gate_c1"]], stringsAsFactors = FALSE, check.names = FALSE)

gene_decision_c1 <- .vdr_read_csv(required_paths[["gene_decision_c1"]], stringsAsFactors = FALSE, check.names = FALSE)

coverage_c1 <- .vdr_read_csv(required_paths[["coverage_c1"]], stringsAsFactors = FALSE, check.names = FALSE)

source_summary_c1 <- .vdr_read_csv(required_paths[["source_summary_c1"]], stringsAsFactors = FALSE, check.names = FALSE)

column_audit_c1 <- .vdr_read_csv(required_paths[["column_audit_c1"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock <- .vdr_read_csv(required_paths[["source_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- gate_c1$status[gate_c1$gate_component == "STEP03B3C1 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_PEAK_ANNOTATION_FIELD_LOCK") {
    stop("STEP03B3C1 does not authorize annotation-field lock. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (nrow(gene_decision_c1) != 2L || !all(gene_decision_c1$column_name == "raw_geneOverlap") || !setequal(gene_decision_c1$source_id, 
    c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME"))) {
    stop("STEP03B3C1 did not identify the expected two raw_geneOverlap fields.")
}

expected_field_rows <- data.frame(source_id = rep(c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME"), each = 3L), 
    column_name = rep(c("raw_geneOverlap", "raw_promoterOverlap", "raw_tssDistance"), times = 2L), stringsAsFactors = FALSE)

field_check <- merge(expected_field_rows, column_audit_c1[c("source_id", "column_name", "nonempty_n")], by = c("source_id", 
    "column_name"), all.x = TRUE, sort = FALSE)

if (nrow(field_check) != 6L || any(is.na(field_check$nonempty_n)) || any(field_check$nonempty_n == 0L)) {
    stop("Expected gene, promoter and TSS-distance annotation fields are incomplete.")
}

if (nrow(source_summary_c1) != 3L || !all(source_summary_c1$primary_gene_field_present[source_summary_c1$source_id %in% 
    c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME")]) || source_summary_c1$primary_gene_field_present[source_summary_c1$source_id == 
    "VDR_SHARED_CISTROME"]) {
    stop("Source-level gene-field summary does not match the reviewed schema.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

primary_sources <- source_lock[source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", 
    , drop = FALSE]

if (nrow(primary_sources) != 3L || !all(primary_sources$source_lock_version == source_lock_version) || !all(primary_sources$source_lock_status == 
    "LOCKED")) {
    stop("Primary VDR cistrome source lock is incomplete.")
}

primary_sources$file_exists <- file.exists(primary_sources$clean_table_path)

primary_sources$observed_md5 <- vapply(primary_sources$clean_table_path, md5_file, character(1))

primary_sources$observed_sha256 <- vapply(primary_sources$clean_table_path, sha256_file, character(1))

primary_sources$identity_pass <- with(primary_sources, file_exists & clean_table_md5 == observed_md5 & clean_table_sha256 == 
    observed_sha256)

if (!all(primary_sources$identity_pass)) {
    stop("One or more locked source tables failed hash revalidation.")
}

source_paths <- setNames(primary_sources$clean_table_path, primary_sources$source_id)

source_tables <- lapply(source_paths, function(path_value) {
    .vdr_fread(path_value, sep = "\t", colClasses = "character", data.table = FALSE, showProgress = FALSE, na.strings = c("", 
        "NA", "N/A", "NULL"))
})

vehicle_result <- build_condition_peak_gene_map(source_tables[["VDR_VEHICLE_CISTROME"]], "vehicle")

calcitriol_result <- build_condition_peak_gene_map(source_tables[["VDR_CALCITRIOL_CISTROME"]], "calcitriol")

shared_result <- inherit_shared_annotations(shared_table = source_tables[["VDR_SHARED_CISTROME"]], vehicle_table = source_tables[["VDR_VEHICLE_CISTROME"]], 
    calcitriol_table = source_tables[["VDR_CALCITRIOL_CISTROME"]])

peak_gene_long <- data.table::rbindlist(list(vehicle_result$long_table, calcitriol_result$long_table, shared_result$long_table), 
    fill = TRUE)

peak_gene_long <- as.data.frame(peak_gene_long, stringsAsFactors = FALSE)

peak_gene_long$peak_gene_relation_id <- paste0(peak_gene_long$source_id, "::", peak_gene_long$coordinate_key, "::", 
    peak_gene_long$gene_symbol_author)

peak_gene_long <- peak_gene_long[!duplicated(peak_gene_long$peak_gene_relation_id), , drop = FALSE]

peak_gene_long <- peak_gene_long[order(peak_gene_long$source_id, peak_gene_long$chromosome, peak_gene_long$start0, 
    peak_gene_long$gene_symbol_author), , drop = FALSE]

promoter_exceptions <- data.table::rbindlist(list(vehicle_result$promoter_exceptions, calcitriol_result$promoter_exceptions), 
    fill = TRUE)

invalid_tokens <- data.table::rbindlist(list(vehicle_result$invalid_tokens, calcitriol_result$invalid_tokens), 
    fill = TRUE)

peak_qc <- data.table::rbindlist(list(vehicle_result$peak_qc, calcitriol_result$peak_qc), fill = TRUE)

promoter_exceptions <- as.data.frame(promoter_exceptions)

invalid_tokens <- as.data.frame(invalid_tokens)

peak_qc <- as.data.frame(peak_qc)

shared_audit <- shared_result$audit

shared_gene_concordance <- mean(shared_audit$gene_set_exact)

shared_promoter_concordance <- mean(shared_audit$promoter_set_exact)

shared_relative_location_concordance <- mean(shared_audit$relative_location_exact)

shared_tss_distance_concordance <- mean(shared_audit$tss_distance_exact)

total_gene_tokens <- nrow(peak_gene_long)

invalid_gene_token_fraction <- if (total_gene_tokens > 0L) {
    nrow(invalid_tokens)/total_gene_tokens
} else {
    NA_real_
}

promoter_subset_pass <- nrow(promoter_exceptions) == 0L

shared_concordance_pass <- all(c(shared_gene_concordance, shared_promoter_concordance) >= minimum_shared_annotation_concordance)

token_shape_pass <- is.finite(invalid_gene_token_fraction) && invalid_gene_token_fraction <= maximum_unparseable_gene_token_fraction

mapping_qc_summary <- data.frame(metric = c("Vehicle nuclear peaks", "Calcitriol nuclear peaks", "Shared nuclear peaks", 
    "Vehicle peak-gene relations", "Calcitriol peak-gene relations", "Shared inherited peak-gene relations", "Unique candidate gene symbols", 
    "Promoter-only relation exceptions", "Invalid gene-token rows", "Invalid gene-token fraction", "Shared gene-set concordance", 
    "Shared promoter-set concordance", "Shared relative-location concordance", "Shared TSS-distance concordance"), 
    value = c(nrow(vehicle_result$peak_source_table), nrow(calcitriol_result$peak_source_table), nrow(shared_audit), 
        sum(peak_gene_long$source_id == "VDR_VEHICLE_CISTROME"), sum(peak_gene_long$source_id == "VDR_CALCITRIOL_CISTROME"), 
        sum(peak_gene_long$source_id == "VDR_SHARED_CISTROME"), length(unique(peak_gene_long$gene_symbol_author)), 
        nrow(promoter_exceptions), nrow(invalid_tokens), invalid_gene_token_fraction, shared_gene_concordance, 
        shared_promoter_concordance, shared_relative_location_concordance, shared_tss_distance_concordance), stringsAsFactors = FALSE)

write_csv_utf8(mapping_qc_summary, file.path(dirs$qc, "STEP03B3C2_MAPPING_QC_SUMMARY.csv"))

write_csv_utf8(shared_audit, file.path(dirs$qc, "STEP03B3C2_SHARED_ANNOTATION_INHERITANCE_AUDIT.csv"))

write_csv_utf8(promoter_exceptions, file.path(dirs$qc, "STEP03B3C2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"))

write_csv_utf8(invalid_tokens, file.path(dirs$qc, "STEP03B3C2_INVALID_GENE_TOKEN_ROWS.csv"))

write_csv_utf8(peak_qc, file.path(dirs$qc, "STEP03B3C2_PEAK_LEVEL_MAPPING_QC.csv"))

peak_gene_long_path <- file.path(map_root, "STEP03B3C2_PEAK_GENE_LONG_CANDIDATE.tsv.gz")

write_tsv_gz(peak_gene_long, peak_gene_long_path)

gene_level_preview <- aggregate(coordinate_key ~ source_id + occupancy_context + gene_symbol_author, data = peak_gene_long, 
    FUN = function(x) {
        length(unique(x))
    })

names(gene_level_preview)[names(gene_level_preview) == "coordinate_key"] <- "supporting_peak_n"

promoter_summary <- aggregate(cbind(author_gene_overlap = as.integer(peak_gene_long$author_gene_overlap %in% TRUE | 
    peak_gene_long$vehicle_gene_overlap %in% TRUE | peak_gene_long$calcitriol_gene_overlap %in% TRUE), author_promoter_overlap = as.integer(peak_gene_long$author_promoter_overlap %in% 
    TRUE | peak_gene_long$vehicle_promoter_overlap %in% TRUE | peak_gene_long$calcitriol_promoter_overlap %in% 
    TRUE)) ~ source_id + occupancy_context + gene_symbol_author, data = peak_gene_long, FUN = max, na.rm = TRUE)

gene_level_preview <- merge(gene_level_preview, promoter_summary, by = c("source_id", "occupancy_context", "gene_symbol_author"), 
    all.x = TRUE, sort = FALSE)

gene_level_preview_path <- file.path(map_root, "STEP03B3C2_GENE_LEVEL_EVIDENCE_PREVIEW.csv")

write_csv_utf8(gene_level_preview, gene_level_preview_path)

field_semantics_lock <- data.frame(field_name = c("raw_geneOverlap", "raw_promoterOverlap", "raw_relativeLocation", 
    "raw_tssDistance", "raw_distance", "VDR_SHARED_CISTROME annotation"), locked_semantics = c("AUTHOR_REPORTED_GENE_OVERLAP_SYMBOL_SET", 
    "AUTHOR_REPORTED_PROMOTER_OVERLAP_SYMBOL_SET", "PEAK_LEVEL_GENOMIC_REGION_CLASS", "PEAK_LEVEL_NUMERIC_TSS_DISTANCE_WITHOUT_AUTHORITATIVE_GENE_IDENTITY", 
    "AUTHOR_NUMERIC_DISTANCE_FIELD_NOT_USED_FOR_GENE_LINKAGE", "INHERITED_BY_EXACT_COORDINATE_MATCH_FROM_VEHICLE_AND_CALCITRIOL"), 
    parsing_rule = c("Trim; treat blank/dot/NA as missing; split on one-or-more whitespace; deduplicate within peak.", 
        "Trim; treat blank/dot/NA as no promoter symbol; split on one-or-more whitespace; deduplicate within peak.", 
        "Retain raw category as peak-level annotation; do not assign it separately to each listed gene.", "Parse numeric and retain at peak level; do not create a nearest-gene relationship.", 
        "Parse numeric and retain at peak level only.", "Require exact coordinate match to both source tables; audit gene and promoter set concordance."), 
    allowed_use = c("Create author gene-overlap peak-gene links.", "Flag promoter-supported genes within each peak.", 
        "Describe peak genomic context.", "Describe peak-level distance only.", "Descriptive provenance only.", 
        "Create transparent shared occupancy links after concordance QC."), prohibited_use = c("Cannot alone prove direct regulation.", 
        "Cannot be merged silently with generic gene overlap.", "Cannot serve as a gene identifier.", "Cannot define a nearest gene because no linked gene identity is supplied.", 
        "Cannot create a peak-gene relationship.", "Cannot inherit silently or override discordant source annotations."), 
    status = "LOCKED_FOR_CANDIDATE_MAP", stringsAsFactors = FALSE)

write_csv_utf8(field_semantics_lock, file.path(dirs$results, "STEP03B3C2_FIELD_SEMANTICS_LOCK.csv"))

candidate_manifest <- data.frame(artifact_id = c("PEAK_GENE_LONG_CANDIDATE", "GENE_LEVEL_EVIDENCE_PREVIEW"), version = c(peak_gene_map_candidate_version, 
    peak_gene_map_candidate_version), path = normalize_slashes(c(peak_gene_long_path, gene_level_preview_path)), 
    size_bytes = as.numeric(file.info(c(peak_gene_long_path, gene_level_preview_path))$size), md5 = vapply(c(peak_gene_long_path, 
        gene_level_preview_path), md5_file, character(1)), sha256 = vapply(c(peak_gene_long_path, gene_level_preview_path), 
        sha256_file, character(1)), lock_status = "CANDIDATE_PENDING_MANUAL_QC", stringsAsFactors = FALSE)

write_csv_utf8(candidate_manifest, file.path(map_root, "STEP03B3C2_CANDIDATE_MAP_MANIFEST.csv"))

gate_status <- if (shared_concordance_pass && promoter_subset_pass && token_shape_pass && nrow(peak_gene_long) > 
    0L) {
    "PROCEED_TO_MANUAL_PEAK_GENE_MAP_LOCK"
} else if (!shared_concordance_pass) {
    "HOLD_SHARED_ANNOTATION_DISCORDANCE"
} else if (!promoter_subset_pass) {
    "HOLD_PROMOTER_GENE_RELATION_EXCEPTIONS"
} else if (!token_shape_pass) {
    "HOLD_GENE_TOKEN_PARSE_QC"
} else {
    "HOLD_EMPTY_PEAK_GENE_MAP"
}

gate_table <- data.frame(gate_component = c("STEP03B3C1 authorization", "Published cistrome source lock", "RNA membership", 
    "ATAC role", "Field semantics", "Vehicle peak-gene map", "Calcitriol peak-gene map", "Shared annotation inheritance", 
    "Promoter relation consistency", "Gene token shape", "Candidate map version", "RNA/CUT&RUN integration", "Definitive direct-target assignment", 
    "Scientific figures", "STEP03B3C2 gate"), status = c("PASS", paste0("IMMUTABLE_", source_lock_version), paste0("IMMUTABLE_", 
    RNA_anchor), ATAC_role, "LOCKED_FOR_CANDIDATE_MAP", "CREATED", "CREATED", ifelse(shared_concordance_pass, "PASS_EXACT", 
    "FAIL"), ifelse(promoter_subset_pass, "PASS", "REVIEW"), ifelse(token_shape_pass, "PASS", "REVIEW"), paste0(peak_gene_map_candidate_version, 
    "_PENDING_MANUAL_QC"), "NOT_STARTED", "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain the field-audit record.", 
    "Do not alter source peaks or coordinates.", "Do not change the 436 RNA interaction candidates.", "Use only as optional accessible context.", 
    "Apply the six locked field interpretations.", "Review relation counts and promoter flags.", "Review relation counts and promoter flags.", 
    "Confirm vehicle/calcitriol annotation concordance for shared peaks.", "Review any promoter-only genes before locking.", 
    "Review malformed tokens before locking.", "Lock only after manual QC.", "Do not integrate until the candidate map is locked.", 
    "Use occupancy-supported candidate terminology only after integration QC.", "No Scientific Figure Skill 3 workflow is triggered.", 
    "Upload candidate manifest and QC tables."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C2_GATE_TABLE.csv"))

write_text(c("# STEP03B3C2 peak-gene map candidate card", "", paste0("- Run ID: ", run_id), paste0("- Source lock: ", 
    source_lock_version, " IMMUTABLE"), paste0("- Candidate version: ", peak_gene_map_candidate_version), paste0("- Peak-gene relation rows: ", 
    nrow(peak_gene_long)), paste0("- Unique author gene symbols: ", length(unique(peak_gene_long$gene_symbol_author))), 
    paste0("- Shared gene-set concordance: ", shared_gene_concordance), paste0("- Shared promoter-set concordance: ", 
        shared_promoter_concordance), paste0("- Promoter-only relation exceptions: ", nrow(promoter_exceptions)), 
    paste0("- Invalid gene-token rows: ", nrow(invalid_tokens)), "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C2_PEAK_GENE_MAP_CANDIDATE_CARD.md"))

write_text(c("# STEP03B3C2 Methods draft", "", "Author-provided `geneOverlap` and `promoterOverlap` fields were interpreted as whitespace-delimited gene-symbol sets after trimming blank, missing and dot placeholder values. One-to-many peak-gene relationships were expanded into a long-form table, while promoter overlap was retained as a gene-specific evidence flag rather than merged with generic gene overlap.", 
    "", "The `relativeLocation`, `tssDistance` and `distance` fields were retained as peak-level annotations. Because the worksheets did not supply an explicit gene identity linked to the TSS-distance field, TSS distance was not used to construct nearest-gene relationships.", 
    "", "The coordinate-only shared cistrome table inherited annotation through exact coordinate matches to both vehicle and calcitriol primary sources. Gene and promoter symbol sets were compared before inheritance. No integration with the frozen RNA program was performed."), 
    file.path(dirs$text, "STEP03B3C2_Methods_draft.md"))

write_text(c("# STEP03B3C2 interpretation constraints", "", "- The long table records author-provided overlap relationships.", 
    "- A gene-overlap row is not equivalent to direct transcriptional regulation.", "- Promoter overlap is retained separately and is stronger positional support than generic overlap, but remains associative.", 
    "- Peak-level TSS distance cannot define a nearest-gene link without an explicit associated gene identity.", 
    "- Shared annotations are inherited only through exact coordinate matches and concordance QC.", "- RNA_PROGRAM_v1.0 membership remains unchanged.", 
    "- No direct-target terminology is authorized before RNA integration and evidence-tier QC."), file.path(dirs$text, 
    "STEP03B3C2_Interpretation_constraints.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, source_lock_version = source_lock_version, peak_gene_map_candidate_version = peak_gene_map_candidate_version, 
    RNA_anchor = RNA_anchor, ATAC_role = ATAC_role, split_rule = "one_or_more_whitespace", missing_gene_tokens = c("", 
        ".", "NA", "N/A", "NULL"), shared_gene_concordance = shared_gene_concordance, shared_promoter_concordance = shared_promoter_concordance, 
    peak_gene_relation_rows = nrow(peak_gene_long), RNA_CUTRUN_integration = FALSE, definitive_direct_target_assignment = FALSE, 
    scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("38_step03B3C2_lock_annotation_semantics_build_peak_gene_map.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "build_peak_gene_map.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C2_FIELD_SEMANTICS_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3C2_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C2_PEAK_GENE_MAP_CANDIDATE_CARD.md"), peak_gene_long_path, 
    gene_level_preview_path, file.path(map_root, "STEP03B3C2_CANDIDATE_MAP_MANIFEST.csv"), file.path(dirs$qc, "STEP03B3C2_MAPPING_QC_SUMMARY.csv"), 
    file.path(dirs$qc, "STEP03B3C2_SHARED_ANNOTATION_INHERITANCE_AUDIT.csv"), file.path(dirs$qc, "STEP03B3C2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv"), 
    file.path(dirs$qc, "STEP03B3C2_INVALID_GENE_TOKEN_ROWS.csv"), file.path(dirs$qc, "STEP03B3C2_PEAK_LEVEL_MAPPING_QC.csv"), 
    file.path(dirs$text, "STEP03B3C2_Methods_draft.md"), file.path(dirs$text, "STEP03B3C2_Interpretation_constraints.md"), 
    file.path(dirs$parameters, "STEP03B3C2_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C2.txt"), 
    file.path(dirs$code, "build_peak_gene_map.R"), log_path)

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
    else if (grepl("Methods|Interpretation", basename(path_value))) {
        "txt"
    }
    else if (grepl("QC|AUDIT|EXCEPTION|INVALID|SUMMARY|CARD|SEMANTICS|MANIFEST", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C2 archive status", "", paste0("- Run ID: ", run_id), paste0("- Candidate map: ", peak_gene_map_candidate_version), 
    paste0("- Peak-gene relations: ", nrow(peak_gene_long)), paste0("- Shared gene concordance: ", shared_gene_concordance), 
    paste0("- Shared promoter concordance: ", shared_promoter_concordance), paste0("- Promoter-only exceptions: ", 
        nrow(promoter_exceptions)), paste0("- Invalid gene-token rows: ", nrow(invalid_tokens)), "- RNA/CUT&RUN integration: NOT STARTED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C2", paste0("- Run ID: ", run_id), paste0("- Published source lock: ", 
    source_lock_version, " IMMUTABLE"), paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE"), paste0("- ATAC role: ", 
    ATAC_role), "- Annotation field semantics: LOCKED_FOR_CANDIDATE_MAP", paste0("- Peak-gene map: ", peak_gene_map_candidate_version, 
    " PENDING_MANUAL_QC"), "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, 
    "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C2 annotation semantics and peak-gene map completed.\n")

cat("Source lock: ", source_lock_version, " IMMUTABLE\n")

cat("Candidate map: ", peak_gene_map_candidate_version, "\n")

cat("Peak-gene relation rows: ", nrow(peak_gene_long), "\n")

cat("Unique author gene symbols: ", length(unique(peak_gene_long$gene_symbol_author)), "\n")

cat("Shared gene-set concordance: ", shared_gene_concordance, "\n")

cat("Shared promoter-set concordance: ", shared_promoter_concordance, "\n")

cat("Promoter-only relation exceptions: ", nrow(promoter_exceptions), "\n")

cat("Invalid gene-token rows: ", nrow(invalid_tokens), "\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C2_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3C2_FIELD_SEMANTICS_LOCK.csv\n")

cat("  08_qc/STEP03B3C2_MAPPING_QC_SUMMARY.csv\n")

cat("  08_qc/STEP03B3C2_SHARED_ANNOTATION_INHERITANCE_AUDIT.csv\n")

cat("  08_qc/STEP03B3C2_PROMOTER_GENE_SUBSET_EXCEPTIONS.csv\n")

cat("  08_qc/STEP03B3C2_INVALID_GENE_TOKEN_ROWS.csv\n")

cat("  06_results/STEP03B3C2_peak_gene_map_candidate/\n")

cat("============================================================\n")

