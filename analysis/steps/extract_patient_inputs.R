#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B02: Extract count archives and sample metadata

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(3600L, getOption("timeout")))

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04B_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "COMBINED_INPUT_SOURCE_v1.0"

sample_design_candidate_version <- "COMBINED_SAMPLE_DESIGN_CANDIDATE_v1.0"

expected_source_rows <- c(GSE250130 = 2L, GSE165004 = 7L, GSE92324 = 2L, GSE111974 = 3L, GSE183837 = 28L)

required_packages <- c("data.table", "digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(step04_root, "00_module_control"), inputs = file.path(step04_root, "01_inputs"), 
    code = file.path(step04_root, "02_code"), parameters = file.path(step04_root, "03_parameters"), environment = file.path(step04_root, 
        "04_environment"), execution = file.path(step04_root, "05_execution"), results = file.path(step04_root, 
        "06_results"), figures = file.path(step04_root, "07_figures_tables"), qc = file.path(step04_root, "08_qc"), 
    text = file.path(step04_root, "09_manuscript_text"), archive_parent = file.path(step04_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

locked_root <- file.path(dirs$results, "STEP04B_locked_combined_inputs", source_lock_version)

design_root <- file.path(dirs$results, "STEP04B_sample_design_candidates")

dir.create(locked_root, recursive = TRUE, showWarnings = FALSE)

dir.create(design_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S4B_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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
        stop("Copy/hash validation failed:\n", source_path)
    }
    TRUE
}

parse_GEO_series_matrix <- function(path) {
    if (!file.exists(path)) {
        stop("Series matrix is missing: ", path)
    }
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    all_lines <- readLines(connection, warn = FALSE, encoding = "UTF-8")
    sample_line_indices <- grep("^!Sample_", all_lines)
    metadata_rows <- list()
    metadata_index <- 1L
    field_occurrence_counter <- list()
    for (line_index in sample_line_indices) {
        fields <- strsplit(all_lines[[line_index]], "\t", fixed = TRUE)[[1]]
        field_name <- sub("^!Sample_", "", fields[[1]])
        occurrence <- if (is.null(field_occurrence_counter[[field_name]])) {
            1L
        }
        else {
            field_occurrence_counter[[field_name]] + 1L
        }
        field_occurrence_counter[[field_name]] <- occurrence
        values <- fields[-1]
        values <- gsub("^\"|\"$", "", values)
        metadata_rows[[metadata_index]] <- data.frame(field_name = field_name, field_occurrence = occurrence, sample_index = seq_along(values), 
            value = values, stringsAsFactors = FALSE)
        metadata_index <- metadata_index + 1L
    }
    metadata_long <- if (length(metadata_rows) > 0L) {
        do.call(rbind, metadata_rows)
    }
    else {
        data.frame(field_name = character(0), field_occurrence = integer(0), sample_index = integer(0), value = character(0), 
            stringsAsFactors = FALSE)
    }
    table_begin <- grep("^!series_matrix_table_begin", all_lines)
    table_end <- grep("^!series_matrix_table_end", all_lines)
    expression_feature_n <- 0L
    expression_sample_headers <- character(0)
    if (length(table_begin) == 1L && length(table_end) == 1L && table_end > table_begin + 1L) {
        expression_feature_n <- table_end - table_begin - 2L
        header_line <- all_lines[table_begin + 1L]
        header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
        header_fields <- gsub("^\"|\"$", "", header_fields)
        expression_sample_headers <- header_fields[-1]
    }
    geo_accessions <- metadata_long$value[metadata_long$field_name == "geo_accession"]
    geo_accessions <- unique(geo_accessions[nzchar(geo_accessions)])
    list(metadata_long = metadata_long, sample_accessions = geo_accessions, sample_n = length(geo_accessions), 
        expression_feature_n = expression_feature_n, expression_sample_headers = expression_sample_headers)
}

metadata_to_sample_sheet <- function(dataset_id, parsed_series) {
    metadata_long <- parsed_series$metadata_long
    if (nrow(metadata_long) == 0L || parsed_series$sample_n == 0L) {
        return(data.frame())
    }
    geo_rows <- metadata_long[metadata_long$field_name == "geo_accession", , drop = FALSE]
    if (nrow(geo_rows) == 0L) {
        stop("No GEO accession row was found for ", dataset_id)
    }
    sample_lookup <- data.frame(sample_index = geo_rows$sample_index, sample_accession = geo_rows$value, stringsAsFactors = FALSE)
    selected_fields <- metadata_long[metadata_long$field_name %in% c("title", "source_name_ch1", "organism_ch1", 
        "characteristics_ch1", "description", "relation"), , drop = FALSE]
    selected_fields$field_key <- paste0(selected_fields$field_name, "__", selected_fields$field_occurrence)
    sample_sheet <- sample_lookup
    for (field_key in unique(selected_fields$field_key)) {
        field_rows <- selected_fields[selected_fields$field_key == field_key, c("sample_index", "value"), drop = FALSE]
        names(field_rows)[names(field_rows) == "value"] <- field_key
        sample_sheet <- merge(sample_sheet, field_rows, by = "sample_index", all.x = TRUE, sort = FALSE)
    }
    sample_sheet <- sample_sheet[order(sample_sheet$sample_index), , drop = FALSE]
    sample_sheet$dataset_id <- dataset_id
    sample_sheet
}

characteristics_to_key_value <- function(dataset_id, sample_sheet) {
    if (nrow(sample_sheet) == 0L) {
        return(data.frame(dataset_id = character(0), sample_accession = character(0), characteristic_field = character(0), 
            raw_value = character(0), parsed_key = character(0), parsed_value = character(0), stringsAsFactors = FALSE))
    }
    characteristic_columns <- grep("^characteristics_ch1__", names(sample_sheet), value = TRUE)
    rows <- list()
    row_index <- 1L
    for (column_name in characteristic_columns) {
        for (sample_index in seq_len(nrow(sample_sheet))) {
            raw_value <- as.character(sample_sheet[sample_index, column_name])
            if (is.na(raw_value) || !nzchar(trimws(raw_value))) {
                next
            }
            colon_position <- regexpr(":", raw_value, fixed = TRUE)[[1]]
            if (colon_position > 0L) {
                parsed_key <- trimws(substr(raw_value, 1L, colon_position - 1L))
                parsed_value <- trimws(substr(raw_value, colon_position + 1L, nchar(raw_value)))
            }
            else {
                parsed_key <- NA_character_
                parsed_value <- raw_value
            }
            rows[[row_index]] <- data.frame(dataset_id = dataset_id, sample_accession = sample_sheet$sample_accession[sample_index], 
                characteristic_field = column_name, raw_value = raw_value, parsed_key = parsed_key, parsed_value = parsed_value, 
                stringsAsFactors = FALSE)
            row_index <- row_index + 1L
        }
    }
    if (length(rows) == 0L) {
        return(data.frame(dataset_id = character(0), sample_accession = character(0), characteristic_field = character(0), 
            raw_value = character(0), parsed_key = character(0), parsed_value = character(0), stringsAsFactors = FALSE))
    }
    do.call(rbind, rows)
}

identify_triplet_component <- function(path) {
    file_name <- tolower(basename(path))
    if (grepl("matrix.*\\.mtx(\\.gz)?$|counts.*\\.mtx(\\.gz)?$", file_name, perl = TRUE)) {
        return("matrix")
    }
    if (grepl("barcodes.*\\.(tsv|txt)(\\.gz)?$", file_name, perl = TRUE)) {
        return("barcodes")
    }
    if (grepl("(features|genes).*\\.(tsv|txt)(\\.gz)?$", file_name, perl = TRUE)) {
        return("features_or_genes")
    }
    "other"
}

extract_GSM <- function(path) {
    match_value <- regmatches(basename(path), regexpr("GSM[0-9]+", basename(path), ignore.case = TRUE, perl = TRUE))
    if (length(match_value) == 0L || !nzchar(match_value)) {
        return(NA_character_)
    }
    toupper(match_value)
}

log_event("Starting", run_id)

step04A_gate_path <- file.path(dirs$results, "STEP04A_v2_GATE_TABLE.csv")

source_candidate_path <- file.path(dirs$results, "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv")

role_lock_path <- file.path(dirs$results, "STEP04A_v2_DATASET_ROLE_LOCK.csv")

slide_lock_path <- file.path(dirs$results, "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv")

required_upstream <- c(step04A_gate_path, source_candidate_path, role_lock_path, slide_lock_path)

missing_upstream <- required_upstream[!file.exists(required_upstream)]

if (length(missing_upstream) > 0L) {
    stop("Required STEP04A v2 output is missing:\n", paste(missing_upstream, collapse = "\n"))
}

step04A_gate <- .vdr_read_csv(step04A_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

source_candidate <- .vdr_read_csv(source_candidate_path, stringsAsFactors = FALSE, check.names = FALSE)

role_lock <- .vdr_read_csv(role_lock_path, stringsAsFactors = FALSE, check.names = FALSE)

slide_lock <- .vdr_read_csv(slide_lock_path, stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- step04A_gate$status[step04A_gate$gate_component == "STEP04A v2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_STEP04_INPUT_SOURCE_LOCK") {
    stop("STEP04A v2 does not authorize source lock.")
}

if (nrow(source_candidate) != sum(expected_source_rows) || any(!source_candidate$file_exists) || any(is.na(source_candidate$md5) | 
    is.na(source_candidate$sha256))) {
    stop("The 42-row source candidate is incomplete.")
}

observed_source_rows <- table(source_candidate$dataset_id)

if (!setequal(names(observed_source_rows), names(expected_source_rows)) || !all(as.integer(observed_source_rows[names(expected_source_rows)]) == 
    as.integer(expected_source_rows))) {
    stop("Dataset-specific source-row counts differ from 2/7/2/3/28.")
}

if (nrow(role_lock) != 5L || role_lock$analysis_role[role_lock$dataset_id == "GSE250130"] != "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT") {
    stop("The corrected five-dataset role lock failed.")
}

if (nrow(slide_lock) != 4L || any(!slide_lock$file_exists)) {
    stop("The four-file GSE165004 slide-QC lock failed.")
}

source_candidate$current_exists <- file.exists(source_candidate$path)

source_candidate$current_size_bytes <- as.numeric(file.info(source_candidate$path)$size)

source_candidate$size_identity_pass <- source_candidate$current_exists & source_candidate$current_size_bytes == 
    source_candidate$size_bytes

source_candidate$hash_recheck_mode <- ifelse(source_candidate$size_bytes <= 500 * 1024^2, "MD5_SHA256_RECHECKED", 
    "LARGE_FILE_SIZE_AND_ARCHIVE_STRUCTURE_RECHECKED")

source_candidate$observed_md5 <- NA_character_

source_candidate$observed_sha256 <- NA_character_

small_file_index <- source_candidate$size_bytes <= 500 * 1024^2

source_candidate$observed_md5[small_file_index] <- vapply(source_candidate$path[small_file_index], md5_file, character(1))

source_candidate$observed_sha256[small_file_index] <- vapply(source_candidate$path[small_file_index], sha256_file, 
    character(1))

source_candidate$hash_identity_pass <- TRUE

source_candidate$hash_identity_pass[small_file_index] <- with(source_candidate[small_file_index, , drop = FALSE], 
    md5 == observed_md5 & sha256 == observed_sha256)

if (!all(source_candidate$size_identity_pass) || !all(source_candidate$hash_identity_pass)) {
    stop("At least one combined source reference failed identity validation.")
}

write_csv_utf8(source_candidate, file.path(dirs$qc, "STEP04B_SOURCE_IDENTITY_RECHECK.csv"))

GSE250130_RAW_path <- source_candidate$path[source_candidate$dataset_id == "GSE250130" & source_candidate$source_role == 
    "RAW_OUTER_TAR_OF_SAMPLE_BUNDLES"]

if (length(GSE250130_RAW_path) != 1L) {
    stop("GSE250130 outer TAR could not be uniquely identified.")
}

GSE250130_extract_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE250130", "raw_counts", 
    "outer_extracted")

GSE250130_sample_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE250130", "raw_counts", 
    "mm_by_sample")

dir.create(GSE250130_extract_root, recursive = TRUE, showWarnings = FALSE)

dir.create(GSE250130_sample_root, recursive = TRUE, showWarnings = FALSE)

outer_marker <- file.path(GSE250130_extract_root, "_OUTER_EXTRACTION_COMPLETE.txt")

if (!file.exists(outer_marker)) {
    log_event("Extracting GSE250130 outer TAR.")
    extraction_status <- tryCatch({
        utils::untar(GSE250130_RAW_path, exdir = GSE250130_extract_root)
        TRUE
    }, error = function(error_object) {
        stop("GSE250130 outer TAR extraction failed: ", conditionMessage(error_object))
    })
    if (isTRUE(extraction_status)) {
        write_text(c(paste0("Run ID: ", run_id), paste0("Source: ", normalize_slashes(GSE250130_RAW_path)), paste0("Completed: ", 
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"))), outer_marker)
    }
}

nested_bundles <- list.files(GSE250130_extract_root, pattern = "\\.(tar|tar\\.gz|tgz)$", recursive = TRUE, full.names = TRUE, 
    ignore.case = TRUE)

if (length(nested_bundles) != 28L) {
    stop("Expected 28 GSE250130 nested sample bundles but found ", length(nested_bundles), ".")
}

GSE250130_bundle_rows <- list()

bundle_row_index <- 1L

for (bundle_path in sort(nested_bundles)) {
    sample_id <- extract_GSM(bundle_path)
    if (is.na(sample_id)) {
        stop("No GSM accession was found in GSE250130 bundle name: ", basename(bundle_path))
    }
    sample_dir <- file.path(GSE250130_sample_root, sample_id)
    marker_path <- file.path(sample_dir, "_EXTRACTION_COMPLETE.txt")
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(marker_path)) {
        log_event("Extracting", sample_id)
        tryCatch({
            utils::untar(bundle_path, exdir = sample_dir)
        }, error = function(error_object) {
            stop("Nested extraction failed for ", sample_id, ": ", conditionMessage(error_object))
        })
        write_text(c(paste0("Run ID: ", run_id), paste0("Bundle: ", normalize_slashes(bundle_path))), marker_path)
    }
    extracted_files <- list.files(sample_dir, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
    components <- vapply(extracted_files, identify_triplet_component, character(1))
    matrix_files <- extracted_files[components == "matrix"]
    barcode_files <- extracted_files[components == "barcodes"]
    feature_files <- extracted_files[components == "features_or_genes"]
    GSE250130_bundle_rows[[bundle_row_index]] <- data.frame(dataset_id = "GSE250130", sample_accession = sample_id, 
        nested_bundle_path = normalize_slashes(bundle_path), matrix_file = if (length(matrix_files) == 1L) {
            normalize_slashes(matrix_files)
        }
        else {
            paste(normalize_slashes(matrix_files), collapse = ";")
        }, barcode_file = if (length(barcode_files) == 1L) {
            normalize_slashes(barcode_files)
        }
        else {
            paste(normalize_slashes(barcode_files), collapse = ";")
        }, feature_file = if (length(feature_files) == 1L) {
            normalize_slashes(feature_files)
        }
        else {
            paste(normalize_slashes(feature_files), collapse = ";")
        }, matrix_file_n = length(matrix_files), barcode_file_n = length(barcode_files), feature_file_n = length(feature_files), 
        complete_triplet = length(matrix_files) == 1L && length(barcode_files) == 1L && length(feature_files) == 
            1L, stringsAsFactors = FALSE)
    bundle_row_index <- bundle_row_index + 1L
}

GSE250130_bundle_map <- do.call(rbind, GSE250130_bundle_rows)

if (nrow(GSE250130_bundle_map) != 28L || anyDuplicated(GSE250130_bundle_map$sample_accession) > 0L || !all(GSE250130_bundle_map$complete_triplet)) {
    stop("GSE250130 did not yield 28 unique complete matrix triplets.")
}

write_csv_utf8(GSE250130_bundle_map, file.path(dirs$results, "STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv"))

GSE183837_triplet_paths <- source_candidate$path[source_candidate$dataset_id == "GSE183837" & source_candidate$source_role == 
    "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER"]

if (length(GSE183837_triplet_paths) != 27L) {
    stop("Expected 27 locked GSE183837 triplet members.")
}

GSE183837_triplet_inventory <- data.frame(dataset_id = "GSE183837", sample_accession = vapply(GSE183837_triplet_paths, 
    extract_GSM, character(1)), component = vapply(GSE183837_triplet_paths, identify_triplet_component, character(1)), 
    path = normalize_slashes(GSE183837_triplet_paths), stringsAsFactors = FALSE)

if (any(is.na(GSE183837_triplet_inventory$sample_accession)) || any(GSE183837_triplet_inventory$component == "other")) {
    stop("GSE183837 triplet filename classification failed.")
}

GSE183837_bundle_summary <- aggregate(component ~ sample_accession, data = GSE183837_triplet_inventory, FUN = function(x) {
    paste(sort(x), collapse = ";")
})

GSE183837_component_count <- aggregate(component ~ sample_accession, data = GSE183837_triplet_inventory, FUN = length)

names(GSE183837_component_count)[names(GSE183837_component_count) == "component"] <- "component_n"

GSE183837_bundle_summary <- merge(GSE183837_bundle_summary, GSE183837_component_count, by = "sample_accession", 
    all = TRUE, sort = FALSE)

GSE183837_bundle_summary$complete_triplet <- GSE183837_bundle_summary$component_n == 3L & GSE183837_bundle_summary$component == 
    paste(sort(c("barcodes", "features_or_genes", "matrix")), collapse = ";")

if (nrow(GSE183837_bundle_summary) != 9L || !all(GSE183837_bundle_summary$complete_triplet)) {
    stop("GSE183837 did not yield nine complete sample triplets.")
}

write_csv_utf8(GSE183837_triplet_inventory, file.path(dirs$qc, "STEP04B_GSE183837_TRIPLET_INVENTORY.csv"))

write_csv_utf8(GSE183837_bundle_summary, file.path(dirs$results, "STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv"))

series_matrix_rows <- source_candidate[grepl("SERIES_MATRIX", source_candidate$source_role), , drop = FALSE]

dataset_ids <- c("GSE250130", "GSE165004", "GSE92324", "GSE111974", "GSE183837")

series_audit_rows <- list()

sample_sheet_list <- list()

characteristics_list <- list()

for (dataset_id in dataset_ids) {
    dataset_series_rows <- series_matrix_rows[series_matrix_rows$dataset_id == dataset_id, , drop = FALSE]
    if (nrow(dataset_series_rows) != 1L) {
        stop("Expected exactly one series matrix for ", dataset_id, " but found ", nrow(dataset_series_rows), ".")
    }
    series_path <- dataset_series_rows$path[[1]]
    parsed_series <- parse_GEO_series_matrix(series_path)
    sample_sheet <- metadata_to_sample_sheet(dataset_id, parsed_series)
    characteristic_long <- characteristics_to_key_value(dataset_id, sample_sheet)
    sample_sheet_list[[dataset_id]] <- sample_sheet
    characteristics_list[[dataset_id]] <- characteristic_long
    series_audit_rows[[dataset_id]] <- data.frame(dataset_id = dataset_id, series_matrix_path = normalize_slashes(series_path), 
        parsed_sample_n = parsed_series$sample_n, expression_feature_n = parsed_series$expression_feature_n, expression_header_sample_n = length(parsed_series$expression_sample_headers), 
        metadata_field_n = length(unique(parsed_series$metadata_long$field_name)), characteristic_key_n = length(unique(characteristic_long$parsed_key[!is.na(characteristic_long$parsed_key)])), 
        sample_design_status = "CANDIDATE_PENDING_MANUAL_LOCK", stringsAsFactors = FALSE)
    write_csv_utf8(parsed_series$metadata_long, file.path(dirs$qc, paste0("STEP04B_", dataset_id, "_SERIES_METADATA_LONG.csv")))
    write_csv_utf8(sample_sheet, file.path(design_root, paste0("STEP04B_", dataset_id, "_SAMPLE_SHEET_CANDIDATE.csv")))
    write_csv_utf8(characteristic_long, file.path(design_root, paste0("STEP04B_", dataset_id, "_CHARACTERISTICS_KEY_VALUE.csv")))
}

series_audit <- do.call(rbind, series_audit_rows)

write_csv_utf8(series_audit, file.path(dirs$results, "STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv"))

GSE250130_geo_samples <- sample_sheet_list[["GSE250130"]]$sample_accession

GSE183837_geo_samples <- sample_sheet_list[["GSE183837"]]$sample_accession

GSE250130_accession_pass <- setequal(GSE250130_geo_samples, GSE250130_bundle_map$sample_accession)

GSE183837_accession_pass <- setequal(GSE183837_geo_samples, GSE183837_bundle_summary$sample_accession)

if (!GSE250130_accession_pass || !GSE183837_accession_pass) {
    stop("Single-cell sample-bundle accessions do not reconcile with GEO metadata.")
}

sheet_165004 <- sample_sheet_list[["GSE165004"]]

sheet_111974 <- sample_sheet_list[["GSE111974"]]

title_column_165004 <- grep("^title__", names(sheet_165004), value = TRUE)

title_column_111974 <- grep("^title__", names(sheet_111974), value = TRUE)

normalize_title <- function(x) {
    tolower(gsub("[^a-z0-9]+", "", as.character(x), perl = TRUE))
}

shared_GSM <- intersect(sheet_165004$sample_accession, sheet_111974$sample_accession)

shared_title_rows <- data.frame(GSE165004_sample = character(0), GSE111974_sample = character(0), normalized_title = character(0), 
    stringsAsFactors = FALSE)

if (length(title_column_165004) >= 1L && length(title_column_111974) >= 1L) {
    titles_165004 <- data.frame(sample_accession = sheet_165004$sample_accession, normalized_title = normalize_title(sheet_165004[[title_column_165004[[1]]]]), 
        stringsAsFactors = FALSE)
    titles_111974 <- data.frame(sample_accession = sheet_111974$sample_accession, normalized_title = normalize_title(sheet_111974[[title_column_111974[[1]]]]), 
        stringsAsFactors = FALSE)
    shared_title_rows <- merge(titles_165004, titles_111974, by = "normalized_title", suffixes = c("_GSE165004", 
        "_GSE111974"), all = FALSE)
    names(shared_title_rows) <- c("normalized_title", "GSE165004_sample", "GSE111974_sample")
    shared_title_rows <- shared_title_rows[nzchar(shared_title_rows$normalized_title), , drop = FALSE]
}

shared_control_audit <- data.frame(metric = c("Exact shared GSM accessions", "Normalized shared sample titles", 
    "GSE165004 parsed samples", "GSE111974 parsed samples", "Shared-control design status"), value = c(length(shared_GSM), 
    nrow(shared_title_rows), nrow(sheet_165004), nrow(sheet_111974), "PENDING_MANUAL_CONTROL_SHARING_LOCK"), stringsAsFactors = FALSE)

write_csv_utf8(shared_control_audit, file.path(dirs$results, "STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv"))

write_csv_utf8(shared_title_rows, file.path(dirs$qc, "STEP04B_GSE165004_GSE111974_SHARED_TITLE_ROWS.csv"))

write_text(shared_GSM, file.path(dirs$qc, "STEP04B_GSE165004_GSE111974_SHARED_GSM.txt"))

locked_source_reference_path <- file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_REFERENCE_LOCK.csv")

write_csv_utf8(source_candidate, locked_source_reference_path)

locked_role_path <- file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_DATASET_ROLE_LOCK.csv")

write_csv_utf8(role_lock, locked_role_path)

locked_GSE250130_bundle_path <- file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_GSE250130_BUNDLE_MAP.csv")

write_csv_utf8(GSE250130_bundle_map, locked_GSE250130_bundle_path)

locked_GSE183837_bundle_path <- file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_GSE183837_BUNDLE_MAP.csv")

write_csv_utf8(GSE183837_bundle_summary, locked_GSE183837_bundle_path)

locked_series_audit_path <- file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_SERIES_AUDIT.csv")

write_csv_utf8(series_audit, locked_series_audit_path)

locked_paths <- c(locked_source_reference_path, locked_role_path, locked_GSE250130_bundle_path, locked_GSE183837_bundle_path, 
    locked_series_audit_path)

locked_manifest <- data.frame(artifact_id = c("COMBINED_SOURCE_REFERENCE_LOCK", "DATASET_ROLE_LOCK", "GSE250130_SAMPLE_BUNDLE_MAP", 
    "GSE183837_SAMPLE_BUNDLE_MAP", "SERIES_MATRIX_SAMPLE_AUDIT"), lock_version = source_lock_version, path = normalize_slashes(locked_paths), 
    size_bytes = as.numeric(file.info(locked_paths)$size), md5 = vapply(locked_paths, md5_file, character(1)), 
    sha256 = vapply(locked_paths, sha256_file, character(1)), lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(locked_manifest, file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv"))

metadata_field_inventory <- do.call(rbind, lapply(dataset_ids, function(dataset_id) {
    characteristic_long <- characteristics_list[[dataset_id]]
    if (nrow(characteristic_long) == 0L) {
        return(data.frame(dataset_id = dataset_id, parsed_key = NA_character_, nonmissing_value_n = 0L, unique_value_n = 0L, 
            example_values = NA_character_, stringsAsFactors = FALSE))
    }
    split_rows <- split(characteristic_long, characteristic_long$parsed_key, drop = TRUE)
    do.call(rbind, lapply(names(split_rows), function(key_name) {
        key_rows <- split_rows[[key_name]]
        values <- unique(key_rows$parsed_value[!is.na(key_rows$parsed_value) & nzchar(trimws(key_rows$parsed_value))])
        data.frame(dataset_id = dataset_id, parsed_key = key_name, nonmissing_value_n = sum(!is.na(key_rows$parsed_value) & 
            nzchar(trimws(key_rows$parsed_value))), unique_value_n = length(values), example_values = paste(head(values, 
            8L), collapse = ";"), stringsAsFactors = FALSE)
    }))
}))

write_csv_utf8(metadata_field_inventory, file.path(dirs$qc, "STEP04B_METADATA_FIELD_INVENTORY.csv"))

gate_status <- "PROCEED_TO_MANUAL_SAMPLE_DESIGN_LOCK_AND_COMBINED_ANALYSIS_BUILD"

gate_table <- data.frame(gate_component = c("STEP04A v2 gate", "Combined source rows", "Source identity recheck", 
    "GSE250130 outer extraction", "GSE250130 complete triplets", "GSE250130 GEO/bundle accession reconciliation", 
    "GSE183837 complete triplets", "GSE183837 GEO/bundle accession reconciliation", "Five series matrices parsed", 
    "GSE165004/GSE111974 shared-control audit", "Combined source lock", "Sample-design candidates", "Biological analysis", 
    "Scientific figures", "STEP04B gate"), status = c("PASS", "LOCKED_42", "PASS", "PASS_28_NESTED_BUNDLES", "PASS_28_OF_28", 
    "PASS", "PASS_9_OF_9", "PASS", "PASS_5_OF_5", paste0("EXACT_GSM_", length(shared_GSM), "_NORMALIZED_TITLE_", 
        nrow(shared_title_rows)), paste0("LOCKED_", source_lock_version), paste0("CREATED_", sample_design_candidate_version), 
    "NOT_STARTED", "NONE_CREATED", gate_status), required_action = c("Retain the acquisition audit.", "Do not change source paths or hashes without a new version.", 
    "Use the locked source-reference table.", "Do not re-extract unless the extraction markers are removed intentionally.", 
    "Use donor-aware single-cell analysis.", "Use the 28-sample GEO design.", "Use nine donor-level sample bundles.", 
    "Use the nine-sample GEO design.", "Review group, timepoint, fertility/RIF and tissue fields.", "Lock shared-control handling before external validation.", 
    "Do not duplicate the 3.6-GB TAR inside results or archives.", "Manually approve sample/group assignments before analysis.", 
    "Build one coordinated biological-analysis package.", "Render only after biological QC.", "Upload the sample audits, metadata inventory, shared-control audit and Gate."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP04B_GATE_TABLE.csv"))

write_text(c("# STEP04B source-lock and sample-design audit card", "", paste0("- Run ID: ", run_id), paste0("- Combined input source: ", 
    source_lock_version, " LOCKED"), "- Source references: 42", "- GSE250130 complete sample triplets: 28/28", 
    "- GSE183837 complete sample triplets: 9/9", paste0("- GSE165004/GSE111974 exact shared GSM: ", length(shared_GSM)), 
    paste0("- GSE165004/GSE111974 normalized shared titles: ", nrow(shared_title_rows)), "- Five series matrices parsed: YES", 
    "- Sample-design candidates: CREATED", "- Biological analysis: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP04B_SOURCE_LOCK_SAMPLE_DESIGN_AUDIT_CARD.md"))

write_text(c("# STEP04B Methods note", "", "The five-dataset source set was locked by reference rather than by duplicating large raw files. All 42 source rows were revalidated by path and file size, with MD5 and SHA-256 recomputed for files no larger than 500 MB. The 3.6-GB GSE250130 outer TAR retained its acquisition hashes and was additionally validated by archive structure and complete extraction.", 
    "", "The 28 GSE250130 nested bundles and nine GSE183837 sample triplets were reconciled with GEO sample accessions. GEO sample metadata, titles and characteristics were parsed into dataset-specific sample-design candidates. No biological association analysis or scientific figure was generated."), 
    file.path(dirs$text, "STEP04B_Methods_note.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, source_lock_version = source_lock_version, source_row_n = nrow(source_candidate), 
    GSE250130_complete_triplet_n = nrow(GSE250130_bundle_map), GSE183837_complete_triplet_n = nrow(GSE183837_bundle_summary), 
    parsed_series_matrix_n = nrow(series_audit), exact_shared_GSM_n = length(shared_GSM), normalized_shared_title_n = nrow(shared_title_rows), 
    biological_analysis_started = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04B_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04B.txt"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    copy_if_different(source_file, file.path(dirs$code, "extract_patient_inputs.R"))
}

archive_artifacts <- c(file.path(dirs$results, "STEP04B_GATE_TABLE.csv"), file.path(dirs$results, "STEP04B_SOURCE_LOCK_SAMPLE_DESIGN_AUDIT_CARD.md"), 
    file.path(dirs$results, "STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv"), file.path(dirs$results, "STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv"), 
    file.path(dirs$results, "STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv"), file.path(dirs$results, "STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv"), 
    file.path(dirs$qc, "STEP04B_SOURCE_IDENTITY_RECHECK.csv"), file.path(dirs$qc, "STEP04B_METADATA_FIELD_INVENTORY.csv"), 
    file.path(locked_root, "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv"), locked_source_reference_path, locked_role_path, 
    locked_GSE250130_bundle_path, locked_GSE183837_bundle_path, locked_series_audit_path, file.path(dirs$text, 
        "STEP04B_Methods_note.md"), file.path(dirs$parameters, "STEP04B_parameters.json"), file.path(dirs$environment, 
        "sessionInfo_POST_RUN_STEP04B.txt"), file.path(dirs$code, "extract_patient_inputs.R"), log_path)

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
    else if (grepl("Methods", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|AUDIT|CARD|MANIFEST|LOCK|MAP|SUMMARY|INVENTORY", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04B_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04B archive status", "", paste0("- Run ID: ", run_id), paste0("- Source lock: ", source_lock_version), 
    "- Source rows: 42", "- GSE250130 complete triplets: 28", "- GSE183837 complete triplets: 9", "- Raw files duplicated into archive: NO", 
    "- Biological analysis: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(archive_root, 
    "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Accepted run: ", run_id), paste0("- Combined input source: ", 
    source_lock_version, " LOCKED"), "- GSE250130: 28 complete sample bundles", "- GSE183837: 9 complete sample bundles", 
    "- Five series matrices: PARSED", "- Sample designs: CANDIDATE_PENDING_MANUAL_LOCK", "- Biological analysis: NOT STARTED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04B combined input lock and sample-design audit completed.\n")

cat("Combined input source: ", source_lock_version, " LOCKED\n")

cat("Source rows: 42\n")

cat("GSE250130 complete sample triplets: 28/28\n")

cat("GSE183837 complete sample triplets: 9/9\n")

cat("GSE165004/GSE111974 exact shared GSM: ", length(shared_GSM), "\n")

cat("GSE165004/GSE111974 normalized shared titles: ", nrow(shared_title_rows), "\n")

cat("Five series matrices parsed: YES\n")

cat("Sample-design candidates: CREATED\n")

cat("Biological analysis: NOT STARTED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP04B_GATE_TABLE.csv\n")

cat("  06_results/STEP04B_SERIES_MATRIX_SAMPLE_AUDIT.csv\n")

cat("  06_results/STEP04B_GSE250130_SAMPLE_BUNDLE_MAP.csv\n")

cat("  06_results/STEP04B_GSE183837_SAMPLE_BUNDLE_SUMMARY.csv\n")

cat("  06_results/STEP04B_GSE165004_GSE111974_SHARED_CONTROL_AUDIT.csv\n")

cat("  08_qc/STEP04B_METADATA_FIELD_INVENTORY.csv\n")

cat("  06_results/STEP04B_sample_design_candidates/\n")

cat("============================================================\n")

