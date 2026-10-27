#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B01: Acquire patient-cohort sources

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(3600L, getOption("timeout")))

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04A_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

download_large_GSE250130_RAW <- TRUE

GSE250130_series_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250130/matrix/GSE250130_series_matrix.txt.gz"

GSE250130_RAW_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE250nnn/GSE250130/suppl/GSE250130_RAW.tar"

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

archive_root <- file.path(dirs$archive_parent, paste0("S4Av2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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
        stop("Destination path exceeds 240 characters:\n", destination_path)
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

download_with_resume <- function(url, destination, large_file = FALSE) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    curl_candidates <- c(Sys.which("curl.exe"), Sys.which("curl"))
    curl_candidates <- curl_candidates[nzchar(curl_candidates)]
    if (length(curl_candidates) > 0L) {
        curl_binary <- curl_candidates[[1]]
        args <- c("--fail", "--location", "--retry", "20", "--retry-all-errors", "--retry-delay", "10", "--connect-timeout", 
            "60", "--continue-at", "-", "--output", shQuote(destination), shQuote(url))
        log_event("Downloading with resume support:", basename(destination))
        status <- system2(command = curl_binary, args = args, stdout = "", stderr = "")
        if (status != 0L || !file.exists(destination)) {
            stop("curl download failed for: ", url)
        }
        return(TRUE)
    }
    if (large_file && file.exists(destination)) {
        stop("A partial large file exists but curl is unavailable for resume: ", destination)
    }
    log_event("Downloading with R download.file:", basename(destination))
    download.file(url = url, destfile = destination, mode = "wb", method = "libcurl", quiet = FALSE)
    if (!file.exists(destination)) {
        stop("download.file failed for: ", url)
    }
    TRUE
}

parse_GEO_series_metadata <- function(path) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    all_lines <- readLines(connection, warn = FALSE, encoding = "UTF-8")
    sample_lines <- all_lines[grepl("^!Sample_", all_lines)]
    if (length(sample_lines) == 0L) {
        return(data.frame())
    }
    parsed_rows <- lapply(sample_lines, function(line_value) {
        fields <- strsplit(line_value, "\t", fixed = TRUE)[[1]]
        field_name <- sub("^!Sample_", "", fields[[1]])
        values <- fields[-1]
        values <- gsub("^\"|\"$", "", values)
        data.frame(field_name = field_name, sample_index = seq_along(values), value = values, stringsAsFactors = FALSE)
    })
    do.call(rbind, parsed_rows)
}

build_source_row <- function(dataset_id, source_role, path, source_status = "CANDIDATE_PENDING_MANUAL_LOCK") {
    data.frame(dataset_id = dataset_id, source_role = source_role, path = normalize_slashes(path), file_name = basename(path), 
        file_exists = file.exists(path), size_bytes = if (file.exists(path)) {
            as.numeric(file.info(path)$size)
        }
        else {
            NA_real_
        }, md5 = if (file.exists(path)) {
            md5_file(path)
        }
        else {
            NA_character_
        }, sha256 = if (file.exists(path)) {
            sha256_file(path)
        }
        else {
            NA_character_
        }, source_status = source_status, stringsAsFactors = FALSE)
}

log_event("Starting", run_id)

v1_gate_path <- file.path(dirs$results, "STEP04A_GATE_TABLE.csv")

v1_readiness_path <- file.path(dirs$results, "STEP04A_INPUT_READINESS_MATRIX.csv")

if (!file.exists(v1_gate_path) || !file.exists(v1_readiness_path)) {
    stop("STEP04A v1 gate or readiness matrix is missing.")
}

v1_gate <- .vdr_read_csv(v1_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

v1_readiness <- .vdr_read_csv(v1_readiness_path, stringsAsFactors = FALSE, check.names = FALSE)

v1_gate_value <- v1_gate$status[v1_gate$gate_component == "STEP04A gate"]

if (length(v1_gate_value) != 1L || v1_gate_value != "HOLD_MISSING_REQUIRED_DOWNSTREAM_DATASET_FILES") {
    stop("Unexpected STEP04A v1 gate: ", paste(v1_gate_value, collapse = ", "))
}

if (v1_readiness$preliminary_readiness[v1_readiness$dataset_id == "GSE250130"] != "MISSING") {
    stop("STEP04A v2 correction is not applicable because GSE250130 was not marked missing.")
}

dataset_role_lock_v2 <- data.frame(dataset_id = c("GSE250130", "GSE165004", "GSE92324", "GSE111974", "GSE183837"), 
    analysis_role = c("TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT", "PRIMARY_PATIENT_PROJECTION", "EXTERNAL_PATIENT_SUPPORT", 
        "EXTERNAL_PATIENT_SUPPORT_WITH_SHARED_CONTROL_CONSTRAINT", "INDEPENDENT_SINGLE_CELL_LOCALIZATION"), priority = c(1L, 
        1L, 2L, 2L, 1L), required_for_main_manuscript = TRUE, primary_input_expected = c("28 sample-level single-cell count bundles plus LH-timepoint, fertility/RIF and donor metadata", 
        "processed endometrial expression matrix plus phenotype and slide metadata", "processed endometrial expression matrix plus phenotype metadata", 
        "processed endometrial expression matrix plus phenotype metadata and explicit shared-control audit", "9 sample-level single-cell count triplets plus Ctrl/RIF and donor metadata"), 
    locked_constraints = c("Use fertile LH+3/+5/+7/+9/+11 samples for temporal support; RIF is available at LH+7. Use donor-aware cell-type or pseudobulk analyses and do not pool cells as independent replicates.", 
        "Use the accepted slide-aware and leave-one-slide robustness framework.", "Use as external supportive validation; do not train a diagnostic classifier.", 
        "Account explicitly for controls shared with GSE165004.", "Use for independent cell-type/state localization; do not infer causality from localization."), 
    version = "STEP04A_DATASET_ROLE_LOCK_v2.0", stringsAsFactors = FALSE)

write_csv_utf8(dataset_role_lock_v2, file.path(dirs$results, "STEP04A_v2_DATASET_ROLE_LOCK.csv"))

GSE250130_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE250130")

GSE250130_matrix_dir <- file.path(GSE250130_root, "series_matrix")

GSE250130_raw_dir <- file.path(GSE250130_root, "raw_counts")

dir.create(GSE250130_matrix_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(GSE250130_raw_dir, recursive = TRUE, showWarnings = FALSE)

GSE250130_series_matrix_path <- file.path(GSE250130_matrix_dir, "GSE250130_series_matrix.txt.gz")

GSE250130_RAW_path <- file.path(GSE250130_raw_dir, "GSE250130_RAW.tar")

if (!file.exists(GSE250130_series_matrix_path)) {
    download_with_resume(GSE250130_series_matrix_url, GSE250130_series_matrix_path, large_file = FALSE)
}

if (download_large_GSE250130_RAW && !file.exists(GSE250130_RAW_path)) {
    download_with_resume(GSE250130_RAW_url, GSE250130_RAW_path, large_file = TRUE)
}

if (download_large_GSE250130_RAW && file.exists(GSE250130_RAW_path) && as.numeric(file.info(GSE250130_RAW_path)$size) < 
    2500000000) {
    stop("GSE250130_RAW.tar is smaller than 2.5 GB and is likely incomplete. Re-run the script to resume.")
}

GSE250130_series_metadata_long <- parse_GEO_series_metadata(GSE250130_series_matrix_path)

if (nrow(GSE250130_series_metadata_long) == 0L) {
    stop("GSE250130 series metadata could not be parsed.")
}

write_csv_utf8(GSE250130_series_metadata_long, file.path(dirs$qc, "STEP04A_v2_GSE250130_SERIES_METADATA_LONG.csv"))

GSE250130_sample_accessions <- GSE250130_series_metadata_long$value[GSE250130_series_metadata_long$field_name == 
    "geo_accession"]

GSE250130_sample_accessions <- unique(GSE250130_sample_accessions[nzchar(GSE250130_sample_accessions)])

if (length(GSE250130_sample_accessions) != 28L) {
    stop("Expected 28 GSE250130 samples but parsed ", length(GSE250130_sample_accessions), ".")
}

tar_member_list <- character(0)

tar_inventory_status <- "RAW_TAR_NOT_PRESENT"

if (file.exists(GSE250130_RAW_path)) {
    tar_member_list <- tryCatch(utils::untar(GSE250130_RAW_path, list = TRUE), error = function(error_object) {
        character(0)
    })
    if (length(tar_member_list) > 0L) {
        tar_inventory_status <- "RAW_TAR_LIST_OK"
    }
    else {
        tar_inventory_status <- "RAW_TAR_PRESENT_BUT_LIST_FAILED"
    }
}

write_text(tar_member_list, file.path(dirs$qc, "STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt"))

tar_summary <- data.frame(metric = c("RAW tar exists", "RAW tar size bytes", "RAW tar member-list status", "RAW tar member count", 
    "Parsed GEO sample count", "Series matrix exists"), value = c(as.character(file.exists(GSE250130_RAW_path)), 
    if (file.exists(GSE250130_RAW_path)) {
        as.character(as.numeric(file.info(GSE250130_RAW_path)$size))
    } else {
        NA_character_
    }, tar_inventory_status, as.character(length(tar_member_list)), as.character(length(GSE250130_sample_accessions)), 
    as.character(file.exists(GSE250130_series_matrix_path))), stringsAsFactors = FALSE)

write_csv_utf8(tar_summary, file.path(dirs$qc, "STEP04A_v2_GSE250130_ACQUISITION_QC.csv"))

step00_locked_stat_root <- file.path(project_root, "06_locked_results", "statistical_outputs", "STEP00")

slide_QC_paths <- file.path(step00_locked_stat_root, c("GSE165004_leave_one_slide_out_metrics.csv", "GSE165004_slide_aggregate_top1000_per_contrast.csv", 
    "GSE165004_slide_by_group_counts.csv", "GSE165004_slide_confounding_summary.csv"))

missing_slide_QC <- slide_QC_paths[!file.exists(slide_QC_paths)]

if (length(missing_slide_QC) > 0L) {
    stop("Accepted GSE165004 slide-QC file is missing:\n", paste(missing_slide_QC, collapse = "\n"))
}

slide_QC_lock <- do.call(rbind, lapply(slide_QC_paths, function(path_value) {
    build_source_row(dataset_id = "GSE165004", source_role = "LOCKED_SLIDE_ROBUSTNESS_STATISTICAL_OUTPUT", path = path_value, 
        source_status = "LOCKED_ACCEPTED_STEP00_SOURCE")
}))

write_csv_utf8(slide_QC_lock, file.path(dirs$results, "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv"))

source_rows <- list()

add_source <- function(dataset_id, source_role, path, status = "CANDIDATE_PENDING_MANUAL_LOCK") {
    source_rows[[length(source_rows) + 1L]] <<- build_source_row(dataset_id = dataset_id, source_role = source_role, 
        path = path, source_status = status)
}

add_source("GSE250130", "SERIES_MATRIX_METADATA", GSE250130_series_matrix_path)

add_source("GSE250130", "RAW_OUTER_TAR_OF_SAMPLE_BUNDLES", GSE250130_RAW_path)

add_source("GSE165004", "SERIES_MATRIX_EXPRESSION", file.path(project_root, "03_data", "raw_external_READ_ONLY", 
    "GSE165004", "series_matrix", "GSE165004_series_matrix.txt.gz"))

add_source("GSE165004", "LOCKED_SERIES_METADATA", file.path(step00_locked_stat_root, "GSE165004_series_metadata_parsed.csv"))

add_source("GSE165004", "LOCKED_GROUP_METADATA", file.path(step00_locked_stat_root, "GSE165004_title_group_audit.csv"))

for (path_value in slide_QC_paths) {
    add_source("GSE165004", "LOCKED_SLIDE_QC", path_value, status = "LOCKED_ACCEPTED_STEP00_SOURCE")
}

add_source("GSE92324", "SERIES_MATRIX_EXPRESSION_AND_METADATA", file.path(project_root, "03_data", "raw_external_READ_ONLY", 
    "GSE92324", "series_matrix", "GSE92324_series_matrix.txt.gz"))

add_source("GSE92324", "PLATFORM_ANNOTATION", file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE92324", 
    "platform_annotation", "GPL10558.soft.gz"))

add_source("GSE111974", "SERIES_MATRIX_EXPRESSION", file.path(project_root, "03_data", "raw_external_READ_ONLY", 
    "GSE111974", "series_matrix", "GSE111974_series_matrix.txt.gz"))

add_source("GSE111974", "LOCKED_SERIES_METADATA", file.path(step00_locked_stat_root, "GSE111974_series_metadata_parsed.csv"))

add_source("GSE111974", "LOCKED_GROUP_METADATA", file.path(step00_locked_stat_root, "GSE111974_title_group_audit.csv"))

GSE183837_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE183837")

GSE183837_triplet_files <- list.files(file.path(GSE183837_root, "raw_counts", "matrix_market_by_donor"), pattern = "\\.(barcodes|genes|matrix)\\.(tsv|mtx)\\.gz$", 
    full.names = TRUE, ignore.case = TRUE)

if (length(GSE183837_triplet_files) != 27L) {
    stop("Expected 27 GSE183837 triplet files but found ", length(GSE183837_triplet_files), ".")
}

for (path_value in sort(GSE183837_triplet_files)) {
    add_source("GSE183837", "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER", path_value)
}

add_source("GSE183837", "SERIES_MATRIX_METADATA", file.path(GSE183837_root, "series_matrix", "GSE183837_series_matrix.txt.gz"))

source_lock_candidate <- do.call(rbind, source_rows)

missing_source_rows <- source_lock_candidate[!source_lock_candidate$file_exists, , drop = FALSE]

write_csv_utf8(source_lock_candidate, file.path(dirs$results, "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv"))

write_csv_utf8(missing_source_rows, file.path(dirs$qc, "STEP04A_v2_MISSING_SOURCE_ROWS.csv"))

source_summary <- aggregate(file_exists ~ dataset_id, data = source_lock_candidate, FUN = function(x) {
    paste0(sum(x), "/", length(x))
})

names(source_summary)[names(source_summary) == "file_exists"] <- "files_present"

source_summary$all_required_files_present <- vapply(source_summary$dataset_id, function(dataset_id) {
    all(source_lock_candidate$file_exists[source_lock_candidate$dataset_id == dataset_id])
}, logical(1))

write_csv_utf8(source_summary, file.path(dirs$results, "STEP04A_v2_DATASET_SOURCE_READINESS.csv"))

all_required_present <- all(source_lock_candidate$file_exists)

GSE250130_RAW_ready <- file.exists(GSE250130_RAW_path) && as.numeric(file.info(GSE250130_RAW_path)$size) >= 2500000000 && 
    tar_inventory_status == "RAW_TAR_LIST_OK"

gate_status <- if (!GSE250130_RAW_ready) {
    "HOLD_GSE250130_RAW_DOWNLOAD_OR_TAR_INVENTORY_INCOMPLETE"
} else if (!all_required_present) {
    "HOLD_OTHER_REQUIRED_SOURCE_FILES_MISSING"
} else {
    "PROCEED_TO_MANUAL_STEP04_INPUT_SOURCE_LOCK"
}

gate_table <- data.frame(gate_component = c("STEP04A v1 hold", "GSE250130 scientific role", "GSE250130 series metadata", 
    "GSE250130 sample count", "GSE250130 RAW tar", "GSE250130 tar inventory", "GSE165004 source candidate", "GSE165004 slide-QC source lock", 
    "GSE92324 source candidate", "GSE111974 source candidate", "GSE183837 triplet source candidate", "All five dataset source rows", 
    "Scientific figures", "STEP04A v2 gate"), status = c("CONFIRMED_ONLY_GSE250130_MISSING", "LOCKED_TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT", 
    ifelse(file.exists(GSE250130_series_matrix_path), "PASS", "MISSING"), paste0("PARSED_", length(GSE250130_sample_accessions)), 
    ifelse(file.exists(GSE250130_RAW_path), paste0("PRESENT_BYTES_", as.numeric(file.info(GSE250130_RAW_path)$size)), 
        "MISSING"), tar_inventory_status, ifelse(all(source_lock_candidate$file_exists[source_lock_candidate$dataset_id == 
        "GSE165004"]), "PASS", "MISSING_SOURCE"), "LOCKED_4_CANONICAL_FILES", ifelse(all(source_lock_candidate$file_exists[source_lock_candidate$dataset_id == 
        "GSE92324"]), "PASS", "MISSING_SOURCE"), ifelse(all(source_lock_candidate$file_exists[source_lock_candidate$dataset_id == 
        "GSE111974"]), "PASS", "MISSING_SOURCE"), ifelse(all(source_lock_candidate$file_exists[source_lock_candidate$dataset_id == 
        "GSE183837"]), "PASS_27_TRIPLET_FILES_PLUS_METADATA", "MISSING_SOURCE"), ifelse(all_required_present, "PASS", 
        paste0("MISSING_ROWS_", nrow(missing_source_rows))), "NONE_CREATED", gate_status), required_action = c("Retain v1 as the discovery audit.", 
    "Use donor-aware temporal single-cell analysis.", "Use for sample metadata only.", "Require exactly 28 GEO samples.", 
    "Resume download until the full outer TAR is present.", "Require a readable TAR member list before source lock.", 
    "Use series matrix plus locked STEP00 metadata and slide evidence.", "Ignore the overinclusive 407-row v1 inventory.", 
    "Use series matrix and GPL10558 annotation.", "Use series matrix and locked metadata with shared-control audit.", 
    "Use nine complete donor triplets.", "Lock only after manual review of this candidate.", "No figure workflow is triggered.", 
    "Upload the source candidate, readiness, GSE250130 QC, narrow slide lock and gate."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP04A_v2_GATE_TABLE.csv"))

write_text(c("# STEP04A v2 acquisition and source-lock candidate card", "", paste0("- Run ID: ", run_id), "- GSE250130 role: TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT", 
    "- GSE250130 design: 28 biopsies across LH+3/+5/+7/+9/+11; RIF at LH+7", paste0("- GSE250130 GEO samples parsed: ", 
        length(GSE250130_sample_accessions)), paste0("- GSE250130 RAW tar present: ", file.exists(GSE250130_RAW_path)), 
    paste0("- GSE250130 RAW tar bytes: ", if (file.exists(GSE250130_RAW_path)) {
        as.numeric(file.info(GSE250130_RAW_path)$size)
    } else {
        NA
    }), paste0("- GSE250130 tar member count: ", length(tar_member_list)), "- GSE165004 slide-QC provenance: 4 canonical locked files", 
    paste0("- Combined source rows: ", nrow(source_lock_candidate)), paste0("- Missing source rows: ", nrow(missing_source_rows)), 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP04A_v2_ACQUISITION_SOURCE_CANDIDATE_CARD.md"))

write_text(c("# STEP04A v2 correction note", "", "STEP04A v1 correctly identified GSE250130 as the only missing mandatory dataset. However, its role was described too generically. GSE250130 is a time-series single-cell endometrial atlas and must be analyzed with donor-aware cell-type or pseudobulk methods.", 
    "", "The v1 GSE165004 slide-QC inventory was also overinclusive because every file under the STEP00 path was collected. The accepted provenance is narrowed to four canonical statistical outputs in the locked STEP00 directory.", 
    "", "No downstream biological analysis or scientific figure is generated in this acquisition step."), file.path(dirs$text, 
    "STEP04A_v2_Correction_note.md"))

parameters <- list(schema_version = "2.0", run_id = run_id, GSE250130_role = "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT", 
    GSE250130_expected_sample_n = 28, GSE250130_parsed_sample_n = length(GSE250130_sample_accessions), GSE250130_RAW_download_enabled = download_large_GSE250130_RAW, 
    GSE250130_RAW_present = file.exists(GSE250130_RAW_path), GSE250130_RAW_size_bytes = if (file.exists(GSE250130_RAW_path)) {
        as.numeric(file.info(GSE250130_RAW_path)$size)
    } else {
        NA
    }, GSE250130_tar_member_n = length(tar_member_list), GSE165004_slide_QC_locked_file_n = nrow(slide_QC_lock), 
    combined_source_row_n = nrow(source_lock_candidate), missing_source_row_n = nrow(missing_source_rows), biological_analysis_started = FALSE, 
    scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04A_v2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04A_v2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("48_step04A_v2_acquire_GSE250130_and_build_source_lock_candidate.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "acquire_patient_data.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP04A_v2_DATASET_ROLE_LOCK.csv"), file.path(dirs$results, "STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv"), 
    file.path(dirs$results, "STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv"), file.path(dirs$results, "STEP04A_v2_DATASET_SOURCE_READINESS.csv"), 
    file.path(dirs$results, "STEP04A_v2_GATE_TABLE.csv"), file.path(dirs$results, "STEP04A_v2_ACQUISITION_SOURCE_CANDIDATE_CARD.md"), 
    file.path(dirs$qc, "STEP04A_v2_GSE250130_SERIES_METADATA_LONG.csv"), file.path(dirs$qc, "STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt"), 
    file.path(dirs$qc, "STEP04A_v2_GSE250130_ACQUISITION_QC.csv"), file.path(dirs$qc, "STEP04A_v2_MISSING_SOURCE_ROWS.csv"), 
    file.path(dirs$text, "STEP04A_v2_Correction_note.md"), file.path(dirs$parameters, "STEP04A_v2_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04A_v2.txt"), file.path(dirs$code, "acquire_patient_data.R"), 
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
    else if (grepl("Correction", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|READINESS|LOCK|CANDIDATE|CARD|MEMBER|MISSING", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04A_v2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04A v2 archive status", "", paste0("- Run ID: ", run_id), "- GSE250130 role corrected: YES", 
    paste0("- GSE250130 RAW present: ", file.exists(GSE250130_RAW_path)), paste0("- GSE250130 GEO samples: ", length(GSE250130_sample_accessions)), 
    "- GSE165004 slide-QC canonical files: 4", paste0("- Missing source rows: ", nrow(missing_source_rows)), "- Biological analysis: NOT STARTED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Active run: ", run_id), "- Dataset role lock: v2 CANDIDATE", 
    "- GSE250130: ACQUISITION_AND_TAR_INVENTORY", "- GSE165004 slide-QC provenance: NARROWED_TO_4_CANONICAL_FILES", 
    "- Five-dataset source lock: CANDIDATE_PENDING_MANUAL_REVIEW", "- Biological analysis: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04A v2 GSE250130 acquisition and source-candidate build completed.\n")

cat("GSE250130 role: TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT\n")

cat("GSE250130 GEO samples parsed: ", length(GSE250130_sample_accessions), "\n")

cat("GSE250130 RAW tar present: ", file.exists(GSE250130_RAW_path), "\n")

cat("GSE250130 RAW tar bytes: ", if (file.exists(GSE250130_RAW_path)) {
    as.numeric(file.info(GSE250130_RAW_path)$size)
} else {
    NA
}, "\n")

cat("GSE250130 tar member count: ", length(tar_member_list), "\n")

cat("GSE165004 slide-QC canonical files: 4\n")

cat("Combined source rows: ", nrow(source_lock_candidate), "\n")

cat("Missing source rows: ", nrow(missing_source_rows), "\n")

cat("Biological analysis: NOT STARTED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP04A_v2_GATE_TABLE.csv\n")

cat("  06_results/STEP04A_v2_COMBINED_INPUT_SOURCE_LOCK_CANDIDATE.csv\n")

cat("  06_results/STEP04A_v2_DATASET_SOURCE_READINESS.csv\n")

cat("  06_results/STEP04A_v2_GSE165004_SLIDE_QC_SOURCE_LOCK.csv\n")

cat("  08_qc/STEP04A_v2_GSE250130_ACQUISITION_QC.csv\n")

cat("  08_qc/STEP04A_v2_GSE250130_RAW_TAR_MEMBER_LIST.txt\n")

cat("============================================================\n")

