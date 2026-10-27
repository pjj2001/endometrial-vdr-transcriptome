#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B03: Prepare accepted cohort designs and platform annotations

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(3600L, getOption("timeout")))

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04C_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

sample_design_lock_version <- "COMBINED_SAMPLE_DESIGN_v1.0"

platform_source_version <- "BULK_PLATFORM_ANNOTATION_SOURCE_v1.0"

GPL16699_self_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16699&targ=self&view=full&form=text"

GPL17077_self_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17077&targ=self&view=full&form=text"

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

locked_root <- file.path(dirs$results, "STEP04C_locked_sample_designs", sample_design_lock_version)

platform_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "platform_annotations")

dir.create(locked_root, recursive = TRUE, showWarnings = FALSE)

dir.create(platform_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S4C_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

try_download_with_retry <- function(url, destination) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(destination) && as.numeric(file.info(destination)$size) > 10000) {
        return(TRUE)
    }
    if (file.exists(destination)) {
        unlink(destination, force = TRUE)
    }
    curl_candidates <- c(Sys.which("curl.exe"), Sys.which("curl"))
    curl_candidates <- curl_candidates[nzchar(curl_candidates)]
    if (length(curl_candidates) > 0L) {
        status <- suppressWarnings(system2(curl_candidates[[1]], args = c("--fail", "--location", "--retry", "10", 
            "--retry-all-errors", "--retry-delay", "5", "--connect-timeout", "60", "--output", shQuote(destination), 
            shQuote(url))))
        return(identical(status, 0L) && file.exists(destination) && as.numeric(file.info(destination)$size) > 10000)
    }
    status <- tryCatch({
        download.file(url, destination, mode = "wb", method = "libcurl", quiet = FALSE)
        TRUE
    }, error = function(error_object) {
        FALSE
    })
    isTRUE(status) && file.exists(destination) && as.numeric(file.info(destination)$size) > 10000
}

read_platform_soft_header <- function(path) {
    if (!file.exists(path)) {
        return(data.frame(read_status = "MISSING", n_columns = NA_integer_, first_columns = NA_character_, stringsAsFactors = FALSE))
    }
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    table_begin_found <- FALSE
    header_line <- NA_character_
    line_counter <- 0L
    repeat {
        line_value <- readLines(connection, n = 1L, warn = FALSE, encoding = "UTF-8")
        if (length(line_value) == 0L) {
            break
        }
        line_counter <- line_counter + 1L
        if (identical(trimws(line_value), "!platform_table_begin")) {
            table_begin_found <- TRUE
            header_value <- readLines(connection, n = 1L, warn = FALSE, encoding = "UTF-8")
            if (length(header_value) == 1L) {
                header_line <- header_value
            }
            break
        }
        if (line_counter > 500000L) {
            break
        }
    }
    if (!table_begin_found || is.na(header_line) || !nzchar(header_line)) {
        return(data.frame(read_status = "PLATFORM_TABLE_NOT_FOUND", n_columns = NA_integer_, first_columns = NA_character_, 
            stringsAsFactors = FALSE))
    }
    header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
    header_fields <- gsub("^\"|\"$", "", header_fields)
    data.frame(read_status = "READ_OK", n_columns = length(header_fields), first_columns = paste(header_fields, 
        collapse = ";"), stringsAsFactors = FALSE)
}

acquire_platform_self <- function(platform_id, source_url, destination_root) {
    destination_path <- file.path(destination_root, paste0(platform_id, "_self_full.soft"))
    download_ok <- try_download_with_retry(source_url, destination_path)
    if (!download_ok) {
        stop("Official GEO self-only platform download failed for ", platform_id, ".")
    }
    header_audit <- read_platform_soft_header(destination_path)
    if (header_audit$read_status != "READ_OK") {
        stop("Downloaded GEO self-only platform file failed table-header validation for ", platform_id, ".")
    }
    data.frame(platform_id = platform_id, source_scope = "SELF_ONLY_PLATFORM_RECORD", source_url = source_url, 
        path = normalize_slashes(destination_path), download_ok = TRUE, read_status = header_audit$read_status, 
        n_columns = header_audit$n_columns, first_columns = header_audit$first_columns, stringsAsFactors = FALSE)
}

strip_prefix <- function(x, prefix) {
    sub(paste0("^", prefix), "", as.character(x))
}

normalize_title <- function(x) {
    tolower(gsub("[^a-z0-9]+", "", as.character(x), perl = TRUE))
}

extract_array_suffix <- function(path) {
    sub("^GSM[0-9]+_", "", basename(path), ignore.case = TRUE)
}

extract_slide_id <- function(path) {
    file_name <- basename(path)
    match_value <- regmatches(file_name, regexpr("SG[0-9]+_[0-9]+", file_name, ignore.case = TRUE, perl = TRUE))
    if (length(match_value) == 0L || !nzchar(match_value)) {
        return(NA_character_)
    }
    toupper(match_value)
}

extract_GSM <- function(path) {
    match_value <- regmatches(basename(path), regexpr("GSM[0-9]+", basename(path), ignore.case = TRUE, perl = TRUE))
    if (length(match_value) == 0L || !nzchar(match_value)) {
        return(NA_character_)
    }
    toupper(match_value)
}

digest_decompressed_gzip <- function(path) {
    connection <- gzfile(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    chunks <- list()
    chunk_index <- 1L
    repeat {
        raw_chunk <- readBin(connection, what = "raw", n = 4L * 1024L * 1024L)
        if (length(raw_chunk) == 0L) {
            break
        }
        chunks[[chunk_index]] <- raw_chunk
        chunk_index <- chunk_index + 1L
    }
    raw_value <- do.call(c, chunks)
    digest::digest(raw_value, algo = "sha256", serialize = FALSE)
}

cramers_v <- function(table_value) {
    chi_value <- suppressWarnings(chisq.test(table_value, correct = FALSE)$statistic)
    n_value <- sum(table_value)
    min_dimension <- min(nrow(table_value) - 1L, ncol(table_value) - 1L)
    if (n_value <= 0L || min_dimension <= 0L) {
        return(NA_real_)
    }
    sqrt(as.numeric(chi_value)/(n_value * min_dimension))
}

log_event("Starting", run_id)

step04B_gate_path <- file.path(dirs$results, "STEP04B_GATE_TABLE.csv")

step04B_locked_root <- file.path(dirs$results, "STEP04B_locked_combined_inputs", "COMBINED_INPUT_SOURCE_v1.0")

step04B_manifest_path <- file.path(step04B_locked_root, "COMBINED_INPUT_SOURCE_v1.0_LOCKED_MANIFEST.csv")

design_candidate_root <- file.path(dirs$results, "STEP04B_sample_design_candidates")

required_upstream <- c(step04B_gate_path, step04B_manifest_path, design_candidate_root)

missing_upstream <- required_upstream[!file.exists(required_upstream) & !dir.exists(required_upstream)]

if (length(missing_upstream) > 0L) {
    stop("Required STEP04B input is missing:\n", paste(missing_upstream, collapse = "\n"))
}

step04B_gate <- .vdr_read_csv(step04B_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

step04B_manifest <- .vdr_read_csv(step04B_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- step04B_gate$status[step04B_gate$gate_component == "STEP04B gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_SAMPLE_DESIGN_LOCK_AND_COMBINED_ANALYSIS_BUILD") {
    stop("STEP04B does not authorize sample-design lock.")
}

if (nrow(step04B_manifest) != 5L || !all(step04B_manifest$lock_status == "LOCKED")) {
    stop("STEP04B locked manifest is incomplete.")
}

step04B_manifest$file_exists <- file.exists(step04B_manifest$path)

step04B_manifest$observed_md5 <- vapply(step04B_manifest$path, md5_file, character(1))

step04B_manifest$observed_sha256 <- vapply(step04B_manifest$path, sha256_file, character(1))

step04B_manifest$identity_pass <- with(step04B_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(step04B_manifest$identity_pass)) {
    stop("A STEP04B locked artifact failed hash validation.")
}

read_design <- function(dataset_id) {
    path <- file.path(design_candidate_root, paste0("STEP04B_", dataset_id, "_SAMPLE_SHEET_CANDIDATE.csv"))
    if (!file.exists(path)) {
        stop("Sample-design candidate is missing: ", path)
    }
    .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

sheet_250130 <- read_design("GSE250130")

sheet_165004 <- read_design("GSE165004")

sheet_92324 <- read_design("GSE92324")

sheet_111974 <- read_design("GSE111974")

sheet_183837 <- read_design("GSE183837")

design_250130 <- data.frame(dataset_id = "GSE250130", sample_accession = sheet_250130$sample_accession, sample_title = sheet_250130$title__1, 
    raw_group = strip_prefix(sheet_250130$characteristics_ch1__3, "group: "), donor_id = sheet_250130$sample_accession, 
    stringsAsFactors = FALSE)

design_250130$clinical_group <- ifelse(design_250130$raw_group == "RIF", "RIF", "FERTILE")

design_250130$LH_day <- ifelse(design_250130$raw_group == "LH3", 3L, ifelse(design_250130$raw_group == "LH5", 5L, 
    ifelse(design_250130$raw_group %in% c("LH7", "LH7_PGT", "RIF"), 7L, ifelse(design_250130$raw_group == "LH9", 
        9L, ifelse(design_250130$raw_group == "LH11", 11L, NA_integer_)))))

design_250130$temporal_fertile_eligible <- design_250130$clinical_group == "FERTILE"

design_250130$LH7_RIF_comparison_eligible <- design_250130$LH_day == 7L

design_250130$LH7_reference_group <- ifelse(design_250130$LH_day != 7L, NA_character_, ifelse(design_250130$clinical_group == 
    "RIF", "RIF", "FERTILE_LH7"))

design_250130$LH7_fertile_subgroup <- ifelse(design_250130$raw_group == "LH7_PGT", "LH7_PGT", ifelse(design_250130$raw_group == 
    "LH7", "LH7_STANDARD", NA_character_))

design_250130$analysis_role <- "TEMPORAL_SINGLE_CELL_DECIDUALIZATION_SUPPORT"

expected_250130_counts <- c(LH3 = 3L, LH5 = 3L, LH7 = 3L, LH7_PGT = 3L, LH9 = 3L, LH11 = 3L, RIF = 10L)

observed_250130_counts <- table(design_250130$raw_group)

if (!setequal(names(observed_250130_counts), names(expected_250130_counts)) || !all(as.integer(observed_250130_counts[names(expected_250130_counts)]) == 
    as.integer(expected_250130_counts)) || sum(design_250130$temporal_fertile_eligible) != 18L || sum(design_250130$LH7_reference_group == 
    "FERTILE_LH7", na.rm = TRUE) != 6L || sum(design_250130$LH7_reference_group == "RIF", na.rm = TRUE) != 10L) {
    stop("GSE250130 design failed 18-fertile/6-LH7/10-RIF reconciliation.")
}

design_165004 <- data.frame(dataset_id = "GSE165004", sample_accession = sheet_165004$sample_accession, sample_title = sheet_165004$title__1, 
    group = strip_prefix(sheet_165004$characteristics_ch1__1, "subject status/group: "), stringsAsFactors = FALSE)

design_165004$group[design_165004$group == "patient with RPL"] <- "RPL"

design_165004$group[design_165004$group == "patient with UIF"] <- "UIF"

design_165004$group[design_165004$group == "Control"] <- "CONTROL"

design_165004$analysis_role <- "PRIMARY_PATIENT_PROJECTION"

step00_raw_root <- file.path(project_root, "05_analysis_steps", "STEP00_dataset_identity_and_technical_audit", 
    "01_inputs", "raw_extracted")

raw_165004 <- list.files(file.path(step00_raw_root, "GSE165004"), pattern = "\\.txt\\.gz$", full.names = TRUE, 
    recursive = FALSE, ignore.case = TRUE)

if (length(raw_165004) != 72L) {
    stop("Expected 72 GSE165004 raw TXT.GZ files but found ", length(raw_165004), ".")
}

raw_165004_map <- data.frame(sample_accession = vapply(raw_165004, extract_GSM, character(1)), raw_file = normalize_slashes(raw_165004), 
    array_suffix = vapply(raw_165004, extract_array_suffix, character(1)), slide_id = vapply(raw_165004, extract_slide_id, 
        character(1)), stringsAsFactors = FALSE)

design_165004 <- merge(design_165004, raw_165004_map, by = "sample_accession", all.x = TRUE, sort = FALSE)

design_165004 <- design_165004[match(sheet_165004$sample_accession, design_165004$sample_accession), , drop = FALSE]

if (any(is.na(design_165004$slide_id)) || length(unique(design_165004$slide_id)) != 9L) {
    stop("GSE165004 slide-ID derivation failed.")
}

slide_group_table <- table(design_165004$slide_id, design_165004$group)

slide_group_nonzero_n <- rowSums(slide_group_table > 0L)

slide_purity_pass <- all(slide_group_nonzero_n == 1L)

slide_cramers_v <- cramers_v(slide_group_table)

if (!slide_purity_pass || abs(slide_cramers_v - 1) > 9.9999999999999998e-13 || any(table(design_165004$group) != 
    24L)) {
    stop("GSE165004 slide/group structure failed the locked phenotype-pure design.")
}

slide_group_counts <- as.data.frame(slide_group_table, stringsAsFactors = FALSE)

names(slide_group_counts) <- c("slide_id", "group", "sample_n")

design_92324 <- data.frame(dataset_id = "GSE92324", sample_accession = sheet_92324$sample_accession, sample_title = sheet_92324$title__1, 
    sample_identifier = strip_prefix(sheet_92324$characteristics_ch1__2, "sample group and identifier: "), age = as.numeric(strip_prefix(sheet_92324$characteristics_ch1__3, 
        "age: ")), collection_day = strip_prefix(sheet_92324$characteristics_ch1__4, "sample collection day: "), 
    stringsAsFactors = FALSE)

design_92324$group <- ifelse(grepl("^control_", design_92324$sample_identifier), "CONTROL", ifelse(grepl("^case_", 
    design_92324$sample_identifier), "IMPLANTATION_FAILURE", NA_character_))

design_92324$analysis_role <- "INDEPENDENT_STIMULATED_CYCLE_VALIDATION"

if (sum(design_92324$group == "CONTROL") != 8L || sum(design_92324$group == "IMPLANTATION_FAILURE") != 10L || any(is.na(design_92324$group))) {
    stop("GSE92324 group assignment failed 8-control/10-case reconciliation.")
}

design_matrix_92324 <- model.matrix(~group + collection_day, data = design_92324)

design_rank_92324 <- qr(design_matrix_92324)$rank

design_full_rank_92324 <- design_rank_92324 == ncol(design_matrix_92324)

if (!design_full_rank_92324) {
    stop("GSE92324 group plus collection-day design matrix is not full rank.")
}

group_day_92324 <- as.data.frame(table(design_92324$group, design_92324$collection_day), stringsAsFactors = FALSE)

names(group_day_92324) <- c("group", "collection_day", "sample_n")

design_111974 <- data.frame(dataset_id = "GSE111974", sample_accession = sheet_111974$sample_accession, sample_title = sheet_111974$title__1, 
    stringsAsFactors = FALSE)

design_111974$group <- ifelse(grepl("_Control_", design_111974$sample_title, fixed = TRUE), "CONTROL", ifelse(grepl("_RIF_", 
    design_111974$sample_title, fixed = TRUE), "RIF", NA_character_))

design_111974$analysis_role <- "SHARED_CONTROL_RIF_DISEASE_CONTEXT_EXTENSION"

if (sum(design_111974$group == "CONTROL") != 24L || sum(design_111974$group == "RIF") != 24L || any(is.na(design_111974$group))) {
    stop("GSE111974 group assignment failed 24-control/24-RIF reconciliation.")
}

raw_111974 <- list.files(file.path(step00_raw_root, "GSE111974"), pattern = "\\.txt\\.gz$", full.names = TRUE, 
    recursive = FALSE, ignore.case = TRUE)

if (length(raw_111974) != 48L) {
    stop("Expected 48 GSE111974 raw TXT.GZ files but found ", length(raw_111974), ".")
}

raw_111974_map <- data.frame(sample_accession = vapply(raw_111974, extract_GSM, character(1)), raw_file = normalize_slashes(raw_111974), 
    array_suffix = vapply(raw_111974, extract_array_suffix, character(1)), slide_id = vapply(raw_111974, extract_slide_id, 
        character(1)), stringsAsFactors = FALSE)

design_111974 <- merge(design_111974, raw_111974_map, by = "sample_accession", all.x = TRUE, sort = FALSE)

design_111974 <- design_111974[match(sheet_111974$sample_accession, design_111974$sample_accession), , drop = FALSE]

control_165004 <- design_165004[design_165004$group == "CONTROL", c("sample_accession", "sample_title", "raw_file", 
    "array_suffix", "slide_id"), drop = FALSE]

control_111974 <- design_111974[design_111974$group == "CONTROL", c("sample_accession", "sample_title", "raw_file", 
    "array_suffix", "slide_id"), drop = FALSE]

control_165004$normalized_title <- normalize_title(control_165004$sample_title)

control_111974$normalized_title <- normalize_title(control_111974$sample_title)

shared_control_pairs <- merge(control_165004, control_111974, by = "normalized_title", suffixes = c("_GSE165004", 
    "_GSE111974"), all = FALSE, sort = TRUE)

if (nrow(shared_control_pairs) != 24L) {
    stop("Expected 24 normalized-title shared control pairs but found ", nrow(shared_control_pairs), ".")
}

shared_control_pairs$array_suffix_match <- shared_control_pairs$array_suffix_GSE165004 == shared_control_pairs$array_suffix_GSE111974

shared_control_pairs$slide_id_match <- shared_control_pairs$slide_id_GSE165004 == shared_control_pairs$slide_id_GSE111974

log_event("Computing decompressed hashes for 24 shared-control array pairs.")

shared_control_pairs$decompressed_sha256_GSE165004 <- vapply(shared_control_pairs$raw_file_GSE165004, digest_decompressed_gzip, 
    character(1))

shared_control_pairs$decompressed_sha256_GSE111974 <- vapply(shared_control_pairs$raw_file_GSE111974, digest_decompressed_gzip, 
    character(1))

shared_control_pairs$decompressed_content_match <- shared_control_pairs$decompressed_sha256_GSE165004 == shared_control_pairs$decompressed_sha256_GSE111974

shared_control_pairs$shared_biological_control_lock <- shared_control_pairs$array_suffix_match & shared_control_pairs$slide_id_match & 
    shared_control_pairs$decompressed_content_match

if (!all(shared_control_pairs$shared_biological_control_lock)) {
    stop("At least one GSE165004/GSE111974 control pair failed exact array-content identity.")
}

design_111974$shared_control_with_GSE165004 <- FALSE

design_111974$paired_GSE165004_sample <- NA_character_

for (row_index in seq_len(nrow(shared_control_pairs))) {
    target_index <- design_111974$sample_accession == shared_control_pairs$sample_accession_GSE111974[row_index]
    design_111974$shared_control_with_GSE165004[target_index] <- TRUE
    design_111974$paired_GSE165004_sample[target_index] <- shared_control_pairs$sample_accession_GSE165004[row_index]
}

design_111974$independent_external_validation <- FALSE

design_111974$validation_scope <- ifelse(design_111974$group == "RIF", "RIF_DISEASE_CONTEXT_EXTENSION_WITH_SHARED_CONTROLS", 
    "SHARED_CONTROL_REUSE_NOT_INDEPENDENT")

design_183837 <- data.frame(dataset_id = "GSE183837", sample_accession = sheet_183837$sample_accession, sample_title = sheet_183837$title__1, 
    group = toupper(strip_prefix(sheet_183837$characteristics_ch1__3, "group: ")), stage = strip_prefix(sheet_183837$characteristics_ch1__2, 
        "Stage: "), donor_id = sheet_183837$sample_accession, stringsAsFactors = FALSE)

design_183837$group[design_183837$group == "CONTROL"] <- "CONTROL"

design_183837$analysis_role <- "INDEPENDENT_SINGLE_CELL_LOCALIZATION"

if (sum(design_183837$group == "CONTROL") != 3L || sum(design_183837$group == "RIF") != 6L || length(unique(design_183837$stage)) != 
    1L) {
    stop("GSE183837 design failed 3-control/6-RIF stage-matched reconciliation.")
}

stale_family_paths <- file.path(platform_root, c("GPL16699_family.soft.gz", "GPL17077_family.soft.gz"))

stale_family_audit <- data.frame(path = normalize_slashes(stale_family_paths), exists = file.exists(stale_family_paths), 
    size_bytes = as.numeric(file.info(stale_family_paths)$size), analysis_use = "NOT_USED_DELETE_AFTER_MANUAL_CONFIRMATION", 
    stringsAsFactors = FALSE)

write_csv_utf8(stale_family_audit, file.path(dirs$qc, "STEP04C_v3_STALE_FAMILY_SOFT_AUDIT.csv"))

GPL16699_self <- acquire_platform_self(platform_id = "GPL16699", source_url = GPL16699_self_url, destination_root = platform_root)

GPL17077_self <- acquire_platform_self(platform_id = "GPL17077", source_url = GPL17077_self_url, destination_root = platform_root)

GPL16699_path <- GPL16699_self$path[[1]]

GPL17077_path <- GPL17077_self$path[[1]]

GPL16699_url <- GPL16699_self$source_url[[1]]

GPL17077_url <- GPL17077_self$source_url[[1]]

platform_download_attempts <- rbind(GPL16699_self, GPL17077_self)

write_csv_utf8(platform_download_attempts, file.path(dirs$qc, "STEP04C_v3_PLATFORM_SELF_DOWNLOAD_AUDIT.csv"))

GPL10558_path <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE92324", "platform_annotation", 
    "GPL10558.soft.gz")

if (!file.exists(GPL10558_path)) {
    stop("Locked GSE92324 GPL10558 annotation is missing.")
}

platform_rows <- list(GPL16699 = cbind(data.frame(platform_id = "GPL16699", source_url = GPL16699_url, path = normalize_slashes(GPL16699_path), 
    stringsAsFactors = FALSE), read_platform_soft_header(GPL16699_path)), GPL17077 = cbind(data.frame(platform_id = "GPL17077", 
    source_url = GPL17077_url, path = normalize_slashes(GPL17077_path), stringsAsFactors = FALSE), read_platform_soft_header(GPL17077_path)), 
    GPL10558 = cbind(data.frame(platform_id = "GPL10558", source_url = "LOCKED_LOCAL_GEO_SOFT_SOURCE", path = normalize_slashes(GPL10558_path), 
        stringsAsFactors = FALSE), read_platform_soft_header(GPL10558_path)))

platform_audit <- do.call(rbind, platform_rows)

platform_audit$size_bytes <- as.numeric(file.info(platform_audit$path)$size)

platform_audit$md5 <- vapply(platform_audit$path, md5_file, character(1))

platform_audit$sha256 <- vapply(platform_audit$path, sha256_file, character(1))

platform_audit$source_version <- platform_source_version

platform_audit$lock_status <- "LOCKED_SOURCE_REFERENCE"

if (any(platform_audit$read_status != "READ_OK") || any(platform_audit$size_bytes < 10000)) {
    stop("At least one platform annotation failed acquisition or header audit.")
}

design_paths <- c(GSE250130 = file.path(locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"), GSE165004 = file.path(locked_root, 
    "COMBINED_SAMPLE_DESIGN_v1.0_GSE165004.csv"), GSE92324 = file.path(locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE92324.csv"), 
    GSE111974 = file.path(locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE111974.csv"), GSE183837 = file.path(locked_root, 
        "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"))

write_csv_utf8(design_250130, design_paths[["GSE250130"]])

write_csv_utf8(design_165004, design_paths[["GSE165004"]])

write_csv_utf8(design_92324, design_paths[["GSE92324"]])

write_csv_utf8(design_111974, design_paths[["GSE111974"]])

write_csv_utf8(design_183837, design_paths[["GSE183837"]])

write_csv_utf8(slide_group_counts, file.path(dirs$results, "STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv"))

write_csv_utf8(group_day_92324, file.path(dirs$results, "STEP04C_GSE92324_GROUP_DAY_COUNTS.csv"))

write_csv_utf8(shared_control_pairs, file.path(dirs$results, "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"))

write_csv_utf8(platform_audit, file.path(dirs$results, "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"))

sample_design_summary <- data.frame(dataset_id = c("GSE250130", "GSE165004", "GSE92324", "GSE111974", "GSE183837"), 
    total_sample_n = c(nrow(design_250130), nrow(design_165004), nrow(design_92324), nrow(design_111974), nrow(design_183837)), 
    locked_design = c("18 fertile temporal; 6 fertile LH7 vs 10 RIF", "24 control; 24 RPL; 24 UIF; 9 phenotype-pure slides", 
        "8 control; 10 implantation failure; HCG day retained", "24 shared controls; 24 RIF; not independent external validation", 
        "3 control; 6 RIF; same post-ovulation stage"), analysis_role = c("TEMPORAL_AND_LH7_RIF_SINGLE_CELL_SUPPORT", 
        "PRIMARY_PATIENT_PROJECTION", "INDEPENDENT_STIMULATED_CYCLE_VALIDATION", "SHARED_CONTROL_RIF_CONTEXT_EXTENSION", 
        "INDEPENDENT_SINGLE_CELL_LOCALIZATION"), lock_version = sample_design_lock_version, stringsAsFactors = FALSE)

write_csv_utf8(sample_design_summary, file.path(dirs$results, "STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv"))

locked_artifact_paths <- c(unname(design_paths), file.path(dirs$results, "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"), 
    file.path(dirs$results, "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"))

locked_manifest <- data.frame(artifact_id = c("GSE250130_SAMPLE_DESIGN", "GSE165004_SAMPLE_DESIGN", "GSE92324_SAMPLE_DESIGN", 
    "GSE111974_SAMPLE_DESIGN", "GSE183837_SAMPLE_DESIGN", "SHARED_CONTROL_IDENTITY_LOCK", "BULK_PLATFORM_ANNOTATION_SOURCE_LOCK"), 
    lock_version = c(rep(sample_design_lock_version, 6L), platform_source_version), path = normalize_slashes(locked_artifact_paths), 
    size_bytes = as.numeric(file.info(locked_artifact_paths)$size), md5 = vapply(locked_artifact_paths, md5_file, 
        character(1)), sha256 = vapply(locked_artifact_paths, sha256_file, character(1)), lock_status = "LOCKED", 
    stringsAsFactors = FALSE)

write_csv_utf8(locked_manifest, file.path(locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv"))

gate_status <- "READY_FOR_COMBINED_BULK_AND_SINGLE_CELL_ANALYSIS_BUILD"

gate_table <- data.frame(gate_component = c("STEP04B gate", "GSE250130 temporal fertile design", "GSE250130 LH7 RIF comparison", 
    "GSE165004 groups", "GSE165004 slide derivation", "GSE165004 slide/group confounding", "GSE92324 groups", "GSE92324 group plus HCG-day model rank", 
    "GSE111974 groups", "GSE111974/GSE165004 shared controls", "GSE111974 external-independence status", "GSE183837 groups and stage", 
    "Bulk platform annotations", "Sample-design lock", "Biological analysis", "Scientific figures", "STEP04C gate"), 
    status = c("PASS", "LOCKED_18_FERTILE_ACROSS_5_TIMEPOINTS", "LOCKED_6_FERTILE_LH7_VS_10_RIF", "LOCKED_24_CONTROL_24_RPL_24_UIF", 
        "PASS_72_SAMPLES_9_SLIDES", paste0("PHENOTYPE_PURE_CRAMERS_V_", format(slide_cramers_v, digits = 6)), "LOCKED_8_CONTROL_10_IF", 
        paste0("PASS_FULL_RANK_", design_rank_92324, "_OF_", ncol(design_matrix_92324)), "LOCKED_24_CONTROL_24_RIF", 
        "LOCKED_24_IDENTICAL_CONTROL_ARRAYS", "NOT_INDEPENDENT_SHARED_CONTROL_EXTENSION", "LOCKED_3_CONTROL_6_RIF_STAGE_MATCHED", 
        "LOCKED_GPL16699_GPL17077_GPL10558", paste0("LOCKED_", sample_design_lock_version), "NOT_STARTED", "NONE_CREATED", 
        gate_status), required_action = c("Retain the STEP04B source lock.", "Use donor-aware temporal single-cell summaries.", 
        "Compare RIF only with fertile LH+7 donors.", "Use separate Control-vs-RPL and Control-vs-UIF projections.", 
        "Use the derived slide IDs for sensitivity analyses.", "Do not include slide as a standard covariate because it is phenotype-pure; use leave-one-slide sensitivity.", 
        "Retain HCG+6/+7 metadata.", "Use group plus collection day when gene-level modeling is required.", "Use as RIF disease-context support.", 
        "Do not double-count controls across datasets.", "Do not call GSE111974 an independent external validation.", 
        "Use donor-aware single-cell summaries.", "Build gene-level matrices only from locked platform sources.", 
        "Do not change assignments without a new version.", "Proceed as one combined analysis package.", "Render only after biological QC.", 
        "Build bulk and single-cell analysis scripts."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP04C_GATE_TABLE.csv"))

write_text(c("# STEP04C sample-design and platform lock card", "", paste0("- Run ID: ", run_id), paste0("- Sample-design lock: ", 
    sample_design_lock_version, " LOCKED"), "- GSE250130: 18 fertile temporal samples; 6 fertile LH+7 versus 10 RIF", 
    "- GSE165004: 24 control, 24 RPL, 24 UIF; 9 phenotype-pure slides", "- GSE92324: 8 control, 10 implantation failure; HCG+6/+7 retained", 
    "- GSE111974: 24 RIF plus 24 controls identical to GSE165004 controls", "- GSE111974 independent external validation: NO", 
    "- GSE183837: 3 control, 6 RIF; stage matched", "- Shared-control decompressed content matches: 24/24", "- Platform sources: GPL16699, GPL17077 and GPL10558 LOCKED", 
    "- Biological analysis: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, 
    "STEP04C_SAMPLE_DESIGN_PLATFORM_LOCK_CARD.md"))

write_text(c("# Combined biological analysis authorization", "", paste0("- Authorized by: ", run_id), paste0("- Sample-design source: ", 
    sample_design_lock_version, " LOCKED"), "- Authorization: YES", "", "## Bulk patient analyses", "", "- GSE165004: primary projection, Control vs UIF and Control vs RPL separately, with leave-one-slide sensitivity.", 
    "- GSE92324: independent stimulated-cycle validation, retaining HCG collection day.", "- GSE111974: RIF disease-context extension sharing all 24 control arrays with GSE165004; never count it as an independent cohort.", 
    "", "## Single-cell analyses", "", "- GSE250130: fertile time series across LH+3/+5/+7/+9/+11 and RIF-vs-fertile comparison restricted to LH+7.", 
    "- GSE183837: independent Ctrl/RIF localization at five days after ovulation.", "- Use donor-aware cell-type pseudobulk or sample-level scores; cells are not independent replicates.", 
    "", "## Scope limits", "", "- Use the locked 284 occupancy-supported RNA candidates as the primary projection set.", 
    "- Retain the 152 RNA-only candidates for secondary sensitivity analysis.", "- Do not use discovery ROC, LASSO, random forest or generic network pharmacology.", 
    "- Do not define definitive direct targets.", "- Do not generate scientific figures until analysis QC passes."), 
    file.path(dirs$control, "COMBINED_BULK_SINGLE_CELL_ANALYSIS_AUTHORIZATION.md"))

write_text(c("# STEP04C Methods note", "", "Sample assignments were locked from GEO metadata and raw array identifiers. For GSE250130, fertile samples at LH+3, LH+5, LH+7, LH+9 and LH+11 were retained for temporal analysis, whereas the RIF comparison was restricted to the six fertile and ten RIF samples collected at LH+7.", 
    "", "For GSE165004, raw Agilent filenames were used to reconstruct nine slide identifiers. Slides were phenotype-pure, so slide was treated through leave-one-slide sensitivity rather than as a conventional adjustment covariate. The 24 control arrays in GSE111974 were confirmed to be the same arrays used in GSE165004 by normalized title, array identifier and decompressed raw-file content.", 
    "", "GSE92324 retained collection day as a design variable, and GSE183837 retained the donor-level 3-control/6-RIF design. Platform annotations were locked before gene-level expression mapping."), 
    file.path(dirs$text, "STEP04C_Methods_note.md"))

parameters <- list(schema_version = "1.0-lock", run_id = run_id, sample_design_lock_version = sample_design_lock_version, 
    platform_source_version = platform_source_version, GSE250130_temporal_fertile_n = 18, GSE250130_fertile_LH7_n = 6, 
    GSE250130_RIF_LH7_n = 10, GSE165004_control_n = 24, GSE165004_RPL_n = 24, GSE165004_UIF_n = 24, GSE165004_slide_n = 9, 
    GSE165004_slide_cramers_v = slide_cramers_v, GSE92324_control_n = 8, GSE92324_IF_n = 10, GSE92324_model_full_rank = design_full_rank_92324, 
    GSE111974_control_n = 24, GSE111974_RIF_n = 24, shared_control_identity_n = sum(shared_control_pairs$shared_biological_control_lock), 
    GSE111974_independent_external_validation = FALSE, GSE183837_control_n = 3, GSE183837_RIF_n = 6, biological_analysis_started = FALSE, 
    scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04C_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04C.txt"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    copy_if_different(source_file, file.path(dirs$code, "prepare_patient_designs.R"))
}

archive_artifacts <- c(file.path(dirs$results, "STEP04C_GATE_TABLE.csv"), file.path(dirs$results, "STEP04C_SAMPLE_DESIGN_PLATFORM_LOCK_CARD.md"), 
    file.path(dirs$results, "STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv"), file.path(dirs$results, "STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv"), 
    file.path(dirs$results, "STEP04C_GSE92324_GROUP_DAY_COUNTS.csv"), file.path(dirs$results, "STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv"), 
    file.path(dirs$results, "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv"), file.path(dirs$qc, "STEP04C_v3_PLATFORM_SELF_DOWNLOAD_AUDIT.csv"), 
    file.path(dirs$qc, "STEP04C_v3_STALE_FAMILY_SOFT_AUDIT.csv"), file.path(locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv"), 
    unname(design_paths), file.path(dirs$control, "COMBINED_BULK_SINGLE_CELL_ANALYSIS_AUTHORIZATION.md"), file.path(dirs$text, 
        "STEP04C_Methods_note.md"), file.path(dirs$parameters, "STEP04C_parameters.json"), file.path(dirs$environment, 
        "sessionInfo_POST_RUN_STEP04C.txt"), file.path(dirs$code, "prepare_patient_designs.R"), log_path)

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
    else if (grepl("Methods|AUTHORIZATION", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("LOCK|CARD|COUNTS|MANIFEST|SUMMARY|AUDIT|GATE", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04C_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04C archive status", "", paste0("- Run ID: ", run_id), paste0("- Sample-design lock: ", sample_design_lock_version), 
    "- Shared control identity: 24/24", "- Platform annotations: 3 locked sources", "- Raw data copied into archive: NO", 
    "- Biological analysis: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(archive_root, 
    "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Accepted run: ", run_id), paste0("- Sample designs: ", 
    sample_design_lock_version, " LOCKED"), paste0("- Platform sources: ", platform_source_version, " LOCKED"), 
    "- GSE111974 role: SHARED_CONTROL_RIF_CONTEXT_EXTENSION", "- Biological analysis: AUTHORIZED_NOT_STARTED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04C v3 sample-design and platform preflight completed.\n")

cat("Sample-design lock: ", sample_design_lock_version, " LOCKED\n")

cat("GSE250130 temporal fertile samples: 18\n")

cat("GSE250130 LH7 comparison: 6 fertile vs 10 RIF\n")

cat("GSE165004: 24 control / 24 RPL / 24 UIF\n")

cat("GSE165004 phenotype-pure slides: 9; Cramer's V = ", format(slide_cramers_v, digits = 6), "\n")

cat("GSE92324: 8 control / 10 implantation failure\n")

cat("GSE92324 group + HCG-day model: FULL RANK\n")

cat("GSE111974: 24 shared controls / 24 RIF\n")

cat("Shared control decompressed identity: 24/24\n")

cat("GSE111974 independent external validation: NO\n")

cat("GSE183837: 3 control / 6 RIF; stage matched\n")

cat("Platform annotations locked: GPL16699 / GPL17077 / GPL10558\n")

cat("Biological analysis: AUTHORIZED_NOT_STARTED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP04C_GATE_TABLE.csv\n")

cat("  06_results/STEP04C_COMBINED_SAMPLE_DESIGN_LOCK_SUMMARY.csv\n")

cat("  06_results/STEP04C_GSE165004_GSE111974_SHARED_CONTROL_IDENTITY_LOCK.csv\n")

cat("  06_results/STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv\n")

cat("  06_results/STEP04C_GSE165004_SLIDE_GROUP_COUNTS.csv\n")

cat("  06_results/STEP04C_GSE92324_GROUP_DAY_COUNTS.csv\n")

cat("============================================================\n")

