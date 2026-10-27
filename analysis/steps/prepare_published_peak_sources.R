#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A07c: Prepare author peak-table sources and IgG context

options(stringsAsFactors = FALSE)

options(timeout = 1800)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3A1_", format(Sys.time(), "%Y%m%d_%H%M%S"))

article_id <- "PMC12642428"

supplement_url <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/", article_id, "/supplementaryFiles")

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

raw_supplement_root <- file.path(project_root, "03_data", "raw_external_READ_ONLY", article_id, "supplementary")

dir.create(raw_supplement_root, recursive = TRUE, showWarnings = FALSE)

supplement_zip <- file.path(raw_supplement_root, paste0(article_id, "_supplementaryFiles.zip"))

supplement_extract_root <- file.path(raw_supplement_root, "extracted")

dir.create(supplement_extract_root, recursive = TRUE, showWarnings = FALSE)

input_reference_root <- file.path(dirs$inputs, "STEP03B3B_published_peak_table_refs")

dir.create(input_reference_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3A1_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

required_packages <- c("digest", "jsonlite", "readxl")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
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

zip_magic_pass <- function(path) {
    if (!file.exists(path)) {
        return(FALSE)
    }
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    magic <- readBin(connection, what = "raw", n = 4L)
    length(magic) == 4L && identical(as.integer(magic)[1:2], c(80L, 75L))
}

download_with_fallback <- function(url, destination) {
    temporary_path <- paste0(destination, ".partial")
    if (file.exists(temporary_path)) {
        unlink(temporary_path)
    }
    methods_to_try <- unique(c("libcurl", if (.Platform$OS.type == "windows") {
        "wininet"
    } else {
        character(0)
    }))
    success <- FALSE
    messages <- character(0)
    for (method_name in methods_to_try) {
        result <- tryCatch({
            utils::download.file(url = url, destfile = temporary_path, method = method_name, mode = "wb", quiet = FALSE)
            TRUE
        }, error = function(error_object) {
            messages <<- c(messages, paste0(method_name, ": ", conditionMessage(error_object)))
            FALSE
        })
        if (isTRUE(result) && file.exists(temporary_path) && file.info(temporary_path)$size > 100000L && zip_magic_pass(temporary_path)) {
            success <- TRUE
            break
        }
        if (file.exists(temporary_path)) {
            unlink(temporary_path)
        }
    }
    if (!success) {
        curl_executable <- Sys.which("curl.exe")
        if (nzchar(curl_executable)) {
            curl_output <- tryCatch(system2(curl_executable, args = c("-L", "--fail", "--retry", "4", "--retry-delay", 
                "3", "--connect-timeout", "30", "--max-time", "1800", "-o", shQuote(normalize_slashes(temporary_path)), 
                shQuote(url)), stdout = TRUE, stderr = TRUE), error = function(error_object) {
                paste0("curl.exe: ", conditionMessage(error_object))
            })
            curl_status <- attr(curl_output, "status")
            if ((is.null(curl_status) || curl_status == 0L) && file.exists(temporary_path) && file.info(temporary_path)$size > 
                100000L && zip_magic_pass(temporary_path)) {
                success <- TRUE
            }
            else {
                messages <- c(messages, paste(curl_output, collapse = " "))
            }
        }
    }
    if (success) {
        file.copy(temporary_path, destination, overwrite = TRUE)
        unlink(temporary_path)
        return(list(success = TRUE, messages = messages))
    }
    if (file.exists(temporary_path)) {
        unlink(temporary_path)
    }
    list(success = FALSE, messages = messages)
}

log_event("Starting", run_id)

required_paths <- c(gate = file.path(dirs$results, "STEP03B3A_GATE_TABLE.csv"), audit = file.path(dirs$results, 
    "STEP03B3A_CUTRUN_BED_integrity_audit.csv"), depth = file.path(dirs$qc, "STEP03B3A_condition_depth_summary.csv"), 
    semantics = file.path(dirs$qc, "STEP03B3A_record_semantics_summary.csv"), IgG_map = file.path(dirs$results, 
        "STEP03B3A_MATCHED_IgG_CONTROL_MAP_LOCKED.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3A output is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate <- .vdr_read_csv(required_paths[["gate"]], stringsAsFactors = FALSE, check.names = FALSE)

audit <- .vdr_read_csv(required_paths[["audit"]], stringsAsFactors = FALSE, check.names = FALSE)

depth <- .vdr_read_csv(required_paths[["depth"]], stringsAsFactors = FALSE, check.names = FALSE)

semantics <- .vdr_read_csv(required_paths[["semantics"]], stringsAsFactors = FALSE, check.names = FALSE)

IgG_map <- .vdr_read_csv(required_paths[["IgG_map"]], stringsAsFactors = FALSE, check.names = FALSE)

step03b3a_status <- gate$status[gate$gate_component == "STEP03B3A gate"]

if (length(step03b3a_status) != 1L || step03b3a_status != "HOLD_CUTRUN_BED_INTEGRITY_OR_DEPTH_QC") {
    stop("Unexpected STEP03B3A gate: ", paste(step03b3a_status, collapse = ", "))
}

if (nrow(audit) != 48L || !all(audit$audit_status == "PASS") || !all(audit$invalid_rows == 0L) || !all(audit$primary_fraction >= 
    0.94999999999999996)) {
    stop("The 48-file BED integrity audit does not support manual release.")
}

if (nrow(semantics) != 1L || semantics$provisional_record_semantics != "FRAGMENT_LIKE_INTERVALS" || semantics$n_files != 
    48L) {
    stop("CUT&RUN record semantics are not uniformly fragment-like.")
}

if (!all(IgG_map$pairing_pass)) {
    stop("Matched IgG routing is incomplete.")
}

depth_review <- depth[depth$cv_rows > 0.02, , drop = FALSE]

expected_low_depth_condition <- paste0("CUTRUN_VDR_cistrome__VDR_overexpression__", "calcitriol__IgG_control")

if (nrow(depth_review) != 1L || depth_review$condition_id_locked != expected_low_depth_condition) {
    stop("The depth exception is not the single prespecified IgG condition.")
}

low_depth_sample <- audit[audit$condition_id_locked == expected_low_depth_condition & audit$valid_rows == min(audit$valid_rows[audit$condition_id_locked == 
    expected_low_depth_condition]), , drop = FALSE]

if (nrow(low_depth_sample) != 1L || low_depth_sample$gsm != "GSM9193665" || low_depth_sample$valid_rows != 6375100L) {
    stop("The expected low-depth control GSM9193665 was not uniquely identified.")
}

low_depth_decision <- data.frame(gsm = low_depth_sample$gsm, condition_id_locked = low_depth_sample$condition_id_locked, 
    valid_rows = low_depth_sample$valid_rows, peer_median_rows = stats::median(audit$valid_rows[audit$condition_id_locked == 
        expected_low_depth_condition & audit$gsm != low_depth_sample$gsm]), depth_ratio_to_peer_median = low_depth_sample$valid_rows/stats::median(audit$valid_rows[audit$condition_id_locked == 
        expected_low_depth_condition & audit$gsm != low_depth_sample$gsm]), invalid_fraction = low_depth_sample$invalid_fraction, 
    primary_fraction = low_depth_sample$primary_fraction, median_width_approx = low_depth_sample$median_width_approx, 
    p95_width_approx = low_depth_sample$p95_width_approx, approximate_duplicate_fraction = low_depth_sample$approximate_duplicate_fraction, 
    integrity_status = "PASS_COMPLETE_FILE", depth_status = "LOWER_THAN_MATCHED_PEERS", sample_decision = "COMPLETE_LOW_DEPTH_CONTROL_RETAINED", 
    sample_excluded = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(low_depth_decision, file.path(dirs$results, "STEP03B3A1_LOW_DEPTH_IgG_DECISION.csv"))

local_peak_sensitivity_rules <- data.frame(rule_id = c("LP01", "LP02", "LP03", "LP04", "LP05", "LP06"), rule = c("Retain GSM9193665 and its matched target replicate.", 
    "Use paired-end fragment semantics for the three-column BED intervals.", "For local peak calling, scale toward the smaller library rather than scaling the small IgG upward.", 
    "Primary condition consensus requires support in at least two of three replicate-specific target-versus-IgG calls.", 
    "GSM9193665 cannot define a primary consensus peak by itself.", "For VDR-overexpression plus calcitriol, compare the 2-of-3 consensus with a replicate-1-and-2-only sensitivity consensus."), 
    status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(local_peak_sensitivity_rules, file.path(dirs$results, "STEP03B3A1_LOCAL_PEAK_CALLING_SENSITIVITY_RULES.csv"))

download_status <- "EXISTING_VALID_ZIP_REUSED"

download_messages <- character(0)

if (!file.exists(supplement_zip) || file.info(supplement_zip)$size <= 100000L || !zip_magic_pass(supplement_zip)) {
    log_event("Downloading supplementary archive from Europe PMC")
    download_result <- download_with_fallback(supplement_url, supplement_zip)
    download_messages <- download_result$messages
    if (!download_result$success) {
        download_status <- "DOWNLOAD_FAILED_MANUAL_ACTION_REQUIRED"
    }
    else {
        download_status <- "DOWNLOADED_AND_VALIDATED"
    }
}

if (download_status == "DOWNLOAD_FAILED_MANUAL_ACTION_REQUIRED") {
    manual_download_note <- c("# Manual supplementary download required", "", paste0("The automated download failed for: `", 
        supplement_url, "`."), "", "Open this Europe PMC endpoint in a web browser and save the downloaded ZIP as:", 
        paste0("`", normalize_slashes(supplement_zip), "`."), "", "Then rerun this script.", "", "Messages:", paste(download_messages, 
            collapse = "\n"))
    write_text(manual_download_note, file.path(dirs$text, "STEP03B3A1_MANUAL_SUPPLEMENT_DOWNLOAD_REQUIRED.md"))
    hold_gate <- data.frame(gate_component = c("Low-depth IgG decision", "Sample exclusion", "Published supplementary archive", 
        "STEP03B3A1 gate"), status = c("COMPLETE_LOW_DEPTH_CONTROL_RETAINED", "NONE", "MANUAL_DOWNLOAD_REQUIRED", 
        "HOLD_SUPPLEMENT_DOWNLOAD_REQUIRED"), required_action = c("Retain GSM9193665.", "Do not remove samples.", 
        "Download the Europe PMC supplementary ZIP to the locked path.", "Rerun STEP03B3A1."), stringsAsFactors = FALSE)
    write_csv_utf8(hold_gate, file.path(dirs$results, "STEP03B3A1_GATE_TABLE.csv"))
    stop("Supplementary archive download failed. See STEP03B3A1_MANUAL_SUPPLEMENT_DOWNLOAD_REQUIRED.md")
}

zip_listing <- utils::unzip(supplement_zip, list = TRUE)

write_csv_utf8(zip_listing, file.path(dirs$inputs, "STEP03B3A1_SUPPLEMENT_ZIP_LISTING.csv"))

utils::unzip(supplement_zip, exdir = supplement_extract_root, overwrite = TRUE)

extracted_files <- list.files(supplement_extract_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

if (length(extracted_files) == 0L) {
    stop("The supplementary ZIP contained no extractable files.")
}

supplement_manifest <- data.frame(file_name = basename(extracted_files), relative_path = substring(normalize_slashes(extracted_files), 
    nchar(normalize_slashes(supplement_extract_root)) + 2L), extension = tolower(tools::file_ext(extracted_files)), 
    size_bytes = as.numeric(file.info(extracted_files)$size), md5 = vapply(extracted_files, md5_file, character(1)), 
    sha256 = vapply(extracted_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(supplement_manifest, file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_SUPPLEMENT_FILE_MANIFEST.csv"))

xlsx_paths <- extracted_files[tolower(tools::file_ext(extracted_files)) == "xlsx"]

if (length(xlsx_paths) == 0L) {
    stop("No XLSX supplementary workbook was found.")
}

workbook_rows <- list()

sheet_rows <- list()

candidate_rows <- list()

peak_keywords <- c("peak", "vdr", "cistrom", "h3k4", "h3k27", "histone", "motif", "deg", "chrom")

for (workbook_index in seq_along(xlsx_paths)) {
    workbook_path <- xlsx_paths[workbook_index]
    sheet_names <- readxl::excel_sheets(workbook_path)
    workbook_rows[[workbook_index]] <- data.frame(workbook_file = basename(workbook_path), workbook_path = normalize_slashes(workbook_path), 
        workbook_size_bytes = as.numeric(file.info(workbook_path)$size), workbook_md5 = md5_file(workbook_path), 
        workbook_sha256 = sha256_file(workbook_path), n_sheets = length(sheet_names), stringsAsFactors = FALSE)
    for (sheet_index in seq_along(sheet_names)) {
        sheet_name <- sheet_names[sheet_index]
        sheet_data <- tryCatch(readxl::read_excel(workbook_path, sheet = sheet_name, .name_repair = "unique_quiet"), 
            error = function(error_object) {
                structure(data.frame(), read_error = conditionMessage(error_object))
            })
        read_error <- attr(sheet_data, "read_error")
        if (is.null(read_error)) {
            read_error <- ""
        }
        column_names <- names(sheet_data)
        searchable_text <- tolower(paste(basename(workbook_path), sheet_name, paste(column_names, collapse = " ")))
        matched_keywords <- peak_keywords[vapply(peak_keywords, function(keyword_value) {
            grepl(keyword_value, searchable_text, fixed = TRUE)
        }, logical(1))]
        sheet_record <- data.frame(workbook_file = basename(workbook_path), sheet_name = sheet_name, n_rows = nrow(sheet_data), 
            n_columns = ncol(sheet_data), column_names = paste(column_names, collapse = " | "), read_error = read_error, 
            candidate_keyword_hits = paste(matched_keywords, collapse = ";"), candidate_for_peak_review = length(matched_keywords) > 
                0L, stringsAsFactors = FALSE)
        sheet_rows[[length(sheet_rows) + 1L]] <- sheet_record
        if (sheet_record$candidate_for_peak_review && !nzchar(read_error)) {
            candidate_rows[[length(candidate_rows) + 1L]] <- sheet_record
        }
    }
}

workbook_inventory <- do.call(rbind, workbook_rows)

sheet_inventory <- do.call(rbind, sheet_rows)

candidate_sheet_inventory <- if (length(candidate_rows) > 0L) {
    do.call(rbind, candidate_rows)
} else {
    sheet_inventory[FALSE, , drop = FALSE]
}

write_csv_utf8(workbook_inventory, file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_WORKBOOK_INVENTORY.csv"))

write_csv_utf8(sheet_inventory, file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_WORKBOOK_SHEET_INVENTORY.csv"))

write_csv_utf8(candidate_sheet_inventory, file.path(dirs$results, "STEP03B3A1_PUBLISHED_PEAK_SHEET_SHORTLIST.csv"))

source_reference_manifest <- data.frame(source_id = c("ARTICLE_SUPPLEMENT_ARCHIVE", paste0("WORKBOOK_", seq_len(nrow(workbook_inventory)))), 
    role = c("Open-access supplementary archive", rep("Author-processed supplementary workbook", nrow(workbook_inventory))), 
    source_path = c(normalize_slashes(supplement_zip), workbook_inventory$workbook_path), size_bytes = c(as.numeric(file.info(supplement_zip)$size), 
        workbook_inventory$workbook_size_bytes), md5 = c(md5_file(supplement_zip), workbook_inventory$workbook_md5), 
    sha256 = c(sha256_file(supplement_zip), workbook_inventory$workbook_sha256), source_status = "LOCKED_READ_ONLY_REFERENCE", 
    stringsAsFactors = FALSE)

write_csv_utf8(source_reference_manifest, file.path(input_reference_root, "STEP03B3A1_PUBLISHED_SOURCE_REFERENCE_MANIFEST.csv"))

candidate_sheet_count <- nrow(candidate_sheet_inventory)

final_gate_status <- if (candidate_sheet_count > 0L) {
    "PROCEED_TO_MANUAL_PUBLISHED_PEAK_TABLE_REVIEW"
} else {
    "HOLD_NO_CANDIDATE_PEAK_SHEETS_IDENTIFIED"
}

gate_table <- data.frame(gate_component = c("STEP03B3A BED integrity", "Record semantics", "Matched IgG routing", 
    "GSM9193665 integrity", "GSM9193665 depth", "Sample exclusion", "Local peak-calling sensitivity rules", "Published supplementary archive", 
    "Published workbook inventory", "Candidate peak-sheet shortlist", "Local peak calling", "RNA/CUT&RUN integration", 
    "STEP03B3A1 gate"), status = c("PASS_48_OF_48", "FRAGMENT_LIKE_INTERVALS_48_OF_48", "PASS_10_OF_10", "PASS_COMPLETE_FILE", 
    "LOW_DEPTH_CONTROL_RETAINED", "NONE", "LOCKED", download_status, paste0("PASS_", nrow(workbook_inventory), 
        "_WORKBOOKS"), paste0("N_", candidate_sheet_count), "NOT_STARTED", "NOT_STARTED", final_gate_status), required_action = c("Retain hashes and audit table.", 
    "Treat three-column records as paired-end fragments.", "Use matched IgG for every target condition.", "Do not classify GSM9193665 as corrupt.", 
    "Report the depth exception in all downstream methods and QC.", "Retain all 48 samples.", "Apply the six locked sensitivity rules if local peak calling is later used.", 
    "Keep the archive and workbooks read-only.", "Use sheet inventory to locate author-processed peak tables.", 
    "Perform manual schema and scientific-role review.", "Do not install or run MACS3 before published-table review.", 
    "Do not integrate until peak-table schema is locked.", "Upload the inventories and candidate-sheet shortlist."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3A1_GATE_TABLE.csv"))

write_text(c("# STEP03B3A1 low-depth IgG lock card", "", paste0("- Run ID: ", run_id), "- Sample: GSM9193665", 
    "- Condition: VDR overexpression + calcitriol + IgG", "- File integrity: PASS", "- Valid intervals: 6,375,100", 
    "- Matched peer intervals: 10,500,100 and 10,500,100", "- Classification: COMPLETE_LOW_DEPTH_CONTROL_RETAINED", 
    "- Sample exclusion: NO", "- Primary replicate consensus rule: at least 2 of 3", "- Replicate-3-only peaks: prohibited from primary consensus", 
    "- Published peak-table audit: primary next action", paste0("- Gate: ", final_gate_status)), file.path(dirs$results, 
    "STEP03B3A1_LOW_DEPTH_IgG_LOCK_CARD.md"))

write_text(c("# STEP03B3B transition", "", "The processed CUT&RUN fragment layer passed structural integrity QC.", 
    "One complete but lower-depth IgG file, GSM9193665, is retained with", "explicit sensitivity rules and without sample deletion.", 
    "", "The next task is to inspect the author-processed supplementary peak", "tables and determine whether they provide condition-specific peak", 
    "coordinates, scores, motifs and transcriptomic overlap suitable for", "the project integration. Local peak calling remains deferred.", 
    "", paste0("Current gate: ", final_gate_status)), file.path(dirs$text, "STEP03B3A1_TRANSITION_TO_PUBLISHED_PEAK_TABLE_REVIEW.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, article_id = article_id, supplement_url = supplement_url, 
    low_depth_sample = "GSM9193665", low_depth_sample_decision = "COMPLETE_LOW_DEPTH_CONTROL_RETAINED", samples_excluded = 0, 
    local_peak_calling = FALSE, RNA_integration = FALSE, candidate_sheet_count = candidate_sheet_count, gate = final_gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3A1_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3A1.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("32_step03B3A1_low_depth_IgG_and_published_peak_table_audit.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "prepare_published_peak_sources.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3A1_LOW_DEPTH_IgG_DECISION.csv"), file.path(dirs$results, 
    "STEP03B3A1_LOCAL_PEAK_CALLING_SENSITIVITY_RULES.csv"), file.path(dirs$results, "STEP03B3A1_PUBLISHED_PEAK_SHEET_SHORTLIST.csv"), 
    file.path(dirs$results, "STEP03B3A1_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3A1_LOW_DEPTH_IgG_LOCK_CARD.md"), 
    file.path(dirs$inputs, "STEP03B3A1_SUPPLEMENT_ZIP_LISTING.csv"), file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_SUPPLEMENT_FILE_MANIFEST.csv"), 
    file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_WORKBOOK_INVENTORY.csv"), file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_WORKBOOK_SHEET_INVENTORY.csv"), 
    file.path(input_reference_root, "STEP03B3A1_PUBLISHED_SOURCE_REFERENCE_MANIFEST.csv"), file.path(dirs$text, 
        "STEP03B3A1_TRANSITION_TO_PUBLISHED_PEAK_TABLE_REVIEW.md"), file.path(dirs$parameters, "STEP03B3A1_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3A1.txt"), file.path(dirs$code, "prepare_published_peak_sources.R"), 
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
    else if (grepl("TRANSITION", basename(path_value))) {
        "txt"
    }
    else if (grepl("DECISION|RULES|LOCK_CARD", basename(path_value))) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3A1_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3A1 archive status", "", paste0("- Run ID: ", run_id), "- GSM9193665: COMPLETE_LOW_DEPTH_CONTROL_RETAINED", 
    "- Sample exclusions: NONE", paste0("- Published workbooks inventoried: ", nrow(workbook_inventory)), paste0("- Candidate peak sheets: ", 
        candidate_sheet_count), "- Local peak calling: NOT_STARTED", paste0("- Gate: ", final_gate_status), paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3A1", paste0("- Run ID: ", run_id), "- RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE", 
    "- ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", "- CUT&RUN fragment integrity: PASS", "- Low-depth IgG: RETAINED WITH SENSITIVITY RULES", 
    "- Sample exclusions: NONE", "- Published peak-table audit: COMPLETED_PENDING_MANUAL_REVIEW", "- Local peak calling: NOT_STARTED", 
    "- RNA/CUT&RUN integration: NOT_STARTED", paste0("- Gate: ", final_gate_status), "- Full STEP03 module: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", final_gate_status)

cat("\n============================================================\n")

cat("STEP03B3A1 low-depth IgG and supplementary audit completed.\n")

cat("GSM9193665: COMPLETE_LOW_DEPTH_CONTROL_RETAINED\n")

cat("Samples excluded: 0\n")

cat("Published workbooks inventoried: ", nrow(workbook_inventory), "\n")

cat("Candidate peak sheets: ", candidate_sheet_count, "\n")

cat("Local peak calling: NOT STARTED\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Gate: ", final_gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3A1_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3A1_LOW_DEPTH_IgG_DECISION.csv\n")

cat("  06_results/STEP03B3A1_PUBLISHED_PEAK_SHEET_SHORTLIST.csv\n")

cat("  01_inputs/STEP03B3A1_PUBLISHED_WORKBOOK_INVENTORY.csv\n")

cat("  01_inputs/STEP03B3A1_PUBLISHED_WORKBOOK_SHEET_INVENTORY.csv\n")

cat("============================================================\n")

