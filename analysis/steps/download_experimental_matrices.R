#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A04: Acquire processed expression and peak sources

options(stringsAsFactors = FALSE)

options(timeout = 7200)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03A_v5_", format(Sys.time(), "%Y%m%d_%H%M%S"))

download_atac_bigwig <- TRUE

max_retries <- 4L

retry_wait_seconds <- 5L

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

invisible(NULL)

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    archive = file.path(module_root, "10_archive_work", run_id))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

progress_path <- file.path(dirs$results, "STEP03A_v5_download_progress.csv")

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

required_pkgs <- c("curl", "digest", "jsonlite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

if (!all(vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required packages could not be loaded.")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

md5_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    unname(tools::md5sum(path))
}

sha256_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

destination_subdir <- function(accession, file_role) {
    if (accession == "GSE254251" && file_role == "processed_expression_matrix") {
        return("expression")
    }
    if (accession == "GSE306094" && grepl("^merged_peak_set_", file_role)) {
        return("merged_peaks")
    }
    if (accession == "GSE306094" && file_role == "replicate_bigwig_same_depth") {
        return("bigwig_same_depth")
    }
    if (accession == "GSE306127" && file_role == "replicate_BED_same_depth") {
        return("bed_same_depth")
    }
    "other"
}

validate_gzip_text <- function(path, expected_type) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    lines <- tryCatch(readLines(con, n = 20L, warn = FALSE), error = function(e) e)
    if (inherits(lines, "error")) {
        return(list(pass = FALSE, detail = paste0("gzip read failed: ", conditionMessage(lines))))
    }
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0L) {
        return(list(pass = FALSE, detail = "gzip file contains no readable nonempty line"))
    }
    if (expected_type == "expression") {
        fields <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
        return(list(pass = length(fields) >= 2L, detail = paste0("header_fields=", length(fields))))
    }
    data_lines <- lines[!grepl("^(#|track\\b|browser\\b)", trimws(lines), ignore.case = TRUE, perl = TRUE)]
    if (length(data_lines) == 0L) {
        return(list(pass = FALSE, detail = "no BED-like data line in first 20 readable lines"))
    }
    fields <- strsplit(data_lines[1], "\t", fixed = TRUE)[[1]]
    numeric_ok <- FALSE
    if (length(fields) >= 3L) {
        start_value <- suppressWarnings(as.numeric(fields[2]))
        end_value <- suppressWarnings(as.numeric(fields[3]))
        numeric_ok <- is.finite(start_value) && is.finite(end_value) && end_value >= start_value
    }
    list(pass = length(fields) >= 3L && numeric_ok, detail = paste0("bed_fields=", length(fields), "; coordinates_numeric=", 
        numeric_ok))
}

validate_bigwig <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    magic <- readBin(con, what = "raw", n = 4L)
    little_endian <- as.raw(c(38, 252, 143, 136))
    big_endian <- rev(little_endian)
    pass <- identical(magic, little_endian) || identical(magic, big_endian)
    list(pass = pass, detail = paste(sprintf("%02X", as.integer(magic)), collapse = " "))
}

validate_download <- function(path, expected_size, file_role) {
    if (!file.exists(path)) {
        return(list(pass = FALSE, size_pass = FALSE, format_pass = FALSE, detail = "file missing"))
    }
    actual_size <- as.numeric(file.info(path)$size)
    size_pass <- is.finite(expected_size) && actual_size == expected_size
    format_result <- if (file_role == "replicate_bigwig_same_depth") {
        validate_bigwig(path)
    }
    else if (file_role == "processed_expression_matrix") {
        validate_gzip_text(path, "expression")
    }
    else {
        validate_gzip_text(path, "bed")
    }
    list(pass = size_pass && format_result$pass, size_pass = size_pass, format_pass = format_result$pass, detail = paste0("actual_size=", 
        actual_size, "; expected_size=", expected_size, "; ", format_result$detail))
}

download_one <- function(url, final_path, expected_size, file_role) {
    existing <- validate_download(final_path, expected_size, file_role)
    if (isTRUE(existing$pass)) {
        return(list(status = "VALID_EXISTING_FILE", attempts = 0L, validation = existing))
    }
    if (file.exists(final_path)) {
        renamed <- paste0(final_path, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        file.rename(final_path, renamed)
    }
    part_path <- paste0(final_path, ".part")
    if (file.exists(part_path)) {
        old_part <- paste0(part_path, ".previous_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        file.rename(part_path, old_part)
    }
    last_error <- ""
    for (attempt in seq_len(max_retries)) {
        log_event("Download attempt", attempt, "of", max_retries, basename(final_path))
        result <- tryCatch({
            handle <- curl::new_handle(followlocation = TRUE, useragent = paste0("Nutrition_Endometrial_VDR_v2/", 
                run_id), timeout = 7200, connecttimeout = 120)
            curl::curl_download(url = url, destfile = part_path, quiet = FALSE, mode = "wb", handle = handle)
            TRUE
        }, error = function(e) {
            last_error <<- conditionMessage(e)
            FALSE
        })
        if (isTRUE(result)) {
            validation <- validate_download(part_path, expected_size, file_role)
            if (isTRUE(validation$pass)) {
                ok <- file.rename(part_path, final_path)
                if (!ok) {
                  stop("Validated part file could not be renamed:\n", final_path)
                }
                return(list(status = "DOWNLOADED_AND_VALIDATED", attempts = attempt, validation = validation))
            }
            last_error <- paste0("validation failed: ", validation$detail)
        }
        if (file.exists(part_path)) {
            unlink(part_path, force = TRUE)
        }
        if (attempt < max_retries) {
            Sys.sleep(retry_wait_seconds * attempt)
        }
    }
    list(status = "FAILED", attempts = max_retries, validation = list(pass = FALSE, size_pass = FALSE, format_pass = FALSE, 
        detail = last_error))
}

remote_manifest_path <- file.path(dirs$results, "STEP03A_v4_processed_file_remote_size_manifest.csv")

gate_path <- file.path(dirs$results, "STEP03A_v4_GATE_TABLE.csv")

design_path <- file.path(dirs$results, "STEP03A_v4_sample_design_LOCKED.csv")

contrast_path <- file.path(dirs$results, "STEP03A_v4_contrast_plan_LOCKED.csv")

for (p in c(remote_manifest_path, gate_path, design_path, contrast_path)) {
    if (!file.exists(p)) {
        stop("Required STEP03A v4 input is missing:\n", p)
    }
}

gate <- .vdr_read_csv(gate_path, stringsAsFactors = FALSE, check.names = FALSE)

gate_value <- gate$status[gate$gate_component == "STEP03A v4 gate"]

if (length(gate_value) != 1L || gate_value != "READY_FOR_TARGETED_PROCESSED_DOWNLOAD") {
    stop("STEP03A v4 did not authorize targeted download. Gate=", paste(gate_value, collapse = ","))
}

manifest <- .vdr_read_csv(remote_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

selected <- manifest[manifest$remote_accessible & manifest$download_decision %in% c("APPROVE_TARGETED_DOWNLOAD", 
    "APPROVE_IF_SIZE_BUDGET_ACCEPTABLE"), , drop = FALSE]

if (!download_atac_bigwig) {
    selected <- selected[selected$file_role != "replicate_bigwig_same_depth", , drop = FALSE]
}

if (download_atac_bigwig && nrow(selected) != 57L) {
    stop("Expected 57 selected files but found ", nrow(selected), ".")
}

if (!download_atac_bigwig && nrow(selected) != 51L) {
    stop("Expected 51 core files but found ", nrow(selected), ".")
}

log_event("Starting", run_id, "selected files=", nrow(selected), "expected compressed GB=", round(sum(selected$size_bytes)/1024^3, 
    3))

rows <- vector("list", nrow(selected))

start_time <- Sys.time()

for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]
    subdir <- destination_subdir(row$accession, row$file_role)
    target_dir <- file.path(project_root, "03_data", "raw_external_READ_ONLY", row$accession, "processed_supplementary", 
        subdir)
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    final_path <- file.path(target_dir, row$file_name)
    log_event(i, "/", nrow(selected), row$file_id, row$file_name, "expected bytes=", row$size_bytes)
    one_start <- Sys.time()
    result <- download_one(url = row$url, final_path = final_path, expected_size = row$size_bytes, file_role = row$file_role)
    one_end <- Sys.time()
    final_exists <- file.exists(final_path)
    rows[[i]] <- data.frame(file_id = row$file_id, accession = row$accession, assay = row$assay, gsm = row$gsm, 
        file_role = row$file_role, file_name = row$file_name, source_url = row$url, destination = final_path, expected_size_bytes = row$size_bytes, 
        actual_size_bytes = if (final_exists) {
            as.numeric(file.info(final_path)$size)
        }
        else {
            NA_real_
        }, status = result$status, attempts = result$attempts, size_validation = result$validation$size_pass, format_validation = result$validation$format_pass, 
        validation_detail = result$validation$detail, md5 = if (final_exists) {
            md5_file(final_path)
        }
        else {
            NA_character_
        }, sha256 = if (final_exists) {
            sha256_file(final_path)
        }
        else {
            NA_character_
        }, runtime_seconds = as.numeric(difftime(one_end, one_start, units = "secs")), stringsAsFactors = FALSE)
    progress <- do.call(rbind, rows[seq_len(i)])
    write_csv_utf8(progress, progress_path)
}

download_manifest <- do.call(rbind, rows)

end_time <- Sys.time()

write_csv_utf8(download_manifest, file.path(dirs$results, "STEP03A_v5_targeted_download_manifest.csv"))

all_pass <- all(download_manifest$status %in% c("VALID_EXISTING_FILE", "DOWNLOADED_AND_VALIDATED") & download_manifest$size_validation & 
    download_manifest$format_validation)

size_summary <- aggregate(list(n_files = download_manifest$file_id, total_bytes = download_manifest$actual_size_bytes), 
    by = list(accession = download_manifest$accession, assay = download_manifest$assay, file_role = download_manifest$file_role, 
        status = download_manifest$status), FUN = function(x) {
        if (is.character(x)) {
            length(x)
        }
        else {
            sum(x, na.rm = TRUE)
        }
    })

size_summary$total_gb <- size_summary$total_bytes/1024^3

write_csv_utf8(size_summary, file.path(dirs$results, "STEP03A_v5_download_size_summary.csv"))

integrity_summary <- data.frame(selected_files = nrow(selected), valid_files = sum(download_manifest$size_validation & 
    download_manifest$format_validation), failed_files = sum(!download_manifest$size_validation | !download_manifest$format_validation), 
    expected_bytes = sum(selected$size_bytes), downloaded_bytes = sum(download_manifest$actual_size_bytes, na.rm = TRUE), 
    expected_gb = sum(selected$size_bytes)/1024^3, downloaded_gb = sum(download_manifest$actual_size_bytes, na.rm = TRUE)/1024^3, 
    runtime_minutes = as.numeric(difftime(end_time, start_time, units = "mins")), all_pass = all_pass, stringsAsFactors = FALSE)

write_csv_utf8(integrity_summary, file.path(dirs$results, "STEP03A_v5_integrity_summary.csv"))

automated_gate <- if (all_pass) {
    "READY_FOR_STEP03B_RNA_ANALYSIS"
} else {
    "HOLD_DOWNLOAD_OR_INTEGRITY_FAILURE"
}

gate_table <- data.frame(gate_component = c("STEP03A v4 design lock", "Selected processed files", "File-size validation", 
    "Format validation", "MD5 and SHA-256 capture", "RAW TAR/FASTQ/BAM", "STEP03A v5 gate"), status = c("PASS", 
    paste0("N=", nrow(selected)), ifelse(all(download_manifest$size_validation), "PASS", "FAIL"), ifelse(all(download_manifest$format_validation), 
        "PASS", "FAIL"), ifelse(all(nzchar(download_manifest$md5) & nzchar(download_manifest$sha256)), "PASS", 
        "FAIL"), "HOLD", automated_gate), required_action = c("Retain locked design and contrast plan.", "Retain exact URL and destination manifest.", 
    "Investigate any size mismatch.", "Investigate unreadable gzip/BED/BigWig.", "Use hashes for all later analysis inputs.", 
    "Do not download without a separate reprocessing decision.", "Begin RNA-seq QC and factorial analysis only after PASS."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03A_v5_GATE_TABLE.csv"))

report <- c("# STEP03A v5 targeted processed-data download report", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: `", 
    automated_gate, "`"), paste0("- Selected files: ", nrow(selected), "."), paste0("- Validated files: ", sum(download_manifest$size_validation & 
    download_manifest$format_validation), "."), paste0("- Downloaded/verified size: ", round(sum(download_manifest$actual_size_bytes, 
    na.rm = TRUE)/1024^3, 3), " GB."), paste0("- Runtime: ", round(as.numeric(difftime(end_time, start_time, units = "mins")), 
    1), " minutes."), "", "## Acquired evidence", "", "- GSE254251 RPKM expression matrix.", "- GSE306094 merged siNT and siVDR peak sets.", 
    if (download_atac_bigwig) {
        "- Six GSE306094 same-depth BigWig files."
    } else {
        "- Six GSE306094 BigWig files were intentionally deferred."
    }, "- Forty-eight GSE306127 same-depth CUT&RUN BED.gz files.", "", "## Preserved restrictions", "", "- DESeq2/edgeR count models are prohibited for the RPKM matrix.", 
    "- BigWig signal is normalized signal, not raw fragment counts.", "- Merged ATAC BED files do not provide replicate uncertainty.", 
    "- CUT&RUN BED files do not replace raw-read QC.", "- RAW TAR, FASTQ, BAM and SRA files remain on hold.", "", 
    "## Next step", "", if (all_pass) {
        "- Proceed to STEP03B1 RNA-seq matrix identity, QC and factorial limma analysis."
    } else {
        "- Do not begin analysis until every failed file is corrected and revalidated."
    }, "", "No scientific figure was generated in this download substep.")

write_text(report, file.path(dirs$results, "STEP03A_v5_TARGETED_DOWNLOAD_REPORT.md"))

write_text(c("# STEP03A v5 figure classification", "", "- Figure class: NONE", "- This substep performs data acquisition and integrity validation.", 
    "- Tables and manifests are authoritative.", "- No scientific or QC figure is required."), file.path(dirs$figures, 
    "STEP03A_v5_FIGURE_CLASSIFICATION.md"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "download_experimental_matrices.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03A_v5", 
    run_id = run_id, download_atac_bigwig = download_atac_bigwig, max_retries = max_retries, selected_file_count = nrow(selected), 
    expected_compressed_gb = sum(selected$size_bytes)/1024^3, raw_archive_policy = "HOLD", analysis_started = FALSE, 
    figure_class = "NONE")

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03A_v5_parameters_v0.1.json"), pretty = TRUE, 
    auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A_v5.txt"))

for (d in c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", "06_results", 
    "07_figures_tables")) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_groups <- list(`01_inputs` = c(remote_manifest_path, gate_path, design_path, contrast_path), `02_code` = list.files(dirs$code, 
    pattern = "07_step03A_v5_", full.names = TRUE), `03_parameters` = list.files(dirs$parameters, pattern = "^STEP03A_v5_", 
    full.names = TRUE), `04_environment` = list.files(dirs$environment, pattern = "STEP03A_v5", full.names = TRUE), 
    `05_execution` = log_path, `06_results` = list.files(dirs$results, pattern = "^STEP03A_v5_", full.names = TRUE), 
    `07_figures_tables` = file.path(dirs$figures, "STEP03A_v5_FIGURE_CLASSIFICATION.md"))

for (subdir in names(copy_groups)) {
    paths <- copy_groups[[subdir]]
    paths <- paths[file.exists(paths)]
    for (p in paths) {
        file.copy(p, file.path(dirs$archive, subdir, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03A_v5_archive_manifest.csv"))

write_text(c("# STEP03A v5 archive status", "", paste0("- Run ID: ", run_id), paste0("- Gate: ", automated_gate), 
    "- Design and contrast plan: LOCKED", "- Targeted processed download: COMPLETE", "- Data-integrity validation: COMPLETE", 
    "- RAW TAR/FASTQ/BAM: HOLD", "- Scientific analysis: NOT_STARTED", "- Figure production: NOT_REQUIRED", "- Full STEP03 module: NOT_LOCKED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03A_v5", paste0("- Run ID: ", run_id), paste0("- Gate: ", 
    automated_gate), "- Experimental design: LOCKED", "- Contrast plan: LOCKED", "- Processed data acquisition: COMPLETE", 
    "- RAW archives: HOLD", "- Next substep: STEP03B1 RNA-seq QC and factorial analysis", "- Full module status: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03A v5 targeted download completed.\n")

cat("Gate: ", automated_gate, "\n")

cat("Files selected: ", nrow(selected), "\n")

cat("Validated: ", sum(download_manifest$size_validation & download_manifest$format_validation), "\n")

cat("Total size: ", round(sum(download_manifest$actual_size_bytes, na.rm = TRUE)/1024^3, 3), " GB\n")

cat("RAW TAR/FASTQ/BAM remain on HOLD.\n")

cat("No scientific figure was generated.\n\n")

cat("Upload these files after completion:\n")

cat("  STEP03A_v5_TARGETED_DOWNLOAD_REPORT.md\n")

cat("  STEP03A_v5_targeted_download_manifest.csv\n")

cat("  STEP03A_v5_download_size_summary.csv\n")

cat("  STEP03A_v5_integrity_summary.csv\n")

cat("  STEP03A_v5_GATE_TABLE.csv\n")

cat("============================================================\n")

