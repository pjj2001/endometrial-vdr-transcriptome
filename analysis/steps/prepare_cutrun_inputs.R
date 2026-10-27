#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A07b: Prepare CUTRUN input inventory and BED sources

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3A_", format(Sys.time(), "%Y%m%d_%H%M%S"))

chunk_lines <- 100000L

width_sample_cap <- 100000L

keys_sample_cap <- 100000L

within_condition_row_count_cv_threshold <- 0.02

minimum_primary_fraction <- 0.94999999999999996

maximum_invalid_fraction <- 0.001

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

required_packages <- c("data.table", "ggplot2", "digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". Install them separately; this script does not update packages.")
}

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

qc_root <- file.path(dirs$figures, "A_QC_only", "STEP03B3A_CUTRUN_BED_audit")

for (d in c(file.path(qc_root, "figures"), file.path(qc_root, "source_data"), file.path(qc_root, "code"), file.path(qc_root, 
    "qc"), file.path(qc_root, "manifests"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

archive_root <- file.path(dirs$archive_parent, paste0("S3B3A_", format(Sys.time(), "%Y%m%d_%H%M%S")))

for (d in c(file.path(archive_root, "A"), file.path(archive_root, "res"), file.path(archive_root, "qc"), file.path(archive_root, 
    "txt"), file.path(archive_root, "code"), file.path(archive_root, "env"), file.path(archive_root, "log"), file.path(archive_root, 
    "ctl"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
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

gzip_magic_pass <- function(path) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    magic <- readBin(con, what = "raw", n = 2L)
    length(magic) == 2L && identical(as.integer(magic), c(31L, 139L))
}

add_named_counts <- function(existing, new_table) {
    if (length(new_table) == 0L) {
        return(existing)
    }
    if (length(existing) == 0L) {
        existing <- integer(0)
    }
    for (name_value in names(new_table)) {
        current_value <- existing[name_value]
        if (length(current_value) == 0L || is.na(current_value)) {
            current_value <- 0L
        }
        existing[name_value] <- as.integer(current_value) + as.integer(new_table[[name_value]])
    }
    existing
}

append_capped_sample <- function(existing, new_values, cap, per_chunk = 2500L) {
    new_values <- new_values[is.finite(new_values)]
    if (length(new_values) == 0L) {
        return(existing)
    }
    take_n <- min(per_chunk, length(new_values))
    selected <- if (length(new_values) > take_n) {
        sample(new_values, size = take_n, replace = FALSE)
    }
    else {
        new_values
    }
    combined <- c(existing, selected)
    if (length(combined) > cap) {
        combined <- sample(combined, size = cap, replace = FALSE)
    }
    combined
}

append_capped_character_sample <- function(existing, new_values, cap, per_chunk = 2500L) {
    new_values <- new_values[!is.na(new_values) & nzchar(new_values)]
    if (length(new_values) == 0L) {
        return(existing)
    }
    take_n <- min(per_chunk, length(new_values))
    selected <- if (length(new_values) > take_n) {
        sample(new_values, size = take_n, replace = FALSE)
    }
    else {
        new_values
    }
    combined <- c(existing, selected)
    if (length(combined) > cap) {
        combined <- sample(combined, size = cap, replace = FALSE)
    }
    combined
}

audit_bed_gz <- function(path, sample_seed) {
    set.seed(sample_seed)
    if (!file.exists(path)) {
        stop("BED.gz file does not exist:\n", path)
    }
    if (!gzip_magic_pass(path)) {
        stop("File does not have a gzip header:\n", path)
    }
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    total_rows <- 0L
    valid_rows <- 0L
    invalid_rows <- 0L
    total_width <- 0
    width_min <- Inf
    width_max <- -Inf
    column_count_min <- Inf
    column_count_max <- -Inf
    chromosome_counts <- integer(0)
    width_sample <- numeric(0)
    key_sample <- character(0)
    first_chunk <- TRUE
    repeat {
        lines <- readLines(con, n = chunk_lines, warn = FALSE)
        if (length(lines) == 0L) {
            break
        }
        lines <- lines[nzchar(lines) & !grepl("^#|^track\\b|^browser\\b", lines, ignore.case = TRUE, perl = TRUE)]
        if (length(lines) == 0L) {
            next
        }
        total_rows <- total_rows + length(lines)
        if (first_chunk) {
            inspect_n <- min(10000L, length(lines))
            observed_columns <- lengths(strsplit(lines[seq_len(inspect_n)], "\t", fixed = TRUE))
            column_count_min <- min(observed_columns)
            column_count_max <- max(observed_columns)
            first_chunk <- FALSE
        }
        fields <- data.table::tstrsplit(lines, "\t", fixed = TRUE, keep = 1:3)
        if (length(fields) < 3L) {
            invalid_rows <- invalid_rows + length(lines)
            next
        }
        chromosome <- as.character(fields[[1]])
        start0 <- suppressWarnings(as.numeric(fields[[2]]))
        end_value <- suppressWarnings(as.numeric(fields[[3]]))
        valid <- !is.na(chromosome) & nzchar(chromosome) & is.finite(start0) & is.finite(end_value) & start0 >= 
            0 & end_value > start0
        valid_rows <- valid_rows + sum(valid)
        invalid_rows <- invalid_rows + sum(!valid)
        if (!any(valid)) {
            next
        }
        valid_chromosome <- chromosome[valid]
        valid_start <- start0[valid]
        valid_end <- end_value[valid]
        widths <- valid_end - valid_start
        chromosome_counts <- add_named_counts(chromosome_counts, table(valid_chromosome))
        total_width <- total_width + sum(widths)
        width_min <- min(width_min, min(widths))
        width_max <- max(width_max, max(widths))
        width_sample <- append_capped_sample(width_sample, widths, cap = width_sample_cap)
        keys <- paste0(valid_chromosome, ":", format(valid_start, scientific = FALSE, trim = TRUE), "-", format(valid_end, 
            scientific = FALSE, trim = TRUE))
        key_sample <- append_capped_character_sample(key_sample, keys, cap = keys_sample_cap)
    }
    if (valid_rows == 0L) {
        stop("No valid intervals were found in:\n", path)
    }
    observed_chromosomes <- names(chromosome_counts)
    has_chr_prefix <- if (length(observed_chromosomes) > 0L) {
        sum(chromosome_counts[grepl("^chr", observed_chromosomes)])/sum(chromosome_counts) >= 0.5
    }
    else {
        FALSE
    }
    primary_chromosomes <- if (has_chr_prefix) {
        paste0("chr", c(1:22, "X", "Y"))
    }
    else {
        as.character(c(1:22, "X", "Y"))
    }
    primary_rows <- sum(chromosome_counts[names(chromosome_counts) %in% primary_chromosomes])
    nonprimary_rows <- valid_rows - primary_rows
    width_median_approx <- stats::median(width_sample)
    width_p95_approx <- as.numeric(stats::quantile(width_sample, probs = 0.94999999999999996, names = FALSE, type = 7))
    approximate_duplicate_fraction <- if (length(key_sample) > 0L) {
        1 - length(unique(key_sample))/length(key_sample)
    }
    else {
        NA_real_
    }
    provisional_semantics <- if (valid_rows >= 100000L && width_median_approx <= 1000) {
        "FRAGMENT_LIKE_INTERVALS"
    }
    else {
        "PEAK_LIKE_OR_MIXED_INTERVALS_REVIEW"
    }
    list(total_rows = total_rows, valid_rows = valid_rows, invalid_rows = invalid_rows, invalid_fraction = invalid_rows/max(total_rows, 
        1L), primary_rows = primary_rows, nonprimary_rows = nonprimary_rows, primary_fraction = primary_rows/valid_rows, 
        chromosome_n = length(chromosome_counts), chromosome_prefix = ifelse(has_chr_prefix, "chr_prefixed", "not_chr_prefixed"), 
        column_count_min = ifelse(is.finite(column_count_min), column_count_min, NA_real_), column_count_max = ifelse(is.finite(column_count_max), 
            column_count_max, NA_real_), total_interval_bp = total_width, mean_width = total_width/valid_rows, 
        minimum_width = ifelse(is.finite(width_min), width_min, NA_real_), maximum_width = ifelse(is.finite(width_max), 
            width_max, NA_real_), median_width_approx = width_median_approx, p95_width_approx = width_p95_approx, 
        interval_sample_n = length(width_sample), coordinate_key_sample_n = length(key_sample), approximate_duplicate_fraction = approximate_duplicate_fraction, 
        provisional_record_semantics = provisional_semantics)
}

save_qc_png <- function(plot_object, path, width_mm, height_mm) {
    ggplot2::ggsave(filename = path, plot = plot_object, width = width_mm, height = height_mm, units = "mm", dpi = 220, 
        bg = "white")
}

log_event("Starting", run_id)

required_paths <- c(step03b2a2_gate = file.path(dirs$results, "STEP03B2A2_GATE_TABLE.csv"), cutrun_authorization = file.path(dirs$control, 
    "STEP03B3_CUTRUN_AUTHORIZATION.md"), locked_design = file.path(dirs$results, "STEP03A_v4_sample_design_LOCKED.csv"), 
    locked_contrast_plan = file.path(dirs$results, "STEP03A_v4_contrast_plan_LOCKED.csv"), download_manifest = file.path(dirs$results, 
        "STEP03A_v5_targeted_download_manifest.csv"), RNA_lock_card = file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3A input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

step03b2a2_gate <- .vdr_read_csv(required_paths[["step03b2a2_gate"]], stringsAsFactors = FALSE, check.names = FALSE)

authorization_text <- readLines(required_paths[["cutrun_authorization"]], warn = FALSE, encoding = "UTF-8")

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- step03b2a2_gate$status[step03b2a2_gate$gate_component == "STEP03B2A2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "READY_FOR_STEP03B3_CUTRUN_WITH_RESTRICTED_ATAC_SUPPORT") {
    stop("STEP03B2A2 gate does not authorize CUT&RUN. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (!any(grepl("Authorization: YES", authorization_text, fixed = TRUE)) || !any(grepl("Primary chromatin mechanism layer: CUT&RUN", 
    authorization_text, fixed = TRUE))) {
    stop("STEP03B3 CUT&RUN authorization is incomplete.")
}

if (!any(grepl("RNA_PROGRAM_v1.0", RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", 
    RNA_lock_text, fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

locked_design <- .vdr_read_csv(required_paths[["locked_design"]], stringsAsFactors = FALSE, check.names = FALSE)

contrast_plan <- .vdr_read_csv(required_paths[["locked_contrast_plan"]], stringsAsFactors = FALSE, check.names = FALSE)

download_manifest <- .vdr_read_csv(required_paths[["download_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

cutrun_design <- locked_design[locked_design$accession == "GSE306127", , drop = FALSE]

required_design_columns <- c("gsm", "experiment_block", "genetic_context", "ligand_status_locked", "assay_target_locked", 
    "replicate_locked", "condition_id_locked")

missing_design_columns <- setdiff(required_design_columns, names(cutrun_design))

if (length(missing_design_columns) > 0L) {
    stop("CUT&RUN design is missing columns: ", paste(missing_design_columns, collapse = ", "))
}

if (nrow(cutrun_design) != 48L || sum(cutrun_design$experiment_block == "CUTRUN_VDR_cistrome") != 24L || sum(cutrun_design$experiment_block == 
    "CUTRUN_histone_marks") != 24L) {
    stop("CUT&RUN design must contain 48 samples split 24/24 between blocks.")
}

condition_counts <- aggregate(gsm ~ condition_id_locked, data = cutrun_design, FUN = length)

if (nrow(condition_counts) != 16L || any(condition_counts$gsm != 3L)) {
    stop("CUT&RUN design must contain 16 locked condition cells with n=3 each.")
}

replicate_complete <- vapply(split(cutrun_design$replicate_locked, cutrun_design$condition_id_locked), function(x) {
    identical(sort(as.integer(x)), 1:3)
}, logical(1))

if (!all(replicate_complete)) {
    stop("At least one CUT&RUN condition does not contain replicates 1, 2 and 3.")
}

cutrun_downloads <- download_manifest[download_manifest$accession == "GSE306127" & download_manifest$file_role == 
    "replicate_BED_same_depth", , drop = FALSE]

if (nrow(cutrun_downloads) != 48L) {
    stop("Expected 48 downloaded CUT&RUN BED.gz files; found ", nrow(cutrun_downloads), ".")
}

cutrun_manifest <- merge(cutrun_design, cutrun_downloads, by = "gsm", all.x = TRUE, sort = FALSE)

cutrun_manifest <- cutrun_manifest[match(cutrun_design$gsm, cutrun_manifest$gsm), , drop = FALSE]

required_download_columns <- c("destination", "file_name", "md5", "sha256")

missing_download_columns <- setdiff(required_download_columns, names(cutrun_manifest))

if (length(missing_download_columns) > 0L) {
    stop("CUT&RUN download manifest is missing columns: ", paste(missing_download_columns, collapse = ", "))
}

cutrun_manifest$file_exists <- file.exists(cutrun_manifest$destination)

cutrun_manifest$observed_md5 <- vapply(cutrun_manifest$destination, md5_file, character(1))

cutrun_manifest$observed_sha256 <- vapply(cutrun_manifest$destination, sha256_file, character(1))

cutrun_manifest$hash_pass <- with(cutrun_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

write_csv_utf8(cutrun_manifest, file.path(dirs$inputs, "STEP03B3A_CUTRUN_input_manifest_PRE_RUN_LOCK.csv"))

if (!all(cutrun_manifest$hash_pass)) {
    stop("One or more CUT&RUN BED.gz files failed existence or hash validation.")
}

cutrun_contrasts <- contrast_plan[grepl("^CUTRUN_", contrast_plan$contrast_id), , drop = FALSE]

write_csv_utf8(cutrun_contrasts, file.path(dirs$results, "STEP03B3A_CUTRUN_contrast_plan_LOCKED.csv"))

target_conditions <- unique(cutrun_design[cutrun_design$assay_target_locked != "IgG_control", c("experiment_block", 
    "genetic_context", "ligand_status_locked", "assay_target_locked", "condition_id_locked"), drop = FALSE])

pairing_rows <- vector("list", nrow(target_conditions))

for (i in seq_len(nrow(target_conditions))) {
    target_row <- target_conditions[i, , drop = FALSE]
    if (target_row$experiment_block == "CUTRUN_VDR_cistrome") {
        control_match <- cutrun_design[cutrun_design$experiment_block == target_row$experiment_block & cutrun_design$genetic_context == 
            target_row$genetic_context & cutrun_design$ligand_status_locked == target_row$ligand_status_locked & 
            cutrun_design$assay_target_locked == "IgG_control", , drop = FALSE]
    }
    else {
        control_match <- cutrun_design[cutrun_design$experiment_block == target_row$experiment_block & cutrun_design$genetic_context == 
            target_row$genetic_context & cutrun_design$assay_target_locked == "IgG_control", , drop = FALSE]
    }
    target_match <- cutrun_design[cutrun_design$condition_id_locked == target_row$condition_id_locked, , drop = FALSE]
    pairing_rows[[i]] <- data.frame(pairing_id = paste0("CUTRUN_PAIR_", sprintf("%02d", i)), experiment_block = target_row$experiment_block, 
        genetic_context = target_row$genetic_context, ligand_status = target_row$ligand_status_locked, target = target_row$assay_target_locked, 
        target_condition_id = target_row$condition_id_locked, target_gsms = paste(target_match$gsm, collapse = ";"), 
        IgG_condition_id = paste(unique(control_match$condition_id_locked), collapse = ";"), IgG_gsms = paste(control_match$gsm, 
            collapse = ";"), target_n = nrow(target_match), IgG_n = nrow(control_match), pairing_pass = (nrow(target_match) == 
            3L && nrow(control_match) == 3L), stringsAsFactors = FALSE)
}

matched_control_map <- do.call(rbind, pairing_rows)

write_csv_utf8(matched_control_map, file.path(dirs$results, "STEP03B3A_MATCHED_IgG_CONTROL_MAP_LOCKED.csv"))

if (!all(matched_control_map$pairing_pass)) {
    stop("One or more CUT&RUN target conditions lack a matched n=3 IgG control.")
}

audit_rows <- vector("list", nrow(cutrun_manifest))

for (i in seq_len(nrow(cutrun_manifest))) {
    log_event("Auditing CUT&RUN BED", i, "/", nrow(cutrun_manifest), cutrun_manifest$gsm[i], cutrun_manifest$file_name[i])
    audit_result <- tryCatch(audit_bed_gz(cutrun_manifest$destination[i], sample_seed = 20260803L + i), error = function(e) {
        list(audit_error = conditionMessage(e))
    })
    if (!is.null(audit_result$audit_error)) {
        audit_rows[[i]] <- data.frame(gsm = cutrun_manifest$gsm[i], audit_status = "FAIL", audit_error = audit_result$audit_error, 
            stringsAsFactors = FALSE)
        next
    }
    audit_rows[[i]] <- data.frame(gsm = cutrun_manifest$gsm[i], audit_status = "PASS", audit_error = "", total_rows = audit_result$total_rows, 
        valid_rows = audit_result$valid_rows, invalid_rows = audit_result$invalid_rows, invalid_fraction = audit_result$invalid_fraction, 
        primary_rows = audit_result$primary_rows, nonprimary_rows = audit_result$nonprimary_rows, primary_fraction = audit_result$primary_fraction, 
        chromosome_n = audit_result$chromosome_n, chromosome_prefix = audit_result$chromosome_prefix, column_count_min = audit_result$column_count_min, 
        column_count_max = audit_result$column_count_max, total_interval_bp = audit_result$total_interval_bp, mean_width = audit_result$mean_width, 
        minimum_width = audit_result$minimum_width, maximum_width = audit_result$maximum_width, median_width_approx = audit_result$median_width_approx, 
        p95_width_approx = audit_result$p95_width_approx, interval_sample_n = audit_result$interval_sample_n, coordinate_key_sample_n = audit_result$coordinate_key_sample_n, 
        approximate_duplicate_fraction = audit_result$approximate_duplicate_fraction, provisional_record_semantics = audit_result$provisional_record_semantics, 
        stringsAsFactors = FALSE)
}

bed_audit <- data.table::rbindlist(audit_rows, fill = TRUE)

bed_audit <- as.data.frame(bed_audit, stringsAsFactors = FALSE)

bed_audit <- merge(cutrun_manifest[c("gsm", "file_name", "destination", "experiment_block", "genetic_context", 
    "ligand_status_locked", "assay_target_locked", "replicate_locked", "condition_id_locked", "md5", "sha256")], 
    bed_audit, by = "gsm", all.x = TRUE, sort = FALSE)

bed_audit <- bed_audit[match(cutrun_manifest$gsm, bed_audit$gsm), , drop = FALSE]

write_csv_utf8(bed_audit, file.path(dirs$results, "STEP03B3A_CUTRUN_BED_integrity_audit.csv"))

parse_pass <- all(bed_audit$audit_status == "PASS")

if (!parse_pass) {
    gate_status <- "HOLD_CUTRUN_BED_PARSE_FAILURE"
} else {
    condition_depth_summary <- aggregate(valid_rows ~ condition_id_locked + experiment_block + genetic_context + 
        ligand_status_locked + assay_target_locked, data = bed_audit, FUN = function(x) {
        c(n = length(x), minimum = min(x), median = stats::median(x), maximum = max(x), mean = mean(x), sd = stats::sd(x), 
            cv = ifelse(mean(x) == 0, NA_real_, stats::sd(x)/mean(x)))
    })
    condition_depth_summary <- data.frame(condition_id_locked = condition_depth_summary$condition_id_locked, experiment_block = condition_depth_summary$experiment_block, 
        genetic_context = condition_depth_summary$genetic_context, ligand_status_locked = condition_depth_summary$ligand_status_locked, 
        assay_target_locked = condition_depth_summary$assay_target_locked, n = condition_depth_summary$valid_rows[, 
            "n"], minimum_rows = condition_depth_summary$valid_rows[, "minimum"], median_rows = condition_depth_summary$valid_rows[, 
            "median"], maximum_rows = condition_depth_summary$valid_rows[, "maximum"], mean_rows = condition_depth_summary$valid_rows[, 
            "mean"], sd_rows = condition_depth_summary$valid_rows[, "sd"], cv_rows = condition_depth_summary$valid_rows[, 
            "cv"], stringsAsFactors = FALSE)
    write_csv_utf8(condition_depth_summary, file.path(dirs$qc, "STEP03B3A_condition_depth_summary.csv"))
    integrity_pass <- all(bed_audit$invalid_fraction <= maximum_invalid_fraction) && all(bed_audit$primary_fraction >= 
        minimum_primary_fraction)
    depth_pass <- all(is.finite(condition_depth_summary$cv_rows) & condition_depth_summary$cv_rows <= within_condition_row_count_cv_threshold)
    semantics_table <- as.data.frame(table(bed_audit$provisional_record_semantics), stringsAsFactors = FALSE)
    names(semantics_table) <- c("provisional_record_semantics", "n_files")
    write_csv_utf8(semantics_table, file.path(dirs$qc, "STEP03B3A_record_semantics_summary.csv"))
    semantics_pass <- all(bed_audit$provisional_record_semantics == "FRAGMENT_LIKE_INTERVALS")
    gate_status <- if (integrity_pass && depth_pass && semantics_pass) {
        "PROCEED_TO_MANUAL_CUTRUN_BED_QC"
    }
    else if (integrity_pass && depth_pass) {
        "HOLD_RECORD_SEMANTICS_REVIEW"
    }
    else {
        "HOLD_CUTRUN_BED_INTEGRITY_OR_DEPTH_QC"
    }
}

if (parse_pass) {
    plot_data <- bed_audit
    plot_data$short_label <- paste0(ifelse(plot_data$experiment_block == "CUTRUN_VDR_cistrome", "CIS", "HIS"), 
        "_", plot_data$assay_target_locked, "_", plot_data$replicate_locked)
    plot_data$sample_order <- seq_len(nrow(plot_data))
    write_csv_utf8(plot_data, file.path(qc_root, "source_data", "SourceData_STEP03B3A_BED_audit_samples.csv"))
    p_rows <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_order, y = valid_rows, color = assay_target_locked, 
        shape = experiment_block)) + ggplot2::geom_point(size = 2.2000000000000002) + ggplot2::facet_wrap(~experiment_block, 
        scales = "free_x", ncol = 1) + ggplot2::labs(title = "CUT&RUN processed BED interval counts", subtitle = "Forty-eight same-depth files; points are biological replicates", 
        x = "Locked sample order", y = "Valid interval rows", color = "Target", shape = "Block") + ggplot2::theme_bw(base_size = 10)
    save_qc_png(p_rows, file.path(qc_root, "figures", "QC_STEP03B3A_interval_counts.png"), width_mm = 165, height_mm = 118)
    p_width <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_order, y = median_width_approx, color = assay_target_locked, 
        shape = experiment_block)) + ggplot2::geom_point(size = 2.2000000000000002) + ggplot2::facet_wrap(~experiment_block, 
        scales = "free_x", ncol = 1) + ggplot2::labs(title = "CUT&RUN interval-width audit", subtitle = "Approximate median from a deterministic capped interval sample", 
        x = "Locked sample order", y = "Approximate median interval width", color = "Target", shape = "Block") + 
        ggplot2::theme_bw(base_size = 10)
    save_qc_png(p_width, file.path(qc_root, "figures", "QC_STEP03B3A_interval_widths.png"), width_mm = 165, height_mm = 118)
    p_primary <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_order, y = primary_fraction, color = assay_target_locked, 
        shape = experiment_block)) + ggplot2::geom_hline(yintercept = minimum_primary_fraction, linetype = 2) + 
        ggplot2::geom_point(size = 2.2000000000000002) + ggplot2::facet_wrap(~experiment_block, scales = "free_x", 
        ncol = 1) + ggplot2::coord_cartesian(ylim = c(0, 1)) + ggplot2::labs(title = "CUT&RUN primary-chromosome fraction", 
        subtitle = "Dashed line marks the prespecified 0.95 review threshold", x = "Locked sample order", y = "Primary-chromosome fraction", 
        color = "Target", shape = "Block") + ggplot2::theme_bw(base_size = 10)
    save_qc_png(p_primary, file.path(qc_root, "figures", "QC_STEP03B3A_primary_chromosome_fraction.png"), width_mm = 165, 
        height_mm = 118)
}

gate_table <- data.frame(gate_component = c("STEP03B2A2 CUT&RUN authorization", "RNA_PROGRAM_v1.0", "ATAC evidence role", 
    "GSE306127 design", "Sixteen condition cells n=3", "Matched IgG control routing", "Forty-eight file/hash identities", 
    "Forty-eight BED.gz parse audits", "Within-condition same-depth consistency", "Record semantics", "Figure class", 
    "CUT&RUN occupancy analysis", "STEP03B3A gate"), status = c("PASS", "IMMUTABLE", "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", 
    "PASS_48_SAMPLES_24_PLUS_24", "PASS_16_OF_16", ifelse(all(matched_control_map$pairing_pass), "PASS", "FAIL"), 
    ifelse(all(cutrun_manifest$hash_pass), "PASS_48_OF_48", "FAIL"), ifelse(parse_pass, "PASS_48_OF_48", "FAIL"), 
    if (parse_pass) {
        ifelse(all(condition_depth_summary$cv_rows <= within_condition_row_count_cv_threshold), "PASS", "REVIEW")
    } else {
        "NOT_EVALUABLE"
    }, if (parse_pass) {
        paste(unique(bed_audit$provisional_record_semantics), collapse = ";")
    } else {
        "NOT_EVALUABLE"
    }, "A_QC_ONLY", "NOT_STARTED", gate_status), required_action = c("Retain authorization record.", "Do not change RNA membership.", 
    "Use ATAC only as accessible context.", "Retain locked sample design.", "Do not collapse or relabel condition cells.", 
    "Use matched IgG in all occupancy interpretation.", "Use only hash-validated BED files.", "Review any parse failure before analysis.", 
    "Review any condition with CV above 0.02.", "Confirm whether records are fragments or called peaks before modeling.", 
    "Do not promote current QC figures.", "Begin only after manual STEP03B3A review.", "Upload audit tables and actual A-class renders."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3A_GATE_TABLE.csv"))

write_text(c("# STEP03B3A CUT&RUN PRE_RUN_LOCK", "", paste0("- Run ID: ", run_id), "- Dataset: GSE306127", "- Samples: 48", 
    "- VDR-cistrome block: 24", "- Histone-mark block: 24", "- Locked condition cells: 16, n=3 each", "- Matched IgG routing: REQUIRED", 
    "- RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE", "- ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", "- Raw-read QC: unavailable from processed BED files", 
    "- Occupancy analysis: NOT STARTED", paste0("- Gate: ", gate_status), "", "No VDR occupancy, histone-mark direction, direct-target or", 
    "decidualization claim is authorized at STEP03B3A."), file.path(dirs$results, "STEP03B3A_PRE_RUN_LOCK_CARD.md"))

write_text(c("# STEP03B3A figure classification", "", "- Figure class: A_QC_ONLY.", "- Interval-count, width and primary-chromosome plots are internal diagnostics.", 
    "- No B-class CUT&RUN candidate is created.", "- Scientific Figure Skill 3 final-release requirements are not triggered.", 
    "- Direct source data and exact code are retained."), file.path(qc_root, "qc", "STEP03B3A_FIGURE_CLASSIFICATION.md"))

if (parse_pass) {
    figure_source_map <- data.frame(figure_file = c("figures/QC_STEP03B3A_interval_counts.png", "figures/QC_STEP03B3A_interval_widths.png", 
        "figures/QC_STEP03B3A_primary_chromosome_fraction.png"), direct_source_data = "source_data/SourceData_STEP03B3A_BED_audit_samples.csv", 
        generator_code = "code/31_step03B3A_CUTRUN_PRE_RUN_LOCK_and_BED_audit.R", figure_class = "A_QC_ONLY", release_status = "INTERNAL_QC_ONLY", 
        stringsAsFactors = FALSE)
    write_csv_utf8(figure_source_map, file.path(qc_root, "manifests", "STEP03B3A_FIGURE_SOURCE_MAP.csv"))
    write_text(c("# STEP03B3A source-data guide", "", "The three QC figures are generated from:", "`SourceData_STEP03B3A_BED_audit_samples.csv`.", 
        "", "The authoritative full audit table remains:", "`06_results/STEP03B3A_CUTRUN_BED_integrity_audit.csv`.", 
        "", "No BED file needs to be rescanned merely to re-render these QC plots."), file.path(qc_root, "source_data", 
        "README_SOURCE_DATA.md"))
}

write_text(c("# STEP03B3A Methods draft", "", "The 48 GSE306127 downsampled CUT&RUN BED files were reconciled to the locked two-block experimental design and verified against the targeted-download manifest using MD5 and SHA-256 hashes. The VDR-cistrome block comprised vector-control and VDR-overexpression cells under vehicle or calcitriol conditions, with VDR and matched IgG profiles in three biological replicates per cell. The histone-mark block comprised siRNA-control and VDR-knockdown cells profiled for H3K4me3, H3K27ac, H3K27me3 and matched IgG in three biological replicates per cell.", 
    "", "Every gzip-compressed BED file was audited in a streaming pass. The audit recorded interval-row counts, coordinate validity, chromosome composition, total interval coverage, approximate interval-width quantiles and approximate duplicate-coordinate burden. This substep was restricted to data integrity and experimental-control routing; no occupancy model, differential-binding test or peak-to-gene annotation was performed."), 
    file.path(dirs$text, "STEP03B3A_Methods_draft.md"))

write_text(c("# STEP03B3A Limitations", "", "- The files are processed/downsampled BED intervals rather than raw FASTQ/BAM.", 
    "- Original-depth library complexity, mapping quality and fragment-length QC cannot be reconstructed.", "- Approximate width and duplicate metrics use capped deterministic samples.", 
    "- Record semantics must be confirmed before choosing an occupancy model.", "- Matched IgG controls are mandatory.", 
    "- No decidualization treatment is present.", "- No direct-target claim is permitted at this stage."), file.path(dirs$text, 
    "STEP03B3A_Limitations.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, dataset = "GSE306127", samples = 48, VDR_cistrome_samples = 24, 
    histone_mark_samples = 24, condition_cells = 16, replicates_per_cell = 3, chunk_lines = chunk_lines, width_sample_cap = width_sample_cap, 
    coordinate_key_sample_cap = keys_sample_cap, maximum_invalid_fraction = maximum_invalid_fraction, minimum_primary_fraction = minimum_primary_fraction, 
    within_condition_row_count_cv_threshold = within_condition_row_count_cv_threshold, occupancy_analysis = FALSE, 
    RNA_integration = FALSE, figure_class = "A_QC_ONLY", gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3A_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3A.txt"))

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

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    safe_file_copy(source_file, .vdr_script_path("31_step03B3A_CUTRUN_PRE_RUN_LOCK_and_BED_audit.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "prepare_cutrun_inputs.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
    safe_file_copy(source_file, file.path(qc_root, "code", "prepare_cutrun_inputs.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

qc_files <- list.files(qc_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

qc_manifest <- data.frame(relative_path = substring(normalize_slashes(qc_files), nchar(normalize_slashes(qc_root)) + 
    2L), size_bytes = as.numeric(file.info(qc_files)$size), md5 = vapply(qc_files, md5_file, character(1)), sha256 = vapply(qc_files, 
    sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(qc_manifest, file.path(qc_root, "manifests", "STEP03B3A_QC_PACKAGE_MANIFEST.csv"))

qc_files <- list.files(qc_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

for (source_path in qc_files) {
    relative_path <- substring(normalize_slashes(source_path), nchar(normalize_slashes(qc_root)) + 2L)
    copy_validate(source_path, file.path(archive_root, "A", relative_path))
}

archive_artifacts <- c(file.path(dirs$inputs, "STEP03B3A_CUTRUN_input_manifest_PRE_RUN_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3A_CUTRUN_contrast_plan_LOCKED.csv"), file.path(dirs$results, "STEP03B3A_MATCHED_IgG_CONTROL_MAP_LOCKED.csv"), 
    file.path(dirs$results, "STEP03B3A_CUTRUN_BED_integrity_audit.csv"), file.path(dirs$results, "STEP03B3A_PRE_RUN_LOCK_CARD.md"), 
    file.path(dirs$results, "STEP03B3A_GATE_TABLE.csv"), file.path(dirs$qc, "STEP03B3A_condition_depth_summary.csv"), 
    file.path(dirs$qc, "STEP03B3A_record_semantics_summary.csv"), file.path(dirs$text, "STEP03B3A_Methods_draft.md"), 
    file.path(dirs$text, "STEP03B3A_Limitations.md"), file.path(dirs$parameters, "STEP03B3A_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3A.txt"), file.path(dirs$code, "prepare_cutrun_inputs.R"), 
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
    else if (grepl("Methods|Limitations", basename(path_value))) {
        "txt"
    }
    else if (grepl("summary|audit|LOCK_CARD", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3A_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3A archive status", "", paste0("- Run ID: ", run_id), paste0("- Gate: ", gate_status), "- Input BED hashes: 48/48 verified", 
    paste0("- BED parse audits: ", sum(bed_audit$audit_status == "PASS", na.rm = TRUE), "/48"), "- Matched IgG routing: LOCKED", 
    "- Occupancy analysis: NOT_STARTED", "- Figure class: A_QC_ONLY", paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3A", paste0("- Run ID: ", run_id), "- RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE", 
    "- ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", "- CUT&RUN design and IgG routing: LOCKED", paste0("- CUT&RUN BED gate: ", 
        gate_status), "- CUT&RUN occupancy analysis: NOT_STARTED", "- CUT&RUN retained figures: A_QC_ONLY", "- Full STEP03 module: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3A CUT&RUN PRE_RUN_LOCK and BED audit completed.\n")

cat("CUT&RUN samples: 48/48\n")

cat("Condition cells: 16; n=3 each\n")

cat("Matched IgG routing: LOCKED\n")

cat("BED parse audits passed: ", sum(bed_audit$audit_status == "PASS", na.rm = TRUE), "/48\n")

cat("RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE\n")

cat("ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY\n")

cat("Occupancy analysis: NOT STARTED\n")

cat("Figure class: A_QC_ONLY\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload after completion:\n")

cat("  06_results/STEP03B3A_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3A_CUTRUN_BED_integrity_audit.csv\n")

cat("  06_results/STEP03B3A_MATCHED_IgG_CONTROL_MAP_LOCKED.csv\n")

cat("  08_qc/STEP03B3A_condition_depth_summary.csv\n")

cat("  08_qc/STEP03B3A_record_semantics_summary.csv\n")

cat("  07_figures_tables/A_QC_only/STEP03B3A_CUTRUN_BED_audit/\n")

cat("============================================================\n")

