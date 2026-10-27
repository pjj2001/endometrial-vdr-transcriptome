#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A08: Read published VDR occupancy tables

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3B_", format(Sys.time(), "%Y%m%d_%H%M%S"))

minimum_coordinate_valid_fraction <- 0.98999999999999999

minimum_source_overlap_fraction <- 0.94999999999999996

required_packages <- c("readxl", "data.table", "digest", "jsonlite", "GenomicRanges", "IRanges")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

source_root <- file.path(dirs$results, "STEP03B3B_published_peak_source_tables")

dir.create(source_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3B_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

read_sheet_all_text <- function(workbook_path, sheet_name, col_names = TRUE) {
    readxl::read_excel(path = workbook_path, sheet = sheet_name, col_names = col_names, col_types = "text", .name_repair = "unique_quiet")
}

find_column_case_insensitive <- function(data_frame, target_name) {
    matching <- which(tolower(names(data_frame)) == tolower(target_name))
    if (length(matching) != 1L) {
        stop("Expected exactly one column named '", target_name, "'; found ", length(matching), ". Available columns: ", 
            paste(names(data_frame), collapse = " | "))
    }
    names(data_frame)[matching]
}

parse_coordinate_text <- function(x) {
    cleaned <- trimws(as.character(x))
    cleaned <- gsub(",", "", cleaned, fixed = TRUE)
    suppressWarnings(as.numeric(cleaned))
}

is_primary_chromosome <- function(chromosome) {
    grepl("^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$", chromosome, ignore.case = FALSE, perl = TRUE)
}

is_standard_chromosome <- function(chromosome) {
    grepl("^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", chromosome, ignore.case = FALSE, perl = TRUE)
}

standardize_coordinate_sheet <- function(data_frame, chrom_column, start_column, end_column, source_id) {
    raw_names <- names(data_frame)
    names(data_frame) <- make.unique(paste0("raw_", raw_names), sep = "__")
    raw_chrom_name <- paste0("raw_", chrom_column)
    raw_start_name <- paste0("raw_", start_column)
    raw_end_name <- paste0("raw_", end_column)
    if (!all(c(raw_chrom_name, raw_start_name, raw_end_name) %in% names(data_frame))) {
        stop("Coordinate columns were not preserved after renaming for ", source_id, ".")
    }
    chromosome <- trimws(as.character(data_frame[[raw_chrom_name]]))
    start0 <- parse_coordinate_text(data_frame[[raw_start_name]])
    end_value <- parse_coordinate_text(data_frame[[raw_end_name]])
    valid_coordinate <- nzchar(chromosome) & is_standard_chromosome(chromosome) & is.finite(start0) & is.finite(end_value) & 
        start0 >= 0 & end_value > start0
    standardized <- cbind(data.frame(source_id = source_id, source_row = seq_len(nrow(data_frame)) + 1L, chromosome = chromosome, 
        start0 = start0, end = end_value, width = end_value - start0, valid_coordinate = valid_coordinate, primary_chromosome = is_primary_chromosome(chromosome), 
        coordinate_key = ifelse(valid_coordinate, paste0(chromosome, ":", format(start0, scientific = FALSE, trim = TRUE), 
            "-", format(end_value, scientific = FALSE, trim = TRUE)), NA_character_), stringsAsFactors = FALSE), 
        data_frame)
    standardized
}

standardize_headerless_overlap_sheet <- function(data_frame, source_id) {
    if (ncol(data_frame) < 4L) {
        stop(source_id, " must contain at least four columns.")
    }
    names(data_frame) <- paste0("raw_col_", seq_len(ncol(data_frame)))
    chromosome <- trimws(as.character(data_frame[["raw_col_2"]]))
    start0 <- parse_coordinate_text(data_frame[["raw_col_3"]])
    end_value <- parse_coordinate_text(data_frame[["raw_col_4"]])
    valid_coordinate <- nzchar(chromosome) & is_standard_chromosome(chromosome) & is.finite(start0) & is.finite(end_value) & 
        start0 >= 0 & end_value > start0
    cbind(data.frame(source_id = source_id, source_row = seq_len(nrow(data_frame)), chromosome = chromosome, start0 = start0, 
        end = end_value, width = end_value - start0, valid_coordinate = valid_coordinate, primary_chromosome = is_primary_chromosome(chromosome), 
        coordinate_key = ifelse(valid_coordinate, paste0(chromosome, ":", format(start0, scientific = FALSE, trim = TRUE), 
            "-", format(end_value, scientific = FALSE, trim = TRUE)), NA_character_), stringsAsFactors = FALSE), 
        data_frame)
}

coordinate_summary <- function(standardized_table, source_role) {
    valid <- standardized_table$valid_coordinate
    valid_table <- standardized_table[valid, , drop = FALSE]
    data.frame(source_id = unique(standardized_table$source_id), source_role = source_role, total_rows = nrow(standardized_table), 
        valid_coordinate_rows = sum(valid), invalid_coordinate_rows = sum(!valid), valid_coordinate_fraction = mean(valid), 
        primary_coordinate_rows = sum(valid_table$primary_chromosome), primary_coordinate_fraction = if (nrow(valid_table) > 
            0L) {
            mean(valid_table$primary_chromosome)
        }
        else {
            NA_real_
        }, unique_coordinate_n = length(unique(valid_table$coordinate_key)), duplicate_coordinate_fraction = if (nrow(valid_table) > 
            0L) {
            1 - length(unique(valid_table$coordinate_key))/nrow(valid_table)
        }
        else {
            NA_real_
        }, median_width = if (nrow(valid_table) > 0L) {
            stats::median(valid_table$width)
        }
        else {
            NA_real_
        }, minimum_width = if (nrow(valid_table) > 0L) {
            min(valid_table$width)
        }
        else {
            NA_real_
        }, maximum_width = if (nrow(valid_table) > 0L) {
            max(valid_table$width)
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
}

to_granges <- function(standardized_table) {
    valid_table <- standardized_table[standardized_table$valid_coordinate & standardized_table$primary_chromosome, 
        , drop = FALSE]
    GenomicRanges::GRanges(seqnames = valid_table$chromosome, ranges = IRanges::IRanges(start = as.integer(valid_table$start0 + 
        1), end = as.integer(valid_table$end)))
}

overlap_audit <- function(query_table, reference_table, audit_id) {
    query_valid <- query_table[query_table$valid_coordinate & query_table$primary_chromosome, , drop = FALSE]
    reference_valid <- reference_table[reference_table$valid_coordinate & reference_table$primary_chromosome, , 
        drop = FALSE]
    query_gr <- to_granges(query_valid)
    reference_gr <- to_granges(reference_valid)
    any_overlap <- GenomicRanges::countOverlaps(query_gr, reference_gr, ignore.strand = TRUE) > 0L
    exact_match <- query_valid$coordinate_key %in% reference_valid$coordinate_key
    data.frame(audit_id = audit_id, query_source = unique(query_valid$source_id), reference_source = unique(reference_valid$source_id), 
        query_primary_regions = length(query_gr), reference_primary_regions = length(reference_gr), any_overlap_n = sum(any_overlap), 
        any_overlap_fraction = mean(any_overlap), exact_match_n = sum(exact_match), exact_match_fraction = mean(exact_match), 
        stringsAsFactors = FALSE)
}

assign_sheet_role <- function(workbook_file, sheet_name) {
    key <- paste(workbook_file, sheet_name, sep = "::")
    role_map <- c(`media-3.xlsx::KeySources` = "ADMIN_KEY_RESOURCES", `media-3.xlsx::MouseDietInfo` = "ADMIN_DIET_INFORMATION", 
        `media-3.xlsx::List of Abbreviations` = "ADMIN_ABBREVIATIONS", `media-4.xlsx::siNTsiVDR_OpenChromatin` = "ATAC_DERIVED_ANNOTATION_RESTRICTED_CONTEXT", 
        `media-4.xlsx::siNT_ATAC` = "ATAC_CONDITION_LEVEL_PEAKS_RESTRICTED_CONTEXT", `media-4.xlsx::siVDR_ATAC` = "ATAC_CONDITION_LEVEL_PEAKS_RESTRICTED_CONTEXT", 
        `media-5.xlsx::DEGs` = "AUTHOR_RNA_REFERENCE_NOT_PROJECT_AUTHORITY", `media-5.xlsx::DEGs_LigandIndependentOverlap` = "AUTHOR_RNA_REFERENCE_NOT_PROJECT_AUTHORITY", 
        `media-5.xlsx::CanonicalPathways` = "AUTHOR_RNA_PATHWAY_REFERENCE", `media-5.xlsx::CanonicalPathways_OverlapDEGs` = "AUTHOR_RNA_PATHWAY_REFERENCE", 
        `media-5.xlsx::UpstreamRegulators` = "AUTHOR_RNA_REGULATOR_REFERENCE", `media-5.xlsx::UpstreamRegulators_OverlapDEGs` = "AUTHOR_RNA_REGULATOR_REFERENCE", 
        `media-6.xlsx::pLentiVeh` = "CORE_VDR_CISTROME_VEHICLE_PEAK_COORDINATES", `media-6.xlsx::pLenti125` = "CORE_VDR_CISTROME_CALCITRIOL_PEAK_COORDINATES", 
        `media-6.xlsx::AnnotationOfPeakDistribution` = "DERIVED_VDR_CISTROME_DISTRIBUTION", `media-6.xlsx::pLentiVeh_125overlap` = "CORE_VDR_CISTROME_SHARED_PEAK_COORDINATES", 
        `media-6.xlsx::HOMER_DeNovo` = "EMPTY_DE_NOVO_MOTIF_SHEET", `media-6.xlsx::HOMER_KnownMotif` = "DERIVED_VDR_CISTROME_KNOWN_MOTIFS", 
        `media-7.xlsx::Hg38_GSE129585_BED` = "EXTERNAL_KIDNEY_VDR_REFERENCE", `media-7.xlsx::AnnotationOfPeakDistribution` = "EXTERNAL_KIDNEY_DISTRIBUTION", 
        `media-7.xlsx::pLentiVeh_KIDNEYoverlap` = "EXTERNAL_KIDNEY_OVERLAP_REFERENCE", `media-7.xlsx::pLenti125_KIDNEYoverlap` = "EXTERNAL_KIDNEY_OVERLAP_REFERENCE", 
        `media-7.xlsx::HOMER_KIDNEY_overlapDeNovo` = "EMPTY_EXTERNAL_DE_NOVO_MOTIF_SHEET", `media-2.xlsx::OverlapAnnGenes_Peaks_DEGs` = "DERIVED_CISTROME_RNA_COMPOSITE_SUMMARY", 
        `media-2.xlsx::OverlapPeaksVeh_ReceptorDEG` = "DERIVED_VEHICLE_CISTROME_RNA_OVERLAP_COORDINATES", `media-2.xlsx::OverlapPeaks125_125DEG` = "DERIVED_CALCITRIOL_CISTROME_RNA_OVERLAP_COORDINATES", 
        `media-2.xlsx::HOMERDeNovo_OverlapPeaksDEGs` = "EMPTY_DE_NOVO_MOTIF_SHEET", `media-2.xlsx::HOMERKnownMotifOverlapPeaksDEGs` = "DERIVED_CISTROME_RNA_KNOWN_MOTIFS")
    if (!key %in% names(role_map)) {
        return("UNMAPPED_REQUIRES_REVIEW")
    }
    unname(role_map[[key]])
}

log_event("Starting", run_id)

required_paths <- c(gate = file.path(dirs$results, "STEP03B3A1_GATE_TABLE.csv"), workbook_inventory = file.path(dirs$inputs, 
    "STEP03B3A1_PUBLISHED_WORKBOOK_INVENTORY.csv"), sheet_inventory = file.path(dirs$inputs, "STEP03B3A1_PUBLISHED_WORKBOOK_SHEET_INVENTORY.csv"), 
    low_depth_decision = file.path(dirs$results, "STEP03B3A1_LOW_DEPTH_IgG_DECISION.csv"), sensitivity_rules = file.path(dirs$results, 
        "STEP03B3A1_LOCAL_PEAK_CALLING_SENSITIVITY_RULES.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3A1 output is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate <- .vdr_read_csv(required_paths[["gate"]], stringsAsFactors = FALSE, check.names = FALSE)

workbook_inventory <- .vdr_read_csv(required_paths[["workbook_inventory"]], stringsAsFactors = FALSE, check.names = FALSE)

sheet_inventory <- .vdr_read_csv(required_paths[["sheet_inventory"]], stringsAsFactors = FALSE, check.names = FALSE)

low_depth_decision <- .vdr_read_csv(required_paths[["low_depth_decision"]], stringsAsFactors = FALSE, check.names = FALSE)

sensitivity_rules <- .vdr_read_csv(required_paths[["sensitivity_rules"]], stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- gate$status[gate$gate_component == "STEP03B3A1 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_PUBLISHED_PEAK_TABLE_REVIEW") {
    stop("STEP03B3A1 does not authorize published-table review. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (nrow(low_depth_decision) != 1L || low_depth_decision$sample_decision != "COMPLETE_LOW_DEPTH_CONTROL_RETAINED" || 
    isTRUE(low_depth_decision$sample_excluded)) {
    stop("Low-depth IgG decision is not locked as retained.")
}

if (nrow(sensitivity_rules) != 6L || !all(sensitivity_rules$status == "LOCKED")) {
    stop("Local peak-calling sensitivity rules are incomplete.")
}

if (nrow(workbook_inventory) != 6L || nrow(sheet_inventory) != 28L) {
    stop("Expected six workbooks and 28 sheets.")
}

workbook_inventory$file_exists <- file.exists(workbook_inventory$workbook_path)

workbook_inventory$observed_size_bytes <- as.numeric(file.info(workbook_inventory$workbook_path)$size)

workbook_inventory$observed_md5 <- vapply(workbook_inventory$workbook_path, md5_file, character(1))

workbook_inventory$observed_sha256 <- vapply(workbook_inventory$workbook_path, sha256_file, character(1))

workbook_inventory$identity_pass <- with(workbook_inventory, file_exists & workbook_size_bytes == observed_size_bytes & 
    workbook_md5 == observed_md5 & workbook_sha256 == observed_sha256)

write_csv_utf8(workbook_inventory, file.path(dirs$inputs, "STEP03B3B_WORKBOOK_IDENTITY_RECHECK.csv"))

if (!all(workbook_inventory$identity_pass)) {
    stop("One or more published workbooks failed path, size or hash recheck.")
}

workbook_path_by_name <- setNames(workbook_inventory$workbook_path, workbook_inventory$workbook_file)

sheet_inventory$scientific_role <- mapply(assign_sheet_role, sheet_inventory$workbook_file, sheet_inventory$sheet_name, 
    USE.NAMES = FALSE)

sheet_inventory$role_status <- ifelse(sheet_inventory$scientific_role == "UNMAPPED_REQUIRES_REVIEW", "HOLD", "LOCKED_BY_CONTENT")

write_csv_utf8(sheet_inventory, file.path(dirs$results, "STEP03B3B_PUBLISHED_WORKBOOK_ROLE_MAP_LOCKED.csv"))

if (any(sheet_inventory$role_status == "HOLD")) {
    stop("At least one workbook sheet remains unmapped.")
}

media6_path <- workbook_path_by_name[["media-6.xlsx"]]

media2_path <- workbook_path_by_name[["media-2.xlsx"]]

log_event("Reading media-6 pLentiVeh as all-text")

vehicle_raw <- read_sheet_all_text(media6_path, "pLentiVeh", col_names = TRUE)

vehicle_chrom <- find_column_case_insensitive(vehicle_raw, "chrom")

vehicle_start <- find_column_case_insensitive(vehicle_raw, "start")

vehicle_end <- find_column_case_insensitive(vehicle_raw, "end")

vehicle <- standardize_coordinate_sheet(vehicle_raw, vehicle_chrom, vehicle_start, vehicle_end, "MEDIA6_pLentiVeh")

log_event("Reading media-6 pLenti125 as all-text")

calcitriol_raw <- read_sheet_all_text(media6_path, "pLenti125", col_names = TRUE)

calcitriol_chrom <- find_column_case_insensitive(calcitriol_raw, "chrom")

calcitriol_start <- find_column_case_insensitive(calcitriol_raw, "start")

calcitriol_end <- find_column_case_insensitive(calcitriol_raw, "end")

calcitriol <- standardize_coordinate_sheet(calcitriol_raw, calcitriol_chrom, calcitriol_start, calcitriol_end, 
    "MEDIA6_pLenti125")

log_event("Reading media-6 pLentiVeh_125overlap as all-text without headers")

shared_raw <- read_sheet_all_text(media6_path, "pLentiVeh_125overlap", col_names = FALSE)

shared <- standardize_headerless_overlap_sheet(shared_raw, "MEDIA6_pLentiVeh_125overlap")

log_event("Reading media-2 vehicle cistrome-RNA overlap as all-text")

vehicle_rna_raw <- read_sheet_all_text(media2_path, "OverlapPeaksVeh_ReceptorDEG", col_names = TRUE)

vehicle_rna_chrom <- find_column_case_insensitive(vehicle_rna_raw, "chrom")

vehicle_rna_start <- find_column_case_insensitive(vehicle_rna_raw, "start")

vehicle_rna_end <- find_column_case_insensitive(vehicle_rna_raw, "end")

vehicle_rna <- standardize_coordinate_sheet(vehicle_rna_raw, vehicle_rna_chrom, vehicle_rna_start, vehicle_rna_end, 
    "MEDIA2_OverlapPeaksVeh_ReceptorDEG")

log_event("Reading media-2 calcitriol cistrome-RNA overlap as all-text")

calcitriol_rna_raw <- read_sheet_all_text(media2_path, "OverlapPeaks125_125DEG", col_names = TRUE)

calcitriol_rna_chrom <- find_column_case_insensitive(calcitriol_rna_raw, "chrom")

calcitriol_rna_start <- find_column_case_insensitive(calcitriol_rna_raw, "start")

calcitriol_rna_end <- find_column_case_insensitive(calcitriol_rna_raw, "end")

calcitriol_rna <- standardize_coordinate_sheet(calcitriol_rna_raw, calcitriol_rna_chrom, calcitriol_rna_start, 
    calcitriol_rna_end, "MEDIA2_OverlapPeaks125_125DEG")

clean_tables <- list(vehicle = vehicle, calcitriol = calcitriol, shared = shared, vehicle_rna = vehicle_rna, calcitriol_rna = calcitriol_rna)

clean_paths <- c(vehicle = file.path(source_root, "MEDIA6_pLentiVeh_clean.tsv.gz"), calcitriol = file.path(source_root, 
    "MEDIA6_pLenti125_clean.tsv.gz"), shared = file.path(source_root, "MEDIA6_pLentiVeh_125overlap_clean.tsv.gz"), 
    vehicle_rna = file.path(source_root, "MEDIA2_OverlapPeaksVeh_ReceptorDEG_clean.tsv.gz"), calcitriol_rna = file.path(source_root, 
        "MEDIA2_OverlapPeaks125_125DEG_clean.tsv.gz"))

for (table_name in names(clean_tables)) {
    write_tsv_gz(clean_tables[[table_name]], clean_paths[[table_name]])
}

coordinate_summaries <- rbind(coordinate_summary(vehicle, "CORE_VDR_CISTROME_VEHICLE"), coordinate_summary(calcitriol, 
    "CORE_VDR_CISTROME_CALCITRIOL"), coordinate_summary(shared, "CORE_VDR_CISTROME_SHARED"), coordinate_summary(vehicle_rna, 
    "DERIVED_VEHICLE_CISTROME_RNA_OVERLAP"), coordinate_summary(calcitriol_rna, "DERIVED_CALCITRIOL_CISTROME_RNA_OVERLAP"))

write_csv_utf8(coordinate_summaries, file.path(dirs$qc, "STEP03B3B_COORDINATE_SCHEMA_SUMMARY.csv"))

overlap_audits <- rbind(overlap_audit(shared, vehicle, "SHARED_TO_VEHICLE"), overlap_audit(shared, calcitriol, 
    "SHARED_TO_CALCITRIOL"), overlap_audit(vehicle_rna, vehicle, "VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME"), overlap_audit(calcitriol_rna, 
    calcitriol, "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME"))

write_csv_utf8(overlap_audits, file.path(dirs$qc, "STEP03B3B_PUBLISHED_TABLE_RECONCILIATION.csv"))

coordinate_pass <- all(coordinate_summaries$valid_coordinate_fraction >= minimum_coordinate_valid_fraction)

source_overlap_pass <- all(overlap_audits$any_overlap_fraction >= minimum_source_overlap_fraction)

source_lock_candidate <- data.frame(source_id = c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME", "VDR_SHARED_CISTROME", 
    "VEHICLE_CISTROME_RNA_OVERLAP", "CALCITRIOL_CISTROME_RNA_OVERLAP"), workbook_file = c("media-6.xlsx", "media-6.xlsx", 
    "media-6.xlsx", "media-2.xlsx", "media-2.xlsx"), sheet_name = c("pLentiVeh", "pLenti125", "pLentiVeh_125overlap", 
    "OverlapPeaksVeh_ReceptorDEG", "OverlapPeaks125_125DEG"), scientific_role = c("PRIMARY_AUTHOR_PROCESSED_VDR_PEAKS_VEHICLE", 
    "PRIMARY_AUTHOR_PROCESSED_VDR_PEAKS_CALCITRIOL", "PRIMARY_AUTHOR_PROCESSED_SHARED_VDR_PEAKS", "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_VEHICLE", 
    "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_CALCITRIOL"), genome_build = "hg38_as_reported_by_source_article", clean_table_path = normalize_slashes(unname(clean_paths)), 
    clean_table_md5 = vapply(unname(clean_paths), md5_file, character(1)), clean_table_sha256 = vapply(unname(clean_paths), 
        sha256_file, character(1)), source_lock_status = if (coordinate_pass && source_overlap_pass) {
        "CANDIDATE_FOR_MANUAL_SOURCE_LOCK"
    } else {
        "HOLD_SCHEMA_OR_RECONCILIATION_FAILURE"
    }, stringsAsFactors = FALSE)

write_csv_utf8(source_lock_candidate, file.path(dirs$results, "STEP03B3B_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"))

role_counts <- as.data.frame(table(sheet_inventory$scientific_role), stringsAsFactors = FALSE)

names(role_counts) <- c("scientific_role", "n_sheets")

write_csv_utf8(role_counts, file.path(dirs$qc, "STEP03B3B_SHEET_ROLE_COUNTS.csv"))

important_noncoordinate_sheets <- data.frame(workbook_file = c("media-6.xlsx", "media-6.xlsx", "media-2.xlsx", 
    "media-2.xlsx"), sheet_name = c("AnnotationOfPeakDistribution", "HOMER_KnownMotif", "OverlapAnnGenes_Peaks_DEGs", 
    "HOMERKnownMotifOverlapPeaksDEGs"), role = c("DERIVED_DISTRIBUTION_SUMMARY", "DERIVED_KNOWN_MOTIF_SUMMARY", 
    "DERIVED_CISTROME_RNA_COMPOSITE_SUMMARY", "DERIVED_CISTROME_RNA_KNOWN_MOTIF_SUMMARY"), coordinate_source = FALSE, 
    future_use = c("descriptive validation only", "motif support after primary VDR source lock", "cross-check only; not primary coordinate source", 
        "motif support after project-specific RNA integration"), stringsAsFactors = FALSE)

write_csv_utf8(important_noncoordinate_sheets, file.path(dirs$results, "STEP03B3B_DERIVED_SHEET_USE_POLICY.csv"))

final_gate <- if (coordinate_pass && source_overlap_pass) {
    "PROCEED_TO_MANUAL_PUBLISHED_PEAK_SOURCE_LOCK"
} else if (coordinate_pass) {
    "HOLD_PUBLISHED_TABLE_RECONCILIATION"
} else {
    "HOLD_PUBLISHED_COORDINATE_SCHEMA"
}

gate_table <- data.frame(gate_component = c("STEP03B3A1 authorization", "Published workbook identities", "Workbook-sheet role map", 
    "Chromosome X/Y preservation", "Primary cistrome coordinate schemas", "Derived RNA-overlap coordinate schemas", 
    "Shared-peak reconciliation", "Cistrome-RNA table reconciliation", "Local peak calling", "RNA/CUT&RUN integration", 
    "Scientific figures", "STEP03B3B gate"), status = c("PASS", "PASS_6_OF_6", "LOCKED_28_OF_28", "PASS_ALL_RELEVANT_SHEETS_READ_AS_TEXT", 
    ifelse(all(coordinate_summaries$valid_coordinate_fraction[coordinate_summaries$source_role %in% c("CORE_VDR_CISTROME_VEHICLE", 
        "CORE_VDR_CISTROME_CALCITRIOL", "CORE_VDR_CISTROME_SHARED")] >= minimum_coordinate_valid_fraction), "PASS", 
        "FAIL"), ifelse(all(coordinate_summaries$valid_coordinate_fraction[grepl("^DERIVED_", coordinate_summaries$source_role)] >= 
        minimum_coordinate_valid_fraction), "PASS", "FAIL"), ifelse(all(overlap_audits$any_overlap_fraction[overlap_audits$audit_id %in% 
        c("SHARED_TO_VEHICLE", "SHARED_TO_CALCITRIOL")] >= minimum_source_overlap_fraction), "PASS", "FAIL"), ifelse(all(overlap_audits$any_overlap_fraction[overlap_audits$audit_id %in% 
        c("VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME", "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME")] >= minimum_source_overlap_fraction), 
        "PASS", "FAIL"), "NOT_STARTED", "NOT_STARTED", "NONE_CREATED", final_gate), required_action = c("Retain low-depth IgG decision and sensitivity rules.", 
    "Use only hash-validated workbooks.", "Use content-defined roles rather than filename order alone.", "Do not permit readxl numeric guessing for chromosome columns.", 
    "Review coordinate counts and duplicates.", "Treat media-2 overlap sheets as derived, not independent discovery.", 
    "Confirm shared table is supported by both primary cistrome sets.", "Confirm derived overlap tables reconcile to media-6 source peaks.", 
    "Do not run MACS3 before source-lock decision.", "Do not integrate until the source lock is manually approved.", 
    "No Skill-3 figure workflow is triggered.", "Upload coordinate summary, reconciliation and source-lock candidate."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3B_GATE_TABLE.csv"))

write_text(c("# STEP03B3B published peak-table audit card", "", paste0("- Run ID: ", run_id), "- Workbook identities: 6/6 PASS", 
    "- Sheet roles: 28/28 content-mapped", "- Core primary source: media-6 pLentiVeh and pLenti125", "- Shared source: media-6 pLentiVeh_125overlap", 
    "- Derived cistrome-RNA sources: two media-2 coordinate sheets", "- Chromosome columns: read as text", "- Genome build: hg38 as reported by the source article", 
    "- Local peak calling: NOT STARTED", "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Gate: ", final_gate), "", "The automatic 21-sheet keyword shortlist is superseded by the", "28-sheet content-defined role map. Keyword hits alone are not used", 
    "to determine scientific authority."), file.path(dirs$results, "STEP03B3B_PUBLISHED_PEAK_TABLE_AUDIT_CARD.md"))

write_text(c("# STEP03B3B Methods draft", "", "All supplementary Excel workbooks were revalidated by file size, MD5 and SHA-256. Workbook identity was determined from sheet content rather than attachment numbering alone. All candidate coordinate sheets were imported with every column forced to text, thereby preserving chromosome X and Y values and avoiding type-guessing losses.", 
    "", "The primary author-processed VDR cistrome tables were defined as the vehicle and calcitriol peak-coordinate sheets in media-6.xlsx. Their reported shared-peak sheet was evaluated against both primary cistrome sets. The vehicle and calcitriol cistrome-transcriptome overlap coordinate sheets in media-2.xlsx were treated as derived tables and reconciled to the corresponding media-6 primary peaks by genomic overlap. Coordinate validity, primary-chromosome composition, duplicate-coordinate burden and interval widths were recorded. No local peak calling or RNA integration was performed."), 
    file.path(dirs$text, "STEP03B3B_Methods_draft.md"))

write_text(c("# STEP03B3B interpretation constraints", "", "- media-6 primary coordinate tables are author-processed VDR peak sources.", 
    "- media-2 cistrome-RNA overlap tables are derived and are not independent validation.", "- media-5 RNA sheets do not replace RNA_PROGRAM_v1.0.", 
    "- media-4 ATAC sheets remain restricted descriptive context.", "- media-7 kidney sheets are external biological comparators, not T-HESC discovery data.", 
    "- known-motif and distribution sheets are secondary annotations.", "- empty de-novo motif sheets provide no positive evidence.", 
    "- no direct-target assignment is permitted before project-specific RNA integration and source lock."), file.path(dirs$text, 
    "STEP03B3B_Interpretation_constraints.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, workbooks = 6, sheets = 28, all_cells_read_as_text = TRUE, 
    genome_build = "hg38_as_reported_by_source_article", minimum_coordinate_valid_fraction = minimum_coordinate_valid_fraction, 
    minimum_source_overlap_fraction = minimum_source_overlap_fraction, local_peak_calling = FALSE, RNA_CUTRUN_integration = FALSE, 
    scientific_figures_created = FALSE, gate = final_gate)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3B_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("33_step03B3B_published_CUTRUN_peak_schema_audit.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "read_published_vdr_peaks.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

source_manifest <- data.frame(source_id = source_lock_candidate$source_id, clean_table_path = source_lock_candidate$clean_table_path, 
    size_bytes = as.numeric(file.info(source_lock_candidate$clean_table_path)$size), md5 = source_lock_candidate$clean_table_md5, 
    sha256 = source_lock_candidate$clean_table_sha256, upstream_workbook = source_lock_candidate$workbook_file, 
    upstream_sheet = source_lock_candidate$sheet_name, current_status = source_lock_candidate$source_lock_status, 
    stringsAsFactors = FALSE)

write_csv_utf8(source_manifest, file.path(source_root, "STEP03B3B_CLEAN_SOURCE_TABLE_MANIFEST.csv"))

archive_artifacts <- c(file.path(dirs$inputs, "STEP03B3B_WORKBOOK_IDENTITY_RECHECK.csv"), file.path(dirs$results, 
    "STEP03B3B_PUBLISHED_WORKBOOK_ROLE_MAP_LOCKED.csv"), file.path(dirs$results, "STEP03B3B_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"), 
    file.path(dirs$results, "STEP03B3B_DERIVED_SHEET_USE_POLICY.csv"), file.path(dirs$results, "STEP03B3B_GATE_TABLE.csv"), 
    file.path(dirs$results, "STEP03B3B_PUBLISHED_PEAK_TABLE_AUDIT_CARD.md"), file.path(dirs$qc, "STEP03B3B_COORDINATE_SCHEMA_SUMMARY.csv"), 
    file.path(dirs$qc, "STEP03B3B_PUBLISHED_TABLE_RECONCILIATION.csv"), file.path(dirs$qc, "STEP03B3B_SHEET_ROLE_COUNTS.csv"), 
    file.path(source_root, "STEP03B3B_CLEAN_SOURCE_TABLE_MANIFEST.csv"), file.path(dirs$text, "STEP03B3B_Methods_draft.md"), 
    file.path(dirs$text, "STEP03B3B_Interpretation_constraints.md"), file.path(dirs$parameters, "STEP03B3B_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B.txt"), file.path(dirs$code, "read_published_vdr_peaks.R"), 
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
    else if (grepl("Methods|Interpretation", basename(path_value))) {
        "txt"
    }
    else if (grepl("SUMMARY|RECONCILIATION|ROLE|POLICY|AUDIT_CARD|MANIFEST", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3B_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3B archive status", "", paste0("- Run ID: ", run_id), "- Workbook identities: 6/6 PASS", 
    "- Sheet role map: 28/28 LOCKED", paste0("- Clean coordinate sources: ", nrow(source_lock_candidate)), paste0("- Gate: ", 
        final_gate), "- Local peak calling: NOT STARTED", "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, "ctl", 
    "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3B", paste0("- Run ID: ", run_id), "- RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE", 
    "- ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", "- CUT&RUN fragment integrity: PASS", "- Low-depth IgG: RETAINED WITH SENSITIVITY RULES", 
    "- Published workbook roles: LOCKED", "- Published peak source lock: PENDING MANUAL REVIEW", "- Local peak calling: NOT STARTED", 
    "- RNA/CUT&RUN integration: NOT STARTED", paste0("- Gate: ", final_gate), "- Full STEP03 module: NOT LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", final_gate)

cat("\n============================================================\n")

cat("STEP03B3B published CUT&RUN peak-table audit completed.\n")

cat("Workbook identities: 6/6 PASS\n")

cat("Sheet roles: 28/28 content-mapped\n")

cat("Coordinate sheets read as text: 5\n")

cat("Local peak calling: NOT STARTED\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", final_gate, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3B_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3B_PUBLISHED_WORKBOOK_ROLE_MAP_LOCKED.csv\n")

cat("  06_results/STEP03B3B_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv\n")

cat("  08_qc/STEP03B3B_COORDINATE_SCHEMA_SUMMARY.csv\n")

cat("  08_qc/STEP03B3B_PUBLISHED_TABLE_RECONCILIATION.csv\n")

cat("============================================================\n")

