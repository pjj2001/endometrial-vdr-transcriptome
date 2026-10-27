#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A09: Harmonize chromosome naming

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3B_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

minimum_coordinate_valid_fraction <- 0.98999999999999999

minimum_source_overlap_fraction <- 0.94999999999999996

required_packages <- c("data.table", "digest", "jsonlite", "GenomicRanges", "IRanges")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(module_root, "00_module_control"), code = file.path(module_root, "02_code"), parameters = file.path(module_root, 
    "03_parameters"), environment = file.path(module_root, "04_environment"), execution = file.path(module_root, 
    "05_execution"), results = file.path(module_root, "06_results"), qc = file.path(module_root, "08_qc"), text = file.path(module_root, 
    "09_manuscript_text"), archive_parent = file.path(module_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

v1_source_root <- file.path(dirs$results, "STEP03B3B_published_peak_source_tables")

v2_source_root <- file.path(dirs$results, "STEP03B3B_v2_published_peak_source_tables_canonical_chr")

dir.create(v2_source_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3Bv2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

canonicalize_chromosome <- function(x) {
    value <- trimws(as.character(x))
    value <- sub("^chr", "", value, ignore.case = TRUE)
    value[value %in% c("MT", "Mt", "mt", "M", "m")] <- "M"
    value[value %in% c("x")] <- "X"
    value[value %in% c("y")] <- "Y"
    allowed <- c(as.character(1:22), "X", "Y", "M")
    canonical <- ifelse(value %in% allowed, paste0("chr", value), NA_character_)
    canonical
}

is_primary_canonical <- function(x) {
    x %in% paste0("chr", c(1:22, "X", "Y"))
}

refresh_table <- function(data_frame) {
    required_columns <- c("source_id", "source_row", "chromosome", "start0", "end", "width", "valid_coordinate", 
        "primary_chromosome", "coordinate_key")
    missing_columns <- setdiff(required_columns, names(data_frame))
    if (length(missing_columns) > 0L) {
        stop("Clean source table is missing columns: ", paste(missing_columns, collapse = ", "))
    }
    data_frame$chromosome_source <- data_frame$chromosome
    data_frame$chromosome <- canonicalize_chromosome(data_frame$chromosome_source)
    data_frame$start0 <- suppressWarnings(as.numeric(data_frame$start0))
    data_frame$end <- suppressWarnings(as.numeric(data_frame$end))
    data_frame$width <- data_frame$end - data_frame$start0
    data_frame$valid_coordinate <- !is.na(data_frame$chromosome) & is.finite(data_frame$start0) & is.finite(data_frame$end) & 
        data_frame$start0 >= 0 & data_frame$end > data_frame$start0
    data_frame$primary_chromosome <- data_frame$valid_coordinate & is_primary_canonical(data_frame$chromosome)
    data_frame$coordinate_key <- ifelse(data_frame$valid_coordinate, paste0(data_frame$chromosome, ":", format(data_frame$start0, 
        scientific = FALSE, trim = TRUE), "-", format(data_frame$end, scientific = FALSE, trim = TRUE)), NA_character_)
    data_frame
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
        }, chromosome_style = "canonical_UCSC_chr", stringsAsFactors = FALSE)
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
        query_seqlevels = paste(sort(unique(as.character(GenomicRanges::seqnames(query_gr)))), collapse = ";"), 
        reference_seqlevels = paste(sort(unique(as.character(GenomicRanges::seqnames(reference_gr)))), collapse = ";"), 
        sequence_level_intersection_n = length(intersect(unique(as.character(GenomicRanges::seqnames(query_gr))), 
            unique(as.character(GenomicRanges::seqnames(reference_gr))))), stringsAsFactors = FALSE)
}

log_event("Starting", run_id)

required_paths <- c(gate_v1 = file.path(dirs$results, "STEP03B3B_GATE_TABLE.csv"), source_lock_v1 = file.path(dirs$results, 
    "STEP03B3B_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"), coordinate_summary_v1 = file.path(dirs$qc, "STEP03B3B_COORDINATE_SCHEMA_SUMMARY.csv"), 
    reconciliation_v1 = file.path(dirs$qc, "STEP03B3B_PUBLISHED_TABLE_RECONCILIATION.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3B v1 output is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v1 <- .vdr_read_csv(required_paths[["gate_v1"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock_v1 <- .vdr_read_csv(required_paths[["source_lock_v1"]], stringsAsFactors = FALSE, check.names = FALSE)

reconciliation_v1 <- .vdr_read_csv(required_paths[["reconciliation_v1"]], stringsAsFactors = FALSE, check.names = FALSE)

v1_gate_status <- gate_v1$status[gate_v1$gate_component == "STEP03B3B gate"]

if (length(v1_gate_status) != 1L || v1_gate_status != "HOLD_PUBLISHED_TABLE_RECONCILIATION") {
    stop("Unexpected STEP03B3B v1 gate: ", paste(v1_gate_status, collapse = ", "))
}

failed_reconciliation_ids <- reconciliation_v1$audit_id[reconciliation_v1$any_overlap_fraction < minimum_source_overlap_fraction]

expected_failed_ids <- c("SHARED_TO_VEHICLE", "SHARED_TO_CALCITRIOL")

if (!setequal(failed_reconciliation_ids, expected_failed_ids)) {
    stop("The v1 reconciliation failures are not restricted to the two shared-peak checks.")
}

source_paths <- setNames(source_lock_v1$clean_table_path, source_lock_v1$source_id)

if (!all(file.exists(source_paths))) {
    stop("One or more v1 clean source tables are missing.")
}

source_id_to_role <- c(VDR_VEHICLE_CISTROME = "CORE_VDR_CISTROME_VEHICLE", VDR_CALCITRIOL_CISTROME = "CORE_VDR_CISTROME_CALCITRIOL", 
    VDR_SHARED_CISTROME = "CORE_VDR_CISTROME_SHARED", VEHICLE_CISTROME_RNA_OVERLAP = "DERIVED_VEHICLE_CISTROME_RNA_OVERLAP", 
    CALCITRIOL_CISTROME_RNA_OVERLAP = "DERIVED_CALCITRIOL_CISTROME_RNA_OVERLAP")

canonical_tables <- list()

canonical_paths <- character(0)

summary_rows <- list()

invalid_rows <- list()

style_rows <- list()

for (source_id in names(source_paths)) {
    log_event("Canonicalizing", source_id)
    table_v1 <- .vdr_fread(source_paths[[source_id]], sep = "\t", data.table = FALSE, showProgress = FALSE)
    table_v2 <- refresh_table(table_v1)
    canonical_tables[[source_id]] <- table_v2
    output_path <- file.path(v2_source_root, paste0(source_id, "_canonical_chr.tsv.gz"))
    write_tsv_gz(table_v2, output_path)
    canonical_paths[source_id] <- output_path
    summary_rows[[source_id]] <- coordinate_summary(table_v2, source_id_to_role[[source_id]])
    invalid_subset <- table_v2[!table_v2$valid_coordinate, , drop = FALSE]
    if (nrow(invalid_subset) > 0L) {
        invalid_rows[[source_id]] <- invalid_subset
    }
    source_chromosomes <- unique(table_v2$chromosome_source)
    canonical_chromosomes <- unique(table_v2$chromosome)
    style_rows[[source_id]] <- data.frame(source_id = source_id, original_chr_prefixed_fraction = mean(grepl("^chr", 
        table_v2$chromosome_source, ignore.case = TRUE), na.rm = TRUE), original_unique_chromosome_n = length(source_chromosomes), 
        canonical_unique_chromosome_n = length(canonical_chromosomes[!is.na(canonical_chromosomes)]), canonicalization_unmapped_rows = sum(is.na(table_v2$chromosome)), 
        stringsAsFactors = FALSE)
}

coordinate_summary_v2 <- do.call(rbind, summary_rows)

chromosome_style_summary <- do.call(rbind, style_rows)

invalid_coordinate_rows <- if (length(invalid_rows) > 0L) {
    data.table::rbindlist(invalid_rows, fill = TRUE, idcol = "source_id_from_list")
} else {
    data.frame(source_id_from_list = character(0), stringsAsFactors = FALSE)
}

write_csv_utf8(coordinate_summary_v2, file.path(dirs$qc, "STEP03B3B_v2_COORDINATE_SCHEMA_SUMMARY.csv"))

write_csv_utf8(chromosome_style_summary, file.path(dirs$qc, "STEP03B3B_v2_CHROMOSOME_STYLE_SUMMARY.csv"))

write_csv_utf8(invalid_coordinate_rows, file.path(dirs$qc, "STEP03B3B_v2_INVALID_COORDINATE_ROWS.csv"))

vehicle <- canonical_tables[["VDR_VEHICLE_CISTROME"]]

calcitriol <- canonical_tables[["VDR_CALCITRIOL_CISTROME"]]

shared <- canonical_tables[["VDR_SHARED_CISTROME"]]

vehicle_rna <- canonical_tables[["VEHICLE_CISTROME_RNA_OVERLAP"]]

calcitriol_rna <- canonical_tables[["CALCITRIOL_CISTROME_RNA_OVERLAP"]]

reconciliation_v2 <- rbind(overlap_audit(shared, vehicle, "SHARED_TO_VEHICLE"), overlap_audit(shared, calcitriol, 
    "SHARED_TO_CALCITRIOL"), overlap_audit(vehicle_rna, vehicle, "VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME"), overlap_audit(calcitriol_rna, 
    calcitriol, "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME"))

write_csv_utf8(reconciliation_v2, file.path(dirs$qc, "STEP03B3B_v2_PUBLISHED_TABLE_RECONCILIATION.csv"))

coordinate_pass <- all(coordinate_summary_v2$valid_coordinate_fraction >= minimum_coordinate_valid_fraction)

source_overlap_pass <- all(reconciliation_v2$any_overlap_fraction >= minimum_source_overlap_fraction)

shared_reconciliation_pass <- all(reconciliation_v2$any_overlap_fraction[reconciliation_v2$audit_id %in% c("SHARED_TO_VEHICLE", 
    "SHARED_TO_CALCITRIOL")] >= minimum_source_overlap_fraction)

derived_reconciliation_pass <- all(reconciliation_v2$any_overlap_fraction[reconciliation_v2$audit_id %in% c("VEHICLE_RNA_OVERLAP_TO_VEHICLE_CISTROME", 
    "CALCITRIOL_RNA_OVERLAP_TO_CALCITRIOL_CISTROME")] >= minimum_source_overlap_fraction)

final_gate <- if (coordinate_pass && source_overlap_pass) {
    "PROCEED_TO_MANUAL_PUBLISHED_PEAK_SOURCE_LOCK"
} else if (coordinate_pass) {
    "HOLD_PUBLISHED_TABLE_RECONCILIATION"
} else {
    "HOLD_PUBLISHED_COORDINATE_SCHEMA"
}

source_lock_v2 <- source_lock_v1

source_lock_v2$clean_table_path <- normalize_slashes(unname(canonical_paths[source_lock_v2$source_id]))

source_lock_v2$clean_table_md5 <- vapply(source_lock_v2$clean_table_path, md5_file, character(1))

source_lock_v2$clean_table_sha256 <- vapply(source_lock_v2$clean_table_path, sha256_file, character(1))

source_lock_v2$chromosome_style <- "canonical_UCSC_chr"

source_lock_v2$source_lock_status <- if (final_gate == "PROCEED_TO_MANUAL_PUBLISHED_PEAK_SOURCE_LOCK") {
    "CANDIDATE_FOR_MANUAL_SOURCE_LOCK"
} else {
    "HOLD_SCHEMA_OR_RECONCILIATION_FAILURE"
}

write_csv_utf8(source_lock_v2, file.path(dirs$results, "STEP03B3B_v2_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"))

canonical_manifest <- data.frame(source_id = names(canonical_paths), canonical_table_path = normalize_slashes(unname(canonical_paths)), 
    size_bytes = as.numeric(file.info(unname(canonical_paths))$size), md5 = vapply(unname(canonical_paths), md5_file, 
        character(1)), sha256 = vapply(unname(canonical_paths), sha256_file, character(1)), upstream_v1_table_path = normalize_slashes(unname(source_paths[names(canonical_paths)])), 
    transformation = "chromosome_style_only_to_canonical_UCSC_chr", coordinate_membership_changed = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(canonical_manifest, file.path(v2_source_root, "STEP03B3B_v2_CANONICAL_SOURCE_TABLE_MANIFEST.csv"))

gate_table <- data.frame(gate_component = c("STEP03B3B v1 hold", "v1 clean-source identities", "Chromosome-style diagnosis", 
    "Canonical chromosome conversion", "Coordinate membership", "Coordinate schema", "Shared-peak reconciliation", 
    "Cistrome-RNA reconciliation", "Local peak calling", "RNA/CUT&RUN integration", "Scientific figures", "STEP03B3B v2 gate"), 
    status = c("CONFIRMED_SHARED_RECONCILIATION_ONLY", "PASS_5_OF_5", "CONFIRMED_CHR_PREFIX_MISMATCH", "PASS_5_OF_5", 
        "UNCHANGED", ifelse(coordinate_pass, "PASS", "FAIL"), ifelse(shared_reconciliation_pass, "PASS", "FAIL"), 
        ifelse(derived_reconciliation_pass, "PASS", "FAIL"), "NOT_STARTED", "NOT_STARTED", "NONE_CREATED", final_gate), 
    required_action = c("Do not interpret v1 zero overlaps as biological absence.", "Retain v1 hashes and paths.", 
        "Record the software-level sequence-style mismatch.", "Use canonical chr names for all later genomic joins.", 
        "Do not add, drop or resize valid peaks.", "Review invalid-coordinate rows separately.", "Manually confirm overlap fractions before source lock.", 
        "Retain media-2 tables as derived, not independent validation.", "Do not run MACS3 before source-lock decision.", 
        "Do not integrate until manual source lock.", "No Skill-3 figure workflow is triggered.", "Upload v2 gate, summaries and reconciliation."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3B_v2_GATE_TABLE.csv"))

write_text(c("# STEP03B3B v2 chromosome-style repair card", "", paste0("- Run ID: ", run_id), "- Root cause: chr-prefixed shared table versus non-prefixed primary tables", 
    "- v1 zero-overlap result: SOFTWARE SEQUENCE-STYLE ARTIFACT", "- Peak membership: UNCHANGED", "- Start/end coordinates: UNCHANGED", 
    "- Canonical style: UCSC chr1...chr22, chrX, chrY, chrM", "- Local peak calling: NOT STARTED", "- RNA/CUT&RUN integration: NOT STARTED", 
    "- Scientific figures: NONE", paste0("- Gate: ", final_gate)), file.path(dirs$results, "STEP03B3B_v2_CHROMOSOME_STYLE_REPAIR_CARD.md"))

write_text(c("# STEP03B3B v2 Methods addendum", "", "The initial genomic reconciliation returned zero overlap for the published shared VDR-peak table because that table used UCSC-style chromosome labels (for example, chr1), whereas the primary vehicle and calcitriol peak tables used non-prefixed labels (for example, 1). The coordinate tables were therefore harmonized to canonical UCSC-style labels without changing interval membership or start/end coordinates. Genomic overlap reconciliation was then repeated using the harmonized labels.", 
    "", "Rows that could not be mapped to canonical chromosomes or valid numeric coordinates were retained in a separate invalid-row audit and excluded only from genomic overlap calculations."), 
    file.path(dirs$text, "STEP03B3B_v2_Methods_addendum.md"))

parameters <- list(schema_version = "2.0", run_id = run_id, root_cause = "chr_prefix_sequence_level_mismatch", 
    source_tables_reused = 5, chromosome_style = "canonical_UCSC_chr", peak_membership_changed = FALSE, interval_coordinates_changed = FALSE, 
    minimum_coordinate_valid_fraction = minimum_coordinate_valid_fraction, minimum_source_overlap_fraction = minimum_source_overlap_fraction, 
    local_peak_calling = FALSE, RNA_CUTRUN_integration = FALSE, scientific_figures_created = FALSE, gate = final_gate)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3B_v2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B_v2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("34_step03B3B_v2_harmonize_chromosome_style_and_reconcile.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "harmonize_peak_chromosomes.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3B_v2_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3B_v2_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv"), 
    file.path(dirs$results, "STEP03B3B_v2_CHROMOSOME_STYLE_REPAIR_CARD.md"), file.path(dirs$qc, "STEP03B3B_v2_COORDINATE_SCHEMA_SUMMARY.csv"), 
    file.path(dirs$qc, "STEP03B3B_v2_CHROMOSOME_STYLE_SUMMARY.csv"), file.path(dirs$qc, "STEP03B3B_v2_INVALID_COORDINATE_ROWS.csv"), 
    file.path(dirs$qc, "STEP03B3B_v2_PUBLISHED_TABLE_RECONCILIATION.csv"), file.path(v2_source_root, "STEP03B3B_v2_CANONICAL_SOURCE_TABLE_MANIFEST.csv"), 
    file.path(dirs$text, "STEP03B3B_v2_Methods_addendum.md"), file.path(dirs$parameters, "STEP03B3B_v2_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3B_v2.txt"), file.path(dirs$code, "harmonize_peak_chromosomes.R"), 
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
    else if (grepl("Methods", basename(path_value))) {
        "txt"
    }
    else if (grepl("SUMMARY|RECONCILIATION|REPAIR_CARD|MANIFEST|INVALID", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3B_v2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3B v2 archive status", "", paste0("- Run ID: ", run_id), "- Root cause: chr-prefix sequence-style mismatch", 
    "- Peak membership changed: NO", "- Coordinate values changed: NO", "- Canonical source tables: 5", paste0("- Gate: ", 
        final_gate), "- Local peak calling: NOT STARTED", "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, "ctl", 
    "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3B v2", paste0("- Run ID: ", run_id), "- RNA anchor: RNA_PROGRAM_v1.0 IMMUTABLE", 
    "- ATAC role: RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY", "- CUT&RUN fragment integrity: PASS", "- Low-depth IgG: RETAINED WITH SENSITIVITY RULES", 
    "- Published workbook roles: LOCKED", "- Chromosome-style mismatch: REPAIRED", "- Peak source lock: PENDING MANUAL REVIEW", 
    "- Local peak calling: NOT STARTED", "- RNA/CUT&RUN integration: NOT STARTED", paste0("- Gate: ", final_gate), 
    "- Full STEP03 module: NOT LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", final_gate)

cat("\n============================================================\n")

cat("STEP03B3B v2 chromosome-style reconciliation completed.\n")

cat("Root cause: chr-prefix sequence-level mismatch\n")

cat("Canonical source tables: 5\n")

cat("Peak membership changed: NO\n")

cat("Start/end coordinates changed: NO\n")

cat("Shared-to-vehicle overlap: ", round(reconciliation_v2$any_overlap_fraction[reconciliation_v2$audit_id == "SHARED_TO_VEHICLE"], 
    4), "\n")

cat("Shared-to-calcitriol overlap: ", round(reconciliation_v2$any_overlap_fraction[reconciliation_v2$audit_id == 
    "SHARED_TO_CALCITRIOL"], 4), "\n")

cat("Local peak calling: NOT STARTED\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", final_gate, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3B_v2_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3B_v2_PUBLISHED_PEAK_SOURCE_LOCK_CANDIDATE.csv\n")

cat("  08_qc/STEP03B3B_v2_COORDINATE_SCHEMA_SUMMARY.csv\n")

cat("  08_qc/STEP03B3B_v2_CHROMOSOME_STYLE_SUMMARY.csv\n")

cat("  08_qc/STEP03B3B_v2_INVALID_COORDINATE_ROWS.csv\n")

cat("  08_qc/STEP03B3B_v2_PUBLISHED_TABLE_RECONCILIATION.csv\n")

cat("============================================================\n")

