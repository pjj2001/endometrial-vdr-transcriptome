#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A16: Resolve the documented symbol aliases

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3A_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

RNA_source_candidate_version <- "RNA_INTERACTION_SOURCE_v1.0_CANDIDATE"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_RNA_gene_n <- 436L

expected_exact_n <- 279L

expected_unmatched_n <- 157L

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

audit_root <- file.path(dirs$results, "STEP03B3C3A_v3_punctuation_version_audit")

dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3Av3_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

clean_symbol <- function(x) {
    value <- trimws(as.character(x))
    value[is.na(value) | !nzchar(value) | toupper(value) %in% c("NA", "N/A", "NULL", "NONE", ".")] <- NA_character_
    value
}

punctuation_key <- function(x) {
    toupper(gsub("[._-]", "", as.character(x), perl = TRUE))
}

empty_candidate_table <- function() {
    data.frame(RNA_row = integer(0), RNA_gene_symbol = character(0), candidate_cistrome_symbol = character(0), 
        candidate_method = character(0), candidate_priority = integer(0), candidate_is_unique_for_method = logical(0), 
        automatic_mapping_applied = logical(0), stringsAsFactors = FALSE)
}

log_event("Starting", run_id)

locked_map_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

v2_audit_root <- file.path(dirs$results, "STEP03B3C3A_v2_RNA_source_and_gene_ID_audit")

required_paths <- c(gate_v2 = file.path(dirs$results, "STEP03B3C3A_v2_GATE_TABLE.csv"), source_candidate = file.path(dirs$results, 
    "STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv"), RNA_gene_audit = file.path(v2_audit_root, "STEP03B3C3A_v2_RNA_GENE_ID_AUDIT.csv"), 
    match_summary = file.path(v2_audit_root, "STEP03B3C3A_v2_RNA_CISTROME_MATCH_SUMMARY.csv"), collision_audit = file.path(v2_audit_root, 
        "STEP03B3C3A_v2_CASE_COLLISION_AUDIT.csv"), unmatched_audit = file.path(v2_audit_root, "STEP03B3C3A_v2_UNMATCHED_RNA_SYMBOL_AUDIT.csv"), 
    map_manifest = file.path(locked_map_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3A v3 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v2 <- .vdr_read_csv(required_paths[["gate_v2"]], stringsAsFactors = FALSE, check.names = FALSE)

source_candidate <- .vdr_read_csv(required_paths[["source_candidate"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_gene_audit <- .vdr_read_csv(required_paths[["RNA_gene_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

match_summary <- .vdr_read_csv(required_paths[["match_summary"]], stringsAsFactors = FALSE, check.names = FALSE)

collision_audit <- .vdr_read_csv(required_paths[["collision_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

unmatched_audit <- .vdr_read_csv(required_paths[["unmatched_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

map_manifest <- .vdr_read_csv(required_paths[["map_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- gate_v2$status[gate_v2$gate_component == "STEP03B3C3A v2 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_RNA_SOURCE_AND_GENE_ID_LOCK") {
    stop("STEP03B3C3A v2 does not authorize manual gene-ID review. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (nrow(source_candidate) != 1L || source_candidate$source_candidate_version != RNA_source_candidate_version || 
    source_candidate$gene_symbol_column != "feature_id" || source_candidate$n_rows != expected_RNA_gene_n || source_candidate$unique_gene_n != 
    expected_RNA_gene_n) {
    stop("RNA source candidate metadata failed validation.")
}

summary_value <- setNames(match_summary$value, match_summary$metric)

required_summary <- c(`RNA rows` = expected_RNA_gene_n, `Unique RNA symbols` = expected_RNA_gene_n, `Valid RNA symbol shape` = expected_RNA_gene_n, 
    `Exact cistrome matches` = expected_exact_n, `Case-only cistrome matches` = 0L, `No cistrome match` = expected_unmatched_n, 
    `RNA case-collision rows` = 0L, `Cistrome case-collision rows` = 0L, `Cistrome symbol universe` = expected_cistrome_gene_n, 
    `Automatic alias mappings` = 0L)

missing_summary_metrics <- setdiff(names(required_summary), names(summary_value))

if (length(missing_summary_metrics) > 0L || !all(as.integer(summary_value[names(required_summary)]) == as.integer(required_summary))) {
    stop("RNA/cistrome match summary differs from the reviewed v2 result.")
}

if (nrow(RNA_gene_audit) != expected_RNA_gene_n || length(unique(RNA_gene_audit$RNA_row)) != expected_RNA_gene_n || 
    length(unique(RNA_gene_audit$RNA_gene_symbol)) != expected_RNA_gene_n || !all(RNA_gene_audit$valid_symbol_shape)) {
    stop("RNA gene-ID audit failed row, uniqueness or symbol-shape checks.")
}

partition_n <- rowSums(RNA_gene_audit[c("exact_cistrome_match", "case_only_cistrome_match", "no_cistrome_match")])

if (!all(partition_n == 1L) || sum(RNA_gene_audit$exact_cistrome_match) != expected_exact_n || sum(RNA_gene_audit$no_cistrome_match) != 
    expected_unmatched_n || nrow(collision_audit) != 0L) {
    stop("RNA exact/case-only/unmatched partition failed validation.")
}

expected_unmatched_symbols <- RNA_gene_audit$RNA_gene_symbol[RNA_gene_audit$no_cistrome_match]

if (nrow(unmatched_audit) != expected_unmatched_n || !identical(as.character(unmatched_audit$RNA_gene_symbol), 
    as.character(expected_unmatched_symbols))) {
    stop("Unmatched-symbol audit is not the exact v2 unmatched subset.")
}

if (nrow(map_manifest) != 2L || !all(map_manifest$lock_version == peak_gene_map_lock_version) || !all(map_manifest$lock_status == 
    "LOCKED")) {
    stop("PEAK_GENE_MAP_v1.1 manifest is incomplete.")
}

map_manifest$file_exists <- file.exists(map_manifest$path)

map_manifest$observed_md5 <- vapply(map_manifest$path, md5_file, character(1))

map_manifest$observed_sha256 <- vapply(map_manifest$path, sha256_file, character(1))

map_manifest$identity_pass <- with(map_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(map_manifest$identity_pass)) {
    stop("A locked PEAK_GENE_MAP_v1.1 artifact failed hash validation.")
}

map_long_path <- map_manifest$path[map_manifest$artifact_id == "PEAK_GENE_LONG"]

if (length(map_long_path) != 1L) {
    stop("Locked peak-gene long map could not be uniquely identified.")
}

cistrome_map <- .vdr_fread(map_long_path, sep = "\t", data.table = FALSE, showProgress = FALSE, select = "gene_symbol_author")

cistrome_symbols <- sort(unique(clean_symbol(cistrome_map$gene_symbol_author)))

cistrome_symbols <- cistrome_symbols[!is.na(cistrome_symbols)]

if (length(cistrome_symbols) != expected_cistrome_gene_n) {
    stop("Locked cistrome symbol universe differs from 10,042.")
}

cistrome_key <- punctuation_key(cistrome_symbols)

candidate_rows <- vector("list", nrow(unmatched_audit))

unmatched_class_rows <- vector("list", nrow(unmatched_audit))

for (row_index in seq_len(nrow(unmatched_audit))) {
    RNA_row <- unmatched_audit$RNA_row[row_index]
    RNA_symbol <- as.character(unmatched_audit$RNA_gene_symbol[row_index])
    row_candidates <- list()
    if (grepl(".", RNA_symbol, fixed = TRUE)) {
        dot_to_hyphen <- gsub(".", "-", RNA_symbol, fixed = TRUE)
        if (dot_to_hyphen %in% cistrome_symbols) {
            row_candidates[[length(row_candidates) + 1L]] <- data.frame(RNA_row = RNA_row, RNA_gene_symbol = RNA_symbol, 
                candidate_cistrome_symbol = dot_to_hyphen, candidate_method = "EXACT_DOT_TO_HYPHEN", candidate_priority = 1L, 
                candidate_is_unique_for_method = TRUE, automatic_mapping_applied = FALSE, stringsAsFactors = FALSE)
        }
    }
    if (grepl("^(AC|AL|AP)[0-9]+\\.[0-9]+$", RNA_symbol, perl = TRUE)) {
        version_stripped <- sub("\\.[0-9]+$", "", RNA_symbol, perl = TRUE)
        if (version_stripped %in% cistrome_symbols) {
            row_candidates[[length(row_candidates) + 1L]] <- data.frame(RNA_row = RNA_row, RNA_gene_symbol = RNA_symbol, 
                candidate_cistrome_symbol = version_stripped, candidate_method = "EXACT_ACCESSION_VERSION_STRIP", 
                candidate_priority = 2L, candidate_is_unique_for_method = TRUE, automatic_mapping_applied = FALSE, 
                stringsAsFactors = FALSE)
        }
    }
    punctuation_candidates <- cistrome_symbols[cistrome_key == punctuation_key(RNA_symbol)]
    punctuation_candidates <- unique(punctuation_candidates)
    if (length(punctuation_candidates) > 0L) {
        row_candidates[[length(row_candidates) + 1L]] <- data.frame(RNA_row = rep(RNA_row, length(punctuation_candidates)), 
            RNA_gene_symbol = rep(RNA_symbol, length(punctuation_candidates)), candidate_cistrome_symbol = punctuation_candidates, 
            candidate_method = "PUNCTUATION_INSENSITIVE_KEY", candidate_priority = 3L, candidate_is_unique_for_method = length(punctuation_candidates) == 
                1L, automatic_mapping_applied = FALSE, stringsAsFactors = FALSE)
    }
    if (length(row_candidates) > 0L) {
        combined_candidates <- do.call(rbind, row_candidates)
        combined_candidates <- combined_candidates[!duplicated(paste(combined_candidates$RNA_gene_symbol, combined_candidates$candidate_cistrome_symbol, 
            combined_candidates$candidate_method, sep = "::")), , drop = FALSE]
        candidate_rows[[row_index]] <- combined_candidates
    }
    candidate_symbols <- unique(unlist(lapply(row_candidates, function(candidate_table) {
        candidate_table$candidate_cistrome_symbol
    }), use.names = FALSE))
    classification <- if (length(candidate_symbols) == 0L) {
        "NO_DETERMINISTIC_PUNCTUATION_OR_VERSION_CANDIDATE"
    }
    else if (length(candidate_symbols) == 1L) {
        "ONE_UNIQUE_DETERMINISTIC_CANDIDATE"
    }
    else {
        "MULTIPLE_DETERMINISTIC_CANDIDATES"
    }
    unmatched_class_rows[[row_index]] <- data.frame(RNA_row = RNA_row, RNA_gene_symbol = RNA_symbol, contains_dot = grepl(".", 
        RNA_symbol, fixed = TRUE), accession_style_versioned = grepl("^(AC|AL|AP)[0-9]+\\.[0-9]+$", RNA_symbol, 
        perl = TRUE), candidate_symbol_n = length(candidate_symbols), candidate_symbols = paste(candidate_symbols, 
        collapse = ";"), review_class = classification, automatic_mapping_applied = FALSE, stringsAsFactors = FALSE)
}

candidate_table <- data.table::rbindlist(candidate_rows, fill = TRUE)

candidate_table <- as.data.frame(candidate_table, stringsAsFactors = FALSE)

if (nrow(candidate_table) == 0L) {
    candidate_table <- empty_candidate_table()
}

unmatched_classification <- do.call(rbind, unmatched_class_rows)

write_csv_utf8(candidate_table, file.path(audit_root, "STEP03B3C3A_v3_PUNCTUATION_VERSION_CANDIDATES.csv"))

write_csv_utf8(unmatched_classification, file.path(audit_root, "STEP03B3C3A_v3_UNMATCHED_SYMBOL_CLASSIFICATION.csv"))

unique_candidate_RNA_n <- sum(unmatched_classification$candidate_symbol_n > 0L)

single_candidate_RNA_n <- sum(unmatched_classification$candidate_symbol_n == 1L)

multiple_candidate_RNA_n <- sum(unmatched_classification$candidate_symbol_n > 1L)

no_candidate_RNA_n <- sum(unmatched_classification$candidate_symbol_n == 0L)

dotted_unmatched_n <- sum(unmatched_classification$contains_dot)

versioned_accession_n <- sum(unmatched_classification$accession_style_versioned)

method_summary <- if (nrow(candidate_table) > 0L) {
    as.data.frame(table(candidate_table$candidate_method), stringsAsFactors = FALSE)
} else {
    data.frame(Var1 = character(0), Freq = integer(0), stringsAsFactors = FALSE)
}

names(method_summary) <- c("candidate_method", "candidate_row_n")

summary_table <- data.frame(metric = c("RNA symbols", "Exact cistrome matches retained", "Unmatched symbols audited", 
    "Dotted unmatched symbols", "Versioned AC/AL/AP unmatched symbols", "Unmatched symbols with at least one deterministic candidate", 
    "Unmatched symbols with exactly one deterministic candidate", "Unmatched symbols with multiple deterministic candidates", 
    "Unmatched symbols with no deterministic candidate", "Automatic mappings applied", "RNA membership changed", 
    "Peak-gene membership changed"), value = c(expected_RNA_gene_n, expected_exact_n, expected_unmatched_n, dotted_unmatched_n, 
    versioned_accession_n, unique_candidate_RNA_n, single_candidate_RNA_n, multiple_candidate_RNA_n, no_candidate_RNA_n, 
    0L, 0L, 0L), stringsAsFactors = FALSE)

write_csv_utf8(summary_table, file.path(dirs$qc, "STEP03B3C3A_v3_PUNCTUATION_VERSION_AUDIT_SUMMARY.csv"))

write_csv_utf8(method_summary, file.path(dirs$qc, "STEP03B3C3A_v3_CANDIDATE_METHOD_SUMMARY.csv"))

gate_status <- if (multiple_candidate_RNA_n > 0L) {
    "HOLD_MULTIPLE_PUNCTUATION_VERSION_CANDIDATES"
} else if (single_candidate_RNA_n > 0L) {
    "PROCEED_TO_MANUAL_PUNCTUATION_VERSION_CANDIDATE_REVIEW"
} else {
    "PROCEED_TO_LOCK_EXACT_ONLY_GENE_ID_POLICY"
}

gate_table <- data.frame(gate_component = c("STEP03B3C3A v2 gate", "RNA source candidate", "RNA gene-symbol field", 
    "RNA exact-match partition", "Peak-gene map lock", "Dotted unmatched symbols", "Versioned accession-style unmatched symbols", 
    "Unique deterministic candidates", "Multiple deterministic candidates", "Automatic punctuation mapping", "Automatic alias mapping", 
    "RNA membership", "Peak-gene membership", "RNA/CUT&RUN integration", "Definitive direct-target assignment", 
    "Scientific figures", "STEP03B3C3A v3 gate"), status = c("PASS", RNA_source_candidate_version, "feature_id", 
    "279_EXACT_157_UNMATCHED", paste0("IMMUTABLE_", peak_gene_map_lock_version), as.character(dotted_unmatched_n), 
    as.character(versioned_accession_n), as.character(single_candidate_RNA_n), as.character(multiple_candidate_RNA_n), 
    "ZERO_NOT_APPLIED", "ZERO_NOT_APPLIED", "UNCHANGED_436", "UNCHANGED_23499", "NOT_STARTED", "NOT_PERFORMED", 
    "NONE_CREATED", gate_status), required_action = c("Retain the accepted v2 source audit.", "Do not lock the combined gene-ID policy until this audit is reviewed.", 
    "Retain feature_id as the RNA identifier candidate.", "Preserve all 436 rows and the exact partition.", "Do not alter cistrome mapping membership.", 
    "Review deterministic punctuation candidates.", "Review version stripping only when a locked cistrome symbol exists.", 
    "Confirm each unique candidate manually.", "Resolve any one-to-many candidate ambiguity before lock.", "Do not apply mappings in this substep.", 
    "Do not apply aliases in this substep.", "Membership changes are prohibited.", "Membership changes are prohibited.", 
    "Do not integrate until gene-ID policy is locked.", "Use occupancy-supported candidate terminology after integration QC.", 
    "No Scientific Figure Skill 3 workflow is triggered.", "Upload candidate and classification tables."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3A_v3_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3A v3 punctuation/version audit card", "", paste0("- Run ID: ", run_id), "- RNA source candidate: STEP03B1_v10_RNA_interaction_program_FROZEN.csv", 
    "- RNA gene field: feature_id", "- Exact cistrome matches retained: 279", "- Unmatched symbols audited: 157", 
    paste0("- Dotted unmatched symbols: ", dotted_unmatched_n), paste0("- Versioned AC/AL/AP unmatched symbols: ", 
        versioned_accession_n), paste0("- Symbols with exactly one deterministic candidate: ", single_candidate_RNA_n), 
    paste0("- Symbols with multiple deterministic candidates: ", multiple_candidate_RNA_n), paste0("- Symbols with no deterministic candidate: ", 
        no_candidate_RNA_n), "- Automatic punctuation mapping: ZERO", "- Automatic alias mapping: ZERO", "- RNA membership changed: NO", 
    "- Peak-gene membership changed: NO", "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C3A_v3_PUNCTUATION_VERSION_AUDIT_CARD.md"))

write_text(c("# STEP03B3C3A v3 Methods draft", "", "The 157 RNA symbols lacking exact cistrome matches were audited for deterministic identifier-format alternatives without changing either locked membership set. Candidate cistrome symbols were identified by exact dot-to-hyphen conversion, restricted version stripping for AC/AL/AP-style accessions and a punctuation-insensitive key that removed periods, underscores and hyphens.", 
    "", "All candidate mappings were retained for manual review. No punctuation, version, alias or capitalization mapping was applied automatically, and RNA/cistrome integration was not performed."), 
    file.path(dirs$text, "STEP03B3C3A_v3_Methods_draft.md"))

parameters <- list(schema_version = "3.0", run_id = run_id, RNA_anchor = RNA_anchor, RNA_source_candidate_version = RNA_source_candidate_version, 
    peak_gene_map_lock_version = peak_gene_map_lock_version, exact_match_n = expected_exact_n, unmatched_n = expected_unmatched_n, 
    dotted_unmatched_n = dotted_unmatched_n, versioned_accession_n = versioned_accession_n, single_candidate_RNA_n = single_candidate_RNA_n, 
    multiple_candidate_RNA_n = multiple_candidate_RNA_n, no_candidate_RNA_n = no_candidate_RNA_n, automatic_punctuation_mapping = FALSE, 
    automatic_alias_mapping = FALSE, RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, RNA_CUTRUN_integration = FALSE, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3A_v3_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v3.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("43_step03B3C3A_v3_punctuation_version_candidate_audit.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "resolve_gene_symbol_aliases.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3A_v3_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3A_v3_PUNCTUATION_VERSION_AUDIT_CARD.md"), 
    file.path(audit_root, "STEP03B3C3A_v3_PUNCTUATION_VERSION_CANDIDATES.csv"), file.path(audit_root, "STEP03B3C3A_v3_UNMATCHED_SYMBOL_CLASSIFICATION.csv"), 
    file.path(dirs$qc, "STEP03B3C3A_v3_PUNCTUATION_VERSION_AUDIT_SUMMARY.csv"), file.path(dirs$qc, "STEP03B3C3A_v3_CANDIDATE_METHOD_SUMMARY.csv"), 
    file.path(dirs$text, "STEP03B3C3A_v3_Methods_draft.md"), file.path(dirs$parameters, "STEP03B3C3A_v3_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v3.txt"), file.path(dirs$code, "resolve_gene_symbol_aliases.R"), 
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
    else if (grepl("Methods", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("AUDIT|SUMMARY|CARD|CANDIDATE|CLASSIFICATION", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3A_v3_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3A v3 archive status", "", paste0("- Run ID: ", run_id), "- Exact matches retained: 279", 
    "- Unmatched audited: 157", paste0("- Unique deterministic candidates: ", single_candidate_RNA_n), paste0("- Multiple deterministic candidates: ", 
        multiple_candidate_RNA_n), "- Automatic mappings applied: 0", "- RNA membership changed: NO", "- Peak-gene membership changed: NO", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C3A v3", paste0("- Run ID: ", run_id), paste0("- RNA anchor: ", 
    RNA_anchor, " IMMUTABLE"), paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE"), "- Exact RNA/cistrome matches: 279", 
    "- Punctuation/version candidates: AUDITED_PENDING_MANUAL_REVIEW", "- Automatic alias/punctuation mapping: PROHIBITED", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3A v3 punctuation/version candidate audit completed.\n")

cat("RNA source: STEP03B1_v10_RNA_interaction_program_FROZEN.csv\n")

cat("Gene field: feature_id\n")

cat("Exact cistrome matches retained: 279\n")

cat("Unmatched symbols audited: 157\n")

cat("Dotted unmatched symbols: ", dotted_unmatched_n, "\n")

cat("Versioned AC/AL/AP symbols: ", versioned_accession_n, "\n")

cat("Exactly one deterministic candidate: ", single_candidate_RNA_n, "\n")

cat("Multiple deterministic candidates: ", multiple_candidate_RNA_n, "\n")

cat("No deterministic candidate: ", no_candidate_RNA_n, "\n")

cat("Automatic punctuation mapping: ZERO\n")

cat("Automatic alias mapping: ZERO\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C3A_v3_GATE_TABLE.csv\n")

cat("  08_qc/STEP03B3C3A_v3_PUNCTUATION_VERSION_AUDIT_SUMMARY.csv\n")

cat("  06_results/STEP03B3C3A_v3_punctuation_version_audit/STEP03B3C3A_v3_PUNCTUATION_VERSION_CANDIDATES.csv\n")

cat("  06_results/STEP03B3C3A_v3_punctuation_version_audit/STEP03B3C3A_v3_UNMATCHED_SYMBOL_CLASSIFICATION.csv\n")

cat("============================================================\n")

