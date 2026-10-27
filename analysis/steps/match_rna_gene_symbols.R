#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A15: Match RNA symbols to author cistrome symbols

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3A_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

RNA_source_candidate_version <- "RNA_INTERACTION_SOURCE_v1.0_CANDIDATE"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_RNA_gene_n <- 436L

expected_attenuated_n <- 431L

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

audit_root <- file.path(dirs$results, "STEP03B3C3A_v2_RNA_source_and_gene_ID_audit")

dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3Av2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

valid_symbol_shape <- function(x) {
    grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", x, perl = TRUE)
}

parse_bool <- function(x) {
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

log_event("Starting", run_id)

locked_map_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

required_paths <- c(previous_gate = file.path(dirs$results, "STEP03B3C3A_GATE_TABLE.csv"), previous_inventory = file.path(dirs$qc, 
    "STEP03B3C3A_RNA_SOURCE_CANDIDATE_INVENTORY.csv"), RNA_source = file.path(dirs$results, "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"), 
    RNA_gate = file.path(dirs$results, "STEP03B1_v10_GATE_TABLE.csv"), RNA_lock_card = file.path(dirs$results, 
        "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"), RNA_count_reconciliation = file.path(dirs$results, "STEP03B1_v10_INTERACTION_COUNT_RECONCILIATION.csv"), 
    map_manifest = file.path(locked_map_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"), authorization = file.path(dirs$control, 
        "STEP03B3C3_AUTHORIZATION.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3A v2 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

previous_gate <- .vdr_read_csv(required_paths[["previous_gate"]], stringsAsFactors = FALSE, check.names = FALSE)

previous_inventory <- .vdr_read_csv(required_paths[["previous_inventory"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_gate <- .vdr_read_csv(required_paths[["RNA_gate"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

RNA_count_reconciliation <- .vdr_read_csv(required_paths[["RNA_count_reconciliation"]], stringsAsFactors = FALSE, 
    check.names = FALSE)

map_manifest <- .vdr_read_csv(required_paths[["map_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

authorization_text <- readLines(required_paths[["authorization"]], warn = FALSE, encoding = "UTF-8")

previous_gate_value <- previous_gate$status[previous_gate$gate_component == "STEP03B3C3A gate"]

if (length(previous_gate_value) != 1L || previous_gate_value != "HOLD_NO_AUTHORITATIVE_RNA_SOURCE_CANDIDATE") {
    stop("Unexpected STEP03B3C3A v1 gate: ", paste(previous_gate_value, collapse = ", "))
}

inventory_row <- previous_inventory[previous_inventory$file_name == "STEP03B1_v10_RNA_interaction_program_FROZEN.csv", 
    , drop = FALSE]

if (nrow(inventory_row) != 1L || inventory_row$read_status != "READ_OK" || inventory_row$n_rows != expected_RNA_gene_n) {
    stop("The v1 inventory did not identify the expected readable 436-row frozen RNA file.")
}

if (!is.na(inventory_row$best_gene_column)) {
    stop("The v1 inventory unexpectedly nominated a gene column; this repair is not applicable.")
}

RNA_gate_value <- RNA_gate$status[grepl("STEP03B1.*gate", RNA_gate$gate_component, ignore.case = TRUE)]

if (length(RNA_gate_value) != 1L || RNA_gate_value != "READY_FOR_STEP03B2_ATAC_ANALYSIS_WITH_FIGURE_RELEASE_DEFERRED") {
    stop("STEP03B1 v10 gate is not the accepted RNA lock gate.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("STEP03B1 v10 lock card does not confirm the frozen RNA program.")
}

required_reconciliation <- c(candidate_total = 436L, summary_total = 436L, candidate_attenuated = 431L, summary_attenuated = 431L)

observed_reconciliation <- setNames(RNA_count_reconciliation$value, RNA_count_reconciliation$metric)

missing_reconciliation <- setdiff(names(required_reconciliation), names(observed_reconciliation))

if (length(missing_reconciliation) > 0L || !all(as.integer(observed_reconciliation[names(required_reconciliation)]) == 
    as.integer(required_reconciliation)) || !all(parse_bool(RNA_count_reconciliation$pass))) {
    stop("STEP03B1 v10 interaction-count reconciliation failed.")
}

if (!any(grepl("Authorization: YES", authorization_text, fixed = TRUE))) {
    stop("STEP03B3C3 authorization is incomplete.")
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

RNA_table <- .vdr_read_csv(required_paths[["RNA_source"]], stringsAsFactors = FALSE, check.names = FALSE)

required_RNA_columns <- c("program_freeze_version", "program_freeze_run", "program_freeze_date", "membership_status", 
    "program_scope", "direct_target_tier", "feature_id", "ligand_logFC_siNT", "ligand_FDR_siNT", "ligand_logFC_siVDR", 
    "ligand_FDR_siVDR", "interaction_logFC", "interaction_FDR", "is_interaction_candidate", "program_class")

missing_RNA_columns <- setdiff(required_RNA_columns, names(RNA_table))

if (length(missing_RNA_columns) > 0L) {
    stop("Frozen RNA source is missing columns: ", paste(missing_RNA_columns, collapse = ", "))
}

RNA_symbols <- clean_symbol(RNA_table$feature_id)

if (nrow(RNA_table) != expected_RNA_gene_n || any(is.na(RNA_symbols)) || length(unique(RNA_symbols)) != expected_RNA_gene_n || 
    !all(valid_symbol_shape(RNA_symbols))) {
    stop("feature_id does not define exactly 436 unique valid RNA symbols.")
}

if (!all(RNA_table$program_freeze_version == RNA_anchor) || !all(RNA_table$membership_status == "FROZEN_AT_RNA_STATISTICAL_LAYER")) {
    stop("Frozen RNA metadata is inconsistent with RNA_PROGRAM_v1.0.")
}

if (!all(parse_bool(RNA_table$is_interaction_candidate)) || !all(is.finite(as.numeric(RNA_table$interaction_logFC))) || 
    !all(is.finite(as.numeric(RNA_table$interaction_FDR))) || !all(as.numeric(RNA_table$interaction_FDR) < 0.050000000000000003) || 
    !all(abs(as.numeric(RNA_table$interaction_logFC)) >= 0.5)) {
    stop("Frozen RNA rows do not all satisfy the locked interaction definition.")
}

attenuated_n <- sum(grepl("^attenuated_ligand_", RNA_table$program_class))

if (attenuated_n != expected_attenuated_n) {
    stop("Frozen RNA attenuated-class count differs from 431.")
}

archived_RNA_copies <- list.files(dirs$archive_parent, pattern = "^STEP03B1_v10_RNA_interaction_program_FROZEN\\.csv$", 
    recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

if (length(archived_RNA_copies) == 0L) {
    stop("No archived STEP03B1 v10 frozen RNA copy was found.")
}

current_md5 <- md5_file(required_paths[["RNA_source"]])

current_sha256 <- sha256_file(required_paths[["RNA_source"]])

archive_identity <- data.frame(archive_path = normalize_slashes(archived_RNA_copies), size_bytes = as.numeric(file.info(archived_RNA_copies)$size), 
    md5 = vapply(archived_RNA_copies, md5_file, character(1)), sha256 = vapply(archived_RNA_copies, sha256_file, 
        character(1)), stringsAsFactors = FALSE)

archive_identity$matches_current <- with(archive_identity, md5 == current_md5 & sha256 == current_sha256)

if (!any(archive_identity$matches_current)) {
    stop("No archived frozen RNA copy matches the current result file.")
}

write_csv_utf8(archive_identity, file.path(dirs$qc, "STEP03B3C3A_v2_RNA_ARCHIVE_IDENTITY_AUDIT.csv"))

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

RNA_upper <- toupper(RNA_symbols)

cistrome_upper <- toupper(cistrome_symbols)

RNA_upper_count <- table(RNA_upper)

cistrome_upper_count <- table(cistrome_upper)

RNA_upper_lookup <- as.integer(RNA_upper_count[RNA_upper])

cistrome_upper_lookup <- as.integer(cistrome_upper_count[RNA_upper])

RNA_upper_lookup[is.na(RNA_upper_lookup)] <- 0L

cistrome_upper_lookup[is.na(cistrome_upper_lookup)] <- 0L

exact_match <- RNA_symbols %in% cistrome_symbols

case_insensitive_match <- RNA_upper %in% cistrome_upper

case_only_match <- !exact_match & case_insensitive_match

no_match <- !case_insensitive_match

cistrome_case_candidates <- vapply(RNA_upper, function(symbol_upper) {
    paste(cistrome_symbols[cistrome_upper == symbol_upper], collapse = ";")
}, character(1))

RNA_gene_ID_audit <- data.frame(RNA_row = seq_len(nrow(RNA_table)), RNA_gene_symbol = RNA_symbols, valid_symbol_shape = valid_symbol_shape(RNA_symbols), 
    exact_cistrome_match = exact_match, case_only_cistrome_match = case_only_match, no_cistrome_match = no_match, 
    RNA_case_collision = RNA_upper_lookup > 1L, cistrome_case_collision = cistrome_upper_lookup > 1L, cistrome_case_candidates = cistrome_case_candidates, 
    automatic_alias_mapping = FALSE, stringsAsFactors = FALSE)

case_collision_audit <- RNA_gene_ID_audit[RNA_gene_ID_audit$RNA_case_collision | RNA_gene_ID_audit$cistrome_case_collision | 
    RNA_gene_ID_audit$case_only_cistrome_match, , drop = FALSE]

unmatched_audit <- RNA_gene_ID_audit[RNA_gene_ID_audit$no_cistrome_match, , drop = FALSE]

match_summary <- data.frame(metric = c("RNA rows", "Unique RNA symbols", "Valid RNA symbol shape", "Exact cistrome matches", 
    "Case-only cistrome matches", "No cistrome match", "RNA case-collision rows", "Cistrome case-collision rows", 
    "Cistrome symbol universe", "Automatic alias mappings"), value = c(nrow(RNA_table), length(unique(RNA_symbols)), 
    sum(RNA_gene_ID_audit$valid_symbol_shape), sum(RNA_gene_ID_audit$exact_cistrome_match), sum(RNA_gene_ID_audit$case_only_cistrome_match), 
    sum(RNA_gene_ID_audit$no_cistrome_match), sum(RNA_gene_ID_audit$RNA_case_collision), sum(RNA_gene_ID_audit$cistrome_case_collision), 
    length(cistrome_symbols), 0L), stringsAsFactors = FALSE)

write_csv_utf8(RNA_gene_ID_audit, file.path(audit_root, "STEP03B3C3A_v2_RNA_GENE_ID_AUDIT.csv"))

write_csv_utf8(match_summary, file.path(audit_root, "STEP03B3C3A_v2_RNA_CISTROME_MATCH_SUMMARY.csv"))

write_csv_utf8(case_collision_audit, file.path(audit_root, "STEP03B3C3A_v2_CASE_COLLISION_AUDIT.csv"))

write_csv_utf8(unmatched_audit, file.path(audit_root, "STEP03B3C3A_v2_UNMATCHED_RNA_SYMBOL_AUDIT.csv"))

RNA_source_candidate <- data.frame(source_candidate_version = RNA_source_candidate_version, RNA_anchor = RNA_anchor, 
    source_file = normalize_slashes(required_paths[["RNA_source"]]), source_file_name = basename(required_paths[["RNA_source"]]), 
    size_bytes = as.numeric(file.info(required_paths[["RNA_source"]])$size), md5 = current_md5, sha256 = current_sha256, 
    matching_archive_copy_n = sum(archive_identity$matches_current), n_rows = nrow(RNA_table), n_columns = ncol(RNA_table), 
    gene_symbol_column = "feature_id", unique_gene_n = length(unique(RNA_symbols)), attenuated_class_n = attenuated_n, 
    source_selection_basis = paste("Exact STEP03B1 v10 frozen output;", "436 unique feature_id values;", "RNA_PROGRAM_v1.0 metadata;", 
        "accepted interaction thresholds;", "matching archived copy"), source_status = "CANDIDATE_PENDING_MANUAL_LOCK", 
    stringsAsFactors = FALSE)

write_csv_utf8(RNA_source_candidate, file.path(dirs$results, "STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv"))

collision_n <- nrow(case_collision_audit)

gate_status <- if (collision_n == 0L) {
    "PROCEED_TO_MANUAL_RNA_SOURCE_AND_GENE_ID_LOCK"
} else {
    "HOLD_RNA_CISTROME_CASE_COLLISION_REVIEW"
}

gate_table <- data.frame(gate_component = c("STEP03B3C3A v1 hold", "Root cause", "RNA source", "RNA source archived identity", 
    "RNA gene-symbol column", "RNA row count", "RNA unique-gene count", "RNA interaction definition", "RNA attenuated class count", 
    "Peak-gene map lock", "RNA/cistrome exact-match audit", "Case-only match audit", "Unmatched-symbol audit", 
    "Automatic alias mapping", "RNA membership", "Peak-gene membership", "ATAC role", "RNA/CUT&RUN integration", 
    "Definitive direct-target assignment", "Scientific figures", "STEP03B3C3A v2 gate"), status = c("CONFIRMED_FEATURE_ID_NOT_DISCOVERED", 
    "GENE_COLUMN_NAME_DETECTION_GAP", "STEP03B1_v10_RNA_interaction_program_FROZEN.csv", paste0("PASS_MATCHING_COPIES_", 
        sum(archive_identity$matches_current)), "LOCK_CANDIDATE_feature_id", "PASS_436", "PASS_436", "PASS_FDR_AND_ABS_LOG2FC", 
    "PASS_431", paste0("IMMUTABLE_", peak_gene_map_lock_version), paste0("COMPLETED_EXACT_", sum(exact_match)), 
    paste0("COMPLETED_CASE_ONLY_", sum(case_only_match)), paste0("COMPLETED_UNMATCHED_", sum(no_match)), "ZERO_PROHIBITED", 
    "UNCHANGED_436", "UNCHANGED_23499", ATAC_role, "NOT_STARTED", "NOT_PERFORMED", "NONE_CREATED", gate_status), 
    required_action = c("Retain the v1 hold as an audit record.", "Add feature_id to future RNA gene-column registries.", 
        "Use only after manual source lock.", "Retain current and archived hashes.", "Lock feature_id only after reviewing this audit.", 
        "Do not add or remove RNA rows.", "Do not deduplicate or alias-map.", "Retain the frozen interaction model.", 
        "Retain the reviewed program-class count.", "Do not alter cistrome mapping membership.", "Review exact matched symbols.", 
        "Review case-only rows without automatic correction.", "Retain unmatched rows explicitly.", "Do not silently update aliases.", 
        "Membership changes are prohibited.", "Membership changes are prohibited.", "Use only as optional accessible context.", 
        "Do not integrate until source and ID lock.", "Use occupancy-supported candidate terminology after integration QC.", 
        "No Scientific Figure Skill 3 workflow is triggered.", "Upload source candidate and all gene-ID audit tables."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3A_v2_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3A v2 source and gene-ID audit card", "", paste0("- Run ID: ", run_id), "- v1 hold cause: feature_id was not recognized as a gene-symbol column", 
    paste0("- RNA source candidate: ", basename(required_paths[["RNA_source"]])), "- Gene-symbol column: feature_id", 
    "- RNA rows: 436", "- Unique RNA symbols: 436", "- Attenuated interaction candidates: 431", paste0("- Matching archived copies: ", 
        sum(archive_identity$matches_current)), paste0("- Exact cistrome matches: ", sum(exact_match)), paste0("- Case-only cistrome matches: ", 
        sum(case_only_match)), paste0("- No cistrome match: ", sum(no_match)), paste0("- Case/collision review rows: ", 
        collision_n), "- Automatic alias mapping: ZERO", "- RNA membership changed: NO", "- Peak-gene membership changed: NO", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C3A_v2_RNA_SOURCE_GENE_ID_AUDIT_CARD.md"))

write_text(c("# STEP03B3C3A v2 correction note", "", "The v1 source-discovery audit correctly inventoried the 436-row frozen RNA interaction table but failed to nominate it because the accepted gene-symbol field is named `feature_id`, not `gene` or `symbol`. The v2 audit therefore used the explicit STEP03B1 v10 frozen result and validated `feature_id` against the original interaction-model definition, freeze metadata, count reconciliation and an identical archived copy.", 
    "", "RNA symbols were compared with the locked cistrome universe using exact and separate case-insensitive checks. No alias or capitalization mapping was performed automatically, and RNA/cistrome integration remained deferred."), 
    file.path(dirs$text, "STEP03B3C3A_v2_Correction_note.md"))

parameters <- list(schema_version = "2.0", run_id = run_id, RNA_anchor = RNA_anchor, RNA_source_candidate_version = RNA_source_candidate_version, 
    RNA_source_file = normalize_slashes(required_paths[["RNA_source"]]), gene_symbol_column = "feature_id", RNA_gene_n = expected_RNA_gene_n, 
    attenuated_n = expected_attenuated_n, peak_gene_map_lock_version = peak_gene_map_lock_version, cistrome_gene_n = expected_cistrome_gene_n, 
    exact_match_n = sum(exact_match), case_only_match_n = sum(case_only_match), unmatched_n = sum(no_match), collision_review_n = collision_n, 
    automatic_alias_mapping = FALSE, RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, RNA_CUTRUN_integration = FALSE, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3A_v2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("42_step03B3C3A_v2_use_locked_RNA_v10_feature_id_and_audit.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "match_rna_gene_symbols.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3A_v2_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv"), 
    file.path(dirs$results, "STEP03B3C3A_v2_RNA_SOURCE_GENE_ID_AUDIT_CARD.md"), file.path(dirs$qc, "STEP03B3C3A_v2_RNA_ARCHIVE_IDENTITY_AUDIT.csv"), 
    file.path(audit_root, "STEP03B3C3A_v2_RNA_GENE_ID_AUDIT.csv"), file.path(audit_root, "STEP03B3C3A_v2_RNA_CISTROME_MATCH_SUMMARY.csv"), 
    file.path(audit_root, "STEP03B3C3A_v2_CASE_COLLISION_AUDIT.csv"), file.path(audit_root, "STEP03B3C3A_v2_UNMATCHED_RNA_SYMBOL_AUDIT.csv"), 
    file.path(dirs$text, "STEP03B3C3A_v2_Correction_note.md"), file.path(dirs$parameters, "STEP03B3C3A_v2_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v2.txt"), file.path(dirs$code, "match_rna_gene_symbols.R"), 
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
    else if (grepl("AUDIT|SUMMARY|CARD|IDENTITY|CANDIDATE", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3A_v2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3A v2 archive status", "", paste0("- Run ID: ", run_id), "- RNA source: STEP03B1_v10_RNA_interaction_program_FROZEN.csv", 
    "- Gene-symbol column: feature_id", "- RNA rows/unique symbols: 436/436", paste0("- Exact cistrome matches: ", 
        sum(exact_match)), paste0("- Case-only matches: ", sum(case_only_match)), paste0("- Unmatched symbols: ", 
        sum(no_match)), "- Automatic alias mapping: NO", "- RNA membership changed: NO", "- Peak-gene membership changed: NO", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C3A v2", paste0("- Run ID: ", run_id), paste0("- RNA anchor: ", 
    RNA_anchor, " IMMUTABLE"), "- RNA source: STEP03B1_v10 frozen interaction program CANDIDATE", "- RNA gene-symbol column: feature_id CANDIDATE", 
    paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE"), "- Automatic alias mapping: PROHIBITED", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3A v2 frozen RNA source and gene-ID audit completed.\n")

cat("Root cause: feature_id was not recognized in v1\n")

cat("RNA source: STEP03B1_v10_RNA_interaction_program_FROZEN.csv\n")

cat("Gene-symbol column: feature_id\n")

cat("RNA rows: 436\n")

cat("Unique RNA symbols: 436\n")

cat("Attenuated interaction candidates: 431\n")

cat("Matching archived copies: ", sum(archive_identity$matches_current), "\n")

cat("Exact cistrome matches: ", sum(exact_match), "\n")

cat("Case-only cistrome matches: ", sum(case_only_match), "\n")

cat("No cistrome match: ", sum(no_match), "\n")

cat("Case/collision review rows: ", collision_n, "\n")

cat("Automatic alias mapping: ZERO\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C3A_v2_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv\n")

cat("  08_qc/STEP03B3C3A_v2_RNA_ARCHIVE_IDENTITY_AUDIT.csv\n")

cat("  06_results/STEP03B3C3A_v2_RNA_source_and_gene_ID_audit/\n")

cat("============================================================\n")

