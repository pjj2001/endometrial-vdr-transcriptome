#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A17: Export the accepted symbol map

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3A_v4_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

RNA_source_lock_version <- "RNA_INTERACTION_SOURCE_v1.0"

gene_ID_policy_lock_version <- "RNA_CISTROME_GENE_ID_POLICY_v1.0"

gene_ID_map_lock_version <- "RNA_CISTROME_GENE_ID_MAP_v1.0"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_RNA_n <- 436L

expected_exact_n <- 279L

expected_manual_mapping_n <- 5L

expected_join_authorized_n <- 284L

expected_unmatched_n <- 152L

expected_peak_gene_relation_n <- 23499L

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

locked_root <- file.path(dirs$results, "STEP03B3C3A_locked_RNA_gene_ID", gene_ID_policy_lock_version)

dir.create(locked_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3Av4_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

parse_bool <- function(x) {
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

clean_symbol <- function(x) {
    value <- trimws(as.character(x))
    value[is.na(value) | !nzchar(value) | toupper(value) %in% c("NA", "N/A", "NULL", "NONE", ".")] <- NA_character_
    value
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

log_event("Starting", run_id)

v2_audit_root <- file.path(dirs$results, "STEP03B3C3A_v2_RNA_source_and_gene_ID_audit")

v3_audit_root <- file.path(dirs$results, "STEP03B3C3A_v3_punctuation_version_audit")

map_locked_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

required_paths <- c(v2_gate = file.path(dirs$results, "STEP03B3C3A_v2_GATE_TABLE.csv"), v2_source_candidate = file.path(dirs$results, 
    "STEP03B3C3A_v2_RNA_INTERACTION_SOURCE_CANDIDATE.csv"), v2_gene_audit = file.path(v2_audit_root, "STEP03B3C3A_v2_RNA_GENE_ID_AUDIT.csv"), 
    v2_match_summary = file.path(v2_audit_root, "STEP03B3C3A_v2_RNA_CISTROME_MATCH_SUMMARY.csv"), v3_gate = file.path(dirs$results, 
        "STEP03B3C3A_v3_GATE_TABLE.csv"), v3_candidates = file.path(v3_audit_root, "STEP03B3C3A_v3_PUNCTUATION_VERSION_CANDIDATES.csv"), 
    v3_classification = file.path(v3_audit_root, "STEP03B3C3A_v3_UNMATCHED_SYMBOL_CLASSIFICATION.csv"), RNA_source = file.path(dirs$results, 
        "STEP03B1_v10_RNA_interaction_program_FROZEN.csv"), map_manifest = file.path(map_locked_root, "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3A v4 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

v2_gate <- .vdr_read_csv(required_paths[["v2_gate"]], stringsAsFactors = FALSE, check.names = FALSE)

source_candidate <- .vdr_read_csv(required_paths[["v2_source_candidate"]], stringsAsFactors = FALSE, check.names = FALSE)

gene_audit <- .vdr_read_csv(required_paths[["v2_gene_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

match_summary <- .vdr_read_csv(required_paths[["v2_match_summary"]], stringsAsFactors = FALSE, check.names = FALSE)

v3_gate <- .vdr_read_csv(required_paths[["v3_gate"]], stringsAsFactors = FALSE, check.names = FALSE)

v3_candidates <- .vdr_read_csv(required_paths[["v3_candidates"]], stringsAsFactors = FALSE, check.names = FALSE)

v3_classification <- .vdr_read_csv(required_paths[["v3_classification"]], stringsAsFactors = FALSE, check.names = FALSE)

map_manifest <- .vdr_read_csv(required_paths[["map_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

v2_gate_value <- v2_gate$status[v2_gate$gate_component == "STEP03B3C3A v2 gate"]

v3_gate_value <- v3_gate$status[v3_gate$gate_component == "STEP03B3C3A v3 gate"]

if (length(v2_gate_value) != 1L || v2_gate_value != "PROCEED_TO_MANUAL_RNA_SOURCE_AND_GENE_ID_LOCK") {
    stop("STEP03B3C3A v2 gate is not accepted.")
}

if (length(v3_gate_value) != 1L || v3_gate_value != "PROCEED_TO_MANUAL_PUNCTUATION_VERSION_CANDIDATE_REVIEW") {
    stop("STEP03B3C3A v3 gate is not accepted.")
}

if (nrow(source_candidate) != 1L || source_candidate$gene_symbol_column != "feature_id" || source_candidate$n_rows != 
    expected_RNA_n || source_candidate$unique_gene_n != expected_RNA_n) {
    stop("RNA source candidate metadata failed.")
}

current_source_md5 <- md5_file(required_paths[["RNA_source"]])

current_source_sha256 <- sha256_file(required_paths[["RNA_source"]])

if (source_candidate$md5 != current_source_md5 || source_candidate$sha256 != current_source_sha256) {
    stop("RNA source candidate hash differs from the frozen source file.")
}

RNA_table <- .vdr_read_csv(required_paths[["RNA_source"]], stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(RNA_table) != expected_RNA_n || !"feature_id" %in% names(RNA_table)) {
    stop("Frozen RNA source failed row or feature_id validation.")
}

RNA_symbols <- clean_symbol(RNA_table$feature_id)

if (any(is.na(RNA_symbols)) || length(unique(RNA_symbols)) != expected_RNA_n) {
    stop("Frozen RNA feature_id values are not 436 unique nonmissing symbols.")
}

if (nrow(gene_audit) != expected_RNA_n || length(unique(gene_audit$RNA_row)) != expected_RNA_n || !identical(as.character(gene_audit$RNA_gene_symbol), 
    as.character(RNA_symbols))) {
    stop("The v2 gene audit does not exactly track the frozen RNA source.")
}

partition_n <- rowSums(gene_audit[c("exact_cistrome_match", "case_only_cistrome_match", "no_cistrome_match")])

if (!all(partition_n == 1L) || sum(parse_bool(gene_audit$exact_cistrome_match)) != expected_exact_n || sum(parse_bool(gene_audit$case_only_cistrome_match)) != 
    0L || sum(parse_bool(gene_audit$no_cistrome_match)) != 157L) {
    stop("The v2 exact/case-only/unmatched partition failed.")
}

manual_decisions <- data.frame(RNA_gene_symbol_original = c("DPP4.DT", "HHIP.AS1", "KRTAP4.8", "MUC12.AS1", "RGMB.AS1"), 
    cistrome_gene_symbol_harmonized = c("DPP4-DT", "HHIP-AS1", "KRTAP4-8", "MUC12-AS1", "RGMB-AS1"), decision = "ACCEPT", 
    mapping_rule = "MANUAL_ACCEPTED_DOT_TO_HYPHEN_OFFICIAL_SYMBOL_EQUIVALENT", local_candidate_evidence = "EXACT_DOT_TO_HYPHEN_AND_UNIQUE_PUNCTUATION_KEY", 
    external_nomenclature_review = "REVIEWED_AGAINST_CURRENT_HGNC_NCBI_SYMBOL_RECORDS", automatic_mapping_applied = FALSE, 
    lock_status = "LOCKED_MANUAL_DECISION", stringsAsFactors = FALSE)

observed_unique_pairs <- unique(v3_candidates[c("RNA_gene_symbol", "candidate_cistrome_symbol")])

observed_unique_pairs <- observed_unique_pairs[order(observed_unique_pairs$RNA_gene_symbol), , drop = FALSE]

expected_pairs <- manual_decisions[c("RNA_gene_symbol_original", "cistrome_gene_symbol_harmonized")]

expected_pairs <- expected_pairs[order(expected_pairs$RNA_gene_symbol_original), , drop = FALSE]

if (nrow(observed_unique_pairs) != expected_manual_mapping_n || !identical(as.character(observed_unique_pairs$RNA_gene_symbol), 
    as.character(expected_pairs$RNA_gene_symbol_original)) || !identical(as.character(observed_unique_pairs$candidate_cistrome_symbol), 
    as.character(expected_pairs$cistrome_gene_symbol_harmonized))) {
    stop("The v3 deterministic candidate pairs differ from the reviewed five mappings.")
}

candidate_method_check <- aggregate(candidate_method ~ RNA_gene_symbol + candidate_cistrome_symbol, data = v3_candidates, 
    FUN = function(x) {
        paste(sort(unique(x)), collapse = ";")
    })

required_method_string <- paste(sort(c("EXACT_DOT_TO_HYPHEN", "PUNCTUATION_INSENSITIVE_KEY")), collapse = ";")

if (nrow(candidate_method_check) != expected_manual_mapping_n || !all(candidate_method_check$candidate_method == 
    required_method_string)) {
    stop("Each reviewed mapping must be supported by both deterministic methods.")
}

reviewed_classification <- v3_classification[v3_classification$RNA_gene_symbol %in% manual_decisions$RNA_gene_symbol_original, 
    , drop = FALSE]

if (nrow(reviewed_classification) != expected_manual_mapping_n || !all(reviewed_classification$candidate_symbol_n == 
    1L) || !all(reviewed_classification$review_class == "ONE_UNIQUE_DETERMINISTIC_CANDIDATE")) {
    stop("The five reviewed mappings are not unique deterministic candidates.")
}

if (nrow(map_manifest) != 2L || !all(map_manifest$lock_version == peak_gene_map_lock_version) || !all(map_manifest$lock_status == 
    "LOCKED")) {
    stop("PEAK_GENE_MAP_v1.1 manifest failed.")
}

map_manifest$file_exists <- file.exists(map_manifest$path)

map_manifest$observed_md5 <- vapply(map_manifest$path, md5_file, character(1))

map_manifest$observed_sha256 <- vapply(map_manifest$path, sha256_file, character(1))

map_manifest$identity_pass <- with(map_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(map_manifest$identity_pass)) {
    stop("A locked peak-gene artifact failed hash validation.")
}

map_long_path <- map_manifest$path[map_manifest$artifact_id == "PEAK_GENE_LONG"]

if (length(map_long_path) != 1L) {
    stop("Locked peak-gene long map could not be uniquely identified.")
}

peak_gene_map <- .vdr_fread(map_long_path, sep = "\t", data.table = FALSE, showProgress = FALSE)

if (nrow(peak_gene_map) != expected_peak_gene_relation_n) {
    stop("Locked peak-gene relation count differs from 23,499.")
}

cistrome_symbols <- sort(unique(clean_symbol(peak_gene_map$gene_symbol_author)))

cistrome_symbols <- cistrome_symbols[!is.na(cistrome_symbols)]

if (length(cistrome_symbols) != expected_cistrome_gene_n) {
    stop("Locked cistrome gene-symbol universe differs from 10,042.")
}

if (!all(manual_decisions$cistrome_gene_symbol_harmonized %in% cistrome_symbols)) {
    stop("At least one accepted harmonized symbol is absent from the locked cistrome universe.")
}

gene_ID_map <- data.frame(RNA_row = gene_audit$RNA_row, RNA_gene_symbol_original = gene_audit$RNA_gene_symbol, 
    RNA_gene_symbol_cistrome_join = NA_character_, mapping_status = NA_character_, mapping_rule = NA_character_, 
    exact_original_match = parse_bool(gene_audit$exact_cistrome_match), manual_mapping_accepted = FALSE, cistrome_join_authorized = FALSE, 
    RNA_membership_status = "FROZEN_IN_RNA_PROGRAM_v1.0", automatic_alias_mapping = FALSE, automatic_punctuation_mapping = FALSE, 
    stringsAsFactors = FALSE)

exact_index <- gene_ID_map$exact_original_match

gene_ID_map$RNA_gene_symbol_cistrome_join[exact_index] <- gene_ID_map$RNA_gene_symbol_original[exact_index]

gene_ID_map$mapping_status[exact_index] <- "EXACT_AUTHOR_SYMBOL_MATCH"

gene_ID_map$mapping_rule[exact_index] <- "IDENTITY"

gene_ID_map$cistrome_join_authorized[exact_index] <- TRUE

for (decision_index in seq_len(nrow(manual_decisions))) {
    original_symbol <- manual_decisions$RNA_gene_symbol_original[decision_index]
    harmonized_symbol <- manual_decisions$cistrome_gene_symbol_harmonized[decision_index]
    map_index <- gene_ID_map$RNA_gene_symbol_original == original_symbol
    if (sum(map_index) != 1L) {
        stop("Manual decision does not map to exactly one RNA row: ", original_symbol)
    }
    gene_ID_map$RNA_gene_symbol_cistrome_join[map_index] <- harmonized_symbol
    gene_ID_map$mapping_status[map_index] <- "MANUAL_ACCEPTED_PUNCTUATION_EQUIVALENT"
    gene_ID_map$mapping_rule[map_index] <- "DOT_TO_HYPHEN_OFFICIAL_SYMBOL_EQUIVALENT"
    gene_ID_map$manual_mapping_accepted[map_index] <- TRUE
    gene_ID_map$cistrome_join_authorized[map_index] <- TRUE
}

unmatched_index <- is.na(gene_ID_map$mapping_status)

gene_ID_map$mapping_status[unmatched_index] <- "NO_CURRENT_CISTROME_SYMBOL_MATCH"

gene_ID_map$mapping_rule[unmatched_index] <- "NONE"

if (nrow(gene_ID_map) != expected_RNA_n || sum(gene_ID_map$mapping_status == "EXACT_AUTHOR_SYMBOL_MATCH") != expected_exact_n || 
    sum(gene_ID_map$manual_mapping_accepted) != expected_manual_mapping_n || sum(gene_ID_map$cistrome_join_authorized) != 
    expected_join_authorized_n || sum(gene_ID_map$mapping_status == "NO_CURRENT_CISTROME_SYMBOL_MATCH") != expected_unmatched_n) {
    stop("The 436-row gene-ID map failed count reconciliation.")
}

authorized_join_symbols <- gene_ID_map$RNA_gene_symbol_cistrome_join[gene_ID_map$cistrome_join_authorized]

if (any(is.na(authorized_join_symbols)) || length(unique(authorized_join_symbols)) != expected_join_authorized_n || 
    !all(authorized_join_symbols %in% cistrome_symbols)) {
    stop("Authorized cistrome join symbols failed completeness, uniqueness or universe validation.")
}

unmatched_locked <- gene_ID_map[!gene_ID_map$cistrome_join_authorized, , drop = FALSE]

RNA_source_lock <- data.frame(lock_version = RNA_source_lock_version, RNA_anchor = RNA_anchor, source_file = normalize_slashes(required_paths[["RNA_source"]]), 
    source_file_name = basename(required_paths[["RNA_source"]]), gene_symbol_column = "feature_id", n_rows = expected_RNA_n, 
    unique_gene_n = expected_RNA_n, size_bytes = as.numeric(file.info(required_paths[["RNA_source"]])$size), md5 = current_source_md5, 
    sha256 = current_source_sha256, membership_status = "IMMUTABLE_436", lock_status = "LOCKED", stringsAsFactors = FALSE)

gene_ID_decision_lock <- manual_decisions

gene_ID_decision_lock$policy_lock_version <- gene_ID_policy_lock_version

gene_ID_map_path <- file.path(locked_root, "RNA_CISTROME_GENE_ID_MAP_v1.0.csv")

unmatched_locked_path <- file.path(locked_root, "RNA_CISTROME_UNMATCHED_152_v1.0.csv")

RNA_source_lock_path <- file.path(locked_root, "RNA_INTERACTION_SOURCE_v1.0_LOCK.csv")

gene_ID_decision_lock_path <- file.path(locked_root, "RNA_CISTROME_GENE_ID_MANUAL_DECISIONS_v1.0.csv")

write_csv_utf8(gene_ID_map, gene_ID_map_path)

write_csv_utf8(unmatched_locked, unmatched_locked_path)

write_csv_utf8(RNA_source_lock, RNA_source_lock_path)

write_csv_utf8(gene_ID_decision_lock, gene_ID_decision_lock_path)

locked_artifact_paths <- c(RNA_source_lock_path, gene_ID_map_path, unmatched_locked_path, gene_ID_decision_lock_path)

locked_manifest <- data.frame(artifact_id = c("RNA_INTERACTION_SOURCE_LOCK", "RNA_CISTROME_GENE_ID_MAP", "RNA_CISTROME_UNMATCHED_SET", 
    "RNA_CISTROME_MANUAL_DECISIONS"), lock_version = c(RNA_source_lock_version, gene_ID_map_lock_version, gene_ID_policy_lock_version, 
    gene_ID_policy_lock_version), path = normalize_slashes(locked_artifact_paths), size_bytes = as.numeric(file.info(locked_artifact_paths)$size), 
    md5 = vapply(locked_artifact_paths, md5_file, character(1)), sha256 = vapply(locked_artifact_paths, sha256_file, 
        character(1)), lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(locked_manifest, file.path(locked_root, "RNA_CISTROME_GENE_ID_POLICY_v1.0_LOCKED_MANIFEST.csv"))

policy_summary <- data.frame(policy_component = c("RNA source", "RNA identifier", "RNA membership", "Exact identity joins", 
    "Manual punctuation joins", "Total cistrome-authorized RNA symbols", "RNA-only unmatched symbols", "Automatic aliases", 
    "Automatic capitalization changes", "Automatic punctuation changes", "Peak-gene membership", "ATAC role", "Direct-target assignment", 
    "Scientific figures"), locked_value = c(RNA_source_lock_version, "feature_id", "436_IMMUTABLE", "279", "5", 
    "284", "152", "0", "0", "0", "23499_IMMUTABLE", ATAC_role, "NOT_PERFORMED", "NONE_CREATED"), lock_status = "LOCKED", 
    stringsAsFactors = FALSE)

write_csv_utf8(policy_summary, file.path(dirs$results, "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK.csv"))

gate_status <- "READY_FOR_STEP03B3C3B_RNA_OCCUPANCY_EVIDENCE_INTEGRATION"

gate_table <- data.frame(gate_component = c("STEP03B3C3A v2 source audit", "STEP03B3C3A v3 punctuation audit", 
    "RNA source identity", "RNA identifier field", "RNA membership", "Exact gene-ID joins", "Manual punctuation decisions", 
    "Manual mapping uniqueness", "Cistrome-authorized RNA symbols", "RNA-only unmatched symbols", "Automatic alias mapping", 
    "Automatic capitalization mapping", "Automatic punctuation mapping", "Peak-gene map", "ATAC role", "RNA/CUT&RUN integration", 
    "Definitive direct-target assignment", "Scientific figures", "STEP03B3C3B authorization", "STEP03B3C3A v4 gate"), 
    status = c("PASS", "PASS", paste0("LOCKED_", RNA_source_lock_version), "LOCKED_feature_id", "IMMUTABLE_436", 
        "LOCKED_279", "LOCKED_5", "PASS_5_OF_5", "LOCKED_284", "LOCKED_152", "ZERO", "ZERO", "ZERO_AUTOMATIC", 
        paste0("IMMUTABLE_", peak_gene_map_lock_version, "_23499"), ATAC_role, "NOT_STARTED", "NOT_PERFORMED", 
        "NONE_CREATED", "AUTHORIZED", gate_status), required_action = c("Retain v2 source audit.", "Retain v3 deterministic candidate audit.", 
        "Use only the locked source and hash.", "Do not replace feature_id in the RNA source.", "Do not add or remove RNA rows.", 
        "Join exact symbols by identity.", "Use only the five locked manual decisions.", "Do not expand to broad punctuation conversion.", 
        "Use the locked harmonized join field.", "Retain as RNA candidates without current occupancy support.", 
        "Broad alias mapping remains prohibited.", "Case-only mapping remains prohibited unless separately reviewed.", 
        "Only five manual mappings are authorized.", "Do not alter peak-gene membership.", "Use only as optional accessible context.", 
        "Begin only in STEP03B3C3B.", "Use occupancy-supported candidate terminology after integration QC.", "No Scientific Figure Skill 3 workflow is triggered.", 
        "Proceed to gene-level occupancy evidence integration.", "STEP03B3C3A is substep-locked."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3A_v4_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3A v4 RNA source and gene-ID policy lock card", "", paste0("- Run ID: ", run_id), paste0("- RNA source: ", 
    RNA_source_lock_version, " LOCKED"), "- RNA identifier: feature_id LOCKED", "- RNA membership: 436 IMMUTABLE", 
    "- Exact cistrome joins: 279", "- Manually accepted punctuation equivalents: 5", "- Total cistrome-authorized RNA symbols: 284", 
    "- RNA-only unmatched symbols: 152", "- Automatic alias mappings: 0", "- Automatic capitalization mappings: 0", 
    "- Automatic punctuation mappings: 0", paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE; 23,499 relations"), 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    "- STEP03B3C3B authorization: YES", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK_CARD.md"))

write_text(c("# STEP03B3C3B authorization", "", paste0("- Authorized by: ", run_id), paste0("- RNA source: ", RNA_source_lock_version, 
    " LOCKED"), paste0("- RNA/cistrome gene-ID policy: ", gene_ID_policy_lock_version, " LOCKED"), paste0("- Gene-ID map: ", 
    gene_ID_map_lock_version, " LOCKED"), paste0("- Peak-gene map: ", peak_gene_map_lock_version, " LOCKED"), "- Authorization: YES", 
    "- Authorized task 1: join the 284 cistrome-authorized RNA candidates to the locked peak-gene map.", "- Authorized task 2: summarize vehicle, calcitriol and shared occupancy separately.", 
    "- Authorized task 3: retain supporting-peak counts and promoter evidence separately.", "- Authorized task 4: retain the 152 unmatched RNA candidates as RNA-only rows.", 
    "- Authorized task 5: construct occupancy-supported candidate evidence tiers without changing RNA membership.", 
    "- Prohibited: broad alias remapping, additional punctuation conversion or nearest-gene invention.", "- Prohibited: definitive direct-target or causal wording.", 
    "- ATAC: optional accessible_context only after CUT&RUN-supported target selection.", "- Scientific figures: only after the integrated evidence table passes content QC."), 
    file.path(dirs$control, "STEP03B3C3B_AUTHORIZATION.md"))

write_text(c("# STEP03B3C3A v4 Methods addendum", "", "The frozen RNA interaction source and `feature_id` identifier field were locked after validation against the accepted RNA program and archived source identity. Of 436 RNA interaction candidates, 279 matched the locked cistrome gene-symbol universe exactly.", 
    "", "Five additional one-to-one punctuation equivalents were accepted after deterministic dot-to-hyphen comparison and manual nomenclature review: DPP4.DT/DPP4-DT, HHIP.AS1/HHIP-AS1, KRTAP4.8/KRTAP4-8, MUC12.AS1/MUC12-AS1 and RGMB.AS1/RGMB-AS1. Original RNA identifiers were retained, and separate harmonized cistrome-join symbols were created. No broad alias or automatic punctuation mapping was applied."), 
    file.path(dirs$text, "STEP03B3C3A_v4_Methods_addendum.md"))

write_text(c("# STEP03B3C3A v4 Results draft", "", "The locked RNA/cistrome gene-ID map retained all 436 RNA interaction candidates. A total of 279 candidates had exact symbol matches and five had manually accepted one-to-one punctuation equivalents, yielding 284 RNA candidates authorized for cistrome joining. The remaining 152 candidates were retained as RNA-only candidates without current author-annotated occupancy support.", 
    "", "No automatic alias, capitalization or broad punctuation mapping was performed, and the locked peak-gene map remained unchanged."), 
    file.path(dirs$text, "STEP03B3C3A_v4_Results_draft.md"))

parameters <- list(schema_version = "4.0-lock", run_id = run_id, RNA_anchor = RNA_anchor, RNA_source_lock_version = RNA_source_lock_version, 
    gene_ID_policy_lock_version = gene_ID_policy_lock_version, gene_ID_map_lock_version = gene_ID_map_lock_version, 
    peak_gene_map_lock_version = peak_gene_map_lock_version, RNA_n = expected_RNA_n, exact_join_n = expected_exact_n, 
    manual_punctuation_join_n = expected_manual_mapping_n, cistrome_join_authorized_n = expected_join_authorized_n, 
    unmatched_RNA_only_n = expected_unmatched_n, automatic_alias_mapping = FALSE, automatic_capitalization_mapping = FALSE, 
    automatic_punctuation_mapping = FALSE, RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, 
    RNA_CUTRUN_integration = FALSE, definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, 
    STEP03B3C3B_authorized = TRUE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3A_v4_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v4.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("44_step03B3C3A_v4_lock_RNA_source_and_gene_ID_policy.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "export_gene_symbol_map.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3C3A_v4_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3A_v4_RNA_SOURCE_GENE_ID_POLICY_LOCK_CARD.md"), 
    RNA_source_lock_path, gene_ID_map_path, unmatched_locked_path, gene_ID_decision_lock_path, file.path(locked_root, 
        "RNA_CISTROME_GENE_ID_POLICY_v1.0_LOCKED_MANIFEST.csv"), file.path(dirs$control, "STEP03B3C3B_AUTHORIZATION.md"), 
    file.path(dirs$text, "STEP03B3C3A_v4_Methods_addendum.md"), file.path(dirs$text, "STEP03B3C3A_v4_Results_draft.md"), 
    file.path(dirs$parameters, "STEP03B3C3A_v4_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A_v4.txt"), 
    file.path(dirs$code, "export_gene_symbol_map.R"), log_path)

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
    else if (grepl("Methods|Results|AUTHORIZATION", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("LOCK|MANIFEST|CARD|DECISION|UNMATCHED", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3A_v4_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3A v4 archive status", "", paste0("- Run ID: ", run_id), paste0("- RNA source: ", RNA_source_lock_version, 
    " LOCKED"), paste0("- Gene-ID policy: ", gene_ID_policy_lock_version, " LOCKED"), "- RNA membership: 436 IMMUTABLE", 
    "- Exact joins: 279", "- Manual punctuation joins: 5", "- Total cistrome-authorized: 284", "- RNA-only unmatched: 152", 
    "- Automatic mappings: 0", "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", "- STEP03B3C3B authorization: YES", 
    paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, "ctl", 
    "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Completed substep: STEP03B3C3A", paste0("- Accepted run: ", run_id), 
    paste0("- RNA source: ", RNA_source_lock_version, " LOCKED"), paste0("- RNA/cistrome gene-ID policy: ", gene_ID_policy_lock_version, 
        " LOCKED"), paste0("- Gene-ID map: ", gene_ID_map_lock_version, " LOCKED"), "- RNA membership: RNA_PROGRAM_v1.0 IMMUTABLE; 436", 
    "- Cistrome-authorized RNA symbols: 284", "- RNA-only unmatched symbols: 152", paste0("- Peak-gene map: ", 
        peak_gene_map_lock_version, " IMMUTABLE"), paste0("- ATAC role: ", ATAC_role), "- RNA/CUT&RUN integration: NOT STARTED", 
    "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", "- Active next substep: STEP03B3C3B", 
    "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3A v4 RNA source and gene-ID policy lock completed.\n")

cat("RNA source: ", RNA_source_lock_version, " LOCKED\n")

cat("Gene-ID policy: ", gene_ID_policy_lock_version, " LOCKED\n")

cat("RNA identifier: feature_id\n")

cat("RNA membership: 436 IMMUTABLE\n")

cat("Exact cistrome joins: 279\n")

cat("Manual punctuation joins: 5\n")

cat("Total cistrome-authorized RNA symbols: 284\n")

cat("RNA-only unmatched symbols: 152\n")

cat("Automatic alias mappings: 0\n")

cat("Automatic punctuation mappings: 0\n")

cat("Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE; 23,499 relations\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("STEP03B3C3B authorization: YES\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n")

cat("============================================================\n")

