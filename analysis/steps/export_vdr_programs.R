#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A19: Export four program definitions and evidence routes

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3B_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

RNA_source_lock_version <- "RNA_INTERACTION_SOURCE_v1.0"

gene_ID_policy_lock_version <- "RNA_CISTROME_GENE_ID_POLICY_v1.0"

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

integration_candidate_version <- "RNA_OCCUPANCY_EVIDENCE_CANDIDATE_v1.0"

integration_lock_version <- "RNA_OCCUPANCY_EVIDENCE_v1.0"

routing_table_version <- "DOWNSTREAM_CANDIDATE_ROUTING_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

expected_counts <- list(RNA_total = 436L, exact_join = 279L, manual_join = 5L, occupancy_supported = 284L, RNA_only = 152L, 
    relation_join = 855L, vehicle_candidates = 77L, calcitriol_candidates = 275L, shared_candidates = 46L, promoter_candidates = 84L, 
    peak_gene_relations_locked = 23499L, attenuation_candidates = 431L)

expected_evidence_groups <- c(RNA_INTERACTION_ONLY_NO_CURRENT_CISTROME_SUPPORT = 152L, RNA_INTERACTION_PLUS_BOTH_CONDITION_OCCUPANCY_NO_SHARED_PEAK = 22L, 
    RNA_INTERACTION_PLUS_CALCITRIOL_OCCUPANCY = 207L, RNA_INTERACTION_PLUS_SHARED_VDR_OCCUPANCY = 46L, RNA_INTERACTION_PLUS_VEHICLE_OCCUPANCY = 9L)

expected_relation_sources <- c(VDR_CALCITRIOL_CISTROME = 694L, VDR_SHARED_CISTROME = 58L, VDR_VEHICLE_CISTROME = 103L)

expected_relation_source_genes <- c(VDR_CALCITRIOL_CISTROME = 275L, VDR_SHARED_CISTROME = 46L, VDR_VEHICLE_CISTROME = 77L)

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

candidate_root <- file.path(dirs$results, "STEP03B3C3B_RNA_occupancy_integration_candidate")

locked_root <- file.path(dirs$results, "STEP03B3C3B_locked_RNA_occupancy_evidence", integration_lock_version)

dir.create(locked_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3Bv2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

log_event("Starting", run_id)

required_paths <- c(gate = file.path(dirs$results, "STEP03B3C3B_GATE_TABLE.csv"), candidate_manifest = file.path(candidate_root, 
    "STEP03B3C3B_INTEGRATION_CANDIDATE_MANIFEST.csv"), integrated_candidate = file.path(candidate_root, "STEP03B3C3B_RNA_OCCUPANCY_INTEGRATED_CANDIDATE.csv"), 
    relation_long = file.path(candidate_root, "STEP03B3C3B_RNA_OCCUPANCY_RELATION_LONG.tsv.gz"), evidence_groups = file.path(candidate_root, 
        "STEP03B3C3B_EVIDENCE_GROUP_COUNTS.csv"), program_class = file.path(candidate_root, "STEP03B3C3B_PROGRAM_CLASS_BY_EVIDENCE_GROUP.csv"), 
    occupancy_source = file.path(candidate_root, "STEP03B3C3B_OCCUPANCY_SOURCE_SUMMARY.csv"), relation_source = file.path(candidate_root, 
        "STEP03B3C3B_RELATION_SOURCE_SUMMARY.csv"), mapping_status = file.path(candidate_root, "STEP03B3C3B_MAPPING_STATUS_BY_OCCUPANCY.csv"), 
    integration_qc = file.path(dirs$qc, "STEP03B3C3B_INTEGRATION_QC_SUMMARY.csv"), manual_mapping_audit = file.path(dirs$qc, 
        "STEP03B3C3B_MANUAL_MAPPING_JOIN_AUDIT.csv"), RNA_only_audit = file.path(dirs$qc, "STEP03B3C3B_RNA_ONLY_152_AUDIT.csv"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3B output is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate <- .vdr_read_csv(required_paths[["gate"]], stringsAsFactors = FALSE, check.names = FALSE)

candidate_manifest <- .vdr_read_csv(required_paths[["candidate_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

integrated <- .vdr_read_csv(required_paths[["integrated_candidate"]], stringsAsFactors = FALSE, check.names = FALSE)

evidence_groups <- .vdr_read_csv(required_paths[["evidence_groups"]], stringsAsFactors = FALSE, check.names = FALSE)

program_class <- .vdr_read_csv(required_paths[["program_class"]], stringsAsFactors = FALSE, check.names = FALSE)

occupancy_source <- .vdr_read_csv(required_paths[["occupancy_source"]], stringsAsFactors = FALSE, check.names = FALSE)

relation_source <- .vdr_read_csv(required_paths[["relation_source"]], stringsAsFactors = FALSE, check.names = FALSE)

mapping_status <- .vdr_read_csv(required_paths[["mapping_status"]], stringsAsFactors = FALSE, check.names = FALSE)

integration_qc <- .vdr_read_csv(required_paths[["integration_qc"]], stringsAsFactors = FALSE, check.names = FALSE)

manual_mapping_audit <- .vdr_read_csv(required_paths[["manual_mapping_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

RNA_only_audit <- .vdr_read_csv(required_paths[["RNA_only_audit"]], stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- gate$status[gate$gate_component == "STEP03B3C3B gate"]

if (length(accepted_gate) != 1L || accepted_gate != "PROCEED_TO_MANUAL_RNA_OCCUPANCY_INTEGRATION_LOCK") {
    stop("STEP03B3C3B does not authorize integration lock. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (nrow(candidate_manifest) != 6L || !all(candidate_manifest$candidate_version == integration_candidate_version) || 
    !all(candidate_manifest$RNA_membership_status == "IMMUTABLE_436") || !all(candidate_manifest$peak_gene_membership_status == 
    "IMMUTABLE_23499")) {
    stop("Integration candidate manifest metadata failed.")
}

candidate_manifest$file_exists <- file.exists(candidate_manifest$path)

candidate_manifest$observed_size_bytes <- as.numeric(file.info(candidate_manifest$path)$size)

candidate_manifest$observed_md5 <- vapply(candidate_manifest$path, md5_file, character(1))

candidate_manifest$observed_sha256 <- vapply(candidate_manifest$path, sha256_file, character(1))

candidate_manifest$identity_pass <- with(candidate_manifest, file_exists & size_bytes == observed_size_bytes & 
    md5 == observed_md5 & sha256 == observed_sha256)

if (!all(candidate_manifest$identity_pass)) {
    stop("One or more integration candidate artifacts failed hash validation.")
}

write_csv_utf8(candidate_manifest, file.path(dirs$qc, "STEP03B3C3B_v2_CANDIDATE_IDENTITY_RECHECK.csv"))

required_integrated_columns <- c("feature_id", "program_class", "interaction_logFC", "interaction_FDR", "mapping_status", 
    "manual_mapping_accepted", "cistrome_join_authorized", "total_unique_peak_n", "total_peak_gene_relation_n", 
    "vehicle_peak_n", "calcitriol_peak_n", "shared_peak_n", "vehicle_nonshared_peak_n", "calcitriol_nonshared_peak_n", 
    "any_promoter_overlap_support", "occupancy_supported", "shared_occupancy_supported", "vehicle_occupancy_supported", 
    "calcitriol_occupancy_supported", "occupancy_evidence_group", "direct_target_assignment", "causal_claim_authorized", 
    "ATAC_evidence_status")

missing_integrated_columns <- setdiff(required_integrated_columns, names(integrated))

if (length(missing_integrated_columns) > 0L) {
    stop("Integrated table is missing columns: ", paste(missing_integrated_columns, collapse = ", "))
}

boolean_columns <- c("manual_mapping_accepted", "cistrome_join_authorized", "any_promoter_overlap_support", "occupancy_supported", 
    "shared_occupancy_supported", "vehicle_occupancy_supported", "calcitriol_occupancy_supported", "causal_claim_authorized")

for (column_name in boolean_columns) {
    integrated[[column_name]] <- parse_bool(integrated[[column_name]])
}

if (nrow(integrated) != expected_counts$RNA_total || length(unique(integrated$feature_id)) != expected_counts$RNA_total || 
    any(is.na(integrated$feature_id) | !nzchar(trimws(integrated$feature_id)))) {
    stop("Integrated table failed 436-row feature_id identity QC.")
}

if (sum(integrated$mapping_status == "EXACT_AUTHOR_SYMBOL_MATCH") != expected_counts$exact_join || sum(integrated$manual_mapping_accepted) != 
    expected_counts$manual_join || sum(integrated$cistrome_join_authorized) != expected_counts$occupancy_supported || 
    sum(integrated$occupancy_supported) != expected_counts$occupancy_supported || sum(!integrated$occupancy_supported) != 
    expected_counts$RNA_only) {
    stop("Integrated gene-ID or occupancy counts failed reconciliation.")
}

if (any(integrated$cistrome_join_authorized != integrated$occupancy_supported) || any(integrated$direct_target_assignment != 
    "NOT_PERFORMED") || any(integrated$causal_claim_authorized) || any(integrated$ATAC_evidence_status != "NOT_ADDED_RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY")) {
    stop("Integrated interpretation safeguards failed.")
}

if (!all(integrated$shared_peak_n <= integrated$vehicle_peak_n) || !all(integrated$shared_peak_n <= integrated$calcitriol_peak_n) || 
    !all(integrated$vehicle_peak_n == integrated$vehicle_nonshared_peak_n + integrated$shared_peak_n) || !all(integrated$calcitriol_peak_n == 
    integrated$calcitriol_nonshared_peak_n + integrated$shared_peak_n)) {
    stop("Shared/nonshared peak-count reconciliation failed.")
}

expected_group_by_row <- ifelse(integrated$shared_peak_n > 0L, "RNA_INTERACTION_PLUS_SHARED_VDR_OCCUPANCY", ifelse(integrated$vehicle_peak_n > 
    0L & integrated$calcitriol_peak_n > 0L, "RNA_INTERACTION_PLUS_BOTH_CONDITION_OCCUPANCY_NO_SHARED_PEAK", ifelse(integrated$calcitriol_peak_n > 
    0L, "RNA_INTERACTION_PLUS_CALCITRIOL_OCCUPANCY", ifelse(integrated$vehicle_peak_n > 0L, "RNA_INTERACTION_PLUS_VEHICLE_OCCUPANCY", 
    "RNA_INTERACTION_ONLY_NO_CURRENT_CISTROME_SUPPORT"))))

if (!identical(as.character(integrated$occupancy_evidence_group), as.character(expected_group_by_row))) {
    stop("At least one occupancy evidence group is inconsistent with peak counts.")
}

observed_evidence_groups <- setNames(evidence_groups$RNA_candidate_n, evidence_groups$occupancy_evidence_group)

if (!setequal(names(observed_evidence_groups), names(expected_evidence_groups)) || !all(as.integer(observed_evidence_groups[names(expected_evidence_groups)]) == 
    as.integer(expected_evidence_groups)) || sum(evidence_groups$RNA_candidate_n) != expected_counts$RNA_total) {
    stop("Evidence-group counts failed reconciliation.")
}

occupancy_value <- setNames(occupancy_source$RNA_candidate_n, occupancy_source$occupancy_metric)

required_occupancy_values <- c(`RNA candidates with any occupancy support` = 284L, `RNA candidates with vehicle occupancy` = 77L, 
    `RNA candidates with calcitriol occupancy` = 275L, `RNA candidates with shared occupancy` = 46L, `RNA candidates with promoter occupancy` = 84L, 
    `RNA-only candidates without current cistrome support` = 152L, `Manually harmonized candidates with occupancy support` = 5L)

if (!all(names(required_occupancy_values) %in% names(occupancy_value)) || !all(as.integer(occupancy_value[names(required_occupancy_values)]) == 
    as.integer(required_occupancy_values))) {
    stop("Occupancy-source summary failed reconciliation.")
}

observed_relation_sources <- setNames(relation_source$joined_peak_gene_relation_n, relation_source$source_id)

observed_relation_source_genes <- setNames(relation_source$joined_RNA_candidate_n, relation_source$source_id)

if (!setequal(names(observed_relation_sources), names(expected_relation_sources)) || !all(as.integer(observed_relation_sources[names(expected_relation_sources)]) == 
    as.integer(expected_relation_sources)) || !all(as.integer(observed_relation_source_genes[names(expected_relation_source_genes)]) == 
    as.integer(expected_relation_source_genes)) || sum(observed_relation_sources) != expected_counts$relation_join) {
    stop("Relation-source summary failed reconciliation.")
}

attenuation_n <- sum(program_class$RNA_candidate_n[grepl("^attenuated_ligand_", program_class$program_class)])

if (sum(program_class$RNA_candidate_n) != expected_counts$RNA_total || attenuation_n != expected_counts$attenuation_candidates) {
    stop("Program-class by evidence-group summary failed.")
}

if (nrow(manual_mapping_audit) != expected_counts$manual_join || any(manual_mapping_audit$total_unique_peak_n < 
    1L) || nrow(RNA_only_audit) != expected_counts$RNA_only) {
    stop("Manual-mapping or RNA-only audit failed.")
}

qc_value <- setNames(integration_qc$value, integration_qc$metric)

required_qc_values <- c(`Frozen RNA rows` = 436L, `Unique frozen RNA feature_id` = 436L, `Exact gene-ID joins` = 279L, 
    `Manual punctuation joins` = 5L, `Cistrome-authorized RNA candidates` = 284L, `RNA-only candidates` = 152L, 
    `Integrated candidate rows` = 436L, `Unique integrated feature_id` = 436L, `Occupancy-supported RNA candidates` = 284L, 
    `RNA-only integrated rows` = 152L, `Relation-level joined rows` = 855L, `Joined unique RNA candidates` = 284L, 
    `Manual mappings successfully joined` = 5L, `Locked peak-gene relation rows` = 23499L, `Peak-gene membership changed` = 0L, 
    `RNA membership changed` = 0L, `Definitive direct-target assignments` = 0L, `Scientific figures` = 0L)

if (!all(names(required_qc_values) %in% names(qc_value)) || !all(as.integer(qc_value[names(required_qc_values)]) == 
    as.integer(required_qc_values))) {
    stop("Integration-QC summary differs from the reviewed values.")
}

routing_table <- integrated[, c("feature_id", "program_class", "interaction_logFC", "interaction_FDR", "RNA_gene_symbol_cistrome_join", 
    "mapping_status", "occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support", "vehicle_peak_n", 
    "calcitriol_peak_n", "shared_peak_n", "total_unique_peak_n", "occupancy_evidence_group"), drop = FALSE]

routing_table$downstream_route_class <- ifelse(routing_table$shared_occupancy_supported & routing_table$any_promoter_overlap_support, 
    "ROUTE_1_SHARED_AND_PROMOTER_SUPPORT", ifelse(routing_table$shared_occupancy_supported, "ROUTE_2_SHARED_NONPROMOTER_SUPPORT", 
        ifelse(routing_table$occupancy_supported & routing_table$any_promoter_overlap_support, "ROUTE_3_NONSHARED_PROMOTER_SUPPORT", 
            ifelse(routing_table$occupancy_supported, "ROUTE_4_OTHER_OCCUPANCY_SUPPORT", "ROUTE_5_RNA_ONLY"))))

routing_table$downstream_analysis_scope <- ifelse(routing_table$occupancy_supported, "PRIMARY_TEMPORAL_PATIENT_AND_SINGLE_CELL_PROJECTION", 
    "SECONDARY_RNA_PROGRAM_RETENTION_AND_SENSITIVITY")

routing_table$routing_table_version <- routing_table_version

routing_table$direct_target_status <- "NOT_ASSIGNED"

routing_counts <- as.data.frame(table(routing_table$downstream_route_class), stringsAsFactors = FALSE)

names(routing_counts) <- c("downstream_route_class", "RNA_candidate_n")

expected_routing_counts <- c(ROUTE_1_SHARED_AND_PROMOTER_SUPPORT = 21L, ROUTE_2_SHARED_NONPROMOTER_SUPPORT = 25L, 
    ROUTE_3_NONSHARED_PROMOTER_SUPPORT = 63L, ROUTE_4_OTHER_OCCUPANCY_SUPPORT = 175L, ROUTE_5_RNA_ONLY = 152L)

observed_routing_counts <- setNames(routing_counts$RNA_candidate_n, routing_counts$downstream_route_class)

if (!setequal(names(observed_routing_counts), names(expected_routing_counts)) || !all(as.integer(observed_routing_counts[names(expected_routing_counts)]) == 
    as.integer(expected_routing_counts))) {
    stop("Downstream routing counts failed reconciliation.")
}

locked_integrated_path <- file.path(locked_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

locked_relation_path <- file.path(locked_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_relation_long_855.tsv.gz")

locked_routing_path <- file.path(locked_root, "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv")

locked_evidence_counts_path <- file.path(locked_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_group_counts.csv")

locked_routing_counts_path <- file.path(locked_root, "DOWNSTREAM_CANDIDATE_ROUTING_v1.0_counts.csv")

copy_validate(required_paths[["integrated_candidate"]], locked_integrated_path)

copy_validate(required_paths[["relation_long"]], locked_relation_path)

write_csv_utf8(routing_table, locked_routing_path)

write_csv_utf8(evidence_groups, locked_evidence_counts_path)

write_csv_utf8(routing_counts, locked_routing_counts_path)

locked_paths <- c(locked_integrated_path, locked_relation_path, locked_routing_path, locked_evidence_counts_path, 
    locked_routing_counts_path)

locked_manifest <- data.frame(artifact_id = c("RNA_OCCUPANCY_INTEGRATED_436", "RNA_OCCUPANCY_RELATION_LONG_855", 
    "DOWNSTREAM_CANDIDATE_ROUTING_436", "OCCUPANCY_EVIDENCE_GROUP_COUNTS", "DOWNSTREAM_ROUTING_COUNTS"), lock_version = c(integration_lock_version, 
    integration_lock_version, routing_table_version, integration_lock_version, routing_table_version), path = normalize_slashes(locked_paths), 
    size_bytes = as.numeric(file.info(locked_paths)$size), md5 = vapply(locked_paths, md5_file, character(1)), 
    sha256 = vapply(locked_paths, sha256_file, character(1)), RNA_membership_status = "IMMUTABLE_436", peak_gene_membership_status = "IMMUTABLE_23499", 
    direct_target_status = "NOT_ASSIGNED", lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(locked_manifest, file.path(locked_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv"))

lock_summary <- data.frame(lock_component = c("RNA source", "RNA membership", "Gene-ID policy", "Peak-gene map", 
    "Integrated evidence table", "Relation-level subset", "Occupancy-supported RNA candidates", "RNA-only candidates", 
    "Vehicle occupancy candidates", "Calcitriol occupancy candidates", "Shared occupancy candidates", "Promoter occupancy candidates", 
    "Direct-target assignments", "ATAC evidence", "Scientific figures"), locked_value = c(RNA_source_lock_version, 
    "436_IMMUTABLE", gene_ID_policy_lock_version, paste0(peak_gene_map_lock_version, "_23499_IMMUTABLE"), paste0(integration_lock_version, 
        "_436_ROWS"), "855_RELATIONS", "284", "152", "77", "275", "46", "84", "0", "NOT_ADDED", "NONE_CREATED"), 
    lock_status = "LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(lock_summary, file.path(dirs$results, "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK.csv"))

gate_status <- "READY_FOR_COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PROJECTION_PACKAGE"

gate_table <- data.frame(gate_component = c("STEP03B3C3B candidate gate", "Candidate artifact hashes", "Integrated RNA rows", 
    "Occupancy-supported RNA candidates", "RNA-only candidates", "Relation-level subset", "Evidence-group reconciliation", 
    "Program-class reconciliation", "Manual punctuation mappings", "Shared/nonshared peak reconciliation", "RNA membership", 
    "Peak-gene membership", "ATAC role", "Definitive direct-target assignment", "Scientific figures", "STEP03 VDR perturbation/cistrome core", 
    "Combined downstream package authorization", "STEP03B3C3B v2 gate"), status = c("PASS", "PASS_6_OF_6", "LOCKED_436", 
    "LOCKED_284", "LOCKED_152", "LOCKED_855", "PASS_436", "PASS_431_ATTENUATION_OF_436", "PASS_5_OF_5", "PASS", 
    "UNCHANGED_436", "UNCHANGED_23499", ATAC_role, "NOT_PERFORMED", "NONE_CREATED", "CORE_MODULE_LOCKED", "AUTHORIZED_AS_ONE_PACKAGE", 
    gate_status), required_action = c("Retain the candidate audit.", "Use only the locked copies and hashes.", 
    "Keep feature_id as RNA row authority.", "Use as the primary downstream projection set.", "Retain for secondary sensitivity analyses.", 
    "Use for traceable peak-level evidence only.", "Do not reinterpret groups as causal tiers.", "Retain all four RNA interaction classes.", 
    "Use only the five locked manual mappings.", "Keep total, shared and nonshared counts separate.", "Membership changes are prohibited.", 
    "Membership changes are prohibited.", "Use selected-locus accessible context only after target selection.", 
    "Use occupancy-supported candidate terminology.", "Begin figure workflow only after downstream biological projection.", 
    "Do not reopen source/gene-ID/peak mapping without a new version.", "Run temporal, patient and single-cell projection as one coordinated package.", 
    "STEP03B3C3B is locked and the core module is closed."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3B_v2_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3B v2 lock card", "", paste0("- Run ID: ", run_id), paste0("- Integrated evidence: ", 
    integration_lock_version, " LOCKED"), "- Frozen RNA candidates: 436", "- Occupancy-supported RNA candidates: 284", 
    "- RNA-only candidates: 152", "- Relation-level joined rows: 855", "- Vehicle occupancy candidates: 77", "- Calcitriol occupancy candidates: 275", 
    "- Shared occupancy candidates: 46", "- Promoter occupancy candidates: 84", "- Attenuation-class RNA candidates: 431", 
    "- Manual punctuation mappings: 5/5", "- RNA membership changed: NO", "- Peak-gene membership changed: NO", 
    "- Direct-target assignments: 0", "- ATAC evidence: NOT ADDED", "- Scientific figures: NONE", "- STEP03 VDR perturbation/cistrome core: LOCKED", 
    "- Next workflow: one combined temporal, patient and single-cell projection package", paste0("- Gate: ", gate_status)), 
    file.path(dirs$results, "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK_CARD.md"))

write_text(c("# Combined downstream package authorization", "", paste0("- Authorized by: ", run_id), paste0("- Integrated evidence source: ", 
    integration_lock_version, " LOCKED"), paste0("- Routing table: ", routing_table_version, " LOCKED"), "- Authorization: YES", 
    "", "## One coordinated analysis package", "", "1. Temporal/decidualization support using GSE250130 or the accepted temporal source.", 
    "2. Patient projection in GSE165004 with the already established leave-one-slide robustness framework.", "3. External patient validation using GSE92324 and GSE111974 under their valid sample-sharing constraints.", 
    "4. Single-cell localization using GSE183837.", "5. Selected-locus ATAC accessible-context review only after integrated biological candidates are narrowed.", 
    "", "## Primary and secondary scopes", "", "- Primary projection set: 284 occupancy-supported RNA interaction candidates.", 
    "- Higher-context routing: 21 shared-plus-promoter and 25 shared-nonpromoter candidates.", "- Promoter-supported nonshared routing: 63 candidates.", 
    "- Other occupancy-supported routing: 175 candidates.", "- Secondary sensitivity set: 152 RNA-only candidates.", 
    "", "## Prohibited", "", "- Do not reopen the frozen RNA, gene-ID or peak-gene memberships.", "- Do not call the integrated candidates definitive direct targets.", 
    "- Do not use discovery ROC, LASSO, random forest or generic network pharmacology.", "- Do not use ATAC as a genome-wide target-ranking layer.", 
    "- Do not split minor file-format corrections into separate scientific substeps."), file.path(dirs$control, 
    "COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PACKAGE_AUTHORIZATION.md"))

write_text(c("# STEP03B3C3B Results closeout draft", "", "Integration of the frozen 436-gene RNA interaction program with the locked VDR peak-gene map identified 284 RNA candidates with author-annotated VDR occupancy support and retained 152 candidates without current cistrome support. The integrated relation-level subset comprised 855 peak-gene relations.", 
    "", "Among the occupancy-supported candidates, 275 had calcitriol-context occupancy, 77 had vehicle-context occupancy and 46 had shared occupancy. Promoter-associated occupancy was present for 84 candidates. The evidence was interpreted as cross-context convergence because the RNA interaction and VDR occupancy experiments used distinct perturbation designs.", 
    "", "No RNA candidate or peak-gene relation was added or removed, no ATAC evidence was incorporated, and no definitive direct-target assignment was made."), 
    file.path(dirs$text, "STEP03B3C3B_v2_Results_closeout.md"))

write_text(c("# STEP03B3C3B Methods closeout", "", "The accepted 436-row RNA/occupancy integration table was revalidated against its manifest hashes, source-specific relation counts, gene-ID mapping states, occupancy evidence groups and promoter-support fields. Evidence groups were reconstructed directly from vehicle, calcitriol and shared peak counts, with shared coordinates separated from nonshared condition-specific coordinates.", 
    "", "The accepted integrated table and its 855-row relation-level subset were copied into a versioned locked directory with MD5 and SHA-256 hashes. A downstream routing table was created for coordinated temporal, patient and single-cell projection; these routing classes were used for analysis organization and were not interpreted as causal or direct-target tiers."), 
    file.path(dirs$text, "STEP03B3C3B_v2_Methods_closeout.md"))

parameters <- list(schema_version = "2.0-lock", run_id = run_id, integration_lock_version = integration_lock_version, 
    routing_table_version = routing_table_version, RNA_anchor = RNA_anchor, RNA_n = 436, occupancy_supported_n = 284, 
    RNA_only_n = 152, relation_join_n = 855, vehicle_candidate_n = 77, calcitriol_candidate_n = 275, shared_candidate_n = 46, 
    promoter_candidate_n = 84, route_1_shared_promoter_n = 21, route_2_shared_nonpromoter_n = 25, route_3_nonshared_promoter_n = 63, 
    route_4_other_occupancy_n = 175, route_5_RNA_only_n = 152, RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, 
    ATAC_evidence_added = FALSE, definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, 
    core_module_locked = TRUE, combined_downstream_package_authorized = TRUE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3B_v2_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3B_v2.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("46_step03B3C3B_v2_lock_RNA_occupancy_and_handoff.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "export_vdr_programs.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK.csv"), file.path(dirs$results, 
    "STEP03B3C3B_v2_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3B_v2_RNA_OCCUPANCY_EVIDENCE_LOCK_CARD.md"), 
    file.path(dirs$qc, "STEP03B3C3B_v2_CANDIDATE_IDENTITY_RECHECK.csv"), file.path(locked_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv"), 
    locked_integrated_path, locked_relation_path, locked_routing_path, locked_evidence_counts_path, locked_routing_counts_path, 
    file.path(dirs$control, "COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PACKAGE_AUTHORIZATION.md"), file.path(dirs$text, 
        "STEP03B3C3B_v2_Results_closeout.md"), file.path(dirs$text, "STEP03B3C3B_v2_Methods_closeout.md"), file.path(dirs$parameters, 
        "STEP03B3C3B_v2_parameters.json"), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3B_v2.txt"), 
    file.path(dirs$code, "export_vdr_programs.R"), log_path)

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
    else if (grepl("QC|RECHECK|MANIFEST|LOCK|CARD|COUNTS|ROUTING", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3B_v2_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3B v2 archive status", "", paste0("- Run ID: ", run_id), paste0("- Integrated evidence: ", 
    integration_lock_version, " LOCKED"), "- RNA candidates: 436", "- Occupancy-supported: 284", "- RNA-only: 152", 
    "- Relation-level subset: 855", "- RNA membership changed: NO", "- Peak-gene membership changed: NO", "- Direct-target assignments: 0", 
    "- Scientific figures: NONE", "- STEP03 core: LOCKED", "- Combined downstream package: AUTHORIZED", paste0("- Gate: ", 
        gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, 
    "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Completed module: STEP03 VDR perturbation and cistrome core", paste0("- Accepted run: ", 
    run_id), paste0("- Integrated evidence: ", integration_lock_version, " LOCKED"), paste0("- Downstream routing: ", 
    routing_table_version, " LOCKED"), "- RNA membership: RNA_PROGRAM_v1.0 IMMUTABLE; 436", "- Occupancy-supported RNA candidates: 284", 
    "- RNA-only candidates: 152", "- Peak-gene map: PEAK_GENE_MAP_v1.1 IMMUTABLE; 23,499", paste0("- ATAC role: ", 
        ATAC_role), "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", "- Active next workflow: COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PROJECTION_PACKAGE", 
    "- STEP03 core module: LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3B v2 RNA occupancy evidence lock completed.\n")

cat("Integrated evidence: ", integration_lock_version, " LOCKED\n")

cat("RNA candidates: 436\n")

cat("Occupancy-supported RNA candidates: 284\n")

cat("RNA-only candidates: 152\n")

cat("Relation-level joined rows: 855\n")

cat("Vehicle occupancy candidates: 77\n")

cat("Calcitriol occupancy candidates: 275\n")

cat("Shared occupancy candidates: 46\n")

cat("Promoter occupancy candidates: 84\n")

cat("Downstream routing classes: 21 / 25 / 63 / 175 / 152\n")

cat("RNA membership changed: NO\n")

cat("Peak-gene membership changed: NO\n")

cat("Direct-target assignments: 0\n")

cat("ATAC evidence added: NO\n")

cat("Scientific figures: NONE\n")

cat("STEP03 VDR perturbation/cistrome core: LOCKED\n")

cat("Combined downstream package: AUTHORIZED\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n")

cat("============================================================\n")

