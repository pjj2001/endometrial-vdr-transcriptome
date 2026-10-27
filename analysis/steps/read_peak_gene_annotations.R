#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A11: Read author peak-gene annotation fields

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C1_", format(Sys.time(), "%Y%m%d_%H%M%S"))

source_lock_version <- "PUBLISHED_VDR_CISTROME_SOURCE_v1.0"

RNA_anchor <- "RNA_PROGRAM_v1.0"

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

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

preview_root <- file.path(dirs$results, "STEP03B3C1_peak_annotation_previews")

dir.create(preview_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C1_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

nonempty_values <- function(x) {
    value <- trimws(as.character(x))
    value[!is.na(value) & nzchar(value) & !tolower(value) %in% c("na", "n/a", "nan", "null", "none")]
}

collapse_examples <- function(x, n = 6L) {
    value <- unique(nonempty_values(x))
    if (length(value) == 0L) {
        return("")
    }
    paste(utils::head(value, n), collapse = " || ")
}

infer_field_role <- function(column_name) {
    name_lower <- tolower(column_name)
    if (name_lower %in% c("source_id", "source_row", "chromosome", "chromosome_source", "start0", "end", "width", 
        "valid_coordinate", "primary_chromosome", "coordinate_key")) {
        return("LOCKED_COORDINATE_OR_PROVENANCE")
    }
    if (grepl("promoter", name_lower, perl = TRUE)) {
        return("CANDIDATE_PROMOTER_FIELD")
    }
    if (grepl("dist.*tss|tss.*dist|nearest.*tss", name_lower, perl = TRUE)) {
        return("CANDIDATE_TSS_DISTANCE_FIELD")
    }
    if (grepl("nearest.*gene|gene.*nearest", name_lower, perl = TRUE)) {
        return("CANDIDATE_NEAREST_GENE_FIELD")
    }
    if (grepl("geneoverlap|gene_overlap|overlap.*gene|gene.*overlap", name_lower, perl = TRUE)) {
        return("CANDIDATE_GENE_OVERLAP_FIELD")
    }
    if (grepl("gene|symbol|entrez|ensembl", name_lower, perl = TRUE)) {
        return("CANDIDATE_GENE_IDENTIFIER_FIELD")
    }
    if (grepl("annot|region|feature|genic|intergenic|intron|exon|utr", name_lower, perl = TRUE)) {
        return("CANDIDATE_GENOMIC_ANNOTATION_FIELD")
    }
    if (grepl("score|signal|fold|enrich|pileup|pvalue|qvalue|fdr|summit", name_lower, perl = TRUE)) {
        return("CANDIDATE_PEAK_SCORE_FIELD")
    }
    if (grepl("motif|homer", name_lower, perl = TRUE)) {
        return("CANDIDATE_MOTIF_FIELD")
    }
    "OTHER_RAW_FIELD"
}

delimiter_audit <- function(x) {
    value <- nonempty_values(x)
    if (length(value) == 0L) {
        return(data.frame(nonempty_n = 0L, semicolon_n = 0L, comma_n = 0L, pipe_n = 0L, slash_n = 0L, whitespace_n = 0L, 
            multi_delimiter_fraction = NA_real_, stringsAsFactors = FALSE))
    }
    semicolon <- grepl(";", value, fixed = TRUE)
    comma <- grepl(",", value, fixed = TRUE)
    pipe <- grepl("|", value, fixed = TRUE)
    slash <- grepl("/", value, fixed = TRUE)
    whitespace <- grepl("\\S+\\s+\\S+", value, perl = TRUE)
    any_delimiter <- semicolon | comma | pipe | slash
    data.frame(nonempty_n = length(value), semicolon_n = sum(semicolon), comma_n = sum(comma), pipe_n = sum(pipe), 
        slash_n = sum(slash), whitespace_n = sum(whitespace), multi_delimiter_fraction = mean(any_delimiter), stringsAsFactors = FALSE)
}

numeric_parse_fraction <- function(x) {
    value <- nonempty_values(x)
    if (length(value) == 0L) {
        return(NA_real_)
    }
    parsed <- suppressWarnings(as.numeric(gsub(",", "", value, fixed = TRUE)))
    mean(is.finite(parsed))
}

build_column_audit <- function(source_id, source_class, data_frame) {
    audit_rows <- vector("list", length(names(data_frame)))
    for (column_index in seq_along(names(data_frame))) {
        column_name <- names(data_frame)[column_index]
        column_value <- data_frame[[column_name]]
        values_nonempty <- nonempty_values(column_value)
        role <- infer_field_role(column_name)
        delimiters <- delimiter_audit(column_value)
        audit_rows[[column_index]] <- data.frame(source_id = source_id, source_class = source_class, column_order = column_index, 
            column_name = column_name, inferred_field_role = role, total_rows = nrow(data_frame), nonempty_n = length(values_nonempty), 
            nonempty_fraction = length(values_nonempty)/max(nrow(data_frame), 1L), unique_nonempty_n = length(unique(values_nonempty)), 
            numeric_parse_fraction = numeric_parse_fraction(column_value), semicolon_n = delimiters$semicolon_n, 
            comma_n = delimiters$comma_n, pipe_n = delimiters$pipe_n, slash_n = delimiters$slash_n, whitespace_n = delimiters$whitespace_n, 
            multi_delimiter_fraction = delimiters$multi_delimiter_fraction, example_values = collapse_examples(column_value), 
            stringsAsFactors = FALSE)
    }
    do.call(rbind, audit_rows)
}

build_presence_matrix <- function(audit_table, role_value) {
    subset_table <- audit_table[audit_table$inferred_field_role == role_value, c("source_id", "column_name", "nonempty_fraction"), 
        drop = FALSE]
    if (nrow(subset_table) == 0L) {
        return(data.frame(source_id = character(0), inferred_field_role = character(0), candidate_columns = character(0), 
            maximum_nonempty_fraction = numeric(0), stringsAsFactors = FALSE))
    }
    source_split <- split(subset_table, subset_table$source_id)
    result <- lapply(source_split, function(source_table) {
        data.frame(source_id = unique(source_table$source_id), inferred_field_role = role_value, candidate_columns = paste(source_table$column_name, 
            collapse = ";"), maximum_nonempty_fraction = max(source_table$nonempty_fraction), stringsAsFactors = FALSE)
    })
    do.call(rbind, result)
}

log_event("Starting", run_id)

required_paths <- c(gate_v3 = file.path(dirs$results, "STEP03B3B_v3_GATE_TABLE.csv"), source_lock = file.path(dirs$results, 
    "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK.csv"), source_lock_card = file.path(dirs$results, "STEP03B3B_v3_PUBLISHED_VDR_CISTROME_SOURCE_LOCK_CARD.md"), 
    authorization = file.path(dirs$control, "STEP03B3C_AUTHORIZATION.md"), RNA_lock_card = file.path(dirs$results, 
        "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C1 input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_v3 <- .vdr_read_csv(required_paths[["gate_v3"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock <- .vdr_read_csv(required_paths[["source_lock"]], stringsAsFactors = FALSE, check.names = FALSE)

source_lock_card_text <- readLines(required_paths[["source_lock_card"]], warn = FALSE, encoding = "UTF-8")

authorization_text <- readLines(required_paths[["authorization"]], warn = FALSE, encoding = "UTF-8")

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- gate_v3$status[gate_v3$gate_component == "STEP03B3B v3 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "READY_FOR_STEP03B3C_PEAK_TO_GENE_AND_RNA_INTEGRATION_PREP") {
    stop("STEP03B3B v3 does not authorize STEP03B3C. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (!any(grepl(source_lock_version, source_lock_card_text, fixed = TRUE)) || !all(source_lock$source_lock_version == 
    source_lock_version) || !all(source_lock$source_lock_status == "LOCKED")) {
    stop("PUBLISHED_VDR_CISTROME_SOURCE_v1.0 is not fully locked.")
}

if (!any(grepl("Authorization: YES", authorization_text, fixed = TRUE))) {
    stop("STEP03B3C authorization is incomplete.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
}

primary_sources <- source_lock[source_lock$source_class == "PRIMARY_AUTHOR_PROCESSED_VDR_CISTROME_COORDINATES", 
    , drop = FALSE]

derived_sources <- source_lock[source_lock$source_class == "DERIVED_AUTHOR_CISTROME_RNA_OVERLAP_REFERENCE", , drop = FALSE]

if (nrow(primary_sources) != 3L || nrow(derived_sources) != 2L) {
    stop("Expected three primary cistrome sources and two derived sources.")
}

primary_sources$file_exists <- file.exists(primary_sources$clean_table_path)

primary_sources$observed_md5 <- vapply(primary_sources$clean_table_path, md5_file, character(1))

primary_sources$observed_sha256 <- vapply(primary_sources$clean_table_path, sha256_file, character(1))

primary_sources$identity_pass <- with(primary_sources, file_exists & clean_table_md5 == observed_md5 & clean_table_sha256 == 
    observed_sha256)

if (!all(primary_sources$identity_pass)) {
    stop("One or more locked primary cistrome tables failed hash revalidation.")
}

write_csv_utf8(primary_sources, file.path(dirs$inputs, "STEP03B3C1_PRIMARY_SOURCE_IDENTITY_RECHECK.csv"))

source_tables <- list()

column_audits <- list()

preview_manifests <- list()

for (source_index in seq_len(nrow(primary_sources))) {
    source_id <- primary_sources$source_id[source_index]
    source_path <- primary_sources$clean_table_path[source_index]
    log_event("Auditing annotation fields for", source_id)
    source_table <- .vdr_fread(source_path, sep = "\t", colClasses = "character", data.table = FALSE, showProgress = FALSE, 
        na.strings = c("", "NA", "N/A", "NULL"))
    source_tables[[source_id]] <- source_table
    audit_table <- build_column_audit(source_id = source_id, source_class = primary_sources$source_class[source_index], 
        data_frame = source_table)
    column_audits[[source_id]] <- audit_table
    candidate_columns <- audit_table$column_name[audit_table$inferred_field_role %in% c("CANDIDATE_PROMOTER_FIELD", 
        "CANDIDATE_TSS_DISTANCE_FIELD", "CANDIDATE_NEAREST_GENE_FIELD", "CANDIDATE_GENE_OVERLAP_FIELD", "CANDIDATE_GENE_IDENTIFIER_FIELD", 
        "CANDIDATE_GENOMIC_ANNOTATION_FIELD", "CANDIDATE_PEAK_SCORE_FIELD") & audit_table$nonempty_n > 0L]
    preview_columns <- unique(c(intersect(c("source_id", "source_row", "chromosome", "start0", "end", "width", 
        "coordinate_key"), names(source_table)), candidate_columns))
    preview_table <- source_table[seq_len(min(100L, nrow(source_table))), preview_columns, drop = FALSE]
    preview_path <- file.path(preview_root, paste0(source_id, "_first100_annotation_preview.csv"))
    write_csv_utf8(preview_table, preview_path)
    preview_manifests[[source_id]] <- data.frame(source_id = source_id, preview_path = normalize_slashes(preview_path), 
        preview_rows = nrow(preview_table), preview_columns = ncol(preview_table), source_total_rows = nrow(source_table), 
        source_total_columns = ncol(source_table), stringsAsFactors = FALSE)
}

column_audit <- do.call(rbind, column_audits)

preview_manifest <- do.call(rbind, preview_manifests)

write_csv_utf8(column_audit, file.path(dirs$qc, "STEP03B3C1_ANNOTATION_COLUMN_AUDIT.csv"))

write_csv_utf8(preview_manifest, file.path(preview_root, "STEP03B3C1_PREVIEW_MANIFEST.csv"))

candidate_roles <- c("CANDIDATE_GENE_OVERLAP_FIELD", "CANDIDATE_GENE_IDENTIFIER_FIELD", "CANDIDATE_NEAREST_GENE_FIELD", 
    "CANDIDATE_PROMOTER_FIELD", "CANDIDATE_TSS_DISTANCE_FIELD", "CANDIDATE_GENOMIC_ANNOTATION_FIELD", "CANDIDATE_PEAK_SCORE_FIELD")

role_coverage_rows <- lapply(candidate_roles, function(role_value) {
    build_presence_matrix(column_audit, role_value)
})

role_coverage <- do.call(rbind, role_coverage_rows)

all_source_role_grid <- expand.grid(source_id = primary_sources$source_id, inferred_field_role = candidate_roles, 
    stringsAsFactors = FALSE)

role_coverage_complete <- merge(all_source_role_grid, role_coverage, by = c("source_id", "inferred_field_role"), 
    all.x = TRUE, sort = FALSE)

role_coverage_complete$candidate_columns[is.na(role_coverage_complete$candidate_columns)] <- ""

role_coverage_complete$maximum_nonempty_fraction[is.na(role_coverage_complete$maximum_nonempty_fraction)] <- 0

role_coverage_complete$field_present <- nzchar(role_coverage_complete$candidate_columns) & role_coverage_complete$maximum_nonempty_fraction > 
    0

write_csv_utf8(role_coverage_complete, file.path(dirs$results, "STEP03B3C1_FIELD_ROLE_COVERAGE_MATRIX.csv"))

gene_field_roles <- c("CANDIDATE_GENE_OVERLAP_FIELD", "CANDIDATE_GENE_IDENTIFIER_FIELD", "CANDIDATE_NEAREST_GENE_FIELD")

gene_field_candidates <- column_audit[column_audit$inferred_field_role %in% gene_field_roles & column_audit$nonempty_n > 
    0L, , drop = FALSE]

if (nrow(gene_field_candidates) > 0L) {
    gene_field_candidates$proposed_field_use <- ifelse(gene_field_candidates$inferred_field_role == "CANDIDATE_GENE_OVERLAP_FIELD", 
        "PRIMARY_AUTHOR_GENE_OVERLAP_CANDIDATE", ifelse(gene_field_candidates$inferred_field_role == "CANDIDATE_NEAREST_GENE_FIELD", 
            "SECONDARY_NEAREST_GENE_CANDIDATE", "GENE_IDENTIFIER_REQUIRES_SEMANTIC_REVIEW"))
    gene_field_candidates$manual_decision <- "PENDING_MANUAL_REVIEW"
    gene_field_candidates$authoritative_mapping_status <- "NOT_YET_LOCKED"
}

write_csv_utf8(gene_field_candidates, file.path(dirs$results, "STEP03B3C1_GENE_FIELD_DECISION_TABLE.csv"))

vehicle_or_calcitriol_ids <- c("VDR_VEHICLE_CISTROME", "VDR_CALCITRIOL_CISTROME")

core_gene_field_presence <- vapply(vehicle_or_calcitriol_ids, function(source_id) {
    any(gene_field_candidates$source_id == source_id)
}, logical(1))

shared_gene_field_presence <- any(gene_field_candidates$source_id == "VDR_SHARED_CISTROME")

annotation_field_presence <- vapply(primary_sources$source_id, function(source_id) {
    any(column_audit$source_id == source_id & column_audit$inferred_field_role == "CANDIDATE_GENOMIC_ANNOTATION_FIELD" & 
        column_audit$nonempty_n > 0L)
}, logical(1))

promoter_field_presence <- vapply(primary_sources$source_id, function(source_id) {
    any(column_audit$source_id == source_id & column_audit$inferred_field_role == "CANDIDATE_PROMOTER_FIELD" & 
        column_audit$nonempty_n > 0L)
}, logical(1))

field_summary <- data.frame(source_id = primary_sources$source_id, primary_gene_field_present = vapply(primary_sources$source_id, 
    function(source_id) {
        any(gene_field_candidates$source_id == source_id)
    }, logical(1)), promoter_field_present = promoter_field_presence[primary_sources$source_id], annotation_field_present = annotation_field_presence[primary_sources$source_id], 
    shared_source_coordinate_only_allowed = primary_sources$source_id == "VDR_SHARED_CISTROME", stringsAsFactors = FALSE)

write_csv_utf8(field_summary, file.path(dirs$qc, "STEP03B3C1_SOURCE_FIELD_SUMMARY.csv"))

author_annotation_sufficient_for_manual_review <- all(core_gene_field_presence)

gate_status <- if (author_annotation_sufficient_for_manual_review) {
    "PROCEED_TO_MANUAL_PEAK_ANNOTATION_FIELD_LOCK"
} else {
    "HOLD_AUTHOR_GENE_ANNOTATION_INSUFFICIENT"
}

mapping_policy <- data.frame(policy_component = c("Primary coordinate source", "Primary author gene field", "Promoter evidence", 
    "Nearest-gene evidence", "Shared peak annotation", "One-to-many peak-gene relation", "RNA membership", "Derived media-2 tables", 
    "ATAC context", "Direct-target terminology"), provisional_policy = c("Use only the three locked media-6 canonical coordinate sources.", 
    "Select only after manual semantic review of candidate gene fields.", "Preserve separately from generic gene overlap.", 
    "Preserve as secondary linkage, not equivalent to promoter overlap.", "May be inherited transparently through exact overlap with vehicle/calcitriol sources if no direct gene field exists.", 
    "Retain peak-by-gene long format before any gene-level aggregation.", "RNA_PROGRAM_v1.0 remains immutable.", 
    "Use only as author-derived cross-checks.", "Optional accessible context only; no differential-accessibility claim.", 
    "Use occupancy-supported candidate terminology only after integration QC."), lock_status = c("LOCKED", "PENDING_MANUAL_REVIEW", 
    "PENDING_MANUAL_REVIEW", "PENDING_MANUAL_REVIEW", "PENDING_MANUAL_REVIEW", "LOCKED", "LOCKED", "LOCKED", "LOCKED", 
    "LOCKED"), stringsAsFactors = FALSE)

write_csv_utf8(mapping_policy, file.path(dirs$results, "STEP03B3C1_PROVISIONAL_PEAK_TO_GENE_MAPPING_POLICY.csv"))

gate_table <- data.frame(gate_component = c("STEP03B3B v3 authorization", "Published source lock", "Primary source identities", 
    "RNA membership", "ATAC role", "Vehicle author gene field", "Calcitriol author gene field", "Shared author gene field", 
    "Peak-to-gene long table", "RNA/CUT&RUN integration", "Direct-target assignment", "Scientific figures", "STEP03B3C1 gate"), 
    status = c("PASS", paste0("LOCKED_", source_lock_version), "PASS_3_OF_3", paste0("IMMUTABLE_", RNA_anchor), 
        ATAC_role, ifelse(core_gene_field_presence[["VDR_VEHICLE_CISTROME"]], "PRESENT_REQUIRES_MANUAL_SEMANTIC_REVIEW", 
            "ABSENT"), ifelse(core_gene_field_presence[["VDR_CALCITRIOL_CISTROME"]], "PRESENT_REQUIRES_MANUAL_SEMANTIC_REVIEW", 
            "ABSENT"), ifelse(shared_gene_field_presence, "PRESENT_REQUIRES_MANUAL_SEMANTIC_REVIEW", "ABSENT_COORDINATE_INHERITANCE_OPTION"), 
        "NOT_CREATED", "NOT_STARTED", "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain authorization.", 
        "Do not modify locked source tables.", "Use only hash-revalidated primary sources.", "Do not alter 436 interaction candidates.", 
        "Use only as optional accessible context.", "Review actual values, delimiters and meaning.", "Review actual values, delimiters and meaning.", 
        "Decide whether exact-overlap inheritance is required.", "Create only after manual field lock.", "Do not integrate until field semantics are locked.", 
        "Use candidate terminology only after later QC.", "No Scientific Figure Skill 3 workflow is triggered.", 
        "Upload field audit, decision table, coverage matrix and previews."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C1_GATE_TABLE.csv"))

write_text(c("# STEP03B3C1 annotation-field audit card", "", paste0("- Run ID: ", run_id), paste0("- Source lock: ", 
    source_lock_version), "- Primary coordinate sources audited: 3", paste0("- Candidate gene-field rows: ", nrow(gene_field_candidates)), 
    paste0("- Vehicle gene field present: ", core_gene_field_presence[["VDR_VEHICLE_CISTROME"]]), paste0("- Calcitriol gene field present: ", 
        core_gene_field_presence[["VDR_CALCITRIOL_CISTROME"]]), paste0("- Shared gene field present: ", shared_gene_field_presence), 
    "- Final peak-to-gene mapping: NOT CREATED", "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C1_ANNOTATION_FIELD_AUDIT_CARD.md"))

write_text(c("# STEP03B3C1 manual review guide", "", "Review the following outputs in order:", "", "1. `STEP03B3C1_ANNOTATION_COLUMN_AUDIT.csv`", 
    "   - confirm which raw columns contain gene symbols, promoter labels,", "     genomic annotations, distances and scores;", 
    "   - inspect example values and delimiter patterns.", "", "2. `STEP03B3C1_GENE_FIELD_DECISION_TABLE.csv`", 
    "   - choose the primary gene-overlap field;", "   - distinguish gene-overlap from nearest-gene fields;", "   - do not merge multiple genes into one value.", 
    "", "3. `STEP03B3C1_FIELD_ROLE_COVERAGE_MATRIX.csv`", "   - compare field availability across vehicle, calcitriol and shared", 
    "     cistrome sources.", "", "4. `STEP03B3C1_peak_annotation_previews/`", "   - inspect the first 100 rows for each source using the original", 
    "     values preserved in the locked clean source table.", "", "The next substep will create a peak-by-gene long table only after", 
    "these field semantics are manually locked."), file.path(dirs$text, "STEP03B3C1_MANUAL_REVIEW_GUIDE.md"))

write_text(c("# STEP03B3C1 Methods draft", "", "The three locked author-processed VDR cistrome coordinate tables were revalidated by MD5 and SHA-256 and imported with all columns represented as character values. Raw fields were audited for completeness, unique-value burden, numeric parseability, delimiter patterns and candidate scientific roles, including gene overlap, nearest gene, promoter overlap, TSS distance, genomic annotation and peak-score fields.", 
    "", "This substep was restricted to annotation-schema assessment. No peak-to-gene relationship was yet treated as authoritative, no peak-by-gene long table was created and no integration with the frozen RNA interaction program was performed."), 
    file.path(dirs$text, "STEP03B3C1_Methods_draft.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, source_lock_version = source_lock_version, primary_sources_audited = 3, 
    derived_sources_used = 0, RNA_anchor = RNA_anchor, ATAC_role = ATAC_role, final_peak_to_gene_mapping_created = FALSE, 
    RNA_CUTRUN_integration = FALSE, direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C1_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C1.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("37_step03B3C1_peak_annotation_field_schema_audit_path_fixed.R"), 
        overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "read_peak_gene_annotations.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$inputs, "STEP03B3C1_PRIMARY_SOURCE_IDENTITY_RECHECK.csv"), file.path(dirs$results, 
    "STEP03B3C1_FIELD_ROLE_COVERAGE_MATRIX.csv"), file.path(dirs$results, "STEP03B3C1_GENE_FIELD_DECISION_TABLE.csv"), 
    file.path(dirs$results, "STEP03B3C1_PROVISIONAL_PEAK_TO_GENE_MAPPING_POLICY.csv"), file.path(dirs$results, 
        "STEP03B3C1_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C1_ANNOTATION_FIELD_AUDIT_CARD.md"), file.path(dirs$qc, 
        "STEP03B3C1_ANNOTATION_COLUMN_AUDIT.csv"), file.path(dirs$qc, "STEP03B3C1_SOURCE_FIELD_SUMMARY.csv"), file.path(preview_root, 
        "STEP03B3C1_PREVIEW_MANIFEST.csv"), file.path(dirs$text, "STEP03B3C1_MANUAL_REVIEW_GUIDE.md"), file.path(dirs$text, 
        "STEP03B3C1_Methods_draft.md"), file.path(dirs$parameters, "STEP03B3C1_parameters.json"), file.path(dirs$environment, 
        "sessionInfo_POST_RUN_STEP03B3C1.txt"), file.path(dirs$code, "read_peak_gene_annotations.R"), log_path)

preview_files <- list.files(preview_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_artifacts <- unique(c(archive_artifacts, preview_files))

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
    else if (grepl("Methods|REVIEW_GUIDE", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("AUDIT|SUMMARY|PREVIEW|DECISION|POLICY|CARD|MATRIX", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C1_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C1 archive status", "", paste0("- Run ID: ", run_id), paste0("- Source lock: ", source_lock_version), 
    "- Primary sources audited: 3", paste0("- Candidate gene-field rows: ", nrow(gene_field_candidates)), "- Final peak-to-gene mapping: NOT_CREATED", 
    "- RNA/CUT&RUN integration: NOT_STARTED", "- Direct-target assignment: NOT_PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C1", paste0("- Run ID: ", run_id), paste0("- Published source lock: ", 
    source_lock_version, " IMMUTABLE"), paste0("- RNA anchor: ", RNA_anchor, " IMMUTABLE"), paste0("- ATAC role: ", 
    ATAC_role), "- Peak annotation fields: AUDITED_PENDING_MANUAL_LOCK", "- Peak-to-gene long table: NOT_CREATED", 
    "- RNA/CUT&RUN integration: NOT_STARTED", "- Direct-target assignment: NOT_PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C1 peak-annotation field audit completed.\n")

cat("Source lock: ", source_lock_version, " IMMUTABLE\n")

cat("Primary sources audited: 3/3\n")

cat("Candidate gene-field rows: ", nrow(gene_field_candidates), "\n")

cat("Vehicle gene field present: ", core_gene_field_presence[["VDR_VEHICLE_CISTROME"]], "\n")

cat("Calcitriol gene field present: ", core_gene_field_presence[["VDR_CALCITRIOL_CISTROME"]], "\n")

cat("Shared gene field present: ", shared_gene_field_presence, "\n")

cat("Peak-to-gene long table: NOT CREATED\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C1_GATE_TABLE.csv\n")

cat("  06_results/STEP03B3C1_GENE_FIELD_DECISION_TABLE.csv\n")

cat("  06_results/STEP03B3C1_FIELD_ROLE_COVERAGE_MATRIX.csv\n")

cat("  08_qc/STEP03B3C1_ANNOTATION_COLUMN_AUDIT.csv\n")

cat("  08_qc/STEP03B3C1_SOURCE_FIELD_SUMMARY.csv\n")

cat("  06_results/STEP03B3C1_peak_annotation_previews/\n")

cat("============================================================\n")

