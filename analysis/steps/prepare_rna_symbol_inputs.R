#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A14b: Locate accepted RNA symbol input tables

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

run_id <- paste0("STEP03B3C3A_", format(Sys.time(), "%Y%m%d_%H%M%S"))

RNA_anchor <- "RNA_PROGRAM_v1.0"

expected_RNA_gene_n <- 436L

peak_gene_map_lock_version <- "PEAK_GENE_MAP_v1.1"

expected_cistrome_gene_n <- 10042L

ATAC_role <- "RESTRICTED_DESCRIPTIVE_SUPPORT_ONLY"

maximum_candidate_file_size_bytes <- 250 * 1024^2

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

audit_root <- file.path(dirs$results, "STEP03B3C3A_RNA_source_and_gene_ID_audit")

dir.create(audit_root, recursive = TRUE, showWarnings = FALSE)

archive_root <- file.path(dirs$archive_parent, paste0("S3B3C3A_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

gene_column_name_score <- function(column_name) {
    name_lower <- tolower(column_name)
    score <- 0
    if (name_lower %in% c("gene_symbol", "genesymbol", "symbol", "hgnc_symbol")) {
        score <- score + 10
    }
    if (grepl("gene.*symbol|symbol.*gene", name_lower, perl = TRUE)) {
        score <- score + 7
    }
    if (grepl("(^|_)gene($|_)", name_lower, perl = TRUE)) {
        score <- score + 4
    }
    if (grepl("ensembl|entrez|probe|transcript", name_lower, perl = TRUE)) {
        score <- score - 4
    }
    score
}

filename_score <- function(path) {
    name_lower <- tolower(basename(path))
    score <- 0
    if (grepl("rna_program", name_lower, fixed = TRUE)) {
        score <- score + 10
    }
    if (grepl("interaction", name_lower, fixed = TRUE)) {
        score <- score + 9
    }
    if (grepl("contrast", name_lower, fixed = TRUE)) {
        score <- score + 5
    }
    if (grepl("candidate", name_lower, fixed = TRUE)) {
        score <- score + 3
    }
    if (grepl("attenuat", name_lower, fixed = TRUE)) {
        score <- score + 3
    }
    if (grepl("gate|manifest|lock|card|qc|summary|figure|plot|mapping|peak|archive", name_lower, perl = TRUE)) {
        score <- score - 8
    }
    score
}

interaction_field_score <- function(column_names) {
    name_lower <- tolower(column_names)
    patterns <- c("interaction", "logfc", "log2fc", "adj.p", "fdr", "p.value", "pvalue", "class", "attenuat", "contrast", 
        "direction")
    sum(vapply(patterns, function(pattern_value) {
        any(grepl(pattern_value, name_lower, fixed = TRUE))
    }, logical(1)))
}

read_candidate_table <- function(path) {
    extension <- tolower(tools::file_ext(path))
    if (extension %in% c("csv", "tsv", "txt")) {
        return(.vdr_fread(path, data.table = FALSE, showProgress = FALSE))
    }
    if (extension == "rds") {
        object_value <- readRDS(path)
        if (is.data.frame(object_value) || data.table::is.data.table(object_value)) {
            return(as.data.frame(object_value, stringsAsFactors = FALSE))
        }
        return(NULL)
    }
    NULL
}

audit_file_candidate <- function(path) {
    file_size <- as.numeric(file.info(path)$size)
    base_result <- data.frame(file_path = normalize_slashes(path), file_name = basename(path), extension = tolower(tools::file_ext(path)), 
        size_bytes = file_size, filename_score = filename_score(path), read_status = "NOT_READ", n_rows = NA_integer_, 
        n_columns = NA_integer_, interaction_field_score = NA_integer_, best_gene_column = NA_character_, best_gene_column_score = NA_real_, 
        best_gene_nonempty_n = NA_integer_, best_gene_unique_n = NA_integer_, best_gene_valid_shape_fraction = NA_real_, 
        exact_436_unique_genes = FALSE, exact_436_rows = FALSE, candidate_score = NA_real_, error_message = NA_character_, 
        stringsAsFactors = FALSE)
    if (!is.finite(file_size) || file_size > maximum_candidate_file_size_bytes) {
        base_result$read_status <- "SKIPPED_FILE_TOO_LARGE"
        return(list(summary = base_result, gene_audit = data.frame(), data = NULL))
    }
    table_value <- tryCatch(read_candidate_table(path), error = function(error_object) {
        structure(list(message = conditionMessage(error_object)), class = "candidate_read_error")
    })
    if (inherits(table_value, "candidate_read_error")) {
        base_result$read_status <- "READ_ERROR"
        base_result$error_message <- table_value$message
        return(list(summary = base_result, gene_audit = data.frame(), data = NULL))
    }
    if (is.null(table_value) || !is.data.frame(table_value)) {
        base_result$read_status <- "UNSUPPORTED_OBJECT"
        return(list(summary = base_result, gene_audit = data.frame(), data = NULL))
    }
    base_result$read_status <- "READ_OK"
    base_result$n_rows <- nrow(table_value)
    base_result$n_columns <- ncol(table_value)
    base_result$interaction_field_score <- interaction_field_score(names(table_value))
    candidate_gene_columns <- names(table_value)[vapply(names(table_value), gene_column_name_score, numeric(1)) > 
        0]
    if (length(candidate_gene_columns) == 0L) {
        base_result$candidate_score <- base_result$filename_score + base_result$interaction_field_score
        return(list(summary = base_result, gene_audit = data.frame(), data = table_value))
    }
    gene_audit_rows <- lapply(candidate_gene_columns, function(column_name) {
        symbols <- clean_symbol(table_value[[column_name]])
        nonempty_symbols <- symbols[!is.na(symbols)]
        unique_symbols <- unique(nonempty_symbols)
        valid_fraction <- if (length(nonempty_symbols) > 0L) {
            mean(valid_symbol_shape(nonempty_symbols))
        }
        else {
            NA_real_
        }
        data.frame(file_path = normalize_slashes(path), column_name = column_name, column_name_score = gene_column_name_score(column_name), 
            nonempty_n = length(nonempty_symbols), unique_n = length(unique_symbols), duplicate_nonempty_n = length(nonempty_symbols) - 
                length(unique_symbols), valid_shape_fraction = valid_fraction, exact_436_unique = length(unique_symbols) == 
                expected_RNA_gene_n, exact_436_rows = nrow(table_value) == expected_RNA_gene_n, column_candidate_score = gene_column_name_score(column_name) + 
                ifelse(length(unique_symbols) == expected_RNA_gene_n, 20, 0) + ifelse(nrow(table_value) == expected_RNA_gene_n, 
                5, 0) + ifelse(is.finite(valid_fraction) && valid_fraction >= 0.98999999999999999, 3, 0), stringsAsFactors = FALSE)
    })
    gene_audit <- do.call(rbind, gene_audit_rows)
    gene_audit <- gene_audit[order(-gene_audit$column_candidate_score, gene_audit$column_name), , drop = FALSE]
    best_row <- gene_audit[1L, , drop = FALSE]
    base_result$best_gene_column <- best_row$column_name
    base_result$best_gene_column_score <- best_row$column_candidate_score
    base_result$best_gene_nonempty_n <- best_row$nonempty_n
    base_result$best_gene_unique_n <- best_row$unique_n
    base_result$best_gene_valid_shape_fraction <- best_row$valid_shape_fraction
    base_result$exact_436_unique_genes <- best_row$exact_436_unique
    base_result$exact_436_rows <- best_row$exact_436_rows
    base_result$candidate_score <- base_result$filename_score + base_result$interaction_field_score + base_result$best_gene_column_score
    list(summary = base_result, gene_audit = gene_audit, data = table_value)
}

log_event("Starting", run_id)

locked_map_root <- file.path(dirs$results, "STEP03B3C2_locked_peak_gene_map", peak_gene_map_lock_version)

required_paths <- c(gate_c2 = file.path(dirs$results, "STEP03B3C2_v3_GATE_TABLE.csv"), map_manifest = file.path(locked_map_root, 
    "PEAK_GENE_MAP_v1.1_LOCKED_MANIFEST.csv"), authorization = file.path(dirs$control, "STEP03B3C3_AUTHORIZATION.md"), 
    RNA_lock_card = file.path(dirs$results, "STEP03B1_v10_RNA_PROGRAM_LOCK_CARD.md"))

missing_required_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_required_paths) > 0L) {
    stop("Required STEP03B3C3A input is missing:\n", paste(missing_required_paths, collapse = "\n"))
}

gate_c2 <- .vdr_read_csv(required_paths[["gate_c2"]], stringsAsFactors = FALSE, check.names = FALSE)

map_manifest <- .vdr_read_csv(required_paths[["map_manifest"]], stringsAsFactors = FALSE, check.names = FALSE)

authorization_text <- readLines(required_paths[["authorization"]], warn = FALSE, encoding = "UTF-8")

RNA_lock_text <- readLines(required_paths[["RNA_lock_card"]], warn = FALSE, encoding = "UTF-8")

accepted_gate <- gate_c2$status[gate_c2$gate_component == "STEP03B3C2 v3 gate"]

if (length(accepted_gate) != 1L || accepted_gate != "READY_FOR_STEP03B3C3_RNA_GENE_ID_HARMONIZATION_AND_OCCUPANCY_INTEGRATION") {
    stop("STEP03B3C2 v3 does not authorize STEP03B3C3. Gate: ", paste(accepted_gate, collapse = ", "))
}

if (!any(grepl("Authorization: YES", authorization_text, fixed = TRUE))) {
    stop("STEP03B3C3 authorization is incomplete.")
}

if (!any(grepl(RNA_anchor, RNA_lock_text, fixed = TRUE)) || !any(grepl("RNA statistical membership: FROZEN", RNA_lock_text, 
    fixed = TRUE))) {
    stop("RNA_PROGRAM_v1.0 is not confirmed as frozen.")
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

search_roots <- c(dirs$inputs, dirs$results)

all_files <- unique(unlist(lapply(search_roots, function(root_path) {
    list.files(root_path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
}), use.names = FALSE))

allowed_extensions <- c("csv", "tsv", "txt", "rds")

candidate_name_pattern <- paste(c("STEP03B1", "RNA_PROGRAM", "interaction", "contrast", "attenuat", "candidate"), 
    collapse = "|")

candidate_files <- all_files[tolower(tools::file_ext(all_files)) %in% allowed_extensions & grepl(candidate_name_pattern, 
    basename(all_files), ignore.case = TRUE, perl = TRUE) & !grepl("STEP03B3C|peak_gene|cistrome|manifest|gate|lock_card|sessionInfo|parameters", 
    all_files, ignore.case = TRUE, perl = TRUE)]

candidate_files <- sort(unique(candidate_files))

if (length(candidate_files) == 0L) {
    stop("No plausible RNA interaction table files were found in 01_inputs or 06_results.")
}

log_event("Candidate files discovered:", length(candidate_files))

candidate_audits <- lapply(candidate_files, audit_file_candidate)

candidate_inventory <- do.call(rbind, lapply(candidate_audits, function(item) item[["summary"]]))

gene_column_audit <- data.table::rbindlist(lapply(candidate_audits, function(item) item[["gene_audit"]]), fill = TRUE)

gene_column_audit <- as.data.frame(gene_column_audit)

candidate_inventory <- candidate_inventory[order(-candidate_inventory$exact_436_unique_genes, -candidate_inventory$candidate_score, 
    candidate_inventory$file_name), , drop = FALSE]

write_csv_utf8(candidate_inventory, file.path(dirs$qc, "STEP03B3C3A_RNA_SOURCE_CANDIDATE_INVENTORY.csv"))

write_csv_utf8(gene_column_audit, file.path(dirs$qc, "STEP03B3C3A_RNA_GENE_COLUMN_AUDIT.csv"))

high_confidence_candidates <- candidate_inventory[candidate_inventory$read_status == "READ_OK" & candidate_inventory$exact_436_unique_genes & 
    candidate_inventory$interaction_field_score >= 2L & candidate_inventory$best_gene_valid_shape_fraction >= 0.98999999999999999, 
    , drop = FALSE]

source_selection_status <- if (nrow(high_confidence_candidates) == 1L) {
    "UNIQUE_HIGH_CONFIDENCE_CANDIDATE"
} else if (nrow(high_confidence_candidates) > 1L) {
    "MULTIPLE_HIGH_CONFIDENCE_CANDIDATES"
} else {
    "NO_HIGH_CONFIDENCE_CANDIDATE"
}

RNA_source_candidate <- data.frame()

RNA_gene_audit <- data.frame()

match_summary <- data.frame()

case_collision_audit <- data.frame()

selected_table <- NULL

selected_gene_column <- NA_character_

selected_path <- NA_character_

if (source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE") {
    selected_path <- high_confidence_candidates$file_path[[1]]
    selected_gene_column <- high_confidence_candidates$best_gene_column[[1]]
    selected_table <- read_candidate_table(selected_path)
    RNA_symbols <- clean_symbol(selected_table[[selected_gene_column]])
    if (length(RNA_symbols) != nrow(selected_table) || any(is.na(RNA_symbols))) {
        stop("Selected RNA candidate contains missing gene symbols.")
    }
    if (length(unique(RNA_symbols)) != expected_RNA_gene_n) {
        stop("Selected RNA candidate does not contain exactly 436 unique symbols.")
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
    cistrome_upper <- toupper(cistrome_symbols)
    RNA_upper <- toupper(RNA_symbols)
    cistrome_upper_count <- table(cistrome_upper)
    RNA_upper_count <- table(RNA_upper)
    exact_match <- RNA_symbols %in% cistrome_symbols
    case_insensitive_match <- RNA_upper %in% cistrome_upper
    case_only_match <- !exact_match & case_insensitive_match
    RNA_case_collision <- as.integer(RNA_upper_count[RNA_upper]) > 1L
    cistrome_case_collision <- as.integer(cistrome_upper_count[RNA_upper]) > 1L
    cistrome_case_candidates <- vapply(RNA_upper, function(symbol_upper) {
        paste(cistrome_symbols[cistrome_upper == symbol_upper], collapse = ";")
    }, character(1))
    RNA_gene_audit <- data.frame(RNA_row = seq_len(nrow(selected_table)), RNA_gene_symbol = RNA_symbols, valid_symbol_shape = valid_symbol_shape(RNA_symbols), 
        exact_cistrome_match = exact_match, case_only_cistrome_match = case_only_match, no_cistrome_match = !case_insensitive_match, 
        RNA_case_collision = RNA_case_collision, cistrome_case_collision = cistrome_case_collision, cistrome_case_candidates = cistrome_case_candidates, 
        automatic_alias_mapping = FALSE, stringsAsFactors = FALSE)
    match_summary <- data.frame(metric = c("RNA rows", "Unique RNA symbols", "Valid RNA symbol shape", "Exact cistrome matches", 
        "Case-only cistrome matches", "No cistrome match", "RNA case-collision rows", "Cistrome case-collision rows", 
        "Cistrome symbol universe", "Automatic alias mappings"), value = c(nrow(selected_table), length(unique(RNA_symbols)), 
        sum(RNA_gene_audit$valid_symbol_shape), sum(RNA_gene_audit$exact_cistrome_match), sum(RNA_gene_audit$case_only_cistrome_match), 
        sum(RNA_gene_audit$no_cistrome_match), sum(RNA_gene_audit$RNA_case_collision), sum(RNA_gene_audit$cistrome_case_collision), 
        length(cistrome_symbols), 0L), stringsAsFactors = FALSE)
    collision_rows <- RNA_gene_audit[RNA_gene_audit$RNA_case_collision | RNA_gene_audit$cistrome_case_collision | 
        RNA_gene_audit$case_only_cistrome_match, , drop = FALSE]
    case_collision_audit <- collision_rows
    RNA_source_candidate <- data.frame(source_candidate_id = "RNA_INTERACTION_436_SOURCE_CANDIDATE", RNA_anchor = RNA_anchor, 
        file_path = normalize_slashes(selected_path), file_name = basename(selected_path), size_bytes = as.numeric(file.info(selected_path)$size), 
        md5 = md5_file(selected_path), sha256 = sha256_file(selected_path), n_rows = nrow(selected_table), n_columns = ncol(selected_table), 
        gene_column = selected_gene_column, unique_gene_n = length(unique(RNA_symbols)), interaction_field_score = high_confidence_candidates$interaction_field_score[[1]], 
        source_status = "CANDIDATE_PENDING_MANUAL_LOCK", stringsAsFactors = FALSE)
    write_csv_utf8(RNA_gene_audit, file.path(audit_root, "STEP03B3C3A_RNA_GENE_ID_AUDIT.csv"))
    write_csv_utf8(match_summary, file.path(audit_root, "STEP03B3C3A_RNA_CISTROME_MATCH_SUMMARY.csv"))
    write_csv_utf8(case_collision_audit, file.path(audit_root, "STEP03B3C3A_CASE_COLLISION_AUDIT.csv"))
    write_csv_utf8(RNA_source_candidate, file.path(dirs$results, "STEP03B3C3A_RNA_INTERACTION_SOURCE_CANDIDATE.csv"))
}

gate_status <- if (source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE") {
    "PROCEED_TO_MANUAL_RNA_SOURCE_AND_GENE_ID_LOCK"
} else if (source_selection_status == "MULTIPLE_HIGH_CONFIDENCE_CANDIDATES") {
    "HOLD_MULTIPLE_RNA_SOURCE_CANDIDATES"
} else {
    "HOLD_NO_AUTHORITATIVE_RNA_SOURCE_CANDIDATE"
}

gate_table <- data.frame(gate_component = c("STEP03B3C2 authorization", "Peak-gene map lock", "RNA program lock", 
    "RNA candidate files inventoried", "High-confidence RNA source candidates", "RNA source-selection status", 
    "RNA unique-gene requirement", "RNA/cistrome exact-match audit", "Automatic alias mapping", "RNA membership", 
    "Peak-gene membership", "ATAC role", "RNA/CUT&RUN integration", "Definitive direct-target assignment", "Scientific figures", 
    "STEP03B3C3A gate"), status = c("PASS", paste0("IMMUTABLE_", peak_gene_map_lock_version), paste0("IMMUTABLE_", 
    RNA_anchor), as.character(length(candidate_files)), as.character(nrow(high_confidence_candidates)), source_selection_status, 
    ifelse(source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE", "PASS_436_UNIQUE", "NOT_LOCKED"), ifelse(source_selection_status == 
        "UNIQUE_HIGH_CONFIDENCE_CANDIDATE", "COMPLETED", "NOT_STARTED"), "ZERO_PROHIBITED", "UNCHANGED_436", "UNCHANGED_23499", 
    ATAC_role, "NOT_STARTED", "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain authorization.", 
    "Do not alter 23,499 peak-gene relations.", "Do not alter the frozen 436 RNA candidates.", "Review the full candidate inventory.", 
    "Require one uniquely supported source before integration.", "Resolve ambiguity manually if present.", "Require exactly 436 nonmissing unique symbols.", 
    "Review exact, case-only and unmatched symbols.", "Do not silently update aliases or capitalization.", "Membership changes are prohibited.", 
    "Membership changes are prohibited.", "Use only as optional accessible context.", "Do not integrate until RNA source and gene IDs are locked.", 
    "Use occupancy-supported candidate terminology only after later QC.", "No Scientific Figure Skill 3 workflow is triggered.", 
    "Upload candidate inventory, source candidate and gene-ID audit."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03B3C3A_GATE_TABLE.csv"))

write_text(c("# STEP03B3C3A RNA source and gene-ID audit card", "", paste0("- Run ID: ", run_id), paste0("- RNA anchor: ", 
    RNA_anchor, " IMMUTABLE"), paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE"), paste0("- Candidate files inventoried: ", 
    length(candidate_files)), paste0("- High-confidence 436-gene candidates: ", nrow(high_confidence_candidates)), 
    paste0("- Source selection: ", source_selection_status), if (source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE") {
        paste0("- Selected candidate: ", basename(selected_path))
    } else {
        "- Selected candidate: NONE"
    }, if (source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE") {
        paste0("- Gene column: ", selected_gene_column)
    } else {
        "- Gene column: NOT_LOCKED"
    }, "- Automatic alias mapping: ZERO", "- RNA/CUT&RUN integration: NOT STARTED", "- Definitive direct-target assignment: NOT PERFORMED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, "STEP03B3C3A_RNA_SOURCE_GENE_ID_AUDIT_CARD.md"))

write_text(c("# STEP03B3C3A Methods draft", "", "Candidate RNA interaction tables were identified from the locked module input and result directories using filename, row-count, gene-column and interaction-field evidence. Candidate tables were required to contain exactly 436 unique nonmissing gene symbols and at least two interaction-model fields before being considered high confidence.", 
    "", "For a uniquely selected candidate, RNA symbols were compared with the locked 10,042-symbol cistrome universe using exact case-sensitive matching followed by a separate case-insensitive audit. No alias, capitalization or historical-symbol mapping was performed automatically. Unmatched and collision-prone symbols were retained explicitly for manual review.", 
    "", "No RNA/cistrome evidence integration was performed in this substep."), file.path(dirs$text, "STEP03B3C3A_Methods_draft.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, RNA_anchor = RNA_anchor, expected_RNA_gene_n = expected_RNA_gene_n, 
    peak_gene_map_lock_version = peak_gene_map_lock_version, expected_cistrome_gene_n = expected_cistrome_gene_n, 
    candidate_files_inventoried = length(candidate_files), high_confidence_candidates = nrow(high_confidence_candidates), 
    source_selection_status = source_selection_status, selected_path = selected_path, selected_gene_column = selected_gene_column, 
    automatic_alias_mapping = FALSE, RNA_membership_changed = FALSE, peak_gene_membership_changed = FALSE, RNA_CUTRUN_integration = FALSE, 
    definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03B3C3A_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("41_step03B3C3A_RNA_source_discovery_gene_ID_audit.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "prepare_rna_symbol_inputs.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP03B3C3A_GATE_TABLE.csv"), file.path(dirs$results, "STEP03B3C3A_RNA_SOURCE_GENE_ID_AUDIT_CARD.md"), 
    file.path(dirs$results, "STEP03B3C3A_RNA_INTERACTION_SOURCE_CANDIDATE.csv"), file.path(dirs$qc, "STEP03B3C3A_RNA_SOURCE_CANDIDATE_INVENTORY.csv"), 
    file.path(dirs$qc, "STEP03B3C3A_RNA_GENE_COLUMN_AUDIT.csv"), file.path(audit_root, "STEP03B3C3A_RNA_GENE_ID_AUDIT.csv"), 
    file.path(audit_root, "STEP03B3C3A_RNA_CISTROME_MATCH_SUMMARY.csv"), file.path(audit_root, "STEP03B3C3A_CASE_COLLISION_AUDIT.csv"), 
    file.path(dirs$text, "STEP03B3C3A_Methods_draft.md"), file.path(dirs$parameters, "STEP03B3C3A_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03B3C3A.txt"), file.path(dirs$code, "prepare_rna_symbol_inputs.R"), 
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
    else if (grepl("AUDIT|SUMMARY|INVENTORY|CARD|CANDIDATE", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP03B3C3A_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP03B3C3A archive status", "", paste0("- Run ID: ", run_id), paste0("- Candidate files: ", length(candidate_files)), 
    paste0("- High-confidence candidates: ", nrow(high_confidence_candidates)), paste0("- Source selection: ", 
        source_selection_status), "- RNA membership changed: NO", "- Peak-gene membership changed: NO", "- Automatic alias mapping: NO", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status), paste0("- Maximum archive path length: ", 
        max(archive_manifest$path_length))), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03B3C3A", paste0("- Run ID: ", run_id), paste0("- RNA anchor: ", 
    RNA_anchor, " IMMUTABLE"), paste0("- Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE"), paste0("- RNA source selection: ", 
    source_selection_status), "- Gene-ID audit: COMPLETED_IF_UNIQUE_SOURCE", "- Automatic alias mapping: PROHIBITED", 
    "- RNA/CUT&RUN integration: NOT STARTED", "- Direct-target assignment: NOT PERFORMED", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status), "- Full STEP03 module: NOT_LOCKED"), file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP03B3C3A RNA source and gene-ID audit completed.\n")

cat("RNA anchor: ", RNA_anchor, " IMMUTABLE\n")

cat("Peak-gene map: ", peak_gene_map_lock_version, " IMMUTABLE\n")

cat("Candidate files inventoried: ", length(candidate_files), "\n")

cat("High-confidence 436-gene candidates: ", nrow(high_confidence_candidates), "\n")

cat("Source selection: ", source_selection_status, "\n")

if (source_selection_status == "UNIQUE_HIGH_CONFIDENCE_CANDIDATE") {
    cat("Selected RNA file: ", basename(selected_path), "\n")
    cat("Gene column: ", selected_gene_column, "\n")
    cat("Exact cistrome matches: ", match_summary$value[match_summary$metric == "Exact cistrome matches"], "\n")
    cat("Case-only matches: ", match_summary$value[match_summary$metric == "Case-only cistrome matches"], "\n")
    cat("No cistrome match: ", match_summary$value[match_summary$metric == "No cistrome match"], "\n")
}

cat("Automatic alias mapping: ZERO\n")

cat("RNA/CUT&RUN integration: NOT STARTED\n")

cat("Definitive direct-target assignment: NOT PERFORMED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP03B3C3A_GATE_TABLE.csv\n")

cat("  08_qc/STEP03B3C3A_RNA_SOURCE_CANDIDATE_INVENTORY.csv\n")

cat("  08_qc/STEP03B3C3A_RNA_GENE_COLUMN_AUDIT.csv\n")

cat("  06_results/STEP03B3C3A_RNA_INTERACTION_SOURCE_CANDIDATE.csv\n")

cat("  06_results/STEP03B3C3A_RNA_source_and_gene_ID_audit/\n")

cat("============================================================\n")

