#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B00: Prepare patient source inventory

options(stringsAsFactors = FALSE)

options(warn = 1)

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04A_", format(Sys.time(), "%Y%m%d_%H%M%S"))

upstream_lock_version <- "RNA_OCCUPANCY_EVIDENCE_v1.0"

routing_lock_version <- "DOWNSTREAM_CANDIDATE_ROUTING_v1.0"

required_RNA_total <- 436L

required_occupancy_supported <- 284L

required_RNA_only <- 152L

required_packages <- c("data.table", "digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(step04_root, "00_module_control"), inputs = file.path(step04_root, "01_inputs"), 
    code = file.path(step04_root, "02_code"), parameters = file.path(step04_root, "03_parameters"), environment = file.path(step04_root, 
        "04_environment"), execution = file.path(step04_root, "05_execution"), results = file.path(step04_root, 
        "06_results"), figures = file.path(step04_root, "07_figures_tables"), qc = file.path(step04_root, "08_qc"), 
    text = file.path(step04_root, "09_manuscript_text"), archive_parent = file.path(step04_root, "10_archive_work"))

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

archive_root <- file.path(dirs$archive_parent, paste0("S4A_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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

classify_file_role <- function(path, n_rows = NA_integer_, n_columns = NA_integer_) {
    path_lower <- tolower(path)
    extension <- tolower(tools::file_ext(path))
    if (grepl("seurat|single.?cell|scrna|sce|h5ad|h5seurat", path_lower, perl = TRUE)) {
        return("SINGLE_CELL_OBJECT_OR_MATRIX_CANDIDATE")
    }
    if (grepl("metadata|phenodata|pheno|clinical|sample.?info|coldata|group", path_lower, perl = TRUE)) {
        return("SAMPLE_METADATA_CANDIDATE")
    }
    if (grepl("expression|expr|matrix|counts|normalized|fpkm|tpm|rpkm|series_matrix", path_lower, perl = TRUE)) {
        return("EXPRESSION_MATRIX_CANDIDATE")
    }
    if (extension %in% c("rds", "rda", "rdata", "h5ad", "h5", "h5seurat", "mtx")) {
        return("BINARY_OR_SPARSE_DATA_OBJECT_CANDIDATE")
    }
    if (is.finite(n_rows) && is.finite(n_columns) && n_rows >= 1000L && n_columns >= 4L) {
        return("TABULAR_EXPRESSION_MATRIX_CANDIDATE")
    }
    if (is.finite(n_rows) && is.finite(n_columns) && n_rows <= 1000L && n_columns >= 3L) {
        return("TABULAR_METADATA_OR_RESULT_CANDIDATE")
    }
    "UNCLASSIFIED_RELEVANT_FILE"
}

lightweight_table_audit <- function(path) {
    extension <- tolower(tools::file_ext(path))
    result <- data.frame(read_status = "NOT_ATTEMPTED", object_class = NA_character_, n_rows = NA_integer_, n_columns = NA_integer_, 
        first_columns = NA_character_, error_message = NA_character_, stringsAsFactors = FALSE)
    file_size <- as.numeric(file.info(path)$size)
    if (!is.finite(file_size) || file_size > 150 * 1024^2) {
        result$read_status <- "SKIPPED_FILE_OVER_150MB"
        return(result)
    }
    if (extension %in% c("csv", "tsv", "txt")) {
        table_value <- tryCatch(.vdr_fread(path, data.table = FALSE, nrows = 200000L, showProgress = FALSE), error = function(error_object) {
            structure(list(message = conditionMessage(error_object)), class = "read_error")
        })
        if (inherits(table_value, "read_error")) {
            result$read_status <- "READ_ERROR"
            result$error_message <- table_value$message
            return(result)
        }
        result$read_status <- if (nrow(table_value) >= 200000L) {
            "READ_OK_CAPPED_AT_200000_ROWS"
        }
        else {
            "READ_OK"
        }
        result$object_class <- paste(class(table_value), collapse = ";")
        result$n_rows <- nrow(table_value)
        result$n_columns <- ncol(table_value)
        result$first_columns <- paste(head(names(table_value), 12L), collapse = ";")
        return(result)
    }
    if (extension == "rds") {
        object_value <- tryCatch(readRDS(path), error = function(error_object) {
            structure(list(message = conditionMessage(error_object)), class = "read_error")
        })
        if (inherits(object_value, "read_error")) {
            result$read_status <- "READ_ERROR"
            result$error_message <- object_value$message
            return(result)
        }
        result$read_status <- "READ_OK"
        result$object_class <- paste(class(object_value), collapse = ";")
        if (is.data.frame(object_value) || is.matrix(object_value)) {
            result$n_rows <- nrow(object_value)
            result$n_columns <- ncol(object_value)
            result$first_columns <- paste(head(colnames(object_value), 12L), collapse = ";")
        }
        return(result)
    }
    result$read_status <- "FORMAT_NOT_READ_IN_LIGHTWEIGHT_AUDIT"
    result
}

log_event("Starting", run_id)

step03_module_root <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons")

step03_gate_path <- file.path(step03_module_root, "06_results", "STEP03B3C3B_v2_GATE_TABLE.csv")

step03_lock_root <- file.path(step03_module_root, "06_results", "STEP03B3C3B_locked_RNA_occupancy_evidence", upstream_lock_version)

step03_manifest_path <- file.path(step03_lock_root, "RNA_OCCUPANCY_EVIDENCE_v1.0_LOCKED_MANIFEST.csv")

step03_routing_path <- file.path(step03_lock_root, "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv")

required_upstream_paths <- c(step03_gate_path, step03_manifest_path, step03_routing_path)

missing_upstream_paths <- required_upstream_paths[!file.exists(required_upstream_paths)]

if (length(missing_upstream_paths) > 0L) {
    stop("Required locked STEP03 artifact is missing:\n", paste(missing_upstream_paths, collapse = "\n"))
}

step03_gate <- .vdr_read_csv(step03_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

step03_manifest <- .vdr_read_csv(step03_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

step03_routing <- .vdr_read_csv(step03_routing_path, stringsAsFactors = FALSE, check.names = FALSE)

accepted_step03_gate <- step03_gate$status[step03_gate$gate_component == "STEP03B3C3B v2 gate"]

if (length(accepted_step03_gate) != 1L || accepted_step03_gate != "READY_FOR_COMBINED_TEMPORAL_PATIENT_SINGLE_CELL_PROJECTION_PACKAGE") {
    stop("STEP03 does not authorize the combined downstream package.")
}

if (nrow(step03_manifest) != 5L || !all(step03_manifest$lock_status == "LOCKED")) {
    stop("STEP03 locked manifest is incomplete.")
}

step03_manifest$file_exists <- file.exists(step03_manifest$path)

step03_manifest$observed_md5 <- vapply(step03_manifest$path, md5_file, character(1))

step03_manifest$observed_sha256 <- vapply(step03_manifest$path, sha256_file, character(1))

step03_manifest$identity_pass <- with(step03_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(step03_manifest$identity_pass)) {
    stop("At least one STEP03 locked artifact failed hash validation.")
}

if (nrow(step03_routing) != required_RNA_total || sum(step03_routing$occupancy_supported) != required_occupancy_supported || 
    sum(!step03_routing$occupancy_supported) != required_RNA_only) {
    stop("STEP03 routing table failed 436/284/152 reconciliation.")
}

dataset_registry <- data.frame(dataset_id = c("GSE250130", "GSE165004", "GSE92324", "GSE111974", "GSE183837"), 
    analysis_role = c("TEMPORAL_DECIDUALIZATION_SUPPORT", "PRIMARY_PATIENT_PROJECTION", "EXTERNAL_PATIENT_SUPPORT", 
        "EXTERNAL_PATIENT_SUPPORT_WITH_SHARED_CONTROL_CONSTRAINT", "SINGLE_CELL_LOCALIZATION"), priority = c(1L, 
        1L, 2L, 2L, 1L), required_for_main_manuscript = TRUE, primary_input_expected = c("processed expression matrix plus time/state metadata", 
        "processed endometrial expression matrix plus phenotype and slide metadata", "processed endometrial expression matrix plus phenotype metadata", 
        "processed endometrial expression matrix plus phenotype metadata and control-sharing record", "processed single-cell object or matrix plus cell metadata"), 
    locked_constraints = c("Use temporal or decidualization support; do not redefine the frozen RNA program.", 
        "Use the established slide-aware and leave-one-slide robustness framework.", "Use as external supportive validation; do not train a diagnostic classifier.", 
        "Account explicitly for controls shared with GSE165004.", "Use for cell-type/state localization; do not infer causality from localization."), 
    stringsAsFactors = FALSE)

write_csv_utf8(dataset_registry, file.path(dirs$results, "STEP04A_DATASET_ROLE_LOCK.csv"))

all_files <- list.files(project_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

all_files <- all_files[!grepl("/10_archive_work/|\\\\10_archive_work\\\\|/\\.git/|\\\\\\.git\\\\|/\\.venv|\\\\\\.venv|/renv/|\\\\renv\\\\", 
    all_files, ignore.case = TRUE, perl = TRUE)]

candidate_rows <- list()

candidate_index <- 1L

for (dataset_id in dataset_registry$dataset_id) {
    dataset_files <- all_files[grepl(dataset_id, all_files, ignore.case = TRUE, fixed = TRUE)]
    if (length(dataset_files) == 0L) {
        candidate_rows[[candidate_index]] <- data.frame(dataset_id = dataset_id, file_path = NA_character_, file_name = NA_character_, 
            extension = NA_character_, size_bytes = NA_real_, md5 = NA_character_, sha256 = NA_character_, read_status = "NO_LOCAL_FILE_FOUND", 
            object_class = NA_character_, n_rows = NA_integer_, n_columns = NA_integer_, first_columns = NA_character_, 
            inferred_file_role = "MISSING_REQUIRED_DATASET_INPUT", error_message = NA_character_, stringsAsFactors = FALSE)
        candidate_index <- candidate_index + 1L
        next
    }
    for (path_value in dataset_files) {
        extension <- tolower(tools::file_ext(path_value))
        audit <- lightweight_table_audit(path_value)
        candidate_rows[[candidate_index]] <- data.frame(dataset_id = dataset_id, file_path = normalize_slashes(path_value), 
            file_name = basename(path_value), extension = extension, size_bytes = as.numeric(file.info(path_value)$size), 
            md5 = md5_file(path_value), sha256 = sha256_file(path_value), read_status = audit$read_status, object_class = audit$object_class, 
            n_rows = audit$n_rows, n_columns = audit$n_columns, first_columns = audit$first_columns, inferred_file_role = classify_file_role(path_value, 
                audit$n_rows, audit$n_columns), error_message = audit$error_message, stringsAsFactors = FALSE)
        candidate_index <- candidate_index + 1L
    }
}

candidate_inventory <- do.call(rbind, candidate_rows)

candidate_inventory$source_candidate_status <- ifelse(candidate_inventory$read_status == "NO_LOCAL_FILE_FOUND", 
    "MISSING", "DISCOVERED_PENDING_MANUAL_SOURCE_LOCK")

write_csv_utf8(candidate_inventory, file.path(dirs$qc, "STEP04A_COMBINED_SOURCE_CANDIDATE_INVENTORY.csv"))

slide_artifacts <- all_files[grepl("STEP00|leave.?one.?slide|slide.?allocation|slide.?robust", all_files, ignore.case = TRUE, 
    perl = TRUE) & grepl("GSE165004|STEP00", all_files, ignore.case = TRUE, perl = TRUE)]

slide_inventory <- data.frame(file_path = normalize_slashes(slide_artifacts), file_name = basename(slide_artifacts), 
    size_bytes = as.numeric(file.info(slide_artifacts)$size), md5 = vapply(slide_artifacts, md5_file, character(1)), 
    sha256 = vapply(slide_artifacts, sha256_file, character(1)), role = "EXISTING_GSE165004_SLIDE_ROBUSTNESS_EVIDENCE", 
    stringsAsFactors = FALSE)

write_csv_utf8(slide_inventory, file.path(dirs$qc, "STEP04A_GSE165004_EXISTING_SLIDE_QC_INVENTORY.csv"))

readiness_rows <- lapply(dataset_registry$dataset_id, function(dataset_id) {
    dataset_candidates <- candidate_inventory[candidate_inventory$dataset_id == dataset_id & candidate_inventory$read_status != 
        "NO_LOCAL_FILE_FOUND", , drop = FALSE]
    role_counts <- table(dataset_candidates$inferred_file_role)
    expression_candidate_n <- sum(names(role_counts) %in% c("EXPRESSION_MATRIX_CANDIDATE", "TABULAR_EXPRESSION_MATRIX_CANDIDATE", 
        "BINARY_OR_SPARSE_DATA_OBJECT_CANDIDATE", "SINGLE_CELL_OBJECT_OR_MATRIX_CANDIDATE") * as.integer(role_counts))
    metadata_candidate_n <- sum(names(role_counts) %in% c("SAMPLE_METADATA_CANDIDATE", "TABULAR_METADATA_OR_RESULT_CANDIDATE") * 
        as.integer(role_counts))
    data.frame(dataset_id = dataset_id, local_file_n = nrow(dataset_candidates), expression_or_object_candidate_n = expression_candidate_n, 
        metadata_candidate_n = metadata_candidate_n, any_read_error = any(dataset_candidates$read_status == "READ_ERROR"), 
        preliminary_readiness = if (nrow(dataset_candidates) == 0L) {
            "MISSING"
        }
        else if (expression_candidate_n >= 1L) {
            "FILES_DISCOVERED_PENDING_MANUAL_SOURCE_LOCK"
        }
        else {
            "FILES_DISCOVERED_BUT_PRIMARY_MATRIX_NOT_YET_IDENTIFIED"
        }, stringsAsFactors = FALSE)
})

readiness_matrix <- do.call(rbind, readiness_rows)

readiness_matrix <- merge(dataset_registry[c("dataset_id", "analysis_role", "required_for_main_manuscript")], readiness_matrix, 
    by = "dataset_id", all.x = TRUE, sort = FALSE)

write_csv_utf8(readiness_matrix, file.path(dirs$results, "STEP04A_INPUT_READINESS_MATRIX.csv"))

execution_plan <- data.frame(package_order = 1:7, workstream = c("UPSTREAM_LOCK_AND_SOURCE_SELECTION", "TEMPORAL_DECIDUALIZATION_SUPPORT", 
    "PRIMARY_PATIENT_PROJECTION", "EXTERNAL_PATIENT_VALIDATION", "SINGLE_CELL_LOCALIZATION", "OPTIONAL_SELECTED_LOCUS_ATAC", 
    "UNIFIED_FIGURE_TABLE_AND_MANUSCRIPT_PACKAGE"), dataset_scope = c("RNA_OCCUPANCY_EVIDENCE_v1.0 plus all five dataset sources", 
    "GSE250130", "GSE165004", "GSE92324 and GSE111974", "GSE183837", "GSE306094 only after candidate narrowing", 
    "all accepted results"), primary_question = c("Are all accepted source files, roles and constraints fixed before analysis?", 
    "Which VDR-occupancy-supported RNA candidates track decidualization or temporal state?", "How does the locked program project into patient endometrium while remaining robust to slide removal?", 
    "Which patient-associated signals replicate externally without diagnostic overclaiming?", "Which endometrial cell types and states express the supported program?", 
    "Do selected final loci have descriptive accessible-chromatin context?", "Can the complete evidence chain be expressed in a stable seven-figure manuscript?"), 
    mandatory_output = c("source lock, sample-role table, hashes and exclusion log", "gene-level temporal statistics, program scores and QC", 
        "slide-aware patient projection, leave-one-slide sensitivity and QC", "external validation tables with shared-control accounting", 
        "cell-type/state localization tables and QC", "selected-locus descriptive table only if scientifically useful", 
        "submission-candidate figures, source data, legends, Methods, Results and Discussion"), direct_target_claim_allowed = FALSE, 
    stringsAsFactors = FALSE)

write_csv_utf8(execution_plan, file.path(dirs$results, "STEP04A_COMBINED_PACKAGE_EXECUTION_PLAN.csv"))

figure_plan <- data.frame(item_id = c("Figure_1", "Figure_2", "Figure_3", "Figure_4", "Figure_5", "Figure_6", "Figure_7", 
    "Table_1", "Table_2", "Supplementary_Table_1", "Supplementary_Table_2"), item_type = c(rep("MAIN_FIGURE", 7L), 
    rep("MAIN_TABLE", 2L), rep("SUPPLEMENTARY_TABLE", 2L)), scientific_question = c("What is the complete perturbation-to-patient evidence chain?", 
    "How does VDR knockdown alter the calcitriol-responsive RNA program?", "Which frozen RNA interaction candidates have VDR occupancy support?", 
    "How is the supported program related to decidualization or temporal state?", "How does the program project into the primary patient cohort and survive slide sensitivity analysis?", 
    "Which patient-associated signals are supported in external cohorts?", "Which endometrial cell types and states localize the final program?", 
    "Which datasets, samples, platforms and roles were used?", "Which final candidates carry which layers of evidence?", 
    "What are the complete 436-row RNA/occupancy results?", "What are the complete 855 peak-gene relations?"), 
    planned_content = c("workflow, datasets, locks and evidence boundaries", "RNA design, interaction contrasts, attenuation summary and selected program views", 
        "436/284/152 integration, evidence groups, promoter/shared occupancy and traceable relations", "GSE250130 temporal or decidualization trajectories and program-level support", 
        "GSE165004 patient scores, group comparisons, slide allocation and leave-one-slide robustness", "GSE92324 and GSE111974 validation with explicit control-sharing constraints", 
        "GSE183837 cell-type/state localization and final candidate distribution", "dataset registry and sample-role summary", 
        "final prioritized candidates with RNA, occupancy, temporal, patient and single-cell evidence", "locked RNA_OCCUPANCY_EVIDENCE_v1.0 table", 
        "locked relation-level subset"), current_status = c("PLANNED", "UPSTREAM_SOURCE_LOCKED_FIGURE_DEFERRED", 
        "SOURCE_DATA_LOCKED_FIGURE_DEFERRED", "AWAITING_ANALYSIS", "AWAITING_ANALYSIS", "AWAITING_ANALYSIS", "AWAITING_ANALYSIS", 
        "PLANNED", "AWAITING_FINAL_PRIORITIZATION", "SOURCE_DATA_LOCKED", "SOURCE_DATA_LOCKED"), stringsAsFactors = FALSE)

write_csv_utf8(figure_plan, file.path(dirs$results, "STEP04A_FIGURE_TABLE_PLAN.csv"))

missing_datasets <- readiness_matrix$dataset_id[readiness_matrix$preliminary_readiness == "MISSING"]

matrix_unresolved_datasets <- readiness_matrix$dataset_id[readiness_matrix$preliminary_readiness == "FILES_DISCOVERED_BUT_PRIMARY_MATRIX_NOT_YET_IDENTIFIED"]

gate_status <- if (length(missing_datasets) > 0L) {
    "HOLD_MISSING_REQUIRED_DOWNSTREAM_DATASET_FILES"
} else {
    "PROCEED_TO_MANUAL_COMBINED_INPUT_SOURCE_LOCK_AND_SINGLE_PACKAGE_ANALYSIS_BUILD"
}

gate_table <- data.frame(gate_component = c("STEP03 core lock", "STEP03 locked hashes", "STEP03 routing reconciliation", 
    "Dataset-role registry", "GSE250130 local files", "GSE165004 local files", "GSE92324 local files", "GSE111974 local files", 
    "GSE183837 local files", "GSE165004 prior slide-QC artifacts", "Combined execution plan", "Figure/table plan", 
    "Scientific figures", "STEP04A gate"), status = c("PASS_LOCKED", "PASS_5_OF_5", "PASS_436_284_152", "LOCKED_5_DATASET_ROLES", 
    readiness_matrix$preliminary_readiness[readiness_matrix$dataset_id == "GSE250130"], readiness_matrix$preliminary_readiness[readiness_matrix$dataset_id == 
        "GSE165004"], readiness_matrix$preliminary_readiness[readiness_matrix$dataset_id == "GSE92324"], readiness_matrix$preliminary_readiness[readiness_matrix$dataset_id == 
        "GSE111974"], readiness_matrix$preliminary_readiness[readiness_matrix$dataset_id == "GSE183837"], paste0("DISCOVERED_", 
        nrow(slide_inventory)), "LOCKED_AS_ONE_COORDINATED_PACKAGE", "PLANNED_NOT_RENDERED", "NONE_CREATED", gate_status), 
    required_action = c("Do not reopen STEP03 without a new version.", "Use only locked upstream artifacts.", "Retain all 436 candidates and 284/152 partition.", 
        "Use the five fixed dataset roles.", "Manually select the authoritative temporal matrix and metadata.", 
        "Manually select matrix, phenotype and slide metadata.", "Manually select the external validation matrix and metadata.", 
        "Manually select matrix and shared-control metadata.", "Manually select the single-cell object and cell metadata.", 
        "Reuse accepted slide-aware QC rather than repeating it.", "Build one analysis package after source lock.", 
        "Render figures only after all biological results are available.", "No figure skill is triggered in STEP04A.", 
        "Upload inventory, readiness matrix, role lock and gate."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP04A_GATE_TABLE.csv"))

write_text(c("# STEP04A combined package audit card", "", paste0("- Run ID: ", run_id), paste0("- Upstream evidence: ", 
    upstream_lock_version, " LOCKED"), "- Upstream routing: 436 total; 284 occupancy-supported; 152 RNA-only", 
    "- Mandatory downstream datasets: 5", paste0("- Missing datasets: ", if (length(missing_datasets) == 0L) {
        "NONE"
    } else {
        paste(missing_datasets, collapse = ";")
    }), paste0("- Datasets with files but unresolved primary matrix: ", if (length(matrix_unresolved_datasets) == 
        0L) {
        "NONE"
    } else {
        paste(matrix_unresolved_datasets, collapse = ";")
    }), paste0("- Existing GSE165004 slide-QC artifacts discovered: ", nrow(slide_inventory)), "- Combined package execution plan: LOCKED", 
    "- Figure/table plan: PLANNED", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$results, 
    "STEP04A_COMBINED_PACKAGE_AUDIT_CARD.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Active run: ", run_id), paste0("- Upstream source: ", 
    upstream_lock_version, " LOCKED"), paste0("- Routing source: ", routing_lock_version, " LOCKED"), "- Dataset roles: LOCKED", 
    "- Source files: PENDING_MANUAL_SELECTION", "- Combined execution plan: LOCKED", "- Figure/table structure: PLANNED", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, "STATUS.md"))

write_text(c("# STEP04A handoff note", "", "This step intentionally combines source discovery, dataset-role locking and manuscript planning. It does not create separate scientific substeps for minor filename or formatting issues.", 
    "", "After the authoritative files for all five datasets are selected, the next script should run the temporal, patient, external-validation and single-cell analyses as one coordinated package with shared inputs, shared candidate definitions and unified QC outputs."), 
    file.path(dirs$text, "STEP04A_Handoff_note.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, upstream_lock_version = upstream_lock_version, routing_lock_version = routing_lock_version, 
    mandatory_dataset_n = 5, missing_datasets = missing_datasets, matrix_unresolved_datasets = matrix_unresolved_datasets, 
    GSE165004_slide_QC_artifact_n = nrow(slide_inventory), execution_plan_locked = TRUE, figure_table_plan_created = TRUE, 
    scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04A_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04A.txt"))

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
    safe_file_copy(source_file, .vdr_script_path("47_step04A_combined_package_input_audit_and_plan.R"), overwrite = TRUE, 
        copy.mode = TRUE, copy.date = TRUE)
    safe_file_copy(source_file, file.path(dirs$code, "prepare_patient_source_inventory.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

archive_artifacts <- c(file.path(dirs$results, "STEP04A_DATASET_ROLE_LOCK.csv"), file.path(dirs$results, "STEP04A_INPUT_READINESS_MATRIX.csv"), 
    file.path(dirs$results, "STEP04A_COMBINED_PACKAGE_EXECUTION_PLAN.csv"), file.path(dirs$results, "STEP04A_FIGURE_TABLE_PLAN.csv"), 
    file.path(dirs$results, "STEP04A_GATE_TABLE.csv"), file.path(dirs$results, "STEP04A_COMBINED_PACKAGE_AUDIT_CARD.md"), 
    file.path(dirs$qc, "STEP04A_COMBINED_SOURCE_CANDIDATE_INVENTORY.csv"), file.path(dirs$qc, "STEP04A_GSE165004_EXISTING_SLIDE_QC_INVENTORY.csv"), 
    file.path(dirs$text, "STEP04A_Handoff_note.md"), file.path(dirs$parameters, "STEP04A_parameters.json"), file.path(dirs$environment, 
        "sessionInfo_POST_RUN_STEP04A.txt"), file.path(dirs$code, "prepare_patient_source_inventory.R"), log_path)

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
    else if (grepl("Handoff", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|INVENTORY|READINESS|AUDIT|CARD|PLAN|LOCK", basename(path_value), ignore.case = TRUE)) {
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

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04A_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04A archive status", "", paste0("- Run ID: ", run_id), paste0("- Upstream: ", upstream_lock_version, 
    " LOCKED"), "- Mandatory datasets: 5", paste0("- Missing datasets: ", length(missing_datasets)), "- Dataset roles: LOCKED", 
    "- Combined execution plan: LOCKED", "- Figure/table plan: CREATED", "- Scientific figures: NONE", paste0("- Gate: ", 
        gate_status), paste0("- Maximum archive path length: ", max(archive_manifest$path_length))), file.path(archive_root, 
    "ctl", "ARCHIVE_STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04A combined package input audit completed.\n")

cat("Upstream evidence: ", upstream_lock_version, " LOCKED\n")

cat("Mandatory downstream datasets: 5\n")

cat("Missing datasets: ", if (length(missing_datasets) == 0L) {
    "NONE"
} else {
    paste(missing_datasets, collapse = ";")
}, "\n")

cat("Datasets with unresolved primary matrix: ", if (length(matrix_unresolved_datasets) == 0L) {
    "NONE"
} else {
    paste(matrix_unresolved_datasets, collapse = ";")
}, "\n")

cat("Existing GSE165004 slide-QC artifacts: ", nrow(slide_inventory), "\n")

cat("Dataset roles: LOCKED\n")

cat("Combined execution plan: LOCKED\n")

cat("Figure/table plan: CREATED\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP04A_GATE_TABLE.csv\n")

cat("  06_results/STEP04A_INPUT_READINESS_MATRIX.csv\n")

cat("  06_results/STEP04A_DATASET_ROLE_LOCK.csv\n")

cat("  08_qc/STEP04A_COMBINED_SOURCE_CANDIDATE_INVENTORY.csv\n")

cat("  08_qc/STEP04A_GSE165004_EXISTING_SLIDE_QC_INVENTORY.csv\n")

cat("============================================================\n")

