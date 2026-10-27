#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# E01: Compare 25 model directions across independent cohorts

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260806)

project_root <- .vdr_project_root()

script_name <- "compare_single_cell_cohorts.R"

script_version <- "STEP04E1P_CROSS_DATASET_SINGLE_CELL_CONCORDANCE_v1.0"

cross_version <- "SCCROSS_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1P_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_path_safe <- function(path, max_chars = 235L) {
    normalized <- normalize_slashes(path)
    path_n <- nchar(normalized, type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Path exceeds the Windows/R safety limit (", path_n, " > ", max_chars, "):\n", normalized)
    }
    invisible(path)
}

write_csv_utf8 <- function(x, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    output <- tryCatch(.vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), 
        error = function(e) NULL)
    if (is.null(output)) {
        output <- .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    output
}

as_logical_safe <- function(x) {
    if (is.logical(x)) 
        return(x)
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

collapse_unique <- function(x, sep = "; ") {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) {
        ""
    }
    else {
        paste(x, collapse = sep)
    }
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

resolve_column <- function(data, aliases, context, required = TRUE) {
    column_names <- colnames(data)
    exact_index <- match(tolower(aliases), tolower(column_names))
    exact_index <- exact_index[!is.na(exact_index)]
    if (length(exact_index)) {
        return(column_names[exact_index[1]])
    }
    if (required) {
        stop("Could not resolve required column for ", context, ". Accepted aliases: ", paste(aliases, collapse = ", "), 
            ". Available columns: ", paste(column_names, collapse = ", "))
    }
    ""
}

extract_optional_numeric <- function(data, aliases, context) {
    column_name <- resolve_column(data = data, aliases = aliases, context = context, required = FALSE)
    if (!nzchar(column_name)) {
        rep(NA_real_, nrow(data))
    }
    else {
        safe_numeric(data[[column_name]])
    }
}

extract_optional_logical <- function(data, aliases, context) {
    column_name <- resolve_column(data = data, aliases = aliases, context = context, required = FALSE)
    if (!nzchar(column_name)) {
        rep(FALSE, nrow(data))
    }
    else {
        as_logical_safe(data[[column_name]])
    }
}

effect_direction <- function(x, tolerance = 9.9999999999999998e-13) {
    output <- rep("NONFINITE", length(x))
    output[is.finite(x) & x > tolerance] <- "RIF_HIGHER"
    output[is.finite(x) & x < -tolerance] <- "RIF_LOWER"
    output[is.finite(x) & abs(x) <= tolerance] <- "ZERO"
    output
}

safe_spearman <- function(x, y) {
    keep <- (is.finite(x) & is.finite(y))
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) {
        return(NA_real_)
    }
    suppressWarnings(stats::cor(x, y, method = "spearman"))
}

save_plot_pair <- function(plot_object, stem, width, height) {
    png_path <- file.path(figure_root, paste0(stem, ".png"))
    pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))
    assert_path_safe(png_path)
    assert_path_safe(pdf_path)
    ggplot2::ggsave(filename = png_path, plot = plot_object, width = width, height = height, units = "in", dpi = 300, 
        bg = "white")
    ggplot2::ggsave(filename = pdf_path, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    c(png_path, pdf_path)
}

log_path <- NULL

log_event <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(line, "\n")
    if (!is.null(log_path)) {
        cat(line, "\n", file = log_path, append = TRUE)
    }
    flush.console()
}

required_packages <- c("data.table", "ggplot2", "scales", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

shared_cell_types <- c("Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells")

shared_endpoints <- c("VDR_expression", "RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")

expected_grid <- expand.grid(broad_cell_type = shared_cell_types, endpoint = shared_endpoints, stringsAsFactors = FALSE)

expected_grid <- expected_grid[order(match(expected_grid$broad_cell_type, shared_cell_types), match(expected_grid$endpoint, 
    shared_endpoints)), , drop = FALSE]

expected_grid$model_key <- paste(expected_grid$broad_cell_type, expected_grid$endpoint, sep = " | ")

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

G250_model_root <- file.path(results_root, "S4E1_sc", "SCI_v1_1_memory_safe", "models", "GSE250130", "VDRMODELS_v1_1")

G250_program_path <- file.path(G250_model_root, "tables", "03_LH7_program_primary.csv")

G250_VDR_path <- file.path(G250_model_root, "tables", "07_LH7_VDR_primary.csv")

G183_model_root <- file.path(results_root, "S4E1_sc", "SCI183_v1_0", "mdl", "G183", "MD183_v1_0")

G183_model_path <- file.path(G183_model_root, "tab", "01_primary_models.csv")

G183_cluster8_path <- file.path(G183_model_root, "tab", "06_cluster8_model_sensitivity.csv")

output_root <- file.path(results_root, "S4E1_sc", "cross_dataset", cross_version)

table_root <- file.path(output_root, "tab")

source_data_root <- file.path(output_root, "src")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "man")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1P_", format(run_start, "%Y%m%d_%H%M%S"), "_crossQC"))

figure_root <- file.path(run_qc_root, "fig")

for (path_value in c(execution_root, parameter_root, environment_root, manuscript_text_root, output_root, table_root, 
    source_data_root, figure_output_root, manifest_root, run_qc_root, figure_root)) {
    dir.create(path_value, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(execution_root, paste0(run_id, "_execution.log"))

run_manifest_path <- file.path(execution_root, paste0(run_id, "_RUN_MANIFEST.csv"))

parameter_path <- file.path(parameter_root, paste0(run_id, "_parameters.json"))

session_path <- file.path(environment_root, paste0("sessionInfo_", run_id, ".txt"))

desktop_candidates <- .vdr_report_root()

desktop_candidates <- .vdr_report_root()

zip_parent <- .vdr_report_root()

handoff_zip <- file.path(zip_parent, paste0(run_id, "_SMALL_REPORT.zip"))

planned_paths <- data.frame(artifact = c("G250_program_models", "G250_VDR_models", "G183_primary_models", "G183_cluster8_sensitivity", 
    "output_root", "run_qc_root", "handoff_zip"), path = normalize_slashes(c(G250_program_path, G250_VDR_path, 
    G183_model_path, G183_cluster8_path, output_root, run_qc_root, handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

if (!all(planned_paths$path_safe)) {
    stop("One or more cross-dataset paths exceed the safe limit.")
}

G250_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

G250_qc_dirs <- G250_qc_dirs[grepl("^STEP04E1G_[0-9]{8}_[0-9]{6}_modelQC$", basename(G250_qc_dirs), perl = TRUE)]

if (!length(G250_qc_dirs)) {
    stop("No GSE250130 STEP04E1G model-QC directory was found.")
}

G250_gate_candidates <- list()

for (directory_path in G250_qc_dirs) {
    gate_path <- file.path(directory_path, "VDR_MODEL_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        G250_gate_candidates[[length(G250_gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(G250_gate_candidates)) {
    stop("No readable GSE250130 donor-model Gate was found.")
}

G250_gate_candidates <- as.data.frame(data.table::rbindlist(G250_gate_candidates, use.names = TRUE, fill = TRUE))

G250_gate_candidates$source_mtime_POSIX <- as.POSIXct(G250_gate_candidates$source_mtime)

G250_gate_candidates <- G250_gate_candidates[order(G250_gate_candidates$source_mtime_POSIX, decreasing = TRUE), 
    , drop = FALSE]

G250_authorized_index <- which(G250_gate_candidates$gate_status == "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE")

if (!length(G250_authorized_index)) {
    stop("No GSE250130 Gate passed donor-level model QC.")
}

G250_gate <- G250_gate_candidates[G250_authorized_index[1], , drop = FALSE]

G183_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

G183_qc_dirs <- G183_qc_dirs[grepl("^STEP04E1O_[0-9]{8}_[0-9]{6}_models$", basename(G183_qc_dirs), perl = TRUE)]

if (!length(G183_qc_dirs)) {
    stop("No GSE183837 STEP04E1O model-QC directory was found.")
}

G183_gate_candidates <- list()

for (directory_path in G183_qc_dirs) {
    gate_path <- file.path(directory_path, "G183_DONOR_MODEL_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        G183_gate_candidates[[length(G183_gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(G183_gate_candidates)) {
    stop("No readable GSE183837 donor-model Gate was found.")
}

G183_gate_candidates <- as.data.frame(data.table::rbindlist(G183_gate_candidates, use.names = TRUE, fill = TRUE))

G183_gate_candidates$source_mtime_POSIX <- as.POSIXct(G183_gate_candidates$source_mtime)

G183_gate_candidates <- G183_gate_candidates[order(G183_gate_candidates$source_mtime_POSIX, decreasing = TRUE), 
    , drop = FALSE]

G183_authorized_index <- which(G183_gate_candidates$gate_status == "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED" & 
    as_logical_safe(G183_gate_candidates$cross_dataset_concordance_review_authorized))

if (!length(G183_authorized_index)) {
    stop("No GSE183837 Gate authorizes cross-dataset concordance review.")
}

G183_gate <- G183_gate_candidates[G183_authorized_index[1], , drop = FALSE]

write_csv_utf8(G250_gate, file.path(run_qc_root, "00_G250_gate.csv"))

write_csv_utf8(G183_gate, file.path(run_qc_root, "00_G183_gate.csv"))

required_input_paths <- c(G250_program_path, G250_VDR_path, G183_model_path, G183_cluster8_path)

missing_input_paths <- required_input_paths[!file.exists(required_input_paths)]

if (length(missing_input_paths)) {
    stop("Missing required cross-dataset input file(s):\n", paste(normalize_slashes(missing_input_paths), collapse = "\n"))
}

G250_program_raw <- safe_read_csv(G250_program_path)

G250_VDR_raw <- safe_read_csv(G250_VDR_path)

G183_raw <- safe_read_csv(G183_model_path)

G183_cluster8 <- safe_read_csv(G183_cluster8_path)

for (data_name in c("G250_program_raw", "G250_VDR_raw", "G183_raw", "G183_cluster8")) {
    data_value <- get(data_name)
    if (is.null(data_value) || !nrow(data_value) || anyDuplicated(colnames(data_value)) > 0L) {
        stop(data_name, " is absent, empty or has duplicate columns.")
    }
}

schema_audit <- data.frame(dataset_table = c("GSE250130_LH7_program_primary", "GSE250130_LH7_VDR_primary", "GSE183837_primary_models", 
    "GSE183837_cluster8_sensitivity"), row_n = c(nrow(G250_program_raw), nrow(G250_VDR_raw), nrow(G183_raw), nrow(G183_cluster8)), 
    column_n = c(ncol(G250_program_raw), ncol(G250_VDR_raw), ncol(G183_raw), ncol(G183_cluster8)), duplicate_column_n = c(sum(duplicated(colnames(G250_program_raw))), 
        sum(duplicated(colnames(G250_VDR_raw))), sum(duplicated(colnames(G183_raw))), sum(duplicated(colnames(G183_cluster8)))), 
    column_names = c(paste(colnames(G250_program_raw), collapse = ";"), paste(colnames(G250_VDR_raw), collapse = ";"), 
        paste(colnames(G183_raw), collapse = ";"), paste(colnames(G183_cluster8), collapse = ";")), stringsAsFactors = FALSE)

write_csv_utf8(schema_audit, file.path(table_root, "00_input_schema_audit.csv"))

log_event("Starting", run_id, "G250_program_rows=", nrow(G250_program_raw), "G250_VDR_rows=", nrow(G250_VDR_raw), 
    "G183_rows=", nrow(G183_raw))

G250_program_cell_column <- resolve_column(data = G250_program_raw, aliases = c("broad_cell_type", "cell_type"), 
    context = "GSE250130 program broad cell type")

G250_program_endpoint_column <- resolve_column(data = G250_program_raw, aliases = c("program_name", "endpoint"), 
    context = "GSE250130 program endpoint")

G250_program_effect_column <- resolve_column(data = G250_program_raw, aliases = c("RIF_minus_fertile", "RIF_minus_reference", 
    "mean_difference", "estimate"), context = "GSE250130 program effect")

G250_program_p_column <- resolve_column(data = G250_program_raw, aliases = c("exact_permutation_p", "exact_permutation_P"), 
    context = "GSE250130 program exact P")

G250_program_q_column <- resolve_column(data = G250_program_raw, aliases = c("FDR_BH", "BH_FDR", "adjusted_p"), 
    context = "GSE250130 program FDR")

G250_program_Hedges_column <- resolve_column(data = G250_program_raw, aliases = c("Hedges_g", "Hedges_g_RIF_minus_Control", 
    "Hedges_g_RIF_minus_fertile"), context = "GSE250130 program Hedges g", required = FALSE)

G250_program_LOO_column <- resolve_column(data = G250_program_raw, aliases = c("LOO_sign_consistency", "LOO_direction_concordance"), 
    context = "GSE250130 program LOO", required = FALSE)

G250_program <- data.frame(dataset = "GSE250130", reference_group = "FERTILE_LH7_ALL", broad_cell_type = as.character(G250_program_raw[[G250_program_cell_column]]), 
    endpoint = as.character(G250_program_raw[[G250_program_endpoint_column]]), endpoint_class = "VDR_PROGRAM", 
    effect_RIF_minus_reference = safe_numeric(G250_program_raw[[G250_program_effect_column]]), exact_P = safe_numeric(G250_program_raw[[G250_program_p_column]]), 
    primary_q = safe_numeric(G250_program_raw[[G250_program_q_column]]), primary_q_scope = "BH_WITHIN_24_LH7_PROGRAM_MODELS", 
    Hedges_g = if (nzchar(G250_program_Hedges_column)) {
        safe_numeric(G250_program_raw[[G250_program_Hedges_column]])
    } else {
        rep(NA_real_, nrow(G250_program_raw))
    }, LOO_direction_concordance = if (nzchar(G250_program_LOO_column)) {
        safe_numeric(G250_program_raw[[G250_program_LOO_column]])
    } else {
        rep(NA_real_, nrow(G250_program_raw))
    }, n_reference = extract_optional_numeric(data = G250_program_raw, aliases = c("fertile_n", "reference_n", 
        "n_fertile", "n_FERTILE"), context = "GSE250130 program reference n"), n_RIF = extract_optional_numeric(data = G250_program_raw, 
        aliases = c("RIF_n", "n_RIF", "comparison_n"), context = "GSE250130 program RIF n"), constant_endpoint = FALSE, 
    stringsAsFactors = FALSE)

G250_VDR_cell_column <- resolve_column(data = G250_VDR_raw, aliases = c("broad_cell_type", "cell_type"), context = "GSE250130 VDR broad cell type")

G250_VDR_effect_column <- resolve_column(data = G250_VDR_raw, aliases = c("RIF_minus_fertile", "RIF_minus_reference", 
    "mean_difference", "estimate"), context = "GSE250130 VDR effect")

G250_VDR_p_column <- resolve_column(data = G250_VDR_raw, aliases = c("exact_permutation_p", "exact_permutation_P"), 
    context = "GSE250130 VDR exact P")

G250_VDR_q_column <- resolve_column(data = G250_VDR_raw, aliases = c("FDR_BH", "BH_FDR", "adjusted_p"), context = "GSE250130 VDR FDR")

G250_VDR_Hedges_column <- resolve_column(data = G250_VDR_raw, aliases = c("Hedges_g", "Hedges_g_RIF_minus_Control", 
    "Hedges_g_RIF_minus_fertile"), context = "GSE250130 VDR Hedges g", required = FALSE)

G250_VDR_LOO_column <- resolve_column(data = G250_VDR_raw, aliases = c("LOO_sign_consistency", "LOO_direction_concordance"), 
    context = "GSE250130 VDR LOO", required = FALSE)

G250_VDR <- data.frame(dataset = "GSE250130", reference_group = "FERTILE_LH7_ALL", broad_cell_type = as.character(G250_VDR_raw[[G250_VDR_cell_column]]), 
    endpoint = "VDR_expression", endpoint_class = "VDR_EXPRESSION", effect_RIF_minus_reference = safe_numeric(G250_VDR_raw[[G250_VDR_effect_column]]), 
    exact_P = safe_numeric(G250_VDR_raw[[G250_VDR_p_column]]), primary_q = safe_numeric(G250_VDR_raw[[G250_VDR_q_column]]), 
    primary_q_scope = "BH_WITHIN_6_LH7_VDR_MODELS", Hedges_g = if (nzchar(G250_VDR_Hedges_column)) {
        safe_numeric(G250_VDR_raw[[G250_VDR_Hedges_column]])
    } else {
        rep(NA_real_, nrow(G250_VDR_raw))
    }, LOO_direction_concordance = if (nzchar(G250_VDR_LOO_column)) {
        safe_numeric(G250_VDR_raw[[G250_VDR_LOO_column]])
    } else {
        rep(NA_real_, nrow(G250_VDR_raw))
    }, n_reference = extract_optional_numeric(data = G250_VDR_raw, aliases = c("fertile_n", "reference_n", "n_fertile", 
        "n_FERTILE"), context = "GSE250130 VDR reference n"), n_RIF = extract_optional_numeric(data = G250_VDR_raw, 
        aliases = c("RIF_n", "n_RIF", "comparison_n"), context = "GSE250130 VDR RIF n"), constant_endpoint = FALSE, 
    stringsAsFactors = FALSE)

G250_standardized <- rbind(G250_VDR, G250_program)

G250_standardized$direction <- effect_direction(G250_standardized$effect_RIF_minus_reference)

G250_standardized$nominal_significant <- is.finite(G250_standardized$exact_P) & G250_standardized$exact_P < 0.050000000000000003

G250_standardized$corrected_significant <- is.finite(G250_standardized$primary_q) & G250_standardized$primary_q < 
    0.050000000000000003

G250_standardized$model_key <- paste(G250_standardized$broad_cell_type, G250_standardized$endpoint, sep = " | ")

if (nrow(G250_standardized) != 30L || anyDuplicated(G250_standardized$model_key) > 0L) {
    stop("GSE250130 standardized model table is not the expected 6 x 5 grid.")
}

G250_common <- G250_standardized[G250_standardized$broad_cell_type %in% shared_cell_types & G250_standardized$endpoint %in% 
    shared_endpoints, , drop = FALSE]

if (nrow(G250_common) != 25L || !setequal(G250_common$model_key, expected_grid$model_key)) {
    stop("GSE250130 common comparison grid is incomplete.")
}

required_G183_columns <- c("broad_cell_type", "endpoint", "endpoint_class", "n_Control", "n_RIF", "mean_difference_RIF_minus_Control", 
    "exact_permutation_P", "BH_global_25", "Hedges_g_RIF_minus_Control", "LOO_direction_concordance", "constant_endpoint")

missing_G183_columns <- setdiff(required_G183_columns, colnames(G183_raw))

if (length(missing_G183_columns)) {
    stop("GSE183837 primary model table is missing: ", paste(missing_G183_columns, collapse = ", "))
}

G183_standardized <- data.frame(dataset = "GSE183837", reference_group = "CONTROL_STAGE_MATCHED", broad_cell_type = as.character(G183_raw$broad_cell_type), 
    endpoint = as.character(G183_raw$endpoint), endpoint_class = as.character(G183_raw$endpoint_class), effect_RIF_minus_reference = safe_numeric(G183_raw$mean_difference_RIF_minus_Control), 
    exact_P = safe_numeric(G183_raw$exact_permutation_P), primary_q = safe_numeric(G183_raw$BH_global_25), primary_q_scope = "GLOBAL_BH_ACROSS_25_PLANNED_MODELS", 
    Hedges_g = safe_numeric(G183_raw$Hedges_g_RIF_minus_Control), LOO_direction_concordance = safe_numeric(G183_raw$LOO_direction_concordance), 
    n_reference = safe_numeric(G183_raw$n_Control), n_RIF = safe_numeric(G183_raw$n_RIF), constant_endpoint = as_logical_safe(G183_raw$constant_endpoint), 
    stringsAsFactors = FALSE)

G183_standardized$direction <- effect_direction(G183_standardized$effect_RIF_minus_reference)

G183_standardized$nominal_significant <- is.finite(G183_standardized$exact_P) & G183_standardized$exact_P < 0.050000000000000003

G183_standardized$corrected_significant <- is.finite(G183_standardized$primary_q) & G183_standardized$primary_q < 
    0.050000000000000003

G183_standardized$model_key <- paste(G183_standardized$broad_cell_type, G183_standardized$endpoint, sep = " | ")

if (nrow(G183_standardized) != 25L || anyDuplicated(G183_standardized$model_key) > 0L || !setequal(G183_standardized$model_key, 
    expected_grid$model_key)) {
    stop("GSE183837 standardized model table is not the expected 5 x 5 grid.")
}

write_csv_utf8(G250_standardized, file.path(table_root, "01_G250_standardized_models.csv"))

write_csv_utf8(G183_standardized, file.path(table_root, "02_G183_standardized_models.csv"))

G250_join <- G250_common[, c("model_key", "broad_cell_type", "endpoint", "endpoint_class", "effect_RIF_minus_reference", 
    "exact_P", "primary_q", "primary_q_scope", "Hedges_g", "LOO_direction_concordance", "n_reference", "n_RIF", 
    "constant_endpoint", "direction", "nominal_significant", "corrected_significant"), drop = FALSE]

colnames(G250_join)[!colnames(G250_join) %in% c("model_key", "broad_cell_type", "endpoint", "endpoint_class")] <- paste0(colnames(G250_join)[!colnames(G250_join) %in% 
    c("model_key", "broad_cell_type", "endpoint", "endpoint_class")], "_G250")

G183_join <- G183_standardized[, c("model_key", "effect_RIF_minus_reference", "exact_P", "primary_q", "primary_q_scope", 
    "Hedges_g", "LOO_direction_concordance", "n_reference", "n_RIF", "constant_endpoint", "direction", "nominal_significant", 
    "corrected_significant"), drop = FALSE]

colnames(G183_join)[colnames(G183_join) != "model_key"] <- paste0(colnames(G183_join)[colnames(G183_join) != "model_key"], 
    "_G183")

cross_models <- merge(G250_join, G183_join, by = "model_key", all = TRUE, sort = FALSE)

cross_models <- merge(expected_grid, cross_models, by = c("model_key", "broad_cell_type", "endpoint"), all.x = TRUE, 
    sort = FALSE)

cross_models <- cross_models[match(expected_grid$model_key, cross_models$model_key), , drop = FALSE]

if (nrow(cross_models) != 25L || anyNA(cross_models$effect_RIF_minus_reference_G250) || anyNA(cross_models$effect_RIF_minus_reference_G183)) {
    stop("Cross-dataset model alignment failed.")
}

cross_models$noninformative_endpoint <- (cross_models$constant_endpoint_G250 | cross_models$constant_endpoint_G183 | 
    cross_models$direction_G250 %in% c("ZERO", "NONFINITE") | cross_models$direction_G183 %in% c("ZERO", "NONFINITE"))

cross_models$direction_concordant <- ifelse(cross_models$noninformative_endpoint, NA, cross_models$direction_G250 == 
    cross_models$direction_G183)

cross_models$direction_status <- ifelse(cross_models$noninformative_endpoint, "NONINFORMATIVE", ifelse(cross_models$direction_concordant, 
    "CONCORDANT", "DISCORDANT"))

cross_models$evidence_class <- ifelse(cross_models$noninformative_endpoint, "NONINFORMATIVE_CONSTANT_OR_ZERO", 
    ifelse(cross_models$direction_concordant & cross_models$corrected_significant_G250 & cross_models$corrected_significant_G183, 
        "REPLICATED_CORRECTED_CONCORDANT", ifelse(cross_models$direction_concordant & (cross_models$corrected_significant_G250 | 
            cross_models$corrected_significant_G183) & (cross_models$nominal_significant_G250 & cross_models$nominal_significant_G183), 
            "CONCORDANT_ONE_CORRECTED_BOTH_NOMINAL", ifelse(cross_models$direction_concordant & cross_models$nominal_significant_G250 & 
                cross_models$nominal_significant_G183, "CONCORDANT_BOTH_NOMINAL", ifelse(cross_models$direction_concordant & 
                (cross_models$nominal_significant_G250 | cross_models$nominal_significant_G183), "CONCORDANT_ONE_NOMINAL_EXPLORATORY", 
                ifelse(cross_models$direction_concordant, "CONCORDANT_NO_STATISTICAL_SUPPORT", ifelse(cross_models$nominal_significant_G250 | 
                  cross_models$nominal_significant_G183 | cross_models$corrected_significant_G250 | cross_models$corrected_significant_G183, 
                  "DISCORDANT_WITH_NOMINAL_OR_CORRECTED_SIGNAL", "DISCORDANT_NO_STATISTICAL_SUPPORT")))))))

cross_models$cross_dataset_interpretation <- ifelse(cross_models$evidence_class == "REPLICATED_CORRECTED_CONCORDANT", 
    "CORRECTED_REPLICATION", ifelse(cross_models$evidence_class %in% c("CONCORDANT_ONE_CORRECTED_BOTH_NOMINAL", 
        "CONCORDANT_BOTH_NOMINAL", "CONCORDANT_ONE_NOMINAL_EXPLORATORY"), "EXPLORATORY_DIRECTIONAL_SUPPORT", ifelse(cross_models$evidence_class == 
        "CONCORDANT_NO_STATISTICAL_SUPPORT", "DIRECTION_ONLY_NO_INFERENTIAL_SUPPORT", ifelse(cross_models$evidence_class == 
        "NONINFORMATIVE_CONSTANT_OR_ZERO", "NONINFORMATIVE", "NO_CROSS_DATASET_SUPPORT"))))

write_csv_utf8(cross_models, file.path(table_root, "03_cross_dataset_25_model_concordance.csv"))

informative_models <- cross_models[!cross_models$noninformative_endpoint, , drop = FALSE]

overall_summary <- data.frame(metric = c("planned_common_model_n", "informative_model_n", "noninformative_model_n", 
    "direction_concordant_n", "direction_discordant_n", "direction_concordance_fraction", "G250_nominal_P_lt_0_05_n", 
    "G250_corrected_q_lt_0_05_n", "G183_nominal_P_lt_0_05_n", "G183_corrected_q_lt_0_05_n", "corrected_replicated_concordant_n", 
    "exploratory_directional_support_n", "cluster8_direction_changed_n"), value = c(nrow(cross_models), nrow(informative_models), 
    sum(cross_models$noninformative_endpoint), sum(informative_models$direction_concordant), sum(!informative_models$direction_concordant), 
    mean(informative_models$direction_concordant), sum(cross_models$nominal_significant_G250), sum(cross_models$corrected_significant_G250), 
    sum(cross_models$nominal_significant_G183), sum(cross_models$corrected_significant_G183), sum(cross_models$evidence_class == 
        "REPLICATED_CORRECTED_CONCORDANT"), sum(cross_models$cross_dataset_interpretation == "EXPLORATORY_DIRECTIONAL_SUPPORT"), 
    sum(G183_cluster8$direction_concordant %in% FALSE, na.rm = TRUE)), stringsAsFactors = FALSE)

by_endpoint_rows <- list()

for (endpoint_index in seq_along(shared_endpoints)) {
    endpoint_value <- shared_endpoints[endpoint_index]
    endpoint_data <- cross_models[cross_models$endpoint == endpoint_value, , drop = FALSE]
    endpoint_informative <- endpoint_data[!endpoint_data$noninformative_endpoint, , drop = FALSE]
    by_endpoint_rows[[endpoint_index]] <- data.frame(endpoint = endpoint_value, planned_model_n = nrow(endpoint_data), 
        informative_model_n = nrow(endpoint_informative), concordant_n = sum(endpoint_informative$direction_concordant), 
        discordant_n = sum(!endpoint_informative$direction_concordant), concordance_fraction = if (nrow(endpoint_informative)) {
            mean(endpoint_informative$direction_concordant)
        }
        else {
            NA_real_
        }, Spearman_effect_rho = safe_spearman(endpoint_informative$effect_RIF_minus_reference_G250, endpoint_informative$effect_RIF_minus_reference_G183), 
        G250_nominal_n = sum(endpoint_data$nominal_significant_G250), G183_nominal_n = sum(endpoint_data$nominal_significant_G183), 
        corrected_replication_n = sum(endpoint_data$evidence_class == "REPLICATED_CORRECTED_CONCORDANT"), stringsAsFactors = FALSE)
}

by_endpoint <- do.call(rbind, by_endpoint_rows)

by_cell_type_rows <- list()

for (cell_index in seq_along(shared_cell_types)) {
    cell_type <- shared_cell_types[cell_index]
    cell_data <- cross_models[cross_models$broad_cell_type == cell_type, , drop = FALSE]
    cell_informative <- cell_data[!cell_data$noninformative_endpoint, , drop = FALSE]
    by_cell_type_rows[[cell_index]] <- data.frame(broad_cell_type = cell_type, planned_model_n = nrow(cell_data), 
        informative_model_n = nrow(cell_informative), concordant_n = sum(cell_informative$direction_concordant), 
        discordant_n = sum(!cell_informative$direction_concordant), concordance_fraction = if (nrow(cell_informative)) {
            mean(cell_informative$direction_concordant)
        }
        else {
            NA_real_
        }, G250_nominal_n = sum(cell_data$nominal_significant_G250), G183_nominal_n = sum(cell_data$nominal_significant_G183), 
        corrected_replication_n = sum(cell_data$evidence_class == "REPLICATED_CORRECTED_CONCORDANT"), stringsAsFactors = FALSE)
}

by_cell_type <- do.call(rbind, by_cell_type_rows)

evidence_class_summary <- as.data.frame(table(evidence_class = cross_models$evidence_class), stringsAsFactors = FALSE)

colnames(evidence_class_summary)[2] <- "model_n"

myeloid_VDR_focus <- cross_models[cross_models$broad_cell_type == "Myeloid cells" & cross_models$endpoint == "VDR_expression", 
    , drop = FALSE]

if (nrow(myeloid_VDR_focus) != 1L) {
    stop("Myeloid VDR focus row was not uniquely recovered.")
}

write_csv_utf8(overall_summary, file.path(table_root, "04_overall_concordance_summary.csv"))

write_csv_utf8(by_endpoint, file.path(table_root, "05_concordance_by_endpoint.csv"))

write_csv_utf8(by_cell_type, file.path(table_root, "06_concordance_by_cell_type.csv"))

write_csv_utf8(evidence_class_summary, file.path(table_root, "07_evidence_class_summary.csv"))

write_csv_utf8(myeloid_VDR_focus, file.path(table_root, "08_myeloid_VDR_focus.csv"))

write_csv_utf8(G183_cluster8, file.path(table_root, "09_cluster8_sensitivity_review.csv"))

corrected_replication_n <- sum(cross_models$evidence_class == "REPLICATED_CORRECTED_CONCORDANT")

exploratory_support_n <- sum(cross_models$cross_dataset_interpretation == "EXPLORATORY_DIRECTIONAL_SUPPORT")

myeloid_VDR_G183_nominal <- (myeloid_VDR_focus$nominal_significant_G183)

myeloid_VDR_G183_corrected <- (myeloid_VDR_focus$corrected_significant_G183)

myeloid_VDR_direction_status <- (myeloid_VDR_focus$direction_status)

single_cell_evidence_state <- if (corrected_replication_n > 0L) {
    "CORRECTED_CROSS_DATASET_REPLICATION_PRESENT"
} else if (exploratory_support_n > 0L) {
    "NO_CORRECTED_REPLICATION_EXPLORATORY_DIRECTIONAL_SUPPORT_ONLY"
} else {
    "NO_CORRECTED_REPLICATION_OR_DIRECTIONAL_SUPPORT"
}

broad_RIF_VDR_program_dysregulation_supported <- (corrected_replication_n > 0L && any(cross_models$endpoint_class == 
    "VDR_PROGRAM" & cross_models$evidence_class == "REPLICATED_CORRECTED_CONCORDANT"))

interpretation_lock <- data.frame(item = c("single_cell_evidence_state", "corrected_cross_dataset_replication_n", 
    "exploratory_directional_support_n", "broad_RIF_VDR_program_dysregulation_supported", "GSE183837_myeloid_VDR_nominal_signal", 
    "GSE183837_myeloid_VDR_corrected_signal", "myeloid_VDR_cross_dataset_direction_status", "cluster8_sensitivity_direction_changed_n", 
    "GSE183837_ciliated_epithelial_inference", "pooled_meta_analysis_performed", "direct_VDR_target_claim_authorized", 
    "causal_mechanism_claim_authorized", "single_cell_component_interpretation_locked", "project_level_result_lock_authorized"), 
    value = c(single_cell_evidence_state, corrected_replication_n, exploratory_support_n, broad_RIF_VDR_program_dysregulation_supported, 
        myeloid_VDR_G183_nominal, myeloid_VDR_G183_corrected, myeloid_VDR_direction_status, sum(G183_cluster8$direction_concordant %in% 
            FALSE, na.rm = TRUE), "DESCRIPTIVE_ONLY", FALSE, FALSE, FALSE, TRUE, FALSE), stringsAsFactors = FALSE)

write_csv_utf8(interpretation_lock, file.path(table_root, "10_single_cell_interpretation_lock.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

scatter_data <- cross_models[!cross_models$noninformative_endpoint, , drop = FALSE]

scatter_data$endpoint <- factor(scatter_data$endpoint, levels = shared_endpoints)

short_cell_label_map <- c(`Endothelial cells` = "Endo", `Myeloid cells` = "Myeloid", `NK/T cells` = "NK/T", `Stromal cells` = "Stromal", 
    `Unciliated epithelial cells` = "Unciliated")

scatter_data$short_cell_label <- unname(short_cell_label_map[scatter_data$broad_cell_type])

if (anyNA(scatter_data$short_cell_label)) {
    stop("One or more broad cell types could not be mapped to a short plot label.")
}

plot_scatter <- ggplot2::ggplot(scatter_data, ggplot2::aes(x = effect_RIF_minus_reference_G250, y = effect_RIF_minus_reference_G183, 
    shape = direction_status)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::geom_hline(yintercept = 0, 
    linetype = "dashed") + ggplot2::geom_point(size = 2.7999999999999998) + ggplot2::geom_text(ggplot2::aes(label = short_cell_label), 
    nudge_y = 0.029999999999999999, size = 2.6000000000000001, check_overlap = TRUE) + ggplot2::facet_wrap(~endpoint, 
    scales = "free", ncol = 3) + ggplot2::labs(title = "Cross-dataset RIF effect directions", subtitle = "GSE250130 compares RIF with fertile LH+7; GSE183837 compares RIF with stage-matched controls", 
    x = "GSE250130 effect: RIF minus fertile", y = "GSE183837 effect: RIF minus Control", shape = "Direction status") + 
    theme_QC

save_plot_pair(plot_scatter, "QC1_cross_dataset_effect_scatter", width = 11, height = 8)

matrix_data <- cross_models

matrix_data$broad_cell_type <- factor(matrix_data$broad_cell_type, levels = rev(shared_cell_types))

matrix_data$endpoint <- factor(matrix_data$endpoint, levels = shared_endpoints)

matrix_data$direction_symbol <- ifelse(matrix_data$direction_status == "CONCORDANT", "C", ifelse(matrix_data$direction_status == 
    "DISCORDANT", "D", "NI"))

plot_direction_matrix <- ggplot2::ggplot(matrix_data, ggplot2::aes(x = endpoint, y = broad_cell_type, fill = direction_status)) + 
    ggplot2::geom_tile(linewidth = 0.40000000000000002) + ggplot2::geom_text(ggplot2::aes(label = direction_symbol), 
    size = 3.5) + ggplot2::scale_fill_manual(values = c(CONCORDANT = "grey35", DISCORDANT = "grey75", NONINFORMATIVE = "white")) + 
    ggplot2::labs(title = "Direction concordance across the 25 shared models", subtitle = "C, concordant; D, discordant; NI, non-informative constant or zero endpoint", 
        x = "Endpoint", y = "Broad cell type", fill = "Status") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_direction_matrix, "QC2_direction_concordance_matrix", width = 10, height = 6)

Hedges_long <- rbind(data.frame(model_key = cross_models$model_key, broad_cell_type = cross_models$broad_cell_type, 
    endpoint = cross_models$endpoint, dataset = "GSE250130", Hedges_g = cross_models$Hedges_g_G250, stringsAsFactors = FALSE), 
    data.frame(model_key = cross_models$model_key, broad_cell_type = cross_models$broad_cell_type, endpoint = cross_models$endpoint, 
        dataset = "GSE183837", Hedges_g = cross_models$Hedges_g_G183, stringsAsFactors = FALSE))

Hedges_long <- Hedges_long[is.finite(Hedges_long$Hedges_g), , drop = FALSE]

Hedges_long$endpoint <- factor(Hedges_long$endpoint, levels = shared_endpoints)

plot_Hedges <- ggplot2::ggplot(Hedges_long, ggplot2::aes(x = dataset, y = Hedges_g, group = broad_cell_type, shape = broad_cell_type)) + 
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + ggplot2::geom_line(alpha = 0.45000000000000001) + 
    ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~endpoint, scales = "free_y", ncol = 3) + ggplot2::labs(title = "Standardized effect-size comparison", 
    subtitle = "Hedges g is oriented as RIF minus reference in both studies", x = "Dataset", y = "Hedges g", shape = "Broad cell type") + 
    theme_QC

save_plot_pair(plot_Hedges, "QC3_Hedges_g_cross_dataset", width = 11, height = 8)

q_long <- rbind(data.frame(broad_cell_type = cross_models$broad_cell_type, endpoint = cross_models$endpoint, dataset = "GSE250130", 
    primary_q = cross_models$primary_q_G250, stringsAsFactors = FALSE), data.frame(broad_cell_type = cross_models$broad_cell_type, 
    endpoint = cross_models$endpoint, dataset = "GSE183837", primary_q = cross_models$primary_q_G183, stringsAsFactors = FALSE))

q_long$minus_log10_q <- -log10(pmax(q_long$primary_q, 9.9999999999999998e-13))

q_long$model_label <- paste(q_long$broad_cell_type, q_long$endpoint, sep = " | ")

q_long$model_label <- factor(q_long$model_label, levels = rev(expected_grid$model_key))

plot_q <- ggplot2::ggplot(q_long, ggplot2::aes(x = dataset, y = model_label, fill = minus_log10_q)) + ggplot2::geom_tile() + 
    ggplot2::scale_fill_viridis_c() + ggplot2::labs(title = "Primary multiple-testing-adjusted evidence", subtitle = "Adjustment scopes differ by the prespecified model families in each study", 
    x = "Dataset", y = "Shared model", fill = "-log10(q)") + theme_QC + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))

save_plot_pair(plot_q, "QC4_adjusted_evidence_matrix", width = 9, height = 10)

myeloid_plot_data <- data.frame(dataset = c("GSE250130", "GSE183837"), effect = c(myeloid_VDR_focus$effect_RIF_minus_reference_G250, 
    myeloid_VDR_focus$effect_RIF_minus_reference_G183), Hedges_g = c(myeloid_VDR_focus$Hedges_g_G250, myeloid_VDR_focus$Hedges_g_G183), 
    exact_P = c(myeloid_VDR_focus$exact_P_G250, myeloid_VDR_focus$exact_P_G183), primary_q = c(myeloid_VDR_focus$primary_q_G250, 
        myeloid_VDR_focus$primary_q_G183), stringsAsFactors = FALSE)

myeloid_plot_data$label <- paste0("P=", format(myeloid_plot_data$exact_P, digits = 3), "\nq=", format(myeloid_plot_data$primary_q, 
    digits = 3))

plot_myeloid <- ggplot2::ggplot(myeloid_plot_data, ggplot2::aes(x = dataset, y = effect, shape = dataset)) + ggplot2::geom_hline(yintercept = 0, 
    linetype = "dashed") + ggplot2::geom_point(size = 3.2000000000000002) + ggplot2::geom_text(ggplot2::aes(label = label), 
    nudge_y = 0.14999999999999999, size = 3) + ggplot2::labs(title = "Myeloid-cell VDR expression across the two studies", 
    subtitle = "The GSE183837 nominal signal remains exploratory unless supported in direction by GSE250130", x = "Dataset", 
    y = "RIF minus reference VDR log2(CPM + 1)", shape = "Dataset") + theme_QC

save_plot_pair(plot_myeloid, "QC5_myeloid_VDR_focus", width = 7, height = 5.5)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

myeloid_G183_exact_P <- myeloid_VDR_focus$exact_P_G183

myeloid_G183_q <- myeloid_VDR_focus$primary_q_G183

cluster8_direction_changed_n <- sum(G183_cluster8$direction_concordant %in% FALSE, na.rm = TRUE)

gate_checks <- data.frame(check = c("GSE250130_model_gate_passes", "GSE183837_model_gate_passes", "GSE250130_program_table_has_24_models", 
    "GSE250130_VDR_table_has_6_models", "GSE250130_common_grid_has_25_models", "GSE183837_grid_has_25_models", 
    "shared_grid_alignment_is_complete", "all_effects_are_finite", "all_exact_P_values_are_finite", "all_primary_q_values_are_finite", 
    "one_GSE183837_constant_endpoint_is_endothelial_VDR", "GSE183837_myeloid_VDR_is_nominal_but_not_corrected", 
    "cluster8_sensitivity_has_5_endpoints", "cluster8_sensitivity_has_no_direction_change", "no_pooled_meta_analysis_performed", 
    "cross_dataset_interpretation_can_be_locked", "project_level_result_lock_remains_disabled"), pass = c(G250_gate$gate_status[[1]] == 
    "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE", G183_gate$gate_status[[1]] == "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED", 
    nrow(G250_program_raw) == 24L, nrow(G250_VDR_raw) == 6L, nrow(G250_common) == 25L, nrow(G183_standardized) == 
        25L, nrow(cross_models) == 25L && !anyNA(cross_models$model_key), all(is.finite(cross_models$effect_RIF_minus_reference_G250)) && 
        all(is.finite(cross_models$effect_RIF_minus_reference_G183)), all(is.finite(cross_models$exact_P_G250)) && 
        all(is.finite(cross_models$exact_P_G183)), all(is.finite(cross_models$primary_q_G250)) && all(is.finite(cross_models$primary_q_G183)), 
    sum(cross_models$constant_endpoint_G183) == 1L && cross_models$broad_cell_type[cross_models$constant_endpoint_G183] == 
        "Endothelial cells" && cross_models$endpoint[cross_models$constant_endpoint_G183] == "VDR_expression", 
    is.finite(myeloid_G183_exact_P) && myeloid_G183_exact_P < 0.050000000000000003 && is.finite(myeloid_G183_q) && 
        myeloid_G183_q >= 0.050000000000000003, nrow(G183_cluster8) == 5L, cluster8_direction_changed_n == 0L, 
    TRUE, TRUE, TRUE), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_SINGLE_CELL_CROSS_DATASET_CONCORDANCE_REVIEW_COMPLETE"
} else {
    "HOLD_SINGLE_CELL_CROSS_DATASET_CONCORDANCE_INCOMPLETE"
}

cross_gate <- data.frame(cross_version = cross_version, gate_status = gate_status, planned_common_model_n = nrow(cross_models), 
    informative_model_n = nrow(informative_models), direction_concordant_n = sum(informative_models$direction_concordant), 
    direction_discordant_n = sum(!informative_models$direction_concordant), direction_concordance_fraction = mean(informative_models$direction_concordant), 
    corrected_replication_n = corrected_replication_n, exploratory_directional_support_n = exploratory_support_n, 
    G250_corrected_signal_n = sum(cross_models$corrected_significant_G250), G183_corrected_signal_n = sum(cross_models$corrected_significant_G183), 
    G183_nominal_signal_n = sum(cross_models$nominal_significant_G183), G183_myeloid_VDR_exact_P = myeloid_G183_exact_P, 
    G183_myeloid_VDR_global_q = myeloid_G183_q, myeloid_VDR_direction_status = myeloid_VDR_direction_status, cluster8_direction_changed_n = cluster8_direction_changed_n, 
    single_cell_evidence_state = single_cell_evidence_state, broad_RIF_VDR_program_dysregulation_supported = broad_RIF_VDR_program_dysregulation_supported, 
    single_cell_component_interpretation_locked = all_gate_checks_pass, project_level_result_lock_authorized = FALSE, 
    gene_level_cross_dataset_DE_authorized = FALSE, pooled_meta_analysis_authorized = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(cross_gate, file.path(run_qc_root, "SC_CROSS_DATASET_GATE.csv"))

card_lines <- c("# STEP04E1P Single-cell Cross-dataset Concordance Card", "", paste0("- Gate: ", gate_status), 
    "- Shared planned models: 25", paste0("- Informative models: ", nrow(informative_models)), paste0("- Direction-concordant models: ", 
        sum(informative_models$direction_concordant)), paste0("- Direction-discordant models: ", sum(!informative_models$direction_concordant)), 
    paste0("- Direction concordance fraction: ", format(mean(informative_models$direction_concordant), digits = 4)), 
    paste0("- Corrected replicated models: ", corrected_replication_n), paste0("- Exploratory directional-support models: ", 
        exploratory_support_n), paste0("- GSE183837 myeloid VDR exact P: ", format(myeloid_G183_exact_P, digits = 4)), 
    paste0("- GSE183837 myeloid VDR global q: ", format(myeloid_G183_q, digits = 4)), paste0("- Myeloid VDR direction status: ", 
        myeloid_VDR_direction_status), paste0("- Single-cell evidence state: ", single_cell_evidence_state), paste0("- Broad RIF VDR-program dysregulation supported: ", 
        broad_RIF_VDR_program_dysregulation_supported), paste0("- Single-cell component interpretation locked: ", 
        all_gate_checks_pass), "- Project-level result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "SC_CROSS_DATASET_CARD.md"))

interpretation_lines <- c("# Single-cell cross-dataset interpretation", "", paste0("- Run ID: `", run_id, "`"), 
    paste0("- Gate: **", gate_status, "**"), "", "## Locked interpretation", "", paste0("- ", nrow(informative_models), 
        " of 25 shared models were directionally informative; ", sum(informative_models$direction_concordant), 
        " were concordant and ", sum(!informative_models$direction_concordant), " were discordant."), paste0("- Corrected cross-dataset replication was observed for ", 
        corrected_replication_n, " model(s)."), paste0("- Exploratory directional support was observed for ", exploratory_support_n, 
        " model(s)."), paste0("- The GSE183837 myeloid VDR-expression signal had exact P=", format(myeloid_G183_exact_P, 
        digits = 4), " and global q=", format(myeloid_G183_q, digits = 4), "; its direction relative to GSE250130 was ", 
        myeloid_VDR_direction_status, "."), paste0("- The final single-cell evidence state is: ", single_cell_evidence_state, 
        "."), "", "## Scientific boundary", "", "- The two studies were not pooled because their reference groups and sampling structures differ.", 
    "- Directional concordance without corrected significance is descriptive rather than confirmatory.", "- The independent GSE183837 cohort does not by itself establish broad RIF-associated VDR-program dysregulation.", 
    "- The myeloid VDR-expression finding remains exploratory because it did not survive global multiplicity correction.", 
    "- No direct-target, VDR-activity or causal mechanism claim is authorized.", "", "## Authorized next operation", 
    "", "- Use this locked single-cell interpretation in the STEP04 closeout.", "- Integrate the temporal single-cell evidence, patient bulk projection and occupancy hierarchy.", 
    "- Finalize the manuscript Figure/Table blueprint and wording QC.", "", "## Not authorized", "", "- No pooled meta-analysis.", 
    "- No additional cell-level testing.", "- No GSE183837 gene-level discovery DE.", "- No project-level final result lock.")

write_text(interpretation_lines, file.path(run_qc_root, "SC_CROSS_DATASET_INTERPRETATION.md"))

methods_lines <- c("# STEP04E1P Methods shell", "", "Cross-dataset consistency was evaluated across the 25 combinations shared by GSE250130 and GSE183837, comprising five broad cell types and five endpoints. Effects were oriented as RIF minus reference in both studies. GSE250130 contrasted RIF with all fertile LH+7 donors, whereas GSE183837 contrasted RIF with stage-matched controls. Direction concordance was classified only for finite non-zero effects. Each study retained its prospectively defined multiplicity correction; no pooled meta-analysis or cross-study combined P value was calculated because the reference groups, sample sizes and cohort structures differed.")

results_lines <- c("# STEP04E1P Results shell", "", paste0("Among ", nrow(informative_models), " directionally informative shared models, ", 
    sum(informative_models$direction_concordant), " showed the same effect direction across studies and ", sum(!informative_models$direction_concordant), 
    " were discordant."), paste0("Corrected cross-dataset replication was observed in ", corrected_replication_n, 
    " model(s)."), paste0("In GSE183837, myeloid-cell VDR expression was lower in RIF at nominal exact P=", format(myeloid_G183_exact_P, 
    digits = 4), ", but the global BH-adjusted q value was ", format(myeloid_G183_q, digits = 4), ". The corresponding effect direction in GSE250130 was ", 
    myeloid_VDR_direction_status, "."), paste0("Retaining the low-confidence cluster 8 changed the unciliated epithelial effect direction for ", 
    cluster8_direction_changed_n, " of five endpoints."), "Together, the two donor-level single-cell analyses did not support a broad, reproducible RIF-associated shift in the locked VDR transcriptional programs.")

discussion_lines <- c("# STEP04E1P Discussion shell", "", "The independent cohort primarily provided a falsification and localization test rather than strong confirmatory power. The absence of multiplicity-corrected replication argues against framing RIF as a uniform VDR-program deficiency across endometrial cell compartments. The nominal myeloid VDR-expression decrease is biologically interesting because its donor-level effect was large and leave-one-donor directionally stable in GSE183837, but it remains exploratory until reproduced in a larger stage-matched cohort. The strong temporal remodeling observed in GSE250130 should therefore be separated from disease-associated dysregulation: normal luteal timing clearly structures the VDR-related transcriptional context, whereas broad RIF-associated disruption was not consistently demonstrated.")

write_text(methods_lines, file.path(manuscript_text_root, paste0(run_id, "_Methods.md")))

write_text(results_lines, file.path(manuscript_text_root, paste0(run_id, "_Results.md")))

write_text(discussion_lines, file.path(manuscript_text_root, paste0(run_id, "_Discussion.md")))

parameter_manifest <- list(script_name = script_name, script_version = script_version, cross_version = cross_version, 
    GSE250130_reference = "ALL_FERTILE_LH7", GSE183837_reference = "STAGE_MATCHED_CONTROL", shared_cell_types = shared_cell_types, 
    shared_endpoints = shared_endpoints, effect_orientation = "RIF_MINUS_REFERENCE", concordance_rule = "SAME_NONZERO_EFFECT_DIRECTION", 
    GSE250130_primary_multiplicity = "PROGRAMS_BH_WITHIN_24_AND_VDR_BH_WITHIN_6", GSE183837_primary_multiplicity = "GLOBAL_BH_ACROSS_25", 
    pooled_meta_analysis = FALSE, corrected_replication_rule = "BOTH_STUDIES_CORRECTED_SIGNIFICANT_AND_DIRECTION_CONCORDANT", 
    single_cell_component_interpretation_locked = all_gate_checks_pass, project_level_result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, cross_version = cross_version, 
    started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, 
        run_start, units = "mins")), 2), gate_status = gate_status, planned_common_model_n = nrow(cross_models), 
    informative_model_n = nrow(informative_models), direction_concordant_n = sum(informative_models$direction_concordant), 
    direction_discordant_n = sum(!informative_models$direction_concordant), corrected_replication_n = corrected_replication_n, 
    exploratory_directional_support_n = exploratory_support_n, myeloid_VDR_direction_status = myeloid_VDR_direction_status, 
    single_cell_component_interpretation_locked = all_gate_checks_pass, project_level_result_lock_authorized = FALSE, 
    output_root = normalize_slashes(output_root), qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

persistent_files <- list.files(output_root, recursive = TRUE, full.names = TRUE)

persistent_files <- persistent_files[file.exists(persistent_files) & !dir.exists(persistent_files)]

persistent_manifest <- data.frame(path = normalize_slashes(persistent_files), size_bytes = file.info(persistent_files)$size, 
    md5 = unname(tools::md5sum(persistent_files)), stringsAsFactors = FALSE)

write_csv_utf8(persistent_manifest, file.path(manifest_root, "persistent_manifest.csv"))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(source_data_root, 
    recursive = TRUE, full.names = TRUE), list.files(figure_output_root, recursive = TRUE, full.names = TRUE), 
    list.files(manifest_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, recursive = TRUE, full.names = TRUE), 
    log_path, run_manifest_path, parameter_path, session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", 
        mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    100 * 1024^2]

artifact_manifest <- data.frame(path = normalize_slashes(small_files), size_bytes = file.info(small_files)$size, 
    md5 = unname(tools::md5sum(small_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

log_event("Completed", run_id, "gate_status=", gate_status, "informative=", nrow(informative_models), "concordant=", 
    sum(informative_models$direction_concordant), "corrected_replication=", corrected_replication_n)

if (!requireNamespace("zip", quietly = TRUE)) {
    try(install.packages("zip", repos = "https://cloud.r-project.org"), silent = TRUE)
}

zip_ok <- FALSE

zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
    zip_ok <- tryCatch({
        zip::zipr(zipfile = handoff_zip, files = small_files, recurse = FALSE, include_directories = FALSE)
        file.exists(handoff_zip)
    }, error = function(e) {
        zip_error <<- conditionMessage(e)
        FALSE
    })
}

cat("\n============================================================\n")

cat("STEP04E1P single-cell cross-dataset concordance completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Shared planned models: ", nrow(cross_models), "\n", sep = "")

cat("Informative models: ", nrow(informative_models), "\n", sep = "")

cat("Direction-concordant models: ", sum(informative_models$direction_concordant), "\n", sep = "")

cat("Direction-discordant models: ", sum(!informative_models$direction_concordant), "\n", sep = "")

cat("Corrected replicated models: ", corrected_replication_n, "\n", sep = "")

cat("Exploratory directional-support models: ", exploratory_support_n, "\n", sep = "")

cat("Myeloid VDR direction status: ", myeloid_VDR_direction_status, "\n", sep = "")

cat("Single-cell evidence state: ", single_cell_evidence_state, "\n", sep = "")

cat("Single-cell component interpretation locked: ", all_gate_checks_pass, "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nCross-dataset outputs are complete in the project folders.\n", sep = "")
}

