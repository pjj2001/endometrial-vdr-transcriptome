#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D07: Fit independent donor-level models

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "fit_gse183837_donor_models.R"

script_version <- "STEP04E1O_GSE183837_DONOR_MODELS_v1.1_TYPE_SAFE_SUMMARY"

integration_version <- "SCI183_v1_0"

pseudobulk_version <- "PB183_v1_0"

model_version <- "MD183_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1O_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

safe_mean <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        mean(x)
    }
}

safe_median <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        median(x)
    }
}

safe_sd <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 2L) {
        NA_real_
    }
    else {
        stats::sd(x)
    }
}

safe_min <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        min(x)
    }
}

safe_max <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) {
        NA_real_
    }
    else {
        max(x)
    }
}

hedges_g <- function(Control_values, RIF_values) {
    Control_values <- as.numeric(Control_values)
    RIF_values <- as.numeric(RIF_values)
    n_Control <- length(Control_values)
    n_RIF <- length(RIF_values)
    if (n_Control < 2L || n_RIF < 2L) {
        return(NA_real_)
    }
    sd_Control <- stats::sd(Control_values)
    sd_RIF <- stats::sd(RIF_values)
    df <- n_Control + n_RIF - 2L
    pooled_variance <- ((n_Control - 1L) * sd_Control^2 + (n_RIF - 1L) * sd_RIF^2)/df
    if (!is.finite(pooled_variance) || pooled_variance <= 0) {
        return(NA_real_)
    }
    Cohen_d <- (mean(RIF_values) - mean(Control_values))/sqrt(pooled_variance)
    correction_J <- 1 - 3/(4 * df - 1)
    correction_J * Cohen_d
}

rank_biserial_RIF_minus_Control <- function(Control_values, RIF_values) {
    Control_values <- as.numeric(Control_values)
    RIF_values <- as.numeric(RIF_values)
    n_Control <- length(Control_values)
    n_RIF <- length(RIF_values)
    if (n_Control < 1L || n_RIF < 1L) {
        return(NA_real_)
    }
    pooled_values <- c(RIF_values, Control_values)
    pooled_ranks <- rank(pooled_values, ties.method = "average")
    RIF_rank_sum <- sum(pooled_ranks[seq_len(n_RIF)])
    U_RIF <- RIF_rank_sum - n_RIF * (n_RIF + 1)/2
    2 * U_RIF/(n_RIF * n_Control) - 1
}

exact_permutation_mean_difference <- function(values, groups, model_id, tolerance = 9.9999999999999998e-13) {
    values <- as.numeric(values)
    groups <- as.character(groups)
    if (length(values) != length(groups) || any(!is.finite(values)) || any(!groups %in% c("CONTROL", "RIF"))) {
        stop("Invalid values or groups for exact permutation model ", model_id, ".")
    }
    Control_index <- which(groups == "CONTROL")
    n_total <- length(values)
    n_Control <- length(Control_index)
    n_RIF <- n_total - n_Control
    if (n_Control < 2L || n_RIF < 2L) {
        stop("Exact permutation requires at least two donors per group for ", model_id, ".")
    }
    observed_difference <- mean(values[groups == "RIF"]) - mean(values[groups == "CONTROL"])
    Control_allocations <- utils::combn(n_total, n_Control)
    permutation_differences <- apply(Control_allocations, 2, function(Control_allocation) {
        allocation_is_Control <- rep(FALSE, n_total)
        allocation_is_Control[Control_allocation] <- TRUE
        mean(values[!allocation_is_Control]) - mean(values[allocation_is_Control])
    })
    exact_P <- mean(abs(permutation_differences) >= abs(observed_difference) - tolerance)
    permutation_table <- data.frame(model_id = model_id, permutation_index = seq_along(permutation_differences), 
        permuted_mean_difference_RIF_minus_Control = as.numeric(permutation_differences), absolute_permuted_difference = abs(permutation_differences), 
        at_least_as_extreme_as_observed = abs(permutation_differences) >= abs(observed_difference) - tolerance, 
        stringsAsFactors = FALSE)
    list(observed_difference = observed_difference, exact_P = exact_P, permutation_n = length(permutation_differences), 
        minimum_resolvable_P = 1/length(permutation_differences), permutation_table = permutation_table)
}

leave_one_donor_out <- function(model_data, full_difference, model_id) {
    output_rows <- list()
    donors <- as.character(model_data$donor_id)
    for (donor_index in seq_along(donors)) {
        excluded_donor <- donors[donor_index]
        retained_data <- model_data[model_data$donor_id != excluded_donor, , drop = FALSE]
        Control_n <- sum(retained_data$analysis_group == "CONTROL")
        RIF_n <- sum(retained_data$analysis_group == "RIF")
        if (Control_n < 2L || RIF_n < 2L) {
            LOO_difference <- NA_real_
        }
        else {
            LOO_difference <- mean(retained_data$value[retained_data$analysis_group == "RIF"]) - mean(retained_data$value[retained_data$analysis_group == 
                "CONTROL"])
        }
        direction_concordant <- if (is.finite(LOO_difference) && is.finite(full_difference) && full_difference != 
            0) {
            sign(LOO_difference) == sign(full_difference)
        }
        else {
            NA
        }
        output_rows[[donor_index]] <- data.frame(model_id = model_id, excluded_donor = excluded_donor, excluded_group = model_data$analysis_group[donor_index], 
            retained_Control_n = Control_n, retained_RIF_n = RIF_n, LOO_mean_difference_RIF_minus_Control = LOO_difference, 
            direction_concordant_with_full = direction_concordant, stringsAsFactors = FALSE)
    }
    do.call(rbind, output_rows)
}

fit_donor_model <- function(model_data, broad_cell_type, endpoint, endpoint_class, matrix_role) {
    required_columns <- c("donor_id", "sample_accession", "analysis_group", "broad_cell_type", "endpoint", "value", 
        "cell_n")
    missing_columns <- setdiff(required_columns, colnames(model_data))
    if (length(missing_columns)) {
        stop("Model data are missing: ", paste(missing_columns, collapse = ", "))
    }
    model_data <- model_data[is.finite(model_data$value), , drop = FALSE]
    model_data <- model_data[order(model_data$analysis_group, model_data$donor_id), , drop = FALSE]
    if (anyDuplicated(model_data$donor_id) > 0L) {
        stop("Duplicate donor rows in ", broad_cell_type, " / ", endpoint, " / ", matrix_role, ".")
    }
    Control_values <- model_data$value[model_data$analysis_group == "CONTROL"]
    RIF_values <- model_data$value[model_data$analysis_group == "RIF"]
    n_Control <- length(Control_values)
    n_RIF <- length(RIF_values)
    model_id <- paste(matrix_role, broad_cell_type, endpoint, sep = "__")
    if (n_Control < 3L || n_RIF < 3L) {
        stop("Model is not authorized by the 3-per-group rule: ", model_id, " (", n_Control, " Control, ", n_RIF, 
            " RIF).")
    }
    constant_endpoint <- (stats::sd(model_data$value) == 0)
    permutation_result <- exact_permutation_mean_difference(values = model_data$value, groups = model_data$analysis_group, 
        model_id = model_id)
    if (constant_endpoint) {
        Welch_P <- NA_real_
        Welch_statistic <- NA_real_
        Welch_df <- NA_real_
        CI_lower <- 0
        CI_upper <- 0
    }
    else {
        Welch_result <- tryCatch(stats::t.test(x = RIF_values, y = Control_values, alternative = "two.sided", var.equal = FALSE, 
            conf.level = 0.94999999999999996), error = function(e) NULL)
        if (is.null(Welch_result)) {
            Welch_P <- NA_real_
            Welch_statistic <- NA_real_
            Welch_df <- NA_real_
            CI_lower <- NA_real_
            CI_upper <- NA_real_
        }
        else {
            Welch_P <- as.numeric(Welch_result$p.value)
            Welch_statistic <- as.numeric(Welch_result$statistic)
            Welch_df <- as.numeric(Welch_result$parameter)
            CI_lower <- as.numeric(Welch_result$conf.int[1])
            CI_upper <- as.numeric(Welch_result$conf.int[2])
        }
    }
    LOO_table <- leave_one_donor_out(model_data = model_data, full_difference = permutation_result$observed_difference, 
        model_id = model_id)
    LOO_valid <- LOO_table[is.finite(LOO_table$LOO_mean_difference_RIF_minus_Control), , drop = FALSE]
    LOO_direction_concordance <- if (nrow(LOO_valid) && permutation_result$observed_difference != 0) {
        mean(LOO_valid$direction_concordant_with_full, na.rm = TRUE)
    }
    else {
        NA_real_
    }
    model_summary <- data.frame(model_id = model_id, matrix_role = matrix_role, broad_cell_type = broad_cell_type, 
        endpoint = endpoint, endpoint_class = endpoint_class, n_Control = n_Control, n_RIF = n_RIF, total_donor_n = n_Control + 
            n_RIF, Control_mean = safe_mean(Control_values), RIF_mean = safe_mean(RIF_values), Control_median = safe_median(Control_values), 
        RIF_median = safe_median(RIF_values), Control_SD = safe_sd(Control_values), RIF_SD = safe_sd(RIF_values), 
        mean_difference_RIF_minus_Control = permutation_result$observed_difference, median_difference_RIF_minus_Control = safe_median(RIF_values) - 
            safe_median(Control_values), exact_permutation_P = permutation_result$exact_P, exact_permutation_n = permutation_result$permutation_n, 
        minimum_resolvable_exact_P = permutation_result$minimum_resolvable_P, Welch_P = Welch_P, Welch_statistic = Welch_statistic, 
        Welch_df = Welch_df, Welch_95CI_lower = CI_lower, Welch_95CI_upper = CI_upper, Hedges_g_RIF_minus_Control = hedges_g(Control_values, 
            RIF_values), rank_biserial_RIF_minus_Control = rank_biserial_RIF_minus_Control(Control_values, RIF_values), 
        constant_endpoint = constant_endpoint, LOO_valid_n = nrow(LOO_valid), LOO_direction_concordance = LOO_direction_concordance, 
        LOO_minimum_difference = safe_min(LOO_valid$LOO_mean_difference_RIF_minus_Control), LOO_maximum_difference = safe_max(LOO_valid$LOO_mean_difference_RIF_minus_Control), 
        effect_direction = if (permutation_result$observed_difference > 0) {
            "RIF_HIGHER"
        }
        else if (permutation_result$observed_difference < 0) {
            "RIF_LOWER"
        }
        else {
            "NO_DIFFERENCE"
        }, biological_interpretation = if (endpoint_class == "VDR_EXPRESSION") {
            "VDR_LOG2_CPM_PLUS_1"
        }
        else {
            "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION"
        }, stringsAsFactors = FALSE)
    input_table <- model_data
    input_table$model_id <- model_id
    input_table$matrix_role <- matrix_role
    list(summary = model_summary, input = input_table, permutation = permutation_result$permutation_table, LOO = LOO_table)
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

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

integration_root <- file.path(results_root, "S4E1_sc", integration_version)

pb_root <- file.path(integration_root, "pb", "G183", pseudobulk_version)

primary_VDR_path <- file.path(pb_root, "tab", "10_primary_VDR_expression.csv")

atlas_VDR_path <- file.path(pb_root, "tab", "11_atlas_VDR_expression.csv")

primary_program_path <- file.path(pb_root, "tab", "12_primary_program_scores.csv")

atlas_program_path <- file.path(pb_root, "tab", "13_atlas_program_scores.csv")

model_eligibility_path <- file.path(pb_root, "tab", "14_model_eligibility.csv")

cluster8_program_sensitivity_path <- file.path(pb_root, "tab", "15_cluster8_program_sensitivity.csv")

cluster8_VDR_sensitivity_path <- file.path(pb_root, "tab", "16_cluster8_VDR_sensitivity.csv")

output_root <- file.path(integration_root, "mdl", "G183", model_version)

table_root <- file.path(output_root, "tab")

source_data_root <- file.path(output_root, "src")

figure_output_root <- file.path(output_root, "fig")

manifest_root <- file.path(output_root, "man")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1O_", format(run_start, "%Y%m%d_%H%M%S"), "_models"))

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

planned_paths <- data.frame(artifact = c("primary_VDR", "atlas_VDR", "primary_program", "atlas_program", "model_eligibility", 
    "output_root", "run_qc_root", "handoff_zip"), path = normalize_slashes(c(primary_VDR_path, atlas_VDR_path, 
    primary_program_path, atlas_program_path, model_eligibility_path, output_root, run_qc_root, handoff_zip)), 
    stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

if (!all(planned_paths$path_safe)) {
    stop("One or more donor-model paths exceed the safe limit.")
}

pb_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

pb_qc_dirs <- pb_qc_dirs[grepl("^STEP04E1N_[0-9]{8}_[0-9]{6}_pbQC$", basename(pb_qc_dirs), perl = TRUE)]

if (!length(pb_qc_dirs)) {
    stop("No STEP04E1N pseudobulk-QC directory was found.")
}

gate_candidates <- list()

for (directory_path in pb_qc_dirs) {
    gate_path <- file.path(directory_path, "G183_PBVDR_GATE.csv")
    gate_value <- safe_read_csv(gate_path)
    if (!is.null(gate_value) && nrow(gate_value) == 1L) {
        gate_value$source_directory <- normalize_slashes(directory_path)
        gate_value$source_mtime <- as.character(file.info(directory_path)$mtime)
        gate_candidates[[length(gate_candidates) + 1L]] <- gate_value
    }
}

if (!length(gate_candidates)) {
    stop("No readable STEP04E1N Gate was found.")
}

gate_candidates <- as.data.frame(data.table::rbindlist(gate_candidates, use.names = TRUE, fill = TRUE))

gate_candidates$source_mtime_POSIX <- as.POSIXct(gate_candidates$source_mtime)

gate_candidates <- gate_candidates[order(gate_candidates$source_mtime_POSIX, decreasing = TRUE), , drop = FALSE]

authorized_index <- which(gate_candidates$gate_status == "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED" & 
    toupper(as.character(gate_candidates$donor_level_VDR_program_models_authorized)) %in% c("TRUE", "T", "1", "YES"))

if (!length(authorized_index)) {
    stop("No pseudobulk Gate authorizes donor-level models.")
}

upstream_gate <- gate_candidates[authorized_index[1], , drop = FALSE]

write_csv_utf8(upstream_gate, file.path(run_qc_root, "00_upstream_PBVDR_gate.csv"))

required_inputs <- c(primary_VDR_path, atlas_VDR_path, primary_program_path, atlas_program_path, model_eligibility_path, 
    cluster8_program_sensitivity_path, cluster8_VDR_sensitivity_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(normalize_slashes(missing_inputs), collapse = "\n"))
}

primary_VDR <- safe_read_csv(primary_VDR_path)

atlas_VDR <- safe_read_csv(atlas_VDR_path)

primary_program <- safe_read_csv(primary_program_path)

atlas_program <- safe_read_csv(atlas_program_path)

model_eligibility <- safe_read_csv(model_eligibility_path)

cluster8_program_sensitivity <- safe_read_csv(cluster8_program_sensitivity_path)

cluster8_VDR_sensitivity <- safe_read_csv(cluster8_VDR_sensitivity_path)

for (data_name in c("primary_VDR", "atlas_VDR", "primary_program", "atlas_program", "model_eligibility")) {
    data_value <- get(data_name)
    if (is.null(data_value) || !nrow(data_value) || anyDuplicated(colnames(data_value)) > 0L) {
        stop(data_name, " is absent, empty or has duplicate columns.")
    }
}

for (column_name in c("minimum_cell_threshold_pass", "donor_level_model_eligible")) {
    if (column_name %in% colnames(primary_VDR)) {
        primary_VDR[[column_name]] <- as_logical_safe(primary_VDR[[column_name]])
    }
    if (column_name %in% colnames(atlas_VDR)) {
        atlas_VDR[[column_name]] <- as_logical_safe(atlas_VDR[[column_name]])
    }
    if (column_name %in% colnames(primary_program)) {
        primary_program[[column_name]] <- as_logical_safe(primary_program[[column_name]])
    }
    if (column_name %in% colnames(atlas_program)) {
        atlas_program[[column_name]] <- as_logical_safe(atlas_program[[column_name]])
    }
    if (column_name %in% colnames(model_eligibility)) {
        model_eligibility[[column_name]] <- as_logical_safe(model_eligibility[[column_name]])
    }
}

eligible_types <- model_eligibility$broad_cell_type[model_eligibility$donor_level_model_eligible]

eligible_types <- sort(as.character(eligible_types))

expected_eligible_types <- sort(c("Endothelial cells", "Myeloid cells", "NK/T cells", "Stromal cells", "Unciliated epithelial cells"))

if (!identical(eligible_types, expected_eligible_types)) {
    stop("Donor-model-eligible broad cell types do not match the locked five-type set.")
}

program_names <- c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", "PROMOTER_OCCUPANCY_84")

if (!all(program_names %in% unique(primary_program$program_name)) || !all(program_names %in% unique(atlas_program$program_name))) {
    stop("One or more locked VDR programs are absent.")
}

log_event("Starting", run_id, "eligible_types=", collapse_unique(eligible_types))

primary_VDR_model_data <- primary_VDR[primary_VDR$broad_cell_type %in% eligible_types & primary_VDR$minimum_cell_threshold_pass, 
    , drop = FALSE]

primary_VDR_model_data$endpoint <- "VDR_expression"

primary_VDR_model_data$endpoint_class <- "VDR_EXPRESSION"

primary_VDR_model_data$value <- as.numeric(primary_VDR_model_data$VDR_log2_CPM_plus_1)

primary_program_model_data <- primary_program[primary_program$broad_cell_type %in% eligible_types & primary_program$minimum_cell_threshold_pass & 
    primary_program$program_name %in% program_names, , drop = FALSE]

primary_program_model_data$endpoint <- as.character(primary_program_model_data$program_name)

primary_program_model_data$endpoint_class <- "VDR_PROGRAM"

primary_program_model_data$value <- as.numeric(primary_program_model_data$score)

model_columns <- c("donor_id", "sample_accession", "analysis_group", "broad_cell_type", "endpoint", "endpoint_class", 
    "value", "cell_n")

primary_model_data <- rbind(primary_VDR_model_data[, model_columns, drop = FALSE], primary_program_model_data[, 
    model_columns, drop = FALSE])

primary_model_data$donor_id <- as.character(primary_model_data$donor_id)

primary_model_data$sample_accession <- as.character(primary_model_data$sample_accession)

primary_model_data$analysis_group <- as.character(primary_model_data$analysis_group)

primary_model_data$broad_cell_type <- as.character(primary_model_data$broad_cell_type)

primary_model_data$endpoint <- as.character(primary_model_data$endpoint)

if (any(!primary_model_data$analysis_group %in% c("CONTROL", "RIF")) || any(!is.finite(primary_model_data$value))) {
    stop("Primary model data contain invalid groups or values.")
}

write_csv_utf8(primary_model_data, file.path(source_data_root, "SD1_primary_model_data.csv"))

endpoint_order <- c("VDR_expression", program_names)

primary_summary_rows <- list()

primary_input_rows <- list()

primary_permutation_rows <- list()

primary_LOO_rows <- list()

model_index <- 0L

for (cell_type in eligible_types) {
    for (endpoint_value in endpoint_order) {
        model_index <- model_index + 1L
        model_data <- primary_model_data[primary_model_data$broad_cell_type == cell_type & primary_model_data$endpoint == 
            endpoint_value, , drop = FALSE]
        endpoint_class_value <- unique(model_data$endpoint_class)
        if (length(endpoint_class_value) != 1L) {
            stop("Endpoint class is not unique for ", cell_type, " / ", endpoint_value, ".")
        }
        log_event("MODEL", model_index, "/25", cell_type, endpoint_value, "donors=", nrow(model_data))
        model_result <- fit_donor_model(model_data = model_data, broad_cell_type = cell_type, endpoint = endpoint_value, 
            endpoint_class = endpoint_class_value, matrix_role = "PRIMARY")
        primary_summary_rows[[model_index]] <- model_result$summary
        primary_input_rows[[model_index]] <- model_result$input
        primary_permutation_rows[[model_index]] <- model_result$permutation
        primary_LOO_rows[[model_index]] <- model_result$LOO
    }
}

primary_models <- do.call(rbind, primary_summary_rows)

raw_primary_models_checkpoint <- file.path(source_data_root, "SD0_primary_models_pre_multiplicity.rds")

saveRDS(primary_models, raw_primary_models_checkpoint, compress = "gzip")

primary_model_inputs <- do.call(rbind, primary_input_rows)

primary_permutations <- do.call(rbind, primary_permutation_rows)

primary_LOO <- do.call(rbind, primary_LOO_rows)

if (nrow(primary_models) != 25L || length(unique(primary_models$broad_cell_type)) != 5L || length(unique(primary_models$endpoint)) != 
    5L || anyDuplicated(primary_models$model_id) > 0L) {
    stop("Primary model grid is not the planned 5 x 5 design.")
}

primary_models$BH_global_25 <- stats::p.adjust(primary_models$exact_permutation_P, method = "BH")

primary_models$BH_within_endpoint <- ave(primary_models$exact_permutation_P, primary_models$endpoint, FUN = function(x) {
    stats::p.adjust(x, method = "BH")
})

primary_models$BH_within_cell_type <- ave(primary_models$exact_permutation_P, primary_models$broad_cell_type, FUN = function(x) {
    stats::p.adjust(x, method = "BH")
})

primary_models$statistical_status <- ifelse(primary_models$constant_endpoint, "CONSTANT_UNTESTABLE_ENDPOINT", ifelse(primary_models$BH_global_25 < 
    0.050000000000000003, "GLOBAL_BH_SIGNIFICANT", ifelse(primary_models$BH_within_endpoint < 0.050000000000000003, 
    "ENDPOINT_FAMILY_BH_SIGNIFICANT", ifelse(primary_models$exact_permutation_P < 0.050000000000000003, "NOMINAL_EXACT_P_LT_0_05_EXPLORATORY", 
        "NO_NOMINAL_EVIDENCE"))))

primary_models$LOO_review_flag <- ifelse(is.finite(primary_models$LOO_direction_concordance) & primary_models$LOO_direction_concordance < 
    0.75, "LOO_DIRECTION_UNSTABLE", "")

primary_models$resolution_review_flag <- ifelse(primary_models$exact_permutation_n < 84, paste0("LIMITED_EXACT_RESOLUTION_MIN_P_", 
    format(primary_models$minimum_resolvable_exact_P, digits = 4)), "")

primary_models$interpretation_restriction <- ifelse(primary_models$endpoint_class == "VDR_PROGRAM", "RELATIVE_ALIGNMENT_SCORE_NOT_DIRECT_VDR_ACTIVITY", 
    "VDR_EXPRESSION_ONLY_NOT_CAUSAL_ACTIVITY")

primary_models <- primary_models[order(match(primary_models$broad_cell_type, eligible_types), match(primary_models$endpoint, 
    endpoint_order)), , drop = FALSE]

write_csv_utf8(primary_models, file.path(table_root, "01_primary_models.csv"))

write_csv_utf8(primary_model_inputs, file.path(source_data_root, "SD2_primary_model_inputs.csv"))

write_csv_utf8(primary_permutations, file.path(source_data_root, "SD3_exact_permutation_distributions.csv"))

write_csv_utf8(primary_LOO, file.path(source_data_root, "SD4_leave_one_donor_out.csv"))

group_summary <- data.table::as.data.table(primary_model_data)[, .(donor_n = as.integer(data.table::uniqueN(donor_id)), 
    mean_value = as.numeric(mean(value)), median_value = as.numeric(median(value)), SD_value = as.numeric(stats::sd(value)), 
    minimum_value = as.numeric(min(value)), maximum_value = as.numeric(max(value)), median_cell_n = as.numeric(median(cell_n)), 
    minimum_cell_n = as.numeric(min(cell_n)), maximum_cell_n = as.numeric(max(cell_n))), by = .(broad_cell_type, 
    endpoint, endpoint_class, analysis_group)]

group_summary <- as.data.frame(group_summary)

group_summary_type_audit <- data.frame(column_name = colnames(group_summary), column_class = vapply(group_summary, 
    function(x) {
        paste(class(x), collapse = ";")
    }, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(group_summary_type_audit, file.path(table_root, "00_group_summary_type_audit.csv"))

write_csv_utf8(group_summary, file.path(table_root, "02_group_summary.csv"))

multiplicity_summary <- data.frame(metric = c("planned_primary_model_n", "constant_endpoint_n", "nominal_exact_P_lt_0_05_n", 
    "global_BH_lt_0_05_n", "endpoint_family_BH_lt_0_05_n", "cell_type_family_BH_lt_0_05_n", "LOO_direction_unstable_n", 
    "unciliated_models_with_35_permutations_n"), value = c(nrow(primary_models), sum(primary_models$constant_endpoint), 
    sum(primary_models$exact_permutation_P < 0.050000000000000003), sum(primary_models$BH_global_25 < 0.050000000000000003), 
    sum(primary_models$BH_within_endpoint < 0.050000000000000003), sum(primary_models$BH_within_cell_type < 0.050000000000000003), 
    sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE"), sum(primary_models$exact_permutation_n == 
        35L)), stringsAsFactors = FALSE)

write_csv_utf8(multiplicity_summary, file.path(table_root, "03_multiplicity_summary.csv"))

review_models <- primary_models[primary_models$statistical_status != "NO_NOMINAL_EVIDENCE" | nzchar(primary_models$LOO_review_flag) | 
    nzchar(primary_models$resolution_review_flag), , drop = FALSE]

write_csv_utf8(review_models, file.path(table_root, "04_models_requiring_review.csv"))

atlas_VDR_model_data <- atlas_VDR[atlas_VDR$broad_cell_type == "Unciliated epithelial cells" & atlas_VDR$minimum_cell_threshold_pass, 
    , drop = FALSE]

atlas_VDR_model_data$endpoint <- "VDR_expression"

atlas_VDR_model_data$endpoint_class <- "VDR_EXPRESSION"

atlas_VDR_model_data$value <- as.numeric(atlas_VDR_model_data$VDR_log2_CPM_plus_1)

atlas_program_model_data <- atlas_program[atlas_program$broad_cell_type == "Unciliated epithelial cells" & atlas_program$minimum_cell_threshold_pass & 
    atlas_program$program_name %in% program_names, , drop = FALSE]

atlas_program_model_data$endpoint <- as.character(atlas_program_model_data$program_name)

atlas_program_model_data$endpoint_class <- "VDR_PROGRAM"

atlas_program_model_data$value <- as.numeric(atlas_program_model_data$score)

atlas_unciliated_data <- rbind(atlas_VDR_model_data[, model_columns, drop = FALSE], atlas_program_model_data[, 
    model_columns, drop = FALSE])

sensitivity_summary_rows <- list()

sensitivity_input_rows <- list()

sensitivity_permutation_rows <- list()

sensitivity_LOO_rows <- list()

for (endpoint_index in seq_along(endpoint_order)) {
    endpoint_value <- endpoint_order[endpoint_index]
    model_data <- atlas_unciliated_data[atlas_unciliated_data$endpoint == endpoint_value, , drop = FALSE]
    endpoint_class_value <- unique(model_data$endpoint_class)
    sensitivity_result <- fit_donor_model(model_data = model_data, broad_cell_type = "Unciliated epithelial cells", 
        endpoint = endpoint_value, endpoint_class = endpoint_class_value, matrix_role = "ATLAS_SENSITIVITY")
    sensitivity_summary_rows[[endpoint_index]] <- sensitivity_result$summary
    sensitivity_input_rows[[endpoint_index]] <- sensitivity_result$input
    sensitivity_permutation_rows[[endpoint_index]] <- sensitivity_result$permutation
    sensitivity_LOO_rows[[endpoint_index]] <- sensitivity_result$LOO
}

atlas_sensitivity_models <- do.call(rbind, sensitivity_summary_rows)

atlas_sensitivity_inputs <- do.call(rbind, sensitivity_input_rows)

atlas_sensitivity_permutations <- do.call(rbind, sensitivity_permutation_rows)

atlas_sensitivity_LOO <- do.call(rbind, sensitivity_LOO_rows)

atlas_sensitivity_models$BH_sensitivity_5 <- stats::p.adjust(atlas_sensitivity_models$exact_permutation_P, method = "BH")

primary_unciliated_models <- primary_models[primary_models$broad_cell_type == "Unciliated epithelial cells", , 
    drop = FALSE]

cluster8_model_sensitivity <- merge(primary_unciliated_models[, c("endpoint", "n_Control", "n_RIF", "mean_difference_RIF_minus_Control", 
    "exact_permutation_P", "minimum_resolvable_exact_P", "Hedges_g_RIF_minus_Control", "LOO_direction_concordance"), 
    drop = FALSE], atlas_sensitivity_models[, c("endpoint", "n_Control", "n_RIF", "mean_difference_RIF_minus_Control", 
    "exact_permutation_P", "minimum_resolvable_exact_P", "Hedges_g_RIF_minus_Control", "LOO_direction_concordance", 
    "BH_sensitivity_5"), drop = FALSE], by = "endpoint", all = TRUE, suffixes = c("_primary", "_atlas"), sort = FALSE)

cluster8_model_sensitivity$effect_delta_atlas_minus_primary <- (cluster8_model_sensitivity$mean_difference_RIF_minus_Control_atlas - 
    cluster8_model_sensitivity$mean_difference_RIF_minus_Control_primary)

cluster8_model_sensitivity$direction_concordant <- ifelse(cluster8_model_sensitivity$mean_difference_RIF_minus_Control_primary == 
    0 | cluster8_model_sensitivity$mean_difference_RIF_minus_Control_atlas == 0, NA, sign(cluster8_model_sensitivity$mean_difference_RIF_minus_Control_primary) == 
    sign(cluster8_model_sensitivity$mean_difference_RIF_minus_Control_atlas))

cluster8_model_sensitivity$sensitivity_interpretation <- ifelse(cluster8_model_sensitivity$direction_concordant %in% 
    TRUE, "DIRECTION_CONCORDANT", ifelse(cluster8_model_sensitivity$direction_concordant %in% FALSE, "DIRECTION_CHANGED_WITH_CLUSTER8", 
    "ZERO_EFFECT_DIRECTION_NOT_CLASSIFIED"))

write_csv_utf8(atlas_sensitivity_models, file.path(table_root, "05_atlas_unciliated_models.csv"))

write_csv_utf8(cluster8_model_sensitivity, file.path(table_root, "06_cluster8_model_sensitivity.csv"))

write_csv_utf8(atlas_sensitivity_inputs, file.path(source_data_root, "SD5_atlas_sensitivity_inputs.csv"))

write_csv_utf8(atlas_sensitivity_permutations, file.path(source_data_root, "SD6_atlas_sensitivity_permutations.csv"))

write_csv_utf8(atlas_sensitivity_LOO, file.path(source_data_root, "SD7_atlas_sensitivity_LOO.csv"))

ciliated_VDR <- primary_VDR[primary_VDR$broad_cell_type == "Ciliated epithelial cells", , drop = FALSE]

ciliated_VDR$endpoint <- "VDR_expression"

ciliated_VDR$value <- ciliated_VDR$VDR_log2_CPM_plus_1

ciliated_program <- primary_program[primary_program$broad_cell_type == "Ciliated epithelial cells" & primary_program$program_name %in% 
    program_names, , drop = FALSE]

ciliated_program$endpoint <- ciliated_program$program_name

ciliated_program$value <- ciliated_program$score

ciliated_columns <- c("donor_id", "sample_accession", "analysis_group", "broad_cell_type", "endpoint", "value", 
    "cell_n", "minimum_cell_threshold_pass")

ciliated_data <- rbind(ciliated_VDR[, ciliated_columns, drop = FALSE], ciliated_program[, ciliated_columns, drop = FALSE])

ciliated_descriptive <- data.table::as.data.table(ciliated_data)[, .(donor_profile_n = as.integer(.N), threshold_pass_profile_n = as.integer(sum(minimum_cell_threshold_pass)), 
    total_cell_n = as.numeric(sum(cell_n)), median_cell_n = as.numeric(median(cell_n)), mean_value = as.numeric(mean(value)), 
    median_value = as.numeric(median(value)), minimum_value = as.numeric(min(value)), maximum_value = as.numeric(max(value))), 
    by = .(endpoint, analysis_group)]

ciliated_descriptive <- as.data.frame(ciliated_descriptive)

write_csv_utf8(ciliated_data, file.path(source_data_root, "SD8_ciliated_descriptive_data.csv"))

write_csv_utf8(ciliated_descriptive, file.path(table_root, "07_ciliated_descriptive_summary.csv"))

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

forest_data <- primary_models

forest_data$broad_cell_type <- factor(forest_data$broad_cell_type, levels = rev(eligible_types))

forest_data$endpoint <- factor(forest_data$endpoint, levels = endpoint_order)

plot_forest <- ggplot2::ggplot(forest_data, ggplot2::aes(x = mean_difference_RIF_minus_Control, y = broad_cell_type, 
    shape = statistical_status)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = Welch_95CI_lower, 
    xmax = Welch_95CI_upper), height = 0.14999999999999999, na.rm = TRUE) + ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~endpoint, 
    scales = "free_x", ncol = 2) + ggplot2::labs(title = "GSE183837 donor-level Control-versus-RIF effects", subtitle = "Effect is RIF minus Control; intervals are secondary Welch 95% confidence intervals", 
    x = "Mean difference", y = "Broad cell type", shape = "Statistical status") + theme_QC

save_plot_pair(plot_forest, "QC1_primary_effect_forest", width = 11, height = 9)

pvalue_data <- primary_models

pvalue_data$minus_log10_exact_P <- -log10(pmax(pvalue_data$exact_permutation_P, 9.9999999999999998e-13))

pvalue_data$minus_log10_global_BH <- -log10(pmax(pvalue_data$BH_global_25, 9.9999999999999998e-13))

plot_pvalues <- ggplot2::ggplot(pvalue_data, ggplot2::aes(x = minus_log10_exact_P, y = minus_log10_global_BH, shape = endpoint)) + 
    ggplot2::geom_vline(xintercept = -log10(0.050000000000000003), linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(0.050000000000000003), 
    linetype = "dashed") + ggplot2::geom_point(size = 2.5, alpha = 0.80000000000000004) + ggplot2::labs(title = "Exact permutation evidence before and after global BH correction", 
    subtitle = "Global multiplicity includes all 25 planned donor-level models", x = "-log10 exact permutation P", 
    y = "-log10 global BH q", shape = "Endpoint") + theme_QC

save_plot_pair(plot_pvalues, "QC2_exact_P_and_global_BH", width = 8, height = 6)

VDR_plot_data <- primary_model_data[primary_model_data$endpoint == "VDR_expression", , drop = FALSE]

VDR_plot_data$broad_cell_type <- factor(VDR_plot_data$broad_cell_type, levels = eligible_types)

plot_VDR_groups <- ggplot2::ggplot(VDR_plot_data, ggplot2::aes(x = analysis_group, y = value, shape = analysis_group)) + 
    ggplot2::geom_point(size = 2.5, alpha = 0.84999999999999998, position = ggplot2::position_jitter(width = 0.080000000000000002, 
        height = 0)) + ggplot2::facet_wrap(~broad_cell_type, scales = "free_y", ncol = 3) + ggplot2::labs(title = "Donor-level VDR expression by group", 
    subtitle = "Each point is one threshold-passing donor pseudobulk", x = "Group", y = "VDR log2(CPM + 1)", shape = "Group") + 
    theme_QC

save_plot_pair(plot_VDR_groups, "QC3_VDR_group_profiles", width = 10, height = 7)

program_plot_data <- primary_model_data[primary_model_data$endpoint_class == "VDR_PROGRAM", , drop = FALSE]

program_plot_data$broad_cell_type <- factor(program_plot_data$broad_cell_type, levels = eligible_types)

program_plot_data$endpoint <- factor(program_plot_data$endpoint, levels = program_names)

plot_program_groups <- ggplot2::ggplot(program_plot_data, ggplot2::aes(x = analysis_group, y = value, shape = analysis_group)) + 
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") + ggplot2::geom_point(size = 1.8999999999999999, alpha = 0.80000000000000004, 
    position = ggplot2::position_jitter(width = 0.070000000000000007, height = 0)) + ggplot2::facet_grid(broad_cell_type ~ 
    endpoint, scales = "free_y") + ggplot2::labs(title = "Locked VDR-program scores by donor group", subtitle = "Positive scores indicate greater relative alignment with the locked interaction direction", 
    x = "Group", y = "Signed program score", shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_program_groups, "QC4_program_group_profiles", width = 13, height = 11)

LOO_plot_data <- primary_models

LOO_plot_data$broad_cell_type <- factor(LOO_plot_data$broad_cell_type, levels = rev(eligible_types))

LOO_plot_data$endpoint <- factor(LOO_plot_data$endpoint, levels = endpoint_order)

plot_LOO <- ggplot2::ggplot(LOO_plot_data, ggplot2::aes(x = mean_difference_RIF_minus_Control, y = broad_cell_type, 
    shape = LOO_review_flag)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = LOO_minimum_difference, 
    xmax = LOO_maximum_difference), height = 0.14999999999999999) + ggplot2::geom_point(size = 2.3999999999999999) + 
    ggplot2::facet_wrap(~endpoint, scales = "free_x", ncol = 2) + ggplot2::labs(title = "Leave-one-donor-out effect stability", 
    subtitle = "Intervals show the minimum and maximum RIF-minus-Control effect after excluding each donor", x = "Mean difference", 
    y = "Broad cell type", shape = "LOO review") + theme_QC

save_plot_pair(plot_LOO, "QC5_leave_one_donor_out", width = 11, height = 9)

sensitivity_plot_data <- cluster8_model_sensitivity

sensitivity_plot_long <- rbind(data.frame(endpoint = sensitivity_plot_data$endpoint, matrix_role = "PRIMARY", mean_difference = sensitivity_plot_data$mean_difference_RIF_minus_Control_primary, 
    exact_P = sensitivity_plot_data$exact_permutation_P_primary, stringsAsFactors = FALSE), data.frame(endpoint = sensitivity_plot_data$endpoint, 
    matrix_role = "ATLAS_SENSITIVITY", mean_difference = sensitivity_plot_data$mean_difference_RIF_minus_Control_atlas, 
    exact_P = sensitivity_plot_data$exact_permutation_P_atlas, stringsAsFactors = FALSE))

sensitivity_plot_long$endpoint <- factor(sensitivity_plot_long$endpoint, levels = endpoint_order)

plot_sensitivity <- ggplot2::ggplot(sensitivity_plot_long, ggplot2::aes(x = matrix_role, y = mean_difference, group = endpoint, 
    shape = endpoint)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + ggplot2::geom_line(alpha = 0.59999999999999998) + 
    ggplot2::geom_point(size = 2.7000000000000002) + ggplot2::labs(title = "Cluster-8 sensitivity of unciliated epithelial effects", 
    subtitle = "Primary excludes cluster 8; atlas sensitivity retains it", x = "Analysis matrix", y = "RIF-minus-Control mean difference", 
    shape = "Endpoint") + theme_QC

save_plot_pair(plot_sensitivity, "QC6_cluster8_model_sensitivity", width = 9, height = 6)

ciliated_plot_data <- ciliated_data

ciliated_plot_data$endpoint <- factor(ciliated_plot_data$endpoint, levels = endpoint_order)

plot_ciliated <- ggplot2::ggplot(ciliated_plot_data, ggplot2::aes(x = analysis_group, y = value, shape = minimum_cell_threshold_pass)) + 
    ggplot2::geom_point(size = 2.3999999999999999, alpha = 0.80000000000000004, position = ggplot2::position_jitter(width = 0.080000000000000002, 
        height = 0)) + ggplot2::facet_wrap(~endpoint, scales = "free_y", ncol = 2) + ggplot2::labs(title = "Ciliated epithelial descriptive profiles", 
    subtitle = "Only two donor profiles pass the 20-cell threshold; no group test is performed", x = "Group", y = "Endpoint value", 
    shape = "20-cell threshold pass") + theme_QC

save_plot_pair(plot_ciliated, "QC7_ciliated_descriptive_only", width = 10, height = 8)

figure_files <- list.files(figure_root, full.names = TRUE)

file.copy(figure_files, figure_output_root, overwrite = TRUE)

nominal_models <- primary_models[primary_models$exact_permutation_P < 0.050000000000000003, , drop = FALSE]

global_BH_models <- primary_models[primary_models$BH_global_25 < 0.050000000000000003, , drop = FALSE]

constant_models <- primary_models[primary_models$constant_endpoint, , drop = FALSE]

direction_changed_sensitivity <- cluster8_model_sensitivity[cluster8_model_sensitivity$direction_concordant %in% 
    FALSE, , drop = FALSE]

gate_checks <- data.frame(check = c("upstream_gate_authorizes_donor_models", "five_locked_broad_cell_types_present", 
    "five_locked_endpoints_present", "planned_primary_model_n_is_25", "all_primary_models_have_at_least_3_donors_per_group", 
    "four_full_design_types_have_84_exact_allocations", "unciliated_models_have_35_exact_allocations", "endothelial_VDR_is_retained_as_constant_endpoint", 
    "global_BH_values_complete", "within_endpoint_BH_values_complete", "within_cell_type_BH_values_complete", "LOO_results_exist_for_all_primary_models", 
    "atlas_unciliated_sensitivity_has_5_models", "atlas_unciliated_uses_3_Control_6_RIF", "cluster8_primary_atlas_comparison_has_5_endpoints", 
    "ciliated_epithelial_has_no_inferential_model", "cell_level_testing_not_performed", "gene_level_DE_not_performed", 
    "independent_validation_review_can_be_authorized"), pass = c(upstream_gate$gate_status[[1]] == "PASS_GSE183837_PSEUDOBULK_VDR_LOCALIZATION_DONOR_MODELS_AUTHORIZED", 
    identical(sort(unique(primary_models$broad_cell_type)), expected_eligible_types), identical(unique(primary_models$endpoint[match(endpoint_order, 
        primary_models$endpoint)]), endpoint_order), nrow(primary_models) == 25L, all(primary_models$n_Control >= 
        3L & primary_models$n_RIF >= 3L), all(primary_models$exact_permutation_n[primary_models$broad_cell_type != 
        "Unciliated epithelial cells"] == 84L), all(primary_models$exact_permutation_n[primary_models$broad_cell_type == 
        "Unciliated epithelial cells"] == 35L), nrow(constant_models) == 1L && constant_models$broad_cell_type == 
        "Endothelial cells" && constant_models$endpoint == "VDR_expression", all(is.finite(primary_models$BH_global_25)), 
    all(is.finite(primary_models$BH_within_endpoint)), all(is.finite(primary_models$BH_within_cell_type)), length(unique(primary_LOO$model_id)) == 
        25L, nrow(atlas_sensitivity_models) == 5L, all(atlas_sensitivity_models$n_Control == 3L & atlas_sensitivity_models$n_RIF == 
        6L), nrow(cluster8_model_sensitivity) == 5L, !("Ciliated epithelial cells" %in% primary_models$broad_cell_type), 
    TRUE, TRUE, TRUE), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_GSE183837_DONOR_MODELS_COMPLETE_INDEPENDENT_VALIDATION_REVIEW_AUTHORIZED"
} else {
    "HOLD_GSE183837_DONOR_MODELS_INCOMPLETE"
}

review_flags <- c("CONTROL_N3_RIF_N6_LOW_POWER", "UNCILIATED_PRIMARY_RIF_N4_EXACT_MIN_P_1_OVER_35", "ENDOTHELIAL_VDR_CONSTANT_ZERO", 
    "CILIATED_EPITHELIAL_DESCRIPTIVE_ONLY", "CLUSTER8_ATLAS_SENSITIVITY_ONLY", if (nrow(nominal_models)) {
        paste0("NOMINAL_EXACT_P_LT_0_05_MODELS_", nrow(nominal_models))
    } else {
        "NO_NOMINAL_EXACT_P_LT_0_05"
    }, if (nrow(global_BH_models)) {
        paste0("GLOBAL_BH_SIGNIFICANT_MODELS_", nrow(global_BH_models))
    } else {
        "NO_GLOBAL_BH_SIGNIFICANT_MODELS"
    }, if (nrow(direction_changed_sensitivity)) {
        paste0("CLUSTER8_DIRECTION_CHANGED_ENDPOINTS_", collapse_unique(direction_changed_sensitivity$endpoint))
    } else {
        "CLUSTER8_DIRECTIONS_CONCORDANT_OR_ZERO"
    })

model_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, pseudobulk_version = pseudobulk_version, 
    model_version = model_version, gate_status = gate_status, primary_model_n = nrow(primary_models), testable_nonconstant_model_n = sum(!primary_models$constant_endpoint), 
    nominal_exact_P_lt_0_05_n = nrow(nominal_models), global_BH_lt_0_05_n = nrow(global_BH_models), endpoint_family_BH_lt_0_05_n = sum(primary_models$BH_within_endpoint < 
        0.050000000000000003), LOO_direction_unstable_n = sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE"), 
    cluster8_direction_changed_n = nrow(direction_changed_sensitivity), nominal_models = collapse_unique(paste(nominal_models$broad_cell_type, 
        nominal_models$endpoint, sep = " | ")), global_BH_models = collapse_unique(paste(global_BH_models$broad_cell_type, 
        global_BH_models$endpoint, sep = " | ")), review_flags = collapse_unique(review_flags), independent_validation_review_authorized = all_gate_checks_pass, 
    cross_dataset_concordance_review_authorized = all_gate_checks_pass, gene_level_pseudobulk_DE_authorized = FALSE, 
    result_lock_authorized = FALSE, stringsAsFactors = FALSE)

write_csv_utf8(model_gate, file.path(run_qc_root, "G183_DONOR_MODEL_GATE.csv"))

card_lines <- c("# STEP04E1O GSE183837 Donor Model Card", "", paste0("- Gate: ", gate_status), "- Primary models: 25", 
    paste0("- Nonconstant testable models: ", sum(!primary_models$constant_endpoint)), paste0("- Nominal exact P < 0.05: ", 
        nrow(nominal_models)), paste0("- Global BH q < 0.05: ", nrow(global_BH_models)), paste0("- Endpoint-family BH q < 0.05: ", 
        sum(primary_models$BH_within_endpoint < 0.050000000000000003)), paste0("- LOO direction-unstable models: ", 
        sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE")), paste0("- Cluster-8 sensitivity direction changes: ", 
        nrow(direction_changed_sensitivity)), paste0("- Independent validation review authorized: ", all_gate_checks_pass), 
    "- Gene-level pseudobulk DE authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "G183_DONOR_MODEL_CARD.md"))

interpretation_lines <- c("# GSE183837 donor-model interpretation", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: **", 
    gate_status, "**"), "", "## Prespecified design", "", "- Twenty-five donor-level models were tested across five broad cell types and five endpoints.", 
    "- Exact permutation tests were primary; Welch tests and confidence intervals were secondary.", "- Global BH correction covered all 25 planned models.", 
    "", "## Interpretation rules", "", "- Nominal exact P < 0.05 without global BH significance is exploratory only.", 
    "- A negative VDR-program effect means lower relative alignment in RIF with the locked calcitriol-dependent interaction direction.", 
    "- Program effects do not establish direct binding, VDR activity or causality.", "- Ciliated epithelial cells remain descriptive only.", 
    "- Cluster 8 is evaluated only through atlas sensitivity analysis.", "", "## Authorized next operation", "", 
    "- Compare the independent GSE183837 effect directions with the locked GSE250130 donor-level results.", "- Determine convergent, discordant and non-informative evidence.", 
    "- Lock the cross-dataset interpretation only after consistency QC.", "", "## Not authorized", "", "- No gene-level pseudobulk DE.", 
    "- No cell-level disease testing.", "- No disease mechanism claim.", "- No final result lock.")

write_text(interpretation_lines, file.path(run_qc_root, "G183_DONOR_MODEL_INTERPRETATION.md"))

nominal_description <- if (nrow(nominal_models)) {
    paste(paste0(nominal_models$broad_cell_type, " / ", nominal_models$endpoint, " (exact P=", format(nominal_models$exact_permutation_P, 
        digits = 4), ", global BH q=", format(nominal_models$BH_global_25, digits = 4), ")"), collapse = "; ")
} else {
    "none"
}

methods_lines <- c("# STEP04E1O Methods shell", "", "For each model-eligible broad cell type, donor-level VDR expression and four prospectively locked signed VDR-program scores were compared between Control and RIF donors. Two-sided exact label-permutation tests based on the RIF-minus-Control mean difference were used as the primary inferential procedure. All possible allocations were enumerated, yielding 84 allocations for the 3-versus-6 design and 35 allocations for the 3-versus-4 unciliated epithelial design. Multiplicity was controlled using the Benjamini-Hochberg procedure across all 25 planned primary models. Welch tests, 95% confidence intervals, Hedges g, rank-biserial effects and leave-one-donor-out analyses were reported as secondary estimates.")

results_lines <- c("# STEP04E1O Results shell", "", paste0("Across 25 planned donor-level models, ", nrow(nominal_models), 
    " met nominal exact P<0.05, whereas ", nrow(global_BH_models), " remained significant after global Benjamini-Hochberg correction."), 
    paste0("The nominal models were: ", nominal_description, "."), paste0("Leave-one-donor-out direction instability was observed in ", 
        sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE"), " models."), paste0("Retaining cluster 8 changed the unciliated epithelial effect direction for ", 
        nrow(direction_changed_sensitivity), " of five endpoints."), "Ciliated epithelial profiles were retained for descriptive localization only because only two donor profiles met the 20-cell threshold.")

write_text(methods_lines, file.path(manuscript_text_root, paste0(run_id, "_Methods.md")))

write_text(results_lines, file.path(manuscript_text_root, paste0(run_id, "_Results.md")))

parameter_manifest <- list(script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_DATA_TABLE_GROUPED_SUMMARY_INTEGER_DOUBLE_TYPE_MISMATCH", 
    grouped_summary_policy = "EXPLICIT_INTEGER_COUNTS_AND_NUMERIC_CONTINUOUS_SUMMARIES", dataset_id = dataset_id, 
    integration_version = integration_version, pseudobulk_version = pseudobulk_version, model_version = model_version, 
    upstream_gate = upstream_gate$gate_status[[1]], eligible_broad_cell_types = eligible_types, endpoints = endpoint_order, 
    primary_test = "TWO_SIDED_EXACT_LABEL_PERMUTATION_OF_RIF_MINUS_CONTROL_MEAN_DIFFERENCE", full_design_permutation_n = 84L, 
    unciliated_primary_permutation_n = 35L, primary_multiplicity = "BH_ACROSS_ALL_25_PLANNED_MODELS", secondary_statistics = c("WELCH_T_TEST", 
        "WELCH_95_PERCENT_CI", "HEDGES_G", "RANK_BISERIAL", "LEAVE_ONE_DONOR_OUT"), cluster8_policy = "ATLAS_SENSITIVITY_ONLY", 
    ciliated_policy = "DESCRIPTIVE_ONLY_NO_GROUP_TEST", independent_validation_review_authorized = all_gate_checks_pass, 
    gene_level_pseudobulk_DE_authorized = FALSE, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), session_path)

run_end <- Sys.time()

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, supersedes_partial_execution = "V1_0_DATA_TABLE_GROUPED_SUMMARY_INTEGER_DOUBLE_TYPE_MISMATCH", 
    dataset_id = dataset_id, integration_version = integration_version, pseudobulk_version = pseudobulk_version, 
    model_version = model_version, started_at = as.character(run_start), completed_at = as.character(run_end), 
    elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 2), gate_status = gate_status, 
    primary_model_n = nrow(primary_models), nominal_exact_P_lt_0_05_n = nrow(nominal_models), global_BH_lt_0_05_n = nrow(global_BH_models), 
    LOO_direction_unstable_n = sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE"), cluster8_direction_changed_n = nrow(direction_changed_sensitivity), 
    independent_validation_review_authorized = all_gate_checks_pass, result_lock_authorized = FALSE, output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

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

log_event("Completed", run_id, "gate_status=", gate_status, "primary_models=", nrow(primary_models), "nominal=", 
    nrow(nominal_models), "global_BH=", nrow(global_BH_models))

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

cat("STEP04E1O GSE183837 donor-level models completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Primary models: ", nrow(primary_models), "\n", sep = "")

cat("Nonconstant testable models: ", sum(!primary_models$constant_endpoint), "\n", sep = "")

cat("Nominal exact P < 0.05: ", nrow(nominal_models), "\n", sep = "")

cat("Global BH q < 0.05: ", nrow(global_BH_models), "\n", sep = "")

cat("LOO direction-unstable models: ", sum(primary_models$LOO_review_flag == "LOO_DIRECTION_UNSTABLE"), "\n", sep = "")

cat("Cluster-8 sensitivity direction changes: ", nrow(direction_changed_sensitivity), "\n", sep = "")

cat("Independent validation review authorized: ", all_gate_checks_pass, "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nDonor-model outputs are complete in the project folders.\n", sep = "")
}

