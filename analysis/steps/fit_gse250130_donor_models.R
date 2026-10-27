#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C07: Fit temporal and LH+7 donor-level models

options(stringsAsFactors = FALSE, warn = 1)

set.seed(20260805)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "fit_gse250130_donor_models.R"

script_version <- "STEP04E1G_GSE250130_DONOR_VDR_MODELS_v1.1_UNBALANCED_FIXED"

input_version <- "SCI_v1_1_memory_safe"

pseudobulk_version <- "PBVDR_v1_0"

model_version <- "VDRMODELS_v1_1"

run_start <- Sys.time()

run_id <- paste0("STEP04E1G_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_windows_path_length <- function(path, max_chars = 245L) {
    normalized <- normalize_slashes(path)
    path_n <- nchar(normalized, type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Output path is too long for safe Windows/R writing (", path_n, " characters; safe limit ", max_chars, 
            "):\n", normalized)
    }
    invisible(path)
}

write_csv_utf8 <- function(x, path) {
    assert_windows_path_length(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    assert_windows_path_length(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    out <- tryCatch(.vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), 
        error = function(e) NULL)
    if (is.null(out)) {
        out <- .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    out
}

safe_package_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) 
        return(NA_character_)
    as.character(utils::packageVersion(pkg))
}

collapse_unique <- function(x, sep = "; ") {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) 
        ""
    else paste(x, collapse = sep)
}

safe_numeric <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}

format_p <- function(x) {
    ifelse(is.na(x), "NA", ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x)))
}

sign_label <- function(effect) {
    ifelse(is.na(effect), "NA", ifelse(effect > 0, "POSITIVE", ifelse(effect < 0, "NEGATIVE", "ZERO")))
}

save_plot_pair <- function(plot_object, stem, width, height) {
    png_path <- file.path(figure_root, paste0(stem, ".png"))
    pdf_path <- file.path(figure_root, paste0(stem, ".pdf"))
    assert_windows_path_length(png_path)
    assert_windows_path_length(pdf_path)
    ggplot2::ggsave(filename = png_path, plot = plot_object, width = width, height = height, units = "in", dpi = 300, 
        bg = "white")
    ggplot2::ggsave(filename = pdf_path, plot = plot_object, width = width, height = height, units = "in", device = grDevices::cairo_pdf, 
        bg = "white")
    c(png_path, pdf_path)
}

oneway_F <- function(y, group) {
    keep <- is.finite(y) & !is.na(group)
    y <- y[keep]
    group <- droplevels(factor(group[keep]))
    n <- length(y)
    k <- nlevels(group)
    if (n <= k || k < 2L) 
        return(NA_real_)
    grand_mean <- mean(y)
    group_means <- tapply(y, group, mean)
    group_n <- table(group)
    ss_between <- sum(as.numeric(group_n) * (as.numeric(group_means) - grand_mean)^2)
    fitted_group_mean <- group_means[as.character(group)]
    ss_within <- sum((y - fitted_group_mean)^2)
    df_between <- k - 1L
    df_within <- n - k
    if (df_within <= 0L || ss_within <= 0) {
        return(NA_real_)
    }
    (ss_between/df_between)/(ss_within/df_within)
}

permutation_oneway <- function(y, group, n_perm = 5000L, seed = 20260805L) {
    keep <- is.finite(y) & !is.na(group)
    y <- as.numeric(y[keep])
    group <- droplevels(factor(group[keep]))
    observed_F <- oneway_F(y, group)
    if (!is.finite(observed_F)) {
        return(list(F = observed_F, p = NA_real_, permutation_n = 0L))
    }
    set.seed(seed)
    perm_F <- replicate(n_perm, oneway_F(sample(y, length(y), replace = FALSE), group))
    perm_F <- perm_F[is.finite(perm_F)]
    p_value <- (1 + sum(perm_F >= observed_F))/(1 + length(perm_F))
    list(F = observed_F, p = p_value, permutation_n = length(perm_F))
}

exact_two_group_test <- function(y, group, reference_level, comparison_level, maximum_exact_combinations = 50000L, 
    Monte_Carlo_n = 20000L, seed = 20260805L) {
    keep <- (is.finite(y) & !is.na(group) & group %in% c(reference_level, comparison_level))
    y <- as.numeric(y[keep])
    group <- as.character(group[keep])
    reference_values <- y[group == reference_level]
    comparison_values <- y[group == comparison_level]
    n_reference <- length(reference_values)
    n_comparison <- length(comparison_values)
    n_total <- length(y)
    if (n_reference < 2L || n_comparison < 2L) {
        return(list(estimate = NA_real_, p = NA_real_, test_mode = "INSUFFICIENT_REPLICATION", permutation_n = 0L))
    }
    observed <- (mean(comparison_values) - mean(reference_values))
    total_sum <- sum(y)
    combination_n <- choose(n_total, n_comparison)
    if (is.finite(combination_n) && combination_n <= maximum_exact_combinations) {
        combinations <- utils::combn(n_total, n_comparison)
        comparison_sums <- colSums(matrix(y[combinations], nrow = n_comparison))
        permuted_differences <- (comparison_sums/n_comparison) - ((total_sum - comparison_sums)/n_reference)
        p_value <- mean(abs(permuted_differences) >= abs(observed) - 9.9999999999999998e-13)
        test_mode <- "EXACT"
        permutation_n <- length(permuted_differences)
    }
    else {
        set.seed(seed)
        permuted_differences <- replicate(Monte_Carlo_n, {
            comparison_index <- sample.int(n_total, n_comparison, replace = FALSE)
            comparison_sum <- sum(y[comparison_index])
            (comparison_sum/n_comparison) - ((total_sum - comparison_sum)/n_reference)
        })
        p_value <- (1 + sum(abs(permuted_differences) >= abs(observed) - 9.9999999999999998e-13))/(1 + length(permuted_differences))
        test_mode <- "MONTE_CARLO"
        permutation_n <- length(permuted_differences)
    }
    list(estimate = observed, p = p_value, test_mode = test_mode, permutation_n = permutation_n)
}

welch_summary <- function(y, group, reference_level, comparison_level) {
    keep <- (is.finite(y) & !is.na(group) & group %in% c(reference_level, comparison_level))
    y <- as.numeric(y[keep])
    group <- as.character(group[keep])
    reference_values <- y[group == reference_level]
    comparison_values <- y[group == comparison_level]
    if (length(reference_values) < 2L || length(comparison_values) < 2L) {
        return(list(estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
    }
    test <- stats::t.test(comparison_values, reference_values, var.equal = FALSE, conf.level = 0.94999999999999996)
    list(estimate = (mean(comparison_values) - mean(reference_values)), ci_low = unname(test$conf.int[1]), ci_high = unname(test$conf.int[2]), 
        p = unname(test$p.value))
}

hedges_g <- function(y, group, reference_level, comparison_level) {
    keep <- (is.finite(y) & !is.na(group) & group %in% c(reference_level, comparison_level))
    y <- as.numeric(y[keep])
    group <- as.character(group[keep])
    x_ref <- y[group == reference_level]
    x_cmp <- y[group == comparison_level]
    n_ref <- length(x_ref)
    n_cmp <- length(x_cmp)
    if (n_ref < 2L || n_cmp < 2L) {
        return(NA_real_)
    }
    pooled_variance <- ((n_ref - 1) * stats::var(x_ref) + (n_cmp - 1) * stats::var(x_cmp))/(n_ref + n_cmp - 2)
    if (!is.finite(pooled_variance) || pooled_variance <= 0) {
        return(NA_real_)
    }
    cohen_d <- (mean(x_cmp) - mean(x_ref))/sqrt(pooled_variance)
    correction <- (1 - 3/(4 * (n_ref + n_cmp) - 9))
    correction * cohen_d
}

leave_one_out_difference <- function(y, group, reference_level, comparison_level) {
    keep <- (is.finite(y) & !is.na(group) & group %in% c(reference_level, comparison_level))
    y <- as.numeric(y[keep])
    group <- as.character(group[keep])
    if (sum(group == reference_level) < 3L || sum(group == comparison_level) < 3L) {
        return(c(loo_min = NA_real_, loo_max = NA_real_, loo_sign_consistency = NA_real_))
    }
    full_effect <- (mean(y[group == comparison_level]) - mean(y[group == reference_level]))
    loo_effects <- vapply(seq_along(y), function(i) {
        y_i <- y[-i]
        group_i <- group[-i]
        (mean(y_i[group_i == comparison_level]) - mean(y_i[group_i == reference_level]))
    }, numeric(1))
    full_sign <- sign(full_effect)
    c(loo_min = min(loo_effects), loo_max = max(loo_effects), loo_sign_consistency = mean(sign(loo_effects) == 
        full_sign))
}

expected_temporal_counts <- function(eligibility_table, cell_type) {
    one <- eligibility_table[eligibility_table$broad_cell_type == cell_type, , drop = FALSE]
    if (nrow(one) != 1L) {
        stop("Eligibility table does not contain exactly one row for: ", cell_type)
    }
    counts <- c(`3` = safe_numeric(one$fertile_LH3_n), `5` = safe_numeric(one$fertile_LH5_n), `7` = safe_numeric(one$fertile_LH7_n), 
        `9` = safe_numeric(one$fertile_LH9_n), `11` = safe_numeric(one$fertile_LH11_n))
    if (any(!is.finite(counts)) || any(counts < 2L)) {
        stop("Temporal eligibility counts are invalid for ", cell_type, ": ", paste(names(counts), counts, sep = "=", 
            collapse = ";"))
    }
    as.integer(counts)
}

group_summary <- function(data, value_column, group_columns) {
    dt <- data.table::as.data.table(data)
    dt[, {
        values <- safe_numeric(get(value_column))
        values <- values[is.finite(values)]
        n <- length(values)
        mean_value <- if (n) 
            mean(values)
        else NA_real_
        sd_value <- if (n >= 2L) 
            stats::sd(values)
        else NA_real_
        se_value <- if (n >= 2L) 
            sd_value/sqrt(n)
        else NA_real_
        critical <- if (n >= 2L) {
            stats::qt(0.97499999999999998, df = n - 1L)
        }
        else {
            NA_real_
        }
        list(n = n, mean = mean_value, sd = sd_value, se = se_value, ci_low = mean_value - critical * se_value, 
            ci_high = mean_value + critical * se_value, median = if (n) stats::median(values) else NA_real_, IQR = if (n) stats::IQR(values) else NA_real_)
    }, by = group_columns]
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

required_packages <- c("data.table", "ggplot2", "jsonlite", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

pb_root <- file.path(results_root, "S4E1_sc", input_version, "pb", dataset_id, pseudobulk_version)

pb_table_root <- file.path(pb_root, "tables")

pb_object_root <- file.path(pb_root, "objects")

program_score_path <- file.path(pb_table_root, "05_VDR_program_scores.csv")

VDR_expression_path <- file.path(pb_table_root, "08_VDR_expression.csv")

eligibility_path <- file.path(pb_table_root, "04_model_eligibility.csv")

pb_metadata_path <- file.path(pb_table_root, "02_pb_metadata.csv")

pb_status_path <- file.path(pb_object_root, "pb_status.csv")

output_root <- file.path(results_root, "S4E1_sc", input_version, "models", dataset_id, model_version)

table_root <- file.path(output_root, "tables")

source_data_root <- file.path(output_root, "source_data")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1G_", format(run_start, "%Y%m%d_%H%M%S"), "_modelQC"))

figure_root <- file.path(run_qc_root, "figures")

for (p in c(execution_root, parameter_root, environment_root, manuscript_text_root, output_root, table_root, source_data_root, 
    run_qc_root, figure_root)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(execution_root, paste0(run_id, "_execution.log"))

run_manifest_path <- file.path(execution_root, paste0(run_id, "_RUN_MANIFEST.csv"))

parameter_path <- file.path(parameter_root, paste0(run_id, "_parameters.json"))

session_path <- file.path(environment_root, paste0("sessionInfo_", run_id, ".txt"))

desktop_candidates <- .vdr_report_root()

desktop_candidates <- .vdr_report_root()

zip_parent <- .vdr_report_root()

handoff_zip <- file.path(zip_parent, paste0(run_id, "_SMALL_REPORT.zip"))

planned_paths <- data.frame(object = c("output_root", "table_root", "source_data_root", "run_qc_root", "figure_root", 
    "log_path", "run_manifest_path", "parameter_path", "session_path", "handoff_zip"), path = normalize_slashes(c(output_root, 
    table_root, source_data_root, run_qc_root, figure_root, log_path, run_manifest_path, parameter_path, session_path, 
    handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

if (.Platform$OS.type == "windows" && any(planned_paths$character_n > 245L)) {
    bad <- planned_paths[planned_paths$character_n > 245L, , drop = FALSE]
    stop("One or more planned paths are too long:\n", paste(paste0(bad$character_n, " chars: ", bad$path), collapse = "\n"))
}

required_inputs <- c(program_score_path, VDR_expression_path, eligibility_path, pb_metadata_path, pb_status_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

pb_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

pb_qc_dirs <- pb_qc_dirs[grepl("^STEP04E1F_[0-9]{8}_[0-9]{6}_pbQC$", basename(pb_qc_dirs), perl = TRUE)]

if (!length(pb_qc_dirs)) {
    stop("No completed STEP04E1F pseudobulk-QC folder found.")
}

latest_pb_qc <- pb_qc_dirs[which.max(file.info(pb_qc_dirs)$mtime)]

pb_gate_path <- file.path(latest_pb_qc, "PBVDR_GATE.csv")

pb_gate <- safe_read_csv(pb_gate_path)

pb_status <- safe_read_csv(pb_status_path)

if (is.null(pb_gate) || nrow(pb_gate) != 1L || pb_gate$gate_status[[1]] != "PASS_PSEUDOBULK_AND_VDR_LOCALIZATION_READY_FOR_DONOR_MODELS" || 
    !isTRUE(pb_gate$donor_level_modeling_authorized[[1]])) {
    stop("STEP04E1F Gate does not authorize donor-level models.")
}

if (is.null(pb_status) || nrow(pb_status) != 1L || pb_status$status[[1]] != "COMPLETE" || pb_status$pseudobulk_n[[1]] != 
    222L || pb_status$aggregated_cell_n[[1]] != 236636L) {
    stop("PBVDR pseudobulk status failed validation.")
}

scores <- .vdr_fread(program_score_path, data.table = FALSE, showProgress = FALSE)

VDR_expression <- .vdr_fread(VDR_expression_path, data.table = FALSE, showProgress = FALSE)

eligibility <- safe_read_csv(eligibility_path)

pb_metadata <- safe_read_csv(pb_metadata_path)

required_score_columns <- c("pb_id", "program_name", "score", "sample_accession", "donor_id", "broad_cell_type", 
    "clinical_group", "LH_day", "temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "LH7_reference_group", 
    "LH7_fertile_subgroup", "minimum_cell_threshold_pass", "score_eligible")

missing_score_columns <- setdiff(required_score_columns, colnames(scores))

if (length(missing_score_columns)) {
    stop("Program-score table is missing columns: ", paste(missing_score_columns, collapse = ", "))
}

required_VDR_columns <- c("pb_id", "sample_accession", "donor_id", "broad_cell_type", "clinical_group", "LH_day", 
    "temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "LH7_reference_group", "LH7_fertile_subgroup", 
    "minimum_cell_threshold_pass", "VDR_log2_CPM_plus_1")

missing_VDR_columns <- setdiff(required_VDR_columns, colnames(VDR_expression))

if (length(missing_VDR_columns)) {
    stop("VDR-expression table is missing columns: ", paste(missing_VDR_columns, collapse = ", "))
}

if (nrow(scores) != 888L || length(unique(scores$pb_id)) != 222L || length(unique(scores$program_name)) != 4L || 
    nrow(VDR_expression) != 222L || nrow(eligibility) != 8L || nrow(pb_metadata) != 222L) {
    stop("PBVDR input identity checks failed.")
}

temporal_types <- as.character(eligibility$broad_cell_type[eligibility$fertile_temporal_model_eligible])

LH7_types <- as.character(eligibility$broad_cell_type[eligibility$LH7_RIF_vs_fertile_model_eligible])

program_names <- sort(unique(as.character(scores$program_name)))

expected_programs <- c("OCCUPANCY_SUPPORTED_284", "PROMOTER_OCCUPANCY_84", "RNA_ALL_436", "SHARED_OCCUPANCY_46")

if (!setequal(program_names, expected_programs)) {
    stop("The four locked VDR programs were not recovered exactly.")
}

if (length(temporal_types) != 5L || length(LH7_types) != 6L) {
    stop("Model-eligible cell-type counts do not match the Gate.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, input_version = input_version, 
    pseudobulk_version = pseudobulk_version, model_version = model_version, supersedes_failed_model_version = "VDRMODELS_v1_0_HARD_CODED_SAMPLE_COUNT_QC_FAILURE", 
    dataset_id = dataset_id, temporal_question = "FERTILE_LH3_LH5_LH7_LH9_LH11", temporal_replication_policy = "USE_ACTUAL_GATE_ELIGIBLE_COUNTS_MINIMUM_2_DONORS_PER_TIMEPOINT", 
    temporal_design_balance = "UNBALANCED_ALLOWED_AND_EXPECTED_FOR_CILIATED_EPITHELIAL_LH11_N2", supersedes_failed_model_version = "VDRMODELS_v1_0_HARD_CODED_SAMPLE_COUNT_QC_FAILURE", 
    temporal_primary_test = "5000_LABEL_PERMUTATION_ONE_WAY_ANOVA", temporal_secondary_test = "CLASSICAL_ONE_WAY_ANOVA", 
    temporal_planned_pairwise = c("LH5_vs_LH3", "LH7_vs_LH5", "LH9_vs_LH7", "LH11_vs_LH9", "LH7_vs_LH3", "LH11_vs_LH3", 
        "LH11_vs_LH7"), LH7_primary_question = "RIF_VS_ALL_FERTILE_LH7", LH7_primary_test = "EXACT_PERMUTATION_MEAN_DIFFERENCE", 
    LH7_secondary_test = "WELCH_T_TEST", LH7_sensitivity_questions = c("RIF_VS_FERTILE_LH7_STANDARD", "RIF_VS_FERTILE_LH7_PGT"), 
    multiple_testing = "BENJAMINI_HOCHBERG_WITHIN_PREDEFINED_MODEL_FAMILY", minimum_cell_policy = "minimum_cell_threshold_pass_TRUE", 
    program_score_role = "RELATIVE_ALIGNMENT_WITH_LOCKED_VDR_INTERACTION_DIRECTION", composition_inference_authorized = FALSE, 
    subtype_inference_authorized = FALSE, random_seed = 20260805L, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

log_event("Starting", run_id, "temporal_types=", collapse_unique(temporal_types), "LH7_types=", collapse_unique(LH7_types))

scores$LH_day <- safe_numeric(scores$LH_day)

VDR_expression$LH_day <- safe_numeric(VDR_expression$LH_day)

scores$score <- safe_numeric(scores$score)

VDR_expression$VDR_log2_CPM_plus_1 <- safe_numeric(VDR_expression$VDR_log2_CPM_plus_1)

scores$minimum_cell_threshold_pass <- as.logical(scores$minimum_cell_threshold_pass)

scores$score_eligible <- as.logical(scores$score_eligible)

VDR_expression$minimum_cell_threshold_pass <- as.logical(VDR_expression$minimum_cell_threshold_pass)

temporal_score_data <- scores[scores$clinical_group == "FERTILE" & scores$temporal_fertile_eligible & scores$minimum_cell_threshold_pass & 
    scores$score_eligible & scores$broad_cell_type %in% temporal_types, , drop = FALSE]

LH7_score_data <- scores[scores$LH7_RIF_comparison_eligible & scores$minimum_cell_threshold_pass & scores$score_eligible & 
    scores$broad_cell_type %in% LH7_types, , drop = FALSE]

temporal_VDR_data <- VDR_expression[VDR_expression$clinical_group == "FERTILE" & VDR_expression$temporal_fertile_eligible & 
    VDR_expression$minimum_cell_threshold_pass & VDR_expression$broad_cell_type %in% temporal_types, , drop = FALSE]

LH7_VDR_data <- VDR_expression[VDR_expression$LH7_RIF_comparison_eligible & VDR_expression$minimum_cell_threshold_pass & 
    VDR_expression$broad_cell_type %in% LH7_types, , drop = FALSE]

input_audit_rows <- list()

for (cell_type in temporal_types) {
    one <- temporal_VDR_data[temporal_VDR_data$broad_cell_type == cell_type, , drop = FALSE]
    observed <- table(factor(one$LH_day, levels = c(3, 5, 7, 9, 11)))
    expected <- expected_temporal_counts(eligibility_table = eligibility, cell_type = cell_type)
    input_audit_rows[[length(input_audit_rows) + 1L]] <- data.frame(model_family = "FERTILE_TEMPORAL", broad_cell_type = cell_type, 
        expected_groups = paste0("LH", names(expected), "=", as.integer(expected), collapse = ";"), observed_groups = paste0("LH", 
            names(observed), "=", as.integer(observed), collapse = ";"), donor_n = length(unique(one$donor_id)), 
        minimum_group_n = min(as.integer(observed)), unbalanced_design = length(unique(as.integer(observed))) > 
            1L, identity_pass = identical(as.integer(observed), as.integer(expected)) && all(as.integer(observed) >= 
            2L), stringsAsFactors = FALSE)
}

for (cell_type in LH7_types) {
    one <- LH7_VDR_data[LH7_VDR_data$broad_cell_type == cell_type, , drop = FALSE]
    observed <- table(factor(one$clinical_group, levels = c("FERTILE", "RIF")))
    input_audit_rows[[length(input_audit_rows) + 1L]] <- data.frame(model_family = "LH7_RIF_VS_FERTILE", broad_cell_type = cell_type, 
        expected_groups = "FERTILE>=3;RIF>=3", observed_groups = paste0(names(observed), "=", as.integer(observed), 
            collapse = ";"), donor_n = length(unique(one$donor_id)), minimum_group_n = min(as.integer(observed)), 
        unbalanced_design = length(unique(as.integer(observed))) > 1L, identity_pass = all(as.integer(observed) >= 
            3L), stringsAsFactors = FALSE)
}

input_audit <- do.call(rbind, input_audit_rows)

write_csv_utf8(input_audit, file.path(run_qc_root, "01_model_input_audit.csv"))

if (!all(input_audit$identity_pass)) {
    stop("One or more donor-level model inputs failed sample-count QC.")
}

temporal_program_rows <- list()

temporal_pairwise_rows <- list()

temporal_summary_rows <- list()

planned_pairs <- data.frame(comparison = c("LH5_vs_LH3", "LH7_vs_LH5", "LH9_vs_LH7", "LH11_vs_LH9", "LH7_vs_LH3", 
    "LH11_vs_LH3", "LH11_vs_LH7"), reference_day = c(3, 5, 7, 9, 3, 3, 7), comparison_day = c(5, 7, 9, 11, 7, 11, 
    11), stringsAsFactors = FALSE)

model_counter <- 0L

for (cell_type in temporal_types) {
    for (program_name in program_names) {
        model_counter <- model_counter + 1L
        one <- temporal_score_data[temporal_score_data$broad_cell_type == cell_type & temporal_score_data$program_name == 
            program_name, , drop = FALSE]
        one$LH_factor <- factor(one$LH_day, levels = c(3, 5, 7, 9, 11))
        expected_counts <- expected_temporal_counts(eligibility_table = eligibility, cell_type = cell_type)
        observed_counts <- as.integer(table(one$LH_factor))
        if (nrow(one) != sum(expected_counts) || !identical(observed_counts, as.integer(expected_counts)) || any(observed_counts < 
            2L)) {
            stop("Temporal score input mismatch for ", cell_type, " / ", program_name, ". Expected ", paste(names(expected_counts), 
                expected_counts, sep = "=", collapse = ";"), "; observed ", paste(levels(one$LH_factor), observed_counts, 
                sep = "=", collapse = ";"))
        }
        permutation <- permutation_oneway(y = one$score, group = one$LH_factor, n_perm = 5000L, seed = 20260805L + 
            model_counter)
        fit <- stats::lm(score ~ LH_factor, data = one)
        anova_table <- stats::anova(fit)
        classical_p <- anova_table$`Pr(>F)`[1]
        classical_F <- anova_table$`F value`[1]
        group_means <- tapply(one$score, one$LH_factor, mean)
        total_ss <- sum((one$score - mean(one$score))^2)
        fitted_ss <- sum((fitted(fit) - mean(one$score))^2)
        eta_squared <- if (is.finite(total_ss) && total_ss > 0) {
            fitted_ss/total_ss
        }
        else {
            NA_real_
        }
        temporal_program_rows[[length(temporal_program_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, 
            program_name = program_name, donor_n = nrow(one), LH3_n = sum(one$LH_day == 3), LH5_n = sum(one$LH_day == 
                5), LH7_n = sum(one$LH_day == 7), LH9_n = sum(one$LH_day == 9), LH11_n = sum(one$LH_day == 11), 
            LH3_mean = unname(group_means["3"]), LH5_mean = unname(group_means["5"]), LH7_mean = unname(group_means["7"]), 
            LH9_mean = unname(group_means["9"]), LH11_mean = unname(group_means["11"]), permutation_F = permutation$F, 
            permutation_p = permutation$p, permutation_n = permutation$permutation_n, classical_F = classical_F, 
            classical_p = classical_p, eta_squared = eta_squared, stringsAsFactors = FALSE)
        summary_one <- group_summary(data = one, value_column = "score", group_columns = c("broad_cell_type", "program_name", 
            "LH_day"))
        temporal_summary_rows[[length(temporal_summary_rows) + 1L]] <- as.data.frame(summary_one)
        for (pair_index in seq_len(nrow(planned_pairs))) {
            pair <- planned_pairs[pair_index, ]
            pair_data <- one[one$LH_day %in% c(pair$reference_day, pair$comparison_day), , drop = FALSE]
            pair_group <- paste0("LH", pair_data$LH_day)
            reference_level <- paste0("LH", pair$reference_day)
            comparison_level <- paste0("LH", pair$comparison_day)
            exact_test <- exact_two_group_test(y = pair_data$score, group = pair_group, reference_level = reference_level, 
                comparison_level = comparison_level, seed = 20260805L + model_counter * 100L + pair_index)
            welch <- welch_summary(y = pair_data$score, group = pair_group, reference_level = reference_level, 
                comparison_level = comparison_level)
            temporal_pairwise_rows[[length(temporal_pairwise_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, 
                program_name = program_name, comparison = pair$comparison, reference_day = pair$reference_day, 
                comparison_day = pair$comparison_day, reference_n = sum(pair_group == reference_level), comparison_n = sum(pair_group == 
                  comparison_level), mean_difference = exact_test$estimate, exact_permutation_p = exact_test$p, 
                exact_test_mode = exact_test$test_mode, exact_permutation_n = exact_test$permutation_n, Welch_ci_low = welch$ci_low, 
                Welch_ci_high = welch$ci_high, Welch_p = welch$p, Hedges_g = hedges_g(y = pair_data$score, group = pair_group, 
                  reference_level = reference_level, comparison_level = comparison_level), stringsAsFactors = FALSE)
        }
        log_event("DONE temporal program model", model_counter, "/", length(temporal_types) * length(program_names), 
            cell_type, program_name, "permutation_p=", format_p(permutation$p))
    }
}

temporal_program <- do.call(rbind, temporal_program_rows)

temporal_program$FDR_BH <- stats::p.adjust(temporal_program$permutation_p, method = "BH")

temporal_program$trajectory_direction <- ifelse(temporal_program$LH11_mean > temporal_program$LH3_mean, "LH11_ABOVE_LH3", 
    ifelse(temporal_program$LH11_mean < temporal_program$LH3_mean, "LH11_BELOW_LH3", "LH11_EQUALS_LH3"))

temporal_pairwise <- do.call(rbind, temporal_pairwise_rows)

temporal_pairwise$FDR_BH_within_comparison <- ave(temporal_pairwise$exact_permutation_p, temporal_pairwise$comparison, 
    FUN = function(p) {
        stats::p.adjust(p, method = "BH")
    })

temporal_pairwise$effect_direction <- sign_label(temporal_pairwise$mean_difference)

temporal_program_summary <- do.call(rbind, temporal_summary_rows)

write_csv_utf8(temporal_program, file.path(table_root, "01_temporal_program_omnibus.csv"))

write_csv_utf8(temporal_pairwise, file.path(table_root, "02_temporal_program_pairwise.csv"))

write_csv_utf8(temporal_program_summary, file.path(source_data_root, "source_temporal_program_summary.csv"))

LH7_primary_rows <- list()

LH7_sensitivity_rows <- list()

model_counter <- 0L

for (cell_type in LH7_types) {
    for (program_name in program_names) {
        model_counter <- model_counter + 1L
        one <- LH7_score_data[LH7_score_data$broad_cell_type == cell_type & LH7_score_data$program_name == program_name, 
            , drop = FALSE]
        one$clinical_group <- factor(one$clinical_group, levels = c("FERTILE", "RIF"))
        exact_test <- exact_two_group_test(y = one$score, group = one$clinical_group, reference_level = "FERTILE", 
            comparison_level = "RIF", seed = 20260805L + model_counter)
        welch <- welch_summary(y = one$score, group = one$clinical_group, reference_level = "FERTILE", comparison_level = "RIF")
        loo <- leave_one_out_difference(y = one$score, group = one$clinical_group, reference_level = "FERTILE", 
            comparison_level = "RIF")
        LH7_primary_rows[[length(LH7_primary_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, program_name = program_name, 
            fertile_n = sum(one$clinical_group == "FERTILE"), RIF_n = sum(one$clinical_group == "RIF"), fertile_mean = mean(one$score[one$clinical_group == 
                "FERTILE"]), RIF_mean = mean(one$score[one$clinical_group == "RIF"]), RIF_minus_fertile = exact_test$estimate, 
            exact_permutation_p = exact_test$p, exact_test_mode = exact_test$test_mode, exact_permutation_n = exact_test$permutation_n, 
            Welch_ci_low = welch$ci_low, Welch_ci_high = welch$ci_high, Welch_p = welch$p, Hedges_g = hedges_g(y = one$score, 
                group = one$clinical_group, reference_level = "FERTILE", comparison_level = "RIF"), LOO_min = unname(loo["loo_min"]), 
            LOO_max = unname(loo["loo_max"]), LOO_sign_consistency = unname(loo["loo_sign_consistency"]), stringsAsFactors = FALSE)
        for (fertile_subgroup in c("LH7_STANDARD", "LH7_PGT")) {
            sensitivity <- one[one$clinical_group == "RIF" | (one$clinical_group == "FERTILE" & one$LH7_fertile_subgroup == 
                fertile_subgroup), , drop = FALSE]
            sensitivity_group <- ifelse(sensitivity$clinical_group == "RIF", "RIF", fertile_subgroup)
            sensitivity_test <- exact_two_group_test(y = sensitivity$score, group = sensitivity_group, reference_level = fertile_subgroup, 
                comparison_level = "RIF", seed = 20260805L + model_counter * 100L + match(fertile_subgroup, c("LH7_STANDARD", 
                  "LH7_PGT")))
            sensitivity_welch <- welch_summary(y = sensitivity$score, group = sensitivity_group, reference_level = fertile_subgroup, 
                comparison_level = "RIF")
            LH7_sensitivity_rows[[length(LH7_sensitivity_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, 
                program_name = program_name, sensitivity_reference = fertile_subgroup, fertile_subgroup_n = sum(sensitivity_group == 
                  fertile_subgroup), RIF_n = sum(sensitivity_group == "RIF"), RIF_minus_reference = sensitivity_test$estimate, 
                exact_permutation_p = sensitivity_test$p, exact_test_mode = sensitivity_test$test_mode, exact_permutation_n = sensitivity_test$permutation_n, 
                Welch_ci_low = sensitivity_welch$ci_low, Welch_ci_high = sensitivity_welch$ci_high, Welch_p = sensitivity_welch$p, 
                Hedges_g = hedges_g(y = sensitivity$score, group = sensitivity_group, reference_level = fertile_subgroup, 
                  comparison_level = "RIF"), stringsAsFactors = FALSE)
        }
        log_event("DONE LH7 program model", model_counter, "/", length(LH7_types) * length(program_names), cell_type, 
            program_name, "effect=", round(exact_test$estimate, 4), "p=", format_p(exact_test$p))
    }
}

LH7_primary <- do.call(rbind, LH7_primary_rows)

LH7_primary$FDR_BH <- stats::p.adjust(LH7_primary$exact_permutation_p, method = "BH")

LH7_primary$effect_direction <- sign_label(LH7_primary$RIF_minus_fertile)

LH7_primary$LOO_stable_sign <- (LH7_primary$LOO_sign_consistency >= 0.80000000000000004)

LH7_sensitivity <- do.call(rbind, LH7_sensitivity_rows)

LH7_sensitivity$FDR_BH_within_reference <- ave(LH7_sensitivity$exact_permutation_p, LH7_sensitivity$sensitivity_reference, 
    FUN = function(p) {
        stats::p.adjust(p, method = "BH")
    })

primary_lookup <- LH7_primary[, c("broad_cell_type", "program_name", "RIF_minus_fertile")]

colnames(primary_lookup)[3] <- "primary_RIF_minus_fertile"

LH7_sensitivity <- merge(LH7_sensitivity, primary_lookup, by = c("broad_cell_type", "program_name"), all.x = TRUE, 
    sort = FALSE)

LH7_sensitivity$sign_concordant_with_primary <- (sign(LH7_sensitivity$RIF_minus_reference) == sign(LH7_sensitivity$primary_RIF_minus_fertile))

write_csv_utf8(LH7_primary, file.path(table_root, "03_LH7_program_primary.csv"))

write_csv_utf8(LH7_sensitivity, file.path(table_root, "04_LH7_program_sensitivity.csv"))

temporal_VDR_rows <- list()

temporal_VDR_pair_rows <- list()

temporal_VDR_summary_rows <- list()

for (type_index in seq_along(temporal_types)) {
    cell_type <- temporal_types[type_index]
    one <- temporal_VDR_data[temporal_VDR_data$broad_cell_type == cell_type, , drop = FALSE]
    one$LH_factor <- factor(one$LH_day, levels = c(3, 5, 7, 9, 11))
    expected_counts <- expected_temporal_counts(eligibility_table = eligibility, cell_type = cell_type)
    observed_counts <- as.integer(table(one$LH_factor))
    if (nrow(one) != sum(expected_counts) || !identical(observed_counts, as.integer(expected_counts)) || any(observed_counts < 
        2L)) {
        stop("Temporal VDR-expression input mismatch for ", cell_type)
    }
    permutation <- permutation_oneway(y = one$VDR_log2_CPM_plus_1, group = one$LH_factor, n_perm = 5000L, seed = 20260805L + 
        10000L + type_index)
    fit <- stats::lm(VDR_log2_CPM_plus_1 ~ LH_factor, data = one)
    anova_table <- stats::anova(fit)
    group_means <- tapply(one$VDR_log2_CPM_plus_1, one$LH_factor, mean)
    total_ss <- sum((one$VDR_log2_CPM_plus_1 - mean(one$VDR_log2_CPM_plus_1))^2)
    fitted_ss <- sum((fitted(fit) - mean(one$VDR_log2_CPM_plus_1))^2)
    temporal_VDR_rows[[length(temporal_VDR_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, donor_n = nrow(one), 
        LH3_mean = unname(group_means["3"]), LH5_mean = unname(group_means["5"]), LH7_mean = unname(group_means["7"]), 
        LH9_mean = unname(group_means["9"]), LH11_mean = unname(group_means["11"]), permutation_F = permutation$F, 
        permutation_p = permutation$p, permutation_n = permutation$permutation_n, classical_F = anova_table$`F value`[1], 
        classical_p = anova_table$`Pr(>F)`[1], eta_squared = if (is.finite(total_ss) && total_ss > 0) {
            fitted_ss/total_ss
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
    summary_one <- group_summary(data = one, value_column = "VDR_log2_CPM_plus_1", group_columns = c("broad_cell_type", 
        "LH_day"))
    temporal_VDR_summary_rows[[length(temporal_VDR_summary_rows) + 1L]] <- as.data.frame(summary_one)
    for (pair_index in seq_len(nrow(planned_pairs))) {
        pair <- planned_pairs[pair_index, ]
        pair_data <- one[one$LH_day %in% c(pair$reference_day, pair$comparison_day), , drop = FALSE]
        pair_group <- paste0("LH", pair_data$LH_day)
        reference_level <- paste0("LH", pair$reference_day)
        comparison_level <- paste0("LH", pair$comparison_day)
        exact_test <- exact_two_group_test(y = pair_data$VDR_log2_CPM_plus_1, group = pair_group, reference_level = reference_level, 
            comparison_level = comparison_level, seed = 20260805L + 11000L + type_index * 100L + pair_index)
        welch <- welch_summary(y = pair_data$VDR_log2_CPM_plus_1, group = pair_group, reference_level = reference_level, 
            comparison_level = comparison_level)
        temporal_VDR_pair_rows[[length(temporal_VDR_pair_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, 
            comparison = pair$comparison, reference_day = pair$reference_day, comparison_day = pair$comparison_day, 
            mean_difference = exact_test$estimate, exact_permutation_p = exact_test$p, exact_test_mode = exact_test$test_mode, 
            exact_permutation_n = exact_test$permutation_n, Welch_ci_low = welch$ci_low, Welch_ci_high = welch$ci_high, 
            Welch_p = welch$p, Hedges_g = hedges_g(y = pair_data$VDR_log2_CPM_plus_1, group = pair_group, reference_level = reference_level, 
                comparison_level = comparison_level), stringsAsFactors = FALSE)
    }
}

temporal_VDR <- do.call(rbind, temporal_VDR_rows)

temporal_VDR$FDR_BH <- stats::p.adjust(temporal_VDR$permutation_p, method = "BH")

temporal_VDR_pairwise <- do.call(rbind, temporal_VDR_pair_rows)

temporal_VDR_pairwise$FDR_BH_within_comparison <- ave(temporal_VDR_pairwise$exact_permutation_p, temporal_VDR_pairwise$comparison, 
    FUN = function(p) {
        stats::p.adjust(p, method = "BH")
    })

temporal_VDR_summary <- do.call(rbind, temporal_VDR_summary_rows)

write_csv_utf8(temporal_VDR, file.path(table_root, "05_temporal_VDR_omnibus.csv"))

write_csv_utf8(temporal_VDR_pairwise, file.path(table_root, "06_temporal_VDR_pairwise.csv"))

write_csv_utf8(temporal_VDR_summary, file.path(source_data_root, "source_temporal_VDR_summary.csv"))

LH7_VDR_primary_rows <- list()

LH7_VDR_sensitivity_rows <- list()

for (type_index in seq_along(LH7_types)) {
    cell_type <- LH7_types[type_index]
    one <- LH7_VDR_data[LH7_VDR_data$broad_cell_type == cell_type, , drop = FALSE]
    one$clinical_group <- factor(one$clinical_group, levels = c("FERTILE", "RIF"))
    exact_test <- exact_two_group_test(y = one$VDR_log2_CPM_plus_1, group = one$clinical_group, reference_level = "FERTILE", 
        comparison_level = "RIF", seed = 20260805L + 12000L + type_index)
    welch <- welch_summary(y = one$VDR_log2_CPM_plus_1, group = one$clinical_group, reference_level = "FERTILE", 
        comparison_level = "RIF")
    loo <- leave_one_out_difference(y = one$VDR_log2_CPM_plus_1, group = one$clinical_group, reference_level = "FERTILE", 
        comparison_level = "RIF")
    LH7_VDR_primary_rows[[length(LH7_VDR_primary_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, fertile_n = sum(one$clinical_group == 
        "FERTILE"), RIF_n = sum(one$clinical_group == "RIF"), fertile_mean = mean(one$VDR_log2_CPM_plus_1[one$clinical_group == 
        "FERTILE"]), RIF_mean = mean(one$VDR_log2_CPM_plus_1[one$clinical_group == "RIF"]), RIF_minus_fertile = exact_test$estimate, 
        exact_permutation_p = exact_test$p, exact_test_mode = exact_test$test_mode, exact_permutation_n = exact_test$permutation_n, 
        Welch_ci_low = welch$ci_low, Welch_ci_high = welch$ci_high, Welch_p = welch$p, Hedges_g = hedges_g(y = one$VDR_log2_CPM_plus_1, 
            group = one$clinical_group, reference_level = "FERTILE", comparison_level = "RIF"), LOO_min = unname(loo["loo_min"]), 
        LOO_max = unname(loo["loo_max"]), LOO_sign_consistency = unname(loo["loo_sign_consistency"]), stringsAsFactors = FALSE)
    for (fertile_subgroup in c("LH7_STANDARD", "LH7_PGT")) {
        sensitivity <- one[one$clinical_group == "RIF" | (one$clinical_group == "FERTILE" & one$LH7_fertile_subgroup == 
            fertile_subgroup), , drop = FALSE]
        sensitivity_group <- ifelse(sensitivity$clinical_group == "RIF", "RIF", fertile_subgroup)
        exact_sensitivity <- exact_two_group_test(y = sensitivity$VDR_log2_CPM_plus_1, group = sensitivity_group, 
            reference_level = fertile_subgroup, comparison_level = "RIF", seed = 20260805L + 13000L + type_index * 
                100L + match(fertile_subgroup, c("LH7_STANDARD", "LH7_PGT")))
        welch_sensitivity <- welch_summary(y = sensitivity$VDR_log2_CPM_plus_1, group = sensitivity_group, reference_level = fertile_subgroup, 
            comparison_level = "RIF")
        LH7_VDR_sensitivity_rows[[length(LH7_VDR_sensitivity_rows) + 1L]] <- data.frame(broad_cell_type = cell_type, 
            sensitivity_reference = fertile_subgroup, fertile_subgroup_n = sum(sensitivity_group == fertile_subgroup), 
            RIF_n = sum(sensitivity_group == "RIF"), RIF_minus_reference = exact_sensitivity$estimate, exact_permutation_p = exact_sensitivity$p, 
            exact_test_mode = exact_sensitivity$test_mode, exact_permutation_n = exact_sensitivity$permutation_n, 
            Welch_ci_low = welch_sensitivity$ci_low, Welch_ci_high = welch_sensitivity$ci_high, Welch_p = welch_sensitivity$p, 
            Hedges_g = hedges_g(y = sensitivity$VDR_log2_CPM_plus_1, group = sensitivity_group, reference_level = fertile_subgroup, 
                comparison_level = "RIF"), stringsAsFactors = FALSE)
    }
}

LH7_VDR_primary <- do.call(rbind, LH7_VDR_primary_rows)

LH7_VDR_primary$FDR_BH <- stats::p.adjust(LH7_VDR_primary$exact_permutation_p, method = "BH")

LH7_VDR_primary$effect_direction <- sign_label(LH7_VDR_primary$RIF_minus_fertile)

LH7_VDR_primary$LOO_stable_sign <- (LH7_VDR_primary$LOO_sign_consistency >= 0.80000000000000004)

LH7_VDR_sensitivity <- do.call(rbind, LH7_VDR_sensitivity_rows)

LH7_VDR_sensitivity$FDR_BH_within_reference <- ave(LH7_VDR_sensitivity$exact_permutation_p, LH7_VDR_sensitivity$sensitivity_reference, 
    FUN = function(p) {
        stats::p.adjust(p, method = "BH")
    })

write_csv_utf8(LH7_VDR_primary, file.path(table_root, "07_LH7_VDR_primary.csv"))

write_csv_utf8(LH7_VDR_sensitivity, file.path(table_root, "08_LH7_VDR_sensitivity.csv"))

significant_temporal_program <- temporal_program[temporal_program$FDR_BH < 0.050000000000000003, , drop = FALSE]

significant_LH7_program <- LH7_primary[LH7_primary$FDR_BH < 0.050000000000000003, , drop = FALSE]

significant_temporal_VDR <- temporal_VDR[temporal_VDR$FDR_BH < 0.050000000000000003, , drop = FALSE]

significant_LH7_VDR <- LH7_VDR_primary[LH7_VDR_primary$FDR_BH < 0.050000000000000003, , drop = FALSE]

cross_model_summary <- data.frame(model_family = c("TEMPORAL_VDR_PROGRAM", "LH7_RIF_VS_FERTILE_VDR_PROGRAM", "TEMPORAL_VDR_EXPRESSION", 
    "LH7_RIF_VS_FERTILE_VDR_EXPRESSION"), tested_hypothesis_n = c(nrow(temporal_program), nrow(LH7_primary), nrow(temporal_VDR), 
    nrow(LH7_VDR_primary)), FDR_lt_0_05_n = c(nrow(significant_temporal_program), nrow(significant_LH7_program), 
    nrow(significant_temporal_VDR), nrow(significant_LH7_VDR)), minimum_raw_p = c(min(temporal_program$permutation_p, 
    na.rm = TRUE), min(LH7_primary$exact_permutation_p, na.rm = TRUE), min(temporal_VDR$permutation_p, na.rm = TRUE), 
    min(LH7_VDR_primary$exact_permutation_p, na.rm = TRUE)), minimum_FDR = c(min(temporal_program$FDR_BH, na.rm = TRUE), 
    min(LH7_primary$FDR_BH, na.rm = TRUE), min(temporal_VDR$FDR_BH, na.rm = TRUE), min(LH7_VDR_primary$FDR_BH, 
        na.rm = TRUE)), stringsAsFactors = FALSE)

write_csv_utf8(cross_model_summary, file.path(table_root, "09_cross_model_summary.csv"))

theme_qc <- ggplot2::theme_bw(base_size = 9) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95"))

temporal_program_summary$LH_day <- factor(temporal_program_summary$LH_day, levels = c(3, 5, 7, 9, 11), labels = c("LH+3", 
    "LH+5", "LH+7", "LH+9", "LH+11"))

plot_temporal_program <- ggplot2::ggplot(temporal_program_summary, ggplot2::aes(x = LH_day, y = mean, group = 1)) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, ymax = ci_high, group = 1), alpha = 0.14999999999999999) + 
    ggplot2::geom_line(linewidth = 0.55000000000000004) + ggplot2::geom_point(size = 1.3) + ggplot2::facet_grid(program_name ~ 
    broad_cell_type, scales = "free_y") + ggplot2::labs(title = "Fertile temporal trajectories of locked VDR programs", 
    subtitle = "Donor-level pseudobulk scores; ribbons show 95% mean confidence intervals", x = "Biopsy timing", 
    y = "Signed VDR-program score") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1), strip.text = ggplot2::element_text(size = 6))

save_plot_pair(plot_temporal_program, "QC1_temporal_programs", width = 16, height = 11)

LH7_primary$FDR_label <- ifelse(LH7_primary$FDR_BH < 0.050000000000000003, "FDR < 0.05", "FDR >= 0.05")

plot_LH7_program <- ggplot2::ggplot(LH7_primary, ggplot2::aes(x = RIF_minus_fertile, y = broad_cell_type)) + ggplot2::geom_vline(xintercept = 0, 
    linetype = "dashed") + ggplot2::geom_segment(ggplot2::aes(x = Welch_ci_low, xend = Welch_ci_high, y = broad_cell_type, 
    yend = broad_cell_type), linewidth = 0.5) + ggplot2::geom_point(ggplot2::aes(shape = FDR_label, size = abs(Hedges_g))) + 
    ggplot2::facet_wrap(~program_name, scales = "free_x", ncol = 2) + ggplot2::labs(title = "LH+7 RIF versus fertile VDR-program effects", 
    subtitle = "Positive values indicate greater alignment in RIF; intervals are Welch 95% CIs", x = "RIF minus fertile score", 
    y = "Broad cell type", shape = "Permutation FDR", size = "|Hedges g|") + theme_qc

save_plot_pair(plot_LH7_program, "QC2_LH7_program_forest", width = 11, height = 9)

temporal_VDR_summary$LH_day <- factor(temporal_VDR_summary$LH_day, levels = c(3, 5, 7, 9, 11), labels = c("LH+3", 
    "LH+5", "LH+7", "LH+9", "LH+11"))

plot_temporal_VDR <- ggplot2::ggplot(temporal_VDR_summary, ggplot2::aes(x = LH_day, y = mean, group = 1)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_low, 
    ymax = ci_high, group = 1), alpha = 0.14999999999999999) + ggplot2::geom_line(linewidth = 0.65000000000000002) + 
    ggplot2::geom_point(size = 1.5) + ggplot2::facet_wrap(~broad_cell_type, scales = "free_y") + ggplot2::labs(title = "Fertile temporal trajectory of VDR expression", 
    subtitle = "Donor-level pseudobulk log2(CPM + 1)", x = "Biopsy timing", y = "VDR log2(CPM + 1)") + theme_qc + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

save_plot_pair(plot_temporal_VDR, "QC3_temporal_VDR", width = 11, height = 7)

LH7_VDR_primary$FDR_label <- ifelse(LH7_VDR_primary$FDR_BH < 0.050000000000000003, "FDR < 0.05", "FDR >= 0.05")

plot_LH7_VDR <- ggplot2::ggplot(LH7_VDR_primary, ggplot2::aes(x = RIF_minus_fertile, y = broad_cell_type)) + ggplot2::geom_vline(xintercept = 0, 
    linetype = "dashed") + ggplot2::geom_segment(ggplot2::aes(x = Welch_ci_low, xend = Welch_ci_high, y = broad_cell_type, 
    yend = broad_cell_type), linewidth = 0.59999999999999998) + ggplot2::geom_point(ggplot2::aes(shape = FDR_label, 
    size = abs(Hedges_g))) + ggplot2::labs(title = "LH+7 RIF versus fertile VDR expression", subtitle = "Positive values indicate higher VDR expression in RIF", 
    x = "RIF minus fertile log2(CPM + 1)", y = "Broad cell type", shape = "Permutation FDR", size = "|Hedges g|") + 
    theme_qc

save_plot_pair(plot_LH7_VDR, "QC4_LH7_VDR_forest", width = 9, height = 6)

sensitivity_plot_data <- LH7_sensitivity

sensitivity_plot_data$comparison_label <- ifelse(sensitivity_plot_data$sensitivity_reference == "LH7_STANDARD", 
    "RIF vs standard fertile", "RIF vs PGT fertile")

plot_sensitivity <- ggplot2::ggplot(sensitivity_plot_data, ggplot2::aes(x = primary_RIF_minus_fertile, y = RIF_minus_reference, 
    shape = comparison_label)) + ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + ggplot2::geom_vline(xintercept = 0, 
    linetype = "dashed") + ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted") + ggplot2::geom_point(size = 2) + 
    ggplot2::facet_wrap(~program_name, scales = "free") + ggplot2::labs(title = "LH+7 fertile-subgroup sensitivity", 
    subtitle = "Concordance of primary and subgroup-specific RIF effects", x = "Primary: RIF minus all fertile", 
    y = "Sensitivity: RIF minus fertile subgroup", shape = "Sensitivity comparison") + theme_qc

save_plot_pair(plot_sensitivity, "QC5_LH7_sensitivity", width = 10, height = 8)

gate_checks <- data.frame(check = c("upstream_gate_authorizes_donor_models", "pseudobulk_status_complete", "temporal_type_n_is_5", 
    "LH7_type_n_is_6", "four_locked_programs_tested", "temporal_program_model_n_is_20", "LH7_program_model_n_is_24", 
    "temporal_VDR_model_n_is_5", "LH7_VDR_model_n_is_6", "all_temporal_program_models_have_permutation_p", "all_LH7_program_models_have_exact_p", 
    "all_primary_LH7_models_have_leave_one_out_QC", "all_FDR_values_generated", "composition_inference_remains_disabled", 
    "subtype_inference_remains_disabled"), pass = c(pb_gate$donor_level_modeling_authorized[[1]], pb_status$status[[1]] == 
    "COMPLETE", length(temporal_types) == 5L, length(LH7_types) == 6L, setequal(program_names, expected_programs), 
    nrow(temporal_program) == 20L, nrow(LH7_primary) == 24L, nrow(temporal_VDR) == 5L, nrow(LH7_VDR_primary) == 
        6L, all(is.finite(temporal_program$permutation_p)), all(is.finite(LH7_primary$exact_permutation_p)), all(is.finite(LH7_primary$LOO_sign_consistency)), 
    all(is.finite(c(temporal_program$FDR_BH, LH7_primary$FDR_BH, temporal_VDR$FDR_BH, LH7_VDR_primary$FDR_BH))), 
    !pb_gate$cell_composition_inference_authorized[[1]], !pb_gate$subtype_modeling_authorized[[1]]), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "02_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE"
} else {
    "HOLD_DONOR_LEVEL_VDR_MODEL_QC_FAILED"
}

model_gate <- data.frame(dataset_id = dataset_id, input_version = input_version, pseudobulk_version = pseudobulk_version, 
    model_version = model_version, gate_status = gate_status, temporal_cell_type_n = length(temporal_types), LH7_cell_type_n = length(LH7_types), 
    program_n = length(program_names), temporal_program_model_n = nrow(temporal_program), temporal_program_FDR_lt_0_05_n = nrow(significant_temporal_program), 
    LH7_program_model_n = nrow(LH7_primary), LH7_program_FDR_lt_0_05_n = nrow(significant_LH7_program), temporal_VDR_model_n = nrow(temporal_VDR), 
    temporal_VDR_FDR_lt_0_05_n = nrow(significant_temporal_VDR), LH7_VDR_model_n = nrow(LH7_VDR_primary), LH7_VDR_FDR_lt_0_05_n = nrow(significant_LH7_VDR), 
    gene_level_pseudobulk_DE_authorized = all_gate_checks_pass, composition_inference_authorized = FALSE, subtype_inference_authorized = FALSE, 
    result_lock_authorized = FALSE, stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "VDR_MODEL_GATE.csv")

write_csv_utf8(model_gate, gate_path)

card_lines <- c("# STEP04E1G Donor-level VDR Model Card", "", paste0("- Dataset: ", dataset_id), paste0("- Pseudobulk version: ", 
    pseudobulk_version), paste0("- Model version: ", model_version), paste0("- Gate: ", gate_status), paste0("- Temporal program tests: ", 
    nrow(temporal_program)), paste0("- Temporal program FDR < 0.05: ", nrow(significant_temporal_program)), paste0("- LH+7 program tests: ", 
    nrow(LH7_primary)), paste0("- LH+7 program FDR < 0.05: ", nrow(significant_LH7_program)), paste0("- Temporal VDR-expression FDR < 0.05: ", 
    nrow(significant_temporal_VDR)), paste0("- LH+7 VDR-expression FDR < 0.05: ", nrow(significant_LH7_VDR)), paste0("- Gene-level pseudobulk DE authorized: ", 
    model_gate$gene_level_pseudobulk_DE_authorized), "- Cell-composition inference authorized: FALSE", "- Subtype inference authorized: FALSE", 
    "- Result lock authorized: FALSE")

card_path <- file.path(run_qc_root, "VDR_MODEL_RESULT_CARD.md")

write_text(card_lines, card_path)

interpretation_lines <- c("# STEP04E1G donor-level VDR model interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), "", "## Model scope", "", "- Each observation represents one donor-by-broad-cell-type pseudobulk.", 
    "- Fertile temporal models compare cross-sectional LH+3, LH+5, LH+7, LH+9 and LH+11 donor groups.", "- The primary LH+7 comparison contrasts RIF with all fertile LH+7 donors.", 
    "- Separate standard-fertile and PGT-fertile comparisons are sensitivity analyses.", "- Program scores indicate relative alignment with the locked VDR interaction direction and are not causal activation scores.", 
    "", "## Statistical scope", "", "- Permutation tests provide the primary p-values.", "- Benjamini-Hochberg correction is applied within predefined model families.", 
    "- Welch confidence intervals and Hedges g are retained as secondary effect-size diagnostics.", "- Leave-one-out sign stability is reported for primary LH+7 comparisons.", 
    "", "## Next step", "", "- Review temporal and LH+7 VDR-program effects.", "- Proceed to gene-level pseudobulk differential expression only after this Gate passes.", 
    "- Integrate gene-level effects with the locked 436/284/46/84 candidate hierarchy.", "- Do not perform cell-composition or subtype-level inference.")

interpretation_path <- file.path(run_qc_root, "VDR_MODEL_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

methods_lines <- c("# STEP04E1G Methods draft", "", paste("Donor-level analyses were performed using one pseudobulk profile per", 
    "donor and broad cell type."), paste("Four frozen VDR RNA programs were examined: the full 436-gene", "interaction program, the 284-gene occupancy-supported program, the", 
    "46-gene shared-occupancy subset and the 84-gene promoter-occupancy", "subset."), paste("Fertile temporal trajectories were evaluated across LH+3, LH+5, LH+7,", 
    "LH+9 and LH+11 using the actual pseudobulk profiles that passed the", "predefined cell-count threshold."), 
    paste("The temporal design was allowed to be unbalanced provided that each", "time point retained at least two independent donors; specifically,", 
        "ciliated epithelial cells included two eligible LH+11 donors whereas", "the remaining eligible temporal cell types retained three LH+11 donors."), 
    paste("Temporal omnibus evidence was assessed using a 5,000-label permutation", "one-way analysis of variance, with classical one-way ANOVA retained as", 
        "a secondary diagnostic."), paste("At LH+7, RIF and fertile donors were compared using an exact", "two-group permutation test of the mean difference whenever the number", 
        "of possible allocations was computationally tractable."), paste("Welch confidence intervals, Hedges standardized effect sizes and", 
        "leave-one-out sign stability were reported as complementary", "diagnostics."), paste("Sensitivity analyses compared RIF donors separately with standard", 
        "fertile LH+7 donors and PGT-associated fertile LH+7 donors."), paste("Benjamini-Hochberg adjustment was applied within each predefined", 
        "model family. Parallel analyses were performed for VDR expression", "measured as log2 counts per million plus one."))

write_text(methods_lines, file.path(manuscript_text_root, "STEP04E1G_VDR_models_Methods_draft.md"))

results_shell <- c("# STEP04E1G Results shell", "", paste0("The donor-level VDR-program analysis included ", length(temporal_types), 
    " broad cell types in the fertile temporal series and ", length(LH7_types), " broad cell types in the LH+7 RIF comparison."), 
    paste0("Across ", nrow(temporal_program), " temporal program tests, ", nrow(significant_temporal_program), 
        " met FDR < 0.05."), paste0("Across ", nrow(LH7_primary), " primary LH+7 program comparisons, ", nrow(significant_LH7_program), 
        " met FDR < 0.05."), paste0("For VDR expression, ", nrow(significant_temporal_VDR), " temporal cell-type model(s) and ", 
        nrow(significant_LH7_VDR), " LH+7 RIF-versus-fertile model(s) met FDR < 0.05."), "", "Detailed biological interpretation will be drafted after review of effect directions, sensitivity analyses and gene-level pseudobulk results.")

write_text(results_shell, file.path(manuscript_text_root, "STEP04E1G_VDR_models_Results_shell.md"))

run_end <- Sys.time()

package_names <- c("data.table", "ggplot2", "jsonlite", "scales")

package_table <- data.frame(package = package_names, version = vapply(package_names, safe_package_version, character(1)), 
    stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- unique(c(list.files(output_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1G_VDR_models_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1G_VDR_models_Results_shell.md")))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, input_version = input_version, 
    pseudobulk_version = pseudobulk_version, model_version = model_version, dataset_id = dataset_id, started_at = as.character(run_start), 
    completed_at = as.character(run_end), elapsed_seconds = round(as.numeric(difftime(run_end, run_start, units = "secs")), 
        2), gate_status = gate_status, temporal_program_model_n = nrow(temporal_program), LH7_program_model_n = nrow(LH7_primary), 
    temporal_VDR_model_n = nrow(temporal_VDR), LH7_VDR_model_n = nrow(LH7_VDR_primary), donor_level_modeling_completed = all_gate_checks_pass, 
    gene_level_pseudobulk_DE_authorized = model_gate$gene_level_pseudobulk_DE_authorized, composition_inference_authorized = FALSE, 
    subtype_inference_authorized = FALSE, result_lock_authorized = FALSE, output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "temporal_program_FDR_lt_0_05_n=", nrow(significant_temporal_program), 
    "LH7_program_FDR_lt_0_05_n=", nrow(significant_LH7_program))

small_files <- unique(c(list.files(output_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, recursive = TRUE, 
    full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1G_VDR_models_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1G_VDR_models_Results_shell.md"), 
    normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & !dir.exists(small_files) & file.info(small_files)$size <= 
    200 * 1024^2]

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

cat("STEP04E1G GSE250130 donor-level VDR models completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Temporal program models: ", nrow(temporal_program), "\n", sep = "")

cat("Temporal program FDR < 0.05: ", nrow(significant_temporal_program), "\n", sep = "")

cat("LH7 program models: ", nrow(LH7_primary), "\n", sep = "")

cat("LH7 program FDR < 0.05: ", nrow(significant_LH7_program), "\n", sep = "")

cat("Temporal VDR-expression FDR < 0.05: ", nrow(significant_temporal_VDR), "\n", sep = "")

cat("LH7 VDR-expression FDR < 0.05: ", nrow(significant_LH7_VDR), "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe model and QC folders are complete.\n", sep = "")
}

