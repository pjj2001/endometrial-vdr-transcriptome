#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C08: Fit gene-level pseudobulk differential-expression models

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260805)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "fit_gse250130_gene_models.R"

script_version <- "STEP04E1H_GSE250130_GENE_PSEUDOBULK_DE_v1.2_PATHSAFE_RESUME"

input_version <- "SCI_v1_1_memory_safe"

pseudobulk_version <- "PBVDR_v1_0"

upstream_model_version <- "VDRMODELS_v1_1"

DE_version <- "DE_v1_2"

run_start <- Sys.time()

run_id <- paste0("STEP04E1H_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

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

as_logical_safe <- function(x) {
    if (is.logical(x)) 
        return(x)
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

slugify <- function(x) {
    out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
    out <- gsub("^_+|_+$", "", out)
    substr(out, 1L, 60L)
}

cell_type_code <- function(cell_type) {
    code_map <- c(`Ciliated epithelial cells` = "CE", `Endothelial cells` = "ENDO", `Myeloid cells` = "MY", `NK/T cells` = "NKT", 
        `Stromal cells` = "STR", `Unciliated epithelial cells` = "UE")
    code <- unname(code_map[as.character(cell_type)])
    if (length(code) != 1L || is.na(code) || !nzchar(code)) {
        stop("No short cell-type code defined for: ", cell_type)
    }
    code
}

contrast_code <- function(contrast_name) {
    code_map <- c(LH5_vs_LH3 = "5v3", LH7_vs_LH5 = "7v5", LH9_vs_LH7 = "9v7", LH11_vs_LH9 = "11v9", LH7_vs_LH3 = "7v3", 
        LH11_vs_LH3 = "11v3", LH11_vs_LH7 = "11v7", RIF_vs_LH7_STANDARD = "STD", RIF_vs_LH7_PGT = "PGT")
    code <- unname(code_map[as.character(contrast_name)])
    if (length(code) != 1L || is.na(code) || !nzchar(code)) {
        stop("No short contrast code defined for: ", contrast_name)
    }
    code
}

sign_label <- function(x) {
    ifelse(is.na(x), "NA", ifelse(x > 0, "UP", ifelse(x < 0, "DOWN", "ZERO")))
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

log_path <- NULL

log_event <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(line, "\n")
    if (!is.null(log_path)) {
        cat(line, "\n", file = log_path, append = TRUE)
    }
    flush.console()
}

required_CRAN_packages <- c("Matrix", "data.table", "ggplot2", "jsonlite", "scales")

required_Bioc_packages <- c("edgeR", "limma")

missing_CRAN <- required_CRAN_packages[!vapply(required_CRAN_packages, requireNamespace, logical(1), quietly = TRUE)]

missing_Bioc <- required_Bioc_packages[!vapply(required_Bioc_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_CRAN) || length(missing_Bioc)) {
    message_lines <- c("Required package(s) are missing.")
    if (length(missing_CRAN)) {
        message_lines <- c(message_lines, paste0("Missing CRAN: ", paste(missing_CRAN, collapse = ", ")), paste0("Run: install.packages(c(", 
            paste(paste0("\"", missing_CRAN, "\""), collapse = ", "), "))"))
    }
    if (length(missing_Bioc)) {
        message_lines <- c(message_lines, paste0("Missing Bioconductor: ", paste(missing_Bioc, collapse = ", ")), 
            "Run:", "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) install.packages(\"BiocManager\")", 
            paste0("BiocManager::install(c(", paste(paste0("\"", missing_Bioc, "\""), collapse = ", "), "), ask = FALSE, update = FALSE)"))
    }
    stop(paste(message_lines, collapse = "\n"))
}

if (!exists("fry", envir = asNamespace("limma"), inherits = FALSE) || !is.function(get("fry", envir = asNamespace("limma"), 
    inherits = FALSE))) {
    stop("limma::fry is unavailable in the installed limma namespace. ", "Reinstall limma with BiocManager, restart R, and rerun this script.")
}

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

execution_root <- file.path(step04_root, "05_execution")

results_root <- file.path(step04_root, "06_results")

qc_base_root <- file.path(step04_root, "08_qc")

parameter_root <- file.path(step04_root, "03_parameters")

environment_root <- file.path(step04_root, "04_environment")

manuscript_text_root <- file.path(step04_root, "09_manuscript_text")

pb_root <- file.path(results_root, "S4E1_sc", input_version, "pb", dataset_id, pseudobulk_version)

pb_counts_path <- file.path(pb_root, "objects", "pb_raw_counts.rds")

pb_status_path <- file.path(pb_root, "objects", "pb_status.csv")

pb_metadata_path <- file.path(pb_root, "tables", "02_pb_metadata.csv")

eligibility_path <- file.path(pb_root, "tables", "04_model_eligibility.csv")

candidate_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

output_root <- file.path(results_root, "S4E1_sc", input_version, "de", dataset_id, DE_version)

previous_partial_root <- file.path(results_root, "S4E1_sc", input_version, "de", dataset_id, "GENEDE_v1_1")

previous_per_model_root <- file.path(previous_partial_root, "per_model")

per_model_root <- file.path(output_root, "m")

full_result_root <- file.path(output_root, "full")

table_root <- file.path(output_root, "tab")

candidate_root <- file.path(output_root, "cand")

source_data_root <- file.path(output_root, "src")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1H_", format(run_start, "%Y%m%d_%H%M%S"), "_deQC"))

figure_root <- file.path(run_qc_root, "figures")

for (p in c(execution_root, parameter_root, environment_root, manuscript_text_root, output_root, per_model_root, 
    full_result_root, table_root, candidate_root, source_data_root, run_qc_root, figure_root)) {
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

planned_paths <- data.frame(object = c("output_root", "per_model_root", "full_result_root", "table_root", "candidate_root", 
    "source_data_root", "run_qc_root", "figure_root", "log_path", "run_manifest_path", "parameter_path", "session_path", 
    "handoff_zip"), path = normalize_slashes(c(output_root, per_model_root, full_result_root, table_root, candidate_root, 
    source_data_root, run_qc_root, figure_root, log_path, run_manifest_path, parameter_path, session_path, handoff_zip)), 
    stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

if (.Platform$OS.type == "windows" && any(planned_paths$character_n > 245L)) {
    bad_paths <- planned_paths[planned_paths$character_n > 245L, , drop = FALSE]
    stop("One or more planned output paths are too long:\n", paste(paste0(bad_paths$character_n, " chars: ", bad_paths$path), 
        collapse = "\n"))
}

required_inputs <- c(pb_counts_path, pb_status_path, pb_metadata_path, eligibility_path, candidate_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Missing required input file(s):\n", paste(missing_inputs, collapse = "\n"))
}

upstream_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

upstream_qc_dirs <- upstream_qc_dirs[grepl("^STEP04E1G_[0-9]{8}_[0-9]{6}_modelQC$", basename(upstream_qc_dirs), 
    perl = TRUE)]

if (!length(upstream_qc_dirs)) {
    stop("No completed STEP04E1G model-QC folder found.")
}

latest_upstream_qc <- upstream_qc_dirs[which.max(file.info(upstream_qc_dirs)$mtime)]

upstream_gate_path <- file.path(latest_upstream_qc, "VDR_MODEL_GATE.csv")

upstream_gate <- safe_read_csv(upstream_gate_path)

pb_status <- safe_read_csv(pb_status_path)

pb_metadata <- safe_read_csv(pb_metadata_path)

eligibility <- safe_read_csv(eligibility_path)

candidate_table <- safe_read_csv(candidate_path)

if (is.null(upstream_gate) || nrow(upstream_gate) != 1L || upstream_gate$gate_status[[1]] != "PASS_DONOR_LEVEL_VDR_MODELS_COMPLETE" || 
    !isTRUE(upstream_gate$gene_level_pseudobulk_DE_authorized[[1]])) {
    stop("STEP04E1G Gate does not authorize gene-level pseudobulk DE.")
}

if (is.null(pb_status) || nrow(pb_status) != 1L || pb_status$status[[1]] != "COMPLETE" || pb_status$gene_n[[1]] != 
    36601L || pb_status$pseudobulk_n[[1]] != 222L || pb_status$aggregated_cell_n[[1]] != 236636L) {
    stop("Pseudobulk status failed identity validation.")
}

required_metadata_columns <- c("pb_id", "sample_accession", "donor_id", "broad_cell_type", "clinical_group", "LH_day", 
    "temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "LH7_fertile_subgroup", "minimum_cell_threshold_pass", 
    "cell_n", "library_size")

missing_metadata_columns <- setdiff(required_metadata_columns, colnames(pb_metadata))

if (length(missing_metadata_columns)) {
    stop("Pseudobulk metadata is missing columns: ", paste(missing_metadata_columns, collapse = ", "))
}

required_candidate_columns <- c("feature_id", "interaction_logFC", "occupancy_supported", "shared_occupancy_supported", 
    "any_promoter_overlap_support")

missing_candidate_columns <- setdiff(required_candidate_columns, colnames(candidate_table))

if (length(missing_candidate_columns)) {
    stop("Candidate table is missing columns: ", paste(missing_candidate_columns, collapse = ", "))
}

if (nrow(pb_metadata) != 222L || nrow(eligibility) != 8L || nrow(candidate_table) != 436L) {
    stop("Input table row-count identity checks failed.")
}

candidate_table$feature_id <- as.character(candidate_table$feature_id)

candidate_table$interaction_logFC <- safe_numeric(candidate_table$interaction_logFC)

candidate_table$occupancy_supported <- as_logical_safe(candidate_table$occupancy_supported)

candidate_table$shared_occupancy_supported <- as_logical_safe(candidate_table$shared_occupancy_supported)

candidate_table$any_promoter_overlap_support <- as_logical_safe(candidate_table$any_promoter_overlap_support)

if (anyDuplicated(candidate_table$feature_id) > 0L || sum(candidate_table$occupancy_supported) != 284L || sum(candidate_table$shared_occupancy_supported) != 
    46L || sum(candidate_table$any_promoter_overlap_support) != 84L || any(!is.finite(candidate_table$interaction_logFC)) || 
    any(candidate_table$interaction_logFC == 0)) {
    stop("Locked VDR candidate hierarchy failed validation.")
}

pb_metadata$LH_day <- safe_numeric(pb_metadata$LH_day)

logical_columns <- c("temporal_fertile_eligible", "LH7_RIF_comparison_eligible", "minimum_cell_threshold_pass")

for (column_name in logical_columns) {
    pb_metadata[[column_name]] <- as_logical_safe(pb_metadata[[column_name]])
}

temporal_types <- as.character(eligibility$broad_cell_type[as_logical_safe(eligibility$fertile_temporal_model_eligible)])

LH7_types <- as.character(eligibility$broad_cell_type[as_logical_safe(eligibility$LH7_RIF_vs_fertile_model_eligible)])

if (length(temporal_types) != 5L || length(LH7_types) != 6L) {
    stop("Eligible cell-type counts do not match the upstream Gate.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, input_version = input_version, 
    pseudobulk_version = pseudobulk_version, upstream_model_version = upstream_model_version, DE_version = DE_version, 
    supersedes_partial_DE_versions = c("GENEDE_v1_0_EDGE_R_NAMESPACE_FRY_FAILURE", "GENEDE_v1_1_DYNAMIC_PATH_LENGTH_FAILURE"), 
    previous_partial_root = normalize_slashes(previous_partial_root), resume_policy = "IMPORT_COMPLETE_PER_MODEL_RDS_FROM_GENEDE_v1_1", 
    dataset_id = dataset_id, primary_gene_level_FDR = "GENOME_WIDE_BH_WITHIN_EACH_CELL_TYPE_AND_CONTRAST", focused_candidate_FDR = "SEPARATE_BH_WITHIN_PROSPECTIVELY_LOCKED_436_284_46_84_HIERARCHIES", 
    reporting_effect_threshold_pairwise = "ABS_LOG2_FOLD_CHANGE_GE_0.5", normalization = "TMM", expression_filter = "edgeR_filterByExpr_DEFAULT_WITH_MODEL_DESIGN", 
    model = "ROBUST_EDGE_R_QUASI_LIKELIHOOD_GLM", temporal_design = "CROSS_SECTIONAL_FACTOR_LH3_LH5_LH7_LH9_LH11_UNBALANCED_ALLOWED", 
    temporal_omnibus = "MULTI_DF_LH5_LH7_LH9_LH11_VERSUS_LH3", LH7_primary = "RIF_VERSUS_ALL_FERTILE_LH7", LH7_sensitivity = c("RIF_VERSUS_LH7_STANDARD", 
        "RIF_VERSUS_LH7_PGT"), locked_program_test = "DIRECTIONALLY_WEIGHTED_EDGE_R_FRY_FOR_ONE_DF_CONTRASTS", 
    temporal_cell_types = temporal_types, LH7_cell_types = LH7_types, composition_inference_authorized = FALSE, 
    subtype_inference_authorized = FALSE, random_seed = 20260805L, result_lock_authorized = FALSE)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

log_event("Starting", run_id, "temporal_types=", collapse_unique(temporal_types), "LH7_types=", collapse_unique(LH7_types))

pb_counts <- readRDS(pb_counts_path)

if (!inherits(pb_counts, "Matrix") && !is.matrix(pb_counts)) {
    stop("pb_raw_counts.rds is not a matrix-like object.")
}

if (nrow(pb_counts) != 36601L || ncol(pb_counts) != 222L || anyDuplicated(rownames(pb_counts)) > 0L || anyDuplicated(colnames(pb_counts)) > 
    0L) {
    stop("Pseudobulk count matrix identity checks failed.")
}

pb_metadata <- pb_metadata[match(colnames(pb_counts), pb_metadata$pb_id), , drop = FALSE]

if (anyNA(pb_metadata$pb_id) || !identical(as.character(pb_metadata$pb_id), colnames(pb_counts))) {
    stop("Pseudobulk metadata did not align with count columns.")
}

candidate_coverage <- data.frame(candidate_hierarchy = c("RNA_ALL_436", "OCCUPANCY_SUPPORTED_284", "SHARED_OCCUPANCY_46", 
    "PROMOTER_OCCUPANCY_84"), locked_gene_n = c(436L, 284L, 46L, 84L), present_in_counts_n = c(sum(candidate_table$feature_id %in% 
    rownames(pb_counts)), sum(candidate_table$feature_id[candidate_table$occupancy_supported] %in% rownames(pb_counts)), 
    sum(candidate_table$feature_id[candidate_table$shared_occupancy_supported] %in% rownames(pb_counts)), sum(candidate_table$feature_id[candidate_table$any_promoter_overlap_support] %in% 
        rownames(pb_counts))), stringsAsFactors = FALSE)

candidate_coverage$coverage_fraction <- (candidate_coverage$present_in_counts_n/candidate_coverage$locked_gene_n)

write_csv_utf8(candidate_coverage, file.path(run_qc_root, "01_candidate_coverage.csv"))

if (any(candidate_coverage$coverage_fraction < 0.90000000000000002)) {
    stop("Locked candidate coverage is below 90% in pseudobulk counts.")
}

candidate_annotation <- candidate_table[, c("feature_id", "interaction_logFC", "occupancy_supported", "shared_occupancy_supported", 
    "any_promoter_overlap_support"), drop = FALSE]

colnames(candidate_annotation)[1] <- "gene"

candidate_annotation$RNA_ALL_436 <- TRUE

add_candidate_annotation <- function(result_table) {
    result_table$gene <- as.character(result_table$gene)
    out <- merge(result_table, candidate_annotation, by = "gene", all.x = TRUE, sort = FALSE)
    out <- out[match(result_table$gene, out$gene), , drop = FALSE]
    logical_flags <- c("RNA_ALL_436", "occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support")
    for (flag in logical_flags) {
        out[[flag]][is.na(out[[flag]])] <- FALSE
    }
    out
}

add_candidate_restricted_FDR <- function(result_table) {
    result_table$FDR_locked436 <- NA_real_
    result_table$FDR_occupancy284 <- NA_real_
    result_table$FDR_shared46 <- NA_real_
    result_table$FDR_promoter84 <- NA_real_
    hierarchy_map <- list(FDR_locked436 = result_table$RNA_ALL_436, FDR_occupancy284 = result_table$occupancy_supported, 
        FDR_shared46 = result_table$shared_occupancy_supported, FDR_promoter84 = result_table$any_promoter_overlap_support)
    for (column_name in names(hierarchy_map)) {
        index <- which(hierarchy_map[[column_name]] & is.finite(result_table$PValue))
        if (length(index)) {
            result_table[[column_name]][index] <- stats::p.adjust(result_table$PValue[index], method = "BH")
        }
    }
    result_table
}

finalize_one_df_result <- function(result_table, cell_type, model_family, contrast_name) {
    result_table <- add_candidate_annotation(result_table)
    result_table <- add_candidate_restricted_FDR(result_table)
    result_table$broad_cell_type <- cell_type
    result_table$model_family <- model_family
    result_table$contrast <- contrast_name
    result_table$effect_direction <- sign_label(result_table$logFC)
    result_table$whole_genome_DE_status <- ifelse(result_table$FDR < 0.050000000000000003 & abs(result_table$logFC) >= 
        0.5, ifelse(result_table$logFC > 0, "UP_FDR_LT_0.05_ABS_LOGFC_GE_0.5", "DOWN_FDR_LT_0.05_ABS_LOGFC_GE_0.5"), 
        "NOT_WHOLE_GENOME_DE")
    result_table$interaction_direction_concordant <- ifelse(result_table$RNA_ALL_436 & is.finite(result_table$interaction_logFC), 
        sign(result_table$logFC) == sign(result_table$interaction_logFC), NA)
    first_columns <- c("gene", "broad_cell_type", "model_family", "contrast", "logFC", "logCPM", "F", "PValue", 
        "FDR", "whole_genome_DE_status", "effect_direction", "RNA_ALL_436", "occupancy_supported", "shared_occupancy_supported", 
        "any_promoter_overlap_support", "interaction_logFC", "interaction_direction_concordant", "FDR_locked436", 
        "FDR_occupancy284", "FDR_shared46", "FDR_promoter84")
    result_table[, c(first_columns, setdiff(colnames(result_table), first_columns)), drop = FALSE]
}

finalize_omnibus_result <- function(result_table, cell_type) {
    result_table <- add_candidate_annotation(result_table)
    result_table <- add_candidate_restricted_FDR(result_table)
    result_table$broad_cell_type <- cell_type
    result_table$model_family <- "FERTILE_TEMPORAL_OMNIBUS"
    result_table$contrast <- "ANY_LH_TIME_EFFECT"
    result_table$whole_genome_DE_status <- ifelse(result_table$FDR < 0.050000000000000003, "TEMPORAL_OMNIBUS_FDR_LT_0.05", 
        "NOT_WHOLE_GENOME_TEMPORAL_DE")
    first_columns <- c("gene", "broad_cell_type", "model_family", "contrast", "logCPM", "F", "PValue", "FDR", "whole_genome_DE_status", 
        "RNA_ALL_436", "occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support", "interaction_logFC", 
        "FDR_locked436", "FDR_occupancy284", "FDR_shared46", "FDR_promoter84")
    result_table[, c(first_columns[first_columns %in% colnames(result_table)], setdiff(colnames(result_table), 
        first_columns)), drop = FALSE]
}

top_table_from_test <- function(test_object) {
    table <- edgeR::topTags(test_object, n = Inf, sort.by = "PValue")$table
    table$gene <- rownames(table)
    rownames(table) <- NULL
    if (!"FDR" %in% colnames(table)) {
        table$FDR <- stats::p.adjust(table$PValue, method = "BH")
    }
    table
}

make_directional_sets <- function(filtered_genes) {
    set_definitions <- list(RNA_ALL_436 = rep(TRUE, nrow(candidate_table)), OCCUPANCY_SUPPORTED_284 = candidate_table$occupancy_supported, 
        SHARED_OCCUPANCY_46 = candidate_table$shared_occupancy_supported, PROMOTER_OCCUPANCY_84 = candidate_table$any_promoter_overlap_support)
    output <- list()
    for (set_name in names(set_definitions)) {
        one <- candidate_table[set_definitions[[set_name]], c("feature_id", "interaction_logFC"), drop = FALSE]
        one$index <- match(one$feature_id, filtered_genes)
        one <- one[!is.na(one$index), , drop = FALSE]
        one$weight <- sign(one$interaction_logFC)
        output[[set_name]] <- data.frame(index = as.integer(one$index), weight = as.numeric(one$weight), stringsAsFactors = FALSE)
    }
    output
}

run_directional_fry <- function(y, design, contrast, cell_type, model_family, contrast_name) {
    sets <- make_directional_sets(rownames(y$counts))
    fry_result <- limma::fry(y = y, index = sets, design = design, contrast = contrast, sort = "none")
    fry_result$program_name <- rownames(fry_result)
    rownames(fry_result) <- NULL
    fry_result$broad_cell_type <- cell_type
    fry_result$model_family <- model_family
    fry_result$contrast <- contrast_name
    first_columns <- c("broad_cell_type", "model_family", "contrast", "program_name")
    fry_result[, c(first_columns, setdiff(colnames(fry_result), first_columns)), drop = FALSE]
}

fit_edgeR_QL <- function(count_matrix, metadata, group_factor) {
    counts_dense <- as.matrix(count_matrix)
    y <- edgeR::DGEList(counts = counts_dense, samples = metadata, group = group_factor)
    design <- stats::model.matrix(~0 + group_factor)
    colnames(design) <- levels(group_factor)
    rownames(design) <- colnames(y$counts)
    keep <- edgeR::filterByExpr(y, design = design)
    filtering_QC <- data.frame(input_gene_n = nrow(y), retained_gene_n = sum(keep), filtered_gene_n = sum(!keep), 
        retained_fraction = mean(keep), stringsAsFactors = FALSE)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y, method = "TMM")
    y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
    normalization_QC <- data.frame(pb_id = colnames(y$counts), group = as.character(group_factor), library_size = y$samples$lib.size, 
        norm_factor = y$samples$norm.factors, effective_library_size = (y$samples$lib.size * y$samples$norm.factors), 
        stringsAsFactors = FALSE)
    log_CPM <- edgeR::cpm(y, log = TRUE, prior.count = 2)
    PCA_fit <- stats::prcomp(t(log_CPM), center = TRUE, scale. = FALSE)
    variance_fraction <- (PCA_fit$sdev^2/sum(PCA_fit$sdev^2))
    PCA_table <- data.frame(pb_id = rownames(PCA_fit$x), PC1 = PCA_fit$x[, 1], PC2 = PCA_fit$x[, 2], PC1_variance = variance_fraction[1], 
        PC2_variance = variance_fraction[2], group = as.character(group_factor), stringsAsFactors = FALSE)
    list(y = y, design = design, fit = fit, filtering_QC = filtering_QC, normalization_QC = normalization_QC, PCA_table = PCA_table)
}

expected_temporal_counts <- function(eligibility_table, cell_type) {
    one <- eligibility_table[eligibility_table$broad_cell_type == cell_type, , drop = FALSE]
    if (nrow(one) != 1L) {
        stop("Eligibility row was not unique for ", cell_type)
    }
    counts <- c(LH3 = safe_numeric(one$fertile_LH3_n), LH5 = safe_numeric(one$fertile_LH5_n), LH7 = safe_numeric(one$fertile_LH7_n), 
        LH9 = safe_numeric(one$fertile_LH9_n), LH11 = safe_numeric(one$fertile_LH11_n))
    if (any(!is.finite(counts)) || any(counts < 2L)) {
        stop("Invalid temporal eligibility counts for ", cell_type)
    }
    as.integer(counts)
}

temporal_levels <- c("LH3", "LH5", "LH7", "LH9", "LH11")

planned_temporal_contrasts <- list(LH5_vs_LH3 = c(LH3 = -1, LH5 = 1, LH7 = 0, LH9 = 0, LH11 = 0), LH7_vs_LH5 = c(LH3 = 0, 
    LH5 = -1, LH7 = 1, LH9 = 0, LH11 = 0), LH9_vs_LH7 = c(LH3 = 0, LH5 = 0, LH7 = -1, LH9 = 1, LH11 = 0), LH11_vs_LH9 = c(LH3 = 0, 
    LH5 = 0, LH7 = 0, LH9 = -1, LH11 = 1), LH7_vs_LH3 = c(LH3 = -1, LH5 = 0, LH7 = 1, LH9 = 0, LH11 = 0), LH11_vs_LH3 = c(LH3 = -1, 
    LH5 = 0, LH7 = 0, LH9 = 0, LH11 = 1), LH11_vs_LH7 = c(LH3 = 0, LH5 = 0, LH7 = -1, LH9 = 0, LH11 = 1))

temporal_omnibus_contrast <- cbind(LH5_vs_LH3 = planned_temporal_contrasts$LH5_vs_LH3, LH7_vs_LH3 = planned_temporal_contrasts$LH7_vs_LH3, 
    LH9_vs_LH3 = c(LH3 = -1, LH5 = 0, LH7 = 0, LH9 = 1, LH11 = 0), LH11_vs_LH3 = planned_temporal_contrasts$LH11_vs_LH3)

dynamic_paths <- character()

for (cell_type in temporal_types) {
    code <- cell_type_code(cell_type)
    dynamic_paths <- c(dynamic_paths, file.path(per_model_root, paste0("T_", code, ".rds")), file.path(full_result_root, 
        paste0("T_", code, "_OMNI.csv")))
    for (contrast_name in names(planned_temporal_contrasts)) {
        dynamic_paths <- c(dynamic_paths, file.path(full_result_root, paste0("T_", code, "_", contrast_code(contrast_name), 
            ".csv")))
    }
}

for (cell_type in LH7_types) {
    code <- cell_type_code(cell_type)
    dynamic_paths <- c(dynamic_paths, file.path(per_model_root, paste0("L7_", code, ".rds")), file.path(full_result_root, 
        paste0("L7_", code, "_PRI.csv")), file.path(full_result_root, paste0("L7_", code, "_STD.csv")), file.path(full_result_root, 
        paste0("L7_", code, "_PGT.csv")))
}

dynamic_path_audit <- data.frame(path = normalize_slashes(dynamic_paths), character_n = nchar(normalize_slashes(dynamic_paths), 
    type = "chars"), stringsAsFactors = FALSE)

if (.Platform$OS.type == "windows" && any(dynamic_path_audit$character_n > 245L)) {
    bad_dynamic_paths <- dynamic_path_audit[dynamic_path_audit$character_n > 245L, , drop = FALSE]
    stop("One or more dynamic output paths remain too long:\n", paste(paste0(bad_dynamic_paths$character_n, " chars: ", 
        bad_dynamic_paths$path), collapse = "\n"))
}

write_csv_utf8(dynamic_path_audit, file.path(run_qc_root, "01_dynamic_path_audit.csv"))

materialize_temporal_outputs <- function(model_object, cell_type) {
    code <- cell_type_code(cell_type)
    write_csv_utf8(model_object$omnibus, file.path(full_result_root, paste0("T_", code, "_OMNI.csv")))
    for (contrast_name in names(model_object$pairwise)) {
        write_csv_utf8(model_object$pairwise[[contrast_name]], file.path(full_result_root, paste0("T_", code, "_", 
            contrast_code(contrast_name), ".csv")))
    }
    invisible(TRUE)
}

materialize_LH7_outputs <- function(model_object, cell_type) {
    code <- cell_type_code(cell_type)
    write_csv_utf8(model_object$primary, file.path(full_result_root, paste0("L7_", code, "_PRI.csv")))
    for (subgroup in names(model_object$sensitivity)) {
        contrast_name <- paste0("RIF_vs_", subgroup)
        write_csv_utf8(model_object$sensitivity[[subgroup]], file.path(full_result_root, paste0("L7_", code, "_", 
            contrast_code(contrast_name), ".csv")))
    }
    invisible(TRUE)
}

temporal_model_objects <- list()

temporal_status_rows <- list()

for (model_index in seq_along(temporal_types)) {
    cell_type <- temporal_types[model_index]
    cell_code <- cell_type_code(cell_type)
    old_cell_slug <- slugify(cell_type)
    model_rds <- file.path(per_model_root, paste0("T_", cell_code, ".rds"))
    old_model_rds <- file.path(previous_per_model_root, paste0("temporal_", old_cell_slug, ".rds"))
    if (!file.exists(model_rds) && file.exists(old_model_rds)) {
        old_object <- tryCatch(readRDS(old_model_rds), error = function(e) NULL)
        if (!is.null(old_object) && identical(old_object$status, "COMPLETE")) {
            saveRDS(old_object, model_rds, compress = "gzip")
            log_event("IMPORTED completed temporal model", cell_type, "from GENEDE_v1_1")
        }
        rm(old_object)
        gc()
    }
    if (file.exists(model_rds)) {
        existing <- tryCatch(readRDS(model_rds), error = function(e) NULL)
        if (!is.null(existing) && identical(existing$status, "COMPLETE")) {
            materialize_temporal_outputs(existing, cell_type)
            temporal_model_objects[[cell_type]] <- existing
            temporal_status_rows[[length(temporal_status_rows) + 1L]] <- existing$status_table
            log_event("SKIP completed temporal model", cell_type)
            next
        }
    }
    model_start <- Sys.time()
    metadata_one <- pb_metadata[pb_metadata$broad_cell_type == cell_type & pb_metadata$clinical_group == "FERTILE" & 
        pb_metadata$temporal_fertile_eligible & pb_metadata$minimum_cell_threshold_pass, , drop = FALSE]
    metadata_one$time_group <- factor(paste0("LH", metadata_one$LH_day), levels = temporal_levels)
    observed_counts <- as.integer(table(metadata_one$time_group))
    expected_counts <- expected_temporal_counts(eligibility, cell_type)
    if (!identical(observed_counts, as.integer(expected_counts)) || any(observed_counts < 2L) || anyDuplicated(metadata_one$donor_id) > 
        0L) {
        stop("Temporal donor identity failed for ", cell_type, ". Observed ", paste(temporal_levels, observed_counts, 
            sep = "=", collapse = ";"))
    }
    column_index <- match(metadata_one$pb_id, colnames(pb_counts))
    if (anyNA(column_index)) {
        stop("Temporal count columns were not found for ", cell_type)
    }
    fit_result <- fit_edgeR_QL(count_matrix = pb_counts[, column_index, drop = FALSE], metadata = metadata_one, 
        group_factor = metadata_one$time_group)
    omnibus_test <- edgeR::glmQLFTest(fit_result$fit, contrast = temporal_omnibus_contrast)
    omnibus_table <- top_table_from_test(omnibus_test)
    omnibus_table <- finalize_omnibus_result(omnibus_table, cell_type)
    omnibus_path <- file.path(full_result_root, paste0("T_", cell_code, "_OMNI.csv"))
    write_csv_utf8(omnibus_table, omnibus_path)
    pairwise_tables <- list()
    fry_tables <- list()
    for (contrast_name in names(planned_temporal_contrasts)) {
        contrast_vector <- planned_temporal_contrasts[[contrast_name]]
        pairwise_test <- edgeR::glmQLFTest(fit_result$fit, contrast = contrast_vector)
        pairwise_table <- top_table_from_test(pairwise_test)
        pairwise_table <- finalize_one_df_result(pairwise_table, cell_type = cell_type, model_family = "FERTILE_TEMPORAL_PAIRWISE", 
            contrast_name = contrast_name)
        pairwise_path <- file.path(full_result_root, paste0("T_", cell_code, "_", contrast_code(contrast_name), 
            ".csv"))
        write_csv_utf8(pairwise_table, pairwise_path)
        pairwise_tables[[contrast_name]] <- pairwise_table
        fry_tables[[contrast_name]] <- run_directional_fry(y = fit_result$y, design = fit_result$design, contrast = contrast_vector, 
            cell_type = cell_type, model_family = "FERTILE_TEMPORAL_PAIRWISE", contrast_name = contrast_name)
    }
    filtering_QC <- fit_result$filtering_QC
    filtering_QC$broad_cell_type <- cell_type
    filtering_QC$model_family <- "FERTILE_TEMPORAL"
    normalization_QC <- fit_result$normalization_QC
    normalization_QC$broad_cell_type <- cell_type
    normalization_QC$model_family <- "FERTILE_TEMPORAL"
    PCA_table <- fit_result$PCA_table
    PCA_table$broad_cell_type <- cell_type
    PCA_table$model_family <- "FERTILE_TEMPORAL"
    PCA_table$model_label <- paste0("Temporal | ", cell_type)
    elapsed_minutes <- round(as.numeric(difftime(Sys.time(), model_start, units = "mins")), 2)
    status_table <- data.frame(model_family = "FERTILE_TEMPORAL", broad_cell_type = cell_type, status = "COMPLETE", 
        donor_n = nrow(metadata_one), minimum_group_n = min(observed_counts), filtered_gene_n = filtering_QC$retained_gene_n, 
        omnibus_genome_FDR_lt_0_05_n = sum(omnibus_table$FDR < 0.050000000000000003), elapsed_minutes = elapsed_minutes, 
        stringsAsFactors = FALSE)
    output_object <- list(status = "COMPLETE", cell_type = cell_type, metadata = metadata_one, filtering_QC = filtering_QC, 
        normalization_QC = normalization_QC, PCA_table = PCA_table, omnibus = omnibus_table, pairwise = pairwise_tables, 
        fry = fry_tables, status_table = status_table)
    saveRDS(output_object, model_rds, compress = "gzip")
    temporal_model_objects[[cell_type]] <- output_object
    temporal_status_rows[[length(temporal_status_rows) + 1L]] <- status_table
    log_event("DONE temporal gene model", model_index, "/", length(temporal_types), cell_type, "retained_genes=", 
        filtering_QC$retained_gene_n, "omnibus_FDR_lt_0.05=", status_table$omnibus_genome_FDR_lt_0_05_n, "elapsed_min=", 
        elapsed_minutes)
    rm(fit_result, omnibus_test, omnibus_table, pairwise_tables, fry_tables, output_object)
    gc()
}

LH7_model_objects <- list()

LH7_status_rows <- list()

for (model_index in seq_along(LH7_types)) {
    cell_type <- LH7_types[model_index]
    cell_code <- cell_type_code(cell_type)
    old_cell_slug <- slugify(cell_type)
    model_rds <- file.path(per_model_root, paste0("L7_", cell_code, ".rds"))
    old_model_rds <- file.path(previous_per_model_root, paste0("LH7_", old_cell_slug, ".rds"))
    if (!file.exists(model_rds) && file.exists(old_model_rds)) {
        old_object <- tryCatch(readRDS(old_model_rds), error = function(e) NULL)
        if (!is.null(old_object) && identical(old_object$status, "COMPLETE")) {
            saveRDS(old_object, model_rds, compress = "gzip")
            log_event("IMPORTED completed LH7 model", cell_type, "from GENEDE_v1_1")
        }
        rm(old_object)
        gc()
    }
    if (file.exists(model_rds)) {
        existing <- tryCatch(readRDS(model_rds), error = function(e) NULL)
        if (!is.null(existing) && identical(existing$status, "COMPLETE")) {
            materialize_LH7_outputs(existing, cell_type)
            LH7_model_objects[[cell_type]] <- existing
            LH7_status_rows[[length(LH7_status_rows) + 1L]] <- existing$status_table
            log_event("SKIP completed LH7 model", cell_type)
            next
        }
    }
    model_start <- Sys.time()
    metadata_primary <- pb_metadata[pb_metadata$broad_cell_type == cell_type & pb_metadata$LH7_RIF_comparison_eligible & 
        pb_metadata$minimum_cell_threshold_pass, , drop = FALSE]
    metadata_primary$group <- factor(metadata_primary$clinical_group, levels = c("FERTILE", "RIF"))
    primary_group_counts <- as.integer(table(metadata_primary$group))
    if (any(primary_group_counts < 3L) || anyDuplicated(metadata_primary$donor_id) > 0L) {
        stop("LH+7 donor identity failed for ", cell_type)
    }
    primary_column_index <- match(metadata_primary$pb_id, colnames(pb_counts))
    if (anyNA(primary_column_index)) {
        stop("LH+7 count columns were not found for ", cell_type)
    }
    primary_fit <- fit_edgeR_QL(count_matrix = pb_counts[, primary_column_index, drop = FALSE], metadata = metadata_primary, 
        group_factor = metadata_primary$group)
    primary_contrast <- c(FERTILE = -1, RIF = 1)
    primary_test <- edgeR::glmQLFTest(primary_fit$fit, contrast = primary_contrast)
    primary_table <- top_table_from_test(primary_test)
    primary_table <- finalize_one_df_result(primary_table, cell_type = cell_type, model_family = "LH7_RIF_VS_ALL_FERTILE", 
        contrast_name = "RIF_vs_all_fertile_LH7")
    primary_path <- file.path(full_result_root, paste0("L7_", cell_code, "_PRI.csv"))
    write_csv_utf8(primary_table, primary_path)
    primary_fry <- run_directional_fry(y = primary_fit$y, design = primary_fit$design, contrast = primary_contrast, 
        cell_type = cell_type, model_family = "LH7_RIF_VS_ALL_FERTILE", contrast_name = "RIF_vs_all_fertile_LH7")
    sensitivity_tables <- list()
    sensitivity_fry <- list()
    sensitivity_filtering <- list()
    sensitivity_normalization <- list()
    sensitivity_PCA <- list()
    for (subgroup in c("LH7_STANDARD", "LH7_PGT")) {
        metadata_sensitivity <- metadata_primary[metadata_primary$clinical_group == "RIF" | (metadata_primary$clinical_group == 
            "FERTILE" & metadata_primary$LH7_fertile_subgroup == subgroup), , drop = FALSE]
        metadata_sensitivity$group <- factor(ifelse(metadata_sensitivity$clinical_group == "RIF", "RIF", subgroup), 
            levels = c(subgroup, "RIF"))
        sensitivity_counts <- as.integer(table(metadata_sensitivity$group))
        if (any(sensitivity_counts < 2L) || anyDuplicated(metadata_sensitivity$donor_id) > 0L) {
            stop("LH+7 sensitivity donor identity failed for ", cell_type, " / ", subgroup)
        }
        sensitivity_columns <- match(metadata_sensitivity$pb_id, colnames(pb_counts))
        sensitivity_fit <- fit_edgeR_QL(count_matrix = pb_counts[, sensitivity_columns, drop = FALSE], metadata = metadata_sensitivity, 
            group_factor = metadata_sensitivity$group)
        sensitivity_contrast <- c(-1, 1)
        names(sensitivity_contrast) <- levels(metadata_sensitivity$group)
        sensitivity_test <- edgeR::glmQLFTest(sensitivity_fit$fit, contrast = sensitivity_contrast)
        sensitivity_table <- top_table_from_test(sensitivity_test)
        sensitivity_table <- finalize_one_df_result(sensitivity_table, cell_type = cell_type, model_family = "LH7_RIF_SENSITIVITY", 
            contrast_name = paste0("RIF_vs_", subgroup))
        sensitivity_path <- file.path(full_result_root, paste0("L7_", cell_code, "_", contrast_code(paste0("RIF_vs_", 
            subgroup)), ".csv"))
        write_csv_utf8(sensitivity_table, sensitivity_path)
        sensitivity_tables[[subgroup]] <- sensitivity_table
        sensitivity_fry[[subgroup]] <- run_directional_fry(y = sensitivity_fit$y, design = sensitivity_fit$design, 
            contrast = sensitivity_contrast, cell_type = cell_type, model_family = "LH7_RIF_SENSITIVITY", contrast_name = paste0("RIF_vs_", 
                subgroup))
        one_filtering <- sensitivity_fit$filtering_QC
        one_filtering$broad_cell_type <- cell_type
        one_filtering$model_family <- "LH7_RIF_SENSITIVITY"
        one_filtering$contrast <- paste0("RIF_vs_", subgroup)
        sensitivity_filtering[[subgroup]] <- one_filtering
        one_normalization <- sensitivity_fit$normalization_QC
        one_normalization$broad_cell_type <- cell_type
        one_normalization$model_family <- "LH7_RIF_SENSITIVITY"
        one_normalization$contrast <- paste0("RIF_vs_", subgroup)
        sensitivity_normalization[[subgroup]] <- one_normalization
        one_PCA <- sensitivity_fit$PCA_table
        one_PCA$broad_cell_type <- cell_type
        one_PCA$model_family <- "LH7_RIF_SENSITIVITY"
        one_PCA$model_label <- paste0("LH7 ", subgroup, " | ", cell_type)
        sensitivity_PCA[[subgroup]] <- one_PCA
        rm(sensitivity_fit, sensitivity_test, sensitivity_table)
        gc()
    }
    primary_filtering <- primary_fit$filtering_QC
    primary_filtering$broad_cell_type <- cell_type
    primary_filtering$model_family <- "LH7_RIF_VS_ALL_FERTILE"
    primary_filtering$contrast <- "RIF_vs_all_fertile_LH7"
    primary_normalization <- primary_fit$normalization_QC
    primary_normalization$broad_cell_type <- cell_type
    primary_normalization$model_family <- "LH7_RIF_VS_ALL_FERTILE"
    primary_normalization$contrast <- "RIF_vs_all_fertile_LH7"
    primary_PCA <- primary_fit$PCA_table
    primary_PCA$broad_cell_type <- cell_type
    primary_PCA$model_family <- "LH7_RIF_VS_ALL_FERTILE"
    primary_PCA$model_label <- paste0("LH7 primary | ", cell_type)
    elapsed_minutes <- round(as.numeric(difftime(Sys.time(), model_start, units = "mins")), 2)
    status_table <- data.frame(model_family = "LH7_RIF_VS_ALL_FERTILE", broad_cell_type = cell_type, status = "COMPLETE", 
        fertile_n = primary_group_counts[1], RIF_n = primary_group_counts[2], filtered_gene_n = primary_filtering$retained_gene_n, 
        primary_genome_DE_n = sum(primary_table$FDR < 0.050000000000000003 & abs(primary_table$logFC) >= 0.5), 
        elapsed_minutes = elapsed_minutes, stringsAsFactors = FALSE)
    output_object <- list(status = "COMPLETE", cell_type = cell_type, metadata_primary = metadata_primary, primary_filtering = primary_filtering, 
        primary_normalization = primary_normalization, primary_PCA = primary_PCA, primary = primary_table, primary_fry = primary_fry, 
        sensitivity = sensitivity_tables, sensitivity_fry = sensitivity_fry, sensitivity_filtering = sensitivity_filtering, 
        sensitivity_normalization = sensitivity_normalization, sensitivity_PCA = sensitivity_PCA, status_table = status_table)
    saveRDS(output_object, model_rds, compress = "gzip")
    LH7_model_objects[[cell_type]] <- output_object
    LH7_status_rows[[length(LH7_status_rows) + 1L]] <- status_table
    log_event("DONE LH7 gene model", model_index, "/", length(LH7_types), cell_type, "retained_genes=", primary_filtering$retained_gene_n, 
        "primary_DE_n=", status_table$primary_genome_DE_n, "elapsed_min=", elapsed_minutes)
    rm(primary_fit, primary_test, primary_table, primary_fry, sensitivity_tables, sensitivity_fry, output_object)
    gc()
}

temporal_status <- do.call(rbind, temporal_status_rows)

LH7_status <- do.call(rbind, LH7_status_rows)

write_csv_utf8(temporal_status, file.path(table_root, "01_temporal_model_status.csv"))

write_csv_utf8(LH7_status, file.path(table_root, "02_LH7_model_status.csv"))

temporal_omnibus_all <- data.table::rbindlist(lapply(temporal_model_objects, function(x) x$omnibus), use.names = TRUE, 
    fill = TRUE)

temporal_pairwise_all <- data.table::rbindlist(unlist(lapply(temporal_model_objects, function(x) x$pairwise), recursive = FALSE), 
    use.names = TRUE, fill = TRUE)

temporal_fry_all <- data.table::rbindlist(unlist(lapply(temporal_model_objects, function(x) x$fry), recursive = FALSE), 
    use.names = TRUE, fill = TRUE)

temporal_filtering_all <- data.table::rbindlist(lapply(temporal_model_objects, function(x) x$filtering_QC), use.names = TRUE, 
    fill = TRUE)

temporal_normalization_all <- data.table::rbindlist(lapply(temporal_model_objects, function(x) x$normalization_QC), 
    use.names = TRUE, fill = TRUE)

temporal_PCA_all <- data.table::rbindlist(lapply(temporal_model_objects, function(x) x$PCA_table), use.names = TRUE, 
    fill = TRUE)

temporal_omnibus_significant <- temporal_omnibus_all[temporal_omnibus_all$FDR < 0.050000000000000003, , drop = FALSE]

temporal_pairwise_significant <- temporal_pairwise_all[temporal_pairwise_all$FDR < 0.050000000000000003 & abs(temporal_pairwise_all$logFC) >= 
    0.5, , drop = FALSE]

write_csv_utf8(as.data.frame(temporal_omnibus_significant), file.path(table_root, "03_temporal_omnibus_significant.csv"))

write_csv_utf8(as.data.frame(temporal_pairwise_significant), file.path(table_root, "04_temporal_pairwise_significant.csv"))

LH7_primary_all <- data.table::rbindlist(lapply(LH7_model_objects, function(x) x$primary), use.names = TRUE, fill = TRUE)

LH7_sensitivity_all <- data.table::rbindlist(unlist(lapply(LH7_model_objects, function(x) x$sensitivity), recursive = FALSE), 
    use.names = TRUE, fill = TRUE)

LH7_fry_primary_all <- data.table::rbindlist(lapply(LH7_model_objects, function(x) x$primary_fry), use.names = TRUE, 
    fill = TRUE)

LH7_fry_sensitivity_all <- data.table::rbindlist(unlist(lapply(LH7_model_objects, function(x) x$sensitivity_fry), 
    recursive = FALSE), use.names = TRUE, fill = TRUE)

LH7_filtering_primary_all <- data.table::rbindlist(lapply(LH7_model_objects, function(x) x$primary_filtering), 
    use.names = TRUE, fill = TRUE)

LH7_normalization_primary_all <- data.table::rbindlist(lapply(LH7_model_objects, function(x) x$primary_normalization), 
    use.names = TRUE, fill = TRUE)

LH7_PCA_primary_all <- data.table::rbindlist(lapply(LH7_model_objects, function(x) x$primary_PCA), use.names = TRUE, 
    fill = TRUE)

LH7_filtering_sensitivity_all <- data.table::rbindlist(unlist(lapply(LH7_model_objects, function(x) x$sensitivity_filtering), 
    recursive = FALSE), use.names = TRUE, fill = TRUE)

LH7_normalization_sensitivity_all <- data.table::rbindlist(unlist(lapply(LH7_model_objects, function(x) x$sensitivity_normalization), 
    recursive = FALSE), use.names = TRUE, fill = TRUE)

LH7_PCA_sensitivity_all <- data.table::rbindlist(unlist(lapply(LH7_model_objects, function(x) x$sensitivity_PCA), 
    recursive = FALSE), use.names = TRUE, fill = TRUE)

LH7_primary_significant <- LH7_primary_all[LH7_primary_all$FDR < 0.050000000000000003 & abs(LH7_primary_all$logFC) >= 
    0.5, , drop = FALSE]

LH7_sensitivity_significant <- LH7_sensitivity_all[LH7_sensitivity_all$FDR < 0.050000000000000003 & abs(LH7_sensitivity_all$logFC) >= 
    0.5, , drop = FALSE]

write_csv_utf8(as.data.frame(LH7_primary_significant), file.path(table_root, "05_LH7_primary_significant.csv"))

write_csv_utf8(as.data.frame(LH7_sensitivity_significant), file.path(table_root, "06_LH7_sensitivity_significant.csv"))

candidate_temporal_omnibus <- temporal_omnibus_all[temporal_omnibus_all$RNA_ALL_436, , drop = FALSE]

candidate_temporal_pairwise <- temporal_pairwise_all[temporal_pairwise_all$RNA_ALL_436, , drop = FALSE]

candidate_LH7_primary <- LH7_primary_all[LH7_primary_all$RNA_ALL_436, , drop = FALSE]

candidate_LH7_sensitivity <- LH7_sensitivity_all[LH7_sensitivity_all$RNA_ALL_436, , drop = FALSE]

write_csv_utf8(as.data.frame(candidate_temporal_omnibus), file.path(candidate_root, "01_locked436_temporal_omnibus.csv"))

write_csv_utf8(as.data.frame(candidate_temporal_pairwise), file.path(candidate_root, "02_locked436_temporal_pairwise.csv"))

write_csv_utf8(as.data.frame(candidate_LH7_primary), file.path(candidate_root, "03_locked436_LH7_primary.csv"))

write_csv_utf8(as.data.frame(candidate_LH7_sensitivity), file.path(candidate_root, "04_locked436_LH7_sensitivity.csv"))

all_result_tables <- list(temporal_omnibus = as.data.frame(temporal_omnibus_all), temporal_pairwise = as.data.frame(temporal_pairwise_all), 
    LH7_primary = as.data.frame(LH7_primary_all), LH7_sensitivity = as.data.frame(LH7_sensitivity_all))

VDR_rows <- data.table::rbindlist(lapply(names(all_result_tables), function(table_name) {
    one <- all_result_tables[[table_name]]
    one <- one[one$gene == "VDR", , drop = FALSE]
    if (nrow(one)) {
        one$result_family <- table_name
    }
    one
}), use.names = TRUE, fill = TRUE)

write_csv_utf8(as.data.frame(VDR_rows), file.path(candidate_root, "05_VDR_gene_all_models.csv"))

fry_all <- data.table::rbindlist(list(temporal_fry_all, LH7_fry_primary_all, LH7_fry_sensitivity_all), use.names = TRUE, 
    fill = TRUE)

write_csv_utf8(as.data.frame(fry_all), file.path(candidate_root, "06_directional_VDR_program_fry.csv"))

temporal_omnibus_summary <- data.table::as.data.table(temporal_omnibus_all)[, .(retained_gene_n = .N, whole_genome_FDR_lt_0_05_n = sum(FDR < 
    0.050000000000000003), locked436_FDR_lt_0_05_n = sum(RNA_ALL_436 & FDR_locked436 < 0.050000000000000003, na.rm = TRUE), 
    occupancy284_FDR_lt_0_05_n = sum(occupancy_supported & FDR_occupancy284 < 0.050000000000000003, na.rm = TRUE), 
    shared46_FDR_lt_0_05_n = sum(shared_occupancy_supported & FDR_shared46 < 0.050000000000000003, na.rm = TRUE), 
    promoter84_FDR_lt_0_05_n = sum(any_promoter_overlap_support & FDR_promoter84 < 0.050000000000000003, na.rm = TRUE)), 
    by = broad_cell_type]

temporal_pairwise_summary <- data.table::as.data.table(temporal_pairwise_all)[, .(whole_genome_up_n = sum(FDR < 
    0.050000000000000003 & logFC >= 0.5), whole_genome_down_n = sum(FDR < 0.050000000000000003 & logFC <= -0.5), 
    locked436_focused_FDR_lt_0_05_n = sum(RNA_ALL_436 & FDR_locked436 < 0.050000000000000003, na.rm = TRUE), locked436_concordant_focused_n = sum(RNA_ALL_436 & 
        FDR_locked436 < 0.050000000000000003 & interaction_direction_concordant, na.rm = TRUE)), by = .(broad_cell_type, 
    contrast)]

LH7_primary_summary <- data.table::as.data.table(LH7_primary_all)[, .(whole_genome_up_n = sum(FDR < 0.050000000000000003 & 
    logFC >= 0.5), whole_genome_down_n = sum(FDR < 0.050000000000000003 & logFC <= -0.5), locked436_focused_FDR_lt_0_05_n = sum(RNA_ALL_436 & 
    FDR_locked436 < 0.050000000000000003, na.rm = TRUE), locked436_concordant_focused_n = sum(RNA_ALL_436 & FDR_locked436 < 
    0.050000000000000003 & interaction_direction_concordant, na.rm = TRUE)), by = broad_cell_type]

write_csv_utf8(as.data.frame(temporal_omnibus_summary), file.path(table_root, "07_temporal_omnibus_summary.csv"))

write_csv_utf8(as.data.frame(temporal_pairwise_summary), file.path(table_root, "08_temporal_pairwise_summary.csv"))

write_csv_utf8(as.data.frame(LH7_primary_summary), file.path(table_root, "09_LH7_primary_summary.csv"))

top_temporal_omnibus <- data.table::as.data.table(temporal_omnibus_all)[order(FDR, -F), head(.SD, 100L), by = broad_cell_type]

top_temporal_pairwise <- data.table::as.data.table(temporal_pairwise_all)[order(FDR, -abs(logFC)), head(.SD, 50L), 
    by = .(broad_cell_type, contrast)]

top_LH7_primary <- data.table::as.data.table(LH7_primary_all)[order(FDR, -abs(logFC)), head(.SD, 100L), by = broad_cell_type]

write_csv_utf8(as.data.frame(top_temporal_omnibus), file.path(table_root, "10_top_temporal_omnibus.csv"))

write_csv_utf8(as.data.frame(top_temporal_pairwise), file.path(table_root, "11_top_temporal_pairwise.csv"))

write_csv_utf8(as.data.frame(top_LH7_primary), file.path(table_root, "12_top_LH7_primary.csv"))

filtering_all <- data.table::rbindlist(list(temporal_filtering_all, LH7_filtering_primary_all, LH7_filtering_sensitivity_all), 
    use.names = TRUE, fill = TRUE)

normalization_all <- data.table::rbindlist(list(temporal_normalization_all, LH7_normalization_primary_all, LH7_normalization_sensitivity_all), 
    use.names = TRUE, fill = TRUE)

PCA_all <- data.table::rbindlist(list(temporal_PCA_all, LH7_PCA_primary_all, LH7_PCA_sensitivity_all), use.names = TRUE, 
    fill = TRUE)

write_csv_utf8(as.data.frame(filtering_all), file.path(run_qc_root, "02_filtering_QC.csv"))

write_csv_utf8(as.data.frame(normalization_all), file.path(run_qc_root, "03_normalization_QC.csv"))

write_csv_utf8(as.data.frame(PCA_all), file.path(source_data_root, "source_model_PCA.csv"))

theme_qc <- ggplot2::theme_bw(base_size = 9) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "grey95"))

PCA_plot_data <- as.data.frame(PCA_all)

plot_PCA <- ggplot2::ggplot(PCA_plot_data, ggplot2::aes(x = PC1, y = PC2, shape = group)) + ggplot2::geom_point(size = 1.7) + 
    ggplot2::facet_wrap(~model_label, scales = "free", ncol = 3) + ggplot2::labs(title = "Gene-level pseudobulk model PCA", 
    subtitle = "PCA calculated after cell-type-specific expression filtering and TMM normalization", x = "PC1", 
    y = "PC2", shape = "Model group") + theme_qc + ggplot2::theme(strip.text = ggplot2::element_text(size = 6))

save_plot_pair(plot_PCA, "QC1_model_PCA", width = 13, height = 12)

temporal_count_plot <- as.data.frame(temporal_omnibus_summary)

temporal_count_plot <- reshape(temporal_count_plot[, c("broad_cell_type", "whole_genome_FDR_lt_0_05_n", "locked436_FDR_lt_0_05_n")], 
    varying = c("whole_genome_FDR_lt_0_05_n", "locked436_FDR_lt_0_05_n"), v.names = "gene_n", timevar = "result_scope", 
    times = c("Genome-wide FDR", "Locked 436 focused FDR"), direction = "long")

plot_temporal_counts <- ggplot2::ggplot(temporal_count_plot, ggplot2::aes(x = broad_cell_type, y = gene_n, fill = result_scope)) + 
    ggplot2::geom_col(position = "dodge") + ggplot2::scale_y_continuous(labels = scales::comma) + ggplot2::labs(title = "Temporal omnibus gene-level findings", 
    subtitle = "Candidate-restricted FDR is secondary and does not replace genome-wide FDR", x = "Broad cell type", 
    y = "Gene count", fill = "Multiplicity scope") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_temporal_counts, "QC2_temporal_DE_counts", width = 10, height = 6)

fry_plot_data <- as.data.frame(fry_all)

direction_column <- if ("Direction" %in% colnames(fry_plot_data)) {
    "Direction"
} else {
    ""
}

p_column <- if ("FDR" %in% colnames(fry_plot_data)) {
    "FDR"
} else if ("PValue" %in% colnames(fry_plot_data)) {
    "PValue"
} else {
    ""
}

if (nzchar(p_column)) {
    fry_plot_data$plot_FDR <- safe_numeric(fry_plot_data[[p_column]])
    fry_plot_data$minus_log10_FDR <- -log10(pmax(fry_plot_data$plot_FDR, .Machine$double.xmin))
    fry_plot_data$direction_for_plot <- if (nzchar(direction_column)) {
        as.character(fry_plot_data[[direction_column]])
    }
    else {
        "NA"
    }
    plot_fry <- ggplot2::ggplot(fry_plot_data, ggplot2::aes(x = contrast, y = broad_cell_type, size = minus_log10_FDR, 
        shape = direction_for_plot)) + ggplot2::geom_point() + ggplot2::facet_wrap(~program_name, scales = "free_x", 
        ncol = 2) + ggplot2::geom_hline(yintercept = NA_real_) + ggplot2::labs(title = "Directionally weighted locked VDR-program tests", 
        subtitle = "limma fry; direction is defined relative to the locked calcitriol-dependent interaction signs", 
        x = "Contrast", y = "Broad cell type", size = "-log10(FDR)", shape = "Direction") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
        vjust = 0.5, hjust = 1, size = 5), strip.text = ggplot2::element_text(size = 7))
    save_plot_pair(plot_fry, "QC3_directional_VDR_fry", width = 16, height = 10)
}

volcano_data <- as.data.frame(LH7_primary_all)

volcano_data$minus_log10_FDR <- -log10(pmax(volcano_data$FDR, .Machine$double.xmin))

volcano_data$plot_class <- ifelse(volcano_data$FDR < 0.050000000000000003 & abs(volcano_data$logFC) >= 0.5, "Whole-genome DE", 
    ifelse(volcano_data$RNA_ALL_436, "Locked VDR candidate", "Other gene"))

plot_LH7_volcano <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = logFC, y = minus_log10_FDR, shape = plot_class)) + 
    ggplot2::geom_point(size = 0.45000000000000001, alpha = 0.55000000000000004) + ggplot2::geom_vline(xintercept = c(-0.5, 
    0.5), linetype = "dashed") + ggplot2::geom_hline(yintercept = -log10(0.050000000000000003), linetype = "dashed") + 
    ggplot2::facet_wrap(~broad_cell_type, scales = "free", ncol = 3) + ggplot2::labs(title = "LH+7 RIF versus fertile gene-level effects", 
    subtitle = "Dashed lines show FDR 0.05 and absolute log2 fold-change 0.5 reporting thresholds", x = "RIF minus fertile log2 fold-change", 
    y = "-log10 genome-wide FDR", shape = "Gene class") + theme_qc

save_plot_pair(plot_LH7_volcano, "QC4_LH7_volcano", width = 12, height = 8)

top_candidate_LH7 <- data.table::as.data.table(candidate_LH7_primary)[order(PValue, -abs(logFC)), head(.SD, 8L), 
    by = broad_cell_type]

top_candidate_LH7$gene_label <- paste0(top_candidate_LH7$gene, " | ", top_candidate_LH7$broad_cell_type)

top_candidate_LH7$gene_label <- factor(top_candidate_LH7$gene_label, levels = rev(unique(top_candidate_LH7$gene_label)))

plot_candidate_LH7 <- ggplot2::ggplot(as.data.frame(top_candidate_LH7), ggplot2::aes(x = logFC, y = gene_label, 
    size = -log10(pmax(PValue, .Machine$double.xmin)), shape = interaction_direction_concordant)) + ggplot2::geom_vline(xintercept = 0, 
    linetype = "dashed") + ggplot2::geom_point() + ggplot2::labs(title = "Top locked VDR candidates in the LH+7 comparison", 
    subtitle = "Ranking uses nominal gene-level evidence; formal claims require genome-wide or explicitly labelled focused FDR", 
    x = "RIF minus fertile log2 fold-change", y = "Candidate gene and broad cell type", size = "-log10(P)", shape = "Concordant with\ninteraction direction") + 
    theme_qc + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5))

save_plot_pair(plot_candidate_LH7, "QC5_LH7_top_locked_candidates", width = 11, height = 13)

LH7_count_plot <- as.data.frame(LH7_primary_summary)

LH7_count_long <- data.table::melt(data.table::as.data.table(LH7_count_plot), id.vars = "broad_cell_type", measure.vars = c("whole_genome_up_n", 
    "whole_genome_down_n", "locked436_focused_FDR_lt_0_05_n"), variable.name = "result_class", value.name = "gene_n")

plot_LH7_counts <- ggplot2::ggplot(as.data.frame(LH7_count_long), ggplot2::aes(x = broad_cell_type, y = gene_n, 
    fill = result_class)) + ggplot2::geom_col(position = "dodge") + ggplot2::scale_y_continuous(labels = scales::comma) + 
    ggplot2::labs(title = "LH+7 RIF versus fertile gene-level summary", subtitle = "Locked-436 focused FDR is displayed separately from genome-wide DE", 
        x = "Broad cell type", y = "Gene count", fill = "Result class") + theme_qc + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, 
    hjust = 1))

save_plot_pair(plot_LH7_counts, "QC6_LH7_DE_counts", width = 10, height = 6)

expected_fry_n <- (length(temporal_types) * length(planned_temporal_contrasts) * 4L) + (length(LH7_types) * 4L) + 
    (length(LH7_types) * 2L * 4L)

gate_checks <- data.frame(check = c("upstream_gate_authorizes_gene_DE", "pseudobulk_matrix_dimensions_locked", 
    "pseudobulk_metadata_aligned", "locked_candidate_hierarchy_valid", "candidate_coverage_above_90_percent", "five_temporal_models_complete", 
    "six_LH7_primary_models_complete", "twelve_LH7_sensitivity_models_present", "all_models_retain_at_least_1000_genes", 
    "all_temporal_omnibus_FDR_values_finite", "all_LH7_primary_FDR_values_finite", "directional_fry_result_count_complete", 
    "composition_inference_remains_disabled", "subtype_inference_remains_disabled"), pass = c(upstream_gate$gene_level_pseudobulk_DE_authorized[[1]], 
    nrow(pb_counts) == 36601L && ncol(pb_counts) == 222L, identical(colnames(pb_counts), as.character(pb_metadata$pb_id)), 
    nrow(candidate_table) == 436L && sum(candidate_table$occupancy_supported) == 284L && sum(candidate_table$shared_occupancy_supported) == 
        46L && sum(candidate_table$any_promoter_overlap_support) == 84L, all(candidate_coverage$coverage_fraction >= 
        0.90000000000000002), nrow(temporal_status) == 5L && all(temporal_status$status == "COMPLETE"), nrow(LH7_status) == 
        6L && all(LH7_status$status == "COMPLETE"), nrow(unique(LH7_sensitivity_all[, c("broad_cell_type", "contrast")])) == 
        12L, all(filtering_all$retained_gene_n >= 1000L), all(is.finite(temporal_omnibus_all$FDR)), all(is.finite(LH7_primary_all$FDR)), 
    nrow(fry_all) == expected_fry_n, !upstream_gate$composition_inference_authorized[[1]], !upstream_gate$subtype_inference_authorized[[1]]), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "04_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass) {
    "PASS_GENE_LEVEL_PSEUDOBULK_DE_COMPLETE"
} else {
    "HOLD_GENE_LEVEL_PSEUDOBULK_DE_QC_FAILED"
}

DE_gate <- data.frame(dataset_id = dataset_id, input_version = input_version, pseudobulk_version = pseudobulk_version, 
    upstream_model_version = upstream_model_version, DE_version = DE_version, gate_status = gate_status, temporal_model_n = nrow(temporal_status), 
    LH7_primary_model_n = nrow(LH7_status), LH7_sensitivity_model_n = nrow(unique(LH7_sensitivity_all[, c("broad_cell_type", 
        "contrast")])), temporal_omnibus_genome_FDR_lt_0_05_n = nrow(temporal_omnibus_significant), temporal_pairwise_genome_DE_n = nrow(temporal_pairwise_significant), 
    LH7_primary_genome_DE_n = nrow(LH7_primary_significant), LH7_sensitivity_genome_DE_n = nrow(LH7_sensitivity_significant), 
    locked_program_fry_test_n = nrow(fry_all), integrated_interpretation_authorized = all_gate_checks_pass, GSE183837_independent_localization_authorized = all_gate_checks_pass, 
    composition_inference_authorized = FALSE, subtype_inference_authorized = FALSE, result_lock_authorized = FALSE, 
    stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "GENE_DE_GATE.csv")

write_csv_utf8(DE_gate, gate_path)

card_lines <- c("# STEP04E1H Gene-level Pseudobulk DE Card", "", paste0("- Dataset: ", dataset_id), paste0("- Pseudobulk version: ", 
    pseudobulk_version), paste0("- DE version: ", DE_version), paste0("- Gate: ", gate_status), paste0("- Temporal cell-type models: ", 
    nrow(temporal_status)), paste0("- LH+7 primary cell-type models: ", nrow(LH7_status)), paste0("- Temporal omnibus genes with genome-wide FDR < 0.05: ", 
    nrow(temporal_omnibus_significant)), paste0("- Temporal planned-pairwise genes meeting genome-wide FDR and effect threshold: ", 
    nrow(temporal_pairwise_significant)), paste0("- LH+7 primary genes meeting genome-wide FDR and effect threshold: ", 
    nrow(LH7_primary_significant)), paste0("- Directional locked-program fry tests: ", nrow(fry_all)), paste0("- Integrated interpretation authorized: ", 
    DE_gate$integrated_interpretation_authorized), paste0("- GSE183837 independent localization authorized: ", 
    DE_gate$GSE183837_independent_localization_authorized), "- Composition inference authorized: FALSE", "- Subtype inference authorized: FALSE", 
    "- Result lock authorized: FALSE")

card_path <- file.path(run_qc_root, "GENE_DE_RESULT_CARD.md")

write_text(card_lines, card_path)

interpretation_lines <- c("# STEP04E1H gene-level pseudobulk interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), "", "## Primary inference", "", "- Genome-wide BH FDR within each cell type and contrast is the primary gene-level multiplicity control.", 
    "- Pairwise DE reporting additionally uses an absolute log2 fold-change threshold of 0.5.", "- The temporal omnibus test asks whether any of the five fertile LH time groups differ.", 
    "- The LH+7 primary contrast compares RIF with all eligible fertile LH+7 donors.", "", "## VDR-focused inference", 
    "", "- The 436/284/46/84 VDR hierarchies were locked before the present gene-level analysis.", "- Candidate-restricted FDR values are therefore reported as focused secondary evidence.", 
    "- Candidate-restricted FDR does not replace genome-wide FDR and must be labelled explicitly.", "- Directionally weighted fry tests evaluate whether a locked gene set changes in the same direction as the calcitriol-dependent VDR interaction program.", 
    "", "## Preserved limitations", "", "- The fertile temporal series is cross-sectional.", "- Cell-composition inference remains outside the current Gate.", 
    "- Candidate subtype labels are not used for differential-expression testing.", "- Results remain candidates until content, statistical and visual QC are reviewed.")

interpretation_path <- file.path(run_qc_root, "GENE_DE_INTERPRETATION.md")

write_text(interpretation_lines, interpretation_path)

methods_lines <- c("# STEP04E1H Methods draft", "", paste("Gene-level differential-expression analyses were performed separately", 
    "for each eligible broad cell type using donor-level raw-count", "pseudobulk profiles."), paste("Genes were filtered using edgeR filterByExpr with the corresponding", 
    "cell-type-specific design matrix, followed by trimmed mean of M-values", "normalization."), paste("Negative-binomial dispersions were estimated with robust empirical", 
    "Bayes procedures, and robust quasi-likelihood generalized linear", "models were fitted."), paste("Fertile temporal effects were tested using a multi-degree-of-freedom", 
    "contrast comparing LH+5, LH+7, LH+9 and LH+11 with LH+3."), paste("Seven planned temporal pairwise contrasts were evaluated separately."), 
    paste("At LH+7, the primary comparison contrasted RIF with all fertile donors;", "sensitivity models contrasted RIF separately with standard fertile and", 
        "PGT-associated fertile donors."), paste("Benjamini-Hochberg correction was applied genome-wide within each cell", 
        "type and contrast. Pairwise genes were summarized as differentially", "expressed when FDR was below 0.05 and the absolute log2 fold-change was", 
        "at least 0.5."), paste("Because the VDR candidate hierarchy was frozen before this analysis,", "separate focused BH values were also calculated within the locked", 
        "436-, 284-, 46- and 84-gene sets; these values were interpreted as", "secondary candidate-focused evidence."), 
    paste("Directionally weighted limma fry tests were used for one-degree-of-", "freedom contrasts, with gene weights defined by the sign of the locked", 
        "calcitriol-dependent VDR interaction log-fold change."))

write_text(methods_lines, file.path(manuscript_text_root, "STEP04E1H_gene_DE_Methods_draft.md"))

results_shell <- c("# STEP04E1H Results shell", "", paste0("Gene-level pseudobulk models were completed for ", 
    nrow(temporal_status), " fertile temporal cell types and ", nrow(LH7_status), " LH+7 RIF-versus-fertile cell types."), 
    paste0("Across the five temporal omnibus models, ", nrow(temporal_omnibus_significant), " cell-type-specific gene results met genome-wide FDR < 0.05."), 
    paste0("Across the seven planned temporal contrasts, ", nrow(temporal_pairwise_significant), " gene results met genome-wide FDR < 0.05 and absolute log2 fold-change >= 0.5."), 
    paste0("In the primary LH+7 RIF-versus-fertile analyses, ", nrow(LH7_primary_significant), " gene results met the corresponding genome-wide FDR and effect-size criteria."), 
    "", "Biological interpretation will be finalized after reviewing cell-type specificity, locked VDR-candidate concordance, directional gene-set tests and sensitivity analyses.")

write_text(results_shell, file.path(manuscript_text_root, "STEP04E1H_gene_DE_Results_shell.md"))

run_end <- Sys.time()

package_names <- c("Matrix", "edgeR", "limma", "data.table", "ggplot2", "jsonlite", "scales")

package_table <- data.frame(package = package_names, version = vapply(package_names, safe_package_version, character(1)), 
    stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

artifact_files <- unique(c(list.files(output_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1H_gene_DE_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1H_gene_DE_Results_shell.md")))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, input_version = input_version, 
    pseudobulk_version = pseudobulk_version, upstream_model_version = upstream_model_version, DE_version = DE_version, 
    supersedes_partial_DE_versions = "GENEDE_v1_0;GENEDE_v1_1", resume_source_root = normalize_slashes(previous_partial_root), 
    dataset_id = dataset_id, started_at = as.character(run_start), completed_at = as.character(run_end), elapsed_minutes = round(as.numeric(difftime(run_end, 
        run_start, units = "mins")), 2), gate_status = gate_status, temporal_model_n = nrow(temporal_status), LH7_primary_model_n = nrow(LH7_status), 
    temporal_omnibus_significant_n = nrow(temporal_omnibus_significant), temporal_pairwise_significant_n = nrow(temporal_pairwise_significant), 
    LH7_primary_significant_n = nrow(LH7_primary_significant), locked_program_fry_test_n = nrow(fry_all), large_Seurat_checkpoint_loaded = FALSE, 
    composition_inference_authorized = FALSE, subtype_inference_authorized = FALSE, result_lock_authorized = FALSE, 
    output_root = normalize_slashes(output_root), qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "temporal_omnibus_significant_n=", nrow(temporal_omnibus_significant), 
    "temporal_pairwise_significant_n=", nrow(temporal_pairwise_significant), "LH7_primary_significant_n=", nrow(LH7_primary_significant))

small_handoff_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(candidate_root, 
    recursive = TRUE, full.names = TRUE), list.files(source_data_root, recursive = TRUE, full.names = TRUE), list.files(run_qc_root, 
    recursive = TRUE, full.names = TRUE), log_path, run_manifest_path, parameter_path, session_path, file.path(manuscript_text_root, 
    "STEP04E1H_gene_DE_Methods_draft.md"), file.path(manuscript_text_root, "STEP04E1H_gene_DE_Results_shell.md"), 
    normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_handoff_files <- small_handoff_files[file.exists(small_handoff_files) & !dir.exists(small_handoff_files) & 
    file.info(small_handoff_files)$size <= 150 * 1024^2]

if (!requireNamespace("zip", quietly = TRUE)) {
    try(install.packages("zip", repos = "https://cloud.r-project.org"), silent = TRUE)
}

zip_ok <- FALSE

zip_error <- ""

if (requireNamespace("zip", quietly = TRUE)) {
    zip_ok <- tryCatch({
        zip::zipr(zipfile = handoff_zip, files = small_handoff_files, recurse = FALSE, include_directories = FALSE)
        file.exists(handoff_zip)
    }, error = function(e) {
        zip_error <<- conditionMessage(e)
        FALSE
    })
}

cat("\n============================================================\n")

cat("STEP04E1H GSE250130 gene-level pseudobulk DE completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Temporal models: ", nrow(temporal_status), "\n", sep = "")

cat("Temporal omnibus FDR < 0.05 results: ", nrow(temporal_omnibus_significant), "\n", sep = "")

cat("Temporal pairwise DE results: ", nrow(temporal_pairwise_significant), "\n", sep = "")

cat("LH7 primary DE results: ", nrow(LH7_primary_significant), "\n", sep = "")

cat("Directional VDR-program fry tests: ", nrow(fry_all), "\n", sep = "")

cat("Output root: ", output_root, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe results and QC folders are complete.\n", sep = "")
}

