#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# D02: Doublet assessment and independent Harmony integration

options(stringsAsFactors = FALSE, warn = 1)

options(timeout = max(14400L, getOption("timeout")))

set.seed(20260806)

project_root <- .vdr_project_root()

dataset_id <- "GSE183837"

script_name <- "integrate_gse183837.R"

script_version <- "STEP04E1J_GSE183837_DOUBLET_HARMONY_INTEGRATION_v1.4_DONOR_GATE_FIXED"

integration_version <- "SCI183_v1_0"

run_start <- Sys.time()

run_id <- paste0("STEP04E1J_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

assert_path_safe <- function(path, max_chars = 235L) {
    path_n <- nchar(normalize_slashes(path), type = "chars")
    if (.Platform$OS.type == "windows" && path_n > max_chars) {
        stop("Path exceeds the Windows/R safety limit (", path_n, " > ", max_chars, "):\n", normalize_slashes(path))
    }
    invisible(path)
}

write_csv_utf8 <- function(x, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    assert_path_safe(path)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Could not create output directory:\n", dirname(path))
    }
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

safe_package_version <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        return(NA_character_)
    }
    as.character(utils::packageVersion(package_name))
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

detect_first_column <- function(column_names, candidates) {
    lower_names <- tolower(column_names)
    for (candidate in tolower(candidates)) {
        hit <- which(lower_names == candidate)
        if (length(hit)) {
            return(column_names[hit[1]])
        }
    }
    ""
}

safe_factor_group <- function(x) {
    values <- toupper(trimws(as.character(x)))
    values[values %in% c("CTRL", "CONTROL", "FERTILE", "NORMAL")] <- "CONTROL"
    values[values %in% c("RIF", "RECURRENT_IMPLANTATION_FAILURE")] <- "RIF"
    factor(values, levels = c("CONTROL", "RIF"))
}

get_free_memory_gb <- function() {
    if (.Platform$OS.type != "windows") {
        return(NA_real_)
    }
    command <- paste("(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory")
    output <- tryCatch(system2("powershell", args = c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, 
        stderr = FALSE), error = function(e) character())
    value <- suppressWarnings(as.numeric(trimws(output[1])))
    if (length(value) != 1L || !is.finite(value)) {
        NA_real_
    }
    else {
        round(value/1024^2, 2)
    }
}

get_free_disk_gb <- function() {
    if (.Platform$OS.type != "windows") {
        return(NA_real_)
    }
    command <- paste("(Get-CimInstance Win32_LogicalDisk", "-Filter \"DeviceID='C:'\").FreeSpace")
    output <- tryCatch(system2("powershell", args = c("-NoProfile", "-Command", shQuote(command)), stdout = TRUE, 
        stderr = FALSE), error = function(e) character())
    value <- suppressWarnings(as.numeric(trimws(output[1])))
    if (length(value) != 1L || !is.finite(value)) {
        NA_real_
    }
    else {
        round(value/1024^3, 2)
    }
}

prepare_unique_UMAP_metadata <- function(metadata, UMAP_embedding, audit_path) {
    if (is.null(rownames(metadata)) || is.null(rownames(UMAP_embedding))) {
        stop("Metadata or UMAP embedding lacks row names.")
    }
    if (!all(rownames(UMAP_embedding) %in% rownames(metadata))) {
        stop("UMAP cells are not fully represented in Seurat metadata.")
    }
    metadata_aligned <- metadata[rownames(UMAP_embedding), , drop = FALSE]
    original_column_names <- colnames(metadata_aligned)
    duplicate_name_table <- data.frame(column_name = unique(original_column_names[duplicated(original_column_names) | 
        duplicated(original_column_names, fromLast = TRUE)]), occurrence_n = vapply(unique(original_column_names[duplicated(original_column_names) | 
        duplicated(original_column_names, fromLast = TRUE)]), function(column_name) {
        sum(original_column_names == column_name)
    }, integer(1)), stringsAsFactors = FALSE)
    reserved_names <- c("cell_barcode", "UMAPH_1", "UMAPH_2")
    reserved_collision <- data.frame(reserved_name = reserved_names, present_in_original_metadata = reserved_names %in% 
        original_column_names, occurrence_n = vapply(reserved_names, function(column_name) {
        sum(original_column_names == column_name)
    }, integer(1)), stringsAsFactors = FALSE)
    metadata_aligned <- metadata_aligned[, !duplicated(colnames(metadata_aligned)), drop = FALSE]
    metadata_aligned <- metadata_aligned[, !colnames(metadata_aligned) %in% reserved_names, drop = FALSE]
    output <- data.frame(cell_barcode = rownames(UMAP_embedding), UMAPH_1 = UMAP_embedding[, 1], UMAPH_2 = UMAP_embedding[, 
        2], metadata_aligned, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(output) <- make.unique(colnames(output), sep = "__duplicate_")
    if (anyDuplicated(colnames(output)) > 0L) {
        stop("QC UMAP metadata still contains duplicate column names after repair.")
    }
    required_columns <- c("cell_barcode", "UMAPH_1", "UMAPH_2", "sample_accession", "donor_id", "analysis_group", 
        "seurat_clusters")
    missing_required <- setdiff(required_columns, colnames(output))
    if (length(missing_required)) {
        stop("QC UMAP metadata is missing required columns: ", paste(missing_required, collapse = ", "))
    }
    audit_summary <- data.frame(metric = c("original_metadata_column_n", "original_duplicate_name_n", "reserved_collision_n", 
        "export_column_n", "export_duplicate_name_n", "required_columns_complete"), value = c(length(original_column_names), 
        nrow(duplicate_name_table), sum(reserved_collision$present_in_original_metadata), ncol(output), anyDuplicated(colnames(output)), 
        !length(missing_required)), stringsAsFactors = FALSE)
    write_csv_utf8(audit_summary, audit_path)
    write_csv_utf8(duplicate_name_table, file.path(dirname(audit_path), "00_metadata_duplicate_names.csv"))
    write_csv_utf8(reserved_collision, file.path(dirname(audit_path), "00_metadata_reserved_collisions.csv"))
    output
}

join_normalized_data_layers_only <- function(object, assay_name = "RNA") {
    layers_before <- SeuratObject::Layers(object, assay = assay_name, search = NA)
    count_layers_before <- grep("^counts", layers_before, value = TRUE)
    data_layers_before <- grep("^data", layers_before, value = TRUE)
    if (!length(data_layers_before)) {
        stop("No normalized data layers were found before JoinLayers.")
    }
    if (length(data_layers_before) > 1L) {
        object <- SeuratObject::JoinLayers(object = object, assay = assay_name, layers = "data", new = "data")
    }
    layers_after <- SeuratObject::Layers(object, assay = assay_name, search = NA)
    count_layers_after <- grep("^counts", layers_after, value = TRUE)
    data_layers_after <- grep("^data", layers_after, value = TRUE)
    if (!identical(count_layers_before, count_layers_after)) {
        stop("Raw count-layer identities changed during normalized-data joining.")
    }
    if (length(data_layers_after) != 1L || data_layers_after[[1]] != "data") {
        stop("Normalized data layers did not resolve to exactly one 'data' layer.")
    }
    attr(object, "STEP04E1J_join_audit") <- list(count_layers_before = count_layers_before, count_layers_after = count_layers_after, 
        data_layers_before = data_layers_before, data_layers_after = data_layers_after)
    object
}

run_harmony_Seurat_compatible <- function(object, group_variable, dimensions = 1:30, reduction_name = "pca", output_reduction = "harmony", 
    verbose = FALSE) {
    if (!requireNamespace("harmony", quietly = TRUE)) {
        stop("The harmony package is not installed.")
    }
    harmony_method <- tryCatch(utils::getS3method("RunHarmony", "Seurat", envir = asNamespace("harmony")), error = function(e) NULL)
    method_formals <- if (is.null(harmony_method)) {
        character()
    }
    else {
        names(formals(harmony_method))
    }
    call_arguments <- list(object = object, group.by.vars = group_variable, reduction = reduction_name, dims.use = dimensions, 
        reduction.save = output_reduction, plot_convergence = FALSE, verbose = verbose)
    output <- do.call(harmony::RunHarmony, call_arguments)
    attr(output, "STEP04E1J_Harmony_method_formals") <- method_formals
    output
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

required_packages <- c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "ggplot2", "jsonlite", "scDblFinder", 
    "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "BiocParallel")

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

design_path <- file.path(results_root, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0", "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv")

preQC_root <- file.path(results_root, "S4D_sc_preQC", "SCP_v1_0", dataset_id)

output_root <- file.path(results_root, "S4E1_sc", integration_version)

sample_checkpoint_root <- file.path(output_root, "sample_ckpt")

checkpoint_root <- file.path(output_root, "checkpoints", dataset_id)

table_root <- file.path(output_root, "tables")

manifest_root <- file.path(output_root, "manifests")

merged_counts_checkpoint <- file.path(checkpoint_root, "merged_singlets.rds")

prejoin_checkpoint_path <- file.path(checkpoint_root, "harmony_prejoin.rds")

final_checkpoint_path <- file.path(checkpoint_root, "integrated.rds")

final_status_path <- file.path(checkpoint_root, "status.csv")

run_qc_root <- file.path(qc_base_root, paste0("STEP04E1J_", format(run_start, "%Y%m%d_%H%M%S"), "_integration"))

figure_root <- file.path(run_qc_root, "figures")

for (path_value in c(execution_root, parameter_root, environment_root, output_root, sample_checkpoint_root, checkpoint_root, 
    table_root, manifest_root, run_qc_root, figure_root)) {
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

planned_paths <- data.frame(artifact = c("output_root", "sample_checkpoint_root", "checkpoint_root", "merged_counts_checkpoint", 
    "prejoin_checkpoint", "final_checkpoint", "final_status", "table_root", "manifest_root", "run_qc_root", "handoff_zip"), 
    path = normalize_slashes(c(output_root, sample_checkpoint_root, checkpoint_root, merged_counts_checkpoint, 
        prejoin_checkpoint_path, final_checkpoint_path, final_status_path, table_root, manifest_root, run_qc_root, 
        handoff_zip)), stringsAsFactors = FALSE)

planned_paths$character_n <- nchar(planned_paths$path, type = "chars")

planned_paths$path_safe <- (planned_paths$character_n <= 235L)

if (!all(planned_paths$path_safe)) {
    stop("One or more planned paths exceed the safety limit:\n", paste(planned_paths$path[!planned_paths$path_safe], 
        collapse = "\n"))
}

write_csv_utf8(planned_paths, file.path(run_qc_root, "00_path_audit.csv"))

preflight_qc_dirs <- list.dirs(qc_base_root, recursive = FALSE, full.names = TRUE)

preflight_qc_dirs <- preflight_qc_dirs[grepl("^STEP04E1I_[0-9]{8}_[0-9]{6}_preflight$", basename(preflight_qc_dirs), 
    perl = TRUE)]

if (!length(preflight_qc_dirs)) {
    stop("No completed STEP04E1I preflight-QC directory found.")
}

latest_preflight_qc <- preflight_qc_dirs[which.max(file.info(preflight_qc_dirs)$mtime)]

preflight_gate_path <- file.path(latest_preflight_qc, "GSE183837_PREFLIGHT_GATE.csv")

preflight_gate <- safe_read_csv(preflight_gate_path)

if (is.null(preflight_gate) || nrow(preflight_gate) != 1L || preflight_gate$gate_status[[1]] != "PASS_GSE183837_INTEGRATION_PREFLIGHT" || 
    !isTRUE(preflight_gate$integration_authorized[[1]])) {
    stop("STEP04E1I Gate does not authorize GSE183837 integration.")
}

if (!file.exists(design_path)) {
    stop("Locked sample design is missing:\n", design_path)
}

design <- safe_read_csv(design_path)

if (is.null(design) || nrow(design) != 9L || !"sample_accession" %in% colnames(design) || anyDuplicated(design$sample_accession) > 
    0L) {
    stop("Locked GSE183837 design failed identity validation.")
}

group_column <- detect_first_column(colnames(design), c("group", "clinical_group", "condition", "phenotype", "case_control_status", 
    "case_control"))

donor_column <- detect_first_column(colnames(design), c("donor_id", "donor", "subject_id", "patient_id"))

if (!nzchar(group_column) || !nzchar(donor_column)) {
    stop("Could not identify group and donor columns in the locked design.")
}

design$analysis_group <- safe_factor_group(design[[group_column]])

design$analysis_donor_id <- as.character(design[[donor_column]])

if (anyNA(design$analysis_group) || anyNA(design$analysis_donor_id) || anyDuplicated(design$analysis_donor_id) > 
    0L || !identical(as.integer(table(design$analysis_group)), c(3L, 6L))) {
    stop("Locked GSE183837 donor/group design is not 3 Control and 6 RIF.")
}

free_memory_before_GB <- get_free_memory_gb()

free_disk_before_GB <- get_free_disk_gb()

if (is.finite(free_memory_before_GB) && free_memory_before_GB < 12) {
    stop("Free physical memory is below 12 GB: ", free_memory_before_GB)
}

if (is.finite(free_disk_before_GB) && free_disk_before_GB < 15) {
    stop("Free C-drive space is below 15 GB: ", free_disk_before_GB)
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, integration_version = integration_version, 
    dataset_id = dataset_id, upstream_gate = preflight_gate$gate_status[[1]], expected_sample_n = 9L, expected_Control_n = 3L, 
    expected_RIF_n = 6L, expected_preQC_cell_n = 50997L, doublet_method = "scDblFinder", doublet_scope = "INDEPENDENT_WITHIN_EACH_DONOR_CAPTURE", 
    scDblFinder_clusters = TRUE, scDblFinder_parallel = "BiocParallel_SerialParam", normalization = "LogNormalize_scale_factor_10000", 
    variable_features = 3000L, PCA_components = 50L, Harmony_batch_variable = "sample_accession", Harmony_dimensions = "1:30", 
    Harmony_API_policy = "MINIMAL_CURRENT_API_NO_ASSAY_USE_NO_PROJECT_DIM", JoinLayers_policy = "LAYERS_DATA_SEARCH_CATEGORY_TO_SINGLE_DATA_LAYER_COUNTS_UNCHANGED", 
    QC_metadata_policy = "DROP_DUPLICATE_OR_RESERVED_METADATA_NAMES_BEFORE_UMAP_EXPORT", donor_gate_policy = "COUNT_UNIQUE_DONOR_GROUP_PAIRS_NOT_CELLS", 
    supersedes_failed_executions = c("V1_0_ASSAY_USE_REJECTED_BY_HARMONY_CHECK_LEGACY_ARGS", "V1_1_EXACT_DATA_LAYER_VECTOR_WITH_SINGLE_NEW_NAME", 
        "V1_2_GGPLOT_DUPLICATE_QC_METADATA_COLUMNS", "V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS"), UMAP_dimensions = "1:30", 
    cluster_resolutions = c(0.20000000000000001, 0.40000000000000002, 0.59999999999999998), primary_resolution = 0.40000000000000002, 
    counts_layer_policy = "RETAIN_NINE_SPLIT_RAW_COUNT_LAYERS", data_layer_policy = "JOIN_NORMALIZED_DATA_LAYERS_ONLY", 
    scale_data_policy = "REMOVE_BEFORE_FINAL_CHECKPOINT_WHEN_SUPPORTED", marker_analysis_authorized = FALSE, disease_testing_authorized = FALSE, 
    result_lock_authorized = FALSE, random_seed = 20260806L)

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

harmony_method <- tryCatch(utils::getS3method("RunHarmony", "Seurat", envir = asNamespace("harmony")), error = function(e) NULL)

harmony_API_audit <- data.frame(metric = c("harmony_version", "Seurat_method_found", "Seurat_method_formals", "assay_use_passed_by_script", 
    "project_dim_passed_by_script", "reduction", "dims_use", "reduction_save"), value = c(safe_package_version("harmony"), 
    !is.null(harmony_method), if (is.null(harmony_method)) {
        ""
    } else {
        paste(names(formals(harmony_method)), collapse = ";")
    }, FALSE, FALSE, "pca", "1:30", "harmony"), stringsAsFactors = FALSE)

write_csv_utf8(harmony_API_audit, file.path(run_qc_root, "00_Harmony_API_audit.csv"))

log_event("Starting", run_id, "free_memory_GB=", free_memory_before_GB, "free_disk_GB=", free_disk_before_GB)

reference_features <- NULL

sample_status_rows <- list()

doublet_metadata_paths <- character()

for (sample_index in seq_len(nrow(design))) {
    sample_accession <- as.character(design$sample_accession[sample_index])
    donor_id <- as.character(design$analysis_donor_id[sample_index])
    analysis_group <- as.character(design$analysis_group[sample_index])
    sample_dir <- file.path(sample_checkpoint_root, sample_accession)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    sample_checkpoint_path <- file.path(sample_dir, "singlet.rds")
    sample_status_path <- file.path(sample_dir, "status.csv")
    sample_doublet_path <- file.path(sample_dir, "doublet_metadata.csv")
    existing_status <- safe_read_csv(sample_status_path)
    if (file.exists(sample_checkpoint_path) && !is.null(existing_status) && nrow(existing_status) == 1L && existing_status$status[[1]] == 
        "COMPLETE") {
        checkpoint_object <- tryCatch(readRDS(sample_checkpoint_path), error = function(e) NULL)
        checkpoint_valid <- (inherits(checkpoint_object, "Seurat") && ncol(checkpoint_object) == existing_status$singlet_cell_n[[1]] && 
            length(unique(as.character(checkpoint_object$sample_accession))) == 1L && unique(as.character(checkpoint_object$sample_accession)) == 
            sample_accession)
        if (checkpoint_valid) {
            sample_status_rows[[sample_index]] <- existing_status
            doublet_metadata_paths <- c(doublet_metadata_paths, sample_doublet_path)
            if (is.null(reference_features)) {
                reference_features <- rownames(checkpoint_object)
            }
            else if (!identical(reference_features, rownames(checkpoint_object))) {
                stop("Feature order differs in existing sample checkpoint: ", sample_accession)
            }
            log_event("SKIP completed donor", sample_index, "/", nrow(design), sample_accession, "singlets=", ncol(checkpoint_object))
            rm(checkpoint_object)
            gc()
            next
        }
    }
    sample_start <- Sys.time()
    input_path <- file.path(preQC_root, paste0(sample_accession, "_preQC.rds"))
    if (!file.exists(input_path)) {
        stop("Pre-QC object is missing:\n", input_path)
    }
    input_object <- readRDS(input_path)
    required_components <- c("counts", "feature_table", "cell_metadata")
    if (!all(required_components %in% names(input_object))) {
        stop("Pre-QC object lacks required components: ", sample_accession)
    }
    counts <- input_object$counts
    cell_metadata <- input_object$cell_metadata
    if ((!inherits(counts, "Matrix") && !is.matrix(counts)) || nrow(counts) != 27984L || ncol(counts) < 1000L || 
        is.null(rownames(counts)) || is.null(colnames(counts)) || anyDuplicated(colnames(counts)) > 0L) {
        stop("Count-matrix identity failed for ", sample_accession)
    }
    if (is.null(reference_features)) {
        reference_features <- rownames(counts)
    }
    else if (identical(reference_features, rownames(counts))) {
        invisible(TRUE)
    }
    else if (setequal(reference_features, rownames(counts))) {
        counts <- counts[reference_features, , drop = FALSE]
    }
    else {
        stop("Feature identities differ across donors at ", sample_accession)
    }
    if (!is.data.frame(cell_metadata) || nrow(cell_metadata) != ncol(counts)) {
        stop("Cell metadata does not align with counts for ", sample_accession)
    }
    if (!is.null(rownames(cell_metadata)) && setequal(rownames(cell_metadata), colnames(counts))) {
        cell_metadata <- cell_metadata[colnames(counts), , drop = FALSE]
    }
    else {
        rownames(cell_metadata) <- colnames(counts)
    }
    original_cell_barcode <- colnames(counts)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    SummarizedExperiment::colData(sce)$sample_accession <- sample_accession
    set.seed(20260806L + sample_index)
    sce <- scDblFinder::scDblFinder(sce = sce, clusters = TRUE, samples = NULL, verbose = FALSE, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE))
    scDblFinder_metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    if (!"scDblFinder.class" %in% colnames(scDblFinder_metadata) || !"scDblFinder.score" %in% colnames(scDblFinder_metadata)) {
        stop("scDblFinder output columns are missing for ", sample_accession)
    }
    predicted_class <- tolower(as.character(scDblFinder_metadata$scDblFinder.class))
    singlet_index <- which(predicted_class == "singlet")
    doublet_index <- which(predicted_class == "doublet")
    input_cell_n <- ncol(counts)
    singlet_cell_n <- length(singlet_index)
    doublet_cell_n <- length(doublet_index)
    if (singlet_cell_n < 1000L || singlet_cell_n + doublet_cell_n != input_cell_n) {
        stop("Invalid singlet/doublet partition for ", sample_accession)
    }
    doublet_fraction <- (doublet_cell_n/input_cell_n)
    doublet_metadata <- data.frame(original_cell_barcode = original_cell_barcode, integrated_cell_barcode = paste0(sample_accession, 
        "__", original_cell_barcode), sample_accession = sample_accession, donor_id = donor_id, analysis_group = analysis_group, 
        scDblFinder_class = predicted_class, scDblFinder_score = as.numeric(scDblFinder_metadata$scDblFinder.score), 
        stringsAsFactors = FALSE)
    optional_scDblFinder_columns <- grep("^scDblFinder\\.", colnames(scDblFinder_metadata), value = TRUE)
    optional_scDblFinder_columns <- setdiff(optional_scDblFinder_columns, c("scDblFinder.class", "scDblFinder.score"))
    for (column_name in optional_scDblFinder_columns) {
        doublet_metadata[[column_name]] <- scDblFinder_metadata[[column_name]]
    }
    write_csv_utf8(doublet_metadata, sample_doublet_path)
    singlet_counts <- counts[, singlet_index, drop = FALSE]
    singlet_metadata <- cell_metadata[singlet_index, , drop = FALSE]
    new_cell_names <- paste0(sample_accession, "__", original_cell_barcode[singlet_index])
    colnames(singlet_counts) <- new_cell_names
    rownames(singlet_metadata) <- new_cell_names
    singlet_metadata$original_cell_barcode <- original_cell_barcode[singlet_index]
    singlet_metadata$sample_accession <- sample_accession
    singlet_metadata$donor_id <- donor_id
    singlet_metadata$analysis_group <- analysis_group
    singlet_metadata$scDblFinder_class <- predicted_class[singlet_index]
    singlet_metadata$scDblFinder_score <- as.numeric(scDblFinder_metadata$scDblFinder.score[singlet_index])
    design_row <- design[design$sample_accession == sample_accession, , drop = FALSE]
    for (design_column in colnames(design_row)) {
        if (!design_column %in% colnames(singlet_metadata)) {
            singlet_metadata[[design_column]] <- rep(as.character(design_row[[design_column]][1]), singlet_cell_n)
        }
    }
    singlet_object <- Seurat::CreateSeuratObject(counts = singlet_counts, project = sample_accession, assay = "RNA", 
        min.cells = 0, min.features = 0, meta.data = singlet_metadata)
    SeuratObject::Project(singlet_object) <- sample_accession
    if (any(grepl("^MT-", rownames(singlet_object), ignore.case = FALSE))) {
        singlet_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(singlet_object, pattern = "^MT-")
    }
    singlet_object@misc$GSE183837_sample_checkpoint <- list(sample_accession = sample_accession, donor_id = donor_id, 
        analysis_group = analysis_group, input_cell_n = input_cell_n, singlet_cell_n = singlet_cell_n, doublet_cell_n = doublet_cell_n, 
        doublet_fraction = doublet_fraction, scDblFinder_version = safe_package_version("scDblFinder"), random_seed = 20260806L + 
            sample_index)
    saveRDS(singlet_object, sample_checkpoint_path, compress = "gzip")
    checkpoint_md5 <- unname(tools::md5sum(sample_checkpoint_path))
    elapsed_minutes <- round(as.numeric(difftime(Sys.time(), sample_start, units = "mins")), 2)
    sample_status <- data.frame(sample_accession = sample_accession, donor_id = donor_id, analysis_group = analysis_group, 
        status = "COMPLETE", input_cell_n = input_cell_n, singlet_cell_n = singlet_cell_n, doublet_cell_n = doublet_cell_n, 
        doublet_fraction = doublet_fraction, high_doublet_review_flag = doublet_fraction > 0.20000000000000001, 
        singlet_checkpoint_size_MB = round(file.info(sample_checkpoint_path)$size/1024^2, 2), singlet_checkpoint_md5 = checkpoint_md5, 
        elapsed_minutes = elapsed_minutes, stringsAsFactors = FALSE)
    write_csv_utf8(sample_status, sample_status_path)
    sample_status_rows[[sample_index]] <- sample_status
    doublet_metadata_paths <- c(doublet_metadata_paths, sample_doublet_path)
    log_event("DONE donor", sample_index, "/", nrow(design), sample_accession, "input=", input_cell_n, "doublets=", 
        doublet_cell_n, "singlets=", singlet_cell_n, "doublet_fraction=", round(doublet_fraction, 4), "elapsed_min=", 
        elapsed_minutes)
    rm(input_object, counts, cell_metadata, sce, scDblFinder_metadata, doublet_metadata, singlet_counts, singlet_metadata, 
        singlet_object)
    gc()
}

sample_status_all <- do.call(rbind, sample_status_rows)

sample_status_all <- sample_status_all[match(design$sample_accession, sample_status_all$sample_accession), , drop = FALSE]

if (nrow(sample_status_all) != 9L || anyNA(sample_status_all$sample_accession) || !all(sample_status_all$status == 
    "COMPLETE") || sum(sample_status_all$input_cell_n) != 50997L) {
    stop("Sample-wise scDblFinder completion audit failed.")
}

write_csv_utf8(sample_status_all, file.path(table_root, "01_doublet_summary.csv"))

doublet_metadata_all <- data.table::rbindlist(lapply(doublet_metadata_paths, function(path_value) {
    .vdr_fread(path_value, data.table = FALSE, showProgress = FALSE)
}), use.names = TRUE, fill = TRUE)

write_csv_utf8(as.data.frame(doublet_metadata_all), file.path(table_root, "02_doublet_metadata_all.csv"))

total_singlets <- sum(sample_status_all$singlet_cell_n)

total_doublets <- sum(sample_status_all$doublet_cell_n)

if (file.exists(final_status_path)) {
    final_status_existing <- safe_read_csv(final_status_path)
} else {
    final_status_existing <- NULL
}

final_checkpoint_complete <- (file.exists(final_checkpoint_path) && !is.null(final_status_existing) && nrow(final_status_existing) == 
    1L && final_status_existing$status[[1]] == "COMPLETE")

if (final_checkpoint_complete) {
    integrated_object <- readRDS(final_checkpoint_path)
    if (!inherits(integrated_object, "Seurat") || ncol(integrated_object) != total_singlets || length(unique(integrated_object$sample_accession)) != 
        9L) {
        stop("Existing final checkpoint failed validation.")
    }
    log_event("SKIP integration; valid final checkpoint found", "cells=", ncol(integrated_object))
} else {
    if (file.exists(merged_counts_checkpoint)) {
        merged_object <- tryCatch(readRDS(merged_counts_checkpoint), error = function(e) NULL)
        merged_checkpoint_valid <- (inherits(merged_object, "Seurat") && ncol(merged_object) == total_singlets && 
            length(unique(merged_object$sample_accession)) == 9L)
        if (!merged_checkpoint_valid) {
            rm(merged_object)
            merged_object <- NULL
        }
    }
    else {
        merged_object <- NULL
    }
    if (is.null(merged_object)) {
        log_event("START merging nine singlet checkpoints")
        object_list <- vector("list", nrow(design))
        names(object_list) <- as.character(design$sample_accession)
        for (sample_index in seq_len(nrow(design))) {
            sample_accession <- as.character(design$sample_accession[sample_index])
            checkpoint_path <- file.path(sample_checkpoint_root, sample_accession, "singlet.rds")
            object_list[[sample_index]] <- readRDS(checkpoint_path)
        }
        merged_object <- merge(x = object_list[[1]], y = object_list[-1], merge.data = FALSE, merge.dr = FALSE, 
            collapse = FALSE, project = dataset_id)
        if (ncol(merged_object) != total_singlets || length(unique(merged_object$sample_accession)) != 9L) {
            stop("Merged singlet checkpoint failed cell/sample identity.")
        }
        saveRDS(merged_object, merged_counts_checkpoint, compress = "gzip")
        log_event("DONE merging singlets", "cells=", ncol(merged_object), "features=", nrow(merged_object), "checkpoint_MB=", 
            round(file.info(merged_counts_checkpoint)$size/1024^2, 2))
        rm(object_list)
        gc()
    }
    else {
        log_event("RESUME from merged singlet checkpoint", "cells=", ncol(merged_object))
    }
    first_preQC_path <- file.path(preQC_root, paste0(as.character(design$sample_accession[1]), "_preQC.rds"))
    first_preQC_object <- readRDS(first_preQC_path)
    original_features <- rownames(first_preQC_object$counts)
    expected_Seurat_features <- gsub("_", "-", original_features, fixed = TRUE)
    expected_Seurat_features <- make.unique(expected_Seurat_features)
    current_features <- rownames(merged_object)
    feature_name_audit <- data.frame(metric = c("original_feature_n", "merged_feature_n", "original_feature_names_with_underscore_n", 
        "expected_transformed_feature_unique_n", "merged_feature_unique_n", "expected_transformation_matches_merged"), 
        value = c(length(original_features), length(current_features), sum(grepl("_", original_features, fixed = TRUE)), 
            length(unique(expected_Seurat_features)), length(unique(current_features)), identical(expected_Seurat_features, 
                current_features)), stringsAsFactors = FALSE)
    write_csv_utf8(feature_name_audit, file.path(table_root, "00_feature_name_audit.csv"))
    if (length(current_features) != length(original_features) || anyDuplicated(current_features) > 0L) {
        stop("Feature-name transformation produced a dimension or uniqueness problem.")
    }
    rm(first_preQC_object, original_features, expected_Seurat_features, current_features)
    gc()
    prejoin_checkpoint_valid <- FALSE
    prejoin_object <- NULL
    if (file.exists(prejoin_checkpoint_path)) {
        prejoin_object <- tryCatch(readRDS(prejoin_checkpoint_path), error = function(e) NULL)
        if (inherits(prejoin_object, "Seurat") && ncol(prejoin_object) == total_singlets && length(unique(prejoin_object$sample_accession)) == 
            9L && all(c("pca", "harmony", "umap.harmony") %in% names(prejoin_object@reductions)) && "harmony_snn_res.0.4" %in% 
            colnames(prejoin_object[[]])) {
            prejoin_checkpoint_valid <- TRUE
        }
    }
    integration_start <- Sys.time()
    if (prejoin_checkpoint_valid) {
        merged_object <- prejoin_object
        primary_cluster_column <- "harmony_snn_res.0.4"
        merged_object$seurat_clusters <- as.factor(merged_object[[primary_cluster_column]][, 1])
        SeuratObject::Idents(merged_object) <- merged_object$seurat_clusters
        log_event("RESUME from Harmony prejoin checkpoint", "cells=", ncol(merged_object), "clusters_res0.4=", 
            length(unique(merged_object$seurat_clusters)))
    }
    else {
        rm(prejoin_object)
        gc()
        SeuratObject::DefaultAssay(merged_object) <- "RNA"
        merged_object <- Seurat::NormalizeData(object = merged_object, assay = "RNA", normalization.method = "LogNormalize", 
            scale.factor = 10000, verbose = FALSE)
        merged_object <- Seurat::FindVariableFeatures(object = merged_object, assay = "RNA", selection.method = "vst", 
            nfeatures = 3000, verbose = FALSE)
        if (length(SeuratObject::VariableFeatures(merged_object)) < 2000L) {
            stop("Too few variable features were recovered.")
        }
        merged_object <- Seurat::ScaleData(object = merged_object, assay = "RNA", features = SeuratObject::VariableFeatures(merged_object), 
            verbose = FALSE)
        merged_object <- Seurat::RunPCA(object = merged_object, assay = "RNA", features = SeuratObject::VariableFeatures(merged_object), 
            npcs = 50, seed.use = 20260806, verbose = FALSE)
        merged_object <- run_harmony_Seurat_compatible(object = merged_object, group_variable = "sample_accession", 
            dimensions = 1:30, reduction_name = "pca", output_reduction = "harmony", verbose = FALSE)
        merged_object <- Seurat::RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony", 
            reduction.key = "UMAPH_", seed.use = 20260806, verbose = FALSE)
        merged_object <- Seurat::FindNeighbors(object = merged_object, reduction = "harmony", dims = 1:30, graph.name = c("harmony_nn", 
            "harmony_snn"), verbose = FALSE)
        merged_object <- Seurat::FindClusters(object = merged_object, graph.name = "harmony_snn", resolution = c(0.20000000000000001, 
            0.40000000000000002, 0.59999999999999998), algorithm = 1, random.seed = 20260806, verbose = FALSE)
        primary_cluster_column <- "harmony_snn_res.0.4"
        if (!primary_cluster_column %in% colnames(merged_object[[]])) {
            candidate_primary_columns <- grep("res\\.0\\.4$", colnames(merged_object[[]]), value = TRUE)
            if (length(candidate_primary_columns) != 1L) {
                stop("Could not resolve the primary 0.4 cluster column.")
            }
            primary_cluster_column <- candidate_primary_columns[1]
        }
        merged_object$seurat_clusters <- as.factor(merged_object[[primary_cluster_column]][, 1])
        SeuratObject::Idents(merged_object) <- merged_object$seurat_clusters
        merged_object@misc$GSE183837_prejoin_checkpoint <- list(script_version = script_version, cell_n = ncol(merged_object), 
            sample_n = length(unique(merged_object$sample_accession)), primary_cluster_column = primary_cluster_column, 
            primary_cluster_n = length(unique(merged_object$seurat_clusters)), reductions = names(merged_object@reductions), 
            checkpoint_role = "POST_HARMONY_UMAP_CLUSTERING_PRE_JOINLAYERS")
        saveRDS(merged_object, prejoin_checkpoint_path, compress = "gzip")
        log_event("SAVED Harmony prejoin checkpoint", "checkpoint_MB=", round(file.info(prejoin_checkpoint_path)$size/1024^2, 
            2))
    }
    data_layers_before_join <- SeuratObject::Layers(merged_object, assay = "RNA", search = "^data")
    count_layers_before_join <- SeuratObject::Layers(merged_object, assay = "RNA", search = "^counts")
    merged_object <- join_normalized_data_layers_only(object = merged_object, assay_name = "RNA")
    data_layers_after_join <- SeuratObject::Layers(merged_object, assay = "RNA", search = "^data")
    count_layers_after_join <- SeuratObject::Layers(merged_object, assay = "RNA", search = "^counts")
    join_audit <- data.frame(metric = c("data_layers_before_n", "data_layers_before", "data_layers_after_n", "data_layers_after", 
        "count_layers_before_n", "count_layers_before", "count_layers_after_n", "count_layers_after", "count_layer_identity_preserved"), 
        value = c(length(data_layers_before_join), collapse_unique(data_layers_before_join), length(data_layers_after_join), 
            collapse_unique(data_layers_after_join), length(count_layers_before_join), collapse_unique(count_layers_before_join), 
            length(count_layers_after_join), collapse_unique(count_layers_after_join), identical(count_layers_before_join, 
                count_layers_after_join)), stringsAsFactors = FALSE)
    write_csv_utf8(join_audit, file.path(table_root, "00_JoinLayers_audit.csv"))
    scale_layers_before <- SeuratObject::Layers(merged_object, assay = "RNA", search = "^scale\\.data")
    scale_data_removed <- FALSE
    if ("scale.data" %in% scale_layers_before) {
        assay_object <- merged_object[["RNA"]]
        scale_data_removed <- tryCatch({
            assay_object$scale.data <- NULL
            merged_object[["RNA"]] <- assay_object
            TRUE
        }, error = function(e) FALSE)
    }
    integrated_object <- merged_object
    integrated_object@misc$GSE183837_integration <- list(integration_version = integration_version, script_version = script_version, 
        Harmony_API_policy = "MINIMAL_CURRENT_API_NO_ASSAY_USE_NO_PROJECT_DIM", JoinLayers_policy = "LAYERS_DATA_TO_SINGLE_DATA_COUNTS_UNCHANGED", 
        QC_metadata_policy = "UNIQUE_COLUMN_NAMES_FOR_EXPORT_AND_GGPLOT", donor_gate_policy = "UNIQUE_DONOR_GROUP_PAIRS", 
        supersedes_failed_executions = c("V1_0_HARMONY_ASSAY_USE_FAILURE", "V1_1_JOINLAYERS_LENGTH_MISMATCH", "V1_2_GGPLOT_DUPLICATE_METADATA_COLUMNS", 
            "V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS"), upstream_preflight_gate = preflight_gate$gate_status[[1]], 
        sample_n = 9L, preQC_cell_n = 50997L, singlet_cell_n = total_singlets, predicted_doublet_n = total_doublets, 
        Harmony_batch_variable = "sample_accession", Harmony_dimensions = 1:30, cluster_resolutions = c(0.20000000000000001, 
            0.40000000000000002, 0.59999999999999998), primary_resolution = 0.40000000000000002, primary_cluster_column = primary_cluster_column, 
        scale_data_removed = scale_data_removed, marker_analysis_authorized = FALSE, result_lock_authorized = FALSE)
    gc()
    saveRDS(integrated_object, final_checkpoint_path, compress = "gzip")
    final_checkpoint_md5 <- unname(tools::md5sum(final_checkpoint_path))
    integration_elapsed_minutes <- round(as.numeric(difftime(Sys.time(), integration_start, units = "mins")), 2)
    final_status <- data.frame(dataset_id = dataset_id, integration_version = integration_version, status = "COMPLETE", 
        feature_n = nrow(integrated_object), singlet_cell_n = ncol(integrated_object), sample_n = length(unique(integrated_object$sample_accession)), 
        donor_n = length(unique(integrated_object$donor_id)), group_n = length(unique(integrated_object$analysis_group)), 
        primary_cluster_column = primary_cluster_column, primary_cluster_n = length(unique(integrated_object$seurat_clusters)), 
        final_checkpoint_size_GB = round(file.info(final_checkpoint_path)$size/1024^3, 3), final_checkpoint_md5 = final_checkpoint_md5, 
        integration_elapsed_minutes = integration_elapsed_minutes, stringsAsFactors = FALSE)
    write_csv_utf8(final_status, final_status_path)
    log_event("DONE Harmony integration", "cells=", ncol(integrated_object), "clusters_res0.4=", length(unique(integrated_object$seurat_clusters)), 
        "elapsed_min=", integration_elapsed_minutes)
    rm(merged_object)
    gc()
}

if (!exists("primary_cluster_column")) {
    primary_cluster_column <- integrated_object@misc$GSE183837_integration$primary_cluster_column
}

RNA_layers <- SeuratObject::Layers(integrated_object, assay = "RNA", search = NA)

count_layers <- grep("^counts", RNA_layers, value = TRUE)

data_layers_final <- grep("^data", RNA_layers, value = TRUE)

scale_layers_final <- grep("^scale\\.data", RNA_layers, value = TRUE)

layer_audit <- data.frame(layer = RNA_layers, layer_type = ifelse(grepl("^counts", RNA_layers), "counts", ifelse(grepl("^data", 
    RNA_layers), "data", ifelse(grepl("^scale\\.data", RNA_layers), "scale.data", "other"))), stringsAsFactors = FALSE)

write_csv_utf8(layer_audit, file.path(table_root, "03_layer_audit.csv"))

metadata <- integrated_object[[]]

UMAP_embedding <- SeuratObject::Embeddings(integrated_object, reduction = "umap.harmony")

UMAP_metadata <- prepare_unique_UMAP_metadata(metadata = metadata, UMAP_embedding = UMAP_embedding, audit_path = file.path(table_root, 
    "00_metadata_name_audit.csv"))

write_csv_utf8(UMAP_metadata, file.path(table_root, "04_UMAP_metadata.csv"))

sample_cell_summary <- aggregate(cell_barcode ~ sample_accession + donor_id + analysis_group, data = UMAP_metadata, 
    FUN = length)

colnames(sample_cell_summary)[colnames(sample_cell_summary) == "cell_barcode"] <- "singlet_cell_n"

sample_cell_summary <- merge(sample_cell_summary, sample_status_all[, c("sample_accession", "input_cell_n", "doublet_cell_n", 
    "doublet_fraction")], by = "sample_accession", all.x = TRUE, sort = FALSE)

sample_cell_summary <- sample_cell_summary[match(design$sample_accession, sample_cell_summary$sample_accession), 
    , drop = FALSE]

write_csv_utf8(sample_cell_summary, file.path(table_root, "05_sample_cell_summary.csv"))

cluster_sample_counts <- aggregate(cell_barcode ~ seurat_clusters + sample_accession + analysis_group, data = UMAP_metadata, 
    FUN = length)

colnames(cluster_sample_counts)[colnames(cluster_sample_counts) == "cell_barcode"] <- "cell_n"

write_csv_utf8(cluster_sample_counts, file.path(table_root, "06_cluster_sample_counts.csv"))

cluster_summary <- aggregate(cell_barcode ~ seurat_clusters, data = UMAP_metadata, FUN = length)

colnames(cluster_summary)[colnames(cluster_summary) == "cell_barcode"] <- "cell_n"

cluster_summary$cell_fraction <- (cluster_summary$cell_n/sum(cluster_summary$cell_n))

cluster_summary <- cluster_summary[order(as.integer(as.character(cluster_summary$seurat_clusters))), , drop = FALSE]

write_csv_utf8(cluster_summary, file.path(table_root, "07_cluster_summary.csv"))

resolution_columns <- grep("harmony_snn_res\\.", colnames(metadata), value = TRUE)

resolution_summary_rows <- list()

for (resolution_column in resolution_columns) {
    resolution_summary_rows[[length(resolution_summary_rows) + 1L]] <- data.frame(resolution_column = resolution_column, 
        cluster_n = length(unique(metadata[[resolution_column]])), minimum_cluster_n = min(table(metadata[[resolution_column]])), 
        maximum_cluster_n = max(table(metadata[[resolution_column]])), stringsAsFactors = FALSE)
}

resolution_summary <- do.call(rbind, resolution_summary_rows)

write_csv_utf8(resolution_summary, file.path(table_root, "08_resolution_summary.csv"))

PCA_embedding <- SeuratObject::Embeddings(integrated_object, reduction = "pca")

Harmony_embedding <- SeuratObject::Embeddings(integrated_object, reduction = "harmony")

reduction_audit <- data.frame(reduction = c("pca", "harmony", "umap.harmony"), cell_n = c(nrow(PCA_embedding), 
    nrow(Harmony_embedding), nrow(UMAP_embedding)), dimension_n = c(ncol(PCA_embedding), ncol(Harmony_embedding), 
    ncol(UMAP_embedding)), stringsAsFactors = FALSE)

write_csv_utf8(reduction_audit, file.path(table_root, "09_reduction_audit.csv"))

if (anyDuplicated(colnames(UMAP_metadata)) > 0L) {
    stop("UMAP plotting metadata contains duplicate column names.")
}

theme_QC <- ggplot2::theme_bw(base_size = 10) + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

sample_status_plot <- sample_status_all[, !duplicated(colnames(sample_status_all)), drop = FALSE]

if (anyDuplicated(colnames(sample_status_plot)) > 0L) {
    stop("Sample-status plotting data contains duplicate column names.")
}

sample_status_plot$sample_accession <- factor(sample_status_plot$sample_accession, levels = design$sample_accession)

plot_doublet <- ggplot2::ggplot(sample_status_plot, ggplot2::aes(x = sample_accession, y = doublet_fraction, shape = analysis_group)) + 
    ggplot2::geom_hline(yintercept = 0.20000000000000001, linetype = "dashed") + ggplot2::geom_point(size = 2.5) + 
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + ggplot2::labs(title = "GSE183837 donor-wise predicted doublet fractions", 
    subtitle = "scDblFinder was run independently within each donor; dashed line is a review threshold", x = "Sample", 
    y = "Predicted doublet fraction", shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_doublet, "QC1_doublet_fraction", width = 9, height = 5)

plot_singlets <- ggplot2::ggplot(sample_status_plot, ggplot2::aes(x = sample_accession, y = singlet_cell_n, shape = analysis_group)) + 
    ggplot2::geom_point(size = 2.5) + ggplot2::scale_y_continuous(labels = scales::comma) + ggplot2::labs(title = "GSE183837 retained singlet cells", 
    x = "Sample", y = "Retained singlet cells", shape = "Group") + theme_QC + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
    hjust = 1))

save_plot_pair(plot_singlets, "QC2_singlet_counts", width = 9, height = 5)

plot_UMAP_sample <- ggplot2::ggplot(UMAP_metadata, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = sample_accession)) + 
    ggplot2::geom_point(size = 0.12, alpha = 0.65000000000000002, stroke = 0) + ggplot2::labs(title = "GSE183837 Harmony UMAP by sample", 
    x = "Harmony UMAP 1", y = "Harmony UMAP 2", color = "Sample") + theme_QC

save_plot_pair(plot_UMAP_sample, "QC3_UMAP_sample", width = 10, height = 7)

plot_UMAP_group <- ggplot2::ggplot(UMAP_metadata, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = analysis_group)) + 
    ggplot2::geom_point(size = 0.12, alpha = 0.65000000000000002, stroke = 0) + ggplot2::labs(title = "GSE183837 Harmony UMAP by clinical group", 
    subtitle = "Visualization only; no cell-level group inference", x = "Harmony UMAP 1", y = "Harmony UMAP 2", 
    color = "Group") + theme_QC

save_plot_pair(plot_UMAP_group, "QC4_UMAP_group", width = 9, height = 7)

cluster_centroids <- aggregate(cbind(UMAPH_1, UMAPH_2) ~ seurat_clusters, data = UMAP_metadata, FUN = median)

plot_UMAP_cluster <- ggplot2::ggplot(UMAP_metadata, ggplot2::aes(x = UMAPH_1, y = UMAPH_2, color = seurat_clusters)) + 
    ggplot2::geom_point(size = 0.12, alpha = 0.65000000000000002, stroke = 0) + ggplot2::geom_text(data = cluster_centroids, 
    ggplot2::aes(x = UMAPH_1, y = UMAPH_2, label = seurat_clusters), inherit.aes = FALSE, size = 3, fontface = "bold") + 
    ggplot2::labs(title = "GSE183837 Harmony UMAP by primary clusters", subtitle = "Primary resolution = 0.4; biological identities are not assigned", 
        x = "Harmony UMAP 1", y = "Harmony UMAP 2", color = "Cluster") + theme_QC

save_plot_pair(plot_UMAP_cluster, "QC5_UMAP_cluster", width = 10, height = 7)

free_memory_after_GB <- get_free_memory_gb()

free_disk_after_GB <- get_free_disk_gb()

group_cell_counts_final <- table(factor(as.character(integrated_object$analysis_group), levels = c("CONTROL", "RIF")))

donor_group_table <- unique(data.frame(donor_id = as.character(integrated_object$donor_id), analysis_group = as.character(integrated_object$analysis_group), 
    stringsAsFactors = FALSE))

donor_group_table <- donor_group_table[order(donor_group_table$donor_id), , drop = FALSE]

if (nrow(donor_group_table) != 9L || anyDuplicated(donor_group_table$donor_id) > 0L || anyNA(donor_group_table$donor_id) || 
    anyNA(donor_group_table$analysis_group)) {
    stop("Final donor/group identity table is not one row per nine unique donors.")
}

group_donor_counts_final <- table(factor(donor_group_table$analysis_group, levels = c("CONTROL", "RIF")))

group_count_audit <- data.frame(metric = c("unique_donor_n", "Control_donor_n", "RIF_donor_n", "Control_cell_n", 
    "RIF_cell_n", "donor_group_design_is_3_Control_6_RIF"), value = c(nrow(donor_group_table), as.integer(group_donor_counts_final["CONTROL"]), 
    as.integer(group_donor_counts_final["RIF"]), as.integer(group_cell_counts_final["CONTROL"]), as.integer(group_cell_counts_final["RIF"]), 
    identical(as.integer(group_donor_counts_final), c(3L, 6L))), stringsAsFactors = FALSE)

write_csv_utf8(donor_group_table, file.path(table_root, "00_donor_group_identity.csv"))

write_csv_utf8(group_count_audit, file.path(table_root, "00_group_count_audit.csv"))

review_flags <- c(if (any(sample_status_all$high_doublet_review_flag)) {
    "SAMPLE_DOUBLET_FRACTION_ABOVE_20_PERCENT"
} else {
    character()
}, if (any(cluster_summary$cell_n < 100L)) {
    "PRIMARY_CLUSTER_BELOW_100_CELLS"
} else {
    character()
}, if (length(scale_layers_final) > 0L) {
    "SCALE_DATA_LAYER_RETAINED"
} else {
    character()
})

review_flags <- unique(review_flags[nzchar(review_flags)])

gate_checks <- data.frame(check = c("upstream_preflight_authorizes_integration", "nine_donor_doublet_jobs_complete", 
    "preQC_cells_partitioned_exactly", "all_donors_retain_at_least_1000_singlets", "final_cells_equal_singlet_sum", 
    "final_sample_n_is_9", "final_donor_n_is_9", "final_group_counts_are_3_Control_6_RIF", "PCA_has_50_dimensions", 
    "Harmony_has_at_least_30_dimensions", "UMAP_has_2_dimensions", "cluster_resolutions_0.2_0.4_0.6_present", "primary_resolution_0.4_present", 
    "nine_split_count_layers_retained", "one_joined_data_layer_present", "QC_UMAP_metadata_columns_unique", "final_checkpoint_exists", 
    "final_checkpoint_md5_available", "free_disk_after_run_at_least_10GB_or_unknown"), pass = c(preflight_gate$integration_authorized[[1]], 
    nrow(sample_status_all) == 9L && all(sample_status_all$status == "COMPLETE"), sum(sample_status_all$input_cell_n) == 
        50997L && sum(sample_status_all$singlet_cell_n + sample_status_all$doublet_cell_n) == 50997L, all(sample_status_all$singlet_cell_n >= 
        1000L), ncol(integrated_object) == total_singlets, length(unique(integrated_object$sample_accession)) == 
        9L, length(unique(integrated_object$donor_id)) == 9L, identical(as.integer(group_donor_counts_final), c(3L, 
        6L)), ncol(PCA_embedding) == 50L, ncol(Harmony_embedding) >= 30L, ncol(UMAP_embedding) == 2L, all(c("harmony_snn_res.0.2", 
        "harmony_snn_res.0.4", "harmony_snn_res.0.6") %in% colnames(metadata)), primary_cluster_column %in% colnames(metadata), 
    length(count_layers) == 9L, length(data_layers_final) == 1L, anyDuplicated(colnames(UMAP_metadata)) == 0L, 
    file.exists(final_checkpoint_path), nzchar(unname(tools::md5sum(final_checkpoint_path))), is.na(free_disk_after_GB) || 
        free_disk_after_GB >= 10), stringsAsFactors = FALSE)

write_csv_utf8(gate_checks, file.path(run_qc_root, "01_gate_checks.csv"))

all_gate_checks_pass <- all(gate_checks$pass)

gate_status <- if (all_gate_checks_pass && !length(review_flags)) {
    "PASS_GSE183837_INTEGRATION_EXECUTION_COMPLETE_QC_REQUIRED"
} else if (all_gate_checks_pass) {
    "PASS_GSE183837_INTEGRATION_EXECUTION_WITH_REVIEW_FLAGS"
} else {
    "HOLD_GSE183837_INTEGRATION_EXECUTION_FAILED"
}

integration_gate <- data.frame(dataset_id = dataset_id, integration_version = integration_version, gate_status = gate_status, 
    preQC_cell_n = 50997L, predicted_doublet_n = total_doublets, retained_singlet_n = total_singlets, predicted_doublet_fraction = total_doublets/50997L, 
    sample_n = 9L, Control_donor_n = as.integer(group_donor_counts_final["CONTROL"]), RIF_donor_n = as.integer(group_donor_counts_final["RIF"]), 
    Control_cell_n = as.integer(group_cell_counts_final["CONTROL"]), RIF_cell_n = as.integer(group_cell_counts_final["RIF"]), 
    feature_n = nrow(integrated_object), primary_cluster_n = length(unique(integrated_object$seurat_clusters)), 
    count_layer_n = length(count_layers), data_layer_n = length(data_layers_final), review_flags = collapse_unique(review_flags), 
    detailed_integration_QC_authorized = all_gate_checks_pass, marker_analysis_authorized = FALSE, biological_annotation_locked = FALSE, 
    disease_pseudobulk_testing_authorized = FALSE, result_lock_authorized = FALSE, stringsAsFactors = FALSE)

gate_path <- file.path(run_qc_root, "GSE183837_INTEGRATION_GATE.csv")

write_csv_utf8(integration_gate, gate_path)

card_lines <- c("# STEP04E1J GSE183837 Integration Execution Card", "", paste0("- Dataset: ", dataset_id), paste0("- Integration version: ", 
    integration_version), paste0("- Gate: ", gate_status), paste0("- Pre-QC cells: ", format(50997L, big.mark = ",")), 
    paste0("- Predicted doublets: ", format(total_doublets, big.mark = ","), " (", round(100 * total_doublets/50997L, 
        2), "%)"), paste0("- Retained singlets: ", format(total_singlets, big.mark = ",")), paste0("- Primary resolution-0.4 clusters: ", 
        length(unique(integrated_object$seurat_clusters))), paste0("- Raw count layers retained: ", length(count_layers)), 
    paste0("- Joined normalized data layers: ", length(data_layers_final)), paste0("- Review flags: ", if (length(review_flags)) {
        collapse_unique(review_flags)
    } else {
        "NONE"
    }), paste0("- Detailed integration QC authorized: ", all_gate_checks_pass), "- Marker analysis authorized: FALSE", 
    "- Biological annotation locked: FALSE", "- Disease pseudobulk testing authorized: FALSE", "- Result lock authorized: FALSE")

write_text(card_lines, file.path(run_qc_root, "GSE183837_INTEGRATION_CARD.md"))

interpretation_lines <- c("# GSE183837 integration execution interpretation", "", paste0("- Run ID: `", run_id, 
    "`"), paste0("- Gate: **", gate_status, "**"), "", "## Completed operations", "", "- scDblFinder was run independently within each donor/capture.", 
    "- Predicted doublets were excluded before integration.", "- Nine donor-level singlet checkpoints were retained for reproducibility and resume support.", 
    "- PCA and Harmony were performed on the merged singlet data.", "- UMAP and clustering were generated from Harmony dimensions 1-30.", 
    "- Resolution 0.4 is the primary technical cluster solution.", "- Nine split raw-count layers were retained and normalized data layers were joined.", 
    "", "## Not authorized yet", "", "- Marker extraction is deferred until detailed integration QC passes.", "- Broad cell-type identities are not assigned or locked.", 
    "- No cell-level Control-versus-RIF inference is permitted.", "- Donor-level disease testing is deferred until annotation and pseudobulk QC pass.", 
    "", "## Next step", "", "- Quantify sample mixing, cluster dominance, rare-cluster composition and group/sample confounding.", 
    "- Review UMAP and cluster-resolution stability.", "- Authorize memory-safe marker extraction only after the detailed QC Gate.")

write_text(interpretation_lines, file.path(run_qc_root, "GSE183837_INTEGRATION_INTERPRETATION.md"))

run_end <- Sys.time()

package_table <- data.frame(package = required_packages, version = vapply(required_packages, safe_package_version, 
    character(1)), stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(run_qc_root, "PACKAGE_VERSIONS.csv"))

resource_summary <- data.frame(metric = c("free_memory_before_GB", "free_memory_after_GB", "free_C_drive_before_GB", 
    "free_C_drive_after_GB", "merged_checkpoint_size_GB", "prejoin_checkpoint_size_GB", "final_checkpoint_size_GB"), 
    value = c(free_memory_before_GB, free_memory_after_GB, free_disk_before_GB, free_disk_after_GB, if (file.exists(merged_counts_checkpoint)) {
        round(file.info(merged_counts_checkpoint)$size/1024^3, 3)
    } else {
        NA_real_
    }, if (file.exists(prejoin_checkpoint_path)) {
        round(file.info(prejoin_checkpoint_path)$size/1024^3, 3)
    } else {
        NA_real_
    }, round(file.info(final_checkpoint_path)$size/1024^3, 3)), stringsAsFactors = FALSE)

write_csv_utf8(resource_summary, file.path(run_qc_root, "02_resource_summary.csv"))

write_text(capture.output(sessionInfo()), session_path)

checkpoint_manifest <- data.frame(artifact = c("merged_singlet_checkpoint", "Harmony_prejoin_checkpoint", "final_integrated_checkpoint"), 
    path = normalize_slashes(c(merged_counts_checkpoint, prejoin_checkpoint_path, final_checkpoint_path)), size_bytes = file.info(c(merged_counts_checkpoint, 
        prejoin_checkpoint_path, final_checkpoint_path))$size, md5 = unname(tools::md5sum(c(merged_counts_checkpoint, 
        prejoin_checkpoint_path, final_checkpoint_path))), stringsAsFactors = FALSE)

write_csv_utf8(checkpoint_manifest, file.path(manifest_root, "checkpoint_manifest.csv"))

artifact_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(manifest_root, 
    recursive = TRUE, full.names = TRUE), list.files(run_qc_root, recursive = TRUE, full.names = TRUE), log_path, 
    run_manifest_path, parameter_path, session_path))

artifact_files <- artifact_files[file.exists(artifact_files) & !dir.exists(artifact_files)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_files), size_bytes = file.info(artifact_files)$size, 
    md5 = unname(tools::md5sum(artifact_files)), stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(run_qc_root, "ARTIFACT_MANIFEST.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, supersedes_failed_executions = "V1_0_HARMONY_ASSAY_USE_FAILURE;V1_1_JOINLAYERS_LENGTH_MISMATCH;V1_2_GGPLOT_DUPLICATE_METADATA_COLUMNS;V1_3_CELL_COUNTS_MISTAKEN_FOR_DONOR_COUNTS", 
    integration_version = integration_version, dataset_id = dataset_id, started_at = as.character(run_start), completed_at = as.character(run_end), 
    elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 2), gate_status = gate_status, 
    preQC_cell_n = 50997L, predicted_doublet_n = total_doublets, retained_singlet_n = total_singlets, sample_n = 9L, 
    primary_cluster_n = length(unique(integrated_object$seurat_clusters)), count_layer_n = length(count_layers), 
    data_layer_n = length(data_layers_final), final_checkpoint_size_GB = round(file.info(final_checkpoint_path)$size/1024^3, 
        3), final_checkpoint_md5 = unname(tools::md5sum(final_checkpoint_path)), detailed_integration_QC_authorized = all_gate_checks_pass, 
    marker_analysis_authorized = FALSE, result_lock_authorized = FALSE, output_root = normalize_slashes(output_root), 
    qc_root = normalize_slashes(run_qc_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

log_event("Completed", run_id, "gate_status=", gate_status, "preQC_cells=", 50997L, "predicted_doublets=", total_doublets, 
    "retained_singlets=", total_singlets, "primary_clusters=", length(unique(integrated_object$seurat_clusters)))

small_files <- unique(c(list.files(table_root, recursive = TRUE, full.names = TRUE), list.files(manifest_root, 
    recursive = TRUE, full.names = TRUE), list.files(run_qc_root, recursive = TRUE, full.names = TRUE), log_path, 
    run_manifest_path, parameter_path, session_path, normalizePath(.vdr_script_path(script_name), winslash = "/", 
        mustWork = FALSE)))

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

cat("STEP04E1J GSE183837 doublet/Harmony integration completed\n")

cat("============================================================\n")

cat("Gate status: ", gate_status, "\n", sep = "")

cat("Pre-QC cells: ", format(50997L, big.mark = ","), "\n", sep = "")

cat("Predicted doublets: ", format(total_doublets, big.mark = ","), " (", round(100 * total_doublets/50997L, 2), 
    "%)\n", sep = "")

cat("Retained singlets: ", format(total_singlets, big.mark = ","), "\n", sep = "")

cat("Primary resolution-0.4 clusters: ", length(unique(integrated_object$seurat_clusters)), "\n", sep = "")

cat("Count layers retained: ", length(count_layers), "\n", sep = "")

cat("Joined data layers: ", length(data_layers_final), "\n", sep = "")

cat("Final checkpoint: ", final_checkpoint_path, "\n", sep = "")

cat("QC root: ", run_qc_root, "\n", sep = "")

if (zip_ok) {
    cat("Upload this ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe integration outputs and QC report are complete.\n", sep = "")
}

