#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C01: Doublet assessment and memory-safe single-cell normalization

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(3600L, getOption("timeout")))

options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260804)

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04E1_", format(Sys.time(), "%Y%m%d_%H%M%S"))

integration_version <- "SINGLE_CELL_INTEGRATION_v1.1_MEMORY_SAFE"

annotation_version <- "PRELIMINARY_BROAD_CELL_ANNOTATION_v1.0"

localization_version <- "SINGLE_CELL_VDR_PROGRAM_LOCALIZATION_v1.0"

required_packages <- c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "digest", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". Install them before rerunning. This script does not install or update packages.")
}

doublet_packages <- c("scDblFinder", "SingleCellExperiment", "SummarizedExperiment")

doublet_mode <- if (all(vapply(doublet_packages, requireNamespace, logical(1), quietly = TRUE))) {
    "SCDBLFINDER"
} else {
    "NOT_RUN_REQUIRED_PACKAGES_MISSING"
}

if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
}

dirs <- list(control = file.path(step04_root, "00_module_control"), code = file.path(step04_root, "02_code"), parameters = file.path(step04_root, 
    "03_parameters"), environment = file.path(step04_root, "04_environment"), execution = file.path(step04_root, 
    "05_execution"), results = file.path(step04_root, "06_results"), figures = file.path(step04_root, "07_figures_tables"), 
    qc = file.path(step04_root, "08_qc"), text = file.path(step04_root, "09_manuscript_text"), archive_parent = file.path(step04_root, 
        "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory_path)) {
        stop("Required directory could not be created: ", directory_path)
    }
}

output_root <- file.path(dirs$results, "S4E1_sc", "SCI_v1_0")

object_root <- file.path(output_root, "objects")

table_root <- file.path(output_root, "tables")

pseudobulk_root <- file.path(output_root, "pseudobulk")

checkpoint_root <- file.path(output_root, "checkpoints")

archive_root <- file.path(dirs$archive_parent, paste0("S4E1_", format(Sys.time(), "%Y%m%d_%H%M%S")))

for (directory_path in c(output_root, object_root, table_root, pseudobulk_root, checkpoint_root, file.path(archive_root, 
    "res"), file.path(archive_root, "qc"), file.path(archive_root, "txt"), file.path(archive_root, "code"), file.path(archive_root, 
    "env"), file.path(archive_root, "log"), file.path(archive_root, "ctl"))) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory_path)) {
        stop("Required output/archive directory could not be created: ", directory_path)
    }
    if (nchar(directory_path) >= 220L) {
        stop("Directory exceeds the 220-character soft limit: ", directory_path)
    }
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    message_text <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(message_text, "\n")
    cat(message_text, "\n", file = log_path, append = TRUE)
}

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
    if (nchar(path) >= 220L) {
        stop("Output path exceeds the 220-character soft limit: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Output directory could not be created: ", dirname(path))
    }
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

write_text <- function(lines, path) {
    if (nchar(path) >= 220L) {
        stop("Output path exceeds the 220-character soft limit: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
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

same_normalized_path <- function(path_a, path_b) {
    if (is.null(path_a) || is.null(path_b) || length(path_a) != 1L || length(path_b) != 1L || is.na(path_a) || 
        is.na(path_b)) {
        return(FALSE)
    }
    normalized_a <- tryCatch(normalizePath(path_a, winslash = "/", mustWork = FALSE), error = function(error_object) {
        normalize_slashes(path_a)
    })
    normalized_b <- tryCatch(normalizePath(path_b, winslash = "/", mustWork = FALSE), error = function(error_object) {
        normalize_slashes(path_b)
    })
    identical(tolower(normalized_a), tolower(normalized_b))
}

copy_if_different <- function(source_path, destination_path) {
    if (same_normalized_path(source_path, destination_path)) {
        return(TRUE)
    }
    dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(source_path, destination_path, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
}

copy_validate <- function(source_path, destination_path) {
    copied <- copy_if_different(source_path, destination_path)
    pass <- isTRUE(copied) && file.exists(destination_path) && identical(as.numeric(file.info(source_path)$size), 
        as.numeric(file.info(destination_path)$size)) && identical(md5_file(source_path), md5_file(destination_path)) && 
        identical(sha256_file(source_path), sha256_file(destination_path))
    if (!pass) {
        stop("Copy/hash validation failed: ", source_path)
    }
    TRUE
}

as_logical_flag <- function(x) {
    if (is.logical(x)) {
        return(x)
    }
    toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

get_layer_matrix <- function(object, layer_name) {
    matrix_value <- tryCatch(SeuratObject::LayerData(object, assay = "RNA", layer = layer_name), error = function(error_object) {
        NULL
    })
    if (is.null(matrix_value)) {
        matrix_value <- tryCatch(Seurat::GetAssayData(object, assay = "RNA", slot = layer_name), error = function(error_object) {
            NULL
        })
    }
    if (is.null(matrix_value)) {
        stop("Could not retrieve RNA ", layer_name, " matrix.")
    }
    matrix_value
}

collapse_duplicate_symbols <- function(counts, gene_symbol) {
    gene_symbol <- trimws(as.character(gene_symbol))
    valid <- !is.na(gene_symbol) & nzchar(gene_symbol) & gene_symbol != "."
    counts <- counts[valid, , drop = FALSE]
    gene_symbol <- gene_symbol[valid]
    unique_symbols <- unique(gene_symbol)
    aggregation_matrix <- Matrix::sparseMatrix(i = seq_along(gene_symbol), j = match(gene_symbol, unique_symbols), 
        x = 1, dims = c(length(gene_symbol), length(unique_symbols)))
    collapsed_counts <- Matrix::t(aggregation_matrix) %*% counts
    collapsed_counts <- methods::as(collapsed_counts, "dgCMatrix")
    rownames(collapsed_counts) <- unique_symbols
    colnames(collapsed_counts) <- colnames(counts)
    collapsed_counts
}

run_doublet_assessment <- function(counts, sample_accession) {
    if (doublet_mode != "SCDBLFINDER") {
        return(list(keep = rep(TRUE, ncol(counts)), class = rep("NOT_RUN", ncol(counts)), score = rep(NA_real_, 
            ncol(counts)), method = "NOT_RUN_REQUIRED_PACKAGES_MISSING"))
    }
    set.seed(20260804 + sum(utf8ToInt(sample_accession)))
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
    sce <- scDblFinder::scDblFinder(sce, verbose = FALSE)
    doublet_metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    class_value <- as.character(doublet_metadata$scDblFinder.class)
    score_value <- as.numeric(doublet_metadata$scDblFinder.score)
    list(keep = class_value == "singlet", class = class_value, score = score_value, method = "SCDBLFINDER")
}

safe_zscore <- function(x) {
    standard_deviation <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(standard_deviation) || standard_deviation == 0) {
        return(rep(0, length(x)))
    }
    as.numeric(scale(x))
}

calculate_module_score <- function(data_matrix, genes) {
    present_genes <- intersect(genes, rownames(data_matrix))
    if (length(present_genes) == 0L) {
        return(list(score = rep(NA_real_, ncol(data_matrix)), present_genes = character(0)))
    }
    list(score = as.numeric(Matrix::colMeans(data_matrix[present_genes, , drop = FALSE])), present_genes = present_genes)
}

calculate_signed_program_score <- function(data_matrix, gene_table) {
    positive_genes <- gene_table$feature_id[gene_table$interaction_logFC > 0]
    negative_genes <- gene_table$feature_id[gene_table$interaction_logFC < 0]
    positive_present <- intersect(positive_genes, rownames(data_matrix))
    negative_present <- intersect(negative_genes, rownames(data_matrix))
    positive_score <- if (length(positive_present) > 0L) {
        as.numeric(Matrix::colMeans(data_matrix[positive_present, , drop = FALSE]))
    }
    else {
        rep(0, ncol(data_matrix))
    }
    negative_score <- if (length(negative_present) > 0L) {
        as.numeric(Matrix::colMeans(data_matrix[negative_present, , drop = FALSE]))
    }
    else {
        rep(0, ncol(data_matrix))
    }
    list(score = positive_score - negative_score, positive_present_n = length(positive_present), negative_present_n = length(negative_present), 
        total_present_n = length(union(positive_present, negative_present)))
}

build_sparse_pseudobulk <- function(counts, metadata, minimum_cells = 50L) {
    group_table <- unique(metadata[, c("sample_accession", "preliminary_broad_cell_type"), drop = FALSE])
    group_table$group_id <- paste(group_table$sample_accession, group_table$preliminary_broad_cell_type, sep = "__")
    cell_group_id <- paste(metadata$sample_accession, metadata$preliminary_broad_cell_type, sep = "__")
    group_levels <- unique(cell_group_id)
    group_factor <- factor(cell_group_id, levels = group_levels)
    membership <- Matrix::sparseMatrix(i = seq_along(group_factor), j = as.integer(group_factor), x = 1, dims = c(length(group_factor), 
        length(group_levels)), dimnames = list(colnames(counts), group_levels))
    cell_n <- as.numeric(Matrix::colSums(membership))
    keep_group <- cell_n >= minimum_cells
    pseudobulk_counts <- counts %*% membership[, keep_group, drop = FALSE]
    pseudobulk_counts <- methods::as(pseudobulk_counts, "dgCMatrix")
    group_metadata <- data.frame(group_id = group_levels, cell_n = cell_n, stringsAsFactors = FALSE)
    group_metadata <- merge(group_metadata, group_table, by = "group_id", all.x = TRUE, sort = FALSE)
    group_metadata$included_in_pseudobulk <- group_metadata$cell_n >= minimum_cells
    list(counts = pseudobulk_counts, metadata = group_metadata)
}

log_event("Starting", run_id, "Doublet mode:", doublet_mode)

qc1_root <- file.path(dirs$results, "S4D_qc_lock", "BVP_v1_3")

qc1_gate_path <- file.path(qc1_root, "STEP04D_QC1_GATE.csv")

candidate_path <- file.path(qc1_root, "cross_context_candidates_v1_3.csv")

watchlist_path <- file.path(qc1_root, "sc_sample_watchlist.csv")

preQC_root <- file.path(dirs$results, "S4D_sc_preQC", "SCP_v1_0")

design_root <- file.path(dirs$results, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0")

required_inputs <- c(qc1_gate_path, candidate_path, watchlist_path, file.path(design_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv"), 
    file.path(design_root, "COMBINED_SAMPLE_DESIGN_v1.0_GSE183837.csv"))

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
    stop("Required STEP04E1 input is missing:\n", paste(missing_inputs, collapse = "\n"))
}

qc1_gate <- .vdr_read_csv(qc1_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

qc1_status <- qc1_gate$status[qc1_gate$gate_component == "STEP04D-QC1 gate"]

if (length(qc1_status) != 1L || qc1_status != "READY_FOR_SINGLE_CELL_INTEGRATION_WITH_GSE183837_SAMPLE_WATCHLIST") {
    stop("STEP04D-QC1 did not authorize single-cell integration.")
}

candidate_table <- .vdr_read_csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)

candidate_table$occupancy_supported <- as_logical_flag(candidate_table$occupancy_supported)

candidate_table$shared_occupancy_supported <- as_logical_flag(candidate_table$shared_occupancy_supported)

candidate_table$any_promoter_overlap_support <- as_logical_flag(candidate_table$any_promoter_overlap_support)

if (nrow(candidate_table) != 436L || length(unique(candidate_table$feature_id)) != 436L) {
    stop("The corrected 436-candidate table failed validation.")
}

watchlist <- .vdr_read_csv(watchlist_path, stringsAsFactors = FALSE, check.names = FALSE)

read_locked_design <- function(dataset_id) {
    .vdr_read_csv(file.path(design_root, paste0("COMBINED_SAMPLE_DESIGN_v1.0_", dataset_id, ".csv")), stringsAsFactors = FALSE, 
        check.names = FALSE)
}

designs <- list(GSE250130 = read_locked_design("GSE250130"), GSE183837 = read_locked_design("GSE183837"))

expected_sample_n <- c(GSE250130 = 28L, GSE183837 = 9L)

for (dataset_id in names(designs)) {
    if (nrow(designs[[dataset_id]]) != expected_sample_n[[dataset_id]] || anyDuplicated(designs[[dataset_id]]$sample_accession)) {
        stop(dataset_id, " locked single-cell design failed sample-count or uniqueness validation.")
    }
}

broad_marker_sets <- list(Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN4", "MUC1"), Stromal_fibroblast = c("COL1A1", 
    "COL1A2", "DCN", "COL3A1", "LUM", "COL6A1", "C7"), Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "ESAM", 
    "RAMP2"), Immune = c("PTPRC", "TYROBP", "LST1", "FCER1G", "HLA-DRA", "CTSS"), Pericyte_SMC = c("RGS5", "CSPG4", 
    "MCAM", "ACTA2", "TAGLN", "MYL9", "COL4A1"), Erythroid = c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA"))

cycling_markers <- c("MKI67", "TOP2A", "UBE2C", "TYMS", "PCNA", "CENPF")

program_definitions <- list(RNA_ALL_436 = candidate_table, OCCUPANCY_SUPPORTED_284 = candidate_table[candidate_table$occupancy_supported, 
    , drop = FALSE], SHARED_OCCUPANCY_46 = candidate_table[candidate_table$shared_occupancy_supported, , drop = FALSE], 
    PROMOTER_OCCUPANCY_84 = candidate_table[candidate_table$any_promoter_overlap_support, , drop = FALSE])

process_dataset <- function(dataset_id, sample_design) {
    log_event("Processing dataset", dataset_id)
    dataset_object_root <- file.path(object_root, dataset_id)
    dataset_table_root <- file.path(table_root, dataset_id)
    dataset_pseudobulk_root <- file.path(pseudobulk_root, dataset_id)
    dataset_checkpoint_root <- file.path(checkpoint_root, dataset_id)
    for (directory_path in c(dataset_object_root, dataset_table_root, dataset_pseudobulk_root, dataset_checkpoint_root)) {
        dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
        if (!dir.exists(directory_path)) {
            stop("Could not create dataset output directory: ", directory_path)
        }
    }
    integrated_checkpoint_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_integrated_checkpoint.rds"))
    checkpoint_status_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_checkpoint_status.csv"))
    merged_counts_checkpoint_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_post_doublet_merged_counts_checkpoint.rds"))
    merged_counts_status_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_post_doublet_merged_counts_status.csv"))
    checkpoint_reusable <- FALSE
    if (file.exists(integrated_checkpoint_path) && file.exists(checkpoint_status_path)) {
        checkpoint_status <- .vdr_read_csv(checkpoint_status_path, stringsAsFactors = FALSE, check.names = FALSE)
        checkpoint_reusable <- nrow(checkpoint_status) == 1L && checkpoint_status$dataset_id[[1]] == dataset_id && 
            checkpoint_status$doublet_mode[[1]] == doublet_mode && checkpoint_status$sample_n[[1]] == nrow(sample_design) && 
            checkpoint_status$checkpoint_status[[1]] == "COMPLETE"
    }
    if (checkpoint_reusable) {
        log_event("Reusing integrated checkpoint for", dataset_id)
        object <- readRDS(integrated_checkpoint_path)
        sample_doublet_qc <- .vdr_read_csv(file.path(dataset_table_root, paste0(dataset_id, "_sample_doublet_qc.csv")), 
            stringsAsFactors = FALSE, check.names = FALSE)
    }
    else {
        merged_counts_reusable <- FALSE
        if (file.exists(merged_counts_checkpoint_path) && file.exists(merged_counts_status_path)) {
            merged_counts_status <- .vdr_read_csv(merged_counts_status_path, stringsAsFactors = FALSE, check.names = FALSE)
            merged_counts_reusable <- nrow(merged_counts_status) == 1L && merged_counts_status$dataset_id[[1]] == 
                dataset_id && merged_counts_status$doublet_mode[[1]] == doublet_mode && merged_counts_status$sample_n[[1]] == 
                nrow(sample_design) && merged_counts_status$checkpoint_status[[1]] == "COMPLETE"
        }
        if (merged_counts_reusable) {
            log_event("Reusing post-doublet merged-count checkpoint for", dataset_id)
            object <- readRDS(merged_counts_checkpoint_path)
            sample_doublet_qc <- .vdr_read_csv(file.path(dataset_table_root, paste0(dataset_id, "_sample_doublet_qc.csv")), 
                stringsAsFactors = FALSE, check.names = FALSE)
        }
        else {
            sample_objects <- list()
            doublet_qc_rows <- list()
            for (sample_index in seq_len(nrow(sample_design))) {
                sample_accession <- sample_design$sample_accession[sample_index]
                log_event(dataset_id, "sample", sample_index, "/", nrow(sample_design), sample_accession)
                preQC_path <- file.path(preQC_root, dataset_id, paste0(sample_accession, "_preQC.rds"))
                if (!file.exists(preQC_path)) {
                  stop("Pre-QC object missing: ", preQC_path)
                }
                preQC_object <- readRDS(preQC_path)
                if (!all(c("counts", "feature_table", "cell_metadata") %in% names(preQC_object))) {
                  stop(dataset_id, " ", sample_accession, " pre-QC object lacks required components.")
                }
                counts <- preQC_object$counts
                if (nrow(preQC_object$feature_table) != nrow(counts)) {
                  stop(dataset_id, " ", sample_accession, " feature table/count matrix row mismatch.")
                }
                counts <- collapse_duplicate_symbols(counts = counts, gene_symbol = preQC_object$feature_table$gene_symbol)
                doublet_result <- run_doublet_assessment(counts = counts, sample_accession = sample_accession)
                if (length(doublet_result$keep) != ncol(counts)) {
                  stop(dataset_id, " ", sample_accession, " doublet-result length mismatch.")
                }
                retained_counts <- counts[, doublet_result$keep, drop = FALSE]
                retained_metadata <- data.frame(row.names = colnames(retained_counts), sample_accession = sample_accession, 
                  dataset_id = dataset_id, doublet_method = doublet_result$method, doublet_class = doublet_result$class[doublet_result$keep], 
                  doublet_score = doublet_result$score[doublet_result$keep], stringsAsFactors = FALSE)
                design_row <- sample_design[sample_design$sample_accession == sample_accession, , drop = FALSE]
                for (column_name in setdiff(names(design_row), c("sample_accession", "dataset_id"))) {
                  retained_metadata[[column_name]] <- design_row[[column_name]][[1]]
                }
                watch_row <- watchlist[watchlist$dataset_id == dataset_id & watchlist$sample_accession == sample_accession, 
                  , drop = FALSE]
                retained_metadata$preQC_watch_tier <- if (nrow(watch_row) == 1L) {
                  watch_row$QC_review_tier[[1]]
                }
                else {
                  "PASS"
                }
                sample_object <- SeuratObject::CreateSeuratObject(counts = retained_counts, project = dataset_id, 
                  min.cells = 0, min.features = 0, meta.data = retained_metadata)
                sample_objects[[sample_accession]] <- sample_object
                doublet_qc_rows[[sample_accession]] <- data.frame(dataset_id = dataset_id, sample_accession = sample_accession, 
                  pre_doublet_cell_n = ncol(counts), post_doublet_cell_n = ncol(retained_counts), predicted_doublet_n = sum(!doublet_result$keep), 
                  predicted_doublet_fraction = mean(!doublet_result$keep), doublet_method = doublet_result$method, 
                  preQC_watch_tier = retained_metadata$preQC_watch_tier[[1]], stringsAsFactors = FALSE)
                rm(preQC_object, counts, retained_counts, sample_object)
                gc()
            }
            sample_doublet_qc <- do.call(rbind, doublet_qc_rows)
            write_csv_utf8(sample_doublet_qc, file.path(dataset_table_root, paste0(dataset_id, "_sample_doublet_qc.csv")))
            log_event("Merging", dataset_id, "sample objects")
            first_sample <- names(sample_objects)[[1]]
            object <- merge(x = sample_objects[[first_sample]], y = sample_objects[setdiff(names(sample_objects), 
                first_sample)], add.cell.ids = NULL, merge.data = FALSE, project = dataset_id)
            rm(sample_objects)
            gc()
            saveRDS(object, merged_counts_checkpoint_path, compress = FALSE)
            write_csv_utf8(data.frame(dataset_id = dataset_id, doublet_mode = doublet_mode, sample_n = nrow(sample_design), 
                cell_n = ncol(object), checkpoint_status = "COMPLETE", checkpoint_path = normalize_slashes(merged_counts_checkpoint_path), 
                stringsAsFactors = FALSE), merged_counts_status_path)
            log_event("Saved post-doublet merged-count checkpoint for", dataset_id)
        }
        log_event("Keeping sample-specific RNA layers through normalization for", dataset_id)
        object <- Seurat::NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
            block.size = 2000, verbose = TRUE)
        object <- Seurat::FindVariableFeatures(object, assay = "RNA", selection.method = "vst", nfeatures = 3000, 
            verbose = TRUE)
        object <- Seurat::ScaleData(object, assay = "RNA", features = SeuratObject::VariableFeatures(object), block.size = 1000, 
            verbose = TRUE)
        object <- Seurat::RunPCA(object, assay = "RNA", features = SeuratObject::VariableFeatures(object), npcs = 50, 
            seed.use = 20260804, verbose = TRUE)
        object <- harmony::RunHarmony(object = object, group.by.vars = "sample_accession", reduction = "pca", dims.use = 1:40, 
            reduction.save = "harmony", plot_convergence = FALSE, verbose = TRUE)
        object <- Seurat::FindNeighbors(object, reduction = "harmony", dims = 1:30, k.param = 30, verbose = TRUE)
        object <- Seurat::FindClusters(object, resolution = 0.20000000000000001, algorithm = 1, random.seed = 20260804, 
            verbose = TRUE)
        object <- Seurat::FindClusters(object, resolution = 0.59999999999999998, algorithm = 1, random.seed = 20260804, 
            verbose = TRUE)
        object <- Seurat::FindClusters(object, resolution = 0.40000000000000002, algorithm = 1, random.seed = 20260804, 
            verbose = TRUE)
        object <- Seurat::RunUMAP(object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", 
            reduction.key = "UMAPH_", seed.use = 20260804, verbose = TRUE)
        reductions_to_keep <- intersect(c("pca", "harmony", "umap_harmony"), SeuratObject::Reductions(object))
        graphs_to_keep <- SeuratObject::Graphs(object)
        RNA_layers_before_join <- SeuratObject::Layers(object[["RNA"]])
        RNA_layers_to_keep <- grep("^(counts|data)(\\.|$)", RNA_layers_before_join, value = TRUE)
        if (length(RNA_layers_to_keep) == 0L) {
            stop(dataset_id, " has no counts/data RNA layers before memory cleanup.")
        }
        object <- Seurat::DietSeurat(object, assays = "RNA", layers = list(RNA = RNA_layers_to_keep), dimreducs = reductions_to_keep, 
            graphs = graphs_to_keep, misc = TRUE)
        gc()
        if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
            log_event("Joining RNA layers after PCA/Harmony/UMAP for", dataset_id)
            object <- SeuratObject::JoinLayers(object, assay = "RNA")
            gc()
        }
        object <- Seurat::DietSeurat(object, assays = "RNA", layers = c("counts", "data"), dimreducs = reductions_to_keep, 
            graphs = graphs_to_keep, misc = TRUE)
        saveRDS(object, integrated_checkpoint_path, compress = FALSE)
        write_csv_utf8(data.frame(dataset_id = dataset_id, doublet_mode = doublet_mode, sample_n = nrow(sample_design), 
            cell_n = ncol(object), checkpoint_status = "COMPLETE", checkpoint_path = normalize_slashes(integrated_checkpoint_path), 
            stringsAsFactors = FALSE), checkpoint_status_path)
    }
    if (!all(sample_design$sample_accession %in% unique(object$sample_accession))) {
        stop(dataset_id, " integrated object does not retain all locked samples.")
    }
    data_matrix <- get_layer_matrix(object, "data")
    module_score_columns <- character(0)
    module_coverage_rows <- list()
    for (marker_name in names(broad_marker_sets)) {
        score_value <- calculate_module_score(data_matrix = data_matrix, genes = broad_marker_sets[[marker_name]])
        column_name <- paste0("lineage_score_", marker_name)
        object[[column_name]] <- safe_zscore(score_value$score)
        module_score_columns <- c(module_score_columns, column_name)
        module_coverage_rows[[marker_name]] <- data.frame(dataset_id = dataset_id, module_name = marker_name, expected_gene_n = length(broad_marker_sets[[marker_name]]), 
            present_gene_n = length(score_value$present_genes), present_genes = paste(score_value$present_genes, 
                collapse = ";"), stringsAsFactors = FALSE)
    }
    cycling_score <- calculate_module_score(data_matrix = data_matrix, genes = cycling_markers)
    object$cycling_score <- safe_zscore(cycling_score$score)
    cluster_column <- "seurat_clusters"
    metadata_dt <- data.table::as.data.table(object@meta.data, keep.rownames = "cell_barcode")
    cluster_score_dt <- metadata_dt[, lapply(.SD, mean, na.rm = TRUE), by = cluster_column, .SDcols = module_score_columns]
    score_matrix <- as.matrix(cluster_score_dt[, module_score_columns, with = FALSE])
    top_index <- apply(score_matrix, 1L, which.max)
    sorted_scores <- t(apply(score_matrix, 1L, sort, decreasing = TRUE))
    top_score <- sorted_scores[, 1L]
    second_score <- if (ncol(sorted_scores) >= 2L) {
        sorted_scores[, 2L]
    }
    else {
        rep(NA_real_, nrow(sorted_scores))
    }
    broad_labels <- names(broad_marker_sets)[top_index]
    label_margin <- top_score - second_score
    preliminary_label <- ifelse(!is.finite(top_score) | top_score < 0 | !is.finite(label_margin) | label_margin < 
        0.14999999999999999, "AMBIGUOUS", broad_labels)
    cluster_annotation <- data.frame(dataset_id = dataset_id, seurat_cluster = as.character(cluster_score_dt[[cluster_column]]), 
        preliminary_broad_cell_type = preliminary_label, top_lineage_score = top_score, second_lineage_score = second_score, 
        lineage_score_margin = label_margin, stringsAsFactors = FALSE)
    cluster_annotation <- cbind(cluster_annotation, as.data.frame(score_matrix, stringsAsFactors = FALSE))
    label_map <- setNames(cluster_annotation$preliminary_broad_cell_type, cluster_annotation$seurat_cluster)
    object$preliminary_broad_cell_type <- unname(label_map[as.character(object$seurat_clusters)])
    object$cycling_state <- ifelse(object$cycling_score >= 1, "CYCLING_HIGH", "CYCLING_LOW")
    write_csv_utf8(do.call(rbind, module_coverage_rows), file.path(dataset_table_root, paste0(dataset_id, "_broad_marker_coverage.csv")))
    write_csv_utf8(cluster_annotation, file.path(dataset_table_root, paste0(dataset_id, "_preliminary_cluster_annotation.csv")))
    SeuratObject::Idents(object) <- "seurat_clusters"
    marker_results <- tryCatch(Seurat::FindAllMarkers(object, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
        test.use = "wilcox", max.cells.per.ident = 1500, random.seed = 20260804, verbose = TRUE), error = function(error_object) {
        log_event(dataset_id, "FindAllMarkers failed:", conditionMessage(error_object))
        data.frame(marker_status = paste0("FAILED_", conditionMessage(error_object)), stringsAsFactors = FALSE)
    })
    write_csv_utf8(marker_results, file.path(dataset_table_root, paste0(dataset_id, "_cluster_markers.csv")))
    if ("cluster" %in% names(marker_results)) {
        fold_change_column <- intersect(c("avg_log2FC", "avg_logFC"), names(marker_results))
        if (length(fold_change_column) > 0L) {
            marker_dt <- data.table::as.data.table(marker_results)
            data.table::setorderv(marker_dt, cols = c("cluster", fold_change_column[[1]]), order = c(1L, -1L))
            top_markers <- marker_dt[, head(.SD, 25L), by = cluster]
            write_csv_utf8(as.data.frame(top_markers), file.path(dataset_table_root, paste0(dataset_id, "_top25_cluster_markers.csv")))
        }
    }
    program_coverage_rows <- list()
    for (program_name in names(program_definitions)) {
        score_value <- calculate_signed_program_score(data_matrix = data_matrix, gene_table = program_definitions[[program_name]])
        score_column <- paste0("VDR_program_", program_name)
        object[[score_column]] <- score_value$score
        program_coverage_rows[[program_name]] <- data.frame(dataset_id = dataset_id, program_name = program_name, 
            frozen_gene_n = nrow(program_definitions[[program_name]]), positive_present_n = score_value$positive_present_n, 
            negative_present_n = score_value$negative_present_n, total_present_n = score_value$total_present_n, 
            stringsAsFactors = FALSE)
    }
    write_csv_utf8(do.call(rbind, program_coverage_rows), file.path(dataset_table_root, paste0(dataset_id, "_VDR_program_coverage.csv")))
    metadata_dt <- data.table::as.data.table(object@meta.data, keep.rownames = "cell_barcode")
    cluster_donor_counts <- metadata_dt[, .(cell_n = .N), by = .(seurat_clusters, sample_accession)]
    cluster_totals <- cluster_donor_counts[, .(cluster_cell_n = sum(cell_n), donor_n = data.table::uniqueN(sample_accession), 
        max_donor_fraction = max(cell_n)/sum(cell_n)), by = seurat_clusters]
    write_csv_utf8(as.data.frame(cluster_donor_counts), file.path(dataset_table_root, paste0(dataset_id, "_cluster_by_donor_counts.csv")))
    write_csv_utf8(as.data.frame(cluster_totals), file.path(dataset_table_root, paste0(dataset_id, "_cluster_donor_mixing_audit.csv")))
    broad_composition <- metadata_dt[, .(cell_n = .N), by = .(sample_accession, preliminary_broad_cell_type)]
    broad_composition[, `:=`(sample_total, sum(cell_n)), by = sample_accession]
    broad_composition[, `:=`(cell_fraction, cell_n/sample_total)]
    write_csv_utf8(as.data.frame(broad_composition), file.path(dataset_table_root, paste0(dataset_id, "_broad_cell_composition.csv")))
    program_columns <- paste0("VDR_program_", names(program_definitions))
    donor_program_summary <- metadata_dt[, c(list(cell_n = .N), lapply(.SD, function(x) {
        mean(x, na.rm = TRUE)
    })), by = .(sample_accession, preliminary_broad_cell_type), .SDcols = program_columns]
    donor_program_summary <- merge(as.data.frame(donor_program_summary), sample_design, by = "sample_accession", 
        all.x = TRUE, sort = FALSE)
    write_csv_utf8(donor_program_summary, file.path(dataset_table_root, paste0(dataset_id, "_donor_celltype_VDR_program_scores.csv")))
    counts_matrix <- get_layer_matrix(object, "counts")
    pseudobulk <- build_sparse_pseudobulk(counts = counts_matrix, metadata = object@meta.data, minimum_cells = 50L)
    pseudobulk_matrix_path <- file.path(dataset_pseudobulk_root, paste0(dataset_id, "_broad_cell_pseudobulk_counts.rds"))
    saveRDS(pseudobulk$counts, pseudobulk_matrix_path, compress = TRUE)
    pseudobulk_metadata <- merge(pseudobulk$metadata, sample_design, by = "sample_accession", all.x = TRUE, sort = FALSE)
    write_csv_utf8(pseudobulk_metadata, file.path(dataset_pseudobulk_root, paste0(dataset_id, "_broad_cell_pseudobulk_metadata.csv")))
    object_path <- file.path(dataset_object_root, paste0(dataset_id, "_integrated_preliminary_annotation.rds"))
    reductions_to_keep <- intersect(c("pca", "harmony", "umap_harmony"), SeuratObject::Reductions(object))
    graphs_to_keep <- SeuratObject::Graphs(object)
    object <- Seurat::DietSeurat(object, assays = "RNA", layers = c("counts", "data"), dimreducs = reductions_to_keep, 
        graphs = graphs_to_keep, misc = TRUE)
    saveRDS(object, object_path, compress = FALSE)
    if (file.exists(integrated_checkpoint_path) && !same_normalized_path(integrated_checkpoint_path, object_path)) {
        unlink(integrated_checkpoint_path, force = TRUE)
    }
    ambiguous_fraction <- mean(object$preliminary_broad_cell_type == "AMBIGUOUS")
    dataset_qc <- data.frame(dataset_id = dataset_id, locked_sample_n = nrow(sample_design), retained_sample_n = data.table::uniqueN(object$sample_accession), 
        retained_cell_n = ncol(object), doublet_mode = doublet_mode, predicted_doublet_n = sum(sample_doublet_qc$predicted_doublet_n), 
        predicted_doublet_fraction = sum(sample_doublet_qc$predicted_doublet_n)/sum(sample_doublet_qc$pre_doublet_cell_n), 
        cluster_n = length(unique(object$seurat_clusters)), preliminary_broad_cell_type_n = length(unique(object$preliminary_broad_cell_type)), 
        ambiguous_cell_fraction = ambiguous_fraction, max_cluster_donor_fraction = max(cluster_totals$max_donor_fraction), 
        pseudobulk_group_n = ncol(pseudobulk$counts), object_path = normalize_slashes(object_path), pseudobulk_matrix_path = normalize_slashes(pseudobulk_matrix_path), 
        stringsAsFactors = FALSE)
    write_csv_utf8(dataset_qc, file.path(dataset_table_root, paste0(dataset_id, "_integration_QC_summary.csv")))
    rm(object, data_matrix, counts_matrix)
    gc()
    list(dataset_qc = dataset_qc, sample_doublet_qc = sample_doublet_qc, cluster_annotation = cluster_annotation, 
        donor_program_summary = donor_program_summary)
}

dataset_results <- list()

for (dataset_id in c("GSE250130", "GSE183837")) {
    dataset_results[[dataset_id]] <- process_dataset(dataset_id = dataset_id, sample_design = designs[[dataset_id]])
}

combined_dataset_qc <- do.call(rbind, lapply(dataset_results, function(x) {
    x$dataset_qc
}))

combined_doublet_qc <- do.call(rbind, lapply(dataset_results, function(x) {
    x$sample_doublet_qc
}))

write_csv_utf8(combined_dataset_qc, file.path(table_root, "STEP04E1_DATASET_QC_SUMMARY.csv"))

write_csv_utf8(combined_doublet_qc, file.path(table_root, "STEP04E1_SAMPLE_DOUBLET_QC.csv"))

all_samples_retained <- all(combined_dataset_qc$locked_sample_n == combined_dataset_qc$retained_sample_n)

ambiguous_pass <- all(combined_dataset_qc$ambiguous_cell_fraction < 0.20000000000000001)

donor_mixing_pass <- all(combined_dataset_qc$max_cluster_donor_fraction < 0.80000000000000004)

doublet_pass <- doublet_mode == "SCDBLFINDER"

gate_status <- if (all_samples_retained && doublet_pass && ambiguous_pass && donor_mixing_pass) {
    "READY_FOR_MANUAL_CELLTYPE_QC_AND_DONOR_AWARE_PSEUDOBULK_REVIEW"
} else if (!doublet_pass) {
    "HOLD_SCDBLFINDER_PACKAGES_MISSING"
} else {
    "HOLD_STEP04E1_BIOLOGICAL_QC_REVIEW_REQUIRED"
}

gate_table <- data.frame(gate_component = c("STEP04D-QC1 authorization", "Dataset integration strategy", "GSE250130 sample retention", 
    "GSE183837 sample retention", "Doublet assessment", "Preliminary broad-cell annotation", "Cluster donor-mixing audit", 
    "Donor-aware VDR-program localization", "Preliminary broad-cell pseudobulk", "GSE183837 watchlist", "Automatic donor exclusion", 
    "Definitive direct-target assignment", "Scientific figures", "STEP04E1 gate"), status = c("PASS", "DATASETS_PROCESSED_SEPARATELY", 
    paste0(combined_dataset_qc$retained_sample_n[combined_dataset_qc$dataset_id == "GSE250130"], "_OF_28"), paste0(combined_dataset_qc$retained_sample_n[combined_dataset_qc$dataset_id == 
        "GSE183837"], "_OF_9"), doublet_mode, ifelse(ambiguous_pass, "PASS_AMBIGUOUS_CELL_FRACTION_LT_0.20", "REVIEW_AMBIGUOUS_CELL_FRACTION"), 
    ifelse(donor_mixing_pass, "PASS_MAX_CLUSTER_DONOR_FRACTION_LT_0.80", "REVIEW_DONOR_DOMINATED_CLUSTER"), "COMPLETED_SAMPLE_BY_BROAD_CELL_TYPE", 
    "BUILT_MINIMUM_50_CELLS_PER_GROUP", "ACTIVE_2_REVIEW_REQUIRED_2_CAUTION", "NONE", "NOT_PERFORMED", "NONE_CREATED", 
    gate_status), required_action = c("Retain corrected v1.3 bulk handoff.", "Do not pool the two studies into one integrated embedding.", 
    "Review temporal-stage representation by donor and cell type.", "Review control/RIF representation and watchlist donors.", 
    "Use scDblFinder singlets for downstream work.", "Inspect cluster markers and revise preliminary labels before locking.", 
    "Inspect any donor-dominated cluster before biological interpretation.", "Use donors, not cells, as biological replicates.", 
    "Do not run final pseudobulk DE until broad-cell labels are locked.", "Inspect GSM5572244 and GSM5572245 most closely; retain GSM5572238 and GSM5572240 as caution donors.", 
    "Do not exclude a donor solely from the pre-QC watchlist.", "Use occupancy-supported candidate terminology.", 
    "Render scientific figures only after cell-type QC and pseudobulk review.", "Proceed only according to the Gate."), 
    stringsAsFactors = FALSE)

gate_path <- file.path(output_root, "STEP04E1_GATE.csv")

write_csv_utf8(gate_table, gate_path)

result_card <- c("# STEP04E1 single-cell integration handoff card", "", paste0("- Run ID: ", run_id), paste0("- Integration version: ", 
    integration_version, " CANDIDATE"), paste0("- Annotation version: ", annotation_version, " PRELIMINARY"), paste0("- Localization version: ", 
    localization_version, " CANDIDATE"), "- Integration strategy: GSE250130 and GSE183837 processed separately", 
    paste0("- Doublet mode: ", doublet_mode), paste0("- GSE250130 retained samples/cells: ", combined_dataset_qc$retained_sample_n[combined_dataset_qc$dataset_id == 
        "GSE250130"], "/", combined_dataset_qc$retained_cell_n[combined_dataset_qc$dataset_id == "GSE250130"]), 
    paste0("- GSE183837 retained samples/cells: ", combined_dataset_qc$retained_sample_n[combined_dataset_qc$dataset_id == 
        "GSE183837"], "/", combined_dataset_qc$retained_cell_n[combined_dataset_qc$dataset_id == "GSE183837"]), 
    paste0("- GSE250130 clusters/ambiguous fraction: ", combined_dataset_qc$cluster_n[combined_dataset_qc$dataset_id == 
        "GSE250130"], "/", format(combined_dataset_qc$ambiguous_cell_fraction[combined_dataset_qc$dataset_id == 
        "GSE250130"], digits = 5)), paste0("- GSE183837 clusters/ambiguous fraction: ", combined_dataset_qc$cluster_n[combined_dataset_qc$dataset_id == 
        "GSE183837"], "/", format(combined_dataset_qc$ambiguous_cell_fraction[combined_dataset_qc$dataset_id == 
        "GSE183837"], digits = 5)), "- Automatic donor exclusions: 0", "- Definitive direct-target assignments: 0", 
    "- Scientific figures: NONE", paste0("- Gate: ", gate_status))

write_text(result_card, file.path(output_root, "STEP04E1_CARD.md"))

write_text(c("# STEP04E1 Methods draft", "", "The two single-cell studies were processed separately because they serve different biological roles: GSE250130 represents a time-series luteal-phase atlas with an LH+7 RIF comparison, whereas GSE183837 provides an independent control-RIF localization cohort. Donor/sample identity was retained throughout.", 
    "", "Within each sample, duplicated gene symbols were collapsed by summing sparse counts. When available, scDblFinder was applied independently to each sample and predicted doublets were removed before merging. Each dataset was log normalized, highly variable genes were selected, and PCA was calculated. Harmony was applied to the PCA representation using sample accession as the batch variable. Shared-nearest-neighbor clustering and UMAP were calculated from the Harmony embedding.", 
    "", "Broad cell compartments were assigned preliminarily from cluster-level canonical marker scores for epithelial, stromal fibroblast, endothelial, immune, pericyte/smooth-muscle and erythroid lineages. Cluster markers, donor mixing and ambiguous-label fractions were retained for manual annotation QC.", 
    "", "Four locked VDR-related transcriptional programs were mapped to single cells using direction-aware signed normalized-expression scores. Scores were summarized at donor-by-broad-cell-type level. Raw counts were aggregated into preliminary donor-by-broad-cell-type pseudobulk matrices for groups containing at least 50 cells. No cell-level significance test was used as evidence of biological replication."), 
    file.path(dirs$text, "STEP04E1_Methods_draft.md"))

write_text(c("# STEP04E1 Results shell", "", "This Results file remains a shell until preliminary broad-cell labels, cluster markers, donor mixing and watchlist-donor behavior have passed manual QC.", 
    "", "No scientific figure, causal claim, diagnostic claim or definitive direct-target statement is authorized at this stage."), 
    file.path(dirs$text, "STEP04E1_Results_shell.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, integration_version = integration_version, annotation_version = annotation_version, 
    localization_version = localization_version, doublet_mode = doublet_mode, variable_feature_n = 3000, PCA_n = 50, 
    harmony_dimension_n = 40, graph_dimension_n = 30, main_cluster_resolution = 0.40000000000000002, pseudobulk_minimum_cell_n = 50, 
    automatic_donor_exclusion_n = 0, definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, 
    gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04E1_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04E1.txt"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) {
    NULL
})

if (!is.null(source_file) && file.exists(source_file)) {
    copy_if_different(source_file, file.path(dirs$code, "53_step04E1_single_cell_integration_and_localization.R"))
}

small_artifacts <- c(gate_path, file.path(output_root, "STEP04E1_CARD.md"), file.path(table_root, "STEP04E1_DATASET_QC_SUMMARY.csv"), 
    file.path(table_root, "STEP04E1_SAMPLE_DOUBLET_QC.csv"), file.path(dirs$text, "STEP04E1_Methods_draft.md"), 
    file.path(dirs$text, "STEP04E1_Results_shell.md"), file.path(dirs$parameters, "STEP04E1_parameters.json"), 
    file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04E1.txt"), file.path(dirs$code, "53_step04E1_single_cell_integration_and_localization.R"), 
    log_path)

dataset_table_artifacts <- list.files(table_root, pattern = "\\.(csv|md)$", recursive = TRUE, full.names = TRUE)

pseudobulk_metadata_artifacts <- list.files(pseudobulk_root, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)

object_artifacts <- list.files(object_root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

pseudobulk_matrix_artifacts <- list.files(pseudobulk_root, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

artifact_paths <- unique(c(small_artifacts, dataset_table_artifacts, pseudobulk_metadata_artifacts, object_artifacts, 
    pseudobulk_matrix_artifacts))

artifact_paths <- artifact_paths[file.exists(artifact_paths)]

artifact_manifest <- data.frame(path = normalize_slashes(artifact_paths), artifact_class = ifelse(artifact_paths %in% 
    object_artifacts, "LARGE_INTEGRATED_OBJECT", ifelse(artifact_paths %in% pseudobulk_matrix_artifacts, "PSEUDOBULK_MATRIX", 
    "SMALL_QC_OR_TEXT_ARTIFACT")), size_bytes = as.numeric(file.info(artifact_paths)$size), md5 = vapply(artifact_paths, 
    md5_file, character(1)), sha256 = vapply(artifact_paths, sha256_file, character(1)), artifact_status = ifelse(gate_status == 
    "READY_FOR_MANUAL_CELLTYPE_QC_AND_DONOR_AWARE_PSEUDOBULK_REVIEW", "CANDIDATE_PENDING_MANUAL_CELLTYPE_QC", "HOLD_PENDING_QC_OR_PACKAGE_COMPLETION"), 
    stringsAsFactors = FALSE)

manifest_path <- file.path(output_root, "STEP04E1_ARTIFACT_MANIFEST.csv")

write_csv_utf8(artifact_manifest, manifest_path)

archive_artifacts <- unique(c(small_artifacts, dataset_table_artifacts, pseudobulk_metadata_artifacts, manifest_path))

archive_artifacts <- archive_artifacts[file.exists(archive_artifacts)]

for (path_value in archive_artifacts) {
    section <- if (grepl("\\.R$", basename(path_value))) {
        "code"
    }
    else if (grepl("sessionInfo|parameters", basename(path_value))) {
        "env"
    }
    else if (grepl("\\.log$", basename(path_value))) {
        "log"
    }
    else if (grepl("\\.md$", basename(path_value))) {
        "txt"
    }
    else if (grepl("GATE|CARD|MANIFEST|QC|MARKER|ANNOTATION|MIXING|COMPOSITION", basename(path_value), ignore.case = TRUE)) {
        "qc"
    }
    else {
        "res"
    }
    destination_path <- file.path(archive_root, section, basename(path_value))
    if (file.exists(destination_path)) {
        destination_path <- file.path(archive_root, section, paste0(basename(dirname(path_value)), "__", basename(path_value)))
    }
    copy_validate(path_value, destination_path)
}

archive_files <- list.files(archive_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(normalize_slashes(archive_files), nchar(normalize_slashes(archive_root)) + 
    2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, character(1)), 
    sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04E1_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04E1 archive status", "", paste0("- Run ID: ", run_id), paste0("- Integration version: ", integration_version, 
    " CANDIDATE"), paste0("- Doublet mode: ", doublet_mode), "- Datasets integrated separately: YES", "- Large Seurat objects copied into short archive: NO", 
    "- Automatic donor exclusions: 0", "- Definitive direct-target assignments: 0", "- Scientific figures: NONE", 
    paste0("- Gate: ", gate_status)), file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Current run: ", run_id), "- Bulk projection: BULK_VDR_PROGRAM_PROJECTION_v1.3 LOCKED", 
    paste0("- Single-cell integration: ", integration_version, " CANDIDATE"), paste0("- Broad-cell annotation: ", 
        annotation_version, " PRELIMINARY"), "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, 
    "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04E1 single-cell integration completed.\n")

cat("Integration version: ", integration_version, " CANDIDATE\n")

cat("Doublet mode: ", doublet_mode, "\n")

for (dataset_id in combined_dataset_qc$dataset_id) {
    dataset_row <- combined_dataset_qc[combined_dataset_qc$dataset_id == dataset_id, , drop = FALSE]
    cat(dataset_id, " samples/cells/clusters: ", dataset_row$retained_sample_n, " / ", dataset_row$retained_cell_n, 
        " / ", dataset_row$cluster_n, "\n", sep = "")
    cat(dataset_id, " ambiguous fraction: ", format(dataset_row$ambiguous_cell_fraction, digits = 6), "\n", sep = "")
}

cat("Automatic donor exclusions: 0\n")

cat("Definitive direct-target assignments: 0\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/S4E1_sc/SCI_v1_0/STEP04E1_GATE.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/STEP04E1_CARD.md\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/STEP04E1_DATASET_QC_SUMMARY.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/STEP04E1_SAMPLE_DOUBLET_QC.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_preliminary_cluster_annotation.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_preliminary_cluster_annotation.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_top25_cluster_markers.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_top25_cluster_markers.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE250130/GSE250130_donor_celltype_VDR_program_scores.csv\n")

cat("  06_results/S4E1_sc/SCI_v1_0/tables/GSE183837/GSE183837_donor_celltype_VDR_program_scores.csv\n")

cat("============================================================\n")

