#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# C02: Continue from the retained post-doublet checkpoint

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(7200L, getOption("timeout")))

options(future.globals.maxSize = 60 * 1024^3)

set.seed(20260804)

project_root <- .vdr_project_root()

dataset_id <- "GSE250130"

script_name <- "integrate_gse250130.R"

script_version <- "STEP04E1B_GSE250130_INTEGRATION_RESUME_v1.0"

integration_version <- "SINGLE_CELL_INTEGRATION_v1.1_MEMORY_SAFE"

output_version <- "SCI_v1_1_memory_safe"

run_start <- Sys.time()

run_id <- paste0("STEP04E1B_", dataset_id, "_", format(run_start, "%Y%m%d_%H%M%S"))

normalize_slashes <- function(x) {
    gsub("\\\\", "/", x)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

safe_read_csv <- function(path) {
    if (!file.exists(path)) 
        return(NULL)
    tryCatch(.vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM"), error = function(e1) {
        .vdr_read_csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    })
}

format_bytes <- function(bytes) {
    if (is.na(bytes)) 
        return("")
    if (bytes < 1024) 
        return(sprintf("%.0f B", bytes))
    if (bytes < 1024^2) 
        return(sprintf("%.1f KB", bytes/1024))
    if (bytes < 1024^3) 
        return(sprintf("%.1f MB", bytes/1024^2))
    sprintf("%.2f GB", bytes/1024^3)
}

safe_package_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) 
        return(NA_character_)
    as.character(utils::packageVersion(pkg))
}

get_free_memory_gb <- function() {
    if (.Platform$OS.type != "windows") 
        return(NA_real_)
    lines <- tryCatch(system(paste0("powershell -NoProfile -Command ", "\"$o=Get-CimInstance Win32_OperatingSystem; ", 
        "[math]::Round($o.FreePhysicalMemory/1MB,2)\""), intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
    if (!length(lines)) 
        return(NA_real_)
    suppressWarnings(as.numeric(trimws(lines[1])))
}

get_free_disk_gb <- function(drive_letter = "C") {
    if (.Platform$OS.type != "windows") 
        return(NA_real_)
    lines <- tryCatch(system(paste0("powershell -NoProfile -Command ", "\"[math]::Round((Get-PSDrive -Name ", drive_letter, 
        ").Free/1GB,2)\""), intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
    if (!length(lines)) 
        return(NA_real_)
    suppressWarnings(as.numeric(trimws(lines[1])))
}

object_size_gb <- function(x) {
    as.numeric(utils::object.size(x))/1024^3
}

md5_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    unname(tools::md5sum(path))
}

log_event <- NULL

run_status <- "STARTED"

failure_stage <- ""

failure_message <- ""

prejoin_saved <- FALSE

data_join_saved <- FALSE

zip_ok <- FALSE

zip_error <- ""

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

dirs <- list(parameters = file.path(step04_root, "03_parameters"), environment = file.path(step04_root, "04_environment"), 
    execution = file.path(step04_root, "05_execution"), results = file.path(step04_root, "06_results"), qc = file.path(step04_root, 
        "08_qc"), archive = file.path(step04_root, "10_archive_work"))

for (p in dirs) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p)) 
        stop("Cannot create required directory: ", p)
}

old_output_root <- file.path(dirs$results, "S4E1_sc", "SCI_v1_0")

new_output_root <- file.path(dirs$results, "S4E1_sc", output_version)

dataset_checkpoint_root <- file.path(new_output_root, "checkpoints", dataset_id)

dataset_table_root <- file.path(new_output_root, "tables", dataset_id)

qc_run_root <- file.path(dirs$qc, paste0(run_id, "_integration_resume"))

for (p in c(new_output_root, dataset_checkpoint_root, dataset_table_root, qc_run_root)) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

run_manifest_path <- file.path(dirs$execution, paste0(run_id, "_RUN_MANIFEST.csv"))

parameter_path <- file.path(dirs$parameters, paste0(run_id, "_parameters.json"))

session_path <- file.path(dirs$environment, paste0("sessionInfo_", run_id, ".txt"))

desktop_candidates <- .vdr_report_root()

desktop_candidates <- .vdr_report_root()

zip_parent <- .vdr_report_root()

handoff_zip <- file.path(zip_parent, paste0(run_id, "_SMALL_REPORT.zip"))

log_event <- function(...) {
    line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(line, "\n")
    cat(line, "\n", file = log_path, append = TRUE)
    flush.console()
}

run_step <- function(label, fun) {
    log_event("START", label)
    t0 <- Sys.time()
    result <- tryCatch(fun(), error = function(e) {
        failure_stage <<- label
        failure_message <<- conditionMessage(e)
        log_event("ERROR", label, failure_message)
        stop(e)
    })
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
    log_event("DONE", label, "elapsed_min=", elapsed)
    result
}

required_packages <- c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "jsonlite")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages)) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
}

source_checkpoint <- file.path(old_output_root, "checkpoints", dataset_id, paste0(dataset_id, "_post_doublet_merged_counts_checkpoint.rds"))

source_status_path <- file.path(old_output_root, "checkpoints", dataset_id, paste0(dataset_id, "_post_doublet_merged_counts_status.csv"))

sample_doublet_qc_path <- file.path(old_output_root, "tables", dataset_id, paste0(dataset_id, "_sample_doublet_qc.csv"))

locked_design_path <- file.path(dirs$results, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0", "COMBINED_SAMPLE_DESIGN_v1.0_GSE250130.csv")

required_inputs <- c(source_checkpoint, source_status_path, sample_doublet_qc_path, locked_design_path)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop("Required input(s) missing: ", paste(missing_inputs, collapse = "; "))
}

source_status <- safe_read_csv(source_status_path)

sample_design <- safe_read_csv(locked_design_path)

sample_doublet_qc <- safe_read_csv(sample_doublet_qc_path)

if (is.null(source_status) || nrow(source_status) != 1L) {
    stop("Invalid source checkpoint status table.")
}

if (source_status$dataset_id[[1]] != dataset_id || source_status$checkpoint_status[[1]] != "COMPLETE" || source_status$sample_n[[1]] != 
    28L || source_status$cell_n[[1]] != 236636L) {
    stop("Source checkpoint status does not match the locked diagnostic facts.")
}

if (is.null(sample_design) || nrow(sample_design) != 28L || !"sample_accession" %in% names(sample_design)) {
    stop("Locked GSE250130 sample design is invalid.")
}

if (is.null(sample_doublet_qc) || nrow(sample_doublet_qc) != 28L) {
    stop("GSE250130 sample doublet QC table is invalid.")
}

free_memory_gb <- get_free_memory_gb()

free_disk_gb <- get_free_disk_gb("C")

resource_snapshot <- data.frame(metric = c("free_physical_memory_GB_before_load", "free_C_drive_GB_before_run", 
    "source_checkpoint_size_GB", "source_checkpoint_md5"), value = c(free_memory_gb, free_disk_gb, round(file.info(source_checkpoint)$size/1024^3, 
    3), md5_file(source_checkpoint)), stringsAsFactors = FALSE)

write_csv_utf8(resource_snapshot, file.path(qc_run_root, "01_RESOURCE_AND_INPUT_SNAPSHOT.csv"))

if (is.finite(free_memory_gb) && free_memory_gb < 18) {
    stop("Free physical memory is below 18 GB. Restart RStudio and rerun in a clean session.")
}

if (is.finite(free_disk_gb) && free_disk_gb < 25) {
    stop("Free C: drive space is below 25 GB. Free disk space before running.")
}

parameter_manifest <- list(script_name = script_name, script_version = script_version, integration_version = integration_version, 
    output_version = output_version, dataset_id = dataset_id, source_checkpoint = normalize_slashes(source_checkpoint), 
    source_checkpoint_status = "COMPLETE", expected_sample_n = 28L, expected_cell_n = 236636L, normalization_method = "LogNormalize", 
    scale_factor = 10000, normalization_block_size = 2000L, variable_feature_method = "vst", variable_feature_n = 3000L, 
    scaling_block_size = 1000L, pca_npcs = 50L, harmony_dims = "1:40", neighbor_dims = "1:30", neighbor_k = 30L, 
    cluster_resolutions = c(0.20000000000000001, 0.59999999999999998, 0.40000000000000002), umap_dims = "1:30", 
    random_seed = 20260804L, layer_join_policy = "JOIN_DATA_ONLY_KEEP_COUNTS_SPLIT", candidate_status = "NOT_LOCKED")

jsonlite::write_json(parameter_manifest, parameter_path, pretty = TRUE, auto_unbox = TRUE)

log_event("Starting", run_id, "source_checkpoint=", source_checkpoint, "free_memory_GB=", free_memory_gb, "free_disk_GB=", 
    free_disk_gb)

prejoin_checkpoint <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_integrated_prejoin_layers_checkpoint.rds"))

prejoin_status_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_integrated_prejoin_layers_status.csv"))

data_join_checkpoint <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_integrated_checkpoint.rds"))

data_join_status_path <- file.path(dataset_checkpoint_root, paste0(dataset_id, "_checkpoint_status.csv"))

for (p in c(prejoin_checkpoint, prejoin_status_path, data_join_checkpoint, data_join_status_path)) {
    if (file.exists(p)) {
        stop("Target already exists. Refusing silent overwrite: ", p)
    }
}

pipeline_error <- NULL

object <- NULL

tryCatch({
    run_status <- "RUNNING"
    object <- run_step("LOAD_POST_DOUBLET_CHECKPOINT", function() readRDS(source_checkpoint))
    if (!inherits(object, "Seurat")) {
        stop("Source checkpoint is not a Seurat object.")
    }
    if (ncol(object) != 236636L) {
        stop("Unexpected cell count after checkpoint load: ", ncol(object))
    }
    if (!"sample_accession" %in% colnames(object@meta.data) || data.table::uniqueN(object$sample_accession) != 
        28L) {
        stop("Checkpoint does not retain all 28 locked samples.")
    }
    SeuratObject::DefaultAssay(object) <- "RNA"
    layer_before <- SeuratObject::Layers(object[["RNA"]])
    write_csv_utf8(data.frame(phase = "SOURCE_POST_DOUBLET", layer = layer_before, stringsAsFactors = FALSE), file.path(qc_run_root, 
        "02_LAYER_MANIFEST_SOURCE.csv"))
    log_event("Loaded source object", "cells=", ncol(object), "features=", nrow(object), "samples=", data.table::uniqueN(object$sample_accession), 
        "object_size_GB=", round(object_size_gb(object), 3), "RNA_layer_n=", length(layer_before))
    object <- run_step("NORMALIZE_DATA_LAYERED", function() {
        Seurat::NormalizeData(object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, 
            block.size = 2000, verbose = TRUE)
    })
    object <- run_step("FIND_VARIABLE_FEATURES", function() {
        Seurat::FindVariableFeatures(object, assay = "RNA", selection.method = "vst", nfeatures = 3000, verbose = TRUE)
    })
    if (length(SeuratObject::VariableFeatures(object)) != 3000L) {
        log_event("WARNING variable feature count=", length(SeuratObject::VariableFeatures(object)))
    }
    object <- run_step("SCALE_VARIABLE_FEATURES", function() {
        Seurat::ScaleData(object, assay = "RNA", features = SeuratObject::VariableFeatures(object), block.size = 1000, 
            verbose = TRUE)
    })
    object <- run_step("RUN_PCA_50", function() {
        Seurat::RunPCA(object, assay = "RNA", features = SeuratObject::VariableFeatures(object), npcs = 50, seed.use = 20260804, 
            verbose = TRUE)
    })
    object <- run_step("RUN_HARMONY_SAMPLE_ACCESSION", function() {
        harmony::RunHarmony(object = object, group.by.vars = "sample_accession", reduction = "pca", dims.use = 1:40, 
            reduction.save = "harmony", plot_convergence = FALSE, verbose = TRUE)
    })
    object <- run_step("FIND_NEIGHBORS_HARMONY", function() {
        Seurat::FindNeighbors(object, reduction = "harmony", dims = 1:30, k.param = 30, verbose = TRUE)
    })
    object <- run_step("FIND_CLUSTERS_RES_0_2", function() {
        Seurat::FindClusters(object, resolution = 0.20000000000000001, algorithm = 1, random.seed = 20260804, verbose = TRUE)
    })
    object <- run_step("FIND_CLUSTERS_RES_0_6", function() {
        Seurat::FindClusters(object, resolution = 0.59999999999999998, algorithm = 1, random.seed = 20260804, verbose = TRUE)
    })
    object <- run_step("FIND_CLUSTERS_RES_0_4", function() {
        Seurat::FindClusters(object, resolution = 0.40000000000000002, algorithm = 1, random.seed = 20260804, verbose = TRUE)
    })
    object <- run_step("RUN_UMAP_HARMONY", function() {
        Seurat::RunUMAP(object, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony", reduction.key = "UMAPH_", 
            seed.use = 20260804, verbose = TRUE)
    })
    reductions_to_keep <- intersect(c("pca", "harmony", "umap_harmony"), SeuratObject::Reductions(object))
    graphs_to_keep <- SeuratObject::Graphs(object)
    layers_before_diet <- SeuratObject::Layers(object[["RNA"]])
    RNA_layers_to_keep <- grep("^(counts|data)(\\.|$)", layers_before_diet, value = TRUE)
    if (!length(RNA_layers_to_keep)) {
        stop("No counts/data RNA layers found before DietSeurat.")
    }
    object <- run_step("DIET_OBJECT_REMOVE_SCALE_DATA", function() {
        Seurat::DietSeurat(object, assays = "RNA", layers = list(RNA = RNA_layers_to_keep), dimreducs = reductions_to_keep, 
            graphs = graphs_to_keep, misc = TRUE)
    })
    gc()
    prejoin_layers <- SeuratObject::Layers(object[["RNA"]])
    cluster_columns <- grep("_snn_res\\.", colnames(object@meta.data), value = TRUE)
    umap_embeddings <- Seurat::Embeddings(object, reduction = "umap_harmony")
    embedding_table <- data.frame(cell_barcode = rownames(umap_embeddings), UMAPH_1 = umap_embeddings[, 1], UMAPH_2 = umap_embeddings[, 
        2], sample_accession = object$sample_accession, seurat_clusters = as.character(object$seurat_clusters), 
        stringsAsFactors = FALSE)
    for (cc in cluster_columns) {
        embedding_table[[cc]] <- object@meta.data[[cc]]
    }
    write_csv_utf8(embedding_table, file.path(dataset_table_root, paste0(dataset_id, "_UMAP_cluster_metadata_prejoin.csv")))
    cluster_sample_counts <- data.table::as.data.table(object@meta.data)[, .(cell_n = .N), by = .(sample_accession, 
        seurat_clusters)]
    write_csv_utf8(as.data.frame(cluster_sample_counts), file.path(dataset_table_root, paste0(dataset_id, "_cluster_by_sample_prejoin.csv")))
    write_csv_utf8(data.frame(phase = "PREJOIN_AFTER_UMAP", layer = prejoin_layers, stringsAsFactors = FALSE), 
        file.path(qc_run_root, "03_LAYER_MANIFEST_PREJOIN.csv"))
    log_event("Saving PRE-JOIN checkpoint", "object_size_GB=", round(object_size_gb(object), 3), "layers=", length(prejoin_layers), 
        "reductions=", paste(reductions_to_keep, collapse = ";"), "graphs=", paste(graphs_to_keep, collapse = ";"))
    run_step("SAVE_PREJOIN_CHECKPOINT", function() {
        saveRDS(object, prejoin_checkpoint, compress = FALSE)
        TRUE
    })
    write_csv_utf8(data.frame(dataset_id = dataset_id, integration_version = integration_version, output_version = output_version, 
        checkpoint_stage = "PREJOIN_AFTER_PCA_HARMONY_CLUSTER_UMAP", checkpoint_status = "COMPLETE", sample_n = data.table::uniqueN(object$sample_accession), 
        cell_n = ncol(object), feature_n = nrow(object), cluster_n = length(unique(object$seurat_clusters)), layer_n = length(prejoin_layers), 
        checkpoint_size_bytes = file.info(prejoin_checkpoint)$size, checkpoint_md5 = md5_file(prejoin_checkpoint), 
        checkpoint_path = normalize_slashes(prejoin_checkpoint), stringsAsFactors = FALSE), prejoin_status_path)
    prejoin_saved <- TRUE
    log_event("PRE-JOIN checkpoint safely saved", prejoin_checkpoint, format_bytes(file.info(prejoin_checkpoint)$size))
    data_layers_before <- SeuratObject::Layers(object[["RNA"]], search = "^data")
    counts_layers_before <- SeuratObject::Layers(object[["RNA"]], search = "^counts")
    if (!length(data_layers_before)) {
        stop("No normalized data layers found for data-only JoinLayers.")
    }
    log_event("Joining DATA layers only", "data_layer_n=", length(data_layers_before), "counts_layer_n_kept_split=", 
        length(counts_layers_before))
    object <- run_step("JOIN_NORMALIZED_DATA_LAYERS_ONLY", function() {
        SeuratObject::JoinLayers(object, assay = "RNA", layers = "data", new = "data")
    })
    gc()
    data_layers_after <- SeuratObject::Layers(object[["RNA"]], search = "^data")
    counts_layers_after <- SeuratObject::Layers(object[["RNA"]], search = "^counts")
    write_csv_utf8(rbind(data.frame(phase = "AFTER_DATA_ONLY_JOIN", layer_type = "data", layer = data_layers_after, 
        stringsAsFactors = FALSE), data.frame(phase = "AFTER_DATA_ONLY_JOIN", layer_type = "counts", layer = counts_layers_after, 
        stringsAsFactors = FALSE)), file.path(qc_run_root, "04_LAYER_MANIFEST_AFTER_DATA_JOIN.csv"))
    if (!identical(data_layers_after, "data")) {
        stop("Data-only JoinLayers did not produce exactly one 'data' layer. Found: ", paste(data_layers_after, 
            collapse = ";"))
    }
    if (!length(counts_layers_after)) {
        stop("Counts layers were unexpectedly lost after data-only join.")
    }
    log_event("Saving DATA-JOINED integrated checkpoint", "object_size_GB=", round(object_size_gb(object), 3), 
        "joined_data_layer=", paste(data_layers_after, collapse = ";"), "counts_layer_n=", length(counts_layers_after))
    run_step("SAVE_DATA_JOINED_INTEGRATED_CHECKPOINT", function() {
        saveRDS(object, data_join_checkpoint, compress = FALSE)
        TRUE
    })
    write_csv_utf8(data.frame(dataset_id = dataset_id, integration_version = integration_version, output_version = output_version, 
        doublet_mode = "SCDBLFINDER_REUSED_FROM_SCI_v1_0", sample_n = data.table::uniqueN(object$sample_accession), 
        cell_n = ncol(object), feature_n = nrow(object), cluster_n = length(unique(object$seurat_clusters)), checkpoint_status = "COMPLETE", 
        join_policy = "DATA_ONLY_COUNTS_REMAIN_SPLIT", joined_data_layer_n = length(data_layers_after), remaining_counts_layer_n = length(counts_layers_after), 
        checkpoint_size_bytes = file.info(data_join_checkpoint)$size, checkpoint_md5 = md5_file(data_join_checkpoint), 
        checkpoint_path = normalize_slashes(data_join_checkpoint), stringsAsFactors = FALSE), data_join_status_path)
    data_join_saved <- TRUE
    run_status <- "PASS_GSE250130_INTEGRATION_AND_DATA_JOIN_COMPLETE"
    log_event("PASS", run_status, data_join_checkpoint, format_bytes(file.info(data_join_checkpoint)$size))
}, error = function(e) {
    pipeline_error <<- e
    if (prejoin_saved && !data_join_saved) {
        run_status <<- "HOLD_PREJOIN_SAVED_DATA_JOIN_FAILED"
    }
    else if (!prejoin_saved) {
        run_status <<- "FAIL_BEFORE_PREJOIN_CHECKPOINT"
    }
    else {
        run_status <<- "FAIL_AFTER_DATA_JOIN_CHECKPOINT"
    }
    log_event("PIPELINE_FINAL_STATUS", run_status, "failure_stage=", failure_stage, "message=", failure_message)
})

run_end <- Sys.time()

package_table <- data.frame(package = c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "jsonlite", 
    "future"), version = vapply(c("Seurat", "SeuratObject", "Matrix", "harmony", "data.table", "jsonlite", "future"), 
    safe_package_version, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(package_table, file.path(qc_run_root, "05_PACKAGE_VERSIONS.csv"))

write_text(capture.output(sessionInfo()), session_path)

final_checkpoint_audit <- data.frame(artifact = c("source_post_doublet_checkpoint", "prejoin_integrated_checkpoint", 
    "data_joined_integrated_checkpoint"), exists = c(file.exists(source_checkpoint), file.exists(prejoin_checkpoint), 
    file.exists(data_join_checkpoint)), size_bytes = c(if (file.exists(source_checkpoint)) file.info(source_checkpoint)$size else NA_real_, 
    if (file.exists(prejoin_checkpoint)) file.info(prejoin_checkpoint)$size else NA_real_, if (file.exists(data_join_checkpoint)) file.info(data_join_checkpoint)$size else NA_real_), 
    size_human = c(if (file.exists(source_checkpoint)) format_bytes(file.info(source_checkpoint)$size) else "", 
        if (file.exists(prejoin_checkpoint)) format_bytes(file.info(prejoin_checkpoint)$size) else "", if (file.exists(data_join_checkpoint)) format_bytes(file.info(data_join_checkpoint)$size) else ""), 
    md5 = c(md5_file(source_checkpoint), md5_file(prejoin_checkpoint), md5_file(data_join_checkpoint)), path = normalize_slashes(c(source_checkpoint, 
        prejoin_checkpoint, data_join_checkpoint)), stringsAsFactors = FALSE)

write_csv_utf8(final_checkpoint_audit, file.path(qc_run_root, "06_FINAL_CHECKPOINT_AUDIT.csv"))

run_manifest <- data.frame(run_id = run_id, script_name = script_name, script_version = script_version, integration_version = integration_version, 
    output_version = output_version, dataset_id = dataset_id, started_at = as.character(run_start), completed_at = as.character(run_end), 
    elapsed_minutes = round(as.numeric(difftime(run_end, run_start, units = "mins")), 2), run_status = run_status, 
    failure_stage = failure_stage, failure_message = failure_message, source_checkpoint_reused = TRUE, sample_doublet_processing_repeated = FALSE, 
    prejoin_checkpoint_saved = prejoin_saved, data_join_checkpoint_saved = data_join_saved, old_output_root_modified = FALSE, 
    new_output_root = normalize_slashes(new_output_root), stringsAsFactors = FALSE)

write_csv_utf8(run_manifest, run_manifest_path)

summary_lines <- c("# STEP04E1B GSE250130 integration resume", "", paste0("- Run ID: `", run_id, "`"), paste0("- Status: **", 
    run_status, "**"), paste0("- Source checkpoint reused: `", normalize_slashes(source_checkpoint), "`"), "- Sample-level doublet detection was not repeated.", 
    paste0("- Pre-join checkpoint saved: ", prejoin_saved), paste0("- Data-joined checkpoint saved: ", data_join_saved), 
    paste0("- Failure stage: ", ifelse(nzchar(failure_stage), failure_stage, "none")), paste0("- Failure message: ", 
        ifelse(nzchar(failure_message), failure_message, "none")), "", "## Interpretation", "", if (identical(run_status, 
        "PASS_GSE250130_INTEGRATION_AND_DATA_JOIN_COMPLETE")) {
        "- GSE250130 normalization, PCA, Harmony, graph, clustering, UMAP, and normalized-data layer joining completed."
    } else if (identical(run_status, "HOLD_PREJOIN_SAVED_DATA_JOIN_FAILED")) {
        "- Integration completed and is safely checkpointed before JoinLayers; only data-layer joining remains unresolved."
    } else {
        "- Integration did not reach a safe pre-join checkpoint; inspect the execution log before rerunning."
    }, "", "## Next step", "", "- Do not start GSE183837 yet.", "- Return the small ZIP and review the layer/checkpoint audit.", 
    "- After PASS, run the separate annotation/VDR-localization/pseudobulk script.", "- Results remain CANDIDATE and are not yet eligible for locked results.")

write_text(summary_lines, file.path(qc_run_root, "00_STEP04E1B_RUN_SUMMARY.md"))

small_handoff_root <- file.path(qc_run_root, "small_handoff")

dir.create(small_handoff_root, recursive = TRUE, showWarnings = FALSE)

small_files <- unique(c(log_path, run_manifest_path, parameter_path, session_path, list.files(qc_run_root, full.names = TRUE, 
    recursive = FALSE), list.files(dataset_table_root, pattern = "\\.csv$", full.names = TRUE), prejoin_status_path, 
    data_join_status_path, normalizePath(.vdr_script_path(script_name), winslash = "/", mustWork = FALSE)))

small_files <- small_files[file.exists(small_files) & file.info(small_files)$size <= 100 * 1024^2]

copy_rows <- list()

for (i in seq_along(small_files)) {
    src <- small_files[i]
    dst <- file.path(small_handoff_root, paste0(sprintf("%03d", i), "_", basename(src)))
    ok <- file.copy(src, dst, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    copy_rows[[i]] <- data.frame(source = normalize_slashes(src), destination = normalize_slashes(dst), copied = isTRUE(ok), 
        stringsAsFactors = FALSE)
}

if (length(copy_rows)) {
    write_csv_utf8(do.call(rbind, copy_rows), file.path(qc_run_root, "07_SMALL_HANDOFF_COPY_MANIFEST.csv"))
}

if (!requireNamespace("zip", quietly = TRUE)) {
    try(install.packages("zip", repos = "https://cloud.r-project.org"), silent = TRUE)
}

if (requireNamespace("zip", quietly = TRUE)) {
    zip_ok <- tryCatch({
        zip::zipr(zipfile = handoff_zip, files = c(file.path(qc_run_root, "00_STEP04E1B_RUN_SUMMARY.md"), file.path(qc_run_root, 
            "01_RESOURCE_AND_INPUT_SNAPSHOT.csv"), file.path(qc_run_root, "02_LAYER_MANIFEST_SOURCE.csv"), file.path(qc_run_root, 
            "03_LAYER_MANIFEST_PREJOIN.csv"), file.path(qc_run_root, "04_LAYER_MANIFEST_AFTER_DATA_JOIN.csv"), 
            file.path(qc_run_root, "05_PACKAGE_VERSIONS.csv"), file.path(qc_run_root, "06_FINAL_CHECKPOINT_AUDIT.csv"), 
            file.path(qc_run_root, "07_SMALL_HANDOFF_COPY_MANIFEST.csv"), small_handoff_root), recurse = TRUE, 
            include_directories = TRUE)
        file.exists(handoff_zip)
    }, error = function(e) {
        zip_error <<- conditionMessage(e)
        FALSE
    })
}

cat("\n============================================================\n")

cat("STEP04E1B GSE250130 resume finished\n")

cat("============================================================\n")

cat("Run status: ", run_status, "\n", sep = "")

cat("Execution log: ", log_path, "\n", sep = "")

cat("Run manifest: ", run_manifest_path, "\n", sep = "")

cat("QC folder: ", qc_run_root, "\n", sep = "")

cat("Pre-join checkpoint saved: ", prejoin_saved, "\n", sep = "")

cat("Data-joined checkpoint saved: ", data_join_saved, "\n", sep = "")

if (zip_ok) {
    cat("Upload this small ZIP:\n", handoff_zip, "\n")
} else {
    cat("ZIP creation failed: ", zip_error, "\nThe QC folder is complete; compress it manually if needed.\n", sep = "")
}

