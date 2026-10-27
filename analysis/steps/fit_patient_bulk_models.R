#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# B04: Build gene matrices, project programs, fit bulk models and prepare single-cell inputs

options(stringsAsFactors = FALSE)

options(warn = 1)

options(timeout = max(3600L, getOption("timeout")))

project_root <- .vdr_project_root()

step04_root <- file.path(project_root, "05_analysis_steps", "STEP04_combined_temporal_patient_single_cell")

run_id <- paste0("STEP04D_v6_", format(Sys.time(), "%Y%m%d_%H%M%S"))

bulk_matrix_version <- "BULK_GENE_MATRIX_v1.2"

bulk_projection_version <- "BULK_VDR_PROGRAM_PROJECTION_v1.2"

sc_preQC_version <- "SINGLE_CELL_PREQC_OBJECTS_v1.0"

GSE92324_profile_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92324/suppl/GSE92324_sample_probe_profile.txt.gz"

GSE92324_profile_path <- file.path(project_root, "03_data", "raw_external_READ_ONLY", "GSE92324", "supplementary", 
    "GSE92324_sample_probe_profile.txt.gz")

required_packages <- c("data.table", "digest", "jsonlite", "limma", "Matrix")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0L) {
    stop("Required package(s) missing: ", paste(missing_packages, collapse = ", "), ". This script does not install or update packages.")
}

dirs <- list(control = file.path(step04_root, "00_module_control"), inputs = file.path(step04_root, "01_inputs"), 
    code = file.path(step04_root, "02_code"), parameters = file.path(step04_root, "03_parameters"), environment = file.path(step04_root, 
        "04_environment"), execution = file.path(step04_root, "05_execution"), results = file.path(step04_root, 
        "06_results"), figures = file.path(step04_root, "07_figures_tables"), qc = file.path(step04_root, "08_qc"), 
    text = file.path(step04_root, "09_manuscript_text"), archive_parent = file.path(step04_root, "10_archive_work"))

for (directory_path in dirs) {
    dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
}

bulk_matrix_root <- file.path(dirs$results, "S4D_bulk_matrix", "BGM_v1_2")

bulk_result_root <- file.path(dirs$results, "S4D_bulk", "BVP_v1_2")

sc_object_root <- file.path(dirs$results, "S4D_sc_preQC", "SCP_v1_0")

previous_bulk_matrix_root <- file.path(dirs$results, "STEP04D_bulk_gene_matrices", "BULK_GENE_MATRIX_v1.1")

for (directory_path in c(bulk_matrix_root, bulk_result_root, sc_object_root)) {
    directory_created <- dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory_path)) {
        stop("Required STEP04D output directory could not be created: ", directory_path)
    }
    if (nchar(directory_path) >= 220L) {
        stop("STEP04D output directory exceeds the 220-character soft limit: ", directory_path)
    }
}

archive_root <- file.path(dirs$archive_parent, paste0("S4D_", format(Sys.time(), "%Y%m%d_%H%M%S")))

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
    if (nchar(path) >= 220L) {
        stop("Output path exceeds the 220-character soft limit: ", path)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(dirname(path))) {
        stop("Output directory does not exist and could not be created: ", dirname(path))
    }
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
        stop("Copy/hash validation failed:\n", source_path)
    }
    TRUE
}

download_with_resume <- function(url, destination) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(destination) && as.numeric(file.info(destination)$size) > 10000) {
        return(TRUE)
    }
    curl_candidates <- c(Sys.which("curl.exe"), Sys.which("curl"))
    curl_candidates <- curl_candidates[nzchar(curl_candidates)]
    if (length(curl_candidates) > 0L) {
        status <- system2(curl_candidates[[1]], args = c("--fail", "--location", "--retry", "15", "--retry-all-errors", 
            "--retry-delay", "5", "--connect-timeout", "60", "--continue-at", "-", "--output", shQuote(destination), 
            shQuote(url)))
        if (status != 0L || !file.exists(destination) || as.numeric(file.info(destination)$size) <= 10000) {
            stop("Download failed: ", url)
        }
        return(TRUE)
    }
    download.file(url = url, destfile = destination, mode = "wb", method = "libcurl", quiet = FALSE)
    if (!file.exists(destination) || as.numeric(file.info(destination)$size) <= 10000) {
        stop("Download failed: ", url)
    }
    TRUE
}

read_tabular_maybe_gz <- function(path) {
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    read.delim(connection, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "", 
        fill = TRUE)
}

extract_soft_table <- function(path, begin_marker, end_marker) {
    if (!file.exists(path)) {
        stop("SOFT/matrix source is missing: ", path)
    }
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    temp_path <- tempfile(fileext = ".tsv")
    output_connection <- file(temp_path, open = "wt")
    output_connection_open <- TRUE
    on.exit({
        if (isTRUE(output_connection_open)) {
            try(close(output_connection), silent = TRUE)
        }
        if (file.exists(temp_path)) {
            unlink(temp_path, force = TRUE)
        }
    }, add = TRUE)
    in_table <- FALSE
    found_begin <- FALSE
    found_end <- FALSE
    repeat {
        lines <- readLines(connection, n = 5000L, warn = FALSE, encoding = "UTF-8")
        if (length(lines) == 0L) {
            break
        }
        if (!in_table) {
            begin_index <- which(trimws(lines) == begin_marker)
            if (length(begin_index) > 0L) {
                found_begin <- TRUE
                in_table <- TRUE
                start_index <- begin_index[[1]] + 1L
                if (start_index <= length(lines)) {
                  lines <- lines[start_index:length(lines)]
                }
                else {
                  lines <- character(0)
                }
            }
            else {
                next
            }
        }
        if (in_table && length(lines) > 0L) {
            end_index <- which(trimws(lines) == end_marker)
            if (length(end_index) > 0L) {
                end_position <- end_index[[1]] - 1L
                if (end_position >= 1L) {
                  writeLines(lines[seq_len(end_position)], output_connection, useBytes = TRUE)
                }
                found_end <- TRUE
                break
            }
            writeLines(lines, output_connection, useBytes = TRUE)
        }
    }
    close(output_connection)
    output_connection_open <- FALSE
    if (!found_begin || !found_end) {
        stop("Could not locate complete table markers in: ", path)
    }
    table_value <- .vdr_fread(temp_path, data.table = FALSE, showProgress = FALSE, check.names = FALSE)
    table_value
}

read_series_expression <- function(path) {
    table_value <- extract_soft_table(path = path, begin_marker = "!series_matrix_table_begin", end_marker = "!series_matrix_table_end")
    if (ncol(table_value) < 3L) {
        stop("Series matrix has fewer than three columns: ", path)
    }
    probe_id <- trimws(gsub("^\"|\"$", "", as.character(table_value[[1]])))
    expression_value <- as.matrix(table_value[, -1, drop = FALSE])
    storage.mode(expression_value) <- "numeric"
    rownames(expression_value) <- probe_id
    sample_names <- gsub("^\"|\"$", "", colnames(expression_value))
    colnames(expression_value) <- sample_names
    list(expression = expression_value, probe_id = probe_id)
}

read_platform_table <- function(path) {
    extract_soft_table(path = path, begin_marker = "!platform_table_begin", end_marker = "!platform_table_end")
}

choose_column <- function(table_value, candidate_names) {
    normalized_names <- toupper(gsub("[^A-Z0-9]+", "_", names(table_value)))
    normalized_candidates <- toupper(gsub("[^A-Z0-9]+", "_", candidate_names))
    match_index <- match(normalized_candidates, normalized_names)
    match_index <- match_index[!is.na(match_index)]
    if (length(match_index) == 0L) {
        return(NA_character_)
    }
    names(table_value)[match_index[[1]]]
}

clean_single_gene_symbol <- function(x) {
    values <- as.character(x)
    vapply(values, function(value) {
        if (is.na(value) || !nzchar(trimws(value)) || trimws(value) %in% c("---", "NA", "N/A", ".")) {
            return(NA_character_)
        }
        tokens <- unlist(strsplit(value, "\\s*(///|//|;|\\||,)\\s*", perl = TRUE))
        tokens <- unique(trimws(tokens))
        tokens <- tokens[nzchar(tokens) & !tokens %in% c("---", "NA", "N/A", ".") & grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", 
            tokens)]
        if (length(tokens) != 1L) {
            return(NA_character_)
        }
        tokens[[1]]
    }, character(1))
}

evaluate_annotation_table <- function(annotation_table, expression_probe_ids, symbol_candidates, source_label) {
    if (nrow(annotation_table) == 0L || ncol(annotation_table) < 2L) {
        stop("Annotation table is empty or has fewer than two columns: ", source_label)
    }
    expression_probe_ids <- trimws(as.character(expression_probe_ids))
    normalized_names <- toupper(gsub("[^A-Z0-9]+", "_", names(annotation_table)))
    id_preferred <- normalized_names %in% c("ID", "ID_REF", "PROBE_ID", "PROBEID", "ARRAY_ADDRESS_ID", "ARRAYADDRESSID")
    id_overlap_n <- vapply(seq_len(ncol(annotation_table)), function(column_index) {
        values <- trimws(gsub("^\"|\"$", "", as.character(annotation_table[[column_index]])))
        length(intersect(unique(values), expression_probe_ids))
    }, integer(1))
    if (any(id_preferred)) {
        preferred_indices <- which(id_preferred)
        preferred_best <- preferred_indices[which.max(id_overlap_n[preferred_indices])]
        if (id_overlap_n[preferred_best] >= max(1000L, floor(0.02 * length(expression_probe_ids)))) {
            id_column_index <- preferred_best
        }
        else {
            id_column_index <- which.max(id_overlap_n)
        }
    }
    else {
        id_column_index <- which.max(id_overlap_n)
    }
    id_overlap_best <- id_overlap_n[id_column_index]
    if (id_overlap_best < max(1000L, floor(0.02 * length(expression_probe_ids)))) {
        stop("No annotation ID column had adequate overlap with expression probes for ", source_label, ". Best overlap: ", 
            id_overlap_best)
    }
    id_column <- names(annotation_table)[id_column_index]
    annotation_probe_ids <- trimws(gsub("^\"|\"$", "", as.character(annotation_table[[id_column_index]])))
    normalized_candidates <- toupper(gsub("[^A-Z0-9]+", "_", symbol_candidates))
    symbol_candidate_index <- which(normalized_names %in% normalized_candidates | grepl("SYMBOL", normalized_names, 
        fixed = TRUE) | normalized_names %in% c("ILMN_GENE", "GENE", "GENE_NAME"))
    symbol_candidate_index <- setdiff(symbol_candidate_index, id_column_index)
    if (length(symbol_candidate_index) == 0L) {
        stop("No plausible gene-symbol column was found for ", source_label, ".")
    }
    expression_match_index <- match(expression_probe_ids, annotation_probe_ids)
    symbol_audit_rows <- list()
    for (column_index in symbol_candidate_index) {
        raw_symbol <- as.character(annotation_table[[column_index]])
        mapped_raw_symbol <- raw_symbol[expression_match_index]
        cleaned_symbol <- clean_single_gene_symbol(mapped_raw_symbol)
        symbol_audit_rows[[as.character(column_index)]] <- data.frame(annotation_source = source_label, id_column = id_column, 
            symbol_column = names(annotation_table)[column_index], probe_overlap_n = id_overlap_best, accepted_probe_n = sum(!is.na(cleaned_symbol)), 
            unique_gene_n = length(unique(cleaned_symbol[!is.na(cleaned_symbol)])), stringsAsFactors = FALSE)
    }
    symbol_audit <- do.call(rbind, symbol_audit_rows)
    best_row <- symbol_audit[order(-symbol_audit$unique_gene_n, -symbol_audit$accepted_probe_n), , drop = FALSE][1L, 
        , drop = FALSE]
    best_symbol_column <- best_row$symbol_column[[1]]
    platform_map <- data.frame(probe_id = annotation_probe_ids, raw_gene_symbol = as.character(annotation_table[[best_symbol_column]]), 
        stringsAsFactors = FALSE)
    platform_map$gene_symbol <- clean_single_gene_symbol(platform_map$raw_gene_symbol)
    platform_map <- platform_map[!is.na(platform_map$probe_id) & nzchar(platform_map$probe_id), , drop = FALSE]
    platform_map <- platform_map[!duplicated(platform_map$probe_id), , drop = FALSE]
    list(platform_map = platform_map, audit = symbol_audit, selected_id_column = id_column, selected_symbol_column = best_symbol_column, 
        selected_unique_gene_n = best_row$unique_gene_n[[1]], selected_accepted_probe_n = best_row$accepted_probe_n[[1]], 
        source_label = source_label, annotation_row_n = nrow(annotation_table))
}

transform_expression_if_needed <- function(expression_value) {
    finite_values <- as.numeric(expression_value[is.finite(expression_value)])
    if (length(finite_values) == 0L) {
        stop("Expression matrix has no finite values.")
    }
    quantiles <- quantile(finite_values, probs = c(0.01, 0.5, 0.98999999999999999), na.rm = TRUE, names = FALSE)
    log2_required <- quantiles[[3]] > 100 || (quantiles[[3]] - quantiles[[1]]) > 50
    transformed <- expression_value
    if (log2_required) {
        transformed[transformed < 0] <- NA_real_
        transformed <- log2(transformed + 1)
    }
    list(expression = transformed, log2_transform_applied = log2_required, q01 = quantiles[[1]], median = quantiles[[2]], 
        q99 = quantiles[[3]])
}

build_gene_matrix <- function(dataset_id, series_path, platform_path, sample_design, symbol_candidates, supplementary_annotation_path = NA_character_) {
    log_event("Building bulk gene matrix for", dataset_id)
    series_value <- read_series_expression(series_path)
    transformed_value <- transform_expression_if_needed(series_value$expression)
    expression_value <- transformed_value$expression
    if (!setequal(colnames(expression_value), sample_design$sample_accession)) {
        stop(dataset_id, " expression sample headers do not match the locked sample design.")
    }
    expression_value <- expression_value[, sample_design$sample_accession, drop = FALSE]
    platform_value <- read_platform_table(platform_path)
    annotation_candidates <- list(GEO_PLATFORM_TABLE = evaluate_annotation_table(annotation_table = platform_value, 
        expression_probe_ids = rownames(expression_value), symbol_candidates = symbol_candidates, source_label = "GEO_PLATFORM_TABLE"))
    if (!is.na(supplementary_annotation_path) && nzchar(supplementary_annotation_path) && file.exists(supplementary_annotation_path)) {
        supplementary_value <- read_tabular_maybe_gz(supplementary_annotation_path)
        annotation_candidates[["AUTHOR_SUPPLEMENTARY_PROFILE"]] <- evaluate_annotation_table(annotation_table = supplementary_value, 
            expression_probe_ids = rownames(expression_value), symbol_candidates = c(symbol_candidates, "Gene Symbol", 
                "Gene_Symbol", "SYMBOL", "ILMN_Gene", "Gene"), source_label = "AUTHOR_SUPPLEMENTARY_PROFILE")
    }
    annotation_candidate_summary <- do.call(rbind, lapply(annotation_candidates, function(candidate_value) {
        data.frame(dataset_id = dataset_id, annotation_source = candidate_value$source_label, annotation_row_n = candidate_value$annotation_row_n, 
            selected_id_column = candidate_value$selected_id_column, selected_symbol_column = candidate_value$selected_symbol_column, 
            accepted_probe_n = candidate_value$selected_accepted_probe_n, unique_gene_n = candidate_value$selected_unique_gene_n, 
            stringsAsFactors = FALSE)
    }))
    selected_source <- annotation_candidate_summary$annotation_source[order(-annotation_candidate_summary$unique_gene_n, 
        -annotation_candidate_summary$accepted_probe_n)][[1]]
    selected_candidate <- annotation_candidates[[selected_source]]
    platform_map <- selected_candidate$platform_map
    probe_map <- platform_map[match(rownames(expression_value), platform_map$probe_id), , drop = FALSE]
    unique_symbol_pass <- !is.na(probe_map$gene_symbol) & nzchar(probe_map$gene_symbol)
    mapped_expression_all <- expression_value[unique_symbol_pass, , drop = FALSE]
    mapped_symbols_all <- probe_map$gene_symbol[unique_symbol_pass]
    complete_probe_pass <- rowSums(is.finite(mapped_expression_all)) == ncol(mapped_expression_all)
    mapped_expression <- mapped_expression_all[complete_probe_pass, , drop = FALSE]
    mapped_symbols <- mapped_symbols_all[complete_probe_pass]
    mapping_audit_path <- file.path(dirs$qc, paste0("STEP04D_v3_", dataset_id, "_ANNOTATION_SOURCE_COMPARISON.csv"))
    write_csv_utf8(annotation_candidate_summary, mapping_audit_path)
    symbol_column_audit <- do.call(rbind, lapply(annotation_candidates, function(candidate_value) {
        candidate_value$audit
    }))
    write_csv_utf8(symbol_column_audit, file.path(dirs$qc, paste0("STEP04D_v3_", dataset_id, "_SYMBOL_COLUMN_AUDIT.csv")))
    if (nrow(mapped_expression) < 5000L || length(unique(mapped_symbols)) < 5000L) {
        stop(dataset_id, " yielded fewer than 5,000 uniquely symbol-mapped probes or genes after evaluating all available annotation sources. Review: ", 
            mapping_audit_path)
    }
    gene_expression_before_model_filter <- limma::avereps(mapped_expression, ID = mapped_symbols)
    gene_complete_pass <- rowSums(is.finite(gene_expression_before_model_filter)) == ncol(gene_expression_before_model_filter)
    gene_variance <- apply(gene_expression_before_model_filter, 1L, stats::var)
    gene_variable_pass <- is.finite(gene_variance) & gene_variance > 0
    gene_model_pass <- gene_complete_pass & gene_variable_pass
    gene_expression <- gene_expression_before_model_filter[gene_model_pass, , drop = FALSE]
    gene_expression <- gene_expression[order(rownames(gene_expression)), , drop = FALSE]
    qc <- data.frame(dataset_id = dataset_id, input_probe_n = nrow(expression_value), selected_annotation_source = selected_source, 
        selected_annotation_row_n = selected_candidate$annotation_row_n, matched_unique_symbol_probe_n_before_complete_filter = nrow(mapped_expression_all), 
        complete_mapped_probe_n = nrow(mapped_expression), incomplete_mapped_probe_removed_n = sum(!complete_probe_pass), 
        gene_n_before_model_filter = nrow(gene_expression_before_model_filter), incomplete_gene_removed_n = sum(!gene_complete_pass), 
        invariant_gene_removed_n = sum(gene_complete_pass & !gene_variable_pass), gene_n = nrow(gene_expression), 
        sample_n = ncol(gene_expression), log2_transform_applied = transformed_value$log2_transform_applied, input_q01 = transformed_value$q01, 
        input_median = transformed_value$median, input_q99 = transformed_value$q99, platform_id_column = selected_candidate$selected_id_column, 
        platform_symbol_column = selected_candidate$selected_symbol_column, stringsAsFactors = FALSE)
    complete_expression_accepted_full <- rep(FALSE, length(unique_symbol_pass))
    complete_expression_accepted_full[unique_symbol_pass] <- complete_probe_pass
    list(gene_expression = gene_expression, qc = qc, probe_map = data.frame(probe_id = rownames(expression_value), 
        gene_symbol = probe_map$gene_symbol, unique_symbol_accepted = unique_symbol_pass, complete_expression_accepted = complete_expression_accepted_full, 
        selected_annotation_source = selected_source, stringsAsFactors = FALSE), annotation_candidate_summary = annotation_candidate_summary, 
        symbol_column_audit = symbol_column_audit)
}

fit_limma_contrasts <- function(expression_value, design_matrix, contrast_matrix, dataset_id) {
    if (ncol(expression_value) != nrow(design_matrix)) {
        stop(dataset_id, " expression/design dimensions do not match.")
    }
    complete_row_pass <- rowSums(is.finite(expression_value)) == ncol(expression_value)
    row_variance <- apply(expression_value, 1L, stats::var)
    variable_row_pass <- is.finite(row_variance) & row_variance > 0
    model_row_pass <- complete_row_pass & variable_row_pass
    model_expression <- expression_value[model_row_pass, , drop = FALSE]
    if (nrow(model_expression) == 0L) {
        stop(dataset_id, " has no complete variable rows for limma modeling.")
    }
    log_event("limma input", dataset_id, "rows:", nrow(expression_value), "retained:", nrow(model_expression), 
        "removed_nonfinite:", sum(!complete_row_pass), "removed_invariant:", sum(complete_row_pass & !variable_row_pass))
    fit <- limma::lmFit(model_expression, design_matrix)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    finite_sigma_pass <- is.finite(fit$sigma) & fit$sigma > 0 & is.finite(fit$df.residual) & fit$df.residual > 
        0
    if (!all(finite_sigma_pass)) {
        log_event("limma post-fit filtering", dataset_id, "removed_invalid_sigma_or_df:", sum(!finite_sigma_pass))
        fit <- fit[finite_sigma_pass, ]
    }
    if (nrow(fit$coefficients) == 0L) {
        stop(dataset_id, " has no valid rows after limma residual-variance QC.")
    }
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
    result_rows <- list()
    for (contrast_name in colnames(contrast_matrix)) {
        table_value <- limma::topTable(fit, coef = contrast_name, number = Inf, sort.by = "none")
        table_value$gene_symbol <- rownames(table_value)
        table_value$dataset_id <- dataset_id
        table_value$contrast <- contrast_name
        result_rows[[contrast_name]] <- table_value
    }
    do.call(rbind, result_rows)
}

row_zscore <- function(expression_value) {
    row_mean <- rowMeans(expression_value, na.rm = TRUE)
    row_sd <- apply(expression_value, 1L, sd, na.rm = TRUE)
    valid <- is.finite(row_sd) & row_sd > 0
    output <- matrix(NA_real_, nrow = nrow(expression_value), ncol = ncol(expression_value), dimnames = dimnames(expression_value))
    output[valid, ] <- (expression_value[valid, , drop = FALSE] - row_mean[valid])/row_sd[valid]
    output
}

calculate_program_scores <- function(gene_expression, candidate_table, program_definitions, dataset_id) {
    z_expression <- row_zscore(gene_expression)
    score_rows <- list()
    coverage_rows <- list()
    for (program_name in names(program_definitions)) {
        program_genes <- program_definitions[[program_name]]
        program_genes <- unique(program_genes[!is.na(program_genes) & nzchar(program_genes)])
        present_genes <- intersect(program_genes, rownames(z_expression))
        candidate_subset <- candidate_table[match(present_genes, candidate_table$feature_id), , drop = FALSE]
        weights <- sign(candidate_subset$interaction_logFC)
        valid_weight <- is.finite(weights) & weights != 0
        present_genes <- present_genes[valid_weight]
        weights <- weights[valid_weight]
        if (length(present_genes) < 10L) {
            stop(dataset_id, " program ", program_name, " has fewer than 10 measurable genes.")
        }
        score_value <- colMeans(sweep(z_expression[present_genes, , drop = FALSE], 1L, weights, "*"), na.rm = TRUE)
        score_rows[[program_name]] <- data.frame(dataset_id = dataset_id, sample_accession = colnames(gene_expression), 
            program_name = program_name, score = as.numeric(score_value), measurable_gene_n = length(present_genes), 
            stringsAsFactors = FALSE)
        coverage_rows[[program_name]] <- data.frame(dataset_id = dataset_id, program_name = program_name, frozen_gene_n = length(program_genes), 
            measurable_gene_n = length(present_genes), coverage_fraction = length(present_genes)/length(program_genes), 
            stringsAsFactors = FALSE)
    }
    list(scores = do.call(rbind, score_rows), coverage = do.call(rbind, coverage_rows))
}

fit_score_model <- function(score_table, sample_design, design_matrix, contrast_matrix, dataset_id) {
    programs <- unique(score_table$program_name)
    result_rows <- list()
    for (program_name in programs) {
        program_scores <- score_table[score_table$program_name == program_name, , drop = FALSE]
        program_scores <- program_scores[match(sample_design$sample_accession, program_scores$sample_accession), 
            , drop = FALSE]
        if (any(is.na(program_scores$score))) {
            stop(dataset_id, " score alignment failed for ", program_name)
        }
        score_matrix <- matrix(program_scores$score, nrow = 1L, dimnames = list(program_name, sample_design$sample_accession))
        fitted <- fit_limma_contrasts(score_matrix, design_matrix, contrast_matrix, dataset_id)
        fitted$program_name <- fitted$gene_symbol
        fitted$gene_symbol <- NULL
        result_rows[[program_name]] <- fitted
    }
    do.call(rbind, result_rows)
}

read_text_table_maybe_gz <- function(path, header = FALSE) {
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    read.delim(connection, header = header, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

read_mtx_maybe_gz <- function(path) {
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rb")
    }
    else {
        file(path, open = "rb")
    }
    on.exit(close(connection), add = TRUE)
    Matrix::readMM(connection)
}

prepare_single_cell_sample <- function(dataset_id, sample_accession, matrix_path, barcode_path, feature_path, output_root) {
    log_event("Single-cell pre-QC:", dataset_id, sample_accession)
    matrix_value <- read_mtx_maybe_gz(matrix_path)
    feature_table <- read_text_table_maybe_gz(feature_path, header = FALSE)
    barcode_table <- read_text_table_maybe_gz(barcode_path, header = FALSE)
    feature_id <- as.character(feature_table[[1]])
    gene_symbol <- if (ncol(feature_table) >= 2L) {
        as.character(feature_table[[2]])
    }
    else {
        feature_id
    }
    gene_symbol[is.na(gene_symbol) | !nzchar(trimws(gene_symbol))] <- feature_id[is.na(gene_symbol) | !nzchar(trimws(gene_symbol))]
    barcodes <- as.character(barcode_table[[1]])
    if (nrow(matrix_value) == length(feature_id) && ncol(matrix_value) == length(barcodes)) {
        orientation <- "FEATURES_BY_CELLS"
    }
    else if (ncol(matrix_value) == length(feature_id) && nrow(matrix_value) == length(barcodes)) {
        matrix_value <- Matrix::t(matrix_value)
        orientation <- "TRANSPOSED_TO_FEATURES_BY_CELLS"
    }
    else {
        stop(dataset_id, " ", sample_accession, " matrix dimensions do not match features and barcodes.")
    }
    matrix_value <- methods::as(matrix_value, "dgCMatrix")
    unique_feature_names <- make.unique(gene_symbol, sep = "__")
    rownames(matrix_value) <- unique_feature_names
    prefixed_barcodes <- paste0(sample_accession, "__", barcodes)
    colnames(matrix_value) <- prefixed_barcodes
    total_counts <- Matrix::colSums(matrix_value)
    detected_features <- diff(matrix_value@p)
    mito_index <- grepl("^MT-", gene_symbol, ignore.case = TRUE)
    mito_counts <- if (any(mito_index)) {
        Matrix::colSums(matrix_value[mito_index, , drop = FALSE])
    }
    else {
        rep(0, ncol(matrix_value))
    }
    percent_mito <- ifelse(total_counts > 0, 100 * mito_counts/total_counts, NA_real_)
    upper_feature_limit <- min(7500, max(2500, median(detected_features, na.rm = TRUE) + 4 * stats::mad(detected_features, 
        na.rm = TRUE)))
    keep_cell <- total_counts >= 500 & detected_features >= 200 & detected_features <= upper_feature_limit & percent_mito <= 
        25
    filtered_matrix <- matrix_value[, keep_cell, drop = FALSE]
    cell_metadata <- data.frame(cell_barcode = prefixed_barcodes, dataset_id = dataset_id, sample_accession = sample_accession, 
        total_counts = as.numeric(total_counts), detected_features = as.numeric(detected_features), percent_mito = as.numeric(percent_mito), 
        preQC_keep = keep_cell, stringsAsFactors = FALSE)
    sample_object <- list(counts = filtered_matrix, feature_table = data.frame(feature_id = feature_id, gene_symbol = gene_symbol, 
        matrix_row_name = unique_feature_names, stringsAsFactors = FALSE), cell_metadata = cell_metadata[keep_cell, 
        , drop = FALSE], sample_metadata = data.frame(dataset_id = dataset_id, sample_accession = sample_accession, 
        stringsAsFactors = FALSE), preQC_parameters = list(minimum_counts = 500, minimum_features = 200, maximum_features = upper_feature_limit, 
        maximum_percent_mito = 25))
    dataset_output_root <- file.path(output_root, dataset_id)
    dir.create(dataset_output_root, recursive = TRUE, showWarnings = FALSE)
    output_path <- file.path(dataset_output_root, paste0(sample_accession, "_preQC.rds"))
    saveRDS(sample_object, output_path, compress = TRUE)
    data.frame(dataset_id = dataset_id, sample_accession = sample_accession, orientation = orientation, feature_n = nrow(matrix_value), 
        raw_cell_n = ncol(matrix_value), retained_cell_n = ncol(filtered_matrix), retained_fraction = ncol(filtered_matrix)/ncol(matrix_value), 
        median_counts = median(total_counts, na.rm = TRUE), median_features = median(detected_features, na.rm = TRUE), 
        median_percent_mito = median(percent_mito, na.rm = TRUE), adaptive_max_features = upper_feature_limit, 
        output_path = normalize_slashes(output_path), stringsAsFactors = FALSE)
}

log_event("Starting", run_id)

planned_path_examples <- c(file.path(bulk_result_root, "loso_gene_results.csv"), file.path(bulk_result_root, "cross_context_candidates.csv"), 
    file.path(sc_object_root, "GSE250130", "GSM0000000_preQC.rds"))

path_length_audit <- data.frame(path = normalize_slashes(planned_path_examples), path_length = nchar(planned_path_examples), 
    under_soft_limit = nchar(planned_path_examples) < 220L, stringsAsFactors = FALSE)

write_csv_utf8(path_length_audit, file.path(dirs$qc, "S4D_path_length_audit.csv"))

if (!all(path_length_audit$under_soft_limit)) {
    stop("At least one planned STEP04D path exceeds the 220-character soft limit.")
}

step04C_gate_path <- file.path(dirs$results, "STEP04C_GATE_TABLE.csv")

step04C_locked_root <- file.path(dirs$results, "STEP04C_locked_sample_designs", "COMBINED_SAMPLE_DESIGN_v1.0")

step04C_manifest_path <- file.path(step04C_locked_root, "COMBINED_SAMPLE_DESIGN_v1.0_LOCKED_MANIFEST.csv")

step04B_locked_root <- file.path(dirs$results, "STEP04B_locked_combined_inputs", "COMBINED_INPUT_SOURCE_v1.0")

source_reference_path <- file.path(step04B_locked_root, "COMBINED_INPUT_SOURCE_v1.0_REFERENCE_LOCK.csv")

GSE250130_bundle_path <- file.path(step04B_locked_root, "COMBINED_INPUT_SOURCE_v1.0_GSE250130_BUNDLE_MAP.csv")

GSE183837_bundle_path <- file.path(step04B_locked_root, "COMBINED_INPUT_SOURCE_v1.0_GSE183837_BUNDLE_MAP.csv")

platform_lock_path <- file.path(dirs$results, "STEP04C_BULK_PLATFORM_ANNOTATION_SOURCE_LOCK.csv")

candidate_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "RNA_OCCUPANCY_EVIDENCE_v1.0_integrated_436.csv")

routing_path <- file.path(project_root, "05_analysis_steps", "STEP03_VDR_perturbation_regulons", "06_results", 
    "STEP03B3C3B_locked_RNA_occupancy_evidence", "RNA_OCCUPANCY_EVIDENCE_v1.0", "DOWNSTREAM_CANDIDATE_ROUTING_v1.0.csv")

required_upstream <- c(step04C_gate_path, step04C_manifest_path, source_reference_path, GSE250130_bundle_path, 
    GSE183837_bundle_path, platform_lock_path, candidate_path, routing_path)

missing_upstream <- required_upstream[!file.exists(required_upstream)]

if (length(missing_upstream) > 0L) {
    stop("Required locked upstream artifact is missing:\n", paste(missing_upstream, collapse = "\n"))
}

step04C_gate <- .vdr_read_csv(step04C_gate_path, stringsAsFactors = FALSE, check.names = FALSE)

accepted_gate <- step04C_gate$status[step04C_gate$gate_component == "STEP04C gate"]

if (length(accepted_gate) != 1L || accepted_gate != "READY_FOR_COMBINED_BULK_AND_SINGLE_CELL_ANALYSIS_BUILD") {
    stop("STEP04C does not authorize the combined analysis build.")
}

step04C_manifest <- .vdr_read_csv(step04C_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

step04C_manifest$file_exists <- file.exists(step04C_manifest$path)

step04C_manifest$observed_md5 <- vapply(step04C_manifest$path, md5_file, character(1))

step04C_manifest$observed_sha256 <- vapply(step04C_manifest$path, sha256_file, character(1))

step04C_manifest$identity_pass <- with(step04C_manifest, file_exists & md5 == observed_md5 & sha256 == observed_sha256)

if (!all(step04C_manifest$identity_pass)) {
    stop("At least one STEP04C locked artifact failed hash validation.")
}

read_locked_design <- function(dataset_id) {
    .vdr_read_csv(file.path(step04C_locked_root, paste0("COMBINED_SAMPLE_DESIGN_v1.0_", dataset_id, ".csv")), stringsAsFactors = FALSE, 
        check.names = FALSE)
}

design_250130 <- read_locked_design("GSE250130")

design_165004 <- read_locked_design("GSE165004")

design_92324 <- read_locked_design("GSE92324")

design_111974 <- read_locked_design("GSE111974")

design_183837 <- read_locked_design("GSE183837")

source_reference <- .vdr_read_csv(source_reference_path, stringsAsFactors = FALSE, check.names = FALSE)

platform_lock <- .vdr_read_csv(platform_lock_path, stringsAsFactors = FALSE, check.names = FALSE)

candidate_table <- .vdr_read_csv(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)

routing_table <- .vdr_read_csv(routing_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(candidate_table) != 436L || length(unique(candidate_table$feature_id)) != 436L || sum(candidate_table$occupancy_supported) != 
    284L) {
    stop("The locked 436/284 RNA occupancy candidate table failed.")
}

get_series_path <- function(dataset_id) {
    rows <- source_reference[source_reference$dataset_id == dataset_id & grepl("SERIES_MATRIX", source_reference$source_role), 
        , drop = FALSE]
    if (nrow(rows) != 1L) {
        stop("Expected one series matrix for ", dataset_id)
    }
    rows$path[[1]]
}

get_platform_path <- function(platform_id) {
    rows <- platform_lock[platform_lock$platform_id == platform_id, , drop = FALSE]
    if (nrow(rows) != 1L) {
        stop("Expected one platform source for ", platform_id)
    }
    rows$path[[1]]
}

download_with_resume(url = GSE92324_profile_url, destination = GSE92324_profile_path)

GSE92324_profile_source_audit <- data.frame(dataset_id = "GSE92324", source_role = "AUTHOR_SUPPLEMENTARY_NORMALIZED_PROFILE_WITH_GENE_ANNOTATION", 
    source_url = GSE92324_profile_url, path = normalize_slashes(GSE92324_profile_path), size_bytes = as.numeric(file.info(GSE92324_profile_path)$size), 
    md5 = md5_file(GSE92324_profile_path), sha256 = sha256_file(GSE92324_profile_path), source_status = "CANDIDATE_FOR_PROBE_GENE_MAPPING", 
    stringsAsFactors = FALSE)

write_csv_utf8(GSE92324_profile_source_audit, file.path(dirs$results, "STEP04D_v3_GSE92324_SUPPLEMENTARY_PROFILE_SOURCE.csv"))

bulk_specs <- list(GSE165004 = list(series_path = get_series_path("GSE165004"), platform_path = get_platform_path("GPL16699"), 
    design = design_165004, symbol_candidates = c("GENE_SYMBOL", "Gene Symbol", "Symbol"), supplementary_annotation_path = NA_character_), 
    GSE92324 = list(series_path = get_series_path("GSE92324"), platform_path = get_platform_path("GPL10558"), design = design_92324, 
        symbol_candidates = c("Symbol", "ILMN_Gene", "GENE_SYMBOL"), supplementary_annotation_path = GSE92324_profile_path), 
    GSE111974 = list(series_path = get_series_path("GSE111974"), platform_path = get_platform_path("GPL17077"), 
        design = design_111974, symbol_candidates = c("GENE_SYMBOL", "Gene Symbol", "Symbol"), supplementary_annotation_path = NA_character_))

bulk_objects <- list()

bulk_qc_rows <- list()

previous_bulk_qc_path <- file.path(dirs$results, "STEP04D_BULK_GENE_MATRIX_QC.csv")

previous_matrix_paths <- file.path(previous_bulk_matrix_root, paste0("BULK_GENE_MATRIX_v1.1_", names(bulk_specs), 
    ".rds"))

reuse_previous_bulk_matrices <- all(file.exists(previous_matrix_paths)) && file.exists(previous_bulk_qc_path)

if (reuse_previous_bulk_matrices) {
    log_event("Reusing three complete v1.1 bulk matrices from the interrupted v5 run.")
    previous_bulk_qc <- .vdr_read_csv(previous_bulk_qc_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!setequal(previous_bulk_qc$dataset_id, names(bulk_specs))) {
        stop("The previous bulk-matrix QC table does not contain all three datasets.")
    }
    for (dataset_id in names(bulk_specs)) {
        previous_path <- file.path(previous_bulk_matrix_root, paste0("BULK_GENE_MATRIX_v1.1_", dataset_id, ".rds"))
        bulk_objects[[dataset_id]] <- readRDS(previous_path)
        if (!is.matrix(bulk_objects[[dataset_id]]) || nrow(bulk_objects[[dataset_id]]) < 5000L) {
            stop("Reusable bulk matrix failed validation: ", dataset_id)
        }
        new_matrix_path <- file.path(bulk_matrix_root, paste0("BGM_v1_2_", dataset_id, ".rds"))
        saveRDS(bulk_objects[[dataset_id]], new_matrix_path, compress = TRUE)
        bulk_qc_rows[[dataset_id]] <- previous_bulk_qc[previous_bulk_qc$dataset_id == dataset_id, , drop = FALSE]
        bulk_qc_rows[[dataset_id]]$matrix_source <- "REUSED_FROM_INTERRUPTED_V5_BULK_GENE_MATRIX_v1.1"
    }
} else {
    log_event("Reusable v1.1 bulk matrices were not all available; rebuilding all three matrices.")
    for (dataset_id in names(bulk_specs)) {
        spec <- bulk_specs[[dataset_id]]
        built <- build_gene_matrix(dataset_id = dataset_id, series_path = spec$series_path, platform_path = spec$platform_path, 
            sample_design = spec$design, symbol_candidates = spec$symbol_candidates, supplementary_annotation_path = spec$supplementary_annotation_path)
        bulk_objects[[dataset_id]] <- built$gene_expression
        built$qc$matrix_source <- "REBUILT_IN_STEP04D_V6"
        bulk_qc_rows[[dataset_id]] <- built$qc
        saveRDS(built$gene_expression, file.path(bulk_matrix_root, paste0("BGM_v1_2_", dataset_id, ".rds")), compress = TRUE)
        write_csv_utf8(built$probe_map, file.path(dirs$qc, paste0("S4D_", dataset_id, "_probe_gene_map.csv")))
    }
}

bulk_matrix_qc <- do.call(rbind, bulk_qc_rows)

write_csv_utf8(bulk_matrix_qc, file.path(dirs$results, "S4D_bulk_matrix_qc.csv"))

program_definitions <- list(RNA_ALL_436 = candidate_table$feature_id, OCCUPANCY_SUPPORTED_284 = candidate_table$feature_id[candidate_table$occupancy_supported], 
    SHARED_OCCUPANCY_46 = candidate_table$feature_id[candidate_table$shared_occupancy_supported], PROMOTER_OCCUPANCY_84 = candidate_table$feature_id[candidate_table$any_promoter_overlap_support])

score_tables <- list()

coverage_tables <- list()

for (dataset_id in names(bulk_objects)) {
    score_value <- calculate_program_scores(gene_expression = bulk_objects[[dataset_id]], candidate_table = candidate_table, 
        program_definitions = program_definitions, dataset_id = dataset_id)
    score_tables[[dataset_id]] <- score_value$scores
    coverage_tables[[dataset_id]] <- score_value$coverage
}

all_scores <- do.call(rbind, score_tables)

all_coverage <- do.call(rbind, coverage_tables)

write_csv_utf8(all_scores, file.path(bulk_result_root, "sample_scores.csv"))

write_csv_utf8(all_coverage, file.path(bulk_result_root, "gene_coverage.csv"))

group_165004 <- factor(design_165004$group, levels = c("CONTROL", "RPL", "UIF"))

design_matrix_165004 <- model.matrix(~0 + group_165004)

colnames(design_matrix_165004) <- levels(group_165004)

contrast_matrix_165004 <- limma::makeContrasts(RPL_vs_CONTROL = RPL - CONTROL, UIF_vs_CONTROL = UIF - CONTROL, 
    levels = design_matrix_165004)

gene_results_165004 <- fit_limma_contrasts(bulk_objects[["GSE165004"]], design_matrix_165004, contrast_matrix_165004, 
    "GSE165004")

score_results_165004 <- fit_score_model(score_tables[["GSE165004"]], design_165004, design_matrix_165004, contrast_matrix_165004, 
    "GSE165004")

group_92324 <- factor(design_92324$group, levels = c("CONTROL", "IMPLANTATION_FAILURE"))

day_92324 <- factor(design_92324$collection_day)

design_matrix_92324 <- model.matrix(~0 + group_92324 + day_92324)

colnames(design_matrix_92324) <- make.names(colnames(design_matrix_92324))

contrast_matrix_92324 <- limma::makeContrasts(IF_vs_CONTROL = group_92324IMPLANTATION_FAILURE - group_92324CONTROL, 
    levels = design_matrix_92324)

gene_results_92324 <- fit_limma_contrasts(bulk_objects[["GSE92324"]], design_matrix_92324, contrast_matrix_92324, 
    "GSE92324")

score_results_92324 <- fit_score_model(score_tables[["GSE92324"]], design_92324, design_matrix_92324, contrast_matrix_92324, 
    "GSE92324")

group_111974 <- factor(design_111974$group, levels = c("CONTROL", "RIF"))

design_matrix_111974 <- model.matrix(~0 + group_111974)

colnames(design_matrix_111974) <- levels(group_111974)

contrast_matrix_111974 <- limma::makeContrasts(RIF_vs_CONTROL = RIF - CONTROL, levels = design_matrix_111974)

gene_results_111974 <- fit_limma_contrasts(bulk_objects[["GSE111974"]], design_matrix_111974, contrast_matrix_111974, 
    "GSE111974")

score_results_111974 <- fit_score_model(score_tables[["GSE111974"]], design_111974, design_matrix_111974, contrast_matrix_111974, 
    "GSE111974")

all_gene_results <- rbind(gene_results_165004, gene_results_92324, gene_results_111974)

all_score_results <- rbind(score_results_165004, score_results_92324, score_results_111974)

write_csv_utf8(all_gene_results, file.path(bulk_result_root, "full_gene_results.csv"))

write_csv_utf8(all_score_results, file.path(bulk_result_root, "score_models.csv"))

candidate_gene_results <- merge(candidate_table, all_gene_results, by.x = "feature_id", by.y = "gene_symbol", all.x = TRUE, 
    sort = FALSE)

write_csv_utf8(candidate_gene_results, file.path(bulk_result_root, "candidate_gene_results.csv"))

full_165004_results <- gene_results_165004[gene_results_165004$gene_symbol %in% candidate_table$feature_id, c("gene_symbol", 
    "contrast", "logFC"), drop = FALSE]

slide_ids <- unique(design_165004$slide_id)

slide_gene_rows <- list()

slide_score_rows <- list()

for (slide_id in slide_ids) {
    keep_samples <- design_165004$slide_id != slide_id
    subset_design <- design_165004[keep_samples, , drop = FALSE]
    subset_group <- factor(subset_design$group, levels = c("CONTROL", "RPL", "UIF"))
    subset_design_matrix <- model.matrix(~0 + subset_group)
    colnames(subset_design_matrix) <- levels(subset_group)
    subset_contrast_matrix <- limma::makeContrasts(RPL_vs_CONTROL = RPL - CONTROL, UIF_vs_CONTROL = UIF - CONTROL, 
        levels = subset_design_matrix)
    subset_gene_result <- fit_limma_contrasts(bulk_objects[["GSE165004"]][, keep_samples, drop = FALSE], subset_design_matrix, 
        subset_contrast_matrix, "GSE165004")
    subset_gene_result <- subset_gene_result[subset_gene_result$gene_symbol %in% candidate_table$feature_id, c("gene_symbol", 
        "contrast", "logFC", "P.Value", "adj.P.Val"), drop = FALSE]
    subset_gene_result$left_out_slide <- slide_id
    slide_gene_rows[[slide_id]] <- subset_gene_result
    subset_scores <- score_tables[["GSE165004"]]
    subset_scores <- subset_scores[subset_scores$sample_accession %in% subset_design$sample_accession, , drop = FALSE]
    subset_score_result <- fit_score_model(subset_scores, subset_design, subset_design_matrix, subset_contrast_matrix, 
        "GSE165004")
    subset_score_result$left_out_slide <- slide_id
    slide_score_rows[[slide_id]] <- subset_score_result
}

slide_gene_results <- do.call(rbind, slide_gene_rows)

slide_score_results <- do.call(rbind, slide_score_rows)

write_csv_utf8(slide_gene_results, file.path(bulk_result_root, "loso_gene_results.csv"))

write_csv_utf8(slide_score_results, file.path(bulk_result_root, "loso_score_results.csv"))

slide_robustness_rows <- list()

for (contrast_name in unique(full_165004_results$contrast)) {
    full_subset <- full_165004_results[full_165004_results$contrast == contrast_name, , drop = FALSE]
    for (slide_id in slide_ids) {
        left_subset <- slide_gene_results[slide_gene_results$contrast == contrast_name & slide_gene_results$left_out_slide == 
            slide_id, , drop = FALSE]
        joined <- merge(full_subset, left_subset, by = c("gene_symbol", "contrast"), suffixes = c("_full", "_leaveout"), 
            all = FALSE)
        slide_robustness_rows[[paste(contrast_name, slide_id, sep = "__")]] <- data.frame(contrast = contrast_name, 
            left_out_slide = slide_id, candidate_gene_n = nrow(joined), spearman_rho = suppressWarnings(cor(joined$logFC_full, 
                joined$logFC_leaveout, method = "spearman", use = "pairwise.complete.obs")), direction_concordance = mean(sign(joined$logFC_full) == 
                sign(joined$logFC_leaveout), na.rm = TRUE), stringsAsFactors = FALSE)
    }
}

slide_robustness <- do.call(rbind, slide_robustness_rows)

write_csv_utf8(slide_robustness, file.path(bulk_result_root, "loso_robustness.csv"))

compact_candidate_results <- candidate_gene_results[, c("feature_id", "program_class", "interaction_logFC", "interaction_FDR", 
    "occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support", "dataset_id", "contrast", 
    "logFC", "P.Value", "adj.P.Val"), drop = FALSE]

compact_candidate_results$context_id <- paste(compact_candidate_results$dataset_id, compact_candidate_results$contrast, 
    sep = "__")

context_ids <- unique(compact_candidate_results$context_id[!is.na(compact_candidate_results$context_id)])

candidate_convergence <- unique(compact_candidate_results[, c("feature_id", "program_class", "interaction_logFC", 
    "interaction_FDR", "occupancy_supported", "shared_occupancy_supported", "any_promoter_overlap_support"), drop = FALSE])

for (context_id in context_ids) {
    context_rows <- compact_candidate_results[compact_candidate_results$context_id == context_id, c("feature_id", 
        "logFC", "P.Value", "adj.P.Val"), drop = FALSE]
    names(context_rows)[-1] <- paste0(context_id, "__", c("logFC", "P.Value", "adj.P.Val"))
    candidate_convergence <- merge(candidate_convergence, context_rows, by = "feature_id", all.x = TRUE, sort = FALSE)
}

logFC_columns <- grep("__logFC$", names(candidate_convergence), value = TRUE)

candidate_convergence$measured_context_n <- rowSums(is.finite(as.matrix(candidate_convergence[, logFC_columns, 
    drop = FALSE])))

candidate_convergence$positive_context_n <- rowSums(as.matrix(candidate_convergence[, logFC_columns, drop = FALSE]) > 
    0, na.rm = TRUE)

candidate_convergence$negative_context_n <- rowSums(as.matrix(candidate_convergence[, logFC_columns, drop = FALSE]) < 
    0, na.rm = TRUE)

candidate_convergence$cross_context_direction_label <- ifelse(candidate_convergence$positive_context_n == candidate_convergence$measured_context_n, 
    "POSITIVE_ALL_MEASURED_CONTEXTS", ifelse(candidate_convergence$negative_context_n == candidate_convergence$measured_context_n, 
        "NEGATIVE_ALL_MEASURED_CONTEXTS", "CONTEXT_DEPENDENT_OR_MIXED_DIRECTION"))

write_csv_utf8(candidate_convergence, file.path(bulk_result_root, "cross_context_candidates.csv"))

GSE250130_bundle_map <- .vdr_read_csv(GSE250130_bundle_path, stringsAsFactors = FALSE, check.names = FALSE)

GSE183837_bundle_map <- .vdr_read_csv(GSE183837_bundle_path, stringsAsFactors = FALSE, check.names = FALSE)

source_183837 <- source_reference[source_reference$dataset_id == "GSE183837" & source_reference$source_role == 
    "SAMPLE_LEVEL_MATRIX_MARKET_TRIPLET_MEMBER", , drop = FALSE]

source_183837$sample_accession <- vapply(source_183837$path, function(path_value) {
    match_value <- regmatches(basename(path_value), regexpr("GSM[0-9]+", basename(path_value), ignore.case = TRUE, 
        perl = TRUE))
    toupper(match_value)
}, character(1))

source_183837$component <- vapply(source_183837$file_name, function(file_name) {
    file_lower <- tolower(file_name)
    if (grepl("matrix.*\\.mtx(\\.gz)?$", file_lower, perl = TRUE)) {
        "matrix"
    }
    else if (grepl("barcodes.*\\.(tsv|txt)(\\.gz)?$", file_lower, perl = TRUE)) {
        "barcodes"
    }
    else if (grepl("(genes|features).*\\.(tsv|txt)(\\.gz)?$", file_lower, perl = TRUE)) {
        "features"
    }
    else {
        "other"
    }
}, character(1))

sc_qc_rows <- list()

sc_qc_index <- 1L

for (row_index in seq_len(nrow(GSE250130_bundle_map))) {
    map_row <- GSE250130_bundle_map[row_index, , drop = FALSE]
    sc_qc_rows[[sc_qc_index]] <- prepare_single_cell_sample(dataset_id = "GSE250130", sample_accession = map_row$sample_accession, 
        matrix_path = map_row$matrix_file, barcode_path = map_row$barcode_file, feature_path = map_row$feature_file, 
        output_root = sc_object_root)
    sc_qc_index <- sc_qc_index + 1L
    gc()
}

for (sample_accession in unique(source_183837$sample_accession)) {
    sample_rows <- source_183837[source_183837$sample_accession == sample_accession, , drop = FALSE]
    matrix_path <- sample_rows$path[sample_rows$component == "matrix"]
    barcode_path <- sample_rows$path[sample_rows$component == "barcodes"]
    feature_path <- sample_rows$path[sample_rows$component == "features"]
    if (length(matrix_path) != 1L || length(barcode_path) != 1L || length(feature_path) != 1L) {
        stop("GSE183837 triplet path reconstruction failed for ", sample_accession)
    }
    sc_qc_rows[[sc_qc_index]] <- prepare_single_cell_sample(dataset_id = "GSE183837", sample_accession = sample_accession, 
        matrix_path = matrix_path, barcode_path = barcode_path, feature_path = feature_path, output_root = sc_object_root)
    sc_qc_index <- sc_qc_index + 1L
    gc()
}

sc_qc_summary <- do.call(rbind, sc_qc_rows)

write_csv_utf8(sc_qc_summary, file.path(dirs$results, "S4D_sc_sample_qc.csv"))

sc_dataset_summary <- aggregate(cbind(raw_cell_n, retained_cell_n) ~ dataset_id, data = sc_qc_summary, FUN = sum)

sc_dataset_summary$sample_n <- vapply(sc_dataset_summary$dataset_id, function(dataset_id) {
    sum(sc_qc_summary$dataset_id == dataset_id)
}, integer(1))

sc_dataset_summary$retained_fraction <- sc_dataset_summary$retained_cell_n/sc_dataset_summary$raw_cell_n

write_csv_utf8(sc_dataset_summary, file.path(dirs$results, "S4D_sc_dataset_qc.csv"))

minimum_candidate_coverage <- min(all_coverage$measurable_gene_n[all_coverage$program_name == "OCCUPANCY_SUPPORTED_284"])

minimum_slide_spearman <- min(slide_robustness$spearman_rho, na.rm = TRUE)

minimum_slide_direction <- min(slide_robustness$direction_concordance, na.rm = TRUE)

expected_sc_samples <- c(GSE250130 = 28L, GSE183837 = 9L)

observed_sc_samples <- table(sc_qc_summary$dataset_id)

sc_sample_pass <- all(as.integer(observed_sc_samples[names(expected_sc_samples)]) == as.integer(expected_sc_samples))

sc_retention_pass <- all(sc_qc_summary$retained_cell_n > 0L)

bulk_matrix_pass <- all(bulk_matrix_qc$gene_n >= 5000L)

candidate_coverage_pass <- minimum_candidate_coverage >= 200L

gate_status <- if (bulk_matrix_pass && candidate_coverage_pass && sc_sample_pass && sc_retention_pass) {
    "READY_FOR_SINGLE_CELL_INTEGRATION_AND_FINAL_MULTILAYER_CANDIDATE_PRIORITIZATION"
} else {
    "HOLD_STEP04D_QC_REVIEW_REQUIRED"
}

gate_table <- data.frame(gate_component = c("STEP04C gate", "Bulk gene matrices", "Occupancy-supported candidate coverage", 
    "Bulk gene-level contrasts", "Program-score models", "GSE165004 leave-one-slide candidate robustness", "GSE250130 donor/sample pre-QC", 
    "GSE183837 donor/sample pre-QC", "Single-cell retained cells", "GSE111974 independence interpretation", "Definitive direct-target assignment", 
    "Scientific figures", "STEP04D gate"), status = c("PASS", ifelse(bulk_matrix_pass, "PASS_3_OF_3", "FAIL"), 
    paste0("MIN_MEASURABLE_", minimum_candidate_coverage, "_OF_284"), "COMPLETED_4_CONTRASTS", "COMPLETED_4_PROGRAMS_ACROSS_4_CONTRASTS", 
    paste0("MIN_SPEARMAN_", format(minimum_slide_spearman, digits = 6), "_MIN_DIRECTION_", format(minimum_slide_direction, 
        digits = 6)), paste0("PASS_", sum(sc_qc_summary$dataset_id == "GSE250130"), "_OF_28"), paste0("PASS_", 
        sum(sc_qc_summary$dataset_id == "GSE183837"), "_OF_9"), ifelse(sc_retention_pass, "PASS_ALL_SAMPLES_RETAIN_CELLS", 
        "FAIL"), "SHARED_CONTROL_CONTEXT_EXTENSION_ONLY", "NOT_PERFORMED", "NONE_CREATED", gate_status), required_action = c("Retain locked designs and source references.", 
    "Use only locked gene-level matrices.", "Review coverage before interpreting absent genes.", "Interpret cross-context direction as supportive or mixed, not causal.", 
    "Do not interpret scores as diagnostic classifiers.", "Retain slide sensitivity in the primary patient result.", 
    "Use donor-aware integration and pseudobulk.", "Use donor-aware integration and pseudobulk.", "Review sample-level QC before clustering.", 
    "Do not count as independent external validation.", "Use occupancy-supported candidate terminology.", "Render figures only after single-cell integration and candidate prioritization.", 
    "Proceed only if the Gate is READY."), stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP04D_GATE_TABLE.csv"))

write_text(c("# STEP04D bulk projection and single-cell pre-QC card", "", paste0("- Run ID: ", run_id), paste0("- Bulk matrix version: ", 
    bulk_matrix_version), paste0("- Bulk projection version: ", bulk_projection_version), paste0("- Single-cell pre-QC version: ", 
    sc_preQC_version), paste0("- Bulk genes per dataset: ", paste(paste0(bulk_matrix_qc$dataset_id, "=", bulk_matrix_qc$gene_n), 
    collapse = "; ")), paste0("- Minimum measurable occupancy-supported candidates: ", minimum_candidate_coverage, 
    "/284"), paste0("- Minimum GSE165004 leave-one-slide Spearman: ", format(minimum_slide_spearman, digits = 6)), 
    paste0("- Minimum leave-one-slide direction concordance: ", format(minimum_slide_direction, digits = 6)), paste0("- GSE250130 raw/retained cells: ", 
        sc_dataset_summary$raw_cell_n[sc_dataset_summary$dataset_id == "GSE250130"], "/", sc_dataset_summary$retained_cell_n[sc_dataset_summary$dataset_id == 
            "GSE250130"]), paste0("- GSE183837 raw/retained cells: ", sc_dataset_summary$raw_cell_n[sc_dataset_summary$dataset_id == 
        "GSE183837"], "/", sc_dataset_summary$retained_cell_n[sc_dataset_summary$dataset_id == "GSE183837"]), "- GSE111974 independent external validation: NO", 
    "- Definitive direct-target assignments: 0", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), 
    file.path(dirs$results, "S4D_result_card.md"))

write_text(c("# STEP04D Methods draft", "", "Processed bulk-expression matrices were mapped to gene symbols using locked GEO platform tables. Only probes with one unambiguous gene symbol were retained, and duplicate probes were averaged with limma. A log2 transformation was applied only when matrix-scale diagnostics indicated that the submitted values were not already log transformed.", 
    "", "The prespecified 436-gene RNA interaction program and its 284-gene occupancy-supported subset were projected into the three bulk cohorts using signed standardized expression scores based on the direction of the locked RNA interaction effect. GSE165004 was analyzed separately for RPL versus control and UIF versus control, with leave-one-slide sensitivity. GSE92324 was adjusted for HCG collection day. GSE111974 was treated as a shared-control RIF context extension.", 
    "", "For GSE250130 and GSE183837, each donor/sample matrix was read independently. Cells with fewer than 500 counts, fewer than 200 detected features, excessive detected features under a sample-adaptive upper threshold, or more than 25% mitochondrial counts were removed. Filtered sparse matrices were retained as donor-level preprocessing objects; cells were not treated as independent biological replicates."), 
    file.path(dirs$text, "STEP04D_Methods_draft.md"))

write_text(c("# STEP04D Results shell", "", "This file is intentionally a result shell. Numeric bulk projection, leave-one-slide and single-cell preprocessing results must be inserted only after the STEP04D Gate and QC tables are reviewed.", 
    "", "No diagnostic, causal or definitive direct-target claim is authorized."), file.path(dirs$text, "STEP04D_Results_shell.md"))

parameters <- list(schema_version = "1.0", run_id = run_id, bulk_matrix_version = bulk_matrix_version, bulk_projection_version = bulk_projection_version, 
    sc_preQC_version = sc_preQC_version, frozen_candidate_n = 436, occupancy_supported_n = 284, minimum_measurable_occupancy_supported_n = minimum_candidate_coverage, 
    minimum_leave_one_slide_spearman = minimum_slide_spearman, minimum_leave_one_slide_direction_concordance = minimum_slide_direction, 
    GSE250130_sample_n = sum(sc_qc_summary$dataset_id == "GSE250130"), GSE183837_sample_n = sum(sc_qc_summary$dataset_id == 
        "GSE183837"), definitive_direct_target_assignment = FALSE, scientific_figures_created = FALSE, gate = gate_status)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP04D_parameters.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP04D.txt"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(error_object) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    copy_if_different(source_file, file.path(dirs$code, "fit_patient_bulk_models.R"))
}

matrix_paths <- list.files(bulk_matrix_root, pattern = "\\.rds$", full.names = TRUE)

result_paths <- c(file.path(dirs$results, "STEP04D_GATE_TABLE.csv"), file.path(dirs$results, "S4D_result_card.md"), 
    file.path(dirs$results, "S4D_bulk_matrix_qc.csv"), file.path(dirs$results, "STEP04D_v3_GSE92324_SUPPLEMENTARY_PROFILE_SOURCE.csv"), 
    file.path(dirs$results, "S4D_sc_sample_qc.csv"), file.path(dirs$results, "S4D_sc_dataset_qc.csv"), file.path(bulk_result_root, 
        "score_models.csv"), file.path(bulk_result_root, "gene_coverage.csv"), file.path(bulk_result_root, "loso_robustness.csv"), 
    file.path(bulk_result_root, "cross_context_candidates.csv"), file.path(dirs$text, "STEP04D_Methods_draft.md"), 
    file.path(dirs$text, "STEP04D_Results_shell.md"), file.path(dirs$parameters, "STEP04D_parameters.json"), file.path(dirs$environment, 
        "sessionInfo_POST_RUN_STEP04D.txt"), file.path(dirs$code, "fit_patient_bulk_models.R"), log_path)

artifact_manifest_paths <- c(matrix_paths, result_paths[file.exists(result_paths)])

artifact_manifest <- data.frame(path = normalize_slashes(artifact_manifest_paths), size_bytes = as.numeric(file.info(artifact_manifest_paths)$size), 
    md5 = vapply(artifact_manifest_paths, md5_file, character(1)), sha256 = vapply(artifact_manifest_paths, sha256_file, 
        character(1)), artifact_status = "CANDIDATE_PENDING_MANUAL_QC", stringsAsFactors = FALSE)

write_csv_utf8(artifact_manifest, file.path(dirs$results, "S4D_artifact_manifest.csv"))

archive_artifacts <- c(result_paths, file.path(dirs$results, "S4D_artifact_manifest.csv"))

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
    else if (grepl("Methods|Results", basename(path_value), ignore.case = TRUE)) {
        "txt"
    }
    else if (grepl("QC|CARD|MANIFEST|ROBUSTNESS|COVERAGE|GATE|SUMMARY", basename(path_value), ignore.case = TRUE)) {
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
    sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(archive_root, "ctl", "STEP04D_SHORT_ARCHIVE_MANIFEST.csv"))

write_text(c("# STEP04D archive status", "", paste0("- Run ID: ", run_id), "- Bulk gene matrices: 3", "- Bulk contrasts: 4", 
    "- Program-score definitions: 4", "- Single-cell sample objects: 37", "- Large single-cell RDS objects copied into archive: NO", 
    "- Definitive direct-target assignments: 0", "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), 
    file.path(archive_root, "ctl", "ARCHIVE_STATUS.md"))

write_text(c("# STEP04 combined package status", "", paste0("- Accepted run candidate: ", run_id), paste0("- Bulk matrices: ", 
    bulk_matrix_version, " CANDIDATE"), paste0("- Bulk projection: ", bulk_projection_version, " CANDIDATE"), paste0("- Single-cell pre-QC: ", 
    sc_preQC_version, " CANDIDATE"), "- Scientific figures: NONE", paste0("- Gate: ", gate_status)), file.path(dirs$control, 
    "STATUS.md"))

log_event("Completed", run_id, "Gate:", gate_status)

cat("\n============================================================\n")

cat("STEP04D v6 bulk projection and single-cell pre-QC completed.\n")

cat("Bulk matrix version: ", bulk_matrix_version, "\n")

cat("Bulk projection version: ", bulk_projection_version, "\n")

cat("Single-cell pre-QC version: ", sc_preQC_version, "\n")

cat("Bulk gene counts: ", paste(paste0(bulk_matrix_qc$dataset_id, "=", bulk_matrix_qc$gene_n), collapse = "; "), 
    "\n")

cat("Minimum measurable occupancy-supported genes: ", minimum_candidate_coverage, "/284\n")

cat("Minimum GSE165004 leave-one-slide Spearman: ", format(minimum_slide_spearman, digits = 6), "\n")

cat("Minimum leave-one-slide direction concordance: ", format(minimum_slide_direction, digits = 6), "\n")

cat("GSE250130 samples processed: ", sum(sc_qc_summary$dataset_id == "GSE250130"), "/28\n")

cat("GSE183837 samples processed: ", sum(sc_qc_summary$dataset_id == "GSE183837"), "/9\n")

cat("Definitive direct-target assignments: 0\n")

cat("Scientific figures: NONE\n")

cat("Gate: ", gate_status, "\n")

cat("Authoritative archive:\n  ", archive_root, "\n\n")

cat("Upload:\n")

cat("  06_results/STEP04D_GATE_TABLE.csv\n")

cat("  06_results/S4D_bulk_matrix_qc.csv\n")

cat("  06_results/S4D_sc_sample_qc.csv\n")

cat("  06_results/S4D_sc_dataset_qc.csv\n")

cat("  06_results/S4D_bulk/BVP_v1_2/score_models.csv\n")

cat("  06_results/S4D_bulk/BVP_v1_2/loso_robustness.csv\n")

cat("  06_results/S4D_bulk/BVP_v1_2/cross_context_candidates.csv\n")

cat("============================================================\n")

