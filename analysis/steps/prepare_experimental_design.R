#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A03: Prepare the accepted experimental sample design

options(stringsAsFactors = FALSE)

options(timeout = 7200)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03A_v4_", format(Sys.time(), "%Y%m%d_%H%M%S"))

if (!dir.exists(project_root)) {
    stop("Project root not found:\n", project_root)
}

invisible(NULL)

dirs <- list(control = file.path(module_root, "00_module_control"), inputs = file.path(module_root, "01_inputs"), 
    code = file.path(module_root, "02_code"), parameters = file.path(module_root, "03_parameters"), environment = file.path(module_root, 
        "04_environment"), execution = file.path(module_root, "05_execution"), results = file.path(module_root, 
        "06_results"), figures = file.path(module_root, "07_figures_tables"), qc = file.path(module_root, "08_qc"), 
    text = file.path(module_root, "09_manuscript_text"), archive = file.path(module_root, "10_archive_work", run_id))

for (d in dirs) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

required_pkgs <- c("curl", "digest", "jsonlite")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

if (!all(vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required CRAN packages could not be loaded.")
}

write_text <- function(lines, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(lines, path, useBytes = TRUE)
}

write_csv_utf8 <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
}

md5_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    unname(tools::md5sum(path))
}

sha256_file <- function(path) {
    if (!file.exists(path)) 
        return(NA_character_)
    digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

extract_last_integer <- function(text) {
    hit <- regexec("([0-9]+)[^0-9]*$", text, perl = TRUE)
    parts <- regmatches(text, hit)[[1]]
    if (length(parts) >= 2L) {
        return(as.integer(parts[2]))
    }
    NA_integer_
}

to_https_ncbi <- function(url) {
    url <- trimws(url)
    sub("^ftp://ftp\\.ncbi\\.nlm\\.nih\\.gov", "https://ftp.ncbi.nlm.nih.gov", url)
}

extract_header_number <- function(header_text, field_pattern) {
    lines <- unlist(strsplit(header_text, "\r\n|\n"))
    lines <- trimws(lines)
    hits <- lines[grepl(field_pattern, lines, ignore.case = TRUE, perl = TRUE)]
    if (length(hits) == 0L) {
        return(NA_real_)
    }
    values <- suppressWarnings(as.numeric(sub("^[^:]+:\\s*", "", hits, perl = TRUE)))
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
        NA_real_
    }
    else {
        tail(values, 1L)
    }
}

extract_content_range_total <- function(header_text) {
    lines <- unlist(strsplit(header_text, "\r\n|\n"))
    hits <- lines[grepl("^content-range:", trimws(lines), ignore.case = TRUE)]
    if (length(hits) == 0L) {
        return(NA_real_)
    }
    total <- suppressWarnings(as.numeric(sub("^.*?/([0-9]+).*$", "\\1", tail(hits, 1L), perl = TRUE)))
    ifelse(is.finite(total), total, NA_real_)
}

remote_size_audit <- function(url) {
    user_agent <- paste0("Nutrition_Endometrial_VDR_v2/", run_id)
    head_result <- tryCatch({
        handle <- curl::new_handle(nobody = TRUE, followlocation = TRUE, useragent = user_agent, timeout = 120)
        response <- curl::curl_fetch_memory(url, handle = handle)
        header_text <- rawToChar(response$headers)
        size <- extract_header_number(header_text, "^content-length:")
        list(status_code = response$status_code, size_bytes = size, method = "HEAD", error = "", header_text = header_text)
    }, error = function(e) {
        list(status_code = NA_integer_, size_bytes = NA_real_, method = "HEAD", error = conditionMessage(e), header_text = "")
    })
    if (is.finite(head_result$status_code) && head_result$status_code < 400 && is.finite(head_result$size_bytes) && 
        head_result$size_bytes > 0) {
        return(head_result)
    }
    range_result <- tryCatch({
        handle <- curl::new_handle(followlocation = TRUE, useragent = user_agent, range = "0-0", timeout = 120)
        response <- curl::curl_fetch_memory(url, handle = handle)
        header_text <- rawToChar(response$headers)
        size <- extract_content_range_total(header_text)
        if (!is.finite(size)) {
            size <- extract_header_number(header_text, "^content-length:")
        }
        list(status_code = response$status_code, size_bytes = size, method = "RANGE_0_0", error = "", header_text = header_text)
    }, error = function(e) {
        list(status_code = NA_integer_, size_bytes = NA_real_, method = "RANGE_0_0", error = conditionMessage(e), 
            header_text = "")
    })
    if (nzchar(head_result$error) && !nzchar(range_result$error)) {
        range_result$error <- paste0("HEAD failed: ", head_result$error)
    }
    else if (nzchar(head_result$error) && nzchar(range_result$error)) {
        range_result$error <- paste0("HEAD: ", head_result$error, " | RANGE: ", range_result$error)
    }
    range_result
}

format_gb <- function(bytes) {
    if (!is.finite(bytes)) 
        return(NA_real_)
    bytes/1024^3
}

log_event("Starting", run_id)

v3_manifest_path <- file.path(dirs$results, "STEP03A_v3_sample_design_manifest_PROVISIONAL.csv")

v3_decision_path <- file.path(dirs$results, "STEP03A_v3_processed_data_decision_matrix.csv")

if (!file.exists(v3_manifest_path)) {
    stop("Missing v3 sample manifest:\n", v3_manifest_path)
}

if (!file.exists(v3_decision_path)) {
    stop("Missing v3 processed-data decision matrix:\n", v3_decision_path)
}

samples <- .vdr_read_csv(v3_manifest_path, stringsAsFactors = FALSE, check.names = FALSE)

prior_decision <- .vdr_read_csv(v3_decision_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(samples) != 66L) {
    stop("Expected 66 total samples but found ", nrow(samples), ".")
}

locked <- samples

locked$experiment_block <- NA_character_

locked$genetic_context <- NA_character_

locked$ligand_status_locked <- NA_character_

locked$decidualization_status_locked <- "not_applied_in_sequencing_dataset"

locked$assay_target_locked <- NA_character_

locked$replicate_locked <- NA_integer_

locked$condition_id_locked <- NA_character_

locked$analysis_role_locked <- NA_character_

locked$manual_review_status <- "LOCKED_BY_TITLE_AND_CHARACTERISTICS"

for (i in seq_len(nrow(locked))) {
    acc <- locked$accession[i]
    title <- tolower(locked$title[i])
    characteristics <- tolower(locked$characteristics_ch1[i])
    if (acc == "GSE254251") {
        locked$experiment_block[i] <- "RNA_2x2_siRNA_by_ligand"
        if (grepl("non-targeting pooled sirna", title, fixed = TRUE)) {
            locked$genetic_context[i] <- "siRNA_control"
        }
        else if (grepl("human vdr sirna", title, fixed = TRUE)) {
            locked$genetic_context[i] <- "VDR_knockdown"
        }
        if (grepl("vehicle", title) || grepl("vehicle", characteristics)) {
            locked$ligand_status_locked[i] <- "vehicle"
        }
        else if (grepl("1,25", title, fixed = TRUE) || grepl("1,25", characteristics, fixed = TRUE)) {
            locked$ligand_status_locked[i] <- "calcitriol_2nM"
        }
        locked$assay_target_locked[i] <- "transcriptome"
        locked$replicate_locked[i] <- extract_last_integer(title)
        locked$analysis_role_locked[i] <- "primary_transcriptomic_perturbation_anchor"
    }
    if (acc == "GSE306094") {
        locked$experiment_block[i] <- "ATAC_siRNA_knockdown"
        if (grepl("non-targeting sirna", title, fixed = TRUE) || grepl("genotype: non-target control", characteristics, 
            fixed = TRUE)) {
            locked$genetic_context[i] <- "siRNA_control"
        }
        else if (grepl("vitamin d receptor sirna", title, fixed = TRUE) || grepl("genotype: vdr knockdown", characteristics, 
            fixed = TRUE)) {
            locked$genetic_context[i] <- "VDR_knockdown"
        }
        locked$ligand_status_locked[i] <- "not_applied_in_dataset"
        locked$assay_target_locked[i] <- "open_chromatin"
        locked$replicate_locked[i] <- extract_last_integer(title)
        locked$analysis_role_locked[i] <- "chromatin_accessibility_support"
    }
    if (acc == "GSE306127") {
        if (grepl("transduced with lentiviral", title, fixed = TRUE)) {
            locked$experiment_block[i] <- "CUTRUN_VDR_cistrome"
            if (grepl("lentiviral gfp", title, fixed = TRUE)) {
                locked$genetic_context[i] <- "vector_control"
            }
            else if (grepl("lentiviral vdr", title, fixed = TRUE)) {
                locked$genetic_context[i] <- "VDR_overexpression"
            }
            if (grepl("1,25-dihydroxylated vitamin d3", characteristics, fixed = TRUE) || grepl("treated with ligand", 
                title, fixed = TRUE)) {
                locked$ligand_status_locked[i] <- "calcitriol"
            }
            else if (grepl("ethyl alcohol", characteristics, fixed = TRUE)) {
                locked$ligand_status_locked[i] <- "vehicle"
            }
        }
        else if (grepl("transfected with pooled non-targeting sirna", title, fixed = TRUE) || grepl("transfected with vdr sirna", 
            title, fixed = TRUE)) {
            locked$experiment_block[i] <- "CUTRUN_histone_marks"
            if (grepl("pooled non-targeting sirna", title, fixed = TRUE)) {
                locked$genetic_context[i] <- "siRNA_control"
            }
            else if (grepl("transfected with vdr sirna", title, fixed = TRUE)) {
                locked$genetic_context[i] <- "VDR_knockdown"
            }
            locked$ligand_status_locked[i] <- "not_applied_in_dataset"
        }
        if (grepl("precipitated with rabbit igg", title, fixed = TRUE) || grepl("precipitated with igg", title, 
            fixed = TRUE)) {
            locked$assay_target_locked[i] <- "IgG_control"
        }
        else if (grepl("h3k4me3", title, fixed = TRUE)) {
            locked$assay_target_locked[i] <- "H3K4me3"
        }
        else if (grepl("h3k27ac", title, fixed = TRUE)) {
            locked$assay_target_locked[i] <- "H3K27ac"
        }
        else if (grepl("h3k27me3", title, fixed = TRUE)) {
            locked$assay_target_locked[i] <- "H3K27me3"
        }
        else if (grepl("precipitated with vdr", title, fixed = TRUE)) {
            locked$assay_target_locked[i] <- "VDR"
        }
        locked$replicate_locked[i] <- extract_last_integer(title)
        locked$analysis_role_locked[i] <- ifelse(locked$experiment_block[i] == "CUTRUN_VDR_cistrome", "direct_VDR_occupancy_support", 
            "histone_modification_support")
    }
}

locked$condition_id_locked <- paste(locked$experiment_block, locked$genetic_context, locked$ligand_status_locked, 
    locked$assay_target_locked, sep = "__")

locked$design_lock_complete <- complete.cases(locked[c("experiment_block", "genetic_context", "ligand_status_locked", 
    "decidualization_status_locked", "assay_target_locked", "replicate_locked", "condition_id_locked", "analysis_role_locked")])

if (!all(locked$design_lock_complete)) {
    write_csv_utf8(locked[!locked$design_lock_complete, , drop = FALSE], file.path(dirs$results, "STEP03A_v4_UNRESOLVED_DESIGN_ROWS.csv"))
    stop("One or more samples remain unresolved. ", "See STEP03A_v4_UNRESOLVED_DESIGN_ROWS.csv.")
}

design_counts <- aggregate(list(n_samples = locked$gsm), by = list(accession = locked$accession, assay = locked$assay, 
    experiment_block = locked$experiment_block, genetic_context = locked$genetic_context, ligand_status = locked$ligand_status_locked, 
    target = locked$assay_target_locked, condition_id = locked$condition_id_locked), FUN = length)

expected_checks <- data.frame(check_id = c("RNA_total", "RNA_four_conditions", "RNA_each_condition_n3", "ATAC_total", 
    "ATAC_two_conditions", "ATAC_each_condition_n3", "CUTRUN_total", "CUTRUN_cistrome_total", "CUTRUN_histone_total", 
    "CUTRUN_each_condition_n3", "All_replicates_1_to_3", "No_decidualization_label"), expected = c("12", "4", "TRUE", 
    "6", "2", "TRUE", "48", "24", "24", "TRUE", "TRUE", "TRUE"), observed = c(sum(locked$accession == "GSE254251"), 
    nrow(design_counts[design_counts$accession == "GSE254251", , drop = FALSE]), all(design_counts$n_samples[design_counts$accession == 
        "GSE254251"] == 3L), sum(locked$accession == "GSE306094"), nrow(design_counts[design_counts$accession == 
        "GSE306094", , drop = FALSE]), all(design_counts$n_samples[design_counts$accession == "GSE306094"] == 3L), 
    sum(locked$accession == "GSE306127"), sum(locked$experiment_block == "CUTRUN_VDR_cistrome"), sum(locked$experiment_block == 
        "CUTRUN_histone_marks"), all(design_counts$n_samples[design_counts$accession == "GSE306127"] == 3L), all(vapply(split(locked$replicate_locked, 
        locked$condition_id_locked), function(x) {
        identical(sort(as.integer(x)), 1:3)
    }, logical(1))), all(locked$decidualization_status_locked == "not_applied_in_sequencing_dataset")), pass = c(sum(locked$accession == 
    "GSE254251") == 12L, nrow(design_counts[design_counts$accession == "GSE254251", , drop = FALSE]) == 4L, all(design_counts$n_samples[design_counts$accession == 
    "GSE254251"] == 3L), sum(locked$accession == "GSE306094") == 6L, nrow(design_counts[design_counts$accession == 
    "GSE306094", , drop = FALSE]) == 2L, all(design_counts$n_samples[design_counts$accession == "GSE306094"] == 
    3L), sum(locked$accession == "GSE306127") == 48L, sum(locked$experiment_block == "CUTRUN_VDR_cistrome") == 
    24L, sum(locked$experiment_block == "CUTRUN_histone_marks") == 24L, all(design_counts$n_samples[design_counts$accession == 
    "GSE306127"] == 3L), all(vapply(split(locked$replicate_locked, locked$condition_id_locked), function(x) {
    identical(sort(as.integer(x)), 1:3)
}, logical(1))), all(locked$decidualization_status_locked == "not_applied_in_sequencing_dataset")), stringsAsFactors = FALSE)

if (!all(expected_checks$pass)) {
    write_csv_utf8(expected_checks, file.path(dirs$results, "STEP03A_v4_design_validation_checks.csv"))
    stop("Experimental design validation failed. ", "Review STEP03A_v4_design_validation_checks.csv.")
}

write_csv_utf8(locked, file.path(dirs$results, "STEP03A_v4_sample_design_LOCKED.csv"))

write_csv_utf8(design_counts, file.path(dirs$results, "STEP03A_v4_design_cell_count_LOCKED.csv"))

write_csv_utf8(expected_checks, file.path(dirs$results, "STEP03A_v4_design_validation_checks.csv"))

contrast_plan <- data.frame(contrast_id = c("RNA_C1_ligand_in_siNT", "RNA_C2_ligand_in_siVDR", "RNA_C3_VDRKD_in_vehicle", 
    "RNA_C4_VDRKD_in_ligand", "RNA_C5_siRNA_by_ligand_interaction", "ATAC_C1_VDRKD_vs_siNT", "CUTRUN_VDR_C1_ligand_in_vector", 
    "CUTRUN_VDR_C2_OE_in_vehicle", "CUTRUN_VDR_C3_OE_in_ligand", "CUTRUN_VDR_C4_construct_by_ligand_interaction", 
    "CUTRUN_H3K4me3_C1_VDRKD_vs_siNT", "CUTRUN_H3K27ac_C1_VDRKD_vs_siNT", "CUTRUN_H3K27me3_C1_VDRKD_vs_siNT"), 
    assay = c(rep("RNA-seq", 5), "ATAC-seq", rep("CUT&RUN-seq", 7)), experiment_block = c(rep("RNA_2x2_siRNA_by_ligand", 
        5), "ATAC_siRNA_knockdown", rep("CUTRUN_VDR_cistrome", 4), rep("CUTRUN_histone_marks", 3)), numerator = c("siRNA_control + calcitriol_2nM", 
        "VDR_knockdown + calcitriol_2nM", "VDR_knockdown + vehicle", "VDR_knockdown + calcitriol_2nM", "(VDRKD_ligand - VDRKD_vehicle)", 
        "VDR_knockdown", "vector_control + calcitriol + VDR", "VDR_overexpression + vehicle + VDR", "VDR_overexpression + calcitriol + VDR", 
        "(OE_ligand - OE_vehicle)", "VDR_knockdown + H3K4me3", "VDR_knockdown + H3K27ac", "VDR_knockdown + H3K27me3"), 
    denominator = c("siRNA_control + vehicle", "VDR_knockdown + vehicle", "siRNA_control + vehicle", "siRNA_control + calcitriol_2nM", 
        "(siNT_ligand - siNT_vehicle)", "siRNA_control", "vector_control + vehicle + VDR", "vector_control + vehicle + VDR", 
        "vector_control + calcitriol + VDR", "(vector_ligand - vector_vehicle)", "siRNA_control + H3K4me3", "siRNA_control + H3K27ac", 
        "siRNA_control + H3K27me3"), primary_or_secondary = c("PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", "PRIMARY", 
        "PRIMARY", "PRIMARY", "SECONDARY", "SECONDARY", "SECONDARY", "PRIMARY", "PRIMARY", "PRIMARY"), planned_method = c(rep("log2(RPKM+pseudocount) + limma factorial model", 
        5), "union peaks + replicate BigWig signal quantification; merged-peak set comparison", rep("VDR-target BED fragments with matched IgG; peak calling/occupancy comparison", 
        4), rep("mark-specific BED fragments with matched IgG; region-level signal/peak comparison", 3)), scientific_role = c("ligand-response program with intact VDR", 
        "residual ligand response after VDR knockdown", "baseline VDR-dependent transcriptome", "VDR dependence under ligand", 
        "formal VDR-dependent ligand-response interaction", "VDR-dependent chromatin accessibility", "ligand-induced VDR occupancy", 
        "VDR abundance effect on occupancy under vehicle", "VDR abundance effect on occupancy under ligand", "construct-by-ligand occupancy interaction", 
        "active-promoter histone response", "active-enhancer/promoter histone response", "repressive-histone response"), 
    notes = c(rep("No decidualization treatment is present in this sequencing dataset.", 5), "No ligand or decidualization factor is present in GSE306094.", 
        rep("Use matched IgG controls; interpret as in-vitro T-HESC cistrome evidence.", 4), rep("No ligand or decidualization factor is present in the histone-mark block.", 
            3)), stringsAsFactors = FALSE)

write_csv_utf8(contrast_plan, file.path(dirs$results, "STEP03A_v4_contrast_plan_LOCKED.csv"))

file_rows <- list()

file_rows[[1]] <- data.frame(accession = "GSE254251", assay = "RNA-seq", gsm = "", file_role = "processed_expression_matrix", 
    file_name = "GSE254251_RNAseq-RPKM.txt.gz", url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", "GSE254nnn/GSE254251/suppl/", 
        "GSE254251_RNAseq-RPKM.txt.gz"), priority = "CORE", analysis_use = paste("sample QC; log2-RPKM limma factorial contrasts;", 
        "frozen VDR and ligand programs"), limitation = paste("RPKM is not raw count data;", "DESeq2/edgeR count models are prohibited"), 
    stringsAsFactors = FALSE)

file_rows[[2]] <- data.frame(accession = "GSE306094", assay = "ATAC-seq", gsm = "", file_role = "merged_peak_set_siNT", 
    file_name = "GSE306094_siNT-merged.bed.gz", url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", "GSE306nnn/GSE306094/suppl/", 
        "GSE306094_siNT-merged.bed.gz"), priority = "CORE", analysis_use = "control consensus open-chromatin peak set", 
    limitation = "condition-level merged peaks; no replicate uncertainty", stringsAsFactors = FALSE)

file_rows[[3]] <- data.frame(accession = "GSE306094", assay = "ATAC-seq", gsm = "", file_role = "merged_peak_set_siVDR", 
    file_name = "GSE306094_siVDR-merged.bed.gz", url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", "GSE306nnn/GSE306094/suppl/", 
        "GSE306094_siVDR-merged.bed.gz"), priority = "CORE", analysis_use = "VDR-knockdown consensus open-chromatin peak set", 
    limitation = "condition-level merged peaks; no replicate uncertainty", stringsAsFactors = FALSE)

atac_samples <- locked[locked$accession == "GSE306094", , drop = FALSE]

for (i in seq_len(nrow(atac_samples))) {
    url <- to_https_ncbi(atac_samples$supplementary_file[i])
    file_rows[[length(file_rows) + 1L]] <- data.frame(accession = "GSE306094", assay = "ATAC-seq", gsm = atac_samples$gsm[i], 
        file_role = "replicate_bigwig_same_depth", file_name = basename(url), url = url, priority = "HIGH_VALUE", 
        analysis_use = paste("replicate-level normalized accessibility signal;", "QC and union-peak quantification"), 
        limitation = paste("normalized signal is not a raw fragment count;", "count-based differential accessibility is prohibited"), 
        stringsAsFactors = FALSE)
}

cutrun_samples <- locked[locked$accession == "GSE306127", , drop = FALSE]

for (i in seq_len(nrow(cutrun_samples))) {
    url <- to_https_ncbi(cutrun_samples$supplementary_file[i])
    file_rows[[length(file_rows) + 1L]] <- data.frame(accession = "GSE306127", assay = "CUT&RUN-seq", gsm = cutrun_samples$gsm[i], 
        file_role = "replicate_BED_same_depth", file_name = basename(url), url = url, priority = "CORE", analysis_use = paste(cutrun_samples$experiment_block[i], 
            cutrun_samples$genetic_context[i], cutrun_samples$ligand_status_locked[i], cutrun_samples$assay_target_locked[i], 
            sep = "; "), limitation = paste("processed/downsampled BED;", "raw-read QC and original-depth analysis are not available"), 
        stringsAsFactors = FALSE)
}

processed_manifest <- do.call(rbind, file_rows)

processed_manifest <- processed_manifest[!is.na(processed_manifest$url) & nzchar(processed_manifest$url), , drop = FALSE]

processed_manifest$file_id <- paste0("STEP03-FILE-", sprintf("%03d", seq_len(nrow(processed_manifest))))

processed_manifest <- processed_manifest[c("file_id", setdiff(names(processed_manifest), "file_id"))]

if (nrow(processed_manifest) != 57L) {
    stop("Expected 57 direct processed files but created ", nrow(processed_manifest), ".")
}

log_event("Auditing remote size/status for", nrow(processed_manifest), "processed files")

remote_rows <- vector("list", nrow(processed_manifest))

for (i in seq_len(nrow(processed_manifest))) {
    audit <- remote_size_audit(processed_manifest$url[i])
    remote_rows[[i]] <- data.frame(file_id = processed_manifest$file_id[i], http_status = audit$status_code, size_bytes = audit$size_bytes, 
        size_gb = format_gb(audit$size_bytes), audit_method = audit$method, audit_error = audit$error, remote_accessible = (is.finite(audit$status_code) && 
            audit$status_code < 400), stringsAsFactors = FALSE)
    log_event(i, "/", nrow(processed_manifest), processed_manifest$file_name[i], "status=", audit$status_code, 
        "size=", audit$size_bytes)
}

remote_audit <- do.call(rbind, remote_rows)

processed_manifest <- merge(processed_manifest, remote_audit, by = "file_id", all.x = TRUE, sort = FALSE)

processed_manifest$download_decision <- ifelse(processed_manifest$priority == "CORE" & processed_manifest$remote_accessible, 
    "APPROVE_TARGETED_DOWNLOAD", ifelse(processed_manifest$priority == "HIGH_VALUE" & processed_manifest$remote_accessible, 
        "APPROVE_IF_SIZE_BUDGET_ACCEPTABLE", "HOLD_REMOTE_AUDIT_OR_ACCESS_FAILURE"))

write_csv_utf8(processed_manifest, file.path(dirs$results, "STEP03A_v4_processed_file_remote_size_manifest.csv"))

size_summary <- aggregate(list(n_files = processed_manifest$file_id, known_size_bytes = processed_manifest$size_bytes), 
    by = list(accession = processed_manifest$accession, assay = processed_manifest$assay, priority = processed_manifest$priority, 
        download_decision = processed_manifest$download_decision), FUN = function(x) {
        if (is.character(x)) {
            length(x)
        }
        else {
            sum(x, na.rm = TRUE)
        }
    })

size_summary$known_size_gb <- size_summary$known_size_bytes/1024^3

write_csv_utf8(size_summary, file.path(dirs$results, "STEP03A_v4_processed_file_size_summary.csv"))

design_pass <- all(expected_checks$pass)

all_core_accessible <- all(processed_manifest$remote_accessible[processed_manifest$priority == "CORE"])

known_core_fraction <- mean(is.finite(processed_manifest$size_bytes[processed_manifest$priority == "CORE"]))

automated_gate <- if (design_pass && all_core_accessible) {
    "READY_FOR_TARGETED_PROCESSED_DOWNLOAD"
} else if (design_pass) {
    "DESIGN_LOCKED_BUT_REMOTE_ACCESS_REVIEW_REQUIRED"
} else {
    "HOLD_DESIGN_LOCK_FAILURE"
}

gate_table <- data.frame(gate_component = c("All 66 samples resolved", "Expected factorial/block counts", "No sequencing decidualization factor", 
    "Contrast plan locked", "Direct processed-file manifest", "All core files remotely accessible", "Known core file-size fraction", 
    "RAW TAR/FASTQ/BAM download", "STEP03A v4 gate"), status = c(ifelse(all(locked$design_lock_complete), "PASS", 
    "FAIL"), ifelse(design_pass, "PASS", "FAIL"), ifelse(all(locked$decidualization_status_locked == "not_applied_in_sequencing_dataset"), 
    "PASS", "FAIL"), "PASS", paste0("N=", nrow(processed_manifest)), ifelse(all_core_accessible, "PASS", "FAIL"), 
    sprintf("%.1f%%", 100 * known_core_fraction), "HOLD", automated_gate), required_action = c("Retain locked sample manifest.", 
    "Retain validation checks.", "Do not label sequencing contrasts as decidualization contrasts.", "Use only prespecified contrasts unless marked exploratory.", 
    "Retain URLs, roles and limitations.", "Review any failed URL before download.", "Use known size summary for storage planning.", 
    "Do not download unless a later raw-reprocessing decision is approved.", "Generate targeted processed-file download script after review."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03A_v4_GATE_TABLE.csv"))

manual_review <- c("# STEP03A v4 manual experimental-design lock", "", paste0("- Run ID: `", run_id, "`"), paste0("- Gate: `", 
    automated_gate, "`"), "", "## Locked dataset designs", "", "### GSE254251 RNA-seq", "", "- 12 T-HESC samples.", 
    "- 2x2 factorial design: siRNA control versus VDR knockdown and vehicle versus 2 nM calcitriol.", "- Three biological wells per condition.", 
    "- No decidualization induction is represented in the sequencing samples.", "", "### GSE306094 ATAC-seq", "", 
    "- Six T-HESC samples.", "- siRNA control versus VDR knockdown, three replicates per group.", "- No ligand or decidualization factor is represented.", 
    "- Two merged peak sets and six sample-level same-depth BigWig files are available.", "", "### GSE306127 CUT&RUN-seq", 
    "", "- 48 T-HESC samples in two separate experimental blocks.", "- VDR-cistrome block: vector control versus VDR overexpression, vehicle versus calcitriol, VDR versus IgG, n=3 per cell.", 
    "- Histone-mark block: siRNA control versus VDR knockdown, H3K4me3/H3K27ac/H3K27me3/IgG, n=3 per cell.", "- No decidualization induction is represented in either sequencing block.", 
    "", "## Scientific consequence", "", "These sequencing datasets can define VDR-dependent transcription, ligand response, chromatin accessibility, VDR occupancy and histone-mark programs. They cannot independently define a decidualization-response program because no sequencing contrast applies a decidualization treatment.", 
    "", "A decidualization program must therefore be obtained from an independent decidualization time course or a prespecified external gene set and kept separate from the VDR/ligand programs.", 
    "", "## Download policy", "", "- Download only the 57 direct processed files after this size audit is reviewed.", 
    "- Keep GSE306094_RAW.tar and GSE306127_RAW.tar on hold.", "- Do not use DESeq2/edgeR on the RPKM matrix.", 
    "- Do not claim replicate-level differential accessibility from merged BED files alone.", "- Use matched IgG controls for CUT&RUN interpretation.")

write_text(manual_review, file.path(dirs$results, "STEP03A_v4_MANUAL_DESIGN_LOCK_REPORT.md"))

write_text(c("# STEP03A v4 figure classification", "", "- Figure class: NONE", "- This substep locks experimental design and audits remote files.", 
    "- Tables and manifests are the authoritative outputs.", "- No QC-only or submission-candidate figure is scientifically necessary."), 
    file.path(dirs$figures, "STEP03A_v4_FIGURE_CLASSIFICATION.md"))

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "prepare_experimental_design.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03A_v4", 
    run_id = run_id, expected_samples = list(GSE254251 = 12, GSE306094 = 6, GSE306127 = 48), decidualization_factor = "not_applied_in_sequencing_datasets", 
    processed_file_count = 57, raw_archive_policy = "HOLD", analysis_started = FALSE, figure_class = "NONE")

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03A_v4_parameters_v0.1.json"), pretty = TRUE, 
    auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A_v4.txt"))

for (d in c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", "06_results", 
    "07_figures_tables", "08_qc")) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_groups <- list(`01_inputs` = c(v3_manifest_path, v3_decision_path), `02_code` = list.files(dirs$code, pattern = "06_step03A_v4_", 
    full.names = TRUE), `03_parameters` = list.files(dirs$parameters, pattern = "^STEP03A_v4_", full.names = TRUE), 
    `04_environment` = list.files(dirs$environment, pattern = "STEP03A_v4", full.names = TRUE), `05_execution` = log_path, 
    `06_results` = list.files(dirs$results, pattern = "^STEP03A_v4_", full.names = TRUE), `07_figures_tables` = file.path(dirs$figures, 
        "STEP03A_v4_FIGURE_CLASSIFICATION.md"))

for (subdir in names(copy_groups)) {
    paths <- copy_groups[[subdir]]
    paths <- paths[file.exists(paths)]
    for (p in paths) {
        file.copy(p, file.path(dirs$archive, subdir, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03A_v4_archive_manifest.csv"))

write_text(c("# STEP03A v4 archive status", "", paste0("- Run ID: ", run_id), paste0("- Gate: ", automated_gate), 
    "- Sample design: LOCKED", "- Contrast plan: LOCKED", "- Processed-file remote audit: COMPLETE", "- Processed-file download: NOT_STARTED", 
    "- RAW TAR/FASTQ/BAM: HOLD", "- Scientific analysis: NOT_STARTED", "- Figure production: NOT_REQUIRED", "- Full STEP03 module: NOT_LOCKED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03A_v4", paste0("- Run ID: ", run_id), paste0("- Gate: ", 
    automated_gate), "- Experimental design: LOCKED", "- Contrast plan: LOCKED", "- Processed data acquisition: PENDING", 
    "- RAW archives: HOLD", "- Perturbation analysis: NOT_STARTED", "- Figure production: NOT_STARTED", "- Full module status: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03A v4 design lock and remote size audit completed.\n")

cat("Gate: ", automated_gate, "\n")

cat("Locked samples: ", nrow(locked), "\n")

cat("Direct processed files audited: ", nrow(processed_manifest), "\n")

cat("No data files were downloaded.\n")

cat("No scientific figure was generated.\n\n")

cat("Upload these files for review:\n")

cat("  STEP03A_v4_MANUAL_DESIGN_LOCK_REPORT.md\n")

cat("  STEP03A_v4_sample_design_LOCKED.csv\n")

cat("  STEP03A_v4_design_cell_count_LOCKED.csv\n")

cat("  STEP03A_v4_design_validation_checks.csv\n")

cat("  STEP03A_v4_contrast_plan_LOCKED.csv\n")

cat("  STEP03A_v4_processed_file_remote_size_manifest.csv\n")

cat("  STEP03A_v4_processed_file_size_summary.csv\n")

cat("  STEP03A_v4_GATE_TABLE.csv\n")

cat("============================================================\n")

