#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A02: Read GEO family SOFT metadata

options(stringsAsFactors = FALSE)

options(timeout = 7200)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03A_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"))

accessions <- c("GSE254251", "GSE306094", "GSE306127")

assay_map <- c(GSE254251 = "RNA-seq", GSE306094 = "ATAC-seq", GSE306127 = "CUT&RUN-seq")

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

for (acc in accessions) {
    dir.create(file.path(project_root, "03_data", "raw_external_READ_ONLY", acc, "metadata"), recursive = TRUE, 
        showWarnings = FALSE)
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

cran_required <- c("curl", "digest", "jsonlite")

missing_pkgs <- cran_required[!vapply(cran_required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

if (!all(vapply(cran_required, requireNamespace, logical(1), quietly = TRUE))) {
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

geo_series_stub <- function(accession) {
    sub("[0-9]{3}$", "nnn", accession)
}

family_soft_urls <- function(accession) {
    stub <- geo_series_stub(accession)
    c(paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", stub, "/", accession, "/soft/", accession, "_family.soft.gz"), 
        paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", stub, "/", accession, "/soft/", accession, "_family.soft.gz"))
}

download_family_soft <- function(accession, destination) {
    urls <- family_soft_urls(accession)
    errors <- character()
    for (url in urls) {
        part <- paste0(destination, ".part")
        if (file.exists(part)) {
            unlink(part, force = TRUE)
        }
        handle <- curl::new_handle(useragent = paste0("Nutrition_Endometrial_VDR_v2/", run_id, " contact=local-research-script"), 
            followlocation = TRUE, timeout = 7200)
        ok <- tryCatch({
            curl::curl_download(url = url, destfile = part, quiet = FALSE, mode = "wb", handle = handle)
            file.exists(part) && file.info(part)$size > 1000
        }, error = function(e) {
            errors <<- c(errors, paste0(url, " :: ", conditionMessage(e)))
            FALSE
        })
        if (isTRUE(ok)) {
            con <- gzfile(part, open = "rt")
            first_lines <- tryCatch(readLines(con, n = 20L, warn = FALSE), error = function(e) character())
            close(con)
            valid_soft <- any(grepl("^\\^SERIES = ", first_lines))
            if (!valid_soft) {
                errors <- c(errors, paste0(url, " :: downloaded file is not a valid GEO family SOFT file"))
                unlink(part, force = TRUE)
                next
            }
            if (file.exists(destination)) {
                old_path <- paste0(destination, ".previous_", format(Sys.time(), "%Y%m%d_%H%M%S"))
                file.rename(destination, old_path)
            }
            renamed <- file.rename(part, destination)
            if (!renamed) {
                stop("Validated family SOFT file could not be renamed:\n", destination)
            }
            return(list(status = "PASS", url = url, error = ""))
        }
    }
    list(status = "FAIL", url = "", error = paste(errors, collapse = " || "))
}

parse_key_value_lines <- function(lines, prefix) {
    pattern <- paste0("^!", prefix, "_([^=]+?) = (.*)$")
    hits <- grep(pattern, lines, value = TRUE, perl = TRUE)
    if (length(hits) == 0L) {
        return(list())
    }
    keys <- sub(pattern, "\\1", hits, perl = TRUE)
    values <- sub(pattern, "\\2", hits, perl = TRUE)
    split_values <- split(values, keys)
    lapply(split_values, function(x) {
        paste(unique(trimws(x)), collapse = " | ")
    })
}

parse_geo_family_soft <- function(path, accession, assay) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE, encoding = "UTF-8")
    if (!any(grepl(paste0("^\\^SERIES = ", accession, "$"), lines))) {
        stop("Family SOFT does not contain expected series: ", accession)
    }
    sample_starts <- grep("^\\^SAMPLE = GSM[0-9]+$", lines)
    if (length(sample_starts) == 0L) {
        stop("No sample records found in family SOFT: ", accession)
    }
    all_record_starts <- grep("^\\^[A-Z]+ = ", lines)
    sample_rows <- vector("list", length(sample_starts))
    for (i in seq_along(sample_starts)) {
        start <- sample_starts[i]
        next_starts <- all_record_starts[all_record_starts > start]
        end <- if (length(next_starts) > 0L) {
            min(next_starts) - 1L
        }
        else {
            length(lines)
        }
        block <- lines[start:end]
        gsm <- sub("^\\^SAMPLE = ", "", block[1])
        meta <- parse_key_value_lines(block, "Sample")
        get_value <- function(key) {
            value <- meta[[key]]
            if (is.null(value)) {
                ""
            }
            else {
                as.character(value)
            }
        }
        sample_rows[[i]] <- data.frame(accession = accession, assay = assay, gsm = gsm, sample_order = i, title = get_value("title"), 
            source_name_ch1 = get_value("source_name_ch1"), organism_ch1 = get_value("organism_ch1"), characteristics_ch1 = get_value("characteristics_ch1"), 
            treatment_protocol_ch1 = get_value("treatment_protocol_ch1"), growth_protocol_ch1 = get_value("growth_protocol_ch1"), 
            extract_protocol_ch1 = get_value("extract_protocol_ch1"), molecule_ch1 = get_value("molecule_ch1"), 
            library_strategy = get_value("library_strategy"), library_source = get_value("library_source"), library_selection = get_value("library_selection"), 
            instrument_model = get_value("instrument_model"), relation = get_value("relation"), supplementary_file = paste(c(get_value("supplementary_file"), 
                get_value("supplementary_file_1"), get_value("supplementary_file_2"), get_value("supplementary_file_3"))[nzchar(c(get_value("supplementary_file"), 
                get_value("supplementary_file_1"), get_value("supplementary_file_2"), get_value("supplementary_file_3")))], 
                collapse = " | "), stringsAsFactors = FALSE)
    }
    series_start <- grep(paste0("^\\^SERIES = ", accession, "$"), lines)[1]
    next_record <- all_record_starts[all_record_starts > series_start]
    series_end <- if (length(next_record) > 0L) {
        min(next_record) - 1L
    }
    else {
        length(lines)
    }
    series_block <- lines[series_start:series_end]
    series_meta <- parse_key_value_lines(series_block, "Series")
    series_df <- data.frame(accession = accession, assay = assay, field = names(series_meta), value = unlist(series_meta, 
        use.names = FALSE), stringsAsFactors = FALSE)
    list(samples = do.call(rbind, sample_rows), series = series_df, n_lines = length(lines))
}

combined_sample_text <- function(df) {
    apply(df[c("title", "source_name_ch1", "characteristics_ch1", "treatment_protocol_ch1", "growth_protocol_ch1", 
        "extract_protocol_ch1")], 1, function(x) {
        tolower(paste(x, collapse = " | "))
    })
}

infer_vdr_manipulation <- function(text) {
    if (grepl("sivdr|si-vdr|vdr knockdown|vdr-kd|knockdown of vdr", text)) {
        return("VDR_knockdown")
    }
    if (grepl("sint|si-nt|nontarget|non-target|non target|scramble|scrambled", text)) {
        return("siRNA_control")
    }
    if (grepl("plenti.?hvdr|plenti.?vdr|vdr overexpress|overexpression of vdr|oe.?vdr|vdr.?oe", text)) {
        return("VDR_overexpression")
    }
    if (grepl("empty vector|vector control|plenti control|mock transfect", text)) {
        return("vector_control")
    }
    if (grepl("\\bvdr\\b", text, perl = TRUE)) {
        return("VDR_related_unresolved")
    }
    "baseline_or_unresolved"
}

infer_ligand <- function(text) {
    if (grepl("1,25|1alpha.?25|calcitriol|dihydroxyvitamin d|1,25d", text)) {
        return("calcitriol")
    }
    if (grepl("vehicle|ethanol control|solvent control|veh\\b", text, perl = TRUE)) {
        return("vehicle_control")
    }
    "unresolved"
}

infer_decidualization <- function(text) {
    if (grepl("undecidual|non.?decidual|without decidual", text)) {
        return("undecidualized")
    }
    if (grepl("decidual|medroxyprogesterone|8-br-camp|8brcamp|\\bmpa\\b|\\bcamp\\b", text, perl = TRUE)) {
        return("decidualization_induced_or_related")
    }
    "unresolved"
}

infer_cutrun_target <- function(text) {
    if (grepl("h3k4me3", text)) {
        return("H3K4me3")
    }
    if (grepl("h3k27ac", text)) {
        return("H3K27ac")
    }
    if (grepl("h3k27me3", text)) {
        return("H3K27me3")
    }
    if (grepl("\\bvdr\\b", text, perl = TRUE)) {
        return("VDR")
    }
    if (grepl("\\bigg\\b", text, perl = TRUE)) {
        return("IgG_control")
    }
    "unresolved"
}

infer_replicate <- function(text) {
    patterns <- c("replicate[ _-]*([0-9]+)", "\\brep[ _-]*([0-9]+)\\b", "[ _-]r([0-9]+)\\b", "[ _-]([0-9]+)$")
    for (pattern in patterns) {
        hit <- regexec(pattern, text, perl = TRUE)
        parts <- regmatches(text, hit)[[1]]
        if (length(parts) >= 2L) {
            return(parts[2])
        }
    }
    ""
}

add_provisional_factors <- function(samples) {
    text <- combined_sample_text(samples)
    samples$vdr_manipulation <- vapply(text, infer_vdr_manipulation, character(1))
    samples$ligand_status <- vapply(text, infer_ligand, character(1))
    samples$decidualization_status <- vapply(text, infer_decidualization, character(1))
    samples$cutrun_target <- ifelse(samples$assay == "CUT&RUN-seq", vapply(text, infer_cutrun_target, character(1)), 
        "not_applicable")
    samples$replicate <- vapply(text, infer_replicate, character(1))
    samples$provisional_condition_id <- paste(samples$vdr_manipulation, samples$ligand_status, samples$decidualization_status, 
        samples$cutrun_target, sep = "__")
    samples$manual_review_required <- (samples$vdr_manipulation %in% c("baseline_or_unresolved", "VDR_related_unresolved")) | 
        samples$ligand_status == "unresolved" | samples$decidualization_status == "unresolved" | (samples$assay == 
        "CUT&RUN-seq" & samples$cutrun_target == "unresolved") | !nzchar(samples$replicate)
    samples
}

make_design_summary <- function(samples) {
    aggregate(list(n_samples = samples$gsm), by = list(accession = samples$accession, assay = samples$assay, vdr_manipulation = samples$vdr_manipulation, 
        ligand_status = samples$ligand_status, decidualization_status = samples$decidualization_status, cutrun_target = samples$cutrun_target, 
        provisional_condition_id = samples$provisional_condition_id), FUN = length)
}

log_event("Starting", run_id)

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "read_geo_soft_metadata.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

pre_run <- data.frame(object_id = paste0(accessions, "_family_SOFT"), accession = accessions, assay = unname(assay_map[accessions]), 
    role = "official_GEO_family_metadata", prior_failure = "GEOquery_getGEO_HTTP_403", recovery_method = "direct_NCBI_GEO_family_SOFT", 
    heavy_data_download = FALSE, status = "PRE_RUN_LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(pre_run, file.path(dirs$inputs, "STEP03A_v3_input_manifest_PRE_RUN_LOCK.csv"))

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03A_v3", 
    run_id = run_id, accessions = accessions, recovery_method = "direct GEO family SOFT download and local parser", 
    heavy_data_download = FALSE, analysis_started = FALSE, figure_class = "NONE", manual_design_lock_required = TRUE)

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03A_v3_parameters_v0.1.json"), pretty = TRUE, 
    auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_PRE_RUN_STEP03A_v3.txt"))

download_rows <- list()

sample_rows <- list()

series_rows <- list()

for (acc in accessions) {
    assay <- unname(assay_map[acc])
    metadata_dir <- file.path(project_root, "03_data", "raw_external_READ_ONLY", acc, "metadata")
    soft_path <- file.path(metadata_dir, paste0(acc, "_family.soft.gz"))
    log_event("Downloading official family SOFT for", acc)
    download_result <- if (file.exists(soft_path) && file.info(soft_path)$size > 1000) {
        list(status = "PASS_EXISTING_FILE", url = "", error = "")
    }
    else {
        download_family_soft(acc, soft_path)
    }
    parse_status <- "NOT_RUN"
    parse_error <- ""
    n_samples <- 0L
    n_soft_lines <- NA_integer_
    if (grepl("^PASS", download_result$status)) {
        parsed <- tryCatch({
            x <- parse_geo_family_soft(soft_path, acc, assay)
            x$samples <- add_provisional_factors(x$samples)
            sample_rows[[acc]] <- x$samples
            series_rows[[acc]] <- x$series
            n_samples <- nrow(x$samples)
            n_soft_lines <- x$n_lines
            parse_status <- "PASS"
            write_csv_utf8(x$samples, file.path(metadata_dir, paste0(acc, "_sample_metadata_from_family_SOFT.csv")))
            write_csv_utf8(x$series, file.path(metadata_dir, paste0(acc, "_series_metadata_from_family_SOFT.csv")))
            TRUE
        }, error = function(e) {
            parse_status <<- "FAIL"
            parse_error <<- conditionMessage(e)
            FALSE
        })
    }
    download_rows[[acc]] <- data.frame(accession = acc, assay = assay, download_status = download_result$status, 
        source_url = download_result$url, download_error = download_result$error, local_path = soft_path, size_bytes = if (file.exists(soft_path)) {
            as.numeric(file.info(soft_path)$size)
        }
        else {
            NA_real_
        }, md5 = md5_file(soft_path), sha256 = sha256_file(soft_path), parse_status = parse_status, parse_error = parse_error, 
        n_samples = n_samples, n_soft_lines = n_soft_lines, stringsAsFactors = FALSE)
}

download_manifest <- do.call(rbind, download_rows)

write_csv_utf8(download_manifest, file.path(dirs$results, "STEP03A_v3_family_SOFT_download_manifest.csv"))

successful_samples <- sample_rows[vapply(sample_rows, nrow, integer(1)) > 0L]

successful_series <- series_rows[vapply(series_rows, nrow, integer(1)) > 0L]

samples <- if (length(successful_samples) > 0L) {
    do.call(rbind, successful_samples)
} else {
    data.frame()
}

series_metadata <- if (length(successful_series) > 0L) {
    do.call(rbind, successful_series)
} else {
    data.frame()
}

write_csv_utf8(samples, file.path(dirs$results, "STEP03A_v3_sample_design_manifest_PROVISIONAL.csv"))

write_csv_utf8(series_metadata, file.path(dirs$results, "STEP03A_v3_series_metadata.csv"))

design_summary <- if (nrow(samples) > 0L) {
    make_design_summary(samples)
} else {
    data.frame()
}

write_csv_utf8(design_summary, file.path(dirs$results, "STEP03A_v3_design_cell_count_summary.csv"))

unresolved <- if (nrow(samples) > 0L) {
    samples[samples$manual_review_required, , drop = FALSE]
} else {
    data.frame()
}

write_csv_utf8(unresolved, file.path(dirs$results, "STEP03A_v3_samples_requiring_manual_design_review.csv"))

cross_assay <- if (nrow(samples) > 0L) {
    as.data.frame.matrix(table(samples$provisional_condition_id, samples$assay), stringsAsFactors = FALSE)
} else {
    data.frame()
}

if (nrow(cross_assay) > 0L) {
    cross_assay$provisional_condition_id <- rownames(cross_assay)
    rownames(cross_assay) <- NULL
    cross_assay <- cross_assay[c("provisional_condition_id", setdiff(names(cross_assay), "provisional_condition_id"))]
}

write_csv_utf8(cross_assay, file.path(dirs$results, "STEP03A_v3_cross_assay_condition_coverage.csv"))

prior_inventory_path <- file.path(dirs$results, "STEP03A_supplementary_file_inventory.csv")

prior_inventory <- if (file.exists(prior_inventory_path)) {
    .vdr_read_csv(prior_inventory_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
    data.frame()
}

if (nrow(prior_inventory) > 0L) {
    decision <- prior_inventory
    decision$provisional_use <- ifelse(decision$file_name == "GSE254251_RNAseq-RPKM.txt.gz", paste("Expression matrix for QC, sample-level visualization,", 
        "limma-style contrasts and program derivation after design lock"), ifelse(grepl("siNT-merged\\.bed", decision$file_name), 
        paste("Control consensus open-chromatin peaks;", "peak-set comparison only"), ifelse(grepl("siVDR-merged\\.bed", 
            decision$file_name), paste("VDR-knockdown consensus open-chromatin peaks;", "peak-set comparison only"), 
            ifelse(decision$accession == "GSE306094" & grepl("_RAW\\.tar$", decision$file_name), paste("Potential replicate-level ATAC reanalysis;", 
                "content and size must be audited first"), ifelse(decision$accession == "GSE306127" & grepl("_RAW\\.tar$", 
                decision$file_name), paste("Required candidate for CUT&RUN replicate/peak analysis;", "content and size must be audited first"), 
                "Manual review")))))
    decision$main_limitation <- ifelse(decision$file_name == "GSE254251_RNAseq-RPKM.txt.gz", paste("RPKM is not an integer raw-count matrix;", 
        "do not use DESeq2/edgeR count models"), ifelse(grepl("merged\\.bed", decision$file_name), paste("Merged condition-level peaks do not support", 
        "replicate QC or differential accessibility"), ifelse(grepl("_RAW\\.tar$", decision$file_name), paste("Unknown content/size until archive listing;", 
        "may require substantial preprocessing"), "Manual review")))
    decision$next_action <- ifelse(decision$file_name == "GSE254251_RNAseq-RPKM.txt.gz", "DOWNLOAD_AFTER_MANUAL_DESIGN_LOCK", 
        ifelse(grepl("merged\\.bed", decision$file_name), "DOWNLOAD_AFTER_MANUAL_DESIGN_LOCK", ifelse(grepl("_RAW\\.tar$", 
            decision$file_name), "HOLD_PENDING_REMOTE_SIZE_AND_TAR_CONTENT_AUDIT", "MANUAL_REVIEW")))
} else {
    decision <- data.frame()
}

write_csv_utf8(decision, file.path(dirs$results, "STEP03A_v3_processed_data_decision_matrix.csv"))

all_download_pass <- all(grepl("^PASS", download_manifest$download_status))

all_parse_pass <- all(download_manifest$parse_status == "PASS")

all_have_samples <- all(download_manifest$n_samples > 0L)

unresolved_n <- nrow(unresolved)

automated_gate <- if (all_download_pass && all_parse_pass && all_have_samples) {
    "PROCEED_TO_MANUAL_DESIGN_LOCK"
} else {
    "HOLD_METADATA_RECOVERY_FAILURE"
}

gate_table <- data.frame(gate_component = c("Official family SOFT downloaded", "All family SOFT records parsed", 
    "Every accession has GSM samples", "Treatment-factor inference", "Processed-file selection", "Heavy RAW archives", 
    "STEP03A v3 gate"), status = c(ifelse(all_download_pass, "PASS", "FAIL"), ifelse(all_parse_pass, "PASS", "FAIL"), 
    ifelse(all_have_samples, "PASS", "FAIL"), paste0("MANUAL_REVIEW_REQUIRED_N=", unresolved_n), "PENDING_MANUAL_DESIGN_LOCK", 
    "HOLD", automated_gate), required_action = c("Retain URL, local path and hashes.", "Review parsed series and sample metadata.", 
    "Verify sample counts and titles against the experiment.", "Resolve every provisional factor and replicate.", 
    "Select only files needed for STEP03B.", "Audit remote size and TAR contents before download.", "Do not begin perturbation analysis before manual design lock."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03A_v3_GATE_TABLE.csv"))

sample_count_text <- paste(paste0(download_manifest$accession, "=", download_manifest$n_samples), collapse = "; ")

report <- c("# STEP03A v3 direct family-SOFT metadata recovery report", "", paste0("- Run ID: `", run_id, "`"), 
    paste0("- Automated gate: `", automated_gate, "`"), "", "## Reason for v3", "", "STEP03A v1 was blocked by GEOquery HTTP 403 errors.", 
    "STEP03A v2 successfully downloaded and parsed the official family SOFT files", "but stopped while classifying five supplementary-file rows because scalar", 
    "logical && was used instead of vectorized & inside ifelse().", "This v3 run reuses the downloaded metadata and corrects that logic.", 
    "", "## Metadata recovery", "", paste0("- Parsed sample counts: ", sample_count_text, "."), paste0("- Samples still requiring manual factor review: ", 
        unresolved_n, "."), "", "## Current file availability from STEP03A v1", "", "- GSE254251 provides one processed RPKM expression file.", 
    "- GSE306094 provides siNT and siVDR merged BED peak files plus a RAW TAR archive.", "- GSE306127 provides a RAW TAR archive and no separate processed peak file in the GEO supplementary inventory.", 
    "", "## Scientific implications", "", "- The RPKM matrix can support expression QC, limma-style contrasts", 
    "  and frozen-program derivation after sample-design lock, but it is", "  not suitable for DESeq2 or edgeR raw-count models.", 
    "- The two merged ATAC BED files can support consensus peak-set", "  comparisons, but not replicate-level differential accessibility.", 
    "- CUT&RUN analysis will require either the RAW archive or a rigorously", "  justified use of paper supplementary peak/overlap tables.", 
    "", "## Next decision", "", if (automated_gate == "PROCEED_TO_MANUAL_DESIGN_LOCK") {
        "- Review and manually lock every sample condition before downloading processed data."
    } else {
        "- Metadata recovery remains incomplete; do not download analysis data."
    }, "", "No differential analysis or scientific figure was generated.")

write_text(report, file.path(dirs$results, "STEP03A_v3_METADATA_RECOVERY_REPORT.md"))

write_text(c("# STEP03A v3 manual design lock checklist", "", "- [ ] Confirm all sample titles and GSM accessions.", 
    "- [ ] Confirm cell line/material and assay for all samples.", "- [ ] Confirm siNT versus siVDR.", "- [ ] Confirm VDR overexpression/vector controls.", 
    "- [ ] Confirm calcitriol versus vehicle.", "- [ ] Confirm decidualized versus nondecidualized state.", "- [ ] Confirm CUT&RUN target and IgG/input controls.", 
    "- [ ] Confirm biological versus technical replicates.", "- [ ] Confirm which conditions are directly comparable across assays.", 
    "- [ ] Lock primary RNA-seq contrasts.", "- [ ] Lock primary ATAC-seq comparisons.", "- [ ] Lock primary CUT&RUN comparisons.", 
    "- [ ] Approve processed-file download list.", "- [ ] Keep large RAW archives on hold until size/content audit.", 
    "", "STEP03B cannot begin until this checklist is completed."), file.path(dirs$qc, "STEP03A_v3_MANUAL_DESIGN_LOCK_CHECKLIST.md"))

write_text(c("# STEP03A v3 figure classification", "", "- Figure class: NONE", "- This is a metadata-recovery and design-audit substep.", 
    "- Tables are the authoritative outputs.", "- No QC-only, submission-candidate or final scientific figure is required."), 
    file.path(dirs$figures, "STEP03A_v3_FIGURE_CLASSIFICATION.md"))

for (d in c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", "06_results", 
    "07_figures_tables", "08_qc")) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_groups <- list(`01_inputs` = list.files(dirs$inputs, pattern = "^STEP03A_v3_", full.names = TRUE), `02_code` = list.files(dirs$code, 
    pattern = "04_step03A_v2_", full.names = TRUE), `03_parameters` = list.files(dirs$parameters, pattern = "^STEP03A_v3_", 
    full.names = TRUE), `04_environment` = list.files(dirs$environment, pattern = "STEP03A_v3", full.names = TRUE), 
    `05_execution` = log_path, `06_results` = list.files(dirs$results, pattern = "^STEP03A_v3_", full.names = TRUE), 
    `07_figures_tables` = file.path(dirs$figures, "STEP03A_v3_FIGURE_CLASSIFICATION.md"), `08_qc` = file.path(dirs$qc, 
        "STEP03A_v3_MANUAL_DESIGN_LOCK_CHECKLIST.md"))

for (subdir in names(copy_groups)) {
    paths <- copy_groups[[subdir]]
    paths <- paths[file.exists(paths)]
    for (p in paths) {
        file.copy(p, file.path(dirs$archive, subdir, basename(p)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A_v3.txt"))

file.copy(file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A_v3.txt"), file.path(dirs$archive, "04_environment", 
    "sessionInfo_POST_RUN_STEP03A_v3.txt"), overwrite = TRUE)

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03A_v3_archive_manifest.csv"))

write_text(c("# STEP03A v3 archive status", "", paste0("- Run ID: ", run_id), "- Archive mode: POST_RUN_CAPTURE", 
    paste0("- Automated gate: ", automated_gate), "- Manual design lock: PENDING", "- Processed data download: NOT_STARTED", 
    "- Heavy RAW download: HOLD", "- Scientific analysis: NOT_STARTED", "- Figure production: NOT_REQUIRED", "- STEP03 module status: NOT_LOCKED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", "- Active substep: STEP03A_v3", paste0("- Accepted metadata-recovery run: ", 
    run_id), paste0("- Automated gate: ", automated_gate), "- Manual design lock: PENDING", "- Processed data acquisition: PENDING", 
    "- Perturbation analysis: NOT_STARTED", "- Figure production: NOT_STARTED", "- Module status: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03A v3 metadata recovery completed.\n")

cat("Automated gate: ", automated_gate, "\n")

cat("No heavy data were downloaded.\n")

cat("No scientific figure was generated.\n\n")

cat("Upload these v2 outputs for review:\n")

cat("  STEP03A_v3_METADATA_RECOVERY_REPORT.md\n")

cat("  STEP03A_v3_family_SOFT_download_manifest.csv\n")

cat("  STEP03A_v3_series_metadata.csv\n")

cat("  STEP03A_v3_sample_design_manifest_PROVISIONAL.csv\n")

cat("  STEP03A_v3_design_cell_count_summary.csv\n")

cat("  STEP03A_v3_samples_requiring_manual_design_review.csv\n")

cat("  STEP03A_v3_cross_assay_condition_coverage.csv\n")

cat("  STEP03A_v3_processed_data_decision_matrix.csv\n")

cat("  STEP03A_v3_GATE_TABLE.csv\n")

cat("============================================================\n")

