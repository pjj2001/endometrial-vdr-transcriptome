#!/usr/bin/env Rscript
# Individual analysis step; paths are supplied by analysis_paths.R.
# Statistical expressions, filters and seeds are retained from the recorded project script.
.vdr_f <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(!length(.vdr_f)) stop("Run this step with Rscript, or use run_analysis.R --step <ID>.")
.vdr_current_step <- normalizePath(sub("^--file=","",.vdr_f[[1]]),winslash="/",mustWork=TRUE)
source(file.path(dirname(dirname(.vdr_current_step)),"analysis_paths.R"),encoding="UTF-8")
# A01: GEO metadata and supplementary listings

options(stringsAsFactors = FALSE)

options(timeout = 7200)

options(warn = 1)

project_root <- .vdr_project_root()

module_name <- "STEP03_VDR_perturbation_regulons"

module_root <- file.path(project_root, "05_analysis_steps", module_name)

run_id <- paste0("STEP03A_", format(Sys.time(), "%Y%m%d_%H%M%S"))

accessions <- c("GSE254251", "GSE306094", "GSE306127")

assay_map <- c(GSE254251 = "RNA-seq", GSE306094 = "ATAC-seq", GSE306127 = "CUT&RUN-seq")

paper_reference <- list(title = paste("Impact of vitamin D deficiency on defective endometrial", "decidualization and the repressive role of vitamin D receptor", 
    "(VDR) in the epigenomic network"), journal = "Pharmacological Research", year = 2026, volume = 227, article = "108162", 
    doi = "10.1016/j.phrs.2026.108162", pmid = "41856382")

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

substep_root <- file.path(dirs$inputs, "STEP03A_GEO_design_audit")

for (d in c("geo_metadata", "supplementary_inventory", "design_manifests", "download_holds")) {
    dir.create(file.path(substep_root, d), recursive = TRUE, showWarnings = FALSE)
}

for (acc in accessions) {
    for (d in c("metadata", "supplementary_inventory", "processed_supplementary_HOLD", "raw_archive_HOLD")) {
        dir.create(file.path(project_root, "03_data", "raw_external_READ_ONLY", acc, d), recursive = TRUE, showWarnings = FALSE)
    }
}

log_path <- file.path(dirs$execution, paste0(run_id, "_execution.log"))

log_event <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "))
    cat(msg, "\n")
    cat(msg, "\n", file = log_path, append = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

cran_required <- c("digest", "jsonlite")

cran_missing <- cran_required[!vapply(cran_required, requireNamespace, logical(1), quietly = TRUE)]

if (length(cran_missing) > 0L) {
    install.packages(cran_missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

required_packages <- c("GEOquery", "digest", "jsonlite")

if (!all(vapply(required_packages, requireNamespace, logical(1), quietly = TRUE))) {
    stop("Required packages could not be loaded.")
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

meta_value <- function(meta_list, key) {
    value <- meta_list[[key]]
    if (is.null(value) || length(value) == 0L) {
        return("")
    }
    value <- unlist(value, use.names = FALSE)
    value <- value[!is.na(value) & nzchar(trimws(as.character(value)))]
    paste(unique(trimws(as.character(value))), collapse = " | ")
}

get_geo_safely <- function(accession, destination) {
    args <- list(GEO = accession, GSEMatrix = FALSE, getGPL = FALSE, destdir = destination)
    valid_args <- intersect(names(args), names(formals(GEOquery::getGEO)))
    do.call(GEOquery::getGEO, args[valid_args])
}

list_supplementary_safely <- function(accession, destination) {
    fn <- GEOquery::getGEOSuppFiles
    args <- list(GEO = accession, makeDirectory = FALSE, baseDir = destination, fetch_files = FALSE)
    valid_args <- intersect(names(args), names(formals(fn)))
    do.call(fn, args[valid_args])
}

normalize_supp_listing <- function(object, accession) {
    if (is.null(object) || length(object) == 0L) {
        return(data.frame(accession = accession, file_reference = character(), file_name = character(), size_bytes = numeric(), 
            mtime = character(), file_class = character(), recommended_action = character(), download_status = character(), 
            notes = character(), stringsAsFactors = FALSE))
    }
    x <- as.data.frame(object, stringsAsFactors = FALSE)
    refs <- rownames(x)
    if (is.null(refs) || all(refs %in% as.character(seq_len(nrow(x))))) {
        candidate_cols <- names(x)[grepl("file|name|url|path", names(x), ignore.case = TRUE)]
        if (length(candidate_cols) > 0L) {
            refs <- as.character(x[[candidate_cols[1]]])
        }
        else {
            refs <- rep("", nrow(x))
        }
    }
    size_candidates <- names(x)[grepl("^size$|size_bytes|length", names(x), ignore.case = TRUE)]
    size_bytes <- if (length(size_candidates) > 0L) {
        suppressWarnings(as.numeric(x[[size_candidates[1]]]))
    }
    else {
        rep(NA_real_, nrow(x))
    }
    mtime_candidates <- names(x)[grepl("mtime|modified|date", names(x), ignore.case = TRUE)]
    mtime <- if (length(mtime_candidates) > 0L) {
        as.character(x[[mtime_candidates[1]]])
    }
    else {
        rep("", nrow(x))
    }
    file_name <- basename(sub("\\?.*$", "", refs))
    lower_name <- tolower(file_name)
    processed_pattern <- paste(c("count", "matrix", "normalized", "expression", "deseq", "edge?r", "peak", "narrowpeak", 
        "broadpeak", "\\.bed", "\\.bw", "bigwig", "\\.csv", "\\.tsv", "\\.txt", "\\.xlsx", "\\.rds", "\\.rdata", 
        "\\.h5"), collapse = "|")
    raw_pattern <- paste(c("_raw\\.tar", "\\.fastq", "\\.fq", "\\.bam", "\\.cram", "\\.sra"), collapse = "|")
    file_class <- ifelse(grepl(raw_pattern, lower_name), "raw_or_large_archive", ifelse(grepl(processed_pattern, 
        lower_name), "processed_or_analysis_ready_candidate", "manual_review_required"))
    recommended_action <- ifelse(file_class == "processed_or_analysis_ready_candidate", "REVIEW_THEN_DOWNLOAD", 
        ifelse(file_class == "raw_or_large_archive", "HOLD_DO_NOT_DOWNLOAD_YET", "MANUAL_REVIEW"))
    data.frame(accession = accession, file_reference = refs, file_name = file_name, size_bytes = size_bytes, mtime = mtime, 
        file_class = file_class, recommended_action = recommended_action, download_status = "NOT_DOWNLOADED_BY_STEP03A", 
        notes = "", stringsAsFactors = FALSE)
}

infer_vdr_manipulation <- function(text) {
    x <- tolower(text)
    if (grepl("sivdr|si-vdr|vdr si|vdr knockdown|knockdown of vdr", x)) {
        return("VDR_knockdown")
    }
    if (grepl("sint|si-nt|non.?target|scrambl|si.?control", x)) {
        return("siRNA_control")
    }
    if (grepl("vdr overexpress|overexpression of vdr|oe.?vdr|vdr.?oe", x)) {
        return("VDR_overexpression")
    }
    if (grepl("empty vector|vector control|control vector|mock transfect", x)) {
        return("vector_control")
    }
    if (grepl("\\bvdr\\b", x, perl = TRUE)) {
        return("VDR_related_unresolved")
    }
    "baseline_or_unresolved"
}

infer_ligand <- function(text) {
    x <- tolower(text)
    if (grepl("1,25|1alpha,25|1alpha25|calcitriol|dihydroxyvitamin d", x)) {
        return("calcitriol")
    }
    if (grepl("vehicle|ethanol control|solvent control", x)) {
        return("vehicle_control")
    }
    "unresolved"
}

infer_decidualization <- function(text) {
    x <- tolower(text)
    if (grepl("undecidual|non.?decidual|without decidual", x)) {
        return("undecidualized")
    }
    if (grepl("decidual|medroxyprogesterone|\\bmpa\\b|8-br-camp|camp", x, perl = TRUE)) {
        return("decidualization_induced_or_related")
    }
    "unresolved"
}

infer_cutrun_target <- function(text) {
    x <- tolower(text)
    if (grepl("h3k27ac", x)) 
        return("H3K27ac")
    if (grepl("h3k4me3", x)) 
        return("H3K4me3")
    if (grepl("h3k27me3", x)) 
        return("H3K27me3")
    if (grepl("\\bvdr\\b", x, perl = TRUE)) 
        return("VDR")
    if (grepl("\\bigg\\b", x, perl = TRUE)) 
        return("IgG_control")
    "unresolved"
}

infer_replicate <- function(text) {
    x <- tolower(text)
    patterns <- c("replicate[ _-]*([0-9]+)", "\\brep[ _-]*([0-9]+)\\b", "[ _-]r([0-9]+)\\b")
    for (pattern in patterns) {
        hit <- regexec(pattern, x, perl = TRUE)
        parts <- regmatches(x, hit)[[1]]
        if (length(parts) >= 2L) {
            return(parts[2])
        }
    }
    ""
}

flatten_gse_samples <- function(gse_object, accession, assay) {
    gsm_list <- GEOquery::GSMList(gse_object)
    if (length(gsm_list) == 0L) {
        stop("No GSM samples found for ", accession)
    }
    out <- lapply(seq_along(gsm_list), function(i) {
        gsm <- gsm_list[[i]]
        meta <- GEOquery::Meta(gsm)
        gsm_id <- names(gsm_list)[i]
        if (is.null(gsm_id) || is.na(gsm_id) || !nzchar(gsm_id)) {
            gsm_id <- meta_value(meta, "geo_accession")
        }
        fields <- c("title", "source_name_ch1", "organism_ch1", "characteristics_ch1", "treatment_protocol_ch1", 
            "growth_protocol_ch1", "extract_protocol_ch1", "molecule_ch1", "library_strategy", "library_source", 
            "library_selection", "instrument_model", "relation", "supplementary_file")
        values <- vapply(fields, function(field) meta_value(meta, field), character(1))
        combined_text <- paste(values, collapse = " | ")
        data.frame(accession = accession, assay = assay, gsm = gsm_id, sample_order = i, title = values["title"], 
            source_name = values["source_name_ch1"], organism = values["organism_ch1"], characteristics = values["characteristics_ch1"], 
            treatment_protocol = values["treatment_protocol_ch1"], growth_protocol = values["growth_protocol_ch1"], 
            extract_protocol = values["extract_protocol_ch1"], molecule = values["molecule_ch1"], library_strategy = values["library_strategy"], 
            library_source = values["library_source"], library_selection = values["library_selection"], instrument_model = values["instrument_model"], 
            relation = values["relation"], supplementary_file = values["supplementary_file"], vdr_manipulation = infer_vdr_manipulation(combined_text), 
            ligand_status = infer_ligand(combined_text), decidualization_status = infer_decidualization(combined_text), 
            cutrun_target = if (assay == "CUT&RUN-seq") {
                infer_cutrun_target(combined_text)
            }
            else {
                "not_applicable"
            }, replicate = infer_replicate(combined_text), provisional_condition_id = paste(infer_vdr_manipulation(combined_text), 
                infer_ligand(combined_text), infer_decidualization(combined_text), if (assay == "CUT&RUN-seq") {
                  infer_cutrun_target(combined_text)
                }
                else {
                  "NA"
                }, sep = "__"), manual_review_required = any(c(infer_vdr_manipulation(combined_text) %in% c("baseline_or_unresolved", 
                "VDR_related_unresolved"), infer_ligand(combined_text) == "unresolved", infer_decidualization(combined_text) == 
                "unresolved", assay == "CUT&RUN-seq" && infer_cutrun_target(combined_text) == "unresolved")), stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

make_design_summary <- function(sample_manifest) {
    aggregate(list(n_samples = sample_manifest$gsm), by = list(accession = sample_manifest$accession, assay = sample_manifest$assay, 
        vdr_manipulation = sample_manifest$vdr_manipulation, ligand_status = sample_manifest$ligand_status, decidualization_status = sample_manifest$decidualization_status, 
        cutrun_target = sample_manifest$cutrun_target, provisional_condition_id = sample_manifest$provisional_condition_id), 
        FUN = length)
}

log_event("Starting", run_id)

source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

if (!is.null(source_file) && file.exists(source_file)) {
    file.copy(source_file, file.path(dirs$code, "prepare_geo_metadata.R"), overwrite = TRUE, copy.mode = TRUE, 
        copy.date = TRUE)
}

pre_run_manifest <- data.frame(object_id = accessions, accession = accessions, assay = unname(assay_map[accessions]), 
    source = paste0("NCBI GEO ", accessions), expected_role = c("transcriptomic perturbation anchor", "chromatin accessibility support", 
        "direct occupancy or histone-mark support"), download_policy = "metadata_and_listing_only_in_STEP03A", 
    status = "PRE_RUN_LOCKED", stringsAsFactors = FALSE)

write_csv_utf8(pre_run_manifest, file.path(dirs$inputs, "STEP03A_input_manifest_PRE_RUN_LOCK.csv"))

parameters <- list(schema_version = "1.0", project = "Nutrition_Endometrial_VDR_v2", module = "STEP03", substep = "STEP03A", 
    run_id = run_id, accessions = accessions, assay_map = as.list(assay_map), paper_reference = paper_reference, 
    data_policy = list(fetch_GEO_metadata = TRUE, list_supplementary_files = TRUE, download_supplementary_files = FALSE, 
        download_FASTQ_BAM_RAW_tar = FALSE), design_policy = list(factors_are_provisional = TRUE, manual_review_required_before_download = TRUE, 
        no_DE_or_peak_analysis_in_STEP03A = TRUE, no_manuscript_figure_in_STEP03A = TRUE))

jsonlite::write_json(parameters, file.path(dirs$parameters, "STEP03A_parameters_v0.1.json"), pretty = TRUE, auto_unbox = TRUE)

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_PRE_RUN_STEP03A.txt"))

all_samples <- list()

all_supp <- list()

retrieval_status <- list()

for (acc in accessions) {
    assay <- unname(assay_map[acc])
    log_event("Retrieving GEO metadata for", acc, assay)
    geo_dir <- file.path(project_root, "03_data", "raw_external_READ_ONLY", acc, "metadata")
    gse_result <- tryCatch({
        gse <- get_geo_safely(acc, geo_dir)
        saveRDS(gse, file.path(geo_dir, paste0(acc, "_GEOquery_full_metadata.rds")), compress = "xz")
        samples <- flatten_gse_samples(gse, acc, assay)
        write_csv_utf8(samples, file.path(geo_dir, paste0(acc, "_sample_metadata_parsed.csv")))
        all_samples[[acc]] <- samples
        list(metadata_status = "PASS", metadata_error = "", n_samples = nrow(samples))
    }, error = function(e) {
        all_samples[[acc]] <<- data.frame()
        list(metadata_status = "FAIL", metadata_error = conditionMessage(e), n_samples = 0L)
    })
    log_event("Listing supplementary files for", acc)
    supp_dir <- file.path(project_root, "03_data", "raw_external_READ_ONLY", acc, "supplementary_inventory")
    supp_result <- tryCatch({
        listing_object <- list_supplementary_safely(acc, supp_dir)
        listing <- normalize_supp_listing(listing_object, acc)
        write_csv_utf8(listing, file.path(supp_dir, paste0(acc, "_supplementary_file_inventory.csv")))
        all_supp[[acc]] <- listing
        list(supplementary_status = "PASS", supplementary_error = "", n_supplementary_files = nrow(listing))
    }, error = function(e) {
        all_supp[[acc]] <<- data.frame()
        list(supplementary_status = "FAIL", supplementary_error = conditionMessage(e), n_supplementary_files = 0L)
    })
    retrieval_status[[acc]] <- data.frame(accession = acc, assay = assay, metadata_status = gse_result$metadata_status, 
        metadata_error = gse_result$metadata_error, n_samples = gse_result$n_samples, supplementary_status = supp_result$supplementary_status, 
        supplementary_error = supp_result$supplementary_error, n_supplementary_files = supp_result$n_supplementary_files, 
        stringsAsFactors = FALSE)
}

retrieval_table <- do.call(rbind, retrieval_status)

sample_manifest <- do.call(rbind, all_samples[vapply(all_samples, nrow, integer(1)) > 0L])

supp_manifest <- do.call(rbind, all_supp[vapply(all_supp, nrow, integer(1)) > 0L])

if (is.null(sample_manifest)) {
    sample_manifest <- data.frame()
}

if (is.null(supp_manifest)) {
    supp_manifest <- data.frame()
}

write_csv_utf8(retrieval_table, file.path(dirs$results, "STEP03A_GEO_retrieval_status.csv"))

write_csv_utf8(sample_manifest, file.path(dirs$results, "STEP03A_sample_design_manifest_PROVISIONAL.csv"))

write_csv_utf8(supp_manifest, file.path(dirs$results, "STEP03A_supplementary_file_inventory.csv"))

if (nrow(sample_manifest) > 0L) {
    design_summary <- make_design_summary(sample_manifest)
} else {
    design_summary <- data.frame()
}

write_csv_utf8(design_summary, file.path(dirs$results, "STEP03A_design_cell_count_summary.csv"))

unresolved_samples <- if (nrow(sample_manifest) > 0L) {
    sample_manifest[sample_manifest$manual_review_required, , drop = FALSE]
} else {
    data.frame()
}

write_csv_utf8(unresolved_samples, file.path(dirs$results, "STEP03A_samples_requiring_manual_design_review.csv"))

download_candidates <- if (nrow(supp_manifest) > 0L) {
    supp_manifest[supp_manifest$recommended_action == "REVIEW_THEN_DOWNLOAD", , drop = FALSE]
} else {
    data.frame()
}

download_holds <- if (nrow(supp_manifest) > 0L) {
    supp_manifest[supp_manifest$recommended_action != "REVIEW_THEN_DOWNLOAD", , drop = FALSE]
} else {
    data.frame()
}

write_csv_utf8(download_candidates, file.path(dirs$results, "STEP03A_processed_file_download_candidates.csv"))

write_csv_utf8(download_holds, file.path(dirs$results, "STEP03A_download_holds_and_manual_review.csv"))

cross_assay_coverage <- if (nrow(sample_manifest) > 0L) {
    as.data.frame.matrix(table(sample_manifest$provisional_condition_id, sample_manifest$assay), stringsAsFactors = FALSE)
} else {
    data.frame()
}

if (nrow(cross_assay_coverage) > 0L) {
    cross_assay_coverage$provisional_condition_id <- rownames(cross_assay_coverage)
    rownames(cross_assay_coverage) <- NULL
    cross_assay_coverage <- cross_assay_coverage[c("provisional_condition_id", setdiff(names(cross_assay_coverage), 
        "provisional_condition_id"))]
}

write_csv_utf8(cross_assay_coverage, file.path(dirs$results, "STEP03A_cross_assay_condition_coverage.csv"))

all_metadata_pass <- all(retrieval_table$metadata_status == "PASS")

all_supp_pass <- all(retrieval_table$supplementary_status == "PASS")

all_have_samples <- all(retrieval_table$n_samples > 0L)

manual_review_n <- nrow(unresolved_samples)

automated_gate <- if (all_metadata_pass && all_supp_pass && all_have_samples) {
    "PROCEED_TO_MANUAL_DESIGN_LOCK"
} else {
    "HOLD_RETRIEVAL_OR_METADATA_FAILURE"
}

gate_table <- data.frame(gate_component = c("Three GEO metadata records retrieved", "Three supplementary inventories retrieved", 
    "Each accession contains sample records", "Treatment-factor parsing", "Heavy-data download", "STEP03A gate"), 
    status = c(ifelse(all_metadata_pass, "PASS", "FAIL"), ifelse(all_supp_pass, "PASS", "FAIL"), ifelse(all_have_samples, 
        "PASS", "FAIL"), ifelse(manual_review_n == 0L, "AUTO_PARSE_COMPLETE", paste0("MANUAL_REVIEW_REQUIRED_N=", 
        manual_review_n)), "NOT_STARTED_BY_DESIGN", automated_gate), required_action = c("Review accession, title and assay.", 
        "Review file names, sizes and classes.", "Confirm sample counts against GEO.", "Resolve treatment, ligand, decidualization and target labels.", 
        "Select processed files only after design lock.", "Do not start differential analysis before manual design lock."), 
    stringsAsFactors = FALSE)

write_csv_utf8(gate_table, file.path(dirs$results, "STEP03A_GATE_TABLE.csv"))

sample_count_text <- paste(paste0(retrieval_table$accession, "=", retrieval_table$n_samples), collapse = "; ")

supp_count_text <- paste(paste0(retrieval_table$accession, "=", retrieval_table$n_supplementary_files), collapse = "; ")

report <- c("# STEP03A GEO metadata and supplementary-file audit", "", paste0("- Run ID: `", run_id, "`"), paste0("- Automated gate: `", 
    automated_gate, "`"), "", "## Scientific role", "", "- GSE254251: RNA-seq perturbation anchor.", "- GSE306094: ATAC-seq chromatin-accessibility support.", 
    "- GSE306127: CUT&RUN direct-occupancy or histone-mark support.", "", "These datasets will define independent experimental programs before", 
    "projection into patient cohorts.", "", "## Retrieval summary", "", paste0("- Sample counts: ", sample_count_text, 
        "."), paste0("- Supplementary-file counts: ", supp_count_text, "."), paste0("- Samples requiring manual treatment-factor review: ", 
        manual_review_n, "."), "", "## Download policy", "", "- STEP03A downloaded GEO metadata only.", "- Supplementary files were inventoried but not downloaded.", 
    "- FASTQ, BAM, SRA and generic RAW archives remain on hold.", "- Processed counts, matrices, peak files and signal tracks must be", 
    "  selected after manual design review.", "", "## Required manual review", "", "- Verify the meaning of every title and characteristics field.", 
    "- Confirm VDR knockdown/control/overexpression labels.", "- Confirm calcitriol versus vehicle labels.", "- Confirm decidualized versus nondecidualized conditions.", 
    "- Confirm CUT&RUN targets and matching controls.", "- Confirm biological replicates and whether libraries are paired.", 
    "- Confirm which conditions are shared across RNA-seq, ATAC-seq and", "  CUT&RUN-seq.", "", "## Gate interpretation", 
    "", if (automated_gate == "PROCEED_TO_MANUAL_DESIGN_LOCK") {
        "- Metadata acquisition succeeded. Proceed to manual design lock and processed-file selection."
    } else {
        "- Retrieval failed for at least one accession. Do not select or download analysis files."
    }, "", "This report does not authorize differential expression, peak calling", "or manuscript-figure production.")

write_text(report, file.path(dirs$results, "STEP03A_DESIGN_AUDIT_REPORT.md"))

checklist <- c("# STEP03A manual design review checklist", "", "- [ ] GEO title and assay verified for all three accessions.", 
    "- [ ] Sample counts verified.", "- [ ] Cell line and biological material verified.", "- [ ] VDR manipulation resolved for every sample.", 
    "- [ ] Ligand/vehicle status resolved for every sample.", "- [ ] Decidualization status resolved for every sample.", 
    "- [ ] CUT&RUN target and control resolved.", "- [ ] Replicate structure resolved.", "- [ ] Cross-assay condition correspondence resolved.", 
    "- [ ] Processed files selected from the supplementary inventory.", "- [ ] Heavy raw files remain on hold unless a justified reprocessing", 
    "      plan is approved.", "- [ ] STEP03A dataset-role and design lock decision recorded.", "", "No STEP03B analysis begins until all mandatory items are checked.")

write_text(checklist, file.path(dirs$qc, "STEP03A_MANUAL_DESIGN_REVIEW_CHECKLIST.md"))

write_text(c("# STEP03A figure classification", "", "- Figure class: NONE", "- Reason: metadata and supplementary-file audit only.", 
    "- No main, supplementary-candidate or final figure is required.", "- Any temporary count plot would be A_QC_only and is omitted because", 
    "  the design tables are more informative."), file.path(dirs$figures, "STEP03A_FIGURE_CLASSIFICATION.md"))

archive_subdirs <- c("00_archive_control", "01_inputs", "02_code", "03_parameters", "04_environment", "05_execution", 
    "06_results", "07_figures_tables", "08_qc")

for (d in archive_subdirs) {
    dir.create(file.path(dirs$archive, d), recursive = TRUE, showWarnings = FALSE)
}

copy_map <- list(`01_inputs` = c(file.path(dirs$inputs, "STEP03A_input_manifest_PRE_RUN_LOCK.csv")), `02_code` = c(file.path(dirs$code, 
    "prepare_geo_metadata.R")), `03_parameters` = c(file.path(dirs$parameters, "STEP03A_parameters_v0.1.json")), 
    `04_environment` = c(file.path(dirs$environment, "sessionInfo_PRE_RUN_STEP03A.txt")), `05_execution` = c(log_path), 
    `06_results` = list.files(dirs$results, full.names = TRUE), `07_figures_tables` = c(file.path(dirs$figures, 
        "STEP03A_FIGURE_CLASSIFICATION.md")), `08_qc` = c(file.path(dirs$qc, "STEP03A_MANUAL_DESIGN_REVIEW_CHECKLIST.md")))

for (target_subdir in names(copy_map)) {
    files <- copy_map[[target_subdir]]
    files <- files[file.exists(files)]
    for (f in files) {
        file.copy(f, file.path(dirs$archive, target_subdir, basename(f)), overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
}

write_text(capture.output(sessionInfo()), file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A.txt"))

file.copy(file.path(dirs$environment, "sessionInfo_POST_RUN_STEP03A.txt"), file.path(dirs$archive, "04_environment", 
    "sessionInfo_POST_RUN_STEP03A.txt"), overwrite = TRUE)

archive_files <- list.files(dirs$archive, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)

archive_manifest <- data.frame(relative_path = substring(gsub("\\\\", "/", archive_files), nchar(gsub("\\\\", "/", 
    dirs$archive)) + 2L), size_bytes = as.numeric(file.info(archive_files)$size), md5 = vapply(archive_files, md5_file, 
    character(1)), sha256 = vapply(archive_files, sha256_file, character(1)), stringsAsFactors = FALSE)

write_csv_utf8(archive_manifest, file.path(dirs$archive, "00_archive_control", "STEP03A_archive_manifest.csv"))

write_text(c("# STEP03A archive status", "", paste0("- Run ID: ", run_id), "- Archive mode: POST_RUN_CAPTURE", 
    paste0("- Automated gate: ", automated_gate), "- Manual design review: PENDING", "- Supplementary-file download: NOT_STARTED", 
    "- Scientific analysis: NOT_STARTED", "- Figure production: NOT_REQUIRED", "- STEP03 module lock: NOT_LOCKED"), 
    file.path(dirs$archive, "00_archive_control", "ARCHIVE_STATUS.md"))

write_text(c("# STEP03 Module Status", "", paste0("- Active substep: STEP03A"), paste0("- Accepted audit run ID: ", 
    run_id), paste0("- Automated gate: ", automated_gate), "- Manual design lock: PENDING", "- Processed data acquisition: PENDING", 
    "- Perturbation analysis: NOT_STARTED", "- Module archive: WORKING_ARCHIVE_CREATED", "- Module status: NOT_LOCKED"), 
    file.path(dirs$control, "STATUS.md"))

log_event("Completed", run_id, "Gate:", automated_gate)

cat("\n============================================================\n")

cat("STEP03A GEO metadata audit completed.\n")

cat("Automated gate: ", automated_gate, "\n")

cat("No heavy supplementary files were downloaded.\n")

cat("No scientific figure was generated.\n\n")

cat("Upload these files for review:\n")

cat("  06_results/STEP03A_DESIGN_AUDIT_REPORT.md\n")

cat("  06_results/STEP03A_GEO_retrieval_status.csv\n")

cat("  06_results/STEP03A_sample_design_manifest_PROVISIONAL.csv\n")

cat("  06_results/STEP03A_design_cell_count_summary.csv\n")

cat("  06_results/STEP03A_samples_requiring_manual_design_review.csv\n")

cat("  06_results/STEP03A_supplementary_file_inventory.csv\n")

cat("  06_results/STEP03A_processed_file_download_candidates.csv\n")

cat("  06_results/STEP03A_cross_assay_condition_coverage.csv\n")

cat("  06_results/STEP03A_GATE_TABLE.csv\n")

cat("============================================================\n")

