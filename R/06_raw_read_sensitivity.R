#!/usr/bin/env Rscript
# Historical analysis module. Use a working copy of the original project data layout.
.vdr_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(.vdr_arg)) { .vdr_module_file <- normalizePath(sub("^--file=","",.vdr_arg[[1]]),winslash="/",mustWork=TRUE) } else {
  .vdr_files <- vapply(sys.frames(),function(f)if(is.null(f$ofile)) "" else as.character(f$ofile),character(1))
  .vdr_files <- .vdr_files[file.exists(.vdr_files)]
  if(!length(.vdr_files)) stop("Source the complete module path or run it with Rscript.")
  .vdr_module_file <- normalizePath(tail(.vdr_files,1),winslash="/",mustWork=TRUE)
}
#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 06_raw_read_sensitivity.R <input_dir> <output_dir>")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  read.delim(file.path(input_dir, name), check.names = FALSE,
             stringsAsFactors = FALSE, quote = "")
}

write_tsv <- function(x, name) {
  write.table(x, file.path(output_dir, name), sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "NA", eol = "\n")
}

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  unname(cor(x[ok], y[ok], method = method))
}

sign_concordance <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & x != 0 & y != 0
  if (!any(ok)) return(NA_real_)
  mean(sign(x[ok]) == sign(y[ok]))
}

required_samples <- c(
  "SRR27731230", "SRR27731229", "SRR27731228",
  "SRR27731227", "SRR27731226", "SRR27731225",
  "SRR27731224", "SRR27731223", "SRR27731222",
  "SRR27731221", "SRR27731220", "SRR27731219"
)
required_groups <- c(
  rep("NT_VEH", 3), rep("NT_CAL", 3),
  rep("KD_VEH", 3), rep("KD_CAL", 3)
)

counts <- read_tsv("Raw_Count_Matrix.tsv")
manifest <- read_tsv("Sample_Manifest.tsv")
gene_map <- read_tsv("Gene_ID_Gene_Name_Mapping.tsv")
program <- read.csv(file.path(input_dir, "Original_Program_Candidates.csv"),
                    check.names = FALSE, stringsAsFactors = FALSE)
original_all <- read.csv(file.path(input_dir, "Original_Interaction_All_Features.csv"),
                         check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(identical(colnames(counts), c("Geneid", required_samples)))
stopifnot(identical(manifest$srr, required_samples))
stopifnot(identical(manifest$group, required_groups))
stopifnot(!anyDuplicated(counts$Geneid))
stopifnot(identical(counts$Geneid, gene_map$gene_id))
stopifnot(nrow(program) == 436L, !anyDuplicated(program$feature_id))

count_matrix <- as.matrix(counts[, required_samples, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
stopifnot(all(is.finite(count_matrix)), all(count_matrix >= 0),
          all(count_matrix == round(count_matrix)))
rownames(count_matrix) <- counts$Geneid

group <- factor(manifest$group,
                levels = c("NT_VEH", "NT_CAL", "KD_VEH", "KD_CAL"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- manifest$srr
stopifnot(qr(design)$rank == ncol(design))

y <- DGEList(counts = count_matrix)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

interaction_test <- glmQLFTest(fit, contrast = c(1, -1, -1, 1))
interaction <- topTags(interaction_test, n = Inf, sort.by = "none")$table
interaction$FDR <- p.adjust(interaction$PValue, method = "BH")
interaction$gene_id <- rownames(interaction)
annotation <- gene_map[match(interaction$gene_id, gene_map$gene_id), ]
interaction_out <- data.frame(
  gene_id = interaction$gene_id,
  gene_id_without_version = annotation$gene_id_without_version,
  gene_name = annotation$gene_name,
  gene_type = annotation$gene_type,
  logFC = interaction$logFC,
  logCPM = interaction$logCPM,
  F = interaction$F,
  PValue = interaction$PValue,
  FDR = interaction$FDR,
  rawcount_candidate = interaction$FDR < 0.05 & abs(interaction$logFC) >= 0.5,
  stringsAsFactors = FALSE
)
write_tsv(interaction_out, "Raw_Count_Interaction_All_Genes.tsv")
write_tsv(interaction_out[interaction_out$rawcount_candidate, ],
          "Raw_Count_Interaction_Candidates.tsv")

symbol_to_ids <- split(gene_map$gene_id, gene_map$gene_name)
map_symbols <- function(symbols) {
  id_lists <- lapply(symbols, function(symbol) unique(symbol_to_ids[[symbol]]))
  status <- vapply(id_lists, function(ids) {
    if (length(ids) == 0L) "UNMAPPED"
    else if (length(ids) == 1L) "EXACT_UNIQUE"
    else "MULTIPLE_GENE_IDS"
  }, character(1))
  mapped_id <- vapply(id_lists, function(ids) {
    if (length(ids) == 1L) ids[[1]] else NA_character_
  }, character(1))
  data.frame(mapping_status = status, mapped_gene_id = mapped_id,
             stringsAsFactors = FALSE)
}

program_mapping <- map_symbols(program$feature_id)
comparison <- cbind(program, program_mapping)
idx <- match(comparison$mapped_gene_id, interaction_out$gene_id)
comparison$count_model_tested <- !is.na(idx)
comparison$comparison_evaluable <-
  comparison$mapping_status == "EXACT_UNIQUE" & comparison$count_model_tested
comparison$raw_interaction_logFC <- interaction_out$logFC[idx]
comparison$raw_interaction_FDR <- interaction_out$FDR[idx]
comparison$raw_threshold_candidate <- with(
  comparison,
  comparison_evaluable & raw_interaction_FDR < 0.05 &
    abs(raw_interaction_logFC) >= 0.5
)
comparison$interaction_sign_concordant <- with(
  comparison,
  comparison_evaluable & sign(interaction_logFC) == sign(raw_interaction_logFC)
)
write_tsv(comparison, "Original_Program_Raw_Count_Comparison.tsv")

evaluable <- comparison[comparison$comparison_evaluable, , drop = FALSE]
original_mapping <- map_symbols(original_all$feature_id)
genomewide <- cbind(original_all, original_mapping)
idx_all <- match(genomewide$mapped_gene_id, interaction_out$gene_id)
genomewide$count_model_tested <- !is.na(idx_all)
genomewide$common_tested <-
  genomewide$mapping_status == "EXACT_UNIQUE" & genomewide$count_model_tested
genomewide$raw_interaction_logFC <- interaction_out$logFC[idx_all]
genomewide$raw_interaction_FDR <- interaction_out$FDR[idx_all]
write_tsv(genomewide, "Genomewide_Interaction_Comparison.tsv")
common <- genomewide[genomewide$common_tested, , drop = FALSE]

summary <- data.frame(
  metric = c(
    "edgeR_version", "input_gene_count", "filtered_tested_gene_count",
    "rawcount_candidate_count", "original_program_count",
    "original_evaluable_count", "candidate_interaction_Pearson",
    "candidate_interaction_Spearman", "candidate_interaction_sign_concordance",
    "candidate_threshold_overlap_count", "genomewide_common_tested_count",
    "genomewide_interaction_Pearson", "genomewide_interaction_Spearman",
    "genomewide_interaction_sign_concordance"
  ),
  value = c(
    as.character(packageVersion("edgeR")), nrow(count_matrix), nrow(y),
    sum(interaction_out$rawcount_candidate), nrow(program), nrow(evaluable),
    safe_cor(evaluable$interaction_logFC, evaluable$raw_interaction_logFC, "pearson"),
    safe_cor(evaluable$interaction_logFC, evaluable$raw_interaction_logFC, "spearman"),
    sign_concordance(evaluable$interaction_logFC, evaluable$raw_interaction_logFC),
    sum(evaluable$raw_threshold_candidate), nrow(common),
    safe_cor(common$logFC, common$raw_interaction_logFC, "pearson"),
    safe_cor(common$logFC, common$raw_interaction_logFC, "spearman"),
    sign_concordance(common$logFC, common$raw_interaction_logFC)
  ),
  stringsAsFactors = FALSE
)
write_tsv(summary, "Raw_Count_Sensitivity_Summary.tsv")
capture.output(sessionInfo(), file = file.path(output_dir, "R_Session_Info.txt"))
