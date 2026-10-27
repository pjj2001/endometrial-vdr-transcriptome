#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
file_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(!length(file_arg)) stop('Run with Rscript run_figures.R <figure> [--check].')
root <- dirname(normalizePath(sub('^--file=','',file_arg[[1]]),winslash='/',mustWork=TRUE))
scripts <- c('2'='figure2_perturbation_occupancy.R','3'='figure3_bulk_projection.R',
 '4'='figure4_fertile_temporal_analysis.R','5'='figure5_donor_level_rif_analysis.R',
 '6'='figure6_cross_dataset_synthesis.R', 'S1'='supplementary_figures.R',
 'S2'='supplementary_figures.R','S3'='supplementary_figures.R','S4'='supplementary_figures.R')
id <- toupper(setdiff(args,'--check')[1])
if(is.na(id) || !id %in% names(scripts)) stop('Choose 2, 3, 4, 5, 6, S1, S2, S3 or S4. Figure 1 uses figures/figure1_workflow.py.')
extra <- c(if(startsWith(id,'S')) id, if('--check' %in% args) '--check')
status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path(root,'figures',scripts[[id]])),extra))
quit(status=status)
