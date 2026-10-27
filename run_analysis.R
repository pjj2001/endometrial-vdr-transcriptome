#!/usr/bin/env Rscript
# Run one selected step; listing and checking never launch an analysis.
args <- commandArgs(TRUE)
fa <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(!length(fa)) stop('Use Rscript run_analysis.R --list or --step <ID>.')
root <- dirname(normalizePath(sub('^--file=','',fa[1]),winslash='/',mustWork=TRUE))
steps <- read.csv(file.path(root,'RUN_ORDER.csv'),stringsAsFactors=FALSE,check.names=FALSE)
if(identical(args,'--list')) {
  for(i in seq_len(nrow(steps))) cat(steps$step[i],steps$script[i],steps$purpose[i],sep=' | ',fill=TRUE)
  quit(status=0)
}
if(identical(args,'--check')) {
  for(p in steps$script) parse(file.path(root,p),encoding='UTF-8')
  parse(file.path(root,'analysis','analysis_paths.R'),encoding='UTF-8')
  cat('All listed step files exist and parse. No analysis was executed.\n')
  data_root <- Sys.getenv('VDR_PROJECT_ROOT',unset='')
  if(!nzchar(data_root)) cat('VDR_PROJECT_ROOT is not set. Set it before running a step.\n')
  else {
    cat('Working data folder:',normalizePath(data_root,winslash='/',mustWork=FALSE),'\n')
    for(p in c('03_data','05_analysis_steps')) cat(p,if(dir.exists(file.path(data_root,p)))'present' else 'missing','\n')
  }
  quit(status=0)
}
if(length(args)!=2L || args[1]!='--step' || !args[2] %in% steps$step) stop('Usage: Rscript run_analysis.R --step <ID>. Use --list for step IDs.')
data_root <- Sys.getenv('VDR_PROJECT_ROOT',unset='')
if(!nzchar(data_root) || !dir.exists(data_root)) stop('Set VDR_PROJECT_ROOT to a writable working copy of the study data. See DATA_ACCESS.md.')
setwd(normalizePath(data_root,winslash='/',mustWork=TRUE))
script <- file.path(root,steps$script[match(args[2],steps$step)])
status <- system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(script)))
quit(status=status)
