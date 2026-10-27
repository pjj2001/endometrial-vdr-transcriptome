#!/usr/bin/env Rscript
# Read-only regression checks. This file never sources an analysis module.
args <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(!length(args)) stop('Run Rscript tests/test_repository.R')
root <- normalizePath(file.path(dirname(sub('^--file=','',args[[1]])),'..'),winslash='/',mustWork=TRUE)
scripts <- list.files(root,pattern='[.]R$',recursive=TRUE,full.names=TRUE)
for(p in scripts) parse(p,encoding='UTF-8')
plot_scripts <- list.files(file.path(root,'figures'),pattern='[.]R$',full.names=TRUE)
for(p in plot_scripts) {
 text <- paste(readLines(p,encoding='UTF-8',warn=FALSE),collapse='\n')
 for(forbidden in c('SKILL.md','FIGURE_STYLE_REFERENCE','LOCAL_R_LIBRARY','OPTIONAL_LOCAL_RUNTIME','USERPROFILE','sys.source(')) {
  if(grepl(forbidden,text,fixed=TRUE)) stop('Private plotting dependency in ',basename(p),': ',forbidden)
 }
}
index <- read.csv(file.path(root,'data/figures/inputs.csv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(nrow(index)>0,!anyDuplicated(paste(index$figure,index$key)))
stopifnot(!any(grepl('(^/|^[A-Za-z]:|(^|/)[.][.](/|$))',index$file)))
stopifnot(all(file.exists(file.path(root,'data/figures',index$file))))
stopifnot(all(c(paste0('Figure',2:6),paste0('S',1:4)) %in% index$figure))
stopifnot(file.exists(file.path(root,'config/Figure1_layout_spec.json')))
ignore <- readLines(file.path(root,'.gitignore'),warn=FALSE)
stopifnot('!data/figures/' %in% ignore,'!data/raw_read_sensitivity/' %in% ignore)
cat(length(scripts),'R scripts parsed;',nrow(index),'figure-input links present; private plotting dependencies absent.\n')
cat('No scientific analysis was executed. Package/font and plotting execution checks are separate.\n')
