# Shared file-location helpers. No statistical analysis is performed here.
.vdr_step_dir <- dirname(.vdr_current_step)
.vdr_code_root <- dirname(dirname(.vdr_step_dir))
.vdr_project_root <- function() {
  p <- Sys.getenv('VDR_PROJECT_ROOT',unset='')
  if(!nzchar(p)) stop('Set VDR_PROJECT_ROOT to a writable working copy of the study data folders. See RUN_ORDER.csv and DATA_ACCESS.md.')
  normalizePath(p,winslash='/',mustWork=TRUE)
}
.vdr_report_root <- function() {
  p <- Sys.getenv('VDR_REPORT_ROOT',unset=file.path(.vdr_project_root(),'reports'))
  dir.create(p,recursive=TRUE,showWarnings=FALSE)
  normalizePath(p,winslash='/',mustWork=TRUE)
}
.vdr_script_path <- function(name) {
  m <- read.csv(file.path(.vdr_code_root,'RUN_ORDER.csv'),stringsAsFactors=FALSE)
  i <- match(basename(name),m$original_script)
  if(is.na(i)) i <- match(basename(name),basename(m$script))
  if(!is.na(i)) return(file.path(.vdr_code_root,m$script[i]))
  # Self-copy and code-availability checks use the public step itself.
  if(basename(name)==basename(.vdr_current_step)) return(.vdr_current_step)
  stop('Source script is not included: ',name)
}
.vdr_rebase <- function(x) {
  if(!is.character(x)) return(x)
  z <- gsub('\\\\','/',x)
  marker <- '/Nutrition_Endometrial_VDR_v2/'
  hit <- !is.na(z) & grepl(marker,z,fixed=TRUE) & grepl('^[A-Za-z]:/',z)
  if(any(hit)) {
    suffix <- substring(z[hit],regexpr(marker,z[hit],fixed=TRUE)+nchar(marker))
    z[hit] <- file.path(.vdr_project_root(),suffix)
  }
  x[hit] <- z[hit]
  x
}
.vdr_rebase_table <- function(x) {
  if(is.data.frame(x)) for(n in names(x)) if(is.character(x[[n]])) x[[n]] <- .vdr_rebase(x[[n]])
  x
}
.vdr_read_csv <- function(file,...) .vdr_rebase_table(utils::read.csv(.vdr_rebase(file),...))
.vdr_fread <- function(input,...){
  if(is.character(input)) input <- .vdr_rebase(input)
  .vdr_rebase_table(data.table::fread(input,...))
}
