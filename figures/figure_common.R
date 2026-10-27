# Shared input, package, font and export handling for manuscript figures.
figure_input <- function(id, key) {
  row <- figure_index[figure_index$figure == id & figure_index$key == key, , drop=FALSE]
  if (nrow(row) != 1L) stop('Unknown or ambiguous figure input: ',id,' / ',key)
  p <- file.path(figure_data_root, row$file)
  if (!file.exists(p)) stop('Missing figure input: ', p)
  p
}

figure_preflight <- function(id, required) {
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop('Install the required R packages: ',paste(missing,collapse=', '))
  rows <- figure_index[figure_index$figure == id,,drop=FALSE]
  if (!nrow(rows)) stop('No input list for ',id)
  invisible(vapply(rows$key,function(key)figure_input(id,key),character(1)))
  match <- systemfonts::match_fonts('Arial')
  if (!nrow(match) || !file.exists(match$path[[1]])) stop('Arial is required to preserve the figure typography.')
  info <- systemfonts::font_info(match$path[[1]],index=match$index[[1]])
  if (!grepl('^Arial$',info$family[[1]],ignore.case=TRUE)) stop('Arial is not installed; font substitution is not permitted.')
  if (.Platform$OS.type=='windows') grDevices::windowsFonts(Arial=grDevices::windowsFont('Arial'))
  invisible(TRUE)
}

# Keep warnings visible. This helper does not change statistical data or fit models.
capture_stage <- function(stage, expr) force(expr)

figure_export <- function(plot, id, width, height) {
  root <- normalizePath(Sys.getenv('VDR_FIGURE_OUTPUT',unset=file.path(repo_root,'results','figures')),winslash='/',mustWork=FALSE)
  if (root==figure_data_root || startsWith(paste0(root,'/'),paste0(figure_data_root,'/'))) stop('Output must not overwrite the input-data directory.')
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  paths <- setNames(file.path(root,paste0(id,c('.pdf','.png','.tiff'))),c('pdf','png','tiff'))
  if (any(file.exists(paths))) stop('Output already exists; choose a new VDR_FIGURE_OUTPUT directory: ',root)
  draw <- function(device) {device();tryCatch(print(plot),finally=grDevices::dev.off())}
  draw(function()grDevices::cairo_pdf(paths[['pdf']],width=width/25.4,height=height/25.4,family='Arial',bg='white',fallback_resolution=600))
  draw(function()ragg::agg_png(paths[['png']],width=width,height=height,units='mm',res=600,background='white'))
  draw(function()ragg::agg_tiff(paths[['tiff']],width=width,height=height,units='mm',res=600,compression='lzw',background='white'))
  message(id,': ',paste(paths,collapse='; '))
  invisible(paths)
}
