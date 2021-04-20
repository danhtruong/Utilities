

MergeATACObjects <- function(seuratlist){
  seuratlist_1 <- seuratlist[[1]]
  
  for(i in 2:length(seuratlist)) {
    seuratlist_1 <- MergeWithRegions(
      object.1 = seuratlist_1,
      object.2 = seuratlist[[i]],
      sep.1 = c(":", "-"),
      sep.2 = c(":", "-"),
      distance = 200
    )
  }
  return(seuratlist_1)
}


CombineFragments <- function(
  filepath,
  output.path,
  compress = TRUE,
  index = TRUE,
  verbose = TRUE,
  ...
){
  require(dplyr)
  fragments.list <- list()
  for(i in 1:length(filepath)) {
    fragments.list[[i]] <- fread(
      file = filepath[i],
      col.names = c('chr', 'start', 'end', 'cell', 'count')
    )
    fragments.list[[i]]$cell <- paste0(
      sapply(strsplit(as.character(fragments.list[[i]]$cell), split="-"), "[[", 1),
      "-", i)
  }
  reads <-  bind_rows(fragments.list)
  reads <- reads[order(chr,start),]
  
  if (verbose) {
    message("Writing output")
  }
  fwrite(
    x = reads,
    #file = output.path,
    file = output.path,
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE,
    sep = '\t'
  )
  
  rm(reads)
  invisible(x = gc())
  if (compress) {
    if (verbose) {
      message("Compressing output")
    }
    outf <- bgzip(file = output.path)
    if (file.exists(outf)) {
      file.remove(output.path)
    }
    if (index) {
      if (verbose) {
        message("Building index")
      }
      index.file <- indexTabix(file = paste0(outf), format = 'bed', zeroBased = TRUE)
    }
  }
}


ATACMetaData <- function(seuratobj){
  seuratobj <- NucleosomeSignal(object = seuratobj)
  seuratobj$pct_reads_in_peaks <- seuratobj$peak_region_fragments / seuratobj$passed_filters * 100
  seuratobj$blacklist_ratio <- seuratobj$blacklist_region_fragments / seuratobj$peak_region_fragments
  return(seuratobj)
}
