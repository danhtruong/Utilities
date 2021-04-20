
path_input <- function(path){
  list_of_files <- list.dirs(path, recursive = FALSE)
  list_of_files <- data.frame('path' = list_of_files)
  list_of_files$path <- sapply(list_of_files$path, function(x) { paste0(x, '/outs')})
  list_of_files$h5 <- sapply(list_of_files$path, function(x) { paste0(x, '/filtered_feature_bc_matrix.h5')})
  list_of_files$library_id <- sapply(strsplit(list_of_files$path, split="/"), "[[", 3)
  list_of_files$sample_type <- sapply(strsplit(list_of_files$library_id, split="-"), "[[", 4)
  list_of_files$metrics_path <- paste0(list_of_files$path,'/metrics_summary.csv')
  return(list_of_files)
}

qc_function <- function(seurat_data, group){
  seurat_data@active.ident <- factor(seurat_data$group)
  seurat_data.sce <- as.SingleCellExperiment(seurat_data)

  discard.mito <- isOutlier(seurat_data.sce$percent.mt,
                            type="higher", batch=seurat_data.sce$group, nmads = 2)
  #mito.plot <- plotColData(seurat_data.sce, x="group", y="percent.mt",
                          # colour_by=I(discard.mito)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + labs(fill = 'Outlier') + ylab('% Mitochondrial Mapping')

  discard.counts <- isOutlier(seurat_data.sce$nCount_RNA, nmads = 2, log = TRUE,
                              type="lower", batch=seurat_data.sce$group)
  #counts.plot <- plotColData(seurat_data.sce, x="group", y="nCount_RNA",
                             #colour_by=I(discard.counts)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Counts')

  discard.features <- isOutlier(seurat_data.sce$nFeature_RNA, nmads = 2, log = TRUE,
                                type="lower", batch=seurat_data.sce$group)
  #gene.plot <- plotColData(seurat_data.sce, x="group", y="nFeature_RNA",
                           #colour_by=I(discard.features)) + theme(axis.text.x = element_text(angle =90)) + xlab('') + ylab('Genes')
  seurat_data <- AddMetaData(seurat_data, discard.counts, col.name = 'discard.counts')
  seurat_data <- AddMetaData(seurat_data, discard.mito, col.name = 'discard.mito')
  seurat_data <- AddMetaData(seurat_data, discard.features, col.name = 'discard.features')
  return(seurat_data)
}

