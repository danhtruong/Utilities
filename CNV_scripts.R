library(data.table)
require(ComplexHeatmap)
require(circlize)


# Read peaks
readSeg <- function(file.path, chromCoords){
  seg <- data.frame(fread(file.path, stringsAsFactors = FALSE, sep="\t"))
  seg <- data.frame(Sample=seg[,1],
                    chromosome=gsub("chr", "", seg[,2]),
                    start=as.numeric(seg[,3]),
                    end=as.numeric(seg[,4]),
                    n_probes=seg[,5],
                    segment_mean=seg[,6])
  seg$start <- as.numeric(lapply(1:nrow(seg),
                                 function(i) subset(chromCoords, chr==seg$chr[i])$start + seg$start[i]))
  seg$end <- as.numeric(lapply(1:nrow(seg),
                               function(i) subset(chromCoords, chr==seg$chr[i])$start + seg$end[i]))
  return(seg[complete.cases(seg),])
}



# Transfer the values from data frame into the matrix
populateCNVMatrix <- function(chromCoords, coordMat, segDf, correction =501225){
  coordMat <- subset(coordMat, chr %in% segDf$chromosome)
  chromCoords <- subset(chromCoords, chr %in% segDf$chromosome)
  samples <- unique(segDf$Sample)
  cnvMat <- matrix(rep(rep(0,nrow(coordMat)),length(samples)),ncol=length(samples))
  colnames(cnvMat) <- samples
  populate.binned.data <- function(start, end, column, value){
    cnvMat[(start/correction):(end/correction), column] <<- value
  }
  for(i in samples) { apply(subset(segDf, Sample == i)[,c(3, 4, 6)], 1 , function(j) populate.binned.data(j[1],j[2],i,j[3]))}
  row.labels <- unlist(lapply(chromCoords$chr, function(i) rep(i, nrow(subset(coordMat, chr==i)))))
  row.names(cnvMat) <- row.labels
  return(cnvMat)
}

#Import segment
import_cnv <- function(file.path = 'OS1-OS31.seg', coords = chromCoordsHg19){
  seg <- readSeg(file.path, coords)
  coordMat <- do.call(rbind, lapply(coords$chr,
                                    function(i) {chr.i <- subset(coords, chr==i);
                                    return(data.frame(row=seq(chr.i$start, chr.i$cumsum, by=501225), chr=i))}))
  cnvMat <- populateCNVMatrix(chromCoords = coords,
                              coordMat = coordMat,
                              segDf = seg)
  coordMat <- subset(coordMat, chr %in% row.names(cnvMat))
  return(list = list('cnvMat' = cnvMat, 'coordMat' = coordMat))
}

cnv_plot <- function(cnvMat, coordMat, rownames = TRUE, title = '', cluster_rows = TRUE, borders = TRUE){
  coordMat$labels <- paste0('chr',coordMat$chr)

  colors = list( chr = setNames(rep(c('black', 'white'), 11), unique(coordMat$labels)))
  chr_annotation_1 = HeatmapAnnotation(chr = coordMat$labels, col = colors, border = TRUE,  show_legend = FALSE, annotation_label = ' ')

  x.center=mean(cnvMat)
  quantiles = quantile(cnvMat[cnvMat != x.center], c(0.01, 0.99))
  # determine max distance from the center.
  delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
  low_threshold = x.center - delta
  high_threshold = x.center + delta
  x.range = c(low_threshold, high_threshold)

  #ha = rowAnnotation(foo = anno_text(row_title, location = 0.25, rot = 0,
                                     #just = "center", gp = gpar(fontsize =9)))

  heatmap <- Heatmap(t(cnvMat), col = colorRamp2(c(x.range[1], x.center , x.range[2]), c("darkblue", "white", "darkred"), transparency = 0, space = "LAB"),
                     column_split = factor(coordMat$labels, levels = unique(coordMat$labels)),
                     #row_split = colnames(cnvMat),
                     top_annotation = chr_annotation_1, #left_annotation  = ha,
                     border = borders,
                     show_row_dend = FALSE,
                     row_title = title,
                     column_title_rot = 90,
                     column_gap = unit(0, 'mm'),
                     row_gap = unit(0, 'mm'),
                     show_column_names = FALSE,
                     column_title_gp = gpar(fontsize = 9),
                     row_title_gp = gpar(fontsize = 9),
                     cluster_columns = FALSE,
                     cluster_rows = cluster_rows,
                     show_row_names = rownames,
                     heatmap_legend_param = list(legend_direction = "horizontal",
                                                 legend_width = unit(2, "cm"),
                                                 title_position='topcenter',
                                                 title = "Expression", border = "black"),
                         #height  = unit(0.25 * dim(cnvMat)[1] , "cm")
                     , show_heatmap_legend = FALSE )
  return(heatmap)

}

load_inferCNV <- function(infercnv_obj, sample_size = 1000, coords = chromCoordsHg38){
  max_sample  <- dim(infercnv_obj@expr.data)[2]
  sample = sample(max_sample,sample_size)
  seg_fil <- cbind(infercnv_obj@expr.data[,sample], infercnv_obj@gene_order)

  seg_fil <- pivot_longer(seg_fil, cols = colnames(infercnv_obj@expr.data[,sample]))
  seg_fil$chr <- factor(seg_fil$chr, unique(seg_fil$chr))
  seg_fil$chr <- sapply(strsplit(as.character(seg_fil$chr) , split = 'chr'), '[[', 2)
  colnames(seg_fil) <- c("chromosome" ,  "start", "stop" , "Sample" , "value")

  seg_fil$start <- as.numeric(seg_fil$start)
  seg_fil$stop <- as.numeric(seg_fil$stop)
  for (i in 1:22){
    seg_fil[seg_fil$chromosome == i,]$start = seg_fil[seg_fil$chromosome == i,]$start + coords$start[coords$chr == i] - 1
    seg_fil[seg_fil$chromosome == i,]$stop = seg_fil[seg_fil$chromosome == i,]$stop + coords$start[coords$chr == i] - 1
  }
  return(list('seg_fil' = seg_fil, 'sample' = sample ))

}

inferCNV_mat <- function(seg_fil, infercnv_obj, bin = 500000, coords = chromCoordsHg38){
  coordMat <- do.call(rbind, lapply(coords$chr,
                                    function(i) {chr.i <- subset(coords, chr==i);
                                    return(data.frame(row=seq(chr.i$start, chr.i$cumsum, by=bin), chr=i))}))
  coordMat <- subset(coordMat, chr %in% seg_fil$chromosome)
  chromCoords <- subset(coords, chr %in% seg_fil$chromosome)

  samples <- unique(seg_fil$Sample)
  x.center=mean(infercnv_obj@expr.data)
  cnvMat <- matrix(rep(rep(x.center,nrow(coordMat)),length(samples)),ncol=length(samples))
  colnames(cnvMat) <- samples
  populate.binned.data <- function(start, end, column, value){
    cnvMat[(start/bin):(end/bin), column] <<- value
  }


  for(i in samples) {
    apply(subset(seg_fil, Sample == i)[,c(2, 3, 5)], 1 , function(j) populate.binned.data(j[1],j[2],i,j[3]))
  }
  row.labels <- unlist(lapply(chromCoords$chr, function(i) rep(i, nrow(subset(coordMat, chr==i)))))
  row.names(cnvMat) <- row.labels
  return(list('cnvMat' = cnvMat, 'coordMat' = coordMat))
}
