queryResolution4Seurat <- function(seurat.obj, k = 10, reduction = 'umap', npc = 20,
                                   min_resl = 0.1, max_resl = 1.5, max_iter = 30, doPCA = F){
  max.dim = ifelse(reduction == 'harmony', npc, 2)
  if(doPCA) {
    seurat.obj <- RunPCA(seurat.obj, npcs = npc, verbose = F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, verbose = F, check_duplicates = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc, verbose = F)
  }


  seurat.obj <- FindNeighbors(seurat.obj, reduction = reduction, verbose = F, dims = 2:max.dim)
  tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, dims = 2:max.dim)@active.ident
  tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, dims = 2:max.dim)@active.ident




  len1 = length(levels(tmp.cluster1))
  len2 = length(levels(tmp.cluster2))

  k1 = k2 = 0
  while(len1 > k ){

    k1 = k1 + 1
    message('min_resl too large, trying to divided it by  2')
    min_resl = min_resl/2
    tmp.cluster1 <- FindClusters(seurat.obj, resolution = min_resl, dims = 2:max.dim)@active.ident
    len1 = length(levels(tmp.cluster1))
    if(k1 == 10) stop('Please specify a much smaller min_res')
  }

  while(len2 < k){
    k2 = k2 + 1
    message('max_resl too small, trying to multiply it by 2')
    max_resl = max_resl * 2
    tmp.cluster2 <- FindClusters(seurat.obj, resolution = max_resl, dims = 2:max.dim)@active.ident
    len2 = length(levels(tmp.cluster2))
    if(k2 == 10) stop('Please specify a much bigger max_res')
  }
  if(len1 == k) {
    return(min_resl)
  }

  if(len2 == k) {
    return(max_resl)
  }

  # repeat in other case

  i = 0
  repeat{
    i = i + 1
    resl0 = min_resl/2 + max_resl/2

    tmp.cluster <- FindClusters(seurat.obj, resolution = resl0, dims = 2:max.dim)@active.ident

    len = length(levels(tmp.cluster))
    if(len == k){
      return(resl0)
    }
    if(len < k){
      min_resl = resl0
      len1 = len
    }
    if(len > k){
      max_resl = resl0
      len2 = len
    }
    if(i == max_iter) break
  }
  return(resl0)
}

#Clustering function using seurat SNN (Seurat v2.3.4)
seuratSNN <- function(matSVD, dims.use = 1:50, ...){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  rownames(mat) <- colnames(matSVD)[1:3]
  obj <- Seurat::CreateSeuratObject(mat, project='scATAC', min.cells=0, min.features=0)
  obj[['pca']] <- Seurat::CreateDimReducObject(embeddings = matSVD)
  obj <- Seurat::FindNeighbors(object = obj, dims  = dims.use)
  obj <- Seurat::FindClusters(object = obj, print.output = TRUE, ...)
  clust <- obj@meta.data[,ncol(obj@meta.data)]
  paste0("Cluster",match(clust, unique(clust)))
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

#Seurat to CDS
seurat2cds <- function(seurat.obj, assay){
  require(monocle3)
  mtx <- seurat.obj@assays[[assay]]@counts
  rnames = data.table('region' = rownames(mtx))
  tmp = tidyr::separate(rnames, col = 'region', into = c('chr', 'start', 'end'))
  rnames = paste0(tmp$chr, ':', tmp$start, '-', tmp$end)
  rownames(mtx) = rnames
  frac_in_cell = 0.05
  mtx = 1 * (mtx > 0)
  rr = Matrix::rowMeans(mtx)
  mtx = mtx[rr >= frac_in_cell, ]
  expression_matrix <- mtx
  gene_annotations <- as.matrix(rownames(expression_matrix))
  colnames(gene_annotations) <- "gene_short_name"
  rownames(gene_annotations) <- rownames(expression_matrix)
  cell_metadata <- seurat.obj@meta.data

  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotations)
  return(cds)
}

#make go plot function
make_go_plots <- function(markers, genes, database){
  require(enrichR)
  require(ggrepel)
  if(is.numeric(genes)){
    genes_for_enrichr <- markers %>%
      #filter(avg_logFC >= log(2)) %>%
      group_by(cluster) %>%
      top_n(genes, avg_log2FC)
  }
  else{
    genes_for_enrichr <- markers %>%
      #filter(avg_logFC >= log(2)) %>%
      group_by(cluster)
  }
  genes_for_enrichr$gene <- as.character(genes_for_enrichr$gene)

  go_terms <- lapply(database, function(x) {
    go_list <- list()
    cluster_length <- length(levels(factor(genes_for_enrichr$cluster)))
    for (i in 1:cluster_length){
      cluster <- levels(factor(genes_for_enrichr$cluster))[i]
      df <- as.data.frame(genes_for_enrichr[genes_for_enrichr$cluster %in% cluster,])
      go_list[i] <-  enrichr(df$gene, x)
      go_list[[i]]$cluster  <- cluster
    }
    go_terms <- dplyr::bind_rows(go_list)
  }
  )
  names(go_terms) <- database
  return(go_terms)
}



col.low = "dodgerblue"
col.mid = "floralwhite"
col.high = "brown1"

colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))

.distinctColorPalette <-function(k) {
  set.seed(123)
  if(packageVersion("scales") >= '1.1.0'){
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  } else {
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=60:100)(2e3))))
  }
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

make_sankey <- function (pred.labels,
                         labels = MTL.integrated$anno,
                         font_size = 4,
                         geom_alluvium_w  =1/16,
                         geom_stratum_w=5/16){
  predictions_df <- as.data.frame(pred.labels$labels)
  rownames(predictions_df) <- pred.labels@rownames
  colnames(predictions_df) <- 'predicted.id'
  predictions_df[,1] <- as.character(predictions_df[,1])
  predictions_df[is.na(predictions_df)] <- 'Unassigned'
  predictions_df$anno <- labels
  predictions_df <- predictions_df %>% group_by(anno) %>% dplyr::count	(predicted.id)

  require(ggalluvial)


  p <- ggplot(as.data.frame(predictions_df),
              aes(y = n, axis1 = anno, axis2 = predicted.id, fill = predicted.id)) +
    geom_alluvium(width = geom_alluvium_w) +
    geom_stratum(width = geom_stratum_w, fill = "white", color = "black") +
    geom_text(stat = 'stratum', infer.label = TRUE, size = font_size) +
    scale_x_discrete(limits = c("anno", "predicted.id"), expand = c(0,0)) +
    #scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle("") + theme_void() + labs(fill= 'Predicted ID') +
    theme(legend.text =element_text(size = 8), legend.title = element_text(size = 10))
  return(p)
}

col <- colorRampPalette(c("#1E90FF", "#FFFFFF", "#CD2626"),
                        space = "Lab"
)(100)

require(circlize)
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("#1E90FF", "#FFFFFF", "#CD2626"))

fig.size <- function (height, width)
{
  options(repr.plot.height = height, repr.plot.width = width)
}



clean_mass_model_object <- function(cm) {
  cm$y <- c()
  #cm$model = c()
  cm$residuals <- c()
  cm$fitted.values <- c()
  cm$effects <- c()
  cm$qr$qr <- c()
  cm$linear.predictors <- c()
  cm$weights <- c()
  cm$prior.weights <- c()
  cm$data <- c()

  cm$family$variance <- c()
  cm$family$dev.resids <- c()
  cm$family$aic <- c()
  cm$family$validmu <- c()
  cm$family$simulate <- c()
  cm$family$control <- c()
  attr(cm$terms,".Environment") <- c()
  attr(cm$formula,".Environment") <- c()
  return(cm)
}

extract_model_status_helper <- function(model){
  if (class(model)[1] == "speedglm") {
    status_str <- ifelse(model$convergence, "OK", "FAIL")
    return (status_str)

  } else if (class(model)[1] == "negbin"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  } else if (class(model) == "zeroinfl"){
    status_str <- ifelse(model$converged, "OK", "FAIL")
    return (status_str)
  }else {
    return("FAIL")
  }
}

cbPalette <- function(n)
{
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  colorRampPalette(cbPalette)(n)
}


flat.cols <- function(n)
{
  fc = c("#34495e", #wet asphalt
         "#9b59b6",  #amythest
         "#3498db",  #peter river
         "#2ecc71",  # emerald
         #"#1abc9c",  #turquiose
         "#f1c40f",  # sunflower
         "#e67e22",   # carrot
         "#e74c3c")   # alizarin
  #"#ecf0f1",   #clouds
  #"#95a5a6")   # concrete
  return(colorRampPalette(fc)(n))
}

#themes
theme_dot_plot <- function() { theme_bw() + theme(axis.title.x = element_blank(),
                                                  axis.title.y = element_blank(),
                                                  axis.text.y = element_text(size = 7),
                                                  axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                                                  aspect.ratio = 2, legend.text = element_text(size = 5),
                                                  legend.title = element_text(size=6), legend.position="right",
                                                  legend.box = 'vertical', legend.box.just = 'left')
                                                  }

require(ggh4x)
box_plot_theme <- theme(panel.background = element_blank(), strip.background = element_blank(), axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"), legend.position = 'none', axis.title = element_blank(), strip.text = element_text(face = 'bold'), text = element_text(size = 8))

#theme for signature plots
signature_theme <- theme(panel.border = element_rect(fill = NA, color = 'black'), legend.position = 'none', panel.background = element_rect(fill= NA), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" ),strip.background = element_rect(fill='gray', colour = 'black'), legend.key=element_blank(), strip.text = element_text(size = 6))

#caclculate effect size
calc_effect_size <- function(df, signature_names)
{
  require(effectsize)
  res <- data.frame()
  for (i in 1:length(signature_names))
  {
    data <- df[df$name %in% signature_names[[i]],]
    m <- lm(value ~ sample_type, data = data)
    lm_res <- as.data.frame(parameters::model_parameters(standardize(m, two_sd = TRUE)))
    lm_res$Coefficient <- round(lm_res$Coefficient,2)
    lm_res <- dplyr::rename(lm_res, Predictor = Parameter)
    lm_res <- lm_res %>% select(Predictor, Coefficient, CI_low, CI_high, p)
    lm_res$Signature <- signature_names[[i]]
    res <- rbind(res, lm_res)

  }
  return(res)
}

#generate a plot based on effect size
effect_size_plot <- function(res){
  p <- ggplot(data= res, aes(y=rev(Predictor), x=Coefficient, xmin=CI_low, xmax=CI_high)) + geom_point() +
    facet_grid(Signature~., scales= "free", space="free") +
    geom_errorbarh(height=.1) +
    scale_x_continuous(limits=c(-3,2.5)) +
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
    theme_minimal() + theme(strip.text = element_text(size = 6)) +
    theme(panel.spacing = unit(1, "lines")) +
    ylab('') + scale_y_discrete(labels =switch(min(1,sum(res$Predictor %in% 'sample_typeCold')) + 1, c('Nuclei', 'Warm'),c('Nuclei', 'Cold', 'Warm') ))

  anno <-  ggplot(data= res, aes(y=rev(Predictor), x=Coefficient, xmin=CI_low, xmax=CI_high)) +
    facet_grid(Signature~., scales= "free", space="free") +
    theme(text = element_blank(), rect = element_blank(), axis.ticks = element_blank()) + xlab(' ') +
    scale_x_continuous(limits=c(0.5,1.5)) +
    geom_text(data = res,
              mapping = aes(x = 0.75, label = paste0(round(CI_low,2), " \u2013 " , round(CI_high,2)),), nudge_y = 0, size = 3) +
    geom_text(data = res,
              mapping = aes(x = 1.25, label = scientific(p, digits = 2)), nudge_y = 0, size = 3)
  plot <-  plot_grid(ggdraw(),ggdraw() ,p, anno, nrow =2, align = 'hv', axis = 'bt', rel_heights = c(1,10)) + annotate('text',x = c(0.65,0.85) , y = 0.95, label = c('CI', 'italic(p)'), parse = TRUE)
  return(plot)
}

make_vp <- function(avg.expression_dissociation, degs){
  require(ggrepel)
  x = degs
  df <- avg.expression_dissociation
  colnames(df) <- c('X','Y')
  df$SYMBOL <- rownames(df)
  clusters = levels(degs$cluster)
  top10 <- x %>% subset(p_val_adj < 0.05 & cluster == clusters[1]) %>% top_n(10, avg_log2FC)
  bot10 <- x %>% subset(p_val_adj < 0.05 & cluster == clusters[2]) %>% top_n(10, avg_log2FC)

  UP <- x %>% subset(p_val_adj < 0.05 & cluster == clusters[1] & avg_log2FC > log2(1.5)) %>% pull(gene)
  DOWN <- x %>% subset(p_val_adj < 0.05 & cluster == clusters[2] & avg_log2FC > log2(1.5))  %>% pull(gene)

  df$DE <- ''
  df[UP,'DE'] <- 'UP'
  df[DOWN, 'DE'] <- 'DOWN'
  df$DE[df$DE == ''] <- 'NS'
  vp_cols <- setNames(c("black", "darkblue", "darkred"), c('NS', 'DOWN', 'UP'))

  df$genelabels <- ""
  df[rownames(rbind(top10,bot10)),'genelabels'] <- TRUE
  gene_ab_plots <- ggplot(df, aes(x= X , y = Y, color = DE)) + geom_point(size = 0.5) +
    geom_label_repel(aes(label = ifelse(genelabels == T, SYMBOL,"")), size = 1.5, max.overlaps = Inf, force = 2,segment.alpha = 4/5, color = 'black') +
    scale_color_manual(values=vp_cols ) + xlab(paste(clusters[2],'avg_log2')) + ylab(paste(clusters[1],'avg_log2')) +
    theme(legend.position = "none", text = element_text(size = 12),
          panel.background=element_blank(),
          panel.border = element_rect(fill = NA))
  return(gene_ab_plots)
}

grid_plot <-  function(x){x + theme(axis.text = element_blank(), axis.title = element_blank(),plot.title = element_text(size=10) , plot.subtitle = element_text(size=8, hjust = 0.5), axis.ticks = element_blank(), panel.border = element_rect(colour = 'black'), plot.margin = unit(c(0, 0, 0, 0), "cm")) + NoLegend()}
