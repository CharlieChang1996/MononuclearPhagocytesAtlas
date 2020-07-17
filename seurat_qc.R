#Seurat QC & Pre-processing
library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)

seurat_qc1 <- function(obj){
  obj[["dataname"]] <- obj@project.name
  ## QC ##
  #QC using mitchondrial metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  #Violin plots
  Idents(obj) <- "dataname"
  violin1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png("violin1.png")
  show(violin1)
  dev.off()
  #Scatter plots
  s.plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  s.plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  png("scatter.png")
  show(s.plot1 + s.plot2)
  dev.off()
  #Filter data
  nfl <- readline(prompt = "Enter lower limit for nFeature")
  nfh <- readline(prompt = "Enter higher cutoff for nFeature")
  mtl <- readline(prompt = "Enter lower limit for percent.mt")
  mth <- readline(prompt = "Enter higher limit for percent.mt")
  nfl <- as.numeric(nfl)
  nfh <- as.numeric(nfh)
  mtl <- as.numeric(mtl)
  mth <- as.numeric(mth)
  obj <- subset(obj, subset = nFeature_RNA > nfl & nFeature_RNA < nfh & 
                  percent.mt > mtl & percent.mt < mth)
  
  return(obj)
}




seurat_qc2 <- function(obj){
  obj[["dataname"]] <- obj@project.name
  ## Normalisation & Feature selection##
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  #ten most variable features & plot them
  top10 <- head(VariableFeatures(obj), 10)
  Fplot <- VariableFeaturePlot(obj)
  Fplot <- LabelPoints(plot = Fplot, points = top10, repel = TRUE)
  png("feature.png")
  show(Fplot)
  dev.off()
  ## Scale the data##
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  
  ## PCA ##
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  #plot results
  VizDimLoadings(obj, dims = 1:2, reduction = "pca")
  png("pcaDim.png")
  show(DimPlot(obj, reduction = "pca"))
  dev.off()
  png("elbow.png")
  show(ElbowPlot(obj,ndims = 50 ))
  dev.off()
  ndim <- readline(prompt = "How many dims for umap/tsne?")
  ndim <- as.integer(ndim)
  obj <- RunUMAP(obj, dims = 1: ndim)
  
  png("umap.png",height = 480,width = 680)
  show(DimPlot(obj, reduction = "umap",group.by = "orig.ident"))
  dev.off()
  return(obj)
}


doublet_rem<- function(obj){
  ##Doublet Removal##
  #pK Identification (no ground-truth)
  ndim <- readline(prompt = "How many dims for doublet removal?")
  ndim <- as.integer(ndim)
  sweep.res.list_obj <- paramSweep_v3(obj, PCs = 1:ndim, sct = FALSE)
  sweep.stats_obj <- summarizeSweep(sweep.res.list_obj, GT = FALSE)
  bcmvn_obj <- find.pK(sweep.stats_obj)
  write.table(bcmvn_obj,"pk.txt")
  drate <- readline(prompt = "Enter doublet formation rate?")
  drate <- as.numeric(drate)
  nExp_poi <- round(drate*length(obj$orig.ident))  ## Assuming doublet formation rate 
  #Find doublets & visualisation
  npk <-readline(prompt = "Enter pK")
  npk <- as.numeric(npk)
  obj <- doubletFinder_v3(obj, PCs = 1:ndim, pN = 0.25, pK = npk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  return(obj)
}


clus_anno <- function(obj){
  png("feature_marker.png",width = 900,height = 900)
  show(FeaturePlot(obj, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
  dev.off()
  png("feature_marker2.png",width = 900,height = 900)
  show(FeaturePlot(obj, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
  dev.off()
  png("feature_marker_epi.png",width = 780,height = 480)
  show(FeaturePlot(obj, features = c("EPCAM", "KRT19")))
  dev.off()
  #SCINA annotation
  exp <- GetAssayData(obj, slot = "scale.data")
  sign_route <- readline(prompt = "Enter route of the marker file")
  signatures <- preprocess.signatures(sign_route)
  results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
                  sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE)
  names(results$cell_labels) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = results$cell_labels,col.name = 'SCINA.type')
  png("umap_SCINA.png",width = 680,height = 680)
  show(DimPlot(obj, reduction = "umap",group.by = 'SCINA.type', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE) )
  dev.off()
  return(obj)
}



#saveRDS(obj,"qc_seurat_object.rds")