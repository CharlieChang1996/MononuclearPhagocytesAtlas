#Seurat QC & Pre-processing
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(umap)
#library(DoubletFinder)
library(SCINA)
library(readxl)
library(FSA)
library(ggplot2)
library(RColorBrewer)
library(projectR)
library(CoGAPS)
library(ggsignif)
library(magrittr)
library(monocle3)
library(topGO)
library(enrichplot)
library(DO.db)
library(clusterProfiler)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(limma)
library(pheatmap)

seurat_qc1 <- function(obj,fname){
  obj[["dataname"]] <- obj@project.name
  ## QC ##
  #QC using mitchondrial metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  #Violin plots
  Idents(obj) <- "dataname"
  violin1 <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  pdf(paste("violin1_",fname,".pdf"))
  show(violin1)
  dev.off()
  #Scatter plots
  s.plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  s.plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  pdf(paste("scatter_",fname,".pdf"))
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
  ## Normalisation & Feature selection##
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000) 
  return(obj)
}




seurat_qc2 <- function(obj,fname){
  
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  #ten most variable features & plot them
  top10 <- head(VariableFeatures(obj), 10)
  Fplot <- VariableFeaturePlot(obj)
  Fplot <- LabelPoints(plot = Fplot, points = top10, repel = TRUE)
  pdf(paste("feature_",fname,".pdf"))
  show(Fplot)
  dev.off()
  ## Scale the data##
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  
  ## PCA ##
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  #plot results
  VizDimLoadings(obj, dims = 1:2, reduction = "pca")
  pdf(paste("pcaDim_",fname,".pdf"))
  show(DimPlot(obj, reduction = "pca"))
  dev.off()
  pdf(paste("elbow",fname,".pdf"))
  show(ElbowPlot(obj,ndims = 50 ))
  dev.off()
  ndim <- readline(prompt = "How many dims for umap/tsne?")
  ndim <- as.integer(ndim)
  obj <- RunUMAP(obj, dims = 1: ndim)
  
  pdf(paste("umap",fname,".pdf"))
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
  
  #SCINA annotation
  exp <- GetAssayData(obj, slot = "scale.data")
  sign_route <- readline(prompt = "Enter route of the marker file")
  signatures <- preprocess.signatures(sign_route)
  results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
                  sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE)
  names(results$cell_labels) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = results$cell_labels,col.name = 'SCINA.type')
  pdf("umap_SCINA.pdf",width = 680,height = 680)
  show(DimPlot(obj, reduction = "umap",group.by = 'SCINA.type', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE) )
  dev.off()
  return(obj)
}

seurat_clus <- function(obj,fname){
  ## Clustering ##
  ndims <-readline(prompt = "Enter dims for clustering")
  ndims <- as.integer(ndims)
  rsl <-readline(prompt = "Enter resolution for clustering")
  rsl <- as.numeric(rsl)
  obj <- FindNeighbors(obj, dims = 1:ndims)
  obj <- FindClusters(obj, resolution = rsl)
  pdf(paste("cluster_umap",fname,".pdf"))
  show(DimPlot(obj, reduction = "umap", label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
  dev.off()
  #Feature Plot
  DefaultAssay(obj) <- "RNA"
  pdf(paste("feature_IMM_",fname,".pdf"))
  show(FeaturePlot(obj, features = c("CD79A", "CD79B","CD3D", "CD8A", "NKG7","LYZ","C1QA","FCN1","CD68")))
  dev.off()
  pdf(paste("feature_other_",fname,".pdf"))
  show(FeaturePlot(obj, features = c("TPSAB1","CPA3","MZB1","IGHA1","COL1A1","EPCAM", "ACTA2","RAMP2","EMCN")))
  dev.off()
  return(obj)
}

seurat_inte <- function(obj.list){
  ndims <-readline(prompt = "Enter dims for Anchor CCA integration")
  ndims <- as.integer(ndims)
  obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:ndims)
  obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:ndims)
  # switch to integrated assay. The variable features of this assay are automatically
  # set during IntegrateData
  DefaultAssay(obj.integrated) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
  obj.integrated <- RunPCA(obj.integrated, npcs = ndims, verbose = FALSE)
  obj.integrated <- RunUMAP(obj.integrated, reduction = "pca", dims = 1:ndims)
  return(obj.integrated)
}

sderr_log <- function(){
  #statistical analysis
  route <-readline(prompt = "Enter route of data: ")
  sheet <-readline(prompt = "Enter sheet: ")
  out <-readline(prompt = "Enter outputname: ")
  data <- read_excel(route, sheet = sheet)
  
  #compute standard error
  sm <- Summarize(
      Freq ~ Health,
            data=data,
            digits=3)
  
  sd1 <- sm$sd[2]
  sd2 <- sm$sd[1]
  #compute P value using Mann-Whitney-Wilcoxon Test 
  w <- wilcox.test(Freq ~ Health, data=data)
  w<- as.character(w)
  v <- log(sd1/sd2,2)
  write(sheet,out,append = TRUE)
  write(w,out,append = TRUE)
  write(v,out,append = TRUE,sep = "\t")
}
chi_log <- function(){
  log1 <- readline(prompt = "Enter log ratio 1: ")
  log2 <- readline(prompt = "Enter log ratio 2: ")
  log1 <- as.numeric(log1)
  log2 <- as.numeric(log2)
  sheet <-readline(prompt = "Enter cluster name: ")
  #compute P value using Mann-Whitney-Wilcoxon Test 
  ch <- chisq.test(c(log1, log2), p = c(1, 1)/2, correct = FALSE)
  ch <- as.character(ch)
 
  write(sheet,"chi-test.txt",append = TRUE)
  write(ch,"chi-test.txt",append = TRUE)
}

de_sub <- function(obj,ct){
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- obj$annotation1
  cto <- subset(obj,idents = ct)
  Idents(cto) <- cto$Health
  i_h.markers <- FindMarkers(cto, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
  n_h.markers <- FindMarkers(cto, ident.1 = "Non-Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
  i_n.markers <- FindMarkers(cto, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD")
  outname1<- readline(prompt = "Enter name for output ivh ")
  outname2<- readline(prompt = "Enter name for output nvh ")
  outname3<- readline(prompt = "Enter name for output ivn ")
  write.csv(i_h.markers,outname1)
  write.csv(n_h.markers,outname2)
  write.csv(i_n.markers,outname3)
}

cds_pre <- function(obj,fname,root){
  cds <- as.cell_data_set(obj)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  pdf(paste("trajectory_",fname,".pdf",sep = ""),width = 9)
  show(plot_cells(cds, label_groups_by_cluster = F, label_leaves = T, 
             label_branch_points = F,trajectory_graph_segment_size = 0.7,group_label_size = 3,cell_size = 0.6))
  dev.off()
  ##pseudotime
  max.avp <- which.max(unlist(FetchData(obj, root)))
  max.avp <- colnames(obj)[max.avp]
  cds <- order_cells(cds,root_cells = max.avp)
  pdf(paste("trajectory_",fname,"_pseudotime.pdf",sep = ""),width = 9)
  show(plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T, 
             label_branch_points = F,trajectory_graph_segment_size = 0.5,group_label_size = 3,cell_size = 1))
  dev.off()
  return(cds)
}
#saveRDS(obj,"qc_seurat_object.rds")