#Seurat QC & Pre-processing
library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
#load the dataset and make seurat object
uc.data <- Read10X(data.dir = "~/Msc project/data/uc",gene.column = 1)
meta.counts <- read.table(file=paste0("~/Msc project/data/uc/all.meta2.txt"),header = T,row.names=1,sep="\t")
uc <- CreateSeuratObject(counts = uc.data, meta.data = meta.counts,project = "UC_Colon",min.cells = 3, min.features = 3)
saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/uc_raw.rds") )
remove(uc.data)
## QC ##
#QC using mitchondrial metrics
uc[["percent.mt"]] <- PercentageFeatureSet(uc, pattern = "^MT-")

#Violin plots
violin1 <- VlnPlot(uc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
png("violin1.png")
show(violin1)
dev.off()
#Scatter plots
s.plot1 <- FeatureScatter(uc, feature1 = "nCount_RNA", feature2 = "percent.mt")
s.plot2 <- FeatureScatter(uc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("scatter.png")
show(s.plot1 + s.plot2)
dev.off()

#Filter data extract healthy cells
uc <- subset(uc, subset =  percent.mt < 50)
Idents(uc)<- "Health"
uc <- subset(uc, idents = c("Healthy"))

## Normalisation & Feature selection##
uc <- NormalizeData(uc, normalization.method = "LogNormalize", scale.factor = 10000)
uc <- FindVariableFeatures(uc, selection.method = "vst", nfeatures = 2000)
#ten most variable features & plot them
top10 <- head(VariableFeatures(uc), 10)
Fplot <- VariableFeaturePlot(uc)
Fplot <- LabelPoints(plot = Fplot, points = top10, repel = TRUE)
png("feature.png")
show(Fplot)
dev.off()
## Scale the data##
all.genes <- rownames(uc)
uc <- ScaleData(uc, features = all.genes)

## PCA ##
uc <- RunPCA(uc, features = VariableFeatures(object = uc))
#plot results
VizDimLoadings(uc, dims = 1:2, reduction = "pca")
DimPlot(uc, reduction = "pca")

#Jackstraw take 4h, too long
#uc <- JackStraw(uc, num.replicate = 100) 
#uc <- ScoreJackStraw(uc, dims = 1:20)
#JackStrawPlot(uc, dims = 1:15)

#elbow plot
png("elbow.png")
show(ElbowPlot(uc,ndims = 50))
dev.off()

#umap/tsne
uc <- RunUMAP(uc, dims = 1:15)
png("ump_sample.png",height = 480,width = 680)
show(DimPlot(uc, reduction = "umap",group.by = "Sample"))
dev.off()

##Doublet Removal##
#pK Identification (no ground-truth)
sweep.res.list_uc <- paramSweep_v3(uc, PCs = 1:15, sct = FALSE)
sweep.stats_uc <- summarizeSweep(sweep.res.list_uc, GT = FALSE)
bcmvn_uc <- find.pK(sweep.stats_uc)
nExp_poi <- round(0.1*length(uc$orig.ident))  ## Assuming doublet formation rate 
#Find doublets & visualisation
uc <- doubletFinder_v3(uc, PCs = 1:15, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
png("ump_doublet.png")
show(DimPlot(uc, reduction = "umap",group.by = "DF.classifications_0.25_0.29_5136"))
dev.off()
uc_doub <- subset(uc, subset = Sample =="N51.LPB")
uc <- subset(uc, subset = DF.classifications_0.25_0.29_5136 =="Singlet" & Sample != "N51.LPB")
uc_doub <-subset(uc_doub, subset = DF.classifications_0.25_0.23_64 =="Singlet")
uc <- merge(uc,y= uc_doub)
##Annotation##
#Feature plot of marker genes
png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(uc, features = c("EPCAM", "KRT19")))
dev.off()
png("umap_celltype.png",width = 780,height = 480)
show(DimPlot(uc, reduction = "umap",group.by = "Cluster",label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()
#Add metadata for integration
uc[["dataname"]] <- uc@project.name

saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/uc_norm_h2.rds") )
 