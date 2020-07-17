library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
uc <- readRDS("~/ds_group/Charlie Msc project data/uc/uc_norm_h2.rds")
ibdp <- readRDS("~/ds_group/Charlie Msc project data/IBDP/ibdp_pre_integration.rds")
col <- readRDS("~/ds_group/Charlie Msc project data/IBD_H/ibdh_col_pre_integration.rds")

# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
col <- FindVariableFeatures(col, selection.method = "vst", nfeatures = 4000)
ibdp <- FindVariableFeatures(ibdp, selection.method = "vst", nfeatures = 4000)
uc <- FindVariableFeatures(uc, selection.method = "vst", nfeatures = 4000)
## option 1: only consider genes that are HVGs in at least 2 studies
colon.HVGs <- c()
colon.HVGs <- Reduce(append, list(VariableFeatures(col), VariableFeatures(ibdp), VariableFeatures(uc)))
colon.HVGs <- unique(colon.HVGs[which(table(colon.HVGs) > 1)]) # 2166 HVGs
colon <- merge(uc, y = c(col,ibdp), project = "colon")
VariableFeatures(colon) <- colon.HVGs

colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(colon)
colon <- ScaleData(colon, features = all.genes)
colon <- RunPCA(colon, features = VariableFeatures(object = colon))
saveRDS(colon,paste0("~/Msc project/data/colon/colon_pre_harmony.rds") )


#Try more covariates
colon <- RunHarmony(colon, "Sample",plot_convergence = TRUE,theta = 3,lambda = 1)
png("elbow.png")
show(ElbowPlot(colon,ndims = 50))
dev.off()
colon <- RunUMAP(colon, dims = 1:20,reduction = "harmony")
png("harmony_colon1.png")
show(DimPlot(object = colon, reduction = "umap", group.by = "dataname"))
dev.off()
png("harmony_colon1_sample.png",width = 780)
show(DimPlot(object = colon, reduction = "umap", group.by = "Sample"))
dev.off()
write.csv(table(colon$dataname),"study.txt")
write.csv(table(colon$Sample),"sample.txt")
png("harmony_ibdh.png")
DimPlot(object = subset(colon, subset = dataname =="IBD_Huang"), reduction = "umap",group.by = "Sample")
dev.off()
png("harmony_ibdp.png")
DimPlot(object = subset(colon, subset = dataname =="IBD_Parikh"), reduction = "umap",group.by = "Sample")
dev.off()
png("harmony_uc.png",width = 780)
DimPlot(object = subset(colon, subset = dataname =="UC_Colon"), reduction = "umap",group.by = "Sample")
dev.off()
png("qc_heatmap.png",width = 780)
show(FeaturePlot(colon, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2))
dev.off()
png("qc_mt_heatmap.png")
show(FeaturePlot(colon, features =  "percent.mt", pt.size = 0.2))
dev.off()
## Clustering ##

colon <- FindNeighbors(colon, dims = 1:20)
colon <- FindClusters(colon, resolution = 1.2)

png("umap_cluster_r1.2.png",width = 680)
show(DimPlot(colon, reduction = "umap",group.by = 'seurat_clusters', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()
#Feature plot
png("BM_marker.png")
show(FeaturePlot(colon, features = c("CD79A", "CD79B","LYZ")))
dev.off()
png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(colon, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(colon, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(colon, features = c("EPCAM", "KRT19")))
dev.off()
##Annotation
#Manual
colon.markers <- FindAllMarkers(colon, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(colon.markers,file = "markers_colon.csv")
new.cluster.ids <- c("Follicular", "Plasma", "Activated Fos-hi T cells", "CD8+ IELs", "Tregs", 
                     "Plasma", "CD8+ LP", "DC2","Cycling T","Monocytes",
                     "GC","Plasma","Plasma","Epithelials","Epithelials",
                     "CD69+ Mast","Plasma","GC","B Cells","NKs","Plasma","DC",
                     "Plasma","ILCs","Plasma","Plasma","DC")
names(new.cluster.ids) <- levels(colon)
colon <- RenameIdents(colon,new.cluster.ids)
colon[["annotation1"]] <- Idents(colon)
png("umap_annotation1.png",width = 680)
show(DimPlot(colon, reduction = "umap",group.by = 'annotation1', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()
#SCINA
exp <- GetAssayData(colon, slot = "scale.data")
signatures <- preprocess.signatures("~/Msc project/data/uc/uc_imm_marker.CSV")
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
                sensitivity_cutoff = 0.9, rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')
names(results$cell_labels) <- colnames(colon)
colon <- AddMetaData(colon, metadata = results$cell_labels,col.name = 'SCINA2.type')
png("umap_SCINA_harmony.png",width = 680)
show(DimPlot(colon, reduction = "umap",group.by = 'SCINA2.type', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE)) 
dev.off()
saveRDS(colon,paste0("~/ds_group/Charlie Msc project data/Colon/colon_ibd_uc_harmony.rds") )
