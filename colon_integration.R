library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)


# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
uc_epi <- subset(uc_epi, subset = Health == "Healthy")

uc_fib <- FindVariableFeatures(uc_fib, selection.method = "vst", nfeatures = 4000)
uc_epi <- FindVariableFeatures(uc_epi, selection.method = "vst", nfeatures = 4000)
## option 1: only consider genes that are HVGs in at least 2 studies
#colon.HVGs <- c()
#colon.HVGs <- Reduce(append, list(VariableFeatures(col), VariableFeatures(ibdp), VariableFeatures(uc)))
#colon.HVGs <- unique(colon.HVGs[which(table(colon.HVGs) > 2)]) # 2166 HVGs
#colon <- merge(uc, y = c(col,ibdp), project = "colon")
 #ibdp <- subset(colon, subset = dataname =="IBD_Parikh")
 #col <- subset(colon, subset = dataname =="IBD_Huang")
 #uc <- subset(colon, subset = dataname =="UC_Colon")

colon.list <- list(colon,uc_fib,uc_epi)
colon.anchors <- FindIntegrationAnchors(object.list = colon.list, dims = 1:30)
colon.integrated <- IntegrateData(anchorset = colon.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(colon.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
colon.integrated <- ScaleData(colon.integrated, verbose = FALSE)
colon.integrated <- RunPCA(colon.integrated, npcs = 30, verbose = FALSE)
colon.integrated <- RunUMAP(colon.integrated, reduction = "pca", dims = 1:30)
png("umap_anchor.png")
show(DimPlot(colon.integrated, reduction = "umap", group.by = "dataname"))
dev.off()
png("anchor_colon1_sample.png",width = 780)
show(DimPlot(object = colon.integrated, reduction = "umap", group.by = "Sample"))
dev.off()
write.csv(table(colon.integrated$dataname),"study.txt")
write.csv(table(colon.integrated$Sample),"sample.txt")
png("anchor_ibdh.png")
DimPlot(object = subset(colon.integrated, subset = dataname =="IBD_HUANG"), reduction = "umap",group.by = "Sample")
dev.off()
png("anchor_ibdp.png")
DimPlot(object = subset(colon.integrated, subset = dataname =="IBD_PARIKH"), reduction = "umap",group.by = "Sample")
dev.off()
png("anchor_uc.png",width = 780)
DimPlot(object = subset(colon.integrated, subset = dataname =="UC_COLON"), reduction = "umap",group.by = "Sample")
dev.off()

##Annotation
DefaultAssay(colon.integrated) <- "RNA"


png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(colon.integrated, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(colon.integrated, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(colon.integrated, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(colon.integrated,paste0("~/Msc project/data/colon/colon_anchor7_15.rds") )
remove(colon,uc_epi,uc_fib)
#SCINA
DefaultAssay(colon.integrated) <- "integrated"
exp <- GetAssayData(colon.integrated, slot = "scale.data")
signatures <- preprocess.signatures("~/Msc project/data/uc/uc_imm_marker.CSV")
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
                sensitivity_cutoff = 0.5, rm_overlap=FALSE, allow_unknown=TRUE)
names(results$cell_labels) <- colnames(colon.integrated)
colon.integrated <- AddMetaData(colon.integrated, metadata = results$cell_labels,col.name = 'SCINA_re.type')
png("umap_SCINA_colon_anchor.png",width = 680)
show(DimPlot(colon.integrated, reduction = "umap",group.by = 'SCINA_re.type', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE)) 
dev.off()

saveRDS(colon.integrated,paste0("~/Msc project/data/colon/colon_anchor7_15.rds") )
