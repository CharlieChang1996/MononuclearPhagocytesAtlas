library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
source('~/Msc project/data/seurat_qc.R')

colon[["Health"]] <- "Healthy"
ibdhs[["Health"]] <- "Inflammed IBD"
ibdpi[["Health"]] <- "Inflammed IBD"
ibdpn[["Health"]] <- "Non-Inflammed IBD"
ibdhb[["Health"]] <- "Inflammed IBD"
uc_i[["Health"]] <- "Inflammed IBD"
uc_n[["Health"]] <- "Non-Inflammed IBD"
gut.list <- list(colon,ibdhs,ibdpi,ibdpn,uc_i,uc_n)
gut.anchors <- FindIntegrationAnchors(object.list = gut.list, dims = 1:30)
gut.integrated <- IntegrateData(anchorset = gut.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(gut.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
gut.integrated <- ScaleData(gut.integrated, verbose = FALSE)
gut.integrated <- RunPCA(gut.integrated, npcs = 30, verbose = FALSE)
gut.integrated <- RunUMAP(gut.integrated, reduction = "pca", dims = 1:30)


png("umap_anchor.png")
show(DimPlot(gut.integrated, reduction = "umap", group.by = "Health"))
dev.off()
png("umap_anchor_study.png")
show(DimPlot(gut.integrated, reduction = "umap", group.by = "dataname"))
dev.off()
png("anchor_gut1_sample.png",width = 780)
show(DimPlot(object = gut.integrated, reduction = "umap", group.by = "Sample"))
dev.off()

write.csv(table(gut.integrated$dataname),"study.txt")
write.csv(table(gut.integrated$Health),"health.txt")

png("anchor_ibdh.png")
DimPlot(object = subset(gut.integrated, subset = dataname =="IBD_Huang"), reduction = "umap",group.by = "Sample")
dev.off()
png("anchor_ibdps.png")
DimPlot(object = subset(gut.integrated, subset = dataname =="s"), reduction = "umap",group.by = "Sample")
dev.off()
png("anchor_uc.png",width = 780)
DimPlot(object = subset(gut.integrated, subset = dataname =="UC_gut"), reduction = "umap",group.by = "Sample")
dev.off()

##Annotation
DefaultAssay(gut.integrated) <- "RNA"


png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(gut.integrated, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(gut.integrated, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(gut.integrated, features = c("EPCAM", "KRT19")))
dev.off()
#SCINA
DefaultAssay(gut.integrated) <- "integrated"
exp <- GetAssayData(gut.integrated, slot = "scale.data")
signatures <- preprocess.signatures("~/Msc project/data/uc/uc_imm_marker.CSV")
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, convergence_rate = 0.999, 
                sensitivity_cutoff = 0.5, rm_overlap=FALSE, allow_unknown=TRUE)
names(results$cell_labels) <- colnames(gut.integrated)
gut.integrated <- AddMetaData(gut.integrated, metadata = results$cell_labels,col.name = 'SCINA_re.type')
png("umap_SCINA_gut_anchor.png",width = 680)
show(DimPlot(gut.integrated, reduction = "umap",group.by = 'SCINA_re.type', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE)) 
dev.off()

saveRDS(gut.integrated,paste0("~/ds_group/Charlie Msc project data/downstream/colon_ibdhs.rds") )

