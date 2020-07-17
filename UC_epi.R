library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
source("~/Msc project/data/seurat_qc.R")

#load the dataset and make seurat object
uc.data <- Read10X(data.dir = "~/ds_group/Charlie Msc project data/uc/raw_epi",gene.column = 1)
meta.counts <- read.table(file=paste0("~/Msc project/data/uc/all.meta2.txt"),header = T,row.names=1,sep="\t")
uc <- CreateSeuratObject(counts = uc.data, meta.data = meta.counts,project = "UC_Colon_epi",min.cells = 3, min.features = 200)
saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/uc_epi_raw.rds") )
remove(uc.data,meta.counts)


uc[["Organ"]]<-"colon"

uc <- seurat_qc1(uc)
uc <- seurat_qc2(uc)

table(uc$Sample)
uc1 <- subset(uc, subset = Sample=="N10.EpiB")
uc2 <- subset(uc, subset = Sample=="N58.LPA2")
uc3 <- subset(uc, subset = Sample== "N15.EpiB")
uc4 <- subset(uc, subset = Sample=="N9.EpiA")
uc5 <- subset(uc, subset = Sample=="N24.EpiA")
uc6 <- subset(uc, subset = Sample=="N51.EpiA")
uc7 <- subset(uc, subset = Sample=="N51.EpiB")
Idents(uc)<- "Sample"
uc8 <- subset(uc, idents = c("N10.EpiB","N58.LPA2","N15.EpiB",
                             "N9.EpiA","N24.EpiA","N51.EpiA","N51.EpiB"),invert = TRUE)

nExp_poi <- round(0.015*length(uc1$orig.ident))
uc1 <- doubletFinder_v3(uc1, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc2$orig.ident))
uc2 <- doubletFinder_v3(uc2, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc3$orig.ident))
uc3 <- doubletFinder_v3(uc3, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc4$orig.ident))
uc4 <- doubletFinder_v3(uc4, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc5$orig.ident))
uc5 <- doubletFinder_v3(uc5, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc6$orig.ident))
uc6 <- doubletFinder_v3(uc6, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.015*length(uc7$orig.ident))
uc7 <- doubletFinder_v3(uc7, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(uc8$orig.ident))
uc8 <- doubletFinder_v3(uc8, PCs = 1:23, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


uc1[["doublet"]]<-uc1$DF.classifications_0.25_0.02_37
uc2[["doublet"]]<-uc2$DF.classifications_0.25_0.02_45
uc3[["doublet"]]<-uc3$DF.classifications_0.25_0.02_38
uc4[["doublet"]]<-uc4$DF.classifications_0.25_0.02_48
uc5[["doublet"]]<-uc5$DF.classifications_0.25_0.02_54
uc6[["doublet"]]<-uc6$DF.classifications_0.25_0.02_41
uc7[["doublet"]]<-uc7$DF.classifications_0.25_0.02_49
uc8[["doublet"]]<-uc8$DF.classifications_0.25_0.02_598

uc <- merge(uc1, y = c(uc2,uc3,uc4,uc5,uc6,uc7,uc8), project = "UC_Colon_noninflammed")

uc <- seurat_qc2(uc)
png("ump_doublet.png")
show(DimPlot(uc, reduction = "umap",group.by = "doublet"))
dev.off()

write.csv(table(uc$doublet),"doublet.txt")
write.csv(table(uc$Health),"health.txt")
uc<- subset(uc,subset = doublet =="Singlet")
remove(uc1,uc2,uc4,uc3,uc5,uc6,uc7,uc8)
uc$pANN_0.25_0.02_37<- NULL
uc$pANN_0.25_0.02_45 <- NULL
uc$pANN_0.25_0.02_38<- NULL
uc$pANN_0.25_0.02_48<- NULL
uc$pANN_0.25_0.02_54<- NULL
uc$pANN_0.25_0.02_41<- NULL
uc$pANN_0.25_0.02_49<- NULL
uc$pANN_0.25_0.02_598<- NULL
uc$DF.classifications_0.25_0.02_37<- NULL
uc$DF.classifications_0.25_0.02_45<- NULL
uc$DF.classifications_0.25_0.02_38<- NULL
uc$DF.classifications_0.25_0.02_48<- NULL
uc$DF.classifications_0.25_0.02_54<- NULL
uc$DF.classifications_0.25_0.02_41<- NULL
uc$DF.classifications_0.25_0.02_49<- NULL
uc$DF.classifications_0.25_0.02_598<- NULL
Idents(uc)<- "Health"
uc_n <- subset(uc, idents = c("Non-inflamed"))
uc_i <- subset(uc, idents = c("Inflamed"))
uc_h <- subset(uc, idents = c("Healthy"))

png("umap_sample.png",width = 780,height = 480)
show(DimPlot(uc, reduction = "umap",group.by = "Sample",pt.size = 0.5,repel = TRUE))
dev.off()
png("umap_health.png",width = 780,height = 480)
show(DimPlot(uc, reduction = "umap",group.by = "Health",pt.size = 0.5,repel = TRUE))
dev.off()
png("umap_celltype.png",width = 780,height = 480)
show(DimPlot(uc_n, reduction = "umap",group.by = "Cluster",label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()

png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(uc, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/epi obj/uc_epi_pre_integration.rds") )
