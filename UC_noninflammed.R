library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
source("~/Msc project/data/seurat_qc.R")
uc[["Organ"]]<-"colon"
Idents(uc)<- "Health"
uc <- subset(uc, idents = c("Non-inflamed"))

uc <- seurat_qc1(uc)
uc <- seurat_qc2(uc)

table(uc$Sample)
uc1 <- subset(uc, subset = Sample=="N49.LPA")
uc2 <- subset(uc, subset = Sample=="N50.LPA")
uc3 <- subset(uc, subset = Sample== "N111.LPA2")
uc4 <- subset(uc, subset = Sample=="N661.LPA1")
uc5 <- subset(uc, subset = Sample=="N661.LPA2")
uc6 <- subset(uc, subset = Sample=="N44.LPA")
Idents(uc)<- "Sample"
uc7 <- subset(uc, idents = c("N111.LPA2","N50.LPA","N661.LPA1",
                             "N661.LPA2","N49.LPA","N44.LPA"),invert = TRUE)

nExp_poi <- round(0.016*length(uc1$orig.ident))
uc1 <- doubletFinder_v3(uc1, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.016*length(uc2$orig.ident))
uc2 <- doubletFinder_v3(uc2, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.016*length(uc3$orig.ident))
uc3 <- doubletFinder_v3(uc3, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.05*length(uc4$orig.ident))
uc4 <- doubletFinder_v3(uc4, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.05*length(uc5$orig.ident))
uc5 <- doubletFinder_v3(uc5, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.03*length(uc6$orig.ident))
uc6 <- doubletFinder_v3(uc6, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(uc7$orig.ident))
uc7 <- doubletFinder_v3(uc7, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

uc1[["doublet"]]<-uc1$DF.classifications_0.25_0.29_56
uc2[["doublet"]]<-uc2$DF.classifications_0.25_0.29_61
uc3[["doublet"]]<-uc3$DF.classifications_0.25_0.29_46
uc4[["doublet"]]<-uc4$DF.classifications_0.25_0.29_538
uc5[["doublet"]]<-uc5$DF.classifications_0.25_0.29_542
uc6[["doublet"]]<-uc6$DF.classifications_0.25_0.29_177
uc7[["doublet"]]<-uc7$DF.classifications_0.25_0.29_289

uc <- merge(uc1, y = c(uc2,uc3,uc4,uc5,uc6,uc7), project = "UC_Colon_noninflammed")

uc <- seurat_qc2(uc)
png("ump_doublet.png")
show(DimPlot(uc, reduction = "umap",group.by = "doublet"))
dev.off()
write.csv(table(uc$doublet),"doublet.txt")
uc<- subset(uc,subset = doublet =="Singlet")
remove(uc1,uc2,uc4,uc3,uc5,uc6,uc7)
uc$pANN_0.25_0.29_56<- NULL
uc$pANN_0.25_0.29_61 <- NULL
uc$pANN_0.25_0.29_46<- NULL
uc$pANN_0.25_0.29_538<- NULL
uc$pANN_0.25_0.29_542<- NULL
uc$pANN_0.25_0.29_177<- NULL
uc$pANN_0.25_0.29_289<- NULL
uc$DF.classifications_0.25_0.29_56<- NULL
uc$DF.classifications_0.25_0.29_61<- NULL
uc$DF.classifications_0.25_0.29_538<- NULL
uc$DF.classifications_0.25_0.29_542<- NULL
uc$DF.classifications_0.25_0.29_177<- NULL
uc$DF.classifications_0.25_0.29_289<- NULL
uc$DF.classifications_0.25_0.29_46<- NULL

png("umap_celltype.png",width = 780,height = 480)
show(DimPlot(uc, reduction = "umap",group.by = "Cluster",label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
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

saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/uc_noninfla_pre_integration.rds") )
saveRDS(uc_infla_pre_integration,paste0("~/ds_group/Charlie Msc project data/uc/uc_infla_pre_integration.rds") )
