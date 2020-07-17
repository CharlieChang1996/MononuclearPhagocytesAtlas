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
uc <- subset(uc, idents = c("Inflamed"))

uc <- seurat_qc1(uc)
uc <- seurat_qc2(uc)


uc1 <- subset(uc, subset = Sample=="N111.LPB1")
uc2 <- subset(uc, subset = Sample=="N661.LPB1")
uc3 <- subset(uc, subset = Sample=="N661.LPB2")
uc4 <- subset(uc, subset = Sample=="N58.LPB1")
Idents(uc)<- "Sample"
uc5 <- subset(uc, idents = c("N111.LPB1","N58.LPB1","N661.LPB1",
              "N661.LPB2"),invert = TRUE)

nExp_poi <- round(0.054*length(uc1$orig.ident))
uc1 <- doubletFinder_v3(uc1, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.04*length(uc2$orig.ident))
uc2 <- doubletFinder_v3(uc2, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.04*length(uc3$orig.ident))
uc3 <- doubletFinder_v3(uc3, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.062*length(uc4$orig.ident))
uc4 <- doubletFinder_v3(uc4, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.02*length(uc5$orig.ident))
uc5 <- doubletFinder_v3(uc5, PCs = 1:20, pN = 0.25, pK = 0.29, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

uc1[["doublet"]]<-uc1$DF.classifications_0.25_0.29_626
uc2[["doublet"]]<-uc2$DF.classifications_0.25_0.29_352
uc3[["doublet"]]<-uc3$DF.classifications_0.25_0.29_339
uc4[["doublet"]]<-uc4$DF.classifications_0.25_0.29_818
uc5[["doublet"]]<-uc5$DF.classifications_0.25_0.29_982

uc <- merge(uc1, y = c(uc2,uc3,uc4,uc5), project = "UC_Colon_inflammed")

uc <- seurat_qc2(uc)
png("ump_doublet.png")
show(DimPlot(uc, reduction = "umap",group.by = "doublet"))
dev.off()
write.csv(table(uc$doublet),"doublet.txt")
uc<- subset(uc,subset = doublet =="Singlet")
remove(uc1,uc2,uc4,uc3,uc5)
uc$pANN_0.25_0.29_626<- NULL
uc$pANN_0.25_0.29_352 <- NULL
uc$pANN_0.25_0.29_339<- NULL
uc$pANN_0.25_0.29_818<- NULL
uc$pANN_0.25_0.29_982<- NULL
uc$DF.classifications_0.25_0.29_626<- NULL
uc$DF.classifications_0.25_0.29_352<- NULL
uc$DF.classifications_0.25_0.29_339<- NULL
uc$DF.classifications_0.25_0.29_818<- NULL
uc$DF.classifications_0.25_0.29_982<- NULL

png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(uc, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(uc, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(uc,paste0("~/ds_group/Charlie Msc project data/uc/uc_infla_pre_integration.rds") )
