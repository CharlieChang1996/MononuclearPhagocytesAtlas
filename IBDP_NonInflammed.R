library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
source("~/Msc project/data/seurat_qc.R")
ibdp[["Organ"]]<-"colon"

ibdp <- subset(ibdp, idents = c("A2","B2","C2"))

ibdp <- seurat_qc1(ibdp)
ibdp <- seurat_qc2(ibdp)

Idents(ibdp) <- "orig.ident"
ibdp1 <- subset(ibdp, idents = "A2")
ibdp2 <- subset(ibdp, idents = "B2")
ibdp3 <- subset(ibdp, idents = "C2")


nExp_poi <- round(0.002*length(ibdp1$orig.ident))
ibdp1 <- doubletFinder_v3(ibdp1, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*length(ibdp2$orig.ident))
ibdp2 <- doubletFinder_v3(ibdp2, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.007*length(ibdp3$orig.ident))
ibdp3 <- doubletFinder_v3(ibdp3, PCs = 1:20, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

ibdp1[["doublet"]]<-ibdp1$DF.classifications_0.25_0.06_1
ibdp2[["doublet"]]<-ibdp2$DF.classifications_0.25_0.06_16
ibdp3[["doublet"]]<-ibdp3$DF.classifications_0.25_0.06_12


ibdp <- merge(ibdp1, y = c(ibdp2,ibdp3), project = "ibdp_Colon_noninflammed")

ibdp <- seurat_qc2(ibdp)

png("ump_doublet.png")
show(DimPlot(ibdp, reduction = "umap",group.by = "doublet"))
dev.off()

write.csv(table(ibdp$doublet),"doublet.txt")

ibdp<- subset(ibdp,subset = doublet =="Singlet")

remove(ibdp1,ibdp2,ibdp3)
ibdp$pANN_0.25_0.06_1<- NULL
ibdp$pANN_0.25_0.06_16 <- NULL
ibdp$pANN_0.25_0.06_12<- NULL
ibdp$DF.classifications_0.25_0.06_1<- NULL
ibdp$DF.classifications_0.25_0.06_16<- NULL
ibdp$DF.classifications_0.25_0.06_12<- NULL


png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(ibdp, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(ibdp, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(ibdp, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(ibdp,paste0("~/ds_group/Charlie Msc project data/IBDP/ibdp_noninfla_pre_integration.rds") )
