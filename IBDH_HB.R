library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
hb[["Sample"]]<-hb$orig.ident
hb[["Organ"]]<-"colon"
hb <- seurat_qc1(hb)
hb <- seurat_qc2(hb)
hb1 <- subset(hb, subset = orig.ident=="1")
hb2 <- subset(hb, subset = orig.ident=="2")
hb3 <- subset(hb, subset = orig.ident=="3")
hb4 <- subset(hb, subset = orig.ident=="4")
hb5 <- subset(hb, subset = orig.ident=="5")


nExp_poi <- round(0.02*length(hb1$orig.ident))
hb1 <- doubletFinder_v3(hb1, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.016*length(hb2$orig.ident))
hb2 <- doubletFinder_v3(hb2, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(hb3$orig.ident))
hb3 <- doubletFinder_v3(hb3, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(hb4$orig.ident))
hb4 <- doubletFinder_v3(hb4, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.004*length(hb5$orig.ident))
hb5 <- doubletFinder_v3(hb5, PCs = 1:17, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

hb1[["doublet"]]<-hb1$DF.classifications_0.25_0.01_82
hb2[["doublet"]]<-hb2$DF.classifications_0.25_0.01_52
hb3[["doublet"]]<-hb3$DF.classifications_0.25_0.01_23
hb4[["doublet"]]<-hb4$DF.classifications_0.25_0.01_30
hb5[["doublet"]]<-hb5$DF.classifications_0.25_0.01_1

hb <- merge(hb1, y = c(hb2,hb3,hb4,hb5), project = "IBD_Huang_HB")

hb <- seurat_qc2(hb)
png("ump_doublet.png")
show(DimPlot(hb, reduction = "umap",group.by = "doublet"))
dev.off()
write.csv(table(hb$doublet),"doublet.txt")
hb<- subset(hb,subset = doublet =="Singlet")
remove(hb1,hb2,hb4,hb3,hb5)
hb$pANN_0.25_0.01_82<- NULL
hb$pANN_0.25_0.01_52 <- NULL
hb$pANN_0.25_0.01_23<- NULL
hb$pANN_0.25_0.01_30<- NULL
hb$pANN_0.25_0.01_1<- NULL
hb$DF.classifications_0.25_0.01_82<- NULL
hb$DF.classifications_0.25_0.01_52<- NULL
hb$DF.classifications_0.25_0.01_23<- NULL
hb$DF.classifications_0.25_0.01_30<- NULL
hb$DF.classifications_0.25_0.01_1<- NULL

png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(hb, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(hb, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(hb, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(hb,paste0("~/ds_group/Charlie Msc project data/IBD_H/ibdh_hb_pre_integration.rds") )
