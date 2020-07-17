library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
s[["Sample"]]<-s$orig.ident
s[["Organ"]]<-"colon"
s <- seurat_qc1(s)
s <- seurat_qc2(s)
s1 <- subset(s, subset = orig.ident=="1")
s2 <- subset(s, subset = orig.ident=="2")
table(s$Sample)


nExp_poi <- round(0.016*length(s1$orig.ident))
s1 <- doubletFinder_v3(s1, PCs = 1:27, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.02*length(s2$orig.ident))
s2 <- doubletFinder_v3(s2, PCs = 1:27, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


s1[["doublet"]]<-s1$DF.classifications_0.25_0.01_54
s2[["doublet"]]<-s2$DF.classifications_0.25_0.01_91


s <- merge(s1, y = s2, project = "IBD_Huang_s")

s <- seurat_qc2(s)
png("ump_doublet.png")
show(DimPlot(s, reduction = "umap",group.by = "doublet"))
dev.off()
write.csv(table(s$doublet),"doublet.txt")
s<- subset(s,subset = doublet =="Singlet")

remove(s1,s2)
s$pANN_0.25_0.01_54<- NULL
s$pANN_0.25_0.01_91 <- NULL
s$DF.classifications_0.25_0.01_54<- NULL
s$DF.classifications_0.25_0.01_91<- NULL


png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(s, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker2.png",width = 900,height = 900)
show(FeaturePlot(s, features = c("CD3D", "CD3E","CD8A", "GZMA", "GZMB", "KLRF1", "TPSAB1","TPSB2","CDC20")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(s, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(s,paste0("~/ds_group/Charlie Msc project data/IBD_H/ibdh_s_pre_integration.rds") )
