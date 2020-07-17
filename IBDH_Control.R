library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
library(DoubletFinder)
library(SCINA)
library(harmony)
col[["Sample"]]<-col$orig.ident
col[["Organ"]]<-"colon"
col <- seurat_qc1(col)
col1 <- subset(col, subset = orig.ident=="1")
col2 <- subset(col, subset = orig.ident=="2")
col3 <- subset(col, subset = orig.ident=="3")
col4 <- subset(col, subset = orig.ident=="4")
col5 <- subset(col, subset = orig.ident=="5")
col6 <- subset(col, subset = orig.ident=="6")

nExp_poi <- round(0.004*length(col1$orig.ident))
col1 <- doubletFinder_v3(col1, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.004*length(col2$orig.ident))
col2 <- doubletFinder_v3(col2, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*lecol2(col3$orig.ident))
col3 <- doubletFinder_v3(col3, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(col4$orig.ident))
col4 <- doubletFinder_v3(col4, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*length(col5$orig.ident))
col5 <- doubletFinder_v3(col5, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*length(col6$orig.ident))
col6 <- doubletFinder_v3(col6, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
col1[["doublet"]]<-col1$DF.classifications_0.25_0.06_2
col2[["doublet"]]<-col2$DF.classifications_0.25_0.06_1
col3[["doublet"]]<-col3$DF.classifications_0.25_0.06_1
col4[["doublet"]]<-col4$DF.classifications_0.25_0.06_20
col5[["doublet"]]<-col5$DF.classifications_0.25_0.06_17
col6[["doublet"]]<-col6$DF.classifications_0.25_0.06_14
col <- merge(col1, y = c(col2,col3,col4,col5,col6), project = "IBD_Parikh")
col <- seurat_qc2(col)
png("ump_sample.png")
show(DimPlot(col, reduction = "umap",group.by = "orig.ident"))
dev.off()
png("ump_doublet.png")
show(DimPlot(col, reduction = "umap",group.by = "doublet"))
dev.off()
col<- subset(col,subset = doublet =="Singlet")
remove(col1,col2,col4,col3,col5,col6)
col$pANN_0.25_0.01_240<- NULL
col$pANN_0.25_0.01_428 <- NULL
col$pANN_0.25_0.01_461<- NULL
col$pANN_0.25_0.01_8<- NULL
col$pANN_0.25_0.01_9<- NULL
col$DF.classifications_0.25_0.01_240<- NULL
col$DF.classifications_0.25_0.01_428<- NULL
col$DF.classifications_0.25_0.01_461<- NULL
col$DF.classifications_0.25_0.01_8<- NULL
col$DF.classifications_0.25_0.01_9<- NULL

png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(col, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
png("feature_marker_epi.png",width = 780,height = 480)
show(FeaturePlot(col, features = c("EPCAM", "KRT19")))
dev.off()

saveRDS(col,paste0("~/ds_group/Charlie Msc project data/IBD_H/ibdh_col_pre_integration.rds") )
