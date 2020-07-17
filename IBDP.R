ibdp[["Sample"]]<-ibdp$orig.ident
ibdp[["Organ"]]<-"colon"
ibdp <- seurat_qc1(ibdp)
ibdp1 <- subset(ibdp, subset = orig.ident=="A1")
ibdp2 <- subset(ibdp, subset = orig.ident=="A2")
ibdp3 <- subset(ibdp, subset = orig.ident=="B1")
ibdp4 <- subset(ibdp, subset = orig.ident=="B2")
ibdp5 <- subset(ibdp, subset = orig.ident=="C1")
ibdp6 <- subset(ibdp, subset = orig.ident=="C2")

nExp_poi <- round(0.004*length(ibdp1$orig.ident))
ibdp1 <- doubletFinder_v3(ibdp1, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.004*length(ibdp2$orig.ident))
ibdp2 <- doubletFinder_v3(ibdp2, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*leibdp2(ibdp3$orig.ident))
ibdp3 <- doubletFinder_v3(ibdp3, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.01*length(ibdp4$orig.ident))
ibdp4 <- doubletFinder_v3(ibdp4, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*length(ibdp5$orig.ident))
ibdp5 <- doubletFinder_v3(ibdp5, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
nExp_poi <- round(0.008*length(ibdp6$orig.ident))
ibdp6 <- doubletFinder_v3(ibdp6, PCs = 1:22, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
ibdp1[["doublet"]]<-ibdp1$DF.classifications_0.25_0.06_2
ibdp2[["doublet"]]<-ibdp2$DF.classifications_0.25_0.06_1
ibdp3[["doublet"]]<-ibdp3$DF.classifications_0.25_0.06_1
ibdp4[["doublet"]]<-ibdp4$DF.classifications_0.25_0.06_20
ibdp5[["doublet"]]<-ibdp5$DF.classifications_0.25_0.06_17
ibdp6[["doublet"]]<-ibdp6$DF.classifications_0.25_0.06_14
ibdp <- merge(ibdp1, y = c(ibdp2,ibdp3,ibdp4,ibdp5,ibdp6), project = "IBD_Parikh")
ibdp <- seurat_qc2(ibdp)
png("ump_sample.png")
show(DimPlot(ibdp, reduction = "umap",group.by = "orig.ident"))
dev.off()
png("ump_doublet.png")
show(DimPlot(ibdp, reduction = "umap",group.by = "doublet"))
dev.off()
ibdp<- subset(ibdp,subset = doublet =="Singlet")
remove(ibdp1,ibdp2,ibdp4,ibdp3,ibdp5,ibdp6)
ibdp$pANN_0.25_0.06_1 <- NULL
ibdp$pANN_0.25_0.06_14<- NULL
ibdp$pANN_0.25_0.06_17<- NULL
ibdp$pANN_0.25_0.06_2<- NULL
ibdp$pANN_0.25_0.06_20<- NULL
ibdp$DF.classifications_0.25_0.06_1<- NULL
ibdp$DF.classifications_0.25_0.06_14<- NULL
ibdp$DF.classifications_0.25_0.06_17<- NULL
ibdp$DF.classifications_0.25_0.06_2<- NULL
ibdp$DF.classifications_0.25_0.06_20<- NULL
png("feature_marker.png",width = 900,height = 900)
show(FeaturePlot(ibdp, features = c("CD79A", "CD79B","CD19", "CD27", "IGHG1", "LILRA4", "IRF7","LYZ","CD68")))
dev.off()
saveRDS(ibdp,paste0("~/ds_group/Charlie Msc project data/IBDP/ibdp_pre_integration.rds") )
