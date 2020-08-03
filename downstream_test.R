source("/home/s1937334/Msc project/Script/seurat_qc.R")

colon <- readRDS("~/ds_group/Charlie Msc project data/Colon/colon_anchor.rds")
uc_n <- readRDS("~/ds_group/Charlie Msc project data/uc/uc_noninfla_pre_integration.rds")
uc_i <- readRDS("~/ds_group/Charlie Msc project data/uc/uc_infla_pre_integration.rds")
ibdhs <- readRDS("~/ds_group/Charlie Msc project data/IBD_H/ibdh_s_disease_preintegration.rds")
ibdhb <- readRDS("~/ds_group/Charlie Msc project data/IBD_H/ibdh_hb_pre_integration.rds")
ibdpn <- readRDS("~/ds_group/Charlie Msc project data/IBDP/ibdp_noninfla_pre_integration.rds")
ibdpi <- readRDS("~/ds_group/Charlie Msc project data/IBDP/ibdp_infla_pre_integration.rds")



colon[["Health"]] <- "Healthy"
ibdhs[["Health"]] <- "Inflammed IBD"
ibdpi[["Health"]] <- "Inflammed IBD"
ibdpn[["Health"]] <- "Non-Inflammed IBD"
ibdhb[["Health"]] <- "Inflammed IBD"
uc_i[["Health"]] <- "Inflammed IBD"
uc_n[["Health"]] <- "Non-Inflammed IBD"

ibdhs[["dataname"]] <- "IBD_Huang_s_Inflammed"
ibdpi[["dataname"]] <- "IBD_Parikh_Inflammed"
ibdpn[["dataname"]] <- "IBD_Parikh_Non-Inflammed"
ibdhb[["dataname"]] <- "IBD_Huang_hb_Inflammed"
uc_i[["dataname"]] <- "UC_Inflammed"
uc_n[["dataname"]] <- "UC_Non-Inflammed"

gut.list <- list(colon,ibdhs,ibdhb,ibdpi,ibdpn,uc_i,uc_n)
#gut,integrated <- seurat_inte(gut.list)
gut.anchors <- FindIntegrationAnchors(object.list = gut.list, dims = 1:20)
gut.integrated <- IntegrateData(anchorset = gut.anchors, dims = 1:20)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(gut.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
gut.integrated <- ScaleData(gut.integrated, verbose = FALSE)
gut.integrated <- RunPCA(gut.integrated, npcs = 20, verbose = FALSE)
gut.integrated <- RunUMAP(gut.integrated, reduction = "pca", dims = 1:20)

write.csv(table(gut.integrated$dataname),"study.txt")
write.csv(table(gut.integrated$Health),"health.txt")

#visualisation
png("umap_anchor.png",width = 780)
show(DimPlot(gut.integrated, reduction = "umap", group.by = "Health"))
dev.off()
png("umap_anchor_study.png",width = 780)
show(DimPlot(gut.integrated, reduction = "umap", group.by = "dataname"))
dev.off()

png("anchor_ibdh.png")
DimPlot(object = subset(gut.integrated, subset = dataname ==c("IBD_Huang_s_Inflammed","IBD_Huang_hb_Inflammed")), reduction = "umap",group.by = "dataname")
dev.off()
png("anchor_ibdp.png")
DimPlot(object = subset(gut.integrated, subset = dataname =="IBD_Parikh_Inflammed"), reduction = "umap",group.by = "dataname")
dev.off()
png("anchor_uc.png",width = 780)
DimPlot(object = subset(gut.integrated, subset = dataname =="UC_Inflammed"), reduction = "umap",group.by = "Health")
dev.off()
png("anchor_healthy.png",width = 780)
DimPlot(object = subset(gut.integrated, subset = Health =="Healthy"), reduction = "umap",group.by = "dataname")
dev.off()
png("anchor_inflammed.png",width = 780)
DimPlot(object = subset(gut.integrated, subset = Health =="Inflammed IBD"), reduction = "umap",group.by = "dataname")
dev.off()
png("anchor_Non-Inflammed.png",width = 780)
DimPlot(object = subset(gut.integrated, subset = Health =="Non-Inflammed IBD"), reduction = "umap",group.by = "dataname")
dev.off()
png("anchor_health.png",width = 980)
show(DimPlot(gut.integrated, reduction = "umap", split.by = "Health"))
dev.off()

##Clustering & Annotation
gut.integrated <- seurat_clus(gut.integrated)# 20dims res=1.8
write.csv(table(gut.integrated$seurat_clusters),"cluster.txt")


#find markers
monocytes.markers <-c(subset(gut.integrated,idents = c(16,21,23,30,43,46,49)), 
                                    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tcell.markers <-FindAllMarkers(subset(gut.integrated,idents = c(3,8,9,10,12,13,14,26,28,29,36,37,42)), 
                               only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rest.markers <-FindAllMarkers(subset(gut.integrated,idents = c(16,21,23,30,43,46,49,3,8,9,10,12,13,14,26,28,29,36,37,42),invert =TRUE), 
                              only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
c26.markers <- FindMarkers(gut.integrated,ident.1 = 26,min.pct = 0.25)
c46.markers <- FindMarkers(gut.integrated,ident.1 = 46,min.pct = 0.25)
c5.markers <- FindMarkers(gut.integrated,ident.1 = 5,min.pct = 0.25)

write.csv(monocytes.markers,"monocytes_markers.csv")
write.csv(tcell.markers,"T_markers.csv")
write.csv(rest.markers,"rest_markers.csv")
write.csv(c26.markers,"c26_markers.csv")
write.csv(c46.markers,"c46_markers.csv")
write.csv(c5.markers,"c5_markers.csv")
#Assign celltypes manually
Idents(gut.integrated) <- gut.integrated$seurat_clusters
new.cluster.ids <- c("Follicular B","Goblet","Plasma","Memory T","Plasma","Follicular B",
                     "Plasma","Plasma","Activated Fos-lo T","CD8+IL17+T","Activated Fos-hi T","Plasma",
                     "Treg1","CD8+LP1 T","CD8+LP2 T","Follicular B","Infla monocytes","Plasma",
                     "CD69+Mast","CD8+LP3 T","GC B","DC2","Plasma","Cycling monocytes",
                     "Plasma","Plasma","CD4 T","Follicular B","Cycling T1","NKs 1",
                     "Macrophages","Treg2","Plasma","Cycling B1","Memory T","Epithelias2",
                     "NKs2","CD8+IEL T","Plasma","B Cells2","Epithelia","Plasma",
                     "ILCs","DC1","CD69-Mast","Cycling T2","DC3","Cycling T3",
                     "Cycling T2","Plasma","Doublet?","Plasma")
names(new.cluster.ids) <- levels(gut.integrated)
gut.integrated <- RenameIdents(gut.integrated,new.cluster.ids)
gut.integrated[["annotation1"]] <- Idents(gut.integrated)

new.cluster.ids2 <- c("B Cells","Epithelias","Plasma","T Cells",
                      "T Cells","T Cells","T Cells","T Cells",
                      "T Cells","T Cells","Myeloid Cells","Mast",
                      "T Cells","B Cells","Myeloid Cells","Myeloid Cells",
                      "T Cells","T Cells","NKs","Myeloid Cells",
                      "T Cells","B Cells","Epithelias","NKs",
                      "T Cells","B Cells","Epithelias","ILCs",
                      "Myeloid Cells","Mast","T Cells","Myeloid Cells",
                      "T Cells","Doublet?")
names(new.cluster.ids2) <- levels(gut.integrated)
gut.integrated <- RenameIdents(gut.integrated,new.cluster.ids2)
gut.integrated[["annotation_major"]] <- Idents(gut.integrated)

png("umap_annotation1.png",width = 1080)
show(DimPlot(gut.integrated, reduction = "umap",group.by = 'annotation1', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()
png("umap_annotation_major.png",width = 1080)
show(DimPlot(gut.integrated, reduction = "umap",group.by = 'annotation_major', label = TRUE,pt.size = 0.5,label.size = 4,repel = TRUE))
dev.off()
write.csv(table(gut.integrated$annotation1),"annotation.txt")

saveRDS(gut.integrated,paste0("~/Msc project/data/downstream/colon_HvD_7_18.rds") )

#SCINA
gut.integrated <- clus_anno(gut.integrated)

saveRDS(gut.integrated,paste0("~/ds_group/Charlie Msc project data/downstream/colon_ibdhs.rds") )

