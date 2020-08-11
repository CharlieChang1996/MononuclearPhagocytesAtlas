source("/home/s1937334/Msc project/Script/seurat_qc.R")
library(RColorBrewer)
library(ggplot2)

gut.integrated <- readRDS("~/Msc project/data/downstream/colon_HvD_7_18.rds")

write.csv(table(gut.integrated$annotation1, gut.integrated$Health),"propotion_cluster.csv")
write.csv(table(gut.integrated$annotation_major, gut.integrated$Health),"propotion_major.csv")
write.csv(table(gut.integrated$annotation_major, gut.integrated$Health,gut.integrated$Sample),"propotion_major_sample.csv")
write.csv(table(gut.integrated$annotation1, gut.integrated$Health,gut.integrated$Sample),"propotion_cluster_sample.csv")
write.csv(table(gut.integrated$annotation1, gut.integrated$dataname,gut.integrated$Sample),"propotion_cluster_dataname.csv")
#Mann-Whitney-Wilcoxon Test
#wilcox.test(Freq ~ Health, data=stat_t)

color1 <- brewer.pal(3,"BrBG")

## visualisation of proportion of cells
png("proportion_celltype1.png",width = 980)
barplot(cbind(pro_h[1:10],pro_inf[1:10],pro_non[1:10]) ~ Celltype[1:10], data = propotion_cluster,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflamed","Non-inflamed"),
       fill = color1)
dev.off()
png("proportion_celltype2.png",width = 1080)
barplot(cbind(pro_h[11:20],pro_inf[11:20],pro_non[11:20]) ~ Celltype[11:20], data = propotion_cluster,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflamed","Non-inflamed"),
       fill = color1)
dev.off()
png("proportion_celltype3.png",width = 980)
barplot(cbind(pro_h[21:30],pro_inf[21:30],pro_non[21:30]) ~ Celltype[21:30], data = propotion_cluster,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflamed","Non-inflamed"),
       fill = color1)
dev.off()
png("proportion_celltype4.png",width = 780)
barplot(cbind(pro_h[31:35],pro_inf[31:35],pro_non[31:35]) ~ Celltype[31:35], data = propotion_cluster,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflamed","Non-inflamed"),
       fill = color1)
dev.off()

#log-ratio plots of major cell types
inf.major.boxplot <- ggplot(propotion_major, aes(x = log_i, y = Celltypes, color = log_i)) +
        geom_point(size = 4) +
        geom_vline(xintercept = 0.4522507, linetype="dashed") +
        geom_segment(aes(x=0, y=7.5, xend=-3, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.91, y=7.5, xend=3.91, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-4, 4)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n inflamed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_major_inf_boxplot.png", width=500,height=800,units="px")
print(inf.major.boxplot)
dev.off()

non.major.boxplot <- ggplot(propotion_major, aes(x = log_n, y = Celltypes, color = log_n)) +
        geom_point(size = 4) +
        geom_vline(xintercept = -0.19171347, linetype="dashed") +
        geom_segment(aes(x=-0.7, y=7.5, xend=-3.7, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.3, y=7.5, xend=3.3, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-4, 4)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n Non-inflamed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_major_non_boxplot.png", width=500,height=800,units="px")
print(non.major.boxplot)
dev.off()

inon.major.boxplot <- ggplot(propotion_major, aes(x = log_in, y = Celltypes, color = log_in)) +
        geom_point(size = 4) +
        geom_vline(xintercept = 0.6439642, linetype="dashed") +
        geom_segment(aes(x=0.3, y=7.5, xend=-2, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=1.28, y=7.5, xend=3.28, yend=7.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-4, 4)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n inflamed/non-inflamed)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_major_innon_boxplot.png", width=500,height=800,units="px")
print(inon.major.boxplot)
dev.off()

#log-ratio plots of subclusters
inf.celltype.boxplot <- ggplot(propotion_cluster, aes(x = log_i, y = Celltype,color = log_i)) +
        geom_point(size = 4) +
        geom_vline(xintercept = 0.4522507, linetype="dashed") +
        geom_segment(aes(x=0, y=12.5, xend=-5, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.91, y=12.5, xend=5.91, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-13, 13)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n inflamed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_celltype_inf_logplot.png", width=500,height=1000,units="px")
print(inf.celltype.boxplot)
dev.off()

non.celltype.boxplot <- ggplot(propotion_cluster, aes(x = log_n, y = Celltype, color = log_n)) +
        geom_point(size = 4) +
        geom_vline(xintercept = -0.19171347, linetype="dashed") +
        geom_segment(aes(x=-0.7, y=12.5, xend=-5.7, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.3, y=12.5, xend=5.3, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-13, 13)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n Non-inflamed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_celltype_non_logplot.png", width=500,height=1000,units="px")
print(non.celltype.boxplot)
dev.off()

inon.celltype.boxplot <- ggplot(propotion_cluster, aes(x = log_in, y = Celltype, color = log_in)) +
        geom_point(size = 4) +
        geom_vline(xintercept = 0.6439642, linetype="dashed") +
        geom_segment(aes(x=0, y=12.5, xend=-5, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=1.28, y=12.5, xend=6.28, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-10, 10)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n inflamed/non-inflamed)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15), 
              axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_celltype_innon_logplot.png", width=500,height=1000,units="px")
print(inon.celltype.boxplot)
dev.off()

#Find DEGs 
DefaultAssay(gut.integrated) <- "RNA"
Idents(gut.integrated) <- gut.integrated$annotation1

monocytes <- subset(gut.integrated,idents = "Infla monocytes")
Idents(monocytes) <- monocytes$Health
c16i_n.markers <- FindMarkers(monocytes, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD", only.pos = TRUE)
c16.markers <- FindMarkers(gut.integrated, ident.1 = "Infla monocytes", ident.2 = NULL, only.pos = TRUE)
write.csv(c16.markers,"c16 markers.csv")
macro.markers <- FindMarkers(gut.integrated, ident.1 = "Cycling monocytes", ident.2 = "Macrophages")
write.csv(macro.markers,"cm vs macro markers.csv")

treg<- subset(gut.integrated,idents = "Treg1")
Idents(treg) <- treg$Health
treg_i_h.markers <- FindMarkers(treg, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
treg_n_h.markers <- FindMarkers(treg, ident.1 = "Non-Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(treg_n_h.markers,"Treg1_n_h markers.csv")
write.csv(treg_i_h.markers,"Treg1_i_h markers.csv")

B.cell1<- subset(gut.integrated,idents = "B Cells1")
Idents(B.cell1) <- B.cell1$Health
B.cell1_i_h.markers <- FindMarkers(B.cell1, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(B.cell1_i_h.markers,"B cell1_i_h markers.csv")

B.cell2<- subset(gut.integrated,idents = "B Cells2")
Idents(B.cell2) <- B.cell2$Health
B.cell2_i_h.markers <- FindMarkers(B.cell2, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(B.cell2_i_h.markers,"B cell2_i_h markers.csv")

ilc<- subset(gut.integrated,idents = "ILCs")
Idents(ilc) <- ilc$Health
ilc_i_h.markers <- FindMarkers(ilc, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(ilc_i_h.markers,"ILC_i_h markers.csv")

mast<- subset(gut.integrated,idents = "CD69-Mast")
Idents(mast) <- mast$Health
mast_i_h.markers <- FindMarkers(mast, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(mast_i_h.markers,"CD69-Mast_i_h markers.csv")

goblet2<- subset(gut.integrated,idents = "Goblet2")
Idents(goblet2) <- goblet2$Health
goblet2_i_h.markers <- FindMarkers(goblet2, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(goblet2_i_h.markers,"Goblet2_i_h markers.csv")

cyclingt3<- subset(gut.integrated,idents = "Cycling T3")
Idents(cyclingt3) <- cyclingt3$Health
cyclingt3_i_h.markers <- FindMarkers(cyclingt3, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
cyclingt3_i_n.markers <- FindMarkers(cyclingt3, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD", only.pos = TRUE)
write.csv(cyclingt3_i_n.markers,"Cycling T3_i_n markers.csv")
write.csv(cyclingt3_i_h.markers,"Cycling T3_i_h markers.csv")

cyclingb<- subset(gut.integrated,idents = "Cycling B1")
Idents(cyclingb) <- cyclingb$Health
cyclingb_i_n.markers <- FindMarkers(cyclingb, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD", only.pos = TRUE)
write.csv(cyclingb_i_n.markers,"Cycling B_i_n markers.csv")

plasma<- subset(gut.integrated,idents = "Plasma")
Idents(plasma) <- plasma$Health
plasma_i_h.markers <- FindMarkers(plasma, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(plasma_i_h.markers,"Plasma_i_h markers.csv")

memory_t<- subset(gut.integrated,idents = "Memory T")
Idents(memory_t) <- memory_t$Health
memory_t_i_h.markers <- FindMarkers(memory_t, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
memory_t_i_n.markers <- FindMarkers(memory_t, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD", only.pos = TRUE)
write.csv(memory_t_i_h.markers,"Memory_T_i_h markers.csv")
write.csv(memory_t_i_n.markers,"Memory_T_i_n markers.csv")

fos_lo_t<- subset(gut.integrated,idents = "Activated Fos-lo T")
Idents(fos_lo_t) <- fos_lo_t$Health
fos_lo_t_i_h.markers <- FindMarkers(fos_lo_t, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(fos_lo_t_i_h.markers,"Active_Fos_loT_i_h markers.csv")

c.monocytes<- subset(gut.integrated,idents = "Cycling monocytes")
Idents(c.monocytes) <- c.monocytes$Health
c.monocytes_i_h.markers <- FindMarkers(c.monocytes, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
c.monocytes_n_h.markers <- FindMarkers(c.monocytes, ident.1 = "Non-Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(c.monocytes_i_h.markers,"Cycling monocytes_i_h markers.csv")
write.csv(c.monocytes_n_h.markers,"Cycling monocytes_n_h markers.csv")

dc1<- subset(gut.integrated,idents = "DC1")
Idents(dc1) <- dc1$Health
dc1_i_h.markers <- FindMarkers(dc1, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(dc1_i_h.markers,"DC1_i_h markers.csv")

dc3<- subset(gut.integrated,idents = "DC3")
Idents(dc3) <- dc3$Health
dc3_i_h.markers <- FindMarkers(dc3, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(dc3_i_h.markers,"DC3_i_h markers.csv")

mast69<- subset(gut.integrated,idents = "CD69+Mast")
Idents(mast69) <- mast69$Health
mast69_i_h.markers <- FindMarkers(mast69, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(mast69_i_h.markers,"CD69+Mast_i_h markers.csv")

gcb<- subset(gut.integrated,idents = "GC B")
Idents(gcb) <- gcb$Health
gcb_i_h.markers <- FindMarkers(gcb, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(gcb_i_h.markers,"GC_B_i_h markers.csv")

lp2<- subset(gut.integrated,idents = "CD8+LP2 T")
Idents(lp2) <- lp2$Health
lp2_i_h.markers <- FindMarkers(lp2, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(lp2_i_h.markers,"CD8+LP2 T_i_h markers.csv")

##Plot plots of marker genes across clusters
t.markergenes <- c("RPL39","RPS29","LTB","GPR183","CREM","YPEL5","FTH1","TNFRSF4",
                  "CTLA4","MAF","ICOS","STMN1","HMGB2","IL32","CD3D","GZMK","CD8A","CD8B","GNLY","NKG7","CCL4",
                  "CST7")
b.markergenes  <-  c("HLA-DRB1","HLA-DRA","CD83","HLA-DPA1","MS4A1","TCL1A","RGS13","HMGB2",
                     "STMN1","KIAA0101","UBE2C","IGLL5","IGLC2","JUN","TPSAB1","CD69","CD19")
m.markergenes <- c("SLC40A1","SEPP1","HSPA1B","S100A8","S100A9","FCN1","LAMP3","HLA-DPB1","IDO1",
                   "TFF1","SPINK4","CADM1","MT-ND3","NEAT1","C1QC","MMP12","LYZ","TNF")
uc.genes <- c("TFF3","HMGB1","ANXA1","MMP9","ADCY7","TNF","TNFRSF4","TNFSF15")
#m1_violin <- VlnPlot(gut.integrated, features = markergenes1, slot = "counts", log = TRUE)



#Heatmaps
DefaultAssay(gut.integrated) <- "integrated"
Idents(gut.integrated) <- gut.integrated$annotation1
tcells <- c("Memory T","Activated Fos-lo T","CD4+PD1+T","Activated Fos-hi T","Treg1",
            "CD8+LP1 T","CD8+LP2 T","CD8+LP3 T","CD4 T","Cycling T1","NKs 1","Treg2",
            "NKs2","CD8+IL17+ T","ILCs","Cycling T2","Cycling T3")
mcells <-c("DC1","DC2","DC3","Infla monocytes","Macrophages","Cycling monocytes")
bcells <- c("Follicular B","B Cells2","GC B","Cycling B1","Plasma","CD69+Mast","CD69-Mast")

t_heatmap <- DoHeatmap(subset(gut.integrated, downsample = 100,idents = tcells),features = t.markergenes, size = 4)
png( "t.markergene_heatmap.png",width=1050,height=650,units="px")
print(t_heatmap)
dev.off()

m_heatmap <- DoHeatmap(subset(gut.integrated, downsample = 100,idents = mcells),features = m.markergenes, size = 5)
png( "m.markergene_heatmap.png",width=1050,height=650,units="px")
print(m_heatmap)
dev.off()

b_heatmap <- DoHeatmap(subset(gut.integrated, downsample = 100,idents = bcells),features = b.markergenes, size = 5)
png( "b.markergene_heatmap.png",width=1050,height=650,units="px")
print(b_heatmap)
dev.off()


#Dotplots for DEGs across clusters
DefaultAssay(gut.integrated) <- "RNA"
Idents(gut.integrated) <- gut.integrated$annotation_major
t_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "T Cells"), features = t.markergenes,group.by = "annotation1") + RotatedAxis()
png( "t.markergene_dotplot.png",width=1050,height=650,units="px")
print(t_dotplot)
dev.off()

b_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = c("B Cells","Mast","Plasma")), features = b.markergenes,group.by = "annotation1") + RotatedAxis()
png( "b.markergene_dotplot.png",width=1050,height=650,units="px")
print(b_dotplot)
dev.off()

m_dotplot <- DotPlot(subset(gut.integrated,idents = "Myeloid Cells"), features = m.markergenes,group.by = "annotation1") + RotatedAxis()
png( "m.markergene_dotplot.png",width=1050,height=650,units="px")
print(m_dotplot)
dev.off()

uc_dotplot <- DotPlot(subset(gut.integrated), features = uc.genes ,group.by = "annotation1") + RotatedAxis()
png( "uc_related_gene_dotplot.png",width=1050,height=650,units="px")
print(uc_dotplot)
dev.off()

Idents(gut.integrated) <- gut.integrated$annotation1
uc2_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Cycling monocytes"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_cyclingmonocytes_dotplot.png",width=1050,height=650,units="px")
print(uc2_dotplot)
dev.off()

uc3_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Cycling B1"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_cyclingB_dotplot.png",width=1050,height=650,units="px")
print(uc3_dotplot)
dev.off()

uc4_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "CD69+Mast"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_cd69Mast_dotplot.png",width=1050,height=650,units="px")
print(uc4_dotplot)
dev.off()

uc5_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Plasma"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_Plasma_dotplot.png",width=1050,height=650,units="px")
print(uc5_dotplot)
dev.off()

uc6_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Cycling T3"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_CyclingT3_dotplot.png",width=1050,height=650,units="px")
print(uc6_dotplot)
dev.off()

uc7_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Treg1"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_Treg1_dotplot.png",width=1050,height=650,units="px")
print(uc7_dotplot)
dev.off()

uc8_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Infla monocytes"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_inflamonocytes_dotplot.png",width=1050,height=650,units="px")
print(uc8_dotplot)
dev.off()

uc9_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "DC1"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_DC1_dotplot.png",width=1050,height=650,units="px")
print(uc9_dotplot)
dev.off()

uc10_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "B Cells2"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_Bcells_dotplot.png",width=1050,height=650,units="px")
print(uc10_dotplot)
dev.off()

uc11_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Memory T"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_memoryT_dotplot.png",width=1050,height=650,units="px")
print(uc11_dotplot)
dev.off()

uc12_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "CD8+LP2 T"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_LP2 T_dotplot.png",width=1050,height=650,units="px")
print(uc12_dotplot)
dev.off()

uc13_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Activated Fos-lo T"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_foslo T_dotplot.png",width=1050,height=650,units="px")
print(uc13_dotplot)
dev.off()

uc0_dotplot <- DotPlot(gut.integrated, features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_health_dotplot.png",width=1050,height=650,units="px")
print(uc0_dotplot)
dev.off()