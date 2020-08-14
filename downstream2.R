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

##Find DEGs between conditions 

de_sub(gut.integrated,"Treg1")
Treg1_i_h markers.csv
Treg1_n_h markers.csv
Treg1_i_n markers.csv

de_sub(gut.integrated,"B Cells1")
B cell1_i_h markers.csv
B cell1_n_h markers.csv
B cell1_i_n markers.csv

de_sub(gut.integrated,"B Cells2")
B cell2_i_h markers.csv
B cell2_n_h markers.csv
B cell2_i_n markers.csv

de_sub(gut.integrated,"ILCs")
ilc_i_h markers.csv
ilc_n_h markers.csv
ilc_i_n markers.csv

de_sub(gut.integrated,"CD69-Mast")
CD69-Mast_i_h markers.csv
CD69-Mast_n_h markers.csv
CD69-Mast_i_n markers.csv

de_sub(gut.integrated,"Goblet")
Goblet_i_h markers.csv
Goblet_n_h markers.csv
Goblet_i_n markers.csv

de_sub(gut.integrated,"CD4+ T2")
Cd4+ T2_i_h markers.csv
Cd4+ T2_n_h markers.csv
Cd4+ T2_i_n markers.csv

de_sub(gut.integrated,"Plasma")
plasma_i_h markers.csv
plasma_n_h markers.csv
plasma_i_n markers.csv

de_sub(gut.integrated,"Memory T")
Memory_T_i_h markers.csv
Memory_T_n_h markers.csv
Memory_T_i_n markers.csv

de_sub(gut.integrated,"Activated Fos-lo T")
Activated Fos-lo T_i_h markers.csv
Activated Fos-lo T_n_h markers.csv
Activated Fos-lo T_i_n markers.csv

de_sub(gut.integrated,"Macrophages1")
Macrophages1_i_h markers.csv
Macrophages1_n_h markers.csv
Macrophages1_i_n markers.csv

de_sub(gut.integrated,"Macrophages2")
Macrophages2_i_h markers.csv
Macrophages2_n_h markers.csv
Macrophages2_i_n markers.csv

de_sub(gut.integrated,"DC1")
DC1_i_h markers.csv
DC1_n_h markers.csv
DC1_i_n markers.csv

de_sub(gut.integrated,"DC3")
DC3_i_h markers.csv
DC3_n_h markers.csv
DC3_i_n markers.csv

de_sub(gut.integrated,"CD69+Mast")
CD69+Mast_i_h markers.csv
CD69+Mast_n_h markers.csv
CD69+Mast_i_n markers.csv

de_sub(gut.integrated,"GC B")
gcb_i_h markers.csv
gcb_n_h markers.csv
gcb_i_n markers.csv

de_sub(gut.integrated,"CD8+LP2 T")
CD8+LP2 T_i_h markers.csv
CD8+LP2 T_n_h markers.csv
CD8+LP2 T_i_n markers.csv



##Plot plots of marker genes across clusters
t.markergenes <- c("RPL39","RPS29","LTB","LST1","GPR183","CREM","YPEL5","FTH1","TNFRSF4",
                  "CTLA4","MAF","ICOS","STMN1","HMGB2","IL17A","CD3D","GZMK","CD8A","CD8B","GNLY","NKG7","CCL4",
                  "CST7")
b.markergenes  <-  c("HLA-DRB1","HLA-DRA","CD83","HLA-DPA1","MS4A1","TCL1A","RGS13","HMGB2",
                     "STMN1","KIAA0101","UBE2C","IGLL5","IGHG1","JUN","MZB1","CD79A")
m.markergenes <- c("SLC40A1","SEPP1","HSPA1B","S100A8","S100A9","FCN1","LAMP3","HLA-DPB1","IDO1",
                   "TFF1","SPINK4","CADM1","MT-ND3","NEAT1","C1QC","MMP12","LYZ","TNF")
ma.markergenes <- c("TPSAB1","CTSG","CD69","KRT1","MAOB")
epi.markergenes <- c("EPCAM","RRM2","STMN1","BIRC5","NEAT1","EGR1","IFI27","AGR2","SPINK4","TFF3","PHGR1","FXYD3","SEPP1")
uc.genes <- c("TFF3","HMGB1","ANXA1","MMP9","ADCY7","TNF","TNFRSF4","TNFSF15")
#m1_violin <- VlnPlot(gut.integrated, features = markergenes1, slot = "counts", log = TRUE)



#Heatmaps
DefaultAssay(gut.integrated) <- "integrated"
Idents(gut.integrated) <- gut.integrated$annotation1
tcells <- c("Memory T","Activated Fos-lo T","Activated Fos-hi T","CD4+ T1","CD4+PD1+T","Treg1","Treg2",
            "CD8+LP1 T","CD8+LP2 T","CD8+LP3 T","CD8+IL17+ T","Cycling T","NKs 1",
            "NKs2","ILCs","CD4+ T2")
mcells <-c("DC1","DC2","Infla monocytes","Macrophages1","Macrophages2","Macrophages3")
bcells <- c("B Cells1","B Cells2","B Cells3","Follicular B","GC B","Plasma")

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

ma_heatmap <- DoHeatmap(subset(gut.integrated, downsample = 100,idents = c("CD69+Mast","CD69-Mast")),features = ma.markergenes, size = 5)
png( "mast.markergene_heatmap.png",width=1050,height=650,units="px")
print(ma_heatmap)
dev.off()

epi_heatmap <- DoHeatmap(subset(gut.integrated, downsample = 100,idents = c("Cycling TA","Goblet","Epithelia","Epithelias2")),features = epi.markergenes, size = 5)
png( "epi.markergene_heatmap.png",width=1050,height=650,units="px")
print(epi_heatmap)
dev.off()
