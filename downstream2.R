source("/home/s1937334/Msc project/Script/seurat_qc.R")
library(RColorBrewer)
library(ggplot2)
gut.integrated <- readRDS("~/Msc project/data/downstream/colon_HvD_7_18.rds")

write.csv(table(gut.integrated$annotation1, gut.integrated$Health),"propotion_cluster.csv")
write.csv(table(gut.integrated$annotation_major, gut.integrated$Health),"propotion_major.csv")
write.csv(table(gut.integrated$annotation_major, gut.integrated$Health,gut.integrated$Sample),"propotion_major_sample.csv")
write.csv(table(gut.integrated$annotation1, gut.integrated$Health,gut.integrated$Sample),"propotion_cluster_sample.csv")
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
        geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
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
        geom_pointrange(aes(xmin=log_n-0.5*sn, xmax=log_n+0.5*sn),size = 0.5) +
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
        geom_pointrange(aes(xmin=log_in-0.5*sin, xmax=log_in+0.5*sin),size = 0.5) +
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
        geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
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
png( "cellstates_celltype_inf_boxplot.png", width=500,height=1000,units="px")
print(inf.celltype.boxplot)
dev.off()

non.celltype.boxplot <- ggplot(propotion_cluster, aes(x = log_n, y = Celltype, color = log_n)) +
        geom_point(size = 4) +
        geom_pointrange(aes(xmin=log_n-0.5*sn, xmax=log_n+0.5*sn),size = 0.5) +
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
png( "cellstates_celltype_non_boxplot.png", width=500,height=1000,units="px")
print(non.celltype.boxplot)
dev.off()

inon.celltype.boxplot <- ggplot(propotion_cluster, aes(x = log_in, y = Celltype, color = log_in)) +
        geom_point(size = 4) +
        geom_pointrange(aes(xmin=log_in-0.5*sin, xmax=log_in+0.5*sin),size = 0.5) +
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
png( "cellstates_celltype_innon_boxplot.png", width=500,height=1000,units="px")
print(inon.celltype.boxplot)
dev.off()

#Find DEGs 
DefaultAssay(gut) <- "RNA"
monocytes <- subset(gut,idents = "Infla monocytes")
Idents(monocytes) <- monocytes$Health
c16i_n.markers <- FindMarkers(monocytes, ident.1 = "Inflammed IBD", ident.2 = "Non-Inflammed IBD", only.pos = TRUE)
c16.markers <- FindMarkers(gut, ident.1 = "Infla monocytes", ident.2 = NULL, only.pos = TRUE)
write.csv(c16.markers,"c16 markers.csv")

treg<- subset(gut,idents = "Treg1")
Idents(treg) <- treg$Health
treg_i_h.markers <- FindMarkers(treg, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(treg_i_h.markers,"treg_i_h markers.csv")

B.cell1<- subset(gut,idents = "B Cells1")
Idents(B.cell1) <- B.cell1$Health
B.cell1_i_h.markers <- FindMarkers(B.cell1, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(B.cell1_i_h.markers,"B cell1_i_h markers.csv")

B.cell2<- subset(gut,idents = "B Cells2")
Idents(B.cell2) <- B.cell2$Health
B.cell2_i_h.markers <- FindMarkers(B.cell2, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(B.cell2_i_h.markers,"B cell2_i_h markers.csv")

ilc<- subset(gut,idents = "ILCs")
Idents(ilc) <- ilc$Health
ilc_i_h.markers <- FindMarkers(ilc, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(ilc_i_h.markers,"ILC_i_h markers.csv")

mast69<- subset(gut,idents = "CD69-Mast")
Idents(mast69) <- mast69$Health
mast69_i_h.markers <- FindMarkers(mast69, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(mast69_i_h.markers,"CD69-Mast_i_h markers.csv")

goblet2<- subset(gut,idents = "Goblet2")
Idents(goblet2) <- goblet2$Health
goblet2_i_h.markers <- FindMarkers(goblet2, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(goblet2_i_h.markers,"Goblet2_i_h markers.csv")

cyclingt3<- subset(gut,idents = "Cycling T3")
Idents(cyclingt3) <- cyclingt3$Health
cyclingt3_i_h.markers <- FindMarkers(cyclingt3, ident.1 = "Inflammed IBD", ident.2 = "Healthy", only.pos = TRUE)
write.csv(cyclingt3_i_h.markers,"Cycling T3_i_h markers.csv")