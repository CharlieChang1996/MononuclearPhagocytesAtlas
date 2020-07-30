source("/home/s1937334/Msc project/Script/seurat_qc.R")
library(RColorBrewer)
library(ggplot2)
gut <- readRDS("~/Msc project/data/downstream/colon_HvD_7_18.rds")

write.csv(table(gut.integrated$annotation1, gut.integrated$Health),"propotion_cluster.csv")
percentge <- read_csv("propotion_cluster.csv")
percentage_celltype <- read_csv("statistics/percentage_celltype.csv")
color1 <- brewer.pal(11,"BrBG")
## visualisation of proportion of cells
png("proportion_celltype1.png",width = 980)
barplot(cbind(ph[1:10],pi[1:10],pn[1:10]) ~ X1[1:10], data = percentage_celltype,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflammed","Non-inflammed"),
       fill = color1)
dev.off()
png("proportion_celltype2.png",width = 1080)
barplot(cbind(ph[11:20],pi[11:20],pn[11:20]) ~ X1[11:20], data = percentage_celltype,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflammed","Non-inflammed"),
       fill = color1)
dev.off()
png("proportion_celltype3.png",width = 980)
barplot(cbind(ph[21:30],pi[21:30],pn[21:30]) ~ X1[21:30], data = percentage_celltype,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflammed","Non-inflammed"),
       fill = color1)
dev.off()
png("proportion_celltype4.png",width = 780)
barplot(cbind(ph[31:35],pi[31:35],pn[31:35]) ~ X1[31:35], data = percentage_celltype,
        xlab = "Cluster",
        ylab = "Percentage(%)",
        ylim = c(0,140),
        col=color1)
legend("topright",
       c("Healthy","Inflammed","Non-inflammed"),
       fill = color1)
dev.off()

inf.cluster.boxplot <- ggplot(percentage_celltype, aes(x = log_i, y = X1, fill = X1)) +
        geom_point() +
        geom_vline(xintercept = 0.45300041, linetype="dashed") +
        geom_segment(aes(x=0, y=12.5, xend=-3, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.91, y=12.5, xend=3.91, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-7.5, 7.5)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n inflammed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12), 
              axis.title.x = element_text(size = 16, vjust = -1),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_boxplot.png", width=500,height=1000,units="px")
print(inf.cluster.boxplot)
dev.off()

non.cluster.boxplot <- ggplot(percentage_celltype, aes(x = log_n, y = X1, fill = X1)) +
        geom_point() +
        geom_vline(xintercept = -0.19096015, linetype="dashed") +
        geom_segment(aes(x=-0.7, y=12.5, xend=-3.7, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) + 
        geom_segment(aes(x=0.3, y=12.5, xend=3.3, yend=12.5), lineend = "round", linejoin = "mitre", arrow=arrow(), size=0.1) +
        scale_x_continuous(limits = c(-7.5, 7.5)) +
        theme_classic() +
        labs(x = "log2(relative proportion of cell health states\n Non-inflammed/healthy)", y = "Cell type", title = " Decreasing | Increasing\n cell health states proportion in different cell types") +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12), 
              axis.title.x = element_text(size = 16, vjust = -1),
              axis.title.y = element_text(size = 16),
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.position="none") 
png( "cellstates_non_boxplot.png", width=500,height=1000,units="px")
print(non.cluster.boxplot)
dev.off()


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