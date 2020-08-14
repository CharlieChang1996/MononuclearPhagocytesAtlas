source("/home/s1937334/Msc project/Script/seurat_qc.R")
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(UpSetR)
gut.integrated <- readRDS("~/Msc project/data/downstream/colon_HvD_7_18.rds")

#Dotplots for DEGs across clusters
DefaultAssay(gut.integrated) <- "RNA"
Idents(gut.integrated) <- gut.integrated$annotation_major

uc_dotplot <- DotPlot(subset(gut.integrated), features = uc.genes ,group.by = "annotation1") + RotatedAxis()
png( "uc_related_gene_dotplot.png",width=1050,height=650,units="px")
print(uc_dotplot)
dev.off()

Idents(gut.integrated) <- gut.integrated$annotation1
uc2_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Macrophages2"), 
                       features = c("HSPA1A","HSPA1B","MMP12","HSP90AA1","IGLC2","RPS26","AIF1","DNAJB1") ,
                       group.by = "Health") + RotatedAxis()
png( "DE_ivh_macrophages2_dotplot.png",width=1050,height=650,units="px")
print(uc2_dotplot)
dev.off()

uc3_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "B Cells1"), 
                       features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_cyclingB_dotplot.png",width=1050,height=650,units="px")
print(uc3_dotplot)
dev.off()

uc4_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "CD69+Mast"), 
                       features = c("CPA3","RPL39","MT-ND3","MT-RNR2","IGHG3") ,
                       group.by = "Health") + RotatedAxis()
png( "DE_ivh_cd69Mast_dotplot.png",width=1050,height=650,units="px")
print(uc4_dotplot)
dev.off()

uc5_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Plasma"), 
                       features = c("IGHG3","IGHG1","IGHG2","MT-ND3","NEAT1") ,
                       group.by = "Health") + RotatedAxis()
png( "DE_ivh_Plasma_dotplot.png",width=1050,height=650,units="px")
print(uc5_dotplot)
dev.off()

uc6_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "CD4+ T2"), 
                       features = c("IGHG1","CD82","BATF","CTLA4","SPOCK2",
                                    "IGHA1","SUPT4H1","IGJ","UBA52","RPL32") ,group.by = "Health") + RotatedAxis()
png( "DE_CD4 T2_dotplot.png",width=1050,height=650,units="px")
print(uc6_dotplot)
dev.off()

uc7_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Treg1"), 
                       features = c("RPS15A","MT-COA1","RPL21","RPS4X","RPL34") ,
                       group.by = "Health") + RotatedAxis()
png( "DE_ivh_Treg1_dotplot.png",width=1050,height=650,units="px")
print(uc7_dotplot)
dev.off()

uc8_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Infla monocytes"), features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_inflamonocytes_dotplot.png",width=1050,height=650,units="px")
print(uc8_dotplot)
dev.off()

uc9_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "DC1"), 
                       features = c("MT-ND3","MEF2C","RPS26","SELL","MALAT1") ,
                       group.by = "Health") + RotatedAxis()
png( "DE_ivh_DC1_dotplot.png",width=1050,height=650,units="px")
print(uc9_dotplot)
dev.off()

uc10_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "B Cells2"), 
                        features = c("CD74","B2M","TMSB4X","HLA-DRA","ACTB") ,
                        group.by = "Health") + RotatedAxis()
png( "DE_ivh_Bcells2_dotplot.png",width=1050,height=650,units="px")
print(uc10_dotplot)
dev.off()

uc11_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Memory T"), 
                        features = c("RPSA","RPL21","RPL12","ACTB","RPS4X",
                                     "RPS26","HLA-B","HLA-C","RPS4Y1","MALAT1") ,
                        group.by = "Health") + RotatedAxis()
png( "DE_memoryT_dotplot.png",width=1050,height=650,units="px")
print(uc11_dotplot)
dev.off()

uc12_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "CD8+LP2 T"), 
                        features = c("RPL39","RPS27","RPS29","RPS15A","MALAT1","MT-ND3") ,
                        group.by = "Health") + RotatedAxis()
png( "DE_ivh_CD8LP2 T_dotplot.png",width=1050,height=650,units="px")
print(uc12_dotplot)
dev.off()

uc13_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "Activated Fos-lo T"), 
                        features = c("RPS15A","RPL34","RPL32","RPL21","MT-CO1") ,group.by = "Health") + RotatedAxis()
png( "DE_ivh_foslo T_dotplot.png",width=1050,height=650,units="px")
print(uc13_dotplot)
dev.off()

uc14_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "GC B"), 
                        features = c("RPS2","TMSB4X","RPL13","RPL3","RPL6") ,group.by = "Health") + RotatedAxis()
png( "DE_ivh_GC b_dotplot.png",width=1050,height=650,units="px")
print(uc14_dotplot)
dev.off()

uc15_dotplot <- DotPlot(subset(gut.integrated, downsample = 100,idents = "B Cells1"), 
                        features = c("CD74","HLA-DRA","CD37","LIMD2","HLA-DPB1",
                                     "CLDN4","MT1E","MT1G","RPS17L","CLDN3") ,group.by = "Health") + RotatedAxis()
png( "DE_Bcell1_dotplot.png",width=1050,height=650,units="px")
print(uc15_dotplot)
dev.off()

uc0_dotplot <- DotPlot(gut.integrated, features = uc.genes ,group.by = "Health") + RotatedAxis()
png( "uc_related_health_dotplot.png",width=1050,height=650,units="px")
print(uc0_dotplot)
dev.off()

##Upset plot
#inflamed vs Healthy
cd4t2 <- read_excel("clusters_i_h markers.xlsx",
sheet = "CD4 T2_i_h markers")

mac1ih <- read_excel("clusters_i_h markers.xlsx",
sheet = "Macrophages1_i_h markers")

mac3 <- read_excel("clusters_i_h markers.xlsx",
sheet = "Macrophages3_i_h markers")

mac2 <- read_excel("clusters_i_h markers.xlsx",
sheet = "Macrophages2_i_h markers")

ilc <- read_excel("clusters_i_h markers.xlsx",
sheet = "ilc_i_h markers")

dc1 <- read_excel("clusters_i_h markers.xlsx",
sheet = "DC1_i_h markers")

lp2<- read_excel("clusters_i_h markers.xlsx",
                  sheet = "CD8+LP2 T_i_h markers")

memt <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "Memory_T_i_h markers")

mas <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "CD69-Mast_i_h markers")

mas69ih <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "CD69+Mast_i_h markers")

b1ih <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "B cell1_i_h markers")

foslo <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "Activated Fos-lo T_i_h markers")

treg1 <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "Treg1_i_h markers")

gcb <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "gcb_i_h markers")

gob <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "Goblet_i_h markers")

b2ih <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "B cell2_i_h markers")

plaih <- read_excel("clusters_i_h markers.xlsx",
                  sheet = "plasma_i_h markers")
deglist <- list(Bcell1=b1$Genes,Bcell2=b2$Genes,GC= gcb$Genes,Plasma = pla$Genes,
                CD4_T2 = cd4t2$Genes,Activated_Fos_loT = foslo$Genes,
                Treg1 = treg1$Genes, ILC = ilc$Genes,CD8LP2T = lp2$Genes,
                Memory_T = memt$Genes,Macrophages1 = mac1$Genes,Macrophages2 = mac2$Genes,
                Macrophages3 = mac3$Genes,DC1 = dc1$Genes,Mast = mas$Genes,
                CD69Mast = mas69$Genes,Goblet = gob$Genes)

png("upset_infla_healthy.png",width=900,height=750,units="px")
print(upset(fromList(deglist),nsets = 17,point.size = 3.5,
            text.scale = c(1.5,2),mainbar.y.label = "DEGs infla vs healthy Intersections", 
            sets.x.label = "DEGs Per Cell type",
            order.by = "freq"))
dev.off()

#Non vs healthy upset
cd4t2 <- read_excel("clusters_n_h markers.xlsx",
                    sheet = "Cd4+ T2_n_h markers")

mac1nh <- read_excel("clusters_n_h markers.xlsx",
                   sheet = "Macrophages1_n_h markers")

#mac3 <- read_excel("clusters_n_h markers.xlsx",sheet = "Macrophages3_n_h markers")

mac2 <- read_excel("clusters_n_h markers.xlsx",
                   sheet = "Macrophages2_n_h markers")

ilc <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "ilc_n_h markers")

dc1 <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "DC1_n_h markers")

lp2<- read_excel("clusters_n_h markers.xlsx",
                 sheet = "CD8+LP2 T_n_h markers")

memt <- read_excel("clusters_n_h markers.xlsx",
                   sheet = "Memory_T_n_h markers")

mas <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "CD69-Mast_n_h markers")

mas69nh <- read_excel("clusters_n_h markers.xlsx",
                    sheet = "CD69+Mast_n_h markers")

b1nh <- read_excel("clusters_n_h markers.xlsx",
                 sheet = "B cell1_n_h markers")

foslo <- read_excel("clusters_n_h markers.xlsx",
                    sheet = "Activated Fos-lo T_n_h markers")

treg1 <- read_excel("clusters_n_h markers.xlsx",
                    sheet = "Treg1_n_h markers")

gcb <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "gcb_n_h markers")

gob <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "Goblet_n_h markers")

b2nh <- read_excel("clusters_n_h markers.xlsx",
                 sheet = "B cell2_n_h markers")

planh <- read_excel("clusters_n_h markers.xlsx",
                  sheet = "plasma_n_h markers")

deglist2 <- list(Bcell1=b1$...1,Bcell2=b2$...1,GC= gcb$...1,Plasma = pla$...1,
                CD4_T2 = cd4t2$...1,Activated_Fos_loT = foslo$...1,
                Treg1 = treg1$...1, ILC = ilc$...1,CD8LP2T = lp2$...1,
                Memory_T = memt$...1,Macrophages1 = mac1$...1,Macrophages2 = mac2$...1,
                DC1 = dc1$...1,Mast = mas$...1,
                CD69Mast = mas69$...1,Goblet = gob$...1)

png("upset_non_healthy.png",width=900,height=750,units="px")
print(upset(fromList(deglist2),nsets = 17,point.size = 3.5,
            text.scale = c(1.5,2),mainbar.y.label = "DEGs Non-infla vs healthy Intersections", 
            sets.x.label = "DEGs Per Cell type",
            order.by = "freq"))
dev.off()

#inf vs non
cd4t2 <- read_excel("clusters_i_n markers.xlsx",
                    sheet = "Cd4+ T2_i_n markers")

mac1in <- read_excel("clusters_i_n markers.xlsx",
                   sheet = "Macrophages1_i_n markers")

mac3 <- read_excel("clusters_i_n markers.xlsx",sheet = "Macrophages3_i_n markers")

mac2 <- read_excel("clusters_i_n markers.xlsx",
                   sheet = "Macrophages2_i_n markers")

ilc <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "ilc_i_n markers")

dc1 <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "DC1_i_n markers")

lp2<- read_excel("clusters_i_n markers.xlsx",
                 sheet = "CD8+LP2 T_i_n markers")

memt <- read_excel("clusters_i_n markers.xlsx",
                   sheet = "Memory_T_i_n markers")

mas <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "CD69-Mast_i_n markers")

mas69 <- read_excel("clusters_i_n markers.xlsx",
                    sheet = "CD69+Mast_i_n markers")

b1 <- read_excel("clusters_i_n markers.xlsx",
                 sheet = "B cell1_i_n markers")

foslo <- read_excel("clusters_i_n markers.xlsx",
                    sheet = "Activated Fos-lo T_i_n markers")

treg1 <- read_excel("clusters_i_n markers.xlsx",
                    sheet = "Treg1_i_n markers")

gcb <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "gcb_i_n markers")

gob <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "Goblet_i_n markers")

b2 <- read_excel("clusters_i_n markers.xlsx",
                 sheet = "B cell2_i_n markers")

pla <- read_excel("clusters_i_n markers.xlsx",
                  sheet = "plasma_i_n markers")

deglist3 <- list(Bcell1=b1$...1,Bcell2=b2$...1,GC= gcb$...1,Plasma = pla$...1,
                 CD4_T2 = cd4t2$...1,Activated_Fos_loT = foslo$...1,
                 Treg1 = treg1$...1, ILC = ilc$...1,CD8LP2T = lp2$...1,
                 Memory_T = memt$...1,Macrophages1 = mac1$...1,Macrophages2 = mac2$...1,
                 Macrophages3 = mac3$...1 ,DC1 = dc1$...1,Mast = mas$...1,
                 CD69Mast = mas69$...1,Goblet = gob$...1)

png("upset_inflam_non.png",width=1050,height=750,units="px")
print(upset(fromList(deglist3),nsets = 17,point.size = 3.5,
            text.scale = c(1.5,2),mainbar.y.label = "DEGs Infla vs Non-infla Intersections", 
            sets.x.label = "DEGs Per Cell type",
            order.by = "freq"))
dev.off()

#upset plot of condtions
maclist <- list(Inflamed = mac1ih$Genes,Non_inflamed = mac1nh$...1)
plalist <- list(Inflamed = plaih$Genes,Non_inflamed = planh$...1)
mastlist <- list(Inflamed = mas69ih$Genes,Non_inflamed = mas69nh$...1)
blist <- list(Inflamed_B1 = b1ih$Genes,Non_inflamed_B1 = b1nh$...1,
              Inflamed_B2 = b2ih$Genes,Non_inflamed_B2 = b2nh$...1)
png("upset_macrophages1.png",width=950,height=650,units="px")
upset(fromList(maclist),point.size = 3.5,
      text.scale = c(1.5,2),mainbar.y.label = "DEGs Infla vs Non-infla Intersections", 
      sets.x.label = "DEGs Per condition in Macrophages1",
      order.by = "freq")
dev.off()

png("upset_plasma.png",width=750,height=650,units="px")
upset(fromList(plalist),point.size = 3.5,
      text.scale = c(1.5,2),mainbar.y.label = "DEGs Infla vs Non-infla Intersections", 
      sets.x.label = "DEGs Per condition in Plasma",
      order.by = "freq")
dev.off()

png("upset_CD69+Mast.png",width=750,height=650,units="px")
upset(fromList(mastlist),point.size = 3.5,
      text.scale = c(1.5,2),mainbar.y.label = "DEGs Infla vs Non-infla Intersections", 
      sets.x.label = "DEGs Per condition in CD69+Mast",
      order.by = "freq")
dev.off()

png("upset_Bcell.png",width=750,height=650,units="px")
upset(fromList(blist),point.size = 3.5,
      text.scale = c(1.5,2),mainbar.y.label = "DEGs Infla vs Non-infla Intersections", 
      sets.x.label = "DEGs Per condition in B Cell",
      order.by = "freq")
dev.off()