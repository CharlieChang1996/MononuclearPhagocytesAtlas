source("/home/s1937334/Msc project/Script/seurat_qc.R")
library(RColorBrewer)
library(ggplot2)

inflamono.boxplot <- ggplot(box_my, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  ylim(0,0.5)+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin inflammatory monocytes") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_inflamonocytes_boxplot.png", width=400,height=600,units="px")
print(inflamono.boxplot)
dev.off()

inflamono.boxplot2 <- ggplot(box_my, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin inflammatory monocytes") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_inflamonocytes_boxplot2.png", width=400,height=600,units="px")
print(inflamono.boxplot2)
dev.off()

macro.boxplot <- ggplot(box_mac, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin Macrophages1") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_macrophage_boxplot.png", width=400,height=600,units="px")
print(macro.boxplot)
dev.off()

cymono.boxplot <- ggplot(box_cm, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin Macrophages2") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_macrophages2_boxplot.png", width=400,height=600,units="px")
print(cymono.boxplot)
dev.off()

dc1.boxplot <- ggplot(box_dc1, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin DC1") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5)+
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_DC1_boxplot.png", width=400,height=600,units="px")
print(dc1.boxplot)
dev.off()

dc2.boxplot <- ggplot(box_dc2, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin DC2") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_DC2_boxplot.png", width=400,height=600,units="px")
print(dc2.boxplot)
dev.off()

mc3.boxplot <- ggplot(box_mc3, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in Myeloid (%)", title = "Box plot\nProportions of health-conditions\nin Macrophages3") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_macophages3_boxplot.png", width=400,height=600,units="px")
print(mc3.boxplot)
dev.off()

#em
ma6.boxplot <- ggplot(box_6ma, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin CD69+ Mast") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_cd69mast_boxplot.png", width=400,height=600,units="px")
print(ma6.boxplot)
dev.off()

mast.boxplot <- ggplot(box_ma, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin CD69- Mast") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_mast_boxplot.png", width=400,height=1000,units="px")
print(mast.boxplot)
dev.off()

mast.boxplot2 <- ggplot(box_ma, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin CD69- Mast") +
  ylim(0,0.3)+
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_mast_boxplot2.png", width=400,height=600,units="px")
print(mast.boxplot2)
dev.off()

gob.boxplot <- ggplot(box_gob, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Goblets") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_goblet_boxplot.png", width=400,height=600,units="px")
print(gob.boxplot)
dev.off()

epi.boxplot <- ggplot(box_epi, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Epithelias") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_epi_boxplot.png", width=400,height=600,units="px")
print(epi.boxplot)
dev.off()

epi2.boxplot <- ggplot(box_epi2, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Epithelias2") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_epi2_boxplot.png", width=400,height=600,units="px")
print(epi2.boxplot)
dev.off()

#b
b.boxplot <- ggplot(box_b, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin B cells") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_bcell_boxplot.png", width=400,height=600,units="px")
print(b.boxplot)
dev.off()

cb.boxplot <- ggplot(box_cb, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Cycling B cells") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_cyclingb_boxplot.png", width=400,height=600,units="px")
print(cb.boxplot)
dev.off()

fb.boxplot <- ggplot(box_fb, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Follicular B cells") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_follicular_boxplot.png", width=400,height=600,units="px")
print(fb.boxplot)
dev.off()

gcb.boxplot <- ggplot(box_gcb, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin GC B cells") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_GC_boxplot.png", width=400,height=600,units="px")
print(gcb.boxplot)
dev.off()

p.boxplot <- ggplot(box_p, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion (%)", title = "Box plot\nProportions of health-conditions\nin Plasma") +
#geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_plasma_boxplot.png", width=400,height=600,units="px")
print(p.boxplot)
dev.off()


#T
treg1.boxplot <- ggplot(box_treg1, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in T cells(%)", title = "Box plot\nProportions of health-conditions\nin Treg1") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_treg1_boxplot.png", width=400,height=600,units="px")
print(treg1.boxplot)
dev.off()

t3.boxplot <- ggplot(box_t3, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in T cells(%)", title = "Box plot\nProportions of health-conditions\nin CD4+ T2") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_CD4T2_boxplot.png", width=400,height=600,units="px")
print(t3.boxplot)
dev.off()

memot.boxplot <- ggplot(box_mt, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in T cells(%)", title = "Box plot\nProportions of health-conditions\nin Memory T") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_Memory T_boxplot.png", width=400,height=600,units="px")
print(memot.boxplot)
dev.off()

lp2.boxplot <- ggplot(box_lp2, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in T cells(%)", title = "Box plot\nProportions of health-conditions\nin CD8+ LP2+ T") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_CD8 LP2 T_boxplot.png", width=400,height=600,units="px")
print(lp2.boxplot)
dev.off()

lo.boxplot <- ggplot(box_t, aes(x = Health, y = Freq,color = Health)) +
  geom_boxplot(outlier.colour = NA) + 
  geom_jitter(size =3, position=position_jitter(0.2))+
  labs( x="condition",y = "Proportion in T cells(%)", title = "Box plot\nProportions of health-conditions\nin Activated Fos-lo T") +
  #geom_pointrange(aes(xmin=log_i-0.5*si, xmax=log_i+0.5*si),size = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position="none") 
png( "cellstates_Activated Fos-lo T_boxplot.png", width=400,height=600,units="px")
print(lo.boxplot)
dev.off()