library(Seurat)
library(dplyr)
library(patchwork)
library(umap)
ibd_counts<-read.table(file=paste0("~/Msc project/data/IBDP/IBDP.txt.gz"),header = T,row.names=1,sep="\t")
row.names(ibd_counts) <- gsub(x = row.names(ibd_counts), pattern = "_", replacement = "-")
names(ibd_counts) <- gsub(x = names(ibd_counts), pattern = "_", replacement = "-")
ibd_counts <- as.sparse(ibd_counts)
ibd <- CreateSeuratObject(counts = ibd_counts, project = "IBD_P",names.field = 2,names.delim = "\\.",min.cells = 3, min.features = 200)
remove(ibd_counts)
saveRDS(ibd,paste0("ibdp_raw.rds") )