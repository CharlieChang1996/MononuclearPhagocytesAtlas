hb_counts<-read.table(file=paste0("~/Msc project/data/IBDH/HB.txt.gz"),header = T,row.names=1,sep="\t")
row.names(hb_counts) <- gsub(x = row.names(hb_counts), pattern = "_", replacement = "-")
names(hb_counts) <- gsub(x = names(hb_counts), pattern = "\\.", replacement = "-")

s_counts<-read.table(file=paste0("~/Msc project/data/IBDH/S.txt.gz"),header = T,row.names=1,sep="\t")
row.names(s_counts) <- gsub(x = row.names(s_counts), pattern = "_", replacement = "-")
names(s_counts) <- gsub(x = names(s_counts), pattern = "\\.", replacement = "-")

col_counts<-read.table(file=paste0("~/Msc project/data/IBDH/Col.txt.gz"),header = T,row.names=1,sep="\t")
row.names(col_counts) <- gsub(x = row.names(col_counts), pattern = "_", replacement = "-")
names(col_counts) <- gsub(x = names(col_counts), pattern = "\\.", replacement = "-")
hb_counts <- as.sparse(hb_counts)
s_counts <- as.sparse(s_counts)
col_counts <- as.sparse(col_counts)

hb <- CreateSeuratObject(counts = hb_counts, project = "hb",names.field = 2,names.delim = "-",min.cells = 3, min.features = 200)
s <- CreateSeuratObject(counts = s_counts, project = "s",names.field = 2,names.delim = "-",min.cells = 3, min.features = 200)
col <- CreateSeuratObject(counts = col_counts, project = "col",names.field = 2,names.delim = "-",min.cells = 3, min.features = 200)
remove(col_counts)
