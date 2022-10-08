library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(NMF)
library(viridis)

setwd("/lustre/scratch117/cellgen/team205/rl20/kidney_final/object")
### tumour cell object
RCC.integrated <- readRDS("RCC_cells_rmdou_by_patient.integrate.SCT.PCA.combineCNV.RDS")

for (Target_tumor in unique(RCC.integrated$patient)) {
     Subset <- subset(x = RCC.integrated, cells = 
                        rownames(RCC.integrated@meta.data)[which(RCC.integrated$patient == Target_tumor)])
     
     DefaultAssay(Subset) <- 'RNA'
     Subset <- NormalizeData(Subset, normalization.method = "LogNormalize", scale.factor = 10000)
     Subset <- FindVariableFeatures(Subset, selection.method = "vst")
     Subset <- ScaleData(object = Subset, vars.to.regress = "percent.mt")
     data <- as.matrix(Subset@assays$RNA@scale.data)
     data[data<0] <- 0
     
     ########## NMF #####
     set.seed(123)
     rk <- 10
     res <- nmf(data, rank = rk, nrun = 1, seed=1)
     s <- extractFeatures(res, 50L)
     data1 <- c()
     for (i in 1:as.numeric(rk)) {
       data1 <- rbind(data1, as.matrix(Subset@assays$RNA@scale.data)[s[[i]],])
       Subset <- AddModuleScore(object = Subset, features = list(rownames(data)[s[[i]]]), name = paste("Pg", i, sep = "_"))
     }
     
     Pg_mat <- Subset@meta.data[ ,paste('Pg',"_",seq(1:10), rep('1',10), sep="")]
     SD <- apply(Pg_mat, 2, sd)
     data2 <- c()
     for (i in which(SD>0.2)) {
       data2 <- rbind(data2, as.matrix(Subset@assays$RNA@scale.data)[s[[i]],])
     }
     
     ### plot a heatmap
     data2[data2>3] <- 3
     data2[data2< -3] <- -3
     x <- pheatmap(data2, cluster_cols=T, cluster_rows=F, color=my_palette, scale="none",
                   fontsize=5, treeheight_row=0, treeheight_col=0, cellheight = 0.8,cellwidth = 0.1,
                   show_rownames=F, show_colnames=F, clustering_method = "ward.D",
                   border_color=NA, 
                   filename = paste(getwd(), paste(Target_tumor, "pdf", sep='.'), sep='/'))
     
     write.table(data2, file=paste(getwd(), paste(Target_tumor, "_program_mat_50L.xls", sep=''), sep='/'), quote=F, sep="\t")
}

#### Meta-programme
jac <- function(x, y) {
  inter <- intersect(x, y)
  total <- union(x, y)
  similarity <- length(inter)/length(total)
  return(similarity)
}

Pgs <- read.table("/home/jovyan/farm/tm/RCC_final/RCC/Program/program_geneset_all.xls", header=T)
Mat <- matrix(0, ncol = ncol(Pgs), nrow = ncol(Pgs))
rownames(Mat) <- colnames(Pgs)
colnames(Mat) <- colnames(Pgs)

for (i in 1:ncol(Pgs)) {
  for (j in 1:ncol(Pgs)) {
    ss <- jac(as.vector(Pgs[,i]), as.vector(Pgs[,j]))
    Mat[i,j] <- ss*100
  }
}

custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
#custom_magma<-colorRampPalette(c("white","blue"))(n=299)

Mat[Mat>40] <- 40
x <- pheatmap(as.matrix(Mat), cluster_cols=T, cluster_rows=T, 
              clustering_distance_rows="euclidean", color=custom_magma, 
              fontsize=8,treeheight_row=0,treeheight_col=30, 
              cellheight = 7,cellwidth = 7,show_rownames=T, 
              show_colnames=F,clustering_method = "ward.D2", border_color = "NA")
