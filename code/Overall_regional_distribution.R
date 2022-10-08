### import packages
packages <- c("Seurat", "dplyr", "Matrix", 
              "ggplot2", "pheatmap", "RColorBrewer", "viridis",
              "reshape2", "ggdendro","cowplot", "patchwork",
              "ggtree", "tidyverse", "aplot"
)
lapply(packages, library, character.only = TRUE)

### Overall object
setwd("/lustre/scratch117/cellgen/team205/rl20/kidney_final")
ALL <- readRDS("allCellsPostQCAnnotated_broadtype_anno.RDS")
DimPlot(object = ALL, reduction = "umap",pt.size = 0.1, group.by = "annotation")

########### Reginal analysis for all cell types
mat <- c()
rown <- c()
for (ct in unique(ALL$annotation)) {
  re <- ALL$summaryDescription[which(ALL$annotation == ct)] %>% table
  mat <- rbind(mat, re)
  rown <- c(rown, ct)
}
colnames(mat) <- ALL$summaryDescription[which(ALL$annotation=='PT')] %>% table %>% names
rownames(mat) <- rown

write.table(mat, file = paste(getwd(), 'region_mat_all_cell_type.xls', sep='/'),
            quote=F, sep='\t')

mat <- mat[, -which(colnames(mat)=='Metastasis' | colnames(mat)=='Thrombus')] 
mat <- mat[-which(rownames(mat)=='Low quality' | rownames(mat)=='Unknown'),]

x <- chisq.test(mat)
expect <- x$expected
data1 <- mat/expect

data2 <- mat
data2[1:11,] <- data2[1:11,]/sum(data2[1:11,])
data2[12:21,] <- data2[12:21,]/sum(data2[12:21,])
data2[22:37,] <- data2[22:37,]/sum(data2[22:37,])
data2[38:48,] <- data2[38:48,]/sum(data2[38:48,])
data2[49:60,] <- data2[49:60,]/sum(data2[49:60,])
data2[61:69,] <- data2[61:69,]/sum(data2[61:69,])
data2[70:87,] <- data2[70:87,]/sum(data2[70:87,])
data2[88:99,] <- data2[88:99,]/sum(data2[88:99,])
data2[100:105,] <- data2[100:105,]/sum(data2[100:105,])

custom_magma <- c(colorRampPalette(c("grey90", rev(viridis::magma(323, begin = 0.15))[1]))(100), 
                  rev(viridis::magma(323, begin = 0.18)))

x <- pheatmap(data1, cluster_cols=F, cluster_rows=F, color = custom_magma,
              clustering_distance_rows="euclidean",
              fontsize=6,treeheight_row=20,treeheight_col=20, 
              cellheight = 7,cellwidth = 14,show_rownames=T, 
              show_colnames=T,clustering_method = "ward.D", border_color = "NA", scale = 'none')

hec_new <- melt(t(data1),value.name = "roe")
hec_new2 <- melt(t(data2),value.name = "prop")
hec_new3 <- melt(t(mat),value.name = "cell_number")
hec_new$prop <- hec_new2$prop
hec_new$num <- hec_new3$cell_number

dotplot <- hec_new %>% filter(prop > 0.001) %>% 
  ggplot(aes(x=Var1, y = Var2, color = roe, size = prop)) + 
  geom_point()+
  scale_size(range = c(0,6))+
  #scale_size_area(max_size=10,guide=F)+
  scale_color_gradientn(colours = custom_magma, limits = c(0,3),
                        oob = scales::squish, name = 'Roe') +
  scale_y_discrete(position = "right")+
  theme_bw()
dotplot


####### Patient by patient analysis
origin_ident <- ALL@active.ident
ALL@active.ident <- as.factor(ALL$patient)

for (Target_tumor in unique(RCC.integrated$patient)) {
     ALL_patient <- subset(ALL, idents = Target_tumor)
     DimPlot(object = ALL_patient, reduction = "umap",pt.size = 0.1, group.by = "patient")
     mat <- c()
     rown <- c()
     for (ct in unique(ALL$annotation)) {
       re <- ALL_patient$region[which(ALL_patient$annotation == ct)] %>% table()
       mat <- rbind(mat, re)
       rown <- c(rown, ct)
     }
     colnames(mat) <- ALL_patient$region[which(ALL_patient$annotation=='PT')] %>% table %>% names
     rownames(mat) <- rown
     
     write.table(mat, file = paste(getwd(), paste(Target_tumor, "regional_mat.xls", sep=''), sep='/'),
                 quote=F, sep='\t')
     
     ########Ro/e patient by patient
     mat <- mat[,c('a','c','d','e')]
     mat <- mat[-which(rownames(mat)=='Low quality' | rownames(mat)=='Unknown'),]
     
     x <- chisq.test(mat)
     expect <- x$expected
     data1 <- mat/expect
     
     data2 <- mat
     for (i in 1:ncol(data2)) {
       data2[,i] <- data2[,i]/sum(data2[,i])
     }
     
     hec_new <- melt(t(data1),value.name = "roe")
     hec_new2 <- melt(t(data2),value.name = "prop")
     hec_new$prop <- hec_new2$prop
     
     ########## heatmap
     ggplot(hec_new,aes(x=Var1,y=Var2,fill=roe))+geom_tile()+
       #geom_text(aes(fill = hec_new$prop, label = round(hec_new$prop, 2)))+
       scale_fill_gradient2(high = "darkblue", 
                            mid = "orange2", 
                            low = "white", 
                            midpoint = 5) + 
       theme(panel.grid.major.x=element_blank(), #no gridlines
             panel.grid.minor.x=element_blank(), 
             panel.grid.major.y=element_blank(), 
             panel.grid.minor.y=element_blank(),
             panel.background=element_rect(fill="white"), # background=white
             axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
             plot.title = element_text(size=20,face="bold"),
             axis.text.y = element_text(size = 12,face = "bold")
       ) + 
       ggtitle("heatmap Plot")+
       theme(legend.title=element_text(face="bold", size=14)) + 
       scale_x_discrete(name="") +
       scale_y_discrete(name="") +
       labs(fill="heatmap legend")
}


########## patient by patient heter combined dotplot

path <- paste(getwd(), 'revision/regional_difference', sep='/')
sample_list <- c('region_mat_PD43824.xls','region_mat_PD45814.xls','region_mat_PD45815.xls',
                 'region_mat_PD45816.xls','region_mat_PD47171.xls','region_mat_PD47172.xls',
                 'region_mat_PD47512.xls')
Data1 <- matrix()
Data2 <- matrix()
Mat <- matrix()
for (sample in sample_list) {
  mat <- read.table(paste(path, sample, sep='/'), sep='\t')
  mat <- mat[,c('a','c','d','e')]
  mat <- mat[-which(rownames(mat)=='Low quality' | rownames(mat)=='Unknown'),]
  colnames(mat) <- paste(sample, colnames(mat), sep='_')
  
  x <- chisq.test(mat)
  expect <- x$expected
  data1 <- mat/expect
  Data1 <- cbind(Data1, data1)
  
  data2 <- mat
  for (i in 1:ncol(data2)) {
    data2[,i] <- data2[,i]/sum(data2[,i])
  }
  Data2 <- cbind(Data2, data2)
  
  Mat <- cbind(Mat, mat)
}

Data1 <- Data1[,-1]
Data2 <- Data2[,-1]
Mat <- Mat[,-1]

hec_new <- melt(t(Data1),value.name = "roe")
hec_new2 <- melt(t(Data2),value.name = "prop")
hec_new$prop <- hec_new2$prop

custom_magma <- c(colorRampPalette(c("grey90", rev(viridis::magma(323, begin = 0.15))[1]))(100), 
                  rev(viridis::magma(323, begin = 0.18)))
dotplot <- hec_new %>% 
  #mutate(Var1 = factor(Var1, levels = clust$labels[clust$order])) %>% 
  ggplot(aes(x=Var1, y = Var2, color = roe, size = prop)) + 
  geom_point()+
  scale_size(range = c(1,9),limits = c(0.001,max(Data2)))+
  #scale_size_area(max_size=10,guide=F)+
  scale_color_gradientn(colours = custom_magma, limits = c(0,3),
                        oob = scales::squish, name = 'Roe') +
  scale_y_discrete(position = "right")+
  theme_bw()
dotplot
