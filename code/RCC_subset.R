library(dplyr)
library(Seurat)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(sctransform)
library(EnhancedVolcano)
library(viridis)

setwd("/lustre/scratch117/cellgen/team205/rl20/kidney_final/object")
RCC.integrated <- readRDS("RCC_cells_rmdou_by_patient.integrate.SCT.PCA.combineCNV.RDS")

RCC.integrated <- FindNeighbors(object = RCC.integrated, dims = 1:30)
RCC.integrated <- FindClusters(object = RCC.integrated, resolution = 0.5)
RCC.integrated <- RunUMAP(object = RCC.integrated, dims = 1:30)

Col_1 <- c("DarkKhaki", "Yellow2","YellowGreen","IndianRed","Thistle2",
           "LightSlateGray","DarkSeaGreen1","CornflowerBlue","LightSlateBlue",
           "DeepSkyBlue","Goldenrod2","Tomato","Honeydew2","SlateGray3",
           "Plum","PaleVioletRed1","Thistle2","SlateGray2","grey71","NavajoWhite1",
           "Cyan2","LightSalmon1","DarkSeaGreen","#6B8E23","grey","grey")
DimPlot(object = RCC.integrated, reduction = "umap",pt.size = 0.4,label=T, cols = Col_1[2:20])
DimPlot(object = RCC.integrated, reduction = "umap",pt.size = 0.4,group.by = "annotation")
DimPlot(object = RCC.integrated, reduction = "umap",pt.size = 0.4, group.by = "summaryDescription", split.by = "summaryDescription")

DefaultAssay(RCC.integrated) <- "SCT"
markers <- FindAllMarkers(object = RCC.integrated, only.pos = TRUE, 
                          min.pct = 0.1, logfc.threshold = 0.25)
x <- markers %>% group_by(cluster) %>% top_n(n =100, wt = avg_log2FC)
write.table(x, file="/home/jovyan/farm/tm/RCC_final/RCC/RCC-DEGs-top100-r0.4.xls",quote=F, sep="\t")
top30 <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DoHeatmap(object = RCC.integrated, features = top20$gene) +
  scale_fill_gradientn(colors = c("steelblue1", "white", "tomato"))

cluster.averages <- AverageExpression(RCC.integrated, return.seurat = TRUE)
DoHeatmap(cluster.averages, features = top30$gene, size = 3, 
          draw.lines = FALSE)+  scale_fill_gradientn(colors = c("steelblue2", "white", "tomato2"))

RCC.integrated@active.ident <- factor(RCC.integrated$annotation)
                                      
features <- c('VEGFA','ANGPTL4','CPT2','PPARA','CPT1A','PRKAA2','PDK2','PRKAB1',
              'CDK2','CDK4','CDK6','BUB1','CCNE1','MKI67','CYP4F3','CYP8B1','NNMT','MGST1',
              'C1S','C1R','CFB','C3')
DotPlot(RCC.integrated, features = features) + RotatedAxis() +
  scale_colour_gradient2(low = "dodgerblue2", mid = "grey90", high = "red",midpoint = 0)

custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
FeaturePlot(object = RCC.integrated, features = features,
            pt.size=0.3,reduction = "umap", ncol=3) &
  scale_colour_gradientn(colours = custom_magma)

############################# current genes for scoring
DefaultAssay(RCC.integrated) <- "SCT"
PT <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_PT.txt')
CD <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_cell_death.txt')
EMT <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_EMT.txt')
Cycling <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_cell_cycle.txt')
MHC <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_MHCII.txt')
stress <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_Stress.txt')

######### Motzer score
F1 <- list(c("CD8A","IFNG","EOMES","PRF1","CD274"))
F2 <- list(c("VEGFA","KDR","ESM1","CD34",'PECAM1','ANGPTL4'))
F3 <- list(c("CPT2","PPARA","CPT1A","PRKAA2","PDK2","PRKAB1"))
F4 <- list(c("CDK2","CDK4","CDK6","BUB1","BUB1B","CCNE1","POLQ","AURKA","MKI67","CCNB2"))
F5 <- list(c("FASN","PARP1","ACACA","G6PD","TKT","TALDO1","PGD"))
F6 <- list(c("FAP","FN1","COL5A1","COL5A2","POSTN",'COL1A1','COL1A2','MMP2'))
F7 <- list(c("CXCL1","CXCL2","CXCL3","CXCL8","IL6",'PTGS2'))
F8 <- list(c("F2","C1S","C9","C1R","CFB",'C3'))
F9 <- list(c("CYP4F3","CYP8B1","NNMT","MGST1","MAOA",'CYP4F11','CYP4F2','CYP4F12'))

########### 16 genes and clearcode 34
Rini16_Va <- list(c("APOLD1","EDNRB","NOS3","PLPP3"))
Rini16_Ce <- list(c("EIF4EBP1","TUBB2A","LMNB1"))
Rini16_Ir <- list(c("CEACAM1","CX3CL1","CCL5"))
Rini16_Ref <- list(c("AAMP","ARF1","ATP5E",'GPX1','RPLP1'))
cc34_ccA <- list(c("MAPT","STK32B","FZD1",'RGS5','GIPC2','PDGFD','EPAS1','MAOB','CDH5',
                   'TCEA3','LEPROTL1','BNIP3L','EHBP1','VCAM1','PHYH','PRKAA2','SLC4A4',
                   'ESD','TLR3','NRP1','C11orf1','ST13','ARNT','SPRYD7'))
cc34_ccB <- list(c("SERPINA3","SLC4A3","MOXD1",'KCNN4','ROR2','FLJ23867','FOXM1','UNG',
                   'GALNT10','GALNT4'))

#################### Scoring
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = PT, name = "PT_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = EMT, name = "EMT_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = CD, name = "CD_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = Cycling, name = "Cycling_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = MHC, name = "MHC_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = stress, name = "Stress_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F1, name = "F1_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F2, name = "F2_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F3, name = "F3_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F4, name = "F4_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F5, name = "F5_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F6, name = "F6_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F7, name = "F7_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F8, name = "F8_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = F9, name = "F9_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = Rini16_Va, name = "Rini16_Va_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = Rini16_Ce, name = "Rini16_Ce_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = Rini16_Ir, name = "Rini16_Ir_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = Rini16_Ref, name = "Rini16_Ref_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = cc34_ccA, name = "ccA_score")
RCC.integrated <- AddModuleScore(object = RCC.integrated, features = cc34_ccB, name = "ccB_score")

FeaturePlot(object = RCC.integrated, features = c("EMT_score1", "PT_score1","CD_score1","Cycling_score1",
                                                  "MHC_score1","Stress_score1",'NMF2_score1'),
            pt.size=0.3,reduction = "umap", ncol=3)

custom_magma <- c(colorRampPalette(c("grey90", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))
FeaturePlot(object = RCC.integrated, features = c("EMT_score1", "PT_score1","CD_score1","Cycling_score1",
                                                  "MHC_score1","Stress_score1",'F1_score1','F2_score1',
                                                  'F3_score1','F4_score1','F5_score1','F6_score1','F7_score1',
                                                  'F8_score1','F9_score1'),
            pt.size=0.3,reduction = "umap", ncol=3) &
  scale_colour_gradientn(colours = custom_magma)

RCC.integrated@active.ident <- factor(RCC.integrated$annotation)
RCC.integrated_sub <- subset(RCC.integrated, idents = c('CD','Cycling','EMT',
                                                        'MHCII','PT','Stress'))

features <- c("CD_score1","Cycling_score1","EMT_score1", "MHC_score1",
              "PT_score1","Stress_score1",'F1_score1','F2_score1',
              'F3_score1','F4_score1','F5_score1','F6_score1','F7_score1',
              'F8_score1','F9_score1','Rini16_Va_score1',
              'Rini16_Ce_score1','Rini16_Ir_score1','ccA_score1','ccB_score1')
#VlnPlot(RCC.integrated, features = features,pt.size = 0.1)
DotPlot(RCC.integrated_sub, features = features) + RotatedAxis() + 
  scale_colour_gradient2(low = "dodgerblue2", mid = "grey90", high = "red",midpoint = 0)

TT <- RCC.integrated@meta.data[,c('EMT_score1','PT_score1','Cycling_score1','Stress_score1','CD_score1','MHC_score1','summaryDescription')]
TT <- TT[order(TT$EMT_score1, decreasing = T),]
tt <- TT[,c(1,2)]

library(pheatmap)
library(RColorBrewer)
my_palette<-colorRampPalette(c("dodgerblue","white","red"))(n=299)
x<-pheatmap(as.matrix(tt),cluster_rows=F, cluster_cols=F,color=my_palette,
            scale="column",fontsize=5,treeheight_row=0,treeheight_col=0, 
            cellheight = 0.02,cellwidth = 20,show_rownames=F,border_color=NA)


rr <- matrix(0, ncol=3, nrow=nrow(TT))
colnames(rr) <- c('interface','core','other')

for (i in 1:nrow(TT)) {
  if (TT[i, 'summaryDescription'] == "Tumour-normal") {rr[i,'interface'] <- 1
  } else if (TT[i,'summaryDescription']=="Tumour") {rr[i,'core'] <- 1
  } else {rr[i,'other'] <- 1
  }
}

my_palette <- c("white","black")
x<-pheatmap(as.matrix(rr),cluster_rows=F,cluster_cols=F,color=my_palette,
            scale="none",fontsize=5,treeheight_row=0,treeheight_col=0, 
            cellheight = 0.02,cellwidth = 20,show_rownames=F,border_color=NA)

