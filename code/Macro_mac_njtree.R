library(phangorn)
library(ape)

cell <- c("Classical-Mono.1", "Classical-Mono.2", "Classical-Mono.3",
          "Classical-Mono.4", "Non-classical-Mono", "TR-Mac.3",
          "An-pres.-TAM", "FN1+-TAM", "GPNMB+-TAM", "TR-Mac.2",
          "TR-Mac.1", "Pro-infla.-TAM", "SPP1+-TAM", "RGS1+-TAM")

setwd("Users/bio-x/Desktop/Teichlab/kidney_cancer/Final/Mye-tree-final")
ff <- file("Mac_Mono-all-rmDEGs-individualCells-annotation-summary.txt")
Lines <- readLines(ff)
close(ff)

mat <- matrix(ncol = 14, nrow = length(Lines), rep(0, length(Lines)*14))
colnames(mat) <- cell

Rnames <- c()
for (n in 1:length(Lines)) {
  Rnames <- c(Rnames, unlist(strsplit(unlist(strsplit(Lines[n],'\t'))[1], "/")))
  cell_type <- unlist(strsplit(unlist(strsplit(Lines[n],'\t'))[2], "/"))
  mat[n, cell_type] <- mat[n, cell_type] + 1
}
rownames(mat) <- Rnames

library(pheatmap)
library(RColorBrewer)

x <- pheatmap(mat, cluster_cols=T, cluster_rows=T, 
              clustering_distance_rows="euclidean", color=c("white", "grey50"), 
              fontsize=8,treeheight_row=30,treeheight_col=30, 
              cellheight = 2,cellwidth = 12,show_rownames=F, 
              show_colnames=T,clustering_method = "ward.D", border_color = "NA")

################ neighbor-joining and boot.phylo in ape #######
dm<-dist(t(mat), method = "binary") ## binary distance
tree <- NJ(dm)
plot(root(tree, outgroup=4))
set.seed(111)
bstrees <- boot.phylo(tree, t(mat), B = 100, FUN = function(x) NJ(dist(x, method = "binary")), trees = TRUE)$trees
treeNJ <- plotBS(root(tree, outgroup=1), bstrees, "phylogram", p=0)

##########hamming bootstraping ######
phy<-phyDat(t(mat),type="USER", levels=c(0,1))
dm <- dist.hamming(phy)
tree <- NJ(dm)
set.seed(111)
NJtrees <- bootstrap.phyDat(phy,
                            FUN=function(x) NJ(dist.hamming(x)), bs=100)
treeNJ <- plotBS(root(tree, outgroup=1), NJtrees, "phylogram", p=0)