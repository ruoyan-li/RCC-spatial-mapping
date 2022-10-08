library(infercnv)
library(Seurat)
library(Matrix)

RCC <- readRDS("/lustre/scratch117/cellgen/team205/rl20/kidney_final/object/NonImmune-for-CNV.RDS")

gene_list <- read.table("/home/jovyan/farm/tm/RCC_final/RCC/CNV/gencode_v21_gen_pos.complete.match.mat.order.intersect.txt")
counts_matrix <- RCC@assays$RNA@counts[as.vector(gene_list[,1]), ]


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="/home/jovyan/farm/tm/RCC_final/RCC/CNV/dataforCNV_meta.txt",
                                    delim="\t",
                                    gene_order_file="/home/jovyan/farm/tm/RCC_final/RCC/CNV/gencode_v21_gen_pos.complete.match.mat.order.intersect.txt",
                                    ref_group_names=c("Epithelial","Endothelial","Fibroblast")) 

infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/lustre/scratch117/cellgen/team205/rl20/kidney_final/infercnv_out_1", 
                             cluster_by_groups=TRUE,
			                 analysis_mode='subclusters', 
                             denoise=TRUE,
                             num_threads=6,
			     HMM=TRUE)