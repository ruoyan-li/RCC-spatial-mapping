library(nichenetr)
library(Seurat)
library(tidyverse)

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

setwd("/lustre/scratch117/cellgen/team205/rl20/kidney_final/object")

Macro_tissue <- readRDS("Macro_tissue.RDS") ### Macrophages
Macro_tissue$celltype <- rep("Macro", nrow(Macro_tissue@meta.data))
Macro_tissue@active.ident <- factor(Macro_tissue$celltype)
expressed_genes_sender <- get_expressed_genes("Macro", Macro_tissue, pct = 0.05, assay_oi = "SCT")

RCC.integrated <- readRDS("RCC_cells_rmdou_by_patient.integrate.SCT.PCA.combineCNV.RDS") ### tumour cells
EMT <- subset(RCC.integrated, cells = rownames(RCC.integrated@meta.data)[which(RCC.integrated$annotation=="EMT")])
EMT$celltype <- rep("EMT", nrow(EMT@meta.data))
EMT@active.ident <- factor(EMT$celltype)
expressed_genes_receiver <- get_expressed_genes("EMT", EMT, pct = 0.05, assay_oi = "SCT")

geneset_oi <- read.table('/home/jovyan/farm/tm/RCC_final/RCC/Program/Sig_EMT.txt')
geneset_oi <- as.vector(geneset_oi[,1])
head(geneset_oi)

background_expressed_genes = intersect(rownames(ligand_target_matrix), expressed_genes_receiver)
length(background_expressed_genes)

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands, expressed_genes_sender)
length(expressed_ligands)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors, expressed_genes_receiver)
length(expressed_receptors)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
head(ligand_activities %>% arrange(-pearson))

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi, 
                                                                 ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

### Plot
pdf(file=paste(getwd(), 'Macro_NichNet.pdf', sep='/'), w=8,h=8)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Macro-ligands","EMT genes in malignant cells", 
                                                                    color = "purple",legend_position = "top", x_axis_position = "top",
                                                                    legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "red", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network

dev.off()

