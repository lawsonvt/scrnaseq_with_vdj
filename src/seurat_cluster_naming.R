library(Seurat)
library(ggplot2)
library(scRepertoire)

out_dir <- "results/seurat_cluster_naming/"
dir.create(out_dir, showWarnings = F)

# cluster markers as supplied by Lukens Lab

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"))

# additional marker genes taken from ssreads list
# source: https://bmblx.bmi.osumc.edu/ssread/help/methods#marker-genes

ssread_cluster_markers <- list("Astrocytes"=c("Gfap", "Aqp4", "Gja1", "Slc1a2","Fgfr3",
                                              "Nkain4", "Agt", "Plxnb1", "Slc1a3"),
                               "Endothelial Cells"=c("Cldn5", "Vwf"),
                               "Neurons"=c("Gls", "Rbfox3", "Camk2a"),
                               "Microglia"=c("P2ry12", "Csf1r", "Cx3cr1", "C3"),
                               "Oligodendrocytes"=c("Mbp", "Mobp", "Plp1", "Myrf", "Mag"),
                               "Oligodendrocyte Precursor Cells"=c("Vcan", "Sox8"),
                               "Pericytes"=c("Ambp", "Higd1b", "Pth1r"),
                               "Excitatory Neurons"=c("Slc17a6", "Slc17a7", "Satb2"),
                               "Inhibitory Neurons"=c("Gad1", "Gad2"))

# read in integrated seurat
int_seu <- LoadSeuratRds("results/integrate_seurat_samples/all_samples.intergrated_seurat.RDS")

analysis_genes <- rownames(int_seu)

# find any issues between marker genes and analysis genes

lab_cluster_markers_missing <- lapply(lab_cluster_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

ssread_cluster_markers_missing <- lapply(ssread_cluster_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})
# missing genes seem to be just missing and not typos


