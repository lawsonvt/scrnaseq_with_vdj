library(Seurat)
library(SeuratObject)
library(ggplot2)

out_dir <- "results/seurat_cluster_microglia_subset/"
dir.create(out_dir, showWarnings = F)


# load in seurat object
int_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

# subset down to microglia
subset_seu <- subset(int_seu, subset = cell_cluster == "Microglia")

# determine the optimal number of dimensions through an elbow plot

subset_seu <- RunPCA(subset_seu, npcs = 50)

# inspect elbow plot
ElbowPlot(subset_seu, ndims=50) + 
  labs(title="Merged Samples (not integrated)") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")






