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
  labs(title="Microglia Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))
ggsave(paste0(out_dir, "microglia_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

# 20 dimensions makes sense for this dataset as well

max_pc_dim <- 20

# cluster the harmonized data
subset_seu <- FindNeighbors(subset_seu, dims = 1:max_pc_dim, reduction = "pca")
subset_seu <- FindClusters(subset_seu, cluster.name = "microglia_clusters")

# create umap
subset_seu <- RunUMAP(subset_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.microglia_pca")

# plot it
DimPlot(subset_seu, reduction="umap.microglia_pca", group.by= c("orig.ident","microglia_clusters"))
ggsave(paste0(out_dir, "microglia_subset_seurat.umap.png"), width=14, height=6)

# save the subset seurat
SaveSeuratRds(subset_seu, paste0(out_dir, "microglia_subset_seurat.RDS"))


