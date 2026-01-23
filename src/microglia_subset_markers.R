library(Seurat)
library(SeuratObject)

out_dir <- "results/microglia_subset_markers/"
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_subset <- LoadSeuratRds("results/seurat_cluster_microglia_subset/microglia_subset_seurat.RDS")

clusters <- unique(Idents(seu_subset))

# subsample to make this go quicker
set.seed(42)

sampled_cells <- sample(colnames(seu_subset), size=12000, replace = F)

seu_subset_sub <- subset(seu_subset, cells=sampled_cells)

cluster_group <- c("1","16","13","2")

cluster_group_markers <- lapply(cluster_group, function(cluster) {

  print(cluster)

  cluster_markers <- FindMarkers(seu_subset_sub,
                                  ident.1 = cluster,
                                  ident.2 = clusters[!clusters %in% cluster_group],
                                test.use="MAST")
  cluster_markers$gene <- rownames(cluster_markers)
  cluster_markers$pct_delta <- cluster_markers$pct.1 -
    cluster_markers$pct.2
  cluster_markers$cluster <- cluster

  return(cluster_markers)
})
names(cluster_group_markers) <- cluster_group

saveRDS(cluster_group_markers, file=paste0(out_dir, "cluster_group_markers.RDS"))
