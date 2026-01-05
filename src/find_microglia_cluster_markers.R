library(Seurat)
library(SeuratObject)

out_dir <- "results/find_microglia_cluster_markers/"
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
seu_subset <- LoadSeuratRds("results/seurat_cluster_microglia_subset/microglia_subset_seurat.RDS")

# clusters
microglia_clusters <- unique(seu_subset@meta.data$microglia_clusters)

Idents(seu_subset) <- "microglia_clusters"

for (cluster in microglia_clusters) {
  
  marker_degs <- FindMarkers(seu_subset,
                             ident.1 = cluster,
                             test.use = "MAST")
  
  marker_degs$gene <- rownames(marker_degs)
  marker_degs$cluster <- cluster
  marker_degs$delta_pct <- marker_degs$pct.1 - marker_degs$pct.2
  
  saveRDS(marker_degs, file=paste0(out_dir, "cluster", cluster, ".markers.RDS"))
}
