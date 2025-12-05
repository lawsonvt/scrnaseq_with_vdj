library(Seurat)

out_dir <- "results/join_and_filter_integrated_seurat/"
dir.create(out_dir, showWarnings = F, recursive = T)

# read in integrated seurat
int_seu <- LoadSeuratRds("results/integrate_seurat_samples/all_samples.integrated_seurat.RDS")

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

# DROP CLUSTER 16 (see cluster_marker_tests and inspect_cluster_marker_test for why)
int_seu <- subset(int_seu, subset = harmony_clusters != 16)

# normalize the data
int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

SaveSeuratRds(int_seu, paste0(out_dir, "filtered.integrated_seurat.RDS"))