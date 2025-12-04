library(Seurat)

out_dir <- "results/cluster_marker_tests/"
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
int_seu <- LoadSeuratRds("results/integrate_seurat_samples/all_samples.integrated_seurat.RDS")

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

# find differences for some clusters

Idents(int_seu) <- "harmony_clusters"

# the one we are interested in
c16_degs <- FindMarkers(int_seu,
                        ident.1 = "16",
                        test.use="wilcox")

# write to file
saveRDS(c16_degs, paste0(out_dir, "cluster16_degs.wilcox.RDS"))

# known ones to double check
c22_degs <- FindMarkers(int_seu,
                        ident.1 = "22",
                        test.use="wilcox")

saveRDS(c22_degs, paste0(out_dir, "cluster22_degs.wilcox.RDS"))


c18_degs <- FindMarkers(int_seu,
                        ident.1 = "18",
                        test.use="wilcox")

saveRDS(c18_degs, paste0(out_dir, "cluster18_degs.wilcox.RDS"))



