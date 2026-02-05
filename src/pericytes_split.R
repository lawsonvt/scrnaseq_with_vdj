library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(scDblFinder)
library(BiocParallel)
library(clustree)

out_dir <- "results/pericytes_split/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming/cell_named.seurat.RDS")

# subset to the pericyte ones
peri_seu <- subset(seu_obj, subset = cell_cluster == "Pericytes")

peri_seu <- RunPCA(peri_seu, npcs = 20)

# inspect elbow plot
ElbowPlot(peri_seu, ndims=20) + 
  labs(title="TCell Subset")
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 20

# cluster the harmonized data
peri_seu <- FindNeighbors(peri_seu, dims = 1:max_pc_dim, reduction = "pca")
peri_seu <- FindClusters(peri_seu, cluster.name = "peri_clusters")

peri_seu <- RunUMAP(peri_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.peri_pca")

DimPlot(peri_seu, reduction="umap.peri_pca", group.by= c("orig.ident","peri_clusters"))

peri_seu <- FindClusters(peri_seu, resolution = c(0.3, 0.5, 0.8, 1.0, 1.2))

clustree(peri_seu, prefix = "RNA_snn_res.")
ggsave(paste0(out_dir, "res_clustering_tree.png"), width=9, height=8)

res_cols <- grep("RNA_snn_res", colnames(peri_seu@meta.data), value=T)

markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                "Monocytes"=c("Ccr2", "Cd44"),
                "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                "T CD"=c("Cd4",
                         "Cd8a",
                         "Cd8b1"),
                "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                "Pericytes"=c("Pdgfrb", "Rgs5"))

# make dot plots for each resolution
res_cols <- grep("RNA_snn_res", colnames(peri_seu@meta.data), value=T)

for (res in res_cols) {
  
  Idents(peri_seu) <- res
  
  dotplot_list <- lapply(names(markers), function(cluster_name) {
    
    markers <- markers[[cluster_name]]
    
    DotPlot(peri_seu, features = markers) + RotatedAxis() + 
      labs(x=NULL, y=NULL, title=cluster_name, color="Avg Exp", size="% Exp") 
    
  })
  
  plot_grid(plotlist = dotplot_list, nrow = 3)
  
  ggsave(paste0(out_dir, res, ".dot_plot.png"), width=16, height=12, bg="white")
  
}


DimPlot(peri_seu, reduction="umap.peri_pca", 
        group.by= res_cols,
        ncol=3, label = T)
ggsave(paste0(out_dir, "res_clustering.umaps.png"), width=13, height=8)

# go with the lowest resolution

DimPlot(peri_seu, reduction="umap.peri_pca", 
        group.by= res_cols[1], label = T)

# cluster marker investigation

Idents(peri_seu) <- "RNA_snn_res.0.3"

cluster2_markers <- FindMarkers(peri_seu,
                                ident.1 = "2",
                                #ident.2 = clusters[!clusters %in% c("5")],
                                test.use="wilcox")
cluster2_markers$gene <- rownames(cluster2_markers)
cluster2_markers$pct_delta <- cluster2_markers$pct.1 -
  cluster2_markers$pct.2

cluster2_markers <- cluster2_markers[order(cluster2_markers$pct_delta, 
                                           decreasing = T),]

DotPlot(peri_seu, features=cluster2_markers$gene[1:10]) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "cluster2_markers.dotplot.png"), width=7, height=4, bg="white")

# NK CELLS

cluster3_markers <- FindMarkers(peri_seu,
                                ident.1 = "3",
                                #ident.3 = clusters[!clusters %in% c("5")],
                                test.use="wilcox")
cluster3_markers$gene <- rownames(cluster3_markers)
cluster3_markers$pct_delta <- cluster3_markers$pct.1 -
  cluster3_markers$pct.2

cluster3_markers <- cluster3_markers[order(cluster3_markers$pct_delta, 
                                           decreasing = T),]

DotPlot(peri_seu, features=cluster3_markers$gene[1:10]) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "cluster3_markers.dotplot.png"), width=7, height=4, bg="white")

cluster4_markers <- FindMarkers(peri_seu,
                                ident.1 = "4",
                                #ident.4 = clusters[!clusters %in% c("5")],
                                test.use="wilcox")
cluster4_markers$gene <- rownames(cluster4_markers)
cluster4_markers$pct_delta <- cluster4_markers$pct.1 -
  cluster4_markers$pct.2

cluster4_markers <- cluster4_markers[order(cluster4_markers$pct_delta, 
                                           decreasing = T),]

DotPlot(peri_seu, features=cluster4_markers$gene[1:10]) + RotatedAxis() +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "cluster4_markers.dotplot.png"), width=7, height=4, bg="white")

# PLATELETS

cluster_xref <- data.frame("RNA_snn_res.0.3"=sort(unique(peri_seu$RNA_snn_res.0.3)),
                           pericyte_type=c("CD8+",
                                           "Pericytes",
                                           "NK",
                                           "Unknown",
                                           "Platelets"))

metadata <- peri_seu@meta.data

metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cluster_xref,
                  by="RNA_snn_res.0.3")
rownames(metadata) <- metadata$cell_id

metadata <- metadata[rownames(peri_seu@meta.data),]

peri_seu$pericyte_type <- metadata$pericyte_type

DimPlot(peri_seu, reduction="umap.peri_pca", 
        group.by= c("RNA_snn_res.0.3","pericyte_type"),
        label = T)
ggsave(paste0(out_dir, "cluster2cell_type.umap.png"), width=11, height=5)

# save Seurat object
SaveSeuratRds(peri_seu, paste0(out_dir, "pericytes_subset.named.seurat.RDS"))

saveRDS(peri_seu@meta.data, file=paste0(out_dir, "pericytes_subset.named.metadata.RDS"))


metadata <- peri_seu@meta.data

ggplot(metadata,
       aes(x=RNA_snn_res.0.3,
           y=nFeature_RNA)) +
  geom_violin() +
  theme_bw()




