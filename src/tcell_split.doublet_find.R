library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(scDblFinder)
library(BiocParallel)
library(clustree)

out_dir <- "results/tcell_split.doublet_find/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming/cell_named.seurat.RDS")

# subset to the tcell ones
tcell_seu <- subset(seu_obj, subset = cell_cluster == "T Cells")

# determine the optimal number of dimensions through an elbow plot

tcell_seu <- RunPCA(tcell_seu, npcs = 20)

# inspect elbow plot
ElbowPlot(tcell_seu, ndims=20) + 
  labs(title="TCell Subset")
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 5

# cluster the harmonized data
tcell_seu <- FindNeighbors(tcell_seu, dims = 1:max_pc_dim, reduction = "pca")
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_clusters")

# run doublet finder
tcell_sce <- as.SingleCellExperiment(tcell_seu)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
tcell_sce <- scDblFinder(tcell_sce, clusters="tcell_clusters", BPPARAM=bp)

# add doublet calls to seurat object
tcell_seu$scDblFinder.class <- tcell_sce$scDblFinder.class

VlnPlot(tcell_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "tcell.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(tcell_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "tcell.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(tcell_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "tcell.doublet_output.RDS"))

# remove doublets, redo clustering
tcell_seu <- subset(tcell_seu, scDblFinder.class == "singlet")

tcell_seu <- RunPCA(tcell_seu, npcs = 20)

# inspect elbow plot
ElbowPlot(tcell_seu, ndims=20) + 
  labs(title="TCell Subset")
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.post_dbl_removal.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 5

# cluster the harmonized data
tcell_seu <- FindNeighbors(tcell_seu, dims = 1:max_pc_dim, reduction = "pca")
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_clusters")

# create umap
tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("orig.ident","tcell_clusters"))

# try different resolutions for clustering
tcell_seu <- FindClusters(tcell_seu, resolution = c(0.3, 0.5, 0.8, 1.0, 1.2))

clustree(tcell_seu, prefix = "RNA_snn_res.")
ggsave(paste0(out_dir, "res_clustering_tree.png"), width=9, height=8)

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
res_cols <- grep("RNA_snn_res", colnames(tcell_seu@meta.data), value=T)

for (res in res_cols) {
  
  Idents(tcell_seu) <- res

  dotplot_list <- lapply(names(markers), function(cluster_name) {
    
    markers <- markers[[cluster_name]]
    
    DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
      theme(legend.position="none")
    
  })
  
  plot_grid(plotlist = dotplot_list, nrow = 1)
  
  ggsave(paste0(out_dir, res, ".dot_plot.png"), width=18, height=7, bg="white")

}


DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= res_cols,
        ncol=3)

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= res_cols,
        ncol=3, label = T)
ggsave(paste0(out_dir, "res_clustering.umaps.png"), width=13, height=8)

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= "RNA_snn_res.0.8",
        label = T)

# proposed cell clusters

cluster_xref <- data.frame("RNA_snn_res.0.8"=sort(unique(tcell_seu$RNA_snn_res.0.8)),
                           tcell_type=c("CD8+",
                                        "CD4+",
                                        "Unknown2",
                                        "CD8+",
                                        "Proliferating Cells",
                                        "CD4+",
                                        "CD4+",
                                        "Macrophages",
                                        "CD8+",
                                        "Unknown9",
                                        "B Cells",
                                        "CD8+",
                                        "Monocytes"))

metadata <- tcell_seu@meta.data

metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cluster_xref,
                  by="RNA_snn_res.0.8")
rownames(metadata) <- metadata$cell_id

metadata <- metadata[rownames(tcell_seu@meta.data),]

tcell_seu$tcell_type <- metadata$tcell_type

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= "tcell_type",
        label = T)

# find markers for clusters 2 and 9

Idents(tcell_seu) <- "RNA_snn_res.0.8"

clusters <- unique(Idents(tcell_seu))

cluster2_markers <- FindMarkers(tcell_seu,
                                ident.1 = "2",
                                ident.2 = clusters[!clusters %in% c("2","9")],
                                test.use="MAST")
cluster2_markers$gene <- rownames(cluster2_markers)
cluster2_markers$pct_delta <- cluster2_markers$pct.1 -
  cluster2_markers$pct.2

cluster9_markers <- FindMarkers(tcell_seu,
                                 ident.1 = "9" ,
                                ident.2 = clusters[!clusters %in% c("2","9")],
                                 test.use="MAST")
cluster9_markers$gene <- rownames(cluster9_markers)
cluster9_markers$pct_delta <- cluster9_markers$pct.1 -
  cluster9_markers$pct.2

cluster_2_9_merged <- merge(cluster2_markers,
                             cluster9_markers,
                             by="gene",
                             suffixes=c("_2", "_9"))

markers_2_9 <- cluster_2_9_merged[cluster_2_9_merged$pct_delta_2 > 0.4 &
                                      cluster_2_9_merged$pct_delta_9 > 0.4,]$gene


cluster_2_9_merged$delta_delta <- cluster_2_9_merged$pct_delta_2 -
  cluster_2_9_merged$pct_delta_9

markers_2 <- cluster_2_9_merged[cluster_2_9_merged$delta_delta> 0.4 &
                                  cluster_2_9_merged$pct_delta_9 > -0.1,]$gene
# drop ensembl names
markers_2 <- markers_2[!grepl("ENSMUSG", markers_2)]

markers_9 <- cluster_2_9_merged[cluster_2_9_merged$delta_delta < -0.4 &
                                  cluster_2_9_merged$pct_delta_2 > -0.2,]$gene

p1 <- DotPlot(tcell_seu, features = markers_2_9) +
  RotatedAxis() + labs(x=NULL, y=NULL, title="Markers 2 & 9") +
  theme(legend.position="none") 

p2 <- DotPlot(tcell_seu, features = markers_2) +
  RotatedAxis() + labs(x=NULL, y=NULL, title="Markers 2") +
  theme(legend.position="none")

p3 <- DotPlot(tcell_seu, features = markers_9) +
  RotatedAxis() + labs(x=NULL, y=NULL, title="Markers 9") +
  theme(legend.position="none")

plot_grid(p1,p2,p3, nrow=1)
ggsave(paste0(out_dir, "marker_2_9.dot_plots.png"), width=12, height=6, bg="white")

# new cluster names
cluster_xref <- data.frame("RNA_snn_res.0.8"=sort(unique(tcell_seu$RNA_snn_res.0.8)),
                           tcell_type=c("CD8+",
                                        "CD4+",
                                        "NK",
                                        "CD8+",
                                        "Proliferating Cells",
                                        "CD4+",
                                        "CD4+",
                                        "Macrophages",
                                        "CD8+",
                                        "NK",
                                        "B Cells",
                                        "CD8+",
                                        "Monocytes"))

metadata <- tcell_seu@meta.data

# clear out what is there
metadata$tcell_type <- NULL

metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cluster_xref,
                  by="RNA_snn_res.0.8")
rownames(metadata) <- metadata$cell_id

metadata <- metadata[rownames(tcell_seu@meta.data),]

tcell_seu$tcell_type <- metadata$tcell_type

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= c("RNA_snn_res.0.8","tcell_type"),
        label = T)
ggsave(paste0(out_dir, "cluster2cell_type.umap.png"), width=11, height=5)

Idents(tcell_seu) <- "tcell_type"

dotplot_list <- lapply(names(markers), function(cluster_name) {
  
  markers <- markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)

ggsave(paste0(out_dir, "cluster_names.dot_plot.png"), width=11, height=8, bg="white")

# save Seurat object
SaveSeuratRds(tcell_seu, paste0(out_dir, "tcell_subset.named.seurat.RDS"))



