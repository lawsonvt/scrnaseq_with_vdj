library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(scDblFinder)
library(BiocParallel)
library(clustree)
library(scRepertoire)
library(dplyr)
library(gtools)
library(tibble)
library(DESeq2)
library(ggrepel)
library(snakecase)
library(openxlsx)

out_dir <- "results/seurat_cluster_naming.microglia_subset.unknown_subset/"
dir.create(out_dir, showWarnings = F)

# load in microglia subset
seu_mg <- LoadSeuratRds("results/seurat_cluster_naming.microglia_subset/cell_named.microglia_subset.seurat.RDS")

# subset down to the unknown cells
seu_unknown <- subset(seu_mg, subset = microglia_cell_name == "MG Unknown")

# recluster them

# determine the optimal number of dimensions through an elbow plot

seu_unknown <- RunPCA(seu_unknown, npcs = 20)

# inspect elbow plot
ElbowPlot(seu_unknown, ndims=20) + 
  labs(title="Unknown Microglia")
ggsave(paste0(out_dir, "mg_unknown.pca_elbow_plot.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 10

# cluster the harmonized data
seu_unknown <- FindNeighbors(seu_unknown, dims = 1:max_pc_dim, reduction = "pca")
seu_unknown <- FindClusters(seu_unknown, cluster.name = "unknown_clusters")

# run doublet finder
sce_unknown <- as.SingleCellExperiment(seu_unknown)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
sce_unknown <- scDblFinder(sce_unknown, clusters="unknown_clusters", BPPARAM=bp)

# add doublet calls to seurat object
seu_unknown$scDblFinder.class <- sce_unknown$scDblFinder.class

VlnPlot(seu_unknown, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "unknown.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(seu_unknown$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "unknown.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(sce_unknown))

saveRDS(dbl_meta, file=paste0(out_dir, "unknown.doublet_output.RDS"))

# remove doublets, redo clustering
seu_unknown <- subset(seu_unknown, scDblFinder.class == "singlet")

seu_unknown <- RunPCA(seu_unknown, npcs = 20)

# inspect elbow plot
ElbowPlot(seu_unknown, ndims=20) + 
  labs(title="Unknown Microglia")
ggsave(paste0(out_dir, "unknown_subset.pca_elbow_plot.post_dbl_removal.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 10

# cluster the harmonized data
seu_unknown <- FindNeighbors(seu_unknown, dims = 1:max_pc_dim, reduction = "pca")
seu_unknown <- FindClusters(seu_unknown, cluster.name = "unknown_clusters")

# create umap
seu_unknown <- RunUMAP(seu_unknown, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.unknown_pca")

DimPlot(seu_unknown, reduction="umap.unknown_pca", group.by= c("orig.ident","unknown_clusters"))

DimPlot(seu_unknown, reduction="umap.unknown_pca", group.by= c("microglia_clusters","unknown_clusters"))


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

Idents(seu_unknown) <- "unknown_clusters"

# cow dot plot for the markers
dotplot_list <- lapply(names(markers), function(cluster_name) {
  
  markers <- markers[[cluster_name]]
  
  DotPlot(seu_unknown, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) 
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)


# try different resolutions for clustering
seu_unknown <- FindClusters(seu_unknown, resolution = c(0.3, 0.5, 0.8, 1.0, 1.2))

res_cols <- grep("RNA_snn_res", colnames(seu_unknown@meta.data), value=T)

DimPlot(seu_unknown, reduction="umap.unknown_pca", 
        group.by= res_cols,
        ncol=3, label=T)


# make dot plots for each resolution

for (res in res_cols) {
  
  Idents(seu_unknown) <- res
  
  dotplot_list <- lapply(names(markers), function(cluster_name) {
    
    markers <- markers[[cluster_name]]
    
    DotPlot(seu_unknown, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
    
  })
  
  plot_grid(plotlist = dotplot_list, nrow = 3)
  
  ggsave(paste0(out_dir, res, ".dot_plot.png"), width=16, height=12, bg="white")
  
  
}

# find markers for the lowest resolution value
Idents(seu_unknown) <- "RNA_snn_res.0.3"

clusters <- sort(unique(seu_unknown@meta.data$RNA_snn_res.0.3))

cluster_markers <- lapply(clusters, function(cluster) {
  
  markers <- FindMarkers(seu_unknown,
                         ident.1 = cluster,
                         test.use="wilcox")
  markers <- rownames_to_column(markers, var="gene")
  
  markers$delta_pct <- markers$pct.1 - markers$pct.2
  markers$cluster <- cluster
  
  markers <- markers[order(markers$delta_pct, decreasing = T),]
  
  return(markers)
})

names(cluster_markers) <- clusters

cluster_markers_df <- bind_rows(cluster_markers)

top_markers <- cluster_markers_df[cluster_markers_df$delta_pct > 0.5,]

DotPlot(seu_unknown, features = unique(top_markers$gene)) + RotatedAxis() +
  labs(x=NULL, y="Cluster")
ggsave(paste0(out_dir, "cluster_markers.unknown_subset.png"), width=8, height=6, bg="white")


DimPlot(seu_unknown, reduction="umap.unknown_pca", group.by = "RNA_snn_res.0.3") +
  labs(x="UMAP 1", y="UMAP 2")




