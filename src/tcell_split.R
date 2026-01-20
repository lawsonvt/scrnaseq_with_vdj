library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)

out_dir <- "results/tcell_split/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming/cell_named.seurat.RDS")

# subset to the tcell ones
tcell_seu <- subset(seu_obj, subset = cell_cluster == "T Cells")

table(tcell_seu@meta.data$orig.ident)

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

# create umap
tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("orig.ident","tcell_clusters"))

grep("Cd[384][a-z]+", rownames(tcell_seu), value=T)

# tcell markers
t_markers <- c("Cd3g",
               "Cd3d",
               "Cd3e",
                "Cd4",
               "Cd8a",
               "Cd8b1")
               # "Ptprc",
               # "Ccr7",
               # "Cd27",
               # "Cd28",
               # "Cd69",
               # "Il2ra",
               # "H2-Eb1",
               # "H2-Eb2")



p1 <- DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= "tcell_clusters",label = T, label.box = T) +
  labs(x="UMAP 1", y="UMAP 2", title='T Cells') + theme(legend.position = "none")

p2 <- DotPlot(tcell_seu, features = t_markers) + RotatedAxis()

plot_grid(p1, p2)
ggsave(paste0(out_dir, "tcell_clusters.umap_and_dotplot.png"), width=10, height=6, bg="white")


markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                "T Cells"=c("Trbc2", "Cd3d", "Lck"))

dotplot_list <- lapply(names(markers), function(cluster_name) {
  
  markers <- markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "b_vs_t_cells.dotplot.png"), width=7, height=5, bg="white")


# find markers for cluster 7, "B Cells"

# the one we are interested in
cluster7_markers <- FindMarkers(tcell_seu,
                        ident.1 = "7",
                        test.use="MAST")
cluster7_markers$gene <- rownames(cluster7_markers)

# find significant Ig genes
sig_ig <- cluster7_markers[cluster7_markers$p_val_adj < 0.01 &
                             grepl("^Ig", cluster7_markers$gene) &
                             cluster7_markers$avg_log2FC > 0 &
                             cluster7_markers$pct.1 > 0.2,]

VlnPlot(tcell_seu, features = sig_ig$gene, ncol=2, raster=F, alpha=0.2)
ggsave(paste0(out_dir, "cluster7_igmarkers.png"), width=10, height=8)


# see what other cell types clusters 2 and 11 can be
lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"))

# cow dot plot for the markers
dotplot_list <- lapply(names(lab_cluster_markers), function(cluster_name) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})
plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "tcell_clusters.full_markers.dotplot.png"), width=16, height=7, bg="white")



DotPlot(tcell_seu, features = c("Trbc2", "Cd3d", "Lck")) + RotatedAxis()

# name the clusters

cluster_xref <- data.frame(tcell_clusters=sort(unique(tcell_seu@meta.data$tcell_clusters)),
                           tcell_type=c("CD8+",
                                        "CD8+",
                                        "T Cell",
                                        "CD4+",
                                        "CD8+",
                                        "CD4+",
                                        "CD8+",
                                        "B Cell",
                                        "CD8+",
                                        "CD4+",
                                        "CD8+",
                                        "T Cell",
                                        "Monocytes"))

metadata <- tcell_seu@meta.data

metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cluster_xref,
                  by="tcell_clusters")
rownames(metadata) <- metadata$cell_id

metadata <- metadata[rownames(tcell_seu@meta.data),]

sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")


metadata$sample <- factor(as.character(metadata$orig.ident),
                          levels=sample_levels)

metadata$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", metadata$orig.ident)

metadata$condition <- factor(as.character(metadata$condition),
                             levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

# add into the seurat object

tcell_seu@meta.data$tcell_type <- metadata$tcell_type
tcell_seu@meta.data$sample <- metadata$sample
tcell_seu@meta.data$condition <- metadata$condition


DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= "tcell_type",) +
  labs(x="UMAP 1", y="UMAP 2", title='T Cells') 
ggsave(paste0(out_dir, "tcell_types.umap.png"), width=5, height=4)

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= "tcell_type",
        split.by = "condition", ncol=2) +
  labs(x="UMAP 1", y="UMAP 2", title='T Cells') 
ggsave(paste0(out_dir, "tcell_types.per_condition.umap.png"), width=7, height=5)

# save Seurat object
SaveSeuratRds(tcell_seu, paste0(out_dir, "tcell_subset.named.seurat.RDS"))



