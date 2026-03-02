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
library(ggsci)

out_dir <- "results/tcell_split.tcell_subset2/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/tcell_split.post_pericytes_split/tcell_cell_named.seurat.RDS")

# filter down to those that are T Cells
tcell_init_cells <- c("Naive T Cells",
                      "CD8+",
                      "CD4+",
                      "Proliferating CD8+")

tcell_seu <- subset(seu_obj, subset = tcell_sub_cell_name %in% tcell_init_cells)

# re subset

# determine the optimal number of dimensions through an elbow plot

tcell_seu <- RunPCA(tcell_seu, npcs = 20)

# inspect elbow plot
ElbowPlot(tcell_seu, ndims=20) + 
  labs(title="TCell Subset")
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 10

# cluster the harmonized data
tcell_seu <- FindNeighbors(tcell_seu, dims = 1:max_pc_dim, reduction = "pca")
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_subset_clusters")

tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_subset_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("orig.ident","tcell_subset_clusters"))
DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("tcell_sub_cell_name","tcell_subset_clusters"), ncol = 1,
        label=T)
ggsave(paste0(out_dir, "init_tcell_v_new_clusters.umap.png"), width=7, height=9)

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("tcell_subset_clusters"), 
        label=T, label.box = T) +
  labs(x="UMAP1", y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "new_clusters.umap.png"), width=7, height=6)

# count up the cells
tcell_meta <- tcell_seu@meta.data

ggplot(tcell_meta,
       aes(x=tcell_subset_clusters,
           fill=orig.ident)) +
  geom_bar() +
  theme_bw()

ggplot(tcell_meta,
       aes(x=tcell_subset_clusters,
           fill=condition)) +
  geom_bar() +
  theme_bw()

# make heat colors for clonal counts
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

ggplot(tcell_meta,
       aes(x=tcell_subset_clusters,
           fill=cloneSize)) +
  facet_wrap(~ condition, ncol=2) +
  geom_bar(color="black") +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x=NULL, y="Cell Count", fill="Clone Size")
ggsave(paste0(out_dir, "clone_count_per_cluster.png"), width=10, height=8)


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

Idents(tcell_seu) <- "tcell_subset_clusters"

dotplot_list <- lapply(names(markers), function(cluster_name) {
  
  markers <- markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)

ggsave(paste0(out_dir, "total_markers_clusters.dot_plot.png"), width=16, height=12, bg="white")


# determine cluster markers

cluster_markers <- lapply(sort(unique(tcell_meta$tcell_subset_clusters)), function(cluster) {
  
  markers <- FindMarkers(tcell_seu,
                         ident.1 = cluster,
                         test.use="wilcox")
  markers <- rownames_to_column(markers, var="gene")
  
  markers$delta_pct <- markers$pct.1 - markers$pct.2
  markers$cluster <- cluster
  
  markers <- markers[order(markers$delta_pct, decreasing = T),]
  
  return(markers)
})
names(cluster_markers) <- sort(unique(tcell_meta$tcell_subset_clusters))

cluster_markers_df <- bind_rows(cluster_markers)

# using a delta pct cutoff doesnt seem to work, top 5 for each?

top_cluster_markers <- lapply(cluster_markers, head)
top_cluster_markers_df <- bind_rows(top_cluster_markers)

# drop the dupes
top_cluster_markers_df <- top_cluster_markers_df[!duplicated(top_cluster_markers_df$gene),]

# get the expression matrix

exp_matrix <- tcell_seu[["RNA"]]$scale.data

sum(top_cluster_markers_df$gene %in% rownames(exp_matrix))

# Tcell markers + markers from paper

tcell_markers <- list("T Cells"=c("Trbc2", "Cd3d", "Lck"),
                      "T CD"=c("Cd4",
                               "Cd8a",
                               "Cd8b1"),
                      "CD4-Foxp3"=c("Cd4",
                                    "Foxp3"),
                      "CD4-Folr4r-Slamf6"=c("Cd4",
                                            "Izumo1r",
                                            "Slamf6"),
                      "CD4-CD40lg"=c("Cd4",
                                     "Cd40lg"),
                      "CD8-Tox-Pdcd1"=c("Cd8a",
                                        "Tox",
                                        "Pdcd1"),
                      "CD8-Xcl"=c("Cd8a",
                                  "Xcl1"),
                      "CD8-Ly6c2"=c("Cd8a",
                                    "Ly6c2"),
                      "CD8-CD7-Gzmb"=c("Cd8a",
                                       "Cd7",
                                       "Gzmb"),
                      "CD8-Isg15"=c("Cd8a",
                                    "Isg15"),
                      "CD8-Egr1"=c("Cd8a",
                                   "Egr1"),
                    "CD8-CD11c-Klre1"=c("Cd8a",
                                          "Itgax",
                                          "Klre1"),
                      "CD8-CCR7-Sell"=c("Cd8a",
                                        "Ccr7",
                                        "Sell"),
                      "CD8-Bcl3-Relb"=c("Cd8a",
                                        "Bcl3",
                                        "Relb"),
                      "CD8-Cx3cr1-Gzma"=c("Cd8a",
                                          "Cx3cr1",
                                          "Gzma"))


analysis_genes <- rownames(tcell_seu)

# find any issues between marker genes and analysis genes

tcell_markers_missing <- lapply(tcell_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

Idents(tcell_seu) <- "tcell_subset_clusters"

dotplot_list <- lapply(names(tcell_markers), function(cluster_name) {
  
  markers <- tcell_markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)

ggsave(paste0(out_dir, "tcell_markers_clusters.dot_plot.png"), width=16, height=12, bg="white")


# make additional plots for clusters 4/6

calc_cluster_markers <- list("T Cells"=c("Trbc2", "Cd3d", "Lck"),
                             "Cluster 4"= head(cluster_markers[["4"]]$gene, n=10),
                             "T CD"=c("Cd4",
                                      "Cd8a",
                                      "Cd8b1"),
                             "Cluster 6"= head(cluster_markers[["6"]]$gene, n=10))


Idents(tcell_seu) <- "tcell_subset_clusters"

dotplot_list <- lapply(names(calc_cluster_markers), function(cluster_name) {
  
  markers <- calc_cluster_markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 2)
ggsave(paste0(out_dir, "cluster_4_6_markers.dotplot.png"), width=10, height=7, bg="white")


# one big dotplot?
top_markers <- unlist(lapply(cluster_markers, function(x) {head(x$gene,n=10)}))

DotPlot(tcell_seu, features=unique(top_markers), 
        split.by = "tcell_subset_clusters",
        cols = pal_d3("category10")(10)) + RotatedAxis() +
  theme(panel.grid = element_line(color = "grey", linewidth = 0.5))

DotPlot(tcell_seu, features=unique(top_markers)) + RotatedAxis() +
  theme(panel.grid = element_line(color = "grey", linewidth = 0.5),
        axis.text.x = element_text(size=9)) +
  labs(x=NULL, y="T Cell Cluster")
ggsave(paste0(out_dir, "top_markers_per_cluster.png"), width=17, height=8, bg='white')

# per cluster

pdf(paste0(out_dir, "top_markers_per_cluster.pdf"), width=8, height=5)

for (cluster_name in names(cluster_markers)) {
  
  top_markers <- head(cluster_markers[[cluster_name]]$gene, n=10)

  print(DotPlot(tcell_seu, features=top_markers) + RotatedAxis() + 
    labs(x=NULL, y="T Cell Cluster", title=paste0("Cluster ", cluster_name)))
  
  
    
}
dev.off()


# assign new categories
tcell_subset_cluster_xref <- data.frame(tcell_subset_clusters=sort(unique(tcell_seu$tcell_subset_clusters)),
                                        tcell_subset_clusters_name=c("CD8_0",
                                                                     "CD8_1",
                                                                     "CD4_2",
                                                                     "CD4_3",
                                                                     "TCELL_4",
                                                                     "CD8_5",
                                                                     "TCELL_6",
                                                                     "CD8_7",
                                                                     "CD8_8",
                                                                     "CD8_9"))

# merge into Seurat object

tcell_meta <- tcell_seu@meta.data
tcell_meta$cell_id <- rownames(tcell_meta)

tcell_meta <- merge(tcell_meta,
                    tcell_subset_cluster_xref,
                    by = "tcell_subset_clusters")


rownames(tcell_meta) <- tcell_meta$cell_id

tcell_meta <- tcell_meta[colnames(tcell_seu),]

tcell_seu$tcell_subset_clusters_name <- tcell_meta$tcell_subset_clusters_name

# save the data
SaveSeuratRds(tcell_seu, file=paste0(out_dir, "tcell_subset2_named.seurat.RDS"))
saveRDS(tcell_meta, file=paste0(out_dir, "tcell_subset2_named.metadata.RDS"))
