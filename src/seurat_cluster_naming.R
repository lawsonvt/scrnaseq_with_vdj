library(Seurat)
library(ggplot2)
library(scRepertoire)
library(snakecase)
library(cowplot)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(openxlsx)

out_dir <- "results/seurat_cluster_naming/"
dir.create(out_dir, showWarnings = F)

# cluster markers as supplied by Lukens Lab

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"))

# read in integrated seurat
int_seu <- LoadSeuratRds("results/join_and_filter_integrated_seurat/filtered.integrated_seurat.RDS")

analysis_genes <- rownames(int_seu)

# find any issues between marker genes and analysis genes

lab_cluster_markers_missing <- lapply(lab_cluster_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

# cow dot plot for the markers
dotplot_list <- lapply(names(lab_cluster_markers), function(cluster_name) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "total_dotplots.png"), width=18, height=7, bg="white")

# The pericytes results are misleading, as they are from just a fraction of the cells

DotPlot(int_seu, features=lab_cluster_markers$Pericytes) + 
  RotatedAxis() + labs(x=NULL, y=NULL, title="Pericytes")
ggsave(paste0(out_dir, "pericytes_markers.png"), width=4, height=7, bg="white")

# pull out cluster 22 markers
cluster22_markers <- readRDS("results/cluster_marker_tests/cluster22_degs.wilcox.RDS")
cluster22_markers$gene <- rownames(cluster22_markers)

cluster22_markers$delta_pct <- cluster22_markers$pct.1 - cluster22_markers$pct.2

cluster22_markers <- cluster22_markers[order(cluster22_markers$delta_pct,
                                             decreasing = T),]

c22 <- cluster22_markers[cluster22_markers$delta_pct > 0.45,]$gene

DotPlot(int_seu, 
        features = c22)
ggsave(paste0(out_dir, "cluster22_markers.dot_plot.png"), width=10, height=7, bg="white")


# make a cross reference after looking over the dot plot

harmony_clusters <- sort(unique(int_seu@meta.data$harmony_clusters))

cluster_xref <- data.frame(harmony_clusters=harmony_clusters,
                           cell_cluster=c("Microglia", 
                                          "Microglia",
                                          "Microglia",
                                          "Microglia",
                                          "Microglia",
                                          "Microglia",
                                          "Microglia",
                                          "Microglia",
                                          "Macrophages",
                                          "Microglia",
                                          "Microglia",
                                          "Macrophages",
                                          "Microglia",
                                          "Microglia",
                                          "T Cells",
                                          "Microglia",
                                          "Neutrophils",
                                          "B Cells",
                                          "Monocytes",
                                          "Proliferating Cells",
                                          "Monocytes",
                                          "Pericytes",
                                          "Macrophages",
                                          "Monocytes",
                                          "Proliferating Cells",
                                          "B Cells",
                                          "Proliferating Cells"))

# merge them into the metadata
metadata <- int_seu@meta.data
# merge in new cluster names

metadata$cell_id <- rownames(metadata)

metadata <- merge(metadata,
                  cluster_xref,
                  by="harmony_clusters")

# fix the order
rownames(metadata) <- metadata$cell_id
metadata <- metadata[rownames(int_seu@meta.data),]

# add back into object
int_seu <- AddMetaData(int_seu, metadata$cell_cluster, col.name="cell_cluster")

# plot this out!
DimPlot(int_seu, reduction="umap.harmony", group.by= c("harmony_clusters", "cell_cluster"),
        label=T)
ggsave(paste0(out_dir, "cluster_cells_umap.with_legend.png"), width=14, height=6)


p1 <- DimPlot(int_seu, reduction="umap.harmony", group.by= "harmony_clusters",
              label=T) + labs(x="UMAP 1", y="UMAP 2", title="Harmony Clusters") + 
  theme(legend.position = "none")

p2 <- DimPlot(int_seu, reduction="umap.harmony", group.by= "cell_cluster",
              label=T) + labs(x="UMAP 1", y="UMAP 2", title="Cell Types") + 
  theme(legend.position = "none")

plot_grid(p1, p2, nrow=1)
ggsave(paste0(out_dir, "cluster_cells_umap.without_legend.png"), width=12, height=6)


DimPlot(int_seu, reduction="umap.harmony", group.by="cell_cluster",
        label=T, split.by = "orig.ident", ncol=4, raster=F) + 
  theme(legend.position = "none") + labs(x="UMAP 1", y= "UMAP 2", title="Cell Types")
ggsave(paste0(out_dir, "cluster_cells_umap.per_sample.png"), width=17, height=10)

# barplots!

metadata$sample <- metadata$orig.ident

# factorize for plotting
sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")

metadata$sample <- factor(as.character(metadata$sample),
                           levels=sample_levels)

ggplot(metadata,
       aes(x=cell_cluster,
           fill=sample)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave(paste0(out_dir, "total_celltypes.cell_count.barplot.png"),
       width=6, height=4)

ggplot(metadata,
       aes(x=sample,
           fill=sample)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave(paste0(out_dir, "total.cell_count.barplot.png"),
       width=6, height=4)

# remove the microglia

ggplot(metadata[metadata$cell_cluster != "Microglia",],
       aes(x=cell_cluster,
           fill=sample)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave(paste0(out_dir, "non_microglia.cell_count.barplot.png"),
       width=6, height=4)


# dot plot of DCs
dc_genes <- c("Xcr1", "Sirpa", "Il3ra", "Ccr7")

DotPlot(int_seu, 
        features=dc_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "dotplot_for_dc_genes.png"), height=6, width=7, bg="white")

Idents(int_seu) <- "cell_cluster"
DotPlot(int_seu, 
        features=dc_genes) + labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "dotplot_for_dc_genes.cell_names.png"), height=4, width=7, bg="white")


# save things
write.xlsx(cluster_xref, paste0(out_dir, "cluster_xref.xlsx"), colWidths="auto")

SaveSeuratRds(int_seu, file=paste0(out_dir, "cell_named.seurat.RDS"))

