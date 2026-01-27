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
library(ggsci)

out_dir <- "results/seurat_cluster_naming.microglia_subset/"
dir.create(out_dir, showWarnings = F)

micro_markers <- list(Homeostatic=c("P2ry12",
                                    "Tmem119",
                                    "Cx3cr1"),
                      DAM=c("Lpl",
                            "Cst7",
                            "Axl",
                            "Apoe",
                            "Clec7a"),
                      "Interferon Responsive"=c("Ifitm3",
                                                "Mx1",
                                                "Irf8",
                                                "Mef2c"),
                      "Antigen Presenting"=c("B2m",
                                             "Ciita",
                                             "Cd80",
                                             "Cd86",
                                             "Tap1",
                                             "H2-K1",
                                             "H2-D1",
                                             "H2-Ab1"),
                      "Inflammatory"=c("Cxcl16",
                                       "Cxcl10",
                                       "Ccl2",
                                       "Il1b",
                                       "Il6",
                                       "Tnf"),
                      "Cross presenting"=c("Wdfy3",
                                           "Wdfy4",
                                           "Sec22b",
                                           "Cd68"),
                      "BAM"=c("Lyve1", 
                              "P2rx7", 
                              "Mrc1", 
                              "Siglec1"))

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"))

# read in integrated seurat
seu_subset <- LoadSeuratRds("results/seurat_cluster_microglia_subset/microglia_subset_seurat.RDS")

analysis_genes <- rownames(seu_subset)

# find any issues between marker genes and analysis genes

micro_markers_missing <- lapply(micro_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

cell_meta <- seu_subset@meta.data

clusters <- sort(unique(cell_meta[cell_meta$cell_cluster == "Microglia",]$microglia_clusters))

# cow dot plot for the markers
dotplot_list <- lapply(names(lab_cluster_markers), function(cluster_name) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  DotPlot(seu_subset, features = markers, idents=clusters) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "lab_cluster_marker_dotplots.png"), width=19, height=7, bg="white")

# dot plot for the other markers
dotplot_list <- lapply(names(micro_markers), function(cluster_name) {
  
  markers <- micro_markers[[cluster_name]]
  
  DotPlot(seu_subset, features = markers, idents=clusters) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)

DimPlot(seu_subset, reduction="umap.microglia_pca", group.by= "microglia_clusters", label=T,
        raster = T) +
  labs(x="UMAP 1", y="UMAP 2")


# find some markers

# possible grouping 5 and 14
clusters <- unique(Idents(seu_subset))

table(Idents(seu_subset))

# subsample

set.seed(42)

sampled_cells <- sample(colnames(seu_subset), size=12000, replace = F)

seu_subset_sub <- subset(seu_subset, cells=sampled_cells)

cluster5_markers <- FindMarkers(seu_subset_sub,
                                ident.1 = "5",
                                ident.2 = clusters[!clusters %in% c("5","14")],
                                test.use="MAST")
cluster5_markers$gene <- rownames(cluster5_markers)
cluster5_markers$pct_delta <- cluster5_markers$pct.1 -
  cluster5_markers$pct.2

cluster14_markers <- FindMarkers(seu_subset_sub,
                                ident.1 = "14" ,
                                ident.2 = clusters[!clusters %in% c("5","14")],
                                test.use="MAST")
cluster14_markers$gene <- rownames(cluster14_markers)
cluster14_markers$pct_delta <- cluster14_markers$pct.1 -
  cluster14_markers$pct.2

cluster_5_14_merged <- merge(cluster5_markers,
                             cluster14_markers,
                             by="gene",
                             suffixes=c("_5","_14"))

markers_5_14 <- cluster_5_14_merged[cluster_5_14_merged$pct_delta_5 > 0.3 &
                                      cluster_5_14_merged$pct_delta_14 > 0.3,]


DotPlot(seu_subset, features = markers_5_14$gene, idents=clusters) + RotatedAxis() + 
  labs(x=NULL, y=NULL, title="Clusters 5,4,12 markers") +
  theme(legend.position="none")


# cluster_group <- c("1","16","13","2")
# 
# cluster1_markers <- FindMarkers(seu_subset_sub,
#                                 ident.1 = "1",
#                                 ident.2 = clusters[!clusters %in% cluster_group],
#                                 test.use="MAST")
# cluster1_markers$gene <- rownames(cluster1_markers)
# cluster1_markers$pct_delta <- cluster1_markers$pct.1 -
#   cluster1_markers$pct.2
# 
# cluster16_markers <- FindMarkers(seu_subset_sub,
#                                 ident.1 = "16",
#                                 ident.2 = clusters[!clusters %in% cluster_group],
#                                 test.use="MAST")
# cluster16_markers$gene <- rownames(cluster16_markers)
# cluster16_markers$pct_delta <- cluster16_markers$pct.1 -
#   cluster16_markers$pct.2

# # try this again, with the other uncategorized
# cluster_group <- c("1","16","13","2")
# 
# cluster_group_markers <- lapply(cluster_group, function(cluster) {
#   
#   print(cluster)
#   
#   cluster_markers <- FindMarkers(seu_subset_sub,
#                                   ident.1 = cluster,
#                                   ident.2 = clusters[!clusters %in% cluster_group],
#                                 test.use="MAST")
#   cluster_markers$gene <- rownames(cluster_markers)
#   cluster_markers$pct_delta <- cluster_markers$pct.1 -
#     cluster_markers$pct.2
#   cluster_markers$cluster <- cluster
#   
#   return(cluster_markers)
# })
# names(cluster_group_markers) <- cluster_group

# cluster groups from another analysis file
cluster_group_markers <- readRDS("results/microglia_subset_markers/cluster_group_markers.RDS")

# combine!
cluster_group_markers_df <- bind_rows(cluster_group_markers)

summary(cluster_group_markers_df)

markers <- unique(cluster_group_markers_df[cluster_group_markers_df$pct_delta > 0.1,]$gene)

DotPlot(seu_subset, features = markers, idents=clusters) + RotatedAxis()




# attempt at naming them
cluster_xref <- data.frame(microglia_clusters=clusters,
                           microglia_cell_name=c("Homeostatic",
                                                 "Microglia",
                                                 "Microglia",
                                                 "Homeostatic",
                                                 "Homeostatic",
                                                 "Cross Presenting",
                                                 "Homeostatic",
                                                 "DAM",
                                                 "DAM",
                                                 "Homeostatic",
                                                 "Homeostatic",
                                                 "Interferon Responsive",
                                                 "Cross Presenting",
                                                 "Microglia",
                                                 "Cross Presenting",
                                                 "Inflammatory",
                                                 "Microglia"))

# merge this into the metadata
cell_meta$cell_id <- rownames(cell_meta)

cell_meta <- merge(cell_meta,
                   cluster_xref, 
                   by="microglia_clusters")

rownames(cell_meta) <- cell_meta$cell_id

cell_meta <- cell_meta[rownames(seu_subset@meta.data),]

# add back into object
seu_subset <- AddMetaData(seu_subset, 
                          cell_meta$microglia_cell_name, 
                          col.name="microglia_cell_name")

DimPlot(seu_subset, reduction="umap.microglia_pca", 
        group.by= c("microglia_clusters",
                    "microglia_cell_name"), label=T,
        raster = F, repel=T)
ggsave(paste0(out_dir, "microglia_subset_seurat.cluster_names.umap.png"), width=14, height=6)

p1 <- DimPlot(seu_subset, reduction="umap.microglia_pca", group.by= "microglia_clusters",
              label=T, raster = T, label.box = T) + labs(x="UMAP 1", y="UMAP 2", title="Microglia Clusters") + 
  theme(legend.position = "none")

p2 <- DimPlot(seu_subset, reduction="umap.microglia_pca", group.by= "microglia_cell_name",
              label=T, raster = T, label.box = T, repel=T) + labs(x="UMAP 1", y="UMAP 2", 
                                                         title="Microglia Cell Names") + 
  theme(legend.position = "none")

plot_grid(p1,p2, ncol=2)
ggsave(paste0(out_dir, "microglia_subset_seurat.cluster_names.umap.png"), width=12, height=6)

write.xlsx(cluster_xref, file=paste0(out_dir, "microglia_subcluster_celltypes.xlsx"), colWidths="auto")



#SaveSeuratRds(seu_subset, file=paste0(out_dir, "cell_named.microglia_subset.seurat.RDS"))

seu_subset <- LoadSeuratRds(file=paste0(out_dir, "cell_named.microglia_subset.seurat.RDS"))

saveRDS(seu_subset@meta.data, file=paste0(out_dir, "cell_named.microglia_subset.cell_metadata.RDS"))

# some basic breakdowns

sample_ids <- unique(sort(cell_meta$orig.ident))

gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", sample_ids)

cell_meta$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", cell_meta$orig.ident)

cell_meta$condition <- factor(as.character(cell_meta$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

seu_subset <- AddMetaData(seu_subset, 
                          cell_meta$condition, 
                          col.name="condition")


DimPlot(seu_subset, reduction="umap.microglia_pca", 
        group.by= "microglia_cell_name",
        split.by = "condition",
        ncol = 2,
        label=T, raster = F, label.box = F, repel=T) + labs(x="UMAP 1", y="UMAP 2", 
                                                         title="Microglia Cell Names")


ggsave(paste0(out_dir, "condition.microglia_cluster_names.umap.png"), width=12, height=9)

# bar plots of cell counts!

cell_meta$sample <- cell_meta$orig.ident

# factorize for plotting
sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")

cell_meta$sample <- factor(as.character(cell_meta$sample),
                           levels=sample_levels)

ggplot(cell_meta,
       aes(x=microglia_cell_name,
           fill=sample)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave(paste0(out_dir, "microglia_subset.cell_count.barplot.png"),
       width=6, height=4)

ggplot(cell_meta,
       aes(x=microglia_clusters,
           fill=sample)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x=NULL, y="Cell Count") +
  theme(axis.text.x = element_text(angle=35, hjust=1))
ggsave(paste0(out_dir, "microglia_subset.cluster_count.barplot.png"),
       width=8, height=4)


