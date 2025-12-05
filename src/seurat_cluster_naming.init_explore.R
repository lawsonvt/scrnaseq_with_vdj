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
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"))

# additional marker genes taken from ssreads list
# source: https://bmblx.bmi.osumc.edu/ssread/help/methods#marker-genes

ssread_cluster_markers <- list("Astrocytes"=c("Gfap", "Aqp4", "Gja1", "Slc1a2","Fgfr3",
                                              "Nkain4", "Agt", "Plxnb1", "Slc1a3"),
                               "Endothelial Cells"=c("Cldn5", "Vwf"),
                               "Neurons"=c("Gls", "Rbfox3", "Camk2a"),
                               "Microglia"=c("P2ry12", "Csf1r", "Cx3cr1", "C3"),
                               "Oligodendrocytes"=c("Mbp", "Mobp", "Plp1", "Myrf", "Mag"),
                               "Oligodendrocyte Precursor Cells"=c("Vcan", "Sox8"),
                               "Pericytes"=c("Ambp", "Higd1b", "Pth1r"),
                               "Excitatory Neurons"=c("Slc17a6", "Slc17a7", "Satb2"),
                               "Inhibitory Neurons"=c("Gad1", "Gad2"))

# read in integrated seurat
int_seu <- LoadSeuratRds("results/integrate_seurat_samples/all_samples.integrated_seurat.RDS")

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

analysis_genes <- rownames(int_seu)

# find any issues between marker genes and analysis genes

lab_cluster_markers_missing <- lapply(lab_cluster_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

ssread_cluster_markers_missing <- lapply(ssread_cluster_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})
# missing genes seem to be just missing and not typos

# lets do some plots!

for (cluster_name in names(lab_cluster_markers)) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  outfile_stump <- paste0(out_dir, to_snake_case(cluster_name))
  
  #RidgePlot(int_seu, features=markers)
  
  VlnPlot(int_seu, features=markers, ncol=1, raster=F)
  ggsave(paste0(outfile_stump, ".violin_plot.png"), width=10, height=8)
  FeaturePlot(int_seu, features=markers, reduction = "umap.harmony", ncol=3)
  ggsave(paste0(outfile_stump, ".feature_umap_plot.png"), width=12, height=4)
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL)
  ggsave(paste0(outfile_stump, ".dotplot.png"), width=6, height=8, bg="white")
  
}

# cow dot plot


dotplot_list <- lapply(names(lab_cluster_markers), function(cluster_name) {
  
  markers <- lab_cluster_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "total_dotplots.png"), width=18, height=7, bg="white")

# convert marker list to dataframe
lab_cluster_markers_df <- bind_rows(lapply(names(lab_cluster_markers), function(cluster_name) {
  
  data.frame(cluster_name=cluster_name,
             marker_gene=lab_cluster_markers[[cluster_name]])
  
}))

DotPlot(int_seu, lab_cluster_markers_df$marker_gene)
ggsave(paste0(out_dir, "simple_total_dotplots.png"), width=18, height=7, bg="white")

# bulk expression heatmap

agg_seu <- AggregateExpression(int_seu, group.by="harmony_clusters")
agg_seu_rna <- t(as.matrix(agg_seu$RNA))
# row scale it
agg_seu_rna_s <- scale(agg_seu_rna)



# make heatmap

subset <- agg_seu_rna_s[,lab_cluster_markers_df$marker_gene]

cluster_colors <- brewer.pal(length(lab_cluster_markers),"Set1")
names(cluster_colors) <- names(lab_cluster_markers)

ha <- HeatmapAnnotation(cell_type=lab_cluster_markers_df$cluster_name,
                        col=list(cell_type=cluster_colors))

Heatmap(subset,
        col = colorRamp2(c(-1, 0, 1), c("blue","white","red")),
        cluster_columns = F,
        bottom_annotation = ha,
        name = "Scaled\nAggregate\nExpression")


# these markers do not work for us
# dotplot_list <- lapply(names(ssread_cluster_markers), function(cluster_name) {
#   
#   markers <- ssread_cluster_markers[[cluster_name]]
#   
#   DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, title=cluster_name) +
#     theme(legend.position="none")
#   
# })
# 
# plot_grid(plotlist = dotplot_list, nrow = 1)

# based on looking over the cluster marker dot plots, here is what I got

cluster_xref <- data.frame(harmony_clusters=0:27,
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
                                          "Unknown",
                                          "Neutrophils",
                                          "B Cells",
                                          "Monocytes",
                                          "Proliferating Cells",
                                          "Monocytes",
                                          "Unknown",
                                          "Macrophages",
                                          "Monocytes",
                                          "Proliferating Cells",
                                          "B Cells",
                                          "Proliferating Cells"))

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
ggsave(paste0(out_dir, "cluster_cells_umap.png"), width=14, height=6)


# save things
write.xlsx(cluster_xref, paste0(out_dir, "cluster_xref.xlsx"), colWidths="auto")

SaveSeuratRds(int_seu, paste0(out_dir, "all_samples.integrated_seurat.cell_cluster_names.RDS"))

table(metadata$orig.ident,
      metadata$cell_cluster)

table(metadata$orig.ident,
      metadata$harmony_clusters)

DimPlot(int_seu, reduction="umap.harmony", group.by="cell_cluster",
        label=T, split.by = "orig.ident", ncol=4, raster=F)
ggsave(paste0(out_dir, "cluster_cells_umap.per_sample.png"), width=17, height=10)


DimPlot(int_seu, reduction="umap.harmony", group.by="harmony_clusters",
        label=T, split.by = "orig.ident", ncol=4, raster=F)
ggsave(paste0(out_dir, "harmony_clusters_umap.per_sample.png"), width=17, height=10)

# additional marker list (see what else is there)

add_markers <- list(Neurons=c("Rbfox3", "Map2"),
                    Oligodendrocytes=c("Mbp", "Olig2"),
                    Astrocytes=c("Gfap", "Aldh1l1"),
                    "Endothelial Cells"=c("Pecam1", "Cdh5"),
                    Pericytes=c("Pdgfrb", "Rgs5"),
                    OPCs=c("Pdgfra", "Cspg4"))

add_markers_missing <- lapply(add_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})
# none missing
# cow dot plot


dotplot_list <- lapply(names(add_markers), function(cluster_name) {
  
  markers <- add_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "extra_markers_dotplots.png"), width=14, height=7, bg="white")

# looks like we have some pericytes too!

# but what about cluster 16?

Idents(int_seu) <- "harmony_clusters"

c16_degs <- FindMarkers(int_seu,
                        ident.1 = "16",
                        test.use="wilcox")






