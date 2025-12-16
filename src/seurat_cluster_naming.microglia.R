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

out_dir <- "results/seurat_cluster_naming.microglia/"
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
                                           "Cd68"))

# read in integrated seurat
int_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

analysis_genes <- rownames(int_seu)

# find any issues between marker genes and analysis genes

micro_markers_missing <- lapply(micro_markers, function(gene_list) {
  
  gene_list[!gene_list %in% analysis_genes]
  
})

# get cell metadata

cell_meta <- int_seu@meta.data

micro_clusters <- unique(cell_meta[cell_meta$cell_cluster == "Microglia",]$harmony_clusters)
rest_clusters <- unique(cell_meta[cell_meta$cell_cluster != "Microglia",]$harmony_clusters)

# cow dot plot for the markers
dotplot_list <- lapply(names(micro_markers), function(cluster_name) {
  
  markers <- micro_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers, idents=micro_clusters) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "microglia_dotplots.png"), width=18, height=7, bg="white")

# cow dot plot for the markers
dotplot_list <- lapply(names(micro_markers), function(cluster_name) {
  
  markers <- micro_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers, idents=rest_clusters) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "non-microglia_dotplots.png"), width=18, height=7, bg="white")

# cow dot plot for the markers
dotplot_list <- lapply(names(micro_markers), function(cluster_name) {
  
  markers <- micro_markers[[cluster_name]]
  
  DotPlot(int_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none")
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "total_dotplots.png"), width=18, height=7, bg="white")


# try out some feature plots
FeaturePlot(int_seu, 
            features = micro_markers[["Cross presenting"]],
            reduction="umap.harmony")

FeaturePlot(int_seu, 
            features = c("Wdfy3","Wdfy4"),
            reduction="umap.harmony", raster = F)

FeaturePlot(int_seu, 
            features = micro_markers[["DAM"]],
            reduction="umap.harmony", raster = F)



