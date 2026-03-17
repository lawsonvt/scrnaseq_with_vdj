library(Seurat)
library(SeuratObject)
library(ggplot2)
library(RColorBrewer)
library(ggsci)

out_dir <- "results/clustering_plots.microglia/"
dir.create(out_dir, showWarnings = F)

# load in seurat objects
total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.just_microglia.seurat.RDS")

total_seu$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", total_seu$orig.ident)

total_seu$condition <- factor(as.character(total_seu$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))



# get microglia sub markers
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
                                                "Mef2c"))


# nothing unique ....

micro_genes <- unlist(micro_markers)
names(micro_genes) <- NULL

DotPlot(total_seu,
        features=micro_genes) +
  RotatedAxis() +labs(x=NULL, y=NULL)

# custom marker plot
total_seu$merged_cell_name <- factor(total_seu$merged_cell_name,
                                     levels=rev(sort(unique(total_seu$merged_cell_name))))

Idents(total_seu) <- "merged_cell_name"

triple_markers <- c("Apoe","Clec7a","Cst7",
                    "P2ry12","Tmem119","Cx3cr1",
                    "Ifitm3","Mx1","Irf7")


DotPlot(total_seu,
        features=triple_markers) +
  RotatedAxis() +labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "microglia_sub.marker_dot_plot.png"), width=7, height=4,
       bg="white")


# determine unique markers for intermediate microglia markers
# Idents(total_seu) <- "merged_cell_name"
# 
# set.seed(42)
# 
# sampled_cells <- sample(colnames(total_seu), size=12000, replace = F)
# 
# sub_seu <- subset(total_seu, cells=sampled_cells)
# 
# im_markers <- FindMarkers(sub_seu,
#                           ident.1 = "Intermediate Microglia",
#                           test.use = "wilcox")
# im_markers$pct_delta <- im_markers$pct.1 - im_markers$pct.2
# im_markers <- im_markers[order(im_markers$pct_delta,
#                                decreasing=T),]
# # nothing unique ....
# 
# micro_genes <- unlist(micro_markers)
# names(micro_genes) <- NULL
# 
# micro_genes <- c(micro_genes, "Hexb")
# 
# # plot all the markers I have
# DotPlot(total_seu,
#         features=micro_genes) +
#   RotatedAxis() +labs(x=NULL, y=NULL)



DimPlot(total_seu, reduction="umap.harmony",
        group.by = "merged_cell_name", raster=F,
        split.by="condition", ncol=2) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "microglia_clusters.per_condition_umap.png"),
       width=8, height=6)




