library(Seurat)
library(ggplot2)

out_dir <- "results/inspect_cluster_marker_tests/"
dir.create(out_dir, showWarnings = F)

# read in integrated seurat
int_seu <- LoadSeuratRds("results/integrate_seurat_samples/all_samples.integrated_seurat.RDS")

all_genes <- rownames(int_seu)

mt_genes <- grep("^mt", all_genes, value=T)

# set assay to RNA and join layers
DefaultAssay(int_seu) <- "RNA"
int_seu <- JoinLayers(int_seu)

int_seu <- NormalizeData(int_seu)
int_seu <- FindVariableFeatures(int_seu)
int_seu <- ScaleData(int_seu)

# pull in cluster DEGs
c16_degs <- readRDS("results/cluster_marker_tests/cluster16_degs.wilcox.RDS")
c16_degs$pct_delta <- c16_degs$pct.1 - c16_degs$pct.2

c16_degs <- c16_degs[order(c16_degs$pct_delta, decreasing = T),]

Idents(int_seu) <- "harmony_clusters"

VlnPlot(int_seu, features = rownames(c16_degs)[1:4], ncol=1, raster=F, alpha=0.2)

DotPlot(int_seu, features = rownames(c16_degs)[1:20]) + RotatedAxis()

VlnPlot(int_seu, features = "mt-Nd6", alpha=0.2) + theme(legend.position = "none") +
  labs(y="Harmony Clusters")
ggsave(paste0(out_dir, "mt_Nd6.violinplot.png"), width=8, height=5)

DotPlot(int_seu, features = mt_genes) + RotatedAxis() + labs(x="Mitochondrial Genes", y="Harmony Clusters")
ggsave(paste0(out_dir, "mitochondria_genes.dotplot.png"), width=8, height=7, bg="white")

c22_degs <- readRDS("results/cluster_marker_tests/cluster22_degs.wilcox.RDS")
c22_degs$pct_delta <- c22_degs$pct.1 - c22_degs$pct.2

c22_degs <- c22_degs[order(c22_degs$pct_delta, decreasing = T),]

peri_markers <- c("Pdgfrb","Rgs5")

c22_degs[peri_markers,]
