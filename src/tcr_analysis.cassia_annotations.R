library(scRepertoire)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(openxlsx)

# load in data
out_dir <- "results/tcr_analysis.cassia_annotations/"

dir.create(out_dir, showWarnings = F)

# read in things
seu_obj <- readRDS("results/tcr_analysis.further_explore/clone_seurat.RDS")

cassia_annot <- read.xlsx("results/cassia_annotate.tcr_clone_cells/results_gemini_lc_summary.mod.xlsx")

combined_tcr <- readRDS("results/tcr_analysis.further_explore/combined_tcr.RDS")

# add annotations to seurat

metadata <- seu_obj@meta.data

cassia_annot_simple <- cassia_annot[,c("Cluster.ID","cell_type")]
colnames(cassia_annot_simple) <- c("clone_clusters","cell_type")

# merge em up
metadata <- merge(metadata,
                  cassia_annot_simple,
                  by="clone_clusters")

# reorder metadata
rownames(metadata) <- metadata$cell_id

metadata <- metadata[colnames(seu_obj),]

# add back into seurat object

seu_obj$cell_type <- metadata$cell_type

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

DimPlot(seu_obj, group.by = c("cell_type","cloneSize"), reduction = "umap.clone_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clone_celltype_size.umap.png"), width=13, height=4)


DimPlot(seu_obj, group.by = c("cloneSize"), reduction = "umap.clone_pca",
        split.by="condition", ncol=2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clone_size.condition_split.umap.png"), width=7, height=5)

DimPlot(seu_obj, group.by = c("cell_type"), reduction = "umap.clone_pca",
        split.by="condition", ncol=2) 
ggsave(paste0(out_dir, "clone_celltype.condition_split.umap.png"), width=8, height=5)

# add cell type into tcr data

combined_tcr <- lapply(combined_tcr, function(data) {
  
  merge(data,
        metadata[,c("cell_id","cell_type",
                    "clonalProportion","clonalFrequency","cloneSize")],
        by.x="barcode", by.y="cell_id")
  
})

clonalProportion(combined_tcr, cloneCall = "gene",
                 group.by = "condition")

clonalProportion(combined_tcr, cloneCall = "gene",
                 group.by = "cell_type")

clonalHomeostasis(combined_tcr,
                  cloneCall = "gene",
                  group.by = "cell_type",
                  cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# make my own plot!

ggplot(metadata,
       aes(y=cell_type,
           fill=cloneSize)) +
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  theme_bw() +
  labs(x="Relative Abundance", y=NULL)
ggsave(paste0(out_dir, "cell_type.clone_size_barplots.png"), width=8, height=5)

ggplot(metadata,
       aes(x=condition,
           fill=cloneSize)) +
  geom_bar(color="black") +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  theme_bw() +
  labs(y="Clone Size Counts", x=NULL)
ggsave(paste0(out_dir, "condition.clone_size_barplots.png"), width=5, height=4)

ggplot(metadata,
       aes(x=condition,
           fill=cloneSize)) +
  geom_bar(color="black") +
  facet_wrap(~ cell_type, ncol=4, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust = 1)) +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x=NULL, y="Clone Size Counts")
ggsave(paste0(out_dir, "cell_type.clone_size_barplots.per_condition.png"), width=12, height=8)


