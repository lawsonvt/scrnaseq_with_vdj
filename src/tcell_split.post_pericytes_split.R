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

out_dir <- "results/tcell_split.post_pericytes_split/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming/cell_named.seurat.RDS")

# merge in TCR clonal data
combined_tcr <- readRDS("results/tcr_analysis.subset_integrated/combined_tcr.RDS")

seu_obj <- combineExpression(combined_tcr,
                             seu_obj,
                               cloneCall="gene",
                               group.by = "sample",
                               proportion = F,
                               cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

total_meta <- seu_obj@meta.data

# get pericytes metadata
peri_meta <- readRDS("results/pericytes_split/pericytes_subset.named.metadata.RDS")
peri_cd8s <- rownames(peri_meta[peri_meta$pericyte_type == "CD8+",])

tcell_cells <- rownames(total_meta[total_meta$cell_cluster == "T Cells",])

# subset to the tcell ones
tcell_seu <- subset(seu_obj, cells = c(peri_cd8s,
                                       tcell_cells))

tcell_seu <- quietTCRgenes(tcell_seu)

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
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_clusters")


# run doublet finder
tcell_sce <- as.SingleCellExperiment(tcell_seu)

bp <- MulticoreParam(2, RNGseed=1234) # equivalent to set seed, for reproducibility
tcell_sce <- scDblFinder(tcell_sce, clusters="tcell_clusters", BPPARAM=bp)

# add doublet calls to seurat object
tcell_seu$scDblFinder.class <- tcell_sce$scDblFinder.class

VlnPlot(tcell_seu, group.by="orig.ident", split.by = "scDblFinder.class",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
ggsave(paste0(out_dir, "tcell.doublet_qc_plots.png"), width=9, height=6)

doublet_counts <- as.data.frame(table(tcell_seu$scDblFinder.class))

ggplot(doublet_counts, 
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", fill="grey", color="black") +
  geom_text(aes(label=Freq), vjust=-0.4) +
  theme_bw() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "tcell.doublet_counts.png"), width=4, height=5)

dbl_meta <- as.data.frame(colData(tcell_sce))

saveRDS(dbl_meta, file=paste0(out_dir, "tcell.doublet_output.RDS"))

# remove doublets, redo clustering
tcell_seu <- subset(tcell_seu, scDblFinder.class == "singlet")

tcell_seu <- RunPCA(tcell_seu, npcs = 20)

# inspect elbow plot
ElbowPlot(tcell_seu, ndims=20) + 
  labs(title="TCell Subset")
ggsave(paste0(out_dir, "tcell_subset.pca_elbow_plot.post_dbl_removal.png"), width=6, height=5, bg="white")

# max dims
max_pc_dim <- 10

# cluster the harmonized data
tcell_seu <- FindNeighbors(tcell_seu, dims = 1:max_pc_dim, reduction = "pca")
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_clusters")

# create umap
tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("orig.ident","tcell_clusters"))

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("cell_cluster","tcell_clusters"))


# try different resolutions for clustering
tcell_seu <- FindClusters(tcell_seu, resolution = c(0.3, 0.5, 0.8, 1.0, 1.2))

clustree(tcell_seu, prefix = "RNA_snn_res.")
ggsave(paste0(out_dir, "res_clustering_tree.png"), width=9, height=8)

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

# make dot plots for each resolution
res_cols <- grep("RNA_snn_res", colnames(tcell_seu@meta.data), value=T)

for (res in res_cols) {
  
  Idents(tcell_seu) <- res
  
  dotplot_list <- lapply(names(markers), function(cluster_name) {
    
    markers <- markers[[cluster_name]]
    
    DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
      
  })
  
  plot_grid(plotlist = dotplot_list, nrow = 3)
  
  ggsave(paste0(out_dir, res, ".dot_plot.png"), width=16, height=12, bg="white")
  
  
}

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= res_cols,
        ncol=3, label=T)
ggsave(paste0(out_dir, "clustering_umaps.png"), width=12, height=8)

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= "RNA_snn_res.0.5",
        label = T)

DimPlot(tcell_seu, group.by = c("RNA_snn_res.0.5","cloneSize"), reduction = "umap.tcell_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

# in between labels

tcell_meta <- tcell_seu@meta.data

# tcell_clusters <- data.frame(RNA_snn_res.0.5=sort(unique(tcell_meta$RNA_snn_res.0.5)),
#                              init_cell_name=c("CD8+",
#                                               "CD4+",
#                                               "Microglia",
#                                               "Unknown3",
#                                               "Unknown4 (T Cells)",
#                                               "Unknown5",
#                                               "Unknown6",
#                                               "B Cells", 
#                                               "Monocytes",
#                                               "CD8+",
#                                               "Microglia"))

clonal_counts <- lapply(res_cols, function(res) {
  
  counts <- as.data.frame(table(tcell_meta[!is.na(tcell_meta$cloneSize),res]))
  colnames(counts) <- c("cluster", "count")
  
  counts$type <- "cells with clones"
  
  no_counts <- as.data.frame(table(tcell_meta[is.na(tcell_meta$cloneSize),res]))
  colnames(no_counts) <- c("cluster", "count")
  
  no_counts$type <- "cells with no clones"
  
  counts <- rbind(counts,
                  no_counts)

  counts$cluster <- as.character(counts$cluster)
  counts$res <- res
  
  return(counts)
})
clonal_counts <- bind_rows(clonal_counts)

clonal_counts$cluster <- factor(clonal_counts$cluster,
                                levels=mixedsort(unique(clonal_counts$cluster)))

ggplot(clonal_counts,
       aes(x=cluster,
           y=count,
           fill=type)) +
  geom_bar(stat="identity") +
  facet_wrap(~ res, scales="free", ncol=3) +
  theme_bw()
ggsave(paste0(out_dir, "cluster_clonal_counts.png"), width=11, height=7)

# clonal size counts
clonal_size_counts <- lapply(res_cols, function(res) {
  
  counts <- as.data.frame(table(tcell_meta[,res],
                                tcell_meta$cloneSize))
  colnames(counts) <- c("cluster", "clone_size","count")
  
  counts$cluster <- as.character(counts$cluster)
  counts$res <- res
  
  return(counts)
})
clonal_size_counts <- bind_rows(clonal_size_counts)

ggplot(clonal_size_counts,
       aes(x=cluster,
           y=count,
           fill=clone_size)) +
  geom_bar(stat="identity", color="black") +
  facet_wrap(~ res, scales="free", ncol=3) +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(3,4,5,7)]))
ggsave(paste0(out_dir, "cluster_clonal_size_counts.png"), width=11, height=7)




# try naming the clusters
tcell_clusters <- data.frame(RNA_snn_res.0.5=sort(unique(tcell_meta$RNA_snn_res.0.5)),
                             init_cell_name=c("CD8+",
                                              "CD4+",
                                              "Unknown2",
                                              "Microglia",
                                              "Unknown4 (T Cells)",
                                              "Unknown5",
                                              "Unknown6",
                                              "B Cells",
                                              "CD8+"))

# merge them in
tcell_meta$cell_id <- rownames(tcell_meta)

tcell_meta <- merge(tcell_meta,
                    tcell_clusters,
                    by="RNA_snn_res.0.5")

rownames(tcell_meta) <- tcell_meta$cell_id

tcell_meta <- tcell_meta[colnames(tcell_seu),]

tcell_seu$init_cell_name <- tcell_meta$init_cell_name

DimPlot(tcell_seu, group.by = c("init_cell_name","cloneSize"), reduction = "umap.tcell_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

DimPlot(tcell_seu, group.by = "init_cell_name", reduction = "umap.tcell_pca", label=T) 


# find cluster markers

Idents(tcell_seu) <- "RNA_snn_res.0.5"

cluster_markers <- lapply(sort(unique(tcell_meta$RNA_snn_res.0.5)), function(cluster) {
  
  markers <- FindMarkers(tcell_seu,
                         ident.1 = cluster,
                         test.use="wilcox")
  markers <- rownames_to_column(markers, var="gene")
  
  markers$delta_pct <- markers$pct.1 - markers$pct.2
  markers$cluster <- cluster
  
  markers <- markers[order(markers$delta_pct, decreasing = T),]
  
  return(markers)
})
names(cluster_markers) <- sort(unique(tcell_meta$RNA_snn_res.0.5))

lapply(cluster_markers, head)


write.xlsx(cluster_markers, paste0(out_dir, "cluster_markers.xlsx"),
           colWidths="auto")

top_num <- 5

dotplot_list <- lapply(cluster_markers, function(data) {
  
  cluster_name <- unique(data$cluster)
  
  markers <- data[1:top_num,]$gene
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)
ggsave(paste0(out_dir, "cluster_marker.dot_plots.png"), width=16, height=12, bg="white")




# microglial markers

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





dotplot_list <- lapply(names(micro_markers), function(cluster_name) {
  
  markers <- micro_markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)

ggsave(paste0(out_dir, res, ".dot_plot.png"), width=16, height=12, bg="white")

# finally put in the new and improved cell names
tcell_clusters <- data.frame(RNA_snn_res.0.5=sort(unique(tcell_meta$RNA_snn_res.0.5)),
                             tcell_sub_cell_name=c("CD8+",
                                                   "CD4+",
                                                   "Unknown2",
                                                   "Microglia",
                                                   "Naive T Cells",
                                                   "Natural Killer",
                                                   "Natural Killer",
                                                   "B Cells",
                                                   "CD8+"))


tcell_meta <- merge(tcell_meta,
                    tcell_clusters,
                    by="RNA_snn_res.0.5")

rownames(tcell_meta) <- tcell_meta$cell_id

tcell_meta <- tcell_meta[colnames(tcell_seu),]

tcell_seu$tcell_sub_cell_name <- tcell_meta$tcell_sub_cell_name

# add in condition as well
tcell_meta$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", tcell_meta$orig.ident)

tcell_meta$condition <- factor(as.character(tcell_meta$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

tcell_seu$condition <- tcell_meta$condition

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= c("RNA_snn_res.0.5",
                    "tcell_sub_cell_name"),
        label = T)
ggsave(paste0(out_dir, "cluster2cell_umap.png"), width=11, height=6)

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= c("RNA_snn_res.0.5"),
        label = T)
ggsave(paste0(out_dir, "cluster_umap.png"), width=6, height=6)

DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= c("tcell_sub_cell_name"),
        label = T)
ggsave(paste0(out_dir, "cell_umap.png"), width=6, height=6)


DimPlot(tcell_seu, reduction="umap.tcell_pca", 
        group.by= "tcell_sub_cell_name",
        split.by = "condition",
        ncol = 2,
        label = T)
ggsave(paste0(out_dir, "cluster2cell_umap.per_condition.png"), width=8, height=6)

# clonal differences


counts <- as.data.frame(table(tcell_meta[!is.na(tcell_meta$cloneSize),]$tcell_sub_cell_name))
colnames(counts) <- c("cell", "count")

counts$type <- "cells with clones"

no_counts <- as.data.frame(table(tcell_meta[is.na(tcell_meta$cloneSize),]$tcell_sub_cell_name))
colnames(no_counts) <- c("cell", "count")

no_counts$type <- "cells with no clones"

counts <- rbind(counts,
                no_counts)


ggplot(counts,
       aes(x=cell,
           y=count,
           fill=type)) +
  geom_bar(stat="identity", color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count", fill=NULL)
ggsave(paste0(out_dir, "cell_names.clonal_counts_bar_plot.png"), width=6, height=4)

# add in condition

counts <- as.data.frame(table(tcell_meta[!is.na(tcell_meta$cloneSize),]$tcell_sub_cell_name,
                              tcell_meta[!is.na(tcell_meta$cloneSize),]$condition))
colnames(counts) <- c("cell", "condition","count")

counts$type <- "cells with clones"

no_counts <- as.data.frame(table(tcell_meta[is.na(tcell_meta$cloneSize),]$tcell_sub_cell_name,
                                 tcell_meta[is.na(tcell_meta$cloneSize),]$condition))
colnames(no_counts) <- c("cell", "condition", "count")

no_counts$type <- "cells with no clones"

counts <- rbind(counts,
                no_counts)

ggplot(counts,
       aes(x=cell,
           y=count,
           fill=type)) +
  geom_bar(stat="identity", color="black") +
  facet_wrap(~ condition, ncol=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count", fill=NULL)
ggsave(paste0(out_dir, "cell_names.clonal_counts_bar_plot.per_condition.png"), width=10, height=8)




counts <- as.data.frame(table(tcell_meta$tcell_sub_cell_name,
                              tcell_meta$cloneSize,
                              tcell_meta$condition))
colnames(counts) <- c("cell", "clone_size","condition","count")


ggplot(counts,
       aes(x=cell,
           y=count,
           fill=clone_size)) +
  facet_wrap(~ condition, ncol=2) +
  geom_bar(stat="identity", color="black") +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(3,4,5,7)])) +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count", fill="Clone Size")
ggsave(paste0(out_dir, "cell_names.clonal_size_counts_bar_plot.per_condition.png"), width=10, height=8)

# one last thing: differential expression analysis for Tcells

tcell_seu$tcell_sub_cell_name <- gsub("\\+", "p", tcell_seu$tcell_sub_cell_name)

# might as well do all of them
cells <- unique(tcell_seu$tcell_sub_cell_name)

sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")


tcell_meta$sample <- factor(as.character(tcell_meta$orig.ident),
                          levels=sample_levels)

# drop microglia
# cells <- cells[cells != "Microglia"]

# aggregate the seurat object

tcell_seu_pseudo <- AggregateExpression(
  tcell_seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("tcell_sub_cell_name","orig.ident")
)

# Extract the count matrix
pseudobulk_matrix <- tcell_seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- tcell_meta %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

# DE workflows for each comparison

# PS19 WT v KO

condition1 <- "PS19-KO"
condition2 <- "PS19-WT"

ps19ko_minus_ps19wt.de_results <- lapply(cells, function(cell) {
  
  cat("Processing:", cell, "\n")
  
  celltype_cols <- grep(paste0("^", cell, "_"), 
                        colnames(pseudobulk_matrix), 
                        value = TRUE)
  
  
  # Extract sample names from column names
  sample_names <- sub(paste0("^", cell, "_"), "", celltype_cols)
  
  # Subset counts
  counts <- pseudobulk_matrix[, celltype_cols]
  colnames(counts) <- sample_names
  
  meta <- sample_metadata[sample_names, , drop = FALSE]
  
  # Filter to only the two conditions being compared
  samples_to_keep <- meta$condition %in% c(condition1, condition2)
  counts <- counts[, samples_to_keep, drop = FALSE]
  meta <- meta[samples_to_keep, , drop = FALSE]
  
  # Keep genes with at least 10 counts in at least 2 samples
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, ]
  
  meta$condition <- factor(as.character(meta$condition), levels = c(condition1, condition2))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ condition
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("condition", condition1, condition2))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(cell_type = cell,
           comparison = paste0(condition1, "_m_", condition2))
  
  # filter out NAs
  res_df <- res_df[!is.na(res_df$padj),]
  
  return(res_df)
})
names(ps19ko_minus_ps19wt.de_results) <- cells

volcano_dir <- paste0(out_dir, "volcano_plots/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

for (cell in names(ps19ko_minus_ps19wt.de_results)) {
  
  
  subset <- ps19ko_minus_ps19wt.de_results[[cell]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.05 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  ggplot(subset,
         aes(x=log2FoldChange,
             y=log_p)) +
    geom_point(alpha=0.4, color="black") +
    geom_hline(yintercept = logp_thresh,
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red") +
    geom_text_repel(data=subset_top,
                    aes(label=gene),
                    color="red", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=cell,
         subtitle="PS19-KO - PS19-WT")
  ggsave(paste0(volcano_dir, to_snake_case(cell), ".volcano_plot.png"), width=5, height=6)
  
  
}

# make a plot of just CD8+ / CD4+

cd_seu <- subset(tcell_seu, subset = tcell_sub_cell_name == "CD8p" |
                   tcell_sub_cell_name == "CD4p")

cd_seu$tcell_sub_cell_name <- gsub("p", "+", cd_seu$tcell_sub_cell_name)

DimPlot(cd_seu, reduction="umap.tcell_pca", 
        group.by= c("tcell_sub_cell_name"),
        label = T) + theme(legend.position = "None") +
  labs(x="UMAP 1", y="UMAP 2", title=NULL) 
ggsave(paste0(out_dir, "cd4_cd8.umap.png"), width=5, height=4)

DimPlot(cd_seu, reduction="umap.tcell_pca", 
        group.by= c("tcell_sub_cell_name"),
        split.by="condition",
        label = T, ncol=2) + theme(legend.position = "None") +
  labs(x="UMAP 1", y="UMAP 2", title=NULL) 
ggsave(paste0(out_dir, "cd4_cd8.umap.per_condition.png"), width=6, height=5)


DimPlot(cd_seu, reduction="umap.tcell_pca", 
        group.by= "cloneSize",
        split.by="condition",
        ncol=2) +scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL) 
ggsave(paste0(out_dir, "cd4_cd8.umap.clonal_size.per_condition.png"), width=7, height=5)









