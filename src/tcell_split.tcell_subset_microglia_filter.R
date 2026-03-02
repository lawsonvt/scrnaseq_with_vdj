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

out_dir <- "results/tcell_split.tcell_subset_microglia_filter/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/tcell_split.post_pericytes_split/tcell_cell_named.seurat.RDS")

# filter down to those that are T Cells
tcell_init_cells <- c("Naive T Cells",
                      "CD8+",
                      "CD4+",
                      "Proliferating CD8+",
                      "Microglia Assoc.")

tcell_seu <- subset(seu_obj, subset = tcell_sub_cell_name %in% tcell_init_cells)

tcell_meta <- tcell_seu@meta.data


hexb_threshold <- 3.5


# ggplot(tcell_meta,
#        aes(x=tcell_sub_cell_name,
#            y=nFeature_RNA)) +
#   geom_boxplot() +
#   geom_jitter(alpha=0.5, width = 0.2) +
#   theme_bw() + labs(x=NULL, y="RNA Count")
# 
VlnPlot(tcell_seu, group.by="tcell_sub_cell_name", 
         features = c("nFeature_RNA", "nCount_RNA"), 
         ncol = 2, pt.size = 1, alpha = 0.4) 
ggsave(paste0(out_dir, "feature_and_counts_violin_plots.png"), width=8, height=5)

VlnPlot(tcell_seu, group.by="tcell_sub_cell_name",
        features=c("Hexb")) + geom_hline(yintercept = hexb_threshold) 
ggsave(paste0(out_dir, "hexb_expression.png"), width=7, height=5)


# get the Hexb exression
hexb_mat <- as.data.frame(FetchData(tcell_seu, vars="Hexb"))

tcell_meta <- data.frame(tcell_meta,
                          Hexb=hexb_mat[rownames(tcell_meta),])


ggplot(tcell_meta,
       aes(y=Hexb,
           x=nFeature_RNA)) +
  facet_wrap(~ tcell_sub_cell_name, ncol=3) +
  geom_hex() +
  theme_bw() +
  geom_hline(yintercept = hexb_threshold) +
  guides(color="none") +
  labs(x="Gene Count", y="Normalized Hexb Expression")
ggsave(paste0(out_dir, "feature_count_v_hexb_expression.hex.png"), width=8, height=6)


ggplot(tcell_meta,
       aes(x=tcell_sub_cell_name,
           fill=Hexb > hexb_threshold)) +
  geom_bar(position="dodge", color="black") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count", fill=paste0("Hexb > ", hexb_threshold)) +
  scale_y_continuous(breaks=seq(0, 1000, 100))
ggsave(paste0(out_dir, "hexb_threshold_counts.png"), width=7, height=5)

# remove the high Hexb expressing things

tcell_meta$hexb_low <- tcell_meta$Hexb < hexb_threshold

table(tcell_meta$hexb_low)


tcell_seu$hexb_low <- tcell_meta[colnames(tcell_seu),]$hexb_low

tcell_seu <- subset(tcell_seu, subset = hexb_low == T)

# recluster

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
tcell_seu <- FindClusters(tcell_seu, cluster.name = "tcell_clusters_mg_filtered")

# create umap
tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("orig.ident","tcell_clusters_mg_filtered"))

DimPlot(tcell_seu, reduction="umap.tcell_pca", group.by= c("tcell_sub_cell_name","tcell_clusters_mg_filtered"))

DimPlot(tcell_seu, reduction="umap.tcell_pca",group.by= c("tcell_clusters_mg_filtered"),
        label=T, label.box = T) + labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered.umap.png"), width=5, height=4)

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

Idents(tcell_seu) <- "tcell_clusters_mg_filtered"

dotplot_list <- lapply(names(markers), function(cluster_name) {
  
  markers <- markers[[cluster_name]]
  
  DotPlot(tcell_seu, features = markers) + RotatedAxis() + labs(x=NULL, y=NULL, title=cluster_name)
  
})

plot_grid(plotlist = dotplot_list, nrow = 3)

ggsave(paste0(out_dir, "tcell_clusters_mg_filtered.dot_plot.png"), width=16, height=12, bg="white")


VlnPlot(tcell_seu, group.by="tcell_clusters_mg_filtered",
        features="Hexb") 
ggsave(paste0(out_dir, "hexb_expression_counts.mg_filtered.png"), width=8, height=6)

# pull out new metadata
tcell_meta <- tcell_seu@meta.data

Idents(tcell_seu) <- "tcell_clusters_mg_filtered"

# find cluster markers
cluster_markers <- lapply(sort(unique(tcell_meta$tcell_clusters_mg_filtered)), function(cluster) {
  
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

top_markers <- unlist(lapply(cluster_markers, function(x) {head(x$gene,n=10)}))


DotPlot(tcell_seu, features=unique(top_markers)) + RotatedAxis() +
  theme(panel.grid = element_line(color = "grey", linewidth = 0.5),
        axis.text.x = element_text(size=9)) +
  labs(x=NULL, y="T Cell Cluster")
ggsave(paste0(out_dir, "top_markers_per_cluster.png"), width=17, height=8, bg='white')

# look at TCR data
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

ggplot(tcell_meta,
       aes(x=tcell_clusters_mg_filtered,
           fill=cloneSize)) +
  facet_wrap(~ condition, ncol=2) +
  geom_bar(color="black") +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x=NULL, y="Cell Count", fill="Clone Size")
ggsave(paste0(out_dir, "clone_count_per_cluster.png"), width=10, height=8)

# assign new categories
tcell_subset_cluster_xref <- data.frame(tcell_clusters_mg_filtered=sort(unique(tcell_seu$tcell_clusters_mg_filtered)),
                                        tcell_clusters_mg_filtered_name=c("CD8-0",
                                                                     "CD4-1",
                                                                     "CD8-2",
                                                                     "TCELL-3",
                                                                     "CD4-4",
                                                                     "CD8-5",
                                                                     "CD8-6",
                                                                     "TCELL-7",
                                                                     "CD8-8",
                                                                     "CD8-9",
                                                                     "PC-10"))


# merge into Seurat object

tcell_meta <- tcell_seu@meta.data
tcell_meta$cell_id <- rownames(tcell_meta)

tcell_meta <- merge(tcell_meta,
                    tcell_subset_cluster_xref,
                    by = "tcell_clusters_mg_filtered")


rownames(tcell_meta) <- tcell_meta$cell_id

tcell_meta <- tcell_meta[colnames(tcell_seu),]

tcell_seu$tcell_clusters_mg_filtered_name <- tcell_meta$tcell_clusters_mg_filtered_name


DimPlot(tcell_seu, reduction="umap.tcell_pca",group.by= c("tcell_clusters_mg_filtered_name"),
        label=T, label.box = T) + labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_name.umap.png"), width=5, height=4)


# one last thing: differential expression analysis for Tcells

# might as well do all of them
cells <- sort(unique(tcell_seu$tcell_clusters_mg_filtered_name))

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
  group.by = c("tcell_clusters_mg_filtered_name","orig.ident")
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


volcano_dir <- paste0(out_dir, "volcano_plots.tcell_clusters_mg_filtered_name/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(ps19ko_minus_ps19wt.de_results), function(cell) {
  
  
  subset <- ps19ko_minus_ps19wt.de_results[[cell]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.05 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  p1 <- ggplot(subset,
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
  p1
  ggsave(paste0(volcano_dir, to_snake_case(cell), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, nrow = 3)
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_name.volcanoes.png"), width=10, height=8, bg="white")


# group by type

tcell_seu$tcell_clusters_mg_filtered_type <- gsub("\\-[0-9]+", "", tcell_seu$tcell_clusters_mg_filtered_name)


cells <- sort(unique(tcell_seu$tcell_clusters_mg_filtered_type))

# aggregate the seurat object

tcell_seu_pseudo <- AggregateExpression(
  tcell_seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("tcell_clusters_mg_filtered_type","orig.ident")
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


volcano_dir <- paste0(out_dir, "volcano_plots.tcell_clusters_mg_filtered_type/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(ps19ko_minus_ps19wt.de_results), function(cell) {
  
  
  subset <- ps19ko_minus_ps19wt.de_results[[cell]]
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.05 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_top <- subset_sig[1:top_genes,]
  
  if (any(is.na(subset_top$gene))) {
    subset_top <- subset_top[!is.na(subset_top$gene),]
  }
  
  logp_thresh <- min(subset_sig$log_p)
  
  p1 <- ggplot(subset,
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
  p1
  ggsave(paste0(volcano_dir, to_snake_case(cell), ".volcano_plot.png"), width=5, height=6)
  
  return(p1)
})


plot_grid(plotlist = volcano_plot_list, nrow = 2)
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_type.volcanoes.png"), width=7, height=6)


DimPlot(tcell_seu, reduction="umap.tcell_pca",group.by= c("tcell_clusters_mg_filtered_type"),
        label=T, label.box = T) + labs(x="UMAP1", y="UMAP2")
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_type.umap.png"), width=5, height=4)

# drop cluster 6

tcell_seu <- subset(tcell_seu, subset = tcell_clusters_mg_filtered != "6")

tcell_meta <- tcell_seu@meta.data

DimPlot(tcell_seu, reduction="umap.tcell_pca",group.by= c("tcell_clusters_mg_filtered_type"),
        label=T, label.box = T) + labs(x="UMAP1", y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_type.drop_cluster6.umap.png"), width=5, height=4)

DimPlot(tcell_seu, reduction="umap.tcell_pca",group.by= c("tcell_clusters_mg_filtered_name"),
        label=T, label.box = T) + labs(x="UMAP1", y="UMAP2", title=NULL)
ggsave(paste0(out_dir, "tcell_clusters_mg_filtered_name.drop_cluster6.umap.png"), width=5, height=4)


SaveSeuratRds(tcell_seu, file=paste0(out_dir, "tcell.microglia_filted_and_named.seurat.RDS"))

saveRDS(tcell_meta, file=paste0(out_dir, "tcell.microglia_filted_and_named.metadata.RDS"))

