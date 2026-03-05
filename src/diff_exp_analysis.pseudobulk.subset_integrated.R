library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)
library(cowplot)

out_dir <- "results/diff_exp_analysis.pseudobulk.subset_integrated/"
dir.create(out_dir, showWarnings = F)


# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.no_low_count.seurat.RDS")

# replace plus sign
seu_obj$merged_cell_name <- gsub("\\+", "p", seu_obj$merged_cell_name)

metadata <- seu_obj@meta.data

sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")


metadata$sample <- factor(as.character(metadata$orig.ident),
                          levels=sample_levels)

metadata$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", metadata$orig.ident)

metadata$condition <- factor(as.character(metadata$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

# Differential expression of cell types ----------------------------------------

# get the cell types

cells <- unique(metadata$merged_cell_name)

# drop microglia
# cells <- cells[cells != "Microglia"]

# aggregate the seurat object

seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("merged_cell_name","orig.ident")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
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

write.xlsx(ps19ko_minus_ps19wt.de_results, 
           file=paste0(out_dir, "ps19ko_minus_ps19wt.cell_types.de_results.xlsx"),
           colWidths="auto")

saveRDS(ps19ko_minus_ps19wt.de_results,
        file=paste0(out_dir, "ps19ko_minus_ps19wt.cell_types.de_results.RDS"))

# Differential expression of cell categories -----------------------------------

cells <- unique(seu_obj$cell_category)

# aggregate the seurat object

seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("cell_category","orig.ident")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

condition1 <- "PS19-KO"
condition2 <- "PS19-WT"

ps19ko_minus_ps19wt.de_results_2 <- lapply(cells, function(cell) {
  
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
names(ps19ko_minus_ps19wt.de_results_2) <- cells


write.xlsx(ps19ko_minus_ps19wt.de_results_2, 
           file=paste0(out_dir, "ps19ko_minus_ps19wt.cell_categories.de_results.xlsx"),
           colWidths="auto")

saveRDS(ps19ko_minus_ps19wt.de_results_2,
        file=paste0(out_dir, "ps19ko_minus_ps19wt.cell_categories.de_results.RDS"))


# Volcano plots ----------------------------------------------------------------

volcano_dir <- paste0(out_dir, "volcano_plots.cell_types/")

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


plot_grid(plotlist = volcano_plot_list, nrow = 4)
ggsave(paste0(out_dir, "cell_types.volcanoes.png"), width=12, height=10, bg="white")

# Cell categories

volcano_dir <- paste0(out_dir, "volcano_plots.cell_categories/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

volcano_plot_list <-  lapply(names(ps19ko_minus_ps19wt.de_results_2), function(cell) {
  
  
  subset <- ps19ko_minus_ps19wt.de_results_2[[cell]]
  
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
ggsave(paste0(out_dir, "cell_categories.volcanoes.png"), width=10, height=8, bg="white")


# DEG counts plots -------------------------------------------------------------

results_df <- bind_rows(ps19ko_minus_ps19wt.de_results)

# sig results
sig_results_df <- results_df[results_df$padj < 0.05 &
                               abs(results_df$log2FoldChange) > 0.5,]

deg_counts <- as.data.frame(table(sig_results_df$cell_type))
deg_counts <- deg_counts[order(deg_counts$Freq, decreasing = T),]
deg_counts$Var1 <- factor(as.character(deg_counts$Var1),
                          levels=as.character(deg_counts$Var1))

ggplot(deg_counts,
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black", fill="grey") +
  geom_text(aes(label=Freq), vjust=-0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="DEG count at FDR < 0.05", title="PS19-KO - PS19-WT")
ggsave(paste0(out_dir, "ps19ko_minus_ps19wt.cell_types.deg_counts.png"), width=8, height=5)

# cell categories

results_df <- bind_rows(ps19ko_minus_ps19wt.de_results_2)

# sig results
sig_results_df <- results_df[results_df$padj < 0.05 &
                               abs(results_df$log2FoldChange) > 0.5,]

deg_counts <- as.data.frame(table(sig_results_df$cell_type))
deg_counts <- deg_counts[order(deg_counts$Freq, decreasing = T),]
deg_counts$Var1 <- factor(as.character(deg_counts$Var1),
                          levels=as.character(deg_counts$Var1))

ggplot(deg_counts,
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black", fill="grey") +
  geom_text(aes(label=Freq), vjust=-0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="DEG count at FDR < 0.05", title="PS19-KO - PS19-WT")
ggsave(paste0(out_dir, "ps19ko_minus_ps19wt.cell_categories.deg_counts.png"), width=8, height=5)

