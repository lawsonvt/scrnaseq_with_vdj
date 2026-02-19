library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)
library(pbayes)
library(reshape2)

out_dir <- "results/tcell_split.post_pericytes.diff_exp_analysis/"
dir.create(out_dir, showWarnings = F)

# load in named seurat object
seu_obj <- readRDS("results/tcell_split.post_pericytes_split/tcell_cell_named.seurat.RDS")

# convert cell names
seu_obj$tcell_sub_cell_name <- gsub("\\+", "p", seu_obj$tcell_sub_cell_name)

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

# start pulling together data
cells <- unique(metadata$tcell_sub_cell_name)

# aggregate the seurat object
seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("tcell_sub_cell_name","orig.ident")
)

# Extract the count matrix
pseudobulk_matrix <- seu_pseudo$RNA

# Extract unique sample-level metadata
sample_metadata <- metadata %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample

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
  
  # calculate posterior probability
  res_df$post_p <- pbayes(res_df$pvalue, n_cores=2)$posterior_prob
  
  return(res_df)
})
names(ps19ko_minus_ps19wt.de_results) <- cells

lapply(ps19ko_minus_ps19wt.de_results, head)

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



# compare between cell types

# WT

# Extract unique sample-level metadata
sample_cell_metadata <- metadata %>%
  select(sample, condition, tcell_sub_cell_name) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_cell_metadata) <- paste0(sample_cell_metadata$tcell_sub_cell_name, "_",
                                         sample_cell_metadata$sample)

all(rownames(sample_cell_metadata) %in% colnames(pseudobulk_matrix))



cell1 <- "Microglia Assoc."
cell2 <- "CD8p"

conditions <- c("PS19-WT", "PS19-KO")

cell_degs <- lapply(conditions, function(condition) {

  meta <- sample_cell_metadata[sample_cell_metadata$condition == condition &
                                 sample_cell_metadata$tcell_sub_cell_name %in% c(cell1,cell2),]
  
  
  # Subset counts
  counts <- pseudobulk_matrix[, rownames(meta)]
  
  # Keep genes with at least 10 counts in at least 2 samples
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, ]
  
  meta$tcell_sub_cell_name <- factor(meta$tcell_sub_cell_name,
                                     levels=c(cell1,cell2))
  
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ tcell_sub_cell_name
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, contrast = c("tcell_sub_cell_name", cell1, cell2))
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(condition = condition,
           comparison = paste0(cell1, "_m_", cell2))
  
  # filter out NAs
  res_df <- res_df[!is.na(res_df$padj),]
  
  # calculate posterior probability
  res_df$post_p <- pbayes(res_df$pvalue, n_cores=2)$posterior_prob
  
  return(res_df)
})
names(cell_degs) <- conditions

# make volcano plots

volcano_dir <- paste0(out_dir, "cell_comp_volcano_plots/")

dir.create(volcano_dir, showWarnings = F)


top_genes <- 25

for (condition in names(cell_degs)) {
  
  
  subset <- cell_degs[[condition]]
  
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
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=condition,
         subtitle=paste0(cell1, " - ", cell2))
  ggsave(paste0(volcano_dir, to_snake_case(condition), ".volcano_plot.png"), width=5, height=6)
  
  
}




# merge them
merged_cell_degs <- merge(cell_degs[[1]],
                          cell_degs[[2]],
                          by=c("gene","comparison"),
                          suffixes=c(".wt", ".ko"))

merged_cell_degs$wt_p <- merged_cell_degs$post_p.wt * 
  (1 - merged_cell_degs$post_p.ko)

merged_cell_degs$ko_p <- merged_cell_degs$post_p.ko * 
  (1 - merged_cell_degs$post_p.wt)

merged_cell_degs$rsi <- merged_cell_degs$post_p.wt * 
  merged_cell_degs$post_p.wt * 
  sign(merged_cell_degs$log2FoldChange.wt) *
  sign(merged_cell_degs$log2FoldChange.ko)


top_ko <- merged_cell_degs[order(merged_cell_degs$ko_p,
                                 decreasing=T),]$gene[1:6] 


# make a subset for plotting

subset_seu <- subset(seu_obj, subset = condition %in% conditions &
                       tcell_sub_cell_name %in% c(cell1, cell2))


VlnPlot(subset_seu, features=top_ko,
        split.by="tcell_sub_cell_name",
        group.by="condition") +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "vlnplot.mga_m_cd8p.top_ko.png"), width=10, height=8)


top_wt <- merged_cell_degs[order(merged_cell_degs$wt_p,
                                 decreasing=T),]$gene[1:6] 


VlnPlot(subset_seu, features=top_wt,
        split.by="tcell_sub_cell_name",
        group.by="condition")+
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "vlnplot.mga_m_cd8p.top_wt.png"), width=10, height=8)

# pull out pseudobulk counts


subset_meta <- sample_cell_metadata[sample_cell_metadata$condition %in% conditions &
                                      sample_cell_metadata$tcell_sub_cell_name %in% c(cell1,cell2),]


counts <- pseudobulk_matrix[,rownames(subset_meta)]


keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]

subset_meta$tcell_sub_cell_name <- factor(subset_meta$tcell_sub_cell_name,
                                   levels=c(cell1,cell2))


# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = subset_meta,
  design = ~ tcell_sub_cell_name
)

# Run DESeq2
dds <- DESeq(dds)

vsd <- assay(vst(dds, blind=F))

vsd_long <- melt(vsd)
colnames(vsd_long) <- c("gene", "sample_cell","value")


# merge in metadata
subset_meta$sample_cell <- rownames(subset_meta)

vsd_long <- merge(vsd_long, subset_meta, by="sample_cell")


top_ko <- merged_cell_degs[order(merged_cell_degs$ko_p,
                                 decreasing=T),]$gene[1:6] 


ggplot(vsd_long[vsd_long$gene %in% top_ko,],
       aes(x=condition,
           y=value,
           color=tcell_sub_cell_name)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge(jitter.width=0.1)) +
  theme_bw() +
  facet_wrap(~ gene, ncol=3, scales="free_y") +
  labs(x=NULL, y="Aggregate, Normalized Value", color="Cell")
ggsave(paste0(out_dir, "boxplot.mga_m_cd8p.top_ko.png"), width=8, height=6)


top_wt <- merged_cell_degs[order(merged_cell_degs$wt_p,
                                 decreasing=T),]$gene[1:6] 


ggplot(vsd_long[vsd_long$gene %in% top_wt,],
       aes(x=condition,
           y=value,
           color=tcell_sub_cell_name)) +
  geom_boxplot() +
  geom_jitter(position = position_jitterdodge(jitter.width=0.1)) +
  theme_bw() +
  facet_wrap(~ gene, ncol=3, scales="free_y") +
  labs(x=NULL, y="Aggregate, Normalized Value", color="Cell")
ggsave(paste0(out_dir, "boxplot.mga_m_cd8p.top_wt.png"), width=8, height=6)
