library(Seurat)
library(SeuratObject)
library(DESeq2)
library(plyr)
library(dplyr)
library(openxlsx)
library(ggplot2)

out_dir <- "results/diff_exp_analysis.pseudobulk/"
dir.create(out_dir, showWarnings = F)


# load in named seurat object
seu_obj <- readRDS("results/seurat_cluster_naming/cell_named.seurat.RDS")

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

# get the cell types

cells <- unique(metadata$cell_cluster)

# drop microglia
cells <- cells[cells != "Microglia"]

# aggregate the seurat object

seu_pseudo <- AggregateExpression(
  seu_obj,
  assays = "RNA",
  slot = "counts",
  group.by = c("cell_cluster","orig.ident")
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
  
  return(res_df)
})
names(ps19ko_minus_ps19wt.de_results) <- cells


# save results

write.xlsx(ps19ko_minus_ps19wt.de_results, 
           file=paste0(out_dir, "ps19ko_minus_ps19wt.de_results.xlsx"),
           colWidths="auto")

saveRDS(ps19ko_minus_ps19wt.de_results,
        file=paste0(out_dir, "ps19ko_minus_ps19wt.de_results.RDS"))


# plots!

results_df <- bind_rows(ps19ko_minus_ps19wt.de_results)

# sig results
sig_results_df <- results_df[results_df$padj < 0.05,]

table(sig_results_df$gene)

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
  labs(x=NULL, y="DEG count at FDR < 0.05")
ggsave(paste0(out_dir, "ps19ko_minus_ps19wt.deg_counts.png"), width=6, height=5)

