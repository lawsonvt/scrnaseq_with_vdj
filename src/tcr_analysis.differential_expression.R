library(scRepertoire)
library(Seurat)
library(ggplot2)
library(cowplot)
library(DESeq2)
library(dplyr)
library(openxlsx)
library(ggrepel)
library(msigdbr)  # for gene sets / pathways
library(fgsea) # GSEA test

# load in data
out_dir <- "results/tcr_analysis.differential_expression/"

dir.create(out_dir, showWarnings = F, recursive = T)

# read in all the samples
contig_files <- list("A-Q637-PS19-KO"="../processed_samples/A-Q637-PS19-KO/vdj_t/filtered_contig_annotations.csv",
                     "B-Q617-WT"="../processed_samples/B-Q617-WT/vdj_t/filtered_contig_annotations.csv",
                     "C-Q635-KO"="../processed_samples/C-Q635-KO/vdj_t/filtered_contig_annotations.csv",
                     "D-Q619-PS19-WT"="../processed_samples/D-Q619-PS19-WT/vdj_t/filtered_contig_annotations.csv",
                     "E-W136-PS19-WT"="../processed_samples/E-W136-PS19-WT/vdj_t/filtered_contig_annotations.csv",
                     "F-W137-KO"="../processed_samples/F-W137-KO/vdj_t/filtered_contig_annotations.csv",
                     "G-W138-WT"="../processed_samples/G-W138-WT/vdj_t/filtered_contig_annotations.csv",
                     "H-E806-PS19-KO"="../processed_samples/H-E806-PS19-KO/vdj_t/filtered_contig_annotations.csv")

# read in contigs
contig_list <- lapply(contig_files, function(file) {
  
  read.csv(file)
  
})

# going with tutorial https://www.borch.dev/uploads/screpertoire/articles/combining_contigs

# combine contigs into clones
combined_tcr <- combineTCR(contig_list,
                           samples=names(contig_list),
                           removeNA = FALSE, # false is default for all of these params
                           removeMulti = FALSE, 
                           filterMulti = FALSE)


# add in condition
combined_tcr <- addVariable(combined_tcr,
                            variable.name = "condition",
                            variables=c("PS19-KO",
                                        "WT",
                                        "KO",
                                        "PS19-WT",
                                        "PS19-WT",
                                        "KO",
                                        "WT",
                                        "PS19-KO"))




# combine with Seurat objects

total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.seurat.RDS")

# fix name
total_seu$merged_cell_name <- gsub("^Microglia$", "Intermediate Microglia", total_seu$merged_cell_name)
total_seu$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", total_seu$orig.ident)

# broad cell type name
cluster_xref <- data.frame(merged_cell_name=c("Intermediate Microglia",
                                              "Cross Presenting",
                                              "Homeostatic",
                                              "Macrophages",
                                              "Monocytes",
                                              "Interferon Responsive",
                                              "Neutrophils",
                                              "DAM",
                                              "Proliferating Cells",
                                              "CD4+",
                                              "CD8+",
                                              "NK",
                                              "B Cells",
                                              "Pericytes",
                                              "Inflammatory"),
                           broad_cell_name=c("Microglia",
                                             "Microglia",
                                             "Microglia",
                                             "Macrophages",
                                             "Monocytes",
                                             "Microglia",
                                             "Neutrophils",
                                             "Microglia",
                                             "Proliferating Cells",
                                             "T Cells",
                                             "T Cells",
                                             "T Cells",
                                             "B Cells",
                                             "Pericytes",
                                             "Microglia"))

total_seu$condition <- factor(as.character(total_seu$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

metadata <- total_seu@meta.data

metadata$cell_id <- rownames(metadata)

# merge in broad names
metadata <- merge(metadata,
                  cluster_xref,
                  by="merged_cell_name")

rownames(metadata) <- metadata$cell_id

metadata <- metadata[colnames(total_seu),]

total_seu$broad_cell_name <- metadata$broad_cell_name


# merge em in
total_seu <- combineExpression(combined_tcr,
                               total_seu,
                               cloneCall="gene",
                               group.by = "sample",
                               proportion = F,
                               cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# filter down to cells with clones
clone_meta <- total_seu@meta.data

clone_seu <- subset(total_seu, subset = !is.na(cloneSize))

clone_meta_f <- clone_seu@meta.data

# pseudobulk differential expression
clone_seu_pseudo <- AggregateExpression(
  clone_seu,
  assays = "RNA",
  slot = "counts",
  group.by = c("condition","orig.ident")
)

# Extract the count matrix
pseudobulk_matrix <- clone_seu_pseudo$RNA

sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")


clone_meta_f$sample <- factor(as.character(clone_meta_f$orig.ident),
                          levels=sample_levels)

# Extract unique sample-level metadata
sample_metadata <- clone_meta_f %>%
  select(sample, condition) %>%
  distinct() %>%
  as.data.frame()
rownames(sample_metadata) <- sample_metadata$sample

condition1 <- "PS19-KO"
condition2 <- "PS19-WT"

celltype_cols <- c(grep(paste0("^", condition1, "_"), 
                      colnames(pseudobulk_matrix), 
                      value = TRUE),
                   grep(paste0("^", condition2, "_"), 
                        colnames(pseudobulk_matrix), 
                        value = TRUE))


# Extract sample names from column names
sample_names <- sub(paste0("^", condition1, "_"), "", celltype_cols)
sample_names <- sub(paste0("^", condition2, "_"), "", sample_names)

# Subset counts
counts <- pseudobulk_matrix[, celltype_cols]
colnames(counts) <- sample_names

meta <- sample_metadata[sample_names, , drop = FALSE]

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
  mutate(comparison = paste0(condition1, "_m_", condition2))

# filter out NAs
res_df <- res_df[!is.na(res_df$padj),]

# save to file
write.xlsx(res_df, file=paste0(out_dir, "clone_degs.PS19-KO_m_PS19-WT.xlsx"),
           colWidths="auto")

res_df$log_p <- -log10(res_df$pvalue)


sig_results_df <- res_df[res_df$padj < 0.05 &
                               abs(res_df$log2FoldChange) > 0.5,]

logp_thresh <- min(sig_results_df$log_p)

ggplot(res_df,
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  geom_hline(yintercept = logp_thresh,
             color="red", linetype=2) +
  geom_vline(xintercept = 0.5,
             color="red", linetype=2) +
  geom_vline(xintercept = -0.5,
             color="red", linetype=2) +
  geom_point(data=sig_results_df,
             color="red") +
  geom_text_repel(data=sig_results_df,
                  aes(label=gene),
                  color="red", size=2.5,
                  max.overlaps = 50) +
  theme_bw() +
  labs(x="Log2 Fold Change", y="-log10(P-Value)", title="Clonally Expanded Cells",
       subtitle="PS19-KO - PS19-WT")
ggsave(paste0(out_dir, "clone_degs.PS19-KO_m_PS19-WT.volcano_plot.png"), width=6, height=6)

# GSEA analysis

# Load in msigdbr and get gene sets --------------------------------------------
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
kegg_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(hallmark=list_convert(hallmark_gene_sets),
                        kegg=list_convert(kegg_gene_sets),
                        reactome=list_convert(reactome_gene_sets),
                        wikipathways=list_convert(wikipath_gene_sets),
                        gobp=list_convert(gobp_gene_sets),
                        gomf=list_convert(gomf_gene_sets))

ranked_genes <- res_df %>%
  filter(!is.na(stat) &
           !is.infinite(stat)) %>%
  arrange(desc(stat)) %>%
  pull(stat, name = gene)

total_gsea_results <- lapply(names(total_gene_sets), function(gs_name) {
  
  print(gs_name)
  
  gene_sets <- total_gene_sets[[gs_name]]
  # run GSEA
  fgsea_results <- fgsea(
    pathways = gene_sets,
    stats = ranked_genes,
    minSize = 15, # Minimal size of a gene set to test
    maxSize = 500 # Maximal size of a gene set to test
  )
  
  return(fgsea_results[order(fgsea_results$pval),])
})
names(total_gsea_results) <- names(total_gene_sets)

write.xlsx(total_gsea_results, 
           file = paste0(out_dir, "clone_degs.gsea_results.PS19-KO_m_PS19-WT.xlsx"),
           colWidths="auto")

