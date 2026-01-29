library(Seurat) # single cell analysis
library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(snakecase) # better file names
library(ggplot2) # plotting
library(cowplot) # combining plots

# Initial setup ----------------------------------------------------------------

# create output directory for results
out_dir <- "results/diff_exp_fgsea.pseudobulk.subset_integrated/"
dir.create(out_dir, showWarnings = F)

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

# Read in data -----------------------------------------------------------------

# can be changes in the future to account for multiple results

results_file <- "results/diff_exp_analysis.pseudobulk.subset_integrated/ps19ko_minus_ps19wt.de_results.RDS"

results_list <- readRDS(results_file)

results_dir <- paste0(out_dir, gsub(".RDS", "", basename(results_file)), "/")

dir.create(results_dir, showWarnings = F)

# Calculate GSEA results per cell type -----------------------------------------

gsea_results <- lapply(results_list, function(degs) {
  
  # Create ranked gene list (by log fold change or stat)
  # Important: remove NAs and sort
  ranked_genes <- degs %>%
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
  
  # return GSEA results as well as ranked list
  return(list(ranked_genes=ranked_genes,
              gsea_results=total_gsea_results,
              degs=degs))
  
})

saveRDS(gsea_results, file=paste0(results_dir, "total.gsea_results.RDS"))

# output to excel files

for (cell in names(gsea_results)) {
  
  cell_results <- gsea_results[[cell]]$gsea_results
  
  write.xlsx(cell_results, 
             file=paste0(results_dir, to_snake_case(cell), ".gsea_results.xlsx"),
             colWidths="auto")
  
  
}

# find common pathways for a big ol heatmap

total_gsea_df <- lapply(names(gsea_results), function(cell) {
  
  cell_results <- gsea_results[[cell]]$gsea_results
  cell_results_df <- bind_rows(cell_results)
  
  cell_results_df$cell <- cell
  
  return(cell_results_df)
})
total_gsea_df <- bind_rows(total_gsea_df)

total_gsea_df$pathway_db <- sapply(total_gsea_df$pathway,
                                   function(x) {
                                     unlist(strsplit(x, "_"))[1]
                                   })

# sigs

sig_gsea_df <- total_gsea_df[total_gsea_df$padj < 0.05,]

table(sig_gsea_df$cell, sig_gsea_df$pathway_db)

sig_pathway_counts <- sort(table(sig_gsea_df$pathway), decreasing=T)

reactome_pathway_counts <- sig_pathway_counts[grep("REACTOME", names(sig_pathway_counts))]
kegg_pathway_counts <- sig_pathway_counts[grep("KEGG", names(sig_pathway_counts))]



