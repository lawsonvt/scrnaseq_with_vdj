library(ggplot2)
library(fgsea)
library(stringr)
library(msigdbr)  # for gene sets / pathways

out_dir <- "results/gsea_pathway_plots/"

dir.create(out_dir, showWarnings = F)

# get microglia category results
gsea_results <- readRDS("results/diff_exp_fgsea.pseudobulk.subset_integrated/ps19ko_minus_ps19wt.cell_categories.de_results/total.gsea_results.RDS")

mg_results <- gsea_results$Microglia

head(mg_results$gsea_results$hallmark)

hallmark_subset <- mg_results$gsea_results$hallmark[1:15,]

hallmark_subset$logp <- -log10(hallmark_subset$pval)
hallmark_subset$pathway_pretty <- sapply(hallmark_subset$pathway, function(p) {
  
  # drop DB name and replace underscores
  p <- paste0(unlist(strsplit(p, "_"))[-1], collapse=" ")
  
  # make title
  p <- str_to_title(p)
  
  # drop WP IDs
  p <- gsub("Wp[0-9]+", "", p)
  p <- trimws(p)
  
  return(p)
  
})

hallmark_subset$pathway_pretty <- factor(hallmark_subset$pathway_pretty,
                                         levels=rev(hallmark_subset$pathway_pretty))

ggplot(hallmark_subset,
       aes(x=logp, y=pathway_pretty,fill=NES)) +
  geom_point(pch=21, size=5) +
  theme_bw() +
  xlim(0, 20) + 
  scale_fill_gradient2(
    low  = "blue",
    mid  = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(x="-log10(PValue)", y=NULL)
ggsave(paste0(out_dir, "microglia_hallmark_top.dot_plot.png"), width=7, height=5)

# enrichment plot
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")

pathways <- hallmark_gene_sets %>%
  split(x = .$gene_symbol, f = .$gs_name)

# get CD8 results


cell_gsea_results <- readRDS("results/diff_exp_fgsea.pseudobulk.subset_integrated/ps19ko_minus_ps19wt.cell_types.de_results/total.gsea_results.RDS")

cd8_results <- cell_gsea_results$CD8p

cd8_hallmark <- cd8_results$gsea_results$hallmark

head(cd8_hallmark[,1:7])

plotEnrichment(
  pathways[["HALLMARK_KRAS_SIGNALING_DN"]],
  cd8_results$ranked_genes
) + labs(title="KRAS Signaling Down")
ggsave(paste0(out_dir, "cd8_kras_signaling.enrichment_plot.png"), width=5, height=3)


plotEnrichment(
  pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  cd8_results$ranked_genes
) + labs(title="IFN-γ Response")
ggsave(paste0(out_dir, "cd8_interferon_gamma.enrichment_plot.png"), width=5, height=3)


# merge up GSEA results
cd8_results_df <- bind_rows(cd8_results$gsea_results)


