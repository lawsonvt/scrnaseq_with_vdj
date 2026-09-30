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
ggsave(paste0(out_dir, "microglia_hallmark_top.dot_plot.png"), width=5.5, height=5)

# same for GO BP

gobp_subset <- mg_results$gsea_results$gobp[1:15,]

gobp_subset$logp <- -log10(gobp_subset$pval)
gobp_subset$pathway_pretty <- sapply(gobp_subset$pathway, function(p) {
  
  # drop DB name and replace underscores
  p <- paste0(unlist(strsplit(p, "_"))[-1], collapse=" ")
  
  # make title
  p <- str_to_title(p)
  
  # drop WP IDs
  p <- gsub("Wp[0-9]+", "", p)
  p <- trimws(p)
  
  return(p)
  
})

gobp_subset$pathway_pretty <- factor(gobp_subset$pathway_pretty,
                                         levels=rev(gobp_subset$pathway_pretty))

ggplot(gobp_subset,
       aes(x=logp, y=pathway_pretty,fill=NES)) +
  geom_point(pch=21, size=5) +
  theme_bw() +
  xlim(0, 30) + 
  scale_fill_gradient2(
    low  = "blue",
    mid  = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(x="-log10(PValue)", y=NULL)
ggsave(paste0(out_dir, "microglia_gobp_top.dot_plot.png"), width=7, height=5)



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
) + labs(title="IFN-γ Response",
         x="Rank", y="Enrichment Score") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
ggsave(paste0(out_dir, "cd8_interferon_gamma.enrichment_plot.png"), width=5, height=3)


adjp <- cd8_hallmark[cd8_hallmark$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$padj
es <- cd8_hallmark[cd8_hallmark$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$ES
nes <- cd8_hallmark[cd8_hallmark$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE",]$NES


plotEnrichment(
  pathways[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  cd8_results$ranked_genes
) + labs(title="IFN-γ Response",
         x="Rank", y="Enrichment Score") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  geom_text(x=1000, y=-0.2, label=paste0("Adj P-Value: ", formatC(adjp, format="e", digits=2)), 
            color="red", size=4, hjust=0) +
  geom_text(x=1000, y=-0.25, label=paste0("ES: ", round(es, digits=3)), color="red", size=4,
            hjust=0) +
  geom_text(x=1000, y=-0.3, label=paste0("NES: ", round(nes, digits=3)), 
            color="red", size=4, hjust = 0)
ggsave(paste0(out_dir, "cd8_interferon_gamma.enrichment_plot.annotated.png"), width=5, height=3)


