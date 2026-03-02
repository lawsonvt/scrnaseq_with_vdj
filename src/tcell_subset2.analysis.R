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

out_dir <- "results/tcell_subset2.analysis/"
dir.create(out_dir, showWarnings = F)

# load in seurat data
tcell_seu <- readRDS("results/tcell_split.tcell_subset2/tcell_subset2_named.seurat.RDS")

# get metadata
tcell_meta <- tcell_seu@meta.data

# make heat colors for clonal counts
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

ggplot(tcell_meta, 
       aes(x=condition,
           fill=cloneSize)) +
  geom_bar(color="black") +
  facet_wrap(~ tcell_subset_clusters_name, nrow=2) +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(y = "Cell Count", x=NULL, fill="Clone Size")


ggplot(tcell_meta[tcell_meta$condition %in% c("PS19-WT",
                                              "PS19-KO"),], 
       aes(x=condition,
           fill=cloneSize)) +
  geom_bar(color="black") +
  facet_wrap(~ tcell_subset_clusters_name, nrow=2) +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(y = "Cell Count", x=NULL, fill="Clone Size")

# load in tcr data
combined_tcr <- readRDS("results/tcr_analysis.subset_integrated/combined_tcr.RDS")
tcell_meta$cell_id <- rownames(tcell_meta)

tcell_tcr <- lapply(combined_tcr, function(data) {
  
  merge(data, tcell_meta[,c("cell_id","tcell_subset_clusters_name")],
        by.x="barcode", by.y="cell_id")
  
})

# does not seem to work for plotting
sample_levels <- c("B-Q617-WT",
                   "G-W138-WT",
                   "C-Q635-KO",
                   "F-W137-KO",
                   "D-Q619-PS19-WT",
                   "E-W136-PS19-WT",
                   "A-Q637-PS19-KO",
                   "H-E806-PS19-KO")

condition_levels <- c("WT","KO","PS19-WT","PS19-KO")
# 
# 
# # factorize the samples and conditions
# tcell_tcr <- lapply(tcell_tcr, function(data) {
#   
#   data$sample <- factor(data$sample,
#                         levels=sample_levels)
#   data$condition <- factor(data$condition,
#                            levels=condition_levels)
#   
#   return(data)
#   
# })

cd_data_s <- clonalDiversity(tcell_tcr,
                cloneCall = "strict",
                x.axis="condition",
                metric="shannon",
                exportTable = T) 

cd_data_g <- clonalDiversity(tcell_tcr,
                cloneCall = "strict",
                x.axis="condition",
                metric="gini",
                exportTable = T) 

cd_data <- rbind(cd_data_s,
                 cd_data_g)

cd_data$Group <- factor(cd_data$Group,
                        levels=sample_levels)
cd_data$condition <- factor(cd_data$condition,
                            levels=condition_levels)

ggplot(cd_data,
       aes(x=condition,
           y=value)) +
  geom_boxplot() +
  geom_jitter(aes(color=Group)) +
  












clonalDiversity(tcell_tcr,
                cloneCall = )

percentGeneUsage(tcell_tcr,
                 genes = "TRBV",
                 group.by = "condition",
                 plot.type = "heatmap",
                 summary.fun = "percent")

percentGeneUsage(tcell_tcr,
                 genes = "TRAV",
                 group.by = "condition",
                 plot.type = "heatmap",
                 summary.fun = "percent")

