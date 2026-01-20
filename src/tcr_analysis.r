library(scRepertoire)
library(Seurat)
library(ggplot2)
library(cowplot)

# load in data
out_dir <- "results/tcr_analysis/"

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


clonalQuant(combined_tcr,
            cloneCall="strict", 
            group.by = "condition",
            chain = "both", 
            scale = TRUE) +
  theme(legend.position = "none") +
  labs(x=NULL) + ylim(0,100)
ggsave(paste0(out_dir, "clonal_quantification.per_condition.png"), width=6, height=4)

clonalAbundance(combined_tcr, 
                cloneCall = "gene", 
                scale = FALSE)


clonalLength(combined_tcr,
             cloneCall = "nt",
             group.by = "condition",
             order.by = c("PS19-WT", "PS19-KO",
                          "WT", "KO"),
             chain="both")


p1 <- clonalLength(subsetClones(combined_tcr,
                          name = "condition",
                          variables=c("WT","KO")),
             cloneCall = "nt",
             group.by = "condition",
             order.by = c("PS19-WT", "PS19-KO",
                          "WT", "KO"),
             chain="both") + scale_x_discrete(breaks = round(seq(0, 200, by = 5), 10)) +
  theme(legend.position = "bottom")

p2 <- clonalLength(subsetClones(combined_tcr,
                          name = "condition",
                          variables=c("PS19-WT","PS19-KO")),
             cloneCall = "nt",
             group.by = "condition",
             order.by = c("PS19-WT", "PS19-KO",
                          "WT", "KO"),
             chain="both") + scale_x_discrete(breaks = round(seq(0, 200, by = 5), 10)) +
  theme(legend.position = "bottom")

plot_grid(p1, p2, ncol=1)
ggsave(paste0(out_dir, "clonal_length.condition_split.png"), width=8, height=6)

# combine with Seurat objects

total_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

test_barcodes <- combined_tcr[[1]]$barcode

# merge em in
total_seu <- combineExpression(combined_tcr,
                             total_seu,
                             cloneCall="gene",
                             group.by = "sample",
                             proportion = F,
                             cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

DimPlot(total_seu, group.by = c("cell_cluster","cloneSize"), reduction = "umap.harmony") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clonal_size.seurat_umap.png"), width=12, height=4)

# DimPlot(total_seu, group.by = "cell_cluster", reduction = "umap.harmony",
#         label=T) + theme(legend.position = "none") + labs(x="UMAP 1", y="UMAP 2")
# 
# DimPlot(total_seu, group.by = "cloneSize", reduction = "umap.harmony") +
#   scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
#   labs()

# same thing, but with the T cell subset

tcell_seu <- LoadSeuratRds("results/tcell_split/tcell_subset.named.seurat.RDS")

tcell_seu <- combineExpression(combined_tcr,
                               tcell_seu,
                               cloneCall="gene",
                               group.by = "sample",
                               proportion = F,
                               cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

DimPlot(tcell_seu, group.by = c("tcell_type","cloneSize"), reduction = "umap.tcell_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) 
ggsave(paste0(out_dir, "clonal_size.tcell_subset.seurat_umap.png"), width=12, height=4)


DimPlot(tcell_seu, group.by = "cloneSize", 
        split.by="condition", reduction = "umap.tcell_pca", ncol=2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x="UMAP 1", y="UMAP 2")
ggsave(paste0(out_dir, "clonal_size.tcell_subset.per_condition.seurat_umap.png"), width=7, height=5)

# look at the metadata

total_metadata <- total_seu@meta.data
tcell_metadata <- tcell_seu@meta.data

# filter down to cells with clones
total_metadata_c <- total_metadata[!is.na(total_metadata$CTnt),]
tcell_metadata_c <- tcell_metadata[!is.na(tcell_metadata$CTnt),]

# cell counts
table(total_metadata_c$cell_cluster)
table(tcell_metadata_c$tcell_type)

# breakdown the counts
total_counts <- as.data.frame(table(total_metadata_c$cell_cluster))
total_counts <- total_counts[order(total_counts$Freq, decreasing=T),]
total_counts$Var1 <- factor(as.character(total_counts$Var1),
                            levels=as.character(total_counts$Var1))


ggplot(total_counts,
       aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black", fill="grey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Clone Count")
ggsave(paste0(out_dir, "clone_counts.per_total_cell.png"), width=5, height=4)


tcell_counts <- as.data.frame(table(tcell_metadata_c$tcell_type,
                                    tcell_metadata_c$condition))
tcell_counts <- tcell_counts[order(tcell_counts$Freq, decreasing=T),]
tcell_counts$Var1 <- factor(as.character(tcell_counts$Var1),
                            levels=unique(as.character(tcell_counts$Var1)))


ggplot(tcell_counts,
       aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Clone Count", fill=NULL)
ggsave(paste0(out_dir, "clone_counts.per_tcell.png"), width=4, height=4)


