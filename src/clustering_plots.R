library(Seurat)
library(SeuratObject)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(scRepertoire)

out_dir <- "results/clustering_plots/"
dir.create(out_dir, showWarnings = F)

# load in seurat objects
total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.no_low_count.seurat.RDS")

# fix name
#total_seu$merged_cell_name <- gsub("^Microglia$", "Intermediate Microglia", total_seu$merged_cell_name)
total_seu$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", total_seu$orig.ident)

total_seu$condition <- factor(as.character(total_seu$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))


# # broad cell type name
# cluster_xref <- data.frame(merged_cell_name=c("Intermediate Microglia",
#                              "Cross Presenting",
#                              "Homeostatic",
#                              "Macrophages",
#                              "Monocytes",
#                              "Interferon Responsive",
#                              "Neutrophils",
#                              "DAM",
#                              "Proliferating Cells",
#                              "CD4+",
#                              "CD8+",
#                              "NK",
#                              "B Cells",
#                              "Pericytes",
#                              "Inflammatory"),
#                            broad_cell_name=c("Microglia",
#                              "Microglia",
#                              "Microglia",
#                              "Macrophages",
#                              "Monocytes",
#                              "Microglia",
#                              "Neutrophils",
#                              "Microglia",
#                              "Proliferating Cells",
#                              "T Cells",
#                              "T Cells",
#                              "T Cells",
#                              "B Cells",
#                              "Pericytes",
#                              "Microglia"))

metadata <- total_seu@meta.data

metadata$cell_id <- rownames(metadata)

merged_cell_colors <- pal_d3("category20")(length(unique(total_seu$merged_cell_name)))


DimPlot(total_seu, reduction="umap.harmony",
        group.by= "cell_category") +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "total_broad_clusters.all_samples_umap.png"),
       width=6, height=4)

DimPlot(total_seu, reduction="umap.harmony", 
        group.by= "merged_cell_name", 
        cols=merged_cell_colors) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "subset_clusters.all_samples_umap.png"),
       width=6, height=4)

# make a dot plot for the broad categories

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"))

lab_cluster_markers <- c("Cd79a","Pf4","P2ry12","Ccr2","Ncr1","Mmp8","Mki67","Cd3d")

total_seu$cell_category <- factor(total_seu$cell_category,
                                 levels=rev(sort(unique(total_seu$cell_category))))

Idents(total_seu) <- "cell_category"


DotPlot(total_seu,
        features=lab_cluster_markers) +
  labs(x=NULL, y=NULL) + RotatedAxis()
ggsave(paste0(out_dir, "cell_category.marker_dot_plot.png"), width=6, height=3,
       bg="white")




DimPlot(total_seu, reduction="umap.harmony", 
        group.by= "merged_cell_name", 
        split.by = "condition",
        raster = T, ncol = 2,
        cols=merged_cell_colors, raster.dpi = c(256,256)) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "cell_types_per_condition.all_samples_umap.png"),
       width=10, height=8)

# microglia subset counts

metadata_mg <- metadata[metadata$cell_category == "Microglia",]

microglia_subset_counts <- as.data.frame(table(metadata_mg$condition,
                                        metadata_mg$merged_cell_name))

microglia_counts <- as.data.frame(table(metadata_mg$condition))

microglia_subset_counts <- merge(microglia_subset_counts,
                                 microglia_counts,
                                 by="Var1")
microglia_subset_counts$frac <- microglia_subset_counts$Freq.x / 
  microglia_subset_counts$Freq.y * 100

microglia_subset_counts <- microglia_subset_counts[order(microglia_subset_counts$frac,
                                                         decreasing = T),]

microglia_subset_counts$Var2 <- factor(microglia_subset_counts$Var2,
                                       levels=unique(microglia_subset_counts$Var2))

ggplot(microglia_subset_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(y=NULL, x="Percentage of Microglia Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_microglia_cells.bar_plot.png"),
       width=5, height=3)

ggplot(microglia_subset_counts,
       aes(x=Var2, y=Freq.x, fill=Var1)) +
  geom_bar(stat="identity", color="black", position="dodge") +
  #scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(y="Cell Count", x=NULL, fill=NULL)
ggsave(paste0(out_dir, "microglia_cell_counts.bar_plot.png"),
       width=5, height=3)


broad_cell_counts <- as.data.frame(table(metadata$condition,
                                         metadata$cell_category))
total_counts <- as.data.frame(table(metadata$condition))

broad_cell_counts <- merge(broad_cell_counts,
                           total_counts,
                           by="Var1")
broad_cell_counts$frac <- broad_cell_counts$Freq.x / 
  broad_cell_counts$Freq.y * 100

broad_cell_counts <- broad_cell_counts[order(broad_cell_counts$frac,
                                             decreasing=T),]

broad_cell_counts$Var2 <- factor(as.character(broad_cell_counts$Var2),
                                 levels=unique(broad_cell_counts$Var2))

ggplot(broad_cell_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y=NULL, x="Percentage of Total Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_total_cells.bar_plot.png"),
       width=5, height=3)

ggplot(broad_cell_counts,
       aes(x=Var2, y=Freq.x, fill=Var1)) +
  geom_bar(stat="identity", color="black", position="dodge") +
  #scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(y="Cell Count", x=NULL, fill=NULL)
ggsave(paste0(out_dir, "total_cell_counts.bar_plot.png"),
       width=5, height=3)

# same percentage, this time without microglia
broad_cell_counts <- as.data.frame(table(metadata[metadata$cell_category != "Microglia",]$condition,
                                         metadata[metadata$cell_category != "Microglia",]$cell_category))
total_counts <- as.data.frame(table(metadata[metadata$cell_category != "Microglia",]$condition))

broad_cell_counts <- merge(broad_cell_counts,
                           total_counts,
                           by="Var1")
broad_cell_counts$frac <- broad_cell_counts$Freq.x / 
  broad_cell_counts$Freq.y * 100

broad_cell_counts <- broad_cell_counts[order(broad_cell_counts$frac,
                                             decreasing=T),]

broad_cell_counts$Var2 <- factor(as.character(broad_cell_counts$Var2),
                                 levels=unique(broad_cell_counts$Var2))

ggplot(broad_cell_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y=NULL, x="Percentage of Non-Microglia Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_non_microglia_cells.bar_plot.png"),
       width=5, height=3)

ggplot(broad_cell_counts,
       aes(x=Var2, y=Freq.x, fill=Var1)) +
  geom_bar(stat="identity", color="black", position="dodge") +
  #scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(y="Cell Count", x=NULL, fill=NULL)
ggsave(paste0(out_dir, "non_microglia_cell_counts.bar_plot.png"),
       width=5, height=3)

# tcell breakdown
t_cell_counts <- as.data.frame(table(metadata[metadata$cell_category == "T Cells",]$condition,
                                     metadata[metadata$cell_category == "T Cells",]$merged_cell_name))
total_counts <- as.data.frame(table(metadata[metadata$cell_category == "T Cells",]$condition))

t_cell_counts <- merge(t_cell_counts,
                           total_counts,
                           by="Var1")
t_cell_counts$frac <- t_cell_counts$Freq.x / 
  t_cell_counts$Freq.y * 100

t_cell_counts <- t_cell_counts[order(t_cell_counts$frac,
                                             decreasing=T),]

t_cell_counts$Var2 <- factor(as.character(t_cell_counts$Var2),
                                 levels=unique(t_cell_counts$Var2))


ggplot(t_cell_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y=NULL, x="Percentage of T Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_t_cells.bar_plot.png"),
       width=5, height=3)

ggplot(t_cell_counts,
       aes(x=Var2, y=Freq.x, fill=Var1)) +
  geom_bar(stat="identity", color="black", position="dodge") +
  #scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(y="Cell Count", x=NULL, fill=NULL)
ggsave(paste0(out_dir, "t_cell_counts.bar_plot.png"),
       width=5, height=3)

# filter down to tcells
tcell_seu <- subset(total_seu, subset = cell_category == "T Cells")

# redo umap
tcell_seu <- RunPCA(tcell_seu, npcs = 20)

max_pc_dim <- 10

tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

DimPlot(tcell_seu, reduction="umap.tcell_pca",
        group.by= "merged_cell_name") +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "tcell_umap.png"),
       width=6, height=4)

DimPlot(tcell_seu, reduction="umap.tcell_pca",
        group.by= "merged_cell_name",
        split.by = "condition", ncol=2) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "tcell_umap.per_condition.png"),
       width=8, height=6)


# integrate clonal data

# merge in TCR clonal data
combined_tcr <- readRDS("results/tcr_analysis.subset_integrated/combined_tcr.RDS")

tcell_seu <- combineExpression(combined_tcr,
                             tcell_seu,
                             cloneCall="gene",
                             group.by = "sample",
                             proportion = F,
                             cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

tcell_meta <- tcell_seu@meta.data

ggplot(tcell_meta, 
       aes(x=condition,
           fill=cloneSize)) +
  geom_bar(color="black") +
  facet_wrap(~ merged_cell_name, nrow=1) +
  scale_fill_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  labs(y = "Cell Count", x=NULL, fill="Clone Size")
ggsave(paste0(out_dir, "tcell_clonal_counts.png"), width=7, height=4)

DimPlot(tcell_seu, reduction="umap.tcell_pca",
        group.by= "cloneSize",
        split.by = "condition", ncol=2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "tcell_umap.clone_size.per_condition.png"),
       width=8, height=6)



