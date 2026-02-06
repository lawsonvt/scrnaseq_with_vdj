library(Seurat)
library(SeuratObject)
library(ggplot2)
library(RColorBrewer)
library(ggsci)

out_dir <- "results/clustering_plots/"
dir.create(out_dir, showWarnings = F)

# load in seurat objects
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

merged_cell_colors <- pal_d3("category20")(length(unique(total_seu$merged_cell_name)))


DimPlot(total_seu, reduction="umap.harmony",
        group.by= "broad_cell_name") +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "total_broad_clusters.all_samples_umap.png"),
       width=8, height=6)

DimPlot(total_seu, reduction="umap.harmony", 
        group.by= "merged_cell_name", 
        cols=merged_cell_colors) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "subset_clusters.all_samples_umap.png"),
       width=8, height=6)

# make a dot plot for the broad categories

lab_cluster_markers <- list("B Cells"=c("Cd19", "Cd79a", "Ms4a1"),
                            "Macrophages"=c("Pf4", "Mrc1", "Ms4a7"),
                            "Microglia"=c("Sall1", "Hexb", "P2ry12"),
                            "Monocytes"=c("Ccr2", "Cd44"),
                            "Neutrophils"=c("Ly6g", "Ngp", "Mmp8"),
                            "Pericytes"=c("Pdgfrb", "Rgs5"),
                            "Proliferating Cells"=c("Mki67", "Ccnb1", "Tpx2"),
                            "T Cells"=c("Trbc2", "Cd3d", "Lck"))

lab_cluster_markers <- c("Cd79a","Pf4","Hexb","Ccr2","Mmp8","Rgs5","Mki67","Cd3d")

total_seu$broad_cell_name <- factor(total_seu$broad_cell_name,
                                 levels=rev(sort(unique(total_seu$broad_cell_name))))

Idents(total_seu) <- "broad_cell_name"


DotPlot(total_seu,
        features=lab_cluster_markers) +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "broad_cluster_names.marker_dot_plot.png"), width=8, height=4,
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

metadata_mg <- metadata[metadata$broad_cell_name == "Microglia",]

microglia_subset_counts <- as.data.frame(table(metadata_mg$condition,
                                        metadata_mg$merged_cell_name))

microglia_counts <- as.data.frame(table(metadata_mg$condition))

microglia_subset_counts <- merge(microglia_counts,
                                 microglia_subset_counts,
                                 by="Var1")
microglia_subset_counts$frac <- microglia_subset_counts$Freq.y / 
  microglia_subset_counts$Freq.x * 100

ggplot(microglia_subset_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(y=NULL, x="Percentage of Microglia Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_microglia_cells.bar_plot.png"),
       width=7, height=4)

broad_cell_counts <- as.data.frame(table(metadata$condition,
                                         metadata$broad_cell_name))
total_counts <- as.data.frame(table(metadata$condition))

broad_cell_counts <- merge(broad_cell_counts,
                           total_counts,
                           by="Var1")
broad_cell_counts$frac <- broad_cell_counts$Freq.x / 
  broad_cell_counts$Freq.y * 100

ggplot(broad_cell_counts,
       aes(x=frac, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", color="black") +
  #scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(y=NULL, x="Percentage of Total Cells", fill=NULL)
ggsave(paste0(out_dir, "percentage_of_total_cells.bar_plot.png"),
       width=7, height=4)



cd_cells <- subset(total_seu, subset = merged_cell_name == "CD4+" | merged_cell_name == "CD8+")
# subset down to PS19 cells
# ps19_seu <- subset(total_seu, subset = condition == "PS19-WT" |
#                      condition == "PS19-KO")

# filter out weird outliers
umap_coords <- Embeddings(cd_cells, reduction = "umap.harmony")

keep_cells <- rownames(umap_coords)[
  umap_coords[,1] < 2.5 & umap_coords[,1] > -5 &
    umap_coords[,2] < -10
]

cd_cells <- subset(cd_cells, cells = keep_cells)

DimPlot(cd_cells, reduction="umap.harmony", 
        group.by= "merged_cell_name") +
  labs(x="UMAP 1", y="UMAP 2", )



