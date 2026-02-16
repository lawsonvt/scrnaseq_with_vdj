library(scRepertoire)
library(Seurat)
library(ggplot2)
library(cowplot)

# load in data
out_dir <- "results/tcr_analysis.subset_integrated/"

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


clonalLength(subsetClones(combined_tcr,
                          name = "condition",
                          variables=c("PS19-WT","PS19-KO")),
             cloneCall = "nt",
             group.by = "condition",
             order.by = c("PS19-WT", "PS19-KO",
                          "WT", "KO"),
             scale=T,
             chain="both")
ggsave(paste0(out_dir, "clonal_length.ps19_density_plot.png"), width=6, height=4)


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
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

total_meta <- total_seu@meta.data


DimPlot(total_seu, group.by = c("broad_cell_name","cloneSize"), reduction = "umap.harmony") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clonal_size.seurat_umap.png"), width=12, height=4)

# meta with contigs
total_meta_contigs <- total_meta[!is.na(total_meta$cloneSize),]




# subset down to T-Cells
tcell_seu <- subset(total_seu, subset = merged_cell_name == "CD4+" | merged_cell_name == "CD8+")

# filter out weird outliers
umap_coords <- Embeddings(tcell_seu, reduction = "umap.harmony")

keep_cells <- rownames(umap_coords)[
  umap_coords[,1] < 2.5 & umap_coords[,1] > -5 &
    umap_coords[,2] < -10
]

tcell_seu <- subset(tcell_seu, cells = keep_cells)

DimPlot(tcell_seu, group.by = c("merged_cell_name",
                                "cloneSize"), reduction = "umap.harmony") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) 
ggsave(paste0(out_dir, "tcell_clonal_size.seurat_umap.png"), width=12, height=4)


DimPlot(tcell_seu, group.by = "cloneSize", 
        split.by = "condition", ncol=2,
        reduction = "umap.harmony",
        raster=F) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)])) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL)
ggsave(paste0(out_dir, "tcell_clonal_umaps.per_condition.png"), 
       width=8, height=5)

tcell_metadata <- tcell_seu@meta.data

table(tcell_metadata$condition,
      tcell_metadata$cloneSize)

# subset by cell type
cd8_meta <- tcell_metadata[tcell_metadata$merged_cell_name == "CD8+",]
cd4_meta <- tcell_metadata[tcell_metadata$merged_cell_name == "CD4+",]

cd8_clone_size <- as.data.frame(table(cd8_meta$condition,
                        cd8_meta$cloneSize))
cd4_clone_size <- as.data.frame(table(cd4_meta$condition,
                                      cd4_meta$cloneSize))

cd8_total <- as.data.frame(table(cd8_meta$condition))
cd4_total <- as.data.frame(table(cd4_meta$condition))

cd8_clone_size <- merge(cd8_clone_size,
                        cd8_total, 
                        by="Var1",
                        suffixes=c(".clone",".total"))

cd8_clone_size$frac <- cd8_clone_size$Freq.clone / cd8_clone_size$Freq.total * 100
cd8_clone_size$type <- "CD8+"

cd4_clone_size <- merge(cd4_clone_size,
                        cd4_total, 
                        by="Var1",
                        suffixes=c(".clone",".total"))

cd4_clone_size$frac <- cd4_clone_size$Freq.clone / cd4_clone_size$Freq.total * 100
cd4_clone_size$type <- "CD4+"

total_clone_size <- rbind(cd8_clone_size,
                          cd4_clone_size)

ggplot(total_clone_size,
       aes(x=Var1,
           y=frac,
           fill=Var2)) +
  facet_wrap(~ type, ncol=2) +
  geom_bar(stat="identity", color="black") +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(3,4,5,7)])) +
  labs(x=NULL, y="Percentage of Celltype", fill=NULL)
ggsave(paste0(out_dir, "clonal_size_per_condition.bar_plots.png"), width=8, height=5)

ggplot(total_clone_size[total_clone_size$Var1 
                        %in% c("PS19-WT","PS19-KO"),],
       aes(x=Var1,
           y=frac,
           fill=Var2)) +
  facet_wrap(~ type, ncol=2) +
  geom_bar(stat="identity", color="black") +
  theme_bw() +
  scale_fill_manual(values=rev(colorblind_vector[c(3,4,5,7)])) +
  labs(x=NULL, y="Percentage of Celltype", fill=NULL)
ggsave(paste0(out_dir, "clonal_size_per_condition.ps19.bar_plots.png"), width=6, height=5)

table(cd4_meta$condition,
      cd4_meta$cloneSize)

