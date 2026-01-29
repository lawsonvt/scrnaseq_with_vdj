library(Seurat)
library(SeuratObject)

out_dir <- "results/seurat_cluster_naming.integrate_subsets/"
dir.create(out_dir, showWarnings = F)

# load in init names
init_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

# load in the metadata from the two subsets
mg_meta <- readRDS("results/seurat_cluster_naming.microglia_subset/cell_named.microglia_subset.cell_metadata.RDS")
tc_meta <- readRDS("results/tcell_split.doublet_find/tcell_subset.named.metadata.RDS")

# merge up the metadata

total_meta <- init_seu@meta.data

total_meta$cell_id <- rownames(total_meta)
mg_meta$cell_id <- rownames(mg_meta)
tc_meta$cell_id <- rownames(tc_meta)

# merge em up

merged_meta <- merge(total_meta,
                     mg_meta[,c("cell_id","microglia_cell_name")],
                     all.x=T,
                     by="cell_id")

merged_meta <- merge(merged_meta,
                     tc_meta[,c("cell_id", "tcell_type")],
                     all.x=T,
                     by="cell_id")

# make sure order is correct
rownames(merged_meta) <- merged_meta$cell_id

merged_meta <- merged_meta[colnames(init_seu),]

# tcell doublets
nrow(merged_meta[is.na(merged_meta$tcell_type) &
              merged_meta$cell_cluster == "T Cells",])

# drop em
merged_meta$tcell_doublet <- "singlet"
merged_meta[is.na(merged_meta$tcell_type) &
              merged_meta$cell_cluster == "T Cells",]$tcell_doublet <- "doublet"


init_seu@meta.data$tcell_doublet <- merged_meta$tcell_doublet

init_seu <- subset(init_seu, subset = tcell_doublet == "singlet")

merged_meta <- merged_meta[colnames(init_seu),]

# create a merged cell name
merged_meta$merged_cell_name <- merged_meta$cell_cluster

merged_meta[!is.na(merged_meta$tcell_type),]$merged_cell_name <-
  merged_meta[!is.na(merged_meta$tcell_type),]$tcell_type

merged_meta[!is.na(merged_meta$microglia_cell_name),]$merged_cell_name <- 
  merged_meta[!is.na(merged_meta$microglia_cell_name),]$microglia_cell_name


table(merged_meta$merged_cell_name)

init_seu@meta.data$merged_cell_name <- merged_meta$merged_cell_name

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.seurat.RDS"))






