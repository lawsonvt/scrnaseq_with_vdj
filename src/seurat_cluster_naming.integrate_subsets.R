library(Seurat)
library(SeuratObject)

out_dir <- "results/seurat_cluster_naming.integrate_subsets/"
dir.create(out_dir, showWarnings = F)

# load in init names
init_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

# load in the metadata from the other subsets
mg_meta <- readRDS("results/seurat_cluster_naming.microglia_subset/cell_named.microglia_subset.cell_metadata.RDS")

peri_meta <- readRDS("results/pericytes_split/pericytes_subset.named.metadata.RDS")

peri_meta$pericyte_type <- gsub("^NK$", "Natural Killer", peri_meta$pericyte_type)
peri_meta$pericyte_type <- gsub("Unknown", "P Unknown", peri_meta$pericyte_type)

tc_meta <- readRDS("results/tcell_split.post_pericytes_split/tcell_cell_named.metadata.RDS")

tc_meta$tcell_sub_cell_name <- gsub("Unknown2", "T Unknown", tc_meta$tcell_sub_cell_name)

# load in doublet metadata
tc_doublets <- readRDS("results/tcell_split.post_pericytes_split/tcell.doublet_output.RDS")

tc_doublet_ids <- rownames(tc_doublets[tc_doublets$scDblFinder.class == "doublet",])

# merge up the metadata

total_meta <- init_seu@meta.data

total_meta$cell_id <- rownames(total_meta)
mg_meta$cell_id <- rownames(mg_meta)
tc_meta$cell_id <- rownames(tc_meta)
peri_meta$cell_id <- rownames(peri_meta)

# merge em up

merged_meta <- merge(total_meta,
                     mg_meta[,c("cell_id","microglia_cell_name")],
                     all.x=T,
                     by="cell_id")

merged_meta <- merge(merged_meta,
                     peri_meta[,c("cell_id", "pericyte_type")],
                     all.x=T,
                     by="cell_id")

merged_meta <- merge(merged_meta,
                     tc_meta[,c("cell_id", "tcell_sub_cell_name")],
                     all.x=T,
                     by="cell_id")



# make sure order is correct
rownames(merged_meta) <- merged_meta$cell_id

merged_meta <- merged_meta[colnames(init_seu),]


# drop doublets

merged_meta <- merged_meta[!merged_meta$cell_id %in% tc_doublet_ids,]

# create a merged cell name
merged_meta$merged_cell_name <- merged_meta$cell_cluster

merged_meta[!is.na(merged_meta$microglia_cell_name),]$merged_cell_name <- 
  merged_meta[!is.na(merged_meta$microglia_cell_name),]$microglia_cell_name

merged_meta[!is.na(merged_meta$pericyte_type),]$merged_cell_name <- 
  merged_meta[!is.na(merged_meta$pericyte_type),]$pericyte_type

merged_meta[!is.na(merged_meta$tcell_sub_cell_name),]$merged_cell_name <-
  merged_meta[!is.na(merged_meta$tcell_sub_cell_name),]$tcell_sub_cell_name


table(merged_meta$merged_cell_name)

# filter down seurat object
init_seu <- subset(init_seu, cells = merged_meta$cell_id)

merged_meta <- merged_meta[colnames(init_seu),]

init_seu@meta.data$merged_cell_name <- merged_meta$merged_cell_name

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.full.seurat.RDS"))

# remove unknown cells


init_seu <- subset(init_seu, subset = merged_cell_name %in% c("P Unknown",
                                                              "T Unknown",
                                                              "MG Unknown"),
                   invert=T)


SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.no_unknown.seurat.RDS"))


# remove low counts
min_count <- 100

cell_counts <- table(init_seu$merged_cell_name)

names(cell_counts)[cell_counts >= min_count]

# filter out unknown and platelets
init_seu <- subset(init_seu, subset = merged_cell_name %in% names(cell_counts)[cell_counts >= min_count])

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.no_unknown_no_low_count.seurat.RDS"))







