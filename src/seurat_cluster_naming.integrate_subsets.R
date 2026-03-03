library(Seurat)
library(SeuratObject)

out_dir <- "results/seurat_cluster_naming.integrate_subsets/"
dir.create(out_dir, showWarnings = F)

# load in init names
init_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

# load in the metadata from the other subsets
mg_meta <- readRDS("results/seurat_cluster_naming.microglia_subset/cell_named.microglia_subset.cell_metadata.RDS")

mg_meta$microglia_cell_name <- gsub("^Microglia$", "Intermediate Microglia", mg_meta$microglia_cell_name)

peri_meta <- readRDS("results/pericytes_split/pericytes_subset.named.metadata.RDS")

peri_meta$pericyte_type <- gsub("^NK$", "Natural Killer", peri_meta$pericyte_type)
peri_meta$pericyte_type <- gsub("Unknown", "P Unknown", peri_meta$pericyte_type)

tc_meta <- readRDS("results/tcell_split.post_pericytes_split/tcell_cell_named.metadata.RDS")

tc_meta$tcell_sub_cell_name <- gsub("Unknown2", "T Unknown", tc_meta$tcell_sub_cell_name)

tc_meta_sub <- readRDS("results/tcell_split.tcell_subset_microglia_filter/tcell.microglia_filted_and_named.metadata.RDS")

# load in doublet metadata
tc_doublets <- readRDS("results/tcell_split.post_pericytes_split/tcell.doublet_output.RDS")

tc_doublet_ids <- rownames(tc_doublets[tc_doublets$scDblFinder.class == "doublet",])

# merge up the metadata

total_meta <- init_seu@meta.data

total_meta$cell_id <- rownames(total_meta)
mg_meta$cell_id <- rownames(mg_meta)
tc_meta$cell_id <- rownames(tc_meta)
tc_meta_sub$cell_id <- rownames(tc_meta_sub)
peri_meta$cell_id <- rownames(peri_meta)

# merge t cell data

unique(tc_meta_sub$tcell_clusters_mg_filtered_type)

tc_meta_sub$tcell_clusters_mg_filtered_type <- gsub("CD8","CD8+",tc_meta_sub$tcell_clusters_mg_filtered_type)
tc_meta_sub$tcell_clusters_mg_filtered_type <- gsub("CD4","CD4+",tc_meta_sub$tcell_clusters_mg_filtered_type)
tc_meta_sub$tcell_clusters_mg_filtered_type <- gsub("TCELL","Naive T Cells",tc_meta_sub$tcell_clusters_mg_filtered_type)
tc_meta_sub$tcell_clusters_mg_filtered_type <- gsub("PC","Proliferating CD8+",tc_meta_sub$tcell_clusters_mg_filtered_type)


tcell_init_cells <- c("Naive T Cells",
                      "CD8+",
                      "CD4+",
                      "Proliferating CD8+",
                      "Microglia Assoc.")

merged_tc_meta <- merge(tc_meta[,c("cell_id","tcell_sub_cell_name")],
                        tc_meta_sub[,c("cell_id","tcell_clusters_mg_filtered_name","tcell_clusters_mg_filtered_type")],
                        all.x=T, by="cell_id")

merged_tc_meta[is.na(merged_tc_meta$tcell_clusters_mg_filtered_name) &
                 merged_tc_meta$tcell_sub_cell_name %in% tcell_init_cells,]$tcell_sub_cell_name <- "T Unknown"

merged_tc_meta[!is.na(merged_tc_meta$tcell_clusters_mg_filtered_type),]$tcell_sub_cell_name <- 
  merged_tc_meta[!is.na(merged_tc_meta$tcell_clusters_mg_filtered_type),]$tcell_clusters_mg_filtered_type

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
                     merged_tc_meta,
                     all.x=T,
                     by="cell_id")

# remove NAs (to make filtering down easier)
merged_meta[is.na(merged_meta$microglia_cell_name),]$microglia_cell_name <- ""
merged_meta[is.na(merged_meta$pericyte_type),]$pericyte_type <- ""
merged_meta[is.na(merged_meta$tcell_sub_cell_name),]$tcell_sub_cell_name <- ""
merged_meta[is.na(merged_meta$tcell_clusters_mg_filtered_name),]$tcell_clusters_mg_filtered_name <- ""
merged_meta[is.na(merged_meta$tcell_clusters_mg_filtered_type),]$tcell_clusters_mg_filtered_type <- ""

# filter out unknowns and other dropped cells

# unknown microglia cells
merged_meta <- merged_meta[merged_meta$microglia_cell_name != "MG Unknown",]

# unknown pericytes
merged_meta <- merged_meta[merged_meta$pericyte_type != "P Unknown",]

# drop filtered out T cells
merged_meta <- merged_meta[merged_meta$tcell_sub_cell_name != "T Unknown",]


# drop doublets

merged_meta <- merged_meta[!merged_meta$cell_id %in% tc_doublet_ids,]

# create a merged cell name
merged_meta$merged_cell_name <- merged_meta$cell_cluster

merged_meta[merged_meta$microglia_cell_name != "",]$merged_cell_name <- 
  merged_meta[merged_meta$microglia_cell_name != "",]$microglia_cell_name

merged_meta[merged_meta$pericyte_type != "",]$merged_cell_name <- 
  merged_meta[merged_meta$pericyte_type != "",]$pericyte_type

merged_meta[merged_meta$tcell_sub_cell_name != "",]$merged_cell_name <-
  merged_meta[merged_meta$tcell_sub_cell_name != "",]$tcell_sub_cell_name


table(merged_meta$merged_cell_name)

# make a cross reference for the broad category
cat_xref <- data.frame(merged_cell_name=c("B Cells",
                                          "CD4+",
                                            "CD8+",
                                            "DAM",
                                          "Homeostatic",
                                          "Interferon Responsive",
                                          "Intermediate Microglia",
                                          "Macrophages",
                                          "Monocytes",
                                          "Naive T Cells",
                                          "Natural Killer",
                                          "Neutrophils",
                                          "Pericytes",
                                          "Platelets",
                                          "Proliferating CD8+",
                                            "Proliferating Cells"),
                       cell_category=c("B Cells",
                                       "T Cells",
                                       "T Cells",
                                       "Microglia",
                                       "Microglia",
                                       "Microglia",
                                       "Microglia",
                                       "Macrophages",
                                       "Monocytes",
                                       "T Cells",
                                       "Natural Killer",
                                       "Neutrophils",
                                       "Pericytes",
                                       "Platelets",
                                       "T Cells",
                                       "Proliferating Cells"))

# merge it in
merged_meta <- merge(merged_meta,
                     cat_xref,
                     by="merged_cell_name")

rownames(merged_meta) <- merged_meta$cell_id

# filter down seurat object
init_seu <- subset(init_seu, cells = merged_meta$cell_id)

merged_meta <- merged_meta[colnames(init_seu),]

init_seu@meta.data$merged_cell_name <- merged_meta$merged_cell_name
init_seu@meta.data$cell_category <- merged_meta$cell_category

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.full.seurat.RDS"))

# # remove unknown cells
# 
# 
# init_seu <- subset(init_seu, subset = merged_cell_name %in% c("P Unknown",
#                                                               "T Unknown",
#                                                               "MG Unknown"),
#                    invert=T)
# 
# 
# SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.no_unknown.seurat.RDS"))


# remove low counts
min_count <- 100

cell_counts <- table(init_seu$merged_cell_name)

names(cell_counts)[cell_counts >= min_count]

# filter out unknown and platelets
init_seu <- subset(init_seu, subset = merged_cell_name %in% names(cell_counts)[cell_counts >= min_count])

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.no_low_count.seurat.RDS"))


# filter down to just microglia

init_seu <- subset(init_seu, subset = cell_category == "Microglia")

SaveSeuratRds(init_seu, file=paste0(out_dir, "subset_names_merged.just_microglia.seurat.RDS"))


