library(Seurat)
library(irGSEA)
library(msigdbr)  # For gene sets
library(dplyr)

# Initial setup ----------------------------------------------------------------

# create output directory for results
out_dir <- "results/aucell_irgsea_workflow/"
dir.create(out_dir, showWarnings = F)

# load in seurat objects
total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.no_low_count.seurat.RDS")

# Load in msigdbr and get gene sets --------------------------------------------
hallmark_gene_sets <- msigdbr(species = "Mus musculus", collection = "H")
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
kegg_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(hallmark=list_convert(hallmark_gene_sets),
                        kegg=list_convert(kegg_gene_sets),
                        reactome=list_convert(reactome_gene_sets),
                        wikipathways=list_convert(wikipath_gene_sets),
                        gobp=list_convert(gobp_gene_sets),
                        gomf=list_convert(gomf_gene_sets))

# Run irGSEA across genesets ---------------------------------------------------

print("Run Hallmark")

total_seu <- irGSEA.score(object = total_seu, 
                     assay = "RNA", 
                     slot = "data", 
                     seeds = 123, 
                     ncores = 4,
                     min.cells = 3, 
                     min.feature = 0,
                     custom = T, 
                     geneset = hallmark_gene_sets, 
                     method = "AUCell",
                     aucell.MaxRank = NULL, 
                     ucell.MaxRank = NULL, 
                     kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_hallmark")

print("Run KEGG")

total_seu <- irGSEA.score(object = total_seu, 
                          assay = "RNA", 
                          slot = "data", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = kegg_gene_sets, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_kegg")

print("Run Reactome")

total_seu <- irGSEA.score(object = total_seu, 
                          assay = "RNA", 
                          slot = "data", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = reactome_gene_sets, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_reactome")

print("Run Wikipathays")

total_seu <- irGSEA.score(object = total_seu, 
                          assay = "RNA", 
                          slot = "data", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = wikipath_gene_sets, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_wikipath")

print("Run GO BP")

total_seu <- irGSEA.score(object = total_seu, 
                          assay = "RNA", 
                          slot = "data", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = gobp_gene_sets, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_gobp")


print("Run GO MF")

total_seu <- irGSEA.score(object = total_seu, 
                          assay = "RNA", 
                          slot = "data", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = gomf_gene_sets, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

total_seu <- RenameAssays(total_seu, "AUCell", "AUCell_gomf")


# Save Seurat object -----------------------------------------------------------

SaveSeuratRds(total_seu, file=paste0(out_dir, "total_seurat.aucell_workflow.RDS"))
