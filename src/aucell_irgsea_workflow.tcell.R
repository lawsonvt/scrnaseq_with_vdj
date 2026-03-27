library(Seurat)
library(irGSEA)
library(msigdbr)  # For gene sets
library(dplyr)

# Initial setup ----------------------------------------------------------------

# create output directory for results
out_dir <- "results/aucell_irgsea_workflow.tcell/"
dir.create(out_dir, showWarnings = F)

# load in seurat objects
total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.no_low_count.seurat.RDS")
tcell_seu <- subset(total_seu, subset = cell_category == "T Cells")

rm(total_seu)
gc()

tcell_seu <- UpdateSeuratObject(tcell_seu)

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

tcell_seu <- irGSEA.score(object = tcell_seu, 
                     assay = "RNA", 
                     slot = "counts", 
                     seeds = 123, 
                     ncores = 4,
                     min.cells = 3, 
                     min.feature = 0,
                     custom = T, 
                     geneset = total_gene_sets$hallmark, 
                     method = "AUCell",
                     aucell.MaxRank = NULL, 
                     ucell.MaxRank = NULL, 
                     kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_hallmark")

print("Run KEGG")

tcell_seu <- irGSEA.score(object = tcell_seu, 
                          assay = "RNA", 
                          slot = "counts", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = total_gene_sets$kegg, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_kegg")

print("Run Reactome")

tcell_seu <- irGSEA.score(object = tcell_seu, 
                          assay = "RNA", 
                          slot = "counts", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = total_gene_sets$reactome, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_reactome")

print("Run Wikipathays")

tcell_seu <- irGSEA.score(object = tcell_seu, 
                          assay = "RNA", 
                          slot = "counts", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = total_gene_sets$wikipathways, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_wikipath")

print("Run GO BP")

tcell_seu <- irGSEA.score(object = tcell_seu, 
                          assay = "RNA", 
                          slot = "counts", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = total_gene_sets$gobp, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_gobp")


print("Run GO MF")

tcell_seu <- irGSEA.score(object = tcell_seu, 
                          assay = "RNA", 
                          slot = "counts", 
                          seeds = 123, 
                          ncores = 4,
                          min.cells = 3, 
                          min.feature = 0,
                          custom = T, 
                          geneset = total_gene_sets$gomf, 
                          method = "AUCell",
                          aucell.MaxRank = NULL, 
                          ucell.MaxRank = NULL, 
                          kcdf = 'Gaussian')

tcell_seu <- RenameAssays(tcell_seu, "AUCell", "AUCell_gomf")


# Save Seurat object -----------------------------------------------------------

SaveSeuratRds(tcell_seu, file=paste0(out_dir, "tcell_seurat.aucell_workflow.RDS"))



# redo umap
tcell_seu <- RunPCA(tcell_seu, npcs = 20)

max_pc_dim <- 10

tcell_seu <- RunUMAP(tcell_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.tcell_pca")

tcell_seu$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", tcell_seu$orig.ident)

tcell_seu$condition <- factor(as.character(tcell_seu$condition),
                              levels=c("WT", "KO", "PS19-WT", "PS19-KO"))

DimPlot(tcell_seu, reduction="umap.tcell_pca",
        group.by= "merged_cell_name",
        split.by = "condition", ncol=2)

irGSEA.density.scatterplot(object = tcell_seu,
                           method = "AUCell_hallmark",
                           show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
                           reduction = "umap.tcell_pca")

FeaturePlot(tcell_seu,
            features="HALLMARK-INFLAMMATORY-RESPONSE",
            reduction = "umap.tcell_pca",
            split.by = "condition", ncol=2)

DefaultAssay(tcell_seu) <- "AUCell_hallmark"
DefaultDimReduc(tcell_seu) <- "umap.tcell_pca"





