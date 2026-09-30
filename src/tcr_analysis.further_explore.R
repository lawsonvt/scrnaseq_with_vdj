library(scRepertoire)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(dplyr)
library(harmony)

# load in data
out_dir <- "results/tcr_analysis.further_explore/"

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

# annotate cell types
combined_tcr <- annotateInvariant(combined_tcr,
                                  type="MAIT",
                                  species = "mouse")

combined_tcr <- annotateInvariant(combined_tcr,
                                  type = "iNKT",
                                  species = "mouse")


combined_tcr_df <- bind_rows(combined_tcr)

clonalQuant(combined_tcr,
            cloneCall="gene", 
            group.by = "condition",
            chain = "both", 
            scale = TRUE) +
  theme(legend.position = "none") +
  labs(x=NULL) + ylim(0,100)

clonalAbundance(combined_tcr,
                cloneCall = "strict",
                scale=F,
                group.by = "condition")


clonalAbundance(combined_tcr,
                cloneCall = "strict",
                scale=T,
                group.by = "condition")


clonalLength(combined_tcr,
             cloneCall = "aa",
             scale=T,
             group.by = "condition")

# create a PS19 subset

combined_tcr_ps19 <- subsetClones(combined_tcr,
                                  name = "condition",
                                  variables=c("PS19-WT","PS19-KO"))

clonalLength(combined_tcr_ps19,
             cloneCall = "nt",
             scale=T,
             group.by = "condition")

clonalLength(combined_tcr_ps19,
             cloneCall = "nt",
             scale=T)


clonalCompare(combined_tcr_ps19,
              top.clones = 10,
              group.by = "condition",
              relabel.clones = F,
              cloneCall = "gene",
              palette = "viridis") 

clonalCompare(combined_tcr_ps19,
              top.clones = 10,
              group.by = "condition",
              relabel.clones = F,
              graph="area",
              cloneCall = "gene") +
  scale_fill_brewer(palette="Set1")

clonalCompare(combined_tcr_ps19,
              top.clones = 5,
              relabel.clones = F,
              graph="area",
              cloneCall = "gene") +
  scale_fill_brewer(palette="Set1")



clonalCompare(combined_tcr_ps19,
              top.clones = 10,
              group.by = "condition",
              relabel.clones = F,
              cloneCall = "aa") 



clonalCompare(combined_tcr_ps19,
              top.clones = 5,
              relabel.clones = F,
              cloneCall = "gene",
              palette = "viridis") 


clonalCompare(combined_tcr_ps19,
              top.clones = 50,
              relabel.clones = T,
              palette = "viridis",
              graph="area") 



clonalCompare(combined_tcr_ps19,
              top.clones = 5,
              group.by = "condition",
              relabel.clones = T,
              graph="area") +
  scale_fill_brewer(palette="Set3")




clonalScatter(combined_tcr_ps19,
              x.axis="PS19-WT",
              y.axis="PS19-KO",
              group.by="condition",
              cloneCall = "gene") +
  scale_fill_brewer(palette="Set1")




clonalScatter(combined_tcr_ps19,
              x.axis="D-Q619-PS19-WT",
              y.axis="A-Q637-PS19-KO") +
  scale_fill_brewer(palette="Set1")


clonalScatter(combined_tcr_ps19,
              x.axis="E-W136-PS19-WT",
              y.axis="H-E806-PS19-KO") +
  scale_fill_brewer(palette="Set1")


clonalHomeostasis(combined_tcr_ps19,
                  cloneCall = "gene",
                  group.by = "condition")

clonalHomeostasis(combined_tcr_ps19,
                  group.by = "sample")


clonalProportion(combined_tcr_ps19)
clonalProportion(combined_tcr_ps19,
                 cloneCall = "gene",
                 group.by="condition")



clonalOverlap(combined_tcr_ps19, method="jaccard", cloneCall = "gene")


clonalDiversity(combined_tcr_ps19,
                cloneCall = "gene",
                metric = "shannon",
                x.axis="condition")


clonalRarefaction(combined_tcr_ps19)

# Integrate Seurat object --------------------------

total_seu <- LoadSeuratRds("results/seurat_cluster_naming/cell_named.seurat.RDS")

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

# filter down to cells with clones

metadata <- total_seu@meta.data
metadata$cell_id <- rownames(metadata)

metadata$condition <- factor(gsub("[A-H]{1}\\-[A-Z0-9]{4}\\-", "", metadata$orig.ident),
                             levels=c("KO","WT","PS19-KO","PS19-WT"))

total_seu$cell_id <- metadata$cell_id
total_seu$condition <- metadata$condition


cell_clones <- metadata[!is.na(metadata$cloneSize),]$cell_id

clone_seu <- subset(total_seu, subset = cell_id %in% cell_clones)

DimPlot(clone_seu, group.by = c("cloneSize"), reduction = "umap.harmony") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))


DimPlot(clone_seu, group.by = c("cloneSize"), reduction = "umap.harmony",
        split.by = "condition", ncol=2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

# recluster
clone_seu <- SCTransform(clone_seu, vars.to.regress = c("percent.mt"), verbose = F)

clone_seu <- RunPCA(clone_seu, npcs = 50)

clone_seu <- RunHarmony(clone_seu, group.by.vars="orig.ident")

ElbowPlot(clone_seu, ndims=50) + 
  labs(title="Clone Subset") +
  scale_x_continuous(breaks=seq(0,50,5)) +
  scale_y_continuous(breaks=seq(0,50,5), limits=c(0,NA))


max_pc_dim <- 25

# cluster the harmonized data
clone_seu <- FindNeighbors(clone_seu, dims = 1:max_pc_dim, reduction = "pca")
clone_seu <- FindClusters(clone_seu, cluster.name = "clone_clusters",
                          resolution = 0.3)

# create umap
clone_seu <- RunUMAP(clone_seu, dims = 1:max_pc_dim, reduction="pca", reduction.name="umap.clone_pca")

DimPlot(clone_seu, group.by = c("cloneSize"), reduction = "umap.clone_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))

DimPlot(clone_seu, group.by = c("clone_clusters","cloneSize"), reduction = "umap.clone_pca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clone_cluster_size.umap.png"), width=12, height=4)

DimPlot(clone_seu, group.by = c("cloneSize"), reduction = "umap.clone_pca",
        split.by="condition", ncol=2) +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
ggsave(paste0(out_dir, "clone_size.condition_split.umap.png"), width=7, height=5)

DimPlot(clone_seu, group.by = c("clone_clusters"), reduction = "umap.clone_pca",
        split.by="condition", ncol=2) 
ggsave(paste0(out_dir, "clone_cluster.condition_split.umap.png"), width=7, height=5)

# find markers

all_markers <- FindAllMarkers(clone_seu)

saveRDS(all_markers, file=paste0(out_dir, "all_markers.RDS"))

SaveSeuratRds(clone_seu, file=paste0(out_dir, "clone_seurat.RDS"))

saveRDS(combined_tcr, file=paste0(out_dir, "combined_tcr.RDS"))

