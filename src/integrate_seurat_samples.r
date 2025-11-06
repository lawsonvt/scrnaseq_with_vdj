library(Seurat)
library(ggplot2)

out_dir <- "results/integrate_seurat_samples/"
dir.create(out_dir, showWarnings = F)


sample_dir <- "results/seurat_qc_samples/"
# read in seurat objects

seurat_files <- list.files(path=sample_dir, 
                           pattern=".qc_processed.seurat.RDS",
                           full.names = T)

# load in files
seurat_list <- lapply(seurat_files, function(file) {
  
  LoadSeuratRds(file)
  
})

# pull in names
names(seurat_list) <- unlist(lapply(seurat_list, function(x) {x@project.name}))

# merge them ...
seurat_merged <- merge(x=seurat_list[[1]],
                       y=seurat_list[-1],
                       add.cell.ids=names(seurat_list),
                       merge.data=T) # merge normalized data as well

# save the merged seurat
SaveSeuratRds(seurat_merged, paste0(out_dir, "all_samples.merged_seurat.RDS"))

# integrate!
# no normalization required, since we merge the normalized data as well
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)

# inspect elbow plot
ElbowPlot(seurat_merged, ndims=30) + 
  labs(title=sample_name) +
  scale_x_continuous(breaks=seq(0,30,5)) +
  scale_y_continuous(breaks=seq(0,30,5), limits=c(0,NA))
ggsave(paste0(out_dir, "merge_seurat.pca_elbow_plot.png"), width=6, height=5, bg="white")







