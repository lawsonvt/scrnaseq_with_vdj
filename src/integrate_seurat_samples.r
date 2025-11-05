library(Seurat)

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
                       add.cell.ids=names(seurat_list))
