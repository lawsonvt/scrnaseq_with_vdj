library(Seurat)
library(SeuratObject)
library(CellChat)
library(snakecase)

out_dir <- "results/cellchat_ps19wt/"

dir.create(out_dir, showWarnings = F)

# load in object
total_seu <- LoadSeuratRds("results/seurat_cluster_naming.integrate_subsets/subset_names_merged.seurat.RDS")

unique(total_seu$orig.ident)

total_seu@meta.data$condition <- gsub("^[A-H]+\\-[A-Z0-9]+\\-", "", total_seu$orig.ident)

# subset the one we want
subset_seu <- subset(total_seu, subset = condition == "PS19-WT")

# drop total seurat object
rm(total_seu)
gc()

# clean up subset
subset_seu <- subset(subset_seu, subset = tcell_doublet == "singlet")

# drop pericytes (too few)
subset_seu <- subset(subset_seu, subset = merged_cell_name != "Pericytes")

# rename microglia
subset_seu$merged_cell_name <- gsub("^Microglia$", "Intermediate Microglia", subset_seu$merged_cell_name)

# prepare the cellchat object
Idents(subset_seu) <- "merged_cell_name"
subset_seu@meta.data$samples <- subset_seu@meta.data$orig.ident

# reorder the cell clusters
Idents(subset_seu) <- factor(as.character(Idents(subset_seu)),
                             levels=sort(unique(subset_seu$merged_cell_name)))


subset_seu@meta.data$samples <- factor(subset_seu@meta.data$orig.ident,
                                       levels=sort(unique(subset_seu@meta.data$orig.ident)))

# create the object
cellChat <- createCellChat(object = subset_seu, group.by = "ident", assay = "RNA")

# set the DB to be mouse
CellChatDB <- CellChatDB.mouse

showDatabaseCategory(CellChatDB)

# subset the Db to remove non-protein signaling
CellChatDB.use <- subsetDB(CellChatDB)

# add it to cell chat object
cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellChat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

options(future.globals.maxSize = 2.5 * 1e9) # increase the max future size from 500MB to 1GB

cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)

# smooth data (replaces projectData function, to apply the results to a PPI)
# https://github.com/jinworks/CellChat/issues/185

cellChat <- smoothData(cellChat, adj=PPI.mouse)

# Compute the communication probability and infer cellular communication network
cellChat <- computeCommunProb(cellChat, type = "triMean")

# filter out low cell count communication
cellChat <- filterCommunication(cellChat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellChat)

# Calculate the aggregated cell-cell communication network
cellChat <- aggregateNet(cellChat)

# make some interaction plots
groupSize <- as.numeric(table(cellChat@idents))


# NOTE: doing the plots combined did NOT work, so I am plotting them individually

pdf(paste0(out_dir, "total.number_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
dev.off()

pdf(paste0(out_dir, "total.weight_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")
dev.off()

# create a folder for the cell plots
cell_dir <- paste0(out_dir, "cell_plots/")

dir.create(cell_dir, showWarnings = F)

mat <- cellChat@net$weight

for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  cell <- rownames(mat)[i]
  pdf(paste0(cell_dir, to_snake_case(cell), ".weight_of_interactions.circle_plot.pdf"), width=8, height=7)
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}

# Infer the cell-cell communication at a signaling pathway level
cellChat <- computeCommunProbPathway(cellChat)

# pathway interactions

pathway_dir <- paste0(out_dir, "pathway_plots/")

dir.create(pathway_dir, showWarnings = F)

sig_pathways <- cellChat@netP$pathways

for (pathway in sig_pathways) {
  
  pdf(paste0(pathway_dir, pathway, ".chord_plot.pdf"), width=8, height=7)
  netVisual_aggregate(cellChat, signaling = pathway, 
                      layout = "chord", vertex.weight = NULL)
  dev.off()
  
  pdf(paste0(pathway_dir, pathway, ".heatmap_plot.pdf"), width=8, height=7)
  print(netVisual_heatmap(cellChat, signaling = pathway, 
                      color.heatmap = 'Reds'))
  dev.off()
  
}

# calculate centrality
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

for (pathway in sig_pathways) {
  
  
  pdf(paste0(pathway_dir, pathway, ".signaling_heatmap_plot.pdf"), width=8, height=7)
  print(netAnalysis_signalingRole_network(cellChat, 
                                  signaling = pathway, 
                                  width = 10, height = 4.5, font.size = 10))
  dev.off()
}



netAnalysis_signalingRole_scatter(cellChat)


# save the cell chat object
saveRDS(cellChat, paste0(out_dir, "cell_chat_object.RDS"))
