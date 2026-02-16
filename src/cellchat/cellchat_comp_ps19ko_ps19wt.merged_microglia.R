library(CellChat)
library(snakecase)
library(ggplot2)
library(ggalluvial)

out_dir <- "results/cellchat_comp_ps19ko_ps19wt.merged_microglia/"
dir.create(out_dir, showWarnings = F)

# load in cell chat data
cc_ps19ko <- readRDS("results/cellchat_ps19ko.merged_microglia/cell_chat_object.RDS")
cc_ps19wt <- readRDS("results/cellchat_ps19wt.merged_microglia/cell_chat_object.RDS")

# merge em up
cc_list <- list("PS19-WT"=cc_ps19wt,
                "PS19-KO"=cc_ps19ko)


cellChat <- mergeCellChat(cc_list, 
                          add.names=names(cc_list))

# compare the total number of interactions
gg1 <- compareInteractions(cellChat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellChat, show.legend = F, group = c(1,2), measure = "weight")

gg1 + gg2
ggsave(paste0(out_dir, "total_interaction_comparison.barplots.png"), width=8, height=5)


# differential number of interactions

# where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
pdf(paste0(out_dir, "total_comparison.number_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_diffInteraction(cellChat, weight.scale = T)
dev.off()

pdf(paste0(out_dir, "total_comparison.weight_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight")
dev.off()


cells <- unique(cellChat@idents$joint)

# create a folder for the cell plots
cell_dir <- paste0(out_dir, "cell_plots/")

dir.create(cell_dir, showWarnings = F)

for (cell in cells) {
  
  pdf(paste0(cell_dir, to_snake_case(cell), "_comparison.weight_of_interactions.circle_plot.pdf"), width=10, height=9)
  netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight",
                            sources.use = cell)
  dev.off()
  
}



gg1 <- netVisual_heatmap(cellChat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellChat, measure = "weight")
#> Do heatmap based on a merged object

pdf(paste0(out_dir, "total_comparison_interactions.heatmap.pdf"), width=10, height=6)
print(gg1 + gg2)
dev.off()


# Identify cell populations with significant changes in sending or receiving signals


num.link <- sapply(cc_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cc_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cc_list[[i]], title = names(cc_list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
ggsave(paste0(out_dir, "total_in_v_out.interaction_strength.scatter.png"), width=10, height=5)

# determine differential signaling
for (cell in cells) {
  
  netAnalysis_signalingChanges_scatter(cellChat, idents.use = cell)
  ggsave(paste0(cell_dir, to_snake_case(cell), ".pathway_signaling_changes.scatter.png"),
         width=7, height=5)
  
  
}

# cellChat <- computeNetSimilarityPairwise(cellChat, type = "functional")
# 
# cellChat <- netEmbedding(cellChat, type = "functional")
# 


# information flow 
gg1 <- rankNet(cellChat, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, 
               stacked = T, do.stat = TRUE) +
  theme(legend.position = "bottom")
gg2 <- rankNet(cellChat, mode = "comparison", measure = "weight", 
               sources.use = NULL, targets.use = NULL, 
               stacked = F, do.stat = TRUE) +
  theme(legend.position = "bottom")

gg1 + gg2
ggsave(paste0(out_dir, "pathway_information_flow_comp.png"), width=10, height=8)



