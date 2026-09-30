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

groupSize <- as.numeric(table(cc_ps19ko@idents)) +
  as.numeric(table(cc_ps19wt@idents))

# where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
pdf(paste0(out_dir, "total_comparison.number_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_diffInteraction(cellChat, weight.scale = T, vertex.weight = groupSize)
dev.off()

pdf(paste0(out_dir, "total_comparison.weight_of_interactions.circle_plot.pdf"), width=10, height=9)
netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight", vertex.weight = groupSize)
dev.off()


cells <- unique(cellChat@idents$joint)

# create a folder for the cell plots
cell_dir <- paste0(out_dir, "cell_plots/")

dir.create(cell_dir, showWarnings = F)

for (cell in cells) {
  
  pdf(paste0(cell_dir, to_snake_case(cell), "_comparison.weight_of_interactions.circle_plot.pdf"), width=7, height=5)
  netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight",
                            sources.use = cell, vertex.weight = groupSize)
  dev.off()
  
  pdf(paste0(cell_dir, to_snake_case(cell), "_comparison.weight_of_interactions.circle_plot.no_labels.pdf"), width=7, height=5)
  netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight",
                            sources.use = cell, vertex.weight = groupSize, vertex.label.cex = 0.000001)
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


pathway_data <- rankNet(cellChat, mode = "comparison", measure = "weight", 
                        sources.use = NULL, targets.use = NULL, 
                        stacked = F, do.stat = TRUE, return.data = T)
pathway_data <- pathway_data$signaling.contribution

top_pathways <- unique(pathway_data[pathway_data$pvalues < 0.05 &
                                      pathway_data$contribution.scaled > 3,]$name)

top_pathway_data <- pathway_data[pathway_data$name %in% top_pathways,]


# add in scaling data

top_pathway_data$total <- 0

for (top_pathway in top_pathways) {
  
  contribution_sum <- sum(top_pathway_data[top_pathway_data$name == top_pathway,]$contribution.scaled)
  
  top_pathway_data[top_pathway_data$name == top_pathway,]$total <- contribution_sum
  
}

top_pathway_data$contribution_rel <- top_pathway_data$contribution.scaled / 
  top_pathway_data$total
top_pathway_data$contribution_rel_d <- top_pathway_data$contribution_rel
top_pathway_data[top_pathway_data$group == "PS19-KO",]$contribution_rel_d <-
  top_pathway_data[top_pathway_data$group == "PS19-KO",]$contribution_rel_d * -1

top_pathway_data <- top_pathway_data[order(top_pathway_data$contribution_rel_d),]
top_pathway_data$name <- factor(as.character(top_pathway_data$name),
                                levels=unique(as.character(top_pathway_data$name)))

ggplot(top_pathway_data,
       aes(y=name,
           x=contribution_rel,
           fill=group)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(y=NULL, x="Relative Information Flow", fill=NULL) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#b2182b","#2166ac"))
ggsave(paste0(out_dir, "pathway_information_flow_comp.relative_filtered.png"), width=7, height=5)

ggplot(top_pathway_data,
       aes(y=name,
           x=contribution.scaled,
           fill=group)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  labs(y=NULL, x="Information Flow", fill=NULL) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#b2182b","#2166ac"))
ggsave(paste0(out_dir, "pathway_information_flow_comp.filtered.png"), width=7, height=5)

# more plots, trying to discern ligand pairs

# CD8+ plots

netVisual_bubble(cellChat,
                 sources.use = "CD8+",
                 comparison=c(1,2),
                 color.text = c("#2166ac","#b2182b"),
                 angle.x = 45)
ggsave(paste0(out_dir, "cd8.total_ligand_bubble_plot.png"), width=11, height=11)



gg1 <- netVisual_bubble(cellChat,
                 sources.use = "CD8+",
                 targets.use = "Neutrophils",
                 comparison=c(1,2),
                 angle.x = 45,
                 max.dataset = 2,
                 color.text = c("#2166ac","#b2182b"),
                 title.name = "Increased signaling in PS19-KO",
                 remove.isolate = T)

gg2 <- netVisual_bubble(cellChat,
                 sources.use = "CD8+",
                 targets.use = "Neutrophils",
                 comparison=c(1,2),
                 angle.x = 45,
                 max.dataset = 1,
                 color.text = c("#2166ac","#b2182b"),
                 title.name = "Decreased signaling in PS19-KO",
                 remove.isolate = T)

gg1 + gg2
ggsave(paste0(out_dir, "cd8_to_neutrophils.ligand_bubble_plot.png"), width=8, height=6)

# same for Naive T Cells

netVisual_bubble(cellChat,
                 sources.use = "Naive T Cells",
                 comparison=c(1,2),
                 color.text = c("#2166ac","#b2182b"),
                 angle.x = 45)
ggsave(paste0(out_dir, "naive_t.total_ligand_bubble_plot.png"), width=11, height=11)


gg1 <- netVisual_bubble(cellChat,
                        sources.use = "Naive T Cells",
                        targets.use = "Neutrophils",
                        comparison=c(1,2),
                        angle.x = 45,
                        max.dataset = 2,
                        color.text = c("#2166ac","#b2182b"),
                        title.name = "Increased signaling in PS19-KO",
                        remove.isolate = T)

gg2 <- netVisual_bubble(cellChat,
                        sources.use = "Naive T Cells",
                        targets.use = "Neutrophils",
                        comparison=c(1,2),
                        angle.x = 45,
                        max.dataset = 1,
                        color.text = c("#2166ac","#b2182b"),
                        title.name = "Decreased signaling in PS19-KO",
                        remove.isolate = T)

gg1 + gg2
ggsave(paste0(out_dir, "naive_t_to_neutrophils.ligand_bubble_plot.png"), width=8, height=6)


gg1 <- netVisual_bubble(cellChat,
                        sources.use = "Naive T Cells",
                        targets.use = "B Cells",
                        comparison=c(1,2),
                        angle.x = 45,
                        max.dataset = 2,
                        color.text = c("#2166ac","#b2182b"),
                        title.name = "Increased signaling in PS19-KO",
                        remove.isolate = T)

gg2 <- netVisual_bubble(cellChat,
                        sources.use = "Naive T Cells",
                        targets.use = "B Cells",
                        comparison=c(1,2),
                        angle.x = 45,
                        max.dataset = 1,
                        color.text = c("#2166ac","#b2182b"),
                        title.name = "Decreased signaling in PS19-KO",
                        remove.isolate = T)

gg1 + gg2
ggsave(paste0(out_dir, "naive_t_to_b_cells.ligand_bubble_plot.png"), width=8, height=6)

# using differential expression analysis


df <- netVisual_bubble(cellChat,
                       sources.use = "Naive T Cells",
                       targets.use = "B Cells",
                       comparison=c(1,2),
                       angle.x = 45,
                       max.dataset = 2,
                       color.text = c("#2166ac","#b2182b"),
                       title.name = "Increased signaling in PS19-KO",
                       remove.isolate = T,
                       return.data = T)




