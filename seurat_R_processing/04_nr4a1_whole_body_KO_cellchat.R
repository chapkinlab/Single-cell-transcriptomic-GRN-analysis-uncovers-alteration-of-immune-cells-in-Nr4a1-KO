library(dplyr)
library(Seurat)
library(ggplot2)
library(Matrix)
library(hdf5r)
library(dittoSeq)
library(SplineDV)
library(presto)
library(ggpubr)
library(tidyr)
library(stringr)
library(tibble)
library(cowplot)
library(openxlsx)
library(patchwork)
library(CellChat)

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))
table(data$CellType)
DefaultAssay(data) <- "RNA"
# large pops: b cells, plasma b cells, T cells, Macrophages, Fibroblasts, Colonocytes, SMCs
# small pops: nk cells, mast cells, EGCs 
# No statistical in nr4a1 ko: NK cells, ILCs,  VECs, LECs
# REMOVE: NK cells, ILCs,  VECs, LECs, Mast cells, EGCs 
# POSSIBLE: b cells, plasma b cells, T cells, Macrophages, Fibroblasts, Colonocytes, SMCs, ENs
# Immune + colonocytes: 'B cells', 'Plasma B cells', 'T cells', 'Macrophages', 'Colonocytes'
#cell_types_tokeep <- c('B cells', 'T cells', 'Macrophages', 'Colonocytes')
data_ko <- subset(data, Condition == 'Nr4a1 KO')
#data_ko <- subset(data_ko, subset = CellType %in% cell_types_tokeep)
data_ko$CellType <- factor(data_ko$CellType)
table(data_ko$CellType)
data_wt <- subset(data, Condition == 'WT')
#data_wt <- subset(data_wt, subset = CellType %in% cell_types_tokeep)
data_wt$CellType <- factor(data_wt$CellType)
table(data_wt$CellType)

# Assign sample labels explicitly and ensure they are factors
data_wt$samples <- factor("WT")
data_ko$samples <- factor("KO")

# Create CellChat objects
cellchat_WT <- createCellChat(object = data_wt, meta = data_wt@meta.data, group.by = "SubCellType", assay = "RNA")
cellchat_KO <- createCellChat(object = data_ko, meta = data_ko@meta.data, group.by = "SubCellType", assay = "RNA")
# cellchat_WT <- createCellChat(object = data_wt, meta = data_wt@meta.data, group.by = "CellType", assay = "RNA")
# cellchat_KO <- createCellChat(object = data_ko, meta = data_ko@meta.data, group.by = "CellType", assay = "RNA")

rm(data_ko, data_wt, data)

# Assign CellChat database
CellChatDB <- CellChatDB.mouse  # Use CellChatDB.human for human data
cellchat_WT@DB <- CellChatDB
cellchat_KO@DB <- CellChatDB

# Set CellType as identity
cellchat_WT <- setIdent(cellchat_WT, ident.use = "SubCellType")
cellchat_KO <- setIdent(cellchat_KO, ident.use = "SubCellType")
# cellchat_WT <- setIdent(cellchat_WT, ident.use = "CellType")
# cellchat_KO <- setIdent(cellchat_KO, ident.use = "CellType")

# Preprocess data
cellchat_WT <- subsetData(cellchat_WT)
#options(future.globals.maxSize = 4 * 1024^3) # Increase to 4 GB (adjust as needed)
#future::plan("multisession", workers = 4) # do parallel
cellchat_WT <- identifyOverExpressedGenes(cellchat_WT)
cellchat_WT <- identifyOverExpressedInteractions(cellchat_WT)
cellchat_WT <- computeCommunProb(cellchat_WT)

cellchat_KO <- subsetData(cellchat_KO)
#options(future.globals.maxSize = 4 * 1024^3) # Increase to 4 GB (adjust as needed)
#future::plan("multisession", workers = 4) # do parallel
cellchat_KO <- identifyOverExpressedGenes(cellchat_KO)
cellchat_KO <- identifyOverExpressedInteractions(cellchat_KO)
cellchat_KO <- computeCommunProb(cellchat_KO)

# Filter communication
cellchat_WT <- filterCommunication(cellchat_WT, min.cells = 50) # Start with a low value
cellchat_KO <- filterCommunication(cellchat_KO, min.cells = 50) # Start with a low value

# Compute at the pathway level
cellchat_WT <- computeCommunProbPathway(cellchat_WT)
cellchat_KO <- computeCommunProbPathway(cellchat_KO)

# Aggregate networks
cellchat_WT <- aggregateNet(cellchat_WT)
cellchat_KO <- aggregateNet(cellchat_KO)

# # Save rds from cell chat objects
# saveRDS(cellchat_WT, file = "cellchat_WT_broad_ct.rds")
# saveRDS(cellchat_KO, file = "cellchat_KO_broad_ct.rds")
# Load
# cellchat_WT <- readRDS("cellchat_WT_broad_ct.rds")
# cellchat_KO <- readRDS("cellchat_KO_broad_ct.rds")

# # Save rds from cell chat objects
# saveRDS(cellchat_WT, file = "cellchat_WT.rds")
# saveRDS(cellchat_KO, file = "cellchat_KO.rds")
# # Load
cellchat_WT <- readRDS("cellchat_WT.rds")
cellchat_KO <- readRDS("cellchat_KO.rds")

# e red(or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
cellchat_merged <- mergeCellChat(list(WT = cellchat_WT, KO = cellchat_KO), add.names = c("WT", "KO"))

# Compare communication networks (individual conditions)
gg1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged, measure = "weight", group = c(1,2))

# Corrected heatmaps for interaction comparison
gg3 <- netVisual_heatmap(cellchat_merged)  # Default "count"
gg4 <- compareInteractions(cellchat_merged, measure = "weight", show.legend = FALSE)
print(gg4)

# Custom visualization of differential interactions (Number of interactions in WT vs KO)
gg_diff <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "count")
print(gg_diff)

# Define your theme
gg3_diff <- netVisual_heatmap(cellchat_merged, font.size = 14, font.size.title = 17)
#gg4_diff <- netVisual_heatmap(cellchat_merged, measure = "weight", font.size = 15, font.size.title = 17)
gg_save <- gg3_diff #+ gg4_diff
print(gg_save)

#png('plots/de_cell_chat_broad_cells.png', height = 6, width = 12, res = 300, units = 'in')
#svg('plots/cellchat_img/de_cell_chat_broad_cells.svg', height = 4, width = 8)
svg('plots/cellchat_img/de_cell_chat_sub_cells.svg', height = 6.5, width = 7)
print(gg_save)
dev.off()

# Define the cell types you want to keep
# Define the specific SubCellTypes you want to keep
desired_subcell_types <- c(
  "TA", "Stem", 
  "B cells", "CD4+ T cells", "CD8+ T cells", "Cycling T cells",
  "DN T cells", "Macrophages",
  "Plasma B cells", "γδ T cells", "Fibroblasts"
)

# Subset the CellChat object using the SubCellType column
cellchat_merged <- subsetCellChat(cellchat_merged, idents.use = desired_subcell_types, group.by = "SubCellType")
# Verify the subsetted object
table(cellchat_merged@meta$CellType)
table(cellchat_merged@meta$SubCellType)

#sources <- c(2, 5, 6, 13, 14)
#targets <- c(2, 4, 5, 6, 8, 13, 14)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_diffInteraction( cellchat_merged, comparison = c(1,2),  weight.scale = T,
                                  #sources.use = sources,
                                  #targets.use = targets,
                                  # --- Font Size Changes ---
                                  vertex.label.cex = 1.2,  # Increase vertex label font size
                                  edge.label.cex = 1.2,    # Increase edge label font size (if label.edge = TRUE)
                                  # --- Radius/Size Changes ---
                                  vertex.size.max = 20,    # Increase maximum vertex size
                                  vertex.weight = 25,    # Optionally adjust base vertex weight
                                  arrow.size = 0.8,
                                  margin = 0.1           # Slightly reduce margin
)
#gg2 <- netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight")
print(gg1)
svg('plots/cellchat_img/de_cell_chat_sub_cells_circ.svg', height = 10, width = 10)
print(gg1)
dev.off()



# Define parameters if not already set
pos.dataset <- "KO"  # Define which dataset represents the positive condition
features.name <- "differential_genes"  # Name for identified genes (optional)

# Run the function on the merged object
cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged,
  group.dataset = "datasets",  # Column that defines WT vs KO
  pos.dataset = pos.dataset,   # Condition to test for overexpression (KO vs WT)
  features.name = features.name,
  only.pos = FALSE,  # Keep both up- and downregulated genes
  thresh.pc = 0.1,   # Minimum proportion of cells expressing gene
  thresh.fc = 0.05,  # 0.05, 0.1 or 0.3 Fold change threshold (adjust as needed)
  thresh.p = 0.05,   # p-value cutoff for significance
  group.DE.combined = FALSE  # Compare datasets separately
)

# Map the results of differential expression analysis onto the inferred cell-cell communications
net <- netMappingDEG(cellchat_merged, features.name = features.name, variable.all = TRUE)
write.csv(net, "plots/cellchat_img/net_mapping_DEG_results_sub_celltypes.csv", row.names = FALSE)

# Extract the ligand-receptor pairs with upregulated ligands in KO (positive condition)
net.up <- subsetCommunication(cellchat_merged, net = net, datasets = "KO", 
                              ligand.logFC = 0.05, receptor.logFC = NULL, 
                              thresh = 0.05, ligand.pvalues = 0.05 )

# Extract the ligand-receptor pairs with downregulated ligands in WT (negative condition)
net.down <- subsetCommunication(cellchat_merged, net = net, datasets = "WT", 
                                ligand.logFC = -0.05, receptor.logFC = NULL,
                                thresh = 0.05, ligand.pvalues = 0.05 )

# Filter small probs
net.up <- net.up[ net.up$prob > 0.01, ]
net.down <- net.down[ net.down$prob > 0.01, ]

# Extract the genes associated with the upregulated and downregulated ligand-receptor pairs
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merged)
print(gene.up)
print(gene.down)

source_name <- levels(cellchat_merged@idents$KO)  # Source names
# sources <- c(1, 2, 3, 6, 7, 8, 9, 10)
# targets <- c(1, 2, 3, 6, 7, 8, 9, 10)
# 
# # Loop for both upregulated and downregulated signaling
# for (source in sources){
#   try({
#     # Generate bubble plot for upregulated signaling
#     pairLR.use.up <- net.up[, "interaction_name", drop = F]
#     gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = source, targets.use = targets, comparison = c(1, 2),
#                             remove.isolate = T, color.heatmap = "viridis",
#                             title.name = paste0("Up-regulated signaling in \n KO (Source ", source_name[source], ")"),
#                             font.size = 19,
#                             font.size.title = 21,
#                             angle.x = 45
#                             )
#     
#     # Save the upregulated plot for each source
#     plot_name_up <- paste0("plots/cellchat_img/upregulated_plot_source_", source_name[source], "_DEup_signaling.svg")
#     ggsave(plot_name_up, plot = gg1, height = 10, width = 9, device = "svg")
#     print(gg1)
#     dev.off()
#   }, silent = TRUE)  # If an error occurs in the upregulated plot, it will silently skip to the downregulated plot
#   
#   # try({
#   #   # Generate bubble plot for downregulated signaling
#   #   pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]
#   #   nrows <- nrow(pairLR.use.up)
#   #   gg2 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.down, sources.use = source, targets.use = targets, comparison = c(1, 2),
#   #                           remove.isolate = T, color.heatmap = "viridis",
#   #                           title.name = paste0("Down-regulated signaling \n in KO (Source ", source_name[source], ")"),
#   #                           font.size = 17,
#   #                           font.size.title = 19, 
#   #                           angle.x = 45
#   #                           )
#   # 
#   #   # Save the downregulated plot for each source
#   #   plot_name_down <- paste0("plots/cellchat_img/downregulated_plot_source_", source_name[source], "_DEdn_signaling.svg")
#   #   if (nrows  <= 4 ){ 
#   #     #png(plot_name_down, height = 4, width = 5, res = 300, units = 'in')
#   #     svg(plot_name_down, height = 4, width = 5)
#   #   } else{#png(plot_name_down, height = 8, width = 10, res = 300, units = 'in')
#   #           svg(plot_name_down, height = 8, width = 10)
#   #   }
#   #   print(gg2)
#   #   dev.off()
#   # }, silent = TRUE)  # If an error occurs in the downregulated plot, it will silently skip to the next source
# }




# Saving only Fibroblasts DEGs comm
source<- 6
print(source_name[source])
# Generate bubble plot for upregulated signaling
pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = source, targets.use = targets, comparison = c(1, 2),
                        remove.isolate = T, color.heatmap = "viridis",
                        title.name = paste0("Up-regulated signaling in \n KO (Source ", source_name[source], ")"),
                        font.size = 19,
                        font.size.title = 21,
                        angle.x = 45
)
# === Add this line to increase legend text and title size ===
gg1 <- gg1 + theme(
  plot.title = element_text(hjust = 0.5), # Centers the main plot title
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16, hjust = 0.5) # Centers the legend titles
)+ guides(
  color = guide_colorbar(title = "Commun.\nProb.")
)

gg1

# Save the upregulated plot for each source
plot_name_up <- paste0("plots/cellchat_img/upregulated_plot_source_", source_name[source], "_DEup_signaling.svg")
ggsave(plot_name_up, plot = gg1, height = 14, width = 9, device = "svg")
print(gg1)
dev.off()


# Saving only CD8+ T cells cells DEGs comm
source<- 3
print(source_name[source])
# Generate bubble plot for upregulated signaling
pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = source, targets.use = targets, comparison = c(1, 2),
                        remove.isolate = T, color.heatmap = "viridis",
                        title.name = paste0("Up-regulated signaling in \n KO (Source ", source_name[source], ")"),
                        font.size = 19,
                        font.size.title = 21,
                        angle.x = 45
)
# === Add this line to increase legend text and title size ===
gg1 <- gg1 + theme(
  plot.title = element_text(hjust = 0.5), # Centers the main plot title
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16, hjust = 0.5) # Centers the legend titles
)+ guides(
  color = guide_colorbar(title = "Commun.\nProb.")
)


gg1
# Save the upregulated plot for each source
plot_name_up <- paste0("plots/cellchat_img/upregulated_plot_source_", source_name[source], "_DEup_signaling.svg")
ggsave(plot_name_up, plot = gg1, height = 12, width = 9, device = "svg")
print(gg1)
dev.off()


# Saving only Fibroblasts DEGs comm
source<- 9
print(source_name[source])
# Generate bubble plot for upregulated signaling
pairLR.use.up <- net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = source, targets.use = targets, comparison = c(1, 2),
                        remove.isolate = T, color.heatmap = "viridis",
                        title.name = paste0("Up-regulated signaling in \n KO (Source ", source_name[source], ")"),
                        font.size = 19,
                        font.size.title = 21,
                        angle.x = 45
)
# === Add this line to increase legend text and title size ===
gg1 <- gg1 + theme(
  plot.title = element_text(hjust = 0.5), # Centers the main plot title
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16, hjust = 0.5) # Centers the legend titles
)+ guides(
  color = guide_colorbar(title = "Commun.\nProb.")
)

gg1

# Save the upregulated plot for each source
plot_name_up <- paste0("plots/cellchat_img/upregulated_plot_source_", source_name[source], "_DEup_signaling.svg")
ggsave(plot_name_up, plot = gg1, height = 13, width = 9, device = "svg")
print(gg1)
dev.off()

