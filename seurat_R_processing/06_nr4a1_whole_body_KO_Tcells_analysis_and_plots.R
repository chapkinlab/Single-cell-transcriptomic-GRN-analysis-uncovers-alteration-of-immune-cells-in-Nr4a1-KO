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
library(patchwork)
library(openxlsx)
library(AUCell)
library(ggsignif)

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))
DefaultAssay(data) <- "RNA"

process_and_extract_cell_type <- function(data, cell_type_name, assay_name = "RNA", num_hvg = 2000, dims_pca = 50, resolution = 1.0) {
  DefaultAssay(data) <- assay_name
  data_cell_type <- subset(data, subset = CellType == cell_type_name)
  data_cell_type <- FindVariableFeatures(data_cell_type, selection.method = "vst", nfeatures = num_hvg)
  data_cell_type <- ScaleData(data_cell_type)
  data_cell_type <- RunPCA(data_cell_type)
  data_cell_type <- RunUMAP(data_cell_type, dims = 1:dims_pca, n.epochs = 500)
  data_cell_type <- FindNeighbors(data_cell_type, dims = 1:dims_pca)
  data_cell_type <- FindClusters(data_cell_type, resolution = resolution)
  return(data_cell_type)
}

# Substract t cells
data_tcells <- process_and_extract_cell_type(data, "T cells", resolution = 1.0)
# Re-factor for plotting 
data_tcells$SubCellType = factor(data_tcells$SubCellType,
                                 levels = c('DN T cells', 'CD4+ T cells', 'CD8+ T cells',
                                            'γδ T cells', 'Cycling T cells'))
# Save t cells rds
saveRDS(data_tcells, file.path(base_dir, "tcells_only.rds"))

# Load this for faster processing
data_tcells <- readRDS(file.path(base_dir, "tcells_only.rds"))

#data_cell_type <- subset(data, subset = CellType == cell_type_name)
Idents(data_tcells) <- "SubCellType"


# UMAP for SubCellType
um_tcells <- DimPlot(object = data_tcells, reduction = "umap", group.by = "SubCellType", 
                     label = T, repel = T, label.size = 7, 
                     label.box = T, alpha = 1, pt.size = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    legend.text = element_text(size = 16),  # Increase legend label size
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) + NoLegend() +
  coord_cartesian(clip = "off")

um_tcells
# Save BatchID UMAP
#png('plots/UMAP_tcells.png', height = 6, width = 8, res = 300, units = 'in')
svg('plots/UMAP_tcells.svg', height = 6, width = 8)
print(um_tcells)
dev.off()

# UMAP for Condition
um_tcells <- DimPlot(object = data_tcells, reduction = "umap", group.by = "Condition", 
                     label = T, repel = T, label.size = 7, 
                     label.box = T, alpha = 1, pt.size = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) + NoLegend() +
  coord_cartesian(clip = "off")

um_tcells
# Save BatchID UMAP
#png('plots/UMAP_tcells_condition.png', height = 6, width = 8, res = 300, units = 'in')
svg('plots/UMAP_tcells_condition.svg', height = 6, width = 8)
print(um_tcells)
dev.off()

##########
# Cell Scores for T cell exhaustion
##########
# Exhaustion
genelists = read.xlsx('Genesets_cell_scores.xlsx')
Exhaustion_chen = list(Exhaustion = intersect(genelists$`Chen.et.al.(2023).T-cell.exhaustion`, rownames(data_tcells)))
data_tcells$exhaustion_chen_score = AUCell_run(GetAssayData(data_tcells),Exhaustion_chen,
                               normAUC = T)@assays@data@listData$AUC[1,]

Exhaustion_cust0 = list(Exhaustion = intersect(genelists$`T-cell.exhaustion.cust0`, rownames(data_tcells)))
data_tcells$exhaustion_cust0_score = AUCell_run(GetAssayData(data_tcells),Exhaustion_cust0,
                                               normAUC = T)@assays@data@listData$AUC[1,]

Exhaustion_cust1 = list(Exhaustion = intersect(genelists$`T-cell.exhaustion.cust1`, rownames(data_tcells)))
data_tcells$exhaustion_cust1_score = AUCell_run(GetAssayData(data_tcells),Exhaustion_cust1,
                                                normAUC = T)@assays@data@listData$AUC[1,]

Exhaustion_cust2 = list(Exhaustion = intersect(genelists$`T-cell.exhaustion.cust2`, rownames(data_tcells)))
data_tcells$exhaustion_cust2_score = AUCell_run(GetAssayData(data_tcells),Exhaustion_cust2,
                                                normAUC = T)@assays@data@listData$AUC[1,]

Exhaustion_cust3 = list(Exhaustion = intersect(genelists$`T-cell.exhaustion.cust3`, rownames(data_tcells)))
data_tcells$exhaustion_cust3_score = AUCell_run(GetAssayData(data_tcells),Exhaustion_cust3,
                                                normAUC = T)@assays@data@listData$AUC[1,]

# T Cell Activation (GO:0042110) from DV
# genes <- c("Lfng","Cd74","Cd4","Cd8a","Trbc1","Kit","Cd7","Trdc","Gata3","Lcp1","Cd247","Cd3e")
# genes = list(genes = intersect(genes, rownames(data_tcells)))
# data_tcells$score_tcell_activation = AUCell_run(GetAssayData(data_tcells),genes,
#                                normAUC = T)@assays@data@listData$AUC[1,]
# 
# # Antigen Receptor-Mediated Signaling Pathway (GO:0050851) from DV
# genes <- c("Ighm", "Gata3", "Thy1", "Cd3e", "Ighg1", "Cd8a", "Inpp5d", "Trbc1",
#            "Ctla4", "Trdc", "Iglc1", "Cd247", "Nfkbid")
# genes = list(genes = intersect(genes, rownames(data_tcells)))
# data_tcells$score_antigen = AUCell_run(GetAssayData(data_tcells),genes,
#                                                 normAUC = T)@assays@data@listData$AUC[1,]

# T Cell Activation Involved in Immune Response (GO:0002286) from DV + DE
genes <- c("Cd74")
genes = list(genes = intersect(genes, rownames(data_tcells)))
data_tcells$score_act = AUCell_run(GetAssayData(data_tcells),genes,
                                       normAUC = T)@assays@data@listData$AUC[1,]

#Positive Regulation of Cytokine-Mediated Signaling Pathway (GO:0001961) from DV + DE
genes <- c("Cd74", "Sting1")
genes = list(genes = intersect(genes, rownames(data_tcells)))
data_tcells$score_pos_reg_cytokine = AUCell_run(GetAssayData(data_tcells),genes,
                                       normAUC = T)@assays@data@listData$AUC[1,]


# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 16), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 13),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 16)
)

p6 = data_tcells[[c('exhaustion_chen_score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = exhaustion_chen_score , fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p6
table(data_tcells$Condition, data_tcells$SubCellType)
#png('plots/Tcell_exhaustion_Condition_subct.png', height = 6, width = 7, res = 300, units = 'in')
svg('plots/Tcell_exhaustion_Condition_subct.svg', height = 3, width = 10)
p6
dev.off()

p = data_tcells[[c('score_act','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = score_act, fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T Cell Activation Involved in Immune Response') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_act_immune_Condition_subct.svg', height = 3, width = 10)
p
dev.off()

p = data_tcells[[c('score_pos_reg_cytokine','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = score_pos_reg_cytokine, fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('Positive Regulation of Cytokine-Mediated Signaling Pathway') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_pos_reg_cyto_Condition_subct.svg', height = 3, width = 10)
p
dev.off()


p = data_tcells[[c('exhaustion_cust0_score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = exhaustion_cust0_score  , fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion custom-0') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_exhaustion_Condition_subct_cust0.svg', height = 3, width = 10)
p
dev.off()

p = data_tcells[[c('exhaustion_cust1_score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = exhaustion_cust1_score , fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion custom-1') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_exhaustion_Condition_subct_cust1.svg', height = 3, width = 10)
p
dev.off()

p = data_tcells[[c('exhaustion_cust2_score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = exhaustion_cust2_score , fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion custom-2') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_exhaustion_Condition_subct_cust2.svg', height = 3, width = 10)
p
dev.off()

p = data_tcells[[c('exhaustion_cust3_score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = exhaustion_cust3_score , fill = Condition))+ 
  geom_violin() + 
  geom_boxplot(color='black', width=0.5, outlier.size = 0.3)+ 
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion custom-3') + ylab('AUCell score')+
  facet_wrap(~ SubCellType, ncol = 5, scales = 'free') + cust_theme + NoLegend()
p
svg('plots/Tcell_exhaustion_Condition_subct_cust3.svg', height = 3, width = 10)
p
dev.off()
# Create DittoSeq plot
cc <- dittoBarPlot(data_tcells, var = 'SubCellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("T cell counts")  # Set title correctly 
 

cc
# Save plot
svg('plots/tcell_counts_subct.svg', height = 4, width = 7)
#png('plots/tcell_counts_subct.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()

# Create DittoSeq plot
cc <- dittoBarPlot(data_tcells, var = 'SubCellType', group.by = 'Condition', scale = 'percent') +
  cust_theme + ggtitle("T cell counts (%)")  # Set title correctly 


cc
# Save plot
svg('plots/tcell_counts_pct_subct.svg', height = 4, width = 7)
#png('plots/tcell_counts_subct.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()


# Markers for annotation (Done at this point)
DotPlot(data_tcells, features = c('Cd8a','Cd8b1','Cd4', 'Cd3e', 'Cd3g','Ncr1','Klrb1b','Trdc', 'Gata3','Ms4a2','Kit','Tpm2'),
        group.by = 'seurat_clusters', dot.min = 0.05, cols = 'RdBu' ) + # Set color range
  coord_flip()

DotPlot(data_tcells, features = c('Cd8a','Cd8b1',
                                  'Cd4', 
                                  'Cd3e', 'Cd3g',
                                  'Ncr1','Klrb1b',
                                  'Batf3','',
                                  'Trdc', 
                                  'Gata3',
                                  'Ms4a2','Kit',
                                  'Tpm2', "Myh11", "Tagln",      
                                  'Pecam1',"Lyve1", "Flt4", 
                                  'Plvap', "Flt1",
                                  "Vip", "Tubb3","Vat1l", 
                                  "Mki67"),
        group.by = 'seurat_clusters', dot.min = 0.05, cols = 'RdBu' ) + # Set color range
  coord_flip()

#> Gamma-delta T cells <- Trdc, Cd3, Cd4-, Cd8-
gene_marker_tcells <- c('Cd3e', 'Cd3g',
                        'Cd4', 
                        'Cd8a','Cd8b1',
                        'Trdc','Trgc1', 
                        "Mki67")

## Canonical marker dotplots ####
p1 = DotPlot(data_tcells, features = gene_marker_tcells, dot.min = 0.05,
             dot.scale = 8, scale = F,
             group.by = 'SubCellType')+
  scale_color_gradientn(colours = c('grey95','firebrick3'))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), # Increase x-axis text size
    axis.text.y = element_text(size = 18),  # Increase y-axis text size
    axis.title.x = element_text(size = 16), # Increase x-axis label size
    axis.title.y = element_text(size = 16)  # Increase y-axis label size
  )
p1
png('plots/Tcells_Dot.png', height = 4, width = 8, res = 300, units = 'in')
#svg('plots/Tcells_Dot.svg', height = 4, width = 8)
p1
dev.off()



########################################
## Differential expression analysis 
########################################

# Ensure `cell_types` is properly extracted
cell_types <- unique(data_tcells$SubCellType)
# Set the grouping variable for Seurat
Idents(data_tcells) <- "Condition"
# Define output directory
output_dir <- "DE_results_genes/tcells"
# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cell in cell_types) {
  # Print the current cell type being processed
  cat("Processing cell type:", cell, "\n")
  
  # Subset the Seurat object for the current cell type
  cell_subset <- subset(data_tcells, SubCellType == cell)
  
  # Perform FindMarkers for this cell type comparing WT vs Nr4a1 KO
  markers <- FindMarkers(cell_subset,
                         ident.1 = "Nr4a1 KO",
                         ident.2 = "WT",
                         test.use = "wilcox",  # Wilcoxon rank-sum test
                         logfc.threshold = 0.1,
                         min.pct = 0.01,
                         min.diff.pct = -Inf,
                         verbose = TRUE)
  
  # Filter significant markers
  #all_markers <- markers[markers$p_val < 0.05, ]
  all_markers <- markers[markers$p_val_adj < 0.05, ]
  
  # Sort by absolute value of log2FC
  all_markers <- all_markers[order(abs(all_markers$avg_log2FC), decreasing = TRUE), ]
  
  upregulated <- all_markers[all_markers$avg_log2FC > 0, ]
  downregulated <- all_markers[all_markers$avg_log2FC < 0, ]
  
  # Add gene names as a column
  all_markers$Gene <- rownames(all_markers)
  upregulated$Gene <- rownames(upregulated)
  downregulated$Gene <- rownames(downregulated)
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add sheets to the workbook
  addWorksheet(wb, "All_Pval_0.05")
  addWorksheet(wb, "Up_regulated")
  addWorksheet(wb, "Down_regulated")
  
  # Write data to respective sheets
  writeData(wb, sheet = "All_Pval_0.05", all_markers)
  writeData(wb, sheet = "Up_regulated", upregulated)
  writeData(wb, sheet = "Down_regulated", downregulated)
  
  # Ensure consistent file naming
  save_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(save_cell_name, "_markers_Wilcoxon.xlsx"))
  
  # Save workbook
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  cat("Saved:", output_file, "\n")
}





library(enrichR)
library(openxlsx)
library(dplyr)

# Define paths
input_dir <- "DE_results_genes/tcells"
output_dir <- "Enrichr_results/tcells"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
dbs <- "GO_Biological_Process_2025"
pval_cutoff <- 0.05
top_n <- 100

# Function to perform Enrichr and process results with waiting time
run_enrichr <- function(genes, direction, dbs, pval_cutoff, wait_time = 2) {
  results_df <- data.frame()
  
  for (db in dbs) {
    if (length(genes) > 0) { # Check if genes vector is not empty
      enrichr_res_list <- enrichr(genes, db)
      if (!is.null(enrichr_res_list) && length(enrichr_res_list) > 0) {
        enrichr_res <- enrichr_res_list[[1]]
        
        if (!is.null(enrichr_res) && nrow(enrichr_res) > 0) {
          enrichr_res <- enrichr_res %>% filter(Adjusted.P.value < pval_cutoff)
          
          if (nrow(enrichr_res) > 0) {
            results_df <- bind_rows(results_df, enrichr_res %>%
                                      transmute(database = db,
                                                Term, Overlap, `Odds.Ratio`,
                                                P.value, Adjusted.P.value,
                                                Genes)) # Use transmute
          }
        }
      }
    } else {
      print(paste("No genes provided for", direction))
    }
    Sys.sleep(wait_time)  # Prevent API overload
  }
  return(results_df)
}

# Assuming 'data_tcells' is a dataframe containing T cell data
# and has a column named 'SubCellType'

# Get unique cell types
if (!exists("data_tcells") || !("SubCellType" %in% colnames(data_tcells))) {
  stop("Error: 'data_tcells' dataframe with 'SubCellType' column not found.")
}
cell_types <- unique(data_tcells$SubCellType)

for (cell in cell_types) {
  print(paste("Working with", cell))
  file_path <- file.path(input_dir, paste0(gsub(" ", "_", cell), "_markers_Wilcoxon.xlsx"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }
  
  wb <- loadWorkbook(file_path)
  upregulated <- read.xlsx(wb, sheet = "Up_regulated")
  downregulated <- read.xlsx(wb, sheet = "Down_regulated")
  
  top_upregulated <- upregulated %>% head(top_n) %>% pull(Gene)
  top_downregulated <- downregulated %>% head(top_n) %>% pull(Gene)
  
  enrichr_results_up <- run_enrichr(top_upregulated, paste0("Up-regulated ", cell), dbs, pval_cutoff)
  enrichr_results_dn <- run_enrichr(top_downregulated, paste0("Down-regulated ", cell), dbs, pval_cutoff)
  
  # Skip file writing if no significant pathways are found
  if (nrow(enrichr_results_up) == 0 && nrow(enrichr_results_dn) == 0) {
    print(paste("Skipping file for", cell, "- No significant pathways found"))
    next
  }
  
  # Write results to Excel
  output_file <- file.path(output_dir, paste0(gsub("[^[:alnum:]_]", "", gsub(" ", "_", cell)), "_enrichr_Wilcoxon_results.xlsx"))
  wb_out <- createWorkbook()
  
  if (nrow(enrichr_results_up) > 0) {
    addWorksheet(wb_out, "Upregulated_Pathways")
    writeData(wb_out, "Upregulated_Pathways", enrichr_results_up)
  }
  
  if (nrow(enrichr_results_dn) > 0) {
    addWorksheet(wb_out, "Downregulated_Pathways")
    writeData(wb_out, "Downregulated_Pathways", enrichr_results_dn)
  }
  
  # Only save workbook if at least one pathway exists
  saveWorkbook(wb_out, output_file, overwrite = TRUE)
  print(paste("Enrichment results saved in", output_file))
}

#genes <- c("Lfng","Cd74","Cd4","Cd8a","Trbc1","Kit","Cd7","Trdc","Gata3","Lcp1","Cd247","Cd3e")
genes <- c("Lfng","Cd74","Cd4","Cd8a","Trbc1","Kit")
genes <- c("Cd7","Trdc","Gata3","Lcp1","Cd247","Cd3e")

VlnPlot(data_tcells, features = genes, group.by = "SubCellType", split.by = "Condition")
