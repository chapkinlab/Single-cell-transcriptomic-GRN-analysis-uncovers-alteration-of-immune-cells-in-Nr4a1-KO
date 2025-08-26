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
library(readxl)

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

# Function to process and extract a specific cell type
process_rna <- function(data, assay_name = "RNA", num_hvg = 2000, dims_pca = 50, resolution = 1.0) {
  DefaultAssay(data) <- assay_name
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = num_hvg)
  data <- ScaleData(data)
  data <- RunPCA(data)
  data <- RunUMAP(data, dims = 1:dims_pca, n.epochs = 500)
  data <- FindNeighbors(data, dims = 1:dims_pca)
  data <- FindClusters(data, resolution = resolution)
  return(data)
}

process_sct_rna <- function(data, assay_name = "SCT", dims_pca = 50, resolution = 1.0) {
  data <- SCTransform(data, vst.flavor = "v2")
  DefaultAssay(data) <- assay_name
  data <- RunPCA(data, assay = "SCT", npcs = dims_pca, verbose = FALSE)
  data <- RunUMAP(data, dims = 1:dims_pca, n.epochs = 500, verbose = FALSE)
  data <- FindNeighbors(data, dims = 1:dims_pca, verbose = FALSE)
  data <- FindClusters(data, resolution = 1.0, verbose = FALSE)
  return(data)
}

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

annotate_cell_type_score <- function(data, cell_type_name, marker_genes, score_threshold = 0.5, use_dynamic_threshold = FALSE, std_mode = 'NA') {
  # Extract the expression data (use the appropriate assay, e.g., "RNA" or "SCT")
  expr_matrix <- GetAssayData(data, assay = "RNA", slot = "data") # Data is scalled and counts not
  
  # Ensure the marker genes are present in the expression data
  marker_genes <- marker_genes[marker_genes %in% rownames(expr_matrix)]
  
  # Calculate the mean expression of marker genes for each cluster
  mean_expression <- sapply(marker_genes, function(gene) {
    # Get expression values for the specific gene
    gene_expression <- expr_matrix[gene, ]
    
    # Group by clusters and calculate mean expression per cluster
    avg_expression <- tapply(gene_expression, data$seurat_clusters, mean)
    
    return(avg_expression)
  })
  
  # Calculate the average expression across all marker genes per cluster
  avg_expression_per_cluster <- rowMeans(mean_expression, na.rm = TRUE)
  
  # Print the average expression per cluster
  print(avg_expression_per_cluster)
  
  # Calculate overall mean and standard deviation of cluster expressions
  exp_mean_all <- mean(avg_expression_per_cluster)
  exp_sd_all <- sd(avg_expression_per_cluster)
  
  print(paste("Mean expression across clusters:", exp_mean_all))
  print(paste("Standard deviation across clusters:", exp_sd_all))
  
  # Initialize the dynamic threshold variable
  dynamic_threshold <- NA
  
  # Calculate the dynamic threshold based on std_mode
  if (use_dynamic_threshold) {
    if (std_mode == "below") {
      dynamic_threshold <- exp_mean_all - (exp_sd_all / 4)
      print(paste("Dynamic threshold (mean - sd/4):", dynamic_threshold))
    } else if (std_mode == "above") {
      dynamic_threshold <- exp_mean_all + (exp_sd_all / 4)
      print(paste("Dynamic threshold (mean + sd/4):", dynamic_threshold))
    } else if (std_mode == "NA") {
      dynamic_threshold <- exp_mean_all
      print(paste("Dynamic threshold (mean):", dynamic_threshold))
    }
  }
  
  # Determine which threshold to apply
  if (use_dynamic_threshold) {
    threshold_to_use <- dynamic_threshold
    print("Using dynamic threshold.")
  } else {
    threshold_to_use <- score_threshold
    print(paste("Using user-defined threshold:", score_threshold))
  }
  
  # Identify clusters with high average expression (above the selected threshold)
  target_clusters <- names(avg_expression_per_cluster[avg_expression_per_cluster > threshold_to_use])
  
  # Initialize the CellType column if it doesn't exist yet
  if (!"CellType" %in% colnames(data@meta.data)) {
    data$CellType <- NA
  }
  
  # Annotate the target clusters with the cell type name
  data$CellType[data$seurat_clusters %in% target_clusters] <- cell_type_name
  
  return(data)
}
# annotate_cell_type_score Another with ucell scores...

options(future.globals.maxSize = 6 * 1024^3) # Set to 6 GB (6 * 1024^3 bytes)
#plan(multisession, globals = FALSE) # No globals

data <- process_sct_rna(data)


#excel_data <- read_excel('cell_potency.xlsx', sheet = 'Sheet1')
#data@meta.data$CellPotency <- excel_data$cell_potency
saveRDS(data, file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))


data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))
plot <- DimPlot(object = data, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1, raster=FALSE) + 
        NoLegend()
plot


DimPlot(object = data, reduction = "umap", group.by = c("BatchID"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1, raster=FALSE)

# Combine all marker gene sets into a list for DotPlot
cell_type_markers <- c(
  "Ms4a1", "Cd19",                           # B cells
  "Ighg1", "Mzb1",                           # Plasma B cells
  "Cd3e", "Cd4", "Cd8a", "Trbc2",            # T cells
  'Ncr1','Klrb1b',                           # NK cells
  "Ms4a2", "Cpa3",                           # Mast cells
  "Gata3",                                   # Innate Lymphoid
  "Cd163", "Csf1r", "Cd68", "Mrc1",  "Aif1", # Macrophages
  "Col1a1", "Col3a1", "Col6a2",              # Fibroblasts
  "Vip", "Tubb3","Vat1l",                    # Enteric neurons
  "Sox10", "Plp1",                           # Enteric glial cells
  "Lyve1", "Flt4", "Pecam1",                 # Lymphatic endothelial
  "Plvap", "Flt1",                           # Vascular endothelial
  "Epcam", "Muc2", "Krt20",                  # Colonocytes
  "Myh11", "Tagln"                           # SMCs
)

# Remove duplicates from the gene list
cell_type_markers <- unique(cell_type_markers)

# Step 1: Ensure row names are unique
rownames(data) <- make.unique(rownames(data), sep = "_")
if (anyDuplicated(rownames(data))) {
  print("Duplicates still exist after make.unique:")
  print(rownames(data)[duplicated(rownames(data))])
} else {
  print("Row names are unique.")
}

# Step 2: Check for missing genes
missing_genes <- setdiff(cell_type_markers, rownames(data))
if (length(missing_genes) > 0) {
  print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
  cell_type_markers <- setdiff(cell_type_markers, missing_genes)  # Remove missing genes
}

# Step 3: Set the identity of the Seurat object to CellType
# unique(data$SubCellType)
# 
# new_subcell_types <- c(
#   "Aqp8+ Abs." = "Abs. Col",
#   "Car1+ Abs." = "Abs. Col",
#   "Cck+ EECs" = "EECs",
#   "EC" = "EECs",
#   "Imm. Goblet" = "Goblet"
# )
# # Rename the cell types in the 'CellType' column
# data$SubCellType <- recode(data$SubCellType, !!!new_subcell_types)
# 
# saveRDS(data, file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))

level_order <- c(
  "B cells", "Plasma B cells", "T cells", "NK cells", "Mast cells", "ILCs",
  "Macrophages", "Fibroblasts", "ENs", "EGCs",
  "LECs", "VECs", "Colonocytes", "SMCs"
)

# Ensure 'CellType' is treated as a factor and reordered
data$CellType <- factor(data$CellType, levels = level_order)

# Create the DotPlot
dot_plot <- DotPlot(object = data, dot.min = 0.05,
                    features = cell_type_markers,
                    group.by = "CellType") +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 13),
    axis.text.y = element_text(size = 15),
    strip.text = element_blank()
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(2, 8))

# Display the plot
print(dot_plot)

# Save the DotPlot
ggsave("plots/CellAtlas_dotplot.svg", plot = dot_plot, height = 6, width = 12, device = "svg")
dot_plot
dev.off()
png('plots/CellAtlas_dotplot.png', height = 6, width = 12, units = "in", res = 300)
dot_plot
dev.off()

# UMAP for BatchID
um_batch <- DimPlot(object = data, reduction = "umap", group.by = "BatchID", 
                   alpha = 1, raster=TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    legend.text = element_text(size = 15), # Increase legend text size (adjust 12 as needed)
    axis.title = element_text(size = 17, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

um_batch
# Save BatchID UMAP
svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_batch)
dev.off()


# UMAP for CellType with larger legend font
um_celltype <- DimPlot(object = data, reduction = "umap", group.by = "CellType",
                       label = TRUE, repel = TRUE, label.size = 5.5, 
                       label.box = TRUE, alpha = 1, raster = TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
    legend.text = element_text(size = 16),  # Increase legend label size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) + NoLegend() +
  coord_cartesian(clip = "off")

um_celltype
# Save CellType UMAP
svg('plots/UMAP_CellType.svg', height = 6, width = 6)
print(um_celltype)
dev.off()

# UMAP for CellType with larger legend font and no title
um_celltype <- DimPlot(object = data, reduction = "umap", group.by = "CellType",
                       alpha = 1, raster = TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
    legend.text = element_text(size = 16),  # Increase legend label size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) + NoLegend() + NoAxes() +
  coord_cartesian(clip = "off")

um_celltype

# Save CellType UMAP
svg('plots/UMAP_CellType_nonames.svg', height = 6, width = 6)
print(um_celltype)
dev.off()

# UMAP for condition
um_batch <- DimPlot(object = data, reduction = "umap", group.by = "Condition",
                    label = TRUE, repel = TRUE, label.size = 5.5, 
                    label.box = TRUE, alpha = 1, raster = TRUE) +
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

# Save Condition UMAP
svg('plots/UMAP_condition.svg', height = 6, width = 6)
print(um_batch)
dev.off()


# Create Feature Plot for Nr4a1 gene, split by Condition
um_feature_condition_list <- FeaturePlot(object = data, features = "Nr4a1",
                                         split.by = "Condition", by.col = F,
                                         cols = c("lightgray", "blue"),
                                         pt.size = 0.3,
                                         alpha = 0.8,
                                         order = F,
                                         ncol = 1,
                                         combine = FALSE) # Important: Set combine = FALSE
# Modify the title size for each individual plot in the list
modified_plots <- lapply(X = um_feature_condition_list,
                         FUN = function(x) x +
                           ggtitle("Nr4a1 gene expression") + # Set title to the condition
                           theme(plot.title = element_text(size = 16, hjust = 0.5)))
# Combine the plots using patchwork
combined_plot <- wrap_plots(modified_plots, ncol = 1) # Arrange in one column as in your original

combined_plot
# Save the feature plot with annotations
svg('plots/FeaturePlot_Nr4a1.svg', height = 6, width = 4)
print(combined_plot)
dev.off()
png('plots/FeaturePlot_Nr4a1.png', height = 6, width = 4, units = "in", res = 300)
print(combined_plot)
dev.off()

# Cell proportions and nr4a1 expression
dev.new()
# Apply custom theme
cust_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 20), # Center and enlarge title
  axis.text.x = element_text(angle = 60, hjust = 1, size = 18),
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 18),
  strip.text = element_text(size = 16)
)

# Create DittoSeq plot
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'BatchID', scale = 'count') +
  cust_theme + ggtitle("Cell counts per \n replicate")  # Set title correctly

cc
# Save plot
svg('plots/cell_counts_replicate.svg', height = 5.5, width = 6.5)
print(cc)
dev.off()

cc <- dittoBarPlot(data, var = 'CellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("Cell counts per \n condition")  # Set title correctly
cc
# Save plot
svg('plots/cell_counts_cond.svg', height = 5.5, width = 6.5)
print(cc)
dev.off()


# Create DittoSeq plot for cell proportions per replicate
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'BatchID', scale = 'percent') +
  cust_theme + ggtitle("Cell proportions per \n replicate (%)")  # Adjust title

cc
# Save plot
svg('plots/cell_proportions_replicate.svg', height = 5.5, width = 6.5)
print(cc)
dev.off()

# Create DittoSeq plot for cell proportions per condition
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'Condition', scale = 'percent') +
  cust_theme + ggtitle("Cell proportions per \n condition (%)")  # Adjust title

# Save plot
svg('plots/cell_proportions_cond.svg', height = 5.5, width = 6.5)
print(cc)
dev.off()



################################
##### Cell potency Plot section
################################

# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 16), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 13),
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 14)
)

# Plotting differentiation potency
# Extract the base score name for the title and filename
score_name <- "CellPotency"
clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name)
score_col <- score_name

# Create the plot
p <- data[[c(score_col, 'Condition', 'CellType')]] %>%
  rename(score = !!sym(score_col)) %>% # Rename the score column for easier use
  ggplot(aes(x = Condition, y = score, fill = Condition)) +
  geom_violin() +
  geom_boxplot(color = 'black', width = 0.5, outlier.size = 0.3) +
  geom_signif(comparisons = list(c("WT", "Nr4a1 KO")),
              map_signif_level = TRUE, step_increase = 0.12, vjust = 2) +
  ggtitle("Differentiation Potency") + ylab('Score') +
  facet_wrap(~CellType,  ncol = 7,  scales = 'free') + cust_theme + NoLegend()

p
# Define the filename for the SVG
filename <- paste0('plots/', gsub("[^[:alnum:]]", "_", score_name), '_Condition_all_celltypes.svg')

# Save the plot as an SVG
#svg(filename, height = 4, width = 10) #5 cols
svg(filename, height = 4, width = 14) #4 cols

print(p) # Important: Use print() to draw the ggplot object to the device
dev.off()

cat("Finished saving SVG plots for each score.\n")


# Load necessary libraries (if not already loaded)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer, pivot_wider
library(ggpubr) # For geom_signif (though not used directly in this specific bar plot style, good to have if needed for error bars)
library(Hmisc) # For mean_cl_normal (if error bars were to be added to this specific bar plot)
library(writexl)

# --- Define the score name ---
score_name <- "CellPotency"
clean_score_name <- gsub("Yuan et al \\(2019\\) ", "", score_name) # Adjust if your score name doesn't have this prefix

# --- 1. Extract and prepare data for CellPotency score ---
# Ensure these columns exist in your data_colo@meta.data
if (!all(c(score_name, 'Condition', 'CellType') %in% colnames(data@meta.data))) {
  stop("One or more required metadata columns (CellPotency, Condition, SubCellType) not found in data_colo.")
}

# Create a data frame from Seurat object metadata
df_cell_potency <- data[[c(score_name, 'Condition', 'SubCellType')]] %>%
  rename(score = !!sym(score_name)) # Rename the score column for easier use

# Ensure Condition is a factor with desired levels
df_cell_potency$Condition <- factor(df_cell_potency$Condition, levels = c("Nr4a1 KO", "WT"))
# Ensure SubCellType is a factor (order can be set later if desired)
df_cell_potency$CellType <- factor(df_cell_potency$SubCellType)


# --- 2. Perform Wilcoxon test and summarize data for CellPotency ---
results_table <- df_cell_potency %>%
  group_by(CellType) %>% # Group by SubCellType, equivalent to CellType in your example
  summarise(
    # Get statistics for KO
    N_KO = sum(Condition == "Nr4a1 KO"),
    Mean_KO = mean(score[Condition == "Nr4a1 KO"], na.rm = TRUE),
    # Get statistics for WT
    N_WT = sum(Condition == "WT"),
    Mean_WT = mean(score[Condition == "WT"], na.rm = TRUE),
    # Perform Wilcoxon test
    p_value = {
      ko_values <- score[Condition == "Nr4a1 KO"]
      wt_values <- score[Condition == "WT"]
      
      # Check if both KO and WT conditions exist and have enough data (at least 2 observations for wilcox.test)
      if (length(ko_values) >= 2 && length(wt_values) >= 2) {
        wilcox_result <- tryCatch(
          wilcox.test(ko_values, wt_values, exact = FALSE)$p.value,
          error = function(e) NA # Handle cases where test still can't be performed (e.g., all values identical)
        )
        wilcox_result
      } else {
        NA # Not enough data for the test (e.g., <2 cells in a group)
      }
    },
    .groups = 'drop' # Drop grouping after summarise
  ) %>%
  ungroup() # Ungroup the data frame

# Display the initial results table
cat("\n--- CellPotency Score Wilcoxon Test Results by SubCellType ---\n")
print(results_table)


# --- 3. Format the p-value for better readability and align with geom_signif ---
results_table_formatted <- results_table %>%
  mutate(p_value_formatted = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    is.na(p_value) ~ "N/A", # For cases where test couldn't be run
    TRUE ~ "NS"
  ))

cat("\n--- Formatted CellPotency Score Wilcoxon Test Results ---\n")
print(results_table_formatted)
write_xlsx(results_table_formatted, "cell_potency_score_stats.xlsx")

# --- If you want to recreate final_summary_table with the correct p-values (as in your example): ---
final_summary_table <- results_table_formatted %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted) %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Mean_"),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N, Mean),
    names_glue = "{.value}_{Condition}"
  ) %>%
  select(SubCellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted)

cat("\n--- Final Summary Table for CellPotency Score ---\n")
print(final_summary_table)






# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 20), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
  axis.title.x = element_blank(),   # element_text(size = 13) 
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 14), 
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 18)
)

# Create and save plot for Nr4a1 gene expression in conditions
df <- as.data.frame(GetAssayData(data)['Nr4a1',])
colnames(df) <- 'Nr4a1_expression'
df$CellType <- data$CellType
df$Condition <- data$Condition

p1 <- df %>% ggplot(aes(y = Nr4a1_expression, x = Condition, fill = Condition))+
  geom_boxplot(width=0.5, outliers = F)+ 
  geom_violin()+ 
  geom_signif(map_signif_level = T, comparisons = list(c('Nr4a1 KO','WT')), vjust = 2,tip_length = 0)+
  facet_wrap(~ CellType, ncol = 7, scales = 'free') + cust_theme + NoLegend()
p1

# Save plot
svg('plots/nr4a1_exp_cond.svg', height = 5, width = 14)
print(p1)
dev.off()


############################ Plots with wilcoxon test for Nr4a1 gene
# Define test batches
test_batches_ko <- c("NR4A1-KO-372", "NR4A1-KO-374", "NR4A1-KO-376", "NR4A1-KO-380")
test_batches_wt <- c("WT-1332", "WT-1353", "WT-1355", "WT-1362")

# Define the list of cell types
cell_types <- unique(data$CellType)

# Define the gene of interest
gl <- "Nr4a1"

# Loop over all cell types for NR4A1-KO test batches
for (cell_type in cell_types) {
  plot_ko <- plot_wilcox_violin_tests(data, cell_type = cell_type, gl = gl,
                                      test_batches = test_batches_ko,
                                      score_name = "Nr4a1 expression KO",
                                      plot_type = "both")
  #png(paste0("plots/", cell_type, "_NR4A1_KO.png"), height = 6, width = 8, res = 300, units = 'in')
  svg(paste0("plots/", cell_type, "_NR4A1_KO.svg"), height = 6, width = 8)
  print(plot_ko)
  dev.off()
}

# Loop over all cell types for WT test batches
for (cell_type in cell_types) {
  plot_wt <- plot_wilcox_violin_tests(data, cell_type = cell_type, gl = gl,
                                      test_batches = test_batches_wt,
                                      score_name = "Nr4a1 expression WT",
                                      plot_type = "both")
  #png(paste0("plots/", cell_type, "_WT.svg"), height = 6, width = 8, res = 300, units = 'in')
  svg(paste0("plots/", cell_type, "_WT.svg"), height = 6, width = 8)
  print(plot_wt)
  dev.off()
}






###########################################
####### Violin plots per cell type of Nr4a1
###########################################
library(dplyr)
library(tidyr) # For pivot_wider if you want a specific output format
library(purrr) # For map and potentially map_dfr

# Apply custom theme
cust_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = 0.5, size = 20), # Center and enlarge title
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
  axis.title.x = element_blank(),   # element_text(size = 13) 
  axis.text.y = element_text(size = 15),
  axis.title.y = element_text(size = 16),
  strip.text = element_text(size = 14), 
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 18)
)

DefaultAssay(data) <- "RNA"
#DefaultAssay(data) <- "RNA"
gene <- "Nr4a1"
gexp <- data@assays$RNA$data[rownames(data)==gene,]
df <- as.data.frame(gexp)
colnames(df) <- 'Nr4a1_expression'
df$CellType <- data$CellType
df$Condition <- data$Condition
df$Batches <- data$BatchID # Assuming this column exists in 'data'
df$Condition <- factor(df$Condition, levels = c("Nr4a1 KO", "WT"))

p1 <- df %>% ggplot(aes(y = Nr4a1_expression, x = Condition, fill = Condition))+
  geom_boxplot(width=0.5, outliers = F)+
  geom_violin(trim = TRUE)+ 
  geom_signif(map_signif_level = T, comparisons = list(c('Nr4a1 KO','WT')), vjust = 2,tip_length = 0, 
              test = "wilcox.test", test.args = list(exact = T)) +
  facet_wrap(~ CellType, ncol = 7, scales = 'free') + cust_theme + NoLegend()
p1

# Save plot
svg('plots/nr4a1_exp_cond_rna_assay.svg', height = 5, width = 14)
print(p1)
dev.off()

###########################################
####### Bar plots per cell type of Nr4a1
###########################################
# --- Perform Wilcoxon test and summarize data ---
results_table <- df %>%
  group_by(CellType) %>%
  summarise(
    # Get statistics for KO
    N_KO = sum(Condition == "Nr4a1 KO"),
    Mean_KO = mean(Nr4a1_expression[Condition == "Nr4a1 KO"], na.rm = TRUE),
    # Get statistics for WT
    N_WT = sum(Condition == "WT"),
    Mean_WT = mean(Nr4a1_expression[Condition == "WT"], na.rm = TRUE),
    # Perform Wilcoxon test
    p_value = {
      # Within a grouped summarise, '.' refers to the current group's data.
      # You can directly use the columns as they are already filtered by group_by.
      ko_values <- Nr4a1_expression[Condition == "Nr4a1 KO"]
      wt_values <- Nr4a1_expression[Condition == "WT"]
      
      # Check if both KO and WT conditions exist and have enough data
      if (length(ko_values) >= 2 && length(wt_values) >= 2) {
        wilcox_result <- tryCatch(
          wilcox.test(ko_values, wt_values, exact = FALSE)$p.value,
          error = function(e) NA # Handle cases where test still can't be performed (e.g., all values identical)
        )
        wilcox_result
      } else {
        NA # Not enough data for the test (e.g., <2 cells in a group)
      }
    }
  ) %>%
  ungroup()

# Display the results table
print(results_table)

# You can also format the p-value for better readability and align with geom_signif
results_table_formatted <- results_table %>%
  mutate(p_value_formatted = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    is.na(p_value) ~ "N/A", # For cases where test couldn't be run
    TRUE ~ "NS"
  ))

print(results_table_formatted)


# If you want to recreate final_summary_table with the correct p-values:
final_summary_table <- results_table_formatted %>%
  select(CellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted) %>%
  pivot_longer(
    cols = starts_with("N_") | starts_with("Mean_"),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)"
  ) %>%
  pivot_wider(
    names_from = Condition,
    values_from = c(N, Mean),
    names_glue = "{.value}_{Condition}"
  ) %>%
  select(CellType, N_KO, N_WT, Mean_KO, Mean_WT, p_value, p_value_formatted)

print(final_summary_table)



library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer

# --- MODIFICATION START ---
# --- MODIFICATION START (from previous response) ---
# Get the desired order of CellTypes from your final_summary_table
desired_celltype_order <- final_summary_table$CellType

# 1. Transform the data into a long format for plotting, carrying N values
plot_data <- final_summary_table %>%
  # Ensure CellType is a factor with the desired order *before* pivoting
  mutate(
    CellType = factor(CellType, levels = desired_celltype_order)
  ) %>%
  pivot_longer(
    cols = c(Mean_KO, Mean_WT, N_KO, N_WT),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)",
    values_drop_na = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("KO", "WT"))
  ) %>%
  rename(Mean_Expression = Mean, N_Cells = N) %>%
  select(CellType, Condition, Mean_Expression, N_Cells)

# Re-prepare data for combined facet labels AND order them
facet_annotation_data <- final_summary_table %>%
  mutate(
    # Ensure CellType is a factor with the desired order here as well
    CellType = factor(CellType, levels = desired_celltype_order),
    Facet_Label = paste0(CellType, "\nW. p-val: ", p_value_formatted),
    max_mean_expression = pmax(Mean_KO, Mean_WT, na.rm = TRUE)
  ) %>%
  # NEW: Convert Facet_Label to a factor with levels ordered by CellType
  mutate(
    Facet_Label = factor(Facet_Label, levels = unique(Facet_Label[order(CellType)]))
  ) %>%
  select(CellType, Facet_Label, max_mean_expression, p_value_formatted)


# Now, crucial for geom_text for p-value stars:
# Make sure plot_data also has the Facet_Label with the correct factor order,
# so that geom_text can inherit it correctly for faceting.
plot_data <- plot_data %>%
  left_join(select(facet_annotation_data, CellType, Facet_Label), by = "CellType") %>%
  # Ensure Facet_Label in plot_data is also a factor with the correct order
  mutate(Facet_Label = factor(Facet_Label, levels = levels(facet_annotation_data$Facet_Label)))

# --- MODIFICATION END ---


# 2. Create the bar plots with mean values and significance stars
p_bar_annotated <- plot_data %>%
  ggplot(aes(x = Condition, y = Mean_Expression, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Facet_Label, ncol = 7, scales = 'free_y') + # Facet by Facet_Label, which now has correctly ordered factor levels
  labs(
    y = "Nr4a1 Mean Expression",
    x = "Condition"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold", margin = margin(b = 10, t = 5), size = 14), # Increased from 10 to 14
    strip.background = element_rect(fill = "white", color = NA),
    
    panel.spacing.y = unit(0.5, "lines"),
    panel.spacing.x = unit(1, "lines"),
    
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    
    axis.title.y = element_text(margin = margin(r = 2), size = 16),
    # Removed the problematic duplicate and incorrect line:

    # Increase font size for y-axis tick labels
    axis.text.y = element_text(size = 14), # Added and increased size
    # Increase font size for x-axis tick labels (KO, WT)
    axis.text.x =  element_blank(), # This is the correct and desired line for x-axis text
    
    legend.position = "top",
    legend.box.margin = margin(b = 2),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  scale_fill_manual(values = c("KO" = scales::alpha("red", 0.6),
                               "WT" = scales::alpha("blue", 0.6))) +
  # Increase y range using expansion to ensure labels are visible
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + # Adjusted to 0.25 or 0.3 if needed
  # Add Mean Value above each bar (THIS WAS KEPT)
  geom_text(
    aes(
      label = sprintf("%.3f", Mean_Expression)
    ),
    vjust = -0.3, # Adjust vertical justification for mean label
    position = position_dodge(width = 0.8),
    size = 5,
    color = "black" # Ensure mean text is visible
  ) 
  # # --- ADDED BACK: geom_text FOR CELL COUNTS (N_Cells) ---
  # geom_text(
  #   aes(
  #     label = paste0("N=", N_Cells), # Format as "N=X"
  #     y = Mean_Expression # Start N at the top of the bar where Mean_Expression is
  #   ),
  #   vjust = 1.5, # Adjust vertical justification relative to Mean_Expression: 1.5 moves it below Mean_Expression
  #   position = position_dodge(width = 0.8),
  #   size = 4, # Slightly smaller size for count
  #   color = "black" # Choose a color that stands out but isn't dominant
  # ) +
  

print(p_bar_annotated)

# Save plot
svg('plots/nr4a1_exp_cond_hist_sct_assay.svg', height = 5, width = 14)
print(p_bar_annotated)
dev.off()


###########################################
####### Bar plots  type of Nr4a1
###########################################
# --- Perform Wilcoxon test and summarize data ---
results_table <- df %>%
  summarise(
    # Get statistics for KO
    N_KO = sum(Condition == "Nr4a1 KO"),
    Mean_KO = mean(Nr4a1_expression[Condition == "Nr4a1 KO"], na.rm = TRUE),
    # Get statistics for WT
    N_WT = sum(Condition == "WT"),
    Mean_WT = mean(Nr4a1_expression[Condition == "WT"], na.rm = TRUE),
    # Perform Wilcoxon test
    p_value = {
      ko_values <- Nr4a1_expression[Condition == "Nr4a1 KO"]
      wt_values <- Nr4a1_expression[Condition == "WT"]
      
      if (length(ko_values) >= 2 && length(wt_values) >= 2) {
        wilcox_result <- tryCatch(
          wilcox.test(ko_values, wt_values, exact = FALSE)$p.value,
          error = function(e) NA # Handle cases where test still can't be performed
        )
        wilcox_result
      } else {
        NA # Not enough data for the test
      }
    }
  ) %>%
  ungroup()

# Format the p-value for better readability and align with common significance stars
results_table_formatted <- results_table %>%
  mutate(p_value_formatted = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    is.na(p_value) ~ "N/A", # For cases where test couldn't be run
    TRUE ~ "NS"
  ))

print("Global Wilcoxon Test Results:")
print(results_table_formatted)

# Prepare data for plotting
plot_data <- results_table_formatted %>%
  pivot_longer(
    cols = c(Mean_KO, Mean_WT, N_KO, N_WT),
    names_to = c(".value", "Condition"),
    names_pattern = "([A-Za-z]+)_(KO|WT)",
    values_drop_na = FALSE
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("KO", "WT")) # Use full labels for consistency
  ) %>%
  rename(Mean_Expression = Mean, N_Cells = N) %>%
  select(Condition, Mean_Expression, N_Cells)

# Prepare data for the global p-value annotation
max_global_mean <- max(plot_data$Mean_Expression, na.rm = TRUE)
# Calculate a Y-position for the annotation, providing some buffer above the highest bar
annotation_y_pos <- max_global_mean * 1.15 

annotation_df <- data.frame(
  x = 1.5, # X-position for the annotation (centered between KO and WT bars)
  y = annotation_y_pos,
  label = results_table_formatted$p_value_formatted[1] # Get the single global p-value formatted string
)

# 2. Create the global bar plot
p_global_bar <- plot_data %>%
  ggplot(aes(x = Condition, y = Mean_Expression, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(
    title = "Global Nr4a1 Expression",
    y = "Nr4a1 Mean Expression",
    x = "Condition"
  ) +
  theme_minimal() +
  theme(
    # Remove strip elements as there are no facets
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    panel.spacing.x = unit(1, "lines"),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
    axis.title.y = element_text(margin = margin(r = 10), size = 16), # Increased right margin
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14), # Show x-axis labels (KO, WT)
    legend.position = "none", # Hide legend as fill is self-explanatory
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold") # Center and style title
  ) +
  scale_fill_manual(values = c("KO" = scales::alpha("red", 0.6),
                               "WT" = scales::alpha("blue", 0.6))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) + # Adjust Y-axis range
  
  # Add Mean Value above each bar
  geom_text(
    aes(
      label = sprintf("%.3f", Mean_Expression)
    ),
    vjust = -0.5, # Adjust vertical justification for mean label
    position = position_dodge(width = 0.8),
    size = 4,
    color = "black"
  ) +
  # Add Cell Counts (N_Cells) below the Mean Value
  geom_text(
    aes(
      label = paste0("N=", N_Cells),
      y = Mean_Expression # Position relative to the top of the bar
    ),
    vjust = 1.5, # Adjust vertical justification to place below mean
    position = position_dodge(width = 0.8),
    size = 4,
    color = "black"
  ) +
  # Add the global p-value annotation
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE, # Crucial: prevents inheriting aesthetics from main plot_data
    size = 6, # Larger text for global p-value
    color = "black",
    vjust = 0 # Align text to the top of its y-position
  )


print(p_global_bar)

# Save plot
svg('plots/global_nr4a1_expression_barplot.svg', height = 4, width = 5) # Adjusted dimensions for a single plot
print(p_global_bar)
dev.off()



###############################
############# Cell numbers 
###############################
# 1. Capture the cross-tabulation table
cell_type_batch_counts <- table(data$CellType, data$BatchID)

# 2. Convert to a data frame (long format) AND RENAME COLUMNS IMMEDIATELY
df_cell_type_batch_long <- as.data.frame(cell_type_batch_counts)
names(df_cell_type_batch_long) <- c("CellType", "BatchID", "Count")

# 3. Calculate percentages per cell type (Percentage of each CellType's total that comes from each batch)
df_percentage_per_celltype <- df_cell_type_batch_long %>%
  group_by(CellType) %>%
  mutate(
    TotalCellsInCellType = sum(Count),
    Percentage = (Count / TotalCellsInCellType) * 100
  ) %>%
  ungroup()

print("df_percentage_per_celltype (long format, % of CellType's total):")
print(head(df_percentage_per_celltype))


# If you want to see it in a wide format (BatchID as columns, percentages as values)
df_percentage_per_celltype_wide <- df_percentage_per_celltype %>%
  select(CellType, BatchID, Percentage) %>%
  pivot_wider(
    names_from = BatchID,
    values_from = Percentage,
    names_prefix = "Percentage_"
  )

print("\ndf_percentage_per_celltype_wide (wide format, % of CellType's total):")
print(head(df_percentage_per_celltype_wide))

# --- YOUR ORIGINAL GROUPED BAR PLOT (THIS WILL WORK) ---
cat("\n--- Grouped Bar Plot (Each CellType Facet Shows % of its total cells per Batch) ---\n")
p_grouped_bar_faceted <- ggplot(df_percentage_per_celltype, aes(x = BatchID, y = Percentage, fill = BatchID)) +
  geom_col(position = "dodge") + # This creates grouped bars
  facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
  labs(
    title = "Distribution of Each Cell Type's Total Cells Across Batches", # Clarified title
    y = "Percentage of Cell Type's Total Cells", # Clarified y-axis label
    x = "Batch ID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

print(p_grouped_bar_faceted)


# --- NOW FOR THE STACKED BAR PLOTS (Cell Type Composition within Each Batch) ---
# Add Condition column for faceting
# Recalculate percentages: Percentage of each CellType *within each Batch*
df_composition_per_batch <- df_cell_type_batch_long %>%
  group_by(BatchID) %>%
  mutate(
    TotalCellsInBatch = sum(Count),
    Percentage_of_Batch = (Count / TotalCellsInBatch) * 100
  ) %>%
  ungroup()

# Add Condition column for faceting (THIS IS CRUCIAL FOR THE FACETED PLOT)
df_composition_per_batch <- df_composition_per_batch %>%
  mutate(Condition = ifelse(grepl("^KO", BatchID), "KO", "WT")) %>%
  mutate(Condition = factor(Condition, levels = c("KO", "WT"))) # Ensure order

# --- Plot: Batches Stacked Up Per Cell Type ---
# Each facet is a CellType. Each bar (one per facet) shows how Batches contribute to that CellType's total.
cat("\n--- Stacked Bar Plot: Batch Contribution to Each Cell Type's Total Cells ---\n")
p_batches_stacked_per_celltype <- ggplot(df_percentage_per_celltype, aes(x = "", y = Percentage, fill = BatchID)) +
  geom_col(position = "stack") + # Use "stack" for stacked bars
  facet_wrap(~ CellType, ncol = 5) + # Facet by CellType
  labs(
    y = "Percentage of Cell Type's Total Cells (across all batches)",
    x = "", # No x-axis label needed as there's only one bar per facet
    fill = "Batch ID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Hide x-axis text
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)), # Increased size and added margin for better spacing
    legend.position = "right",
    legend.title = element_text(size = 14), # Make legend title larger and bold
    legend.text = element_text(size = 12),          
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 

print(p_batches_stacked_per_celltype)

# Save plot
svg('plots/cell_proportions_stacked_per_celltype.svg', height = 6, width = 10)
print(p_batches_stacked_per_celltype)
dev.off()




library(writexl)
library(dplyr) # Ensure dplyr is loaded for data manipulation
library(tidyr) # Ensure tidyr is loaded for data manipulation

# --- START TUNING FOR KO VS WT BATCHES ---
# Define KO and WT batch names
ko_batches <- c("NR4A1-KO-372", "NR4A1-KO-374", "NR4A1-KO-376", "NR4A1-KO-380")
wt_batches <- c("WT-1332", "WT-1353", "WT-1355", "WT-1362")

# --- A. Analysis for KO Batches ---

# 1. Subset raw counts for KO batches
ko_cell_type_batch_counts <- cell_type_batch_counts[, ko_batches]
print("KO Batch Raw Counts:")
print(ko_cell_type_batch_counts)

# 2. Run Chi-squared test for KO batches (Overall independence test)
cat("\n--- Chi-squared Test for KO Batches (Overall Cell Type Distribution) ---\n")
chi_sq_ko_result <- chisq.test(ko_cell_type_batch_counts, simulate.p.value = T, B = 10000)
print(chi_sq_ko_result)

# New Section: Test for significant differences in abundance for each cell type across KO batches
cat("\n--- Fisher's Exact Test for Each Cell Type Abundance in KO Batches ---\n")

# Calculate total cells per KO batch (needed for 'other cells' counts)
total_ko_batch_counts <- colSums(ko_cell_type_batch_counts)

# Initialize a list to store p-values for each cell type
ko_p_values_list <- list()

# Loop through each cell type
for (cell_type_name in rownames(ko_cell_type_batch_counts)) {
  # Get counts for the specific cell type across all KO batches
  specific_cell_counts <- ko_cell_type_batch_counts[cell_type_name, ]
  
  # Calculate counts for 'all other cells' in each KO batch
  other_cells_counts <- total_ko_batch_counts - specific_cell_counts
  
  # Create the 2xN contingency table for Fisher's Exact Test
  # Rows: [Specific Cell Type], [All Other Cell Types]
  # Columns: Each KO Batch
  contingency_table <- rbind(specific_cell_counts, other_cells_counts)
  rownames(contingency_table) <- c(cell_type_name, "Other_cells")
  
  # Run Fisher's Exact Test with Monte Carlo simulation for p-value
  # This is robust for tables larger than 2x2 and with low counts
  test_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
  ko_p_values_list[[cell_type_name]] <- test_result$p.value
}

# Convert results to a data frame
ko_p_values_df <- data.frame(
  CellType = names(ko_p_values_list),
  P_Value = unlist(ko_p_values_list)
)
# Then access CellType via rownames(ko_p_values_df)

# Apply False Discovery Rate (FDR) adjustment for multiple comparisons
ko_p_values_df$Adjusted_P_Value_FDR <- p.adjust(ko_p_values_df$P_Value, method = "BH") # Benjamini-Hochberg

# Print the results, sorted by adjusted p-value
print("P-values for individual cell type abundance differences in KO Batches (with FDR adjustment):")
print(ko_p_values_df[order(ko_p_values_df$Adjusted_P_Value_FDR), ])


# 3. Subset wide percentage table for KO batches
# If you want to see it in a wide format (BatchID as columns, percentages as values)
# df_percentage_per_celltype_wide <- df_percentage_per_celltype %>%
#   select(CellType, BatchID, Percentage) %>%
#   pivot_wider(
#     names_from = BatchID,
#     values_from = Percentage,
#     names_prefix = "Percentage_"
#   )
# print("\ndf_percentage_per_celltype_wide (wide format, % of CellType's total):")
# print(head(df_percentage_per_celltype_wide))

df_percentage_ko_wide <- df_percentage_per_celltype_wide %>%
  dplyr::select(CellType, starts_with("Percentage_NR4A1-KO-")) # Selects all KO percentage columns
print("\nKO Batch Percentages (Wide Format):")
print(df_percentage_ko_wide)

# 4. Calculate CV for cell types within KO batches
cat("\n--- Coefficient of Variation for KO Batches ---\n")
ko_cell_type_cv <- df_percentage_ko_wide %>%
  rowwise() %>%
  mutate(
    Mean_Percentage = mean(c_across(starts_with("Percentage_NR4A1-ko-")), na.rm = TRUE),
    SD_Percentage = sd(c_across(starts_with("Percentage_NR4A1-KO-")), na.rm = TRUE),
    CV_Percentage = (SD_Percentage / Mean_Percentage) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(CellType, Mean_Percentage, SD_Percentage, CV_Percentage)

print(ko_cell_type_cv)

# --- B. Analysis for WT Batches ---
# 1. Subset raw counts for WT batches
wt_cell_type_batch_counts <- cell_type_batch_counts[, wt_batches]
print("\nWT Batch Raw Counts:")
print(wt_cell_type_batch_counts)

# 2. Run Chi-squared test for WT batches (Overall independence test)
cat("\n--- Chi-squared Test for WT Batches (Overall Cell Type Distribution) ---\n")
chi_sq_wt_result <- chisq.test(wt_cell_type_batch_counts, simulate.p.value = T, B = 10000)
print(chi_sq_wt_result)

# New Section: Test for significant differences in abundance for each cell type across WT batches
cat("\n--- Fisher's Exact Test for Each Cell Type Abundance in WT Batches ---\n")

# Calculate total cells per WT batch (needed for 'other cells' counts)
total_wt_batch_counts <- colSums(wt_cell_type_batch_counts)

# Initialize a list to store p-values for each cell type
wt_p_values_list <- list()

# Loop through each cell type
for (cell_type_name in rownames(wt_cell_type_batch_counts)) {
  # Get counts for the specific cell type across all WT batches
  specific_cell_counts <- wt_cell_type_batch_counts[cell_type_name, ]
  
  # Calculate counts for 'all other cells' in each WT batch
  other_cells_counts <- total_wt_batch_counts - specific_cell_counts
  
  # Create the 2xN contingency table for Fisher's Exact Test
  # Rows: [Specific Cell Type], [All Other Cell Types]
  # Columns: Each WT Batch
  contingency_table <- rbind(specific_cell_counts, other_cells_counts)
  rownames(contingency_table) <- c(cell_type_name, "Other_cells")
  
  # Run Fisher's Exact Test with Monte Carlo simulation for p-value
  test_result <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
  wt_p_values_list[[cell_type_name]] <- test_result$p.value
}

# Convert results to a data frame
wt_p_values_df <- data.frame(
  CellType = names(wt_p_values_list),
  P_Value = unlist(wt_p_values_list)
)

# Apply False Discovery Rate (FDR) adjustment for multiple comparisons
wt_p_values_df$Adjusted_P_Value_FDR <- p.adjust(wt_p_values_df$P_Value, method = "BH") # Benjamini-Hochberg

# Print the results, sorted by adjusted p-value
print("P-values for individual cell type abundance differences in WT Batches (with FDR adjustment):")
print(wt_p_values_df[order(wt_p_values_df$Adjusted_P_Value_FDR), ])


# 3. Subset wide percentage table for WT batches
df_percentage_wt_wide <- df_percentage_per_celltype_wide %>%
  dplyr::select(CellType, starts_with("Percentage_WT-")) # Selects all WT percentage columns
print("\nWT Batch Percentages (Wide Format):")
print(df_percentage_wt_wide)

# 4. Calculate CV for cell types within WT batches
cat("\n--- Coefficient of Variation for WT Batches ---\n")
wt_cell_type_cv <- df_percentage_wt_wide %>%
  rowwise() %>%
  mutate(
    Mean_Percentage = mean(c_across(starts_with("Percentage_WT-")), na.rm = TRUE),
    SD_Percentage = sd(c_across(starts_with("Percentage_WT-")), na.rm = TRUE),
    CV_Percentage = (SD_Percentage / Mean_Percentage) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(CellType, Mean_Percentage, SD_Percentage, CV_Percentage)

print(wt_cell_type_cv)

# --- Writing all results to an Excel spreadsheet ---

# Prepare the list of data frames to be written
# Note: For chi-squared test results, you'll need to extract key stats into a data frame.
# For simplicity, let's create data frames for the p-value and relevant statistics.

# Extract Chi-squared test results into data frames
# For KO batches
chi_sq_ko_df <- data.frame(
  Statistic = chi_sq_ko_result$statistic,
  df = chi_sq_ko_result$parameter,
  P_Value = chi_sq_ko_result$p.value,
  Method = "Chi-squared Test (simulated p-value)",
  Group = "KO Batches"
)

# For WT batches
chi_sq_wt_df <- data.frame(
  Statistic = chi_sq_wt_result$statistic,
  df = chi_sq_wt_result$parameter,
  P_Value = chi_sq_wt_result$p.value,
  Method = "Chi-squared Test (simulated p-value)",
  Group = "WT Batches"
)


# List of data frames to write to Excel
output_list <- list(
  "KO_Raw_Counts" = as.data.frame(ko_cell_type_batch_counts), # Convert matrix to dataframe
  "KO_ChiSq_Overall" = chi_sq_ko_df,
  "KO_Fisher_Per_CellType" = ko_p_values_df[order(ko_p_values_df$Adjusted_P_Value_FDR), ],
  "KO_Percentages_Wide" = df_percentage_ko_wide,
  "KO_CV_Per_CellType" = ko_cell_type_cv,
  "WT_Raw_Counts" = as.data.frame(wt_cell_type_batch_counts), # Convert matrix to dataframe
  "WT_ChiSq_Overall" = chi_sq_wt_df,
  "WT_Fisher_Per_CellType" = wt_p_values_df[order(wt_p_values_df$Adjusted_P_Value_FDR), ],
  "WT_Percentages_Wide" = df_percentage_wt_wide,
  "WT_CV_Per_CellType" = wt_cell_type_cv
)

# Define the output file path
output_file <- "Cell_Proportion_Batch_Analysis_Results.xlsx"

# Write the list of data frames to an Excel file, each as a separate sheet
writexl::write_xlsx(output_list, path = output_file)

cat(paste0("\nAll results have been successfully written to: ", output_file, "\n"))
