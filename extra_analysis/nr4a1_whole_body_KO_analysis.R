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

base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "AllSamplesMerged_OldGenesRemoved_StringentFiltered_DoubletsRemoved_AmbientRNARemoved_AllGenesUMAP_FixedKOvsWT_Annotated_SampleSplit.rds"))

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

plot_wilcox_violin_tests <- function(data, cell_type, gl, test_batches = NULL, score_name = "Score", plot_type = "both") {
  
  # Subset metadata first to keep only relevant cell types
  cell_indices <- which(data$CellType == cell_type)
  
  if (length(cell_indices) == 0) {
    stop("No cells found for the specified cell type.")
  }
  
  # Subset assay data to relevant cells and genes
  assay_data <- GetAssayData(data, slot = 'data')[gl, cell_indices, drop = FALSE]
  
  # Convert to full matrix only for this subset
  if (inherits(assay_data, "dgCMatrix")) {
    assay_data_full <- as.matrix(assay_data)
  } else {
    assay_data_full <- assay_data
  }
  
  # Process gene expression values
  if (length(gl) == 1) {
    df <- as.data.frame(t(assay_data_full))
    colnames(df) <- paste0(score_name)
  } else {
    gene_data <- assay_data_full
    cell_scores <- rowMeans(gene_data, na.rm = TRUE)
    df <- data.frame(cell_scores)
    colnames(df) <- paste0(score_name)
  }
  
  # Clean up to free memory
  rm(assay_data_full)
  
  # Add metadata columns
  df$CellType <- factor(data$CellType[cell_indices])
  df$BatchID <- data$BatchID[cell_indices]
  df$condition <- data$Condition[cell_indices]
  
  # Remove rows with NA values
  df <- df[!is.na(df[[score_name]]) & !is.na(df$CellType), ]
  
  # If specific batches are provided, subset further
  if (!is.null(test_batches)) {
    df <- df[df$BatchID %in% test_batches, ]
    df$BatchID <- factor(df$BatchID)
  }
  
  # Debugging checks
  print("Subset data for cell type:")
  print(head(df))
  
  print("Distribution of BatchID:")
  print(table(df$BatchID))
  
  # Perform pairwise Wilcoxon test
  wilcox_test <- pairwise.wilcox.test(df[[score_name]], df$BatchID, p.adjust.method = "bonferroni")
  
  if (all(is.na(wilcox_test$p.value))) {
    stop("No valid comparisons found in the Wilcoxon test.")
  }
  
  # Convert Wilcoxon results into a long format data frame
  comparison_df <- as.data.frame(wilcox_test$p.value)
  comparison_df_long <- comparison_df %>%
    rownames_to_column("Group1") %>%
    pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_value")
  
  # Apply p-value labels
  p_value_labels <- comparison_df_long %>%
    mutate(label = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      P_value < 0.1   ~ ".",
      TRUE            ~ "ns"
    ))
  
  # Filter significant comparisons
  comparisons_list <- p_value_labels %>%
    filter(label != "ns") %>%
    select(Group1, Group2) %>%
    rowwise() %>%
    mutate(pairwise = list(c(as.character(Group1), as.character(Group2)))) %>%
    pull(pairwise)
  
  # Calculate max value for proper geom_signif placement
  y_max <- max(df[[score_name]], na.rm = TRUE)
  #y_max <- quantile(df[[score_name]], probs = 0.75, na.rm = TRUE)
  
  # Adjust significance bar height based on plot type
  y_offset <- if (plot_type == "violin") {
    0.05 * y_max  # Slightly above violin plot
  } else if (plot_type == "boxplot") {
    0.05 * y_max  # Higher above boxplot
  } else {
    0.05 * y_max # Balanced for both
  }
  
  # Create the plot
  plot <- ggplot(df, aes(x = BatchID, y = .data[[score_name]], fill = BatchID))
  
  # Add plots based on plot_type argument
  if (plot_type == "violin") {
    plot <- plot + geom_violin(alpha = 0.6, width = 1)
  } else if (plot_type == "boxplot") {
    plot <- plot + geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6, color = "black")
  } else if (plot_type == "both") {
    plot <- plot +
      geom_violin(alpha = 0.6, width = 1) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4, color = "black")
  }
  
  # Add significance annotations at adjusted y position
  plot <- plot +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    ggtitle(paste(score_name, "for", cell_type)) +
    geom_signif(
      comparisons = comparisons_list,
      map_signif_level = TRUE,
      step_increase = y_offset / y_max,  # Dynamically scale step increase
      y_position = y_max + y_offset      # Place significance bars above max data point
    )
  
  
  return(plot)
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


data <- process_rna(data, resolution = 1.0)
Idents(data) <- "CellType"

DimPlot(object = data, reduction = "umap", group.by = c("BatchID","CellType"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1)


saveRDS(data, file.path(base_dir, "AllSamplesMerged_OldGenesRemoved_StringentFiltered_DoubletsRemoved_AmbientRNARemoved_AllGenesUMAP_FixedKOvsWT_Annotated_SampleSplit.rds"))


data <- readRDS(file.path(base_dir, "AllSamplesMerged_OldGenesRemoved_StringentFiltered_DoubletsRemoved_AmbientRNARemoved_AllGenesUMAP_FixedKOvsWT_Annotated_SampleSplit.rds"))

data$Condition <- factor(ifelse(grepl("WT", data$BatchID), "WT", "Nr4a1 KO"))


# Combine all marker gene sets into a list for DotPlot
cell_type_markers <- c(
  "Ms4a1", "Cd19",                           # B cells
  "Ighg1", "Mzb1",                           # Plasma B cells
  "Cd3e", "Cd4", "Cd8a", "Trbc2",            # T cells
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
Idents(data) <- "CellType"

unique(data$CellType)

# Step 4: Define the level order you want for plotting
level_order <- c(
  "B cells", "Plasma B cells", "T cells", "Macrophages", "Fibroblasts",
  "Enteric neurons", "Enteric glial cells", "Lymphatic endothelial",
  "Vascular endothelial", "Colonocytes", "SMCs"
)

# Ensure 'CellType' is treated as a factor and reordered
data$CellType <- factor(data$CellType, levels = level_order)

# Create the DotPlot
dot_plot <- DotPlot(object = data, dot.min = 0.05,
                    features = cell_type_markers,
                    group.by = "CellType") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 16),
    strip.text = element_blank()
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(range = c(2, 8))

# Display the plot
print(dot_plot)


# Save the DotPlot
png('plots/CellAtlas_dotplot.png', height = 6, width = 12, res = 300, units = 'in')
#pdf('plots/CellAtlas_dotplot.pdf', height = 6, width = 12)
print(dot_plot)
dev.off()

# UMAP for BatchID
um_batch <- DimPlot(object = data, reduction = "umap", group.by = "BatchID", label = TRUE,
                    repel = TRUE, label.size = 5, label.box = TRUE, alpha = 1, raster=TRUE) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

# Save BatchID UMAP
png('plots/UMAP_BatchID.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_batch)
dev.off()

# UMAP for CellType with larger legend font
um_celltype <- DimPlot(object = data, reduction = "umap", group.by = "CellType",
                       label = TRUE, repel = TRUE, label.size = 5, label.box = TRUE, alpha = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 10, color = "black"),  # Keep axis titles
    legend.text = element_text(size = 16),  # Increase legend label size
    legend.title = element_text(size = 14)  # Increase legend title size
  ) +
  coord_cartesian(clip = "off")

# Save CellType UMAP
png('plots/UMAP_CellType.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_CellType.svg', height = 6, width = 8)
print(um_celltype)
dev.off()

# UMAP for condition
um_batch <- DimPlot(object = data, reduction = "umap", group.by = "Condition",
                    label = TRUE, repel = TRUE, label.size = 5, label.box = TRUE, alpha = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    legend.text = element_text(size = 16),  # Increase legend label size
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

# Save BatchID UMAP
png('plots/UMAP_condition.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_condition.svg', height = 6, width = 8)
print(um_batch)
dev.off()

# Create Feature Plot for Nr4a1 gene, split by Condition
um_feature_condition <- FeaturePlot(object = data, features = "Nr4a1",
                                    split.by = "Condition",
                                    cols = c("lightgray", "blue"),
                                    pt.size = 0.15,
                                    alpha = 0.5,
                                    label.size = 3,
                                    order = T
                                    )
# Save the feature plot with annotations
png('plots/FeaturePlot_Nr4a1.png', height = 6, width = 10, res = 300, units = 'in')
#svg('plots/FeaturePlot_Nr4a1.svg', height = 6, width = 10)
print(um_feature_condition)
dev.off()







# Cell proportions and nr4a1 expression
dev.new()
# Apply custom theme
cust_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 24), # Center and enlarge title
  axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 18),
  strip.text = element_text(size = 16)
)

# Create DittoSeq plot
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'BatchID', scale = 'count') +
  cust_theme + ggtitle("Cell counts per replicate")  # Set title correctly

# Save plot
png('plots/cell_counts_replicate.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()


cc <- dittoBarPlot(data, var = 'CellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("Cell counts per condition")  # Set title correctly

# Save plot
png('plots/cell_counts_cond.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()


# Create DittoSeq plot for cell proportions per replicate
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'BatchID', scale = 'percent') +
  cust_theme + ggtitle("Cell proportions per replicate (%)")  # Adjust title

# Save plot
png('plots/cell_proportions_replicate.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)
dev.off()

# Create DittoSeq plot for cell proportions per condition
cc <- dittoBarPlot(data, var = 'CellType', group.by = 'Condition', scale = 'percent') +
  cust_theme + ggtitle("Cell proportions per condition (%)")  # Adjust title

# Save plot
png('plots/cell_proportions_cond.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)
dev.off()


# Apply custom theme
cust_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 24), # Center and enlarge title
  axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 18),
  axis.title.y = element_text(size = 18),
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 18),
  strip.text = element_text(size = 16)
)

# Create and save plot for Nr4a1 gene expression in replicates with 3 panels per row
gexp <- dittoPlot(data, var = 'Nr4a1', plots = c('boxplot', 'vlnplot'), group.by = 'BatchID',
                  split.by = 'CellType', boxplot.show.outliers = F, boxplot.width = 0.8) +
  cust_theme +
  ggtitle("Nr4a1 gene expressions in replicates") +
  facet_wrap(~ CellType, ncol = 4)  # Set 3 panels per row

# Save plot
png('plots/nr4a1_exp_rep.png', height = 12, width = 12, res = 300, units = 'in')
print(gexp)
dev.off()

# Create and save plot for Nr4a1 gene expression in conditions
gexp_cond <- dittoPlot(data, var = 'Nr4a1', plots = c('boxplot', 'vlnplot'), group.by = 'Condition',
                       split.by = 'CellType', boxplot.show.outliers = F, boxplot.width = 0.8) + cust_theme +
  ggtitle("Nr4a1 gene expressions in conditions")  # Adjust title

# Save plot
png('plots/nr4a1_exp_cond.png', height = 12, width = 12, res = 300, units = 'in')
print(gexp_cond)
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
  png(paste0("plots/", cell_type, "_NR4A1_KO.png"), height = 6, width = 8, res = 300, units = 'in')
  print(plot_ko)
  dev.off()
}

# Loop over all cell types for WT test batches
for (cell_type in cell_types) {
  plot_wt <- plot_wilcox_violin_tests(data, cell_type = cell_type, gl = gl,
                                      test_batches = test_batches_wt,
                                      score_name = "Nr4a1 expression WT",
                                      plot_type = "both")
  png(paste0("plots/", cell_type, "_WT.png"), height = 6, width = 8, res = 300, units = 'in')
  print(plot_wt)
  dev.off()
}



library(openxlsx)
# Ensure `cell_types` is properly extracted
cell_types <- unique(data$CellType)
# Set the grouping variable for Seurat
Idents(data) <- "Condition"
# Define output directory
output_dir <- "DE_results"
# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cell in cell_types) {
  # Print the current cell type being processed
  cat("Processing cell type:", cell, "\n")

  # Subset the Seurat object for the current cell type
  cell_subset <- subset(data, CellType == cell)

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
  safe_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(safe_cell_name, "_markers_Wilcoxon_updated.xlsx"))

  # Save workbook
  saveWorkbook(wb, file = output_file, overwrite = TRUE)

  cat("Saved:", output_file, "\n")
}





library(enrichR)
library(openxlsx)
library(dplyr)

# Define paths
input_dir <- "DE_results"
output_dir <- "Enrichr_results"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
dbs <- "GO_Biological_Process_2023"
pval_cutoff <- 0.05
top_n <- 100

# Function to perform Enrichr and process results with waiting time
run_enrichr <- function(genes, direction, dbs, pval_cutoff, wait_time = 2) {
  results_df <- data.frame()

  for (db in dbs) {
    enrichr_res <- enrichr(genes, db)[[1]]

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
  Sys.sleep(wait_time)  # Prevent API overload
  return(results_df)
}

# Get unique cell types (assuming 'data' is defined elsewhere)
cell_types <- unique(data$CellType)

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

  # **Skip file writing if no significant pathways are found**
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

  # **Only save workbook if at least one pathway exists**
  saveWorkbook(wb_out, output_file, overwrite = TRUE)
  print(paste("Enrichment results saved in", output_file))
}









library(openxlsx)

# Ensure `cell_types` is properly extracted
cell_types <- unique(data$CellType)

# Define output directory for DV results
output_dir <- "DV_results"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (cell in cell_types) {
  # Print the current cell type being processed
  cat("Processing cell type:", cell, "\n")

  # Subset the Seurat object for the current cell type
  cell_subset <- subset(data, CellType == cell)

  # Further subset for KO and WT conditions
  data_ko <- subset(cell_subset, Condition == 'Nr4a1 KO')
  data_wt <- subset(cell_subset, Condition == 'WT')

  # Extract counts for KO and WT
  data_ko <- GetAssayData(data_ko, layer = 'counts')  # slot = 'counts' instead of layers
  data_wt <- GetAssayData(data_wt, layer = 'counts')

  # Compute Differential Variability (DV)
  DV_res <- data.frame(splineDV(X = data_ko, Y = data_wt))

  # Filter for significant results
  DV_res <- DV_res[DV_res$Pval < 0.05,]

  # Create a new workbook for DV results
  wb <- createWorkbook()

  # Add a worksheet to the workbook for the overall DV results
  addWorksheet(wb, "DV_Results")
  writeData(wb, sheet = "DV_Results", DV_res)

  # Create additional sheets for variability up and down
  DV_up <- DV_res %>% filter(Direction == 1)  # Variability up
  DV_down <- DV_res %>% filter(Direction == -1)  # Variability down

  # Add the Variability Up sheet if there are any upregulated genes
  if (nrow(DV_up) > 0) {
    addWorksheet(wb, "Variability_Up")
    writeData(wb, sheet = "Variability_Up", DV_up)
  } else {
    cat("No variability up genes for", cell, "\n")
  }

  # Add the Variability Down sheet if there are any downregulated genes
  if (nrow(DV_down) > 0) {
    addWorksheet(wb, "Variability_Down")
    writeData(wb, sheet = "Variability_Down", DV_down)
  } else {
    cat("No variability down genes for", cell, "\n")
  }

  # Ensure consistent file naming and save the workbook
  safe_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(safe_cell_name, "_DV_results.xlsx"))
  saveWorkbook(wb, file = output_file, overwrite = TRUE)

  cat("Saved DV results for:", cell, "\n")
}







library(enrichR)
library(openxlsx)
library(dplyr)

# Define paths
input_dir <- "DV_results"
output_dir <- "Enrichr_results"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
dbs <- "GO_Biological_Process_2023"
pval_cutoff <- 0.05
top_n <- 100

# Get unique cell types (assuming 'data' is defined elsewhere)
cell_types <- unique(data$CellType)
for (cell in cell_types) {
  print(paste("Working with", cell))
  file_path <- file.path(input_dir, paste0(gsub(" ", "_", cell), "_DV_results.xlsx"))

  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    next
  }

  wb <- loadWorkbook(file_path)
  dv_genes_up <- read.xlsx(wb, sheet = "Variability_Up")
  dv_genes_down <- read.xlsx(wb, sheet = "Variability_Down")
  
  # Subsample top_n genes for upregulated (Direction == 1)
  top_dv_genes_up <- dv_genes_up %>%head(top_n) %>% pull(genes)

  # Subsample top_n genes for downregulated (Direction == -1)
  top_dv_genes_down <- dv_genes_down %>% head(top_n) %>% pull(genes)

  # Run Enrichr for upregulated and downregulated genes
  enrichr_results_up <- run_enrichr(top_dv_genes_up, paste0("Increased variability in KO ", cell), dbs, pval_cutoff)
  enrichr_results_down <- run_enrichr(top_dv_genes_down, paste0("Decreased variability in KO ", cell), dbs, pval_cutoff)

  # Write results to Excel
  output_file <- file.path(output_dir, paste0(gsub("[^[:alnum:]_]", "", gsub(" ", "_", cell)), "_enrichr_DV_KO_var_results.xlsx"))
  wb_out <- createWorkbook()

  # Add upregulated results to the workbook if available
  if (nrow(enrichr_results_up) > 0) {
    addWorksheet(wb_out, "Upregulated_DV_Pathways")
    writeData(wb_out, "Upregulated_DV_Pathways", enrichr_results_up)
  } else {
    print(paste("No upregulated DV pathways for", cell))
  }

  # Add downregulated results to the workbook if available
  if (nrow(enrichr_results_down) > 0) {
    addWorksheet(wb_out, "Downregulated_DV_Pathways")
    writeData(wb_out, "Downregulated_DV_Pathways", enrichr_results_down)
  } else {
    print(paste("No downregulated DV pathways for", cell))
  }

  # Save the workbook after both sheets have been added
  saveWorkbook(wb_out, output_file, overwrite = TRUE)
  print(paste("Enrichment results saved in", output_file))
}









library(CellChat)
library(patchwork)
data_ko <- subset(data, Condition == 'Nr4a1 KO')
table(data_ko)
data_wt <- subset(data, Condition == 'WT')
table(data_ko)


# Assign sample labels explicitly and ensure they are factors
data_wt$samples <- factor("WT")
data_ko$samples <- factor("KO")

# Create CellChat objects
cellchat_WT <- createCellChat(object = data_wt, meta = data_wt@meta.data, group.by = "CellType", assay = "RNA")
cellchat_KO <- createCellChat(object = data_ko, meta = data_ko@meta.data, group.by = "CellType", assay = "RNA")

rm(data_ko, data_wt)

# Assign CellChat database
CellChatDB <- CellChatDB.mouse  # Use CellChatDB.human for human data
cellchat_WT@DB <- CellChatDB
cellchat_KO@DB <- CellChatDB

# Set CellType as identity
cellchat_WT <- setIdent(cellchat_WT, ident.use = "CellType")
cellchat_KO <- setIdent(cellchat_KO, ident.use = "CellType")

# Preprocess data
cellchat_WT <- subsetData(cellchat_WT)
cellchat_WT <- identifyOverExpressedGenes(cellchat_WT)
cellchat_WT <- identifyOverExpressedInteractions(cellchat_WT)
cellchat_WT <- computeCommunProb(cellchat_WT)

cellchat_KO <- subsetData(cellchat_KO)
cellchat_KO <- identifyOverExpressedGenes(cellchat_KO)
cellchat_KO <- identifyOverExpressedInteractions(cellchat_KO)
cellchat_KO <- computeCommunProb(cellchat_KO)

# Filter communication
cellchat_WT <- filterCommunication(cellchat_WT, min.cells = 10) # Start with a low value
cellchat_KO <- filterCommunication(cellchat_KO, min.cells = 10) # Start with a low value

# Compute at the pathway level
cellchat_WT <- computeCommunProbPathway(cellchat_WT)
cellchat_KO <- computeCommunProbPathway(cellchat_KO)

# Aggregate networks
cellchat_WT <- aggregateNet(cellchat_WT)
cellchat_KO <- aggregateNet(cellchat_KO)



# Individual plot show!
groupSize <- as.numeric(table(cellchat_WT@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_WT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_WT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_WT@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_WT, color.heatmap = "Reds")


groupSize <- as.numeric(table(cellchat_KO@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_KO@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_KO@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat_KO@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

par(mfrow=c(1,1))
netVisual_heatmap(cellchat_KO, color.heatmap = "Reds")


# Merge datasets for comparison
cellchat_merged <- mergeCellChat(list(WT = cellchat_WT, KO = cellchat_KO), add.names = c("KO", "WT"))


# Load CellChat objects
cellchat_merged <- readRDS("cellchat_merged.rds")
cellchat_WT <- readRDS("cellchat_WT.rds")
cellchat_KO <- readRDS("cellchat_KO.rds")


# Compare communication networks (individual conditions)
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE)
gg2 <- compareInteractions(cellchat_merged, measure = "weight")

# Corrected heatmaps for interaction comparison
gg3 <- netVisual_heatmap(cellchat_merged)  # Default "count"
gg4 <- compareInteractions(cellchat_merged, measure = "weight", show.legend = FALSE)
print(gg4)

# Custom visualization of differential interactions (Number of interactions in WT vs KO)
gg_diff <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "count")
print(gg_diff)


# Plot differential interactions (number of interactions)
gg3_diff <- netVisual_heatmap(cellchat_merged)
gg4_diff <- netVisual_heatmap(cellchat_merged, measure = "weight")

# Re-draw the heatmap with updated settings
draw(gg3_diff+gg4_diff)

png('plots/de_cell_chat_broad_cells.png', height = 6, width = 12, res = 300, units = 'in')
print(gg3_diff+gg4_diff)
dev.off()


# # Plot the joint interactions (both WT and KO)
# gg_joint <- compareInteractions(cellchat_merged, show.legend = FALSE, measure = "count")
# print(gg_joint)
#
# # Rank network changes (differentially expressed interactions)
# ggrank <- rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
# print(ggrank)
# png('plots/information_flow_cell_chat_broad_cells.png', height = 14, width = 6, res = 300, units = 'in')
# print(ggrank)
# dev.off()

# Signaling role scatter plot
gg0 <- netAnalysis_signalingRole_scatter(cellchat_merged)
print(gg0)
png('plots/signaling_role_cell_chat_broad_cells.png', height = 8, width = 8, res = 300, units = 'in')
print(gg0)
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
  thresh.fc = 0.05,  # Fold change threshold (adjust as needed)
  thresh.p = 0.05,   # p-value cutoff for significance
  group.DE.combined = FALSE  # Compare datasets separately
)

# Map the results of differential expression analysis onto the inferred cell-cell communications
net <- netMappingDEG(cellchat_merged, features.name = features.name, variable.all = TRUE)
write.csv(net, "net_mapping_DEG_results.csv", row.names = FALSE)

# Extract the ligand-receptor pairs with upregulated ligands in KO (positive condition)
net.up <- subsetCommunication(cellchat_merged, net = net, datasets = "KO", ligand.logFC = 0.5, receptor.logFC = NULL)

# Extract the ligand-receptor pairs with downregulated ligands in WT (negative condition)
net.down <- subsetCommunication(cellchat_merged, net = net, datasets = "WT", ligand.logFC = -0.5, receptor.logFC = NULL)

# Extract the genes associated with the upregulated and downregulated ligand-receptor pairs
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merged)


sources <- c(1:11)  # Sources you want to loop over
source_name <- levels(cellchat_merged@idents$KO)  # Source names

# Loop for both upregulated and downregulated signaling
for (source in sources){
  try({
    # Generate bubble plot for upregulated signaling
    pairLR.use.up <- net.up[, "interaction_name", drop = FALSE]
    gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = source, targets.use = c(1:11), comparison = c(1, 2),
                            remove.isolate = TRUE, color.heatmap = "Spectral",
                            title.name = paste0("Up-regulated signaling in KO (Source ", source_name[source], ")"))

    # Manually rotate x-axis labels for upregulated plot
    gg1 <- gg1 + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 12),
                       axis.text.y = element_text(size = 9))  # Adjust label size

    # Save the upregulated plot for each source
    plot_name_up <- paste0("plots/cellchat_img/upregulated_plot_source_", source_name[source], "_DEup_signaling.png")
    png(plot_name_up, height = 8, width = 10, res = 300, units = 'in')
    print(gg1)
    dev.off()
  }, silent = TRUE)  # If an error occurs in the upregulated plot, it will silently skip to the downregulated plot

  try({
    # Generate bubble plot for downregulated signaling
    pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]
    gg2 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.down, sources.use = source, targets.use = c(1:11), comparison = c(1, 2),
                            remove.isolate = TRUE, color.heatmap = "Spectral",
                            title.name = paste0("Down-regulated signaling in KO (Source ", source_name[source], ")"))

    # Manually rotate x-axis labels for downregulated plot
    gg2 <- gg2 + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 12),
                       axis.text.y = element_text(size = 9))  # Adjust label size

    # Save the downregulated plot for each source
    plot_name_down <- paste0("plots/cellchat_img/downregulated_plot_source_", source_name[source], "_DEdn_signaling.png")
    png(plot_name_down, height = 8, width = 10, res = 300, units = 'in')
    print(gg2)
    dev.off()
  }, silent = TRUE)  # If an error occurs in the downregulated plot, it will silently skip to the next source
}



# Save CellChat objects for future use
saveRDS(cellchat_merged, file = "cellchat_merged.rds")
saveRDS(cellchat_WT, file = "cellchat_WT.rds")
saveRDS(cellchat_KO, file = "cellchat_KO.rds")





# Annotate T cells
data <- readRDS(file.path(base_dir, "AllSamplesMerged_OldGenesRemoved_StringentFiltered_DoubletsRemoved_AmbientRNARemoved_AllGenesUMAP_FixedKOvsWT_Annotated_SampleSplit.rds"))

data_tcells <- process_and_extract_cell_type(data, "T cells", resolution = 1.0)

Idents(data_tcells) <- "CellType"

DimPlot(object = data_tcells, reduction = "umap", group.by = c("seurat_clusters"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1)

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

#> Mast cells <- Kit, Ms4a2, Cd3-
#> Innate Lymphoid cells <- Gata3, Cd3-
#> Gamma-delta T cells <- Trdc, Cd3, Cd4-, Cd8-
#> NK cells <- Ncr1, Klrb1b, Cd3-
#> SMCs <- Tpm2
DimPlot(object = data_tcells, reduction = "umap", group.by = c("seurat_clusters"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1,pt.size = 2)

# Find markers for cluster 9 ONLY
# cluster_markers <- FindMarkers(data_tcells,
#                                 ident.1 = 9, only.pos = TRUE,
#                                 group.by = "seurat_clusters",
#                                 logfc.threshold = 0.5, # Adjust as needed
#                                 min.pct = 0.1) # Adjust as needed
# top10_genes <- rownames(head(cluster_markers[order(cluster_markers$p_val_adj), ], 20))
# print(top10_genes)

data_tcells$new_celltypes = recode_factor(data_tcells$seurat_clusters,
                                           '0' = 'Innate Lymphoid',
                                           '1' = 'T cells',
                                           '2' = 'Cd4+ T cells',
                                           '3' = 'Cd4+ T cells',
                                           '4' = 'NK cells',
                                           '5' = 'Gamma-delta T cells',
                                           '6' = 'Cd8+ T cells',
                                           '7' = 'Cd4+ T cells',
                                           '8' = 'Innate Lymphoid',
                                           '9' = 'Cycling T cells',
                                           '10' = 'NK cells',
                                           '11' = 'SMCs',
                                           '12' = 'Mast cells',
                                           '13' = 'NK cells',
                                           '14' = 'Vascular endothelial',
                                          )

DimPlot(object = data_tcells, reduction = "umap", group.by = 'new_celltypes')

data_tcells$new_celltypes = factor(data_tcells$new_celltypes,
                                 levels = c( 'T cells','Cd4+ T cells',
                                             'Cd8+ T cells','NK cells',
                                             'Gamma-delta T cells',
                                             'Innate Lymphoid',
                                             'Mast cells','SMCs')
                                 )

DotPlot(data_tcells, features = c('Cd3e', 'Cd3g', 'Cd4', 'Cd8a','Cd8b1','Ncr1','Klrb1b','Trdc', 'Trgc1', 'Gata3','Ms4a2','Kit','Tpm2', 'Myh11'),
        group.by = 'new_celltypes', dot.min = 0.05, cols = 'RdBu' ) 


# Import cell types that do not belong to T cells
str_ntc <- c('SMCs', 'NK cells','Innate Lymphoid', 'Mast cells', 'Vascular endothelial')
idx <- data_tcells$new_celltypes %in% str_ntc
barcodes_ntc <- colnames(data_tcells)[idx]  # Works only if colnames(data_tcells) are barcodes
idx2 <- which(colnames(data) %in% barcodes_ntc)
data$CellType[idx2] <- as.character(data_tcells$new_celltypes[idx])


DimPlot(object = data, reduction = "umap", group.by = c("CellType"),
        label = TRUE, repel = TRUE, label.size = 3, label.box = TRUE, alpha = 1)
DimPlot(object = data, reduction = "umap", split.by = c("CellType"), group.by = c("CellType"),
        alpha = 1, ncol = 4)

# Create sub cell type annotation layer
data$SubCellType <- data$CellType  # Preserve broad annotation
common_cells_tcell <- intersect(rownames(data@meta.data), rownames(data_tcells@meta.data))
new_cell_types_subset <- as.character(data_tcells@meta.data[common_cells_tcell, "new_celltypes"])
data@meta.data[common_cells_tcell, "SubCellType"] <- new_cell_types_subset

DimPlot(object = data, reduction = "umap", group.by = c("SubCellType"),
        repel = TRUE, label.size = 3, alpha = 1)

saveRDS(data, file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated.rds"))





# Importing Colonocytes sub annotation
cell_ids <- read.csv("cell_ids_export.csv", header = FALSE)
colnames(data) <- cell_ids$V1

data_colo <- process_and_extract_cell_type(data, "Colonocytes", resolution = 1.5)

ctdf <- read.csv("shreyan_analysis/data/R/celltypes_colonocytes.csv")

# 1. Ensure barcodes are character vectors (important for matching)
colnames(data_colo) <- as.character(colnames(data_colo))  # Convert data_colo colnames to character
ctdf$X <- as.character(ctdf$X) # Convert ctdf$X to character

# 2. Match barcodes and get corresponding cell types (corrected)
cell_idx <- match(colnames(data_colo), ctdf$X)  # Correct order: data_colo -> ctdf$X
matched_celltypes <- ctdf$x[cell_idx]

# 3. Handle unmatched barcodes (essential)
unmatched_barcodes <- is.na(matched_celltypes)

if (any(unmatched_barcodes)) {
  warning(paste("Found", sum(unmatched_barcodes), "unmatched barcodes."))
  # b. Replace NA with "Unknown":
  matched_celltypes[unmatched_barcodes] <- "Unknown"
  # c. Keep NAs and handle them later.
}

# 4. Add the cell types to data_colo
data_colo$new_celltypes <- matched_celltypes

DimPlot(data_colo, reduction = "umap", group.by = 'new_celltypes',
        label = T, label.box = T)+ NoLegend()

data_colo$new_celltypes = factor(data_colo$new_celltypes,
                                   levels = c('Stem','TA','DCS','Imm Goblet',
                                              'Goblet','Aqp8+ Absorptive',
                                              'Car1+ Absorptive',
                                              'Tuft','Enterochromaffin',
                                              'Cck+ Enteroendocrine'))

## Canonical marker dotplots ####
p1 = DotPlot(data_colo, features = c('Lgr5','Smoc2','Mki67','Kcnq1','Cftr',
                                'Reg4','Ccl6','Spdef','Fam3d',
                                'Tff3','Sytl2','Aqp8','Ceacam20','Car1','Aqp4',
                                'Trpm5','Dclk1','Tph1','Chga',
                                'Cck'),dot.min = 0.15,
             dot.scale = 8, scale = T,
             group.by = 'new_celltypes')+
  scale_color_gradientn(colours = c('dodgerblue3','grey95','firebrick3'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p1
png('plots/Colonocytes_Dot.png', height = 4, width = 10, res = 300, units = 'in')
p1
dev.off()

# Create sub cell type annotation layer to origianl data
common_cells <- intersect(rownames(data@meta.data), rownames(data_colo@meta.data))
new_cell_types_subset <- as.character(data_colo@meta.data[common_cells, "new_celltypes"])
data@meta.data[common_cells, "SubCellType"] <- new_cell_types_subset

DimPlot(object = data, reduction = "umap", group.by = c("SubCellType"),
        alpha = 1)



# Colonocytes UMAP plot
data_colo$SubCellType = factor(data_colo$new_celltypes,
                                 levels = c('Stem','TA','DCS','Imm Goblet',
                                            'Goblet','Aqp8+ Absorptive',
                                            'Car1+ Absorptive',
                                            'Tuft','Enterochromaffin',
                                            'Cck+ Enteroendocrine'))

# UMAP for SubCellType
um_colo <- DimPlot(object = data_colo, reduction = "umap", group.by = "SubCellType", 
                     repel = TRUE, alpha = 1, pt.size = 0.1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

um_colo
# Save BatchID UMAP
png('plots/UMAP_colonocytes.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_colo)
dev.off()

# UMAP for SubCellType
um_colo <- DimPlot(object = data_colo, reduction = "umap", group.by = "Condition", 
                   repel = TRUE, alpha = 1, pt.size = 0.1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

um_colo
# Save BatchID UMAP
png('plots/UMAP_colonocytes_condition.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_colo)
dev.off()

#data_colo$Condition <- factor(ifelse(grepl("WT", data_colo$BatchID), "WT", "Nr4a1 KO"))
# Create DittoSeq plot
cc <- dittoBarPlot(data_colo, var = 'SubCellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("Colonocyte cell counts")  # Set title correctly

cc
# Save plot
png('plots/colonocytes_counts_subct.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()




saveRDS(data, file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated2.rds"))


data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated2.rds"))

##########
# Cell Scores for T cell exhaustion
##########
library(openxlsx)
library(AUCell)
library(ggsignif)

genelists = read.xlsx('Genesets_cell_scores.xlsx')

data_tcells <- process_and_extract_cell_type(data, "T cells", resolution = 1.0)
#data_cell_type <- subset(data, subset = CellType == cell_type_name)
Idents(data_tcells) <- "SubCellType"

DimPlot(object = data_tcells, reduction = "umap", group.by = c("SubCellType"),
        repel = TRUE, label.size = 3, alpha = 1, pt.size = 1)

# UMAP for SubCellType
um_tcells <- DimPlot(object = data_tcells, reduction = "umap", group.by = "SubCellType", 
                     repel = TRUE, alpha = 1, pt.size = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

um_tcells
# Save BatchID UMAP
png('plots/UMAP_tcells.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_tcells)
dev.off()

# UMAP for Condition
um_tcells <- DimPlot(object = data_tcells, reduction = "umap", group.by = "Condition", 
                     repel = TRUE, alpha = 1, pt.size = 1) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.text = element_blank(),  # Remove axis numbers
    axis.title = element_text(size = 14, color = "black")  # Keep axis titles
  ) +
  coord_cartesian(clip = "off")

um_tcells
# Save BatchID UMAP
png('plots/UMAP_tcells_condition.png', height = 6, width = 8, res = 300, units = 'in')
#svg('plots/UMAP_BatchID.svg', height = 6, width = 8)
print(um_tcells)
dev.off()


Exhaustion = list(Exhaustion = intersect(genelists$`Chen.et.al.(2023).T-cell.exhaustion`, rownames(data_tcells)))

data_tcells$score = AUCell_run(GetAssayData(data_tcells),Exhaustion,
                              normAUC = T)@assays@data@listData$AUC[1,]

p6 = data_tcells[[c('score','Condition','SubCellType')]] %>% group_by(SubCellType) %>%
  ggplot( aes(x = Condition, y = score, fill = Condition))+ theme_classic2() +
  geom_boxplot(color='black', outlier.size = 0.3) +
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion - Condition')+ylab('AUCell score')+
  facet_wrap(~SubCellType, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1), 
        plot.title = element_text(hjust = 0.5))
p6

table(data_tcells$Condition, data_tcells$SubCellType)

png('plots/Tcell_exhaustion_Condition_subct.png', height = 6, width = 7, res = 300, units = 'in')
p6
dev.off()

p6 = data_tcells[[c('score','Condition','CellType')]] %>% group_by(Condition) %>%
  ggplot( aes(x = Condition, y = score, fill = Condition))+ theme_classic2() +
  geom_boxplot(color='black', outlier.size = 0.3) +
  geom_signif(comparisons = list(c("WT","Nr4a1 KO")),
              map_signif_level = T, step_increase = 0.12, vjust = 2)+
  ggtitle('T cell exhaustion - Condition')+ylab('AUCell score')+
  facet_wrap(~CellType, scales = 'free')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1), 
        plot.title = element_text(hjust = 0.5))
p6

table(data_tcells$Condition, data_tcells$SubCellType)

png('plots/Tcell_exhaustion_Condition_ct.png', height = 6, width = 7, res = 300, units = 'in')
p6
dev.off()

# Create DittoSeq plot
cc <- dittoBarPlot(data_tcells, var = 'SubCellType', group.by = 'Condition', scale = 'count') +
  cust_theme + ggtitle("T cell counts")  # Set title correctly

cc
# Save plot
png('plots/tcell_counts_subct.png', height = 6, width = 8, res = 300, units = 'in')
print(cc)  # Use 'cc' instead of 'um_batch'
dev.off()

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
                        'Trdc', 
                        "Mki67")

# Colonocytes UMAP plot
data_tcells$SubCellType = factor(data_tcells$SubCellType,
                               levels = c('T cells', 'Cd4+ T cells', 'Cd8+ T cells',
                                          'Gamma-delta T cells', 'Cycling T cells'))

## Canonical marker dotplots ####
p1 = DotPlot(data_tcells, features = gene_marker_tcells, dot.min = 0.05,
             dot.scale = 8, scale = F,
             group.by = 'SubCellType')+
  scale_color_gradientn(colours = c('grey95','firebrick3'))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), # Increase x-axis text size
    axis.text.y = element_text(size = 15),  # Increase y-axis text size
    axis.title.x = element_text(size = 16), # Increase x-axis label size
    axis.title.y = element_text(size = 16)  # Increase y-axis label size
  )
p1
png('plots/Tcells_Dot.png', height = 4, width = 8, res = 300, units = 'in')
p1
dev.off()




library(ggplot2)
library(dplyr)
library(tidyr)
library(densityClust) # For density estimation


plot_wilcox_violin_tests_ct <- function(data, cell_type, gl, score_name = "Score", plot_type = "both") {
  
  # Subset metadata first to keep only relevant cell types
  cell_indices <- which(data$CellType == cell_type)
  if (length(cell_indices) == 0) {
    print(paste("No cells found for the specified cell type:", cell_type, ". Skipping."))
    return(ggplot() + ggtitle(paste("No cells found for", cell_type))) # Return empty plot
  }
  
  # Check if gl (gene) exists
  if (!all(gl %in% rownames(data))) {  # Check if ALL genes in gl are in the data
    missing_genes <- gl[!gl %in% rownames(data)]
    print(paste("Gene(s) not found in data:", paste(missing_genes, collapse = ", "), ". Skipping."))
    return(ggplot() + ggtitle(paste("Gene(s) not found:", paste(missing_genes, collapse = ", ")))) # Return empty plot
  }
  
  if (any(cell_indices > ncol(data))) { # Check if any index is too big
    print("cell_indices out of bounds!")
    # ... handle the error
  }
  
  # Subset assay data to relevant cells and genes
  assay_data <- GetAssayData(data, slot = 'data')[gl, cell_indices, drop = FALSE]
  
  # Convert to full matrix only for this subset
  if (inherits(assay_data, "dgCMatrix")) {
    assay_data_full <- as.matrix(assay_data)
  } else {
    assay_data_full <- assay_data
  }
  
  # Process gene expression values
  if (length(gl) == 1) {
    df <- as.data.frame(t(assay_data_full))
    colnames(df) <- paste0(score_name)
  } else {
    gene_data <- assay_data_full
    cell_scores <- rowMeans(gene_data, na.rm = TRUE)
    df <- data.frame(cell_scores)
    colnames(df) <- paste0(score_name)
  }
  
  # Clean up to free memory
  rm(assay_data_full)
  
  # Add metadata columns
  df$CellType <- factor(data$CellType[cell_indices])
  df$condition <- factor(data$Condition[cell_indices]) # Make condition a factor
  
  # Remove rows with NA values
  df <- df[!is.na(df[[score_name]]) & !is.na(df$CellType) & !is.na(df$condition), ]
  
  # Perform pairwise Wilcoxon test based on condition
  wilcox_test <- pairwise.wilcox.test(df[[score_name]], df$condition, p.adjust.method = "bonferroni")
  
  if (all(is.na(wilcox_test$p.value))) {
    print("No valid comparisons found in the Wilcoxon test.  Check that you have at least two conditions present")
    return(ggplot() + ggtitle("No valid comparisons found")) # Return an empty plot with a message
  }
  
  # Convert Wilcoxon results into a long format data frame
  comparison_df <- as.data.frame(wilcox_test$p.value)
  comparison_df_long <- comparison_df %>%
    rownames_to_column("Group1") %>%
    pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_value")
  
  # Apply p-value labels
  p_value_labels <- comparison_df_long %>%
    mutate(label = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      P_value < 0.1   ~ ".",
      TRUE            ~ "ns"
    ))
  
  # *** REMOVE THE FILTERING OF SIGNIFICANT COMPARISONS ***
  comparisons_list <- p_value_labels %>% # No more filtering!
    select(Group1, Group2) %>%
    rowwise() %>%
    mutate(pairwise = list(c(as.character(Group1), as.character(Group2)))) %>%
    pull(pairwise)
  
  # Calculate max value for proper geom_signif placement
  y_max <- max(df[[score_name]], na.rm = TRUE)
  
  # Adjust significance bar height based on plot type
  y_offset <- 0.05 * y_max
  
  # Create the plot
  plot <- ggplot(df, aes(x = condition, y = .data[[score_name]], fill = condition)) # Condition on x-axis
  
  # Add plots based on plot_type argument
  if (plot_type == "violin") {
    plot <- plot + geom_violin(alpha = 0.6, width = 1)
  } else if (plot_type == "boxplot") {
    plot <- plot + geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.6, color = "black")
  } else if (plot_type == "both") {
    plot <- plot +
      geom_violin(alpha = 0.6, width = 1) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4, color = "black")
  }
  
  # Add significance annotations at adjusted y position
  plot <- plot +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    stat_summary(fun=mean, colour="black", geom="point", 
                 shape=18, size=3, show.legend=FALSE) + 
    ggtitle(paste(score_name, "for", cell_type)) 
  
  plot <- plot + geom_signif(
    comparisons = comparisons_list,
    map_signif_level = TRUE,
    step_increase = y_offset / y_max,
    y_position = y_max + y_offset
  )
  return(plot)
}


library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl) # For reading Excel files (pd.read is Python)

# Create directory if it doesn't exist
store <- "plots/immune_genes"
if (!dir.exists(store)) {
  dir.create(store, recursive = TRUE) # recursive = TRUE creates parent directories if needed
}

# Read the Excel file.  Assumes the first column is the gene names.  Adjust as necessary.
genes_df <- read_excel("Immune_cell_markers_NR4A1_project.xlsx")
genes <- genes_df[[2]] # Extract the gene names from the first column.

cell_types <- c("Macrophages", "B cells", "Plasma B cells", "T cells")
# Loop over all cell types and genes
for (cell_type in cell_types) {
  for (gene in genes) {
    plot_wt <- plot_wilcox_violin_tests_ct(data, cell_type = cell_type, gl = gene,
                                           score_name = paste0(gene, " expression"),
                                           plot_type = "both")
    
    # Construct the filename carefully to avoid issues
    filename <- paste0(store, "/", cell_type, "_", gene, ".png") # Include gene name in file
    
    png(filename, height = 6, width = 6, res = 300, units = 'in')
    print(plot_wt)
    dev.off()
    
    # Print a message to track progress (optional but helpful)
    print(paste("Plot saved:", filename))
  }
}




# SGS results
HVG_res = HVG_splinefit(GetAssayData(pbmc,layer = 'counts'), diptest=TRUE,
                        dropout.filter=TRUE, nHVGs=500)

# HVG results
head(HVG_res[HVG_res$HVG == TRUE,])

# Highly Variable Genes that are known TFs 
HVG_res[(HVG_res$HVG == TRUE & HVG_res$TF == TRUE),][1:5,]

SGS_res = scSGS(GetAssayData(pbmc, layer='counts'), 'S100A9', nHVG=500,
                rm.mt=T, rm.rp=T,
                calcHVG=TRUE, verbose=TRUE)

# scSGS DE results
head(SGS_res$DE)

