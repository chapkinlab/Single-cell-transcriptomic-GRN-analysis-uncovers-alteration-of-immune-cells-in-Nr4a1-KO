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
library(enrichR)


base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
#base_dir <- "/mnt/SCDC/Optimus/selim_working_dir/2023_nr4a1_colon/results"
setwd(base_dir)
data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated4.rds"))


# Ensure `cell_types` is properly extracted
cell_types <- unique(data$CellType)
# Set the grouping variable for Seurat
Idents(data) <- "Condition"
# Define output directory
output_dir <- "DE_results_genes"
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
  addWorksheet(wb, "All_Pval_adj_0.05")
  addWorksheet(wb, "Up_regulated")
  addWorksheet(wb, "Down_regulated")
  
  # Write data to respective sheets
  writeData(wb, sheet = "All_Pval_adj_0.05", all_markers)
  writeData(wb, sheet = "Up_regulated", upregulated)
  writeData(wb, sheet = "Down_regulated", downregulated)
  
  # Ensure consistent file naming
  save_cell_name <- gsub(" ", "_", cell)  # Replace spaces with underscores *once*
  output_file <- file.path(output_dir, paste0(save_cell_name, "_markers_Wilcoxon.xlsx"))
  
  # Save workbook
  saveWorkbook(wb, file = output_file, overwrite = TRUE)
  
  cat("Saved:", output_file, "\n")
}






# Define paths
input_dir <- "DE_results_genes"
output_dir <- "Enrichr_results"

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set Enrichr database and parameters
setEnrichrSite('Enrichr')
dbs <- "GO_Biological_Process_2025"
pval_cutoff <- 0.05
top_n <- 500

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




library(data.table)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(Seurat) # Ensure Seurat is loaded for your data object and FindMarkers
library(openxlsx) # For your DE results saving

# --- Your initial setup code ---
base_dir <- "Z:\\selim_working_dir\\2023_nr4a1_colon\\results"
setwd(base_dir)
#data <- readRDS(file.path(base_dir, "SamplesSplit_AmbientRNARemoved_FixedKOvsWT_SubAnnotated3.rds"))

DefaultAssay(data) <- "RNA"
cell_types <- unique(data$CellType)
Idents(data) <- "Condition"

# Define primary output directories
output_dir_de <- "DE_results_all"
output_dir_gsea <- "GSEA_results" # Main GSEA results directory

# Create main output directories if they don't exist
if (!dir.exists(output_dir_de)) {
  dir.create(output_dir_de, recursive = TRUE)
}
if (!dir.exists(output_dir_gsea)) {
  dir.create(output_dir_gsea, recursive = TRUE)
}

# --- Load MSigDB gene sets ---
# GOBP = C5
msigdbr_df <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
pathways_gobp <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# --- Your Loop for Cell Types ---
for (cell in cell_types) {
  cat("Processing cell type:", cell, "\n")
  
  # Sanitize cell name for directory and file names
  save_cell_name <- gsub(" ", "_", cell) # Replace spaces with underscores
  # You might want to remove other problematic characters too, e.g., using gsub("[^a-zA-Z0-9_]", "", cell)
  
  # Define and create the cell-specific subdirectory for GSEA results
  cell_gsea_subdir <- file.path(output_dir_gsea, save_cell_name)
  if (!dir.exists(cell_gsea_subdir)) {
    dir.create(cell_gsea_subdir, recursive = TRUE)
  }
  
  cell_subset <- subset(data, CellType == cell)
  
  markers <- FindMarkers(cell_subset,
                         ident.1 = "Nr4a1 KO",
                         ident.2 = "WT",
                         test.use = "wilcox",
                         logfc.threshold = 0.1,
                         min.pct = 0.01,
                         min.diff.pct = -Inf,
                         verbose = TRUE)

  # --- Differential Expression (DE) results saving (your original code) ---
  # This part also benefits from subdirectories
  cell_de_subdir <- file.path(output_dir_de, save_cell_name)
  if (!dir.exists(cell_de_subdir)) {
    dir.create(cell_de_subdir, recursive = TRUE)
  }
  
  all_markers <- markers[markers$p_val_adj < 0.05, ]
  all_markers <- all_markers[order(abs(all_markers$avg_log2FC), decreasing = TRUE), ]
  
  upregulated <- all_markers[all_markers$avg_log2FC > 0, ]
  downregulated <- all_markers[all_markers$avg_log2FC < 0, ]
  
  all_markers$Gene <- rownames(all_markers)
  upregulated$Gene <- rownames(upregulated)
  downregulated$Gene <- rownames(downregulated)
  
  wb <- createWorkbook()
  addWorksheet(wb, "All_Pval_adj_0.05")
  addWorksheet(wb, "Up_regulated")
  addWorksheet(wb, "Down_regulated")
  writeData(wb, sheet = "All_Pval_adj_0.05", all_markers)
  writeData(wb, sheet = "Up_regulated", upregulated)
  writeData(wb, sheet = "Down_regulated", downregulated)
  
  # Save DE workbook into cell-specific DE subdirectory
  output_file_de <- file.path(cell_de_subdir, paste0(save_cell_name, "_markers_Wilcoxon.xlsx"))
  saveWorkbook(wb, file = output_file_de, overwrite = TRUE)
  cat("Saved DE results:", output_file_de, "\n")
  
  # --- fgsea analysis starts here ---
  
  ranks <- markers$avg_log2FC
  names(ranks) <- rownames(markers)
  
  ranks <- ranks[!is.na(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  
  ranks <- sort(ranks, decreasing = TRUE)
  
  fgseaRes <- fgsea(pathways = pathways_gobp,
                    stats = ranks,
                    minSize = 10,
                    maxSize = 1000,
                    nPermSimple = 10000)
  
  fgseaRes$padj <- p.adjust(fgseaRes$pval, method = "BH")
  fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
  
  # Save fgsea results into the cell-specific GSEA subdirectory
  output_file_gsea <- file.path(cell_gsea_subdir, paste0(save_cell_name, "_fgsea_GOBP.tsv"))
  fwrite(fgseaRes, file = output_file_gsea, sep = "\t", row.names = FALSE)
  cat("Saved GSEA results:", output_file_gsea, "\n")
  
  # Optional: Generate and save enrichment plots for top pathways
  top_pathways_up <- fgseaRes[NES > 0 & padj < 0.05][head(order(padj), n=5), pathway]
  top_pathways_down <- fgseaRes[NES < 0 & padj < 0.05][head(order(padj), n=5), pathway]
  
  selected_pathways <- unique(c(top_pathways_up, top_pathways_down))
  
  if (length(selected_pathways) > 0) {
    for (pathway_name in selected_pathways) {
      gsea_plot <- plotEnrichment(pathways_gobp[[pathway_name]], ranks) +
        labs(title = pathway_name) +
        theme_minimal()
      
      # Save GSEA plots into the cell-specific GSEA subdirectory
      plot_file <- file.path(cell_gsea_subdir, paste0(save_cell_name, "_GSEA_plot_", gsub(" ", "_", pathway_name), ".png"))
      ggsave(plot_file, plot = gsea_plot, width = 8, height = 6, dpi = 300)
      cat("Saved GSEA plot:", plot_file, "\n")
    }
  } else {
    cat("No significant GOBP pathways found for plotting in", cell, "\n")
  }
}

